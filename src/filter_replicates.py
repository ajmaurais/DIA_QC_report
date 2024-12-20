
import sys
import argparse
import sqlite3
import os
from datetime import datetime

from .submodules.dia_db_utils import update_meta_value
from .submodules.dia_db_utils import check_schema_version
from .submodules.dia_db_utils import mark_all_reps_included, mark_reps_skipped
from .submodules.dia_db_utils import update_command_log
from .submodules.dia_db_utils import is_normalized, is_imputed
from .submodules.logger import LOGGER
from . import __version__ as PROGRAM_VERSION

COMMAND_DESCRIPTION = 'Filter replicates in database.'
CONFLICTING_ARG_MESSAGE = 'exclude Rep/Project and --include_all arguments conflict!'

def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION + '''
Replicates or projects to exclude from imputation and normalization.
The replicates.includeRep value will simply be set to FALSE the replicate 
will not be deleted from the database.''')

    parser.add_argument('-x', '--excludeRep', action='append', default=[],
                        help='Add replicate to exclude.')
    parser.add_argument('-p', '--excludeProject', action='append', default=[],
                        help='Exclude project.')
    parser.add_argument('-a', '--includeAll', action='store_true', default=False,
                        dest='include_all', help='Include all replicates. All '
                        'replicates.includeRep values are set to TRUE.')
    parser.add_argument('db', help='Path to sqlite batch/qc database.')

    return parser.parse_args(argv)


def _main(args):
    '''
    Actual main method. `args` Should be an initialized argparse namespace.
    '''

    exclude_reps = sum([len(args.excludeRep), len(args.excludeProject)]) > 0
    if exclude_reps and args.include_all:
        LOGGER.error(CONFLICTING_ARG_MESSAGE)
        sys.exit(1)

    if os.path.isfile(args.db):
        conn = sqlite3.connect(args.db)
    else:
        LOGGER.error(f'Database file ({args.db}) does not exist!')
        sys.exit(1)

    # check database version
    if not check_schema_version(conn):
        sys.exit(1)

    # remove replicates if applicable
    if exclude_reps:
        if not mark_reps_skipped(conn, reps=args.excludeRep,
                                 projects=args.excludeProject):
            sys.exit(1)

    if args.include_all:
        mark_all_reps_included(conn)

    # unset imputed values if applicable
    # if is_imputed(conn):
    #     cur = conn.cursor()
    #     LOGGER.info('Setting existing precursor normalizedArea values to NULL.')
    #     cur.execute('UPDATE precursors SET totalAreaFragment = NULL WHERE isImputed == 1')
    #     cur.execute('UPDATE precursors SET isImputed = 0')
    #     conn.commit()

    # unset normalized values if applicable
    if is_normalized(conn):
        cur = conn.cursor()
        cur.execute('UPDATE precursors SET normalizedArea = NULL')
        cur.execute('UPDATE proteinQuants SET normalizedAbundance = NULL')
        conn.commit()

    # update commandLog
    update_command_log(conn, sys.argv, os.getcwd())

    # Update normalization method in metadata
    LOGGER.info('Updating metadata...')
    metadata = {'is_normalized': 'False'}
                # 'is_imputed': 'False'} # set to false if we are changing replicates
    for key, value in metadata.items():
        conn = update_meta_value(conn, key, value)

    conn.close()


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc normalize" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()

