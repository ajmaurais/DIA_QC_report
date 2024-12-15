
import sys
import argparse
import sqlite3
import os
from datetime import datetime

from .submodules.dia_db_utils import update_meta_value
from .submodules.dia_db_utils import check_schema_version
from .submodules.dia_db_utils import mark_all_reps_includced, mark_reps_skipped
from .submodules.logger import LOGGER
from . import __version__ as PROGRAM_VERSION

COMMAND_DESCRIPTION = 'Filter replicates in database.'

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
    Actual main method. `args` Should be initialized argparse namespace.
    '''

    exclude_reps = sum([len(args.excludeRep), len(args.excludeProject)]) > 0
    if exclude_reps and args.includeAll:
        LOGGER.error('exclude Rep/Project and --includeAll arguments conflict!')
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
        mark_all_reps_includced(conn)

    # get commands previously run on db
    cur = conn.cursor()
    cur.execute('SELECT value FROM metadata WHERE key == "command_log"')
    previous_commands = cur.fetchall()
    if len(previous_commands) == 0:
        LOGGER.warning('Missing command_log metadata entry!')
        previous_commands = 'MISSING_COMMAND_LOG\n'
    else:
        previous_commands = previous_commands[0][0] + '\n'
    current_command = ' '.join(sys.argv)

    # Update normalization method in metadata
    LOGGER.info('Updating metadata...')
    metadata = {'is_normalized': 'False',
                'command_log': previous_commands + current_command}
    for key, value in metadata.items():
        conn = update_meta_value(conn, key, value)

    conn.close()


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc normalize" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()

