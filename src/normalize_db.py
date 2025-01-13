
import sys
import argparse
import sqlite3
import os
from datetime import datetime

from .submodules.normalization import MedianNormalizer, DirectlfqNormalizer
from .submodules.normalization import DirectlfqNormalizer, NORMALIZATION_METHODS
from .submodules.dia_db_utils import METADATA_TIME_FORMAT
from .submodules.dia_db_utils import update_meta_value
from .submodules.dia_db_utils import check_schema_version
from .submodules.dia_db_utils import IS_NORMALIZED, PRECURSOR_NORM_METHOD, PROTEIN_NORM_METHOD
from .submodules.dia_db_utils import update_command_log
from .submodules.logger import LOGGER
from . import __version__ as PROGRAM_VERSION

COMMAND_DESCRIPTION = 'Perform DirectLFQ or median normalization on batch database.'


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)

    norm_settings = parser.add_argument_group('Normalization settings')
    norm_settings.add_argument('-m', '--method', choices=NORMALIZATION_METHODS, default='DirectLFQ',
                               help='Normalization method to use. Default is "DirectLFQ"')
    norm_settings.add_argument('--keepMissing', default=False, action='store_true', dest='keep_missing',
                               help="Don't exclude precursors which are missing in 1 or more "
                                    "replicates from normalization. This option is only compatible "
                                    "with median normalization.")

    parser.add_argument('db', help='Path to sqlite batch/qc database.')

    return parser.parse_args(argv)


def _main(args):
    '''
    Actual main method. `args` Should be initialized argparse namespace.
    '''

    if os.path.isfile(args.db):
        conn = sqlite3.connect(args.db)
    else:
        LOGGER.error(f'Database file ({args.db}) does not exist!')
        sys.exit(1)

    # check database version
    if not check_schema_version(conn):
        sys.exit(1)

    # init normalization manager
    if args.method == 'median':
        norm_manager = MedianNormalizer(conn, keep_na=args.keep_missing)
        precursor_normalization_method = 'median'
        protein_normalization_method = 'median'
    elif args.method == 'DirectLFQ':
        if args.keep_missing:
            LOGGER.error('--keepMissing option not compatible with DirectLFQ')
            sys.exit(1)
        norm_manager = DirectlfqNormalizer(conn)
        precursor_normalization_method = 'median'
        protein_normalization_method = 'DirectLFQ'
    else:
        LOGGER.error(f'"{args.method}" is an unknown normalization method!')
        sys.exit(1)

    if not norm_manager.normalize():
        sys.exit(1)

    # delete existing normalized values
    cur = conn.cursor()
    LOGGER.info('Setting existing precursor normalizedArea values to NULL.')
    cur.execute('UPDATE precursors SET normalizedArea = NULL')
    LOGGER.info('Setting existing protein normalizedAbundance values to NULL.')
    cur.execute('UPDATE proteinQuants SET normalizedAbundance = NULL')
    conn.commit()

    # update precursor rows
    LOGGER.info('Updating precursor normalizedArea values...')
    cur = conn.cursor()
    prec_data = [(row.normalizedArea,
                  row.replicateId,
                  row.peptideId,
                  row.precursorCharge) for row in norm_manager.precursors.itertuples()]
    cur.executemany('''
                    UPDATE precursors
                        SET normalizedArea = ?
                    WHERE replicateId == ? AND
                          peptideId == ? AND
                          precursorCharge == ? ;''',
                    prec_data)
    conn.commit()
    LOGGER.info('Done updating precursor normalizedArea values.')

    # insert or update proteinQuant rows
    protein_query = '''
        INSERT INTO proteinQuants (replicateId, proteinId, normalizedAbundance)
        VALUES (?, ?, ?)
        ON CONFLICT(replicateId, proteinId) DO UPDATE SET normalizedAbundance = ?
        '''
    LOGGER.info('Updating protein normalizedAbundance values...')
    cur = conn.cursor()
    for row in norm_manager.proteins.itertuples():
        cur.execute(protein_query, (row.replicateId, row.proteinId, row.normalizedAbundance,
                                    row.normalizedAbundance))
    conn.commit()
    LOGGER.info('Done updating protein normalizedAbundance values.')

    # update commandLog
    update_command_log(conn, sys.argv, os.getcwd())

    # Update normalization method in metadata
    LOGGER.info('Updating metadata...')
    metadata = {PRECURSOR_NORM_METHOD: precursor_normalization_method,
                PROTEIN_NORM_METHOD: protein_normalization_method,
                IS_NORMALIZED: 'True'}
    for key, value in metadata.items():
        conn = update_meta_value(conn, key, value)

    conn.close()


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc normalize" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()

