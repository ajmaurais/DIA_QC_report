
import sys
import argparse
import sqlite3
import os
from datetime import datetime

from .submodules.normalization import MedianNormalizer, DirectlfqNormalizer
from .submodules.dia_db_utils import METADATA_TIME_FORMAT
from .submodules.dia_db_utils import update_meta_value
from .submodules.dia_db_utils import check_schema_version
from .submodules.dia_db_utils import mark_all_reps_includced, mark_reps_skipped
from .submodules.logger import LOGGER


def main():
    parser = argparse.ArgumentParser(description='Perform DirectLFQ or median normalization on batch database.')
    parser.add_argument('-m', '--method', choices=['DirectLFQ', 'median'], default='median',
                        help='Normalization method to use. ')
    exclude_args = parser.add_argument_group('Filter replicates',
                                             'Add replicates or projects to exclude from normalization. '
                                             'The replicates.includeRep value will simply be set to FALSE '
                                             'the replicate will not be deleted from the database.')
    exclude_args.add_argument('-x', '--excludeRep', action='append', default=[],
                              help='Add replicate to exclude from normalization.')
    exclude_args.add_argument('-p', '--excludeProject', action='append', default=[],
                              help='Exclude project from normalization.')
    exclude_args.add_argument('-a', '--useAll', action='store_true', default=False,
                              help='Use all replicates for normalization and set all '
                                    'replicates.includeRep values to TRUE.')
    parser.add_argument('db', help='Path to sqlite batch database.')

    args = parser.parse_args()

    exclude_reps = sum([len(args.excludeRep), len(args.excludeProject)]) > 0
    if exclude_reps and args.useAll:
        LOGGER.error('exclude Rep/Project and --useAll arguments conflict!')
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

    if args.useAll:
        mark_all_reps_includced(conn)

    # init normalization manager
    if args.method == 'median':
        norm_manager = MedianNormalizer(conn, keep_na=False)
        precursor_normalization_method = 'median'
        protein_normalization_method = 'median'
    elif args.method == 'DirectLFQ':
        norm_manager = DirectlfqNormalizer(conn)
        precursor_normalization_method = 'median'
        protein_normalization_method = 'DirectLFQ'
    else:
        LOGGER.error(f'"{args.method}" is an unknown normalization method!')
        sys.exit(1)

    norm_manager.normalize()

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
                    WHERE replicateId = ? AND
                          peptideId = ? AND
                          precursorCharge = ? ;''',
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
    metadata = {'Normalization time': datetime.now().strftime(METADATA_TIME_FORMAT),
                'precursor_normalization_method': precursor_normalization_method,
                'protein_normalization_method': protein_normalization_method,
                'is_normalized': 'True',
                'Normalization command': current_command,
                'command_log': previous_commands + current_command}
    for key, value in metadata.items():
        conn = update_meta_value(conn, key, value)

    conn.close()


if __name__ == '__main__':
    main()

