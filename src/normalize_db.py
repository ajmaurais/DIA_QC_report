
import sys
import argparse
import sqlite3
import os
from datetime import datetime
from multiprocessing import cpu_count

import numpy as np
import pandas as pd
import directlfq.utils as lfqutils
from directlfq.normalization import NormalizationManagerSamplesOnSelectedProteins as dlfq_norm
import directlfq.protein_intensity_estimation as lfqprot_estimation

from .submodules.dia_db_utils import METADATA_TIME_FORMAT
from .submodules.dia_db_utils import update_meta_value
from .submodules.dia_db_utils import check_schema_version
from .submodules.dia_db_utils import update_acquired_ranks
from .submodules.dia_db_utils import mark_all_reps_includced, mark_reps_skipped
from .submodules.logger import LOGGER


def median_normalize(df, key_cols, value_col):
    '''
    Median normalize peak areas in long formatted datafraeme.

    Parameters
    ----------
    df: pd.DataFrame
        Long formatted data frame
    key_cols: list
        The names of column(s) which uniquely identify each replicate.
    value_col: str
        The name of the column with peak areas.
        This function will log2 transform this column so the areas should be in linear space.
    '''
    cols = df.columns.tolist()
    cols.append('normalizedArea')

    # log2 transform
    df['log2Area'] = np.log2(df[value_col].replace(0, np.nan))

    # determine value to shift log2Area for each replicate
    medians = df.groupby(key_cols)['log2Area'].median()
    median_shifts = medians - np.mean(medians)
    median_shifts = median_shifts.reset_index()
    median_shifts = {row.replicateId: row.log2Area for row in median_shifts.itertuples()}

    # shift log2Areas to normalize medians
    df['log2NormalizedArea'] = df['log2Area'] - df['replicateId'].apply(lambda x: median_shifts[x])

    # convert back to linear space
    df['normalizedArea'] = np.power(2, df['log2NormalizedArea'])
    df['normalizedArea'] = df['normalizedArea'].replace(np.nan, 0)

    return df[cols]


def main():

    parser = argparse.ArgumentParser(description='Perform DirectLFQ or median normalization on batch database.')
    # parser.add_argument('-m', '--method', choices=['DirectLFQ', 'median'], default='DirectLFQ',
    #                     help='Normalization method to use.')
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

    # get precursor table from db
    LOGGER.info('Reading precursors from database...')
    df = pd.read_sql('''
        SELECT
            p.replicateId,
            p.peptideId,
            ptp.proteinId as protein,
            p.modifiedSequence,
            p.precursorCharge,
            p.totalAreaFragment
        FROM precursors p
        LEFT JOIN peptideToProtein ptp ON ptp.peptideId == p.peptideId
        LEFT JOIN replicates r ON p.replicateId == r.replicateId
        WHERE r.includeRep == TRUE; ''',
        conn)
    LOGGER.info('Finished reading precursors.')

    if len(df.index) == 0:
        LOGGER.error('All replicates in database have been excluded!')
        sys.exit(1)

    df['ion'] = df['modifiedSequence'] + '_' + df['precursorCharge'].astype(str)
    precursor_ids = df[['ion', 'peptideId', 'modifiedSequence', 'precursorCharge']].drop_duplicates()
    precursor_ids = {x.ion: (x.peptideId,
                             x.modifiedSequence,
                             x.precursorCharge) for x in precursor_ids.itertuples()}

    REP_COLUMN_NAME = 'replicateId'

    # add 1 to all precursor areas so the log2 of zeros is finite
    df['totalAreaFragment'] = df['totalAreaFragment'] + 1
    df = df.reset_index()

    # pivot wider
    input_df = df.pivot(index=['protein', 'ion'], columns=REP_COLUMN_NAME, values='totalAreaFragment')

    # remove missing rows
    n_not_missing = len(input_df.index)
    input_df = input_df[~input_df.apply(lambda x: any(x.isnull().values), axis = 1)]
    keep_precursors = set([x[1] for x in input_df.index.tolist()])
    if len(input_df.index) < n_not_missing:
        LOGGER.warning(f'Removed {n_not_missing - len(input_df.index):,} of {n_not_missing:,} precursors with missing values.')
        LOGGER.warning(f'{len(input_df.index):,} precursors without missing values remain.')

    # format input_df for DirectLFQ normalization
    input_df.columns.name = None
    input_df = input_df.reset_index()
    input_df = lfqutils.index_and_log_transform_input_df(input_df)
    input_df = lfqutils.remove_allnan_rows_input_df(input_df)
    input_df = dlfq_norm(input_df, num_samples_quadratic=10).complete_dataframe

    # Do DirectLFQ normalization
    LOGGER.info('Running DirectLFQ normalization.')
    protein_df, _ = lfqprot_estimation.estimate_protein_intensities(input_df,
                                                                    min_nonan=1,
                                                                    num_samples_quadratic=10,
                                                                    num_cores=cpu_count())
                                                                    # num_cores = 1)
    LOGGER.info('Finished with DirectLFQ normalization.')

    # pivot protein_df longer
    protein_df = protein_df.melt(id_vars='protein', value_name='normalizedArea', var_name='replicateId')


    # remove precursors with missing values
    df = df[df['ion'].apply(lambda x: x in keep_precursors)]

    # median normalize precursor quants
    ion_df = median_normalize(df[[REP_COLUMN_NAME, 'ion', 'totalAreaFragment']].drop_duplicates(),
                              [REP_COLUMN_NAME], 'totalAreaFragment')
    ion_df['peptideId'] = ion_df['ion'].apply(lambda x: precursor_ids[x][0])
    ion_df['modifiedSequence'] = ion_df['ion'].apply(lambda x: precursor_ids[x][1])
    ion_df['precursorCharge'] = ion_df['ion'].apply(lambda x: precursor_ids[x][2])

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
    ion_data = [(row.normalizedArea,
                 getattr(row, REP_COLUMN_NAME),
                 row.peptideId,
                 row.precursorCharge) for row in ion_df.itertuples()]
    cur.executemany('''
                    UPDATE precursors
                        SET normalizedArea = ?
                    WHERE replicateId = ? AND
                          peptideId = ? AND
                          precursorCharge = ? ;''',
                    ion_data)
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
    for row in protein_df.itertuples():
        cur.execute(protein_query, (row.replicateId, row.protein, row.normalizedArea,
                                    row.normalizedArea))
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
                'precursor_normalization_method': 'median',
                'protein_normalization_method': 'DirectLFQ',
                'is_normalized': 'True',
                'Normalization command': current_command,
                'command_log': previous_commands + current_command}
    for key, value in metadata.items():
        conn = update_meta_value(conn, key, value)

    conn.close()


if __name__ == '__main__':
    main()

