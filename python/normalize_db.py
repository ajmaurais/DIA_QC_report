
import sys
import argparse
import sqlite3
import os
import logging
from datetime import datetime
import re
from multiprocessing import cpu_count

import numpy as np
import pandas as pd
import directlfq.utils as lfqutils
from directlfq.normalization import NormalizationManagerSamplesOnSelectedProteins as dlfq_norm
import directlfq.protein_intensity_estimation as lfqprot_estimation

logging.basicConfig(
    level=logging.INFO, format='%(asctime)s - %(filename)s %(funcName)s - %(levelname)s: %(message)s'
)
LOGGER = logging.getLogger()

TIME_FORMAT = '%m/%d/%Y %H:%M:%S'


def update_meta_value(conn, key, value):
    '''
    Add or update value in metadata table.
    
    Parameters
    ----------
    conn: sqlite3.Connection:
        Database connection.
    key: str
        The metadata key
    value: str
        The metadata value
    '''
    cur = conn.cursor()
    cur.execute('''
        INSERT INTO metadata
            (key, value) VALUES (?, ?)
        ON CONFLICT(key) DO UPDATE SET value = ? ''',
                (key, value, value))
    conn.commit()

    return conn


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

    parser = argparse.ArgumentParser(description='Perform DirectLFQ or median normalization on precursor DB.')
    # parser.add_argument('-m', '--method', choices=['DirectLFQ', 'median'], default='DirectLFQ',
    #                     help='Normalization method to use.')
    parser.add_argument('db')

    args = parser.parse_args()

    conn = sqlite3.connect(args.db)

    # get precursor table from db
    LOGGER.info('Reading precursors from database...')
    df = pd.read_sql('''
        SELECT
            p.replicateId,
            ptp.proteinId as protein,
            p.modifiedSequence,
            p.precursorCharge,
            p.totalAreaFragment
        FROM precursors p
        LEFT JOIN peptideToProtein ptp ON ptp.modifiedSequence == p.modifiedSequence; ''',
        conn)

    df['ion'] = df['modifiedSequence'] + '_' + df['precursorCharge'].astype(str)
    precursor_ids = df[['ion', 'modifiedSequence', 'precursorCharge']].drop_duplicates()
    precursor_ids = {x.ion: (x.modifiedSequence, x.precursorCharge) for x in precursor_ids.itertuples()}

    REP_COLUMN_NAME = 'replicateId'

    # pivot wider
    input_df = df.pivot(index=['protein', 'ion'], columns=REP_COLUMN_NAME, values='totalAreaFragment')

    # format input_df for DirectLFQ normalization
    input_df.columns.name = None
    input_df = input_df.reset_index()
    input_df = lfqutils.index_and_log_transform_input_df(input_df)
    input_df = lfqutils.remove_allnan_rows_input_df(input_df)
    input_df = dlfq_norm(input_df, num_samples_quadratic=10).complete_dataframe

    # input_df.to_csv('/home/ajm/code/DIA_QC_report/testData/normalize_db/input_df.tsv', sep='\t', index=True)
     
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
    
    # protein_df.to_csv('/home/ajm/code/DIA_QC_report/testData/normalize_db/protein_df.tsv', sep='\t', index=False)
    
    # median normalize precursor quants
    ion_df = median_normalize(df[[REP_COLUMN_NAME, 'ion', 'totalAreaFragment']].drop_duplicates(),
                              [REP_COLUMN_NAME], 'totalAreaFragment')
    ion_df['modifiedSequence'] = ion_df['ion'].apply(lambda x: precursor_ids[x][0])
    ion_df['precursorCharge'] = ion_df['ion'].apply(lambda x: precursor_ids[x][1])

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
                 row.modifiedSequence,
                 row.precursorCharge) for row in ion_df.itertuples()]
    cur.executemany('''
                    UPDATE precursors
                        SET normalizedArea = ?
                    WHERE replicateId = ? AND modifiedSequence = ? AND precursorCharge = ? ''',
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
        cur.execute(protein_query, (row.replicateId, row.protein, row.normalizedArea, row.normalizedArea))
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

    # update normalization method in metadata
    LOGGER.info('Updating metadata...')
    conn = update_meta_value(conn, 'Normalization time', datetime.now().strftime(TIME_FORMAT))
    conn = update_meta_value(conn, 'precursor_normalization_method', 'median')
    conn = update_meta_value(conn, 'protein_normalization_method', 'DirectLFQ')
    conn = update_meta_value(conn, 'is_normalized', True)
    conn = update_meta_value(conn, 'Normalization command', current_command)
    conn = update_meta_value(conn, 'command_log', previous_commands + current_command)
     
    conn.close()


if __name__ == '__main__':
    main()

