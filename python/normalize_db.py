
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
import directlfq.normalization as lfqnorm
import directlfq.protein_intensity_estimation as lfqprot_estimation

logging.basicConfig(
    level=logging.INFO, format='%(asctime)s: %(levelname)s: %(message)s'
)
LOGGER = logging.getLogger()

# def run_lfq(df):

def update_meta_value(conn, key, value):
    cur = conn.cursor()
    cur.execute(f'INSERT INTO metadata (key, value) VALUES ("{key}", "{value}") ON CONFLICT(key) DO UPDATE SET value = "{value}"')
    conn.commit()
    return conn

def median_normalize(df, key_cols, value_col):
    cols = df.columns.tolist()
    cols.append('normalizedArea')
    df['log2Area'] = np.log2(df[value_col])
    df['log2NormalizedArea'] = df.groupby(key_cols)['log2Area'].transform(lambda x: x / x.median())
    df['log2NormalizedArea'] = df['log2NormalizedArea'] * np.mean(df['log2Area'])
    df['normalizedArea'] = np.power(2, df['log2NormalizedArea'])
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
    df = pd.read_sql('SELECT replicateId, proteinId as protein, modifiedSequence, precursorCharge, totalAreaFragment FROM precursors', conn)

    df['ion'] = df['modifiedSequence'] + '_' + df['precursorCharge'].astype(str)
    precursor_ids = df[['ion', 'modifiedSequence', 'precursorCharge']].drop_duplicates()
    precursor_ids = {x.ion: (x.modifiedSequence, x.precursorCharge) for x in precursor_ids.itertuples()}

    # set zero areas to the mininum non-zero value
    df.loc[df['totalAreaFragment'].apply(lambda x: not np.isfinite(x) or x == 0), 'totalAreaFragment'] = min(df[df['totalAreaFragment'] > 0]['totalAreaFragment'])

    REP_COLUMN_NAME = 'replicateId'

    # pivot wider for normalization and log2 transform
    input_df = df.pivot(index=['protein', 'ion'], columns=REP_COLUMN_NAME, values='totalAreaFragment')

    input_df.columns.name = None
    input_df = input_df.reset_index()

    # input_df = input_df.drop_duplicates(subset='ion')
    input_df = lfqutils.index_and_log_transform_input_df(input_df)
    input_df = lfqutils.remove_allnan_rows_input_df(input_df)

    input_df = lfqnorm.NormalizationManagerSamplesOnSelectedProteins(input_df, num_samples_quadratic=10).complete_dataframe

    # input_df.to_csv('/home/ajm/code/DIA_QC_report/testData/normalize_db/input_df.tsv', sep='\t', index=True)
     
    LOGGER.info('Running DirectLFQ normalization.')
    protein_df, _ = lfqprot_estimation.estimate_protein_intensities(input_df,
                                                                    min_nonan=1,
                                                                    num_samples_quadratic=10,
                                                                    num_cores=cpu_count())
                                                                    # num_cores = 1)
    LOGGER.info('Finished with DirectLFQ normalization.')

    protein_df = protein_df.melt(id_vars='protein', value_name='normalizedArea', var_name='replicateId')
    # ion_df = ion_df.reset_index().melt(id_vars=['protein', 'ion'], var_name='replicateId', value_name='normalizedArea')
    # ion_df['modifiedSequence'] = ion_df['ion'].apply(lambda x: precursor_ids[x][0])
    # ion_df['precursorCharge'] = ion_df['ion'].apply(lambda x: precursor_ids[x][1])
    
    # protein_df.to_csv('/home/ajm/code/DIA_QC_report/testData/normalize_db/protein_df.tsv', sep='\t', index=False)
    # ion_df.to_cion_dfsv('/home/ajm/code/DIA_QC_report/testData/normalize_db/ion_df.tsv', sep='\t' , index=False)
    
    # median normalize precursor quants
    precursor_to_protein = df[['protein', 'ion']].drop_duplicates()
    ion_df = median_normalize(df[[REP_COLUMN_NAME, 'ion', 'totalAreaFragment']].drop_duplicates(), [REP_COLUMN_NAME], 'totalAreaFragment')
    ion_df['modifiedSequence'] = ion_df['ion'].apply(lambda x: precursor_ids[x][0])
    ion_df['precursorCharge'] = ion_df['ion'].apply(lambda x: precursor_ids[x][1])

    # join back in protein IDs
    ion_df = precursor_to_protein.merge(ion_df, how='left', left_index=False, right_index=False)

    # delete existing normalized values
    LOGGER.info('Setting existing normalizedArea values to NULL.')
    cur = conn.cursor()
    cur.execute(f'UPDATE precursors SET normalizedArea = NULL')
    cur.execute(f'UPDATE proteinQuants SET abundance = NULL')
    conn.commit()

    # update precursor rows
    LOGGER.info('Updating precursor normalizedArea values...')
    cur = conn.cursor()
    ion_data = [(row.normalizedArea, getattr(row, REP_COLUMN_NAME), row.protein, row.modifiedSequence, row.precursorCharge) for row in ion_df.itertuples()]
    cur.executemany("UPDATE precursors SET normalizedArea = ? WHERE replicateId = ? AND proteinId = ? AND modifiedSequence = ? AND precursorCharge = ?", ion_data)
    conn.commit()
    LOGGER.info('Done updating precursor normalizedArea values.')

    # insert or update proteinQuant rows
    protein_query = '''
        INSERT INTO proteinQuants (replicateId, proteinId, abundance)
        VALUES (?, ?, ?)
        ON CONFLICT(replicateId, proteinId) DO UPDATE SET abundance = ?
        '''
    LOGGER.info('Updating protein abundance values...')
    cur = conn.cursor()
    for row in protein_df.itertuples():
        cur.execute(protein_query, (row.replicateId, row.protein, row.normalizedArea, row.normalizedArea))
    conn.commit()
    LOGGER.info('Done updating protein abundance values.')

    # update normalization method in metadata
    LOGGER.info('Updating metadata...')
    conn = update_meta_value(conn, 'precursor_normalization_method', 'median')
    conn = update_meta_value(conn, 'protein_normalization_method', 'DirectLFQ')
    conn = update_meta_value(conn, 'is_normalized', True)
     
    conn.close()

if __name__ == '__main__':
    main()

