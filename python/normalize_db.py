
import sys
import argparse
import sqlite3
import os
import logging
from datetime import datetime
import re

import numpy as np
import pandas as pd
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


# def main():

parser = argparse.ArgumentParser(description='Perform DirectLFQ or median normalization on precursor DB.')
parser.add_argument('-m', '--method', choices=['DirectLFQ', 'median'], default='DirectLFQ',
                    help='Normalization method to use.')
parser.add_argument('db')

args = parser.parse_args()

conn = sqlite3.connect(args.db)

# get precursor table from db
df = pd.read_sql('SELECT replicateId, proteinId as protein, modifiedSequence, precursorCharge, totalAreaFragment FROM precursors', conn)
df['ion'] = df['modifiedSequence'] + '_' + df['precursorCharge'].astype(str)
precursor_ids = df[['ion', 'modifiedSequence', 'precursorCharge']].drop_duplicates()
precursor_ids = {x.ion: (x.modifiedSequence, x.precursorCharge) for x in precursor_ids.itertuples()}

# set zero areas to the mininum non-zero value
df.loc[df['totalAreaFragment'].apply(lambda x: not np.isfinite(x) or x == 0), 'totalAreaFragment'] = min(df[df['totalAreaFragment'] > 0]['totalAreaFragment'])

# pivot wider for normalization and log2 transform
input_df = df.pivot(index=['protein', 'ion'], columns='replicateId', values='totalAreaFragment')
input_df = np.log2(input_df.replace(0, np.nan))
input_df = input_df.dropna(axis = 0, how = 'all')

print("Estimating lfq intensities.")
protein_df, ion_df = lfqprot_estimation.estimate_protein_intensities(input_df,
                                                                     min_nonan=1,
                                                                     num_samples_quadratic=10,
                                                                     num_cores = 16)

# protein_df.to_csv('/home/ajm/code/DIA_QC_report/testData/normalize_db/protein_df.tsv', sep='\t', index=False)
# ion_df.to_csv('/home/ajm/code/DIA_QC_report/testData/normalize_db/ion_df.tsv', sep='\t' , index=False)

protein_df = protein_df.melt(id_vars='protein', value_name='normalizedArea')
ion_df = ion_df.reset_index().melt(id_vars=['protein', 'ion'], value_name='normalizedArea')
ion_df['modifiedSequence'] = ion_df['ion'].apply(lambda x: precursor_ids[x][0])
ion_df['precursorCharge'] = ion_df['ion'].apply(lambda x: precursor_ids[x][1])


# delete existing normalized values
cur = conn.cursor()
cur.execute('UPDATE precursors SET normalizedArea = -1')
cur.execute('UPDATE proteinQuants SET abundance = -1')
conn.commit()

# update precursor rows
cur = conn.cursor()
ion_data = [(row.normalizedArea, row.replicateId, row.protein, row.modifiedSequence, row.precursorCharge) for row in ion_df.itertuples()]
cur.executemany("UPDATE precursors SET normalizedArea = ? WHERE replicateId = ? AND proteinId = ? AND modifiedSequence = ? AND precursorCharge = ?", ion_data)
conn.commit()
# cur.executemany('INSERT INTO proteinQuants (replicateId, proteinId, abundance) VALUES (?, ?, ?) ON CONFLICT(replicateId, proteinId) DO UPDATE SET value = ')

# insert or update proteinQuant rows
protein_query = '''
    INSERT INTO proteinQuants (replicateId, proteinId, abundance)
    VALUES (?, ?, ?)
    ON CONFLICT(replicateId, proteinId) DO UPDATE SET value = ?
    '''
cur = conn.cursor()
for row in protein_df.itertuples():
    cur.execute(query, (row.replicateId, row.protein, row.normalizedArea, row.normalizedArea))
conn.commit()

# update normalization method in metadata
conn = update_meta_value(conn, 'normalization_method', args.method)
conn = update_meta_value(conn, 'is_normalized', True)

conn.close()

# if __name__ == '__main__':
#     main()

