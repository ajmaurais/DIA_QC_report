
import argparse
import sqlite3
import json
from typing import TypeVar
import pandas as pd
import os
import logging
from datetime import datetime
from csv import DictReader as CsvDictReader

PandasDataFrame = TypeVar('pandas.core.frame.DataFrame')

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s: %(levelname)s: %(message)s"
)
LOGGER = logging.getLogger()

schema = [
'''
CREATE TABLE replicates (
    replicateId INTEGER PRIMARY KEY,
    replicate TEXT UNIQUE,
    acquiredTime BLOB NOT NULL,
    acquiredRank INTEGER NOT NULL,
    ticArea REAL NOT NULL
)''',
'''
CREATE TABLE precursors (
    replicateId INTEGER NOT NULL,
    proteinName VARCHAR(50),
    proteinAccession VARCHAR(25),
    peptide VARCHAR(50),
    modifiedSequence VARCHAR(200) NOT NULL,
    precursorCharge INTEGER NOT NULL,
    precursorMz REAL,
    averageMassErrorPPM REAL,
    totalAreaFragment REAL,
    totalAreaMs1 REAL,
    normalizedArea REAL,
    rt REAL,
    minStartTime REAL,
    maxEndTime REAL,
    maxFwhm REAL,
    libraryDotProduct REAL,
    isotopeDotProduct REAL,
    PRIMARY KEY (replicateId, proteinName, modifiedSequence, precursorCharge),
    FOREIGN KEY (replicateId) REFERENCES replicates(replicateId)
)''',
'''
CREATE TABLE metadata (
    replicateId INTEGER NOT NULL,
    annotationKey TEXT NOT NULL,
    annotationValue TEXT,
    FOREIGN KEY (replicateId) REFERENCES replicates(replicateId)
)
'''
]

PRECURSOR_QUALITY_COLNAMES = {'ReplicateName': 'replicateName',
                              'ProteinAccession': 'proteinAccession',
                              'ProteinName': 'proteinName',
                              'ModifiedSequence': 'modifiedSequence',
                              'Peptide': 'peptide',
                              'ModifiedSequence': 'modifiedSequence',
                              'PrecursorCharge': 'precursorCharge',
                              'PrecursorMz': 'precursorMz',
                              'AverageMassErrorPPM': 'averageMassErrorPPM',
                              'TotalAreaFragment': 'totalAreaFragment',
                              'TotalAreaMs1': 'totalAreaMs1',
                              'NormalizedArea': 'normalizedArea',
                              'BestRetentionTime': 'rt',
                              'MinStartTime': 'minStartTime',
                              'MaxEndTime': 'maxEndTime',
                              'MaxFwhm': 'maxFwhm',
                              'LibraryDotProduct': 'libraryDotProduct',
                              'IsotopeDotProduct': 'isotopeDotProduct'}

PRECURSOR_QUALITY_REQUIRED_COLUMNS = ['replicateId', 'proteinAccession', 'proteinName',
                                      'peptide', 'modifiedSequence', 'precursorCharge', 'precursorMz',
                                      'averageMassErrorPPM', 'totalAreaFragment', 'totalAreaMs1', 'normalizedArea',
                                      'rt', 'minStartTime', 'maxEndTime', 'maxFwhm',
                                      'libraryDotProduct', 'isotopeDotProduct']

REPLICATE_QUALITY_REQUIRED_COLUMNS = {'Replicate': 'replicate',
                                      'AcquiredTime': 'acquiredTime',
                                      'TicArea': 'ticArea'}

def _initialize(fname):
    conn = sqlite3.connect(fname)
    cur = conn.cursor()
    for command in schema:
        cur.execute(command)
    conn.commit()
    return conn


def read_metadata(fname, metadata_format=None):
    if metadata_format:
        _format = metadata_format
    else:
        _format = os.path.splitext(fname)[1][1:]

    if _format == 'tsv':
        with open(fname, 'r') as inF:
            data = list(CsvDictReader(inF, delimiter='\t'))
    elif _format == 'json':
        with open(fname, 'r') as inF:
            data = json.load(inF)
    else:
        raise RuntimeError(f'Unknown metadata file format: {_format}')

    df = pd.DataFrame(data)
    df['file_name'] = df["file_name"].apply(lambda x: os.path.splitext(x)[0])
    df = pd.melt(df, id_vars=['file_name'], var_name='annotationKey', value_name='annotationValue',
                 value_vars=[x for x in list(df.columns) if x != "file_name"])

    return df


def write_db(fname: str,
             replicates: PandasDataFrame,
             precursors: PandasDataFrame,
             metadata: PandasDataFrame=None):

    if os.path.isfile(fname):
        LOGGER.warning(f'{fname} already exists. Overwriting...')
        os.remove(fname)

    conn = _initialize(fname)
    replicates.to_sql('replicates', conn, if_exists='append', index=True, index_label='replicateId')
    precursors.to_sql('precursors', conn, index=False, if_exists='append')
    if metadata is not None:
        query='SELECT replicate, replicateId FROM replicates'
        rep_ids = {fname: rep_id for fname, rep_id in conn.execute(query).fetchall()}
        metadata['replicateId'] = metadata['file_name'].apply(lambda x: rep_ids[x])
        metadata = metadata[['replicateId', 'annotationKey', 'annotationValue']]
        metadata.to_sql('metadata', conn, index=False, if_exists='append')


def main():
    parser = argparse.ArgumentParser(description='Generate QC_metrics input database from '
                                                 'precursor_quality and replicate_quality reports.')
    parser.add_argument('-m', '--metadata', default=None,
                        help='Annotations corresponding to each file.')
    parser.add_argument('--metadataFormat', default=None,
                        help='Specify metadata file format. '
                             'By default the format is infered from the file extension.')
    parser.add_argument('-o', '--ofname', default='data.db3',
                        help='Output file name. Default is ./data.db3')
    parser.add_argument('replicates', help='Skyline replicate_report')
    parser.add_argument('precursors', help='Skyline precursor_report')
    args = parser.parse_args()

    # read annotatiuons if applicable
    metadata = None
    if args.metadata is not None:
        metadata = read_metadata(args.metadata, args.metadataFormat)

    # read replicates
    replicates = pd.read_csv(args.replicates, sep='\t')
    replicates = replicates[REPLICATE_QUALITY_REQUIRED_COLUMNS.keys()]
    replicates = replicates.rename(columns=REPLICATE_QUALITY_REQUIRED_COLUMNS)

    # parse acquired times and add acquiredRank
    replicates['acquiredTime'] = replicates['acquiredTime'].apply(lambda x: datetime.strptime(x, '%m/%d/%Y %H:%M:%S'))
    ranks = [(rank, i) for rank, i in enumerate(replicates["acquiredTime"].sort_values().index)]
    replicates['acquiredRank'] = [x[0] for x in sorted(ranks, key=lambda x: x[1])]
    repIndex = {r: i for i, r in zip(replicates["replicate"].index, replicates["replicate"])}

    # read precursor quality report
    precursors = pd.read_csv(args.precursors, sep='\t')
    precursors = precursors.rename(columns=PRECURSOR_QUALITY_COLNAMES)

    # add replicateId column which links the replicate and precursor tables
    precursors['replicateId'] = precursors['replicateName'].apply(lambda x: repIndex[x])
    precursors = precursors[PRECURSOR_QUALITY_REQUIRED_COLUMNS].drop_duplicates()

    write_db(args.ofname, replicates, precursors, metadata=metadata)


if __name__ == '__main__':
    main()

