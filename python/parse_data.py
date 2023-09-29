
import sys
import argparse
import sqlite3
import json
from typing import TypeVar
import os
import logging
from datetime import datetime
from csv import DictReader as CsvDictReader

import pandas as pd
from jsonschema import validate, ValidationError


logging.basicConfig(
    level=logging.INFO, format='%(asctime)s: %(levelname)s: %(message)s'
)
LOGGER = logging.getLogger()

METADATA_SCHEMA = {
    'type': 'object',
    'additionalProperties': {
        'type': 'object',
        'additionalProperties': {
            'oneOf': [ {'type': 'string'}, {'type': 'number'} ]
        },
        'minProperties': 1
    },
    'minProperties': 1
}


SCHEMA = [
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
    proteinId INTEGER NOT NULL,
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
    PRIMARY KEY (replicateId, proteinId, modifiedSequence, precursorCharge),
    FOREIGN KEY (replicateId) REFERENCES replicates(replicateId),
    FOREIGN KEY (proteinId) REFERENCES proteins(proteinId)
)''',
'''
CREATE TABLE metadata (
    replicateId INTEGER NOT NULL,
    annotationKey TEXT NOT NULL,
    annotationValue TEXT,
    FOREIGN KEY (replicateId) REFERENCES replicates(replicateId)
)
''',
'''
CREATE TABLE proteins (
    proteinId INTEGER PRIMARY KEY,
    accession VARCHAR(25) UNIQUE,
    name VARCHAR(50),
    description VARCHAR(200)
)
''',
'''
CREATE TABLE proteinQuants (
    replicateId INTEGER NOT NULL,
    proteinId INTEGER NOT NULL,
    abundance REAL,
    PRIMARY KEY (replicateId, proteinId),
    FOREIGN KEY (replicateId) REFERENCES replicates(replicateId),
    FOREIGN KEY (proteinId) REFERENCES proteins(proteinId)
)
'''
]

PRECURSOR_QUALITY_REQUIRED_COLUMNS = {'ReplicateName': 'replicateName',
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

PRECURSOR_QUALITY_COLUMNS = ['replicateId', 'proteinId',
                             'peptide', 'modifiedSequence', 'precursorCharge', 'precursorMz',
                             'averageMassErrorPPM', 'totalAreaFragment', 'totalAreaMs1', 'normalizedArea',
                             'rt', 'minStartTime', 'maxEndTime', 'maxFwhm',
                             'libraryDotProduct', 'isotopeDotProduct']

REPLICATE_QUALITY_REQUIRED_COLUMNS = {'Replicate': 'replicate',
                                      'AcquiredTime': 'acquiredTime',
                                      'TicArea': 'ticArea'}

PROTEIN_QUANTS_REQUIRED_COLUMNS = {'ProteinAccession': 'accession',
                                   'Protein': 'name',
                                   'ProteinDescription': 'description',
                                   'ReplicateName': 'replicateName',
                                   'ProteinAbundance': 'abundance'}

def _initialize(fname):
    conn = sqlite3.connect(fname)
    cur = conn.cursor()
    for command in SCHEMA:
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
            try:
                validate(data, METADATA_SCHEMA)
            except ValidationError as e:
                raise ValidationError(f'Invalid metadata format:\n{e.message}')

            # reshape data to so it can be converted into a DataFrame
            data = [{'Replicate':k} | v for k, v in data.items()]
    else:
        raise ValueError(f'Unknown metadata file format: {_format}')

    df = pd.DataFrame(data)
    df['Replicate'] = df['Replicate'].apply(lambda x: os.path.splitext(x)[0])
    df = pd.melt(df, id_vars=['Replicate'], var_name='annotationKey', value_name='annotationValue',
                 value_vars=[x for x in list(df.columns) if x != 'Replicate'])

    return df


def write_db(fname, replicates, precursors, protein_quants=None, metadata=None):
    '''
    Write reports to QC sqlite database.

    Parameters
    ----------
    fname: str
        The name of the database file.
    replicates: pd.DataFrame
        Replicates dataframe
    precursors: pd.DataFrame
        Precursors dataframe
    protein_quants: pd.DataFrame
        Long formated protein quants dataframe (optional)
    metadata: pd.DataFrame
        Metadata dataframe (optional)

    Returns
    -------
    sucess: bool
        True if all operations were sucessful, false otherwise.
    '''

    if os.path.isfile(fname):
        LOGGER.warning(f'{fname} already exists. Overwriting...')
        os.remove(fname)

    # make dict for replicateId column which links the replicate and precursor tables
    repIndex = {r: i for i, r in zip(replicates['replicate'].index, replicates['replicate'])}

    if protein_quants is not None:
        proteins = protein_quants[['name', 'accession', 'description']].drop_duplicates().reset_index(drop=True)
        proteinIndex = {r: i for i, r in zip(proteins['accession'].index, proteins['accession'])}

        # check that protein accessions in protein and precursor table match
        if set(proteins['accession'].to_list()) != set(precursors['proteinAccession'].to_list()):
            LOGGER.error('Protein and precursor ProteinAccessions differ!')
            return False

        protein_quants['replicateId'] = protein_quants['replicateName'].apply(lambda x: repIndex[x])
        protein_quants['proteinId'] = protein_quants['accession'].apply(lambda x: proteinIndex[x])
        protein_quants = protein_quants[['replicateId', 'proteinId', 'abundance']]


    else:
        proteins = precursors[['protein', 'accession']].drop_duplicates().reset_index(drop=True)
        proteins['description'] = None
        proteinIndex = {r: i for i, r in zip(proteins['accession'].index, replicates['accession'])}

    precursors['replicateId'] = precursors['replicateName'].apply(lambda x: repIndex[x])
    precursors['proteinId'] = precursors['proteinAccession'].apply(lambda x: proteinIndex[x])
    precursors = precursors[PRECURSOR_QUALITY_COLUMNS].drop_duplicates()

    conn = _initialize(fname)
    replicates.to_sql('replicates', conn, if_exists='append', index=True, index_label='replicateId')
    proteins.to_sql('proteins', conn, if_exists='append', index=True, index_label='proteinId')
    precursors.to_sql('precursors', conn, index=False, if_exists='append')

    if protein_quants is not None:
        protein_quants.to_sql('proteinQuants', conn, index=False, if_exists='append')

    if metadata is not None:
        query='SELECT replicate, replicateId FROM replicates'
        rep_ids = {fname: rep_id for fname, rep_id in conn.execute(query).fetchall()}
        metadata['replicateId'] = metadata['Replicate'].apply(lambda x: rep_ids[x])
        metadata = metadata[['replicateId', 'annotationKey', 'annotationValue']]
        metadata.to_sql('metadata', conn, index=False, if_exists='append')

    conn.close()
    return True


def main():
    parser = argparse.ArgumentParser(description='Generate QC_metrics input database from '
                                                 'precursor_quality and replicate_quality reports.')
    parser.add_argument('-m', '--metadata', default=None,
                        help='Annotations corresponding to each file.')
    parser.add_argument('--metadataFormat', default=None, choices=('json', 'tsv'),
                        help='Specify metadata file format. '
                             'By default the format is infered from the file extension.')
    parser.add_argument('--proteins', default=None,
                        help='Long formated protein abundance report.')
    parser.add_argument('-o', '--ofname', default='data.db3',
                        help='Output file name. Default is ./data.db3')
    parser.add_argument('replicates', help='Skyline replicate_report')
    parser.add_argument('precursors', help='Skyline precursor_report')
    args = parser.parse_args()

    # read annotatiuons if applicable
    metadata = None
    if args.metadata is not None:
        try:
            metadata = read_metadata(args.metadata, args.metadataFormat)
        except (ValueError, ValidationError) as e:
            LOGGER.error(e)
            sys.exit(1)

    # read replicates
    replicates = pd.read_csv(args.replicates, sep='\t')
    replicates = replicates[REPLICATE_QUALITY_REQUIRED_COLUMNS.keys()]
    replicates = replicates.rename(columns=REPLICATE_QUALITY_REQUIRED_COLUMNS)
    LOGGER.info('Done reading replicates table...')

    # parse acquired times and add acquiredRank
    replicates['acquiredTime'] = replicates['acquiredTime'].apply(lambda x: datetime.strptime(x, '%m/%d/%Y %H:%M:%S'))
    ranks = [(rank, i) for rank, i in enumerate(replicates['acquiredTime'].sort_values().index)]
    replicates['acquiredRank'] = [x[0] for x in sorted(ranks, key=lambda x: x[1])]

    # read precursor quality report
    precursors = pd.read_csv(args.precursors, sep='\t')
    precursors = precursors.rename(columns=PRECURSOR_QUALITY_REQUIRED_COLUMNS)
    LOGGER.info('Done reading precursors table...')

    # read protein quants
    protein_quants = None
    if args.proteins:
        protein_quants = pd.read_csv(args.proteins, sep='\t')
        protein_quants = protein_quants.rename(columns=PROTEIN_QUANTS_REQUIRED_COLUMNS)
        protein_quants = protein_quants[PROTEIN_QUANTS_REQUIRED_COLUMNS.values()]
        LOGGER.info('Done reading proteins table...')

    # build database
    LOGGER.info('Writing database...')
    if not write_db(args.ofname, replicates, precursors, protein_quants=protein_quants, metadata=metadata):
        LOGGER.error('Failed to create database!')
        sys.exit(1)
    LOGGER.info('Done writing database...')


if __name__ == '__main__':
    main()

