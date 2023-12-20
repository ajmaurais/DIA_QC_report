
import sys
import argparse
import sqlite3
import json
import os
import logging
from datetime import datetime
from csv import DictReader as CsvDictReader
import re

import pandas as pd
from jsonschema import validate, ValidationError

logging.basicConfig(
    level=logging.INFO, format='%(asctime)s: %(levelname)s: %(message)s'
)
LOGGER = logging.getLogger()

# sample metadata json schema
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


TIME_FORMAT = '%m/%d/%Y %H:%M:%S'

SCHEMA = [
'''
CREATE TABLE replicates (
    replicateId INTEGER PRIMARY KEY,
    replicate TEXT NOT NULL,
    project TEXT NOT NULL,
    acquiredTime BLOB NOT NULL,
    acquiredRank INTEGER NOT NULL,
    ticArea REAL NOT NULL,
    UNIQUE(replicate, project) ON CONFLICT FAIL
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
CREATE TABLE sampleMetadata (
    replicateId INTEGER NOT NULL,
    annotationKey TEXT NOT NULL,
    annotationValue TEXT,
    annotationType VARCHAR(6) CHECK( annotationType IN ('INT', 'FLOAT', 'BOOL', 'STRING')) NOT NULL DEFAULT 'STRING',
    FOREIGN KEY (replicateId) REFERENCES replicates(replicateId)
)
''',
'''
CREATE TABLE metadata (
    key TEXT NOT NULL,
    value TEXT,
    PRIMARY KEY (key)
)
''',
'''
CREATE TABLE proteins (
    proteinId INTEGER PRIMARY KEY,
    accession VARCHAR(25),
    name VARCHAR(50) UNIQUE,
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
    ''' Initialize empty database with SCHEMA '''
    conn = sqlite3.connect(fname)
    cur = conn.cursor()
    for command in SCHEMA:
        cur.execute(command)
    conn.commit()
    return conn


def pd_type_to_db_type(pdType):
    ''' Convert pandas data types into string '''
    if re.search('^int[0-9]+$', pdType):
        return 'INT'
    if re.search('^float[0-9]+$', pdType):
        return 'FLOAT'
    if pdType == 'bool':
        return 'BOOL'
    else:
        return 'STRING'


def read_metadata(fname, metadata_format=None):
    '''
    Read sample metadata file and format dataframe to be added to sampleMetadata table.

    Parameters
    ----------
    fname: str
        The path to the metadata file.
    metadata_format: str
        One of ('tsv', 'json')

    Returns
    -------
    df: pd.DataFrame
        The sample metadata dataframe.
    '''
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
    types = {k: pd_type_to_db_type(v.name) for k, v in df.dtypes.to_dict().items()}
    df = pd.melt(df, id_vars=['Replicate'], var_name='annotationKey', value_name='annotationValue',
                 value_vars=[x for x in list(df.columns) if x != 'Replicate'])
    df['annotationType'] = df['annotationKey'].apply(lambda x: types[x])

    return df


def insert_program_metadata(conn, metadata):
    '''
    Insert multiple metadata key, value pairs into the metadata table.
    If the key already exists it is overwritten.

    conn: sqlite3.Connection:
        Database connection.
    metadata: dict
        A dict with key, value pairs.
    '''
    cur = conn.cursor()
    for k, v in metadata.items():
        cur.execute(f'DELETE FROM metadata WHERE key = "{k}";')
        cur.execute(f'INSERT INTO metadata (key, value) VALUES ("{k}", "{v}")')
    conn.commit()
    return conn


def write_db(fname, replicates, precursors, protein_quants=None, sample_metadata=None,
             projectName=None, overwriteMode='error'):
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
    sample_metadata: pd.DataFrame
        Metadata dataframe (optional)
    projectName: str
        The project name to use in the replicates table.
    overwriteMode: str
        Behavior when the output file already exists.
        Choices are 'error', 'overwrite', or 'append'

    Returns
    -------
    sucess: bool
        True if all operations were sucessful, false otherwise.
    '''

    conn = None
    append = False
    if os.path.isfile(fname):
        if overwriteMode == 'error':
            LOGGER.error(f'{fname} already exists. Use the --overwriteMode option to append or overwrite')
            return False
        elif overwriteMode == 'overwrite':
            LOGGER.warning(f'{fname} already exists. Overwriting...')
            os.remove(fname)
            conn = _initialize(fname)
        elif overwriteMode == 'append':
            conn = sqlite3.connect(fname)
            append = True
        else:
            LOGGER.error(f'"{overwriteMode}" is an unknown overwriteMode!')
            return False
    else:
        conn = _initialize(fname)
            
    if protein_quants is not None:
        proteins = protein_quants[['name', 'accession', 'description']].drop_duplicates().reset_index(drop=True)
        proteinIndex = {r: i for i, r in zip(proteins['name'].index, proteins['name'])}

        # check that protein accessions in protein and precursor table match
        if set(proteins['name'].to_list()) != set(precursors['proteinName'].to_list()):
            LOGGER.error('Protein and precursor ProteinAccessions differ!')
            return False

    else:
        proteins = precursors[['proteinName', 'proteinAccession']].drop_duplicates().reset_index(drop=True)
        proteins = proteins.rename(columns={'proteinName': 'name', 'proteinAccession': 'accession'})
        proteins['description'] = None

    if projectName is None:
        projects = pd.read_sql('SELECT DISTINCT project FROM replicates;', conn)
        projectName = f'project_{len(projects.index) + 1}'
    replicates['project'] = projectName

    # log_metadata = {'parse_data git info': get_git_version()}
    log_metadata = dict()
    log_metadata[f'Add {projectName}'] = datetime.now().strftime(TIME_FORMAT)
    log_metadata[f'is_normalized'] = False


    # deal with existing proteins and replicates
    if append:
        conn = sqlite3.connect(fname)

        cur = conn.cursor()
        cur.execute(f'SELECT DISTINCT project FROM replicates WHERE project = "{projectName}"')
        existing_project = cur.fetchall()
        if len(existing_project) > 0:
            LOGGER.error(f'{projectName} already exists in db!')
            conn.close()
            return False

        # deal with existing replicates
        curReplicates = pd.read_sql('SELECT * FROM replicates;', conn)
        replicates.index = pd.RangeIndex(start=len(curReplicates.index), stop=len(curReplicates.index) + len(replicates.index), step=1)
        repIndex = {r: i for i, r in zip(replicates['replicate'].index, replicates['replicate'])}

        # deal with existing proteins
        curProteins = pd.read_sql('SELECT * FROM proteins;', conn)
        curProteinNames = set(curProteins['name'].to_list())
        proteins = proteins[proteins['name'].apply(lambda x: x not in curProteinNames)]
        proteins.index = pd.RangeIndex(start=len(curProteins.index), stop=len(curProteins.index) + len(proteins.index), step=1)
        proteinIndex = {r: i for i, r in zip(curProteins['proteinId'], curProteins['name'])}
        proteinIndex = {**proteinIndex, **{r: i for i, r in zip(proteins['name'].index, proteins['name'])}}

    else:
        # make dict for replicateId column which links the replicate and precursor tables
        repIndex = {r: i for i, r in zip(replicates['replicate'].index, replicates['replicate'])}

        # make dict for proteinId column which links the protein and proteinQuants and precursor tables
        proteinIndex = {r: i for i, r in zip(proteins['name'].index, proteins['name'])}

    conn = insert_program_metadata(conn, log_metadata)

    precursors['replicateId'] = precursors['replicateName'].apply(lambda x: repIndex[x])
    precursors['proteinId'] = precursors['proteinName'].apply(lambda x: proteinIndex[x])
    precursors = precursors[PRECURSOR_QUALITY_COLUMNS].drop_duplicates()

    if protein_quants is not None:
        protein_quants['replicateId'] = protein_quants['replicateName'].apply(lambda x: repIndex[x])
        protein_quants['proteinId'] = protein_quants['name'].apply(lambda x: proteinIndex[x])
        protein_quants = protein_quants[['replicateId', 'proteinId', 'abundance']]

    proteins.to_sql('proteins', conn, if_exists='append', index=True, index_label='proteinId')
    replicates.to_sql('replicates', conn, if_exists='append', index=True, index_label='replicateId')
    precursors.to_sql('precursors', conn, index=False, if_exists='append')

    if protein_quants is not None:
        protein_quants.to_sql('proteinQuants', conn, index=False, if_exists='append')

    if sample_metadata is not None:
        missing_metadata = [x for x in sample_metadata['Replicate'].drop_duplicates().to_list() if x not in repIndex]
        for rep in missing_metadata:
            LOGGER.warning(f'Metadata row: \"{rep}\" does not exist in replicate report!')

        sample_metadata = sample_metadata.loc[sample_metadata['Replicate'].apply(lambda x: x in repIndex),]

        sample_metadata['replicateId'] = sample_metadata['Replicate'].apply(lambda x: repIndex[x])
        sample_metadata = sample_metadata[['replicateId', 'annotationKey', 'annotationValue', 'annotationType']]
        sample_metadata.to_sql('sample_metadata', conn, index=False, if_exists='append')

    conn.close()
    return True


def main():
    parser = argparse.ArgumentParser(description='Generate QC_metrics input database from Skyline '
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
    parser.add_argument('--overwriteMode', choices=['error', 'overwrite', 'append'], default='error',
                        help='Behavior if output file already exists. '
                             'By default the script will exit with an error if the file already exists.')
    parser.add_argument('--projectName', default=None, help='Project name to use in replicates table.')
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
    replicates['acquiredTime'] = replicates['acquiredTime'].apply(lambda x: datetime.strptime(x, TIME_FORMAT))
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
    if not write_db(args.ofname, replicates, precursors,
                    protein_quants=protein_quants, sample_metadata=metadata,
                    projectName=args.projectName, overwriteMode=args.overwriteMode):
        LOGGER.error('Failed to create database!')
        sys.exit(1)
    LOGGER.info('Done writing database...')


if __name__ == '__main__':
    main()

