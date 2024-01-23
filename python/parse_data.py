
import sys
import argparse
import sqlite3
import json
import os
import logging
from datetime import datetime
from csv import DictReader as CsvDictReader
import re
from enum import Enum

import pandas as pd
from jsonschema import validate, ValidationError

logging.basicConfig(
    level=logging.INFO, format='%(asctime)s - %(filename)s %(funcName)s - %(levelname)s: %(message)s'
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

PRECURSOR_KEY_COLS = ('replicateId', 'modifiedSequence', 'precursorCharge')

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
f'''
CREATE TABLE precursors (
    replicateId INTEGER NOT NULL,
    peptide VARCHAR(60),
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
    PRIMARY KEY ({', '.join(PRECURSOR_KEY_COLS)}),
    FOREIGN KEY (replicateId) REFERENCES replicates(replicateId)
)''',
'''
CREATE TABLE sampleMetadata (
    replicateId INTEGER NOT NULL,
    annotationKey TEXT NOT NULL,
    annotationValue TEXT,
    annotationType VARCHAR(6) CHECK( annotationType IN ('BOOL', 'INT', 'FLOAT', 'STRING')) NOT NULL DEFAULT 'STRING',
    PRIMARY KEY (replicateId, annotationKey),
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
    normalizedAbundance REAL,
    PRIMARY KEY (replicateId, proteinId),
    FOREIGN KEY (replicateId) REFERENCES replicates(replicateId),
    FOREIGN KEY (proteinId) REFERENCES proteins(proteinId)
)
''',
'''
CREATE TABLE peptideToProtein (
    proteinId INTEGER NOT NULL,
    modifiedSequence VARCHAR(200) NOT NULL,
    PRIMARY KEY (modifiedSequence, proteinId),
    FOREIGN KEY (proteinId) REFERENCES proteins(proteinId),
    FOREIGN KEY (modifiedSequence) REFERENCES precursors(modifiedSequence)
)''']

PRECURSOR_QUALITY_REQUIRED_COLUMNS = {'ReplicateName': 'replicateName',
                                      'ProteinAccession': 'proteinAccession',
                                      'ProteinName': 'proteinName',
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

PRECURSOR_QUALITY_NUMERIC_COLUMNS = ['precursorMz', 'averageMassErrorPPM', 'totalAreaFragment',
                                     'totalAreaMs1', 'normalizedArea', 'rt', 'minStartTime',
                                     'maxEndTime', 'maxFwhm', 'libraryDotProduct', 'isotopeDotProduct']

PRECURSOR_QUALITY_COLUMNS = list(PRECURSOR_KEY_COLS) + PRECURSOR_QUALITY_NUMERIC_COLUMNS

REPLICATE_QUALITY_REQUIRED_COLUMNS = {'Replicate': 'replicate',
                                      'AcquiredTime': 'acquiredTime',
                                      'TicArea': 'ticArea'}

PROTEIN_QUANTS_REQUIRED_COLUMNS = {'ProteinAccession': 'accession',
                                   'Protein': 'name',
                                   'ProteinDescription': 'description',
                                   'ReplicateName': 'replicateName',
                                   'ProteinAbundance': 'abundance'}

DUPLICATE_PRECURSOR_CHOICES = ('e', 'm', 'f', 'i')

def _initialize(fname):
    ''' Initialize empty database with SCHEMA '''
    conn = sqlite3.connect(fname)
    cur = conn.cursor()
    for command in SCHEMA:
        try:
            cur.execute(command)
        except sqlite3.OperationalError as e:
            raise RuntimeError(f'Error running command: {command}')

    conn.commit()
    return conn


NA_RE = re.compile('^(NA|NULL|#N/A|NaN)$')
BOOL_RE = re.compile('^(true|True|TRUE|false|FALSE|False)$')
INT_RE = re.compile('^[+\-]?\d+$')
FLOAT_RE = re.compile('^[+\-]?(\d+(.\d*)?|.\d+|\d+(.\d*)?[eE][+\-]?\d+)$')


class Dtype(Enum):
    NULL=0
    BOOL=1
    INT=2
    FLOAT=3
    STRING=4

    def __str__(self):
        return self.name

    def __lt__(self, rhs):
        if type(rhs) == type(self):
            return self.value < rhs.value

        raise ValueError(f'Cannot compare {type(self)} to {type(rhs)}!')


    def __ge__(self, rhs):
        if type(rhs) == type(self):
            return self.value >= rhs.value

        raise ValueError(f'Cannot compare {type(self)} to {type(rhs)}!')

    @staticmethod
    def infer_type(s):
        '''
        Infer datatype for string.

        Parameters
        ----------
        s: str
            The string to test.

        Returns
        -------
        Dtype object
        '''
        if s == '' or NA_RE.search(s):
            return Dtype.NULL
        if BOOL_RE.search(s):
            return Dtype.BOOL
        if INT_RE.search(s):
            return Dtype.INT
        if FLOAT_RE.search(s):
            return Dtype.FLOAT
        return Dtype.STRING


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

    # pivot longer
    df = pd.melt(df, id_vars=['Replicate'], var_name='annotationKey', value_name='annotationValue',
                 value_vars=[x for x in list(df.columns) if x != 'Replicate'])

    # determine annotationValue types
    df['annotationType'] = df['annotationValue'].apply(lambda x: Dtype.infer_type(x))
    df['annotationValue'] = df.apply(lambda x: pd.NA if x.annotationType is Dtype.NULL else x.annotationValue, axis=1)
    types=dict()
    for annotationKey, group in df.groupby('annotationKey'):
        thisType = max(group['annotationType'].tolist())
        types[annotationKey] = str(Dtype.BOOL if thisType is Dtype.NULL else thisType)

    df['annotationType'] = df['annotationKey'].apply(lambda x: types[x])

    return df


def insert_program_metadata(conn, metadata):
    '''
    Insert multiple metadata key, value pairs into the metadata table.
    If the key already exists it is overwritten.

    Parameters
    ----------
    conn: sqlite3.Connection:
        Database connection.
    metadata: dict
        A dict with key, value pairs.
    '''
    cur = conn.cursor()
    for key, value in metadata.items():
        cur.execute(f'''
            INSERT INTO metadata
                (key, value) VALUES ("{key}", "{value}")
            ON CONFLICT(key) DO UPDATE SET value = "{value}" ''')
    conn.commit()
    return conn


def update_acquired_ranks(conn):
    '''
    Populate acquiredRank column in replicates table.

    Parameters
    ----------
    conn: sqlite3.Connection:
        Database connection.
    '''

    replicates = pd.read_sql('SELECT replicateId, acquiredTime FROM replicates;', conn)

    # parse acquired times and add acquiredRank
    replicates['acquiredTime'] = replicates['acquiredTime'].apply(lambda x: datetime.strptime(x, TIME_FORMAT))
    ranks = [(rank, i) for rank, i in enumerate(replicates['acquiredTime'].sort_values().index)]
    replicates['acquiredRank'] = [x[0] for x in sorted(ranks, key=lambda x: x[1])]

    acquired_ranks = [(row.acquiredRank, row.replicateId) for row in replicates.itertuples()]
    cur = conn.cursor()
    cur.executemany('UPDATE replicates SET acquiredRank = ? WHERE replicateId = ?', acquired_ranks)
    conn.commit()

    insert_program_metadata(conn, {'replicates.acquiredRank updated': True})

    return conn


def update_metadata_dtypes(conn):
    '''
    Update metadata annotationType column to fix cases where
    two projects have a different annotationTypes for the same
    annotationKey. This function will consolidate conflicting
    types using the order in the Dtype Enum class.

    Parameters
    ----------
    conn: sqlite3.Connection:
        Database connection.
    '''

    # Consolidate differing annotationTypes
    types = pd.read_sql('SELECT DISTINCT annotationKey, annotationType FROM sampleMetadata;', conn)
    types['dtype'] = types['annotationType'].apply(lambda x: Dtype[x])
    types = types.groupby('annotationKey')['dtype'].max().apply(lambda x: str(x)).reset_index()
    types = [(row.dtype, row.annotationKey) for row in types.itertuples()]

    # Update database
    cur = conn.cursor()
    cur.executemany('UPDATE sampleMetadata SET annotationType = ? WHERE annotationKey = ?', types)
    conn.commit()

    return conn


def _add_index_column(col, index, df_name=None, index_name=None):
    '''
    Add a column to dataframe from a dict index.
    Also handle KeyErrors when values in key_col do not exist as a key in index.

    Parameters
    ----------
    col: pd.Series
        A Series with key values mapping to index.
    index: dict
        A dict mapping key_col to values in index.
    df_name: str
        The name of the df to use in any error messages.
    index_name: str
        The name of the index to use in any error messages.

    Returns
    -------
    pd.Series: The new column to add to df, or None if there was a missing value in index.
    '''

    # Check that all values in col exist in index
    table_names = '' if df_name is None or index_name is None else f' for {df_name} in {index_name}'
    all_good = True
    for key in col.drop_duplicates().tolist():
        if key not in index:
            LOGGER.error(f'Missing required value: {key}{table_names}!')
            all_good = False

    # Exit now if there are one or more missing values
    if not all_good:
        return None

    return col.apply(lambda x: index[x])


def write_db(fname, replicates, precursors, protein_quants=None, sample_metadata=None,
             projectName=None, overwriteMode='error'):
    '''
    Write reports to precursor sqlite database.

    Parameters
    ----------
    fname: str
        The name of the database file.
    replicates: pd.DataFrame
        Replicates dataframe
    precursors: pd.DataFrame
        Precursors dataframe
    protein_quants: pd.DataFrame
        Long formatted protein quants dataframe (optional)
    sample_metadata: pd.DataFrame
        Metadata dataframe (optional)
    projectName: str
        The project name to use in the replicates table.
    overwriteMode: str
        Behavior when the output file already exists.
        Choices are 'error', 'overwrite', or 'append'

    Returns
    -------
    success: bool
        True if all operations were successful, false otherwise.
    '''

    # Metadata to add to db
    current_command = ' '.join(sys.argv)
    log_metadata = {'command_log': current_command}

    # Initialize database if applicable and get database connection
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

            # get commands previously run on database
            cur = conn.cursor()
            cur.execute('SELECT value FROM metadata WHERE key == "command_log"')
            log_metadata['command_log'] = f'{cur.fetchall()[0][0]}\n{log_metadata["command_log"]}'
        else:
            LOGGER.error(f'"{overwriteMode}" is an unknown overwriteMode!')
            return False
    else:
        conn = _initialize(fname)

    # create proteins table
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

    # if projectName is not specified, set it to project_(n + 1)
    if projectName is None:
        projects = pd.read_sql('SELECT DISTINCT project FROM replicates;', conn)
        projectName = f'project_{len(projects.index) + 1}'
    replicates['project'] = projectName

    # populate some metadata values not that we have the projectName
    log_metadata[f'Add {projectName} time'] = datetime.now().strftime(TIME_FORMAT)
    log_metadata[f'Add {projectName} command'] = current_command
    log_metadata['replicates.acquiredRank updated'] = False
    log_metadata[f'is_normalized'] = False # set this to False because we are adding unnormalized data

    # deal with existing replicates, proteins, and protein to precursor pairs
    if append:
        conn = sqlite3.connect(fname)

        cur = conn.cursor()
        cur.execute('SELECT DISTINCT project FROM replicates WHERE project = ?', (projectName,))
        existing_project = cur.fetchall()
        if len(existing_project) > 0:
            LOGGER.error(f'{projectName} already exists in db!')
            conn.close()
            return False

        # deal with existing replicates
        curReplicates = pd.read_sql('SELECT * FROM replicates;', conn)
        replicates.index = pd.RangeIndex(start=len(curReplicates.index),
                                         stop=len(curReplicates.index) + len(replicates.index),
                                         step=1)
        repIndex = {r: i for i, r in zip(replicates['replicate'].index, replicates['replicate'])}

        # deal with existing proteins
        curProteins = pd.read_sql('SELECT proteinId, name FROM proteins;', conn)
        curProteinNames = set(curProteins['name'].to_list())
        proteins = proteins[proteins['name'].apply(lambda x: x not in curProteinNames)]
        proteins.index = pd.RangeIndex(start=len(curProteins.index),
                                       stop=len(curProteins.index) + len(proteins.index),
                                       step=1)
        proteinIndex = {r: i for i, r in zip(curProteins['proteinId'], curProteins['name'])}
        proteinIndex = {**proteinIndex, **{r: i for i, r in zip(proteins['name'].index, proteins['name'])}}

        # deal with existing protein to peptide pairs
        curPairs = pd.read_sql('''SELECT
                                    prot.name as proteinName, p.modifiedSequence
                                FROM peptideToProtein p
                                LEFT JOIN proteins prot ON prot.proteinId = p.proteinId;''',
                                conn)
        newPairs = pd.merge(curPairs, precursors[['proteinName', 'modifiedSequence']].drop_duplicates(),
                            how='right', indicator='exists')
        newPairs = newPairs[newPairs['exists'] == 'right_only']

    else:
        # make dict for replicateId column which links the replicate and precursor tables
        repIndex = {r: i for i, r in zip(replicates['replicate'].index, replicates['replicate'])}

        # make dict for proteinId column which links the protein and proteinQuants and precursor tables
        proteinIndex = {r: i for i, r in zip(proteins['name'].index, proteins['name'])}

        newPairs = precursors[['proteinName', 'modifiedSequence']].drop_duplicates()

    if protein_quants is None:
        protein_quants = precursors.groupby(['replicateName', 'proteinName'])['totalAreaFragment'].sum().reset_index()
        protein_quants = protein_quants.rename(columns={'proteinName': 'name',
                                                        'totalAreaFragment': 'abundance'})

    # add replicateId and proteinId columns to protein_quants
    if (rep_id_col := _add_index_column(protein_quants['replicateName'], repIndex,
                                        'protein_quants', 'repIndex')) is None:
        return False
    protein_quants['replicateId'] = rep_id_col

    if (prot_id_col := _add_index_column(protein_quants['name'], proteinIndex,
                                         'protein_quants', 'proteinIndex')) is None:
        return False
    protein_quants['proteinId'] = prot_id_col

    protein_quants = protein_quants[['replicateId', 'proteinId', 'abundance']]

    # add replicateId and proteinId columns to precursors and drop duplicate sequences
    if (rep_id_col := _add_index_column(precursors['replicateName'], repIndex,
                                        'precursors', 'repIndex')) is None:
        return False
    precursors['replicateId'] = rep_id_col

    if (prot_id_col := _add_index_column(precursors['proteinName'], proteinIndex,
                                         'precursors', 'proteinIndex')) is None:
        return False
    precursors['proteinId'] = prot_id_col

    precursors = precursors[PRECURSOR_QUALITY_COLUMNS].drop_duplicates()

    # add proteinId column and make peptideToProtein df
    if (prot_id_col := _add_index_column(newPairs['proteinName'], proteinIndex,
                                         'peptideToProtein', 'proteinIndex')) is None:
        return False
    newPairs['proteinId'] = prot_id_col
    peptideToProtein = newPairs[['proteinId', 'modifiedSequence']]

    peptideToProtein.to_sql('peptideToProtein', conn, if_exists='append', index=False)
    proteins.to_sql('proteins', conn, if_exists='append', index=True, index_label='proteinId')
    protein_quants.to_sql('proteinQuants', conn, index=False, if_exists='append')
    replicates.to_sql('replicates', conn, if_exists='append', index=True, index_label='replicateId')
    precursors.to_sql('precursors', conn, index=False, if_exists='append')

    if sample_metadata is not None:
        missing_metadata = [x for x in sample_metadata['Replicate'].drop_duplicates().to_list() if x not in repIndex]
        for rep in missing_metadata:
            LOGGER.warning(f'Metadata row: \"{rep}\" does not exist in replicate report!')

        sample_metadata = sample_metadata.loc[sample_metadata['Replicate'].apply(lambda x: x in repIndex),]

        if (rep_id_col := _add_index_column(sample_metadata['Replicate'], repIndex,
                                            'sample_metadata', 'repIndex')) is None:
            return False
        sample_metadata['replicateId'] = rep_id_col
        sample_metadata = sample_metadata[['replicateId', 'annotationKey', 'annotationValue', 'annotationType']]
        sample_metadata.to_sql('sampleMetadata', conn, index=False, if_exists='append')

    conn = insert_program_metadata(conn, log_metadata)
    conn = update_acquired_ranks(conn)

    if append:
        conn = update_metadata_dtypes(conn)

    conn.close()
    return True


def check_duplicate_precursors(precursors, mode):
    '''
    Check that all precursors are unique.
    The same precursor could be assigned to multiple proteins in Skyline.
    When you manually adjust the peak integration boundries for a precursor, Skyline
    will only adjust the integration boundries for the precursor you are viewing.
    If the precursor is assigned to more than one protein the boundries will not be synced.
    It is not possible to store duplicate precursors in the database which have different peak
    areas. We must ask the user to specify how to deal with duplicate precursors if they exist.

    Parameters
    ----------
    precursors: pd.DataFrame
        Precursor dataframe.
    mode: str
        One of ['e', 'm', 'f', 'i'] indicating how to handle duplicates.

    Returns
    -------
    precursors: pd.DataFrame
        precursors df with duplicates removed.

    '''

    if mode not in DUPLICATE_PRECURSOR_CHOICES:
        raise ValueError(f'"{mode}" is an invlaid mode!')

    LOGGER.info('Checking that all precursors are unique...')

    key_cols = [PRECURSOR_QUALITY_REQUIRED_COLUMNS[k] for k in ('ReplicateName', 'ModifiedSequence', 'PrecursorCharge')]
    n_unique = len(precursors[key_cols].drop_duplicates().index)

    unique_precursors = precursors.drop_duplicates(subset=key_cols + ['totalAreaFragment'])
    n_unique_areas = len(unique_precursors.index)

    all_cols = [v for v in PRECURSOR_QUALITY_REQUIRED_COLUMNS.values() if v not in ['proteinAccession', 'proteinName']]
    n_unique_rows = len(precursors[all_cols].drop_duplicates().index)

    if n_unique_areas != n_unique_rows:
        LOGGER.error(f'There are {n_unique_rows - n_unique} precursor rows which are not unique!')
        return None

    if n_unique == n_unique_areas:
        return precursors

    # If there duplicate precursors return None
    if mode == 'e':
        LOGGER.error(f'There are {n_unique_areas - n_unique} non-unique precursors!')
        return None

    LOGGER.warning(f'There are {n_unique_areas - n_unique} non-unique precursors!')

    precursors = precursors.set_index(keys=key_cols)
    keep_cols = [x for x in precursors.columns if x not in PRECURSOR_QUALITY_NUMERIC_COLUMNS]

    # Select first precursor
    if mode == 'f':
        LOGGER.warning(f'Selecting first occurance of duplicate precursors...')
        areas = precursors[PRECURSOR_QUALITY_NUMERIC_COLUMNS].loc[~precursors.index.duplicated()]

        ret = precursors[keep_cols].merge(areas, how='left', left_index=True, right_index=True)
        return ret.reset_index()

    counts = unique_precursors.groupby(key_cols).apply(lambda x: len(x))
    non_unique = precursors[precursors.index.isin(counts[counts > 1].index)]
    unique = precursors[precursors.index.isin(counts[counts == 1].index)]

    # Use precursor peak areas that were manually set.
    if mode == 'm':
        # check UserSetTotal column
        if 'UserSetTotal' not in precursors.columns:
            LOGGER.error('Need "UserSetTotal" column to select precursors with user set peak boundries!')
            return None
        LOGGER.info('Found "UserSetTotal" column.')
        if precursors.dtypes['UserSetTotal'] != 'bool':
            LOGGER.error('"UserSetTotal" column must be type bool!')
            return None

        # check that each group has at least 1 precursor area that was manually set
        manually_adj = non_unique.groupby(key_cols)['UserSetTotal'].agg(lambda x: any(x))
        if not all(manually_adj):
            LOGGER.error(f'{sum(~manually_adj)} precursor groups have no user set peak boundries!')
            return None

        # remove precursors which were not manually set
        selections = non_unique[non_unique['UserSetTotal']]
        non_unique_len = len(selections.index)
        selections = selections[~selections.index.duplicated(keep='first')]

        if non_unique_len != len(selections.index):
            LOGGER.warning(f'After selecing precursors with user set peak boundries, '
                           f'{non_unique_len - len(selections.index)} non-unique precursors remain.')

    # Interactive selection mode
    if mode == 'i':
        have_user_set_col = 'UserSetTotal' in non_unique.columns
        groups = non_unique.groupby(key_cols)
        n_groups = len(groups)
        int_re = re.compile('[0-9]+')
        selections = []
        for group_i, (name, group) in enumerate(groups):
            valid_choice = False
            sele_int = None
            while not valid_choice:
                # format table of precursor options
                options = [['Selection', 'Protein', 'Sequence', 'Charge', 'Area']]
                if have_user_set_col:
                    options[-1].append('UserSet')
                for i, row in enumerate(group.itertuples()):
                    options.append([f'{i})', row.proteinName, name[1], str(name[2]), str(row.totalAreaFragment)])
                    if have_user_set_col:
                        options[-1].append(str(row.UserSetTotal))
                col_widths = [max([len(options[row][col]) for row in range(len(options))]) for col in range(len(options[0]))]

                # write precursor options
                sys.stdout.write(f'Choose precursor ({group_i + 1} of {n_groups})\n')
                sys.stdout.write(f'Replicate: "{name[0]}"\n')
                for row in options:
                    sys.stdout.write('  {}\n'.format('  '.join(row[col].ljust(col_widths[col]) for col in range(len(row)))))

                selection = input('Select: ')
                if int_re.search(selection):
                    sele_int = int(selection)
                    if sele_int >= 0 and sele_int < len(group.index):
                        valid_choice = True
                        break
                sys.stdout.write('Invalid choice!\n')
            selections.append(group.iloc[sele_int])

        selections = pd.DataFrame(selections)
        selections.index = selections.index.set_names(key_cols)

    # overwrite duplicate precursors with selections
    non_unique = non_unique[keep_cols].merge(selections[PRECURSOR_QUALITY_NUMERIC_COLUMNS],
                                             how='left', left_index=True, right_index=True)

    # Merge duplicates back into main dataframe
    ret = pd.concat([unique, non_unique])
    ret = ret.reset_index()

    return ret


def check_df_columns(df, required_cols, df_name=None):
    ''' Check that df has all of required_cols. '''
    all_good = True
    df_cols = set(df.columns.to_list())
    for col in required_cols:
        if col not in df_cols:
            LOGGER.error(f'Missing required column: "{col}"' + f' in {df_name}' if df_name else '')
            all_good = False
    return all_good


def main():
    parser = argparse.ArgumentParser(description='Generate QC and batch correction database from '
                                                 'Skyline precursor_quality and replicate_quality reports.')
    parser.add_argument('-m', '--metadata', default=None,
                        help='Annotations corresponding to each file.')
    parser.add_argument('--metadataFormat', default=None, choices=('json', 'tsv'),
                        help='Specify metadata file format. '
                             'By default the format is inferred from the file extension.')
    parser.add_argument('--proteins', default=None,
                        help='Long formatted protein abundance report. '
                              'If no protein report is given, proteins are quantified in the database '
                              'by summing all the precursors belonging to that protein.')
    parser.add_argument('-o', '--ofname', default='data.db3',
                        help='Output file name. Default is ./data.db3')
    parser.add_argument('--overwriteMode', choices=['error', 'overwrite', 'append'], default='error',
                        help='Behavior if output file already exists. '
                             'By default the script will exit with an error if the file already exists.')
    parser.add_argument('-d', '--duplicatePrecursors', default='e', choices=DUPLICATE_PRECURSOR_CHOICES,
                        help="How to handle precursors with the same sequence and charge, but different area. "
                             "'e' for error, 'm' to use the peak area with user adjusted integration boundries, "
                             "'f' to select first occurance of duplicate precursors, "
                             "and 'i' to interactively choose which peak are to use. 'e' is the default.")
    parser.add_argument('--projectName', default=None, help='Project name to use in replicates table.')
    parser.add_argument('replicates', help='Skyline replicate_report')
    parser.add_argument('precursors', help='Skyline precursor_report')
    args = parser.parse_args()

    # read annotations if applicable
    metadata = None
    if args.metadata is not None:
        try:
            metadata = read_metadata(args.metadata, args.metadataFormat)
        except (ValueError, ValidationError) as e:
            LOGGER.error(e)
            sys.exit(1)

    # read replicates
    replicates = pd.read_csv(args.replicates, sep='\t')
    replicates = replicates.rename(columns=REPLICATE_QUALITY_REQUIRED_COLUMNS)
    if not check_df_columns(replicates, REPLICATE_QUALITY_REQUIRED_COLUMNS.values(), 'replicates'):
        sys.exit(1)

    replicates = replicates[REPLICATE_QUALITY_REQUIRED_COLUMNS.values()]
    replicates['acquiredRank'] = -1 # this will be updated later
    LOGGER.info('Done reading replicates table...')

    # read precursor quality report
    precursors = pd.read_csv(args.precursors, sep='\t')
    LOGGER.info('Done reading precursors table...')
    precursors = precursors.rename(columns=PRECURSOR_QUALITY_REQUIRED_COLUMNS)
    if not check_df_columns(precursors, PRECURSOR_QUALITY_REQUIRED_COLUMNS.values(), 'precursors'):
        sys.exit(1)

    if (precursors := check_duplicate_precursors(precursors, args.duplicatePrecursors)) is None:
        sys.exit(1)

    # read protein quants
    protein_quants = None
    if args.proteins:
        protein_quants = pd.read_csv(args.proteins, sep='\t')
        protein_quants = protein_quants.rename(columns=PROTEIN_QUANTS_REQUIRED_COLUMNS)
        if not check_df_columns(protein_quants, PROTEIN_QUANTS_REQUIRED_COLUMNS, 'protein_quants'):
            sys.exit(1)

        protein_quants = protein_quants[PROTEIN_QUANTS_REQUIRED_COLUMNS.keys()]
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

