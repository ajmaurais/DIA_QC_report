
import re
from datetime import datetime
from collections import Counter

import pandas as pd

from .dtype import Dtype
from .logger import LOGGER
from .. import __version__ as PROGRAM_VERSION

METADATA_TIME_FORMAT = '%m/%d/%Y %H:%M:%S'

PRECURSOR_KEY_COLS = ('replicateId', 'peptideId', 'precursorCharge')

SCHEMA_VERSION = '1.18'

SCHEMA = ['PRAGMA foreign_keys = ON',
'''
CREATE TABLE replicates (
    id INTEGER PRIMARY KEY,
    replicate TEXT NOT NULL,
    project TEXT NOT NULL,
    includeRep BOOL NOT NULL DEFAULT TRUE,
    acquiredTime BLOB NOT NULL,
    acquiredRank INTEGER NOT NULL,
    ticArea REAL NOT NULL,
    UNIQUE(replicate, project) ON CONFLICT FAIL
)''',
f'''
CREATE TABLE precursors (
    replicateId INTEGER NOT NULL,
    peptideId INTEGER NOT NULL,             -- Unique for every peptide sequence and project
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
    FOREIGN KEY (replicateId) REFERENCES replicates(id) ON DELETE CASCADE
)''',
'CREATE INDEX precProjId ON precursors (peptideId, replicateId, modifiedSequence)',
'''
CREATE TABLE sampleMetadataTypes (
    annotationKey TEXT NOT NULL,
    annotationType VARCHAR(6) CHECK(
        annotationType IN ('NULL', 'BOOL', 'INT', 'FLOAT', 'STRING')
    ) NOT NULL DEFAULT 'STRING',
    PRIMARY KEY (annotationKey)
)''',
'''
CREATE TABLE sampleMetadata (
    replicateId INTEGER NOT NULL,
    annotationKey TEXT NOT NULL,
    annotationValue TEXT,
    PRIMARY KEY (replicateId, annotationKey),
    FOREIGN KEY (replicateId) REFERENCES replicates(id) ON DELETE CASCADE,
    FOREIGN KEY (annotationKey) REFERENCES sampleMetadataTypes(annotationKey) ON DELETE CASCADE
)''',
'''
CREATE TABLE metadata (
    key TEXT NOT NULL,
    value TEXT,
    PRIMARY KEY (key)
)''',
'''
CREATE TABLE proteins (
    proteinId INTEGER PRIMARY KEY,
    accession VARCHAR(25),
    name VARCHAR(50) UNIQUE,
    description VARCHAR(200)
)''',
'''
CREATE TABLE proteinQuants (
    replicateId INTEGER NOT NULL,
    proteinId INTEGER NOT NULL,
    abundance REAL,
    normalizedAbundance REAL,
    PRIMARY KEY (replicateId, proteinId),
    FOREIGN KEY (replicateId) REFERENCES replicates(id) ON DELETE CASCADE,
    FOREIGN KEY (proteinId) REFERENCES proteins(proteinId) ON DELETE CASCADE
)''',
'''
CREATE TABLE peptideToProtein (
    proteinId INTEGER NOT NULL,
    peptideId INTEGER NOT NULL,
    PRIMARY KEY (peptideId, proteinId),
    FOREIGN KEY (proteinId) REFERENCES proteins(proteinId) ON DELETE CASCADE
)''']


def get_meta_value(conn, key):
    ''' Get the value for a key from the metadata table '''
    cur = conn.cursor()
    cur.execute('SELECT value FROM metadata WHERE key == ?', (key,))
    value = cur.fetchall()
    if len(value) == 1:
        return value[0][0]
    LOGGER.error(f"Could not get key '{key}' from metadata table!")
    return None


def get_meta_value_bool(conn, key):
    ''' Get the value of a boolean variable in metadata table '''
    value = get_meta_value(conn, key)
    if value is None:
        return None

    # convert this way at some point in the future
    # value = Dtype.convert(Dtype.BOOL, value[0][0])
    value = value.lower()
    if value in ('true', '1'):
        return True
    return False


def is_normalized(conn):
    ''' Determine if metadata.is_normalized is True '''
    normalized = get_meta_value_bool(conn, 'is_normalized')
    if normalized is None or normalized == False:
        return False
    return True


def check_schema_version(conn):
    '''
    Check the database schema version and the version of dia_qc used to build the database.

    If the version of dia_qc does not match the current verision, a warning is printed to the log.
    If the SCHEMA_VERSION is different an error is printed to the log and the function returns False.

    Parameters
    ----------
    conn: sqlite.Connection
        A connection to the database.

    Returns
    -------
    schema_matches: bool
    '''

    # check the version of dia_qc used to build database
    db_program_version = get_meta_value(conn, 'dia_qc version')
    if db_program_version != PROGRAM_VERSION:
        LOGGER.warning('The database was created with dia_qc version %s but the current version is %s',
                       db_program_version, PROGRAM_VERSION)

    # check the database schema version
    db_version = get_meta_value(conn, 'schema_version')
    if db_version is None or db_version != SCHEMA_VERSION:
        LOGGER.error(f'Database schema version (%s) does not match program (%s)',
                     db_version, SCHEMA_VERSION)
        return False

    return True


def update_meta_value(conn, key, value):
    '''
    Add or update value in metadata table.

    Parameters
    ----------
    conn: sqlite3.Connection
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


def update_acquired_ranks(conn):
    '''
    Populate acquiredRank column in replicates table.

    Parameters
    ----------
    conn: sqlite3.Connection
        Database connection.
    '''

    replicates = pd.read_sql('SELECT id, acquiredTime FROM replicates WHERE includeRep == TRUE;', conn)

    # parse acquired times and add acquiredRank
    replicates['acquiredTime'] = replicates['acquiredTime'].apply(lambda x: datetime.strptime(x, METADATA_TIME_FORMAT))
    ranks = list(enumerate(replicates['acquiredTime'].sort_values().index))
    replicates['acquiredRank'] = [x[0] for x in sorted(ranks, key=lambda x: x[1])]

    acquired_ranks = [(row.acquiredRank, row.id) for row in replicates.itertuples()]
    cur = conn.cursor()
    cur.executemany('UPDATE replicates SET acquiredRank = ? WHERE id = ?', acquired_ranks)
    conn.commit()

    return update_meta_value(conn, 'replicates.acquiredRank updated', True)


def update_metadata_dtypes(conn, new_types):
    '''
    Update metadata annotationType column to fix cases where
    two projects have a different annotationTypes for the same
    annotationKey. This function will consolidate conflicting
    types using the order in the Dtype Enum class.

    Parameters
    ----------
    conn: sqlite3.Connection
        Database connection.
    new_types: dict
        A dictionary of new annotationKey, annotationType pairs.
    '''

    # Consolidate differing annotationTypes
    cur = conn.cursor()
    cur.execute('SELECT annotationKey, annotationType FROM sampleMetadataTypes;')
    existing_types = {x[0]: Dtype[x[1]] for x in cur.fetchall()}

    # consolidate new and existing data types
    for key, value in new_types.items():
        if key not in existing_types:
            existing_types[key] = value
            continue
        existing_types[key] = max(existing_types[key], value)

    # Update database
    insert_query = '''
        INSERT INTO sampleMetadataTypes (annotationKey, annotationType)
        VALUES(?, ?)
        ON CONFLICT(annotationKey) DO UPDATE SET annotationType = ?
    '''
    cur = conn.cursor()
    for annotationKey, dtype in existing_types.items():
        annotationType = str(dtype)
        cur.execute(insert_query, (annotationKey, annotationType, annotationType))
    conn.commit()

    return conn


def mark_reps_skipped(conn, reps=None, projects=None):
    '''
    Mark replicates and/or projects skipped.
    '''

    reps = reps if reps is not None else []
    projects = projects if projects is not None else []

    # retrieve data from replicates table
    cur = conn.cursor()
    cur.execute('SELECT id, replicate, project, includeRep FROM replicates;')
    db_reps = cur.fetchall()

    db_rep_index = dict()
    for i, rep in enumerate(db_reps):
        if rep in db_rep_index:
            db_rep_index[rep[1]][rep[2]] = i
        else:
            db_rep_index[rep[1]] = {rep[2]: i}

    def check_missing(var_name, missing_vals):
        if len(missing_vals) > 0:
            for rep in reps:
                LOGGER.error(f"{var_name} '{rep}' is not in database!")
            return True
        return False

    # make sure all reps and projects specified exist in db
    missing_reps = [rep for rep in reps if rep not in db_rep_index]
    all_projects = {x for xs in [list(rep.keys()) for rep in db_rep_index.values()] for x in xs}
    missing_projects = [proj for proj in projects if proj not in all_projects]
    if check_missing('Replicate', missing_reps) or check_missing('Project', missing_projects):
        return False

    rep_index_to_false = []
    for rep in reps:
        for i in db_rep_index[rep].values():
            rep_index_to_false.append(i)

    for rep, rep_projects in db_rep_index.items():
        for proj in projects:
            if proj in rep_projects:
                rep_index_to_false.append(rep_projects[proj])

    rep_index_to_false = Counter(rep_index_to_false)
    for rep_i, count in rep_index_to_false.items():
        if count > 1:
            LOGGER.warning(f"Replicate '{db_reps[rep_i][1]}' was set to be excluded {count} times!")

    LOGGER.info(f'Excluding {len(rep_index_to_false)} replicates.')
    cur = conn.cursor()
    cur.executemany('UPDATE replicates SET includeRep = FALSE WHERE id = ?;',
                    [(db_reps[rep_i][0],) for rep_i in rep_index_to_false])
    conn.commit()
    conn = update_acquired_ranks(conn)

    return True


def mark_all_reps_includced(conn):
    '''
    Set all replicates.includeRep values to TRUE and update replicates.acquiredRank if necissary.
    '''
    cur = conn.cursor()
    cur.execute('SELECT includeRep, COUNT(includeRep) FROM replicates GROUP BY includeRep;')
    include_rep_counts = Counter({x[0]: x[1] for x in cur.fetchall()})

    if 0 in include_rep_counts:
        LOGGER.info(f'Setting {include_rep_counts[0]} includeRep values to TRUE.')
        cur.execute('UPDATE replicates SET includeRep = TRUE;')
        conn.commit()
        conn = update_acquired_ranks(conn)
    else:
        LOGGER.warning(f'All replicates are already included.')


def read_wide_metadata(conn, meta_vars=None, read_all=True):
    '''
    Read metadata and replicates from database.

    Parameters
    ----------
    conn: sqlite3.Connection
    meta_vars: list
        A list of metadata variables to include.
        If None either none or all variables are included depending on the `read_all` parameter.
    read_all: bool
        If `meta_vars` is None, should all metadata variables be included or
        only replicateId and acquiredRank? Default is True.

    Returns
    -------
    metadata: pd.DataFrame
        A wide formatted Pandas DataFrame.
    '''

    acquired_ranks = pd.read_sql('''
        SELECT
            replicate,
            id as replicateId,
            project,
            acquiredRank
        FROM replicates; ''', conn)

    if meta_vars or read_all:
        meta_query = '''
        SELECT
            replicateId,
            annotationKey as key,
            annotationValue as value
        FROM sampleMetadata '''

        cur = conn.cursor()
        if meta_vars:
            meta_query += f'''\nWHERE key IN ({', '.join('?' * len(meta_vars))});'''
            cur.execute(meta_query, tuple(meta_vars))
        else:
            cur.execute(meta_query)
        metadata = cur.fetchall()

        keys = ['replicateId', 'key', 'value']
        data = {key: list() for key in keys}
        for row in metadata:
            for i, key in enumerate(keys):
                data[key].append(row[i])

        metadata = pd.DataFrame(data=data).pivot(index='replicateId', columns='key', values='value')

        cur = conn.cursor()
        cur.execute('SELECT annotationKey as key, annotationType as type FROM sampleMetadataTypes')
        types = {var: Dtype[t] for var, t in cur.fetchall()}

        for col in metadata.columns:
            metadata[col] = metadata[col].apply(types[col].convert)
            if types[col] is Dtype.STRING:
                metadata.loc[metadata[col] == '', col] = pd.NA

        acquired_ranks = acquired_ranks.join(metadata)

    return acquired_ranks


def validate_bit_mask(mask, n_options=3, n_digits=2):
    '''
    Validate bit mask command line argument

    Parameters
    ----------
    mask: str
        Bit mask.
    n_options: tuple
        Number of options. (1-3)
    n_digits: int
        Expected number of digits in mask.
    '''

    assert n_options in range(4)

    max_values = (0, 1, 3, 7)

    max_value = max_values[n_options]
    if not re.search(f"^[0-{str(max_value)}]+$", mask):
        LOGGER.error(f'Bit mask digits must be between 0 and {str(max_value)}!')
        return False

    if len(mask) != n_digits:
        LOGGER.error(f'Bit mask must be {n_digits} digits!')
        return False

    return True


def parse_bitmask_options(mask, digit_names, options):
    '''
    Parse table option bitmask.

    Parameters
    ----------
    mask: str
        Bit mask.
    digit_names tuple:
        A tuple with the name to use for each digit in the mask.
    options: tuple
        A tuple with a length of 3 mapping bit values to options.

    Return
    ------
    ret: dict
        A dictionary where the keys are the digit_names and values are
        dictionaries mapping options to booleans if their bit was set.
    '''

    assert len(digit_names) == len(mask)
    n_options = len(options)

    def _parse_bitmask(mask):
        mask_a = [int(c) for c in mask]
        for digit in mask_a:
            yield [bool(digit & (1 << i)) for i in range(n_options)]

    ret = dict()
    for key, value in zip(digit_names, _parse_bitmask(mask)):
        ret[key] = dict(zip(options, value))
    return ret


