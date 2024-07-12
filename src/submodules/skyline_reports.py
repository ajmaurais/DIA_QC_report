
from csv import Sniffer
from datetime import datetime
from numpy import float64, int8
import pandas as pd

from .logger import LOGGER
from .dia_db_utils import METADATA_TIME_FORMAT


PRECURSOR_QUALITY_REQUIRED_COLUMNS = {'ReplicateName': 'replicateName',
                                      'ProteinAccession': 'proteinAccession',
                                      'ProteinName': 'proteinName',
                                      'ModifiedSequence': 'modifiedSequence',
                                      'PrecursorCharge': 'precursorCharge',
                                      'PrecursorMz': 'precursorMz',
                                      'AverageMassErrorPPM': 'averageMassErrorPPM',
                                      'TotalAreaFragment': 'totalAreaFragment',
                                      'TotalAreaMs1': 'totalAreaMs1',
                                      'BestRetentionTime': 'rt',
                                      'MinStartTime': 'minStartTime',
                                      'MaxEndTime': 'maxEndTime',
                                      'MaxFwhm': 'maxFwhm',
                                      'LibraryDotProduct': 'libraryDotProduct',
                                      'IsotopeDotProduct': 'isotopeDotProduct'}

PRECURSOR_QUALITY_NUMERIC_COLUMNS = ['precursorMz', 'averageMassErrorPPM', 'totalAreaFragment',
                                     'totalAreaMs1', 'rt', 'minStartTime', 'maxEndTime',
                                     'maxFwhm', 'libraryDotProduct', 'isotopeDotProduct']


REPLICATE_QUALITY_REQUIRED_COLUMNS = {'Replicate': 'replicate',
                                      'AcquiredTime': 'acquiredTime',
                                      'TicArea': 'ticArea'}

# PROTEIN_QUANTS_REQUIRED_COLUMNS = {'ProteinAccession': 'accession',
#                                    'Protein': 'name',
#                                    'ProteinDescription': 'description',
#                                    'ReplicateName': 'replicateName',
#                                    'ProteinAbundance': 'abundance'}

LANGUAGES = ('English',)

REPLICATE_LANGUAGE_TO_INVARIANT = {'Replicate': {LANGUAGES[0]: 'Replicate'},
                                   'FileName': {LANGUAGES[0]: 'File Name'},
                                   'AcquiredTime': {LANGUAGES[0]: 'Acquired Time'},
                                   'TicArea': {LANGUAGES[0]: 'Total Ion Current Area'}}

PRECURSOR_LANGUAGE_TO_INVARIANT = {'ProteinAccession': {LANGUAGES[0]: 'Protein Accession'},
                                  'ProteinName': {LANGUAGES[0]: 'Protein Name'},
                                  'ProteinGene': {LANGUAGES[0]: 'Protein Gene'},
                                  'Precursor': {LANGUAGES[0]: 'Precursor'},
                                  'Peptide': {LANGUAGES[0]: 'Peptide'},
                                  'PrecursorCharge': {LANGUAGES[0]: 'Precursor Charge'},
                                  'PrecursorMz': {LANGUAGES[0]: 'Precursor Mz'},
                                  'ModifiedSequence': {LANGUAGES[0]: 'Modified Sequence'},
                                  'ReplicateName': {LANGUAGES[0]: 'Replicate Name'},
                                  'AverageMassErrorPPM': {LANGUAGES[0]: 'Average Mass Error PPM'},
                                  'TotalAreaFragment': {LANGUAGES[0]: 'Total Area Fragment'},
                                  'TotalAreaMs1': {LANGUAGES[0]: 'Total Area MS1'},
                                  'BestRetentionTime': {LANGUAGES[0]: 'Best Retention Time'},
                                  'MinStartTime': {LANGUAGES[0]: 'Min Start Time'},
                                  'MaxEndTime': {LANGUAGES[0]: 'Max End Time'},
                                  'UserSetTotal': {LANGUAGES[0]: 'User Set Total'},
                                  'MaxFwhm': {LANGUAGES[0]: 'Max Fwhm'},
                                  'NormalizedArea': {LANGUAGES[0]: 'Normalized Area'},
                                  'LibraryDotProduct': {LANGUAGES[0]: 'Library Dot Product'},
                                  'IsotopeDotProduct': {LANGUAGES[0]: 'Isotope Dot Product'}}

def detect_delim(file):
    '''
    Detect csv/tsv file delimiter.

    Parameters
    ----------
    file: file handle

    Returns
    -------
    delimiter: str
    '''
    first_line = next(file).strip()
    file.seek(0)
    sniffer = Sniffer()
    dialect = sniffer.sniff(first_line)
    return dialect.delimiter


def check_df_columns(df, required_cols, df_name=None):
    ''' Check that df has all of required_cols. '''
    all_good = True
    df_cols = set(df.columns.to_list())
    for col in required_cols:
        if col not in df_cols:
            LOGGER.error(f'Missing required column: "{col}"' + f' in {df_name}' if df_name else '')
            all_good = False
    return all_good


def _detect_language(df, col_dict):
    # first check invariant
    col_languages = {col: set() for col in df.columns}
    for col in col_languages:
        if col in col_dict:
            col_languages[col].add('invariant')

    # iterate through every column in df
    for col in col_languages:
        # iterate through every language
        for lang in LANGUAGES:
            matches = []
            for tran_col in col_dict.values():
                if col == tran_col[lang]:
                    matches.append(tran_col)

            # populate set with every match for column
            if len(matches) == 1:
                col_languages[col].add(lang)
                continue
            if len(matches) > 1:
                raise RuntimeError('Ambigious column header: {col}')

    # remove cols with no matches
    col_languages = {col: langs for col, langs in col_languages.items() if len(langs) > 0}

    if len(col_languages) > 0:
        # If multiple language matches, report invariant or first match
        if all('invariant' in langs for langs in col_languages.values()):
            return 'invariant'
        for lang in LANGUAGES:
            if all(lang in langs for langs in col_languages.values()):
                return lang

    return None


def read_replicate_report(fname):
    with open(fname, 'r') as inF:
        df = pd.read_csv(inF, sep=detect_delim(inF))

    # Translate report into invariant format
    report_language = _detect_language(df, REPLICATE_LANGUAGE_TO_INVARIANT)
    LOGGER.info(f'Found {report_language} replicate report...')
    if report_language != 'invariant':
        col_dict = {}
        for inv_col, langs in REPLICATE_LANGUAGE_TO_INVARIANT.items():
            col_dict[langs[report_language]] = inv_col
        df = df.rename(columns=col_dict)

        # convert acquired time from 12 hr into 24 hr format
        if report_language == 'English':
            df['AcquiredTime'] = df['AcquiredTime'].apply(lambda x: datetime.strptime(x, '%m/%d/%Y %I:%M:%S %p'))
            df['AcquiredTime'] = df['AcquiredTime'].apply(lambda x: x.strftime(METADATA_TIME_FORMAT))

    df = df.rename(columns=REPLICATE_QUALITY_REQUIRED_COLUMNS)
    if not check_df_columns(df, REPLICATE_QUALITY_REQUIRED_COLUMNS.values(), 'replicates'):
        return None
    LOGGER.info('Done reading replicates table...')

    return df


def read_precursor_report(fname, by_gene=False):
    # create dict of precursor report column types
    col_types = dict()
    for sky_col, my_col in PRECURSOR_QUALITY_REQUIRED_COLUMNS.items():
        col_types[sky_col] = float64 if my_col in PRECURSOR_QUALITY_NUMERIC_COLUMNS else str
    col_types['PrecursorCharge'] = int8
    col_types['UserSetTotal'] = bool
    for col in ['ProteinGene', 'Precursor', 'Peptide']:
        col_types[col] = str

    # read report df
    with open(fname, 'r') as inF:
        df = pd.read_csv(inF, sep=detect_delim(inF), dtype=col_types)

    # Translate report into invariant format
    report_language = _detect_language(df, PRECURSOR_LANGUAGE_TO_INVARIANT)
    LOGGER.info(f'Found {report_language} precursor report...')
    if report_language != 'invariant':
        col_dict = {}
        for inv_col, langs in PRECURSOR_LANGUAGE_TO_INVARIANT.items():
            col_dict[langs[report_language]] = inv_col
        df = df.rename(columns=col_dict)

        for col in df.columns:
            if col in col_types:
                df[col] = df[col].astype(col_types[col])

    # rename columns
    df = df.rename(columns=PRECURSOR_QUALITY_REQUIRED_COLUMNS)

    precursor_cols = list(PRECURSOR_QUALITY_REQUIRED_COLUMNS.values())
    if by_gene == 'gene':
        precursor_cols.append('ProteinGene')

    if not check_df_columns(df, precursor_cols, 'precursors'):
        return None

    LOGGER.info('Done reading precursors table...')

    return df
