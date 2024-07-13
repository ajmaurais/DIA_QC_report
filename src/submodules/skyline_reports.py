
from abc import ABC, abstractmethod
from csv import Sniffer
from datetime import datetime
from numpy import float64, int8
import pandas as pd

from .logger import LOGGER
from .dia_db_utils import METADATA_TIME_FORMAT


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


class ReportColumn():
    '''
    Represents a Skyline report column.

    Attributes
    ----------
    name: str
        The name of the column. This is the same as the header of the column in the database.
    skyline_name: str
        The 'invariant' name of the column in Skyline
    is_numeric: bool
        Is the column numeric?
    is_required: bool
        Is the column required to create the database?
    skyline_aliases: dict
        A dictionary of alternate names for the column used in skyline for reports in different
        languages. The key should be the language of the report, the value should be column
        name in that language.
    '''

    def __init__(self, name, skyline_name,
                 is_numeric=False, is_required=True):
        self.name = name
        self.skyline_name = skyline_name
        self.is_numeric = is_numeric
        self.is_required = is_required
        self.skyline_aliases = dict()


class SkylineReport(ABC):
    '''
    Abstract base class for Skyline reports

    Attributes
    ----------
    report_name: str
        The name of the report.
    '''

    def __init__(self, report_name=None):
        self._columns = list()
        self.report_name = report_name
        self._languages = set()


    def set_columns(self, columns):
        '''
        Set the _columns and _languages attrubutes
        '''
        self._columns = {col.name: col for col in columns}
        if len(self._columns) < len(columns):
            raise RuntimeError('Ambiguous column names in report!')

        self._languages = {lang for col in self.columns() for lang in col.skyline_aliases}


    def required_columns(self):
        '''
        Yield a generator to required Skyline ReportColumn(s)
        '''
        for column in self._columns.values():
            if column.is_required:
                yield column


    def numeric_columns(self):
        '''
        Yield a generator to numeric Skyline ReportColumn(s)
        '''
        for column in self._columns.values():
            if column.is_numeric:
                yield column


    def columns(self):
        '''
        Yield a generator to all Skyline ReportColumn(s)
        '''
        yield from self._columns.values()


    def check_df_columns(self, df):
        ''' Check that df has all of required_cols. '''
        all_good = True
        df_cols = set(df.columns.to_list())
        for col in self.required_columns():
            if col.name not in df_cols:
                LOGGER.error(f'Missing required column: "{col}"' + f' in {self.report_name}' if self.report_name else '')
                all_good = False
        return all_good


    def _detect_language(self, df):
        # first check invariant
        df_languages = {col: set() for col in df.columns}
        invariant_cols = {col.skyline_name for col in self.columns()}
        for col in df_languages:
            if col in invariant_cols:
                df_languages[col].add('invariant')

        # iterate through every column in df
        for df_col in df_languages:
            # iterate through every language
            for lang in self._languages:
                matches = []
                for col in self.columns():
                    if df_col == col.skyline_aliases[lang]:
                        matches.append(col.skyline_aliases[lang])

                # populate set with every match for column
                if len(matches) == 1:
                    df_languages[df_col].add(lang)
                    continue
                if len(matches) > 1:
                    raise RuntimeError('Ambigious column header: {df_col}')

        # remove cols with no matches
        df_languages = {col: langs for col, langs in df_languages.items() if len(langs) > 0}

        if len(df_languages) > 0:
            # If multiple language matches, report invariant or first match
            if all('invariant' in langs for langs in df_languages.values()):
                return 'invariant'
            for lang in LANGUAGES:
                if all(lang in langs for langs in df_languages.values()):
                    return lang

        return None


    @abstractmethod
    def read_report(self, fname):
        pass



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

PRECURSOR_QUALITY_OPTIONAL_COLUMNS = {'NormalizedArea': 'normalizedArea',
                                      }

PRECURSOR_QUALITY_NUMERIC_COLUMNS = ['precursorMz', 'averageMassErrorPPM', 'totalAreaFragment',
                                     'totalAreaMs1', 'rt', 'minStartTime', 'maxEndTime',
                                     'maxFwhm', 'libraryDotProduct', 'isotopeDotProduct']

# PROTEIN_QUANTS_REQUIRED_COLUMNS = {'ProteinAccession': 'accession',
#                                    'Protein': 'name',
#                                    'ProteinDescription': 'description',
#                                    'ReplicateName': 'replicateName',
#                                    'ProteinAbundance': 'abundance'}

LANGUAGES = ('English',)

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


class ReplicateReport(SkylineReport):
    def __init__(self):
        super().__init__(report_name='replicate')

        columns = [ReportColumn('replicate', 'Replicate'),
                   ReportColumn('acquiredTime', 'AcquiredTime'),
                   ReportColumn('ticArea', 'TicArea', is_numeric=True),
                   ReportColumn('fileName', 'FileName', is_required=False)]

        english_col_names = {'replicate': 'Replicate',
                             'fileName': 'File Name',
                             'acquiredTime': 'Acquired Time',
                             'ticArea': 'Total Ion Current Area'}

        for i in range(len(columns)):
            columns[i].skyline_aliases['English'] = english_col_names[columns[i].name]

        self.set_columns(columns)


    def read_report(self, fname):
        with open(fname, 'r') as inF:
            df = pd.read_csv(inF, sep=detect_delim(inF))

        # Translate report into invariant format
        report_language = self._detect_language(df)
        LOGGER.info(f'Found {report_language} {self.report_name} report...')
        if report_language != 'invariant':
            col_dict = {}
            for col in self.columns():
                col_dict[col.skyline_aliases[report_language]] = col.skyline_name
            df = df.rename(columns=col_dict)

            # convert acquired time from 12 hr into 24 hr format
            if report_language == 'English':
                df['AcquiredTime'] = df['AcquiredTime'].apply(lambda x: datetime.strptime(x, '%m/%d/%Y %I:%M:%S %p'))
                df['AcquiredTime'] = df['AcquiredTime'].apply(lambda x: x.strftime(METADATA_TIME_FORMAT))

        df = df.rename(columns={col.skyline_name: col.name for col in self.columns()})
        if not self.check_df_columns(df):
            return None
        LOGGER.info('Done reading replicates table...')

        return df


class PrecursorReport(SkylineReport):
    def __init__(self):
        super().__init__(report_name='precursor')

    def read_report(self, fname, by_gene=False):
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
        LOGGER.info(f'Found {report_language} {self.report_name} report...')
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

        if not self.check_df_columns(df):
            return None

        LOGGER.info('Done reading precursors table...')

        return df
