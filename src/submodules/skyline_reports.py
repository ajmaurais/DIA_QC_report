
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


    def get_alias(self, alias_language):
        '''
        Get the column alias in the specified language.

        Parameters
        ----------
        alias_language: str:
            The language of the alias.

        Return
        ------
        alias: str

        Raises
        ------
        KeyError
            If `alias_language` is not a language.
        '''

        if alias_language == 'invariant':
            return self.skyline_name
        alias = self.skyline_aliases.get(alias_language)
        if alias is None:
            raise KeyError(f'Unknown language: {alias_language} for column: {self.name}')
        return alias


class SkylineReport(ABC):
    '''
    Abstract base class for Skyline reports

    Attributes
    ----------
    report_name: str
        The name of the report.
    quiet: bool
        Should log information be printed?
    '''

    def __init__(self, report_name=None, quiet=False):
        self._columns = list()
        self.report_name = report_name
        self._languages = set()
        self.quiet = quiet


    def set_columns(self, columns):
        '''
        Set the _columns and _languages attributes
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


    def check_df_columns(self, headers, language=None, quiet=False):
        ''' Check that df has all of required_cols. '''
        all_good = True
        df_cols = set(headers)
        for col in self.required_columns():
            col_name = col.get_alias(language) if language else col.name
            if col_name not in df_cols:
                if quiet:
                    return False
                self._error(f'Missing required column: "{col_name}"' + f' in {self.report_name}' if self.report_name else '')
                all_good = False
        return all_good


    def _error(self, message):
        if self.quiet:
            self.RuntimeError(message)
        else:
            LOGGER.error(message)


    def _info(self, message):
        if not self.quiet:
            LOGGER.info(message)


    @staticmethod
    def read_headers(file, delim=None):
        # read the first line from the file
        if delim is None:
            delim = detect_delim(file)
        headers = next(file).strip().split(delim)
        file.seek(0)
        return headers


    def detect_language(self, headers):
        # first check invariant
        df_languages = {col: set() for col in headers}
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
            for lang in self._languages:
                if all(lang in langs for langs in df_languages.values()):
                    return lang

        return None


    @abstractmethod
    def read_report(self, fname):
        pass


class ReplicateReport(SkylineReport):
    def __init__(self, quiet=False):
        super().__init__(report_name='replicate', quiet=quiet)

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


    def read_report(self, fname, return_invariant=False, remove_unknown_cols=False, quiet=False):
        if fname.endswith('.parquet'):
            df = pd.read_parquet(fname)
            report_language = self.detect_language(df.columns)

        else:
            with open(fname, 'r') as inF:
                delim = detect_delim(inF)
                report_language = self.detect_language(self.read_headers(inF, delim=delim))
                df = pd.read_csv(inF, sep=delim)

        self._info(f'Found {report_language} {self.report_name} report...')

        # Translate report into invariant format
        if report_language != 'invariant':
            col_dict = {}
            for col in self.columns():
                col_dict[col.skyline_aliases[report_language]] = col.skyline_name
            df = df.rename(columns=col_dict)

            # convert acquired time from 12 hr into 24 hr format
            if report_language == 'English':
                df['AcquiredTime'] = df['AcquiredTime'].apply(lambda x: datetime.strptime(x, '%m/%d/%Y %I:%M:%S %p'))
                df['AcquiredTime'] = df['AcquiredTime'].apply(lambda x: x.strftime(METADATA_TIME_FORMAT))

        if remove_unknown_cols:
            known_invariant_cols = {col.skyline_name for col in self.columns()}
            df = df[[col for col in df.columns if col in known_invariant_cols]]

        if return_invariant:
            return df

        df = df.rename(columns={col.skyline_name: col.name for col in self.columns()})
        if not self.check_df_columns(df.columns):
            return None
        self._info('Done reading replicates table...')

        return df


class PrecursorReport(SkylineReport):
    def __init__(self, quiet=False):
        super().__init__(report_name='precursor', quiet=quiet)

        columns = [ReportColumn('replicateName', 'ReplicateName'),
                   ReportColumn('proteinAccession', 'ProteinAccession', is_required=False),
                   ReportColumn('proteinName', 'ProteinName'),
                   ReportColumn('proteinGene', 'ProteinGene', is_required=False),
                   ReportColumn('modifiedSequence', 'ModifiedSequence'),
                   ReportColumn('precursorCharge', 'PrecursorCharge', is_numeric=False),
                   ReportColumn('precursorMz', 'PrecursorMz', is_numeric=True),
                   ReportColumn('averageMassErrorPPM', 'AverageMassErrorPPM', is_numeric=True),
                   ReportColumn('totalAreaFragment', 'TotalAreaFragment', is_numeric=True),
                   ReportColumn('userSetTotal', 'UserSetTotal', is_required=False),
                   ReportColumn('totalAreaMs1', 'TotalAreaMs1', is_numeric=True),
                   ReportColumn('rt', 'BestRetentionTime', is_numeric=True),
                   ReportColumn('minStartTime', 'MinStartTime', is_numeric=True),
                   ReportColumn('maxEndTime', 'MaxEndTime', is_numeric=True),
                   ReportColumn('maxFwhm', 'MaxFwhm', is_numeric=True),
                   ReportColumn('libraryDotProduct', 'LibraryDotProduct', is_numeric=True),
                   ReportColumn('isotopeDotProduct', 'IsotopeDotProduct', is_numeric=True)]

        english_col_names = {'proteinAccession': 'Protein Accession',
                             'proteinName': 'Protein Name',
                             'proteinGene': 'Protein Gene',
                             'precursor': 'Precursor',
                             'peptide': 'Peptide',
                             'precursorCharge': 'Precursor Charge',
                             'precursorMz': 'Precursor Mz',
                             'modifiedSequence': 'Modified Sequence',
                             'replicateName': 'Replicate Name',
                             'averageMassErrorPPM': 'Average Mass Error PPM',
                             'totalAreaFragment': 'Total Area Fragment',
                             'totalAreaMs1': 'Total Area MS1',
                             'rt': 'Best Retention Time',
                             'minStartTime': 'Min Start Time',
                             'maxEndTime': 'Max End Time',
                             'userSetTotal': 'User Set Total',
                             'maxFwhm': 'Max Fwhm',
                             'libraryDotProduct': 'Library Dot Product',
                             'isotopeDotProduct': 'Isotope Dot Product'}

        for i in range(len(columns)):
            columns[i].skyline_aliases['English'] = english_col_names[columns[i].name]

        self.set_columns(columns)


    def read_report(self, fname, by_gene=False,
                    remove_unknown_cols=False, return_invariant=False,
                    quiet=False):
        # read report df
        if fname.endswith('.parquet'):
            df = pd.read_parquet(fname)
            report_language = self.detect_language(df.columns)
            self._info(f'Found {report_language} {self.report_name} report...')

            taf_col = self._columns['totalAreaFragment'].get_alias(report_language)
            df[taf_col] = df[taf_col].astype(float64)
            tams1_col = self._columns['totalAreaMs1'].get_alias(report_language)
            df[tams1_col] = df[tams1_col].astype(float64)

        else:
            with open(fname, 'r') as inF:
                # detect report delim
                delim = detect_delim(inF)

                # detect report language
                report_language = self.detect_language(self.read_headers(inF, delim=delim))
                self._info(f'Found {report_language} {self.report_name} report...')

                # create dtypes dict for detected language
                col_types = dict()
                for col in self.columns():
                    col_types[col.get_alias(report_language)] = float64 if col.is_numeric else str
                col_types[self._columns['precursorCharge'].get_alias(report_language)] = int8
                col_types[self._columns['userSetTotal'].get_alias(report_language)] = bool

                df = pd.read_csv(inF, sep=detect_delim(inF), dtype=col_types, low_memory=False)

        # Translate report into invariant format
        if report_language != 'invariant':
            col_dict = {}
            for col in self.columns():
                col_dict[col.skyline_aliases[report_language]] = col.skyline_name
            df = df.rename(columns=col_dict)

        if remove_unknown_cols:
            known_invariant_cols = {col.skyline_name for col in self.columns()}
            df = df[[col for col in df.columns if col in known_invariant_cols]]

        if return_invariant:
            return df

        # rename columns
        df = df.rename(columns={col.skyline_name: col.name for col in self.columns()})

        if by_gene:
            def protein_uid(row):
                return row.proteinName if pd.isna(row.proteinGene) else row.proteinGene
            df['proteinName'] = df.apply(protein_uid, axis=1)

        if not self.check_df_columns(df.columns):
            return None

        self._info('Done reading precursors table...')

        return df
