
import unittest
from io import StringIO
import pandas as pd
from pandas.testing import assert_frame_equal

from setup_functions import TEST_DIR, AbstractTestsBase

import DIA_QC_report.submodules.skyline_reports as skyline_reports


class TestDetectDelim(unittest.TestCase):
    def setUp(self):
        self.data = pd.DataFrame(data={'a':['x'], 'b':['y'], 'c':[0.1], 'd':[1]})


    @unittest.expectedFailure
    def test_pd_assertion(self):
        df = pd.DataFrame(data={'a':['x'], 'b':['y'], 'c':[0.1], 'd':[0]})
        assert_frame_equal(self.data, df)


    def do_tests(self, test_string, expected_delim):
        # make sure delim is correct
        test_file_stream = StringIO(test_string)
        delim = skyline_reports.detect_delim(test_file_stream)
        self.assertEqual(delim, expected_delim)

        # make sure file stream is unmodified
        self.assertEqual(test_file_stream.read(), test_string)

        # make sure df is as expected
        test_file_stream.seek(0)
        assert_frame_equal(self.data, pd.read_csv(test_file_stream, sep=delim))


    def test_tsv(self):
        self.do_tests('a\tb\tc\td\nx\ty\t0.1\t1', '\t')
        self.do_tests('a\tb\tc\td\r\nx\ty\t0.1\t1', '\t')


    def test_csv(self):
        self.do_tests('a,b,c,d\nx,y,0.1,1', ',')
        self.do_tests('a,b,c,d\r\nx,y,0.1,1', ',')


class TestDetectLanguage(unittest.TestCase):
    def test_invariant_replicates(self):
        df = pd.read_csv(f'{TEST_DIR}/data/skyline_reports/Strap_replicate_quality.tsv', sep='\t')
        lang = skyline_reports._detect_language(df, skyline_reports.REPLICATE_LANGUAGE_TO_INVARIANT)
        self.assertEqual(lang, 'invariant')


    def test_english_replicates(self):
        df = pd.read_csv(f'{TEST_DIR}/data/skyline_reports/Strap_replicate_quality_english.tsv', sep='\t')
        lang = skyline_reports._detect_language(df, skyline_reports.REPLICATE_LANGUAGE_TO_INVARIANT)
        self.assertEqual(lang, 'English')


    def test_invariant_precursors(self):
        df = pd.read_csv(f'{TEST_DIR}/data/skyline_reports/Strap_by_protein_precursor_quality.tsv', sep='\t')
        lang = skyline_reports._detect_language(df, skyline_reports.PRECURSOR_LANGUAGE_TO_INVARIANT)
        self.assertEqual(lang, 'invariant')


    def test_english_precursors(self):
        df = pd.read_csv(f'{TEST_DIR}/data/skyline_reports/Strap_by_protein_precursor_quality_english.tsv', sep='\t')
        lang = skyline_reports._detect_language(df, skyline_reports.PRECURSOR_LANGUAGE_TO_INVARIANT)
        self.assertEqual(lang, 'English')


    def test_no_matching_cols(self):
        df = pd.read_csv(f'{TEST_DIR}/data/skyline_reports/Strap_replicate_quality_english.tsv', sep='\t')

        col_dict = {'dummy0': {skyline_reports.LANGUAGES[0]: 'Dummy1'},
                    'dummy1': {skyline_reports.LANGUAGES[0]: 'Dummy2'},
                    'dummy2': {skyline_reports.LANGUAGES[0]: 'Dummy3'},
                    'dummy3': {skyline_reports.LANGUAGES[0]: 'Dummy4'}}

        lang = skyline_reports._detect_language(df, col_dict)
        self.assertIsNone(lang)


class TestReadReportsBase(AbstractTestsBase):
    def __init__(self):
        self.report_base = None
        self.report_type = None
        self.data = None
        self.df_keys = None
        self.df_values = None


    @staticmethod
    def df_to_dict(df, key_cols, value_cols):
        ret = dict()
        for row in df.itertuples():
            ret['_'.join([str(getattr(row, key)) for key in key_cols])] = {col: getattr(row, col) for col in value_cols}

        return ret


    def do_test(self, language, ext, places=6, col_deltas=None):
        suffix = ''
        if language == 'English':
            suffix = '_english'

        with self.assertLogs(skyline_reports.LOGGER) as cm:
            df = self.read_report(f'{self.report_base}_quality{suffix}.{ext}')

        self.assertIsNotNone(df)
        self.assertTrue(f'Found {language} {self.report_type} report...' in cm.output[0])
        self.assertTrue(f'Done reading {self.report_type}s table...' in cm.output[1])

        test_dict = self.df_to_dict(df, self.df_keys, self.df_values)
        self.assertDataDictEqual(test_dict, self.data, places=places, col_deltas=col_deltas)


class TestReadReplicateReport(unittest.TestCase, TestReadReportsBase):
    TEST_PROJECT = 'Strap'

    def setUp(self):
        # file vars
        self.report_type = 'replicate'
        self.report_base = f'{TEST_DIR}/data/skyline_reports/{self.TEST_PROJECT}_{self.report_type}'

        # report column vars
        self.df_keys = ['FileName']
        self.df_values = list(skyline_reports.REPLICATE_QUALITY_REQUIRED_COLUMNS.values())

        # setup gt data
        df = pd.read_csv(f'{TEST_DIR}/data/intermediate_files/{self.TEST_PROJECT}_replicates_df.tsv', sep='\t')
        self.data = self.df_to_dict(df, self.df_keys, self.df_values)

        # setup read_report fxn
        self.read_report = skyline_reports.read_replicate_report


    def test_read_invariant_tsv(self):
        self.do_test('invariant', 'tsv', places=6)


    def test_read_english_tsv(self):
        self.do_test('English', 'tsv', places=6)


    def test_read_invariant_csv(self):
        self.do_test('invariant', 'csv', places=6)


    def test_read_english_csv(self):
        self.do_test('English', 'csv', places=6)


class TestReadPrecursorReport(unittest.TestCase, TestReadReportsBase):
    TEST_PROJECT = 'Strap'

    COL_DIFFS = {'precursorMz': 0.0001,
                 'averageMassErrorPPM': 0.1,
                 'totalAreaFragment': 0.5,
                 'totalAreaMs1': 0.5,
                 'normalizedArea': 0.5,
                 'rt': 0.01,
                 'minStartTime': 0.01,
                 'maxEndTime': 0.01,
                 'maxFwhm': 0.1,
                 'libraryDotProduct': 0.0001,
                 'isotopeDotProduct': 0.0001}

    def setUp(self):
        # file vars
        self.report_type = 'precursor'
        self.report_base = f'{TEST_DIR}/data/skyline_reports/{self.TEST_PROJECT}_by_protein_{self.report_type}'

        # report column vars
        self.df_keys = ['replicateName', 'proteinAccession', 'modifiedSequence', 'precursorCharge']
        self.df_values = list(skyline_reports.PRECURSOR_QUALITY_NUMERIC_COLUMNS)

        # setup gt data
        df = pd.read_csv(f'{TEST_DIR}/data/intermediate_files/{self.TEST_PROJECT}_precursors_df.tsv', sep='\t')

        self.data = self.df_to_dict(df, self.df_keys, self.df_values)

        # setup read_report fxn
        self.read_report = skyline_reports.read_precursor_report


    def test_read_invariant_tsv(self):
        self.do_test('invariant', 'tsv', places=6)


    def test_read_english_tsv(self):
        self.do_test('English', 'tsv', col_deltas=self.COL_DIFFS)


    def test_read_invariant_csv(self):
        self.do_test('invariant', 'csv', places=6)


    def test_read_english_csv(self):
        self.do_test('English', 'csv', col_deltas=self.COL_DIFFS)


if __name__ == '__main__':
    unittest.main()
