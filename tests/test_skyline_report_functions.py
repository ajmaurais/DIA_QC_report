
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


class TestReadReportsBase(AbstractTestsBase):
    def __init__(self):
        self.report_base = None
        self.report_type = None


    @staticmethod
    def df_to_dict(df, key_cols, value_cols):
        ret = dict()
        for row in df.itertuples():
            ret['_'.join([getattr(row, key) for key in key_cols])] = {col: getattr(row, col) for col in value_cols}

        return ret


    def test_read_invariant_tsv(self):
        with self.assertLogs(skyline_reports.LOGGER) as cm:
            df = self.read_report(f'{self.report_base}_quality.tsv')

        self.assertIsNotNone(df)
        self.assertTrue(f'Found invariant {self.report_type} report...' in cm.output[0])
        self.assertTrue(f'Done reading {self.report_type}s table...' in cm.output[1])

        test_dict = self.df_to_dict(df, self.df_keys, self.df_values)
        self.assertDictEqual(test_dict, self.data)


    def test_read_english_tsv(self):
        with self.assertLogs(skyline_reports.LOGGER) as cm:
            df = self.read_report(f'{self.report_base}_quality_english.tsv')

        self.assertIsNotNone(df)
        self.assertTrue(f'Found English {self.report_type} report...' in cm.output[0])
        self.assertTrue(f'Done reading {self.report_type}s table...' in cm.output[1])

        test_dict = self.df_to_dict(df, self.df_keys, self.df_values)
        self.assertDictEqual(test_dict, self.data)


    def test_read_invariant_csv(self):
        with self.assertLogs(skyline_reports.LOGGER) as cm:
            df = self.read_report(f'{self.report_base}_quality.csv')

        self.assertIsNotNone(df)
        self.assertTrue(f'Found invariant {self.report_type} report...' in cm.output[0])
        self.assertTrue(f'Done reading {self.report_type}s table...' in cm.output[1])

        test_dict = self.df_to_dict(df, self.df_keys, self.df_values)
        self.assertDictEqual(test_dict, self.data)


    def test_read_english_csv(self):
        with self.assertLogs(skyline_reports.LOGGER) as cm:
            df = self.read_report(f'{self.report_base}_quality_english.csv')

        self.assertIsNotNone(df)
        self.assertTrue(f'Found English {self.report_type} report...' in cm.output[0])
        self.assertTrue(f'Done reading {self.report_type}s table...' in cm.output[1])

        test_dict = self.df_to_dict(df, self.df_keys, self.df_values)
        self.assertDictEqual(test_dict, self.data)


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


