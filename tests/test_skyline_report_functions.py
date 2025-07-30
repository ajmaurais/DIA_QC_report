
import unittest
from io import StringIO
import pandas as pd

from setup_functions import TEST_DIR, AbstractTestsBase

import DIA_QC_report.submodules.skyline_reports as skyline_reports


class TestDetectDelim(unittest.TestCase, AbstractTestsBase):
    def setUp(self):
        self.data = pd.DataFrame(data={'a':['x'], 'b':['y'], 'c':[0.1], 'd':[1]})


    def do_tests(self, test_string, expected_delim):
        # make sure delim is correct
        test_file_stream = StringIO(test_string)
        delim = skyline_reports.detect_delim(test_file_stream)
        self.assertEqual(delim, expected_delim)

        # make sure file stream is unmodified
        self.assertEqual(test_file_stream.read(), test_string)

        # make sure df is as expected
        test_file_stream.seek(0)
        self.assertDataFrameEqual(self.data, pd.read_csv(test_file_stream, sep=delim))


    def test_tsv(self):
        self.do_tests('a\tb\tc\td\nx\ty\t0.1\t1', '\t')
        self.do_tests('a\tb\tc\td\r\nx\ty\t0.1\t1', '\t')


    def test_csv(self):
        self.do_tests('a,b,c,d\nx,y,0.1,1', ',')
        self.do_tests('a,b,c,d\r\nx,y,0.1,1', ',')


class TestDetectLanguage(unittest.TestCase):
    def test_invariant_replicates(self):
        with open(f'{TEST_DIR}/data/skyline_reports/Strap_replicate_quality.tsv', 'r') as inF:
            headers = skyline_reports.ReplicateReport.read_headers(inF)
            lang = skyline_reports.ReplicateReport().detect_language(headers)
        self.assertEqual(lang, 'invariant')


    def test_english_replicates(self):
        with open(f'{TEST_DIR}/data/skyline_reports/Strap_replicate_quality_english.tsv', 'r') as inF:
            headers = skyline_reports.ReplicateReport.read_headers(inF)
            lang = skyline_reports.ReplicateReport().detect_language(headers)
        self.assertEqual(lang, 'English')


    def test_invariant_precursors(self):
        with open(f'{TEST_DIR}/data/skyline_reports/Strap_by_protein_precursor_quality.tsv', 'r') as inF:
            headers = skyline_reports.ReplicateReport.read_headers(inF)
            lang = skyline_reports.PrecursorReport().detect_language(headers)
        self.assertEqual(lang, 'invariant')


    def test_english_precursors(self):
        with open(f'{TEST_DIR}/data/skyline_reports/Strap_by_protein_precursor_quality_english.tsv', 'r') as inF:
            headers = skyline_reports.ReplicateReport.read_headers(inF)
            lang = skyline_reports.PrecursorReport().detect_language(headers)
        self.assertEqual(lang, 'English')


    def test_no_matching_cols(self):
        ''' Test that df with no matching headers returns None '''

        # setup test fp with dummy headers
        test_string = 'Dummy1\tDummy2\tDummy3\tDummy4\n'
        with open(f'{TEST_DIR}/data/skyline_reports/Strap_replicate_quality.tsv', 'r') as inF:
            next(inF)
            for line in inF:
                test_string += line

        test_ss = StringIO(test_string)
        lang = skyline_reports.ReplicateReport().detect_language(test_ss)

        self.assertIsNone(lang)


class TestReadReportsBase(AbstractTestsBase):
    def __init__(self):
        self.report_basename = None
        self.report_type = None
        self.data = None
        self.df_keys = None
        self.df_values = None
        self.report_class = None


    @staticmethod
    def df_to_dict(df, key_cols, value_cols):
        ret = dict()
        for row in df.itertuples():
            key = '_'.join([str(getattr(row, key)) for key in key_cols])
            ret[key] = {col: getattr(row, col) for col in value_cols}

        return ret


    def do_test(self, language, ext, places=6, col_deltas=None, test_report=None):
        suffix = ''
        if language == 'English':
            suffix = '_english'

        if test_report is None:
            _test_report = f'{self.report_basename}_quality{suffix}.{ext}'
        else:
            _test_report = test_report

        with self.assertLogs(skyline_reports.LOGGER, level='INFO') as cm:
            df = self.report_class.read_report(_test_report)

        self.assertIsNotNone(df)
        self.assertIn(f'Found {language} {self.report_type} report...', cm.output[0])
        self.assertIn(f'Done reading {self.report_type}s table...', cm.output[1])

        test_dict = self.df_to_dict(df, self.df_keys, self.df_values)
        self.assertDataDictEqual(test_dict, self.data, places=places, col_deltas=col_deltas)


    def do_quiet_test(self, language, ext, places=6, col_deltas=None):
        suffix = ''
        if language == 'English':
            suffix = '_english'

        with self.assertNoLogs(skyline_reports.LOGGER) as cm:
            old_quiet = self.report_class.quiet
            self.report_class.quiet = True
            df = self.report_class.read_report(f'{self.report_basename}_quality{suffix}.{ext}')
            self.report_class.quiet = old_quiet

        self.assertIsNotNone(df)

        test_dict = self.df_to_dict(df, self.df_keys, self.df_values)
        self.assertDataDictEqual(test_dict, self.data, places=places, col_deltas=col_deltas)


class TestReadReplicateReport(unittest.TestCase, TestReadReportsBase):
    TEST_PROJECT = 'Strap'

    def setUp(self):
        # file vars
        self.report_type = 'replicate'
        self.report_basename = f'{TEST_DIR}/data/skyline_reports/{self.TEST_PROJECT}_{self.report_type}'

        # setup SkylineReport class
        self.report_class = skyline_reports.ReplicateReport()

        # report column vars
        self.df_keys = ['fileName']
        self.df_values = [col.name for col in self.report_class.required_columns()]

        # setup gt data
        df = pd.read_csv(f'{TEST_DIR}/data/intermediate_files/{self.TEST_PROJECT}_replicates_df.tsv', sep='\t')
        self.data = self.df_to_dict(df, self.df_keys, self.df_values)


    def test_read_invariant_tsv(self):
        self.do_test('invariant', 'tsv', places=6)


    def test_read_english_tsv(self):
        self.do_test('English', 'tsv', places=6)


    def test_read_invariant_csv(self):
        self.do_test('invariant', 'csv', places=6)


    def test_read_english_csv(self):
        self.do_test('English', 'csv', places=6)


    def test_read_invariant_parquet(self):
        self.do_test('invariant', 'parquet', places=6)


    def test_read_english_parquet(self):
        self.do_test('English', 'parquet', places=6)


    def test_read_quiet(self):
        self.do_quiet_test('invariant', 'tsv', places=6)
        self.do_quiet_test('invariant', 'csv', places=6)
        self.do_quiet_test('invariant', 'parquet', places=6)
        self.do_quiet_test('English', 'csv', places=6)
        self.do_quiet_test('English', 'parquet', places=6)
        self.do_quiet_test('English', 'tsv', places=6)


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
        self.report_basename = f'{TEST_DIR}/data/skyline_reports/{self.TEST_PROJECT}_by_protein_{self.report_type}'

        # setup SkylineReport class
        self.report_class = skyline_reports.PrecursorReport()

        # report column vars
        self.df_keys = ['replicateName', 'proteinName', 'modifiedSequence', 'precursorCharge']
        self.df_values = [col.name for col in self.report_class.required_columns()]

        # setup gt data
        df = pd.read_csv(f'{TEST_DIR}/data/intermediate_files/{self.TEST_PROJECT}_precursors_df.tsv', sep='\t')

        self.data = self.df_to_dict(df, self.df_keys, self.df_values)


    def test_read_invariant_tsv(self):
        self.do_test('invariant', 'tsv', places=6)


    def test_read_english_tsv(self):
        self.do_test('English', 'tsv', col_deltas=self.COL_DIFFS)


    def test_read_invariant_csv(self):
        self.do_test('invariant', 'csv', places=6)


    def test_read_english_csv(self):
        self.do_test('English', 'csv', col_deltas=self.COL_DIFFS)


    def test_read_invariant_parquet(self):
        self.do_test('invariant', 'parquet', places=6)


    def test_read_english_parquet(self):
        self.do_test('English', 'parquet', col_deltas=self.COL_DIFFS)


    def test_read_missing_report_rows(self):
        test_tsv = f'{TEST_DIR}/data/invalid_reports/Strap_by_protein_precursor_quality_empty_row.tsv'
        with self.assertLogs(skyline_reports.LOGGER, level='INFO') as cm:
            df = self.report_class.read_report(test_tsv)

        self.assertIsNotNone(df)
        self.assertIn('Found invariant precursor report...', cm.output[0])
        self.assertIn('Removing 1 row(s) with missing UserSetTotal value', cm.output[1])
        self.assertIn('Done reading precursors table...', cm.output[2])


    def test_read_quiet(self):
        self.do_quiet_test('invariant', 'tsv', places=6)
        self.do_quiet_test('invariant', 'csv', places=6)
        self.do_quiet_test('invariant', 'parquet', places=6)
        self.do_quiet_test('English', 'csv', col_deltas=self.COL_DIFFS)
        self.do_quiet_test('English', 'parquet', col_deltas=self.COL_DIFFS)
        self.do_quiet_test('English', 'tsv', col_deltas=self.COL_DIFFS)