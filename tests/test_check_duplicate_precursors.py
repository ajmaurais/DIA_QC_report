
import unittest
import re

import setup_functions
from setup_functions import TEST_DIR

from DIA_QC_report.submodules import skyline_reports
from DIA_QC_report import parse_data


class TestCheckDuplicatePrecursorsFxn(unittest.TestCase):
    PROJECT = 'Sp3'

    @staticmethod
    def precursor_db_to_data(df, use_user_set_total):
        data = dict()
        for (rep, seq, charge), group in df.groupby(['replicateName',
                                                     'modifiedSequence',
                                                     'precursorCharge']):
            precursor = f'{seq}_{charge}'
            if precursor not in data:
                data[precursor] = dict()

            areas = group['totalAreaFragment'].to_list()
            if all(areas[0] == area for area in areas):
                data[precursor][rep] = areas[0]
                continue

            for row in group.itertuples():
                assert rep not in data[precursor]

                if not use_user_set_total:
                    data[precursor][rep] = row.totalAreaFragment
                    break

                if row.userSetTotal:
                    data[precursor][rep] = row.totalAreaFragment
                    break

        return data


    def test_invalid_mode_fails(self):
        fname = f'{TEST_DIR}/data/invalid_reports/{self.PROJECT}_by_protein_duplicate_precursor_quality.tsv'
        with self.assertLogs(skyline_reports.LOGGER) as cm:
            precursors = skyline_reports.PrecursorReport().read_report(fname)

        with self.assertRaises(ValueError):
            parse_data.check_duplicate_precursors(precursors, 'x')


    def test_invalid_no_user_set_report_fails(self):
        fname = f'{TEST_DIR}/data/invalid_reports/{self.PROJECT}_by_protein_invalid_no_user_set_precursor_quality.tsv'
        with self.assertLogs(skyline_reports.LOGGER) as cm:
            precursors = skyline_reports.PrecursorReport().read_report(fname)

        with self.assertLogs(parse_data.LOGGER, level='ERROR') as cm:
            self.assertIsNone(parse_data.check_duplicate_precursors(precursors, 'm'))
        output_re = re.compile(r'[0-9]+ precursor groups have no user set peak boundaries!')
        self.assertTrue(any(output_re.search(out) is not None for out in cm.output))


    def test_error_option_fails(self):
        fname = f'{TEST_DIR}/data/invalid_reports/{self.PROJECT}_by_protein_duplicate_precursor_quality.tsv'
        with self.assertLogs(skyline_reports.LOGGER) as cm:
            precursors = skyline_reports.PrecursorReport().read_report(fname)

        with self.assertLogs(parse_data.LOGGER, level='ERROR') as cm:
            self.assertIsNone(parse_data.check_duplicate_precursors(precursors, 'e'))
        output_re = re.compile(r'There are [0-9]+ non-unique precursor areas!')
        self.assertTrue(any(output_re.search(out) is not None for out in cm.output))


    def test_invalid_other_diff_report_fails(self):
        fname = f'{TEST_DIR}/data/invalid_reports/{self.PROJECT}_by_protein_invalid_other_diff_precursor_quality.tsv'
        with self.assertLogs(skyline_reports.LOGGER) as cm:
            precursors = skyline_reports.PrecursorReport().read_report(fname)

        with self.assertLogs(parse_data.LOGGER, level='ERROR') as cm:
            self.assertIsNone(parse_data.check_duplicate_precursors(precursors, 'm'))
        output_re = re.compile(r'There are [0-9]+ precursor rows which are not unique!')
        self.assertTrue(any(output_re.search(out) is not None for out in cm.output))


    def test_use_user_set_total_option(self):
        fname = f'{TEST_DIR}/data/invalid_reports/{self.PROJECT}_by_protein_duplicate_precursor_quality.tsv'
        with self.assertLogs(skyline_reports.LOGGER) as cm:
            precursors = skyline_reports.PrecursorReport().read_report(fname)

        with self.assertLogs(parse_data.LOGGER) as cm:
            unique_precursors = parse_data.check_duplicate_precursors(precursors, 'm')

        self.assertTrue(any(re.search(r'There are [0-9]+ non-unique precursor areas!', out) for out in cm.output))
        self.assertTrue(any(re.search(r'After selecting precursors with user set peak boundaries, [0-9]+', out) for out in cm.output))
        self.assertTrue(any('Found "UserSetTotal" column.' in out for out in cm.output))

        test_data = self.precursor_db_to_data(unique_precursors, use_user_set_total=False)
        gt_data = self.precursor_db_to_data(precursors, use_user_set_total=True)

        self.assertEqual(len(gt_data), len(test_data))
        for precursor in gt_data:
            self.assertEqual(len(gt_data[precursor]), len(test_data[precursor]))
            for rep in gt_data[precursor]:
                self.assertAlmostEqual(gt_data[precursor][rep], test_data[precursor][rep])


    def test_first_option(self):
        fname = f'{TEST_DIR}/data/invalid_reports/{self.PROJECT}_by_protein_duplicate_precursor_quality.tsv'
        with self.assertLogs(skyline_reports.LOGGER) as cm:
            precursors = skyline_reports.PrecursorReport().read_report(fname)

        with self.assertLogs(parse_data.LOGGER) as cm:
            unique_precursors = parse_data.check_duplicate_precursors(precursors, 'f')

        self.assertTrue(any(re.search(r'There are [0-9]+ non-unique precursor areas!', out) for out in cm.output))
        self.assertFalse(any('Found "UserSetTotal" column.' in out for out in cm.output))

        test_data = self.precursor_db_to_data(unique_precursors, use_user_set_total=False)
        gt_data = self.precursor_db_to_data(precursors, use_user_set_total=False)

        self.assertEqual(len(gt_data), len(test_data))
        for precursor in gt_data:
            self.assertEqual(len(gt_data[precursor]), len(test_data[precursor]))
            for rep in gt_data[precursor]:
                self.assertAlmostEqual(gt_data[precursor][rep], test_data[precursor][rep])


    def test_first_option_no_user_set(self):
        fname = f'{TEST_DIR}/data/invalid_reports/{self.PROJECT}_by_protein_invalid_no_user_set_precursor_quality.tsv'
        with self.assertLogs(skyline_reports.LOGGER) as cm:
            precursors = skyline_reports.PrecursorReport().read_report(fname)

        with self.assertLogs(parse_data.LOGGER) as cm:
            unique_precursors = parse_data.check_duplicate_precursors(precursors, 'f')

        self.assertTrue(any(re.search(r'There are [0-9]+ non-unique precursor areas!', out) for out in cm.output))
        self.assertFalse(any('Found "UserSetTotal" column.' in out for out in cm.output))

        test_data = self.precursor_db_to_data(unique_precursors, use_user_set_total=False)
        gt_data = self.precursor_db_to_data(precursors, use_user_set_total=False)

        self.assertEqual(len(gt_data), len(test_data))
        for precursor in gt_data:
            self.assertEqual(len(gt_data[precursor]), len(test_data[precursor]))
            for rep in gt_data[precursor]:
                self.assertAlmostEqual(gt_data[precursor][rep], test_data[precursor][rep])
