
import unittest
import os
import re
import pandas as pd
from pyarrow.lib import ArrowInvalid

import setup_functions
from setup_functions import TEST_DIR

from DIA_QC_report import skyline_report_convert
from DIA_QC_report.submodules import skyline_reports


class TestConvertToParquet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.skyline_report_dir = f'{TEST_DIR}/data/skyline_reports'
        cls.metadata_dir = f'{TEST_DIR}/data/metadata'
        cls.work_dir = f'{TEST_DIR}/work/test_report_convert'
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)


    def assertIsParquet(self, file_path):
        try:
            with open(file_path, 'rb') as parquet_file:
                pd.read_parquet(parquet_file)
        except FileNotFoundError:
            self.fail(f"File does not exist: {file_path}")
        except ArrowInvalid:
            self.fail(f'Not a valid Parquet file: {file_path}')


    def test_replicate_report(self):
        report_base = 'Sp3_replicate_quality'
        report_ext = 'tsv'

        command = ['dia_qc', 'report_convert',
                   f'{self.skyline_report_dir}/{report_base}.{report_ext}']
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(0, result.returncode, result.stderr)

        target_path = f'{self.work_dir}/{report_base}.parquet'
        self.assertIsParquet(target_path)
        replicateReport = skyline_reports.ReplicateReport()
        with self.assertLogs(skyline_reports.LOGGER, level='INFO') as log:
            self.assertIsNotNone(replicateReport.read_report(target_path))
        self.assertIn('Found invariant replicate report', log.output[0])


    def test_precursor_report_unknown_cols(self):
        report_base = 'Sp3_DiaNN_precursor_quality'
        report_ext = 'tsv'

        command = ['dia_qc', 'report_convert', '--remove-unknown-cols',
                   f'{self.skyline_report_dir}/{report_base}.{report_ext}']
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(0, result.returncode, result.stderr)

        target_path = f'{self.work_dir}/{report_base}.parquet'
        self.assertIsParquet(target_path)
        precursorReport = skyline_reports.PrecursorReport()
        with self.assertLogs(skyline_reports.LOGGER, level='INFO') as log:
            self.assertIsNotNone(precursorReport.read_report(target_path))
        self.assertIn('Found invariant precursor report', log.output[0])

        precursorReport.quiet = True
        df_in = precursorReport.read_report(f'{self.skyline_report_dir}/{report_base}.{report_ext}')
        all_headers = set(df_in.columns)

        df_out = precursorReport.read_report(target_path)
        df_headers = set(df_out.columns)

        # check that there were columns removed from orriginal tsv
        diff_cols = all_headers - df_headers
        self.assertGreater(len(diff_cols), 0)

        # check that only unknown columns were removed
        known_cols = {col.name for col in precursorReport.columns()}
        for col in diff_cols:
            self.assertNotIn(col, known_cols)


    def test_precursor_report(self):
        report_base = 'Sp3_by_protein_precursor_quality'
        report_ext = 'tsv'

        command = ['dia_qc', 'report_convert',
                   f'{self.skyline_report_dir}/{report_base}.{report_ext}']
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(0, result.returncode, result.stderr)

        target_path = f'{self.work_dir}/{report_base}.parquet'
        self.assertIsParquet(target_path)
        precursorReport = skyline_reports.PrecursorReport()
        with self.assertLogs(skyline_reports.LOGGER, level='INFO') as log:
            self.assertIsNotNone(precursorReport.read_report(target_path))
        self.assertIn('Found invariant precursor report', log.output[0])


    def test_generic_report(self):
        report_base = 'acquired_ranks'
        report_ext = 'tsv'

        command = ['dia_qc', 'report_convert',
                   f'{self.metadata_dir}/{report_base}.{report_ext}']
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(0, result.returncode, result.stderr)

        self.assertRegex(result.stderr, fr'WARNING: Reading .*?{report_base}\.{report_ext} as a generic Skyline report')

        target_path = f'{self.work_dir}/{report_base}.parquet'
        self.assertIsParquet(target_path)