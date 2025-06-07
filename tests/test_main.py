
import unittest
from unittest.mock import patch
import io
import contextlib

import DIA_QC_report
from DIA_QC_report.main import Main


class TestMain(unittest.TestCase):
    def _run_cmd(self, argv, exit_code=0, stdout_contains=None, stderr_contains=None):
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()
        with patch('sys.argv', argv):
            with self.assertRaises(SystemExit) as cm:
                with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
                    Main()
        self.assertEqual(cm.exception.code, exit_code)
        if stdout_contains:
            self.assertIn(stdout_contains, stdout_buf.getvalue())
        if stderr_contains:
            self.assertIn(stderr_contains, stderr_buf.getvalue())


    def test_main_help(self):
        self._run_cmd(
            ['dia_qc', '--help'], exit_code=0,
            stdout_contains='dia_qc <command> [<args>]'
        )


    def test_main_version(self):
        self._run_cmd(
            ['dia_qc', '--version'], exit_code=0,
            stdout_contains=f'dia_qc version {DIA_QC_report.__version__}'
        )


    def test_main_invalid_subcommand(self):
        self._run_cmd(
            ['dia_qc', 'not_a_command'], exit_code=2,
            stderr_contains="dia_qc: 'not_a_command' is not a valid command!"
        )


    def test_parse_command(self):
        self._run_cmd(
            ['dia_qc', 'parse', '--help'], exit_code=0,
            stdout_contains=DIA_QC_report.parse_data.COMMAND_DESCRIPTION
        )


    def test_report_convert_command(self):
        self._run_cmd(
            ['dia_qc', 'report_convert', '--help'], exit_code=0,
            stdout_contains=DIA_QC_report.skyline_report_convert.COMMAND_DESCRIPTION
        )


    def test_filter_command(self):
        self._run_cmd(
            ['dia_qc', 'filter', '--help'], exit_code=0,
            stdout_contains=DIA_QC_report.filter_replicates.COMMAND_DESCRIPTION
        )


    def test_metadata_convert_command(self):
        self._run_cmd(
            ['dia_qc', 'metadata_convert', '--help'], exit_code=0,
            stdout_contains=DIA_QC_report.metadata_convert.COMMAND_DESCRIPTION
        )


    def test_validate_command(self):
        self._run_cmd(
            ['dia_qc', 'validate', '--help'], exit_code=0,
            stdout_contains=DIA_QC_report.validate_pipeline_params.COMMAND_DESCRIPTION
        )


    def test_impute_command(self):
        self._run_cmd(
            ['dia_qc', 'impute', '--help'], exit_code=0,
            stdout_contains=DIA_QC_report.impute_missing.COMMAND_DESCRIPTION
        )

    def test_normalize_command(self):
        self._run_cmd(
            ['dia_qc', 'normalize', '--help'], exit_code=0,
            stdout_contains=DIA_QC_report.normalize_db.COMMAND_DESCRIPTION
        )


    def test_qc_qmd_command(self):
        self._run_cmd(
            ['dia_qc', 'qc_qmd', '--help'], exit_code=0,
            stdout_contains=DIA_QC_report.generate_qc_qmd.COMMAND_DESCRIPTION
        )

    def test_batch_rmd_command(self):
        self._run_cmd(
            ['dia_qc', 'batch_rmd', '--help'], exit_code=0,
            stdout_contains=DIA_QC_report.generate_batch_rmd.COMMAND_DESCRIPTION
        )


    def test_export_gene_matrix_command(self):
        self._run_cmd(
            ['dia_qc', 'export_gene_matrix', '--help'], exit_code=0,
            stdout_contains=DIA_QC_report.export_gene_matrix.COMMAND_DESCRIPTION
        )


    def test_db_export_command(self):
        self._run_cmd(
            ['dia_qc', 'db_export', '--help'], exit_code=0,
            stdout_contains=DIA_QC_report.export_tables.COMMAND_DESCRIPTION
        )


    def test_coverage(self):
        main_methods = [
            method for method in dir(Main)
            if callable(getattr(Main, method)) and not method.startswith("__")
        ]

        self.assertEqual(set(main_methods), DIA_QC_report.main.SUBCOMMANDS,
                         msg="Main class methods do not match expected methods.")

        for method in main_methods:
            with self.subTest(subcommand=method):
                self.assertTrue(
                    hasattr(self, f'test_{method}_command'),
                    f"Test for {method} command is missing."
                )