
import unittest
import os

import setup_functions

from DIA_QC_report import normalize_db
from DIA_QC_report import export_tables

class TestExportTables(unittest.TestCase):
    TEST_PROJECT = 'Strap'
    RENDER_QMD = False

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_export_tables/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        # remove tables subdirectory in work_dir if necissary
        subdirectories = ['norm_tables', 'no_norm_tables', 'parquet_tables']
        for sub_dir in subdirectories:
            d = f'{cls.work_dir}/{sub_dir}'
            if os.path.isdir(d):
                for file in os.listdir(d):
                    os.remove(f'{d}/{file}')
                os.rmdir(d)

        cls.parse_result = setup_functions.setup_single_db(
            cls.data_dir, cls.work_dir, cls.TEST_PROJECT,
            clear_dir=True, group_by_gene=False
        )

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')


    def test_valid_table_options(self):
        self.assertEqual(self.parse_result.returncode, 0)

        export_command = [
            '--outputDir=no_norm_tables',
            '--precursorTables=33', '--proteinTables=33',
            '--metadataTables=11', self.db_path
        ]
        result = setup_functions.run_main(
            export_tables._main, export_command, self.work_dir,
            prefix='export_no_norm_tables', prog='dia_qc db_export'
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        table_names = [
            'proteins_long.tsv', 'precursors_long.tsv',
            'proteins_wide_unnormalized.tsv', 'precursors_wide_unnormalized.tsv',
            'metadata_long.tsv', 'metadata_wide.tsv'
        ]
        for table in table_names:
            test_path = f'{self.work_dir}/no_norm_tables/{table}'
            self.assertTrue(os.path.isfile(test_path), test_path)

        norm_result = setup_functions.run_main(
            normalize_db._main, [self.db_path], self.work_dir,
            prefix='normalize_db', prog='dia_qc normalize'
        )
        self.assertEqual(norm_result.returncode, 0, result.stderr)

        export_command[0] = '--outputDir=norm_tables'
        result = setup_functions.run_main(
            export_tables._main, export_command, self.work_dir,
            prefix='export_norm_tables', prog='dia_qc db_export'
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        table_names += ['proteins_wide_normalized.tsv', 'proteins_wide_normalized.tsv']
        for table in table_names:
            test_path = f'{self.work_dir}/norm_tables/{table}'
            self.assertTrue(os.path.isfile(test_path), test_path)


    def test_parquet_output(self):
        self.assertEqual(self.parse_result.returncode, 0)

        export_command = [
            '--outputDir=parquet_tables',
            '--precursorTables=11', '--proteinTables=11',
            '--metadataTables=11', '--outputFormat=parquet',
            self.db_path
        ]
        result = setup_functions.run_main(
            export_tables._main, export_command, self.work_dir,
            prefix='export_parquet', prog='dia_qc db_export'
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        table_names = [
            'proteins_long.parquet', 'precursors_long.parquet',
            'proteins_wide_unnormalized.parquet', 'precursors_wide_unnormalized.parquet',
            'metadata_long.parquet', 'metadata_wide.parquet'
        ]
        for table in table_names:
            test_path = f'{self.work_dir}/parquet_tables/{table}'
            self.assertTrue(os.path.isfile(test_path), test_path)

        norm_tables = ['proteins_wide_normalized.parquet', 'precursors_wide_normalized.parquet']
        for table in norm_tables:
            test_path = f'{self.work_dir}/parquet_tables/{table}'
            self.assertFalse(os.path.isfile(test_path), test_path)


if __name__ == '__main__':
    unittest.main()
