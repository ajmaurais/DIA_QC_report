
import unittest
import os

import setup_functions


class OptionStep():
    def __init__(self, flag):
        self.flag = flag


class TestExportTables(unittest.TestCase):
    TEST_PROJECT = 'Strap'
    RENDER_QMD = False

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_export_tables/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        # remove tables subdirectory in work_dir if necissary
        for d in (f'{cls.work_dir}/norm_tables', f'{cls.work_dir}/no_norm_tables'):
            if os.path.isdir(d):
                for file in os.listdir(d):
                    os.remove(f'{d}/{file}')
                os.rmdir(d)

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           clear_dir=True,
                                                           group_by_gene=False)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')


    def test_valid_table_options(self):
        self.assertEqual(self.parse_result.returncode, 0)

        export_command = ['dia_qc', 'db_export', '--outputDir=no_norm_tables',
                          '--precursorTables=33', '--proteinTables=33',
                          '--metadataTables=11', self.db_path]
        result = setup_functions.run_command(export_command, self.work_dir)
        self.assertEqual(result.returncode, 0)

        table_names = ['proteins_long.tsv', 'precursors_long.tsv',
                       'proteins_wide_unnormalized.tsv', 'precursors_wide_unnormalized.tsv',
                       'metadata_long.tsv', 'metadata_wide.tsv']
        for table in table_names:
            test_path = f'{self.work_dir}/no_norm_tables/{table}'
            self.assertTrue(os.path.isfile(test_path), test_path)

        normalize_command = ['normalize_db', self.db_path]
        norm_result = setup_functions.run_command(normalize_command, self.work_dir)
        self.assertEqual(norm_result.returncode, 0)

        export_command[2] = '--outputDir=norm_tables'
        result = setup_functions.run_command(export_command, self.work_dir)
        self.assertEqual(result.returncode, 0)

        table_names += ['proteins_wide_normalized.tsv', 'proteins_wide_normalized.tsv']
        for table in table_names:
            test_path = f'{self.work_dir}/norm_tables/{table}'
            self.assertTrue(os.path.isfile(test_path), test_path)


if __name__ == '__main__':
    unittest.main()
