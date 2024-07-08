
import unittest
import os
import pandas as pd

import setup_functions

class TestExportGeneMatrix(unittest.TestCase):
    TEST_PROJECT = 'Strap'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_export_gene_matrix/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'
        cls.gene_id_path = f'{cls.data_dir}/metadata/prhuman2gene_2023_05_24_subset.csv'

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           clear_dir=True,
                                                           group_by_gene=True)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')


    def assertIsFile(self, path):
        if not os.path.isfile(path):
            raise AssertionError(f'File does not exist: {path}')


    def test_is_successful(self):
        prefix='test_is_sucessful'
        command = ['dia_qc', 'export_gene_matrix', f'--prefix={prefix}',
                    self.gene_id_path, self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 0)

        # make sure expected reports are generated
        for file in ('precursors', 'proteins'):
            self.assertIsFile(f'{self.work_dir}/{prefix}_{file}_unnormalized.tsv')


    def test_gene_table_tsv(self):
        prefix='test_tsv'
        command = ['dia_qc', 'export_gene_matrix',
                   f'--prefix={prefix}',
                   f'{self.data_dir}/metadata/prhuman2gene_2023_05_24_subset.tsv',
                   self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 0)

        # make sure expected reports are generated
        for file in ('precursors', 'proteins'):
            self.assertIsFile(f'{self.work_dir}/{prefix}_{file}_unnormalized.tsv')


    def test_use_gene_hash_option(self):
        # make sure id table without hash fails
        command = ['dia_qc', 'export_gene_matrix', '--addGeneUuid', self.gene_id_path, self.db_path]
        result = setup_functions.run_command(command, self.work_dir,
                                             prefix='failed_gene_id_hash_table')
        self.assertEqual(result.returncode, 1)

        # make sure correct table succeeds
        prefix='test_uuid'
        gene_id_hash_path = f'{self.data_dir}/metadata/prhuman2gene_gene_uuid_2023_05_24_subset.csv'
        command = ['dia_qc', 'export_gene_matrix', f'--prefix={prefix}', '--addGeneUuid',
                   gene_id_hash_path, self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 0)

        # make sure expected reports are generated
        for file in ('precursors', 'proteins'):
            report_path = f'{self.work_dir}/{prefix}_{file}_unnormalized.tsv'
            self.assertIsFile(report_path)

            # make sure gene_uuid column is in reports
            df = pd.read_csv(report_path, sep='\t')
            self.assertTrue('gene_uuid' in df.columns)


if __name__ == '__main__':
    unittest.main()
