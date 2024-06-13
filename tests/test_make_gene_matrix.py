
import unittest
import os
import sqlite3
import pandas as pd

import setup_functions

class TestMakeGeneMatrix(unittest.TestCase):
    TEST_PROJECT = 'Strap'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_make_gene_matrix/'
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


    def test_is_successful(self):
        command = ['make_gene_matrix', self.gene_id_path, self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 0)


    def test_use_gene_hash_option(self):
        # make sure id table without hash fails
        command = ['make_gene_matrix', '--addGeneUuid', self.gene_id_path, self.db_path]
        result = setup_functions.run_command(command, self.work_dir,
                                             prefix='failed_gene_id_hash_table')
        self.assertEqual(result.returncode, 1)

        # make sure correct table succedes
        gene_id_hash_path = f'{self.data_dir}/metadata/prhuman2gene_gene_uuid_2023_05_24_subset.csv'
        command = ['make_gene_matrix', '--addGeneUuid', gene_id_hash_path, self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 0)


if __name__ == '__main__':
    unittest.main()
