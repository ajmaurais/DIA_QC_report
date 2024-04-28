
import unittest
import os
import sqlite3
import pandas as pd

import setup_functions


TEST_DIR = os.path.dirname(os.path.abspath(__file__))

class TestMakeGeneMatrix(unittest.TestCase):
    TEST_PROJECT = 'Sp3'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_make_gene_matrix/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{TEST_DIR}/data/'
        cls.gene_id_path = f'{data_dir}/metadata/prhuman2gene_2023_05_24_subset.csv'

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           clear_dir=True,
                                                           group_by_gene=True)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')


    def test_is_successful(self):
        self.assertEqual(self.parse_result.returncode, 0)

        command = ['make_gene_matrix', self.db_path]
        cls.normalize_result = setup_functions.run_command(normalize_command,
                                                           cls.work_dir,
                                                           prefix='normalize_single_proj')
        # self.assertEqual(self.normalize_result.returncode, 0)


if __name__ == '__main__':
    unittest.main()