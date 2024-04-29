
import unittest
from unittest import mock
import argparse
import sys
import os
import sqlite3

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

class TestMakeQCrmd(unittest.TestCase):
    RENDER_RMD = False

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_generate_batch_rmd/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{TEST_DIR}/data/'

        # remove tables subdirectory in work_dir if necissary
        if os.path.isdir(f'{cls.work_dir}/tables'):
            for file in os.listdir(f'{cls.work_dir}/tables'):
                os.remove(f'{cls.work_dir}/tables/{file}')
            os.rmdir(f'{cls.work_dir}/tables')

        cls.parse_results = setup_functions.setup_multi_db(cls.data_dir,
                                                          cls.work_dir,
                                                          clear_dir=True)

        if not all(result.returncode == 0 for result in cls.parse_results):
            raise RuntimeError('Setup of test db failed!')

        normalize_command = ['normalize_db', cls.db_path]
        cls.normalize_result = setup_functions.run_command(normalize_command,
                                                           cls.work_dir,
                                                           prefix='normalize_single_proj')

        cls.conn = None
        if cls.normalize_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_is_successful(self):
        self.assertTrue(self.conn is not None)
        self.assertTrue(db_utils.is_normalized(self.conn))

        rmd_name = 'basic_test'
        command = ['generate_batch_rmd',
                   '-o', f'{rmd_name}.rmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.rmd'))

        if self.RENDER_RMD:
            render_command = ['Rscript', '-e', f"rmarkdown::render('{rmd_name}.rmd')"]
            render_result = setup_functions.run_command(render_command, self.work_dir)
            self.assertEqual(render_result.returncode, 0)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{rmd_name}.html'))


    def test_is_normalized_check(self):
        self.assertTrue(self.conn is not None)
        self.assertTrue(db_utils.is_normalized(self.conn))

        # set is_normalized metadata entry to False
        self.conn = db_utils.update_meta_value(self.conn, 'is_normalized', 'False')

        # make sure generate_batch_rmd fails
        command = ['generate_batch_rmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 1)
        self.assertTrue('Database file it not normalized!' in result.stderr.decode('utf-8'))

        # reset is_normalized metadata entry to True
        self.conn = db_utils.update_meta_value(self.conn, 'is_normalized', 'True')


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_single_batch_fails(self):
        self.assertTrue(self.conn is not None)

        # set 'Strap' project replicates skipped
        self.assertTrue(db_utils.mark_reps_skipped(self.conn, projects=('Strap',)))

        # make sure generate_batch_rmd fails
        command = ['generate_batch_rmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 1)
        self.assertTrue('Only 1 project in replicates!' in result.stderr.decode('utf-8'))

        # reset is_normalized metadata entry to True
        self.conn = db_utils.mark_all_reps_includced(self.conn)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Tests for generate_batch_rmd')
    parser.add_argument('-r', '--render', action='store_true', default=False,
                        help='Also test if rmd file can be rendered?')
    parser.add_argument('unittest_args', nargs='*')
    args = parser.parse_args()

    TestMakeQCrmd.RENDER_RMD = args.render

    unittest_args = args.unittest_args
    unittest_args.insert(0, sys.argv[0])
    unittest.main(argv=unittest_args, verbosity=2)
