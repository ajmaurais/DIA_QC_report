
import unittest
import argparse
import sys
import os

import setup_functions


TEST_DIR = os.path.dirname(os.path.abspath(__file__))

class TestMakeQCqmd(unittest.TestCase):
    TEST_PROJECT = 'Strap'
    RENDER_QMD = False

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_generate_qc_qmd/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{TEST_DIR}/data/'

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           clear_dir=True,
                                                           group_by_gene=False)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')


    def test_is_successful(self):
        self.assertEqual(self.parse_result.returncode, 0)

        command = ['generate_qc_qmd',
                   '-a', 'iRT', '-a', 'sp|P00924|ENO1_YEAST',
                   self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 0)

        if self.RENDER_QMD:
            render_command = ['quarto', 'render', 'qc_report.qmd', '--to', 'html']
            render_result = setup_functions.run_command(render_command, self.work_dir)
            self.assertEqual(render_result.returncode, 0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Tests for generate_qc_qmd')
    parser.add_argument('-r', '--render', action='store_true', default=False,
                        help='Also test if qmd file can be rendered?')
    parser.add_argument('unittest_args', nargs='*')
    args = parser.parse_args()

    TestMakeQCqmd.RENDER_QMD = args.render

    sys.argv[1:] = args.unittest_args
    unittest.main()
