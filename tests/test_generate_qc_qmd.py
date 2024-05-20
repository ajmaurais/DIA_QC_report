
import unittest
import argparse
import sys
import os

import setup_functions


class TestMakeQCqmd(unittest.TestCase):
    TEST_PROJECT = 'Strap'
    RENDER_QMD = False

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_generate_qc_qmd/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           clear_dir=True,
                                                           group_by_gene=False)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')


    def test_is_successful(self):
        self.assertEqual(self.parse_result.returncode, 0)

        qmd_name = 'basic_test'
        command = ['generate_qc_qmd',
                   '-a', 'iRT', '-a', 'sp|P00924|ENO1_YEAST',
                   '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))

        if self.RENDER_QMD:
            render_command = ['quarto', 'render', f'{qmd_name}.qmd', '--to', 'html']
            render_result = setup_functions.run_command(render_command, self.work_dir)
            self.assertEqual(render_result.returncode, 0)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.html'))


    def test_missing_std_protein_fails(self):
        self.assertEqual(self.parse_result.returncode, 0)

        qmd_name = 'failing_test'
        command = ['generate_qc_qmd', '-a', 'NOT_A_PROTEIN',
                   '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 1)
        self.assertTrue('ERROR: Missing standard protein: "NOT_A_PROTEIN"' in result.stderr)
        self.assertFalse(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))


    def test_missing_color_var_fails(self):
        qmd_name = 'failing_test'
        command = ['generate_qc_qmd', '-c', 'NOT_A_VAR',
                   '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 1)
        self.assertTrue('Missing annotationKey: "NOT_A_VAR"' in result.stderr)
        self.assertFalse(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Tests for generate_qc_qmd')
    parser.add_argument('-r', '--render', action='store_true', default=False,
                        help='Also test if qmd file can be rendered?')
    parser.add_argument('unittest_args', nargs='*')
    args = parser.parse_args()

    TestMakeQCqmd.RENDER_QMD = args.render

    unittest_args = args.unittest_args
    unittest_args.insert(0, sys.argv[0])
    unittest.main(argv=unittest_args, verbosity=2)
