
import unittest
import argparse
import sys
import os
import sqlite3
import random
import pandas as pd
from numpy import isnan

import setup_functions

from DIA_QC_report.submodules.pca_plot import convert_string_cols


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
        command = ['dia_qc', 'qc_qmd',
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
        command = ['dia_qc', 'qc_qmd', '-a', 'NOT_A_PROTEIN',
                   '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 1)
        self.assertTrue('ERROR: Missing standard protein: "NOT_A_PROTEIN"' in result.stderr)
        self.assertFalse(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))


    def test_missing_color_var_fails(self):
        qmd_name = 'failing_test'
        command = ['dia_qc', 'qc_qmd', '-c', 'NOT_A_VAR',
                   '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 1)
        self.assertTrue('Missing annotationKey: "NOT_A_VAR"' in result.stderr)
        self.assertFalse(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))


class TestMissingMetadata(unittest.TestCase):
    TEST_PROJECT = 'Strap'
    RENDER_QMD = False

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_qc_report_missing_metadata/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           metadata_suffix='_missing_multi_var_metadata.tsv',
                                                           clear_dir=True)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        if os.path.isfile(cls.db_path):
            cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def do_query(self, keys):
        query = '''SELECT replicateId,
                 m.annotationKey as key,
                 m.annotationValue as value,
                 t.annotationType as type
             FROM sampleMetadata m
             LEFT JOIN sampleMetadataTypes t ON t.annotationKey == m.annotationKey
             WHERE m.annotationKey IN ("{}")'''.format('", "'.join(keys))

        return pd.read_sql(query, self.conn)


    def convert_strings_test(self, var):
        self.assertIsNotNone(self.conn)
        df = self.do_query([var])
        na_reps = df.loc[(df['key'] == var) & pd.isna(df['value']), 'replicateId'].to_list()
        metadata = convert_string_cols(df)
        metadata = metadata.set_index('replicateId')

        for rep in na_reps:
            self.assertTrue(isnan(metadata.loc[rep, var]))


    def test_missing_float(self):
        self.convert_strings_test('float_var')


    def test_missing_string(self):
        self.convert_strings_test('string_var')


    def test_missing_int(self):
        self.convert_strings_test('int_var')


    def test_missing_bool(self):
        self.convert_strings_test('bool_var')


    def test_is_successful(self):
        self.assertEqual(self.parse_result.returncode, 0)

        qmd_name = 'basic_test'
        command = ['dia_qc', 'qc_qmd',
                   '-a', 'iRT', '-a', 'sp|P00924|ENO1_YEAST',
                   '-c=string_var', '-c=bool_var', '-c=int_var', '-c=float_var',
                   '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))

        if self.RENDER_QMD:
            render_command = ['quarto', 'render', f'{qmd_name}.qmd', '--to', 'html']
            render_result = setup_functions.run_command(render_command, self.work_dir)
            self.assertEqual(render_result.returncode, 0)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.html'))


class TestBadMetadataHeaders(unittest.TestCase):
    TEST_PROJECT = 'Strap'
    RENDER_QMD = False

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_qc_report_bad_metadata_headers/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           metadata_suffix='_bad_headers_metadata.tsv',
                                                           clear_dir=True)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        if os.path.isfile(cls.db_path):
            cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_is_successful(self):
        self.assertEqual(self.parse_result.returncode, 0)

        qmd_name = 'bad_header_test'
        command = ['dia_qc', 'qc_qmd',
                   '-a', 'iRT', '-a', 'sp|P00924|ENO1_YEAST',
                   '-c', 'string var', '-c', 'bool var', '-c', 'int var', '-c', 'float var',
                   '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))

        if self.RENDER_QMD:
            render_command = ['quarto', 'render', f'{qmd_name}.qmd', '--to', 'html']
            render_result = setup_functions.run_command(render_command, self.work_dir)
            self.assertEqual(render_result.returncode, 0)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.html'))


class TestAllPrecursorsMissing(unittest.TestCase):
    TEST_PROJECT = 'GFP'
    RENDER_QMD = False

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_qc_report_missing_precursors/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        command = ['dia_qc', 'parse', f'--projectName={cls.TEST_PROJECT}',
                   f'{cls.data_dir}/skyline_reports/{cls.TEST_PROJECT}_replicate_quality.tsv',
                   f'{cls.data_dir}/skyline_reports/{cls.TEST_PROJECT}_precursor_quality.tsv']

        setup_functions.make_work_dir(cls.work_dir, True)
        cls.parse_result = setup_functions.run_command(command, cls.work_dir,
                                                       prefix='parse_missing_precursors')

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        if os.path.isfile(cls.db_path):
            cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_is_successful(self):
        self.assertEqual(self.parse_result.returncode, 0)

        # Generate unnormalized qmd
        unorm_qmd_name = 'unnormalized_test'
        generate_qmd_command = ['dia_qc', 'qc_qmd', '-o', f'{unorm_qmd_name}.qmd', self.db_path]
        unorm_result = setup_functions.run_command(generate_qmd_command, self.work_dir,
                                                   prefix=unorm_qmd_name)

        self.assertEqual(unorm_result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{unorm_qmd_name}.qmd'))

        if self.RENDER_QMD:
            render_command = ['quarto', 'render', f'{unorm_qmd_name}.qmd', '--to', 'html']
            render_result = setup_functions.run_command(render_command, self.work_dir)
            self.assertEqual(render_result.returncode, 0)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{unorm_qmd_name}.html'))

        # Normalize database
        normalize_command = ['dia_qc', 'normalize', '-m=median', self.db_path]
        norm_db_result = setup_functions.run_command(normalize_command, self.work_dir,
                                                  prefix='normalize')
        self.assertEqual(norm_db_result.returncode, 0)

        # Generate normalized qmd
        norm_db_result = setup_functions.run_command(generate_qmd_command, self.work_dir,
                                                     prefix=unorm_qmd_name)
        self.assertEqual(norm_db_result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{unorm_qmd_name}.qmd'))

        if self.RENDER_QMD:
            render_command = ['quarto', 'render', f'{unorm_qmd_name}.qmd', '--to', 'html']
            render_result = setup_functions.run_command(render_command, self.work_dir)
            self.assertEqual(render_result.returncode, 0)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{unorm_qmd_name}.html'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Tests for generate_qc_qmd')
    parser.add_argument('-r', '--render', action='store_true', default=False,
                        help='Also test if qmd file can be rendered?')
    parser.add_argument('unittest_args', nargs='*')
    args = parser.parse_args()

    TestMakeQCqmd.RENDER_QMD = args.render
    TestMissingMetadata.RENDER_QMD = args.render
    TestAllPrecursorsMissing.RENDER_QMD = args.render

    unittest_args = args.unittest_args
    unittest_args.insert(0, sys.argv[0])
    unittest.main(argv=unittest_args, verbosity=2)
