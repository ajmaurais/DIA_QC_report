
import unittest
import os
import sqlite3
import pandas as pd
from numpy import isnan

import setup_functions

from DIA_QC_report import parse_data
from DIA_QC_report import normalize_db
from DIA_QC_report.submodules.dia_db_utils import read_wide_metadata
from DIA_QC_report import generate_qc_qmd


class TestMakeQCqmd(unittest.TestCase):
    TEST_PROJECT = 'Strap'

    @classmethod
    def setUpClass(cls):
        cls.render_qmd = os.getenv('RENDER_QMD', 'False').lower() == 'true'

        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_generate_qc_qmd/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        setup_functions.remove_test_dir(f'{cls.work_dir}/pdf_tab_basic_test_files', recursive=True)
        cls.parse_result = setup_functions.setup_single_db(
            cls.data_dir, cls.work_dir, cls.TEST_PROJECT,
            clear_dir=True, group_by_gene=False, subprocess=True
        )

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')


    def render_test(self, report_format, pca_format='tab'):
        self.assertEqual(self.parse_result.returncode, 0)

        qmd_name = f'{report_format}_{pca_format}_basic_test'
        command = ['-a', 'iRT', '-a', 'sp|P00924|ENO1_YEAST', '--pcaFormat', pca_format,
                   '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_main(
            generate_qc_qmd._main, command, self.work_dir,
            prefix=qmd_name, prog='dia_qc qc_qmd'
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))

        if self.render_qmd:
            render_command = ['quarto', 'render', f'{qmd_name}.qmd', '--to', report_format]
            render_result = setup_functions.run_command(
                render_command, self.work_dir, prefix=f'render_{qmd_name}'
            )
            self.assertEqual(render_result.returncode, 0, render_result.stderr)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.{report_format}'))


    def test_pdf(self):
        self.render_test('pdf')


    def test_html_stacked(self):
        self.render_test('html', 'stack')


    def test_html_tabbed(self):
        self.render_test('html', 'tab')


    def test_missing_std_protein_fails(self):
        self.assertEqual(self.parse_result.returncode, 0)

        qmd_name = 'failing_test'
        command = ['-a', 'NOT_A_PROTEIN', '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_main(
            generate_qc_qmd._main, command, self.work_dir, prog='dia_qc qc_qmd'
        )
        self.assertEqual(result.returncode, 1)
        self.assertTrue('Missing standard protein: "NOT_A_PROTEIN"' in result.stderr, result.stderr)
        self.assertFalse(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))


    def test_missing_color_var_fails(self):
        qmd_name = 'failing_test'
        command = ['-c', 'NOT_A_VAR', '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_main(
            generate_qc_qmd._main, command, self.work_dir, prog='dia_qc qc_qmd'
        )

        self.assertEqual(result.returncode, 1)
        self.assertTrue('Missing annotationKey: "NOT_A_VAR"' in result.stderr)
        self.assertFalse(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))


class TestMissingMetadata(unittest.TestCase):
    TEST_PROJECT = 'Strap'

    @classmethod
    def setUpClass(cls):
        cls.render_qmd = os.getenv('RENDER_QMD', 'False').lower() == 'true'

        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_qc_report_missing_metadata/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        setup_functions.remove_test_dir(f'{cls.work_dir}/pdf_tab_test_files', recursive=True)
        cls.parse_result = setup_functions.setup_single_db(
            cls.data_dir, cls.work_dir, cls.TEST_PROJECT,
            metadata_suffix='_missing_multi_var_metadata.tsv',
            clear_dir=True, subprocess=True
        )

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

        metadata = read_wide_metadata(self.conn)

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


    def render_test(self, report_format, pca_format='tab'):
        self.assertEqual(self.parse_result.returncode, 0)

        qmd_name = f'{report_format}_{pca_format}_test'
        command = ['-a', 'iRT', '-a', 'sp|P00924|ENO1_YEAST',
                   '-c=string_var', '-c=bool_var', '-c=int_var', '-c=float_var', '--pcaFormat', pca_format,
                   '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_main(
            generate_qc_qmd._main, command, self.work_dir,
            prefix=qmd_name, prog='dia_qc qc_qmd'
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))

        if self.render_qmd:
            render_command = ['quarto', 'render', f'{qmd_name}.qmd', '--to', report_format]
            render_result = setup_functions.run_command(
                render_command, self.work_dir, prefix=f'render_{qmd_name}'
            )
            self.assertEqual(render_result.returncode, 0, render_result.stderr)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.{report_format}'))


    def test_pdf(self):
        self.render_test('pdf')


    def test_html_stacked(self):
        self.render_test('html', 'stack')


    def test_html_tabbed(self):
        self.render_test('html', 'tab')


class TestBadMetadataHeaders(unittest.TestCase):
    TEST_PROJECT = 'Strap'

    @classmethod
    def setUpClass(cls):
        cls.render_qmd = os.getenv('RENDER_QMD', 'False').lower() == 'true'

        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_qc_report_bad_metadata_headers/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        setup_functions.remove_test_dir(f'{cls.work_dir}/bad_header_test_files', recursive=True)
        cls.parse_result = setup_functions.setup_single_db(
            cls.data_dir, cls.work_dir, cls.TEST_PROJECT,
            metadata_suffix='_bad_headers_metadata.tsv',
            clear_dir=True, subprocess=True
        )

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
        command = ['-a', 'iRT', '-a', 'sp|P00924|ENO1_YEAST',
                   '-c', 'string var', '-c', 'bool var', '-c', 'int var', '-c', 'float var',
                   '-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_main(
            generate_qc_qmd._main, command, self.work_dir,
            prefix=qmd_name, prog='dia_qc qc_qmd'
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))

        if self.render_qmd:
            render_command = ['quarto', 'render', f'{qmd_name}.qmd', '--to', 'html']
            render_result = setup_functions.run_command(
                render_command, self.work_dir, prefix='render_bad_header_qmd'
            )
            self.assertEqual(render_result.returncode, 0, render_result.stderr)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.html'))


class TestSingleReplicate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.render_qmd = os.getenv('RENDER_QMD', 'False').lower() == 'true'

        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_qc_report_single_replicate/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        parse_command = [f'{cls.data_dir}/invalid_reports/Sp3_single_replicate_replicate_quality.tsv',
                         f'{cls.data_dir}/invalid_reports/Sp3_by_protein_single_replicate_precursor_quality.tsv']

        setup_functions.make_work_dir(cls.work_dir, True)
        cls.parse_result = setup_functions.run_main(
            parse_data._main, parse_command, cls.work_dir, prog='dia_qc parse'
        )

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

        qmd_name = 'one_replicate_test'
        command = ['-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_main(
            generate_qc_qmd._main, command, self.work_dir,
            prefix='generate_qc_qmd', prog='dia_qc qc_qmd'
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))

        if self.render_qmd:
            render_command = ['quarto', 'render', f'{qmd_name}.qmd', '--to', 'html']
            render_result = setup_functions.run_command(
                render_command, self.work_dir, prefix=f'render_{qmd_name}'
            )
            self.assertEqual(render_result.returncode, 0, render_result.stderr)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.html'))


class TestSingleReplicate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.render_qmd = os.getenv('RENDER_QMD', 'False').lower() == 'true'

        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_qc_report_single_replicate/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        parse_command = [f'{cls.data_dir}/invalid_reports/Sp3_single_replicate_replicate_quality.tsv',
                         f'{cls.data_dir}/invalid_reports/Sp3_by_protein_single_replicate_precursor_quality.tsv']

        setup_functions.make_work_dir(cls.work_dir, True)
        cls.parse_result = setup_functions.run_main(
            parse_data._main, parse_command, cls.work_dir, prog='dia_qc parse'
        )
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

        qmd_name = 'one_replicate_test'
        command = ['-o', f'{qmd_name}.qmd', self.db_path]
        result = setup_functions.run_main(
            generate_qc_qmd._main, command, self.work_dir, prog='dia_qc qc_qmd'
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.qmd'))

        if self.render_qmd:
            render_command = ['quarto', 'render', f'{qmd_name}.qmd', '--to', 'html']
            render_result = setup_functions.run_command(render_command, self.work_dir)
            self.assertEqual(render_result.returncode, 0, render_result.stderr)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{qmd_name}.html'))


class TestAllPrecursorsMissing(unittest.TestCase):
    TEST_PROJECT = 'GPF'

    @classmethod
    def setUpClass(cls):
        cls.render_qmd = os.getenv('RENDER_QMD', 'False').lower() == 'true'

        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_qc_report_missing_precursors/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        command = [f'--projectName={cls.TEST_PROJECT}',
                   f'{cls.data_dir}/skyline_reports/{cls.TEST_PROJECT}_replicate_quality.tsv',
                   f'{cls.data_dir}/skyline_reports/{cls.TEST_PROJECT}_precursor_quality.tsv']

        setup_functions.remove_test_dir(f'{cls.work_dir}/unnormalized_test_files', recursive=True)
        setup_functions.make_work_dir(cls.work_dir, True)
        cls.parse_result = setup_functions.run_main(
            parse_data._main, command, cls.work_dir,
            prefix='parse_missing_precursors', prog='dia_qc parse'
        )

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
        prog = 'dia_qc qc_qmd'

        # Generate unnormalized qmd
        unorm_qmd_name = 'unnormalized_test'
        generate_qmd_command = ['-o', f'{unorm_qmd_name}.qmd', self.db_path]
        unorm_result = setup_functions.run_main(
            generate_qc_qmd._main, generate_qmd_command, self.work_dir,
            prefix=unorm_qmd_name, prog=prog
        )
        self.assertEqual(unorm_result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{unorm_qmd_name}.qmd'))

        if self.render_qmd:
            render_command = ['quarto', 'render', f'{unorm_qmd_name}.qmd', '--to', 'html']
            render_result = setup_functions.run_command(
                render_command, self.work_dir, prefix=f'render_{unorm_qmd_name}'
            )
            self.assertEqual(render_result.returncode, 0, render_result.stderr)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{unorm_qmd_name}.html'))

        # Normalize database
        norm_db_result = setup_functions.run_main(
            normalize_db._main, ['-m=median', self.db_path], self.work_dir,
            prefix='normalize', prog='dia_qc normalize'
        )
        self.assertEqual(norm_db_result.returncode, 0)

        # Generate normalized qmd
        norm_qmd_name = 'normalized_test'
        generate_qmd_command = ['-o', f'{norm_qmd_name}.qmd', self.db_path]
        norm_db_result = setup_functions.run_main(
            generate_qc_qmd._main, generate_qmd_command, self.work_dir,
            prefix=norm_qmd_name, prog=prog
        )
        self.assertEqual(norm_db_result.returncode, 0)
        self.assertTrue(os.path.isfile(f'{self.work_dir}/{norm_qmd_name}.qmd'))

        if self.render_qmd:
            render_command = ['quarto', 'render', f'{norm_qmd_name}.qmd', '--to', 'html']
            render_result = setup_functions.run_command(
                render_command, self.work_dir, prefix=f'render_{norm_qmd_name}'
            )
            self.assertEqual(render_result.returncode, 0, render_result.stderr)
            self.assertTrue(os.path.isfile(f'{self.work_dir}/{norm_qmd_name}.html'))


class TestProjectOption(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.render_qmd = os.getenv('RENDER_QMD', 'False').lower() == 'true'
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_qc_report_project_option'
        setup_results = setup_functions.setup_multi_db(f'{setup_functions.TEST_DIR}/data',
                                                       cls.work_dir, clear_dir=True, normalize=True)

        if not all(r.returncode == 0 for r in setup_results):
            raise RuntimeError('Setup of test db failed!')

        cls.db_path = f'{cls.work_dir}/data.db3'


    def test_missing_project_fails(self):
        dummy_project = 'NOT_A_PROJECT'
        command = [f'--project={dummy_project}', self.db_path]
        result = setup_functions.run_main(
            generate_qc_qmd._main, command, self.work_dir, prog='dia_qc qc_qmd'
        )
        self.assertEqual(result.returncode, 1)
        self.assertTrue(f"Project '{dummy_project}' does not exist!" in result.stderr,
                        result.stderr)


    def do_project_test(self, project):
        command = ['-a=iRT', f'--project={project}', self.db_path]
        result = setup_functions.run_main(
            generate_qc_qmd._main, command, self.work_dir,
            prefix=f'test_{project}', prog='dia_qc qc_qmd'
        )
        self.assertEqual(result.returncode, 0, result.stderr)

        qmd_name = f'{project}_{generate_qc_qmd.DEFAULT_OFNAME}'
        qmd_path = f'{self.work_dir}/{qmd_name}'
        self.assertTrue(os.path.isfile(qmd_path), f'{qmd_path} does not exist!')

        if self.render_qmd:
            render_command = ['quarto', 'render', qmd_path, '--to', 'html']
            render_result = setup_functions.run_command(
                render_command, self.work_dir, prefix=f'render_{qmd_name}'
            )
            self.assertEqual(render_result.returncode, 0, render_result.stderr)
            html_path = f'{self.work_dir}/{os.path.splitext(qmd_name)[0]}.html'
            self.assertTrue(os.path.isfile(html_path), f'{html_path} Does not exist!')


    def test_Sp3(self):
        self.do_project_test('Sp3')


    def test_Strap(self):
        self.do_project_test('Strap')