
import os
import re
import unittest
from unittest.mock import patch
import sqlite3
from io import StringIO
from itertools import product

import pandas as pd

import setup_functions

from DIA_QC_report.submodules import dia_db_utils as db_utils
from DIA_QC_report.submodules import imputation
from DIA_QC_report import impute_missing


class TestGetManager(unittest.TestCase):
    ERROR_PATTERNS = {'value':  r"Cannot convert value <[a-z]+> to <class '[a-z]+'>",
                      'choice': r"Value '[a-z]+' must be in \('[a-z]",
                      'range':  r'Value must be [<>]=?',
                      'unknown_arg': r"Unknown option: '[a-z_]+'",
                      'unparsable': r"Could not parse option argument: '[a-z_]+'"}


    def do_valid_arg_test(self, method, method_args, method_kwargs, all_kwargs):
        manager, metadata = impute_missing.get_manager(method, method_args, **method_kwargs)

        impute_data = {p: True if f'impute_{p}s' not in all_kwargs else all_kwargs.pop(f'impute_{p}s')
                       for p in ('protein', 'precursor')}

        for level, impute in impute_data.items():
            if impute:
                for key, value in all_kwargs.items():
                    metadata_key = f'{level}_imputation_{key}'
                    self.assertIn(metadata_key, metadata)
                    self.assertEqual(value, metadata[metadata_key])

                    self.assertTrue(hasattr(manager, key))
                    self.assertEqual(getattr(manager, key), value)
                    self.assertTrue(getattr(manager, f'impute_{level}s'))

            else:
                self.assertFalse(getattr(manager, f'impute_{level}s'))


    def do_invalid_arg_test(self, method, method_args, error_re):
        with self.assertLogs(impute_missing.LOGGER, level='ERROR') as log_cm:
            with self.assertRaises(SystemExit) as error_cm:
                impute_missing.get_manager(method, method_args)

        self.assertEqual(error_cm.exception.code, 2)
        self.assertTrue(any(re.search(error_re, m) for m in log_cm.output))


    def do_invalid_kwarg_test(self, method, kwargs):
        with self.assertRaises(TypeError) as cm:
            impute_missing.get_manager(method, [], **kwargs)

        error_msg = 'ImputationManagerBase.__init__()'
        self.assertTrue(error_msg in str(cm.exception))


    def test_inavlid_kwargs(self):
        args = [{'dummy': 0}, {'level': 0, 'another_dummy': 1}]

        for arg in args:
            self.do_invalid_kwarg_test('KNN', arg)


    def test_default_args(self):
        '''
        Test that impute_missing.get_manager returns manager with default arguments
        when method_args is None.
        '''
        method = 'KNN'
        manager, metadata = impute_missing.get_manager(method, None)

        # metadata should have at least as many entrires for precursors and proteins
        # as the number of KNNImputer.__init__ keyword arguments
        args = dict(setup_functions.function_kwargs(imputation.KNNImputer.__init__))
        args.pop('conn')
        target_metadata = [db_utils.PRECURSOR_IMPUTE_METHOD, db_utils.PROTEIN_IMPUTE_METHOD]
        for level, value in product(['precursor', 'protein'], args):
            target_metadata.append(f'{level}_imputation_{value}')

        self.assertEqual(len(metadata), len(target_metadata))

        default_manager = imputation.KNNImputer()

        self.assertEqual(metadata.pop(db_utils.PRECURSOR_IMPUTE_METHOD), method)
        self.assertEqual(metadata.pop(db_utils.PROTEIN_IMPUTE_METHOD), method)

        for k, v in metadata.items():
            key_arg_name = re.sub(r'^(precursor|protein)_imputation_', '', k)
            self.assertTrue(hasattr(default_manager, key_arg_name))
            self.assertTrue(hasattr(manager, key_arg_name))
            self.assertEqual(getattr(default_manager, key_arg_name), v)
            self.assertEqual(getattr(manager, key_arg_name), v)


    def test_invalid_method(self):
        with self.assertRaises(ValueError) as cm:
            impute_missing.get_manager('Dummy', [], method_help=True)

        self.assertTrue("Unknown imputation method: 'Dummy'" in str(cm.exception))


    def test_method_help(self):
        with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
            with self.assertRaises(SystemExit) as cm:
                impute_missing.get_manager('KNN', [], method_help=True)

        self.assertEqual(cm.exception.code, 0)
        mock_stdout.seek(0)
        self.assertTrue('Options for KNNImputer' in mock_stdout.read())


    def test_KNN_invalid_args(self):
        errors = [('choice', 'weights=none'),
                  ('range', 'n_neighbors=0'),
                  ('value', 'n_neighbors=wdummy'),
                  ('unknown_arg', 'not_a_argument=0'),
                  ('unparsable', 'weights')]

        for err_type, arg in errors:
            self.do_invalid_arg_test('KNN', [arg], self.ERROR_PATTERNS[err_type])


    def test_KNN_valid_args(self):
        args = [[['n_neighbors=6', "weights='uniform'"], {},
                 {'n_neighbors': 6, 'weights': 'uniform'}],
                [['n_neighbors=8', "weights='uniform'"],
                 {'impute_proteins': False, 'missing_threshold': 0.6},
                 {'n_neighbors': 8, 'weights': 'uniform'}],
                [[], {'impute_precursors': False, 'missing_threshold': 0.75}, {}]]

        for method_args, method_kwargs, target_args in args:
            all_args = method_kwargs.copy()
            all_args.update(target_args)

            self.do_valid_arg_test('KNN', method_args, method_kwargs, all_args)


class CommonTests(setup_functions.AbstractTestsBase):
    def __init__(self):
        self.work_dir = None
        self.db_path = None
        self.data_dir = None
        self.conn = None


    def get_not_imputed(self):
        df_prec = pd.read_sql('''
            SELECT
                r.project, p.replicateId, p.peptideId, p.precursorCharge,
                p.totalAreaFragment, p.normalizedArea
            FROM precursors p
            LEFT JOIN replicates r ON r.id == p.replicateId
            WHERE isImputed = 0 and
                  totalAreaFragment is not NULL and
                  normalizedArea is not NULL; ''', self.conn)

        df_prot = pd.read_sql('''
            SELECT
                r.project, q.replicateId, q.proteinId,
                q.abundance, q.normalizedAbundance
            FROM proteinQuants q
            LEFT JOIN replicates r ON r.id == q.replicateId
            WHERE isImputed = 0 and
                  abundance is not NULL and
                  normalizedAbundance is not NULL; ''', self.conn)

        return df_prec, df_prot


    def get_all(self):
        df_prec = pd.read_sql('''
            SELECT
                r.project, p.replicateId, p.peptideId, p.precursorCharge,
                p.isImputed, p.totalAreaFragment as quant, p.normalizedArea as normQuant
            FROM precursors p
            LEFT JOIN replicates r ON r.id == p.replicateId; ''', self.conn)

        df_prec['id'] = df_prec['peptideId'].astype(str) + '_' + df_prec['precursorCharge'].astype(str)
        df_prec = df_prec.drop(columns=['peptideId', 'precursorCharge'])

        df_prot = pd.read_sql('''
            SELECT
                r.project, q.replicateId, q.proteinId as id,
                q.isImputed, q.abundance as quant, q.normalizedAbundance as normQuant
            FROM proteinQuants q
            LEFT JOIN replicates r ON r.id == q.replicateId ; ''', self.conn)

        return df_prec, df_prot


    def get_quant_level(self, normalized=False):
        '''
        Get all the precursor and protein quantities which are unnormalized or normalized.

        Parameters
        ----------
        normalized: bool
            Should the quant value be normalized? Default is False.
        '''
        df_prec = pd.read_sql(f'''
            SELECT
                r.project, p.replicateId, p.peptideId, p.precursorCharge,
                p.{'normalizedArea' if normalized else 'totalAreaFragment'}
            FROM precursors p
            LEFT JOIN replicates r ON r.id == p.replicateId; ''', self.conn)

        df_prot = pd.read_sql(f'''
            SELECT
                r.project, q.replicateId, q.proteinId,
                q.{'normalizedAbundance' if normalized else 'abundance'}
            FROM proteinQuants q
            LEFT JOIN replicates r ON r.id == q.replicateId; ''', self.conn)

        return df_prec, df_prot


    def do_imputation_criteria_test(self, df,
                                    group_by_project=True,
                                    missing_threshold=0.5,
                                    normalized=False):
        cur = self.conn.cursor()
        cur.execute('SELECT project, COUNT(id) FROM replicates GROUP BY project;')
        projects = {row[0]: row[1] for row in cur.fetchall()}
        value_col = 'normQuant' if normalized else 'quant'

        def group_na_count(row):
            if any(row['isImputed']):
                return row[value_col].isna().sum()
            return 0

        for project, n_reps in projects.items():
            project_df = df[df['project'] == project]
            na_counts = project_df.groupby('id').apply(group_na_count, include_groups=False)
            self.assertTrue(all(na_counts / n_reps < missing_threshold))


    def do_imputed_value_count_test(self, expect_prot=False, expect_prec=False):
        cur = self.conn.cursor()

        cur.execute('SELECT isImputed, COUNT(isImputed) FROM precursors GROUP BY isImputed;')
        prec_imputed = dict(cur.fetchall())

        cur.execute('SELECT isImputed, COUNT(isImputed) FROM proteinQuants GROUP BY isImputed;')
        prot_imputed = dict(cur.fetchall())

        if expect_prec:
            self.assertGreater(prec_imputed.get(1, 0), 0)
        else:
            self.assertEqual(prec_imputed.get(1, 0), 0)
        if expect_prot:
            self.assertGreater(prot_imputed.get(1, 0), 0)
        else:
            self.assertEqual(prot_imputed.get(1, 0), 0)
        self.assertGreater(prec_imputed.get(0, 0), 0)
        self.assertGreater(prot_imputed.get(0, 0), 0)


class TestSingleImputation(unittest.TestCase, CommonTests):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_impute_missing_single/'
        cls.base_db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'
        cls.conn = None

        parse_command = ['dia_qc', 'parse', '-n=Sp3',
                         f'{cls.data_dir}/skyline_reports/Sp3_replicate_quality.tsv',
                         f'{cls.data_dir}/skyline_reports/Sp3_DiaNN_precursor_quality.tsv']
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        parse_result = setup_functions.run_command(parse_command, cls.work_dir, prefix='parse')

        if parse_result.returncode != 0:
            raise RuntimeError('Error creating test DB!')


    def setUp(self):
        if not os.path.isfile(self.base_db_path):
            raise RuntimeError('Database does not exist!')

        # copy test db to unique db for test
        self.db_path = f'{self.work_dir}/{self._testMethodName}.db3'
        conn = sqlite3.connect(self.base_db_path)
        self.conn = sqlite3.connect(self.db_path)
        conn.backup(self.conn)
        conn.close()


    def tearDown(self):
        if self.conn is not None:
            self.conn.close()


    def do_imputation_test(self, command,
                           impute_precursors=True, impute_proteins=True,
                           weights='uniform', missing_threshold=0.5,
                           level=0):

        self.assertIsNotNone(self.conn)
        self.assertFalse(db_utils.any_imputed(self.conn))
        self.do_imputed_value_count_test()

        # test that initially there is no imputation metadata
        cur = self.conn.cursor()
        cur.execute('SELECT key, value FROM metadata;')
        metadata_i = dict(cur.fetchall())
        self.assertNotIn(db_utils.PRECURSOR_IMPUTE_METHOD, metadata_i)
        self.assertNotIn(db_utils.PROTEIN_IMPUTE_METHOD, metadata_i)

        prec_lhs, prot_lhs = self.get_not_imputed()
        prec_quant_lhs, prot_quant_lhs = self.get_quant_level(normalized=level==0)

        impute_result = setup_functions.run_command(command, self.work_dir)

        # test that imputation was sucessful and there are now 1 or more imputed values
        self.assertEqual(impute_result.returncode, 0)
        self.do_imputed_value_count_test(expect_prot=impute_proteins,
                                         expect_prec=impute_precursors)

        # test that metadata table is updated
        cur = self.conn.cursor()
        cur.execute('SELECT key, value FROM metadata;')
        metadata = dict(cur.fetchall())

        if impute_precursors:
            self.assertIn(db_utils.PRECURSOR_IMPUTE_METHOD, metadata)
            self.assertIn('precursor_imputation_weights', metadata)
            self.assertIn('precursor_imputation_level', metadata)
            self.assertEqual(metadata['precursor_imputation_level'], str(level))
            self.assertEqual(metadata['precursor_imputation_weights'], weights)
        else:
            self.assertNotIn(db_utils.PRECURSOR_IMPUTE_METHOD, metadata)
            self.assertNotIn('precursor_imputation_weights', metadata)
            self.assertNotIn('precursor_imputation_level', metadata)

        if impute_proteins:
            self.assertIn(db_utils.PROTEIN_IMPUTE_METHOD, metadata)
            self.assertIn('protein_imputation_weights', metadata)
            self.assertIn('protein_imputation_level', metadata)
            self.assertEqual(metadata['protein_imputation_level'], str(level))
            self.assertEqual(metadata['protein_imputation_weights'], weights)
        else:
            self.assertNotIn(db_utils.PROTEIN_IMPUTE_METHOD, metadata)
            self.assertNotIn('protein_imputation_weights', metadata)
            self.assertNotIn('protein_imputation_level', metadata)

        # test that missing values meeting imputation criteria are updated
        prec_i, prot_i = self.get_all()
        if impute_precursors:
            self.do_imputation_criteria_test(prec_i, normalized=level==1,
                                             missing_threshold=missing_threshold)
        if impute_proteins:
            self.do_imputation_criteria_test(prot_i, normalized=level==1,
                                             missing_threshold=missing_threshold)

        # test that non-missing values are not modified.
        prec_rhs, prot_rhs = self.get_not_imputed()
        self.assertDataFrameEqual(prec_lhs, prec_rhs)
        self.assertDataFrameEqual(prot_lhs, prot_rhs)

        # test that opposite level (unnormalized or normalized) are not modified
        prec_quant_rhs, prot_quant_rhs = self.get_quant_level(normalized=level==0)
        self.assertDataFrameEqual(prec_quant_lhs, prec_quant_rhs)
        self.assertDataFrameEqual(prot_quant_lhs, prot_quant_rhs)


    def test_weights(self):
        impute_command = ['dia_qc', 'impute', '-t=0.5', '-s', 'weights=weighted', self.db_path]
        self.do_imputation_test(impute_command, weights='weighted')


    def test_normalized(self):
        normalize_command = ['dia_qc', 'normalize', '-m=median', '--keepMissing', self.db_path]
        normalize_result = setup_functions.run_command(normalize_command,
                                                       self.work_dir,
                                                       prefix='normalize')
        self.assertEqual(normalize_result.returncode, 0)
        self.assertTrue(db_utils.is_normalized(self.conn))

        impute_command = ['dia_qc', 'impute', '-t=0.5', '-l=1', self.db_path]
        self.do_imputation_test(impute_command, level=1)


    def test_na_threshold(self):
        impute_command = ['dia_qc', 'impute', '-t=1', self.db_path]
        self.do_imputation_test(impute_command, missing_threshold=1)


    def test_protein_only(self):
        impute_command = ['dia_qc', 'impute', '-i=2', self.db_path]
        self.do_imputation_test(impute_command,
                                impute_precursors=False,
                                impute_proteins=True)


    def test_precursor_only(self):
        impute_command = ['dia_qc', 'impute', '-i=1', self.db_path]
        self.do_imputation_test(impute_command,
                                impute_precursors=True,
                                impute_proteins=False)


    def test_method_help(self):
        command = ['dia_qc', 'impute', '-m=KNN', '--methodHelp']
        help_result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(help_result.returncode, 0)
        self.assertIn('Options for KNNImputer', help_result.stdout)


    def test_invalid_method(self):
        command = ['dia_qc', 'impute', '-m=dummy', self.db_path]
        error_result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(error_result.returncode, 2)
        message = "invalid choice: 'dummy' (choose from '{}')".format("', '".join(imputation.IMPUTATION_METHODS))
        self.assertIn(message, error_result.stderr)


class TestMultiImputation(unittest.TestCase, CommonTests):
    def test_group_by(self):
        pass
