
import os
import re
import unittest
from unittest.mock import patch
import sqlite3
from io import StringIO
from itertools import product

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

        self.assertEqual(error_cm.exception.code, 1)
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


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


class TestSingleImputation(unittest.TestCase, CommonTests):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_impute_missing_single/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        parse_command = ['dia_qc', 'parse', '-n=Sp3',
                         f'{cls.data_dir}/skyline_reports/Sp3_replicate_quality.tsv',
                         f'{cls.data_dir}/skyline_reports/Sp3_DiaNN_precursor_quality.tsv']
        setup_functions.make_work_dir(cls.work_dir, clear_dir=True)
        parse_result = setup_functions.run_command(parse_command, cls.work_dir, prefix='parse')

        normalize_command = ['dia_qc', 'normalize', '-m=median', '--keepMissing', cls.db_path]
        normalize_result = setup_functions.run_command(normalize_command,
                                                       cls.work_dir,
                                                       prefix='normalize')

        if sum((parse_result.returncode, normalize_result.returncode)) > 0:
            raise RuntimeError('Error creating test DB!')

        if os.path.isfile(cls.db_path):
            cls.conn = sqlite3.connect(cls.db_path)


    def test_imputation(self):
        self.assertIsNotNone(self.conn)
        self.assertTrue(db_utils.is_normalized(self.conn))

        impute_command = ['dia_qc', 'impute', self.db_path]
        impute_result = setup_functions.run_command(impute_command, self.work_dir)

        self.assertEqual(impute_result.returncode, 0)

        # test that metadata table is updated

        # test that missing values meeting imputation criteria are updated

        # test that non-missing values are not modified.


class TestMultiImputation(unittest.TestCase, CommonTests):
    pass

