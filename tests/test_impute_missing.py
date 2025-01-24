
import os
import re
import unittest
from unittest.mock import patch
import sqlite3
from io import StringIO

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils
import DIA_QC_report.submodules.imputation as imputation
import DIA_QC_report.impute_missing as impute_missing


class TestGetManager(unittest.TestCase):
    ERROR_PATTERNS = {'value':  r"Cannot convert value <[a-z]+> to <class '[a-z]+'>",
                      'choice': r"Value '[a-z]+' must be in \('[a-z]",
                      'range':  r'Value must be [<>]=?',
                      'unknown_arg': r"Unknown option: '[a-z_]+'",
                      'unparsable': r"Could not parse option argument: '[a-z_]+'"}


    def do_valid_arg_test(self, method, method_args, method_kwargs, all_kwargs):
        manager, metadata = impute_missing.get_manager(method, method_args, **method_kwargs)

        for key, value in all_kwargs.items():
            self.assertTrue(key in metadata)
            self.assertEqual(value, metadata[key])
            self.assertTrue(hasattr(manager, key))
            self.assertEqual(getattr(manager, key), value)


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
        manager, metadata = impute_missing.get_manager('KNN', None)

        # metadata should have at least as many entrires as the number of
        # KNNImputer.__init__ keyword arguments
        args = dict(setup_functions.function_kwargs(imputation.KNNImputer.__init__))
        args.pop('conn')
        self.assertEqual(len(metadata), len(args))

        default_manager = imputation.KNNImputer()

        for k, v in metadata.items():
            self.assertTrue(hasattr(default_manager, k))
            self.assertTrue(hasattr(manager, k))
            self.assertEqual(getattr(default_manager, k), v)
            self.assertEqual(getattr(manager, k), v)


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
                 {'impute_data': 'precursors', 'missing_threshold': 0.6},
                 {'n_neighbors': 8, 'weights': 'uniform'}],
                [[], {'impute_data': 'proteins', 'missing_threshold': 0.75}, {}]]

        for method_args, method_kwargs, target_args in args:
            all_args = method_kwargs.copy()
            all_args.update(target_args)

            self.do_valid_arg_test('KNN', method_args, method_kwargs, all_args)


    def test_KNN_inavlid_kwargs(self):
        args = [{'dummy': 0}, {'level': 0, 'another_dummy': 1}]

        for arg in args:
            self.do_invalid_kwarg_test('KNN', arg)
