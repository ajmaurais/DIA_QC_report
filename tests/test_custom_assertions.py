
import unittest
from numpy import nan

from setup_functions import AbstractTestsBase

class TestCustomAssertion(unittest.TestCase, AbstractTestsBase):
    def setUp(self):
        self.gt_data = {'rep_1': {'conc': 1.0, 'type': 'pool', 'rep': 1, 'is_std': True, 'na_var': None},
                        'rep_2': {'conc': 0.5, 'type': 'test', 'rep': 2, 'is_std': False, 'na_var': None},
                        'rep_3': {'conc': 0.25, 'type': 'test', 'rep': 3, 'is_std': False, 'na_var': None}}


    @unittest.expectedFailure
    def test_invalid_dtype_fails(self):
        self.assertDataDictEqual(self.gt_data, ['dummy'])


    @unittest.expectedFailure
    def test_missing_rep_fails(self):
        test_data = {k: v for i, (k, v) in enumerate(self.gt_data.items()) if i != 0}
        self.assertDataDictEqual(self.gt_data, test_data)


    @unittest.expectedFailure
    def test_missing_var_fails(self):
        test_data = dict()
        for i, rep in enumerate(self.gt_data):
            test_data[rep] = dict()
            for j, (key, value) in enumerate(self.gt_data[rep].items()):
                if i != 0 and j != 0:
                    test_data[rep][key] = value

        self.assertDataDictEqual(self.gt_data, test_data)


    @unittest.expectedFailure
    def test_different_var_type_fails(self):
        test_data = {rep: {k: str(v) if isinstance(v, float) else v for k, v in row.items()}
                        for rep, row in self.gt_data.items()}
        self.assertDataDictEqual(self.gt_data, test_data)


    def test_both_nan_equal(self):
        self.gt_data['rep_1']['conc'] = nan
        self.assertDataDictEqual(self.gt_data, self.gt_data)


    @unittest.expectedFailure
    def test_not_almost_equal_float_fails(self):
        test_data = {rep: {k: v + 1 if isinstance(v, float) else v for k, v in row.items()}
                        for rep, row in self.gt_data.items()}
        self.assertDataDictEqual(self.gt_data, test_data)


    @unittest.expectedFailure
    def test_not_almost_equal_col_delta_fails(self):
        test_data = {rep: {k: v + 1 if isinstance(v, float) else v for k, v in row.items()}
                        for rep, row in self.gt_data.items()}
        self.assertDataDictEqual(self.gt_data, test_data, col_deltas={'conc': 0.001})


    def test_almost_equal_float_succedes(self):
        test_data = {rep: {k: v + 1e-7 if isinstance(v, float) else v for k, v in row.items()}
                        for rep, row in self.gt_data.items()}
        self.assertDataDictEqual(self.gt_data, test_data)


    def test_almost_equal_col_delta_succedes(self):
        test_data = {rep: {k: v + 1e-7 if isinstance(v, float) else v for k, v in row.items()}
                        for rep, row in self.gt_data.items()}
        self.assertDataDictEqual(self.gt_data, test_data, col_deltas={'conc': 0.001})


    @unittest.expectedFailure
    def test_not_equal_int_fails(self):
        test_data = {rep: {k: v + 1 if isinstance(v, int) else v for k, v in row.items()}
                        for rep, row in self.gt_data.items()}
        self.assertDataDictEqual(self.gt_data, test_data)


