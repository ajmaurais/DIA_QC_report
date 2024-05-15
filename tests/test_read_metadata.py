
import unittest
from unittest import mock
import json
from jsonschema import validate, ValidationError

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils
from DIA_QC_report.submodules.metadata import Dtype


def df_to_dict(df, types):
    ret = dict()
    for row in df.itertuples():
        if row.Replicate not in ret:
            ret[row.Replicate] = {}

        new_value = types[row.annotationKey].convert(row.annotationValue)
        ret[row.Replicate][row.annotationKey] = new_value

    return ret


class TestReadMetadataBase(unittest.TestCase):

    META_TYPES = {'string_var': Dtype.STRING,
                  'bool_var': Dtype.BOOL,
                  'int_var': Dtype.INT,
                  'float_var': Dtype.FLOAT,
                  'na_var': Dtype.NULL}

    @classmethod
    def setUpClass(cls):
        cls.metadata_dir = f'{setup_functions.TEST_DIR}/data/metadata'

        with open(f'{cls.metadata_dir}/HeLa_metadata.json') as inF:
            cls.gt_data = json.load(inF)
            validate(cls.gt_data, db_utils.METADATA_SCHEMA)


    @staticmethod
    def remove_null_data(d):
        na_types = {k for k, v in TestParseSkylineAnnotations.ALL_META_TYPES.items() if v is Dtype.NULL}
        return {rep: {k: v for k, v in row.items() if k not in na_types} for rep, row in d.items()}


    @staticmethod
    def remove_null_types(d):
        return {k: v for k, v in d.items() if v is not Dtype.NULL}


    def assertDataDictEqual(self, lhs, rhs):
        lhs_keys = set(lhs.keys())
        rhs_keys = set(rhs.keys())

        self.assertEqual(lhs_keys, rhs_keys)

        for rep in lhs_keys:
            lhs_vars = set(lhs[rep].keys())
            rhs_vars = set(rhs[rep].keys())

            self.assertEqual(lhs_vars, rhs_vars)
            for var in lhs_vars:
                self.assertEqual(type(lhs[rep][var]), type(rhs[rep][var]))
                if isinstance(lhs[rep][var], float):
                    self.assertAlmostEqual(lhs[rep][var], rhs[rep][var], places=6)
                else:
                    self.assertEqual(lhs[rep][var], rhs[rep][var])


class TestCustomAssertion(TestReadMetadataBase):
    @unittest.expectedFailure
    def test_missing_rep_fails(self):
        test_data = {k: v for i, (k, v) in enumerate(self.gt_data.items()) if i != 0}
        self.assertDataDictEqual(self.gt_data, test_data)


    @unittest.expectedFailure
    def test_missing_var_fails(self):
        test_data = self.remove_null_data(self.gt_data)
        self.assertDataDictEqual(self.gt_data, test_data)


    @unittest.expectedFailure
    def test_different_var_type_fails(self):
        test_data = {rep: {k: str(v) if isinstance(v, float) else v for k, v in row.items()}
                        for rep, row in self.gt_data.items()}
        self.assertDataDictEqual(self.gt_data, test_data)


    @unittest.expectedFailure
    def test_not_almost_equal_float_fails(self):
        test_data = {rep: {k: v + 1 if isinstance(v, float) else v for k, v in row.items()}
                        for rep, row in self.gt_data.items()}
        self.assertDataDictEqual(self.gt_data, test_data)


    def test_almost_equal_float_succedes(self):
        test_data = {rep: {k: v + 1e-7 if isinstance(v, float) else v for k, v in row.items()}
                        for rep, row in self.gt_data.items()}
        self.assertDataDictEqual(self.gt_data, test_data)


    @unittest.expectedFailure
    def test_not_equal_int_fails(self):
        test_data = {rep: {k: v + 1 if isinstance(v, int) else v for k, v in row.items()}
                        for rep, row in self.gt_data.items()}
        self.assertDataDictEqual(self.gt_data, test_data)


class TestParseMetadata(TestReadMetadataBase, unittest.TestCase):
    def test_json(self):
        data, types = db_utils.read_metadata(f'{self.metadata_dir}/HeLa_metadata.json')

        test_df = df_to_dict(data, self.META_TYPES)

        self.assertDataDictEqual(self.gt_data, test_df)
        self.assertDictEqual(self.META_TYPES, types)


    def test_csv(self):
        data, types = db_utils.read_metadata(f'{self.metadata_dir}/HeLa_metadata.csv')

        test_df = df_to_dict(data, self.META_TYPES)

        self.assertDataDictEqual(self.gt_data, test_df)
        self.assertDictEqual(self.META_TYPES, types)


    def test_tsv(self):
        data, types = db_utils.read_metadata(f'{self.metadata_dir}/HeLa_metadata.tsv')

        test_df = df_to_dict(data, self.META_TYPES)

        self.assertDataDictEqual(self.gt_data, test_df)
        self.assertDictEqual(self.META_TYPES, types)


class TestParseSkylineAnnotations(TestReadMetadataBase, unittest.TestCase):

    ALL_META_TYPES = {'BioReplicate': Dtype.NULL,
                      'bool_var': Dtype.BOOL,
                      'Concentration': Dtype.NULL,
                      'float_var': Dtype.FLOAT,
                      'int_var': Dtype.INT,
                      'na_var': Dtype.NULL,
                      'string_var': Dtype.STRING}


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_skyline_csv(self):
        data, types = db_utils.read_metadata(f'{self.metadata_dir}/HeLa_annotations.csv')

        test_df = df_to_dict(data, self.META_TYPES)

        gt_data = self.remove_null_data(self.gt_data)
        gt_types = self.remove_null_types(self.META_TYPES)

        self.assertDataDictEqual(gt_data, test_df)
        self.assertDictEqual(gt_types, types)


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_full_skyline_csv(self):
        data, types = db_utils.read_metadata(f'{self.metadata_dir}/HeLa_all_annotations.csv')

        test_df = df_to_dict(data, self.ALL_META_TYPES)

        gt_data = self.remove_null_data(self.gt_data)
        gt_types = self.remove_null_types(self.ALL_META_TYPES)

        self.assertDataDictEqual(gt_data, test_df)
        self.assertDictEqual(gt_types, types)


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_full_skyline_csv_include_null(self):
        data, types = db_utils.read_metadata(f'{self.metadata_dir}/HeLa_all_annotations.csv',
                                             exclude_null_from_skyline=False)

        test_df = df_to_dict(data, self.ALL_META_TYPES)

        self.assertDataDictEqual(self.remove_null_data(self.gt_data),
                                 self.remove_null_data(test_df))
        self.assertDictEqual(self.ALL_META_TYPES, types)


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_skyline_csv_includ_null(self):
        data, types = db_utils.read_metadata(f'{self.metadata_dir}/HeLa_annotations.csv',
                                             exclude_null_from_skyline=False)

        test_df = df_to_dict(data, self.META_TYPES)

        self.assertDataDictEqual(self.remove_null_data(self.gt_data),
                                 self.remove_null_data(test_df))
        self.assertDictEqual(self.META_TYPES, types)



if __name__ == '__main__':
    unittest.main()
