
import unittest
from unittest import mock
import json
from jsonschema import validate, ValidationError

import setup_functions

from DIA_QC_report.submodules.read_metadata import JSON_SCHEMA, Metadata
from DIA_QC_report.submodules.dtype import Dtype


def df_to_dict(df, types):
    ret = dict()
    for row in df.itertuples():
        if row.Replicate not in ret:
            ret[row.Replicate] = {}

        new_value = types[row.annotationKey].convert(row.annotationValue)
        ret[row.Replicate][row.annotationKey] = new_value

    return ret


class TestReadMetadataBase(unittest.TestCase, setup_functions.AbstractTestsBase):

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
            validate(cls.gt_data, JSON_SCHEMA)


    def setUp(self):
        self.meta_reader = Metadata()


    @staticmethod
    def remove_null_types(d):
        return {k: v for k, v in d.items() if v is not Dtype.NULL}


class TestParseMetadata(TestReadMetadataBase):
    def test_json(self):
        self.assertTrue(self.meta_reader.read(f'{self.metadata_dir}/HeLa_metadata.json'))
        self.assertEqual('json', self.meta_reader.input_format)

        test_df = df_to_dict(self.meta_reader.df, self.META_TYPES)

        self.assertDictEqual(self.META_TYPES, self.meta_reader.types)
        self.assertDataDictEqual(self.gt_data, test_df)


    def test_csv(self):
        self.assertTrue(self.meta_reader.read(f'{self.metadata_dir}/HeLa_metadata.csv'))
        self.assertEqual('csv', self.meta_reader.input_format)

        test_df = df_to_dict(self.meta_reader.df, self.META_TYPES)

        self.assertDictEqual(self.META_TYPES, self.meta_reader.types)
        self.assertDataDictEqual(self.gt_data, test_df)


    def test_tsv(self):
        self.assertTrue(self.meta_reader.read(f'{self.metadata_dir}/HeLa_metadata.tsv'))
        self.assertEqual('tsv', self.meta_reader.input_format)

        test_df = df_to_dict(self.meta_reader.df, self.META_TYPES)

        self.assertDictEqual(self.META_TYPES, self.meta_reader.types)
        self.assertDataDictEqual(self.gt_data, test_df)


class TestParseMissingMetadata(TestReadMetadataBase):
    def test_tsv(self):
        self.assertTrue(self.meta_reader.read(f'{self.metadata_dir}/Strap_missing_multi_var_metadata.tsv'))
        self.assertDictEqual(self.META_TYPES, self.meta_reader.types)


    def test_json(self):
        self.assertTrue(self.meta_reader.read(f'{self.metadata_dir}/Strap_missing_multi_var_metadata.json'))
        self.assertDictEqual(self.META_TYPES, self.meta_reader.types)


class TestParseSkylineAnnotations(TestReadMetadataBase):
    ALL_META_TYPES = {'BioReplicate': Dtype.NULL,
                      'bool_var': Dtype.BOOL,
                      'Concentration': Dtype.NULL,
                      'float_var': Dtype.FLOAT,
                      'int_var': Dtype.INT,
                      'na_var': Dtype.NULL,
                      'string_var': Dtype.STRING}


    @staticmethod
    def remove_null_data(d):
        na_types = {k for k, v in TestParseSkylineAnnotations.ALL_META_TYPES.items() if v is Dtype.NULL}
        return {rep: {k: v for k, v in row.items() if k not in na_types} for rep, row in d.items()}


    @mock.patch('DIA_QC_report.submodules.read_metadata.LOGGER', mock.Mock())
    def test_skyline_csv(self):
        self.assertTrue(self.meta_reader.read(f'{self.metadata_dir}/HeLa_annotations.csv'))
        self.assertEqual('skyline', self.meta_reader.input_format)

        test_df = df_to_dict(self.meta_reader.df, self.META_TYPES)

        gt_data = self.remove_null_data(self.gt_data)
        gt_types = self.remove_null_types(self.META_TYPES)

        self.assertDataDictEqual(gt_data, test_df)
        self.assertDictEqual(gt_types, self.meta_reader.types)


    @mock.patch('DIA_QC_report.submodules.read_metadata.LOGGER', mock.Mock())
    def test_full_skyline_csv(self):
        self.assertTrue(self.meta_reader.read(f'{self.metadata_dir}/HeLa_all_annotations.csv'))
        self.assertEqual('skyline', self.meta_reader.input_format)

        test_df = df_to_dict(self.meta_reader.df, self.ALL_META_TYPES)

        gt_data = self.remove_null_data(self.gt_data)
        gt_types = self.remove_null_types(self.ALL_META_TYPES)

        self.assertDataDictEqual(gt_data, test_df)
        self.assertDictEqual(gt_types, self.meta_reader.types)


    @mock.patch('DIA_QC_report.submodules.read_metadata.LOGGER', mock.Mock())
    def test_full_skyline_csv_include_null(self):
        self.assertTrue(self.meta_reader.read(f'{self.metadata_dir}/HeLa_all_annotations.csv',
                                              exclude_null_from_skyline=False))

        test_df = df_to_dict(self.meta_reader.df, self.ALL_META_TYPES)

        self.assertDataDictEqual(self.remove_null_data(self.gt_data),
                                 self.remove_null_data(test_df))
        self.assertDictEqual(self.ALL_META_TYPES, self.meta_reader.types)


    @mock.patch('DIA_QC_report.submodules.read_metadata.LOGGER', mock.Mock())
    def test_skyline_csv_includ_null(self):
        self.assertTrue(self.meta_reader.read(f'{self.metadata_dir}/HeLa_annotations.csv',
                                              exclude_null_from_skyline=False))

        test_df = df_to_dict(self.meta_reader.df, self.META_TYPES)

        self.assertDataDictEqual(self.remove_null_data(self.gt_data),
                                 self.remove_null_data(test_df))
        self.assertDictEqual(self.META_TYPES, self.meta_reader.types)


if __name__ == '__main__':
    unittest.main()
