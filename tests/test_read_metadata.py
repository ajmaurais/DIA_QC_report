
import unittest
from io import StringIO
from itertools import product
import json
import re
import csv
import random

from jsonschema import validate
import pandas as pd

import setup_functions

from DIA_QC_report.submodules import logger
from DIA_QC_report.submodules.read_metadata import JSON_SCHEMA, Metadata
from DIA_QC_report.submodules.dtype import Dtype


class TestMetadataBase(unittest.TestCase, setup_functions.AbstractTestsBase):
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


class TestMetadataClassMethods(TestMetadataBase):
    def test_parse_from_fp_error(self):
        meta_reader = Metadata()
        with open(f'{self.metadata_dir}/HeLa_metadata.json') as inF:
            with self.assertRaises(RuntimeError):
                meta_reader.read(inF)


    def test_validate_duplicate_keys(self):
        meta_reader = Metadata()
        self.assertTrue(meta_reader.read(f'{self.metadata_dir}/HeLa_metadata.json'))
        self.assertEqual('json', meta_reader.input_format)

        meta_reader.df = pd.concat([meta_reader.df,
                                    meta_reader.df.loc[0:2]]).reset_index(drop=True)

        with self.assertLogs(logger.LOGGER, level='ERROR') as cm:
            self.assertFalse(meta_reader.validate())
        self.assertTrue("Duplicate key: '", cm.output[-1])


    def test_validate_missing_type_key(self):
        meta_reader = Metadata()
        self.assertTrue(meta_reader.read(f'{self.metadata_dir}/HeLa_metadata.json'))
        self.assertEqual('json', meta_reader.input_format)

        key, _ = meta_reader.types.popitem()

        with self.assertLogs(logger.LOGGER, level='ERROR') as cm:
            self.assertFalse(meta_reader.validate())
        self.assertTrue(f"Missing key: '{key}' in Metadata.types!", cm.output[-1])


class TestMetadataEqMethod(TestMetadataBase):
    ''' Test custom Metadata __eq__ method '''

    FILES = ['tsv', 'csv', 'json', 'skyline']

    def read_ext(self, ext):
        ''' Read metadata file with specified `ext` '''
        suffix = 'metadata'
        if ext == 'skyline':
            ext = 'csv'
            suffix = 'annotations'
        ret = Metadata(quiet=True)
        self.assertTrue(ret.read(f'{self.metadata_dir}/HeLa_{suffix}.{ext}',
                                 exclude_null_from_skyline=False))

        return ret


    def test_eq(self):
        ''' Test that each permutation of file extension are equal '''
        for l_ext, r_ext in product(self.FILES, self.FILES):
            lhs = self.read_ext(l_ext)
            rhs = self.read_ext(r_ext)

            self.assertEqual(lhs, rhs)


    def test_neq_different_type(self):
        ''' Test different types are not equal '''
        for l_ext, r_ext in product(self.FILES, self.FILES):
            rhs = self.read_ext(r_ext)

            lhs = self.read_ext(l_ext)
            lhs.types['int_var'] = Dtype.BOOL

            self.assertNotEqual(lhs, rhs)

        for l_ext, r_ext in product(self.FILES, self.FILES):
            lhs = self.read_ext(r_ext)

            rhs = self.read_ext(r_ext)
            rhs.types['int_var'] = Dtype.BOOL

            self.assertNotEqual(lhs, rhs)


    def test_neq_different_key(self):
        ''' Test different annotationKeys are not equal '''
        for l_ext, r_ext in product(self.FILES, self.FILES):
            rhs = self.read_ext(r_ext)

            lhs = self.read_ext(l_ext)
            key, t = lhs.types.popitem()
            lhs.types['other_key'] = t
            lhs.df.loc[lhs.df['annotationKey'] == key, 'annotationKey'] = 'other_key'

            self.assertNotEqual(lhs, rhs)

        for l_ext, r_ext in product(self.FILES, self.FILES):
            lhs = self.read_ext(l_ext)

            rhs = self.read_ext(r_ext)
            key, t = rhs.types.popitem()
            rhs.types['other_key'] = t
            rhs.df.loc[rhs.df['annotationKey'] == key, 'annotationKey'] = 'other_key'

            self.assertNotEqual(lhs, rhs)


    def test_neq_different_replicate(self):
        ''' Test different replicates are not equal '''
        for l_ext, r_ext in product(self.FILES, self.FILES):
            rhs = self.read_ext(r_ext)

            lhs = self.read_ext(l_ext)
            rep = lhs.df['Replicate'].drop_duplicates().to_list()[0]
            lhs.df.loc[lhs.df['Replicate'] == rep, 'Replicate'] = 'other_rep'

            self.assertNotEqual(lhs, rhs)

        for l_ext, r_ext in product(self.FILES, self.FILES):
            lhs = self.read_ext(l_ext)

            rhs = self.read_ext(l_ext)
            rep = rhs.df['Replicate'].drop_duplicates().to_list()[0]
            rhs.df.loc[rhs.df['Replicate'] == rep, 'Replicate'] = 'other_rep'

            self.assertNotEqual(lhs, rhs)


    def test_neq_different_value(self):
        ''' Test different annotationValues are not equal '''
        for l_ext, r_ext in product(self.FILES, self.FILES):
            rhs = self.read_ext(r_ext)

            lhs = self.read_ext(l_ext)
            lhs.df.loc[lhs.df['annotationKey'] == 'string_var', 'annotationValue'] = 'other_value'

            self.assertNotEqual(lhs, rhs)

        for l_ext, r_ext in product(self.FILES, self.FILES):
            lhs = self.read_ext(l_ext)

            rhs = self.read_ext(r_ext)
            rhs.df.loc[rhs.df['annotationKey'] == 'string_var', 'annotationValue'] = 'other_value'

            self.assertNotEqual(lhs, rhs)


class TestWithExtension(TestMetadataBase):
    '''
    Test that results are equal regardless of whether replicate names have a file extension.
    '''

    def run_table_tests(self, ext, delim, replicate_ext):
        ext = 'tsv'
        delim = '\t'
        no_ext = Metadata()
        with_ext = Metadata()
        with open(f'{self.metadata_dir}/HeLa_metadata.{ext}') as inF:
            self.assertTrue(no_ext.read(inF, metadata_format=ext))

            # read fp into csv reader and add file extensions to metadata
            inF.seek(0)
            reader = csv.reader(inF, delimiter=delim)
            rows = list()
            rows.append(next(reader))
            for row in reader:
                rows.append(row)
                rows[-1][0] += f'.{replicate_ext}'

        string_stream = StringIO('\n'.join(delim.join(row) for row in rows))
        self.assertTrue(with_ext.read(string_stream, metadata_format=ext))

        self.assertEqual(no_ext, with_ext)
        self.assertDataDictEqual(with_ext.df_to_dict(), self.gt_data)


    def test_tsv(self):
        self.run_table_tests('tsv', '\t', 'mzML')
        self.run_table_tests('tsv', '\t', 'raw')


    def test_csv(self):
        self.run_table_tests('csv', ',', 'mzML')
        self.run_table_tests('csv', ',', 'raw')


    def test_json(self):
        ext = 'json'
        no_ext = Metadata()
        with_ext = Metadata()
        with open(f'{self.metadata_dir}/HeLa_metadata.{ext}') as inF:
            self.assertTrue(no_ext.read(inF, metadata_format=ext))

            # read fp into dict
            inF.seek(0)
            data = json.load(inF)

        for rep_ext in ('raw', 'mzML'):
            test_data = {key + f'.{rep_ext}': value for key, value in data.items()}

            string_stream = StringIO(json.dumps(test_data))
            self.assertTrue(with_ext.read(string_stream, metadata_format=ext))

            self.assertEqual(no_ext, with_ext)
            self.assertDataDictEqual(with_ext.df_to_dict(), self.gt_data)


class TestReadMetadata(TestMetadataBase):
    def do_test_extension(self, ext, skyline=False):
        path_reader = Metadata()
        fp_reader = Metadata()

        file_suffix = '_annotations' if skyline else '_metadata'

        # read fp
        with open(f'{self.metadata_dir}/HeLa{file_suffix}.{ext}') as inF:
            fp_reader.read(inF, metadata_format=ext,
                           exclude_null_from_skyline=not skyline)
        self.assertDictEqual(self.META_TYPES, fp_reader.types)
        self.assertDataDictEqual(self.gt_data, fp_reader.df_to_dict())

        # read file
        self.assertTrue(path_reader.read(f'{self.metadata_dir}/HeLa{file_suffix}.{ext}',
                                         exclude_null_from_skyline=not skyline))
        self.assertEqual('skyline' if skyline else ext,
                         path_reader.input_format)
        self.assertDictEqual(self.META_TYPES, path_reader.types)
        self.assertDataDictEqual(self.gt_data, path_reader.df_to_dict())

        # test equality
        self.assertEqual(path_reader, fp_reader)


    def test_json(self):
        self.do_test_extension('json')


    def test_csv(self):
        self.do_test_extension('csv')


    def test_tsv(self):
        self.do_test_extension('tsv')


    def test_skyline_csv(self):
        with self.assertLogs(logger.LOGGER) as cm:
            self.do_test_extension('csv', skyline=True)
        self.assertTrue('Found Skyline annotations csv.' in cm.output[0])


class TestReadMissingMetadata(TestMetadataBase):
    def test_tsv(self):
        meta_reader = Metadata()
        self.assertTrue(meta_reader.read(f'{self.metadata_dir}/Strap_missing_multi_var_metadata.tsv'))
        self.assertDictEqual(self.META_TYPES, meta_reader.types)


    def test_json(self):
        meta_reader = Metadata()
        self.assertTrue(meta_reader.read(f'{self.metadata_dir}/Strap_missing_multi_var_metadata.json'))
        self.assertDictEqual(self.META_TYPES, meta_reader.types)


class TestParseSkylineAnnotations(TestMetadataBase):
    ALL_META_TYPES = {'BioReplicate': Dtype.NULL,
                      'bool_var': Dtype.BOOL,
                      'Concentration': Dtype.NULL,
                      'float_var': Dtype.FLOAT,
                      'int_var': Dtype.INT,
                      'na_var': Dtype.NULL,
                      'string_var': Dtype.STRING}


    @staticmethod
    def remove_null_types(d):
        return {k: v for k, v in d.items() if v is not Dtype.NULL}


    @staticmethod
    def remove_null_data(d):
        na_types = {k for k, v in TestParseSkylineAnnotations.ALL_META_TYPES.items() if v is Dtype.NULL}
        return {rep: {k: v for k, v in row.items() if k not in na_types} for rep, row in d.items()}


    def test_skyline_csv(self):
        meta_reader = Metadata(quiet=True)
        self.assertTrue(meta_reader.read(f'{self.metadata_dir}/HeLa_annotations.csv'))
        self.assertEqual('skyline', meta_reader.input_format)

        gt_data = self.remove_null_data(self.gt_data)
        gt_types = self.remove_null_types(self.META_TYPES)

        self.assertDataDictEqual(gt_data, meta_reader.df_to_dict())
        self.assertDictEqual(gt_types, meta_reader.types)


    def test_full_skyline_csv(self):
        meta_reader = Metadata(quiet=True)
        self.assertTrue(meta_reader.read(f'{self.metadata_dir}/HeLa_all_annotations.csv'))
        self.assertEqual('skyline', meta_reader.input_format)

        gt_data = self.remove_null_data(self.gt_data)
        gt_types = self.remove_null_types(self.ALL_META_TYPES)

        self.assertDataDictEqual(gt_data, meta_reader.df_to_dict())
        self.assertDictEqual(gt_types, meta_reader.types)


    def test_full_skyline_csv_include_null(self):
        meta_reader = Metadata(quiet=True)
        self.assertTrue(meta_reader.read(f'{self.metadata_dir}/HeLa_all_annotations.csv',
                                         exclude_null_from_skyline=False))

        self.assertDataDictEqual(self.remove_null_data(self.gt_data),
                                 self.remove_null_data(meta_reader.df_to_dict()))
        self.assertDictEqual(self.ALL_META_TYPES, meta_reader.types)


    def test_skyline_csv_include_null(self):
        meta_reader = Metadata(quiet=True)
        self.assertTrue(meta_reader.read(f'{self.metadata_dir}/HeLa_annotations.csv',
                                         exclude_null_from_skyline=False))

        self.assertDataDictEqual(self.remove_null_data(self.gt_data),
                                 self.remove_null_data(meta_reader.df_to_dict()))
        self.assertDictEqual(self.META_TYPES, meta_reader.types)


class TestMetadataWriteMethods(unittest.TestCase):
    FILES = ['tsv', 'csv', 'json', 'skyline']

    @classmethod
    def setUpClass(cls):
        cls.metadata_dir = f'{setup_functions.TEST_DIR}/data/metadata'


    def do_test(self, write_ext,
                prefix='HeLa', files=FILES):

        write_method = 'to_skyline_annotations' if write_ext == 'skyline' else f'to_{write_ext}'

        for ext in files:
            suffix = f'metadata.{ext}' if ext != 'skyline' else 'annotations.csv'

            writer = Metadata(quiet=True)
            self.assertTrue(writer.read(f'{self.metadata_dir}/{prefix}_{suffix}'))

            exclude_null_from_skyline = True
            if write_ext == 'skyline':
                exclude_null_from_skyline = ext == 'skyline'

            reader = Metadata(quiet=True)
            with StringIO() as sstream:
                getattr(writer, write_method)(sstream)
                sstream.seek(0)
                self.assertTrue(reader.read(sstream,
                    metadata_format=write_ext if write_ext != 'skyline' else 'csv',
                    exclude_null_from_skyline=exclude_null_from_skyline))

            self.assertEqual(writer, reader)


    def test_to_tsv(self):
        self.do_test('tsv')


    def test_to_csv(self):
        self.do_test('csv')


    def test_to_json(self):
        self.do_test('json')


    def test_to_skyline_annotations(self):
        self.do_test('skyline')


    def test_to_tsv_missing(self):
        self.do_test('tsv', prefix='Strap_missing_multi_var',
                     files=['tsv', 'json'])


    def test_to_csv_missing(self):
        self.do_test('csv', prefix='Strap_missing_multi_var',
                     files=['tsv', 'json'])


    def test_to_json_missing(self):
        self.do_test('json', prefix='Strap_missing_multi_var',
                     files=['tsv', 'json'])


    # def test_to_skyline_annotations_missing(self):
    #     self.do_test('skyline', prefix='Strap_missing_multi_var',
    #                  files=['tsv', 'json'])


class TestToSkylineAnnotations(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metadata_dir = f'{setup_functions.TEST_DIR}/data/metadata'


    def test_missing(self):
        reader = Metadata()
        self.assertTrue(reader.read(f'{self.metadata_dir}/Strap_missing_multi_var_metadata.tsv'))

        with StringIO() as sstream:
            reader.to_skyline_annotations(sstream)
            sstream.seek(0)
            rows = [row for row in csv.DictReader(sstream)]

        for row in rows:
            self.assertEqual('', row['annotation_na_var'])

            if not re.search(r'[0-9]+', row['annotation_int_var']):
                self.assertEqual('', row['annotation_int_var'])


class TestToSkylineDefinitions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.metadata_dir = f'{setup_functions.TEST_DIR}/data/metadata'


    def do_test(self, ext, prefix='HeLa'):
        command_re = re.compile(r'^--annotation-name="([\w\.\-]+)" --annotation-targets=replicate --annotation-type=([a-z_]+)$')

        suffix = f'metadata.{ext}' if ext != 'skyline' else 'annotations.csv'

        writer = Metadata(quiet=True)
        self.assertTrue(writer.read(f'{self.metadata_dir}/{prefix}_{suffix}'))

        sstream = StringIO()
        writer.to_skyline_definitions(sstream)
        sstream.seek(0)

        for line in sstream:
            self.assertTrue((match := command_re.search(line)) is not None)
            self.assertTrue(match[1] in writer.types)
            self.assertEqual(Dtype.to_sky_type(writer.types[match[1]]), match[2])


    def test_tsv(self):
        self.do_test('tsv')


    def test_csv(self):
        self.do_test('csv')


    def test_json(self):
        self.do_test('json')


    def test_json_skyline(self):
        self.do_test('skyline')


    def test_csv_missing(self):
        self.do_test('tsv', prefix='Strap_missing_multi_var')


    def test_json_missing(self):
        self.do_test('json', prefix='Strap_missing_multi_var')


class TestJsonMissingKeys(unittest.TestCase, setup_functions.AbstractTestsBase):
    def setUp(self):
        self.reader = Metadata()
        with open(f'{setup_functions.TEST_DIR}/data/metadata/HeLa_metadata.json') as inF:
            if not self.reader.read(inF, metadata_format='json'):
                raise RuntimeError('Failed to read HeLa metadata!')

        random.seed(12)
        self.bad_replicate = random.choice(self.reader.df['Replicate'].drop_duplicates().to_list())
        self.missing_key = random.choice(list(self.reader.types.keys()))
        self.reader.df = self.reader.df[~((self.reader.df['Replicate'] == self.bad_replicate) &
                                           (self.reader.df['annotationKey'] == self.missing_key))].reset_index(drop=True)

        row_data = {rep: {} for rep in self.reader.df['Replicate'].drop_duplicates().to_list()}
        for row in self.reader.df.itertuples():
            row_data[row.Replicate][row.annotationKey] = self.reader.types[row.annotationKey].convert(row.annotationValue)

        self.data = json.dumps(row_data)


    def test_read(self):
        reader = Metadata()
        with self.assertLogs(logger.LOGGER, level='WARNING') as cm:
            self.assertTrue(reader.read(StringIO(self.data), metadata_format='json'))

        self.assertInLog(
            f"annotationKey '{self.missing_key}' is missing for replicate '{self.bad_replicate}'!", cm
        )


    def test_write(self):
        sstream = StringIO()
        self.reader.to_json(sstream)
        sstream.seek(0)

        out_data = json.load(sstream)

        for row in out_data.values():
            for key in self.reader.types:
                self.assertIn(key, row)

        self.assertIsNone(out_data[self.bad_replicate][self.missing_key])


class TestAdd(unittest.TestCase):
    def do_add_test(self, datas, target_types):
        readers = []
        for data in datas:
            readers.append(Metadata())
            readers[-1].read(StringIO(json.dumps(data)), metadata_format='json')

        combined = Metadata()
        for reader in readers:
            combined.add(reader)

        with self.subTest('df_structure'):
            self.assertEqual(len(combined.df['Replicate'].drop_duplicates()), sum(len(data) for data in datas))
            self.assertIsInstance(combined, Metadata)
            self.assertDictEqual(combined.types, target_types)

        with self.subTest('df_to_dict'):
            d = combined.df_to_dict()
            self.assertIsInstance(d, dict)
            self.assertEqual(len(d), len(combined.df['Replicate'].drop_duplicates()))
            for row in d.values():
                self.assertIsInstance(row, dict)
                self.assertEqual(set(row.keys()), set(target_types.keys()))

        with self.subTest('get_wide_data'):
            wide_data = combined.get_wide_data()
            self.assertIsInstance(wide_data, pd.DataFrame)
            self.assertEqual(len(wide_data), len(combined.df['Replicate'].drop_duplicates()))
            self.assertTrue(set(target_types.keys()).issubset(set(wide_data.columns)))


    def test_add(self):
        data_1 = {
           "Rep1": { "string_var": "value1", "int_var": 1, "float_var": 1.0, "bool_var": True },
           "Rep2": { "string_var": "value2", "int_var": 2, "float_var": 2.0, "bool_var": False }
        }
        data_2 = {
           "Rep3": { "string_var": "value3", "int_var": 3, "float_var": 3.0, "bool_var": True },
           "Rep4": { "string_var": "value4", "int_var": 4, "float_var": 4.0, "bool_var": False }
        }
        target_types = {
            "string_var": Dtype.STRING,
            "int_var": Dtype.INT,
            "float_var": Dtype.FLOAT,
            "bool_var": Dtype.BOOL
        }
        self.do_add_test([data_1, data_2], target_types)


    def test_add_disperate_types(self):
        data_1 = {'Rep1': {'string_var': 'value1', 'int_var': 1, 'float_var': 1.0, 'bool_var': True}}
        data_2 = {'Rep2': {'string_var': 'value2', 'int_var': 1.0, 'float_var': '1', 'bool_var': 'True'}}
        target_types = {
            'string_var': Dtype.STRING,
            'int_var': Dtype.FLOAT,     # coerced to FLOAT
            'float_var': Dtype.STRING,  # coerced to STRING
            'bool_var': Dtype.STRING    # coerced to STRING
        }
        self.do_add_test([data_1, data_2], target_types)


    def test_add_missing_string(self):
        data_1 = {'Rep1': {'string_var': 'value1', 'int_var': 1, 'float_var': 1.0, 'bool_var': True}}
        data_2 = {'Rep2': {'int_var': 2, 'float_var': 2.0, 'bool_var': False}}
        target_types = {
            'string_var': Dtype.STRING,  # remains STRING
            'int_var': Dtype.INT,
            'float_var': Dtype.FLOAT,
            'bool_var': Dtype.BOOL
        }
        self.do_add_test([data_1, data_2], target_types)


    def test_add_missing_int(self):
        data_1 = {'Rep1': {'string_var': 'value1', 'int_var': 1, 'float_var': 1.0, 'bool_var': True}}
        data_2 = {'Rep2': {'string_var': 'value2', 'float_var': 2.0, 'bool_var': False}}
        target_types = {
            'string_var': Dtype.STRING,
            'int_var': Dtype.INT,
            'float_var': Dtype.FLOAT,
            'bool_var': Dtype.BOOL
        }
        self.do_add_test([data_1, data_2], target_types)


    def test_add_missing_float(self):
        data_1 = {'Rep1': {'string_var': 'value1', 'int_var': 1, 'float_var': 1.0, 'bool_var': True}}
        data_2 = {'Rep2': {'string_var': 'value2', 'int_var': 2, 'bool_var': False}}
        target_types = {
            'string_var': Dtype.STRING,
            'int_var': Dtype.INT,
            'float_var': Dtype.FLOAT,
            'bool_var': Dtype.BOOL
        }
        self.do_add_test([data_1, data_2], target_types)


    def test_add_missing_bool(self):
        data_1 = {'Rep1': {'string_var': 'value1', 'int_var': 1, 'float_var': 1.0, 'bool_var': True}}
        data_2 = {'Rep2': {'string_var': 'value2', 'int_var': 2, 'float_var': 2.0}}
        target_types = {
            'string_var': Dtype.STRING,
            'int_var': Dtype.INT,
            'float_var': Dtype.FLOAT,
            'bool_var': Dtype.BOOL
        }
        self.do_add_test([data_1, data_2], target_types)


    def test_add_missing_na(self):
        data_1 = {'Rep1': {'string_var': 'value1', 'int_var': 1, 'float_var': 1.0, 'bool_var': True}}
        data_2 = {'Rep2': {'string_var': 'value2', 'int_var': 2, 'float_var': 2.0, 'bool_var': False, 'na_var': None}}
        target_types = {
            'string_var': Dtype.STRING,
            'int_var': Dtype.INT,
            'float_var': Dtype.FLOAT,
            'bool_var': Dtype.BOOL,
            'na_var': Dtype.NULL
        }
        self.do_add_test([data_1, data_2], target_types)


    def test_add_empty(self):
        data_1 = {'Rep1': {'string_var': 'value1'}}
        reader_1 = Metadata()
        reader_1.read(StringIO(json.dumps(data_1)), metadata_format='json')
        reader_2 = Metadata()

        combined = Metadata()
        combined.add(reader_1)
        combined.add(reader_2)
        self.assertEqual(len(combined.df['Replicate'].drop_duplicates()), 1)
        self.assertDictEqual(combined.types, reader_1.types)


    def test_both_empty(self):
        reader_1 = Metadata()
        reader_2 = Metadata()
        combined = Metadata()

        combined.add(reader_1)
        combined.add(reader_2)
        self.assertIsInstance(combined, Metadata)
        self.assertTrue(combined.df.empty)
        self.assertDictEqual(combined.types, {})


    def test_add_duplicate_reps(self):
        ''' Adding 2 metadata files with the same replicate should raise an error '''
        data_1 = {'Rep1': {'string_var': 'value1'}}
        data_2 = {'Rep1': {'string_var': 'value2'}}
        reader_1 = Metadata()
        reader_1.read(StringIO(json.dumps(data_1)), metadata_format='json')
        reader_2 = Metadata()
        reader_2.read(StringIO(json.dumps(data_2)), metadata_format='json')

        with self.assertRaises(ValueError):
            reader_1.add(reader_2)