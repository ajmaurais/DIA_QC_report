
import unittest
import random
from unittest import mock
import sys
import os
import sqlite3
from itertools import product

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils
from DIA_QC_report.submodules.metadata import Dtype

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


class TestDBHelperFunctions(unittest.TestCase):
    TEST_PROJECT = 'Strap'
    N_REPS = 20

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{TEST_DIR}/work/test_metadata_fxns/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{TEST_DIR}/data/'

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           metadata_suffix='_multi_var_metadata.tsv',
                                                           clear_dir=True)

        if cls.parse_result.returncode != 0:
            raise RuntimeError('Setup of test db failed!')

        cls.conn = None
        if os.path.isfile(cls.db_path):
            cls.conn = sqlite3.connect(cls.db_path)


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_mark_reps_skipped_errors(self):
        self.assertTrue(self.conn is not None)
        self.assertFalse(db_utils.mark_reps_skipped(self.conn, reps=('dummy',)))


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_mark_reps_skipped_included(self):
        self.assertTrue(self.conn is not None)

        def excluded_reps():
            cur = self.conn.cursor()
            cur.execute('SELECT replicate FROM replicates WHERE includeRep == FALSE;')
            return [x[0] for x in cur.fetchall()]

        def included_reps():
            cur = self.conn.cursor()
            cur.execute('SELECT replicate FROM replicates WHERE includeRep == TRUE;')
            return [x[0] for x in cur.fetchall()]

        included = included_reps()
        self.assertTrue(len(included), self.N_REPS)

        # Select single rep to exclude
        random.seed(1)
        single_rep = included.pop(random.randint(0, self.N_REPS - 1))
        self.assertTrue(db_utils.mark_reps_skipped(self.conn, reps=[single_rep]))
        self.assertEqual(len(excluded_reps()), 1)
        self.assertEqual(len(included_reps()), self.N_REPS - 1)
        self.assertTrue(single_rep in excluded_reps())

        # Select multiple reps to exclude
        excluded = [included[i] for i in {random.randint(0, self.N_REPS - 2) for _ in range(10)}]
        included = [x for x in included if x not in excluded]
        self.assertTrue(db_utils.mark_reps_skipped(self.conn, reps=excluded))
        excluded.append(single_rep)
        db_excluded = excluded_reps()
        self.assertEqual(len(db_excluded), len(excluded))
        self.assertEqual(len(included_reps()), len(included))
        self.assertTrue(all(rep in db_excluded for rep in excluded + [single_rep]))

        # Set all reps included
        db_utils.mark_all_reps_includced(self.conn)
        included = included_reps()
        self.assertEqual(len(included), self.N_REPS)
        self.assertEqual(len(excluded_reps()), 0)

        # Exclude project
        self.assertTrue(db_utils.mark_reps_skipped(self.conn, projects=[self.TEST_PROJECT]))
        self.assertEqual(len(excluded_reps()), self.N_REPS)
        self.assertEqual(len(included_reps()), 0)

        # Set all reps included again
        db_utils.mark_all_reps_includced(self.conn)
        self.assertEqual(len(included_reps()), self.N_REPS)
        self.assertEqual(len(excluded_reps()), 0)


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_metadata_functions(self):
        self.assertTrue(self.conn is not None)

        # validate current version
        self.assertTrue(db_utils.check_schema_version(self.conn))

        # test get_meta_value
        schema_version = db_utils.get_meta_value(self.conn, 'schema_version')
        self.assertEqual(db_utils.SCHEMA_VERSION, schema_version)

        # set current version to 'NULL' and make sure validation fails
        self.conn = db_utils.update_meta_value(self.conn, 'schema_version', 'NULL')
        self.assertFalse(db_utils.check_schema_version(self.conn))
        self.assertEqual(db_utils.get_meta_value(self.conn, 'schema_version'), 'NULL')

        # set schema back to correct version
        self.conn = db_utils.update_meta_value(self.conn, 'schema_version', schema_version)
        self.assertEqual(db_utils.get_meta_value(self.conn, 'schema_version'), schema_version)
        self.assertTrue(db_utils.check_schema_version(self.conn))


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_is_normalized(self):
        self.assertTrue(self.conn is not None)

        self.assertFalse(db_utils.is_normalized(self.conn))
        normalize_command = ['normalize_db', self.db_path]
        normalize_result = setup_functions.run_command(normalize_command,
                                                       self.work_dir,
                                                       prefix='normalize_single_proj')

        self.assertEqual(normalize_result.returncode, 0)
        self.assertTrue(db_utils.is_normalized(self.conn))


    def test_update_metadata_types(self):
        self.assertTrue(self.conn is not None)

        def get_dtypes():
            cur = self.conn.cursor()
            cur.execute('SELECT annotationKey, annotationType FROM sampleMetadataTypes;')
            return dict((x[0], x[1]) for x in cur.fetchall())

        current_types = {'string_var': 'STRING',
                         'bool_var': 'BOOL',
                         'int_var': 'INT',
                         'float_var': 'FLOAT',
                         'na_var': 'BOOL'}
        db_current_types = get_dtypes()
        self.assertDictEqual(current_types, db_current_types)

        new_types = {'string_var': Dtype.INT,
                     'bool_var': Dtype.STRING,
                     'int_var': Dtype.FLOAT,
                     'float_var': Dtype.INT,
                     'na_var': Dtype.NULL}
        self.conn = db_utils.update_metadata_dtypes(self.conn, new_types)

        result_types = {'string_var': 'STRING',
                        'bool_var': 'STRING',
                        'int_var': 'FLOAT',
                        'float_var': 'FLOAT',
                        'na_var': 'BOOL'}
        db_new_types = get_dtypes()
        self.assertDictEqual(result_types, db_new_types)


    def test_make_gene_matrix_by_protein_fails(self):
        self.assertTrue(self.conn is not None)

        self.assertEqual('protein', db_utils.get_meta_value(self.conn, 'group_precursors_by'))

        gene_id_path = f'{self.data_dir}/metadata/prhuman2gene_2023_05_24_subset.csv'
        command = ['make_gene_matrix', gene_id_path, self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 1)
        self.assertTrue('Precursors in database must be grouped by gene!' in result.stderr.decode('utf-8'))


    # Not sure if I am going to include this one
    # def test_update_acquired_ranks(self):
    #     pass


if __name__ == '__main__':
    unittest.main()
