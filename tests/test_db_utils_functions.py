
import unittest
import random
from unittest import mock
import os
import sqlite3
import pandas as pd

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils
from DIA_QC_report.submodules.metadata.dtype import Dtype


class TestDBHelperFunctions(unittest.TestCase):
    TEST_PROJECT = 'Strap'
    N_REPS = 20

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_metadata_fxns/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

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
                         'na_var': 'NULL'}
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
                        'na_var': 'NULL'}
        db_new_types = get_dtypes()
        self.assertDictEqual(result_types, db_new_types)


    def test_make_gene_matrix_by_protein_fails(self):
        self.assertTrue(self.conn is not None)

        self.assertEqual('protein', db_utils.get_meta_value(self.conn, 'group_precursors_by'))

        gene_id_path = f'{self.data_dir}/metadata/prhuman2gene_2023_05_24_subset.csv'
        command = ['make_gene_matrix', gene_id_path, self.db_path]
        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 1)
        self.assertTrue('Precursors in database must be grouped by gene!' in result.stderr)


class TestUpdateAcquiredRanks(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        dfs = []
        for project in ('Sp3', 'Strap'):
            fname = f'{setup_functions.TEST_DIR}/data/skyline_reports/{project}_replicate_quality.tsv'
            dfs.append(pd.read_csv(fname, sep='\t'))
            dfs[-1]['project'] = project

        cls.replicates = pd.concat(dfs, ignore_index=True)
        cls.replicates = cls.replicates.rename(columns={'Replicate': 'replicate',
                                                        'AcquiredTime': 'acquiredTime',
                                                        'TicArea': 'ticArea'})
        cls.replicates['acquiredRank'] = -1
        cls.replicates = cls.replicates[['replicate', 'project',
                                         'acquiredTime', 'acquiredRank',
                                         'ticArea']]

        cls.acquired_ranks = pd.read_csv(f'{setup_functions.TEST_DIR}/data/metadata/acquired_ranks.tsv',
                                         sep='\t')
        cls.acquired_ranks = {row.replicate: row.acquiredRank for row in cls.acquired_ranks.itertuples()}


    def setUp(self):
        self.conn = sqlite3.connect('file::memory:?cache=shared', uri=True)

        cur = self.conn.cursor()
        cur.execute('''
            CREATE TABLE replicates (
                replicateId INTEGER PRIMARY KEY,
                replicate TEXT NOT NULL,
                project TEXT NOT NULL,
                includeRep BOOL NOT NULL DEFAULT TRUE,
                acquiredTime BLOB NOT NULL,
                acquiredRank INTEGER NOT NULL,
                ticArea REAL NOT NULL,
                UNIQUE(replicate, project) ON CONFLICT FAIL
            );''')
        cur.execute('''
            CREATE TABLE metadata (
                key TEXT NOT NULL,
                value TEXT,
                PRIMARY KEY (key)
            ); ''')
        self.conn.commit()


    def tearDown(self):
        self.conn.close()


    @staticmethod
    def get_db_ranks(conn):
        cur = conn.cursor()
        cur.execute(''' SELECT replicate, acquiredRank FROM replicates
                        WHERE includeRep == TRUE; ''')
        return {row[0]: row[1] for row in cur.fetchall()}


    @staticmethod
    def update_ground_truth_ranks(reps):
        proj_ranks = {rep: TestUpdateAcquiredRanks.acquired_ranks[rep] for rep in reps}
        return {k: i for i, (k, v) in enumerate(sorted(proj_ranks.items(), key=lambda x: x[1]))}


    def test_update_acquired_ranks(self):
        self.replicates.to_sql('replicates', self.conn, if_exists='append',
                                index=True, index_label='replicateId')
        db_utils.update_acquired_ranks(self.conn)

        db_ranks = self.get_db_ranks(self.conn)
        self.assertDictEqual(self.acquired_ranks, db_ranks)


    def test_add_project(self):

        # add first project
        self.replicates[self.replicates['project'] == 'Sp3'].to_sql('replicates', self.conn,
                                                                    if_exists='append', index=True,
                                                                    index_label='replicateId')
        db_utils.update_acquired_ranks(self.conn)

        db_ranks = self.get_db_ranks(self.conn)
        proj_ranks = self.update_ground_truth_ranks(db_ranks.keys())
        self.assertDictEqual(proj_ranks, db_ranks)

        # add second project
        self.replicates[self.replicates['project'] == 'Strap'].to_sql('replicates', self.conn,
                                                                      if_exists='append', index=True,
                                                                      index_label='replicateId')

        db_utils.update_acquired_ranks(self.conn)
        db_ranks = self.get_db_ranks(self.conn)
        self.assertDictEqual(self.acquired_ranks, db_ranks)


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_mark_project_skipped(self):
        self.replicates.to_sql('replicates', self.conn, if_exists='append',
                                index=True, index_label='replicateId')
        db_utils.update_acquired_ranks(self.conn)

        db_ranks = self.get_db_ranks(self.conn)
        self.assertDictEqual(self.acquired_ranks, db_ranks)

        db_utils.mark_reps_skipped(self.conn, projects=('Strap',))
        db_utils.update_acquired_ranks(self.conn)

        db_ranks = self.get_db_ranks(self.conn)
        proj_ranks = self.update_ground_truth_ranks(db_ranks.keys())
        self.assertDictEqual(proj_ranks, db_ranks)


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_mark_replicates_skipped(self):
        self.replicates.to_sql('replicates', self.conn, if_exists='append',
                                index=True, index_label='replicateId')
        db_utils.update_acquired_ranks(self.conn)
        random.seed(1)
        max_rep_index = len(self.replicates.index) - 1
        max_subset_len = 15

        for _ in range(50):
            skip_reps = list()
            for i in random.sample(range(max_rep_index), random.randint(1, max_subset_len)):
                skip_reps.append(self.replicates.loc[i, 'replicate'])

            self.assertTrue(db_utils.mark_reps_skipped(self.conn, reps=skip_reps))
            db_utils.update_acquired_ranks(self.conn)

            db_ranks = self.get_db_ranks(self.conn)
            proj_ranks = self.update_ground_truth_ranks(db_ranks.keys())
            self.assertDictEqual(proj_ranks, db_ranks)
            db_utils.mark_all_reps_includced(self.conn)


if __name__ == '__main__':
    unittest.main()
