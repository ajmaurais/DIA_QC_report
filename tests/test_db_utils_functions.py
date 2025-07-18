
import unittest
import random
import re
from itertools import combinations
from unittest import mock
import os
import sqlite3
from socket import gethostname
from datetime import datetime

import pandas as pd

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils
from DIA_QC_report.submodules import logger
from DIA_QC_report.submodules.dtype import Dtype
from DIA_QC_report.submodules.read_metadata import Metadata
from DIA_QC_report.submodules.skyline_reports import ReplicateReport
from DIA_QC_report import __version__ as PROGRAM_VERSION
from DIA_QC_report import normalize_db
from DIA_QC_report import export_gene_matrix


class TestDBFunctionsBase(setup_functions.AbstractTestsBase):
    def mark_reps_skipped(self, *args, **kwargs):
        with self.assertLogs(logger.LOGGER, level='INFO') as cm:
            ret=db_utils.mark_reps_skipped(*args, **kwargs)

        self.assertEqual(len(cm.output), 1)
        self.assertIsNotNone(re.search('Excluding [0-9]+ replicates.', cm.output[0]))

        return ret


    def mark_all_reps_included(self, *args, **kwargs):
        with self.assertLogs(logger.LOGGER, level='INFO') as cm:
            db_utils.mark_all_reps_included(*args, **kwargs)

        self.assertEqual(len(cm.output), 1)
        messages = [re.compile(m) for m in ('All replicates are already included.',
                                            'Setting [0-9]+ includeRep values to TRUE')]
        self.assertTrue(any(p.search(cm.output[0]) is not None for p in messages))


class TestDBHelperFunctions(unittest.TestCase, TestDBFunctionsBase):
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


    def test_project_exists(self):
        self.assertIsNotNone(self.conn)
        self.assertTrue(db_utils.project_exists(self.conn, self.TEST_PROJECT))
        self.assertFalse(db_utils.project_exists(self.conn, 'NOT_A_PROJECT'))


    def test_mark_reps_skipped_errors(self):
        self.assertTrue(self.conn is not None)
        with self.assertLogs(logger.LOGGER, level='ERROR') as cm:
            self.assertFalse(db_utils.mark_reps_skipped(self.conn, reps=('dummy',)))

        self.assertEqual(len(cm.output), 1)
        self.assertIsNotNone("Replicate 'dummy' is not in database!" in cm.output[0])


    def test_mark_reps_skipped_included(self):
        self.assertTrue(self.conn is not None)

        def excluded_reps(conn):
            cur = conn.cursor()
            cur.execute('SELECT replicate FROM replicates WHERE includeRep == FALSE;')
            return [x[0] for x in cur.fetchall()]

        def included_reps(conn):
            cur = conn.cursor()
            cur.execute('SELECT replicate FROM replicates WHERE includeRep == TRUE;')
            return [x[0] for x in cur.fetchall()]

        test_conn = sqlite3.connect(':memory:')
        try:
            self.conn.backup(test_conn)

            included = included_reps(test_conn)
            self.assertTrue(len(included), self.N_REPS)

            # Select single rep to exclude
            random.seed(1)
            single_rep = included.pop(random.randint(0, self.N_REPS - 1))
            self.assertTrue(self.mark_reps_skipped(test_conn, reps=[single_rep]))
            self.assertEqual(len(excluded_reps(test_conn)), 1)
            self.assertEqual(len(included_reps(test_conn)), self.N_REPS - 1)
            self.assertTrue(single_rep in excluded_reps(test_conn))

            # Select multiple reps to exclude
            excluded = [included[i] for i in {random.randint(0, self.N_REPS - 2) for _ in range(10)}]
            included = [x for x in included if x not in excluded]
            self.assertTrue(self.mark_reps_skipped(test_conn, reps=excluded))
            excluded.append(single_rep)
            db_excluded = excluded_reps(test_conn)
            self.assertEqual(len(db_excluded), len(excluded))
            self.assertEqual(len(included_reps(test_conn)), len(included))
            self.assertTrue(all(rep in db_excluded for rep in excluded + [single_rep]))

            # Set all reps included
            self.mark_all_reps_included(test_conn)
            included = included_reps(test_conn)
            self.assertEqual(len(included), self.N_REPS)
            self.assertEqual(len(excluded_reps(test_conn)), 0)

            # Exclude project
            self.assertTrue(self.mark_reps_skipped(test_conn, projects=[self.TEST_PROJECT]))
            self.assertEqual(len(excluded_reps(test_conn)), self.N_REPS)
            self.assertEqual(len(included_reps(test_conn)), 0)

            # Set all reps included again
            self.mark_all_reps_included(test_conn)
            self.assertEqual(len(included_reps(test_conn)), self.N_REPS)
            self.assertEqual(len(excluded_reps(test_conn)), 0)

        finally:
            test_conn.close()


    def test_metadata_functions(self):
        self.assertTrue(self.conn is not None)

        # validate current version
        self.assertTrue(db_utils.check_schema_version(self.conn))

        # test get_meta_value
        schema_version = db_utils.get_meta_value(self.conn, 'schema_version')
        self.assertEqual(db_utils.SCHEMA_VERSION, schema_version)
        dia_qc_version = db_utils.get_meta_value(self.conn, 'dia_qc version')
        self.assertEqual(PROGRAM_VERSION, dia_qc_version)

        test_conn = sqlite3.connect(':memory:')
        try:
            self.conn.backup(test_conn)

            # set dia_qc version to NULL and check that we get a warning.
            db_utils.update_meta_value(test_conn, 'dia_qc version', 'NULL')
            with self.assertLogs(db_utils.LOGGER, level='WARNING') as cm:
                self.assertTrue(db_utils.check_schema_version(test_conn))
            message = f'The database was created with dia_qc version NULL but the current version is {PROGRAM_VERSION}'
            self.assertEqual(len(cm.output), 1)
            self.assertTrue(message in cm.output[0])
            self.assertEqual(db_utils.get_meta_value(test_conn, 'dia_qc version'), 'NULL')
            self.assertEqual(len(cm.output), 1)

            # set dia_qc version back to correct version
            db_utils.update_meta_value(test_conn, 'dia_qc version', dia_qc_version)
            self.assertEqual(db_utils.get_meta_value(test_conn, 'dia_qc version'), dia_qc_version)
            with self.assertNoLogs(db_utils.LOGGER) as cm:
                self.assertTrue(db_utils.check_schema_version(test_conn))

            # set current version to 'NULL' and make sure validation fails
            db_utils.update_meta_value(test_conn, 'schema_version', 'NULL')
            with self.assertLogs(db_utils.LOGGER, level='ERROR') as cm:
                self.assertFalse(db_utils.check_schema_version(test_conn))
            self.assertEqual(db_utils.get_meta_value(test_conn, 'schema_version'), 'NULL')
            self.assertEqual(len(cm.output), 1)
            self.assertTrue(f'Database schema version (NULL) does not match program ({db_utils.SCHEMA_VERSION})' in cm.output[0])

        finally:
            test_conn.close()


    def test_is_normalized(self):
        self.assertTrue(self.conn is not None)

        self.assertFalse(db_utils.is_normalized(self.conn))
        normalize_result = setup_functions.run_main(
            normalize_db._main, [self.db_path], self.work_dir,
            prefix='normalize_single_proj', prog='dia_qc normalize'
        )
        self.assertEqual(normalize_result.returncode, 0, normalize_result.stderr)
        self.assertTrue(db_utils.is_normalized(self.conn))


    def test_update_metadata_types(self):
        self.assertTrue(self.conn is not None)

        def get_dtypes():
            cur = test_conn.cursor()
            cur.execute('SELECT annotationKey, annotationType FROM sampleMetadataTypes;')
            return dict((x[0], x[1]) for x in cur.fetchall())

        test_conn = sqlite3.connect(':memory:')
        try:
            self.conn.backup(test_conn)


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
            db_utils.update_metadata_dtypes(test_conn, new_types)

            result_types = {'string_var': 'STRING',
                            'bool_var': 'STRING',
                            'int_var': 'FLOAT',
                            'float_var': 'FLOAT',
                            'na_var': 'NULL'}
            db_new_types = get_dtypes()

        finally:
            test_conn.close()

        self.assertDictEqual(result_types, db_new_types)


    def test_make_gene_matrix_by_protein_fails(self):
        self.assertTrue(self.conn is not None)

        self.assertEqual('protein', db_utils.get_meta_value(self.conn, 'group_precursors_by'))

        gene_id_path = f'{self.data_dir}/metadata/prhuman2gene_2023_05_24_subset.csv'
        result = setup_functions.run_main(
            export_gene_matrix._main, [gene_id_path, self.db_path], self.work_dir,
            prog='dia_qc export_gene_matrix'
        )
        self.assertEqual(result.returncode, 1)
        self.assertTrue('Precursors in database must be grouped by gene!' in result.stderr)


    def test_reset_imputed_values(self):
        test_conn = sqlite3.connect(':memory:')

        random.seed(7)
        try:
            self.conn.backup(test_conn)

            # count number of distinct proteins in db
            cur = test_conn.cursor()
            cur.execute('SELECT DISTINCT proteinId FROM proteins;')
            n_proteins = len(cur.fetchall())

            # add fake imputation metadata
            impute_metadata = {db_utils.PRECURSOR_IMPUTE_METHOD: 'KNN',
                               db_utils.PROTEIN_IMPUTE_METHOD: 'KNN'}
            for k, v in impute_metadata.items():
                db_utils.update_meta_value(test_conn, k, v)

            # make sure metadata was updated
            cur = test_conn.cursor()
            cur.execute('SELECT key FROM metadata;')
            metadata = {x[0] for x in cur.fetchall()}
            for key in impute_metadata:
                self.assertIn(key, metadata)

            # randomly select proteins and precursors
            df = pd.read_sql('''
                SELECT p.replicateId, p.peptideId, ptp.proteinId, p.precursorCharge
                FROM precursors p
                LEFT JOIN peptideToProtein ptp ON ptp.peptideId = p.peptideId;''', test_conn)

            null_proteins = random.sample(df['proteinId'].drop_duplicates().to_list(),
                                          int(n_proteins * 0.2))
            null_reps = random.sample(df['replicateId'].drop_duplicates().to_list(),
                                      self.N_REPS - 1)

            null_protein_ids = {(r, p) for r in null_reps for p in null_proteins}
            mask = df.apply(lambda x: x['proteinId'] in null_proteins and
                                      x['replicateId'] in null_reps,
                            axis=1)
            null_prec_ids = {(row.replicateId, row.peptideId) for row in df[mask].itertuples()}

            # update fake isImputed values in database
            cur = test_conn.cursor()
            cur.executemany('''
                UPDATE precursors
                    SET isImputed = TRUE
                WHERE replicateId = ? AND peptideId = ?''', null_prec_ids)
            cur.executemany('''
                UPDATE proteinQuants
                    SET isImputed = TRUE
                WHERE replicateId = ? AND proteinId = ?''', null_protein_ids)
            test_conn.commit()

            # make sure the fake imputed values were set
            cur.execute('SELECT replicateId, peptideId, isImputed, totalAreaFragment FROM precursors;')
            precursor_data = {(r, p): (i, a) for r, p, i, a in cur.fetchall()}
            for key in null_prec_ids:
                self.assertEqual(precursor_data[key][0], 1)
                self.assertIsNotNone(precursor_data[key][1])
            cur.execute('SELECT replicateId, proteinId, isImputed, abundance FROM proteinQuants;')
            protein_data = {(r, p): (i, a) for r, p, i, a in cur.fetchall()}
            for key in null_protein_ids:
                self.assertEqual(protein_data[key][0], 1)
                self.assertIsNotNone(protein_data[key][1])

            # run reset_imputed_values on test db
            with self.assertLogs(logger.LOGGER, level='INFO') as cm:
                db_utils.reset_imputed_values(test_conn)
            self.assertTrue(all(re.search(r'Setting imputed (protein|precursor) quantities to NULL\.', out)
                                for out in cm.output))

            # make sure the metadata was reset
            cur = test_conn.cursor()
            cur.execute('SELECT key FROM metadata;')
            metadata = {x[0] for x in cur.fetchall()}
            for key in impute_metadata:
                self.assertNotIn(key, metadata)

            # make sure the fake imputed values were reset
            cur.execute('SELECT replicateId, peptideId, isImputed, totalAreaFragment FROM precursors;')
            precursor_data = {(r, p): (i, a) for r, p, i, a in cur.fetchall()}
            for key in null_prec_ids:
                self.assertEqual(precursor_data[key][0], 0)
                self.assertIsNone(precursor_data[key][1])
            cur.execute('SELECT replicateId, proteinId, isImputed, abundance FROM proteinQuants;')
            protein_data = {(r, p): (i, a) for r, p, i, a in cur.fetchall()}
            for key in null_protein_ids:
                self.assertEqual(protein_data[key][0], 0)
                self.assertIsNone(protein_data[key][1])

        finally:
            test_conn.close()


class TestCommandLog(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_tables = ['commandLog']
        cls.test_schema = list()
        for command in db_utils.SCHEMA:
            for table in cls.test_tables:
                if command.strip().startswith(f'CREATE TABLE {table}'):
                    cls.test_schema.append(command)
                    continue


    def setUp(self):
        self.assertEqual(len(self.test_tables), len(self.test_schema))
        self.conn = sqlite3.connect('file::memory:?cache=shared', uri=True)
        cur = self.conn.cursor()
        for command in self.test_schema:
            cur.execute(command)
        self.conn.commit()


    def tearDown(self):
        self.conn.close()


    def test_add_entry(self):
        command = ['command', '-a', '123']
        db_utils.update_command_log(self.conn, command, os.getcwd())

        cur = self.conn.cursor()
        cur.execute('''SELECT commandNumber, command, version, workingDirectory, time, user, hostname
                       FROM commandLog;''')
        index, db_command, version, wd, time, user, host = cur.fetchall()[0]

        self.assertEqual(index, 1)
        self.assertEqual(db_command, ' '.join(command))
        self.assertEqual(version, PROGRAM_VERSION)
        self.assertEqual(user, os.environ.get('USER', os.environ.get('USERNAME', 'UNKNOWN')))
        self.assertEqual(host, gethostname())
        self.assertTrue(os.path.isdir(wd))
        self.assertLess(datetime.strptime(time, db_utils.METADATA_TIME_FORMAT),
                        datetime.now())


    def test_get_last_command(self):
        commands = ['abc', 'bca', 'cab', 'aabc']
        for command in commands:
            db_utils.update_command_log(self.conn, [command], os.getcwd())
            self.assertEqual(db_utils.get_last_command(self.conn), [command])


class TestUpdateAcquiredRanks(unittest.TestCase, TestDBFunctionsBase):
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

        # find replicate and metadata table commands in db_utils.SCHEMA
        cls.test_schema = list()
        cls.test_tables = ['replicates', 'metadata']
        for command in db_utils.SCHEMA:
            for table in cls.test_tables:
                if command.strip().startswith(f'CREATE TABLE {table}'):
                    cls.test_schema.append(command)
                    continue


    def setUp(self):
        self.assertEqual(len(self.test_tables), len(self.test_schema))

        # create in memory database
        self.conn = sqlite3.connect('file::memory:?cache=shared', uri=True)
        cur = self.conn.cursor()
        for command in self.test_schema:
            cur.execute(command)

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
                                index=True, index_label='id')
        db_utils.update_acquired_ranks(self.conn)

        db_ranks = self.get_db_ranks(self.conn)
        self.assertDictEqual(self.acquired_ranks, db_ranks)


    def test_add_project(self):

        # add first project
        self.replicates[self.replicates['project'] == 'Sp3'].to_sql('replicates', self.conn,
                                                                    if_exists='append', index=True,
                                                                    index_label='id')
        db_utils.update_acquired_ranks(self.conn)

        db_ranks = self.get_db_ranks(self.conn)
        proj_ranks = self.update_ground_truth_ranks(db_ranks.keys())
        self.assertDictEqual(proj_ranks, db_ranks)

        # add second project
        self.replicates[self.replicates['project'] == 'Strap'].to_sql('replicates', self.conn,
                                                                      if_exists='append', index=True,
                                                                      index_label='id')

        db_utils.update_acquired_ranks(self.conn)
        db_ranks = self.get_db_ranks(self.conn)
        self.assertDictEqual(self.acquired_ranks, db_ranks)


    def test_mark_project_skipped(self):
        self.replicates.to_sql('replicates', self.conn, if_exists='append',
                                index=True, index_label='id')
        db_utils.update_acquired_ranks(self.conn)

        db_ranks = self.get_db_ranks(self.conn)
        self.assertDictEqual(self.acquired_ranks, db_ranks)

        self.mark_reps_skipped(self.conn, projects=('Strap',))
        db_utils.update_acquired_ranks(self.conn)

        db_ranks = self.get_db_ranks(self.conn)
        proj_ranks = self.update_ground_truth_ranks(db_ranks.keys())
        self.assertDictEqual(proj_ranks, db_ranks)


    def test_mark_replicates_skipped(self):
        self.replicates.to_sql('replicates', self.conn, if_exists='append',
                                index=True, index_label='id')
        db_utils.update_acquired_ranks(self.conn)
        random.seed(1)
        max_rep_index = len(self.replicates.index) - 1
        max_subset_len = 15

        for _ in range(50):
            skip_reps = list()
            for i in random.sample(range(max_rep_index), random.randint(1, max_subset_len)):
                skip_reps.append(self.replicates.loc[i, 'replicate'])

            self.assertTrue(self.mark_reps_skipped(self.conn, reps=skip_reps))

            db_utils.update_acquired_ranks(self.conn)

            db_ranks = self.get_db_ranks(self.conn)
            proj_ranks = self.update_ground_truth_ranks(db_ranks.keys())
            self.assertDictEqual(proj_ranks, db_ranks)

            with self.assertLogs(logger.LOGGER, level='INFO'):
                db_utils.mark_all_reps_included(self.conn)


class TestReadWideMetadata(unittest.TestCase, setup_functions.AbstractTestsBase):
    TEST_PROJECT = 'HeLa'

    @staticmethod
    def df_to_dict(df, cols):
        ret = {row['replicate']: {var: row[var] for var in cols}
               for _, row in df.iterrows()}
        return ret

    @classmethod
    def setUpClass(cls):
        # read replicates
        rep_path = f'{setup_functions.TEST_DIR}/data/skyline_reports/{cls.TEST_PROJECT}_replicate_quality.tsv'
        cls.replicates = ReplicateReport(quiet=True).read_report(rep_path)
        cls.replicates = cls.replicates[[col.name for col in ReplicateReport().required_columns()]]
        cls.replicates['project'] = cls.TEST_PROJECT
        cls.replicates['acquiredRank'] = -1
        rep_index = {r: i for i, r in zip(cls.replicates['replicate'].index,
                                          cls.replicates['replicate'])}

        # read metadata
        cls.metadata_path = f'{setup_functions.TEST_DIR}/data/metadata/{cls.TEST_PROJECT}_metadata.tsv'
        metadata_df = pd.read_csv(cls.metadata_path, sep='\t')

        # add replicate data to cls.metadat_df
        metadata_df = metadata_df.rename(columns={'Replicate': 'replicate'})
        metadata_df['replicateId'] = metadata_df['replicate'].apply(lambda x: rep_index[x])
        metadata_df = metadata_df.set_index('replicateId')
        metadata_df = metadata_df.join(cls.replicates[['project']])

        for col in metadata_df.columns:
            if all(pd.isna(metadata_df[col])):
                metadata_df[col] = None

        metadata_df = metadata_df.reset_index()

        keep_cols = [col for col in metadata_df.columns if col != 'replicate']
        cls.metadata_dict = cls.df_to_dict(metadata_df, keep_cols)

        # construct test db schema
        cls.test_tables = ['replicates', 'sampleMetadata ', 'sampleMetadataTypes', 'metadata']
        cls.test_schema = list()
        for command in db_utils.SCHEMA:
            for table in cls.test_tables:
                if command.strip().startswith(f'CREATE TABLE {table}'):
                    cls.test_schema.append(command)
                    continue


    def setUp(self):
        # read metadata
        meta_reader = Metadata()
        self.assertTrue(meta_reader.read(self.metadata_path))
        metadata = meta_reader.df

        # add replicateId column
        rep_index = {r: i for i, r in zip(self.replicates['replicate'].index,
                                          self.replicates['replicate'])}
        metadata['replicateId'] = metadata['Replicate'].apply(lambda x: rep_index[x])
        metadata = metadata[['replicateId', 'annotationKey', 'annotationValue']]

        # create in memory database
        self.assertEqual(len(self.test_tables), len(self.test_schema))
        self.conn = sqlite3.connect('file::memory:?cache=shared', uri=True)
        cur = self.conn.cursor()
        for command in self.test_schema:
            cur.execute(command)
        self.conn.commit()

        # write data to database
        self.replicates.to_sql('replicates', self.conn, if_exists='append', index=True, index_label='id')
        metadata.to_sql('sampleMetadata', self.conn, if_exists='append', index=False)
        db_utils.update_metadata_dtypes(self.conn, meta_reader.types)
        db_utils.update_acquired_ranks(self.conn)


    def tearDown(self):
        self.conn.close()


    @staticmethod
    def subset_data_dict(d, selection):
        ret = dict()
        for rep in d:
            ret[rep] = {k: v for k, v in d[rep].items() if k in selection}

        return ret


    def test_read_all(self):
        metadata = db_utils.read_wide_metadata(self.conn)

        skip_cols = {'replicate', 'acquiredRank'}
        keep_cols = [col for col in metadata.columns if col not in skip_cols]

        metadata_dict = self.df_to_dict(metadata, keep_cols)

        self.assertDataDictEqual(metadata_dict, self.metadata_dict)


    def test_only_acquired_ranks(self):
        metadata = db_utils.read_wide_metadata(self.conn, read_all=False)

        keep_cols = {'replicateId', 'project'}
        metadata_dict = self.df_to_dict(metadata, keep_cols)

        # select keep_cols from ground truth dict
        test_df = self.subset_data_dict(self.metadata_dict, keep_cols)

        self.assertDataDictEqual(metadata_dict, test_df)


    def test_select_vars(self):
        keep_cols = set()
        for d in self.metadata_dict.values():
            keep_cols |= set(d.keys())

        selections = list()
        for n in range(1, len(keep_cols) + 1):
            for sele in combinations(keep_cols, n):
                selections.append(sele)

        for sele in selections:
            metadata = db_utils.read_wide_metadata(self.conn, meta_vars=sele, read_all=False)

            metadata_dict = self.df_to_dict(metadata, sele)

            self.assertDataDictEqual(metadata_dict,
                                     self.subset_data_dict(self.metadata_dict, sele))
