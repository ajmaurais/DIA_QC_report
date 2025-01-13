
import unittest
import os
import sqlite3
import random
import pandas as pd

import setup_functions

import DIA_QC_report.submodules.imputation as imputation
import DIA_QC_report.submodules.dia_db_utils as db_utils
from DIA_QC_report.submodules.skyline_reports import PrecursorReport


def create_missing_values(df, key_cols,
                          replicate_col='replicate', value_col='area',
                          min_missing=0, max_missing=1, seed=1):
    '''
    Add NAs to random replicates in df

    Parameters
    ----------
    df: pd.DataFrame
        A long formated dataframe with peptide/protein quantities
    key_cols: list
        A list of key columns which uniquely identify each peptide or protein
    replicate_col: str
        The name of the replicate colum.
    value_col: str
        The name of the value column.
    min_missing: int
        Mininum number of missing values.
    max_missing: int
        Maximum number of missing values.
    seed: int
        The seed for the random number generator.

    Returns
    -------
    df_na: pd.DataFrame
        A long formated dataframe the same shape as `df` with random values in
        `value_col` set to pd.NA.
    '''
    replicates = df['replicate'].drop_duplicates().to_list()

    def add_nas(d):
        n = random.randint(min_missing, max_missing)
        subset = set(random.sample(replicates, n))
        d.loc[d[replicate_col].apply(lambda x: x in subset), value_col] = pd.NA
        return d.reset_index(drop=True)

    random.seed(seed)
    df_s = df.copy()
    df_s['area'] = df_s.groupby(key_cols).apply(add_nas, include_groups=False).reset_index(drop=True)['area']
    return df_s


class TestKNNImputeDF(unittest.TestCase, setup_functions.AbstractTestsBase):
    PRECURSOR_KEY_COLS = ['modifiedSequence', 'precursorCharge']
    def setUp(self):
        fname = f'{setup_functions.TEST_DIR}/data/skyline_reports/Strap_by_protein_precursor_quality.csv'
        df = PrecursorReport(quiet=True).read_report(fname)
        df = df[self.PRECURSOR_KEY_COLS + ['replicateName', 'totalAreaFragment']].drop_duplicates()
        self.df = df.rename(columns={'replicateName': 'replicate', 'totalAreaFragment': 'area'})
        self.n_reps = len(self.df['replicate'].drop_duplicates())


    def test_no_imputation_when_no_missing(self):
        df = self.df.copy()
        df = imputation.knn_impute_df(df, self.PRECURSOR_KEY_COLS,
                                      'area', replicate_col='replicate')
        self.assertTrue(all(~df['isImputed']))
        df = df.drop(columns=['isImputed'])
        self.assertEqual(set(self.df.columns), set(df))
        self.assertDataFrameEqual(self.df, df[self.df.columns.to_list()])


    def test_create_missing_values(self):
        key_cols = self.PRECURSOR_KEY_COLS

        tests = [(1, self.n_reps - 1), (0, self.n_reps - 1), (5, 5), (5, self.n_reps - 1)]
        for i, test in enumerate(tests):
            df = create_missing_values(self.df, self.PRECURSOR_KEY_COLS,
                                       value_col='area', seed=i,
                                       min_missing=test[0], max_missing=test[1])
            n_nas = df.groupby(key_cols)['area'].apply(lambda x: x.isna().sum())

            self.assertGreaterEqual(min(n_nas), test[0])
            self.assertLessEqual(max(n_nas), test[1])


    def test_imputed_values_set(self):
        key_cols = self.PRECURSOR_KEY_COLS

        df = create_missing_values(self.df, key_cols,
                                   value_col='area', seed=1,
                                   min_missing=0, max_missing=self.n_reps - 1)

        self.assertTrue(any(pd.isna(df['area'])))

        df_i = imputation.knn_impute_df(df, key_cols, 'area',
                                        replicate_col='replicate', max_missing=None)

        df_i = df_i.set_index(['replicate'] + key_cols)
        df = df.set_index(['replicate'] + key_cols)
        df_j = df.join(df_i, rsuffix='_imputed')

        for _, group in df_j.groupby(key_cols):
            if any(group['area'].isna()):
                self.assertTrue(all(group['isImputed'] == group['area'].isna()))
                self.assertFalse(any(group['area_imputed'].isna()))
            else:
                self.assertFalse(any(group['isImputed']))


    def test_threshold_NAs_skipped(self):
        key_cols = self.PRECURSOR_KEY_COLS

        df = create_missing_values(self.df, key_cols,
                                   value_col='area', seed=1,
                                   min_missing=0, max_missing=self.n_reps - 1)

        self.assertTrue(any(pd.isna(df['area'])))

        for t in range(1, self.n_reps):
            df_i = imputation.knn_impute_df(df, key_cols, 'area',
                                            replicate_col='replicate', max_missing=t)

            df_i = df_i.set_index(['replicate'] + key_cols)
            df_j = df.set_index(['replicate'] + key_cols)
            df_j = df_j.join(df_i, rsuffix='_imputed')

            for _, group in df_j.groupby(key_cols):
                self.assertEqual(group.shape[0], self.n_reps)
                n_na = sum(group['area'].isna())

                if n_na == 0:
                    self.assertFalse(any(group['isImputed']))
                    self.assertTrue(all(group['area'] == group['area_imputed']))
                elif n_na >= t:
                    self.assertFalse(any(group['isImputed']))
                    self.assertEqual(sum(group['area_imputed'].isna()), n_na)
                else:
                    self.assertFalse(any(group['area_imputed'].isna()))


class TestImputationBase(setup_functions.AbstractTestsBase):
    def __init__(self):
        self.work_dir = None
        self.db_path = None
        self.data_dir = None
        self.conn = None


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    @staticmethod
    def set_random_nulls(conn, min_missing=0, max_missing=1, seed=1):
        '''
        Set quantities for random precursors and proteins in conn to NA.

        Parameters
        ----------
        conn: sqlite3.Connection
            A connection to a batch database.
        min_missing: int
            Mininum number of missing values.
        max_missing: int
            Maximum number of missing values.
        seed: int
            The seed for the random number generator.
        '''
        random.seed(seed)

        cur = conn.cursor()
        cur.execute('SELECT replicateId, peptideId, precursorCharge FROM precursors;')
        precursors = dict()
        for rep, peptide, charge in cur.fetchall():
            p = (peptide, charge)
            if p not in precursors:
                precursors[p] = list()
            precursors[p].append(rep)

        precursor_nas = list()
        for precursor, reps in precursors.items():
            n = len(reps) - random.randint(min_missing, max_missing)
            for r in random.sample(reps, n):
                precursor_nas.append((r, precursor[0], precursor[1]))

        cur.executemany('''
            UPDATE precursors
                SET totalAreaFragment = NULL, normalizedArea = NULL
            WHERE replicateId = ? AND
                  peptideId = ? AND
                  precursorCharge = ? ;''',
                        precursor_nas)
        conn.commit()

        df_pre = pd.read_sql('''
            SELECT
                p.replicateId, ptp.proteinId, p.peptideId, p.precursorCharge, p.totalAreaFragment as area
            FROM precursors p
            LEFT JOIN peptideToProtein ptp ON ptp.peptideId = p.peptideId;''', conn)

        df_prot = df_pre.groupby(['replicateId', 'proteinId'])['area'].sum(min_count=1).reset_index()

        protein_nas = [(row.replicateId, row.proteinId)
                       for row in df_prot[df_prot['area'].isna()].itertuples()]

        cur = conn.cursor()
        cur.executemany('''
            UPDATE proteinQuants
                SET abundance = NULL, normalizedAbundance = NULL
            WHERE replicateId = ? AND
                  proteinId = ? ;''',
                        protein_nas)
        conn.commit()


    def test_abc_instantiation_fails(self):
        self.assertIsNotNone(self.conn)

        with self.assertRaises(TypeError):
            imputation.ImputationManagerBase(self.conn)


    def test_db_has_missing_values(self):
        '''
        Test that test database was initialized properly.
        There should be missing precursor and protein quantities.
        '''
        self.assertIsNotNone(self.conn)

        df_pre = pd.read_sql('''
            SELECT
                p.replicateId, ptp.proteinId, p.peptideId, p.precursorCharge, p.totalAreaFragment as area
            FROM precursors p
            LEFT JOIN peptideToProtein ptp ON ptp.peptideId = p.peptideId;''', self.conn)

        protein_keys = ['replicateId', 'proteinId']
        na_precursors = df_pre.groupby(protein_keys)['area'].apply(lambda x: all(pd.isna(x)))
        na_precursors.name = 'all_na'
        self.assertTrue(any(na_precursors))

        df_prot = pd.read_sql('''
            SELECT
                q.replicateId, q.proteinId, q.abundance
            FROM proteinQuants q; ''', self.conn)

        self.assertTrue(any(df_prot['abundance'].isna()))

        df_prot = df_prot.set_index(protein_keys).join(na_precursors)
        self.assertFalse(any(df_prot['all_na'].isna()))
        self.assertSeriesEqual(df_prot['all_na'], df_prot['abundance'].isna(),
                               check_names=False)


    def do_level_test(self, level, df_precursor, df_protein):
        manager = imputation.KNNImputer(self.conn, level=level)
        with self.assertNoLogs(imputation.LOGGER, level='WARNING') as cm:
            self.assertTrue(manager.impute())

        db_pre = manager.precursors
        db_pre = db_pre[~db_pre['isImputed']]
        db_pre = db_pre.set_index(df_precursor.index.names)

        db_pro = manager.proteins
        db_pro = db_pro[~db_pro['isImputed']]
        db_pro = db_pro.set_index(df_protein.index.names)

        db_pre = db_pre.join(df_precursor)
        db_pro = db_pro.join(df_protein)

        self.assertSeriesEqual(db_pre['area'],
                               db_pre['quant' if level == 0 else 'normQuant'],
                               check_names=False, check_exact=True)
        self.assertSeriesEqual(db_pro['abundance'],
                               db_pro['quant' if level == 0 else 'normQuant'],
                               check_names=False, check_exact=True)


    def test_imputation_level(self):
        self.assertIsNotNone(self.conn)

        df_pre = pd.read_sql('''SELECT replicateId, peptideId, precursorCharge,
                                totalAreaFragment as quant, normalizedArea as normQuant
                                FROM precursors;''', self.conn)
        df_pre = df_pre.set_index(['replicateId', 'peptideId', 'precursorCharge'])

        df_pro = pd.read_sql('''SELECT replicateId, proteinId,
                                abundance as quant, normalizedAbundance as normQuant
                                FROM proteinQuants;''', self.conn)
        df_pro = df_pro.set_index(['replicateId', 'proteinId'])

        self.do_level_test(0, df_pre, df_pro)
        self.do_level_test(1, df_pre, df_pro)


    def test_all_reps_skipped(self):
        self.assertIsNotNone(self.conn)

        manager = imputation.KNNImputer(self.conn)
        try:
            cur = self.conn.cursor()
            cur.execute('UPDATE replicates SET includeRep = FALSE;')
            self.conn.commit()

            with self.assertLogs(imputation.LOGGER, level='ERROR') as cm:
                self.assertFalse(manager._read_precursors())
            self.assertTrue(any('All replicates in database have been excluded!' in entry for entry in cm.output))

            with self.assertLogs(imputation.LOGGER, level='ERROR') as cm:
                self.assertFalse(manager.impute())
            self.assertTrue(any('All replicates in database have been excluded!' in entry for entry in cm.output))

        finally:
            db_utils.mark_all_reps_included(self.conn, quiet=True)


class TestSingleImputation(unittest.TestCase, TestImputationBase):
    TEST_PROJECT = 'Sp3'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_imputation_functions_single/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                       cls.work_dir,
                                                       cls.TEST_PROJECT,
                                                       clear_dir=True)

        cls.conn = None
        if parse_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                conn = sqlite3.connect(cls.db_path)

                cur = conn.cursor()
                cur.execute('SELECT id FROM replicates;')
                n_reps = len(cur.fetchall())

                # set db quantities randomly to NULL
                TestImputationBase.set_random_nulls(conn, min_missing=0, max_missing=n_reps - 1)
                conn.close()

        else:
            return

        normalize_command = ['dia_qc', 'normalize', '-m=median', '--keepMissing', cls.db_path]
        normalize_result = setup_functions.run_command(normalize_command,
                                                       cls.work_dir,
                                                       prefix='normalize')

        if normalize_result.returncode == 0:
            if os.path.isfile(cls.db_path):
                cls.conn = sqlite3.connect(cls.db_path)


    def test_read_precursors(self):
        self.assertIsNotNone(self.conn)

        manager = imputation.KNNImputer(self.conn)

        with self.assertNoLogs(imputation.LOGGER, level='WARNING') as cm:
            self.assertTrue(manager._read_precursors())


    def test_read_proteins(self):
        self.assertIsNotNone(self.conn)

        manager = imputation.KNNImputer(self.conn)

        with self.assertNoLogs(imputation.LOGGER, level='WARNING') as cm:
            self.assertTrue(manager._read_proteins())

