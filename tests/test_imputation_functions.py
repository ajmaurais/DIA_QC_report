
import unittest
import os
import sqlite3
import pandas as pd
import random

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
    key_cols:
    '''
    replicates = df['replicate'].drop_duplicates().to_list()

    def add_nas(d):
        subset = set(random.sample(replicates, random.randint(min_missing, max_missing)))
        d.loc[d[replicate_col].apply(lambda x: x in subset), value_col] = pd.NA
        return d.reset_index(drop=True)

    random.seed(seed)
    df_s = df.copy()
    df_s = df_s.groupby(key_cols).apply(add_nas, include_groups=False)
    return df_s.reset_index()



class TestKNNImputeDF(unittest.TestCase, setup_functions.AbstractTestsBase):
    def setUp(self):
        fname = f'{setup_functions.TEST_DIR}/data/skyline_reports/Strap_by_protein_precursor_quality.csv'
        df = PrecursorReport(quiet=True).read_report(fname)
        df = df[['replicateName', 'modifiedSequence', 'precursorCharge', 'totalAreaFragment']].drop_duplicates()
        self.df = df.rename(columns={'replicateName': 'replicate', 'totalAreaFragment': 'area'})


    def test_no_imputation_when_no_missing(self):
        df = self.df.copy()
        df = imputation.knn_impute_df(df, ['modifiedSequence', 'precursorCharge'],
                                      'area', replicate_col='replicate')
        self.assertTrue(all(~df['isImputed']))
        df = df.drop(columns=['isImputed'])
        self.assertEqual(set(self.df.columns), set(df))
        self.assertDataFrameEqual(self.df, df[self.df.columns.to_list()])


    def test_imputed_values_set(self):
        replicates = self.df['replicate'].drop_duplicates().to_list()

        df = create_missing_values(self.df, ['modifiedSequence', 'precursorCharge'],
                                   value_col='area', seed=1,
                                   min_missing=1, max_missing=len(replicates) - 1)

        self.assertTrue(any(pd.isna(df['area'])))

        df_i = imputation.knn_impute_df(df, ['modifiedSequence', 'precursorCharge'],
                                        'area', replicate_col='replicate', missing_threshold=None)

        df_i = df_i.set_index(['replicate', 'modifiedSequence', 'precursorCharge'])
        df = df.set_index(['replicate', 'modifiedSequence', 'precursorCharge'])
        df_j = df.join(df_i, rsuffix='_imputed')

        for row in df_j.itertuples():
            if pd.isna(row.area):
                self.assertTrue(row.isImputed)
                self.assertFalse(pd.isna(row.area_imputed))
            else:
                self.assertFalse(row.isImputed)


    def test_threshold_NAs_skipped(self):
        pass


class TestImputationBase(setup_functions.AbstractTestsBase):
    def __init__(self):
        self.work_dir = None
        self.db_path = None
        self.data_dir = None
        self.parse_result = None
        self.conn = None


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    def test_abc_instantiation_fails(self):
        self.assertIsNotNone(self.conn)

        with self.assertRaises(TypeError):
            imputation.ImputationManagerBase(self.conn)


class TestSingleImputation(unittest.TestCase, TestImputationBase):
    TEST_PROJECT = 'Sp3'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_imputation_functions/'
        cls.db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        cls.parse_result = setup_functions.setup_single_db(cls.data_dir,
                                                           cls.work_dir,
                                                           cls.TEST_PROJECT,
                                                           clear_dir=True)

        cls.conn = None
        if cls.parse_result.returncode == 0:
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

