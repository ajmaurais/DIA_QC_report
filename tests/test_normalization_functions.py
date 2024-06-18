
import unittest
from unittest import mock
from abc import ABC, abstractmethod
import os
import sqlite3
import pandas as pd
from numpy import log2

import setup_functions

import DIA_QC_report.submodules.normalization as normalization
import DIA_QC_report.submodules.dia_db_utils as db_utils


class TestCammelCase(unittest.TestCase):
    def test(self):
        self.assertEqual(normalization.cammel_case('log2', 'area'), 'log2Area')
        norm_value = normalization.cammel_case('normalized', 'area')
        self.assertEqual(norm_value, 'normalizedArea')
        self.assertEqual(normalization.cammel_case('log2', norm_value), 'log2NormalizedArea')


class TestNormalizationBase(ABC):
    def __init__(self):
        self.work_dir = None
        self.db_path = None
        self.data_dir = None
        self.parse_result = None
        self.conn = None


    @classmethod
    @abstractmethod
    def setUpClass(cls):
        pass


    @classmethod
    def tearDownClass(cls):
        if cls.conn is not None:
            cls.conn.close()


    @abstractmethod
    def assertIsNotNone(self, expr):
        pass


    @abstractmethod
    def assertTrue(self, expr):
        pass


    @abstractmethod
    def assertFalse(self, expr):
        pass


    @abstractmethod
    def assertRaises(self, exception):
        pass


    @abstractmethod
    def assertLogs(self, logger, level=None):
        pass


    def test_abc_instantiation_fails(self):
        self.assertIsNotNone(self.conn)

        with self.assertRaises(TypeError):
            normalization.NormalizationManagerBase(self.conn)


    def test_median_normalization(self):
        self.assertIsNotNone(self.conn)

        # import pudb
        # pudb.set_trace()

        manager = normalization.MedianNormalizer(self.conn)
        with self.assertLogs(normalization.LOGGER) as cm:
            manager.normalize()

        precursors, proteins = manager.get_long_tables()

        def medians_equal(df, quant_col):
            norm_col = normalization.cammel_case('normalized', quant_col)
            log2_norm_col = normalization.cammel_case('log2', norm_col)

            medians = df.groupby('replicateId')[log2_norm_col].median()
            self.assertTrue(all(medians.iloc[0] == medians))

        medians_equal(precursors, 'area')
        medians_equal(proteins, 'abundance')


    def test_read_precursors(self):
        self.assertIsNotNone(self.conn)

        manager = normalization.MedianNormalizer(self.conn)

        with self.assertLogs(normalization.LOGGER) as cm:
            self.assertTrue(manager._read_precursors())
        self.assertTrue(any('Reading precursors from database...' in entry for entry in cm.output))
        self.assertTrue(any('Finished reading precursors.' in entry for entry in cm.output))


    @mock.patch('DIA_QC_report.submodules.dia_db_utils.LOGGER', mock.Mock())
    def test_all_reps_skipped(self):
        self.assertIsNotNone(self.conn)

        manager = normalization.MedianNormalizer(self.conn)

        try:
            cur = self.conn.cursor()
            cur.execute('UPDATE replicates SET includeRep = FALSE;')
            self.conn.commit()

            with self.assertLogs(normalization.LOGGER, level='ERROR') as cm:
                self.assertFalse(manager._read_precursors())
            self.assertTrue(any('All replicates in database have been excluded!' in entry for entry in cm.output))

            with self.assertLogs(normalization.LOGGER, level='ERROR') as cm:
                self.assertFalse(manager.normalize())
            self.assertTrue(any('All replicates in database have been excluded!' in entry for entry in cm.output))

        finally:
            db_utils.mark_all_reps_includced(self.conn)


class TestNormalization(unittest.TestCase, TestNormalizationBase):
    TEST_PROJECT = 'Sp3'

    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_normalization_functions/'
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


