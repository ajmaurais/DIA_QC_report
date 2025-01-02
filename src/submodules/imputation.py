
from abc import abstractmethod

import numpy as np
import pandas as pd

from .logger import LOGGER
from .transformation import TransformationManagerBase

IMPUTATION_METHODS = ('k-means',)


class ImputationManagerBase(TransformationManagerBase):
    '''
    ImputationManager abastract base class.

    Child classes must implement `impute` method.

    Attributes
    ----------
    conn: sqlite.Connection
        A connection to a precursor database
    precursors: pd.DataFrame
        Long formatted dataframe of precursor quantities.
    proteins: pd.DataFrame
        Long formatted dataframe of protein quantities.
    '''

    def __init__(self, conn,
                 level=0, missing_threshold=0.5,
                 impute_data='both', group_by_project=True):
        '''
        Parameters
        ----------
        conn: sqlite3.Connection
            Initialized connection to batch database.
        level: int
            Data level to impute. 0 for unnormalized, 1 for normalized.
            Default is 0.
        missing_threshold: float
            Mininum percent of non-missing values to impute value.
            Default is 0.5.
        impute_data: str
            Impute 'precursors', 'proteins', or 'both'. Default is 'both'.
        group_by_project: bool
            Only use replicates in the same project to 
            Default is 'project'.
        '''
        super().__init__(conn)

        if level in [0, 1]:
            self.level = level
        else:
            raise RuntimeError('level must be 0 or 1')
        self.missing_threshold = 0.5
        if impute_data in ['precursors', 'proteins', 'both']:
            self.impute_data = impute_data
        else:
            raise RuntimeError("impute_data must be 'precursors', 'proteins', or 'both'")
        

    def _read_precursors(self):
        # get precursor table from db
        LOGGER.info('Reading precursors from database...')
        df = pd.read_sql(f'''
            SELECT
                r.project,
                p.replicateId,
                p.peptideId,
                p.modifiedSequence,
                p.precursorCharge,
                p.isImputed,
                p.{'normalizedArea' if self.level==1 else 'totalAreaFragment'} as area
            FROM precursors p
            LEFT JOIN replicates r ON p.replicateId == r.id
            WHERE r.includeRep == TRUE; ''',
            self.conn)
        LOGGER.info('Finished reading precursors.')

        if len(df.index) == 0:
            LOGGER.error('All replicates in database have been excluded!')
            return False

        df['isImputed'] = df['isImputed'].map({0: False, 1: True})

        self.precursors = df
        return True


    def _read_proteins(self):
        LOGGER.info('Reading precursors from database...')
        df = pd.read_sql(f'''
            SELECT
                r.project,
                q.replicateId,
                q.proteinId,
                q.isImputed,
                q.{'normalizedAbundance' if self.level==1 else 'abundance'} as abundance
            FROM proteinQuants q
            LEFT JOIN replicates r ON q.replicateId == r.id
            WHERE r.includeRep == TRUE; ''',
            self.conn)
        LOGGER.info('Finished reading precursors.')

        if len(df.index) == 0:
            LOGGER.error('All replicates in database have been excluded!')
            return False
        
        df['isImputed'] = df['isImputed'].map({0: False, 1: True})
        # n_missing = df.group_by(['proteinId'])['abundance'].size()

        self.proteins = df
        return True


    @abstractmethod
    def impute(self):
        pass


class KMeansNormalizer(ImputationManagerBase):
    def __init__(self, conn, **kwargs):
        super().__init__(conn, **kwargs)


    def impute(self):
        if self.impute_data in ('precursors', 'both'):
            if not self._read_precursors():
                return False
        if self.impute_data in ('proteins', 'both'):
            if not self._read_proteins():
                return False
