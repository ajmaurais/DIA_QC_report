
from abc import abstractmethod

import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer as skl_KNNImputer

from .logger import LOGGER
from .transformation import cammel_case
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

        self.missing_threshold = missing_threshold
        self.group_by_project = group_by_project

        if impute_data in ('precursors', 'proteins', 'both'):
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
        LOGGER.info('Reading proteins from database...')
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

        self.proteins = df
        return True


    @abstractmethod
    def impute(self):
        pass


def knn_impute_df(df, key_cols, value_col, replicate_col='replicateId',
                  is_imputed_col='isImputed', max_missing=None,
                  n_neighbors=5, weights='uniform'):
    '''
    KNN impute peak areas in a long formated dataframe.

    Parameters
    ----------
    df: pd.DataFrame
        Long formated dataframe
    key_cols: list
        The names of column(s) which uniquely identify each row.
    replicate_col: str
        The name of the column that uniquely identify each replicate.
    value_col: str
        The name of the column with peak areas.
    is_imputed_col: str
        The name of the boolean column which indicates whether the value_col is imputed.
        If the column does not already exist in df it will be added.
        If the column does exist all the existing imputed values will be set to NA and re-imputed.
        Default is 'isImputed'.
    max_missing: int
        Maximum number of missing values; above which imputation is skipped.
        If None, no threshold is used. Default is None.
    n_neighbors: int
        n_neighbors passed to sklearn.impute.KNNImputer. Default is 5
    weights: int
        weights passed to sklearn.impute.KNNImputer. Default is 'uniform'
    '''

    imputed_value_col = cammel_case('imputed', value_col)

    # reset existing imputed values
    df = df.copy()
    if is_imputed_col in df.columns:
        df.loc[df[is_imputed_col]][value_col] = pd.NA
    else:
        df[is_imputed_col] = False

    # pivot df wider
    df_w = df.pivot(index=key_cols, columns=replicate_col, values=value_col)

    # drop rows with more than missing_threshold missing values
    if max_missing is not None:
        mask = df_w.apply(lambda x: pd.isna(x).sum(), axis=1) < max_missing
        df_w = df_w[mask]

    # transpose so peptides are columns
    df_w = df_w.T

    # impute missing values
    imputer = skl_KNNImputer(n_neighbors=n_neighbors, weights=weights,
                             keep_empty_features=True)
    df_w = df_w.T

    df_i = pd.DataFrame(imputer.fit_transform(df_w), columns=df_w.columns, index=df_w.index)
    df_i = df_i.T
    df_i = df_i.melt(value_name=imputed_value_col, ignore_index=False).reset_index()

    # Join imputed values to df and clean up
    df = df.set_index(key_cols + [replicate_col])
    df_i = df_i.set_index(key_cols + [replicate_col])
    df = df.join(df_i)
    df[is_imputed_col] = pd.isna(df[value_col]) & ~pd.isna(df[imputed_value_col])
    df[value_col] = df.apply(lambda x: x[imputed_value_col] if x[is_imputed_col] else x[value_col], axis=1)
    df = df.drop(columns=[imputed_value_col]).reset_index()

    return df


class KNNImputer(ImputationManagerBase):
    '''
    KNN imputation manager.

    Inherets from ImputationManagerBase

    Attributes
    ----------
    n_neighbors: int
        Number of nearest neigbors. Passed to KNN imputaion algorithm.
    '''

    def __init__(self, conn, n_neighbors=5, **kwargs):
        '''
        Parameters
        ----------
        conn: sqlite3.Connection
        n_neighbors: int
        **kwargs: dict
            Additional kwargs passed to ImputationManagerBase
        '''
        super().__init__(conn, **kwargs)
        self.n_neighbors = n_neighbors


    def impute(self):
        if self.impute_data in ('precursors', 'both'):
            if not self._read_precursors():
                return False
            replicates = self.precursors['replicateId'].drop_duplicates().to_list()
            self.precursors = knn_impute_df(self.precursors, ['peptideId', 'precursorCharge'], 'area',
                                            max_missing=int(len(replicates) * self.missing_threshold),
                                            n_neighbors=self.n_neighbors)

        if self.impute_data in ('proteins', 'both'):
            if not self._read_proteins():
                return False
            replicates = self.proteins['replicateId'].drop_duplicates().to_list()
            self.proteins = knn_impute_df(self.proteins, ['proteinId'], 'abundance',
                                          max_missing=int(len(replicates) * self.missing_threshold),
                                          n_neighbors=self.n_neighbors)

        return True
