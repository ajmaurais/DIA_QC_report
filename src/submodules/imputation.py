
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


def knn_impute_df(df, key_cols, value_col,
                  replicate_col='replicateId', is_imputed_col='isImputed',
                  max_missing=None, skip_all_missing=True,
                  n_neighbors=5, weights='uniform'):
    '''
    KNN impute peak areas in a long formated dataframe.

    Parameters
    ----------
    df: pd.DataFrame
        Long formated dataframe
    key_cols: list, str
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
    skip_all_missing: bool
        If True, don't impute values for groups that are all missing.
        If False, values for compleetly missing groups are set to 0.
        Default is True.
    n_neighbors: int
        n_neighbors passed to sklearn.impute.KNNImputer. Default is 5
    weights: int
        weights passed to sklearn.impute.KNNImputer. Default is 'uniform'

    Returns
    -------
    df_imputed: pd.DataFrame
        A dataframe the same shape as `df` where the `value_col` and `is_imputed_col` have been updated.
    '''

    imputed_value_col = cammel_case('imputed', value_col)
    _key_cols = [key_cols] if isinstance(key_cols, str) else key_cols

    # reset existing imputed values
    df = df.copy()
    if is_imputed_col in df.columns:
        df.loc[df[is_imputed_col]][value_col] = pd.NA
    else:
        df[is_imputed_col] = False

    # pivot df wider
    df_w = df.pivot(index=_key_cols, columns=replicate_col, values=value_col)

    # drop rows with all missing values
    if skip_all_missing:
        mask = df_w.apply(lambda x: all(pd.isna(x)), axis=1)
        df_w = df_w[~mask]

    # drop rows with more than max_missing values
    if max_missing is not None:
        mask = df_w.apply(lambda x: pd.isna(x).sum(), axis=1) < max_missing
        df_w = df_w[mask]

    # Return if there are no rows remaining after applying max_missing and skip_all_missing filters.
    if len(df_w.index) == 0:
        LOGGER.warning('Not enough values for imputation!')
        return df

    # transpose so peptides are columns
    df_w = df_w.T

    # impute missing values
    imputer = skl_KNNImputer(n_neighbors=n_neighbors, weights=weights,
                             keep_empty_features=True)
    df_i = pd.DataFrame(imputer.fit_transform(df_w), columns=df_w.columns, index=df_w.index)
    df_i = df_i.T
    df_i = df_i.melt(value_name=imputed_value_col, ignore_index=False).reset_index()

    # Join imputed values to df and clean up
    df = df.set_index(_key_cols + [replicate_col])
    df_i = df_i.set_index(_key_cols + [replicate_col])
    df = df.join(df_i)
    df[is_imputed_col] = pd.isna(df[value_col]) & ~pd.isna(df[imputed_value_col])
    df[value_col] = df.apply(lambda x: x[imputed_value_col] if x[is_imputed_col] else x[value_col], axis=1)
    df = df.drop(columns=[imputed_value_col]).reset_index()

    return df


def split_by_project(df, project_col='project', drop=False):
    '''
    Split dataframe into a dict where the keys are the project name
    and the values are a subset dataframe for the project.

    Parameters
    ----------
    df: pd.DataFrame
        The dataframe to split.
    project_col: str
        The name of the column to split on. Default is 'project'
    drop: bool
        Should the `project_col` be dropped? Default is False.

    Returns
    -------
    df_dict: Dictionary
        A dictionary of key, pd.DataFrame pairs.
    '''
    projects = dict()
    for project in df[project_col].drop_duplicates().to_list():
        projects[project] = df[df[project_col] == project]
        if drop:
            projects[project] = projects[project].drop(columns=[project_col])
    return projects


def bind_rows(df_dict, key_name=None):
    '''
    Take a dictionary of pd.DataFrame(s) and return a single dataframe where all the dataframes
    in have been stacked on top of each other.

    Parameters
    ----------
    df_dict: Dictionary
        A dictionary of key, pd.DataFrame pairs.
    key_name: str
        The name of the column to add for the dictionary key.
        If None, a column is not added for the key. Default is None.

    Returns
    -------
    df: pd.DataFrame
        A single combined dataframe.
    '''
    merged_dfs = []
    for key, df in df_dict.items():
        if key_name:
            df = df.copy()
            df[key_name] = key
        merged_dfs.append(df)

    # Combine all DataFrames
    return pd.concat(merged_dfs, ignore_index=True)


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


    def _impute(self, df, key_cols, value_col, group_by_project):
        replicates = df['replicateId'].drop_duplicates().to_list()

        # split df by project if applicable
        if group_by_project:
            projects = split_by_project(df)
            for project in projects:
                max_missing = int(len(replicates) * self.missing_threshold)
                projects[project] = knn_impute_df(projects[project],
                                                  key_cols, value_col,
                                                  max_missing=max_missing,
                                                  n_neighbors=self.n_neighbors)
            return bind_rows(projects)

        return knn_impute_df(df, key_cols, value_col,
                             max_missing=int(len(replicates) * self.missing_threshold),
                             n_neighbors=self.n_neighbors)


    def impute(self):

        cur = self.conn.cursor()
        cur.execute('SELECT DISTINCT project FROM replicates;')
        n_projects = len(cur.fetchall())
        group_by_project = False if n_projects == 1 else self.group_by_project

        if self.impute_data in ('precursors', 'both'):
            if not self._read_precursors():
                return False
            self.precursors = self._impute(self.precursors, ['peptideId', 'precursorCharge'], 'area',
                                           group_by_project)

        if self.impute_data in ('proteins', 'both'):
            if not self._read_proteins():
                return False
            self.proteins = self._impute(self.proteins, 'proteinId', 'abundance',
                                         group_by_project)

        return True
