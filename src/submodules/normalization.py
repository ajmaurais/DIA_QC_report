
from multiprocessing import cpu_count
from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
import directlfq.utils as lfqutils
from directlfq.normalization import NormalizationManagerSamplesOnSelectedProteins as dlfq_norm
import directlfq.protein_intensity_estimation as lfqprot_estimation

from .logger import LOGGER


class NormalizationManagerBase(ABC):
    def __init__(self, conn, keep_na=False):
        self.conn = conn
        self.precursors = None
        self.proteins = None
        self.keep_na = keep_na


    def _read_precursors(self):
        # get precursor table from db
        LOGGER.info('Reading precursors from database...')
        df = pd.read_sql('''
            SELECT
                p.replicateId,
                p.peptideId,
                ptp.proteinId,
                p.modifiedSequence,
                p.precursorCharge,
                p.totalAreaFragment as area
            FROM precursors p
            LEFT JOIN peptideToProtein ptp ON ptp.peptideId == p.peptideId
            LEFT JOIN replicates r ON p.replicateId == r.id
            WHERE r.includeRep == TRUE; ''',
            self.conn)
        LOGGER.info('Finished reading precursors.')

        if len(df.index) == 0:
            LOGGER.error('All replicates in database have been excluded!')
            return False

        if not self.keep_na:
            df = df.set_index(['modifiedSequence', 'precursorCharge'])
            n_reps = len(df['replicateId'].drop_duplicates().index)
            na_counts = df.groupby(['modifiedSequence', 'precursorCharge'])['area'].apply(lambda x: len(x[~pd.isna(x)]))
            na_counts = n_reps - na_counts
            n_precursors = len(na_counts.index)
            na_counts = na_counts[na_counts == 0]

            n_not_missing = len(na_counts.index)
            n_missing = n_precursors - n_not_missing
            if n_missing > 0:
                LOGGER.warning(f'Removed {n_missing:,} of {n_precursors:,} precursors with missing values.')
                LOGGER.warning(f'{n_not_missing:,} precursors without missing values remain.')

                df = df[df.index.isin(na_counts.index)]
            df = df.reset_index()

        self.precursors = df

        return True


    @abstractmethod
    def normalize(self):
        pass


    def get_long_tables(self, use_db_ids=False):
        '''
        Get long formated precursor and protein tables.
        '''
        return self.precursors, self.proteins


    def get_wide_tables(self, normalized=True, use_db_ids=False):
        pass


def cammel_case(prefix, suffix):
    return f'{prefix}{suffix[0].upper()}{suffix[1:]}'


def median_normalize_df(df, key_cols, value_col):
    '''
    Median normalize peak areas in long formatted datafraeme.

    Parameters
    ----------
    df: pd.DataFrame
        Long formatted data frame
    key_cols: list
        The names of column(s) which uniquely identify each replicate.
    value_col: str
        The name of the column with peak areas.
        This function will log2 transform this column so the areas should be in linear space.
    '''

    log2_value_col = cammel_case('log2', value_col)
    norm_value_col = cammel_case('normalized', value_col)
    log2_norm_value_col = cammel_case('log2', norm_value_col)

    cols = df.columns.tolist()
    cols.append(norm_value_col)

    # log2 transform
    df[log2_value_col] = np.log2(df[value_col] + 1)

    # determine value to shift log2Area for each replicate
    medians = df.groupby(key_cols)[log2_value_col].median()
    median_shifts = medians - np.mean(medians)
    median_shifts = median_shifts.reset_index()
    median_shifts = {row.replicateId: getattr(row, log2_value_col) for row in median_shifts.itertuples()}

    # shift log2Areas to normalize medians
    df[log2_norm_value_col] = df[log2_value_col] - df['replicateId'].apply(lambda x: median_shifts[x])

    # convert back to linear space
    df[norm_value_col] = np.power(2, df[log2_norm_value_col]) - 1

    return df


class MedianNormalizer(NormalizationManagerBase):

    def normalize(self):
        if not self._read_precursors():
            return False

        # median normalize proteins
        self.proteins = self.precursors.groupby(['replicateId', 'proteinId'])['area'].sum()
        self.proteins.name = 'abundance'
        self.proteins = self.proteins.reset_index()
        self.proteins = median_normalize_df(self.proteins, ['replicateId'], 'abundance')

        # median normalize precursors
        self.precursors = self.precursors[['replicateId',
                                           'peptideId',
                                           'modifiedSequence',
                                           'precursorCharge',
                                           'area']].drop_duplicates()
        self.precursors = median_normalize_df(self.precursors, ['replicateId'], 'area')

# def directlfq(conn):
#
#     df['ion'] = df['modifiedSequence'] + '_' + df['precursorCharge'].astype(str)
#
#     # add 1 to all precursor areas so the log2 of zeros is finite
#     df['totalAreaFragment'] = df['totalAreaFragment'] + 1
#     df = df.reset_index()
#
#     # pivot wider
#     input_df = df.pivot(index=['protein', 'ion'], columns='replicateId', values='totalAreaFragment')
#
#     # remove missing rows
#     n_not_missing = len(input_df.index)
#     input_df = input_df[~input_df.apply(lambda x: any(x.isnull().values), axis = 1)]
#     keep_precursors = set([x[1] for x in input_df.index.tolist()])
#     if len(input_df.index) < n_not_missing:
#         LOGGER.warning(f'Removed {n_not_missing - len(input_df.index):,} of {n_not_missing:,} precursors with missing values.')
#         LOGGER.warning(f'{len(input_df.index):,} precursors without missing values remain.')
#
#     # format input_df for DirectLFQ normalization
#     input_df.columns.name = None
#     input_df = input_df.reset_index()
#     input_df = lfqutils.index_and_log_transform_input_df(input_df)
#     input_df = lfqutils.remove_allnan_rows_input_df(input_df)
#     input_df = dlfq_norm(input_df, num_samples_quadratic=10).complete_dataframe
#
#     # Do DirectLFQ normalization
#     LOGGER.info('Running DirectLFQ normalization.')
#     protein_df, _ = lfqprot_estimation.estimate_protein_intensities(input_df,
#                                                                     min_nonan=1,
#                                                                     num_samples_quadratic=10,
#                                                                     num_cores=cpu_count())
#     LOGGER.info('Finished with DirectLFQ normalization.')
#
#     # pivot protein_df longer
#     protein_df = protein_df.melt(id_vars='protein', value_name='normalizedArea', var_name='replicateId')
#
#
#     # remove precursors with missing values
#     df = df[df['ion'].apply(lambda x: x in keep_precursors)]
#
#     # median normalize precursor quants
#     ion_df = median_normalize(df[['replicateId', 'peptideId', 'precursorCharge', 'totalAreaFragment']].drop_duplicates(),
#                               ['replicateId'], 'totalAreaFragment')
#
#     return ion_df, protein_df
#
