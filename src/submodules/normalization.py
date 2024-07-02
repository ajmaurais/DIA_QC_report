
from multiprocessing import cpu_count
from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
import directlfq.utils as lfqutils
from directlfq.normalization import NormalizationManagerSamplesOnSelectedProteins as dlfq_norm
import directlfq.protein_intensity_estimation as lfqprot_estimation

from .logger import LOGGER


class NormalizationManagerBase(ABC):
    '''
    NormalizationManager abastract base class.

    Child classes must implement `normalize` method.

    Attributes
    ----------
    conn: sqlite.Connection
        A connection to a precursor database
    precursors: pd.DataFrame
        Long formatted dataframe of precursor quantities.
    proteins: pd.DataFrame
        Long formatted dataframe of protein quantities.
    keep_na: bool
        Should values missing in 1 or more batches be kept? Default is False.
    '''

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
            df = df.set_index(['proteinId', 'modifiedSequence', 'precursorCharge'])
            n_reps = len(df['replicateId'].drop_duplicates().index)
            na_counts = df.groupby(level=[0, 1, 2])['area'].apply(lambda x: len(x[~pd.isna(x)]))
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
        Get long formatted precursor and protein tables.

        Parameters
        ----------
        use_db_ids: bool
            Should dataframe have database protein and replicate IDs?
        '''

        if use_db_ids:
            return self.precursors, self.proteins

        proteins = self.proteins.copy()
        precursors = self.precursors.copy()

        # add replicate column
        cur = self.conn.cursor()
        cur.execute('SELECT id, replicate FROM replicates;')
        rep_ids = {int(x[0]): x[1] for x in cur.fetchall()}
        precursors['replicate'] = precursors['replicateId'].apply(lambda x: rep_ids[x])
        proteins['replicate'] = proteins['replicateId'].apply(lambda x: rep_ids[x])

        # add protein name column
        cur.execute('SELECT proteinId, name FROM proteins;')
        prot_ids = {int(x[0]): x[1] for x in cur.fetchall()}
        proteins['protein'] = proteins['proteinId'].apply(lambda x: prot_ids[x])

        return precursors, proteins


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
        '''
        Populate self.proteins and self.precursors with median normalized values.
        '''
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

        return True


class DirectlfqNormalizer(NormalizationManagerBase):
    def __init__(self, conn):
        NormalizationManagerBase.__init__(self, conn, keep_na=False)


    def normalize(self):
        '''
        Populate self.proteins with DirectLFQ normalized values and
        self.precursors with median normalized values.
        '''
        if not self._read_precursors():
            return False

        if len(self.precursors.index) == 0:
            LOGGER.error('Can not perform DirectLFQ normalization with 0 precursors!')
            return False

        input_df = self.precursors.rename(columns={'proteinId': 'protein'})
        input_df['ion'] = input_df['modifiedSequence'] + '_' + input_df['precursorCharge'].astype(str)
        input_df = input_df.pivot(index=['protein', 'ion'],
                                  columns='replicateId', values='area')

        # format input_df for DirectLFQ normalization
        input_df.columns.name = None
        input_df = input_df.reset_index()
        input_df = lfqutils.index_and_log_transform_input_df(input_df)
        input_df = lfqutils.remove_allnan_rows_input_df(input_df)
        input_df = dlfq_norm(input_df, num_samples_quadratic=10).complete_dataframe

        # Do DirectLFQ normalization
        LOGGER.info('Running DirectLFQ normalization.')
        protein_df, _ = lfqprot_estimation.estimate_protein_intensities(input_df,
                                                                        min_nonan=1,
                                                                        num_samples_quadratic=10,
                                                                        num_cores=cpu_count())
        LOGGER.info('Finished with DirectLFQ normalization.')

        # pivot protein_df longer
        self.proteins = protein_df.melt(id_vars='protein',
                                        value_name='normalizedAbundance',
                                        var_name='replicateId')
        self.proteins = self.proteins.rename(columns={'protein': 'proteinId'})

        # median normalize precursors
        self.precursors = self.precursors[['replicateId',
                                           'peptideId',
                                           'modifiedSequence',
                                           'precursorCharge',
                                           'area']].drop_duplicates()
        self.precursors = median_normalize_df(self.precursors, ['replicateId'], 'area')

        return True
