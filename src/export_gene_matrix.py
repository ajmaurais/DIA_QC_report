
import sys
import os
import argparse
import re
from collections import namedtuple, Counter
import sqlite3

import pandas as pd

from .submodules.dia_db_utils import get_meta_value, is_normalized
from .submodules.dia_db_utils import check_schema_version
from .submodules.logger import LOGGER

COMMAND_DESCRIPTION = 'Export PDC gene tables from batch database.'
TABLE_TYPES = ('split', 'drop')

SPLIT_RE = re.compile(r'\s/\s')

GENE_ID_TABLE_COLS = {'protein': 'accession',
                      'accession': 'Gene',
                      'geneid': 'NCBIGeneID',
                      'authid': 'Authority',
                      'description': 'Description',
                      'chromosome': 'Chromosome',
                      'locus': 'Locus'}

PROTEIN_QUERY = '''
SELECT
    prot.accession,
	r.replicate,
    a.aliquot_run_metadata_id,
    s.aliquot_submitter_id,
	q.abundance,
	q.normalizedAbundance
FROM proteinQuants q
LEFT JOIN proteins prot ON prot.proteinId == q.proteinId
LEFT JOIN replicates r ON r.id == q.replicateId
LEFT JOIN (
    SELECT
        replicateId, annotationValue as aliquot_run_metadata_id
    FROM sampleMetadata
    WHERE annotationKey == 'aliquot_run_metadata_id'
) a ON a.replicateId == q.replicateId
LEFT JOIN (
    SELECT
        replicateId, annotationValue as aliquot_submitter_id
    FROM sampleMetadata
    WHERE annotationKey == 'aliquot_submitter_id'
) s ON s.replicateId == q.replicateId
WHERE prot.accession IS NOT NULL AND r.includeRep == TRUE;'''

PRECURSOR_QUERY = '''
SELECT
	r.replicate,
    a.aliquot_run_metadata_id,
    s.aliquot_submitter_id,
	prot.accession as accession,
	p.modifiedSequence,
	p.precursorCharge,
	p.totalAreaFragment as area,
	p.normalizedArea
FROM precursors p
LEFT JOIN peptideToProtein ptp ON ptp.peptideId == p.peptideId
LEFT JOIN proteins prot ON prot.proteinId == ptp.proteinId
LEFT JOIN replicates r ON r.id == p.replicateId
LEFT JOIN (
    SELECT
        replicateId, annotationValue as aliquot_run_metadata_id
    FROM sampleMetadata
    WHERE annotationKey == 'aliquot_run_metadata_id'
) a ON a.replicateId == p.replicateId
LEFT JOIN (
    SELECT
        replicateId, annotationValue as aliquot_submitter_id
    FROM sampleMetadata
    WHERE annotationKey == 'aliquot_submitter_id'
) s ON s.replicateId == p.replicateId
WHERE prot.accession IS NOT NULL AND r.includeRep == TRUE;
'''


def check_aliquot_id(df, df_name):
    # check that aliquot_run_metadata_id column exists in df
    arm_id = 'aliquot_run_metadata_id'
    as_id = 'aliquot_submitter_id'
    for col in (arm_id, as_id):
        if col not in df.columns:
            LOGGER.error('In %s df: %s column does not exist!', df_name, col)
            return False

    # make sure that no aliquot_ids are NA
    rep_ids = df[['replicate', arm_id, as_id]].drop_duplicates()
    all_good = True
    for col in (arm_id, as_id):
        na_counts = Counter(pd.isna(rep_ids[col]))
        if True in na_counts:
            LOGGER.error(
                'Missing %s for %d replicate%s!',
                col, na_counts[True], 's' if na_counts[True] > 1 else ''
            )
            all_good = False
    if not all_good:
        return False

    # Make sure that all aliquot_ids are unique
    ids = {f'{row[as_id]}:{row[arm_id]}' for _, row in rep_ids.iterrows()}
    id_diff = len(set(rep_ids['replicate'].to_list())) - len(ids)
    if id_diff != 0:
        LOGGER.error('Not all aliquot_ids are unique!')
        return False

    return True


def unambigious_match(choices, substr):
    matches = [choice for choice in choices if choice.startswith(substr)]
    if len(matches) == 1:
        return matches[0]
    return None


def check_df_columns(df, required_cols, df_name=None):
    ''' Check that df has all of required_cols. '''
    all_good = True
    df_cols = set(df.columns.to_list())
    for col in required_cols:
        if col not in df_cols:
            LOGGER.error(f'Missing required column: "{col}"' + (f' in {df_name}' if df_name else ''))
            all_good = False
    return all_good


def pivot_data_wider(dat, quant_col, index,
                     group_method='split', use_arm_id_headers=False):

    if group_method not in TABLE_TYPES:
        LOGGER.warning(f"group_method '{group_method}' unknown!")

    rep_id_col = 'replicate'
    if use_arm_id_headers:
        rep_id_col = 'as_id:arm_id'
        dat['as_id:arm_id'] = dat.apply(
            lambda row: f"{row['aliquot_submitter_id']}:{row['aliquot_run_metadata_id']}", axis=1
        )

    ret = dat.pivot(
        columns=rep_id_col, index=index, values=quant_col
    ).reset_index()
    ret.columns.name = None

    ret['accession'] = ret['accession'].apply(lambda x: tuple(SPLIT_RE.split(x)))
    ret['gene_group'] = ret['accession']

    if group_method == 'drop':
        ret = ret[ret['gene_group'].apply(lambda x: len(x) == 1)]

    ret = ret.explode('accession')

    return ret


def concat_gene_data(accession_set, gene_data, sep=' / ', gene_uuid=False):
    '''
    Make a lookup table for gene data for each unique protein accession group.

    Parameters
    ----------
    accession_set: set
        A set where each element in the set is a tuple of protein accessions that were
        observed as part of a gene group in Skyline.
    gene_data: pd.DataFrame
        A data frame with the columns:
        ('accession', 'Gene', 'NCBIGeneID', 'Authority', 'Description', 'Chromosome', 'Locus')
    sep: str
        A string used to separate gene group data. Default is ' / '
    gene_uuid: bool
        Also include gene_uuid column?

    Returns
    -------
    gene_data_lookup: DataFrame
        A dataframe with gene data that can be safely joined to protein or precursor table
        by accession tuple.

    missing_accessions: set
        A set of accession groups that were not found in gene_data
    '''

    gene_data_sets = {}
    gene_data_keys = ['Gene', 'NCBIGeneID', 'Authority', 'Description', 'Chromosome', 'Locus']

    if gene_uuid:
        gene_data_keys.append('gene_uuid')

    Gene = namedtuple('Gene', ' '.join(gene_data_keys))
    gene_data_dict = {}
    for row in gene_data.itertuples():
        gene_data = [row.Gene, row.NCBIGeneID, row.Authority,
                     row.Description, row.Chromosome, row.Locus]

        if gene_uuid:
            gene_data.append(row.gene_uuid)

        gene_data_dict[row.accession] = Gene(*gene_data)

    def iter_accession_group(group):
        if isinstance(group, (list, tuple)):
            for i in group:
                yield i
        elif isinstance(group, str):
            yield group
        else:
            raise ValueError('Unknown group type!')

    missing_accessions = set()
    for accession_group in accession_set:
        row_initiated = False
        for accession in iter_accession_group(accession_group):
            if accession not in gene_data_dict:
                missing_accessions.add(accession)
                continue

            # init empty accession group in gene_data_sets
            if not row_initiated:
                gene_data_sets[accession_group] = set()
                row_initiated = True

            gene_data_sets[accession_group].add(gene_data_dict[accession])

    # convert gene_data_sets into a DataFrame
    ret = {key: [] for key in ['accession'] + list(gene_data_keys)}
    for accession_group, gene_data_group in gene_data_sets.items():
        ret['accession'].append(accession_group)
        for gene_data_key in gene_data_keys:
            ret[gene_data_key].append(sep.join([str(getattr(x, gene_data_key)) for x in gene_data_group]))
    ret = pd.DataFrame(ret)

    return ret, missing_accessions


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)

    gene_group_args = parser.add_argument_group(
        'Gene grouping',
        description="Choose how to display gene groups in output tables. "
                    # "'combined': Show all genes in group in the same row separated by --sep. "
                    "'split': Show each gene in group on a separate row. "
                    "'skip': Don't include gene groups with more than one gene in table."
    )
    gene_group_args.add_argument(
        '-g', '--groupMethod', default='split', type=str, dest='group_method',
        help=f''' How to display gene groups in tables. Grouping method should unambiguously match
                  one of ['{"', '".join(TABLE_TYPES)}']'''
    )
    gene_group_args.add_argument(
        '-s', '--sep', default='; ', type=str, dest='gene_group_sep',
        help="Gene group separator. Default is '; '"
    )
    parser.add_argument('--prefix', default=None, help='Prefix to add to output file names.')
    parser.add_argument(
        '--useAliquotId', default=False, action='store_true', dest='use_aliquot_id',
        help='Use aliquot_run_metadata_id as column headers instead of replicate name.'
    )
    parser.add_argument('--addGeneUuid', default=False, action='store_true',
                        dest='add_gene_uuid', help='Add column for gene id hash.')
    parser.add_argument('gene_table', help='A tsv with gene data.')
    parser.add_argument('database', help='Path to sqlite batch/qc database.')

    return parser.parse_args(argv)


def _main(argv, prog=None):
    ''' Actual main method.  '''
    args = parse_args(argv, prog=prog)

    if (gene_group_method := unambigious_match(TABLE_TYPES, args.group_method)) is None:
        LOGGER.error(f"Could not unambiguously determine protein table type: '{args.group_method}'\n")
        sys.exit(1)

    if os.path.isfile(args.database):
        conn = sqlite3.connect(args.database)
    else:
        LOGGER.error('Database file does not exist!')
        sys.exit(1)

    # check database version
    if not check_schema_version(conn):
        sys.exit(1)

    # check that existing precursors in db were grouped by current method
    if get_meta_value(conn, 'group_precursors_by') != 'gene':
        LOGGER.error('Precursors in database must be grouped by gene!')
        sys.exit(1)

    if not (db_normalized := is_normalized(conn)):
        LOGGER.warning('metadata.is_normalized is False. Only using unnormalized values.')

    # read proteins
    LOGGER.info('Reading protein table from database...')
    dat_protein = pd.read_sql(PROTEIN_QUERY, conn)
    LOGGER.info('Done reading protein table!')

    # read precursors
    LOGGER.info('Reading precursor table from database...')
    dat_precursor = pd.read_sql(PRECURSOR_QUERY, conn)
    LOGGER.info('Done reading precursor table!')

    conn.close()

    # check aliquot_run_metadata_id column if necissary
    rep_id_col = 'replicate'
    if args.use_aliquot_id:
        if not check_aliquot_id(dat_protein, 'protein') or not check_aliquot_id(dat_precursor, 'precursor'):
            sys.exit(1)

    gene_id_table_cols = GENE_ID_TABLE_COLS.copy()
    if args.add_gene_uuid:
        gene_id_table_cols['gene_uuid'] = 'gene_uuid'

    # read gene ID lookup table
    gene_ids = pd.read_csv(args.gene_table,
                           sep=('\t' if os.path.splitext(args.gene_table)[1] == '.tsv' else ','))
    if not check_df_columns(gene_ids, gene_id_table_cols):
        sys.exit(1)
    gene_ids = gene_ids[list(gene_id_table_cols.keys())].rename(columns=gene_id_table_cols)

    proteins = {'proteins_unnormalized': 'abundance'}
    precursors = {'precursors_unnormalized': 'area'}
    if db_normalized:
        proteins['proteins_normalized'] = 'normalizedAbundance'
        precursors['precursors_normalized'] = 'normalizedArea'
    else:
        LOGGER.warning('Database is not normalized. Only writing unnormalized reports.')

    # pivot proteins wider
    accessions = set()
    for method in proteins:
        proteins[method] = pivot_data_wider(
            dat_protein, proteins[method], 'accession',
            group_method=gene_group_method, use_arm_id_headers=args.use_aliquot_id
        )
        accessions = accessions | set(proteins[method]['accession'].drop_duplicates().to_list())

    # pivot precursors wider
    for method in precursors:
        precursors[method] = pivot_data_wider(
            dat_precursor, precursors[method], ['accession', 'modifiedSequence', 'precursorCharge'],
            group_method=gene_group_method, use_arm_id_headers=args.use_aliquot_id
        )
        accessions = accessions | set(precursors[method]['accession'].drop_duplicates().to_list())

    data = proteins | precursors

    # make accession gene data lookup
    gene_ids = gene_ids.drop_duplicates(subset='accession')
    LOGGER.info('Looking up gene data for protein accessions.')
    gene_data_lookup, missing_accessions = concat_gene_data(
        accessions, gene_ids, sep=args.gene_group_sep, gene_uuid=args.add_gene_uuid
    )
    LOGGER.info(f'Found gene data for {len(gene_data_lookup)} protein accessions.')
    if len(missing_accessions) > 0:
        LOGGER.warning(f'Missing gene data for {len(missing_accessions)} protein accessions!')

    # make gene group lookup
    gene_id_lookup = {row.accession: row.Gene for row in gene_data_lookup.itertuples()}
    gene_groups = {}
    for group in data['proteins_unnormalized']['gene_group'].drop_duplicates().to_list():
        gene_id_group = set()
        for prot_id in group:
            if prot_id in gene_id_lookup:
                gene_id_group.add(gene_id_lookup[prot_id])
        gene_groups[group] = args.gene_group_sep.join(gene_id_group)

    # join gene data to tables
    for method in data:
        data[method] = data[method].set_index('accession').join(
            gene_data_lookup.set_index('accession'),
            how='inner', on='accession'
        )
        data[method] = data[method].reset_index(drop=True)

        data[method]['GeneGroup'] = data[method]['gene_group'].apply(lambda x: gene_groups[x])

        columns = list(data[method].columns)
        columns.pop(columns.index('gene_group'))
        columns.insert(0, columns.pop(columns.index('Gene')))
        columns.insert(1, columns.pop(columns.index('GeneGroup')))
        data[method] = data[method].loc[:, columns]

    # write tables
    for method in data:
        fname = f'{"" if args.prefix is None else args.prefix + "_"}{method}.tsv'
        LOGGER.info(f'Writing {fname}...')
        data[method].to_csv(fname, sep='\t', index=False)
        LOGGER.info(f'Done writing {fname}')


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc export_gene_matrix" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()