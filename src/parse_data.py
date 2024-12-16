
import sys
import argparse
import sqlite3
import os
from datetime import datetime
import re

import pandas as pd

from .submodules.logger import LOGGER
from .submodules.dia_db_utils import SCHEMA, SCHEMA_VERSION
from .submodules.skyline_reports import ReplicateReport, PrecursorReport
from .submodules.dia_db_utils import PRECURSOR_KEY_COLS, METADATA_TIME_FORMAT
from .submodules.dia_db_utils import update_metadata_dtypes, update_acquired_ranks
from .submodules.dia_db_utils import update_meta_value
from .submodules.dia_db_utils import get_meta_value
from .submodules.dia_db_utils import check_schema_version
from .submodules.dia_db_utils import update_command_log
from .submodules.read_metadata import Metadata
from . import __version__ as PROGRAM_VERSION

COMMAND_DESCRIPTION = 'Generate QC and batch correction database from Skyline reports.'
DUPLICATE_PRECURSOR_CHOICES = ('e', 'm', 'f', 'i')


def _initialize(fname):
    ''' Initialize empty database with SCHEMA '''
    conn = sqlite3.connect(fname)
    cur = conn.cursor()
    for command in SCHEMA:
        try:
            cur.execute(command)
        except sqlite3.OperationalError as e:
            raise RuntimeError(f'Error running command: {command}') from e

    conn.commit()
    return conn


def _add_index_column(col, index, df_name=None, index_name=None):
    '''
    Add a column to dataframe from a dict index.
    Also handle KeyErrors when values in key_col do not exist as a key in index.

    Parameters
    ----------
    col: pd.Series
        A Series with key values mapping to index.
    index: dict
        A dict mapping key_col to values in index.
    df_name: str
        The name of the df to use in any error messages.
    index_name: str
        The name of the index to use in any error messages.

    Returns
    -------
    pd.Series: The new column to add to df, or None if there was a missing value in index.
    '''

    ret = col.copy()

    # Check that all values in col exist in index
    table_names = '' if df_name is None or index_name is None else f' for {df_name} in {index_name}'
    all_good = True
    for key in ret.drop_duplicates().tolist():
        if key not in index:
            LOGGER.error(f'Missing required value: {key}{table_names}!')
            all_good = False

    # Exit now if there are one or more missing values
    if not all_good:
        return None

    return ret.apply(lambda x: index[x])


def write_db(fname, replicates, precursors, protein_quants=None,
             sample_metadata=None, sample_metadata_types=None,
             project_name=None, overwriteMode='error',
             group_precursors_by='protein'):
    '''
    Write reports to precursor sqlite database.

    Parameters
    ----------
    fname: str
        The name of the database file.
    replicates: pd.DataFrame
        Replicates dataframe
    precursors: pd.DataFrame
        Precursors dataframe
    protein_quants: pd.DataFrame
        Long formatted protein quants dataframe (optional)
    sample_metadata: pd.DataFrame
        Metadata values dataframe (optional)
    sample_metadata_types: dict
        Dict mapping metadata annotationKey to annotationType (optional)
    project_name: str
        The project name to use in the replicates table.
    overwriteMode: str
        Behavior when the output file already exists.
        Choices are 'error', 'overwrite', or 'append'
    group_precursors_by: str
        How are precursors grouped?
        The only effect this argument has is to add an entry in the
        database metadata table.

    Returns
    -------
    success: bool
        True if all operations were successful, false otherwise.
    '''

    # Metadata to add to db
    log_metadata = dict()
    log_metadata['group_precursors_by'] = group_precursors_by
    log_metadata['schema_version'] = SCHEMA_VERSION
    log_metadata['dia_qc version'] = PROGRAM_VERSION

    if sample_metadata is not None and sample_metadata_types is None:
        LOGGER.error('Must specify both sample_metadata and sample_metadata_types!')
        return False

    # Initialize database if applicable and get database connection
    conn = None
    append = False
    if os.path.isfile(fname):
        if overwriteMode == 'error':
            LOGGER.error(f'{fname} already exists. Use the --overwriteMode option to append or overwrite')
            return False

        if overwriteMode == 'overwrite':
            LOGGER.warning(f'{fname} already exists. Overwriting...')
            os.remove(fname)
            conn = _initialize(fname)

        elif overwriteMode == 'append':
            conn = sqlite3.connect(fname)
            append = True

            # check database version
            if not check_schema_version(conn):
                return False

            # check that existing precursors in db were grouped by current method
            if (current_group_by := get_meta_value(conn, 'group_precursors_by')) is None:
                return False
            if current_group_by != group_precursors_by:
                LOGGER.error(f"Database grouping method '{current_group_by}' does not match '{group_precursors_by}'")
                return False

        else:
            LOGGER.error(f'"{overwriteMode}" is an unknown overwriteMode!')
            return False

    else:
        conn = _initialize(fname)

    # create proteins table
    if protein_quants is not None:
        proteins = protein_quants[['name', 'accession', 'description']].drop_duplicates().reset_index(drop=True)
        proteinIndex = {r: i for i, r in zip(proteins['name'].index, proteins['name'])}

        # check that protein accessions in protein and precursor table match
        if set(proteins['name'].to_list()) != set(precursors['proteinName'].to_list()):
            LOGGER.error('Protein and precursor ProteinAccessions differ!')
            return False
    else:
        proteins = precursors[['proteinName', 'proteinAccession']].drop_duplicates().reset_index(drop=True)
        proteins = proteins.rename(columns={'proteinName': 'name', 'proteinAccession': 'accession'})
        proteins['description'] = None

    # if project_name is not specified, set it to project_(n + 1)
    if project_name is None:
        projects = pd.read_sql('SELECT DISTINCT project FROM replicates;', conn)
        project_name = f'project_{len(projects.index) + 1}'
    replicates['project'] = project_name

    # populate some metadata values now that we have the project_name
    log_metadata['replicates.acquiredRank updated'] = False
    log_metadata[f'is_normalized'] = False # set this to False because we are adding unnormalized data

    # deal with existing replicates, proteins, and protein to precursor pairs
    if append:
        conn = sqlite3.connect(fname)

        cur = conn.cursor()
        cur.execute('SELECT DISTINCT project FROM replicates WHERE project = ?', (project_name,))
        existing_project = cur.fetchall()
        if len(existing_project) > 0:
            LOGGER.error(f'{project_name} already exists in db!')
            conn.close()
            return False

        # deal with existing replicates
        curReplicates = pd.read_sql('SELECT * FROM replicates;', conn)
        replicates.index = pd.RangeIndex(start=len(curReplicates.index),
                                         stop=len(curReplicates.index) + len(replicates.index),
                                         step=1)
        repIndex = {r: i for i, r in zip(replicates['replicate'].index, replicates['replicate'])}

        # deal with existing proteins
        curProteins = pd.read_sql('SELECT proteinId, name FROM proteins;', conn)
        curProteinNames = set(curProteins['name'].to_list())
        proteins = proteins[proteins['name'].apply(lambda x: x not in curProteinNames)]
        proteins.index = pd.RangeIndex(start=len(curProteins.index),
                                       stop=len(curProteins.index) + len(proteins.index),
                                       step=1)
        proteinIndex = {r: i for i, r in zip(curProteins['proteinId'], curProteins['name'])}
        proteinIndex = {**proteinIndex, **{r: i for i, r in zip(proteins['name'].index, proteins['name'])}}

        # Get initial peptideId
        cur = conn.cursor()
        cur.execute('SELECT MAX(peptideId) FROM peptideToProtein;')
        initial_peptide_id = cur.fetchall()
        if len(initial_peptide_id) != 1:
            LOGGER.error('Could not determine initial peptideId!')
            return False
        initial_peptide_id = initial_peptide_id[0][0] + 1

    else:
        # make dict for replicateId column which links the replicate and precursor tables
        repIndex = {r: i for i, r in zip(replicates['replicate'].index, replicates['replicate'])}

        # make dict for proteinId column which links the protein and proteinQuants and precursor tables
        proteinIndex = {r: i for i, r in zip(proteins['name'].index, proteins['name'])}

        initial_peptide_id = 0

    if protein_quants is None:
        protein_quants = precursors.groupby(['replicateName', 'proteinName'])['totalAreaFragment'].sum().reset_index()
        protein_quants = protein_quants.rename(columns={'proteinName': 'name',
                                                        'totalAreaFragment': 'abundance'})

    peptide_id_index = {s: i + initial_peptide_id for i, s in enumerate(precursors['modifiedSequence'].drop_duplicates())}
    peptide_to_protein = precursors[['proteinName', 'modifiedSequence']].drop_duplicates()

    # add replicateId column to protein_quants
    if (rep_id_col := _add_index_column(protein_quants['replicateName'], repIndex,
                                        'protein_quants', 'repIndex')) is None:
        return False
    protein_quants['replicateId'] = rep_id_col

    # Add proteinId column to protein_quants
    if (prot_id_col := _add_index_column(protein_quants['name'], proteinIndex,
                                         'protein_quants', 'proteinIndex')) is None:
        return False
    protein_quants['proteinId'] = prot_id_col

    protein_quants = protein_quants[['replicateId', 'proteinId', 'abundance']]

    # add replicateId column to precursors
    if (rep_id_col := _add_index_column(precursors['replicateName'], repIndex,
                                        'precursors', 'repIndex')) is None:
        return False
    precursors['replicateId'] = rep_id_col

    # add peptideId column to precursors
    if (pep_id_col := _add_index_column(precursors['modifiedSequence'], peptide_id_index,
                                        'precursors', 'peptideId')) is None:
        return False
    precursors['peptideId'] = pep_id_col

    # Drop duplicate precursors here
    precursor_quality_columns = list(PRECURSOR_KEY_COLS)
    precursor_quality_columns.append('modifiedSequence')
    precursor_quality_columns += [col.name for col in PrecursorReport().numeric_columns()]

    precursors = precursors[precursor_quality_columns].drop_duplicates()

    # add proteinId column to peptide_to_protein
    if (prot_id_col := _add_index_column(peptide_to_protein['proteinName'], proteinIndex,
                                        'peptideToProtein', 'proteinId')) is None:
        return False
    peptide_to_protein['proteinId'] = prot_id_col

    # add peptideId column to peptide_to_protein
    if (pep_id_col := _add_index_column(peptide_to_protein['modifiedSequence'], peptide_id_index,
                                        'peptideToProtein', 'peptideId')) is None:
        return False
    peptide_to_protein['peptideId'] = pep_id_col
    peptide_to_protein = peptide_to_protein[['proteinId', 'peptideId']]

    replicates.to_sql('replicates', conn, if_exists='append', index=True, index_label='id')
    precursors.to_sql('precursors', conn, index=False, if_exists='append')
    proteins.to_sql('proteins', conn, if_exists='append', index=True, index_label='proteinId')
    protein_quants.to_sql('proteinQuants', conn, index=False, if_exists='append')
    peptide_to_protein.to_sql('peptideToProtein', conn, if_exists='append', index=False)

    if sample_metadata is not None:
        missing_metadata = [x for x in sample_metadata['Replicate'].drop_duplicates().to_list() if x not in repIndex]
        for rep in missing_metadata:
            LOGGER.warning(f'Metadata row: \"{rep}\" does not exist in replicate report!')

        sample_metadata = sample_metadata.loc[sample_metadata['Replicate'].apply(lambda x: x in repIndex),]

        if (rep_id_col := _add_index_column(sample_metadata['Replicate'], repIndex,
                                            'sample_metadata', 'repIndex')) is None:
            return False

        if append:
            conn = update_metadata_dtypes(conn, sample_metadata_types)
        else:
            insert_query = 'INSERT INTO sampleMetadataTypes (annotationKey, annotationType) VALUES (?, ?)'
            cur = conn.cursor()
            cur.executemany(insert_query, [(k, str(v)) for k, v in sample_metadata_types.items()])
            conn.commit()

        sample_metadata['replicateId'] = rep_id_col
        sample_metadata = sample_metadata[['replicateId', 'annotationKey', 'annotationValue']]
        sample_metadata.to_sql('sampleMetadata', conn, index=False, if_exists='append')

    # update metadata tablr
    for key, value in log_metadata.items():
        conn = update_meta_value(conn, key, value)
    conn = update_acquired_ranks(conn)

    # update commandLog
    update_command_log(conn, sys.argv, os.getcwd())

    conn.close()
    return True


def check_duplicate_precursors(precursors, mode):
    '''
    Check that all precursors are unique.
    The same precursor could be assigned to multiple proteins in Skyline.
    When you manually adjust the peak integration boundaries for a precursor, Skyline
    will only adjust the integration boundaries for the precursor you are viewing.
    If the precursor is assigned to more than one protein the boundaries will not be synced.
    It is not possible to store duplicate precursors in the database which have different peak
    areas. We must ask the user to specify how to deal with duplicate precursors if they exist.

    Parameters
    ----------
    precursors: pd.DataFrame
        Precursor dataframe.
    mode: str
        One of ['e', 'm', 'f', 'i'] indicating how to handle duplicates.

    Returns
    -------
    precursors: pd.DataFrame
        precursors df with duplicates removed.

    '''

    if mode not in DUPLICATE_PRECURSOR_CHOICES:
        raise ValueError(f'"{mode}" is an invalid mode!')

    LOGGER.info('Checking that all precursors are unique...')

    precursor_cols = list(PrecursorReport().columns())

    key_cols = [col.name for col in precursor_cols
                if col.skyline_name in ('ReplicateName', 'ModifiedSequence', 'PrecursorCharge')]
    n_unique = len(precursors[key_cols].drop_duplicates().index)

    unique_precursors = precursors.drop_duplicates(subset=key_cols + ['totalAreaFragment'])
    n_unique_areas = len(unique_precursors.index)

    all_cols = [col.name for col in precursor_cols
                if col.skyline_name not in ['ProteinAccession', 'ProteinName', 'ProteinGene']]
    all_cols = [col for col in all_cols if col in precursors.columns]
    n_unique_rows = len(precursors[all_cols].drop_duplicates().index)

    # This should never happen.
    # If this is true it means that there are duplicate precursos with the same sequence, charge, and area,
    # but some other floating point column has a non-unique value.
    if n_unique_areas != n_unique_rows:
        LOGGER.error(f'There are {n_unique_rows - n_unique_areas} precursor rows which are not unique!')
        return None

    # If there are no precursors with non-unique areas, simply return precursor df
    if n_unique == n_unique_areas:
        return precursors

    # If there duplicate precursors return None
    if mode == 'e':
        LOGGER.error(f'There are {n_unique_areas - n_unique} non-unique precursor areas!')
        return None

    LOGGER.warning(f'There are {n_unique_areas - n_unique} non-unique precursor areas!')

    precursors = precursors.set_index(keys=key_cols)
    keep_cols = [col.name for col in precursor_cols
                 if col.name in precursors.columns and not col.is_numeric]

    # Select first precursor
    if mode == 'f':
        LOGGER.warning(f'Selecting first occurrence of duplicate precursors...')
        areas = precursors[[col.name for col in precursor_cols if col.is_numeric]].loc[~precursors.index.duplicated()]

        ret = precursors[keep_cols].merge(areas, how='left', left_index=True, right_index=True)
        return ret.reset_index()

    counts = unique_precursors.groupby(key_cols).apply(len, include_groups=False)
    non_unique = precursors[precursors.index.isin(counts[counts > 1].index)]
    unique = precursors[precursors.index.isin(counts[counts == 1].index)]

    # Use precursor peak areas that were manually set.
    if mode == 'm':
        # check UserSetTotal column
        if 'userSetTotal' not in precursors.columns:
            LOGGER.error('Need "UserSetTotal" column to select precursors with user set peak boundaries!')
            return None
        LOGGER.info('Found "UserSetTotal" column.')
        if precursors.dtypes['userSetTotal'] != 'bool':
            LOGGER.error('"UserSetTotal" column must be type bool!')
            return None

        # check that each group has at least 1 precursor area that was manually set
        manually_adj = non_unique.groupby(key_cols)['userSetTotal'].agg(lambda x: any(x))
        if not all(manually_adj):
            LOGGER.error(f'{sum(~manually_adj)} precursor groups have no user set peak boundaries!')
            return None

        # remove precursors which were not manually set
        selections = non_unique[non_unique['userSetTotal']]
        non_unique_len = len(selections.index)
        selections = selections[~selections.index.duplicated(keep='first')]

        if non_unique_len != len(selections.index):
            LOGGER.warning(f'After selecting precursors with user set peak boundaries, '
                           f'{non_unique_len - len(selections.index)} non-unique precursors remain.')

    # Interactive selection mode
    if mode == 'i':
        have_user_set_col = 'userSetTotal' in non_unique.columns
        groups = non_unique.groupby(key_cols)
        n_groups = len(groups)
        int_re = re.compile('[0-9]+')
        selections = []
        for group_i, (name, group) in enumerate(groups):
            valid_choice = False
            sele_int = None
            while not valid_choice:
                # format table of precursor options
                options = [['Selection', 'Protein', 'Sequence', 'Charge', 'Area']]
                if have_user_set_col:
                    options[-1].append('UserSet')
                for i, row in enumerate(group.itertuples()):
                    options.append([f'{i})', row.proteinName, name[1], str(name[2]), str(row.totalAreaFragment)])
                    if have_user_set_col:
                        options[-1].append(str(row.userSetTotal))
                col_widths = [max([len(options[row][col]) for row in range(len(options))]) for col in range(len(options[0]))]

                # write precursor options
                sys.stdout.write(f'Choose precursor ({group_i + 1} of {n_groups})\n')
                sys.stdout.write(f'Replicate: "{name[0]}"\n')
                for row in options:
                    sys.stdout.write('  {}\n'.format('  '.join(row[col].ljust(col_widths[col]) for col in range(len(row)))))

                selection = input('Select: ')
                if int_re.search(selection):
                    sele_int = int(selection)
                    if sele_int >= 0 and sele_int < len(group.index):
                        valid_choice = True
                        break
                sys.stdout.write('Invalid choice!\n')
            selections.append(group.iloc[sele_int])

        selections = pd.DataFrame(selections)
        selections.index = selections.index.set_names(key_cols)

    # overwrite duplicate precursors with selections
    non_unique = non_unique[keep_cols].merge(selections[[col.name for col in precursor_cols if col.is_numeric]],
                                             how='left', left_index=True, right_index=True)

    # Merge duplicates back into main dataframe
    ret = pd.concat([unique, non_unique])
    ret = ret.reset_index()

    return ret


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)
    parser.add_argument('-m', '--metadata', default=None,
                        help='Annotations corresponding to each file.')
    parser.add_argument('--metadataFormat', default=None, choices=('json', 'csv', 'tsv'),
                        help='Specify metadata file format. '
                             'By default the format is inferred from the file extension.')
    # parser.add_argument('--proteins', default=None,
    #                     help='Long formatted protein abundance report. '
    #                           'If no protein report is given, proteins are quantified in the database '
    #                           'by summing all the precursors belonging to that protein.')
    parser.add_argument('-o', '--ofname', default='data.db3',
                        help='Output database name. Default is ./data.db3')
    parser.add_argument('-w', '--overwriteMode', choices=['error', 'overwrite', 'append'], default='error',
                        help='Behavior if output database already exists. '
                             'By default the script will exit with an error if the database file already exists.')
    parser.add_argument('-d', '--duplicatePrecursors', default='e', choices=DUPLICATE_PRECURSOR_CHOICES,
                        help="How to handle precursors with the same sequence and charge, but different area. "
                             "'e' for error, 'm' to use the peak area with user adjusted integration boundaries, "
                             "'f' to select first occurrence of duplicate precursors, "
                             "and 'i' to interactively choose which peak area to use. 'e' is the default.")
    parser.add_argument('--groupBy', choices=['protein', 'gene'], default='protein', dest='group_precursors_by',
                        help="Choose whether to group precursors by gene or protein. Default is 'protein'")
    parser.add_argument('-n', '--projectName', default=None, dest='project_name',
                        help='Project name to use in replicates table.')
    parser.add_argument('replicates', help='Skyline replicate_quality report')
    parser.add_argument('precursors', help='Skyline precursor_quality report')
    return parser.parse_args(argv)


def _main(args):
    '''
    Actual main method. `args` Should be initialized argparse namespace.
    '''

    # read annotations if applicable
    metadata = None
    metadata_types = None
    if args.metadata is not None:
        metadata_reader = Metadata()
        if not metadata_reader.read(args.metadata, args.metadataFormat):
            sys.exit(1)
        metadata = metadata_reader.df
        metadata_types = metadata_reader.types

    # read replicates
    replicates = ReplicateReport().read_report(args.replicates)
    if replicates is None:
        sys.exit(1)

    # read precursor quality report
    precursors = PrecursorReport().read_report(args.precursors, by_gene=(args.group_precursors_by=='gene'))
    if precursors is None:
        sys.exit(1)

    if (precursors := check_duplicate_precursors(precursors, args.duplicatePrecursors)) is None:
        sys.exit(1)

    precursor_reps = set(precursors['replicateName'].drop_duplicates().tolist())
    if not all(k in precursor_reps for k in replicates['replicate'].tolist()):
        LOGGER.warning('Not all replicates in replicate report are in precursor report.')

    # read protein quants
    protein_quants = None
    # if args.proteins:
    #     protein_quants = pd.read_csv(args.proteins, sep='\t')
    #     protein_quants = protein_quants.rename(columns=PROTEIN_QUANTS_REQUIRED_COLUMNS)
    #     if not check_df_columns(protein_quants, PROTEIN_QUANTS_REQUIRED_COLUMNS, 'protein_quants'):
    #         sys.exit(1)

    #     protein_reps = set(protein_quants['replicateName'].drop_duplicates().tolist())
    #     if not all(k in protein_reps for k in replicates['replicate'].tolist()):
    #         LOGGER.warning('Not all replicates in replicate report are in protein quants.')

    #     protein_quants = protein_quants[PROTEIN_QUANTS_REQUIRED_COLUMNS.keys()]
    #     LOGGER.info('Done reading proteins table...')

    # check if file names match replicate names
    if 'fileName' in replicates.columns:
        if any(replicates['replicate'] != replicates['fileName'].apply(lambda x: os.path.splitext(x)[0])):
            LOGGER.warning('Not all file names match replicate names, using FileName instead.')

            # fix names in replicates
            replicates['fileName'] = replicates['fileName'].apply(lambda x: os.path.splitext(x)[0])
            rep_name_map = {row.replicate: row.fileName for row in replicates.itertuples()}
            replicates['replicate'] = replicates['fileName']

            # fix names in precursors
            try:
                precursors['replicateName'] = precursors['replicateName'].apply(lambda x: rep_name_map[x])
                if protein_quants is not None:
                    protein_quants['replicateName'] = protein_quants['replicateName'].apply(lambda x: rep_name_map[x])
            except KeyError as e:
                LOGGER.error(f'KeyError: {str(e)}')
                LOGGER.error('Replicate and precursor report replicate names do not match!')
                sys.exit(1)

    else:
        LOGGER.warning('"FileName" column not found in replciate report, could not check that replicate names match!')

    replicates = replicates[[col.name for col in ReplicateReport().required_columns()]]
    replicates['acquiredRank'] = -1 # this will be populated later

    # build database
    LOGGER.info('Writing database...')
    if not write_db(args.ofname, replicates, precursors,
                    protein_quants=protein_quants,
                    sample_metadata=metadata, sample_metadata_types=metadata_types,
                    project_name=args.project_name, overwriteMode=args.overwriteMode,
                    group_precursors_by=args.group_precursors_by):
        LOGGER.error(f'Failed to {"create" if args.overwriteMode != "append" else "append to"} database!')
        sys.exit(1)
    LOGGER.info('Done writing database...')


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc parse" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()

