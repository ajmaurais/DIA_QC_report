
import sys
import argparse
import sqlite3
import os
from shutil import get_terminal_size

from .submodules.dia_db_utils import validate_bit_mask, parse_bitmask_options
from .submodules.dia_db_utils import check_schema_version
from .submodules.dia_db_utils import PRECURSOR_IMPUTE_METHOD, PROTEIN_IMPUTE_METHOD
from .submodules.dia_db_utils import update_meta_value
from .submodules.dia_db_utils import reset_imputed_values
from .submodules.dia_db_utils import update_command_log
from .submodules.transformation import MethodOptions
from .submodules.imputation import IMPUTATION_METHODS
from .submodules.imputation import KNNImputer
from .submodules.logger import LOGGER

COMMAND_DESCRIPTION = 'Impute missing precursor and/or protein values.'


def get_manager(method, method_args, method_help=False,
                impute_precursors=True, impute_proteins=True,
                **kwargs):
    '''
    Get imputation manager and data processing metadata for imputation method.

    Parameters
    ----------
    method: str
        The name of the imputation method.
    method_args: list
        The list of method specific options strings.
    method_help:
        Show imputation method help and exit?
    impute_precursors: bool
        Impute precursors? Default is True.
    impute_proteins: bool
        Impute proteins? Default is True.
    **kwargs: dict
        Additional arguments passed to all imputation manager for all methods.

    Raises
    ------
    SystemExit:
        If method_help is True show method help and exit with code 0.
        If method_args are invalid exit with code 1.
    '''

    if not (impute_precursors or impute_proteins):
        LOGGER.warning('Both precursor and protein imputation is skipped.')
        sys.exit(0)

    metadata = kwargs.copy()
    parser = MethodOptions()

    if method == 'KNN':
        parser.add_option('n_neighbors', default=5, min_value=1, min_inclusive=True, dtype=int,
                          help_str='Number of neighboring replicates to use for imputation.')
        parser.add_option('weights', choices=('uniform', 'weighted'), default='uniform', dtype=str,
                          help_str='How to weigh the contribution of each neighbor on the '
                                   'imputed value.')
        Manager = KNNImputer
    else:
        raise ValueError(f"Unknown imputation method: '{method}'")

    parser.description = f'Options for {Manager.__name__}'
    if method_help:
        terminal_width = get_terminal_size().columns
        parser.get_help(sys.stdout, max_width=terminal_width)
        sys.exit(0)

    if method_args is not None:
        if not parser.parse_strings(method_args):
            sys.exit(2)

    manager_args = parser.get_option_dict()
    manager = Manager(impute_precursors=impute_precursors, impute_proteins=impute_proteins,
                      **manager_args, **kwargs)
    metadata.update(manager_args)

    ret_metadata = dict()
    if impute_precursors:
        for key, value in metadata.items():
            ret_metadata[f'precursor_imputation_{key}'] = value
            ret_metadata[PRECURSOR_IMPUTE_METHOD] = method
    if impute_proteins:
        for key, value in metadata.items():
            ret_metadata[f'protein_imputation_{key}'] = value
            ret_metadata[PROTEIN_IMPUTE_METHOD] = method

    return manager, ret_metadata


def update_db(conn, imputation_manager):
    '''
    Update database with imputed values in imputation manager.

    Parameters
    ----------
    conn: sqlite3.Connection
    imputation_manager: imputation.ImputationManagerBase
        A child class of ImputationManagerBase where the impute method has already been called.
    '''
    LOGGER.info('Writing imputed values to database...')

    # add imputed values to db
    if imputation_manager.impute_precursors:
        df = imputation_manager.precursors[imputation_manager.precursors['isImputed']]
        prec_data = [(row.area,
                      row.replicateId,
                      row.peptideId,
                      row.precursorCharge) for row in df.itertuples()]
        cur = conn.cursor()
        cur.executemany(f'''
                        UPDATE precursors SET
                            {'totalAreaFragment' if imputation_manager.level == 0 else 'normalizedArea'} = ?,
                            isImputed = 1
                        WHERE replicateId == ? AND
                            peptideId == ? AND
                            precursorCharge == ? ;''', prec_data)
        conn.commit()
        LOGGER.info('Finished writing imputed precursors!')

    if imputation_manager.impute_proteins:
        df = imputation_manager.proteins[imputation_manager.proteins['isImputed']]
        prot_data = [(row.abundance,
                      row.replicateId,
                      row.proteinId) for row in df.itertuples()]
        cur = conn.cursor()
        cur.executemany(f'''
                        UPDATE proteinQuants SET
                            {'abundance' if imputation_manager.level == 0 else 'normalizedAbundance'} = ?,
                            isImputed = 1
                        WHERE replicateId = ? AND proteinId = ? ;''', prot_data)
        conn.commit()
        LOGGER.info('Finished writing imputed proteins!')


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)

    impute_settings = parser.add_argument_group('Global imputation settings')
    impute_settings.add_argument('-i', '--imputeData', default='3', dest='impute_data',
                                 help='One digit bit mask. 1 for precursor imputation, 2 for protein '
                                      'imputation, 3 for both precursor and protein imputation.')
    impute_settings.add_argument('-m', '--method', choices=IMPUTATION_METHODS, default='KNN',
                                 help='Normalization method to use. Default is "KNN"')
    impute_settings.add_argument('-l', '--level', choices=(0, 1), default=0, type=int,
                                 help='Which values to use for imputation. '
                                 '0 for unnormalized, 1 for normalized. Default is 0')
    impute_settings.add_argument('-t', '--missingThreshold', default=0.5,
                                 type=float, dest='missing_threshold',
                                 help='The fraction of missing values above which imputation will be '
                                      'skipped. 0.5 is the default, meaning that if more than 50%% of '
                                      'values are missing in a group, the missing quantities for that '
                                      'peptide/precursor will not be imputed.')
    impute_settings.add_argument('--groupBy', choices=('all', 'project'),
                                 default='project', dest='group_by',
                                 help="How to group data for imputation. "
                                      "If 'project', imputation is performed separately for each project. "
                                      "If 'all' imputation is performed across all replicates. "
                                      "'project' is the default.")

    method_settings = parser.add_argument_group('Imputation method specific settings')
    method_settings.add_argument('-s', '--methodSetting',
                                 default=None, dest='method_settings', action='append',
                                 help='Add method specific settings having the form: <key>=<value>')
    method_settings.add_argument('--methodHelp', action='store_true', default=False,
                                 dest='method_help', help='See imputation method specific settings.')

    parser.add_argument('db', nargs='?', default=None,
                        help='Path to sqlite batch/qc database.')

    args = parser.parse_args(argv)
    if args.db is None and not args.method_help:
        parser.print_usage()
        parser.exit(1, f'{prog}: error: the following arguments are required: db\n')

    return args


def _main(args):
    '''
    Actual main method. `args` Should be an initialized argparse namespace.
    '''

    if not validate_bit_mask(args.impute_data, 2, 1):
        LOGGER.error('Error parsing --imputeData')
        sys.exit(1)

    # parse bitmask options
    impute_data = parse_bitmask_options(args.impute_data, ('impute',), ('precursors', 'proteins'))
    impute_precursors = impute_data['impute']['precursors']
    impute_proteins = impute_data['impute']['proteins']

    # initialize imputation manager
    imputation_manager, metadata = get_manager(args.method, args.method_settings,
                                               method_help=args.method_help,
                                               impute_precursors=impute_precursors,
                                               impute_proteins=impute_proteins,
                                               level=args.level,
                                               missing_threshold=args.missing_threshold,
                                               group_by_project=args.group_by == 'project')

    if os.path.isfile(args.db):
        conn = sqlite3.connect(args.db)
    else:
        LOGGER.error('Database file (%s) does not exist!', args.db)
        sys.exit(1)

    # check database version
    if not check_schema_version(conn):
        sys.exit(1)

    # delete existing imputed values in db
    reset_imputed_values(conn, precursors=impute_precursors, proteins=impute_proteins)

    # add db connection to imputation_manager and impute
    imputation_manager.conn = conn
    if not imputation_manager.impute():
        sys.exit(1)

    # add imputed values from manager to db
    update_db(conn, imputation_manager)

    # update commandLog
    update_command_log(conn, sys.argv, os.getcwd())

    # Update normalization method in metadata
    LOGGER.info('Updating metadata...')
    for key, value in metadata.items():
        update_meta_value(conn, key, value)

    conn.close()


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc impute" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()
