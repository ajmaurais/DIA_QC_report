
import sys
import argparse
import sqlite3
import os

from .submodules.dia_db_utils import check_schema_version
from .submodules.dia_db_utils import METADATA_TIME_FORMAT
from .submodules.dia_db_utils import IS_IMPUTED, PRECURSOR_IMPUTE_METHOD, PROTEIN_IMPUTE_METHOD
from .submodules.dia_db_utils import update_meta_value
from .submodules.dia_db_utils import validate_bit_mask, parse_bitmask_options
from .submodules.transformation import MethodOptions
from .submodules.imputation import IMPUTATION_METHODS
from .submodules.imputation import KNNImputer
from .submodules.logger import LOGGER
from . import __version__ as PROGRAM_VERSION

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
        parser.get_help(sys.stdout, max_width=os.get_terminal_size().columns)
        sys.exit(0)

    if method_args is not None:
        if not parser.parse_strings(method_args):
            sys.exit(1)

    manager_args = parser.get_option_dict()
    manager = Manager(impute_precursors=impute_precursors, impute_proteins=impute_proteins,
                      **manager_args, **kwargs)
    metadata.update(manager_args)

    ret_metadata = dict()
    if impute_precursors:
        for key, value in metadata.items():
            ret_metadata[f'precursor_imputation_{key}'] = value
    if impute_proteins:
        for key, value in metadata.items():
            ret_metadata[f'protein_imputation_{key}'] = value

    return manager, ret_metadata


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)

    impute_settings = parser.add_argument_group('Global imputation settings')
    impute_settings.add_argument('-i', '--imputeData', default='3', dest='impute_data',
                                 help='One digit bit mask. 1 for precursor imputation, 2 for protein '
                                      'imputation, 3 for both precursor and protein imputation.')
    impute_settings.add_argument('-m', '--method', choices=IMPUTATION_METHODS, default='KNN',
                                 help='Normalization method to use. Default is "k-means"')
    impute_settings.add_argument('-l', '--level', choices=(0, 1), default=0,
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

    if not validate_bit_mask(args.interactive, 2, 1):
        LOGGER.error('Error parsing --imputeData')
        sys.exit(1)

    # parse bitmask options
    impute_data = parse_bitmask_options(args.interactive, ('impute',), ('precursors', 'proteins'))

    # initialize imputation manager
    imputation_manager, metadata = get_manager(args.method, args.method_settings,
                                               method_help=args.method_help,
                                               impute_precursors=impute_data['impute']['precursors'],
                                               impute_proteins=impute_data['impute']['proteins'],
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

    # add db connection to imputation_manager
    imputation_manager.conn = conn

    imputation_manager.impute()


    # get commands previously run on db

    # Update normalization method in metadata
    LOGGER.info('Updating metadata...')
    # metadata = {PRECURSOR_IMPUTE_METHOD: precursor_normalization_method,
    #             PROTEIN_IMPUTE_METHOD: protein_normalization_method,
    #             IS_IMPUTED: 'True'}
    for key, value in metadata.items():
        conn = update_meta_value(conn, key, value)

    conn.close()


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc impute" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()

