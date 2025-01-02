
import sys
import argparse
import sqlite3
import os

from .submodules.dia_db_utils import check_schema_version
from .submodules.dia_db_utils import METADATA_TIME_FORMAT
from .submodules.dia_db_utils import update_meta_value
from .submodules.imputation import IMPUTATION_METHODS
from .submodules.logger import LOGGER
from . import __version__ as PROGRAM_VERSION

COMMAND_DESCRIPTION = 'Impute missing precursor and/or protein values.'


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)

    norm_settings = parser.add_argument_group('Imputation settings')
    norm_settings.add_argument('-m', '--method', choices=IMPUTATION_METHODS, default='k-means',
                               help='Normalization method to use. Default is "k-means"')
    norm_settings.add_argument('-l', '--level', choices=(0, 1), default=0,
                               help='Which values to use for imputation. '
                                    '0 for unnormalized, 1 for normalized.')
    norm_settings.add_argument('--methodSettings', default=None, dest='method_settings',
                               help='Imputation method specific settings having the form <key>=<value>.')
    norm_settings.add_argument('--methodHelp', action='store_true', default=False,
                               help='See imputation method specific settings.')

    parser.add_argument('db', help='Path to sqlite batch/qc database.')

    return parser.parse_args(argv)


def _main(args):
    '''
    Actual main method. `args` Should be an initialized argparse namespace.
    '''

    # init imputation manager

    # get commands previously run on db

    # Update normalization method in metadata
    LOGGER.info('Updating metadata...')
    metadata = {'precursor_imputation_method': precursor_normalization_method,
                'protein_imputation_method': protein_normalization_method,
                'is_imputed': 'True'}
    for key, value in metadata.items():
        conn = update_meta_value(conn, key, value)

    conn.close()


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc impute" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()

