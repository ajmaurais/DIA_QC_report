
import sys
import os
import argparse
import sqlite3

import pandas as pd

from .submodules.logger import LOGGER
from .submodules.dia_db_utils import is_normalized
from .submodules.dia_db_utils import check_schema_version
from .submodules.dia_db_utils import validate_bit_mask, parse_bitmask_options

COMMAND_DESCRIPTION = 'Export selected table(s) from precursor database.'

METHOD_NAMES = ['unnormalized', 'normalized']
PRECURSOR_METHOD_NAMES = dict(zip(METHOD_NAMES, ['totalAreaFragment', 'normalizedArea']))
PROTEIN_METHOD_NAMES = dict(zip(METHOD_NAMES, ['abundance', 'normalizedAbundance']))
TABLE_TYPES = ('wide', 'long')


def write_precusor_tables(conn, dest, tables):
    '''
    Write precursor table(s)

    Parameters
    ----------
    conn: sqlite3.Connection
        An initialized connection to the precursor database.
    dest: str
        The Parent directory
    tables: dict
        A nested dictionary where the top level dictionary is the table directions, 'wide' or 'long'
        and the child dictionaries are the method names, 'unnormalized' or 'normalized'
    '''

    df = pd.read_sql('''
    SELECT
        r.replicate,
        prot.name as protein,
        p.modifiedSequence,
        p.precursorCharge,
        p.totalAreaFragment,
        p.normalizedArea
    FROM precursors p
    LEFT JOIN replicates r ON r.id == p.replicateId
    LEFT JOIN peptideToProtein ptp ON ptp.peptideId == p.peptideId
    LEFT JOIN proteins prot ON prot.proteinId == ptp.proteinId;''', conn)

    # long tables
    if any(v for v in tables['long'].values()):
        long_headers = ['replicate', 'protein', 'modifiedSequence', 'precursorCharge']
        long_headers += [PRECURSOR_METHOD_NAMES[k] for k, v in tables['long'].items() if v]

        df[long_headers].to_csv(f'{dest}/precursors_long.tsv', sep='\t', index=False)

    # wide tables, if any
    for method, write in tables['wide'].items():
        if write:
            df.pivot(index=['protein', 'modifiedSequence', 'precursorCharge'],
                     columns='replicate',
                     values=PRECURSOR_METHOD_NAMES[method]).to_csv(f'{dest}/precursors_wide_{method}.tsv',
                                                                       sep='\t', index=True)


def write_protein_tables(conn, dest, tables):
    '''
    Write precursor table(s)

    Parameters
    ----------
    conn: sqlite3.Connection
        An initialized connection to the precursor database.
    dest: str
        The Parent directory
    tables: dict
        A nested dictionary where the top level dictionary is the table directions, 'wide' or 'long'
        and the child dictionaries are the method names, 'unnormalized' or 'normalized'
    '''

    df = pd.read_sql('''
    SELECT
        r.replicate,
        prot.name as protein,
        q.abundance,
        q.normalizedAbundance
    FROM proteinQuants q
    LEFT JOIN replicates r ON r.id == q.replicateId
    LEFT JOIN proteins prot ON prot.proteinId == q.proteinId;''', conn)

    # long tables
    if any(v for v in tables['long'].values()):
        long_headers = ['replicate', 'protein']
        long_headers += [PROTEIN_METHOD_NAMES[k] for k, v in tables['long'].items() if v]

        df[long_headers].to_csv(f'{dest}/proteins_long.tsv', sep='\t', index=False)

    # wide tables, if any
    for method, write in tables['wide'].items():
        if write:
            df.pivot(index='protein', columns='replicate',
                     values=PROTEIN_METHOD_NAMES[method]).to_csv(f'{dest}/proteins_wide_{method}.tsv',
                                                                 sep='\t', index=True)


def write_metadata_tables(conn, dest, tables):
    df = pd.read_sql('''
    SELECT
        r.replicate,
        m.annotationKey as key,
        m.annotationValue as value
    FROM sampleMetadata m
    LEFT JOIN replicates r ON r.id == m.replicateId;''', conn)

    if tables['long']['write']:
        df.to_csv(f'{dest}/metadata_long.tsv', sep='\t', index=False)

    if tables['wide']['write']:
        df.pivot(index='replicate',
                columns='key',
                values='value').to_csv(f'{dest}/metadata_wide.tsv', sep='\t', index=False)


def _any_tables(table_opts):
    return any(table_opts[d][m] for d in table_opts for m in table_opts[d])


def parse_args(argv):
    parser = argparse.ArgumentParser(description=COMMAND_DESCRIPTION)
    parser.add_argument('-o', '--outputDir', default=None, dest='output_dir',
                        help=f'Output directory. Default is the current working directory.')

    table_args = parser.add_argument_group('Output tables',
                     'The tsv files to write are specified by a 2 digit bit mask. '
                     'The first digit is for the wide formatted report, and the second digit is for '
                     'the long formatted report. An integer between 0 and 2 is assigned for each stage in '
                     'the normalization process. 0 for no report, 1 is for unnormalized, '
                     'and 2 is for normalized.')

    table_args.add_argument('-p', '--precursorTables', default='20',
                        help='Tables to write for precursors. 20 is the default')
    table_args.add_argument('-r', '--proteinTables', default='20',
                            help='Tables to write for proteins. 20 is the default')
    table_args.add_argument('-m', '--metadataTables', default='00',
                            help='Tables to write for metadata. Only 0 or 1 are supported. '
                                 '0 for false, 1 for true. 00 is the default')

    parser.add_argument('db', help='Path to precursor quality database.')
    return parser.parse_args(argv)


def _main(args):
    '''
    Actual main method. `args` Should be initialized argparse namespace.
    '''

    # check args
    if not validate_bit_mask(args.precursorTables, 2, 2):
        LOGGER.error('Error parsing --precursorTables')
        sys.exit(1)
    if not validate_bit_mask(args.proteinTables, 2, 2):
        LOGGER.error('Error parsing --proteinTables')
        sys.exit(1)
    if not validate_bit_mask(args.metadataTables, 1, 2):
        LOGGER.error('Error parsing --metadataTables')
        sys.exit(1)

    if os.path.isfile(args.db):
        conn = sqlite3.connect(args.db)
    else:
        LOGGER.error(f"Database file: '{args.db}' does not exist!")
        sys.exit(1)

    # check database version
    if not check_schema_version(conn):
        sys.exit(1)

    # Check if metadata.is_normalized is set to True
    db_is_normalized = is_normalized(conn)

    # create output directory if applicable
    if args.output_dir:
        if not os.path.isdir(args.output_dir):
            output_dir = args.output_dir
            try:
                os.mkdir(args.output_dir)
            except (FileExistsError, FileNotFoundError) as e:
                LOGGER.error(f'Could not create output directory: {e}')
                sys.exit(1)
            LOGGER.info(f'Created tables output directory: {output_dir}')
        LOGGER.info(f'Output directory already exists.')
    else:
        output_dir = os.getcwd()
        print(output_dir)

    # parse table_args
    protein_tables = parse_bitmask_options(args.proteinTables, TABLE_TYPES, METHOD_NAMES)
    precursor_tables = parse_bitmask_options(args.precursorTables, TABLE_TYPES, METHOD_NAMES)
    metadata_tables = parse_bitmask_options(args.metadataTables, TABLE_TYPES, ('write',))

    for table_type, table, in [['protein', protein_tables], ['precursor', precursor_tables]]:
        if any(table[direction]['normalized'] for direction in TABLE_TYPES) and not db_is_normalized:
            LOGGER.warning(f'Database is not normalized. Skipping normalized {table_type} table(s).')
            for direction in TABLE_TYPES:
                table[direction]['normalized'] = False

    if _any_tables(precursor_tables):
        LOGGER.info('Writing precursor tables..')
        write_precusor_tables(conn, output_dir, precursor_tables)
        LOGGER.info('Done writing precursor tables.')

    if _any_tables(protein_tables):
        LOGGER.info('Writing protein tables..')
        write_protein_tables(conn, output_dir, protein_tables)
        LOGGER.info('Done writing protein tables.')

    if _any_tables(metadata_tables):
        LOGGER.info('Writing metadata tables..')
        write_metadata_tables(conn, output_dir, metadata_tables)
        LOGGER.info('Done writing metadata tables.')

    conn.close()


def main():
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()
