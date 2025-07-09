
import sys
import argparse
from os.path import basename, splitext
from pandas import isna

from .submodules.read_metadata import Metadata
from .submodules.logger import LOGGER

COMMAND_DESCRIPTION = 'Convert metadata annotations to specified file format.'
IN_FORMATS = ('json', 'csv', 'tsv')
OUT_FORMATS = ('json', 'skyline', 'csv', 'tsv')


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)
    parser.add_argument(
        '-i', '--in', default=None, choices=IN_FORMATS, dest='input_format',
        help='Specify metadata file format. '
             'By default the format is inferred from the file extension.'
    )
    parser.add_argument(
        '-o', '--out', default='skyline', choices=OUT_FORMATS, dest='output_format',
        help="Specify metadata output format. 'skyline' is the default."
    )
    parser.add_argument(
        '-p', '--prefix', default=None, dest='out_prefix',
        help='Output file prefix. By default the input file basename is used.'
    )
    parser.add_argument(
        'metadata_file', nargs='+',
        help='The metadata file(s) to convert. If multiple files are given, they will be '
             'concatenated into a single output file.'
    )

    return parser.parse_args(argv)


def _warnings_per_file(combined, metadata_files):
    var_na_counter = {var: 0 for var in combined.types}
    for file_name, reader in metadata_files.items():
        for var, dtype in combined.types.items():
            if var not in reader.types:
                LOGGER.warning(
                    "Variable '%s' found in combined metadata but not in file '%s'",
                    var, file_name
                )
                continue

            if reader.types[var] != dtype:
                LOGGER.warning(
                    "Variable '%s' in file '%s' has type '%s', after merging it has type '%s'",
                    var, file_name, reader.types[var], dtype
                )

        for rep in reader.df.itertuples():
            if isna(rep.annotationValue) or rep.annotationValue == '':
                var_na_counter[rep.annotationKey] += 1

    n_reps = len(combined.df['Replicate'].drop_duplicates())
    for var, count in var_na_counter.items():
        if count > 0:
            LOGGER.warning(
                "Variable '%s' has missing values in %d of %d replicates.",
                var, count, n_reps
            )


def _main(argv, prog=None):
    ''' Actual main method. '''
    args = parse_args(argv, prog=prog)

    if len(args.metadata_file) == 1:
        output_prefix = args.out_prefix if args.out_prefix else splitext(basename(args.metadata_file[0]))[0]
    else:
        output_prefix = args.out_prefix if args.out_prefix else 'combined_metadata'

    # read metadata file(s)
    metadata_files = {}
    for metadata_file in args.metadata_file:
        metadata_reader = Metadata()
        if not metadata_reader.read(metadata_file, args.input_format):
            sys.exit(1)
        metadata_files[metadata_file] = metadata_reader

    # combine metadata files
    combined_metadata = Metadata()
    for reader in metadata_files.values():
        combined_metadata += reader

    _warnings_per_file(combined_metadata, metadata_files)
    for var, dtype in combined_metadata.types.items():
        LOGGER.info("Variable '%s' has type '%s'", var, dtype)

    if args.output_format in ('tsv', 'json', 'csv'):
        with open(f'{output_prefix}.{args.output_format}', 'w') as outF:
            getattr(combined_metadata, f'to_{args.output_format}')(outF)
        LOGGER.info('Metadata written to %s.%s', output_prefix, args.output_format)

    elif args.output_format == 'skyline':
        # write skyline annotation csv
        with open(f'{output_prefix}.annotations.csv', 'w') as outF:
            combined_metadata.to_skyline_annotations(outF)
        LOGGER.info('Skyline annotations written to %s.annotations.csv', output_prefix)

        # write annotation definition batch file
        with open(f'{output_prefix}.definitions.bat', 'w') as outF:
            combined_metadata.to_skyline_definitions(outF)
        LOGGER.info('Skyline annotation definitions written to %s.definitions.bat', output_prefix)

    else:
        LOGGER.error('Unknown output format!')
        sys.exit(1)


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc metadata_convert" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()
