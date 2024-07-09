
import sys
import argparse
from os.path import basename, splitext
from csv import DictReader

from .submodules.dtype import Dtype
from .submodules.read_metadata import Metadata
from .submodules.logger import LOGGER

COMMAND_DESCRIPTION = 'Convert metadata annotations to specified file format.'
IN_FORMATS = ('json', 'csv', 'tsv')
OUT_FORMATS = ('json', 'skyline', 'csv', 'tsv')


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)
    parser.add_argument('-i', '--in', default=None,
                        choices=IN_FORMATS, dest='input_format',
                        help='Specify metadata file format. '
                             'By default the format is inferred from the file extension.')
    parser.add_argument('-o', '--out', default='skyline',
                        choices=OUT_FORMATS, dest='output_format',
                        help="Specify metadata output format. 'skyline' is the default.")
    parser.add_argument('-p', '--prefix', default=None, dest='out_prefix',
                        help='Output file prefix. By default the input file basename is used.')
    parser.add_argument('metadata_file', help='The metadata file to convert.')

    return parser.parse_args(argv)


def _main(args):
    '''
    Actual main method. `args` Should be initialized argparse namespace.
    '''

    output_prefix = args.out_prefix if args.out_prefix else splitext(basename(args.metadata_file))[0]

    metadata_reader = Metadata()
    if not metadata_reader.read(args.metadata_file, args.input_format):
        sys.exit(1)

    if args.output_format in ('tsv', 'json', 'csv'):
        with open(f'{output_prefix}.{args.output_format}', 'w') as outF:
            getattr(metadata_reader, f'to_{args.output_format}')(outF)

    elif args.output_format == 'skyline':
        # write skyline annotation csv
        with open(f'{output_prefix}.annotations.csv', 'w') as outF:
            metadata_reader.to_skyline_annotations(outF)

        # write annotation definition batch file
        with open(f'{output_prefix}.definitions.bat', 'w') as outF:
            metadata_reader.to_skyline_definitions(outF)

    else:
        LOGGER.error('Unknown output format!')
        sys.exit(1)


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc metadata_convert" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()
