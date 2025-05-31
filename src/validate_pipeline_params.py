
import sys
import argparse

from .submodules.read_metadata import Metadata
from .submodules.logger import LOGGER

COMMAND_DESCRIPTION = 'Validate Nextflow pipeline params for the nf-skyline-dia-ms pipeline.'


def parse_input_files(input_files, file_glob='*.raw', api_key=None):
    '''
    Parse quant_spectra_dir or chromatogram_library_spectra_dir parameters
    and return a list or dict of files.

    Parameters
    ----------
    input_files : str, dict, or list
        Input files or directories to parse.
    file_glob : str, optional
        Glob pattern to match files in directories (default: '*.raw').
    api_key : str, optional
        API key for accessing files on Panorama.

    Returns
    -------
    multi_batch_mode : bool
        True if directories have multiple batches, False otherwise.
    files : list or dict
        List of files or in multi-batch mode, a dict mapping batch names to files.
    '''
    pass


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)

    return parser.parse_args(argv)


def _main(args):
    '''
    Actual main method. `args` Should be initialized argparse namespace.
    '''

    metadata_reader = Metadata()
    if not metadata_reader.read(args.metadata_file, args.input_format):
        sys.exit(1)


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc metadata_convert" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()
