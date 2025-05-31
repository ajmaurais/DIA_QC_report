
import sys
import os
import argparse
from fnmatch import fnmatch
from pathlib import PurePosixPath
from urllib.parse import quote


from .submodules import nextflow_pipeline_config
from .submodules.panorama import list_panorama_files, download_webdav_file
from .submodules.read_metadata import Metadata
from .submodules.logger import LOGGER

COMMAND_DESCRIPTION = 'Validate Nextflow pipeline params for the nf-skyline-dia-ms pipeline.'

DEFAULT_PIPELINE = 'mriffle/nf-skyline-dia-ms'
DEFAULT_PIPELINE_REVISION = 'main'

def _process_directorty(directory, file_glob='*.raw', api_key=None):
    if directory.startswith('https://panoramaweb.org'):
        files = list_panorama_files(directory, api_key=api_key, panorama_pubilc=True)
    else:
        files = os.listdir(directory)

    files = fnmatch.filter(files, file_glob)
    if len(files) == 0:
        LOGGER.error(f'No files found in {directory} matching {file_glob}.')
    return files


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
    if isinstance(input_files, str):
        multi_batch_mode = False
        files = _process_directorty(input_files, file_glob, api_key)
    elif isinstance(input_files, list):
        multi_batch_mode = False
        files = []
        for item in input_files:
            files.extend(_process_directorty(item, file_glob, api_key))
    elif isinstance(input_files, dict):
        multi_batch_mode = True
        files = {}
        for batch_name, batch_dir in input_files.items():
            files[batch_name] = _process_directorty(batch_dir, file_glob, api_key)
    else:
        LOGGER.error(f'Unsupported input_files type: {type(input_files)}')
        raise TypeError('Parameter must be a str, list, or dict.')

    return multi_batch_mode, files


def generate_schema_url(repo: str, revision: str, filename='nextflow_schema.json') -> str:
    '''
    Build the permanent raw-file download URL for filename hosted on GitHub.

    Parameters
    ----------
    repo : str
        Repository name formated as "user/repo" (e.g., "mriffle/nf-skyline-dia-ms").
    revision : str
        Any valid ref—branch name ("main"), tag ("v1.2.0"), or full / short commit SHA.
    filename : str, optional
        Path to the file in the repository (default: 'nextflow_schema.json').

    Returns
    -------
    str
        Direct URL of the form
    '''
    # Normalise and percent-encode every path component except the slashes
    clean_path = PurePosixPath(filename).as_posix()
    encoded_parts = (quote(part, safe="") for part in clean_path.split("/"))
    encoded_path = "/".join(encoded_parts)

    return f"https://raw.githubusercontent.com/{quote(repo, safe='')}/{quote(revision, safe='')}/{encoded_path}"


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
