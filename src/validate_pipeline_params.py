
import sys
import os
import argparse
import re
import copy
from pathlib import PurePosixPath
from urllib.parse import quote


from .submodules import nextflow_pipeline_config
from .submodules.panorama import PANORA_PUBLIC_KEY
from .submodules.panorama import list_panorama_files, download_webdav_file
from .submodules.read_metadata import Metadata
from .submodules.logger import LOGGER

COMMAND_DESCRIPTION = 'Validate Nextflow pipeline params for the nf-skyline-dia-ms pipeline.'

DEFAULT_PIPELINE = 'mriffle/nf-skyline-dia-ms'
DEFAULT_PIPELINE_REVISION = 'main'


def glob_to_regex(file_glob: str) -> str:
    '''
    Convert a Nextflow-style glob (only * is a wildcard) into an anchored-regex string

    Parameters
    ----------
    file_glob : str
        Pattern such as "*.txt" or "Replicate_*_R1.mzML"

    Returns
    -------
    str
        Regex string beginning with ^ and ending with $.
    '''

    # Any metacharacter except * should be taken literally.
    def escape_regex(s: str) -> str:
        return re.sub(
            r'([.\^$+?{}\[\]\\|()])',     # chars to escape 1-for-1
            lambda m: r'\{}'.format(m.group(1)),
            s,
        )

    # build the final pattern exactly as in the workflow
    regex_body = escape_regex(file_glob).replace('*', '.*')
    return f'^{regex_body}$'


def _process_directorty(directory, file_regex, api_key=None):
    if directory.startswith('https://panoramaweb.org'):
        files = list_panorama_files(directory, api_key=api_key)
    else:
        files = os.listdir(directory)

    files = [file for file in files if re.search(file_regex, file)]
    if len(files) == 0:
        LOGGER.error(f'No files found in {directory} matching {file_regex}')
    return files


def parse_input_files(input_files, file_glob=None, file_regex=None, api_key=None):
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
    if file_regex is not None and file_glob is not None:
        raise ValueError('Only one of file_regex or file_glob should be provided.')

    _file_regex = file_regex or glob_to_regex(file_glob)

    if isinstance(input_files, str):
        input_files = [line.strip() for line in input_files.split('\n') if line.strip()]

    if isinstance(input_files, list):
        multi_batch_mode = False
        files = []
        for item in input_files:
            files.extend(_process_directorty(item, _file_regex, api_key))
    elif isinstance(input_files, dict):
        multi_batch_mode = True
        files = {}
        for batch_name, batch_dir in input_files.items():
            if isinstance(batch_dir, str):
                batch_dirs = [line.strip() for line in batch_dir.split('\n') if line.strip()]
            elif isinstance(batch_dir, list):
                batch_dirs = batch_dir
            else:
                raise TypeError('Batch directories must be a str or list.')

            files[batch_name] = []
            for dir in batch_dirs:
                files[batch_name] += _process_directorty(dir, _file_regex, api_key)
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
