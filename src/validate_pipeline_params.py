
import sys
import os
import argparse
import re
import copy
from pathlib import PurePosixPath
from urllib.parse import quote
import json
from shutil import which
import subprocess
from io import StringIO
from collections import Counter
from types import SimpleNamespace

from jsonschema import validate, ValidationError
from requests import HTTPError
from pandas import isna as pd_isna

from .submodules.nextflow_pipeline_config import parse_params, namespace_to_dict
from .submodules.panorama import PANORA_PUBLIC_KEY
from .submodules.panorama import list_panorama_files, get_webdav_file
from .submodules.panorama import get_http_file
from .submodules.read_metadata import Metadata
from .submodules.dtype import Dtype
from .submodules.logger import LOGGER

COMMAND_DESCRIPTION = 'Validate Nextflow pipeline params for the nf-skyline-dia-ms pipeline.'

DEFAULT_PIPELINE = 'mriffle/nf-skyline-dia-ms'
DEFAULT_PIPELINE_REVISION = 'main'

QUANT_SPECTRA_SCHEMA = {
    "oneOf": [
        {
            "type": "string",
            'minProperties': 1
        },
        {
            "type": "array",
            "items": {"type": "string"},
            'minProperties': 1
        },
        {
            "type": "object",
            "additionalProperties": {
                "oneOf": [
                    {
                        "type": "string",
                        'minProperties': 1
                    },
                    {
                        "type": "array",
                        "items": {"type": "string"},
                        'minProperties': 1
                    },
                ]
            },
            'minProperties': 1
        },
    ],
}


def merge_params(lhs, rhs):
    """
    Recursively “add” two params objects together.

    Both *lhs* and *rhs* may be either dicts or SimpleNamespace
    objects produced by `parse_params`.
    """

    if not isinstance(lhs, (dict, SimpleNamespace)) \
       or not isinstance(rhs, (dict, SimpleNamespace)):
        raise TypeError('both arguments must be dict- or namespace-like')

    def _to_mapping(obj):
        """Return the underlying dict for a mapping-like object."""
        return vars(obj) if isinstance(obj, SimpleNamespace) else obj

    def _wrap(template, mapping):
        """
        Wrap *mapping* back into the same container type as *template*
        (namespace or plain dict).
        """
        return (SimpleNamespace(**mapping)
                if isinstance(template, SimpleNamespace)
                else mapping)

    def _merge(a, b):
        if isinstance(a, (dict, SimpleNamespace)) \
           and isinstance(b, (dict, SimpleNamespace)):
            map_a, map_b = _to_mapping(a), _to_mapping(b)

            merged = {}
            for key in map_a.keys() | map_b.keys():          # union of keys
                if key in map_a and key in map_b:
                    merged[key] = _merge(map_a[key], map_b[key])
                elif key in map_a:
                    merged[key] = copy.deepcopy(map_a[key])
                else:
                    merged[key] = copy.deepcopy(map_b[key])

            # If either side was a namespace, keep the result a namespace
            if isinstance(a, SimpleNamespace) or isinstance(b, SimpleNamespace):
                return SimpleNamespace(**merged)
            return merged
        else:
            # scalars / lists / literal maps → rhs wins (deep-copied)
            return copy.deepcopy(b)

    return _merge(lhs, rhs)


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
        return re.sub( r'([.\^$+?{}\[\]\\|()])', lambda m: r'\{}'.format(m.group(1)), s)

    # build the final pattern exactly as in the workflow
    regex_body = escape_regex(file_glob).replace('*', '.*')
    return f'^{regex_body}$'


def get_api_key_from_nextflow_secrets(name='PANORAMA_API_KEY'):
    '''
    Get Panorama API key from Nextflow secrets manager.

    Returns
    -------
    str or None
        Panorama API key if available, None otherwise.
    '''
    nextflow_exe = which('nextflow')
    if nextflow_exe is None:
        LOGGER.error('Nextflow executable not found in PATH.')
        return None

    result = subprocess.run(
        [nextflow_exe, 'secrets', 'get', name],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        LOGGER.error(f"Nextflow secret retrieval failed: {result.stderr.strip()}")
        return None

    api_key = result.stdout.strip()
    if api_key is None or api_key == 'null':
        LOGGER.error("Key '%s' not found in Nextflow secrets.", name)
        return None

    return api_key


def get_file(file_path, log_name='file', api_key=None,
             dest_path=None, return_text=True):
    '''
    Return the text content of a local or remote file.

    Parameters
    ----------
    file_path : str
        Path to the file or Panorama URL.
    log_name : str, optional
        Name to use in log messages (default: 'file').
    api_key : str
        API key for accessing files on Panorama.
    dest_path : str, optional
        If the file is downloaded, save it to this path.
        If None, the base name of the file in the URL will be used.
    return_text : bool, optional
        If True, return the text content of the file. If False, return the file path (default: True).

    Returns
    -------
    str or None
        Text content of the file if it exists, None otherwise.
    '''
    if file_path.startswith('https://panoramaweb.org'):
        try:
            text = get_webdav_file(file_path, api_key=api_key, dest_path=dest_path, return_text=return_text)
        except HTTPError as e:
            LOGGER.error("Failed to download %s '%s': %s", log_name, file_path, e)
            return None

        if not return_text: # return file path
            return text

    elif file_path.startswith('http://') or file_path.startswith('https://'):
        try:
            text = get_http_file(file_path, dest_path=dest_path, return_text=return_text)
        except HTTPError as e:
            LOGGER.error("Failed to download %s '%s': %s", log_name, file_path, e)
            return None

    else:
        file_path = os.path.abspath(file_path)
        if not os.path.exists(file_path):
            LOGGER.error("%s '%s' does not exist.", log_name[0].upper() + log_name[1:], file_path)
            return None
        if return_text:
            with open(file_path, 'r') as schema_file:
                text = schema_file.read()
        else:
            return file_path

    return text


def _process_directorty(directory, file_regex, api_key=None):
    if directory.startswith('https://panoramaweb.org'):
        files = list_panorama_files(directory, api_key=api_key)
    else:
        files = os.listdir(directory)

    files = [file for file in files if re.search(file_regex, file)]
    if len(files) == 0:
        raise RuntimeError(f'No files found in {directory} matching {file_regex}')
    return files


def parse_input_files(input_files, file_glob=None, file_regex=None, api_key=None, strict=True):
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
    strict : bool, optional
        If True, exit with error if any listing any directory fails.

    Returns
    -------
    multi_batch_mode : bool
        True if directories have multiple batches, False otherwise.
    files : list or dict
        List of files or in multi-batch mode, a dict mapping batch names to files.

    Raises
    ------
    SystemExit
        If strict is True and any directory listing fails, the program exits with code 1.
    '''
    if file_regex is not None and file_glob is not None:
        raise ValueError('Only one of file_regex or file_glob should be provided.')

    _file_regex = file_regex or glob_to_regex(file_glob)

    if isinstance(input_files, str):
        input_files = [line.strip() for line in input_files.split('\n') if line.strip()]

    missing_dirs = 0
    if isinstance(input_files, list):
        multi_batch_mode = False
        files = []
        for item in input_files:
            try:
                files.extend(_process_directorty(item, _file_regex, api_key))
            except RuntimeError as e:
                LOGGER.error(str(e))
                missing_dirs += 1

        n_files = len(files)

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
                try:
                    files[batch_name] += _process_directorty(dir, _file_regex, api_key)
                except RuntimeError as e:
                    LOGGER.error(str(e))
                    missing_dirs += 1

        n_files = sum(len(batch_files) for batch_files in files.values())

    else:
        raise TypeError(f'Parameter must be a str, list, or dict. Got: {type(input_files)}')

    if missing_dirs > 0:
        if strict:
            sys.exit(1)
        else:
            dir_str = 'directory' if missing_dirs == 1 else 'directories'
            LOGGER.warning('No MS files could be found in %d %s. Continuing with %d available files.',
                           missing_dirs, dir_str, n_files)

    return multi_batch_mode, files


def generate_git_url(repo: str, revision: str, filename='nextflow_schema.json') -> str:
    '''
    Build the permanent raw-file download URL for filename hosted on GitHub.

    Parameters
    ----------
    repo : str
        Repository name formatted as "user/repo" (e.g., "mriffle/nf-skyline-dia-ms").
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

    base_url = 'https://raw.githubusercontent.com'
    return f"{base_url}/{quote(repo, safe='')}/{quote(revision, safe='')}/{encoded_path}"


def _remove_none_from_param_dict(data):
    ''' Recursively remove any keys in dictionaries whose value is None.  '''
    if isinstance(data, dict):
        cleaned = {}
        for key, value in data.items():
            if value is None:
                continue
            if isinstance(value, dict):
                nested = _remove_none_from_param_dict(value)
                if nested:
                    cleaned[key] = nested
            else:
                cleaned[key] = value
        return cleaned
    return data


def validate_config_files(config_paths, schema_path, api_key=None):
    '''
    Read, merge, and validate Nextflow pipeline config files against the schema.

    Parameters
    ----------
    config_files : list of str
        List of local file paths to Nextflow pipeline config files.
    schema_path : str
        URL or local path of the JSON schema to validate the config files against.

    Returns
    -------
    bool
        True pipeline config is valid, False otherwise.
    dict
        Merged config parameters from all config files.
    '''
    config_data = {}
    all_good = True
    if isinstance(config_paths, str):
        config_paths = [config_paths]

    for config_path in config_paths:
        config_text = get_file(
            config_path, log_name='pipeline config file', api_key=api_key
        )
        if config_text is None:
            all_good = False
            continue
        this_config = parse_params(text=config_text)
        config_data = merge_params(config_data, this_config)

    if not all_good:
        return False, {}

    schema_text = get_file(
        schema_path, log_name='pipeline config schema file', api_key=api_key
    )
    if schema_text is None:
        return False, config_data
    schema = json.loads(schema_text)

    try:
        param_dict = _remove_none_from_param_dict(namespace_to_dict(config_data))
        validate(param_dict, schema)
    except ValidationError as e:
        LOGGER.error(f'Pipeline config validation failed: {e.message}')
        return False, config_data

    return True, config_data


class Replicate:
    '''
    Class representing a single replicate in the metadata.
    '''
    def __init__(self, name, batch=None):
        self.name = name
        self.batch = batch
        self.metadata = {}

    def __repr__(self):
        return f'Replicate(name={self.name}, metadata={self.metadata})'


def _write_metadata_report(meta_params, metadata_types, rep_na_counts, n_reps,
                           output_path='metadata_validation_report.json'):
    '''
    Write a report of metadata parameters and types to a file.

    Parameters
    ----------
    meta_params : list
        List of tuples containing metadata variable and parameter names.
    metadata_types : dict
        Dictionary of metadata variable names to their dtypes.
    rep_na_counts : dict
        Dictionary mapping variables to the number of missing values in replicates.
    n_reps : int
        Total number of replicates.
    output_path : str, optional
        Path to the output report file. Default is 'metadata_validation_report.json'.

    Raises
    ------
    ValueError
        If the output report format is not 'json' or 'tsv'.
    '''
    if len(meta_params) == 0:
        LOGGER.warning('No metadata parameters to report.')
        return

    report_data = []
    for var, name in meta_params:
        if var not in metadata_types:
            raise ValueError(f"Missing metadata type for variable '{var}'.")
        if var not in rep_na_counts:
            raise ValueError(f"Missing NA counts for variable '{var}'.")

        report_data.append({
            'variable': var,
            'parameter': name,
            'type': str(metadata_types[var]),
            'missing_in_n_replicates': rep_na_counts[var],
            'found_in_n_replicates': n_reps - rep_na_counts[var]
        })

    report_ext = os.path.splitext(output_path)[1].lower()
    if report_ext == '.json':
        with open(output_path, 'w') as outF:
            json.dump(report_data, outF, indent=4)
    elif report_ext == '.tsv':
        with open(output_path, 'w') as outF:
            outF.write('\t'.join(report_data[0].keys()))
            outF.write('\n')
            for row in report_data:
                outF.write('\t'.join(str(value) for value in row.values()))
                outF.write('\n')
    else:
        raise ValueError("Output report format must be either 'json' or 'tsv'.")


def _write_replicate_report(replicates, output_path='replicate_validation_report.json'):
    '''
    Write a report of replicate names and linked metadata.

    Parameters
    ----------
    replicates : dict
        Dictionary mapping replicate names to Replicate objects.
    output_path : str, optional
        Path to the output report file. Default is 'replicate_validation_report.json'.
    '''
    metadata_names = set(key for rep in replicates.values() for key in rep.metadata.keys())
    if 'ParameterBatch' in metadata_names:
        i = 1
        while f'ParameterBatch_{i}' in metadata_names:
            i += 1

        LOGGER.warning("'ParameterBatch' is a reserved column header name. "
                       f"It will be changed to 'ParameterBatch_{i}' in the report.")
        for rep in replicates.values():
            rep.metadata[f'ParameterBatch_{i}'] = rep.metadata.pop('ParameterBatch', rep.batch)

    report_data = []
    for rep_name, rep in replicates.items():
        report_data.append({'ParameterBatch': rep.batch, 'Replicate': rep_name})
        for key, value in rep.metadata.items():
            report_data[-1][key] = value

    report_ext = os.path.splitext(output_path)[1].lower()
    if report_ext == '.json':
        with open(output_path, 'w') as outF:
            json.dump(report_data, outF, indent=4)
    elif report_ext == '.tsv':
        with open(output_path, 'w') as outF:
            outF.write('\t'.join(report_data[0].keys()))
            outF.write('\n')
            for row in report_data:
                outF.write('\t'.join(str(value) for value in row.values()))
                outF.write('\n')
    else:
        raise ValueError("Output report format must be either 'json' or 'tsv'.")


def _log_warn_error(message, *args, warning=True):
    if warning:
        LOGGER.warning(message, *args)
    else:
        LOGGER.error(message, *args)


def validate_metadata(
    ms_files, metadata_df, metadata_types,
    color_vars=None, batch1=None, batch2=None,
    control_key=None, control_values=None,
    write_replicate_report=False, write_metadata_report=False,
    report_format='json', report_dir=None, strict=True
):
    '''
    Validate metadata against the provided parameters.

    Parameters
    ----------
    ms_files : dict
        Dictionary mapping batch names to lists of MS file names.
    metadata_df : pandas.DataFrame
        DataFrame containing metadata for the MS files.
    metadata_types : dict
        Dictionary mapping metadata variable names to their dtypes.
    color_vars : list of str, optional
        List of metadata variables to use for coloring PCA plots.
    batch1 : str, optional
        Metadata variable to use for batch 1.
    batch2 : str, optional
        Metadata variable to use for batch 2.
    control_key : str, optional
        Metadata key indicating control samples.
    control_values : list of str, optional
        List of metadata values indicating control samples.
    write_replicate_report : bool, optional
        Write a report of replicate annotations to a file.
    write_metadata_report : bool, optional
        Write a summary of metadata parameters to a file.
    report_format : str, optional
        Format of the report files. Either 'tsv', or 'json'. Default is 'json'.
    strict : bool, optional
        If True, any validation problems which will not cause the pipeline to fail,
        but are likely to result in unintended behavior, will be logged as errors.

    Returns
    -------
    bool
        True if metadata is valid, False otherwise.
    '''

    if (control_key is None) != (control_values is None):
        LOGGER.error("Both 'control_key' and 'color_vars' must be specified or both omitted.")
        return False

    ###### Validate metadata_df ######

    # Check that metadata parameters match in replicate metadata
    meta_params = []
    if color_vars is not None:
        for var in color_vars:
            meta_params.append((var, 'color_vars'))
    other_vars = {batch1: 'batch1', batch2: 'batch2', control_key: 'control_key'}
    for var, name in other_vars.items():
        if var is not None:
            meta_params.append((var, name))

    if metadata_df is None:
        if len(meta_params) > 0:
            LOGGER.error('Replicate metadata is None, but metadata parameters are specified.')
            for var, name in meta_params:
                if var is not None:
                    LOGGER.error("Without a metadata file, parameter %s = '%s' will cause an error.", name, var)
            return False

    else:
        if metadata_types is None:
            raise ValueError('metadata_types must be provided if metadata_df is not None.')

        all_meta_vars_good = True
        for var, name in meta_params:
            if var not in metadata_types:
                LOGGER.error("Metadata variable '%s' from '%s' parameter not found in metadata.", var, name)
                all_meta_vars_good = False
                continue
            if metadata_types[var] is Dtype.NULL:
                LOGGER.warning("All values for metadata variable '%s' from '%s' parameter are NULL.", var, name)

        if control_values:
            df_control_values = set(metadata_df[metadata_df['annotationKey'] == control_key]['annotationValue'].to_list())
            for param_value in control_values:
                if param_value not in df_control_values:
                    LOGGER.error(f"Control value '{param_value}' for key '{control_key}' not found in metadata.")
                    all_meta_vars_good = False

        if not all_meta_vars_good:
            return False
        LOGGER.info('All metadata parameters are present in replicate metadata.')

    ###### Validate replicates ######

    # make sure that there are not duplicated replicates
    rep_counts = Counter(rep for batch_files in ms_files.values() for rep in batch_files)
    if any(count > 1 for count in rep_counts.values()):
        for rep, count in rep_counts.items():
            if count > 1:
                LOGGER.error("Replicate '%s' is duplicated %d time%s in quant_spectra_dir.",
                             rep, count - 1, 's' if count > 2 else '')
        return False

    replicates = {}
    for batch_name, files in ms_files.items():
        for file_name in files:
            rep_name = os.path.splitext(file_name)[0]
            replicates[rep_name] = Replicate(file_name, batch=batch_name)

    metadata_reps = set(metadata_df['Replicate'].to_list())
    ms_file_reps = set(replicates.keys())

    # Warn if there are extra replicates in the metadata
    if len(metadata_reps - ms_file_reps) > 0:
        LOGGER.warning('There are %d replicates in the metadata that are not in quant_spectra_dir',
                       len(metadata_reps - ms_file_reps))

    # Error if there are missing replicates in the metadata
    bad_reps = ms_file_reps - metadata_reps
    for bad_rep in bad_reps:
        _log_warn_error(
            "Replicate '%s' in quant_spectra_dir is not present in the metadata.",
            bad_rep, warning=not strict
        )

    if len(bad_reps) > 0 and strict:
        return False

    # Add metadata to replicates
    meta_vars_unique = True
    for row in metadata_df.itertuples():
        if row.Replicate in replicates:
            if row.annotationKey in replicates[row.Replicate].metadata:
                _log_warn_error(
                    "Replicate '%s' already has metadata key '%s'. "
                    "Overwriting with value '%s'.",
                    row.Replicate, row.annotationKey, row.annotationValue,
                    warning=not strict
                )
                meta_vars_unique = False

            replicates[row.Replicate].metadata[row.annotationKey] = row.annotationValue

    if not meta_vars_unique and strict:
        LOGGER.error('Duplicate metadata keys found in replicates.')
        return False

    # Check that all replicates have the required metadata
    all_reps_good = True
    rep_na_counts = {var: 0 for var, _ in meta_params}
    for rep_name, rep in replicates.items():
        for var, name in meta_params:
            if var not in rep.metadata:
                _log_warn_error(
                    "Replicate '%s' is missing metadata key '%s' from '%s' parameter.",
                    rep_name, var, name, warning=not strict
                )
                all_reps_good = False
                continue

            value_empty = pd_isna(rep.metadata[var]) or rep.metadata[var] is None or rep.metadata[var] == ''
            if value_empty and metadata_types[var] is not Dtype.NULL:
                rep_na_counts[var] += 1

    if not all_reps_good and strict:
        return False

    # Log missing metadata values
    for var, na_count in rep_na_counts.items():
        if na_count > 0:
            LOGGER.warning("%s variable '%s' is missing in %d of %d replicates.",
                           metadata_types[var].name, var, na_count, len(replicates))
        else:
            LOGGER.info("%s variable '%s' is present in all %d replicates.",
                        metadata_types[var].name, var, len(replicates))

    # Write metadata report
    if write_metadata_report and len(meta_params) > 0:
        output_path = f'{report_dir + '/' if report_dir else ''}metadata_validation_report.{report_format}'
        _write_metadata_report(
            meta_params, metadata_types, rep_na_counts,
            len(replicates), output_path=output_path
        )
        LOGGER.info(f'Metadata report written to {output_path}')

    # Write replicate report
    if write_replicate_report:
        output_path = f'{report_dir + '/' if report_dir else ''}replicate_validation_report.{report_format}'
        _write_replicate_report(replicates, output_path=output_path)
        LOGGER.info(f'Replicate report written to {output_path}')

    return True


def parse_args(argv, prog=None):

    ###### Common subcommand arguments ######
    common_subcommand_args = argparse.ArgumentParser(add_help=False)
    panorama_args = common_subcommand_args.add_mutually_exclusive_group()
    panorama_args.add_argument(
        '-k', '--api-key', dest='api_key', default=None,
        help='API key to use for authentication.'
    )
    panorama_args.add_argument(
        '-n', '--nextflow-key', dest='nextflow_key', action='store_true', default=False,
        help="Get Panorama API key from Nextflow secrets manager. "
             "A Nextflow secret with the name 'PANORAMA_API_KEY' must be set up to use this option."
    )
    panorama_args.add_argument(
        '--panorama-public', dest='panorama_public', action='store_true',
        help='Use Panorama Public API key instead of user API key.'
    )
    common_subcommand_args.add_argument(
        '--report-format', dest='report_format', default=argparse.SUPPRESS, choices=['json', 'tsv'],
        help='Write metadata and MS file validation reports in specified format. '
             'By default no reports are written.'
    )
    common_subcommand_args.add_argument(
        '--metadata-output-path', dest='metadata_output_path', default=None,
        help='Path to write replicate metadata file if it is downloaded.'
    )

    ####### Top level parser #######
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)
    subparsers = parser.add_subparsers(required=True, dest='subcommand', description='Validation modes')

    ####### Config subcommand #######
    config_description = 'Validate nf-skyline-dia-ms pipeline parameters.'
    config_command_args = subparsers.add_parser(
        'config',
        parents=[common_subcommand_args],
        help=config_description,
        description=config_description
    )
    config_command_args.add_argument(
        'pipeline_config', nargs='+',
        help='Path to pipeline config file(s). If multiple config files are provided, '
             'they will be merged from left to right.'
    )

    schema_args = config_command_args.add_argument_group(
        'Options for specifying pipeline schema.'
        'By default, the pipeline config scheaa is downloaded from the latest '
        f'revision of "{DEFAULT_PIPELINE}"'
    )
    schema_args.add_argument(
        '--pipeline', default=argparse.SUPPRESS,
        help=f'Remote pipeline repository to use for validation. Default is "{DEFAULT_PIPELINE}".'
    )
    schema_args.add_argument(
        '-r', '--revision', default=argparse.SUPPRESS,
        help="Remote pipeline branch or tag to use for validation. "
             f"Default is '{DEFAULT_PIPELINE_REVISION}'"
    )
    schema_args.add_argument(
        '-s', '--schema', default=None, help='Path to local pipeline config schema file.'
    )

    ###### Params subcommand #######
    params_description = 'Validate parameters passed as arguments.'
    params_command_args = subparsers.add_parser(
        'params',
        parents=[common_subcommand_args],
        description=f"{params_description} This subcommand is intended for use in computational pipelines. "
                     "Human users should use the 'config' subcommand instead.",
        help=params_description
    )

    metadata_args = params_command_args.add_argument_group(
        'Metadata options', 'Options to manually specify metadata parameters.'
    )
    metadata_args.add_argument(
        '--bcMethod', choices=('limma', 'combat'), default='combat',
        help='Batch correction method. Default is "combat".'
    )
    metadata_args.add_argument(
        '--batch1', default=None, help='sampleMetadata variable to use for batch 1.'
    )
    metadata_args.add_argument(
        '--batch2', default=None, help='sampleMetadata variable to use for batch 2.'
    )
    metadata_args.add_argument(
        '--addCovariate', action='append', dest='covariate_vars',
        help='Add a sampleMetadata annotationKey to use as a covariate for batch correction.'
    )
    metadata_args.add_argument(
        '--addColor', action='append', dest='color_vars',
        help='Add a sampleMetadata annotationKey to use to color PCA plots.'
    )
    metadata_args.add_argument(
        '--controlKey', default=None, dest='control_key',
        help='sampleMetadata annotationKey that has variable indication whether a replicate is a control.'
    )
    metadata_args.add_argument(
        '--addControlValue', action='append', dest='control_values',
        help='Add sampleMetadata annotationValue(s) which indicate whether a replicate is a control.'
    )

    input_args = params_command_args.add_argument_group('Input files options')
    input_args.add_argument(
        '-m', '--metadata', help='Replicate metadata file. Can be a local file or Panorama URL.'
    )
    input_args.add_argument(
        '--chrom-lib-dir', dest='chrom_lib_spectra_dir', action='append', default=None,
        help='Add chromatogram library spectra directory.'
    )
    quant_dir_args = input_args.add_mutually_exclusive_group(required=True)
    quant_dir_args.add_argument(
        '-q', '--quant-file-dir', dest='quant_spectra_dir', action='append', default=None,
        help='Add quantative spectra directory.'
    )
    quant_dir_args.add_argument(
        '--quant-file-param', dest='quant_spectra_param',
        help='JSON file with quant_spectra_dir parameter'
    )
    quant_spectra_patterns = input_args.add_mutually_exclusive_group(required=True)
    quant_spectra_patterns.add_argument(
        '--quant-spectra-glob', dest='quant_spectra_glob',
        help='Glob pattern to match quant_spectra_dir files.'
    )
    quant_spectra_patterns.add_argument(
        '--quant-spectra-regex', dest='quant_spectra_regex',
        help='Regex pattern to match quant_spectra_dir files.'
    )
    chrom_spectra_patterns = input_args.add_mutually_exclusive_group(required=False)
    chrom_spectra_patterns.add_argument(
        '--chromatogram-library-spectra-glob', dest='chrom_lib_spectra_glob',
        help='Glob pattern to match chromatogram_library_spectra_dir files.'
    )
    chrom_spectra_patterns.add_argument(
        '--chromatogram-library-spectra-regex', dest='chrom_lib_spectra_regex',
        help='Regex pattern to match chromatogram_library_spectra_dir files.'
    )

    ###### Argument validation ######
    args = parser.parse_args(argv)
    if args.subcommand == 'config':
        if (hasattr(args, 'pipeline') and not hasattr(args, 'schema')) or \
           (hasattr(args, 'schema') and not hasattr(args, 'pipeline')):
            sys.stderr.write('--pipeline and --schema options conflict.\n')
            sys.exit(2)

        if not hasattr(args, 'pipeline'):
            args.pipeline = DEFAULT_PIPELINE

        revision = DEFAULT_PIPELINE_REVISION if not hasattr(args, 'revision') else args.revision
        args.schema = generate_git_url(
            args.pipeline, revision, filename='nextflow_schema.json'
        ) if not hasattr(args, 'schema') else args.schema

    if args.subcommand == 'params':
        have_pattern = args.chrom_lib_spectra_glob is not None or args.chrom_lib_spectra_regex is not None
        if args.chrom_lib_spectra_dir is not None and not have_pattern:
            sys.stderr.write('Chromatogram library spectra glob or regex pattern must be provided.\n')
            sys.exit(2)

    return args


def _get_config_path(data, vars, default=None):
    if len(vars) == 1:
        return data.get(vars[0], default)

    if vars[0] not in data:
        return default

    _node = data[vars[0]]
    return _get_config_path(_node, vars[1:], default=default)


def _main(argv, prog=None):
    args = parse_args(
        sys.argv[1:] if argv is None else argv,
        prog=prog
    )

    # Input files
    quant_spectra_dir = None
    quant_spectra_param = None
    quant_spectra_glob = None
    quant_spectra_regex = None
    chromatogram_library_spectra_param = None
    chromatogram_library_spectra_dir = None
    chromatogram_library_spectra_glob = None
    chromatogram_library_spectra_regex = None

    # Metadata settings
    color_vars = None
    batch1 = None
    batch2 = None
    control_key = None
    control_values = None
    metadata_path = None
    replicate_metadata = None

    # Determine API key
    api_key = None
    if args.nextflow_key:
        api_key = get_api_key_from_nextflow_secrets()
        if api_key is None:
            sys.exit(1)
    elif args.api_key:
        api_key = args.api_key
    elif args.panorama_public:
        api_key = PANORA_PUBLIC_KEY

    write_reports = False
    report_format = None
    if hasattr(args, 'report_format'):
        write_reports = True
        report_format = args.report_format

    if args.subcommand == 'config':
        if not os.path.exists(args.pipeline_config):
            LOGGER.error(f'Pipeline config file {args.pipeline_config} does not exist.')
            sys.exit(1)

        # Read and validate pipeline config file(s)
        LOGGER.info('Reading pipeline config files...')
        success, config_data = validate_config_files(args.config_paths, args.schema)
        if not success:
            LOGGER.error('Pipeline config validation failed.')
            sys.exit(1)
        LOGGER.info('Pipeline config validation succeeded.')

        color_vars = _get_config_path(config_data, ['qc_report', 'color_vars'])
        batch1 = _get_config_path(config_data, ['batch_report', 'batch1'])
        batch2 = _get_config_path(config_data, ['batch_report', 'batch2'])
        control_key = _get_config_path(config_data, ['qc_report', 'control_key'])
        control_values = _get_config_path(config_data, ['qc_report', 'control_values'])

        quant_spectra_param = config_data.get('quant_spectra_dir')
        chromatogram_library_spectra_param = config_data.get('chromatogram_library_spectra_dir')
        metadata_path = config_data.get('metadata_file')

        quant_spectra_regex = config_data.get('quant_spectra_regex')
        quant_spectra_glob = config_data.get('quant_spectra_glob')
        chromatogram_library_spectra_regex = config_data.get('chromatogram_library_spectra_regex')
        chromatogram_library_spectra_glob = config_data.get('chromatogram_library_spectra_glob')

    elif args.subcommand == 'params':
        color_vars = args.color_vars
        batch1 = args.batch1
        batch2 = args.batch2
        control_key = args.control_key
        control_values = args.control_values

        quant_spectra_glob = args.quant_spectra_glob
        quant_spectra_regex = args.quant_spectra_regex
        chromatogram_library_spectra_regex = args.chrom_lib_spectra_regex
        chromatogram_library_spectra_glob = args.chrom_lib_spectra_glob
        metadata_path = args.metadata

        if args.quant_spectra_dir is not None:
            quant_spectra_param = args.quant_spectra_dir
        elif args.quant_spectra_param is not None:
            quant_spectra_param = json.load(args.quant_spectra_dir)
            try:
                validate(quant_spectra_param, QUANT_SPECTRA_SCHEMA)
            except ValidationError as e:
                LOGGER.error(f'Quant spectra parameter validation failed: {e.message}')
                sys.exit(1)
        else:
            raise RuntimeError('quant_spectra_dir or quant_spectra_param must be provided.')

        if args.chrom_lib_spectra_dir is not None:
            chromatogram_library_spectra_param = args.chrom_lib_spectra_dir

    else:
        raise RuntimeError(f'Unknown subcommand: {args.subcommand}')

    quant_spectra_dir = parse_input_files(
        quant_spectra_param, api_key=api_key,
        file_glob=quant_spectra_glob, file_regex=quant_spectra_regex,
    )
    if chromatogram_library_spectra_param is not None:
        chromatogram_library_spectra_dir = parse_input_files(
            chromatogram_library_spectra_param, api_key=api_key,
            file_glob=chromatogram_library_spectra_glob,
            file_regex=chromatogram_library_spectra_regex,
        )

        if write_reports:
            chrom_reps = {os.path.splitext(file_name)[0]: Replicate(file_name)
                          for file_name in chromatogram_library_spectra_dir}
            _write_replicate_report(
                chrom_reps, output_path=f'chromatogram_library_replicate_report.{report_format}'
            )

    # read metadata
    if metadata_path:
        download_metadata = args.metadata_output_path is not None
        metadata = get_file(
            metadata_path, log_name='replicate metadata', api_key=api_key,
            dest_path=args.metadata_output_path, return_text=not download_metadata
        )
        if metadata is None:
            sys.exit(1)

        metadata_reader = Metadata()
        if download_metadata:
            with open(args.metadata_output_path, 'r') as inF:
                if not metadata_reader.read(inF):
                    sys.exit(1)
        else:
            metadata_stream = StringIO(metadata)
            metadata_stream.seek(0)
            if not metadata_reader.read(metadata_stream):
                sys.exit(1)

        replicate_metadata = metadata_reader.df

    LOGGER.info('Validating metadata...')
    success = validate_metadata(
        quant_spectra_dir, replicate_metadata,
        color_vars=color_vars, batch1=batch1, batch2=batch2,
        control_key=control_key, control_values=control_values,
        write_replicate_report=write_reports, write_metadata_report=write_reports,
        report_format=report_format
    )
    if not success:
        LOGGER.error('Metadata validation failed.')
        sys.exit(1)
    LOGGER.info('Metadata validation succeeded.')


def main():
    LOGGER.warning("Calling this script directly is deprecated. Use 'dia_qc validate' instead.")
    _main(sys.argv[1:])


if __name__ == '__main__':
    main()