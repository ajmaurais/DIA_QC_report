
import os
import subprocess
import re
from shlex import join as join_shell
import inspect
from abc import ABC, abstractmethod
import io
import contextlib
import traceback
import logging

from DIA_QC_report import parse_data
from DIA_QC_report import normalize_db
from pandas import testing as pd_testing
from numpy import isnan

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def function_kwargs(f):
    '''
    Yield pairs of argument names and default values for all the named kwargs in function `f`
    '''
    signature = inspect.signature(f)
    for param in signature.parameters.values():
        if param.default is not inspect.Parameter.empty:
            yield param.name, param.default


class AbstractTestsBase(ABC):
    '''
    Base class for abstract shared unittest base classes.
    '''
    @classmethod
    @abstractmethod
    def setUpClass(cls):
        pass


    @classmethod
    @abstractmethod
    def tearDownClass(cls):
        pass


    @abstractmethod
    def assertIsNotNone(self, expr):
        pass


    @abstractmethod
    def assertTrue(self, expr):
        pass


    @abstractmethod
    def assertEqual(self, lhs, rhs):
        pass


    @abstractmethod
    def assertFalse(self, expr):
        pass


    @abstractmethod
    def assertRaises(self, exception):
        pass


    @abstractmethod
    def assertLogs(self, logger, level=None):
        pass


    @abstractmethod
    def assertAlmostEqual(self, lhs, rhs, places=None, delta=None):
        pass


    @abstractmethod
    def assertDictEqual(self, lhs, rhs):
        pass


    @abstractmethod
    def assertIsInstance(self, obj, obj_type):
        pass


    @abstractmethod
    def assertNoLogs(self, logger, level=None):
        pass


    @abstractmethod
    def assertGreater(self, lhs, rhs):
        pass


    @abstractmethod
    def fail(self, msg):
        pass


    def assertDataFrameEqual(self, a, b, **kwargs):
        try:
            pd_testing.assert_frame_equal(a, b, **kwargs)
        except AssertionError as e:
            self.fail(str(e))


    def assertSeriesEqual(self, a, b, **kwargs):
        try:
            pd_testing.assert_series_equal(a, b, **kwargs)
        except AssertionError as e:
            self.fail(str(e))


    def assertInLog(self, message, cm):
        self.assertTrue(any(message in m for m in cm.output),
                        f"Expected log message '{message}' not found in {cm.output}")


    def assertNotInLog(self, message, cm):
        self.assertFalse(any(message in m for m in cm.output),
                         f"Expected log message '{message}' not found in {cm.output}")


    def assertDataDictEqual(self, lhs, rhs, places=6, col_deltas=None):
        '''
        Check that 2 dictionaries having the format {'rep_name': {<var_name>: <value>, ...}, ...}
        are almost equal. <value>(s) with type `float` will be compared using the `places` or
        `col_deltas` arguments using unittest.assertAlmostEqual. All other types will be tested
        for equality.

        Parameters
        ----------
        lhs: dict
        rhs: dict
        places: int
            Round to given number of decimal places.
        col_deltas: dict
            Dictionary of delta(s) to pass to unittest.assertAlmostEqual for each <var_name>
        '''
        self.assertIsInstance(lhs, dict)
        self.assertIsInstance(rhs, dict)

        lhs_keys = set(lhs.keys())
        rhs_keys = set(rhs.keys())

        self.assertEqual(lhs_keys, rhs_keys)

        for rep in lhs_keys:
            lhs_vars = set(lhs[rep].keys())
            rhs_vars = set(rhs[rep].keys())

            self.assertEqual(lhs_vars, rhs_vars)
            for var in lhs_vars:
                if type(lhs[rep][var]) != type(rhs[rep][var]):
                    raise AssertionError(f"Data types differ in column '{var}'")
                if isinstance(lhs[rep][var], float):
                    if isnan(lhs[rep][var]) and isnan(rhs[rep][var]):
                        continue
                    try:
                        if col_deltas is None:
                            self.assertAlmostEqual(lhs[rep][var], rhs[rep][var], places=places)
                        else:
                            delta = col_deltas.get(var, None)
                            self.assertAlmostEqual(lhs[rep][var], rhs[rep][var], delta=delta)
                    except AssertionError as e:
                        raise AssertionError(f"In column '{var}': {str(e)}") from e
                else:
                    if lhs[rep][var] != rhs[rep][var]:
                        raise AssertionError("Values differ in column '{var}'. {lhs[rep][var]} != {rhs[rep][var]}")


def remove_test_dir(test_dir, recursive=False):
    """Remove test directory if it exists. If recursive=True, delete all nested files/dirs."""
    if not os.path.isdir(test_dir):
        return

    if recursive:
        # walk bottom‐up so files and subdirs are removed before their parents
        for root, dirs, files in os.walk(test_dir, topdown=False):
            for fname in files:
                os.remove(os.path.join(root, fname))
            for dname in dirs:
                os.rmdir(os.path.join(root, dname))
        os.rmdir(test_dir)
    else:
        # only delete files in the top‐level directory
        for entry in os.listdir(test_dir):
            path = os.path.join(test_dir, entry)
            if os.path.isfile(path):
                os.remove(path)
        os.rmdir(test_dir)


def make_work_dir(work_dir, clear_dir=False):
    '''
    Setup work directory for test.

    Parameters
    ----------
    clear_dir: bool
        If the directory already exists, should the files already in directory be deleted?
        Will not work recursively or delete directories.
    '''
    if not os.path.isdir(work_dir):
        if os.path.isfile(work_dir):
            raise RuntimeError('Cannot create work directory!')
        os.makedirs(work_dir)
    else:
        if clear_dir:
            for file in os.listdir(work_dir):
                os.remove(f'{work_dir}/{file}')


def run_command(command, wd, prefix=None):
    '''
    Run command in subprocess and write stdout, stderr, return code and command to
    textfiles in specified directory.

    Parameters
    ----------
    command: list
        The command to run. Each argument should be a separate list element.
    wd: str
        The directory to run the command from.
    prefix: str
        A prefix to add to stdout, stderr, rc and command files.
        If None, the name of the calling function is used as the prefix.
    '''
    result = subprocess.run(command, cwd=wd,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            shell=False, check=False)

    prefix_path = f'{wd}/{prefix if prefix else inspect.stack()[1][3]}'

    result.stderr = result.stderr.decode('utf-8')
    result.stdout = result.stdout.decode('utf-8')

    with open(f'{prefix_path}.command.txt', 'w') as outF:
        outF.write(f"{join_shell(command)}\n")
    with open(f'{prefix_path}.stdout.txt', 'w') as outF:
        outF.write(f"{result.stdout}\n")
    with open(f'{prefix_path}.stderr.txt', 'w') as outF:
        outF.write(f"{result.stderr}\n")
    with open(f'{prefix_path}.rc.txt', 'w') as outF:
        outF.write(f'{str(result.returncode)}\n')

    return result


class ProcessResult:
    def __init__(self, stdout, stderr, returncode):
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode


    def __str__(self):
        return f"ProcessResult(stdout={self.stdout!r}, stderr={self.stderr!r}, returncode={self.returncode})"


def run_main(main_fxn, argv, wd, prog=None, prefix=None):
    '''
    Run a module main function and capture its output.

    Parameters
    ----------
    main_fxn: callable
        The main function to run. Should accept `argv` and `prog` as arguments.
    argv: list
        The command line arguments to pass to the main function.
    wd: str
        The working directory to run the main function in.
    prog: str, optional
        The program name to use in the output files. If None, the name of the calling function is used.
    prefix: str, optional
        A prefix to add to stdout, stderr, rc and command files.

    Returns
    -------
    ProcessResult
        An object containing the captured stdout, stderr, and return code.
    '''

    # ensure working dir exists
    if not os.path.isdir(wd):
        os.makedirs(wd)
    prefix_path = f'{wd}/{prefix if prefix else inspect.stack()[1][3]}'

    # capture
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
        for handler in logging.root.handlers:
            if hasattr(handler, 'stream'):
                handler.stream = stderr_buf
        try:
            curwd = os.getcwd()
            os.chdir(wd)
            main_fxn(argv, prog=prog)
            rc = 0
        except SystemExit as e:
            rc = e.code if isinstance(e.code, int) else 0
        except Exception:
            traceback.print_exc(file=stderr_buf)
            rc = 1
            os.chdir(curwd)

    out_str = stdout_buf.getvalue()
    err_str = stderr_buf.getvalue()

    # write files
    command = f'{main_fxn.__module__}.{main_fxn.__name__}' if prog is None else re.split(r'\s+', prog)
    with open(f"{prefix_path}.command.txt", "w") as outF:
        outF.write(f"{join_shell(command + argv)}\n")
    with open(f"{prefix_path}.stdout.txt", "w") as outF:
        outF.write(out_str)
    with open(f"{prefix_path}.stderr.txt", "w") as outF:
        outF.write(err_str)
    with open(f"{prefix_path}.rc.txt", "w") as outF:
        outF.write(f"{rc}\n")

    return ProcessResult(out_str, err_str, rc)


def setup_single_db(data_dir, output_dir, project,
                    metadata_suffix='_metadata.tsv', output_prefix=None, subprocess=False,
                    overwrite_mode='error', group_by_gene=False, clear_dir=False):
    make_work_dir(output_dir, clear_dir)
    grouping = 'by_gene' if group_by_gene else 'by_protein'

    command = [f'--projectName={project}',
               f'--overwriteMode={overwrite_mode}',
               '-m', f'{data_dir}/metadata/{project}{metadata_suffix}',
               f'{data_dir}/skyline_reports/{project}_replicate_quality.tsv',
               f'{data_dir}/skyline_reports/{project}_{grouping}_precursor_quality.tsv']

    if group_by_gene:
        command.insert(2, '--groupBy=gene')

    if subprocess:
        command = ['dia_qc', 'parse'] + command
        return run_command(command, output_dir, prefix=output_prefix)

    return run_main(parse_data._main, command, output_dir, prog='dia_qc parse', prefix=output_prefix)


def setup_multi_db(data_dir, output_dir,
                   group_by_gene=False, clear_dir=False,
                   normalize=False, subprocess=False,
                   metadata_suffix='_metadata.tsv'):
    make_work_dir(output_dir, clear_dir)
    grouping = 'by_gene' if group_by_gene else 'by_protein'

    db_path = f'{output_dir}/data.db3'

    commands = [['--projectName=Sp3',
                 '-m', f'{data_dir}/metadata/Sp3{metadata_suffix}',
                 f'{data_dir}/skyline_reports/Sp3_replicate_quality.tsv',
                 f'{data_dir}/skyline_reports/Sp3_{grouping}_precursor_quality.tsv'],
                ['--overwriteMode=append', '--projectName=Strap',
                 '-m', f'{data_dir}/metadata/Strap{metadata_suffix}',
                 f'{data_dir}/skyline_reports/Strap_replicate_quality.tsv',
                 f'{data_dir}/skyline_reports/Strap_{grouping}_precursor_quality.tsv']]

    if os.path.isfile(db_path):
        commands[0].insert(0, '--overwriteMode=overwrite')

    if group_by_gene:
        for i in range(len(commands)):
            commands[i].insert(0, '--groupBy=gene')

    results = list()
    for i, command in enumerate(commands):
        if subprocess:
            command = ['dia_qc', 'parse'] + command
            results.append(run_command(command, output_dir, prefix=f'add_project_{i}'))

        else:
            results.append(run_main(
                parse_data._main, command, output_dir,
                prefix=f'add_project_{i}', prog='dia_qc parse'
            ))

    if normalize:
        normalize_command = ['-m=median', db_path]
        if subprocess:
            results.append(run_command(
                ['dia_qc', 'normalize'] + normalize_command, output_dir, prefix='normalize_db'
            ))
        else:
            results.append(run_main(
                normalize_db._main, normalize_command, output_dir,
                prefix='normalize_db', prog='dia_qc normalize'
            ))

    return results