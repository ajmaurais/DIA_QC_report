
import os
import subprocess
import inspect
from abc import ABC, abstractmethod
from numpy import isnan

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


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


    def assertDataDictEqual(self, lhs, rhs, places=6, col_deltas=None):
        '''
        Check that 2 dictionaries having format
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


def make_work_dir(work_dir, clear_dir=False):
    if not os.path.isdir(work_dir):
        if os.path.isfile(work_dir):
            raise RuntimeError('Cannot create work directory!')
        os.makedirs(work_dir)
    else:
        if clear_dir:
            for file in os.listdir(work_dir):
                os.remove(f'{work_dir}/{file}')


def run_command(command, wd, prefix=None):
    result = subprocess.run(command, cwd=wd,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            shell=False, check=False)

    prefix_path = f'{wd}/{prefix if prefix else inspect.stack()[1][3]}'

    result.stderr = result.stderr.decode('utf-8')
    result.stdout = result.stdout.decode('utf-8')

    with open(f'{prefix_path}.command.txt', 'w') as outF:
        outF.write(f"{' '.join(command)}\n")
    with open(f'{prefix_path}.stdout.txt', 'w') as outF:
        outF.write(f"{result.stdout}\n")
    with open(f'{prefix_path}.stderr.txt', 'w') as outF:
        outF.write(f"{result.stderr}\n")
    with open(f'{prefix_path}.rc.txt', 'w') as outF:
        outF.write(f'{str(result.returncode)}\n')

    return result


def setup_single_db(data_dir, output_dir, project,
                    metadata_suffix='_metadata.tsv', output_prefix=None,
                    overwrite_mode='error', group_by_gene=False, clear_dir=False):
    make_work_dir(output_dir, clear_dir)
    grouping = 'by_gene' if group_by_gene else 'by_protein'

    command = ['parse_data', f'--projectName={project}',
               f'--overwriteMode={overwrite_mode}',
               '-m', f'{data_dir}/metadata/{project}{metadata_suffix}',
               f'{data_dir}/skyline_reports/{project}_replicate_quality.tsv',
               f'{data_dir}/skyline_reports/{project}_{grouping}_precursor_quality.tsv']

    if group_by_gene:
        command.insert(1, '--groupBy=gene')

    return run_command(command, output_dir, prefix=output_prefix)


def setup_multi_db(data_dir, output_dir, group_by_gene=False, clear_dir=False):
    make_work_dir(output_dir, clear_dir)
    grouping = 'by_gene' if group_by_gene else 'by_protein'

    commands = [['parse_data', '--projectName=Sp3',
                 '-m', f'{data_dir}/metadata/Sp3_metadata.tsv',
                 f'{data_dir}/skyline_reports/Sp3_replicate_quality.tsv',
                 f'{data_dir}/skyline_reports/Sp3_{grouping}_precursor_quality.tsv'],
                ['parse_data', '--overwriteMode=append', '--projectName=Strap',
                 '-m', f'{data_dir}/metadata/Strap_metadata.tsv',
                 f'{data_dir}/skyline_reports/Strap_replicate_quality.tsv',
                 f'{data_dir}/skyline_reports/Strap_{grouping}_precursor_quality.tsv']]

    if os.path.isfile(f'{output_dir}/data.db3'):
        commands[0].insert(1, '--overwriteMode=overwrite')

    if group_by_gene:
        for i in range(len(commands)):
            commands[i].insert(1, '--groupBy=gene')

    results = list()
    for command in commands:
        results.append(run_command(command, output_dir))

    return results

