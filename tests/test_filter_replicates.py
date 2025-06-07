
import unittest
import os
import sqlite3
import random

import setup_functions

import DIA_QC_report.submodules.dia_db_utils as db_utils
import DIA_QC_report.filter_replicates as filter_replicates


class TestFilterReplicates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work_dir = f'{setup_functions.TEST_DIR}/work/test_filter_replicates/'
        cls.base_db_path = f'{cls.work_dir}/data.db3'
        cls.data_dir = f'{setup_functions.TEST_DIR}/data/'

        # remove tables subdirectory in work_dir if necissary
        for d in (f'{cls.work_dir}/norm_tables', f'{cls.work_dir}/no_norm_tables'):
            if os.path.isdir(d):
                for file in os.listdir(d):
                    os.remove(f'{d}/{file}')
                os.rmdir(d)

        cls.parse_results = setup_functions.setup_multi_db(cls.data_dir,
                                                           cls.work_dir,
                                                           clear_dir=True,
                                                           group_by_gene=False)

        conn = None
        if all(cmd.returncode == 0 for cmd in cls.parse_results):
            if os.path.isfile(cls.base_db_path):
                conn = sqlite3.connect(cls.base_db_path)

                cur = conn.cursor()
                cur.execute('SELECT replicate FROM replicates;')
                cls.replicates = [x[0] for x in cur.fetchall()]
                cur.execute('SELECT DISTINCT project FROM replicates;')
                cls.projects = [x[0] for x in cur.fetchall()]
                conn.close()


    def setUp(self):
        if not os.path.isfile(self.base_db_path):
            raise RuntimeError('Database does not exist!')

        # copy test db to unique db for test
        self.db_path = f'{self.work_dir}/{self._testMethodName}.db3'
        conn = sqlite3.connect(self.base_db_path)
        self.conn = sqlite3.connect(self.db_path)
        conn.backup(self.conn)
        conn.close()


    def tearDown(self):
        if self.conn is not None:
            self.conn.close()


    def do_command_log_test(self, command):
        self.assertIsNotNone(self.conn)

        last_command = db_utils.get_last_command(self.conn)
        last_command[0] = os.path.basename(last_command[0])

        self.assertEqual(last_command, command)


    def test_conflicting_rep_options_fail(self):
        self.assertIsNotNone(self.conn)
        command = ['dia_qc', 'filter',
                   '--excludeRep', self.replicates[0],
                   '--includeAll', self.db_path]

        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 1)
        self.assertTrue(filter_replicates.CONFLICTING_ARG_MESSAGE in result.stderr)

        last_command = db_utils.get_last_command(self.conn)
        last_command[0] = os.path.basename(last_command[0])
        self.assertNotEqual(last_command, command)


    def test_conflicting_project_options_fail(self):
        self.assertIsNotNone(self.conn)
        command = ['dia_qc', 'filter',
                   '--excludeProject', self.projects[0],
                   '--includeAll', self.db_path]

        result = setup_functions.run_command(command, self.work_dir)

        self.assertEqual(result.returncode, 1)
        self.assertTrue(filter_replicates.CONFLICTING_ARG_MESSAGE in result.stderr)

        last_command = db_utils.get_last_command(self.conn)
        last_command[0] = os.path.basename(last_command[0])
        self.assertNotEqual(last_command, command)


    def test_include_all(self):
        ''' Test --includeAll option. '''
        self.assertIsNotNone(self.conn)
        command = ['dia_qc', 'filter',
                   '--excludeProject', self.projects[0],
                   self.db_path]

        result = setup_functions.run_command(command, self.work_dir,
                                             prefix='test_include_all_exclude')
        self.assertEqual(result.returncode, 0)
        self.do_command_log_test(command)

        cur = self.conn.cursor()
        cur.execute('SELECT replicate FROM replicates WHERE includeRep == 0;')
        self.assertGreater(len(cur.fetchall()), 0)

        command = ['dia_qc', 'filter', '--includeAll', self.db_path]

        result = setup_functions.run_command(command, self.work_dir)
        self.assertEqual(result.returncode, 0)
        self.do_command_log_test(command)
        cur.execute('SELECT replicate FROM replicates WHERE includeRep == 0;')
        self.assertEqual(len(cur.fetchall()), 0)


    def test_exclude_rep_option(self):
        ''' Test --excludeRep option. '''
        self.assertIsNotNone(self.conn)
        n_itterations = 5
        random.seed(1)

        max_subset_len = int(len(self.replicates) * 0.75)
        cur = self.conn.cursor()

        for i in range(n_itterations):
            skip_reps = random.sample(self.replicates, random.randint(1, max_subset_len))
            command = ['dia_qc', 'filter'] + [f'-x={rep}' for rep in skip_reps]
            command.append(self.db_path)

            result = setup_functions.run_command(command, self.work_dir,
                                                 prefix=f'filter_replicates_{i + 1}')

            self.assertEqual(result.returncode, 0)
            self.do_command_log_test(command)
            self.assertTrue(f'Excluding {len(skip_reps)} replicates.' in result.stdout)

            cur.execute('SELECT replicate, includeRep FROM replicates;')
            db_skip_reps = cur.fetchall()
            self.assertEqual(len([include for _, include in db_skip_reps if include == 0]),
                             len(skip_reps))

            for rep, included in db_skip_reps:
                if included == 1:
                    self.assertFalse(rep in skip_reps)
                else:
                    self.assertEqual(included, 0)
                    self.assertTrue(rep in skip_reps)

            db_utils.mark_all_reps_included(self.conn, quiet=True)


    def test_exclude_project_option(self):
        ''' Test --excludeProject option. '''
        self.assertIsNotNone(self.conn)

        cur = self.conn.cursor()
        cur.execute('SELECT DISTINCT project FROM replicates;')
        projects = [x[0] for x in cur.fetchall()]

        for i, project in enumerate(projects):
            cur.execute('SELECT replicate FROM replicates WHERE project == ?;', (project,))
            project_replicates = [x[0] for x in cur.fetchall()]

            command = ['dia_qc', 'filter', f'-p={project}', self.db_path]
            result = setup_functions.run_command(command, self.work_dir,
                                                 prefix=f'filter_projects_{i + 1}')

            self.assertEqual(result.returncode, 0)
            self.do_command_log_test(command)
            self.assertTrue(f'Excluding {len(project_replicates)} replicates.' in result.stdout)

            cur.execute('SELECT replicate, includeRep FROM replicates;')
            db_skip_reps = cur.fetchall()
            self.assertEqual(len([include for _, include in db_skip_reps if include == 0]),
                             len(project_replicates))

            for rep, included in db_skip_reps:
                if included == 1:
                    self.assertFalse(rep in project_replicates)
                else:
                    self.assertEqual(included, 0)
                    self.assertTrue(rep in project_replicates)

            db_utils.mark_all_reps_included(self.conn, quiet=True)
