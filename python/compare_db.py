
import sys
import argparse
import sqlite3
import os

import pandas as pd

from pyDIAUtils.dia_db_utils import check_schema_version

def open_db(path):
    if os.path.isfile(path):
        conn = sqlite3.connect(path)
    else:
        sys.stderr.write(f'Database file (path) does not exist!\n')
        return None

    # check database version
    if not check_schema_version(conn):
        return None

    return conn


def compare_metadata(lhs, rhs):
    pass


class DtypeDiff():
    def __init__(self, name):
        self.name = name
        self.n_eq = 0
        self.n_neq = 0

    def is_numeric(self):
        return False

    def add(self, lhs, rhs):
        if lhs == rhs:
            self.n_eq += 1
        else:
            self.n_neq += 1

    def __str__(self):
        return f'{self.name}: object, n_eq={self.n_eq}, n_neq={self.n_neq}'

    def to_dict(self):
        return {'name': self.name, 'n_eq': self.n_eq, 'n_neq': self.n_neq}


class NumericDiff(DtypeDiff):
    def __init__(self, name):
        DtypeDiff.__init__(self, name)
        self.max_diff = 0

    def is_numeric(self):
        return True

    def add(self, lhs, rhs):
        DtypeDiff.add(self, lhs, rhs)
        self.max_diff = max((self.max_diff, abs(lhs) - abs(rhs)))

    def __str__(self):
        return f'{self.name}: numeric, n_eq={self.n_eq}, n_neq={self.n_neq}, max_diff={self.max_diff}'

    def to_dict(self):
        return DtypeDiff.to_dict(self) | {'max_diff': self.max_diff}


def get_diff_type(type_str, name):
    if type_str.startswith('int') or type_str.startswith('float'):
        return NumericDiff(name)
    return DtypeDiff(name)


def write_diff_table(diff_counts):

    table_data = [['name', 'numeric', 'n_eq', 'n_neq']]

    any_numeric = False
    for diff in diff_counts:
        if diff.is_numeric():
            any_numeric = True
            break

    if any_numeric:
        table_data[0].append('max_diff')


def compare_replicates(lhs, rhs):

    def read_table(conn):
        reps = pd.read_sql('SELECT * FROM replicates', conn)
        reps['replicateKey'] = reps['replicate'] + ';' + reps['project']
        ids = {row.replicateKey for row in reps.itertuples()}

        return reps, ids

    lhs_reps, lhs_ids = read_table(lhs)
    rhs_reps, rhs_ids = read_table(rhs)

    overlap_ids = rhs_ids & lhs_ids

    sys.stdout.write('Checking replicate tables...\n')

    sys.stdout.write(f'\t{len(overlap_ids)} replicates common to both.\n')
    if rhs_ids == lhs_ids:
        sys.stdout.write('\tReplicate keys are identical\n')

    else:
        sys.stdout.write(f'\t{len(lhs_ids - rhs_ids)} unique to lhs\n')
        sys.stdout.write(f'\t{len(rhs_ids - lhs_ids)} unique to rhs\n')

    # for rep_key in (rhs_ids & lhs_ids):
    #     print(rep_key)

    rhs_ids = {row.replicateId: row.replicateKey for row in rhs_reps.itertuples() if row.replicateKey in overlap_ids}
    lhs_ids = {row.replicateId: row.replicateKey for row in lhs_reps.itertuples() if row.replicateKey in overlap_ids}

    # rhs_reps

    types = {k: str(v) for k, v in rhs_reps.dtypes.items() if k != 'replicateKey'}

    rep_data = {k: {} for k in types}
    for row in lhs_reps.itertuples():
        if row.replicateKey in overlap_ids:
            for col in types:
                rep_data[col][row.replicateKey] = [getattr(row, col)]

    for row in rhs_reps.itertuples():
        if row.replicateKey in overlap_ids:
            for col in types:
                rep_data[col][row.replicateKey].append(getattr(row, col))

    # df_j = lhs_reps.join(rhs_reps, on='replicateKey', how='inner')
    # df_j = lhs_reps.join(rhs_reps, on='replicateKey', how='inner', lsuffix='lhs_', rsuffix='rhs_')

    diff_counts = {col: get_diff_type(col_type, col) for col, col_type in types.items()}
    for col in types:
        for pair in rep_data[col].values():
            diff_counts[col].add(*pair)

    # for key, value in diff_counts.items():
    #     sys.stdout.write(f'\t{str(value)}\n')

    

    return lhs_ids, rhs_ids


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-l', '--level', default='replicate',
                        help='Level of comparison. Default is replicate.')

    parser.add_argument('lhs')
    parser.add_argument('rhs')

    args = parser.parse_args()

    lhs = open_db(args.lhs)
    rhs = open_db(args.rhs)

    if lhs is None or rhs is None:
        sys.exit(1)

    compare_replicates(lhs, rhs)


if __name__ == '__main__':
    main()
