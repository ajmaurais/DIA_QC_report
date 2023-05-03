
import re
import warnings
import sqlite3
from statistics import stdev
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def pc_matrix(df):
    df_s = StandardScaler().fit_transform(df.transpose())

    pca = PCA()
    pca.fit(df_s)

    df_s = pd.DataFrame(df_s.transpose(), index=df.index, columns=df.columns)
    pc = np.matmul(df_s.transpose(), pca.components_.transpose())
    pc_var = (pca.singular_values_ ** 2) / sum(pca.singular_values_ ** 2) * 100

    return pc, pc_var


def pca_plot(df, label_col, pc_var, label_type='discrete',
             fname=None, dpi=250, x_axis_pc=0, y_axis_pc=1):

    cmap = plt.get_cmap('viridis')
    colors = {label: color for color, label in zip(cmap(np.linspace(0, 1, len(pc[label_col].drop_duplicates()))),
                                                   pc[label_col].drop_duplicates())}
    pc['color'] = pc[label_col].apply(lambda x: colors[x])

    if label_type == 'discrete':
        fig = plt.figure(figsize = (6, 4), dpi=dpi)
        ax = fig.add_axes([0.1, 0.15, 0.65, 0.75])

        for label in sorted(pc[label_col].drop_duplicates()):
            sele = pc[label_col] == label
            ax.scatter(pc[sele][x_axis_pc], pc[sele][y_axis_pc], color=colors[label], label=label)
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), title=label_col, alignment='left')

    elif label_type == 'continuous':
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        points=ax.scatter(pc[x_axis_pc], pc[y_axis_pc], c=pc[label_col], cmap='viridis')
        fig.colorbar(points, label=label_col, use_gridspec=False,
                     ticks=MaxNLocator(integer=True) if pd.api.types.is_integer_dtype(pc[label_col]) else None)

    # axis labels
    ax.set_xlabel(f'PC {x_axis_pc + 1} {pc_var[x_axis_pc]:.1f}% var')
    ax.set_ylabel(f'PC {y_axis_pc + 1} {pc_var[y_axis_pc]:.1f}% var')

    if fname is None:
        plt.show()
    else:
        plt.savefig(fname)
    plt.close()


def convert_string_cols(df):
    float_re = re.compile(r'^[+-]?[0-9]+\.[0-9]?$')
    int_re = re.compile(r'^[+-]?[0-9]+')

    def rank_types(x):
        if int_re.search(x):
            return 1
        if float_re.search(x):
            return 2
        return 3

    ret = df
    if len(ret.index) > 0:
        for col in df.columns:
            if isinstance(df[col][0], str):
                type_rank = min(df[col].apply(rank_types))
                if type_rank == 1:
                    df[col] = df[col].apply(int)
                elif type_rank == 2:
                    df[col] = df[col].apply(float)
    return ret


QUERY = '''
SELECT
    p.replicateId,
    r.acquiredRank,
    p.modifiedSequence,
    p.precursorCharge,
    p.totalAreaFragment
FROM precursors p
LEFT JOIN replicates r
ON p.replicateId = r.replicateId; '''

# conn = sqlite3.connect("/home/ajm/code/DIA_QC_report/testData/csf_data.db3")
# conn = sqlite3.connect("/home/ajm/code/DIA_QC_report/testData/full_data.db3")
conn = sqlite3.connect("/home/ajm/code/DIA_QC_report/testData/data.db3")

df_test = pd.read_sql(QUERY, conn)

# set zero areas to the mininum non-zero value
df_test.loc[df_test['totalAreaFragment'] == 0, 'totalAreaFragment'] = min(df_test[df_test['totalAreaFragment'] > 0]['totalAreaFragment'])

df_test['log2TotalAreaFragment'] = np.log10(df_test['totalAreaFragment'])
df_test['zScore'] = df_test.groupby(["acquiredRank"])['log2TotalAreaFragment'].apply(lambda x: (x - np.mean(x)) / stdev(x))

df_wide = df_test.pivot_table(index=['modifiedSequence', 'precursorCharge'],
                         columns="replicateId", values='zScore')

# actually do pc analysis
pc, pc_var = pc_matrix(df_wide)

meta_values = {'gender': 'discrete', 'vital_status': 'discrete',
               'tumor_grade': 'discrete', 'year_of_birth': 'continuous'}

# check if metadata table exists in db
cur = conn.cursor()
cur.execute('SELECT name FROM sqlite_master WHERE type="table" AND name="metadata"')
table = cur.fetchall()

acquired_ranks = df_test[["replicateId", "acquiredRank"]].drop_duplicates()
acquired_ranks = pd.DataFrame(acquired_ranks['acquiredRank'], index=acquired_ranks['replicateId'])

# add metadata to 
if len(table) > 0:
    METADATA_QUERY = '''
    SELECT
        replicateId, annotationKey as columns, annotationValue
    FROM metadata
    WHERE annotationKey IN ("%s")''' % '", "'.join(meta_values.keys())

    # get metadata labels for pca plot
    metadata_df = pd.read_sql(METADATA_QUERY, conn)
    metadata = metadata_df.pivot(index="replicateId", columns="columns", values="annotationValue")
    metadata = acquired_ranks.join(metadata)
    meta_values['acquiredRank'] = 'continuous'
    metadata = convert_string_cols(metadata)

    # join meatadata to pc matrix
    pc = pc.join(metadata)
    for label_name, label_type in meta_values.items():
        if label_type == 'discrete':
            if any(pc[label_name].apply(pd.isna)):
                warnings.warn('Missing label values!', Warning)
            pc[label_name] = pc[label_name].apply(str)
        elif label_type == 'continuous':
            if any(pc[label_name].apply(pd.isna)):
                raise RuntimeError('Cannot have missing label values in continious scale!')
        else:
            raise RuntimeError(f'"{label_type}" is an unknown label_type!')

else:
    meta_values = {'acquiredRank': 'continuous'}
    pc = pc.join(acquired_ranks)

conn.close()


# labels = metadata
# label_names = meta_values
# x_axis_pc = 0
# y_axis_pc = 1
# dpi=250
# df = df_wide
# fname='fig/pca_test.pdf'

for label_name, label_type in meta_values.items():
    ofname = f'fig/pc_{label_name}.pdf'
    pca_plot(pc, label_name, pc_var, label_type=label_type, fname=ofname)

