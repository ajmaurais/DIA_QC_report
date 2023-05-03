
import sqlite3
from statistics import stdev
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def pca_plot(df, fname=None, x_axis_pc=0, y_axis_pc=1):

    df_s = StandardScaler().fit_transform(df.transpose())

    pca = PCA()
    pca.fit(df_s)

    df_s = pd.DataFrame(df_s.transpose(), index=df.index, columns=df.columns)

    dpi=250
    pc = np.matmul(df.transpose(), pca.components_.transpose())

    cmap = plt.get_cmap('viridis')
    colors = {label: color for color, label in zip(cmap(np.linspace(0, 1, len(pc.index))), pc.index)}
    pc['color'] = [colors[x] for x in pc.index]

    # fig = plt.figure(figsize = (8, 4), dpi=dpi)
    # ax = fig.add_axes([0.1, 0.15, 0.58, 0.75])
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    # for label in pc.index:
    #     sele = pc['treatment'] == treatment
    #     ax.scatter(pc[sele][x_axis_pc], pc[sele][y_axis_pc], color=colors[treatment], label=treatment)

    points=ax.scatter(pc[x_axis_pc], pc[y_axis_pc], c=pc.index, cmap='viridis')

    fig.colorbar(points, label='Acquisition number', use_gridspec=False)

    # axis labels
    pc_var = (pca.singular_values_ ** 2) / sum(pca.singular_values_ ** 2) * 100
    ax.set_xlabel(f'PC {x_axis_pc + 1} {pc_var[x_axis_pc]:.1f}% var')
    ax.set_ylabel(f'PC {y_axis_pc + 1} {pc_var[y_axis_pc]:.1f}% var')

    # ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))

    if fname is None:
        plt.show()
    else:
        plt.savefig('fig/pca_test.pdf')
    plt.close()


query = '''
SELECT
    p.replicateId,
    r.acquiredRank,
    p.modifiedSequence,
    p.precursorCharge,
    p.totalAreaFragment,
FROM precursors p
LEFT JOIN replicates r
ON p.replicateId = r.replicateId;'''

# conn = sqlite3.connect("/home/ajm/code/DIA_QC_report/testData/csf_data.db3")
# conn = sqlite3.connect("/home/ajm/code/DIA_QC_report/testData/full_data.db3")
conn = sqlite3.connect("/home/ajm/code/DIA_QC_report/testData/data.db3")

df = pd.read_sql(query, conn)
metadata = pd.read_sql('SELECT r.acquiredRank', conn)
conn.close()

# set zero areas to the mininum non-zero value
df.loc[df['totalAreaFragment'] == 0, 'totalAreaFragment'] = min(df[df['totalAreaFragment'] > 0]['totalAreaFragment'])

df['log2TotalAreaFragment'] = np.log10(df['totalAreaFragment'])
df['zScore'] = df.groupby(["acquiredRank"])['log2TotalAreaFragment'].apply(lambda x: (x - np.mean(x)) / stdev(x))

df_wide = df.pivot_table(index=['modifiedSequence', 'precursorCharge'],
                         columns="acquiredRank", values='zScore')

pca_plot(df_wide, fname='fig/pca_test.pdf')
# pca_plot(df_wide, fname='fig/csf_pca_test.pdf')

