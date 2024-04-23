
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from .metadata import Dtype

def pc_matrix(df):
    df_s = StandardScaler().fit_transform(df.transpose())

    pca = PCA()
    pca.fit(df_s)

    df_s = pd.DataFrame(df_s.transpose(), index=df.index, columns=df.columns)
    pc = np.matmul(df_s.transpose(), pca.components_.transpose())
    pc_var = (pca.singular_values_ ** 2) / sum(pca.singular_values_ ** 2) * 100

    return pc, pc_var


def pca_plot(pc, label_col, pc_var, label_type='discrete',
             fname=None, dpi=250, x_axis_pc=0, y_axis_pc=1, add_title=True):

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
        fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=dpi)
        points=ax.scatter(pc[x_axis_pc], pc[y_axis_pc], c=pc[label_col], cmap='viridis')
        fig.colorbar(points, label=label_col, use_gridspec=False,
                     ticks=MaxNLocator(integer=True) if pd.api.types.is_integer_dtype(pc[label_col]) else None)

    # axis labels
    ax.set_xlabel(f'PC {x_axis_pc + 1} {pc_var[x_axis_pc]:.1f}% var')
    ax.set_ylabel(f'PC {y_axis_pc + 1} {pc_var[y_axis_pc]:.1f}% var')

    if add_title:
        plt.title(f'Colored by {label_col}')

    if fname is None:
        plt.show()
    else:
        plt.savefig(fname)
    plt.close()


def convert_string_cols(df):
    '''
    Convert string annotation key columns in DataFrame to annotationType
    '''

    types = {row.key: Dtype[row.type] for row in df[['key', 'type']].drop_duplicates().itertuples()}
    ret = df.pivot(index="replicateId", columns="key", values="value")
    for column in ret.columns:
        ret[column] = ret[column].apply(lambda x: types[column].convert(x))

    return ret.rename_axis(columns=None).reset_index()

