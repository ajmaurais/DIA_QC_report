
import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from .dtype import Dtype

def pc_matrix(df):
    df_s = StandardScaler().fit_transform(df.transpose())

    pca = PCA()
    pca.fit(df_s)

    df_s = pd.DataFrame(df_s.transpose(), index=df.index, columns=df.columns)
    pc = np.matmul(df_s.transpose(), pca.components_.transpose())
    pc_var = (pca.singular_values_ ** 2) / sum(pca.singular_values_ ** 2) * 100

    return pc, pc_var


def pca_plot(pc_unorm, label_col, pc_norm=None, label_type='discrete',
             fname=None, dpi=250, x_axis_pc=0, y_axis_pc=1, add_title=True):
    
    assert isinstance(pc_unorm, tuple)
    (pc_unorm, pc_var_unorm) = pc_unorm

    if pc_norm is not None:
        assert isinstance(pc_norm, tuple)
        (pc_norm, pc_var_norm) = pc_norm

    # cmap = plt.get_cmap('viridis')

    # if label_type == 'discrete':
    #     colors = {label: color for color, label in zip(cmap(np.linspace(0, 1, len(pc[label_col].drop_duplicates()))),
    #                                                    pc[label_col].drop_duplicates())}

    #     fig = plt.figure(figsize = (6, 4), dpi=dpi)
    #     ax = fig.add_axes([0.1, 0.15, 0.65, 0.75])

    #     for label in sorted(pc[label_col].drop_duplicates()):
    #         sele = pc[label_col] == label
    #         ax.scatter(pc[sele][x_axis_pc], pc[sele][y_axis_pc], color=colors[label], label=label)
    #     ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), title=label_col, alignment='left')

    # elif label_type == 'continuous':
    #     fig, ax = plt.subplots(1, 1, figsize=(6, 4), dpi=dpi)
    #     cmap.set_bad(color='grey')
    #     points=ax.scatter(pc[x_axis_pc], pc[y_axis_pc], c=pc[label_col], cmap=cmap, plotnonfinite=True)
    #     fig.colorbar(points, label=label_col, use_gridspec=False,
    #                  ticks=MaxNLocator(integer=True) if pd.api.types.is_integer_dtype(pc[label_col]) else None)

    # axis labels
    # ax.set_xlabel(f'PC {x_axis_pc + 1} {pc_var[x_axis_pc]:.1f}% var')
    # ax.set_ylabel(f'PC {y_axis_pc + 1} {pc_var[y_axis_pc]:.1f}% var')

    # fig = px.scatter(x=pc[x_axis_pc], y=pc[y_axis_pc], color=pc[label_col])
    fig = make_subplots(rows=1, cols=(1 if pc_norm is None else 2))

    fig.add_trace(
        go.Scatter(x=pc_unorm[x_axis_pc], y=pc_unorm[y_axis_pc],
                   mode='markers', name='Unnormalized'),
        row=1, col=1
    )

    if pc_norm is not None:
        fig.add_trace(
            go.Scatter(x=pc_norm[x_axis_pc], y=pc_norm[y_axis_pc],
                       mode='markers', name='Normalized'),
            row=1, col=2
        )
    
    fig.update_layout(plot_bgcolor='white')
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False)

    fig.show()
    
    # if add_title:
    #     plt.title(f'Colored by {label_col}')

    # if fname is None:
    #     plt.show()
    # else:
    #     plt.savefig(fname)
    # plt.close()


def convert_string_cols(df):
    '''
    Convert string annotation key columns in DataFrame to annotationType
    '''
    types = {row.key: Dtype[row.type] for row in df[['key', 'type']].drop_duplicates().itertuples()}
    ret = df.pivot(index="replicateId", columns="key", values="value")
    for column in ret.columns:
        ret[column] = ret[column].apply(types[column].convert)
        if types[column] is Dtype.STRING:
            ret.loc[ret[column] == '', column] = pd.NA

    return ret.rename_axis(columns=None).reset_index()
