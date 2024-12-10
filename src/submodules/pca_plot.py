
from os.path import splitext
from colorsys import hsv_to_rgb

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from plotly.subplots import make_subplots
import plotly.graph_objects as go

from .dtype import Dtype


NA_CATEGORY_NAME = '<NA>'
GREY_RGB = (0.36, 0.36, 0.36)


def calc_pc_matrix(df):
    df_s = StandardScaler().fit_transform(df.transpose())

    pca = PCA()
    pca.fit(df_s)

    df_s = pd.DataFrame(df_s.transpose(), index=df.index, columns=df.columns)
    pc = np.matmul(df_s.transpose(), pca.components_.transpose())
    pc_var = (pca.singular_values_ ** 2) / sum(pca.singular_values_ ** 2) * 100

    return pc, pc_var


def generate_hcl_colors(n):
    ''' Generate `n` evenly spaced colors around the HCL color wheel. '''

    colors = []
    for i in range(n):
        hue = i / n  # Evenly spaced around the color wheel
        rgb = hsv_to_rgb(hue, 0.76, 0.82)  # Using HSV as a proxy for HCL
        colors.append(rgb)
    return colors


def map_discrete_colors(replicates, plotly_color=False):

    assert isinstance(replicates, dict)
    assert all(isinstance(key, (str, int)) for key in replicates)
    assert all(isinstance(value, (str, float, int, bool, type(None), type(pd.NA)))
               for value in replicates.values())

    categories = dict()
    for key in replicates:
        if replicates[key] is None or pd.isna(replicates[key]):
            categories[key] = NA_CATEGORY_NAME
        else:
            categories[key] = str(replicates[key])

    unique_categories = sorted(set(categories.values()))
    non_na = [c for c in unique_categories if c != NA_CATEGORY_NAME]

    colors = dict(zip(non_na, generate_hcl_colors(len(non_na))))

    if len(non_na) != len(unique_categories):
        colors[NA_CATEGORY_NAME] = GREY_RGB

    if plotly_color:
        for category, color in colors.items():
            colors[category] = f"rgb({', '.join([str(int(v * 255)) for v in color])})"

    return categories, colors


def mpl_pca_plot(pc_data, metadata, label_col, label_type='discrete',
                 fname=None, x_axis_pc=0, y_axis_pc=1, add_title=True, dpi=250):
    '''
    Generate a Matplotlib PCA plot for a PC matrix.

    Parameters
    ----------
    pc_data: dict
        A dictionary where the keys are the column title ('Unnormalized', 'Normalized' etc.),
        and the values are a tuple where the first element is the PC matrix, and the second
        element are the variances explained by each PC.

    metadata: pd.DataFrame
        A pandas dataframe with replicate metadata.
        Metadata is joined to PC data by the replicateId column.

    label_col: str
        A variable in the metadata to color pannles.

    label_type: str
        The label type of `label_col` ('discrete' or 'continious')

    fname: str
        The name of the file to write figure to. If None, the plt.show method is called.

    x_axis_pc: int
        The X axis PC index. Default is 0.

    y_axis_pc: int
        The Y axis PC index. Default is 1.

    add_title: bool
        Add title to plot? Default is False.

    dpi: int
        The plot dpi if figure is being written to a file.
    '''

    data_levels = list(pc_data.keys())

    fig = plt.figure(figsize=(len(pc_data) * 4 + 1, 4), dpi=dpi, layout='constrained')
    axs = fig.subplots(nrows=1, ncols=len(pc_data))

    axs = np.atleast_1d(axs)

    if label_type == 'discrete':
        categories = metadata[['replicateId', label_col]].set_index('replicateId')
        categories = categories[label_col].to_dict()
        local_metadata, colors = map_discrete_colors(categories)

        for i, level in enumerate(data_levels):

            plot_df = pc_data[level][0]
            plot_df[label_col] = plot_df.index.map(lambda x: local_metadata[x])

            for label, color in colors.items():
                sele = plot_df[label_col] == label
                axs[i].scatter(plot_df[sele][x_axis_pc], plot_df[sele][y_axis_pc],
                               color=color, label=label)

            axs[-1].legend(loc='upper left', bbox_to_anchor=(1.05, 1),
                           title=label_col, alignment='left')

    elif label_type == 'continuous':
        cmap = plt.get_cmap('viridis')
        cmap.set_bad(color='grey')

        for i, level in enumerate(data_levels):
            plot_df = pc_data[level][0]
            local_metadata = metadata[['replicateId', label_col]].set_index('replicateId')
            plot_df = plot_df.join(local_metadata)

            points = axs[i].scatter(plot_df[x_axis_pc], plot_df[y_axis_pc],
                                    c=plot_df[label_col], cmap=cmap, plotnonfinite=True)

        ticks = MaxNLocator(integer=True) if pd.api.types.is_integer_dtype(plot_df[label_col]) else None
        fig.colorbar(points, ax=axs, label=label_col,
                     use_gridspec=False,
                     shrink=0.8,
                     # pad=0.1,
                     ticks=ticks)

    else:
        raise ValueError('Unknown label_type!')

    # axis labels
    for i, _ in enumerate(data_levels):
        axs[i].set_xlabel(f'PC {x_axis_pc + 1} {pc_data[data_levels[i]][1][x_axis_pc]:.1f}% var')
        axs[i].set_ylabel(f'PC {y_axis_pc + 1} {pc_data[data_levels[i]][1][y_axis_pc]:.1f}% var')
        axs[i].set_title(data_levels[i])

    if add_title:
        plt.suptitle(f'Colored by {label_col}')
    if fname is None:
        plt.show()
    else:
        plt.savefig(fname)
    plt.close()


def plotly_pca_plot(pc_data, metadata, label_col, label_type='discrete',
                    fname=None, x_axis_pc=0, y_axis_pc=1, add_title=False):
    '''
    Generate a plotly PCA plot for a PC matrix.

    Parameters
    ----------
    pc_data: dict
        A dictionary where the keys are the column title ('Unnormalized', 'Normalized' etc.),
        and the values are a tuple where the first element is the PC matrix, and the second
        element are the variances explained by each PC.

    metadata: pd.DataFrame
        A pandas dataframe with replicate metadata.
        Metadata is joined to PC data by the replicateId column.

    label_col: str
        A variable in the metadata to color pannles.

    label_type: str
        The label type of `label_col` ('discrete' or 'continious')

    fname: str
        The name of the file to write figure to. If None, the plot object is returned.

    x_axis_pc: int
        The X axis PC index. Default is 0.

    y_axis_pc: int
        The Y axis PC index. Default is 1.

    add_title: bool
        Add title to plot? Default is False.
    '''

    assert isinstance(pc_data, dict)
    for k, v in pc_data.items():
        assert isinstance(v, tuple)

    fig = make_subplots(rows=1, cols=len(pc_data),
                        subplot_titles=list(pc_data.keys()))

    if label_type == 'discrete':
        categories = metadata[['replicateId', label_col]].set_index('replicateId')
        categories = categories[label_col].to_dict()
        local_metadata, colors = map_discrete_colors(categories, plotly_color=True)

        for i, level in enumerate(pc_data):
            show_legend = i >= len(pc_data) - 1

            plot_df = pc_data[level][0].join(metadata[['replicateId', 'replicate']].set_index('replicateId'))
            plot_df[label_col] = plot_df.index.map(lambda x: local_metadata[x])

            for label, color in colors.items():
                mask = plot_df[label_col] == label
                fig.add_trace(
                    go.Scatter(x=plot_df.loc[mask, x_axis_pc],
                               y=plot_df.loc[mask, y_axis_pc],
                               mode='markers', name=label,
                               hovertext=plot_df['replicate'],
                               marker={'color': color},
                               showlegend=show_legend),
                    row=1, col=i + 1
                )

            fig.update_xaxes(title_text=f'PC {x_axis_pc + 1} {pc_data[level][1][x_axis_pc]:.1f}% var', row=1, col=i + 1)
            fig.update_yaxes(title_text=f'PC {y_axis_pc + 1} {pc_data[level][1][y_axis_pc]:.1f}% var', row=1, col=i + 1)

    elif label_type == 'continuous':
        for i, level in enumerate(pc_data):
            show_legend = i >= len(pc_data) - 1

            plot_df = pc_data[level][0]
            local_metadata = metadata[['replicateId', 'replicate', label_col]].set_index('replicateId')
            plot_df = plot_df.join(local_metadata)

            marker = go.scatter.Marker(color=plot_df[label_col], colorscale='Viridis',
                                       showscale=show_legend,
                                       colorbar={'title': label_col} if show_legend else None)

            fig.add_trace(
                go.Scatter(x=plot_df[x_axis_pc], y=plot_df[y_axis_pc],
                           mode='markers', name=level,
                           hovertext=plot_df['replicate'],
                           marker=marker, showlegend=False),
                row=1, col=i + 1
            )

            fig.update_xaxes(title_text=f'PC {x_axis_pc + 1} {pc_data[level][1][x_axis_pc]:.1f}% var', row=1, col=i + 1)
            fig.update_yaxes(title_text=f'PC {y_axis_pc + 1} {pc_data[level][1][y_axis_pc]:.1f}% var', row=1, col=i + 1)

    line_width = None
    axis_format_args = {'showgrid': False,
                        'zeroline': False,
                        'showline': True,
                        'linecolor': 'black',
                        'linewidth': line_width,
                        'mirror': True,
                        'ticks': 'outside',
                        'tickcolor': 'black',
                        'tickwidth': line_width}

    fig.update_xaxes(**axis_format_args)
    fig.update_yaxes(**axis_format_args)

    fig.update_layout(plot_bgcolor='white',
                      height=480, width=(450 * len(pc_data) + 150),
                      legend={'title': label_col})

    fig.update_traces(marker={'size': 8})

    if add_title:
        fig.update_layout(title=f'Colored by {label_col}')

    if fname is None:
        fig.show()
    else:
        ext = splitext(fname)[1]
        if ext == '.html':
            fig.write_html(fname)
        else:
            fig.write_image(fname)
