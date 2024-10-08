
from os.path import splitext
from colorsys import hsv_to_rgb

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from .dtype import Dtype

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
        rgb = hsv_to_rgb(hue, 0.6, 0.8)  # Using HSV as a proxy for HCL
        colors.append(f'rgb({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)})')
    return colors


def pca_plot(pc_data, label_col, label_type='discrete',
             fname=None, x_axis_pc=0, y_axis_pc=1, add_title=False):
    '''
    Generate a plotly PCA plot for a PC matrix.

    Parameters
    ----------
    pc_data: dict
        A dictionary where the keys are the column title ('Unnormalized', 'Normalized' etc.),
        and the values are a tuple wher the first element is the PC matrix, and the second
        element are the variances explained by each PC.

    label_col: str
        A column name in each PC matrix to color pannles.

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
        unique_categories = list()
        for data, _ in pc_data.values():
            unique_categories += data[label_col].drop_duplicates().to_list()
        unique_categories = sorted(set(unique_categories))

        category_colors = dict(zip(unique_categories, generate_hcl_colors(len(unique_categories))))

        for i, (title, (pannel_data, pannel_var)) in enumerate(pc_data.items()):
            if pannel_data is not None:
                show_legend = i >= len(pc_data) - 1

                for category in unique_categories:
                    mask = pannel_data[label_col] == category
                    fig.add_trace(
                        go.Scatter(x=pannel_data.loc[mask, x_axis_pc],
                                   y=pannel_data.loc[mask, y_axis_pc],
                                   mode='markers', name=category,
                                   hovertext=pannel_data['replicate'],
                                   marker={'color': category_colors[category]},
                                   showlegend=show_legend),
                        row=1, col=i + 1
                    )

                fig.update_xaxes(title_text=f'PC {x_axis_pc + 1} {pannel_var[x_axis_pc]:.1f}% var', row=1, col=i + 1)
                fig.update_yaxes(title_text=f'PC {y_axis_pc + 1} {pannel_var[y_axis_pc]:.1f}% var', row=1, col=i + 1)

    elif label_type == 'continuous':
        for i, (title, (pannel_data, pannel_var)) in enumerate(pc_data.items()):
            if pannel_data is not None:
                show_legend = i >= len(pc_data) - 1

                marker = go.scatter.Marker(color=pannel_data[label_col], colorscale='Viridis',
                                           showscale=show_legend,
                                           colorbar={'title': label_col} if show_legend else None)

                fig.add_trace(
                    go.Scatter(x=pannel_data[x_axis_pc], y=pannel_data[y_axis_pc],
                               mode='markers', name=title,
                               hovertext=pannel_data['replicate'],
                               marker=marker, showlegend=False),
                    row=1, col=i + 1
                )

                fig.update_xaxes(title_text=f'PC {x_axis_pc + 1} {pannel_var[x_axis_pc]:.1f}% var', row=1, col=i + 1)
                fig.update_yaxes(title_text=f'PC {y_axis_pc + 1} {pannel_var[y_axis_pc]:.1f}% var', row=1, col=i + 1)

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

    # legend_width = estimate_legend_width(fig)
    legend_width = 120
    fig.update_layout(plot_bgcolor='white',
                      height=450, width=(400 * len(pc_data) + 150),
                      legend={'title': label_col})

    if add_title:
        fig.update_layout(title=f'Colored by {label_col}')

    if fname is None:
        return fig
    else:
        ext = splitext(fname)[1]
        if ext == '.html':
            fig.write_html(fname)
        else:
            fig.write_image(fname)


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
