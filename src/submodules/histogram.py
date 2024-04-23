
import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np


def histogram(dat, xlab, ylab='Count', limits=None, fname=None, dpi=250):
    fig, ax = plt.subplots(1, 1, figsize = (6, 4), dpi=dpi)
    ax.hist(dat, bins = 100)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    if limits:
        plt.xlim(limits)

    if fname:
        plt.savefig(fname)
    else:
        plt.show()
    plt.close()


def box_plot(data, ylab, xlab='Acquisition number', fname=None, hline=None, limits=None, dpi=250):

    # init plot
    fig, ax = plt.subplots(1, 1, figsize = (10, 4), dpi=dpi)

    # plot outliers as scatter
    for i, c in enumerate(data.columns):
        y = data[c]
        q3, q1 = np.percentile(y, [75, 25])
        iqr = q3 - q1
        outliers = y[(y > q3 + (1.5 * iqr)) | (y < q1 - (1.5 * iqr))]
        x = np.random.normal(i + 1, 0.04, size=len(outliers))
        ax.scatter(x, outliers, c = 'black', s = 0.3, alpha = 0.5)

    if hline is not None:
        plt.axhline(y=hline, color='black', linestyle=':', linewidth=1)

    # plot boxes
    ax.boxplot(data, showfliers=False)

    if limits:
        plt.ylim(limits)

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    n_cols = len(data.columns)
    plt.xticks(np.linspace(1, n_cols, num=min(6, n_cols), dtype=int),
               np.linspace(0, n_cols - 1, num=min(6, n_cols), dtype=int))

    if fname:
        plt.savefig(fname)
    else:
        plt.show()
    plt.close()


def multi_boxplot(dats, data_levels,
                  xlab='Acquisition number', ylab='Log2(Area)',
                  fname=None, dpi=250):
    '''
    Draw box plots with a row for each level in data_levels.

    Parameters
    ----------
    dats: dict
        A dict where the key is the data level name and the value is a
        wide formated pd.DataFrame
    data_levels: dict
        A dict where the key is the data level name and the value is the
        index of the level.
    xlab: str
        X axis label
    ylab: str
        Y axis label
    fname: str
        The name of the file to save the plot.
    dpi: int
        250 is the default.
    '''

    fig, axs = plt.subplots(len(data_levels), 1,
                            figsize = (10, len(data_levels) * 3),
                            dpi=dpi)

    if len(data_levels) == 1:
        axs = [axs]

    # iterate through each level
    for level, index in data_levels.items():

        # plot outliers as scatter
        for i, c in enumerate(dats[level].columns):
            y = dats[level][c]
            q3, q1 = np.percentile(y, [75, 25])
            iqr = q3 - q1
            outliers = y[(y > q3 + (1.5 * iqr)) | (y < q1 - (1.5 * iqr))]
            x = np.random.normal(i + 1, 0.04, size=len(outliers))
            axs[index].scatter(x, outliers, c = 'black', s = 0.3, alpha = 0.5)

        axs[index].boxplot(dats[level], showfliers=False) #, vert=False)

        if index == len(data_levels) - 1:
            axs[index].set_xlabel(xlab)

        axs[index].set_ylabel(ylab)

        n_cols = len(dats[level].columns)
        axs[index].set_xticks(np.linspace(1, n_cols, num=min(6, n_cols), dtype=int),
                              np.linspace(0, n_cols - 1, num=min(6, n_cols), dtype=int))
        axs[index].set_title(level)

    if fname:
        plt.savefig(fname)
    else:
        plt.show()
    plt.close()

