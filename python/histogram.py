
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


def box_plot(data, ylab, xlab='Acquired rank', fname=None, hline=None, limits=None, dpi=250):

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


# # read precursors from db
# query = '''
# SELECT 
#     r.acquiredRank,
#     p.modifiedSequence,
#     p.precursorCharge,
#     p.totalAreaFragment,
#     p.totalAreaMs1,
#     p.rt,
#     p.maxFwhm,
#     p.averageMassErrorPPM as massError,
#     p.libraryDotProduct,
#     p.isotopeDotProduct
# FROM precursors p
# LEFT JOIN replicates r
#     ON p.replicateId = r.replicateId;
# '''
# 
# # conn = sqlite3.connect('/home/ajm/code/DIA_QC_report/testData/data.db3')
# conn = sqlite3.connect('/home/ajm/code/DIA_QC_report/testData/full_data.db3')
# df = pd.read_sql(query, conn)
# df = df.drop_duplicates()
# 
# # peptide charge histogram
# 
# # mass accuracy boxplot
# # ylab='Mass error (ppm)'
# # fname='mass_error.png'
# ext = 'png'
# # ext = 'pdf'
# 
# # replicate RT cor histogram
# agg = df.groupby(['modifiedSequence', 'precursorCharge'])['rt'].agg(np.std)
# histogram(agg, 'Replicate RT SD (min)',
#           limits = (agg.quantile(0.05) * -3, agg.quantile(0.95) * 3),
#           fname = f'fig/replicate_rt_cor.{ext}')
# 
# data = df.pivot_table(index=['modifiedSequence', 'precursorCharge'],
#          columns="acquiredRank", values='massError', aggfunc=sum)
# box_plot(data, 'Mass error (ppm)', fname=f'fig/mass_error.{ext}', hline=0)
# 
# data = df.pivot_table(index=['modifiedSequence', 'precursorCharge'],
#          columns="acquiredRank", values='maxFwhm', aggfunc=sum)
# box_plot(data, 'Peak FWHM', fname=f'fig/peak_fwhm.{ext}')
# 
# data = df.pivot_table(index=['modifiedSequence', 'precursorCharge'],
#          columns="acquiredRank", values='libraryDotProduct', aggfunc=sum)
# box_plot(data, 'Library dot product', fname=f'fig/library_dotp.{ext}')
# 
# data = df.pivot_table(index=['modifiedSequence', 'precursorCharge'],
#          columns="acquiredRank", values='isotopeDotProduct', aggfunc=sum)
# box_plot(data, 'Isotope dot product', fname=f'fig/isotope_dotp.{ext}')
# 
# df['ms2_ms1_ratio'] = df['totalAreaFragment'] / df['totalAreaMs1']
# df['ms2_ms1_ratio'] = df['ms2_ms1_ratio'].apply(lambda x: x if np.isfinite(x) else 0)
# data = df.pivot_table(index=['modifiedSequence', 'precursorCharge'],
#          columns="acquiredRank", values='ms2_ms1_ratio', aggfunc=sum)
# box_plot(data, 'Transition / precursor ratio', fname=f'fig/ms2_ms1_ratio.{ext}',
#          limits = (0, df['ms2_ms1_ratio'].quantile(0.95) * 3))


