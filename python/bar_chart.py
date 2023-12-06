
import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re


def bar_chart(dat, ylab, xlab='Acquisition order', legend_title=None, fname=None, dpi=250):

    '''
    Make bar chart from aggregated bar chart.
    '''

    data = {c: dat[c].to_numpy() for c in dat.columns}
    labels = dat.index.to_numpy()
    
    fig = plt.figure(figsize = (6, 3), dpi=dpi)
    
    ax = fig.add_axes([0.1, 0.15, 0.86, 0.65])
    
    bottom = np.zeros(len(dat))
    for n_missed, counts in data.items():
        ax.bar(labels, counts, 0.5, label = n_missed, bottom = bottom)
        bottom += counts
    
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    n_cols = len(dat.index)
    plt.xticks(np.linspace(0, n_cols - 1, num=min(6, n_cols), dtype=int),
               np.linspace(0, n_cols - 1, num=min(6, n_cols), dtype=int))


    if legend_title:
        ax.legend(title=legend_title,
                  bbox_to_anchor=(0., 1.05, 1.0, 0), loc='lower right',
                  ncols=2, borderaxespad=-0.1)

    if fname:
        plt.savefig(fname)
    else:
        plt.show()
    plt.close()


# read precursors from db
# query = '''
# SELECT 
#     r.acquiredRank,
#     p.peptide,
#     p.modifiedSequence,
#     p.precursorCharge
# FROM precursors p
# LEFT JOIN replicates r
#     ON p.replicateId = r.replicateId
# WHERE p.totalAreaFragment > 0;
# '''

# # conn = sqlite3.connect('/home/ajm/code/DIA_QC_report/testData/data.db3')
# conn = sqlite3.connect('/home/ajm/code/DIA_QC_report/testData/full_data.db3')
# df = pd.read_sql(query, conn)
# df = df.drop_duplicates()
# 
# # precursor bar chart colored by missed cleavages
# trypsin_re = re.compile('([RK])(?=[^P])')
# df['nMissed'] = df['peptide'].apply(lambda x: len(trypsin_re.findall(x)))
# agg = df.groupby(["acquiredRank", "nMissed"])["nMissed"].agg(["count"]).reset_index()
# agg = agg.pivot_table(index='acquiredRank', columns='nMissed', values='count')
# bar_chart(agg, 'Number of precursors', legend_title='Missed cleavages', fname='fig/n_missed_dist.pdf')
# 
# # precursor bar chart colored by precursorCharge
# agg = df.groupby(["acquiredRank", "precursorCharge"])["precursorCharge"].agg(["count"]).reset_index()
# agg = agg.pivot_table(index='acquiredRank', columns='precursorCharge', values='count')
# bar_chart(agg, 'Number of precursors', legend_title='Missed cleavages', fname='fig/peptide_charge_dist.pdf')
# 
# # replicate tic bar chart
# tic = pd.read_sql('SELECT acquiredRank, ticArea FROM replicates', conn)
# bar_chart(tic, 'TIC area', fname='fig/replicate_tics.pdf')

