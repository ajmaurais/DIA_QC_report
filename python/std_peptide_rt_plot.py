
import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.ticker import MaxNLocator
from numpy import mean


def peptide_rt_plot(protein_id, conn, fname=None):
    query = '''
    SELECT
        r.acquiredRank,
        p.modifiedSequence,
        p.precursorCharge,
        p.rt,
        p.minStartTime,
        p.maxEndTime
    FROM precursors p
    LEFT JOIN replicates r
        ON p.replicateId = r.replicateId
    WHERE p.proteinName = "%s";
    ''' % protein_id

    df = pd.read_sql(query, conn)

    # rank peptides by mean RT across replicates
    df['meanRT'] = df.groupby('modifiedSequence')['rt'].transform('mean')
    df['peptideLabel'] = df[["modifiedSequence", "precursorCharge"]].apply(lambda x: f'{x[0]}{"+" * x[1]}', axis = 1)
    ranks = {row['peptideLabel']: row['meanRT'] for _, row in df[['meanRT', 'peptideLabel']].drop_duplicates().iterrows()}

    # generate color scale
    cmap = plt.get_cmap('viridis')
    pallet = cmap([i/len(ranks) for i in range(len(ranks))])
    colors = {peptide: color for color, (peptide, rank) in zip(pallet, sorted(ranks.items(), key=lambda x: x[1]))}
    
    # fig, ax = plt.subplots(1, 1)
    fig = plt.figure(figsize = (8, 4), dpi=1000)
    ax = fig.add_axes([0.1, 0.15, 0.58, 0.75])

    line_width = 0.4
    for i, replicate in enumerate(df['acquiredRank'].drop_duplicates()):
        for _, row in df[df['acquiredRank'] == replicate].iterrows():
            ax.add_patch(Rectangle((i - (line_width / 2), row['minStartTime']),
                                    line_width, row['maxEndTime'] - row['minStartTime'],
                                    color=colors[row['peptideLabel']]))

    ax.legend(handles = [Patch(color=color, label=label) for label, color in colors.items()],
              bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    ax.set_xlabel('Acquisition order')
    ax.set_ylabel('RT (min)')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    min_rt = min(df['minStartTime'])
    max_rt = max(df['maxEndTime'])
    rt_range = max_rt - min_rt
    y_padding = rt_range * 0.05

    n_reps = len(df['acquiredRank'].drop_duplicates())
    x_padding = n_reps * 0.02

    plt.xlim([0 - 0.5 - x_padding, n_reps - 0.5 + x_padding])
    plt.ylim([min_rt - y_padding, max_rt + y_padding])
    plt.title(protein_id)

    if fname:
        plt.savefig(fname)
    else:
        plt.show()
    plt.close()

# std_proteins = ['iRT']
# # conn = sqlite3.connect('/home/ajm/code/DIA_QC_report/testData/data.db3')
# conn = sqlite3.connect('/home/ajm/code/DIA_QC_report/testData/full_data.db3')
# 
# peptide_rt_plot(std_proteins[0], conn, 'fig/rt_plot.pdf')
# conn.close()

