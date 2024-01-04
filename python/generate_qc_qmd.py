
import sys
import os
import argparse
import re
from inspect import stack

PYTHON_DIR = os.path.dirname(os.path.abspath(__file__))

DEFAULT_OFNAME = 'qc_report.qmd'
DEFAULT_EXT = 'html'
DEFAULT_TITLE = 'DIA QC report'
DEFAULT_DPI = 250

PEPTIDE_QUERY = '''
query = \'\'\'SELECT 
    r.acquiredRank,
    p.peptide,
    p.modifiedSequence,
    p.precursorCharge
FROM precursors p
LEFT JOIN replicates r
    ON p.replicateId = r.replicateId
WHERE p.totalAreaFragment > 0;\'\'\'

df = pd.read_sql(query, conn)
df = df.drop_duplicates()
'''

PRECURSOR_QUERY = '''
query = \'\'\'SELECT 
    p.replicateId,
    r.acquiredRank,
    p.modifiedSequence,
    p.precursorCharge,
    p.totalAreaFragment,
    p.totalAreaMs1,
    p.rt,
    p.maxFwhm,
    p.averageMassErrorPPM as massError,
    p.libraryDotProduct,
    p.isotopeDotProduct
FROM precursors p
LEFT JOIN replicates r
    ON p.replicateId = r.replicateId;\'\'\'

df = pd.read_sql(query, conn)
df = df.drop_duplicates()
'''


def doc_header(command, title=DEFAULT_TITLE):
    header='''<!-- This document was automatically generated. Do not edit. -->
<!-- %s -->\n
---
title: "%s"
toc: true
format:
    html:
        code-fold: true
        grid:
            body-width: 1000px
            sidebar-width: 0px
            margin-width: 300px
        self-contained: true
    pdf:
        fig-pos: 'H'
        extra_dependencies: ["float"]
        fig-format: 'png'
jupyter: python3
---\n\n''' % (command, title)
    return header


def add_header(text, level=1):
    return f'{"#" * level} {text}\n'


def python_block_header(label):
    text = '''```{python}
#| label: %s
#| echo: false
#| warning: false
#| fig-align: "left"''' % label
    return text


def doc_initialize(db_path, python_dir):
    text='''\n%s\n
import sys
import re
import sqlite3
from statistics import stdev
from scipy.stats import zscore
import numpy as np
import pandas as pd

sys.path.append('%s')
from pyDIAUtils.std_peptide_rt_plot import peptide_rt_plot
from pyDIAUtils.bar_chart import bar_chart
from pyDIAUtils.histogram import histogram, box_plot
from pyDIAUtils.pca_plot import pc_matrix, pca_plot, convert_string_cols

conn = sqlite3.connect('%s')
```\n\n''' % (python_block_header(stack()[0][3]), python_dir, db_path)

    return text


def doc_finalize():
    text='''\n%s\n
    conn.close()
```\n\n''' % (python_block_header(stack()[0][3]))

    return text


# plots

def replicate_tic_areas(dpi=DEFAULT_DPI):
    text='''\n%s\n
# replicate tic bar chart
tic = pd.read_sql('SELECT acquiredRank, ticArea FROM replicates', conn)
bar_chart(tic, 'TIC area', dpi=%s)
```\n\n''' % (python_block_header(stack()[0][3]), dpi)
    return text


def std_rt_dist(std_proteins, dpi=DEFAULT_DPI):
    text = '''\n%s\n
for std in [%s]:
    peptide_rt_plot(std, conn, dpi=%i)
```\n\n''' % (python_block_header(stack()[0][3]),
              "'{}'".format("', '".join(std_proteins)) if len(std_proteins) > 0 else '',
              dpi)
    return text


def missed_cleavages(do_query=True, dpi=DEFAULT_DPI):
    text = '''\n%s\n%s
# precursor bar chart colored by missed cleavages
trypsin_re = re.compile('([RK])(?=[^P])')
df['nMissed'] = df['peptide'].apply(lambda x: len(trypsin_re.findall(x)))
agg = df.groupby(["acquiredRank", "nMissed"])["nMissed"].agg(["count"]).reset_index()
agg = agg.pivot_table(index='acquiredRank', columns='nMissed', values='count')
bar_chart(agg, 'Number of precursors', legend_title='Missed cleavages', dpi=%i)
```\n\n''' % (python_block_header(stack()[0][3]), PEPTIDE_QUERY if do_query else '\n', dpi)

    return text


def precursor_charges(do_query=True, dpi=DEFAULT_DPI):
    text = '''\n%s\n%s
# precursor bar chart colored by precursorCharge
agg = df.groupby(["acquiredRank", "precursorCharge"])["precursorCharge"].agg(["count"]).reset_index()
agg = agg.pivot_table(index='acquiredRank', columns='precursorCharge', values='count')
bar_chart(agg, 'Number of precursors', legend_title='Precursor charge', dpi=%i)
```\n\n''' % (python_block_header(stack()[0][3]), PEPTIDE_QUERY if do_query else '\n', dpi)

    return text


def rep_rt_sd(do_query=True, dpi=DEFAULT_DPI):
    text = '''\n%s\n%s\n
# replicate RT cor histogram
agg = df.groupby(['modifiedSequence', 'precursorCharge'])['rt'].agg(np.std)
histogram(agg, 'Replicate RT SD (min)', dpi=%i,
          limits=(agg.quantile(0.05) * -3, agg.quantile(0.95) * 3))
```\n\n''' % (python_block_header(stack()[0][3]), PRECURSOR_QUERY if do_query else '\n', dpi)

    return text


def _pivot_box_plot(values, title=None, dpi=DEFAULT_DPI):
    if title is None:
        title = ""
        for i, c in enumerate(values):
            if i == 0:
                title += c.upper()
                continue
            if c.isupper():
                title += f" {c.lower()}"
            else:
                title += c
    
    text = '''
data = df.pivot_table(index=['modifiedSequence', 'precursorCharge'],
    columns="acquiredRank", values='%s', aggfunc=sum)
box_plot(data, '%s', dpi=%i)\n''' % (values, title, dpi)

    return text


def cammel_to_underscore(s):
    ret = "".join("_" + c.lower() if c.isupper() else c for c in s).strip("_")
    ret = re.sub("(.)([A-Z])", r"\1_\2", ret).lower()
    return ret


def pivot_box_plot(values, do_query=True, **kwargs):
    text = '\n{}\n{}{}```\n\n'.format(python_block_header(cammel_to_underscore(values)),
                                      PRECURSOR_QUERY if do_query else '',
                                      _pivot_box_plot(values, **kwargs))
    return text

def ms1_ms2_ratio(do_query=True, dpi=DEFAULT_DPI):
    text = '''\n%s\n%s
df['ms2_ms1_ratio'] = df['totalAreaFragment'] / df['totalAreaMs1']
df['ms2_ms1_ratio'] = df['ms2_ms1_ratio'].apply(lambda x: x if np.isfinite(x) else 0)
data = df.pivot_table(index=['modifiedSequence', 'precursorCharge'],
                 columns="acquiredRank", values='ms2_ms1_ratio', aggfunc=sum)
box_plot(data, 'Transition / precursor ratio',
         limits = (0, df['ms2_ms1_ratio'].quantile(0.95) * 3), dpi=%i)
```\n\n''' % (python_block_header(stack()[0][3]), PRECURSOR_QUERY if do_query else '\n', dpi)

    return text


def pc_analysis(do_query, dpi):
    text = '''\n%s\n%s
df_pc = df[['replicateId', 'modifiedSequence', 'precursorCharge', 'totalAreaFragment']]
df_pc['acquisition_number'] = df['acquiredRank']

# set zero areas to the minimum non-zero value
df_pc.loc[df_pc['totalAreaFragment'].apply(lambda x: not np.isfinite(x) or x == 0), 'totalAreaFragment'] = min(df_pc[df_pc['totalAreaFragment'] > 0]['totalAreaFragment'])

df_pc['log2TotalAreaFragment'] = np.log10(df_pc['totalAreaFragment'])
df_pc['zScore'] = df_pc.groupby('acquisition_number')['log2TotalAreaFragment'].transform(lambda x: np.abs(zscore(x)))

df_wide = df_pc.pivot_table(index=['modifiedSequence', 'precursorCharge'],
                         columns="replicateId", values='zScore')

# actually do pc analysis
pc, pc_var = pc_matrix(df_wide)

# meta_values = {'gender': 'discrete', 'vital_status': 'discrete',
#                'tumor_grade': 'discrete', 'year_of_birth': 'continuous'}

# check if metadata table exists in db
# cur = conn.cursor()
# cur.execute('SELECT name FROM sqlite_master WHERE type="table" AND name="metadata"')
# table = cur.fetchall()
# 
# metadata_length = 0
# if len(table) > 0:
#     cur = conn.cursor()
#     cur.execute('SELECT COUNT(*) FROM sampleMetadata')
#     metadata_length = cur.fetchall()[0][0]

acquired_ranks = df_pc[["replicateId", "acquisition_number"]].drop_duplicates()
acquired_ranks = pd.DataFrame(acquired_ranks['acquisition_number'], index=acquired_ranks['replicateId'])

# add metadata to pc matrix
# if metadata_length > 0:
#     METADATA_QUERY = 'SELECT replicateId, annotationKey as columns, annotationValue FROM sampleMetadata WHERE annotationKey IN ("{}")'.format('", "'.join(meta_values.keys()))
# 
#     # get metadata labels for pca plot
#     metadata_df = pd.read_sql(METADATA_QUERY, conn)
#     metadata = metadata_df.pivot(index="replicateId", columns="columns", values="annotationValue")
#     metadata = acquired_ranks.join(metadata)
#     meta_values['acquisition_number'] = 'continuous'
#     metadata = convert_string_cols(metadata)
# 
#     # join metadata to pc matrix
#     pc = pc.join(metadata)
#     for label_name, label_type in meta_values.items():
#         if label_type == 'discrete':
#             if any(pc[label_name].apply(pd.isna)):
#                 warnings.warn('Missing label values!', Warning)
#             pc[label_name] = pc[label_name].apply(str)
#         elif label_type == 'continuous':
#             if any(pc[label_name].apply(pd.isna)):
#                 raise RuntimeError('Cannot have missing label values in continuous scale!')
#         else:
#             raise RuntimeError(f'"{label_type}" is an unknown label_type!')
# 
# else:
meta_values = {'acquisition_number': 'continuous'}
pc = pc.join(acquired_ranks)

for label_name, label_type in sorted(meta_values.items(), key=lambda x: x[0]):
    pca_plot(pc, label_name, pc_var, label_type=label_type, dpi=%i)
```\n\n''' % (python_block_header(stack()[0][3]), PRECURSOR_QUERY if do_query else '\n', dpi)

    return text


def main():
    parser = argparse.ArgumentParser(description='Generate qmd report.')
    parser.add_argument('--dpi', default=DEFAULT_DPI, type=int,
                        help=f'Figure DPI in report. {DEFAULT_DPI} is the default.')
    parser.add_argument('-o', '--ofname', default=f'{DEFAULT_OFNAME}',
                        help=f'Output file basename. Default is "{DEFAULT_OFNAME}"')
    # parser.add_argument('-f', '--format', default=f'{DEFAULT_EXT}',
    #                     help=f'Output file format. Default is "{DEFAULT_EXT}"')
    parser.add_argument('-a', '--addStdProtein', action='append', default=[],
                        help='Add standard protein name for retention time plot.')
    parser.add_argument('--title', default=DEFAULT_TITLE,
                        help=f'Report title. Default is "{DEFAULT_TITLE}"')
    parser.add_argument('db', help='Path to precursor quality database.')
    args = parser.parse_args()

    # ofname = f'{args.ofname}.{args.format}'

    with open(args.ofname, 'w') as outF:
        outF.write(doc_header('{} {}'.format(os.path.abspath(__file__),
                              ' '.join(sys.argv[1:])), title=args.title))
        outF.write(add_header('Peptide independent metrics', level=1))
        outF.write(doc_initialize(args.db, PYTHON_DIR))

        outF.write(add_header('Replicate TIC areas', level=2))
        outF.write(replicate_tic_areas(dpi=args.dpi))

        outF.write(add_header('Peptide dependent metrics\n', level=1))
        outF.write(add_header('Standard retention times across replicates', level=2))
        outF.write(std_rt_dist(args.addStdProtein, dpi=args.dpi))

        outF.write(add_header('Number of missed cleavages', level=2))
        outF.write(missed_cleavages(do_query=True, dpi=args.dpi))

        outF.write(add_header('Precursor charges', level=2))
        outF.write(precursor_charges(do_query=False, dpi=args.dpi))

        outF.write(add_header('Mass error distribution', level=2))
        outF.write(pivot_box_plot('massError', title='Mass error (ppm)',
                                  do_query=True, dpi=args.dpi))

        outF.write(add_header('Peak width distribution', level=2))
        outF.write(pivot_box_plot('maxFwhm', title='Peak FWHM',
                                  do_query=False, dpi=args.dpi))

        outF.write(add_header('Library dot product distribution', level=2))
        outF.write(pivot_box_plot('libraryDotProduct', do_query=False, dpi=args.dpi))

        outF.write(add_header('MS1 isotopic window dot product distribution', level=2))
        outF.write(pivot_box_plot('isotopeDotProduct', do_query=False, dpi=args.dpi))

        outF.write(add_header('MS1 isotopic window dot product distribution', level=2))
        outF.write(ms1_ms2_ratio(do_query=False, dpi=args.dpi))

        outF.write(add_header('Batch effects', level=2))
        outF.write(pc_analysis(do_query=False, dpi=args.dpi))

        outF.write(doc_finalize())

if __name__ == '__main__':
    main()

