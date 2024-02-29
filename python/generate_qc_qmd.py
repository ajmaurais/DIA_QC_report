import sys
import os
import argparse
import re
import sqlite3
from inspect import stack

from enum import Enum

from pyDIAUtils.metadata import Dtype
from pyDIAUtils.logger import LOGGER
from pyDIAUtils.dia_db_utils import is_normalized

DEFAULT_OFNAME = 'qc_report.qmd'
DEFAULT_EXT = 'html'
DEFAULT_TITLE = 'DIA QC report'
DEFAULT_DPI = 250

PEPTIDE_QUERY = '''
query = \'\'\'SELECT
    r.acquiredRank,
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


def open_pannel_tabset():
    return '\n::: { .panel-tabset }\n\n'


def close_pannel_tabset():
    return '\n:::\n'


def python_block_header(label):
    text = '''```{python}
#| label: %s
#| echo: false
#| warning: false
#| fig-align: "left"''' % label
    return text


def doc_initialize(db_path):
    text='''\n%s\n
import sys
import re
import sqlite3
from statistics import stdev
from scipy.stats import zscore
import numpy as np
import pandas as pd

from pyDIAUtils.std_peptide_rt_plot import peptide_rt_plot
from pyDIAUtils.bar_chart import bar_chart
from pyDIAUtils.histogram import histogram, box_plot, multi_boxplot
from pyDIAUtils.pca_plot import pc_matrix, pca_plot, convert_string_cols

conn = sqlite3.connect('%s')
```\n\n''' % (python_block_header(stack()[0][3]), db_path)

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
trypsin_re = re.compile(r'([RK])(?=[^P])')
mod_re = re.compile(r'\[[\w+]\]')
df['peptide'] = df['modifiedSequence'].apply(lambda x: mod_re.sub('', x))
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
    ret = re.sub(r"(.)([A-Z])", r"\1_\2", ret).lower()
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


def get_meta_key_types(db_path, keys):
    '''
    Determine whether each metadata key is discrete or continious.
    '''

    if keys is None or len(keys) == 0:
        return {}

    def continious_discrete(var_type):
        if var_type is Dtype['FLOAT'] or var_type is Dtype['INT']:
            return 'continuous'
        return 'discrete'

    with sqlite3.connect(db_path) as conn:
        cur = conn.cursor()
        cur.execute('SELECT annotationKey as key, annotationType as type FROM sampleMetadataTypes;')
        db_types = [(x[0], Dtype[x[1]]) for x in cur.fetchall() if x[0] in keys]

    types = {key: continious_discrete(var_type) for key, var_type in db_types}

    return types


def precursor_areas(quant_cols):

    quant_col_query = ',\n'.join([f'{key} as {value}' for key, value in quant_cols.items()])

    Stack = stack()
    text = f"""\n{python_block_header(Stack[0][3])}

query = '''SELECT
    p.replicateId,
    r.acquiredRank as acquisitionNumber,
    p.modifiedSequence,
    p.precursorCharge,
    {quant_col_query}
FROM precursors p
LEFT JOIN replicates r ON p.replicateId = r.replicateId;'''

df_areas = pd.read_sql(query, conn)\n"""

    data_levels = {v: i for i, v in enumerate(quant_cols.values())}
    text += f'data_levels = {str(data_levels)}\n'

    text += '''\ndf_areas = df_areas.melt(id_vars=['acquisitionNumber', 'modifiedSequence', 'precursorCharge'],
                         value_vars=list(data_levels.keys()),
                         var_name='method', value_name='area')

df_areas['log2Area'] = df_areas['area'].apply(np.log2)

dats = dict()
for level in data_levels:
    dats[level] = df_areas[df_areas['method'] == level].pivot_table(index=['modifiedSequence', 'precursorCharge'],
                                                                    columns='acquisitionNumber', values='log2Area')
    dats[level] = dats[level].dropna()

multi_boxplot(dats, data_levels, dpi=250)\n\n```\n\n'''

    return text


def pc_metadata(meta_keys=None):

    text = f'''\n{python_block_header(stack()[0][3])}\n'''

    if meta_keys is None or len(meta_keys) == 0:
        text += '\nmeta_values = {}\n'
    else:
        text += f'\nmeta_values = {str(meta_keys)}\n'
        text += """\nMETADATA_QUERY = '''SELECT replicateId,
                         m.annotationKey as key,
                         m.annotationValue as value,
                         t.annotationType as type
                     FROM sampleMetadata m
                     LEFT JOIN sampleMetadataTypes t ON t.annotationKey == m.annotationKey
                     WHERE m.annotationKey IN ("{}")'''.format('", "'.join(meta_values.keys()))

# get metadata labels for pca plot
metadata = pd.read_sql(METADATA_QUERY, conn)
metadata = convert_string_cols(metadata)\n"""

    text += """\nmeta_values['acquisition_number'] = 'continuous'

def join_metadata(pc, acquired_ranks, meta_values, metadata=None):
    ''' join metadata to pc matrix '''

    if metadata is None:
        this_metadata = acquired_ranks
    else:
        this_metadata = acquired_ranks.join(metadata)
    pc = pc.join(this_metadata)

    for label_name, label_type in meta_values.items():
        if label_type == 'discrete':
            if any(pc[label_name].apply(pd.isna)):
                warnings.warn('Missing label values!', Warning)
            pc[label_name] = pc[label_name].apply(str)
        elif label_type == 'continuous':
            if any(pc[label_name].apply(pd.isna)):
                raise RuntimeError('Cannot have missing label values in continuous scale!')
        else:
            raise RuntimeError(f'"{label_type}" is an unknown label_type!')
    return pc
\n```\n\n"""

    return text


def pc_analysis(do_query, dpi, quant_col='totalAreaFragment', set_zero_to_min=True, have_color_vars=False):

    text = '\n{}\n'.format(python_block_header(f'{stack()[0][3]}_{quant_col}'))

    if do_query:
        text += f"""\nquery = '''SELECT
    p.replicateId,
    r.acquiredRank as acquisition_number,
    p.modifiedSequence,
    p.precursorCharge,
    p.{quant_col}
FROM PRECURSORS p
LEFT JOIN replicates r ON p.replicateId = r.replicateId
WHERE p.{quant_col} IS NOT NULL;'''

df_pc = pd.read_sql(query, conn)
"""
    else:
        text += f'''
df_pc = df[['replicateId', 'modifiedSequence', 'precursorCharge', '{quant_col}']]
df_pc['acquisition_number'] = df['acquiredRank']\n'''

    if set_zero_to_min:
        text += f'''\n# set zero areas to the minimum non-zero value
sele = df_pc['{quant_col}'].apply(lambda x: not np.isfinite(x) or x == 0)
df_pc.loc[sele, '{quant_col}'] = min(df_pc[df_pc['{quant_col}'] > 0]['{quant_col}'])\n'''

    text += f'''
df_pc['log2TotalAreaFragment'] = np.log2(df_pc['{quant_col}'])
# df_pc['zScore'] = df_pc.groupby('acquisition_number')['log2TotalAreaFragment'].transform(lambda x: np.abs(zscore(x)))

df_wide = df_pc.pivot_table(index=['modifiedSequence', 'precursorCharge'],
                            columns="replicateId", values='{quant_col}')
df_wide = df_wide.dropna()

# do pc analysis
pc, pc_var = pc_matrix(df_wide)

acquired_ranks = df_pc[['replicateId', 'acquisition_number']].drop_duplicates().set_index('replicateId')

pc = join_metadata(pc, acquired_ranks, meta_values, {'metadata=metadata' if have_color_vars else ''})\n'''

    text += f'''
for label_name, label_type in sorted(meta_values.items(), key=lambda x: x[0]):
    pca_plot(pc, label_name, pc_var, label_type=label_type, dpi={dpi})
```\n\n'''

    return text


def check_meta_keys_exist(conn, keys):
    '''
    Check that each metadata key exists in sampleMetadata table.
    '''
    all_good = True
    cur = conn.cursor()
    cur.execute('SELECT DISTINCT annotationKey FROM sampleMetadata;')
    metadata_keys = {x[0] for x in cur.fetchall()}

    for key in keys:
        if key not in metadata_keys:
            all_good = False
            LOGGER.error(f'Missing annotationKey: "{key}"')

    return all_good


def main():
    parser = argparse.ArgumentParser(description='Generate qmd report.')
    parser.add_argument('--dpi', default=DEFAULT_DPI, type=int,
                        help=f'Figure DPI in report. {DEFAULT_DPI} is the default.')
    parser.add_argument('-o', '--ofname', default=f'{DEFAULT_OFNAME}',
                        help=f'Output file basename. Default is "{DEFAULT_OFNAME}"')
    parser.add_argument('-a', '--addStdProtein', action='append', default=[],
                        help='Add standard protein name for retention time plot.')
    parser.add_argument('-c', '--addColorVar', action='append', default=None, dest='color_vars',
                        help='Add a annotationKey to color PCA plots.')
    parser.add_argument('--title', default=DEFAULT_TITLE,
                        help=f'Report title. Default is "{DEFAULT_TITLE}"')
    parser.add_argument('db', help='Path to precursor quality database.')
    args = parser.parse_args()

    if os.path.isfile(args.db):
        conn = sqlite3.connect(args.db)
    else:
        LOGGER.error(f"Database file: '{args.db}' does not exist!")
        sys.exit(1)

    # check that metadata annotatioKeys exist
    if args.color_vars:
        check_meta_keys_exist(conn, args.color_vars)

    # Check if metadata.is_normalized is set to True
    db_is_normalized = is_normalized(conn)
    conn.close()


    with open(args.ofname, 'w') as outF:
        outF.write(doc_header('{} {}'.format(os.path.abspath(__file__),
                              ' '.join(sys.argv[1:])), title=args.title))
        outF.write(add_header('Peptide independent metrics', level=1))
        outF.write(doc_initialize(args.db))

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

        outF.write(add_header('MS1 to MS2 ratio', level=2))
        outF.write(ms1_ms2_ratio(do_query=False, dpi=args.dpi))

        # Precursor areas section
        quant_cols = {'totalAreaFragment': 'Unnormalized'}
        if db_is_normalized:
            quant_cols['normalizedArea'] = 'Normalized'

        outF.write(add_header('Precursor areas', level=2))
        outF.write(precursor_areas(quant_cols))

        # Batch effects section
        outF.write(add_header('Batch effects', level=2))

        outF.write(pc_metadata(meta_keys=get_meta_key_types(args.db, args.color_vars)))

        outF.write(open_pannel_tabset())

        outF.write(add_header('Unnormalized', level=3))
        outF.write(pc_analysis(do_query=False, dpi=args.dpi,
                               have_color_vars=args.color_vars is not None))

        if db_is_normalized:
            outF.write(add_header('Normalized', level=3))
            outF.write(pc_analysis(do_query=True, dpi=args.dpi,
                                   quant_col='normalizedArea', set_zero_to_min=False,
                                   have_color_vars=args.color_vars is not None))

        outF.write(close_pannel_tabset())

        outF.write(doc_finalize())

if __name__ == '__main__':
    main()

