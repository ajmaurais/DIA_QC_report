import sys
import os
import argparse
import re
import sqlite3
from inspect import stack

from .submodules.dtype import Dtype
from .submodules.logger import LOGGER
from .submodules.dia_db_utils import is_normalized
from .submodules.dia_db_utils import check_schema_version

COMMAND_DESCRIPTION = 'Generate QC qmd report.'

DEFAULT_OFNAME = 'qc_report.qmd'
DEFAULT_EXT = 'html'
DEFAULT_TITLE = 'DIA QC report'
DEFAULT_DPI = 250

METHOD_NAMES = ['unnormalized', 'normalized']
PRECURSOR_METHOD_NAMES = dict(zip(METHOD_NAMES, ['totalAreaFragment', 'normalizedArea']))
PROTEIN_METHOD_NAMES = dict(zip(METHOD_NAMES, ['abundance', 'normalizedAbundance']))
TABLE_TYPES = ('wide', 'long')

PEPTIDE_QUERY = '''
query = \'\'\'SELECT
    r.acquiredRank,
    p.modifiedSequence,
    p.precursorCharge
FROM precursors p
LEFT JOIN replicates r
    ON p.replicateId = r.id
WHERE p.totalAreaFragment > 0 AND r.includeRep == TRUE;\'\'\'

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
    p.rt * 60 as rt,
    p.maxFwhm * 60 as maxFwhm,
    p.averageMassErrorPPM as massError,
    p.libraryDotProduct,
    p.isotopeDotProduct
FROM precursors p
LEFT JOIN replicates r
    ON p.replicateId = r.id
WHERE r.includeRep == TRUE;\'\'\'

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
            body-width: 1200px
            sidebar-width: 0px
            margin-width: 200px
        self-contained: true
    pdf:
        fig-pos: 'H'
        extra_dependencies: ["float"]
        fig-format: 'png'
jupyter: python3
---\n\n''' % (command, title)
    return header


def figure_style():
    return '''<style>
img, figure {
    max-width: 100%;
    height: auto;
}

pre, code {
    white-space: pre-wrap;
    word-wrap: break-word;
}
</style>
'''


def add_header(text, level=1):
    return f'{"#" * level} {text}\n'


def open_pannel_tabset():
    return '\n::: { .panel-tabset }\n'


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
import os
import re
import sqlite3
import warnings
from statistics import stdev
from scipy.stats import zscore
import numpy as np
import pandas as pd

from DIA_QC_report.submodules.std_peptide_rt_plot import peptide_rt_plot
from DIA_QC_report.submodules.bar_chart import bar_chart
from DIA_QC_report.submodules.histogram import histogram, box_plot, multi_boxplot
from DIA_QC_report.submodules.pca_plot import calc_pc_matrix
from DIA_QC_report.submodules.pca_plot import mpl_pca_plot, plotly_pca_plot
from DIA_QC_report.submodules.dia_db_utils import read_wide_metadata

conn = sqlite3.connect('%s')
```\n\n''' % (python_block_header(stack()[0][3]), db_path)

    return text


def doc_finalize():
    text=f'\n{python_block_header(stack()[0][3])}\n\nconn.close()\n```\n\n'

    return text


# plots

def replicate_tic_areas(dpi=DEFAULT_DPI):
    text='''\n%s\n
# replicate tic bar chart
tic = pd.read_sql('SELECT acquiredRank, ticArea FROM replicates WHERE includeRep == TRUE', conn)
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
    text = '''\n%s\n%s
# replicate RT cor histogram
agg = df.groupby(['modifiedSequence', 'precursorCharge'])['rt'].agg('std')

if all(np.isnan(agg)):
    print('Can not plot histogram! 0 precursors found in all replicates.')
else:
    histogram(agg, 'Replicate RT SD (seconds)', dpi=%i,
              limits=(0, agg.quantile(0.95) * 3))
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
    columns="acquiredRank", values='%s')
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
                      columns="acquiredRank", values='ms2_ms1_ratio')
box_plot(data, 'Transition / precursor ratio',
         limits = (0, df['ms2_ms1_ratio'].quantile(0.95) * 3), dpi=%i)
```\n\n''' % (python_block_header(stack()[0][3]), PRECURSOR_QUERY if do_query else '\n', dpi)

    return text


def get_meta_key_types(db_path, keys):
    '''
    Determine whether each metadata key is discrete or continuous.
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
LEFT JOIN replicates r ON p.replicateId = r.id
WHERE r.includeRep == TRUE;'''\n"""

    data_levels = {v: i for i, v in enumerate(quant_cols.values())}
    text += f'\ndata_levels = {str(data_levels)}\n\n'

    text += '''df_areas = pd.read_sql(query, conn)
df_areas = df_areas.melt(id_vars=['acquisitionNumber', 'modifiedSequence', 'precursorCharge'],
                    value_vars=list(data_levels.keys()),
                    var_name='method', value_name='area')
df_areas = df_areas.dropna()
df_areas = df_areas.astype({'area': 'float64'})
df_areas['log2Area'] = np.log2(df_areas['area'] + 1)\n\n'''

    text += '''\ndats = dict()
for level in data_levels:
    dats[level] = df_areas[df_areas['method'] == level].pivot_table(index=['modifiedSequence', 'precursorCharge'],
                                                                    columns='acquisitionNumber', values='log2Area')

multi_boxplot(dats, data_levels, dpi=250)\n\n```\n\n'''

    return text


def pc_analysis(show_normalized=False, meta_keys=None):

    text = f'''\n{python_block_header(stack()[0][3])}\n'''

    if meta_keys is None or len(meta_keys) == 0:
        text += '\nmeta_values = {}\n'
    else:
        text += f"""\nmeta_values = {str(meta_keys)}"""

    text += """\nmeta_values['acquisition_number'] = 'continuous'

metadata = read_wide_metadata(conn, meta_vars=meta_values.keys())
metadata = metadata.rename(columns={'acquiredRank': 'acquisition_number'})

query = '''SELECT
    p.replicateId,
    r.replicate,
    p.modifiedSequence,
    p.precursorCharge,
    p.totalAreaFragment as Unnormalized"""

    if show_normalized:
        text += ',\n\tp.normalizedArea as Normalized'

    text += """\nFROM precursors p
LEFT JOIN replicates r ON p.replicateId = r.id
WHERE r.includeRep == TRUE;'''"""

    text += f'''\n\ndf_pc = pd.read_sql(query, conn)

pc_data = dict()
for col in ('Unnormalized', {"'Normalized'" if show_normalized else ''}):
    pc_data[col] = df_pc.pivot_table(index=['modifiedSequence', 'precursorCharge'],
                                     columns="replicateId", values=col).dropna()
    pc_data[col] = np.log2(pc_data[col] + 1)

skip_pca_plot = False
if max(len(data.index) for data in pc_data.values()) == 0:
    print('Can not perform PC analysis! 0 precursors found in all replicates.')
    skip_pca_plot = True

else:
    pc = dict()
    pc_var = dict()
    for col in pc_data:
        # do pc analysis
        if len(pc_data[col].index) > 0:
            pc[col] = calc_pc_matrix(pc_data[col])
        else:
            pc[col] = (None, None)
```\n\n'''

    return text


def normalize_var_names(meta_vars):
    ''' Normalize metadata variable name so it can be used in code block label '''

    fixed_vars = dict()
    new_vars_count = dict()
    for var in meta_vars:
        fixed_var = var

        # replace leading character that is not alpha
        fixed_var = re.sub(r'^[^A-Za-z_]+', '_', fixed_var)

        # replace non alpha numeric characters
        fixed_var = re.sub(r'[^A-Za-z0-9_]+', '_', fixed_var)

        fixed_vars[var] = fixed_var

        if fixed_var in new_vars_count:
            new_vars_count[fixed_var] += 1
        else:
            new_vars_count[fixed_var] = 1

    # Fix duplicated fixed var names
    for var, count in new_vars_count.items():
        if count > 1:
            i = 1
            for var_name in fixed_vars:
                if fixed_vars[var_name] == var:
                    fixed_vars[var_name] += f'_{i}'
                    i += 1

    assert len(fixed_vars) == len(set(fixed_vars.values()))

    return fixed_vars


def pca_tab_layout(meta_vars, plot_fxn):
    '''
    Layout PCA plots so each variable is in a seperate tab in the final report.

    A seperate code block is generated for each variable and the plot_fxn is called in each block.

    Parameters
    ----------
    meta_vars: list
        The list of metadata variables.
    plot_fxn: str
        The name of the PCA plotting function to use.

    Returns
    -------
    text: str
        The string of the formated pannel-tabset and code blocks.
    '''

    text = open_pannel_tabset()

    fixed_vars = normalize_var_names(meta_vars)

    for var, block_name in fixed_vars.items():
        text += f'\n{add_header(var, level=3)}'

        label = f'{block_name}_pca'
        text += f'\n{python_block_header(label)}\n\n'

        text += f'''if not skip_pca_plot:
    {plot_fxn}(pc, metadata, '{var}', label_type=meta_values['{var}'], add_title=False)\n```\n'''

    return text + close_pannel_tabset()


def pca_stacked_layout(plot_fxn, label=None):
    '''
    Layout PCA plots so they are stacked in the final report.

    The plot_fxn is called in a loop when the final report is rendered to
    generate a plot for each metadata variable.

    Parameters
    ----------
    plot_fxn: str
        The name of the PCA plotting function to use.
    lable: str
        The code block label. If None the function name is used. Default is None.

    Returns
    -------
    text: str
        The string of the formated code block.
    '''

    block_label = stack()[0][3] if label is None else label
    text = f'''\n{python_block_header(block_label)}

if not skip_pca_plot:
    for color_var, label_type in meta_values.items():
        {plot_fxn}(pc, metadata, color_var, label_type=label_type)
```\n\n'''

    return text


def pc_plots(meta_vars, pca_format):
    text = '\n::: {.content-visible when-format="html"}\n'

    if pca_format == 'tab':
        text += pca_tab_layout(meta_vars, 'plotly_pca_plot')
    elif pca_format == 'stack':
        text += pca_stacked_layout('plotly_pca_plot', label='html_pca_plots')

    text += '\n:::\n\n::: {.content-visible when-format="pdf"}\n'

    text += pca_stacked_layout('mpl_pca_plot', label='pdf_pca_plots')

    return text + ':::\n'


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


def check_std_proteins_exist(conn, proteins):
    '''
    Check that each standard protein exists in database.
    '''
    all_good = True

    cur = conn.cursor()
    cur.execute('SELECT name FROM proteins;')
    db_proteins = {x[0] for x in cur.fetchall()}

    for protein in proteins:
        if protein not in db_proteins:
            all_good = False
            LOGGER.error(f'Missing standard protein: "{protein}"')

    return all_good


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)
    parser.add_argument('--dpi', default=DEFAULT_DPI, type=int,
                        help=f'Figure DPI in report. {DEFAULT_DPI} is the default.')
    parser.add_argument('-o', '--ofname', default=f'{DEFAULT_OFNAME}',
                        help=f'Output file basename. Default is "{DEFAULT_OFNAME}"')
    parser.add_argument('-a', '--addStdProtein', action='append', default=None, dest='std_proteins',
                        help='Add standard protein name for retention time plot.')
    parser.add_argument('-c', '--addColorVar', action='append', default=None, dest='color_vars',
                        help='Add a annotationKey to color PCA plots.')
    parser.add_argument('--pcaFormat', choices=['tab', 'stack'], default='tab', dest='pca_format',
                        help='How to arrange PCA plots for each color variable in html report. '
                             'Either in seperate tabs or stacked. Plots are always stacked in pdf report.')
    parser.add_argument('--title', default=DEFAULT_TITLE,
                        help=f'Report title. Default is "{DEFAULT_TITLE}"')

    parser.add_argument('db', help='Path to sqlite qc database.')
    return parser.parse_args(argv)


def _main(args):
    '''
    Actual main method. `args` Should be initialized argparse namespace.
    '''

    # check args
    if os.path.isfile(args.db):
        conn = sqlite3.connect(args.db)
    else:
        LOGGER.error(f"Database file: '{args.db}' does not exist!")
        sys.exit(1)

    # check database version
    if not check_schema_version(conn):
        sys.exit(1)

    # check that metadata annotatioKeys exist
    if args.color_vars:
        if not check_meta_keys_exist(conn, args.color_vars):
            sys.exit(1)

    # check that metadata annotatioKeys exist
    if args.std_proteins:
        if not check_std_proteins_exist(conn, args.std_proteins):
            sys.exit(1)

    # Check if metadata.is_normalized is set to True
    db_is_normalized = is_normalized(conn)
    conn.close()

    with open(args.ofname, 'w') as outF:
        outF.write(doc_header('{} {}'.format(os.path.abspath(__file__),
                              ' '.join(sys.argv[1:])), title=args.title))
        outF.write(add_header('Peptide independent metrics', level=1))
        outF.write(doc_initialize(args.db))

        outF.write(figure_style())

        outF.write(add_header('Replicate TIC areas', level=2))
        outF.write(replicate_tic_areas(dpi=args.dpi))

        outF.write(add_header('Peptide dependent metrics\n', level=1))

        if args.std_proteins:
            outF.write(add_header('Standard retention times across replicates', level=2))
            outF.write(std_rt_dist(args.std_proteins, dpi=args.dpi))

        outF.write(add_header('Standard deviation of precursor RT across replicates', level=2))
        outF.write(rep_rt_sd(dpi=args.dpi, do_query=True))

        outF.write(add_header('Number of missed cleavages', level=2))
        outF.write(missed_cleavages(do_query=False, dpi=args.dpi))

        outF.write(add_header('Precursor charges', level=2))
        outF.write(precursor_charges(do_query=False, dpi=args.dpi))

        outF.write(add_header('Mass error distribution', level=2))
        outF.write(pivot_box_plot('massError', title='Mass error (ppm)',
                                  do_query=True, dpi=args.dpi))

        outF.write(add_header('Peak width distribution', level=2))
        outF.write(pivot_box_plot('maxFwhm', title='Peak FWHM (seconds)',
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

        meta_keys = get_meta_key_types(args.db, args.color_vars)
        outF.write(pc_analysis(meta_keys=meta_keys, show_normalized=db_is_normalized))

        outF.write(pc_plots(['acquisition_number'] + list(meta_keys.keys()), args.pca_format))

        outF.write(doc_finalize())


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc qc_qmd" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()

