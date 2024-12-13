
import sys
import os
import argparse
import sqlite3

from .submodules.logger import LOGGER
from .submodules.normalization import NORMALIZATION_METHODS
from .submodules.dia_db_utils import check_schema_version
from .submodules.dia_db_utils import is_normalized
from .submodules.dia_db_utils import validate_bit_mask, parse_bitmask_options

COMMAND_DESCRIPTION = 'Generate batch correction rmd report.'

DEFAULT_OFNAME = 'bc_report.rmd'
DEFAULT_EXT = 'html'
DEFAULT_TITLE = 'DIA Batch Correction Report'

# Column names for values at each stage of processing
PYTHON_METHOD_NAMES = ('unnormalized', 'normalized', 'batch_corrected')
PRECURSOR_METHOD_NAMES = ('totalAreaFragment', 'normalizedArea', 'normalizedArea.bc')
PROTEIN_METHOD_NAMES = ('abundance', 'normalizedAbundance', 'normalizedAbundance.bc')

# The method names used in plot headers
NORMALIZATION_METHOD_NAMES = {'DirectLFQ': 'DirectLFQ', 'median': 'Median'}

def doc_header(command, title=DEFAULT_TITLE):
    header='''<!-- This document was automatically generated. Command used to generate document: -->
<!-- %s -->\n
---
title: "%s"
input:
    html_document:
        toc: true
        toc_float: true
---\n\n\n''' % (command, title)
    return header


def add_header(text, level=1):
    return f'{"#" * level} {text}\n'


def r_block_header(label, echo=False, warning=False, message=False,
                   fig_width=None, fig_height=None):

    def bool_to_str(boo):
        return 'TRUE' if boo else 'FALSE'

    text = '```{' + f'r {label}, echo={bool_to_str(echo)}, warning={bool_to_str(warning)}, message={bool_to_str(message)}'

    if fig_width:
        text += f', fig.width={fig_width}'
    if fig_height:
        text += f', fig.height={fig_height}'

    return text + '}'


def ggsave(plot_path, plot_name, dim):
    '''
    Add ggsave line.

    plot_path: str
        Plot file path.
    plot_name: str
        The name of the plot variable.
    dim: tuple
        Plot dimensions. Tuple of 2 integers. (height, width)
    '''
    assert isinstance(dim, tuple) and len(dim) == 2 and all(isinstance(x, (int, float)) for x in dim)
    return f"ggsave('{plot_path}', {plot_name}, create.dir=T, height={dim[0]}, width={dim[1]})\n"


def skip_batch_correction(conn, batch1, batch2):
    # If no batch vars, test if we can use replicates.project
    if batch1 is None and batch2 is None:
        LOGGER.info('No batch variables specified. Attempting to use replicates.project as batch variable.')
        cur = conn.cursor()
        cur.execute('SELECT DISTINCT project FROM replicates WHERE includeRep == TRUE;')
        n_projects = len([x[0] for x in cur.fetchall()])
        if n_projects <= 1:
            LOGGER.warning('Only 1 project in replicates! Skipping batch correction.')
            return True
    return False


def norm_method_found(conn):
    '''
    Test whether protein_normalization_method and
    precursor_normalization_method keys exist in metadata table.
    '''

    cur = conn.cursor()
    cur.execute('SELECT * FROM metadata;')
    metadata = cur.fetchall()

    keys = ['protein_normalization_method', 'precursor_normalization_method']
    norm_methods = [value for key, value in metadata if key in keys]
    if len(norm_methods) != len(keys):
        LOGGER.error('Missing normalization methods in metadata table!')
        return False

    all_good = True
    for method in norm_methods:
        if method not in NORMALIZATION_METHODS:
            LOGGER.error(f"Unknown normalization method: '{method}'")
            all_good = False

    return all_good


def test_metadata_variables(conn, batch1=None, batch2=None,
                            covariate_vars=None, color_vars=None,
                            control_key=None, control_values=None):
    '''
    Test whether it is possible to build a valid query to generate the dat.metadata
    dataframe in the rmd report from the user specified parameters.

    Parameters
    ----------
    conn: sqlite.connection
        Initialized connection to precursor database.
    batch1: str
        The batch 1 annotationKey in sampleMetadata.
    batch2: str
        The batch 2 annotationKey in sampleMetadata.
    covariate_vars: list
        A list of sampleMetadata annotationKey(s) to use as covariate_vars.
    color_vars: list
        A list of sampleMetadata annotationKey(s) to color PCA plots.
    control_key: str
        A sampleMetadata annotationKey corresponding to whether a replicate is a control.
    control_values: list
        A list of sampleMetadata annotationValues for each control type.

    Returns
    -------
    all_good: boolean
        True if a valid query can be constructed, False if not.
    '''

    # get all metadata variable names
    cur = conn.cursor()
    cur.execute('SELECT annotationKey FROM sampleMetadata;')
    metadata_variables = set([x[0] for x in cur.fetchall()])

    def check_meta_key(color_cov, variable):
        if variable not in metadata_variables:
            LOGGER.error(f'Missing {color_cov} variable: "{variable}" in sampleMetadata table!')
            return False
        return True

    # Test that batch variables exist in metadata tables
    all_good = True
    if batch1:
        all_good = check_meta_key('batch1', batch1)
    if batch2:
        all_good = check_meta_key('batch2', batch2) and all_good

    # Test that covariate variables exist in metadata table
    all_good = (True if covariate_vars is None else all([check_meta_key('Covariate', v) for v in covariate_vars])) and all_good

    # Test that color variables exist in metadata table
    all_good = (True if color_vars is None else all([check_meta_key('Color', v) for v in color_vars])) and all_good

    # Test that control_key exists in metadata table
    if control_key is not None:
        if (control_key is not None) != (control_values is not None):
            LOGGER.error('Must specify both control_key and control_values')
            return False

        if check_meta_key('Control', control_key):
            cur = conn.cursor()
            cur.execute('SELECT annotationValue FROM sampleMetadata WHERE annotationKey == ?', (control_key,))

            db_control_values = set([x[0] for x in cur.fetchall()])

            for value in control_values:
                if value not in db_control_values:
                    LOGGER.error(f'Missing annotationValue "{value}" for annotationKey '
                                  '"{control_key}" in sampleMetadata table!')
                    all_good = False

        else:
            all_good = False

    return all_good


def pivot_df_longer(p):
    text = f'''
# pivot dat.{p} longer
dat.{p}.l <- dat.{p} %>%
    tidyr::pivot_longer(dplyr::all_of({p}.methods),
                        names_to='method', values_to='value') %>%
    dplyr::select(replicate, {p}, method, value)\n'''

    return text


def doc_initialize(db_path, batch1, batch2=None,
                   skip_bc=False, remove_missing=True,
                   color_vars=None, covariate_vars=None,
                   control_key=None, filter_metadata=False):
    '''
    Add initial block to load libraries, query batch db, construct metadata df,
    and setup precursor and protein dataframes.
    '''

    text = r_block_header('setup')

    text += f'''\n
library(rDIAUtils)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DBI)

# names for columns at different stages in analysis
precursor.methods <- c('{"', '".join(PRECURSOR_METHOD_NAMES[:2] if skip_bc else PRECURSOR_METHOD_NAMES)}')
protein.methods <- c('{"', '".join(PROTEIN_METHOD_NAMES[:2] if skip_bc else PROTEIN_METHOD_NAMES)}')
'''

    if batch1:
        text += f"\n# variables specified by --batch1\nbatch1 = '{batch1}'"
    if batch2:
        text += f"\n\n# variables specified by --batch2\nbatch2 = '{batch2}'"

    if color_vars:
        text += '''\n\n# variables specified by --addColor
color.vars <- c('{}')'''.format("', '".join(color_vars))

    if covariate_vars:
        text += '''\n\n# variables specified by --addCoveriate
color.vars <- c('{}')'''.format("', '".join(covariate_vars))

    if control_key:
        text += f"\n\n# control sample key specified by --controlKey\ncontrol.key = '{control_key}'"

    precursor_filter = ''
    protein_filter = ''
    if remove_missing:
        precursor_filter = ' AND normalizedArea IS NOT NULL'
        protein_filter = ' AND normalizedAbundance IS NOT NULL'

    text += f'''\n
# Load necissary tables from database
conn <- DBI::dbConnect(RSQLite::SQLite(), '{db_path}')
dat.precursor <- DBI::dbGetQuery(conn, 'SELECT
                                      p.replicateId,
                                      p.peptideId,
                                      p.modifiedSequence,
                                      p.precursorCharge,
                                      p.totalAreaFragment,
                                      p.normalizedArea
                                   FROM precursors p
                                   LEFT JOIN replicates r ON r.id == p.replicateId
                                   WHERE r.includeRep == TRUE{precursor_filter};')
peptideToProtein <- DBI::dbGetQuery(conn, 'SELECT
                                            prot.name as protein,
                                            p.modifiedSequence,
                                            p.replicateId
                                           FROM peptideToProtein ptp
                                           LEFT JOIN proteins prot ON prot.proteinId == ptp.proteinId
                                           LEFT JOIN (
                                            SELECT DISTINCT
                                            peptideId, replicateId, modifiedSequence
                                            FROM precursors
                                           ) p ON p.peptideId == ptp.peptideID')
dat.protein <- DBI::dbGetQuery(conn, 'SELECT
                                    q.replicateId,
                                    p.name as protein,
                                    q.abundance,
                                    q.normalizedAbundance
                                FROM proteinQuants q
                                LEFT JOIN proteins p ON p.proteinId = q.proteinId
                                LEFT JOIN replicates r ON r.id == q.replicateId
                                WHERE r.includeRep == TRUE{protein_filter};')
norm.methods <- DBI::dbGetQuery(conn, "SELECT * FROM metadata
                               WHERE key IN ('precursor_normalization_method', 'protein_normalization_method')")
dat.metadata <- rDIAUtils::readWideMetadata(conn)
DBI::dbDisconnect(conn)\n'''

    text += '''
# Format normalization method strings
norm.methods$value <- factor(norm.methods$value, levels=c('median', 'DirectLFQ'),
                            labels=c('Median', 'DirectLFQ'))
norm.methods <- setNames(norm.methods$value, sub('_normalization_method$', '', norm.methods$key))

# fix special characters in replicate names so they can be R headers
dat.metadata$replicate <- make.names(dat.metadata$replicate)
dat.precursor <- dplyr::left_join(dat.precursor,
                                  dplyr::select(dat.metadata, replicate, replicateId),
                                  by='replicateId')
dat.protein <- dplyr::left_join(dat.protein,
                                dplyr::select(dat.metadata, replicate, replicateId),
                                by='replicateId')
peptideToProtein <- dplyr::left_join(peptideToProtein,
                                     dplyr::select(dat.metadata, replicate, replicateId),
                                     by='replicateId') %>%
    dplyr::select(-replicateId)

dat.metadata[[batch1]] <- factor(dat.metadata[[batch1]])\n'''

    if batch2:
        text += 'dat.metadata[[batch2]] <- factor(dat.metadata[[batch2]])\n'

    text += '''
# combine sequence and precursorCharge cols
dat.precursor$precursor <- with(dat.precursor, paste(modifiedSequence, precursorCharge, sep='_'))

# log2 transform abundances
dat.precursor$totalAreaFragment <- log2(dat.precursor$totalAreaFragment + 1)
dat.precursor$normalizedArea <- log2(dat.precursor$normalizedArea + 1)
dat.protein$abundance <- log2(dat.protein$abundance + 1)
dat.protein$normalizedAbundance <- log2(dat.protein$normalizedAbundance + 1)\n'''

    if skip_bc:
        text += f"\n{pivot_df_longer('precursor')}\n{pivot_df_longer('protein')}\n"

    text += '```\n\n\n'

    return text


def batch_correction(p, batch1, batch2=None, covariate_vars=None, bc_method='combat'):
    '''
    Add batch correction section for precursor or protein
    '''

    meta_cols = [s for v, s in zip([batch1, batch2, covariate_vars],
                                   ['batch1', 'batch2', 'covariate.vars']) if v is not None]

    a = 'Area' if p != 'protein' else 'Abundance'

    text = r_block_header(f'{p}_batch_correction')
    text += f'''\n
# batch correction
dat.{p}.bc <- dat.{p} %>%
    dplyr::left_join(dplyr::select(dat.metadata, replicateId, all_of(c({', '.join(meta_cols)}))),
                     by='replicateId') %>%
    rDIAUtils::batchCorrection('normalized{a}', batch1=batch1, {'batch2=batch2,' if batch2 else ''}
                               {'covariate.cols=covariate.vars, ' if covariate_vars else ''}bc.method='{bc_method}', rowsName='{p}')

# Join batch corrected areas to initial df
dat.{p} <- dat.{p} %>%
    dplyr::left_join(dplyr::select(dat.{p}.bc, replicate, {p},
                     normalized{a}.bc=normalized{a}),
                     by=c('replicate', '{p}'))

    {pivot_df_longer(p)}```\n\n\n'''

    return text


def precursor_norm_plot(plot_file_path=None, skip_bc=False):
    '''
    Generate precursor normalization box plot.
    '''

    fig_width = 8
    fig_height = 10

    text = '\n' + r_block_header('area_dist_plot', fig_width=fig_width, fig_height=fig_height)
    text += f'''\n
dat.p <- dat.precursor.l %>%
    dplyr::left_join(dplyr::select(dat.metadata, replicate, acquiredRank, all_of(batch1)))
dat.p$method <- factor(dat.p$method, levels=precursor.methods,
                       labels=c('Unnormalized', 'Median Normalized'{'' if skip_bc else ", 'Normalized, Batch Corrected'"}))
dat.p$acquiredRank <- as.numeric(as.factor(rank(-dat.p$acquiredRank)))

p.norm <- ggplot(dat.p, aes(x=acquiredRank, y=value, group=acquiredRank{'' if skip_bc else ', color=get(batch1)'})) +
    facet_wrap(~method, nrow = 1, scales = 'free_x', strip.position = 'top') +
    geom_boxplot(outlier.size = 0.5) +'''

    if not skip_bc:
        text += '''
    scale_color_discrete(name='Batch') +'''

    text += '''
    coord_flip() +
    ylab('log2(Area)') +
    xlab('Acquisition number') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'top')\n'''

    if plot_file_path:
        text += ggsave(plot_file_path, 'p.norm', (fig_height, fig_width))

    text += '\np.norm\n```\n\n\n'
    return text


def cv_plot(control_values=None, plot_file_path=None, skip_bc=False):

    fig_width = 6.5
    fig_height = 1 + 3 * len(control_values) if control_values else 4
    text = '\n' + r_block_header('cv_dist_plot', fig_width=fig_width, fig_height=fig_height) + '\n'

    # additional text to add if faceting by control samples
    filter_text = ''
    group_by_text = ''
    facet_text = ''

    if control_values:
        text += '''\n# control sample values specified by --controlValues
control.values <- c('{}')'''.format("', '".join(control_values))
        text += '''
control.reps <- dat.metadata[dat.metadata[[control.key]] %in% control.values,] %>%
    dplyr::select(replicate, all_of(control.key))
control.reps[[control.key]] <- factor(control.reps[[control.key]])\n'''

        filter_text = '''dplyr::filter(replicate %in% control.reps$replicate) %>%
    dplyr::left_join(control.reps, by='replicate') %>%\n'''
        group_by_text = ', !!rlang::sym(control.key)'
        facet_text = "facet_wrap(as.formula(paste0('~ ', '`', control.key, '`')), ncol=1) + \n\t"

    text += f'''
# calculate CV for each precursor
dat.cv <- dat.precursor.l %>%{filter_text}
    dplyr::mutate(value=2^value) %>% # Convert back to linear space for CV calculation
    dplyr::group_by(precursor, method{group_by_text}) %>%
    dplyr::summarize(cv = (sd(value) / mean(value)) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::filter(cv < 100)

# Convert method to factor
dat.cv$method <- factor(dat.cv$method, levels=precursor.methods,
                       labels=c('Unnormalized', 'Median Normalized'{'' if skip_bc else ", 'Batch Corrected'"}))

# CV distribution plot
p.cv <- ggplot(dat.cv, aes(x=cv, fill=method, color=method)) +
    {facet_text}geom_density(alpha=0.2) +
    ylab('Density') +
    xlab('CV (%)') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title=element_blank(),
          legend.direction='horizontal',
          legend.position = 'top')\n'''

    if plot_file_path:
        text += ggsave(plot_file_path, 'p.cv', (fig_height, fig_width))

    text += '\np.cv\n```\n\n\n'

    return text


def pca_plot(p, color_vars=None, plot_file_path=None, skip_bc=False):

    n_methods = 2 if skip_bc else 3

    fig_width = n_methods * 3 + 2
    fig_height = ((len(color_vars) if color_vars else 0) + (1 if skip_bc else 2)) * 3 + 0.4

    text = '\n' + r_block_header(f'{p}_pca', fig_width=fig_width, fig_height=fig_height)

    color_vars_text = '' if color_vars is None else ', color.vars'
    batch_var_text = '' if skip_bc else ", 'Batch'=batch1"

    text += f'\n\n# {p} PCAs\npcs.{p} <- list()\nfor(method in {p}.methods)'
    text += '{' + f'''
    pcs.{p}[[method]] <- rDIAUtils::pcAnalysis(dat.{p}.l[dat.{p}.l$method == method,],
                                                     'value', rowsName='{p}', columnsName='replicate')'''
    text += '\n}\n\n'
    text += f'''names({p}.methods) <- c('Unnormalized', paste(norm.methods['{p}'], 'Normalized'){'' if skip_bc else ", 'Batch corrected'"})
pca.{p} <- rDIAUtils::arrangePlots(pcs.{p}, row.cols={p}.methods,
    color.cols=c('Acquisition\\nnumber'='acquiredRank'{batch_var_text}{color_vars_text}),
    dat.metadata=dat.metadata)\n'''

    if plot_file_path:
        text += ggsave(plot_file_path, f'pca.{p}', (fig_height, fig_width))

    text += f'\npca.{p}\n```\n\n\n'

    return text


def wide_table(p, py_name, r_name, df_name):
    row_names = ['protein']
    if p == 'precursor':
        row_names.append('modifiedSequence')
        row_names.append('precursorCharge')

    text = f'''\nmessage('Writing {p}s {py_name} to: "{p}s_{py_name}_wide.tsv"')
write.table(rDIAUtils::pivotLonger({df_name}, valuesFrom='{r_name}',
                                   rowsName=c('{"', '".join(row_names)}'),
                                   columnsName='replicate'),
            file='{p}s_{py_name}_wide.tsv', sep='\\t', row.names=F, quote=F)\n'''

    return text


def long_table(p, df_name, unnormalized=True, normalized=True, batch_corrected=True):
    row_identifers = ['replicate', 'protein']
    quant_value_index = [i for i, value in zip((0, 1, 2), (unnormalized, normalized, batch_corrected)) if value]
    if p == 'precursor':
        row_identifers += ['modifiedSequence', 'precursorCharge']
        quant_values = [PRECURSOR_METHOD_NAMES[i] for i in quant_value_index]
    else:
        quant_values = [PROTEIN_METHOD_NAMES[i] for i in quant_value_index]

    text = f'''\nmessage('Writing {p}s long to: "{p}s_long.tsv"')
write.table(dplyr::select({df_name}, {', '.join(row_identifers)},
                          {', '.join(quant_values)}),
            file='{p}s_long.tsv', sep='\\t', row.names=F, quote=F)\n'''

    return text


def write_tables_section(precursor_tables, protein_tables, metadata_tables):

    text = add_header('TSV files generated:', level=1)

    text += '\n' + r_block_header('write_tables', message=True)
    text += '\n\n'

    if any(t for d in precursor_tables.values() for t in d.values()):
        text += '''dat.precursor.j <- dplyr::inner_join(peptideToProtein, dat.precursor,
                                by=c('replicate', 'modifiedSequence'), relationship='many-to-many')\n'''

    # precursor wide tables
    for (py_name, write), r_name in zip(precursor_tables['wide'].items(), PRECURSOR_METHOD_NAMES):
        if write:
            text += wide_table('precursor', py_name, r_name, 'dat.precursor.j')

    # precursor long tables
    if any(precursor_tables['long'].values()):
        text += long_table('precursor', 'dat.precursor.j', **precursor_tables['long'])

    # protein wide tables
    for (py_name, write), r_name in zip(protein_tables['wide'].items(), PROTEIN_METHOD_NAMES):
        if write:
            text += wide_table('protein', py_name, r_name, 'dat.protein')

    # protein long tables
    if any(protein_tables['long'].values()):
        text += long_table('protein', 'dat.protein', **protein_tables['long'])

    # metadata tables
    if metadata_tables['long']['write']:
        text += '''
message('Writing metadata long to "metadata_long.tsv"')
dat.metadata.l <- dplyr::left_join(dat.meta.l, dplyr::select(dat.metadata, replicateId, replicate),
                                   by='replicateId') %>%
    dplyr::select(replicate, annotationKey, annotationValue, annotationType)
write.table(dat.metadata.l, file='metadata_long.tsv', sep='\\t', row.names=F, quote=F)\n'''

    if metadata_tables['wide']['write']:
        text += '''
message('Writing metadata wide to "metadata_wide.tsv"')
write.table(dplyr::select(dat.metadata, -replicateId),
            file='metadata_wide.tsv', sep='\\t', row.names=F, quote=F)\n'''

    text += '\n```\n\n'

    return text


def remove_bc_tables(table, name=None):
    for direction in ('wide', 'long'):
        if table[direction]['batch_corrected']:
            LOGGER.warning(f"{name + ' b' if name else 'B'}atch corrected {direction} table not available when batch correction is skipped!")
            table[direction]['batch_corrected'] = False

    return table


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)

    file_settings = parser.add_argument_group('R markdown file settings')
    file_settings.add_argument('-o', '--ofname', default=f'{DEFAULT_OFNAME}',
                               help=f'Output file basename. Default is "{DEFAULT_OFNAME}"')
    file_settings.add_argument('--title', default=DEFAULT_TITLE,
                               help=f'Report title. Default is "{DEFAULT_TITLE}"')
    file_settings.add_argument('-f', '--filterMetadata', default=False,
                               action='store_true', dest='filter_metadata',
                               help='Filter metadata table to only include batch, color, and '
                               'covariate variables.')
    file_settings.add_argument('-m', '--bcMethod', choices=['limma', 'combat'], default='combat',
                               help='Batch correction method. Default is "combat".')
    file_settings.add_argument('--allPrecursors', default=False, action='store_true', dest='all_precursors',
                               help="Don't remove precursors and proteins with missing values. "
                                    "Setting this option could cause an error when the rmd renders.")
    file_settings.add_argument('--savePlots', default=None, dest='plot_ext',
                               help='Save all plots to file with specified extension.')

    batch_vars = parser.add_argument_group('Batch variables',
                     'Batch variables are optional. By default the project column in the replicates '
                     'table is used as the batch variable. The project column corresponds to the '
                     'Skyline document which each sample was in. There must be at least 2 different '
                     'values in the project column to it as the batch variable. '
                     'If batch variables are specified, they must exist in the sampleMetadata table.')

    batch_vars.add_argument('-1', '--batch1', default=None,
                            help='sampleMetadata variable to use for batch 1.')
    batch_vars.add_argument('-2', '--batch2', default=None,
                            help='sampleMetadata variable to use for batch 2. Using 2 batch variables '
                            'is only available with limma as the batch correction method. '
                            '--batch1 must be also be specified to use --batch2.')

    color_vars = parser.add_argument_group('Color and covariate variables',
                    'These options can be specified multiple times to add multiple variables.')
    color_vars.add_argument('--addCovariate', action='append', dest='covariate_vars',
                            help='Add a sampleMetadata annotationKey to use as a covariate for '
                                 'batch correction.')
    color_vars.add_argument('-c', '--addColor', action='append', dest='color_vars',
                            help='Add a sampleMetadata annotationKey to use to color PCA plots.')

    control_vars = parser.add_argument_group('Control sample variables',
                        'Specify a metadata variable to indicate whether a replicate is a control. '
                        'These key, value pairs will be used to separate the panels of the control CV plot. '
                        'If --controlKey is specified, at least one value must be specified with --addControlValue')

    control_vars.add_argument('--controlKey', default=None, dest='control_key',
                              help='sampleMetadata annotationKey that has variable indication whether a '
                                   'replicate is a control.')
    control_vars.add_argument('--addControlValue', action='append', dest='control_values',
                              help='Add sampleMetadata annotationValue(s) which indicate whether '
                                   'a replicate is a control.')

    parser.add_argument('-s', '--skipTests', default=False, action='store_true',
                        help='Skip tests that batch, color, and covariate variables exist in '
                             'sampleMetadata table?')

    table_args = parser.add_argument_group('Output tables',
                     'The tsv files to write are specified by a 2 digit bit mask. '
                     'The first digit is for the wide formatted report, and the second digit is for '
                     'the long formatted report. An integer between 0 and 4 is assigned for each stage in '
                     'the normalization process. 0 for no report, 1 is for unnormalized, 2 is for '
                     'normalized, and 4 is for batch corrected.')

    # Some of these examples might be good to include in the help, but is is already a lot of text.
    # For example, a mask of 42 would produce a wide batch corrected and long normalized tsv file.
    # A mask of 30 (1+2=3) would produce both a unnormalized and normalized wide tsv file.
    # A mask of 70 (1+2+4=7) would produce a unnormalized, normalized, and batch
    # corrected wide tsv file.

    table_args.add_argument('-p', '--precursorTables', default='40',
                        help='Tables to write for precursors. 40 is the default')
    table_args.add_argument('-r', '--proteinTables', default='40',
                            help='Tables to write for proteins. 40 is the default')
    table_args.add_argument('--metadataTables', default='10',
                            help='Tables to write for metadata. Only 0 or 1 are supported. '
                                 '0 for false, 1 for true. 10 is the default')

    parser.add_argument('db', help='Path to sqlite batch database.')
    return parser.parse_args(argv)


def _main(args):
    '''
    Actual main method. `args` Should be initialized argparse namespace.
    '''

    # string variable args
    string_args = ['batch1', 'batch2', 'covariate_vars', 'color_vars', 'control_key']
    string_args = {x: getattr(args, x) for x in string_args}

    # check table_args
    if not validate_bit_mask(args.precursorTables, 3, 2):
        LOGGER.error('Error parsing --precursorTables')
        sys.exit(1)
    if not validate_bit_mask(args.proteinTables, 3, 2):
        LOGGER.error('Error parsing --proteinTables')
        sys.exit(1)
    if not validate_bit_mask(args.metadataTables, 1, 2):
        LOGGER.error('Error parsing --metadataTables')

    # parse table_args
    protein_tables = parse_bitmask_options(args.proteinTables, ('wide', 'long'), PYTHON_METHOD_NAMES)
    precursor_tables = parse_bitmask_options(args.precursorTables, ('wide', 'long'), PYTHON_METHOD_NAMES)
    metadata_tables = parse_bitmask_options(args.metadataTables, ('wide', 'long'), ('write',))

    try:
        # Initialize db connection
        conn = None
        if os.path.isfile(args.db):
            conn = sqlite3.connect(args.db)
        else:
            raise FileNotFoundError(f'Database file ({args.db}) does not exist!')

        # Determine if batch correction should be skipped
        skip_bc = skip_batch_correction(conn, args.batch1, args.batch2)

        if not args.skipTests:
            # Check database version
            if not check_schema_version(conn):
                raise RuntimeError

            if not is_normalized(conn):
                raise RuntimeError('Database file it not normalized!')

            if not norm_method_found(conn):
                raise RuntimeError

            # Check control_vars args
            if args.control_key is None and args.control_values is not None:
                raise RuntimeError('No control key specified!')
            if args.control_key is not None and args.control_values is None:
                raise RuntimeError('No control value(s) specified!')

            if not test_metadata_variables(conn, **string_args,
                                           control_values=args.control_values):
                raise RuntimeError

    except (RuntimeError, FileNotFoundError) as e:
        message = str(e)
        if message != '':
            LOGGER.error(message)
        sys.exit(1)

    finally:
        if conn is not None:
            conn.close()

    # make sure batch corrected tables aren't written when batch correction is skipped.
    if skip_bc:
        protein_tables = remove_bc_tables(protein_tables, name='Protein')
        precursor_tables = remove_bc_tables(precursor_tables, name='Precursor')

    # set batch1 to default if it is None
    if string_args['batch1'] is None:
        string_args['batch1'] = 'project'

    # generate rmd
    with open(args.ofname, 'w') as outF:
        outF.write(doc_header('{} {}'.format(os.path.basename(sys.argv[0]), ' '.join(sys.argv[1:])),
                              title=args.title))

        # setup
        outF.write(doc_initialize(args.db, **string_args,
                                  skip_bc=skip_bc,
                                  filter_metadata=args.filter_metadata,
                                  remove_missing=not args.all_precursors))

        if not skip_bc:
            # precursor batch correction
            outF.write(batch_correction('precursor', string_args['batch1'],
                                        batch2=args.batch2,
                                        covariate_vars=args.covariate_vars,
                                        bc_method=args.bcMethod))

            # protein batch correction
            outF.write(batch_correction('protein', string_args['batch1'],
                                        batch2=args.batch2,
                                        covariate_vars=args.covariate_vars,
                                        bc_method=args.bcMethod))

        # precursor normalization plot
        outF.write(add_header('Precursor normalization', level=1))
        outF.write(precursor_norm_plot(plot_file_path='plots/precursor_normalization.tiff' if args.plot_ext else None,
                                       skip_bc=skip_bc))

        # CV distribution plot
        outF.write(add_header(f"{'Control ' if args.control_values else ''}CV distribution", level=1))
        cv_plot_file_path = None
        if args.plot_ext:
            cv_plot_file_path = f'plots/{"control_" if args.control_values else ""}cv_dist.{args.plot_ext}'
        outF.write(cv_plot(control_values=args.control_values,
                           plot_file_path=cv_plot_file_path,
                           skip_bc=skip_bc))

        # precursor PCA plot
        outF.write(add_header('Precursor batch correction PCA', level=1))
        outF.write(pca_plot('precursor',
                            color_vars=args.color_vars,
                            skip_bc=skip_bc,
                            plot_file_path=f'plots/precursor_pca.{args.plot_ext}' if args.plot_ext else None))

        # protein PCA plot
        outF.write(add_header('Protein batch correction PCA', level=1))
        outF.write(pca_plot('protein',
                            color_vars=args.color_vars,
                            skip_bc=skip_bc,
                            plot_file_path=f'plots/protein_pca.{args.plot_ext}' if args.plot_ext else None))

        # Optional output tables
        # Check if there is at least 1 table to be written.
        if sum(int(arg) for arg in (args.proteinTables, args.precursorTables, args.metadataTables)) > 0:
            outF.write(write_tables_section(precursor_tables, protein_tables, metadata_tables))


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc batch_rmd" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()

