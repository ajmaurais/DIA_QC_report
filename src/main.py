
import sys
import argparse

from . import parse_data
from . import skyline_report_convert
from . import metadata_convert
from . import validate_pipeline_params
from . import filter_replicates
from . import impute_missing
from . import normalize_db
from . import generate_qc_qmd
from . import generate_batch_rmd
from . import export_gene_matrix
from . import export_tables
from .submodules.logger import LOGGER, DEBUG

from . import __version__ as PACKAGE_VERSION

SUBCOMMANDS = {'parse', 'report_convert', 'metadata_convert',
               'filter', 'impute', 'normalize', 'validate',
               'qc_qmd', 'batch_rmd', 'export_gene_matrix', 'db_export'}


def _first_subcommand(argv):
    for i in range(1, len(argv)):
        if argv[i] in SUBCOMMANDS:
            return i
    return len(argv)


class VersionAction(argparse.Action):
    def __init__(self, option_strings,
                 dest=argparse.SUPPRESS, default=argparse.SUPPRESS,
                 help='Show program version and exit.', **kwargs):
        super().__init__(option_strings=option_strings, dest=dest, default=default,
                         nargs=0, help=help, **kwargs)
        self.version = PACKAGE_VERSION


    def __call__(self, parser, namespace, values, option_string=None):
        if self.version is None:
            sys.exit(1)
        sys.stdout.write(f'dia_qc version {self.version}\n')
        sys.exit(0)


class Main():
    '''
    A class to parse subcommands.
    Inspired by this blog post: https://chase-seibert.github.io/blog/2014/03/21/python-multilevel-argparse.html
    '''

    def __init__(self):
        parser = argparse.ArgumentParser(description='Tools to generate QC and batch reports from DIA proteomics data',
                                         add_help=False,
                                         usage = f'''dia_qc <command> [<args>]

Available commands:
   parse                {parse_data.COMMAND_DESCRIPTION}
   report_convert       {skyline_report_convert.COMMAND_DESCRIPTION}
   metadata_convert     {metadata_convert.COMMAND_DESCRIPTION}
   validate             {validate_pipeline_params.COMMAND_DESCRIPTION}
   filter               {filter_replicates.COMMAND_DESCRIPTION}
   impute               {impute_missing.COMMAND_DESCRIPTION}
   normalize            {normalize_db.COMMAND_DESCRIPTION}
   qc_qmd               {generate_qc_qmd.COMMAND_DESCRIPTION}
   batch_rmd            {generate_batch_rmd.COMMAND_DESCRIPTION}
   export_gene_matrix   {export_gene_matrix.COMMAND_DESCRIPTION}
   db_export            {export_tables.COMMAND_DESCRIPTION}''')

        parser.add_argument('-h', '--help', action='help',
                            help='Show this help message and exit.')
        parser.add_argument('-v', '--version', action=VersionAction)
        parser.add_argument(
            '--debug', action='store_true', default=False,
            help='Print function names and line numbers in log messages.'
        )
        parser.add_argument('command', help='Subcommand to run.')

        subcommand_start = _first_subcommand(sys.argv)
        args = parser.parse_args(sys.argv[1:(subcommand_start + 1)])
        argv = sys.argv[subcommand_start + 1:]

        # Set up logging
        if args.debug:
            LOGGER.setLevel(DEBUG)
            LOGGER.show_date = True
            LOGGER.set_debug(True)

        if not hasattr(self, args.command):
            sys.stderr.write(f"dia_qc: '{args.command}' is not a valid command!\n")
            sys.exit(2)

        getattr(self, args.command)(argv)


    def parse(self, argv):
        parse_data._main(argv, prog='dia_qc parse')


    def report_convert(self, argv):
        skyline_report_convert._main(argv, prog='dia_qc report_convert')


    def metadata_convert(self, argv):
        metadata_convert._main(argv, prog='dia_qc metadata_convert')


    def validate(self, argv):
        validate_pipeline_params._main(argv, prog='dia_qc validate')


    def filter(self, argv):
        filter_replicates._main(argv, prog='dia_qc filter')


    def impute(self, argv):
        impute_missing._main(argv, prog='dia_qc impute')


    def normalize(self, argv):
        normalize_db._main(argv, prog='dia_qc normalize')


    def qc_qmd(self, argv):
        generate_qc_qmd._main(argv, prog='dia_qc qc_qmd')


    def batch_rmd(self, argv):
        generate_batch_rmd._main(argv, prog='dia_qc batch_rmd')


    def export_gene_matrix(self, argv):
        export_gene_matrix._main(argv, prog='dia_qc export_gene_matrix')


    def db_export(self, argv):
        export_tables._main(argv, prog='dia_qc db_export')


def main():
    _ = Main()


if __name__ == '__main__':
    main()