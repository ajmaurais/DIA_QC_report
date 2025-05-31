
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
        parser.add_argument('command', help='Subcommand to run.')

        subcommand_start = _first_subcommand(sys.argv)
        args = parser.parse_args(sys.argv[1:(subcommand_start + 1)])
        argv = sys.argv[subcommand_start + 1:]

        if not hasattr(self, args.command):
            sys.stderr.write(f"dia_qc: '{args.command}' is not a valid command!\n")
            sys.exit(1)

        getattr(self, args.command)(argv)


    def parse(self, argv):
        args = parse_data.parse_args(argv, prog='dia_qc parse')
        parse_data._main(args)


    def report_convert(self, argv):
        args = skyline_report_convert.parse_args(argv, prog='dia_qc report_convert')
        skyline_report_convert._main(args)


    def metadata_convert(self, argv):
        args = metadata_convert.parse_args(argv, prog='dia_qc metadata_convert')
        metadata_convert._main(args)


    def validate(self, argv):
        args = validate_pipeline_params.parse_args(argv, prog='dia_qc validate')
        validate_pipeline_params._main(args)


    def filter(self, argv):
        args = filter_replicates.parse_args(argv, prog='dia_qc filter')
        filter_replicates._main(args)


    def impute(self, argv):
        args = impute_missing.parse_args(argv, prog='dia_qc impute')
        impute_missing._main(args)


    def normalize(self, argv):
        args = normalize_db.parse_args(argv, prog='dia_qc normalize')
        normalize_db._main(args)


    def qc_qmd(self, argv):
        args = generate_qc_qmd.parse_args(argv, prog='dia_qc qc_qmd')
        generate_qc_qmd._main(args)


    def batch_rmd(self, argv):
        args = generate_batch_rmd.parse_args(argv, prog='dia_qc batch_rmd')
        generate_batch_rmd._main(args)


    def export_gene_matrix(self, argv):
        args = export_gene_matrix.parse_args(argv, prog='dia_qc export_gene_matrix')
        export_gene_matrix._main(args)


    def db_export(self, argv):
        args = export_tables.parse_args(argv, prog='dia_qc db_export')
        export_tables._main(args)


def main():
    _ = Main()


if __name__ == '__main__':
    main()