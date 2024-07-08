
import sys
import argparse

from . import parse_data
from . import metadata_convert
from . import normalize_db
from . import generate_qc_qmd
from . import generate_batch_rmd
from . import export_gene_matrix
from . import export_tables

# SUBCOMMANDS = ['parse', 'db_export']


class Main():
    '''
    A class to parse subcommands.
    Inspired by this blog post: https://chase-seibert.github.io/blog/2014/03/21/python-multilevel-argparse.html
    '''


    def __init__(self):
        parser = argparse.ArgumentParser(description='Tools to generate QC and batch reports from DIA proteomics data',
                                         usage = f'''dia_qc <command> [<args>]

Available commands:
   parse                {parse_data.COMMAND_DESCRIPTION}
   metadata_convert     {metadata_convert.COMMAND_DESCRIPTION}
   normalize            {normalize_db.COMMAND_DESCRIPTION}
   qc_qmd               {generate_qc_qmd.COMMAND_DESCRIPTION}
   batch_rmd            {generate_batch_rmd.COMMAND_DESCRIPTION}
   export_gene_matrix   {export_gene_matrix.COMMAND_DESCRIPTION}
   db_export            {export_tables.COMMAND_DESCRIPTION}''')

        parser.add_argument('command', help='Subcommand to run.')

        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            sys.stderr.write(f"dia_qc: '{args.command}' is not a valid command!\n")
            sys.exit(1)

        getattr(self, args.command)()


    def parse(self):
        parse_data._main(parse_data.parse_args(sys.argv[2:]))


    def metadata_convert(self):
        metadata_convert._main(metadata_convert.parse_args(sys.argv[2:]))


    def normalize(self):
        normalize_db._main(normalize_db.parse_args(sys.argv[2:]))


    def qc_qmd(self):
        generate_qc_qmd._main(generate_qc_qmd.parse_args(sys.argv[2:]))


    def batch_rmd(self):
        generate_batch_rmd._main(generate_batch_rmd.parse_args(sys.argv[2:]))


    def export_gene_matrix(self):
        export_gene_matrix._main(export_gene_matrix.parse_args(sys.argv[2:]))


    def db_export(self):
        export_tables._main(export_tables.parse_args(sys.argv[2:]))


def main():
    _ = Main()


if __name__ == '__main__':
    main()
