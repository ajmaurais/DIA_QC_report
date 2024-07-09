
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
        args = parse_data.parse_args(sys.argv[2:], prog='dia_qc parse')
        parse_data._main(args)


    def metadata_convert(self):
        args = metadata_convert.parse_args(sys.argv[2:], prog='dia_qc metadata_convert')
        metadata_convert._main(args)


    def normalize(self):
        args = normalize_db.parse_args(sys.argv[2:], prog='dia_qc normalize')
        normalize_db._main(args)


    def qc_qmd(self):
        args = generate_qc_qmd.parse_args(sys.argv[2:], prog='dia_qc qc_qmd')
        generate_qc_qmd._main(args)


    def batch_rmd(self):
        args = generate_batch_rmd.parse_args(sys.argv[2:], prog='dia_qc batch_rmd')
        generate_batch_rmd._main(args)


    def export_gene_matrix(self):
        args = export_gene_matrix.parse_args(sys.argv[2:], prog='dia_qc export_gene_matrix')
        export_gene_matrix._main(args)


    def db_export(self):
        args = export_tables.parse_args(sys.argv[2:], prog='dia_qc db_export')
        export_tables._main(args)


def main():
    _ = Main()


if __name__ == '__main__':
    main()
