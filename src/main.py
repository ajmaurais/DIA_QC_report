
import sys
import argparse

from . import parse_data
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
   parse        {parse_data.COMMAND_DESCRIPTION}
   db_export    {export_tables.COMMAND_DESCRIPTION}''')

        parser.add_argument('command', help='Subcommand to run.')

        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            sys.stderr.write(f"dia_qc: '{args.command}' is not a valid command!\n")
            sys.exit(1)

        getattr(self, args.command)()


    def parse(self):
        parse_data._main(parse_data.parse_args(sys.argv[2:]))


    def db_export(self):
        export_tables._main(export_tables.parse_args(sys.argv[2:]))


def main():
    _ = Main()


if __name__ == '__main__':
    main()
