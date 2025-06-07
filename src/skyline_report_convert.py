
import sys
from os.path import splitext, basename
import argparse
import pandas as pd

from .submodules.skyline_reports import SkylineReport, ReplicateReport, PrecursorReport
from .submodules.skyline_reports import detect_delim
from .submodules.logger import LOGGER

COMMAND_DESCRIPTION = 'Convert Skyline report(s) to parquet format.'


def parse_args(argv, prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=COMMAND_DESCRIPTION)
    parser.add_argument('-r', '--remove-unknown-cols', default=False, action='store_true',
                        help='Remove unknown columns from the report before writing to parquet.')
    parser.add_argument('reports', nargs='+', help='Path(s) to Skylime report(s) to convert.')

    return parser.parse_args(argv)


def _main(argv, prog=None):
    ''' Actual main method. '''

    args = parse_args(argv, prog=prog)
    report_parsers = [PrecursorReport(), ReplicateReport()]

    for report_path in args.reports:
        with open(report_path, 'r') as report_file:
            headers = SkylineReport.read_headers(report_file)

        report_written = False
        for parser in report_parsers:
            if (language := parser.detect_language(headers)):
                if parser.check_df_columns(headers, language=language, quiet=True):
                    LOGGER.info('Found %s report: %s', parser.report_name, report_path)
                    df = parser.read_report(report_path, return_invariant=True,
                                            remove_unknown_cols=args.remove_unknown_cols)
                    parquet_path = f'{splitext(basename(report_path))[0]}.parquet'
                    LOGGER.info('Writing %s report to %s', parser.report_name, parquet_path)
                    df.to_parquet(parquet_path, index=False)

                    report_written = True
                    break

        if not report_written:
            LOGGER.warning('Reading %s as a generic Skyline report', report_path)
            with open(report_path, 'r') as report_file:
                df = pd.read_csv(report_path, sep=detect_delim(report_file))
            parquet_path = f'{splitext(basename(report_path))[0]}.parquet'
            LOGGER.info('Writing generic Skyline report to %s', parquet_path)
            df.to_parquet(parquet_path, index=False)


def main():
    LOGGER.warning('Calling this script directly is deprecated. Use "dia_qc report_convert" instead.')
    _main(parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()
