
import sys
import argparse
from csv import DictReader

from .submodules.dtype import Dtype
from .submodules.read_metadata import Metadata

def write_csv_row(elems, out):
    out.write('"{}"\n'.format('","'.join(elems)))


def get_sky_type(dtype):
    if dtype is Dtype.BOOL:
        return 'true_false'
    if dtype is Dtype.INT or dtype is Dtype.FLOAT:
        return 'number'
    return 'text'


IN_FORMATS = ('json', 'csv', 'tsv')
OUT_FORMATS = ('json', 'skyline', 'csv', 'tsv')


def main():
    parser = argparse.ArgumentParser(description='Convert metadata annotation tsv to skyline '
                                                 'annotation csv and batch file to add annotation '
                                                 'definition to a skyline document.')
    parser.add_argument('-i', '--in', default=None,
                        choices=IN_FORMATS, dest='input_format',
                        help='Specify metadata file format. '
                             'By default the format is inferred from the file extension.')
    parser.add_argument('metadata')

    args = parser.parse_args()

    metadata_reader = Metadata()
    if not metadata_reader.read(args.metadata, args.input_format):
        sys.exit(1)

    data = metadata_reader.df.pivot(index='Replicate',
                                    columns='annotationKey',
                                    values='annotationValue').reset_index()

    # convert NULL columns to empty string
    for col, col_type in metadata_reader.types.items():
        if col_type is Dtype.NULL:
            data[col] = ''

    annotation_headers = [x for x in data.columns if x != 'Replicate']
    with open('sky_annotations.csv', 'w') as outF:
        # write header
        write_csv_row(['ElementLocator'] + [f'annotation_{x}' for x in annotation_headers], outF)

        for _, row in data.iterrows():
            line = [f'Replicate:/{row.Replicate}']
            for header in annotation_headers:
                line.append(row[header])
            write_csv_row(line, outF)

    # write commands to add annotationd definitions to skyline file
    with open('sky_annotation_definitions.bat', 'w') as outF:
        for name, dtype in metadata_reader.types.items():
            outF.write(f'--annotation-name="{name}" --annotation-targets=replicate')
            outF.write(f' --annotation-type={get_sky_type(dtype)}\n')


if __name__ == '__main__':
    main()
