
import argparse
from csv import DictReader

from .submodules.metadata import Dtype


def write_csv_row(elems, out):
    out.write('"{}"\n'.format('","'.join(elems)))


def get_sky_type(dtype):
    if dtype is Dtype.BOOL:
        return 'true_false'
    if dtype is Dtype.INT or dtype is Dtype.FLOAT:
        return 'number'
    return 'text'


def main():
    parser = argparse.ArgumentParser(description='Convert metadata annotation tsv to skyline '
                                                 'annotation csv and batch file to add annotation '
                                                 'definition to a skyline document.')
    parser.add_argument('annotation_tsv')

    args = parser.parse_args()

    with open(args.annotation_tsv, 'r') as inF:
        data = list(DictReader(inF, delimiter='\t'))

    annotation_headers = [x for x in data[0].keys() if x != 'Replicate']
    with open('sky_annotations.csv', 'w') as outF:
        # write header
        write_csv_row(['ElementLocator'] + [f'annotation_{x}' for x in annotation_headers], outF)

        for line in data:
            row = [f'Replicate:/{line["Replicate"]}']
            for header in annotation_headers:
                row.append(line[header])
            write_csv_row(row, outF)

    types = dict()
    for header in annotation_headers:
        types[header] = max(Dtype.infer_type(row[header]) for row in data)

    # write commands to add annotationd definitions to skyline file
    with open('sky_annotation_definitions.bat', 'w') as outF:
        for name, dtype in types.items():
            outF.write(f'--annotation-name="{name}" --annotation-targets=replicate')
            outF.write(f' --annotation-type={get_sky_type(dtype)}\n')


if __name__ == '__main__':
    main()
