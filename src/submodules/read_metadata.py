
import re
import json
from csv import reader as csv_reader, DictReader
from os.path import splitext

from jsonschema import validate, ValidationError
import pandas as pd

from .dtype import Dtype
from .logger import LOGGER


# sample metadata json schema
JSON_SCHEMA = {
    'type': 'object',
    'additionalProperties': {
        'type': 'object',
        'additionalProperties': {
            'oneOf': [ {'type': 'null'},
                       {'type': 'boolean'},
                       {'type': 'number'},
                       {'type': 'string'} ]
        },
        'minProperties': 1
    },
    'minProperties': 1
}


class Metadata():
    def __init__(self):
        self.df = None
        self.types = None
        self.input_format = None
        self.input_file = None


    @staticmethod
    def _metadata_data_to_df(data):
        df = pd.DataFrame(data)
        df['Replicate'] = df['Replicate'].apply(lambda x: splitext(x)[0])

        # pivot longer
        df = pd.melt(df, id_vars=['Replicate'], var_name='annotationKey', value_name='annotationValue',
                     value_vars=[x for x in list(df.columns) if x != 'Replicate'])
        df['annotationValue'] = df['annotationValue'].astype(str)
        return df


    @staticmethod
    def _infer_types(df):
        ''' Determine annotationValue types '''
        types = dict()
        for row in df.itertuples():
            if row.annotationKey in types:
                types[row.annotationKey] = max(types[row.annotationKey],
                                               Dtype.infer_type(row.annotationValue))
            else:
                types[row.annotationKey] = Dtype.infer_type(row.annotationValue)
        return types


    def _read_metadata_skyline_csv(self, rows, exclude_null=True):
        self.input_format = 'skyline'

        headers = ['Replicate'] + [re.sub(r'^annotation_', '', h) for h in rows[0][1:]]
        rows = rows[1:]
        annotation_objects = set()

        # iterate through rows and remove all annotation objects that are not Replicate
        new_rows = list()
        for row in rows:
            annotation_match = re.search(f'^([A-Z][A-Za-z]+):/', row[0])
            if annotation_match is None:
                LOGGER.warning(f'Found unknown annotation object: {row[0]}')
                continue
            annotation_object = annotation_match.group(1)
            annotation_objects.add(annotation_object)

            if annotation_object == 'Replicate':
                row[0] = re.sub('^Replicate:/', '', row[0])
                new_rows.append(row)

        if 'Replicate' not in annotation_objects:
            LOGGER.warning('No Replicate annotations in csv!')
            return False

        for obj in annotation_objects:
            if obj != 'Replicate':
                LOGGER.warning(f'Found annotations for: {obj} in csv!')

        data = [dict(zip(headers, row)) for row in new_rows]
        self.df = self._metadata_data_to_df(data)
        self.types = self._infer_types(self.df)

        # set values of Dtype.NULL to NA
        self.df['annotationValue'] = self.df.apply(lambda x: pd.NA if self.types[x.annotationKey] is Dtype.NULL
                                                                   else x.annotationValue, axis=1)

        if exclude_null:
            self.df = self.df[~pd.isna(self.df['annotationValue'])].reset_index()
            self.types = {k: v for k, v in self.types.items() if v is not Dtype.NULL}

        # This is an annoying edge case with metadata annotations exported from Skyline
        # Boolean annotations are either True or blank instead of False
        # Therefore we will set all blank boolean annotations to False
        for key, var_type in self.types.items():
            if var_type is Dtype.BOOL:
                sele = self.df['annotationKey'] == key
                self.df.loc[sele, 'annotationValue'] = self.df[sele]['annotationValue'].apply(lambda x: 'False' if x == '' else x)

        return True


    def _read_metadata_json(self, fname):
        with open(fname, 'r') as inF:
            data = json.load(inF)
            try:
                validate(data, JSON_SCHEMA)
            except ValidationError as e:
                LOGGER.error(f'Invalid metadata format:\n{e.message}')
                return False

        # determine metadata types
        self.types = dict()
        for row in data.values():
            for k, v in row.items():
                if k in self.types:
                    self.types[k] = max(Dtype.var_to_type(v), self.types[k])
                else:
                    self.types[k] = Dtype.var_to_type(v)

        # reshape data to so it can be converted into a DataFrame
        data = [{'Replicate':k} | v for k, v in data.items()]
        self.df = self._metadata_data_to_df(data)

        return True


    def read(self, fname, metadata_format=None, exclude_null_from_skyline=True):
        '''
        Read sample metadata file and format dataframe to be added to sampleMetadata table.

        Parameters
        ----------
        fname: str
            The path to the metadata file.
        metadata_format: str
            One of ('tsv', 'csv', 'json')
        exclude_null_from_skyline: bool
            Exclude NULL annotations from skyline annotation csv?
            Default is True.

        Returns
        -------
        df: pd.DataFrame
            The sample metadata dataframe.
        types: dict
            A dictionary of metadata keys and Dtypes
        '''

        self.input_format = metadata_format
        self.input_file = fname

        if metadata_format:
            self.input_format = metadata_format
        else:
            self.input_format = splitext(fname)[1][1:]

        if self.input_format == 'tsv':
            with open(fname, 'r') as inF:
                data = list(DictReader(inF, delimiter='\t'))
        elif self.input_format == 'csv':
            with open(fname, 'r') as inF:
                rows = list(csv_reader(inF))

            # check if file is skyline annotations csv
            if rows[0][0] == 'ElementLocator':
                LOGGER.info('Found Skyline annotations csv.')
                self.input_format = 'skyline'
                return self._read_metadata_skyline_csv(rows, exclude_null=exclude_null_from_skyline)

            # Otherwise it is a normal csv
            if rows[0][0] == 'Replicate':
                headers = rows[0]
                rows = rows[1:]
                data = [dict(zip(headers, row)) for row in rows]
            else:
                LOGGER.error('Invalid metadata format!')
                return False

        elif self.input_format == 'json':
            return self._read_metadata_json(fname)
        else:
            LOGGER.error(f'Unknown metadata file format: {self.input_format}')
            return False

        self.df = self._metadata_data_to_df(data)
        self.types = self._infer_types(self.df)

        # set values of Dtype.NULL to NA
        self.df['annotationValue'] = self.df.apply(lambda x: pd.NA if self.types[x.annotationKey] is Dtype.NULL
                                                                   else x.annotationValue, axis=1)

        return True

