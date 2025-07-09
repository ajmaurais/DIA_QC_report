
import re
import json
from csv import reader as csv_reader, DictReader
from os.path import splitext
from difflib import get_close_matches

from jsonschema import validate, ValidationError
import pandas as pd

from .dtype import Dtype, NA_RE
from .logger import quiet_log_info, quiet_log_error, quiet_log_warning


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


def closest_match(var: str, candidates: list[str]) -> str:
    '''
    Return closets string match to *var* from *candidates*.

    Returns
    -------
    str or None
        Closest match if found, None otherwise.
    '''
    if var in candidates:
        return var

    # Try to find a close match (diff ratio >= 0.6)
    close = get_close_matches(var, candidates, n=1, cutoff=0.5)
    if close is None or len(close) == 0:
        return None
    return close[0]


class Metadata():
    '''
    A class to read and process replicate metadata.

    Attributes
    ----------
    df: pd.DataFrame
        The sample metadata dataframe.
    types: dict
        A dictionary of metadata keys and Dtypes
    input_format: str
        The input file extension.
    self.input_file: str
        The path of the input file
    '''

    def __init__(self, quiet=False):
        self.df = pd.DataFrame({ 'Replicate': [], 'annotationKey': [], 'annotationValue': [] })
        self.types = {}
        self.input_format = None
        self.input_file = None
        self.quiet = quiet


    def copy(self):
        ret = Metadata()
        ret.df = self.df.copy() if self.df is not None else None
        ret.types = self.types.copy()
        ret.input_format = self.input_format
        ret.input_file = self.input_file
        return ret


    @staticmethod
    def combine_types(lhs, rhs):
        combined = lhs.copy()
        for key, dtype in rhs.items():
            if key in combined:
                combined[key] = max(dtype, combined[key])
            else:
                combined[key] = dtype
        return combined


    def add(self, rhs):
        if not isinstance(rhs, Metadata):
            raise TypeError(f'Cannot combine {type(rhs)} with Metadata!')

        if rhs.df.empty:
            return self

        duplicate_reps = set(self.df['Replicate']) & set(rhs.df['Replicate'])
        if len(duplicate_reps) > 0:
            raise ValueError(f"Cannot combine Metadata with duplicate replicates: {', '.join(duplicate_reps)}")

        self.df = pd.concat([self.df, rhs.df], ignore_index=True)
        self.types = self.combine_types(self.types, rhs.types)

        self._complete_missing_keys()
        self._coerce_null_numerics()
        return self


    def _coerce_null_numerics(self):
        new_string_types = set()
        for var in self.types:
            # Coerce INT, FLOAT, or BOOL Dtype to string if there are any missing values
            if self.types[var] in (Dtype.INT, Dtype.FLOAT, Dtype.BOOL):
                for row in self.df.itertuples():
                    if row.annotationKey == var:
                        if pd.isna(row.annotationValue) or row.annotationValue is None:
                            self.df.at[row.Index, 'annotationValue'] = ''
                            new_string_types.add(var)


    def _complete_missing_keys(self):
        '''
        Ensure that all replicates have all keys in self.types.
        If a key is missing, it will be added with a value of None.
        '''
        for rep in self.df['Replicate'].drop_duplicates():
            for key in self.types:
                if not ((self.df['Replicate'] == rep) & (self.df['annotationKey'] == key)).any():
                    self.df = pd.concat([self.df, pd.DataFrame({
                        'Replicate': [rep],
                        'annotationKey': [key],
                        'annotationValue': [None]
                    })], ignore_index=True)


    def log_info(self, message, *args):
        quiet_log_info(self.quiet, message, *args)


    def log_error(self, message, *args):
        quiet_log_error(self.quiet, message, *args)


    def log_warning(self, message, *args):
        quiet_log_warning(self.quiet, message, *args)


    def validate(self):
        data = dict()
        for row in self.df.itertuples():
            if row.Replicate not in data:
                data[row.Replicate] = dict()
            if row.annotationKey in data[row.Replicate]:
                self.log_error("Duplicate key: '%s' in '%s'!", row.annotationKey, row.Replicate)
                return False
            if row.annotationKey not in self.types:
                self.log_error("Missing key: '%s' in Metadata.types!", row.annotationKey)
                return False
            data[row.Replicate][row.annotationKey] = self.types[row.annotationKey]

        return True


    def _normalize_nas(self):
        '''
        Convert values in self.data which match dtype.NA_RE to pd.NA
        '''

        def normalize_na_fxn(value):
            if value is None or NA_RE.search(value):
                return ''
            return value

        self.df['annotationValue'] = self.df['annotationValue'].apply(normalize_na_fxn)


    def df_to_dict(self):
        if not self.validate():
            return None

        ret = {rep: {} for rep in self.df['Replicate'].drop_duplicates().to_list()}
        for row in self.df.itertuples():
            ret[row.Replicate][row.annotationKey] = self.types[row.annotationKey].convert(row.annotationValue)

        # ensure all replicates have all keys
        for rep in ret:
            for key in self.types:
                if key not in ret[rep]:
                    ret[rep][key] = None

        return ret


    def __eq__(self, rhs):
        if not isinstance(rhs, Metadata):
            return False

        if len(self.types) != len(rhs.types):
            return False

        for l_key, l_type in self.types.items():
            if l_key not in rhs.types:
                return False
            if l_type is not rhs.types[l_key]:
                return False

        l_data = self.df_to_dict()
        r_data = rhs.df_to_dict()

        if len(l_data) != len(r_data):
            return False

        for l_rep in l_data:
            if l_rep not in r_data:
                return False
            for l_key, l_value in l_data[l_rep].items():
                if l_key not in r_data[l_rep]:
                    return False
                if (n_na := sum(pd.isna(v) for v in (l_value, r_data[l_rep][l_key]))) > 0:
                    if n_na != 2:
                        return False
                    continue

                if l_value != r_data[l_rep][l_key]:
                    return False

        return True


    def __ne__(self, rhs):
        return not self.__eq__(rhs)


    def __repr__(self):
        return f'Metadata(input_format={self.input_format})'


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
                types[row.annotationKey] = max(
                    types[row.annotationKey], Dtype.infer_type(row.annotationValue)
                )
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
                self.log_warning("Found unknown annotation object: %s", row[0])
                continue
            annotation_object = annotation_match.group(1)
            annotation_objects.add(annotation_object)

            if annotation_object == 'Replicate':
                row[0] = re.sub('^Replicate:/', '', row[0])
                new_rows.append(row)

        if 'Replicate' not in annotation_objects:
            self.log_warning('No Replicate annotations in csv!')
            return False

        for obj in annotation_objects:
            if obj != 'Replicate':
                self.log_warning("Found annotations for: %s in csv!", obj)

        data = [dict(zip(headers, row)) for row in new_rows]
        self.df = self._metadata_data_to_df(data)
        self.types = self._infer_types(self.df)

        # set values of Dtype.NULL to NA
        self._normalize_nas()
        self.df['annotationValue'] = self.df.apply(lambda x: pd.NA if self.types[x.annotationKey] is Dtype.NULL
                                                                   else x.annotationValue, axis=1)

        if exclude_null:
            self.df = self.df[~pd.isna(self.df['annotationValue'])].reset_index(drop=True)
            self.types = {k: v for k, v in self.types.items() if v is not Dtype.NULL}

        # This is an annoying edge case with metadata annotations exported from Skyline.
        # Boolean annotations are either True or blank instead of False.
        # Therefore, we will set all blank boolean annotations to False.
        for key, var_type in self.types.items():
            if var_type is Dtype.BOOL:
                sele = self.df['annotationKey'] == key
                self.df.loc[sele, 'annotationValue'] = self.df[sele]['annotationValue'].apply(lambda x: 'False' if x == '' else x)

        return True


    def _read_metadata_json(self, fp):
        data = json.load(fp)
        try:
            validate(data, JSON_SCHEMA)
        except ValidationError as e:
            self.log_error("Invalid metadata format:\n%s", e.message)
            return False

        # determine metadata types
        self.types = dict()
        for row in data.values():
            for k, v in row.items():
                if k in self.types:
                    self.types[k] = max(Dtype.var_to_type(v), self.types[k])
                else:
                    self.types[k] = Dtype.var_to_type(v)

        for var in self.types:
            # Coerce INT Dtype to string if there are any missing values
            if self.types[var] is Dtype.INT or self.types[var] is Dtype.FLOAT:
                for rep in data:
                    if var not in data[rep]:
                        self.log_warning("annotationKey '%s' is missing for replicate '%s'!", var, rep)
                        data[rep][var] = ''

                    if data[rep][var] is None:
                        data[rep][var] = ''
                    else:
                        data[rep][var] = str(data[rep][var])

        # reshape data to so it can be converted into a DataFrame
        data = [{'Replicate':k} | v for k, v in data.items()]
        self.df = self._metadata_data_to_df(data)

        self.df['annotationValue'] = self.df.apply(lambda x: pd.NA if self.types[x.annotationKey] is Dtype.NULL
                                                                   else x.annotationValue, axis=1)

        return True


    def read(self, input_file, metadata_format=None, exclude_null_from_skyline=True):
        '''
        Read sample metadata file and format dataframe to be added to sampleMetadata table.

        Parameters
        ----------
        input_file: str or file
            Input file path or file pointer.
        metadata_format: str
            One of ('tsv', 'csv', 'json')
        exclude_null_from_skyline: bool
            Exclude NULL annotations from skyline annotation csv?
            Default is True.

        Returns
        -------
        success: bool
            True if successful, False if not
        '''

        read_from_path = False
        try:
            if isinstance(input_file, str):
                inF = open(input_file, 'r')

                read_from_path = True
                self.input_file = input_file

                if metadata_format:
                    self.input_format = metadata_format
                else:
                    self.input_format = splitext(input_file)[1][1:]
            else:
                if not hasattr(input_file, 'read'):
                    raise TypeError(f'input_file must be a string or file-like object, got {type(input_file)}!')
                inF = input_file

                # check input_format
                if metadata_format is None:
                    raise RuntimeError('Must specify metadata_format when input_file is not a string!')
                self.input_format = metadata_format

            if self.input_format == 'tsv':
                data = list(DictReader(inF, delimiter='\t'))
            elif self.input_format == 'csv':
                rows = list(csv_reader(inF))

                # check if file is skyline annotations csv
                if rows[0][0] == 'ElementLocator':
                    quiet_log_info(self.quiet, 'Found Skyline annotations csv.')
                    self.input_format = 'skyline'
                    return self._read_metadata_skyline_csv(rows, exclude_null=exclude_null_from_skyline)

                # Otherwise it is a normal csv
                if rows[0][0] == 'Replicate':
                    headers = rows[0]
                    rows = rows[1:]
                    data = [dict(zip(headers, row)) for row in rows]
                else:
                    self.log_error('Invalid metadata format!')
                    return False

            elif self.input_format == 'json':
                return self._read_metadata_json(inF)
            else:
                self.log_error("Unknown metadata file format: %s", self.input_format)
                return False

            self.df = self._metadata_data_to_df(data)
            self.types = self._infer_types(self.df)

            # set values of Dtype.NULL to NA
            self._normalize_nas()
            self.df['annotationValue'] = self.df.apply(
                lambda x: pd.NA if self.types[x.annotationKey] is Dtype.NULL
                                else x.annotationValue,
                axis=1
            )

        finally:
            if read_from_path:
                inF.close()

        return True


    def get_wide_data(self):
        df_wide =  self.df.pivot(
            index='Replicate', columns='annotationKey', values='annotationValue'
        ).reset_index().rename_axis(None, axis=1)

        for var, dtype in self.types.items():
            if dtype is Dtype.NULL:
                continue

            col_type = dtype.to_pd_type()

            if dtype is Dtype.BOOL:
                if any(df_wide[var] == ''):
                    col_type = str
                else:
                    df_wide[var] = df_wide[var].map({'True': True, 'False': False})
                    continue

            if dtype is Dtype.INT:
                if any(df_wide[var] == ''):
                    col_type = str
            elif dtype is Dtype.FLOAT:
                if any(df_wide[var] == ''):
                    df_wide.loc[df_wide[var] == '', var] = None

            try:
                df_wide = df_wide.astype({var: col_type})
            except ValueError as e:
                raise ValueError(f'Error converting column {var} to type {col_type}: {e}')

        return df_wide


    def to_tsv(self, out):
        data = self.get_wide_data()
        data.to_csv(out, sep='\t', index=False)


    def to_csv(self, out):
        data = self.get_wide_data()
        data.to_csv(out, index=False)


    def to_json(self, out):
        data = self.df_to_dict()
        json.dump(data, out, indent='\t')


    def to_skyline_annotations(self, out):
        data = self.df.pivot(
            index='Replicate', columns='annotationKey', values='annotationValue'
        ).reset_index()

        def write_csv_row(elems):
            if any(not isinstance(x, str) for x in elems):
                print(f'Invalid elements in row: {elems}')
            out.write('"{}"\n'.format('","'.join(elems)))

        # convert NULL columns to empty string
        for col, col_type in self.types.items():
            if col_type is Dtype.NULL:
                data[col] = ''

        annotation_headers = [x for x in data.columns if x != 'Replicate']
        # write header
        write_csv_row(['ElementLocator'] + [f'annotation_{x}' for x in annotation_headers])

        for _, row in data.iterrows():
            line = [f'Replicate:/{row.Replicate}']
            for header in annotation_headers:
                line.append(row[header])
            write_csv_row(line)


    def to_skyline_definitions(self, out):
        for name, dtype in self.types.items():
            out.write(f'--annotation-name="{name}" --annotation-targets=replicate')
            out.write(f' --annotation-type={dtype.to_sky_type()}\n')