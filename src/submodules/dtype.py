
import re
from enum import Enum


NA_RE = re.compile(r'^(NA|NULL|#N/A|NaN)$')
BOOL_RE = re.compile(r'^(true|True|TRUE|false|FALSE|False)$')
INT_RE = re.compile(r'^[+\-]?\d+$')
FLOAT_RE = re.compile(r'^[+\-]?(\d+(\.\d*)?|\.\d+|\d+(\.\d*)?[eE][+\-]?\d+)$')


class Dtype(Enum):
    NULL=0
    BOOL=1
    INT=2
    FLOAT=3
    STRING=4

    def __str__(self):
        return self.name

    def __lt__(self, rhs):
        if isinstance(rhs, Dtype):
            return self.value < rhs.value

        raise ValueError(f'Cannot compare {type(self)} to {type(rhs)}!')


    def __ge__(self, rhs):
        if isinstance(rhs, Dtype):
            return self.value >= rhs.value

        raise ValueError(f'Cannot compare {type(self)} to {type(rhs)}!')

    def convert(self, val):
        ''' Convert val to the same type as self. '''
        if self is Dtype.NULL:
            return None
        if self is Dtype.BOOL:
            if val.lower() in ('true', '1', 't'):
                return True
            if val.lower() in ('false', '0', 'f'):
                return False
            raise ValueError(f'Invalid conversion to boolean: "{val}"')
        if self is Dtype.INT:
            return int(val)
        if self is Dtype.FLOAT:
            return float(val)
        return val


    @staticmethod
    def var_to_type(v):
        ''' Convert variable with compatible type to Dtype. '''
        if v is None:
            return Dtype.NULL
        if isinstance(v, bool):
            return Dtype.BOOL
        if isinstance(v, int):
            return Dtype.INT
        if isinstance(v, float):
            return Dtype.FLOAT
        if isinstance(v, str):
            return Dtype.STRING
        raise ValueError(f'Unknown type: {str(v)}')


    @staticmethod
    def infer_type(s):
        '''
        Infer datatype for string.

        Parameters
        ----------
        s: str
            The string to test.

        Returns
        -------
        Dtype object
        '''
        if s == '' or NA_RE.search(s) is not None:
            return Dtype.NULL
        if BOOL_RE.search(s):
            return Dtype.BOOL
        if INT_RE.search(s):
            return Dtype.INT
        if FLOAT_RE.search(s):
            return Dtype.FLOAT
        return Dtype.STRING
