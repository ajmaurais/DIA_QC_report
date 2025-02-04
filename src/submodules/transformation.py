
from abc import ABC
import re
from textwrap import fill as fill_text

from .logger import LOGGER
from .logger import quiet_log_info, quiet_log_warning, quiet_log_error

HELP_TAB = '    '


def cammel_case(prefix, suffix):
    return f'{prefix}{suffix[0].upper()}{suffix[1:]}'


class Option:
    '''
    Individual method options.

    Attributes
    ----------
    name: str
        The name or key of the option.
    dtype: type
        The option value type.
    default: (str, int, bool, float)
        The default value for the option.
    choices: tuple
        The valid choices for option value.
    min_value: (int, float)
        For numeric dtypes the minimum value.
    min_inclusive: bool
        Is min_value inclusive?
    max_value: (int, float)
        For numeric dtype, the maxiumn value.
    max_inclusive: bool
        Is max_value inclusive?
    help_str: str
        The help string describing the option.
    '''

    def __init__(self, name,
                 default=None, dtype=int,
                 choices=None,
                 min_value=None, max_value=None,
                 min_inclusive=False, max_inclusive=False,
                 help_str=None):
        '''
        Parameters
        ----------
        name: str
            The name or key of the option.
        default: (str, int, bool, float)
            Default is None.
        dtype: type
            The option value type. Default is str.
        choices: tuple
            If this option is set min_value and max_value can not be set.
        min_value: (float, int)
             Default is None.
        min_inclusive: bool
            Default False.
        max_value: (float, int)
             Default is None.
        max_inclusive: bool
            Default False.
        help_str: str
            Default is None.

        Raises
        ------
        ValueError:
            If min_value or max_value can not be validly set.
        '''
        self.name = name
        self.dtype = dtype
        self.set_choices(choices)
        self.help_str = help_str

        self.max_value = None
        self.min_value = None
        self.max_inclusive = False
        self.min_inclusive = False

        if max_value is not None or min_value is not None:
            if self.choices is not None or not (self.dtype is int or self.dtype is float):
                raise ValueError('Cannot set min/max value for option `dtype` or `choices`!')

            if max_value is not None:
                self.max_value = max_value
                self.max_inclusive = max_inclusive

            if min_value is not None:
                self.min_value = min_value
                self.min_inclusive = min_inclusive

        if default is None:
            self.default = None
            return

        if self.is_valid(default, quiet=True):
            self.default = default
        else:
            raise ValueError(f"Invalid default argument: '{default}'!")


    def get_help(self, max_width=120):
        ''' Return help message for option '''

        msg = [f'{self.name}:']
        if self.choices is not None:
            msg.append(str(self.choices))
        else:
            msg.append(self.dtype.__name__)

        if self.min_value is not None or self.max_value is not None:
            msg.append(self._format_range(err_prefix=False))

        ret = ' '.join(msg)

        if self.help_str is not None or self.default is not None:
            second_line = ''
            if self.help_str:
                second_line += f'{self.help_str}'
            if self.default:
                second_line += f"{' ' if self.help_str else ''}Default is: "
                if self.dtype is str:
                    second_line += f"'{self.default}'"
                else:
                    second_line += str(self.default)
            ret += '\n' + fill_text(second_line, width=max_width,
                                    initial_indent=HELP_TAB, subsequent_indent=HELP_TAB)

        return ret


    def set_choices(self, choices):
        self.choices = None

        if choices is None:
            return

        # check that choices is tuple
        if not isinstance(choices, tuple):
            raise ValueError('choices must be an instance of tuple!')

        # check that all choices are of type self.dtype
        for choice in choices:
            if not isinstance(choice, self.dtype):
                raise ValueError(f"Choice '{choice}' is not of type {self.dtype}")

        self.choices = choices
        self.min_value = None
        self.max_value = None


    def _convert_to_dtype(self, value):
        if self.dtype is bool:
            if isinstance(value, bool):
                return value
            if value.lower() in ('true', 't', '1', 'on'):
                return True
            if value.lower() in ('false', 'f', '0', 'off'):
                return False
            raise ValueError(f"Cannot convert '{value}' to boolean!")

        try:
            ret = self.dtype(value)
        except ValueError as e:
            raise ValueError(f"Cannot convert '{value}' to {self.dtype}!") from e

        return ret


    def valid_choice(self, value):
        '''
        Test whether value is one of accepted option choices.
        '''
        if self.choices is not None:
            if value not in self.choices:
                return False
        return True


    def _format_range(self, err_prefix=False):
        err_l = list()

        if self.min_value is not None:
            err_l.append(f">{'=' if self.min_inclusive else ''} {self.min_value}")
        if self.max_value is not None:
            err_l.append(f"<{'=' if self.max_inclusive else ''} {self.max_value}")

        if len(err_l) > 0:
            if err_prefix:
                return 'Value must be ' + ' and '.join(err_l) + '!'
            return ' and '.join(err_l)

        return ''


    def in_range(self, value):
        '''
        Test whether value is between min/max constraints.
        Value should be self.dtype
        '''
        if not (self.dtype is float or self.dtype is int):
            return True

        if self.min_value is not None:
            if self.min_inclusive:
                if not value >= self.min_value:
                    return False
            elif not value > self.min_value:
                return False

        if self.max_value is not None:
            if self.max_inclusive:
                if not value <= self.max_value:
                    return False
            elif not value < self.max_value:
                return False

        return True


    def can_convert(self, value):
        '''
        Test whether value can be converted to option dtype
        '''
        if self.dtype is bool:
            if value in (True, False):
                return True
            if isinstance(value, str) and value.lower() in ('true', 't', '1', 'on', 'false', 'f', '0', 'off'):
                return True
            if isinstance(value, int):
                return True
            return False

        try:
            self.dtype(value)
        except ValueError:
            return False

        return True


    def is_valid(self, value, quiet=False):
        '''
        Test whether value is valid for option given choice, or min/max constraints,
        and whether value can be converted to option dtype.
        '''
        if not self.can_convert(value):
            quiet_log_error(quiet, 'Cannot convert value <%s> to %s!', value, str(self.dtype))
            return False

        _value = self._convert_to_dtype(value)

        if not self.valid_choice(_value):
            quiet_log_error(quiet, "Value '%s' must be in %s!", _value, self.choices)
            return False

        if not self.in_range(_value):
            quiet_log_error(quiet, self._format_range(err_prefix=True))
            return False

        return True


    def parse(self, value, quiet=False):
        '''
        Convert option value to option dtype.

        Parameters
        ----------
        value: str
            The string to attempt to convert.

        Returns
        -------
        value:
            The value converted to dtype and in option constraints.
            None if value is not valid.
        '''
        if not self.is_valid(value, quiet):
            return None

        return self._convert_to_dtype(value)


class MethodOptions:
    '''
    A collection of keyword arguments passed to transformer classes.

    Attributes
    ----------
    options: dict
        A dictionary of Option(s) objects for the transformer class.
    values: dict
        A dictionary of parsed and validated option names and values.
    description: str
        The description for the transformer to show in the help message.
    OPTION_RE: re.Pattern
        Regex used to parse option string arguments.
    '''

    OPTION_RE = re.compile(r'''[ \t]*([a-zA-Z\-_]+)[ \t]*[= \t]{1}[ \t]*((['"])[a-zA-Z0-9\-_ \.]+\3|[a-zA-Z0-9\.\-_]+)''')

    def __init__(self, description=None):
        self.options = dict()
        self.values = dict()
        self.description = description


    def add_option(self, name, **kwargs):
        '''
        Add option for method.

        Parameters
        ----------
        name: str
            Name of option.
        kwargs: dict
            Additional arguments passed to Option constructor.
        '''
        self.options[name] = Option(name, **kwargs)


    def get_help(self, out, max_width=120):
        '''
        Print description and help message for options.

        Parameters
        ----------
        out: ostream
            Output stream to write to.
        '''
        if self.description:
            out.write(f'{self.description}\n\n')

        for i, opt in enumerate(self.options.values()):
            if i > 0:
                out.write('\n')
            out.write(opt.get_help(max_width=max_width))
        out.write('\n')


    def get_option_dict(self):
        ''' Get a dictionary of validated names and values with default arguments added. '''
        defaults = {n: o.default for n, o in self.options.items() if o.default is not None}

        # first add explicit options
        ret = self.values.copy()

        # add default values for arguments that have defaults but are not in ret
        for name, default in defaults.items():
            if name not in ret:
                ret[name] = default

        return ret


    @staticmethod
    def _get_re_key_value_str(match):
        if len(match.groups()) != 3:
            raise ValueError('Option match must have 3 groups!')
        value = match.group(2) if match.group(3) is None else match.group(2).strip(match.group(3))
        return match.group(1), value


    def parse_strings(self, args, quiet=False):
        '''
        Parse list of string arguments and populate values attribute
        with validated option names and values.

        args: list
            A list of string arguments where each element has the form '<name>=<value>'
        quiet: bool
            Write error message for in valid args to LOGGER?

        Returns
        -------
        success: bool
            True if all args could be parsed, False if not.
        '''

        all_good = True
        for arg in args:
            match = self.OPTION_RE.search(arg)
            if match is None:
                quiet_log_error(quiet, "Could not parse option argument: '%s'", arg)
                all_good = False
                continue

            key, value = self._get_re_key_value_str(match)
            if key not in self.options:
                quiet_log_error(quiet, "Unknown option: '%s'", key)
                all_good = False
                continue

            _value = self.options[key].parse(value)
            if _value is None:
                all_good = False
                continue

            self.values[key] = _value

        return all_good


class TransformationManagerBase(ABC):
    '''
    TransformationManager abastract base class.

    Attributes
    ----------
    conn: sqlite.Connection
        A connection to a precursor database
    precursors: pd.DataFrame
        Long formatted dataframe of precursor quantities.
    proteins: pd.DataFrame
        Long formatted dataframe of protein quantities.
    '''

    def __init__(self, conn=None):
        self.conn = conn
        self.precursors = None
        self.proteins = None


    def get_long_tables(self, use_db_ids=False):
        '''
        Get long formatted precursor and protein tables.

        Parameters
        ----------
        use_db_ids: bool
            Should dataframe have database protein and replicate IDs?
        '''

        if use_db_ids:
            return self.precursors, self.proteins

        proteins = self.proteins.copy()
        precursors = self.precursors.copy()

        # add replicate column
        cur = self.conn.cursor()
        cur.execute('SELECT id, replicate FROM replicates;')
        rep_ids = {int(x[0]): x[1] for x in cur.fetchall()}
        precursors['replicate'] = precursors['replicateId'].apply(lambda x: rep_ids[x])
        proteins['replicate'] = proteins['replicateId'].apply(lambda x: rep_ids[x])

        # add protein name column
        cur.execute('SELECT proteinId, name FROM proteins;')
        prot_ids = {int(x[0]): x[1] for x in cur.fetchall()}
        proteins['protein'] = proteins['proteinId'].apply(lambda x: prot_ids[x])

        return precursors, proteins


    # def get_wide_tables(self, normalized=True, use_db_ids=False):
    #     pass
