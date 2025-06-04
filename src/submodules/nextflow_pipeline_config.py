
from groovy_parser.parser import parse_and_digest_groovy_content
import ast, pathlib, re
from typing import Any, Dict, List

QUANT_SPECTRA_SCHEMA = {
    "oneOf": [
        {
            "type": "string",
            'minProperties': 1
        },
        {
            "type": "array",
            "items": {"type": "string"},
            'minProperties': 1
        },
        {
            "type": "object",
            "additionalProperties": {
                "oneOf": [
                    {
                        "type": "string",
                        'minProperties': 1
                    },
                    {
                        "type": "array",
                        "items": {"type": "string"},
                        'minProperties': 1
                    },
                ]
            },
            'minProperties': 1
        },
    ],
}

CHROM_SPECTRA_SCHEMA = {
    "oneOf": [
        {
            "type": "string",
            'minProperties': 1
        },
        {
            "type": "array",
            "items": {"type": "string"},
            'minProperties': 1
        }
    ]
}

_GROOVY_BOOL_NULL = {
    re.compile(r'\btrue\b'):  'True',
    re.compile(r'\bfalse\b'): 'False',
    re.compile(r'\bnull\b'):  'None',
}

#  Low-level tree helpers
def _flatten_identifiers(node) -> str:
    '''Join IDENTIFIER leaves into a dotted name:  qc_report.color_vars'''
    parts: List[str] = []
    stack = [node]
    while stack:
        n = stack.pop()
        if isinstance(n, dict):
            if n.get('leaf') == 'IDENTIFIER':
                parts.append(n['value'])
            stack.extend(reversed(n.get('children', [])))
    return '.'.join(parts)


def _src(node) -> str:
    '''Re-assemble source code for any subtree.'''
    if isinstance(node, dict):
        if 'leaf' in node:
            return node['value']
        return ''.join(_src(c) for c in node.get('children', []))
    return ''


def _split_top(expr: str, delimiter: str):
    '''
    Split *expr* on *delimiter* only when the delimiter is at the top level
    (outside any quotes or nested brackets).
    '''
    parts, buf = [], []
    depth = 0
    in_single = in_double = False
    for ch in expr:
        if ch == "'" and not in_double:
            in_single = not in_single
        elif ch == '"' and not in_single:
            in_double = not in_double
        elif not in_single and not in_double:
            if ch in '([{':
                depth += 1
            elif ch in ')]}':
                depth -= 1
            elif ch == delimiter and depth == 0:
                parts.append(''.join(buf))
                buf = []
                continue
        buf.append(ch)
    parts.append(''.join(buf))
    return parts


def _has_top_level_colon(text: str) -> bool:
    '''True if *text* contains a colon that is outside quotes/brackets.'''
    depth = 0
    in_single = in_double = False
    for ch in text:
        if ch == "'" and not in_double:
            in_single = not in_single
        elif ch == '"' and not in_single:
            in_double = not in_double
        elif not in_single and not in_double:
            if ch in '([{':
                depth += 1
            elif ch in ')]}':
                depth -= 1
            elif ch == ':' and depth == 0:
                return True
    return False


def _split_key_value(item: str):
    '''
    Split a *single* map entry at the first top-level colon
    and return (key_part, value_part).
    '''
    depth = 0
    in_single = in_double = False
    buf_key, buf_val = [], []
    target = buf_key
    for ch in item:
        if ch == "'" and not in_double:
            in_single = not in_single
        elif ch == '"' and not in_single:
            in_double = not in_double
        elif not in_single and not in_double:
            if ch in '([{':
                depth += 1
            elif ch in ')]}':
                depth -= 1
            elif ch == ':' and depth == 0:
                target = buf_val
                continue
        target.append(ch)
    return ''.join(buf_key).strip(), ''.join(buf_val).strip()


def _parse_map(body: str):
    '''Parse a Groovy map body into a Python dict.'''
    mapping = {}
    for item in _split_top(body, ','):
        item = item.strip()
        if not item:
            continue
        key_part, val_part = _split_key_value(item)

        # key: quoted -> eval, else verbatim
        key = ast.literal_eval(key_part) if key_part and key_part[0] in '\'"' else key_part
        mapping[key] = _g_to_py(val_part)   # recurse for value
    return mapping


def _g_to_py(src: str):
    '''
    Convert a Groovy literal to its Python equivalent.
    Falls back to the raw text when not valid Python after conversion.
    '''
    src = src.strip()

    # 0. Braced literal: decide between map and brace-list
    if (src.startswith('{') and src.endswith('}')) or \
       (src.startswith('[') and src.endswith(']')):

        body = src[1:-1]

        if _has_top_level_colon(body):          # ← new robust test
            return _parse_map(body)             # treat as MAP

        # otherwise treat as brace-list {a,b,c}
        items = [s.strip() for s in _split_top(body, ',') if s.strip()]
        return [_g_to_py(item) for item in items]

    # 1. Replace Groovy booleans / null
    for pat, repl in _GROOVY_BOOL_NULL.items():
        src = re.sub(pat, repl, src)

    # 2. Evaluate if now valid Python, else return raw string
    try:
        return ast.literal_eval(src)
    except Exception:
        return src


def _insert_nested(target: dict, dotted: str, value):
    '''
    Insert *value* into *target* following a dotted path,
    e.g. dotted='panorama.upload'  ->  target['panorama']['upload'] = value
    '''
    parts = dotted.split('.')
    curr = target
    for key in parts[:-1]:
        curr = curr.setdefault(key, {})
    curr[parts[-1]] = value


def _collect(node: dict, into: Dict[str, Any]):
    if not isinstance(node, dict):
        return
    children = node.get('children', [])
    for i, tok in enumerate(children):
        if isinstance(tok, dict) and tok.get('leaf') == 'ASSIGN':
            lhs = _flatten_identifiers(children[i - 1])
            rhs = _g_to_py(_src(children[i + 1]))
            _insert_nested(into, lhs, rhs)
        _collect(tok, into)


def _node_has_ident(node, name: str) -> bool:
    '''True if *node* (or any descendant) is IDENTIFIER==name.'''
    if not isinstance(node, dict):
        return False
    if node.get('leaf') == 'IDENTIFIER' and node.get('value') == name:
        return True
    return any(_node_has_ident(c, name) for c in node.get('children', []))


def _contains_leaf(node, leaf_value: str) -> bool:
    '''True if a descendant token has .leaf == *leaf_value*.'''
    if not isinstance(node, dict):
        return False
    if node.get('leaf') == leaf_value:
        return True
    return any(_contains_leaf(c, leaf_value) for c in node.get('children', []))


def _find_params_block(tree) -> dict:
    '''
    Return the subtree that holds the braces of the top-level

        params { ... }

    call. We look for any node whose *first* child contains the identifier
    'params' and whose *second* child contains a '{' token (the closure).
    '''
    stack = [tree]
    while stack:
        node = stack.pop()
        children = node.get('children', [])
        if len(children) >= 2 and _node_has_ident(children[0], 'params') \
                          and _contains_leaf(children[1], 'LBRACE'):
            # children[1] is the whole argument_list / closure node;
            # that's good enough for the walker that gathers assignments.
            return children[1]
        stack.extend(reversed(children))
    raise ValueError('No `params { ... }` block found in the config.')


def parse_params(file=None, text=None) -> Dict[str, Any]:
    '''
    Parses a Nextflow pipeline configuration file and extracts parameters from the 'params' block.

    Parameters
    ----------
    file : str, optional
        The path to the Nextflow configuration file. If not provided, *text* must be given.
    text : str, optional
        The content of the Nextflow configuration file as a string. If not provided, *file* must be given.

    Returns:
        Dict[str, Any]: A dictionary containing parameter names and their corresponding values extracted from the 'params' block.

    Raises:
        FileNotFoundError: If the specified configuration file does not exist.
        ValueError: If the 'params' block cannot be found or parsed in the configuration file.
    '''
    if file is not None:
        text = pathlib.Path(file).read_text()
    if text is None:
        raise ValueError('Either `file` or `text` must be provided.')

    tree = parse_and_digest_groovy_content(text)
    params_block = _find_params_block(tree)
    out: Dict[str, Any] = {}
    _collect(params_block, out)
    return out


def param_to_list(param_variable):
    """
    Convert *param_variable* to a list.

    If it is already a list, return it unchanged.
    If it is a string, split on new-lines, strip whitespace, and drop blank lines.
    Otherwise, wrap it in a one-element list.

    This is the python translation of the function the pipeline uses to process list parameters.

    Parameters
    ----------
    param_variable : Any
        A single value, a multi-line string, or a list.

    Returns
    -------
    list
        The parameter represented as a list.
    """
    if isinstance(param_variable, list):
        return param_variable

    if isinstance(param_variable, str):
        return [line.strip() for line in param_variable.split('\n') if line.strip()]

    return [param_variable]