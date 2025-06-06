
import ast, pathlib, re
from typing import Any, Dict, List
from types import SimpleNamespace

from groovy_parser.parser import parse_and_digest_groovy_content

_GROOVY_BOOL_NULL = {
    re.compile(r'\btrue\b'):  'True',
    re.compile(r'\bfalse\b'): 'False',
    re.compile(r'\bnull\b'):  'None',
}

_BLOCK_DICT = '__block_dict__'

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
    ''' Walk the AST subtree rooted at *node* and rebuild the Groovy source text. '''
    if not isinstance(node, dict):
        return ''

    # leaf token
    if 'leaf' in node:
        leaf_type = node['leaf']
        value     = node['value']

        if leaf_type == 'STRING_LITERAL':
            if '\n' in value:
                return "'''{}'''".format(value)
            return "'" + value + "'"

        return value
    return ''.join(_src(child) for child in node.get('children', []))


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

        key = ast.literal_eval(key_part) if key_part and key_part[0] in '\'"' else key_part
        mapping[key] = _g_to_py(val_part)
    return mapping


def _g_to_py(src: str):
    """
    Convert a Groovy literal string to its Python equivalent.

    * Square-bracket text is a MAP iff it has a *top-level* colon,
      otherwise it is treated as a LIST.
    * Groovy booleans / null are mapped to Python.
    * Everything else falls back to ast.literal_eval or raw string.
    """
    src = src.strip()

    if src.startswith('[') and src.endswith(']'):
        body = src[1:-1]
        if _has_top_level_colon(body):
            return _parse_map(body)

        parts = [p.strip() for p in _split_top(body, ',') if p.strip()]
        return [_g_to_py(p) for p in parts]

    if src.startswith('{') and src.endswith('}') and _has_top_level_colon(src[1:-1]):
        return _parse_map(src[1:-1])

    for pat, repl in _GROOVY_BOOL_NULL.items():
        src = re.sub(pat, repl, src)

    try:
        return ast.literal_eval(src)
    except Exception:
        return src


def _insert_nested(target: dict, dotted: str, value):
    '''
    Insert *value* under *dotted* path.

    Any dict we create along that path (i.e. because that level did not
    exist yet) is tagged with _BLOCK_DICT so we know it represents a
    section, not a literal map.
    '''
    parts = dotted.split('.')
    curr = target
    for key in parts[:-1]:
        if key not in curr:
            curr[key] = {_BLOCK_DICT: True}
        curr = curr[key]
    curr[parts[-1]] = value


def _leftmost_leaf(node):
    '''Return the first leaf token name inside *node* (skips whitespace).'''
    if not isinstance(node, dict):
        return None
    if 'leaf' in node:
        t = node['leaf']
        return None if t in {'NLS', 'SEMI'} else t
    for child in node.get('children', []):
        t = _leftmost_leaf(child)
        if t is not None:
            return t
    return None


def _is_block_header(left, right) -> bool:
    '''
    True when *left* starts with an identifier and *right*'s first significant token is LBRACE
    '''
    return _first_identifier(left) is not None and _leftmost_leaf(right) == 'LBRACE'


def _collect(node, into, prefix: str = ''):
    ''' DFS through *node*. '''

    if not isinstance(node, dict):
        return

    children = node.get('children', [])
    i = 0
    while i < len(children):
        left = children[i]

        if i + 1 < len(children) and _is_block_header(left, children[i + 1]):
            ident = _first_identifier(left)
            block = children[i + 1]
            _collect(block, into, f'{prefix}{ident}.')
            i += 2
            continue

        if isinstance(left, dict) and left.get('leaf') == 'ASSIGN':
            lhs = _flatten_identifiers(children[i - 1])
            rhs = _g_to_py(_src(children[i + 1]))
            key = f'{prefix}{lhs}' if prefix else lhs
            _insert_nested(into, key, rhs)

        _collect(left, into, prefix)
        i += 1


def _node_has_ident(node, name: str) -> bool:
    '''True if *node* (or any descendant) is IDENTIFIER==name.'''
    if not isinstance(node, dict):
        return False
    if node.get('leaf') == 'IDENTIFIER' and node.get('value') == name:
        return True
    return any(_node_has_ident(c, name) for c in node.get('children', []))


def _first_identifier(node):
    '''Return the first IDENTIFIER leaf under *node*, or None.'''
    if not isinstance(node, dict):
        return None
    if node.get('leaf') == 'IDENTIFIER':
        return node.get('value')
    for child in node.get('children', []):
        ident = _first_identifier(child)
        if ident is not None:
            return ident
    return None


def _desc_with_lbrace(node):
    '''Return the first descendant node that owns an LBRACE token.'''
    if not isinstance(node, dict):
        return None
    if any(isinstance(c, dict) and c.get('leaf') == 'LBRACE'
           for c in node.get('children', [])):
        return node
    for child in node.get('children', []):
        found = _desc_with_lbrace(child)
        if found is not None:
            return found
    return None


def _top_level_statements(tree):
    '''
    Yield the direct children that represent file-level statements,
    regardless of whether the parser inserted a script_statements wrapper.
    '''
    if not isinstance(tree, dict):
        return
    rule = tree.get('rule', [])
    if rule and rule[0] in ('script_statements', 'compilation_unit'):
        # dive one level down
        for child in tree.get('children', []):
            yield from _top_level_statements(child)
        return
    # otherwise this node itself is a statement candidate
    yield tree


def _find_params_block(tree):
    """
    Return the subtree that encloses the FIRST top-level

        params { ... }

    block found in the file.
    """
    for stmt in _top_level_statements(tree):
        if _first_identifier(stmt) == 'params':
            block = _desc_with_lbrace(stmt)
            if block is not None:
                return block
    raise ValueError('No top-level params { ... } block found in the config')


def _dict_to_ns(obj):
    '''
    Recursively walk *obj*.

    If we meet a dict tagged with _BLOCK_DICT, convert it (minus the tag)
    into SimpleNamespace(**children).
    Untagged dicts (literal maps) are left unchanged.
    '''
    if isinstance(obj, dict):
        # recurse into children first
        converted = {k: _dict_to_ns(v) for k, v in obj.items()
                     if k != _BLOCK_DICT}
        if obj.get(_BLOCK_DICT):
            return SimpleNamespace(**converted)
        return converted
    return obj


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
    return SimpleNamespace(**_dict_to_ns(out))


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


def _quote_str(s: str) -> str:
    """Return *s* as a Groovy string literal (single or triple quotes)."""
    if '\n' in s:
        return f"'''{s}'''"
    return f"'{s}'"


def _render_value(val, indent: str, lvl: int) -> str:
    """
    Convert *val* (scalar, list, dict) to Groovy syntax.
    Lists & maps are rendered inline; nested namespaces are handled by caller.
    """
    if isinstance(val, SimpleNamespace):
        raise ValueError('Nested namespace must be rendered by caller')

    if isinstance(val, dict):                                 # literal map
        pairs = [f'{k}: {_render_value(v, indent, lvl)}' for k, v in val.items()]
        return '[ ' + ', '.join(pairs) + ' ]'

    if isinstance(val, list):                                 # list
        items = ', '.join(_render_value(x, indent, lvl) for x in val)
        return f'[ {items} ]'

    if isinstance(val, str):
        return _quote_str(val)

    if val is None:
        return 'null'

    if isinstance(val, bool):
        return 'true' if val else 'false'

    # numbers or anything else repr-able
    return str(val)


def _dump_ns(ns: SimpleNamespace, fh, indent: str, lvl: int):
    """
    Recursively write *ns* to *fh* as Groovy config with indentation.
    """
    pad = indent * lvl
    for key, val in vars(ns).items():
        if isinstance(val, SimpleNamespace):
            fh.write(f'{pad}{key} {{\n')
            _dump_ns(val, fh, indent, lvl + 1)
            fh.write(f'{pad}}}\n')
        else:
            groovy_val = _render_value(val, indent, lvl)
            fh.write(f'{pad}{key} = {groovy_val}\n')


def write_params_namespace(obj: SimpleNamespace, file: str, indent: str = '  '):
    '''
    '''

    _needs_close = isinstance(file, str)
    fh = open(file, 'w', encoding='utf-8') if _needs_close else file

    try:
        _dump_ns(obj, fh, indent, 0)
    finally:
        if _needs_close:
            fh.close()