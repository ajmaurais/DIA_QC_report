
# groovy_params.py ----------------------------------------------------------
from groovy_parser.parser import parse_and_digest_groovy_content
import ast, pathlib, re
from typing import Any, Dict, Union, List

# ---------------------------------------------------------------------------
#  Low-level tree helpers
# ---------------------------------------------------------------------------
def _flatten_identifiers(node) -> str:
    """Join IDENTIFIER leaves into a dotted name:  qc_report.color_vars"""
    parts: List[str] = []
    stack = [node]
    while stack:
        n = stack.pop()
        if isinstance(n, dict):
            if n.get("leaf") == "IDENTIFIER":
                parts.append(n["value"])
            stack.extend(reversed(n.get("children", [])))
    return ".".join(parts)

def _src(node) -> str:
    """Re-assemble source code for any subtree."""
    if isinstance(node, dict):
        if "leaf" in node:        # a token
            return node["value"]
        return "".join(_src(c) for c in node.get("children", []))
    return ""

# ---------------------------------------------------------------------------
#  Groovy-literal → Python-object
# ---------------------------------------------------------------------------

_GROOVY_BOOL_NULL = {
    re.compile(r"\btrue\b"):  "True",
    re.compile(r"\bfalse\b"): "False",
    re.compile(r"\bnull\b"):  "None",
}
_KEY_UNQUOTED = re.compile(r"([a-zA-Z_]\w*)\s*:(?!//)")

def _g_to_py(src: str):
    """
    Convert a Groovy literal to the matching Python object.
    Falls back to the raw string if it cannot be evaluated.
    """
    src = src.strip()

    # 0. Groovy 'set' / brace-list:  {a,b,c}  ->  ['a', 'b', 'c']
    if src.startswith("{") and src.endswith("}") and ":" not in src:
        inner = src[1:-1]
        items = [part.strip() for part in inner.split(",") if part.strip()]
        return items

    # 1. Groovy map  [foo:'bar'] -> {'foo': 'bar'}
    if src.startswith("[") and ":" in src:
        src = _KEY_UNQUOTED.sub(r"'\1':", src)
        src = "{" + src[1:-1] + "}"

    # 2. Booleans / null
    for pat, repl in _GROOVY_BOOL_NULL.items():
        src = re.sub(pat, repl, src)

    # 3. Evaluate if valid Python, else return as-is
    try:
        return ast.literal_eval(src)
    except (ValueError, SyntaxError):
        return src

# ---------------------------------------------------------------------------
#  Main walker
# ---------------------------------------------------------------------------

def _insert_nested(target: dict, dotted: str, value):
    """
    Insert *value* into *target* following a dotted path,
    e.g. dotted='panorama.upload'  ->  target['panorama']['upload'] = value
    """
    parts = dotted.split(".")
    curr = target
    for key in parts[:-1]:
        curr = curr.setdefault(key, {})
    curr[parts[-1]] = value


def _collect(node: dict, into: Dict[str, Any]):
    if not isinstance(node, dict):
        return
    kids = node.get("children", [])
    for i, tok in enumerate(kids):
        if isinstance(tok, dict) and tok.get("leaf") == "ASSIGN":
            lhs = _flatten_identifiers(kids[i - 1])
            rhs = _g_to_py(_src(kids[i + 1]))
            _insert_nested(into, lhs, rhs)            # ← nested insert
        _collect(tok, into)


def _node_has_ident(node, name: str) -> bool:
    """True if *node* (or any descendant) is IDENTIFIER==name."""
    if not isinstance(node, dict):
        return False
    if node.get("leaf") == "IDENTIFIER" and node.get("value") == name:
        return True
    return any(_node_has_ident(c, name) for c in node.get("children", []))


def _contains_leaf(node, leaf_value: str) -> bool:
    """True if a descendant token has .leaf == *leaf_value*."""
    if not isinstance(node, dict):
        return False
    if node.get("leaf") == leaf_value:
        return True
    return any(_contains_leaf(c, leaf_value) for c in node.get("children", []))


def _find_params_block(tree) -> dict:
    """
    Return the subtree that holds the braces of the top-level

        params { ... }

    call. We look for any node whose *first* child contains the identifier
    'params' and whose *second* child contains a '{' token (the closure).
    """
    stack = [tree]
    while stack:
        node = stack.pop()
        kids = node.get("children", [])
        if len(kids) >= 2 and _node_has_ident(kids[0], "params") \
                          and _contains_leaf(kids[1], "LBRACE"):
            # kids[1] is the whole argument_list / closure node;
            # that's good enough for the walker that gathers assignments.
            return kids[1]
        stack.extend(reversed(kids))
    raise ValueError("No `params { ... }` block found in the config.")

# ---------------------------------------------------------------------------
#  Public API
# ---------------------------------------------------------------------------
def parse_params(path="pipeline.config") -> Dict[str, Any]:
    text = pathlib.Path(path).read_text()
    tree = parse_and_digest_groovy_content(text)
    params_block = _find_params_block(tree)
    out: Dict[str, Any] = {}
    _collect(params_block, out)
    return out

def get_param(path: str, dotted: str, default: Any = None) -> Any:
    return parse_params(path).get(dotted, default)

# ---------------------------------------------------------------------------

cfg = parse_params("pipeline.config")
for k, v in cfg.items():
    print(f"{k:30s} → {v!r}")


