
import requests
import xml.etree.ElementTree as ET
from urllib.parse import urlparse, unquote
from base64 import b64encode


PANORA_PUBLIC_KEY = '7d503a4147133c448c6eaf83bc9b8bc22ace4b7f6d36ca61c9d1ca836c510d10'


def list_panorama_files(
    url: str,
    *,
    api_key: str | None = None,
    verify_ssl: bool = True,
    depth: int = 1
) -> list[str]:
    r'''
    List file names in a Panorama/LabKey WebDAV directory.

    The helper automatically disables ``requests`` environment helpers
    (``~/.netrc``, proxy variables, custom CA bundles) whenever an
    ``api_key`` is supplied, ensuring that the explicit credentials you
    pass are the only ones used for authentication.

    Parameters
    ----------
    url : str
        WebDAV directory URL. A trailing ``'/'`` is appended if missing.
    api_key : str, optional
        Personal or public Panorama API-key.  If *None*, the function
        falls back to whatever credentials ``requests`` finds in
        ``~/.netrc`` or proceeds anonymously.
    verify_ssl : bool, default ``True``
        Set to ``False`` to ignore TLS-certificate validation errors
        (e.g., with self-signed certs on dev servers).
    depth : int, default ``1``
        WebDAV *Depth* header.  ``0`` = directory itself, ``1`` =
        immediate children, etc.

    Returns
    -------
    list[str]
        File names contained *directly* in the specified directory
        (sub-directories are excluded).

    Raises
    ------
    requests.HTTPError
        Propagated after printing the server's response body if the HTTP
        status is 400 or higher.
    '''
    # Build request
    if not url.endswith('/'):
        url += '/'

    headers: dict[str, str] = {
        'Depth': str(depth),
        'Content-Type': 'text/xml; charset=utf-8',
    }

    if api_key:
        token = b64encode(f'apikey:{api_key}'.encode()).decode()
        headers['Authorization'] = f'Basic {token}'
    else:
        token = b64encode(f'apikey:{PANORA_PUBLIC_KEY}'.encode()).decode()
        headers['Authorization'] = f'Basic {token}'

    body = (
        '<?xml version=\'1.0\' encoding=\'utf-8\'?>'
        '<d:propfind xmlns:d=\"DAV:\">'
        '<d:prop><d:displayname/><d:resourcetype/></d:prop>'
        '</d:propfind>'
    )

    # Send request
    if api_key:
        session = requests.Session()
        session.trust_env = False          # ignore ~/.netrc, proxies, etc.
        request_fn = session.request
        dummy_auth = ()                    # prevents fallback auth
    else:
        request_fn = requests.request
        dummy_auth = None

    resp = request_fn(
        'PROPFIND',
        url,
        headers=headers,
        data=body,
        verify=verify_ssl,
        auth=dummy_auth,
    )

    if not resp.ok:
        print(f'HTTP {resp.status_code} error')
        print('----- response body -----')
        print(resp.text)
        print('-------------------------')
        resp.raise_for_status()

    # Parse XML response
    ns = {'d': 'DAV:'}
    root = ET.fromstring(resp.text)
    base_path = urlparse(url).path.rstrip('/') + '/'

    files: list[str] = []
    for node in root.findall('d:response', ns):
        href = node.find('d:href', ns).text
        full_path = unquote(urlparse(href).path)

        if full_path == base_path:
            continue   # skip the directory node itself

        is_dir = node.find('.//d:collection', ns) is not None
        if not is_dir:
            files.append(full_path.rsplit('/', 1)[-1])

    return files