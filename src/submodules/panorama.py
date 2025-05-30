
import requests
import xml.etree.ElementTree as ET
from urllib.parse import urlparse, unquote
from base64 import b64encode


PANORA_PUBLIC_KEY = '7d503a4147133c448c6eaf83bc9b8bc22ace4b7f6d36ca61c9d1ca836c510d10'


def list_panorama_files(
    url: str,
    api_key: str | None = None,
    verify_ssl: bool = True,
    depth: int = 1,
    full_url: bool = False,
) -> list[str]:
    r'''
    List file names (or full URLs) located **directly** in a LabKey WebDAV
    directory.

    Parameters
    ----------
    url : str
        WebDAV directory URL. A trailing '/' is appended if missing.
    api_key : str, optional
        LabKey API key.
    verify_ssl : bool, default True
        Set to False to skip TLS-certificate validation.
    depth : int, default 1
        WebDAV *Depth* header.
    full_url : bool, default False
        If True, return the complete WebDAV URL for each file rather
        than just the basename.

    Returns
    -------
    list[str]
        Either basenames or full URLs of the files immediately contained
        in *url*.  Sub-directories are excluded.

    Raises
    ------
    requests.HTTPError
        Propagated after printing the server's response body when an
        HTTP error is encountered.
    '''
    # Prepare request
    if not url.endswith('/'):
        url += '/'

    parsed = urlparse(url)
    headers = {
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
        '<?xml version="1.0" encoding="utf-8"?>'
        '<d:propfind xmlns:d="DAV:">'
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
        error_message = (
            f'HTTP {resp.status_code} error\n'
            '----- response body -----\n'
            f'{resp.text}\n'
            '-------------------------'
        )
        raise requests.HTTPError(error_message, response=resp)

    # Parse XML response
    ns = {'d': 'DAV:'}
    root = ET.fromstring(resp.text)
    base_path = parsed.path.rstrip('/') + '/'

    results = list()
    for node in root.findall('d:response', ns):
        href = node.find('d:href', ns).text
        full_path = unquote(urlparse(href).path)

        if full_path == base_path:
            continue

        is_dir = node.find('.//d:collection', ns) is not None
        if is_dir:
            continue

        name = full_path.rsplit('/', 1)[-1]
        results.append(
            f'{parsed.scheme}://{parsed.netloc}{full_path}' if full_url else name
        )

    return results