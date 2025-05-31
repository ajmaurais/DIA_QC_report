
import os
import xml.etree.ElementTree as ET
from urllib.parse import urlparse, unquote
from base64 import b64encode

import requests


PANORA_PUBLIC_KEY = '7d503a4147133c448c6eaf83bc9b8bc22ace4b7f6d36ca61c9d1ca836c510d10'


def list_panorama_files(
    url: str,
    api_key: str | None = None,
    panorama_pubilc: bool = False,
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
    panorama_pubilc : bool, default False
        If True, use the Panorama Public API key.
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
    elif panorama_pubilc:
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


def download_webdav_file(
    url: str,
    dest_path: str | None = None,
    api_key: str | None = None,
    verify_ssl: bool = True,
    chunk_size: int = 8192,
) -> str:
    r'''
    Download a single file from a LabKey WebDAV URL.

    Parameters
    ----------
    url : str
        Direct WebDAV URL of the file to fetch.
    dest_path : str, optional
        Local path to write the file.  If None, the basename of url
        is used in the current working directory.
    api_key : str, optional
        LabKey API key.
    verify_ssl : bool, default True
        False disables TLS-certificate validation.
    chunk_size : int, default 8192
        Number of bytes written per iteration.

    Returns
    -------
    str
        Absolute path of the downloaded file.

    Raises
    ------
    requests.HTTPError
        Propagated after printing the server's response body if the HTTP
        status code indicates an error.
    '''
    if dest_path is None:
        dest_path = unquote(urlparse(url).path.rsplit('/', 1)[-1])
    dest_path = os.path.abspath(dest_path)

    headers: dict[str, str] = {}
    if api_key:
        token = b64encode(f'apikey:{api_key}'.encode()).decode()
        headers['Authorization'] = f'Basic {token}'

    if api_key:
        session = requests.Session()
        session.trust_env = False
        get_fn = session.get
    else:
        get_fn = requests.get

    resp = get_fn(url, headers=headers, stream=True, verify=verify_ssl)
    if not resp.ok:
        error_message = (
            f'HTTP {resp.status_code} error while downloading {url}\n'
            '----- response body -----\n'
            f'{resp.text}\n'
            '-------------------------'
        )
        raise requests.HTTPError(error_message, response=resp)

    os.makedirs(os.path.dirname(dest_path) or '.', exist_ok=True)
    with open(dest_path, 'wb') as fh:
        for chunk in resp.iter_content(chunk_size=chunk_size):
            if chunk:                      # skip keep-alive chunks
                fh.write(chunk)

    return dest_path