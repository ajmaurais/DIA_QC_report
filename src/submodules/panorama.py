
import os
import xml.etree.ElementTree as ET
from urllib.parse import urlparse, unquote
from base64 import b64encode
from io import BytesIO
import socket
from contextlib import closing

import requests

PANORAMA_URL = 'https://panoramaweb.org'
PANORAMA_PUBLIC_KEY = '7d503a4147133c448c6eaf83bc9b8bc22ace4b7f6d36ca61c9d1ca836c510d10'

_LIST_BODY = (
    '<?xml version="1.0" encoding="utf-8"?>'
    '<d:propfind xmlns:d="DAV:">'
    '<d:prop><d:displayname/><d:resourcetype/></d:prop>'
    '</d:propfind>'
)

def _encode_api_key(api_key: str) -> str:
    '''
    Encode the LabKey API key for use in the HTTP Basic Authorization header.

    Parameters
    ----------
    api_key : str
        LabKey API key.

    Returns
    -------
    str
        Base64-encoded string for the HTTP Basic Authorization header.
    '''
    return b64encode(f'apikey:{api_key}'.encode()).decode()


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
        headers['Authorization'] = f'Basic {_encode_api_key(api_key)}'

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
        data=_LIST_BODY,
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
        full_path = urlparse(href).path

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


def get_webdav_file(
    url: str,
    dest_path: str | None = None,
    api_key: str | None = None,
    verify_ssl: bool = True,
    chunk_size: int = 8192,
    return_text: bool = False
) -> str:
    r'''
    Download a single file from a LabKey WebDAV URL.

    Parameters
    ----------
    url : str
        Direct WebDAV URL of the file to fetch.
    dest_path : str, optional
        Local path to write the file.
        If None, the file is written to the current working directory.
    api_key : str, optional
        LabKey API key.
    return_text : bool, default False
        When True, the function streams the payload into memory and
        returns it as a string.

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
    if return_text:
        dest_fh = BytesIO()
        close_after = False
    else:
        if dest_path is None:
            dest_path = unquote(urlparse(url).path.rsplit('/', 1)[-1])
        dest_path = os.path.abspath(dest_path)
        os.makedirs(os.path.dirname(dest_path) or '.', exist_ok=True)
        dest_fh = open(dest_path, 'wb')
        close_after = True

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

    for chunk in resp.iter_content(chunk_size=chunk_size):
        if chunk:                      # skip keep-alive chunks
            dest_fh.write(chunk)

    if return_text:
        dest_fh.seek(0)
        text = dest_fh.read().decode('utf-8')
        return text

    if close_after:
        dest_fh.close()

    return dest_path


def webdav_file_exists(
    url: str,
    *,
    api_key: str | None = None,
    verify_ssl: bool = True,
    timeout: int | float = 10,
) -> bool:
    r'''
    Check whether a single WebDAV resource exists on a Panorama/LabKey server.

    Parameters
    ----------
    url : str
        Full WebDAV URL pointing directly to the file (not a directory).
    api_key : str, optional
        LabKey API key.
    verify_ssl : bool, default True
        Set to False to ignore TLS certificate validation.
    timeout : int or float, default 10
        Seconds to wait before the request is aborted.

    Returns
    -------
    bool
        True if the server returns HTTP 200, 204, or 207 (success
        for WebDAV resources); False if it returns HTTP 404.

    Raises
    ------
    requests.HTTPError
        Re-raised for any status code that is not 200, 204, 207, or 404
        (e.g. 401 Unauthorized, 403 Forbidden, 500 Server Error)
    '''
    # Build headers
    headers = dict()
    if api_key:
        headers['Authorization'] = f'Basic {_encode_api_key(api_key)}'

    # Choose transport
    if api_key:
        session = requests.Session()
        session.trust_env = False        # ignore ~/.netrc, proxies, etc.
        head_fn = session.head
    else:
        head_fn = requests.head

    # Send HEAD request
    resp = head_fn(url, headers=headers, allow_redirects=True,
                   verify=verify_ssl, timeout=timeout)

    if resp.status_code in (200, 204, 207):           # success codes
        return True
    if resp.status_code == 404:                       # not found
        return False

    # For any other code, show the server's message then raise.
    error_message = (
        f'HTTP {resp.status_code} while checking {url}\n'
        '----- response body -----\n'
        f'{resp.text}\n'
        '-------------------------'
    )
    raise requests.HTTPError(error_message, response=resp)


def have_internet(target=("8.8.8.8", 53), timeout=3.0):
    ''' Return True if a TCP handshake to target succeeds within timeout.'''
    try:
        with closing(socket.create_connection(target, timeout=timeout)):
            return True
    except OSError:
        return False


def url_exists(url: str, timeout: float = 5.0) -> bool:
    '''
    Return True if an HTTP HEAD request to *url* succeeds (status-code 200).

    Works with GitHub raw-file URLs and any well-behaved HTTP server.
    Falls back to a zero-byte GET if HEAD is not allowed.
    '''
    try:
        r = requests.head(url, allow_redirects=True, timeout=timeout)
        if r.status_code == 405: # HEAD not supported
            r = requests.get(url, headers={"Range": "bytes=0-0"},
                             stream=True, timeout=timeout)
        return r.status_code == 200
    except requests.RequestException:
        return False


def get_http_file(url: str, dest_path: str | None = None, return_text: bool = False,
                  timeout: float = 10.0, chunk_size=8192, verify_ssl: bool = True) -> str:
    '''
    Download a file from the given URL and save it to a local path or return its
    contents as a string.

    Parameters
    ----------
    url : str
        The HTTP(S) URL pointing to the text file.
    timeout : float
        How many seconds to wait for the server to send data before giving up.
    dest_path : str, optional
        Local path to write the file.
        If None, the file is written to the current working directory.
    return_text : bool, default False
        If True, the function streams the payload into memory and
        returns it as a string instead of writing to a file.

    Returns
    -------
    str
        Absolute path of the downloaded file, or the text content if
        `return_text` is True.
    '''

    if return_text:
        dest_fh = BytesIO()
        close_after = False
    else:
        if dest_path is None:
            dest_path = unquote(urlparse(url).path.rsplit('/', 1)[-1])
        dest_path = os.path.abspath(dest_path)
        os.makedirs(os.path.dirname(dest_path) or '.', exist_ok=True)
        dest_fh = open(dest_path, 'wb')
        close_after = True

    response = requests.get(url, stream=True, timeout=timeout, verify=verify_ssl)
    response.raise_for_status()

    for chunk in response.iter_content(chunk_size=chunk_size):
        if chunk:                      # skip keep-alive chunks
            dest_fh.write(chunk)

    if return_text:
        dest_fh.seek(0)
        text = dest_fh.read().decode('utf-8')
        return text

    if close_after:
        dest_fh.close()

    return dest_path