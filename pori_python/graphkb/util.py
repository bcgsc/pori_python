import requests
from requests.adapters import HTTPAdapter

import hashlib
import json
import logging
import os
import re
import time
from datetime import datetime
from requests_cache import CacheMixin
from typing import Any, Dict, Iterable, List, Optional, Union, cast
from urllib3.util.retry import Retry
from urllib.parse import urlsplit

from pori_python.types import ParsedVariant, PositionalVariant, Record

from .constants import DEFAULT_LIMIT, TYPES_TO_NOTATION, AA_3to1_MAPPING

# name the logger after the package to make it simple to disable for packages using this one as a dependency
# https://stackoverflow.com/questions/11029717/how-do-i-disable-log-messages-from-the-requests-library

logger = logging.getLogger('graphkb')


def convert_to_rid_list(records: Iterable[Record]) -> List[str]:
    """Given a list of records or record id strings, return their record IDs."""
    result = []
    for record in records:
        if isinstance(record, str):
            result.append(record)  # assume an @rid string
        else:
            result.append(record['@rid'])
    return result


class FeatureNotFoundError(Exception):
    pass


def looks_like_rid(rid: str) -> bool:
    """Check if an input string looks like a GraphKB ID."""
    if re.match(r'^#-?\d+:-?\d+$', rid):
        return True
    return False


def convert_aa_3to1(three_letter_notation: str) -> str:
    """Convert an Input string from 3 letter AA notation to 1 letter AA notation."""
    result = []

    if ':' in three_letter_notation:
        # do not include the feature/gene in replacements
        pos = three_letter_notation.index(':')
        result.append(three_letter_notation[: pos + 1])
        three_letter_notation = three_letter_notation[pos + 1 :]

    last_match_end = 0  # exclusive interval [ )

    for match in re.finditer(r'[A-Z][a-z][a-z]', three_letter_notation):
        # add the in-between string
        result.append(three_letter_notation[last_match_end : match.start()])
        text = three_letter_notation[match.start() : match.end()]
        result.append(AA_3to1_MAPPING.get(text, text))
        last_match_end = match.end()

    result.append(three_letter_notation[last_match_end:])
    return ''.join(result)


def join_url(base_url: str, *parts) -> str:
    """Join parts of a URL into a full URL."""
    if not parts:
        return base_url

    url = [base_url.rstrip('/')] + [part.strip('/') for part in parts]

    return '/'.join(url)


def millis_interval(start: datetime, end: datetime) -> int:
    """Millisecond time from start and end datetime instances."""
    diff = end - start
    millis = diff.days * 24 * 60 * 60 * 1000
    millis += diff.seconds * 1000
    millis += diff.microseconds // 1000
    return millis


class CustomSession(CacheMixin, requests.Session):
    pass


class GraphKBConnection:
    def __init__(
        self,
        url: str = os.environ.get('GRAPHKB_URL'),
        username: str = '',
        password: str = '',
        use_global_cache: bool = True,
        cache_name: str = '',
    ):
        """
        Docstring for __init__

        Args:
        - use_global_cache: cache requests across all requests to GKB
        - cache_db: Path or connection URL to the database which stors the requests cache. see https://requests-cache.readthedocs.io/en/v0.6.4/user_guide.html#cache-name
        """
        if use_global_cache:
            if not cache_name:
                self.http = CustomSession(
                    backend='memory',
                    cache_control=True,
                    allowable_methods=['GET', 'POST'],
                    ignored_parameters=['Authorization'],
                )
            else:
                self.http = CustomSession(
                    cache_name,
                    cache_control=True,
                    allowable_methods=['GET', 'POST'],
                    ignored_parameters=['Authorization'],
                )
        else:
            self.http = requests.Session()
        retries = Retry(
            total=100,
            connect=5,
            status=5,
            backoff_factor=5,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        self.http.mount('https://', HTTPAdapter(max_retries=retries))
        self.token = ''
        self.token_kc = ''
        self.url = url
        self.username = username
        self.password = password
        self.headers = {
            'Accept': 'application/json',
            'Content-Type': 'application/json',
        }
        self.request_count = 0
        self.first_request: Optional[datetime] = None
        self.last_request: Optional[datetime] = None
        if username and password:
            self.login(username=username, password=password)

        # URL check
        if not self.url:
            raise ValueError('URL to a GraphKB API instance is required')

    @property
    def load(self) -> Optional[float]:
        if self.first_request and self.last_request:
            msec = millis_interval(self.first_request, self.last_request)
            if msec:
                return self.request_count * 1000 / msec
        return None

    def request(
        self,
        endpoint: str,
        method: str = 'GET',
        headers: Optional[dict[str, str]] = None,
        **kwargs,
    ) -> Dict:
        """Request wrapper to handle adding common headers and logging.

        Args:
            endpoint (string): api endpoint, excluding the base uri
            method (str, optional): the http method. Defaults to 'GET'.

        Returns:
            dict: the json response as a python dict
        """
        url = join_url(self.url, endpoint)
        self.request_count += 1
        connect_timeout = 7
        read_timeout = 61

        # don't want to use a read timeout if the request is not idempotent
        # otherwise you may wind up making unintended changes
        timeout = None
        if endpoint in ['query', 'parse']:
            timeout = (connect_timeout, read_timeout)

        request_headers = {}
        request_headers.update(self.headers)
        if headers is not None:
            request_headers.update(headers)

        start_time = datetime.now()

        if not self.first_request:
            self.first_request = start_time
        self.last_request = start_time

        # using a manual retry as well as using the requests Retry() object because
        # a ConnectionError or OSError might be thrown and we still want to retry in those cases.
        # about catching OSError as well as ConnectionError:
        # https://stackoverflow.com/questions/74253820
        attempts = range(15)
        need_refresh_login = False
        for attempt in attempts:
            if attempt > 0:
                time.sleep(2)  # wait between retries
            try:
                if need_refresh_login:
                    self.refresh_login()
                    need_refresh_login = False

                self.request_count += 1
                resp = self.http.request(
                    method, url, headers=request_headers, timeout=timeout, **kwargs
                )
                if resp.status_code == 401 or resp.status_code == 403:
                    logger.debug(f'/{endpoint} - {resp.status_code} - retrying')
                    # try to re-login if the token expired
                    need_refresh_login = True
                    continue
                else:
                    break
            except (requests.exceptions.ConnectionError, OSError) as err:
                if attempt < len(attempts) - 1:
                    logger.debug(f'/{endpoint} - {str(err)} - retrying')
                    continue
                raise err
            except Exception as err2:
                raise err2

        timing = millis_interval(start_time, datetime.now())
        logger.debug(f'/{endpoint} - {resp.status_code} - {timing} ms')  # type: ignore

        try:
            resp.raise_for_status()
        except requests.exceptions.HTTPError as err:
            # try to get more error details
            message = str(err)
            try:
                message += ' ' + resp.json()['message']
            except Exception:
                pass

            raise requests.exceptions.HTTPError(message)

        return resp.json()

    def post(self, uri: str, data: Dict = {}, **kwargs) -> Dict:
        """Convenience method for making post requests."""
        return self.request(uri, method='POST', data=json.dumps(data), **kwargs)

    def login_demo(self) -> None:
        """
        KBDEV-1328
        GraphKB login in GSC's PORI online demo is done in 2 steps:
        1. get a first token from KeyCloak using username and password; self.login_demo()
        2. get a second token from the GraphKB API using keyCloakToken; self.login()
        """
        url_parts = urlsplit(self.url)
        base_url = f'{url_parts.scheme}://{url_parts.netloc}'

        try:
            resp = requests.request(
                url=f'{base_url}/auth/realms/PORI/protocol/openid-connect/token',
                method='POST',
                data={
                    'client_id': 'GraphKB',
                    'grant_type': 'password',
                    'password': self.password,
                    'username': self.username,
                },
            )
        except Exception as err:
            logger.debug(f'unable to fetch a token from KeyCloak: {err}')
            raise err

        resp.raise_for_status()
        content = resp.json()
        self.token_kc = content['access_token']

    def login(self, username: str, password: str, pori_demo: bool = False) -> None:
        self.username = username
        self.password = password
        connect_timeout = 7
        read_timeout = 61

        # KBDEV-1328. Alt. GraphKB login for GSC's PORI online demo
        if pori_demo or 'pori-demo' in self.url:
            logger.warning('login demo')
            self.login_demo()

        # use requests package directly to avoid recursion loop on login failure
        attempts = range(10)
        for attempt in attempts:
            if attempt > 0:
                time.sleep(2)  # wait between retries
            try:
                self.request_count += 1
                resp = requests.request(
                    url=f'{self.url}/token',
                    method='POST',
                    headers=self.headers,
                    timeout=(connect_timeout, read_timeout),
                    data=json.dumps(
                        # KBDEV-1328. Alt. GraphKB login for GSC's PORI online demo
                        {'keyCloakToken': self.token_kc}
                        if self.token_kc
                        else {'username': username, 'password': password}
                    ),
                )
                break
            except (requests.exceptions.ConnectionError, OSError) as err:
                if attempt < len(attempts) - 1:
                    logger.debug(f'/login - {str(err)} - retrying')
                    continue
                raise err
            except Exception as err2:
                raise err2
        resp.raise_for_status()
        content = resp.json()
        self.token = content['kbToken']
        self.headers['Authorization'] = self.token

    def refresh_login(self) -> None:
        self.login(self.username, self.password)

    def query(
        self,
        request_body: Dict = {},
        paginate: bool = True,
        ignore_cache: bool = False,
        force_refresh: bool = False,
        limit: int = DEFAULT_LIMIT,
    ) -> List[Record]:
        """
        Query GraphKB
        """
        headers = {}
        if ignore_cache or force_refresh:
            headers = {'Cache-Control': 'no-cache'}

        result: List[Record] = []
        while True:
            content = self.post(
                'query',
                data={**request_body, 'limit': limit, 'skip': len(result)},
                headers=headers,
            )
            records = content['result']
            result.extend(records)
            if len(records) < limit or not paginate:
                break

        return result

    def parse(self, hgvs_string: str, requireFeatures: bool = False) -> ParsedVariant:
        content = self.post(
            'parse', data={'content': hgvs_string, 'requireFeatures': requireFeatures}
        )
        return cast(ParsedVariant, content['result'])

    def get_records_by_id(self, record_ids: List[str]) -> List[Record]:
        if not record_ids:
            return []
        result = self.query({'target': record_ids})
        if len(record_ids) != len(result):
            raise AssertionError(
                f'The number of Ids given ({len(record_ids)}) does not match the number of records fetched ({len(result)})'
            )
        return result

    def get_record_by_id(self, record_id: str) -> Record:
        result = self.get_records_by_id([record_id])
        return result[0]

    def get_source(self, name: str) -> Record:
        source = self.query({'target': 'Source', 'filters': {'name': name}})
        if len(source) != 1:
            raise AssertionError(f'Unable to unqiuely identify source with name {name}')
        return source[0]


def get_rid(conn: GraphKBConnection, target: str, name: str) -> str:
    """
    Retrieve a record by name and target

    Args:
        conn: GraphKBConnection
        target: record type to query
        name: the name of the record to retrieve

    Returns:
        str: @rid of the record

    Raises:
        AssertionError: if the term was not found or more than 1 match was found (expected to be unique)
    """
    result = conn.query(
        {'target': target, 'filters': {'name': name}, 'returnProperties': ['@rid']},
        ignore_cache=False,
    )
    assert len(result) == 1, f"unable to find unique '{target}' ID for '{name}'"

    return result[0]['@rid']


def stripParentheses(breakRepr: str) -> str:
    match = re.search(r'^([a-z])\.\((.+)\)$', breakRepr)

    if match:
        return f'{match.group(1)}.{match.group(2)}'
    return breakRepr


def stripRefSeq(breakRepr: str) -> str:
    # 1 leading RefSeq
    match = re.search(r'^([a-z])\.([A-Z]*|\?)([0-9]*[A-Z]*)$', breakRepr)
    if match:
        return f'{match.group(1)}.{match.group(3)}'

    # TODO: Deal with cases like "p.?889_?890", "chr4:g.55593604_55593605delGGinsTT", ...

    return breakRepr


def stripDisplayName(displayName: str, withRef: bool = True, withRefSeq: bool = True) -> str:
    match = re.search(r'^(.*)(\:)(.*)$', displayName)
    if match and not withRef:
        if withRefSeq:
            return match.group(3)
        displayName = match.group(2) + match.group(3)

    match = re.search(r'^(.*\:)([a-z]\.)(.*)$', displayName)
    if match and not withRefSeq:
        ref: str = match.group(1) if match.group(1) != ':' else ''
        prefix: str = match.group(2)
        rest: str = match.group(3)
        new_matches: Union[bool, object] = True

        # refSeq before position
        while new_matches:
            new_matches = re.search(r'(.*)([A-Z]|\?)([0-9]+)(.*)', rest)
            if new_matches:
                rest = new_matches.group(1) + new_matches.group(3) + new_matches.group(4)

        # refSeq before '>'
        new_matches = re.search(r'^([0-9]*)([A-Z]*|\?)(\>)(.*)$', rest)
        if new_matches:
            rest = new_matches.group(1) + new_matches.group(3) + new_matches.group(4)

        displayName = ref + prefix + rest

    # TODO: Deal with more complex cases like fusion, cds with offset (ex. 'VHL:c.464-2G>A')
    # and other complex cases (ex. 'VHL:c.330_331delCAinsTT')

    return displayName


def stringifyVariant(
    variant: Union[PositionalVariant, ParsedVariant], withRef: bool = True, withRefSeq: bool = True
) -> str:
    """
    Convert variant record to a string representation (displayName/hgvs)

    Args:
        variant: the input variant
        withRef (bool, optional): include the reference part
        withRefSeq (bool, optional): include the reference sequence in the variant part

    Returns:
        str: The string representation
    """

    displayName: str = variant.get('displayName') or ''  # type: ignore

    # If variant is a PositionalVariant (i.e. variant with a displayName) and
    # we already have the appropriate string representation,
    # then return it right away
    if displayName != '' and (withRef and withRefSeq):
        return displayName

    # If variant is a PositionalVariant (i.e. variant with a displayName) and
    # we DO NOT have the appropriate string representation,
    # then strip unwanted features, then return it right away
    if displayName != '':
        return stripDisplayName(displayName, withRef, withRefSeq)

    # If variant is a ParsedVariant (i.e. variant without a displayName yet),
    # the following will return a stringify representation (displayName/hgvs) of that variant
    # based on: https://github.com/bcgsc/pori_graphkb_parser/blob/ae3738842a4c208ab30f58c08ae987594d632504/src/variant.ts#L206-L292

    parsed = variant
    result: List[str] = []

    # Extracting parsed values into individual variables
    break1Repr: str = str(parsed.get('break1Repr', ''))
    break2Repr: str = str(parsed.get('break2Repr', ''))
    multiFeature: bool = bool(parsed.get('multiFeature'))
    noFeatures: bool = bool(parsed.get('noFeatures'))
    notationType: str = str(parsed.get('notationType', ''))
    reference1: str = ''
    if ref1 := parsed.get('reference1'):
        if isinstance(ref1, str):
            reference1 = ref1
        else:
            reference1 = ref1.get('displayName', str(ref1))
    reference2: str = ''
    if ref2 := parsed.get('reference2'):
        if isinstance(ref2, str):
            reference2 = ref2
        else:
            reference2 = ref2.get('displayName', str(ref2))
    refSeq: str = parsed.get('refSeq') or ''
    truncation: int = parsed.get('truncation') or 0  # type: ignore
    variantType: str = parsed.get('type', '')
    untemplatedSeq: str = parsed.get('untemplatedSeq') or ''
    untemplatedSeqSize: int = parsed.get('untemplatedSeqSize') or 0

    # formating notationType
    if not notationType:
        notationType = TYPES_TO_NOTATION.get(variantType, re.sub(r'\s', '-', variantType))

    # If multiFeature
    if multiFeature or (reference2 != '' and reference1 != reference2):
        if withRef and not noFeatures:
            result.append(f'({reference1}:{reference2})')
        result.append(notationType)
        if withRefSeq:
            break1Repr_noParentheses = stripParentheses(break1Repr)
            break2Repr_noParentheses = stripParentheses(break2Repr)
            result.append(f'({break1Repr_noParentheses},{break2Repr_noParentheses})')
        else:
            break1Repr_noParentheses_noRefSeq = stripRefSeq(stripParentheses(break1Repr))
            break2Repr_noParentheses_noRefSeq = stripRefSeq(stripParentheses(break2Repr))
            result.append(
                f'({break1Repr_noParentheses_noRefSeq},{break2Repr_noParentheses_noRefSeq})'
            )
        if untemplatedSeq != '':
            result.append(untemplatedSeq)
        elif untemplatedSeqSize:
            result.append(str(untemplatedSeqSize))
        return ''.join(result)

    # Continuous notation...

    # Reference
    if withRef and not noFeatures:
        result.append(f'{reference1}:')

    # BreakRep
    if withRefSeq:
        result.append(break1Repr)
        if break2Repr != '':
            result.append(f'_{break2Repr[2:]}')
    else:
        result.append(stripRefSeq(break1Repr))
        if break2Repr != '':
            result.append(f'_{stripRefSeq(break2Repr)[2:]}')

    # refSeq, truncation, notationType, untemplatedSeq, untemplatedSeqSize
    if any(i in notationType for i in ['ext', 'fs']) or (
        notationType == '>' and break1Repr.startswith('p.')
    ):
        result.append(untemplatedSeq)
    if notationType == 'mis' and break1Repr.startswith('p.'):
        result.append(untemplatedSeq)
    elif notationType != '>':
        if notationType == 'delins':
            if withRefSeq:
                result.append(f'del{refSeq}ins')
            else:
                result.append('delins')
        else:
            result.append(notationType)
        if truncation and truncation != 1:
            if truncation < 0:
                result.append(str(truncation))
            else:
                result.append(f'*{truncation}')
        if any(i in notationType for i in ['dup', 'del', 'inv']):
            if withRefSeq:
                result.append(refSeq)
        if any(i in notationType for i in ['ins', 'delins']):
            if untemplatedSeq != '':
                result.append(untemplatedSeq)
            elif untemplatedSeqSize:
                result.append(str(untemplatedSeqSize))
    elif not break1Repr.startswith('p.'):
        if withRefSeq:
            refSeq = refSeq if refSeq != '' else '?'
        else:
            refSeq = ''
        untemplatedSeq = untemplatedSeq if untemplatedSeq != '' else '?'
        result.append(f'{refSeq}{notationType}{untemplatedSeq}')

    # TODO: Deal with more complexes cases like 'MED12:p.(?34_?68)mut'

    return ''.join(result)
