import io
import json
import os
import re

import pytest
import requests
import urllib3
from requests.adapters import HTTPAdapter
from requests_ratelimiter import LimiterAdapter

from pori_python.graphkb import GraphKBConnection, util

EXCLUDE_BCGSC_TESTS = os.environ.get('EXCLUDE_BCGSC_TESTS') == '1'


class CountingAdapter(HTTPAdapter):
    """Test transport adapter that returns a canned JSON response without any
    real network I/O, and counts how many times it was actually invoked (i.e.
    how many requests were NOT served from cache).
    """

    def __init__(self, *args, **kwargs):
        self.calls = 0
        super().__init__(*args, **kwargs)

    def send(self, request, **kwargs):
        self.calls += 1
        body = json.dumps({'call': self.calls}).encode('utf-8')
        raw = urllib3.HTTPResponse(
            body=io.BytesIO(body),
            status=200,
            preload_content=False,
            headers={},
            request_url=request.url,
        )
        return self.build_response(request, raw)


class OntologyTerm:
    def __init__(self, name, sourceId, displayName):
        self.name = name
        self.sourceId = sourceId
        self.displayName = displayName


@pytest.fixture(scope='module')
def conn() -> GraphKBConnection:
    conn = GraphKBConnection(url=os.environ['GRAPHKB_URL'])
    conn.login(os.environ['GRAPHKB_USER'], os.environ['GRAPHKB_PASS'])
    return conn


class TestLooksLikeRid:
    @pytest.mark.parametrize('rid', ['#3:4', '#50:04', '#-3:4', '#-3:-4', '#3:-4'])
    def test_valid(self, rid):
        assert util.looks_like_rid(rid)

    @pytest.mark.parametrize('rid', ['-3:4', 'KRAS'])
    def test_invalid(self, rid):
        assert not util.looks_like_rid(rid)


@pytest.mark.parametrize(
    'input,result',
    [
        ['GP5:p.Leu113His', 'GP5:p.L113H'],
        ['GP5:p.Lys113His', 'GP5:p.K113H'],
        ['CDK11A:p.Arg536Gln', 'CDK11A:p.R536Q'],
        ['APC:p.Cys1405*', 'APC:p.C1405*'],
        ['ApcTer:p.Cys1405*', 'ApcTer:p.C1405*'],
        ['GP5:p.Leu113_His114insLys', 'GP5:p.L113_H114insK'],
        ['NP_003997.1:p.Lys23_Val25del', 'NP_003997.1:p.K23_V25del'],
        ['LRG_199p1:p.Val7del', 'LRG_199p1:p.V7del'],
    ],
)
def test_convert_aa_3to1(input, result):
    assert util.convert_aa_3to1(input) == result


class TestStripParentheses:
    @pytest.mark.parametrize(
        'breakRepr,StrippedBreakRepr',
        [
            ['p.(E2015_Q2114)', 'p.E2015_Q2114'],
            ['p.(?572_?630)', 'p.?572_?630'],
            ['g.178916854', 'g.178916854'],
            ['e.10', 'e.10'],
        ],
    )
    def test_stripParentheses(self, breakRepr, StrippedBreakRepr):
        assert util.stripParentheses(breakRepr) == StrippedBreakRepr


class TestStripRefSeq:
    @pytest.mark.parametrize(
        'breakRepr,StrippedBreakRepr',
        [
            ['p.L2209', 'p.2209'],
            ['p.?891', 'p.891'],
            # TODO: ['p.?572_?630', 'p.572_630'],
        ],
    )
    def test_stripRefSeq(self, breakRepr, StrippedBreakRepr):
        assert util.stripRefSeq(breakRepr) == StrippedBreakRepr


class TestStripDisplayName:
    @pytest.mark.parametrize(
        'opt,stripDisplayName',
        [
            [{'displayName': 'ABL1:p.T315I', 'withRef': True, 'withRefSeq': True}, 'ABL1:p.T315I'],
            [{'displayName': 'ABL1:p.T315I', 'withRef': False, 'withRefSeq': True}, 'p.T315I'],
            [{'displayName': 'ABL1:p.T315I', 'withRef': True, 'withRefSeq': False}, 'ABL1:p.315I'],
            [{'displayName': 'ABL1:p.T315I', 'withRef': False, 'withRefSeq': False}, 'p.315I'],
            [
                {'displayName': 'chr3:g.41266125C>T', 'withRef': False, 'withRefSeq': False},
                'g.41266125>T',
            ],
            [
                {
                    'displayName': 'chrX:g.99662504_99662505insG',
                    'withRef': False,
                    'withRefSeq': False,
                },
                'g.99662504_99662505insG',
            ],
            [
                {
                    'displayName': 'chrX:g.99662504_99662505dup',
                    'withRef': False,
                    'withRefSeq': False,
                },
                'g.99662504_99662505dup',
            ],
            # TODO: [{'displayName': 'VHL:c.330_331delCAinsTT', 'withRef': False, 'withRefSeq': False}, 'c.330_331delinsTT'],
            # TODO: [{'displayName': 'VHL:c.464-2G>A', 'withRef': False, 'withRefSeq': False}, 'c.464-2>A'],
        ],
    )
    def test_stripDisplayName(self, opt, stripDisplayName):
        assert util.stripDisplayName(**opt) == stripDisplayName


class TestStringifyVariant:
    @pytest.mark.parametrize(
        'hgvs_string,opt,stringifiedVariant',
        [
            ['VHL:c.345C>G', {'withRef': True, 'withRefSeq': True}, 'VHL:c.345C>G'],
            ['VHL:c.345C>G', {'withRef': False, 'withRefSeq': True}, 'c.345C>G'],
            ['VHL:c.345C>G', {'withRef': True, 'withRefSeq': False}, 'VHL:c.345>G'],
            ['VHL:c.345C>G', {'withRef': False, 'withRefSeq': False}, 'c.345>G'],
            [
                '(LMNA,NTRK1):fusion(e.10,e.12)',
                {'withRef': False, 'withRefSeq': False},
                'fusion(e.10,e.12)',
            ],
            ['ABCA12:p.N1671Ifs*4', {'withRef': False, 'withRefSeq': False}, 'p.1671Ifs*4'],
            ['x:y.p22.33copyloss', {'withRef': False, 'withRefSeq': False}, 'y.p22.33copyloss'],
            # TODO: ['MED12:p.(?34_?68)mut', {'withRef': False, 'withRefSeq': False}, 'p.(34_68)mut'],
            # TODO: ['FLT3:p.(?572_?630)_(?572_?630)ins', {'withRef': False, 'withRefSeq': False}, 'p.(572_630)_(572_630)ins'],
        ],
    )
    def test_stringifyVariant_parsed(self, conn, hgvs_string, opt, stringifiedVariant):
        opt['variant'] = conn.parse(hgvs_string)
        assert util.stringifyVariant(**opt) == stringifiedVariant

    # Based on the assumption that these variants are in the database.
    # createdAt date help avoiding errors if assumption tuns to be false
    @pytest.mark.parametrize(
        'rid,createdAt,stringifiedVariant',
        [
            ['#157:0', 1565627324397, 'p.315I'],
            ['#157:79', 1565627683602, 'p.776_777insVGC'],
            ['#158:35317', 1652734056311, 'c.1>G'],
        ],
    )
    @pytest.mark.skipif(EXCLUDE_BCGSC_TESTS, reason='db-dependent rids')
    def test_stringifyVariant_positional(self, conn, rid, createdAt, stringifiedVariant):
        opt = {'withRef': False, 'withRefSeq': False}
        variant = conn.get_record_by_id(rid)
        if variant and variant.get('createdAt', None) == createdAt:
            assert util.stringifyVariant(variant=variant, **opt) == stringifiedVariant


class TestGraphKBConnection:
    def test_version(self, conn):
        version = conn.version
        assert version['db'] in [
            'production',
            'production-sync-dev',
            'production-sync-staging',
        ]
        SEMANTIC_VERSIONING_REGEX = re.compile(r'^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)$')
        assert SEMANTIC_VERSIONING_REGEX.match(version['api'])
        assert SEMANTIC_VERSIONING_REGEX.match(version['parser'])
        assert SEMANTIC_VERSIONING_REGEX.match(version['schema'])

    def test_get_related_records(self, conn):
        base = util.convert_to_rid_list(
            conn.query({'target': 'Vocabulary', 'filters': {'name': 'missense'}})
        )
        records = conn.get_related_records(
            base=base,
            ontology='Vocabulary',
            subgraphType='similar',
            returnProperties=['displayName'],
        )
        assert 'missense mutation' in list(map(lambda x: x['displayName'], records.values()))

    def test_get_related_terms(self, conn):
        # with defaults
        vocab_terms = conn.get_related_terms(
            terms='missense',
        )
        assert 'missense mutation' in vocab_terms

        # overriding ontology & subgraphType defaults
        disease_terms = conn.get_related_terms(
            terms='all solid tumors',
            ontology='Disease',
            subgraphType='parents',
        )
        assert 'cancer' in disease_terms


class TestRateLimitingOptIn:
    """Rate limiting is opt-in per-connection via the `limiter_kwargs` argument.

    Passing a truthy `limiter_kwargs` dict (e.g. `{'per_second': 3}`) mounts a
    `LimiterAdapter`; leaving it unset (None/default) or empty leaves the plain
    `HTTPAdapter` in place, i.e. rate limiting disabled. Regardless of
    `limiter_kwargs`, rate limiting is always forced off while running under
    pytest (`PYTEST_CURRENT_TEST` set) or when `only_if_cached=True`.
    """

    def test_default_has_no_rate_limiting(self):
        conn = GraphKBConnection(url='http://localhost:8080')
        assert isinstance(conn.http.adapters['https://'], HTTPAdapter)
        assert not isinstance(conn.http.adapters['https://'], LimiterAdapter)

    def test_empty_limiter_kwargs_has_no_rate_limiting(self):
        conn = GraphKBConnection(url='http://localhost:8080', limiter_kwargs={})
        assert not isinstance(conn.http.adapters['https://'], LimiterAdapter)

    def test_limiter_kwargs_forced_off_while_under_pytest(self, caplog):
        # PYTEST_CURRENT_TEST is set while running under pytest, so rate
        # limiting should stay disabled even when limiter_kwargs is provided.
        assert 'PYTEST_CURRENT_TEST' in os.environ
        conn = GraphKBConnection(url='http://localhost:8080', limiter_kwargs={'per_second': 3})
        assert not isinstance(conn.http.adapters['https://'], LimiterAdapter)
        assert 'rate limiting is by default turned off' in caplog.text

    def test_limiter_kwargs_enables_rate_limiting_outside_pytest(self, monkeypatch):
        monkeypatch.delenv('PYTEST_CURRENT_TEST', raising=False)
        conn = GraphKBConnection(url='http://localhost:8080', limiter_kwargs={'per_second': 3})
        assert isinstance(conn.http.adapters['https://'], LimiterAdapter)

    def test_only_if_cached_forces_rate_limiting_off(self, monkeypatch, caplog):
        monkeypatch.delenv('PYTEST_CURRENT_TEST', raising=False)
        conn = GraphKBConnection(
            url='http://localhost:8080', only_if_cached=True, limiter_kwargs={'per_second': 3}
        )
        assert conn.only_if_cached is True
        assert 'rate limiting is by default turned off' in caplog.text

    def test_limiter_kwargs_with_custom_session_raises(self):
        with pytest.raises(NotImplementedError):
            GraphKBConnection(
                url='http://localhost:8080',
                session=requests.Session(),
                limiter_kwargs={'per_second': 3},
            )


class TestCacheFilter:
    """Unit tests for util.cache_filter, which decides which responses requests-cache
    is allowed to store: all GETs, but only POSTs to the /query endpoint (other POSTs
    create content and must never be cached).
    """

    class _FakeRequest:
        def __init__(self, method, url):
            self.method = method
            self.url = url

    class _FakeResponse:
        def __init__(self, method, url):
            self.request = TestCacheFilter._FakeRequest(method, url)

    def test_get_requests_are_always_cacheable(self):
        response = self._FakeResponse('GET', 'http://fake/api/statement')
        assert util.cache_filter(response) is True

    @pytest.mark.parametrize('url', ['http://fake/api/query', 'http://fake/api/query/'])
    def test_post_to_query_endpoint_is_cacheable(self, url):
        response = self._FakeResponse('POST', url)
        assert util.cache_filter(response) is True

    @pytest.mark.parametrize('url', ['http://fake/api/statement', 'http://fake/api/query/similar'])
    def test_post_to_other_endpoints_is_not_cacheable(self, url):
        response = self._FakeResponse('POST', url)
        assert util.cache_filter(response) is False


class TestCachingBehavior:
    """End-to-end tests that mount a counting fake adapter under GraphKBConnection's
    session so we can assert whether a real request was made (adapter.calls
    incremented) or the response was served from the requests-cache layer instead.
    """

    def test_repeated_get_is_served_from_cache(self):
        conn = GraphKBConnection(url='http://fake-graphkb', use_global_cache=True)
        adapter = CountingAdapter()
        conn.http.mount('http://', adapter)

        first = conn.http.get('http://fake-graphkb/api/version')
        second = conn.http.get('http://fake-graphkb/api/version')

        assert adapter.calls == 1
        assert first.json() == second.json()
        assert getattr(second, 'from_cache', False) is True

    def test_no_cache_header_bypasses_cache(self):
        conn = GraphKBConnection(url='http://fake-graphkb', use_global_cache=True)
        adapter = CountingAdapter()
        conn.http.mount('http://', adapter)

        conn.http.get('http://fake-graphkb/api/version')
        conn.http.get('http://fake-graphkb/api/version', headers={'Cache-Control': 'no-cache'})

        assert adapter.calls == 2

    def test_only_if_cached_returns_504_when_not_cached(self):
        conn = GraphKBConnection(url='http://fake-graphkb', use_global_cache=True)
        adapter = CountingAdapter()
        conn.http.mount('http://', adapter)

        resp = conn.http.get(
            'http://fake-graphkb/api/never-requested',
            headers={'Cache-Control': 'only-if-cached'},
        )

        assert resp.status_code == 504
        assert adapter.calls == 0

    def test_only_if_cached_returns_cached_value_without_new_request(self):
        conn = GraphKBConnection(url='http://fake-graphkb', use_global_cache=True)
        adapter = CountingAdapter()
        conn.http.mount('http://', adapter)

        conn.http.get('http://fake-graphkb/api/version')
        resp = conn.http.get(
            'http://fake-graphkb/api/version', headers={'Cache-Control': 'only-if-cached'}
        )

        assert resp.status_code == 200
        assert adapter.calls == 1

    def test_post_to_query_endpoint_is_cached(self):
        conn = GraphKBConnection(url='http://fake-graphkb', use_global_cache=True)
        adapter = CountingAdapter()
        conn.http.mount('http://', adapter)

        conn.http.post('http://fake-graphkb/api/query', data='{}')
        conn.http.post('http://fake-graphkb/api/query', data='{}')

        assert adapter.calls == 1

    def test_post_to_non_query_endpoint_is_never_cached(self):
        conn = GraphKBConnection(url='http://fake-graphkb', use_global_cache=True)
        adapter = CountingAdapter()
        conn.http.mount('http://', adapter)

        conn.http.post('http://fake-graphkb/api/statement', data='{}')
        conn.http.post('http://fake-graphkb/api/statement', data='{}')

        assert adapter.calls == 2

    def test_use_global_cache_false_hits_adapter_every_time(self):
        conn = GraphKBConnection(url='http://fake-graphkb', use_global_cache=False)
        adapter = CountingAdapter()
        conn.http.mount('http://', adapter)

        conn.http.get('http://fake-graphkb/api/version')
        conn.http.get('http://fake-graphkb/api/version')

        assert adapter.calls == 2

    def test_use_global_cache_false_raises_with_cache_name(self):
        with pytest.raises(NotImplementedError):
            GraphKBConnection(
                url='http://fake-graphkb', use_global_cache=False, cache_name='somewhere'
            )
