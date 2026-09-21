import os
import re

import pytest
from requests_ratelimiter import LimiterAdapter

from pori_python.graphkb import GraphKBConnection, util

EXCLUDE_BCGSC_TESTS = os.environ.get('EXCLUDE_BCGSC_TESTS') == '1'


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
    """Rate limiting must stay off unless GRAPHKB_RATE_LIMIT is explicitly truthy."""

    @pytest.fixture(autouse=True)
    def clear_env_var(self, monkeypatch):
        monkeypatch.delenv(util.RATE_LIMIT_ENV_VAR, raising=False)
        monkeypatch.delenv(util.RATE_LIMIT_PER_SECOND_ENV_VAR, raising=False)

    @pytest.mark.parametrize('value', ['1', 'true', 'TRUE', 'yes', 'on'])
    def test_env_var_truthy_values_enable(self, monkeypatch, value):
        monkeypatch.setenv(util.RATE_LIMIT_ENV_VAR, value)
        assert util.rate_limiting_enabled() is True

    @pytest.mark.parametrize('value', ['0', 'false', 'no', 'off', ''])
    def test_env_var_falsy_values_disable(self, monkeypatch, value):
        monkeypatch.setenv(util.RATE_LIMIT_ENV_VAR, value)
        assert util.rate_limiting_enabled() is False

    def test_env_var_unset_defaults_to_disabled(self):
        assert util.RATE_LIMIT_ENV_VAR not in os.environ
        assert util.rate_limiting_enabled() is False

    def test_connection_default_has_rate_limiting_disabled(self):
        conn = GraphKBConnection(url='http://localhost:8080')
        assert conn.rate_limiting_enabled is False

    def test_connection_env_var_enables_rate_limiting(self, monkeypatch):
        monkeypatch.setenv(util.RATE_LIMIT_ENV_VAR, '1')
        conn = GraphKBConnection(url='http://localhost:8080')
        assert conn.rate_limiting_enabled is True

    def test_explicit_none_overrides_env_var(self, monkeypatch):
        monkeypatch.setenv(util.RATE_LIMIT_ENV_VAR, '1')
        conn = GraphKBConnection(url='http://localhost:8080', limiter=None)
        assert conn.rate_limiting_enabled is False

    def test_explicit_limiter_overrides_env_var(self):
        custom_limiter = LimiterAdapter(per_second=1)
        conn = GraphKBConnection(url='http://localhost:8080', limiter=custom_limiter)
        assert conn.rate_limiting_enabled is True

    def test_use_global_cache_false_no_longer_raises_by_default(self):
        # limiter defaults to disabled now, so this combination should be valid
        conn = GraphKBConnection(url='http://localhost:8080', use_global_cache=False)
        assert conn.rate_limiting_enabled is False

    def test_rate_limit_per_second_defaults_when_unset(self):
        assert util.rate_limit_per_second() == util.DEFAULT_RATE_LIMIT_PER_SECOND

    def test_rate_limit_per_second_reads_override(self, monkeypatch):
        monkeypatch.setenv(util.RATE_LIMIT_PER_SECOND_ENV_VAR, '25')
        assert util.rate_limit_per_second() == 25

    @pytest.mark.parametrize('value', ['not-a-number', '0', '-5'])
    def test_rate_limit_per_second_rejects_invalid_values(self, monkeypatch, value):
        monkeypatch.setenv(util.RATE_LIMIT_PER_SECOND_ENV_VAR, value)
        with pytest.raises(ValueError):
            util.rate_limit_per_second()

    def test_connection_env_var_uses_overridden_rate(self, monkeypatch):
        monkeypatch.setenv(util.RATE_LIMIT_ENV_VAR, '1')
        monkeypatch.setenv(util.RATE_LIMIT_PER_SECOND_ENV_VAR, '25')
        conn = GraphKBConnection(url='http://localhost:8080')
        assert conn.rate_limiting_enabled is True
