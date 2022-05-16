import pytest
from graphkb import statement as gkb_statement
from graphkb import vocab as gkb_vocab
from graphkb.types import Statement
from unittest.mock import Mock

from ipr.ipr import convert_statements_to_alterations, filter_kb_matches, germline_kb_matches

DISEASE_RIDS = ['#138:12', '#138:13']
APPROVED_EVIDENCE_RIDS = ['approved1', 'approved2']
GERMLINE_VARIANTS = [
    {
        'altSeq': 'T',
        'chromosome': 'chr9',
        'endPosition': 84286011,
        'gene': 'SLC28A3',
        'germline': True,
        'hgvsCds': 'SLC28A3:c.1381C>T',
        'hgvsGenomic': 'chr9:g.84286011G>A',
        'hgvsProtein': 'SLC28A3:p.L461L',
        'key': '584ffc37bfc37efc41a53a221a93a1f3',
        'ncbiBuild': 'GRCh38',
        'normalAltCount': 37,
        'normalDepth': 37,
        'normalRefCount': 0,
        'proteinChange': 'p.L461L',
        'refSeq': 'C',
        'rnaAltCount': '',
        'rnaDepth': '',
        'rnaRefCount': '',
        'startPosition': 84286011,
        'transcript': 'ENST00000376238',
        'tumourAltCount': '',
        'tumourDepth': '',
        'tumourRefCount': '',
        'variant': 'SLC28A3:p.L461L',
        'variantType': 'mut',
        'zygosity': '',
    }
]

SOMATIC_VARIANTS = [
    {
        'altSeq': 'T',
        'chromosome': 'chr9',
        'endPosition': 84286011,
        'gene': 'SLC28A3',
        'germline': False,
        'hgvsCds': 'SLC28A3:c.1381C>T',
        'hgvsGenomic': 'chr9:g.84286011G>A',
        'hgvsProtein': 'SLC28A3:p.L461L',
        'key': '584ffc37bfc37efc41a53a221a93a1f3',
        'ncbiBuild': 'GRCh38',
        'normalAltCount': 0,
        'normalDepth': 37,
        'normalRefCount': 37,
        'proteinChange': 'p.L461L',
        'refSeq': 'C',
        'rnaAltCount': '',
        'rnaDepth': '',
        'rnaRefCount': '',
        'startPosition': 84286011,
        'transcript': 'ENST00000376238',
        'tumourAltCount': 37,
        'tumourDepth': 37,
        'tumourRefCount': 0,
        'variant': 'SLC28A3:p.L461L',
        'variantType': 'mut',
        'zygosity': '',
    }
]

PCP_KB_MATCHES = [
    {
        'approvedTherapy': False,
        'category': 'pharmacogenomic',
        'context': 'anthracyclines',
        'disease': '',
        'evidenceLevel': '',
        'externalSource': None,
        'externalStatementId': None,
        'kbContextId': '#122:20944',
        'kbRelevanceId': '#147:38',
        'kbStatementId': '#154:13387',
        'kbVariant': 'SLC28A3:c.1381C>T',
        'kbVariantId': '#159:5426',
        'matchedCancer': False,
        'reference': 'PMID: 27197003',
        'relevance': 'decreased toxicity',
        'reviewStatus': 'initial',
        'variant': '584ffc37bfc37efc41a53a221a93a1f3',
        'variantType': 'mut',
    },
    {
        'approvedTherapy': False,
        'category': 'cancer predisposition',
        'context': 'prostate adenocarcinoma [PRAD]',
        'disease': 'prostate adenocarcinoma [PRAD]',
        'evidenceLevel': 'MOAlmanac FDA-Approved',
        'externalSource': 'MOAlmanac',
        'externalStatementId': '621',
        'kbContextId': '#135:8764',
        'kbRelevanceId': '#147:32',
        'kbStatementId': '#155:13511',
        'kbVariant': 'BRCA1 mutation',
        'kbVariantId': '#161:938',
        'matchedCancer': False,
        'reference': 'MOAlmanac FDA-56',
        'relevance': 'pathogenic',
        'reviewStatus': None,
        'variant': '7158e8931eda66a7d0cf9e0313d82561',
        'variantType': 'mut',
    },
]


@pytest.fixture
def graphkb_conn():
    class QueryMock:
        return_values = [
            # get approved evidence levels
            [{'@rid': v} for v in APPROVED_EVIDENCE_RIDS],
        ]
        index = -1

        def __call__(self, *args, **kwargs):
            self.index += 1
            ret_val = self.return_values[self.index] if self.index < len(self.return_values) else []
            return ret_val

    conn = Mock(query=QueryMock(), cache={})

    return conn


def base_graphkb_statement(disease_id: str = 'disease', relevance_rid: str = 'other') -> Statement:
    statement = Statement(
        {
            'conditions': [
                {'@class': 'Disease', '@rid': disease_id, 'displayName': 'disease_display_name'},
                {
                    '@class': 'CategoryVariant',
                    '@rid': 'variant_rid',
                    'displayName': 'KRAS increased expression',
                },
            ],
            'evidence': [],
            'subject': None,
            'source': None,
            'sourceId': None,
            'relevance': {'@rid': relevance_rid, 'displayName': 'relevance_display_name'},
            '@rid': 'statement_rid',
        }
    )
    return statement


@pytest.fixture(autouse=True)
def mock_get_term_tree(monkeypatch):
    def mock_func(*pos, **kwargs):
        return [{'@rid': d} for d in DISEASE_RIDS]

    monkeypatch.setattr(gkb_vocab, 'get_term_tree', mock_func)


@pytest.fixture(autouse=True)
def mock_categorize_relevance(monkeypatch):
    def mock_func(_, relevance_id):
        return relevance_id

    monkeypatch.setattr(gkb_statement, 'categorize_relevance', mock_func)


class TestConvertStatementsToAlterations:
    def test_disease_match(self, graphkb_conn, mock_get_term_tree) -> None:
        statement = base_graphkb_statement(DISEASE_RIDS[0])
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )

        assert len(result) == 1
        row = result[0]
        assert row['kbVariantId'] == 'variant_rid'
        assert row['kbStatementId'] == 'statement_rid'
        assert row['matchedCancer']
        assert row['kbVariant'] == 'KRAS increased expression'
        assert row['relevance'] == 'relevance_display_name'

    def test_no_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement('other')
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )

        assert len(result) == 1
        row = result[0]
        assert not row['matchedCancer']

    def test_multiple_disease_not_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement('disease')
        statement['conditions'].append(
            {'@class': 'Disease', '@rid': 'other', 'displayName': 'disease_display_name'}
        )
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )

        assert len(result) == 1
        row = result[0]
        assert not row['matchedCancer']

    def test_biological(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'biological'

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 1
        row = result[0]
        assert row['category'] == 'biological'

    def test_prognostic_no_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'prognostic'

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 0

    def test_prognostic_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement(DISEASE_RIDS[0])
        statement['relevance']['@rid'] = 'prognostic'

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 1
        row = result[0]
        assert row['category'] == 'prognostic'

    def test_diagnostic(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'diagnostic'

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 1
        row = result[0]
        assert row['category'] == 'diagnostic'

    def test_unapproved_therapeutic(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'therapeutic'
        statement['evidenceLevel'] = [{'@rid': 'other', 'displayName': 'level'}]

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 1
        row = result[0]
        assert row['category'] == 'therapeutic'

    def test_approved_therapeutic(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'therapeutic'
        statement['evidenceLevel'] = [{'@rid': APPROVED_EVIDENCE_RIDS[0], 'displayName': 'level'}]

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 1
        row = result[0]
        assert row['category'] == 'therapeutic'


class TestKbmatchFilters:
    def test_germline_kb_matches(self):
        assert germline_kb_matches(
            PCP_KB_MATCHES, GERMLINE_VARIANTS
        ), "Germline variant improperly excluded by germline_kb_matches"
        assert not germline_kb_matches(
            PCP_KB_MATCHES, SOMATIC_VARIANTS
        ), "Somatic variant matched to KB pharmacogenomic by germline_kb_matches"

    def test_filter_kb_matches(self):
        assert filter_kb_matches(PCP_KB_MATCHES, []), "filter_kb_matches no filter returned nothing"
        # filter from GERO-238 - only report cancer predisposition matches from CGL source.
        kb_match_filters = [
            {'category': (True, ['pharmacogenomic'])},
            {'category': (True, ['cancer predisposition']), 'externalSource': (False, ['CGL'])},
        ]
        assert not filter_kb_matches(
            PCP_KB_MATCHES, kb_match_filters
        ), "filter_kb_matches returned non-CGL filter"
        kb_match_filters = [
            {'category': (True, ['pharmacogenomic'])},
            {
                'category': (True, ['cancer predisposition']),
                'externalSource': (False, ['MOAlmanac']),
            },
        ]
        assert filter_kb_matches(
            PCP_KB_MATCHES, kb_match_filters
        ), "filter_kb_matches failed externalSource filter"
