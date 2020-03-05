from unittest.mock import Mock

import pytest

from genomic_report.ipr import convert_statements_to_alterations


@pytest.fixture()
def graphkb_conn():
    def make_rid_list(*values):
        return [{'@rid': v} for v in values]

    return_values = [
        make_rid_list('disease'),  # disease call 1
        [],  # disease call 2
        make_rid_list('approved1', 'approved2'),
        make_rid_list('ther1'),
        make_rid_list('diag1'),
        make_rid_list('prog1'),
        make_rid_list('bio1', 'bio2'),
    ]
    query_mock = Mock()
    query_mock.side_effect = return_values
    conn = Mock(query=query_mock)

    return conn


def base_graphkb_statement(disease_id='disease'):
    statement = {
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
        'relevance': {'@rid': 'relevance_rid', 'displayName': 'relevance_display_name'},
        '@rid': 'statement_rid',
    }
    return statement


class TestConvertStatementsToAlterations:
    def test_disease_match(self, graphkb_conn):
        statement = base_graphkb_statement('disease')
        result = convert_statements_to_alterations(graphkb_conn, [statement], 'disease')

        assert len(result) == 1
        row = result[0]
        assert row['kb_event_key'] == 'variant_rid'
        assert row['kb_entry_key'] == 'statement_rid'
        assert row['matched_cancer']
        assert row['kbVariant'] == 'KRAS increased expression'
        assert row['association'] == 'relevance_display_name'

    def test_no_disease_match(self, graphkb_conn):
        statement = base_graphkb_statement('other')
        result = convert_statements_to_alterations(graphkb_conn, [statement], 'disease')

        assert len(result) == 1
        row = result[0]
        assert not row['matched_cancer']

    def test_multiple_disease_not_match(self, graphkb_conn):
        statement = base_graphkb_statement('disease')
        statement['conditions'].append(
            {'@class': 'Disease', '@rid': 'other', 'displayName': 'disease_display_name'}
        )
        result = convert_statements_to_alterations(graphkb_conn, [statement], 'disease')

        assert len(result) == 1
        row = result[0]
        assert not row['matched_cancer']

    def test_biological_statement(self, graphkb_conn):
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'bio1'

        result = convert_statements_to_alterations(graphkb_conn, [statement], 'disease')
        assert len(result) == 1
        row = result[0]
        assert row['kb_entry_type'] == 'biological'
        assert row['alterationType'] == 'biological'

    def test_prognostic_statement(self, graphkb_conn):
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'prog1'

        result = convert_statements_to_alterations(graphkb_conn, [statement], 'disease')
        assert len(result) == 1
        row = result[0]
        assert row['kb_entry_type'] == 'prognostic'
        assert row['alterationType'] == 'prognostic'

    def test_diagnostic_statement(self, graphkb_conn):
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'diag1'

        result = convert_statements_to_alterations(graphkb_conn, [statement], 'disease')
        assert len(result) == 1
        row = result[0]
        assert row['kb_entry_type'] == 'diagnostic'
        assert row['alterationType'] == 'diagnostic'

    def test_unapproved_therapeutic_statement(self, graphkb_conn):
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'ther1'
        statement['evidenceLevel'] = [{'@rid': 'other', 'displayName': 'level'}]

        result = convert_statements_to_alterations(graphkb_conn, [statement], 'disease')
        assert len(result) == 1
        row = result[0]
        assert row['kb_entry_type'] == 'therapeutic'
        assert row['alterationType'] == 'therapeutic'

    def test_approved_therapeutic_statement(self, graphkb_conn):
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'ther1'
        statement['evidenceLevel'] = [{'@rid': 'approved1', 'displayName': 'level'}]

        result = convert_statements_to_alterations(graphkb_conn, [statement], 'disease')
        assert len(result) == 1
        row = result[0]
        assert row['kb_entry_type'] == 'therapeutic'
        assert row['alterationType'] == 'therapeutic'
