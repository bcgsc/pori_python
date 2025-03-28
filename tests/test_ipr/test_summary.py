from unittest.mock import MagicMock
from copy import copy

from pori_python.ipr.summary import (
    GRAPHKB_GUI,
    get_preferred_drug_representation,
    substitute_sentence_template,
    get_ipr_analyst_comments,
)


class TestGetPreferredDrugRepresentation:
    def test_prefers_non_alias(self):
        api = MagicMock(
            query=MagicMock(
                side_effect=[
                    [],
                    [
                        {"sourceId": "1", "alias": False, "source": "source", "name": "name"},
                        {"sourceId": "2", "alias": True, "source": "source", "name": "name"},
                    ],
                ]
            )
        )
        rec = get_preferred_drug_representation(api, "anything")
        assert rec["sourceId"] == "1"

    def test_prefers_non_deprecated(self):
        api = MagicMock(
            query=MagicMock(
                side_effect=[
                    [],
                    [
                        {"sourceId": "1", "deprecated": True, "source": "source", "name": "name"},
                        {"sourceId": "2", "deprecated": False, "source": "source", "name": "name"},
                    ],
                ]
            )
        )
        rec = get_preferred_drug_representation(api, "anything")
        assert rec["sourceId"] == "2"

    def test_prefers_lower_sort_source(self):
        api = MagicMock(
            query=MagicMock(
                side_effect=[
                    [{"@rid": "source2", "sort": 0}, {"@rid": "source1", "sort": 1}],
                    [
                        {"sourceId": "1", "deprecated": False, "source": "source1", "name": "name"},
                        {"sourceId": "2", "deprecated": False, "source": "source2", "name": "name"},
                    ],
                ]
            )
        )
        rec = get_preferred_drug_representation(api, "anything")
        assert rec["sourceId"] == "2"

    def test_prefers_newer_version(self):
        api = MagicMock(
            query=MagicMock(
                side_effect=[
                    [],
                    [
                        {
                            "sourceId": "2",
                            "deprecated": True,
                            "source": "source",
                            "name": "name",
                            "sourceIdVersion": "1",
                        },
                        {
                            "sourceId": "2",
                            "deprecated": True,
                            "source": "source",
                            "name": "name",
                            "sourceIdVersion": "2",
                        },
                    ],
                ]
            )
        )
        rec = get_preferred_drug_representation(api, "anything")
        assert rec["sourceIdVersion"] == "1"


class TestSubstituteSentenceTemplate:
    def test_multiple_diseases_no_matches(self):
        template = "{conditions:variant} is associated with {relevance} to {subject} in {conditions:disease} ({evidence})"
        relevance = {"displayName": "senitivity"}
        disease_matches = {"1"}
        diseases = [
            {"@class": "Disease", "@rid": "2", "displayName": "disease 1"},
            {"@class": "Disease", "@rid": "3", "displayName": "disease 2"},
        ]
        variants = [
            {
                "@class": "CategoryVariant",
                "displayName": "KRAS increased RNA expression",
                "@rid": "4",
            }
        ]
        subjects = [{"@class": "Therapy", "displayName": "some drug", "@rid": "5"}]
        sentence = substitute_sentence_template(
            template, diseases + variants, subjects, relevance, [], ["6", "7"], disease_matches
        )
        assert (
            sentence
            == f'KRAS increased RNA expression is associated with senitivity to some drug in other disease types (<a href="{GRAPHKB_GUI}/data/table?complex=eyJ0YXJnZXQiOiBbIjYiLCAiNyJdfQ%3D%3D&%40class=Statement" target="_blank" rel="noopener"></a>)'
        )

    def test_multiple_diseases_some_matches(self):
        template = "{conditions:variant} is associated with {relevance} to {subject} in {conditions:disease} ({evidence})"
        relevance = {"displayName": "senitivity"}
        disease_matches = {"1"}
        diseases = [
            {"@class": "Disease", "@rid": "2", "displayName": "disease 2"},
            {"@class": "Disease", "@rid": "1", "displayName": "disease 1"},
            {"@class": "Disease", "@rid": "3", "displayName": "disease 3"},
        ]
        variants = [
            {
                "@class": "CategoryVariant",
                "displayName": "KRAS increased RNA expression",
                "@rid": "4",
            }
        ]
        subjects = [{"@class": "Therapy", "displayName": "some drug", "@rid": "5"}]
        sentence = substitute_sentence_template(
            template, diseases + variants, subjects, relevance, [], ["6", "7"], disease_matches
        )
        assert (
            sentence
            == f'KRAS increased RNA expression is associated with senitivity to some drug in disease 1, and other disease types (<a href="{GRAPHKB_GUI}/data/table?complex=eyJ0YXJnZXQiOiBbIjYiLCAiNyJdfQ%3D%3D&%40class=Statement" target="_blank" rel="noopener"></a>)'
        )

    def test_multiple_diseases_only_matches(self):
        template = "{conditions:variant} is associated with {relevance} to {subject} in {conditions:disease} ({evidence})"
        relevance = {"displayName": "senitivity"}
        disease_matches = {"1", "2", "3"}
        diseases = [
            {"@class": "Disease", "@rid": "2", "displayName": "disease 2"},
            {"@class": "Disease", "@rid": "1", "displayName": "disease 1"},
            {"@class": "Disease", "@rid": "3", "displayName": "disease 3"},
        ]
        variants = [
            {
                "@class": "CategoryVariant",
                "displayName": "KRAS increased RNA expression",
                "@rid": "4",
            }
        ]
        subjects = [{"@class": "Therapy", "displayName": "some drug", "@rid": "5"}]
        sentence = substitute_sentence_template(
            template, diseases + variants, subjects, relevance, [], ["6", "7"], disease_matches
        )
        assert (
            sentence
            == f'KRAS increased RNA expression is associated with senitivity to some drug in disease 1, disease 2, and disease 3 (<a href="{GRAPHKB_GUI}/data/table?complex=eyJ0YXJnZXQiOiBbIjYiLCAiNyJdfQ%3D%3D&%40class=Statement" target="_blank" rel="noopener"></a>)'
        )


mock_ipr_results = [
    [
        {
            'text': '<p>no cancerType</p>',
            'variantName': 'ERBB2 amplification',
            'cancerType': [],
            'template': {'name': 'test3'},
            'project': {'name': 'test2'},
        },
        {
            'text': '<p>normal</p>',
            'variantName': 'ERBB2 amplification',
            'cancerType': ['test1', 'test'],
            'template': {'name': 'test3'},
            'project': {'name': 'test2'},
        },
        {
            'text': '<p>no project</p>',
            'variantName': 'ERBB2 amplification',
            'cancerType': ['test1', 'test'],
            'template': {'name': 'test3'},
        },
        {
            'text': '<p>no template</p>',
            'variantName': 'ERBB2 amplification',
            'cancerType': ['test1', 'test'],
            'project': {'name': 'test2'},
        },
    ],
    [
        {
            'text': '<p>normal, second variant</p>',
            'variantName': 'second variant',
            'cancerType': ['test1', 'test'],
            'template': {'name': 'test3'},
            'project': {'name': 'test2'},
        },
    ],
]

no_comments_found_output = 'No comments found in IPR for variants in this report'


class TestVariantTextFromIPR:
    def test_gets_fully_matched_output_when_possible(self):
        ipr_conn = MagicMock(get=MagicMock(side_effect=copy(mock_ipr_results)))
        matches = [{'kbVariant': 'ERBB2 amplification'}]
        ipr_summary = get_ipr_analyst_comments(
            ipr_conn,
            matches=matches,
            disease_name='test1',
            project_name='test2',
            report_type='test3',
            include_nonspecific_project=False,
            include_nonspecific_disease=True,
            include_nonspecific_template=True,
        )
        summary_lines = ipr_summary.split('\n')
        assert summary_lines[1] == '<h2>ERBB2 amplification (test1,test)</h2>'
        assert summary_lines[2] == '<p><p>normal</p></p>'
        assert len(summary_lines) == 3

    def test_omits_nonspecific_project_matches_when_specified(self):
        ipr_conn = MagicMock(get=MagicMock(side_effect=copy(mock_ipr_results)))
        matches = [{'kbVariant': 'ERBB2 amplification'}]
        ipr_summary = get_ipr_analyst_comments(
            ipr_conn,
            matches=matches,
            disease_name='test1',
            project_name='notfound',
            report_type='test3',
            include_nonspecific_project=False,
            include_nonspecific_disease=True,
            include_nonspecific_template=True,
        )
        assert ipr_summary == no_comments_found_output

    def test_omits_nonspecific_template_matches_when_specified(self):
        ipr_conn = MagicMock(get=MagicMock(side_effect=copy(mock_ipr_results)))
        matches = [{'kbVariant': 'ERBB2 amplification'}]
        ipr_summary = get_ipr_analyst_comments(
            ipr_conn,
            matches=matches,
            disease_name='test1',
            project_name='test2',
            report_type='notfound',
            include_nonspecific_project=True,
            include_nonspecific_disease=True,
            include_nonspecific_template=False,
        )
        assert ipr_summary == no_comments_found_output

    def test_omits_nonspecific_disease_matches_when_specified(self):
        ipr_conn = MagicMock(get=MagicMock(side_effect=copy(mock_ipr_results)))
        matches = [{'kbVariant': 'ERBB2 amplification'}]
        ipr_summary = get_ipr_analyst_comments(
            ipr_conn,
            matches=matches,
            disease_name='notfound',
            project_name='test2',
            report_type='test3',
            include_nonspecific_project=True,
            include_nonspecific_disease=False,
            include_nonspecific_template=True,
        )
        assert ipr_summary == no_comments_found_output

    def test_includes_nonspecific_project_matches_when_specified(self):
        ipr_conn = MagicMock(get=MagicMock(side_effect=copy(mock_ipr_results)))
        matches = [{'kbVariant': 'ERBB2 amplification'}]
        ipr_summary = get_ipr_analyst_comments(
            ipr_conn,
            matches=matches,
            disease_name='test1',
            project_name='notfound',
            report_type='test3',
            include_nonspecific_project=True,
            include_nonspecific_disease=False,
            include_nonspecific_template=False,
        )
        summary_lines = ipr_summary.split('\n')
        assert summary_lines[2] == '<p><p>no project</p></p>'
        assert len(summary_lines) == 3

    def test_includes_nonspecific_template_matches_when_specified(self):
        ipr_conn = MagicMock(get=MagicMock(side_effect=copy(mock_ipr_results)))
        matches = [{'kbVariant': 'ERBB2 amplification'}]
        ipr_summary = get_ipr_analyst_comments(
            ipr_conn,
            matches=matches,
            disease_name='test1',
            project_name='test2',
            report_type='notfound',
            include_nonspecific_project=False,
            include_nonspecific_disease=False,
            include_nonspecific_template=True,
        )
        summary_lines = ipr_summary.split('\n')
        assert summary_lines[2] == '<p><p>no template</p></p>'
        assert len(summary_lines) == 3

    def test_includes_nonspecific_disease_matches_when_specified(self):
        ipr_conn = MagicMock(get=MagicMock(side_effect=copy(mock_ipr_results)))
        matches = [{'kbVariant': 'ERBB2 amplification'}]
        ipr_summary = get_ipr_analyst_comments(
            ipr_conn,
            matches=matches,
            disease_name='notfound',
            project_name='test2',
            report_type='test3',
            include_nonspecific_project=False,
            include_nonspecific_disease=True,
            include_nonspecific_template=False,
        )
        summary_lines = ipr_summary.split('\n')
        assert summary_lines[1] == '<h2>ERBB2 amplification (no specific cancer types)</h2>'
        assert summary_lines[2] == '<p><p>no cancerType</p></p>'
        assert len(summary_lines) == 3

    def test_prepare_section_for_multiple_variants(self):
        ipr_conn = MagicMock(get=MagicMock(side_effect=copy(mock_ipr_results)))
        # NB this test relies on matches being processed in this order
        matches = [{'kbVariant': 'ERBB2 amplification'}, {'kbVariant': 'second variant'}]
        ipr_summary = get_ipr_analyst_comments(
            ipr_conn,
            matches=matches,
            disease_name='test1',
            project_name='test2',
            report_type='test3',
            include_nonspecific_project=False,
            include_nonspecific_disease=False,
            include_nonspecific_template=False,
        )
        summary_lines = ipr_summary.split('\n')
        assert len(summary_lines) == 5
        assert (
            '\n'.join(summary_lines[1:])
            == '<h2>ERBB2 amplification (test1,test)</h2>\n<p><p>normal</p></p>\n<h2>second variant (test1,test)</h2>\n<p><p>normal, second variant</p></p>'
        )

    def test_empty_section_when_no_variant_match(self):
        ipr_conn = MagicMock(get=MagicMock(side_effect=[[], []]))
        matches = [{'kbVariant': 'notfound1'}, {'kbVariant': 'notfound2'}]
        ipr_summary = get_ipr_analyst_comments(
            ipr_conn,
            matches=matches,
            disease_name='test1',
            project_name='test2',
            report_type='test3',
            include_nonspecific_project=False,
            include_nonspecific_disease=False,
            include_nonspecific_template=False,
        )
        assert ipr_summary == no_comments_found_output
