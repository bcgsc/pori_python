import pytest
from unittest.mock import Mock, patch

from pori_python.graphkb import statement as gkb_statement
from pori_python.graphkb import vocab as gkb_vocab
from pori_python.ipr.ipr import (
    convert_statements_to_alterations,
    germline_kb_matches,
    get_kb_matched_statements,
    get_kb_statement_matched_conditions,
    get_kb_variants,
)
from pori_python.types import Statement

DISEASE_RIDS = ["#138:12", "#138:13"]
APPROVED_EVIDENCE_RIDS = ["approved1", "approved2"]
GERMLINE_VARIANTS = [
    {
        "key": "1",
        "germline": True,
        "hgvsCds": "SLC28A3:c.1381C>T",
        "hgvsGenomic": "chr9:g.84286011G>A",
        "hgvsProtein": "SLC28A3:p.L461L",
        "ncbiBuild": "GRCh38",
        "normalAltCount": 37,
        "normalDepth": 37,
        "normalRefCount": 0,
        "proteinChange": "p.L461L",
        "rnaAltCount": "",
        "rnaDepth": "",
        "rnaRefCount": "",
        "startPosition": 84286011,
        "transcript": "ENST00000376238",
        "tumourAltCount": "",
        "tumourDepth": "",
        "tumourRefCount": "",
        "variant": "SLC28A3:p.L461L",
        "variantType": "mut",
        "zygosity": "",
    },
    {
        "key": "2",
        "germline": True,
        "hgvsCds": "BRCA1:c.4837A>",
        "hgvsGenomic": "chr17:g.43071077T>C",
        "hgvsProtein": "BRCA1:p.S1613G",
        "normalAltCount": 33,
        "normalDepth": 33,
        "normalRefCount": 0,
        "tumourAltCount": 37,
        "tumourDepth": 37,
        "tumourRefCount": 0,
    },
]

SOMATIC_VARIANTS = [
    {
        "key": "1",
        "gene": "SLC28A3",
        "germline": False,
        "hgvsCds": "SLC28A3:c.1381C>T",
        "hgvsGenomic": "chr9:g.84286011G>A",
        "hgvsProtein": "SLC28A3:p.L461L",
        "ncbiBuild": "GRCh38",
        "normalAltCount": 0,
        "normalDepth": 37,
        "normalRefCount": 37,
        "tumourAltCount": 37,
        "tumourDepth": 37,
        "tumourRefCount": 0,
        "variant": "SLC28A3:p.L461L",
        "variantType": "mut",
        "zygosity": "",
    },
    {
        "key": "2",
        "germline": False,
        "hgvsCds": "BRCA1:c.4837A>",
        "hgvsGenomic": "chr17:g.43071077T>C",
        "hgvsProtein": "BRCA1:p.S1613G",
        "normalAltCount": 1,
        "normalDepth": 33,
        "normalRefCount": 32,
        "tumourAltCount": 37,
        "tumourDepth": 37,
        "tumourRefCount": 0,
    },
]

GERMLINE_KB_MATCHES = [
    {
        "variant": "1",
        "approvedTherapy": False,
        "category": "pharmacogenomic",
        "context": "anthracyclines",
        "kbContextId": "#122:20944",
        "kbRelevanceId": "#147:38",
        "kbStatementId": "#154:13387",
        "kbVariant": "SLC28A3:c.1381C>T",
        "kbVariantId": "#159:5426",
        "matchedCancer": False,
        "reference": "PMID: 27197003",
        "relevance": "decreased toxicity",
        "reviewStatus": "initial",
    },
    {
        "variant": "2",
        "approvedTherapy": True,
        "category": "cancer predisposition",
        "kbContextId": "#135:8764",
        "kbRelevanceId": "#147:32",
        "kbStatementId": "#155:13511",
        "kbVariant": "BRCA1 mutation",
        "kbVariantId": "#161:938",
        "matchedCancer": False,
        "reference": "MOAlmanac FDA-56",
        "relevance": "therapy",
        "reviewStatus": None,
    },
]

SOMATIC_KB_MATCHES = [
    {
        "variant": "1",
        "approvedTherapy": False,
        "category": "prognostic",
        "kbContextId": "somatic_test",
        "kbRelevanceId": "#147:38",
        "kbStatementId": "#154:13387",
        "kbVariant": "SLC28A3:c.1381C>T",
        "kbVariantId": "#159:5426",
        "relevance": "prognostic",
        "reviewStatus": "initial",
    },
    {
        "variant": "2",
        "approvedTherapy": True,
        "category": "therapy",
        "kbContextId": "#135:8764",
        "kbRelevanceId": "#147:32",
        "kbStatementId": "#155:13511",
        "kbVariant": "BRCA1 mutation",
        "kbVariantId": "#161:938",
        "matchedCancer": False,
        "reference": "MOAlmanac FDA-56",
        "relevance": "therapy",
        "reviewStatus": None,
    },
]

KB_MATCHES_STATEMENTS = [
    {
        "@rid": SOMATIC_KB_MATCHES[0]["kbStatementId"],
        "conditions": [
            {
                "@class": "PositionalVariant",
                "@rid": SOMATIC_KB_MATCHES[0]["kbVariantId"],
            },
            {"@class": "CategoryVariant", "@rid": SOMATIC_KB_MATCHES[1]["kbVariantId"]},
            {"@class": "Disease", "@rid": ""},
        ],
    },
    {
        "@rid": SOMATIC_KB_MATCHES[1]["kbStatementId"],
        "conditions": [
            {"@class": "CategoryVariant", "@rid": SOMATIC_KB_MATCHES[1]["kbVariantId"]},
            {"@class": "PositionalVariant", "@rid": "157:0", "type": "#999:99"},
        ],
    },
]


@pytest.fixture
def graphkb_conn():
    class QueryMock:
        return_values = [
            # get approved evidence levels
            [{"@rid": v} for v in APPROVED_EVIDENCE_RIDS]
        ]
        index = -1

        def __call__(self, *args, **kwargs):
            self.index += 1
            ret_val = self.return_values[self.index] if self.index < len(self.return_values) else []
            return ret_val

    class PostMock:
        def __call__(self, *args, **kwargs):
            # custom return tailored for multi_variant_filtering() testing
            return {"result": KB_MATCHES_STATEMENTS}

    def mock_get_source(source):
        return {"@rid": 0}

    conn = Mock(query=QueryMock(), cache={}, get_source=mock_get_source, post=PostMock())

    return conn


def base_graphkb_statement(disease_id: str = "disease", relevance_rid: str = "other") -> Statement:
    statement = Statement(  # type: ignore
        {
            "conditions": [
                {
                    "@class": "Disease",
                    "@rid": disease_id,
                    "displayName": "disease_display_name",
                },
                {
                    "@class": "CategoryVariant",
                    "@rid": "variant_rid",
                    "displayName": "KRAS increased expression",
                },
            ],
            "evidence": [],
            "subject": {
                "@class": "dummy_value",
                "@rid": "101:010",
                "displayName": "dummy_display_name",
            },
            "source": None,
            "sourceId": None,
            "relevance": {
                "@rid": relevance_rid,
                "displayName": "relevance_display_name",
                "name": "relevance_name",
            },
            "@rid": "statement_rid",
        }
    )
    return statement


@pytest.fixture(autouse=True)
def mock_get_term_tree(monkeypatch):
    def mock_func(*pos, **kwargs):
        return [{"@rid": d} for d in DISEASE_RIDS]

    monkeypatch.setattr(gkb_vocab, "get_term_tree", mock_func)


@pytest.fixture(autouse=True)
def get_terms_set(monkeypatch):
    def mock_func(*pos, **kwargs):
        return {"#999:99"}

    monkeypatch.setattr(gkb_vocab, "get_terms_set", mock_func)


@pytest.fixture(autouse=True)
def mock_categorize_relevance(monkeypatch):
    def mock_func(_, relevance_id):
        return relevance_id

    monkeypatch.setattr(gkb_statement, "categorize_relevance", mock_func)


class TestConvertStatementsToAlterations:
    def test_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement(DISEASE_RIDS[0])
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], DISEASE_RIDS, {"variant_rid"}
        )

        assert len(result) == 1
        row = result[0]
        assert row["kbVariantId"] == "variant_rid"
        assert row["kbStatementId"] == "statement_rid"
        assert row["matchedCancer"]
        assert row["kbVariant"] == "KRAS increased expression"
        assert row["relevance"] == "relevance_display_name"

    def test_no_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement("other")
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], DISEASE_RIDS, {"variant_rid"}
        )

        assert len(result) == 1
        row = result[0]
        assert not row["matchedCancer"]

    def test_multiple_disease_not_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement("disease")
        statement["conditions"].append(
            {"@class": "Disease", "@rid": "other", "displayName": "disease_display_name"}  # type: ignore
        )
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], DISEASE_RIDS, {"variant_rid"}
        )

        assert len(result) == 1
        row = result[0]
        assert not row["matchedCancer"]

    def test_biological(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement["relevance"]["@rid"] = "biological"

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], DISEASE_RIDS, {"variant_rid"}
        )
        assert len(result) == 1
        row = result[0]
        assert row["category"] == "biological"

    def test_prognostic_no_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement["relevance"]["@rid"] = "prognostic"

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], DISEASE_RIDS, {"variant_rid"}
        )
        assert len(result) == 0

    def test_prognostic_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement(DISEASE_RIDS[0])
        statement["relevance"]["@rid"] = "prognostic"

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], DISEASE_RIDS, {"variant_rid"}
        )
        assert len(result) == 1
        row = result[0]
        assert row["category"] == "prognostic"

    def test_diagnostic(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement["relevance"]["@rid"] = "diagnostic"

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], DISEASE_RIDS, {"variant_rid"}
        )
        assert len(result) == 1
        row = result[0]
        assert row["category"] == "diagnostic"

    @patch("pori_python.ipr.ipr.get_evidencelevel_mapping")
    def test_unapproved_therapeutic(self, mock_get_evidencelevel_mapping, graphkb_conn) -> None:
        mock_get_evidencelevel_mapping.return_value = {"other": "test"}

        statement = base_graphkb_statement()
        statement["relevance"]["@rid"] = "therapeutic"
        statement["evidenceLevel"] = [{"@rid": "other", "displayName": "level"}]  # type: ignore

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], DISEASE_RIDS, {"variant_rid"}
        )
        assert len(result) == 1
        row = result[0]
        assert row["category"] == "therapeutic"

    @patch("pori_python.ipr.ipr.get_evidencelevel_mapping")
    def test_approved_therapeutic(self, mock_get_evidencelevel_mapping, graphkb_conn) -> None:
        mock_get_evidencelevel_mapping.return_value = {APPROVED_EVIDENCE_RIDS[0]: "test"}

        statement = base_graphkb_statement()
        statement["relevance"]["@rid"] = "therapeutic"
        statement["evidenceLevel"] = [{"@rid": APPROVED_EVIDENCE_RIDS[0], "displayName": "level"}]  # type: ignore

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], DISEASE_RIDS, {"variant_rid"}
        )
        assert len(result) == 1
        row = result[0]
        assert row["category"] == "therapeutic"


class TestKbmatchFilters:
    def test_germline_kb_matches(self):
        assert len(germline_kb_matches(GERMLINE_KB_MATCHES, GERMLINE_VARIANTS)) == len(
            GERMLINE_KB_MATCHES
        ), "Germline variant not matched to germline KB statement."
        assert not germline_kb_matches(
            GERMLINE_KB_MATCHES, SOMATIC_VARIANTS
        ), "Somatic variant matched to KB germline statement."
        assert len(germline_kb_matches(SOMATIC_KB_MATCHES, SOMATIC_VARIANTS)) == len(
            SOMATIC_KB_MATCHES
        ), "Somatic variant not matched to somatic KB statement."
        assert not germline_kb_matches(
            SOMATIC_KB_MATCHES, GERMLINE_VARIANTS
        ), "Germline variant matched to KB somatic statement."


GKB_MATCHES = [
    {
        "variant": "1",
        "approvedTherapy": False,
        "category": "prognostic",
        "kbContextId": "somatic_test",
        "kbRelevanceId": "#147:38",
        "kbStatementId": "#154:13387",
        "requiredKbMatches": ["#159:5426"],
        "kbVariant": "SLC28A3:c.1381C>T",
        "kbVariantId": "#159:5426",
        "relevance": "prognostic",
        "variantType": "mut",
        "reviewStatus": "initial",
    },
    {
        "variant": "2",
        "approvedTherapy": True,
        "category": "therapy",
        "kbContextId": "#135:8764",
        "kbRelevanceId": "#147:32",
        "kbStatementId": "#155:13511",
        "requiredKbMatches": ["#161:938"],
        "kbVariant": "BRCA1 mutation",
        "kbVariantId": "#161:938",
        "matchedCancer": False,
        "reference": "MOAlmanac FDA-56",
        "relevance": "therapy",
        "variantType": "mut",
        "reviewStatus": None,
    },
    {
        "variant": "3",
        "approvedTherapy": True,
        "category": "therapy",
        "kbContextId": "#135:8764",
        "kbRelevanceId": "#147:32",
        "kbStatementId": "#155:13511",
        "requiredKbMatches": ["#161:938"],
        "kbVariant": "BRCA1 mutation",
        "kbVariantId": "#161:938",
        "matchedCancer": False,
        "reference": "MOAlmanac FDA-56",
        "relevance": "therapy",
        "variantType": "mut",
        "reviewStatus": None,
    },
    {
        "variant": "4",
        "approvedTherapy": True,
        "category": "therapy",
        "kbContextId": "#135:8764",
        "kbRelevanceId": "#147:32",
        "kbStatementId": "#155:13511",
        "requiredKbMatches": ["#159:5426", "#161:938"],
        "kbVariant": "BRCA1 mutation",
        "kbVariantId": "#161:938",
        "matchedCancer": False,
        "reference": "MOAlmanac FDA-56",
        "relevance": "therapy",
        "variantType": "mut",
        "reviewStatus": None,
    },
]

BASIC_GKB_MATCH = {
    "approvedTherapy": False,
    "category": "test",
    "context": "test",
    "kbContextId": "#124:24761",
    "disease": "test",
    "evidenceLevel": "test",
    "iprEvidenceLevel": "test",
    "matchedCancer": False,
    "reference": "test",
    "relevance": "test",
    "kbRelevanceId": "#148:31",
    "externalSource": "",
    "externalStatementId": "",
    "reviewStatus": "passed",
    "kbData": {},
}


def create_gkb_matches(input_fields):
    matches = []
    for item in input_fields:
        temp = BASIC_GKB_MATCH.copy()
        temp.update(item)
        matches.append(temp)
    return matches


def get_condition_set_string_rep(condition_set):
    for item in condition_set:
        item["observedKeysStrs"] = [
            "-".join([elem["kbVariantId"], elem["observedVariantKey"]])
            for elem in item["matchedConditions"]
        ]
        item["observedKeysStrs"].sort()
        item["observedKeysStr"] = ",".join(item["observedKeysStrs"])
    condition_set = [f"{item['kbStatementId']},{item['observedKeysStr']}" for item in condition_set]
    condition_set.sort()
    return condition_set


class TestKbMatchSectionPrep:
    def test_matched_variant_pairs_extracted_only_once_for_multiple_statements(self):
        input_fields = [
            {"variant": "A", "kbVariantId": "test1", "kbStatementId": "test1"},
            {
                "variant": "A",
                "kbVariantId": "test1",
                "kbStatementId": "test2",
            },  # diff statement
        ]
        for item in input_fields:  # we don't care about these for this test
            item["variantType"] = "test"
            item["kbVariant"] = "test"
            item["requiredKbMatches"] = ["test1", "test2"]
        gkb_matches = create_gkb_matches(input_fields)
        kb_variants = get_kb_variants(gkb_matches)
        found_variants = [f"{item['variant']},{item['kbVariantId']}" for item in kb_variants]
        found_variants.sort()
        assert found_variants == ["A,test1"]

    def test_all_distinct_observed_and_matched_variant_pairs_extracted(self):
        input_fields = [
            {"variant": "A", "kbVariantId": "test1", "kbStatementId": "test1"},
            {"variant": "A", "kbVariantId": "test2", "kbStatementId": "test1"},
            {"variant": "B", "kbVariantId": "test1", "kbStatementId": "test1"},
            {"variant": "B", "kbVariantId": "test1", "kbStatementId": "test1"},
        ]
        for item in input_fields:  # we don't care about these for this test
            item["variantType"] = "test"
            item["kbVariant"] = "test"
            item["requiredKbMatches"] = ["test1", "test2"]
        gkb_matches = create_gkb_matches(input_fields)
        kb_variants = get_kb_variants(gkb_matches)
        found_variants = [f"{item['variant']},{item['kbVariantId']}" for item in kb_variants]
        found_variants.sort()
        assert found_variants == ["A,test1", "A,test2", "B,test1"]

    def test_statements_extracted_only_once(self):
        input_fields = [
            {"variant": "A", "kbVariantId": "test1", "kbStatementId": "X"},
            {"variant": "A", "kbVariantId": "test2", "kbStatementId": "X"},
            {"variant": "B", "kbVariantId": "test1", "kbStatementId": "X"},
            {"variant": "B", "kbVariantId": "test2", "kbStatementId": "X"},
            {"variant": "C", "kbVariantId": "test1", "kbStatementId": "Y"},
        ]
        for item in input_fields:  # we don't care about these for this test
            item["variantType"] = "test"
            item["kbVariant"] = "test"
            item["requiredKbMatches"] = ["test1", "test2"]
        gkb_matches = create_gkb_matches(input_fields)
        kb_stmts = get_kb_matched_statements(gkb_matches)
        kb_stmts = [item["kbStatementId"] for item in kb_stmts]
        kb_stmts.sort()
        assert kb_stmts == ["X", "Y"]

    def test_singlevar_statements_with_multiple_satisfying_condition_sets(self):
        input_fields = [
            {"variant": "A", "kbVariantId": "test1", "kbStatementId": "X"},
            {"variant": "B", "kbVariantId": "test1", "kbStatementId": "X"},
            {"variant": "C", "kbVariantId": "test1", "kbStatementId": "X"},
        ]
        for item in input_fields:  # we don't care about these for this test
            item["variantType"] = "test"
            item["kbVariant"] = "test"
            item["requiredKbMatches"] = ["test1"]
        gkb_matches = create_gkb_matches(input_fields)
        kbcs = get_kb_statement_matched_conditions(gkb_matches)
        kbcs_string_rep = get_condition_set_string_rep(kbcs)
        assert kbcs_string_rep == ["X,test1-A", "X,test1-B", "X,test1-C"]

    def test_multivar_statements_with_multiple_satisfying_condition_sets(self):
        input_fields = [
            {"variant": "A", "kbVariantId": "test1", "kbStatementId": "X"},
            {"variant": "B", "kbVariantId": "test2", "kbStatementId": "X"},
            {"variant": "C", "kbVariantId": "test2", "kbStatementId": "X"},
        ]
        for item in input_fields:  # we don't care about these for this test
            item["variantType"] = "test"
            item["kbVariant"] = "test"
            item["requiredKbMatches"] = ["test1", "test2"]
        gkb_matches = create_gkb_matches(input_fields)
        kbcs = get_kb_statement_matched_conditions(gkb_matches)
        kbcs_string_rep = get_condition_set_string_rep(kbcs)
        assert kbcs_string_rep == ["X,test1-A,test2-B", "X,test1-A,test2-C"]

    def test_do_not_infer_possible_matches(self):
        """edge case - when infer_possible_matches is false, do not allow var/kbvar
        pairs to satisfy conditions for statements they are not explicitly linked
        to in the input"""
        input_fields = [
            {"variant": "A", "kbVariantId": "test1", "kbStatementId": "X"},
            {"variant": "B", "kbVariantId": "test1", "kbStatementId": "Y"},
        ]
        for item in input_fields:  # we don't care about these for this test
            item["variantType"] = "test"
            item["kbVariant"] = "test"
            item["requiredKbMatches"] = ["test1"]
        gkb_matches = create_gkb_matches(input_fields)
        kbcs = get_kb_statement_matched_conditions(gkb_matches)
        kbcs_string_rep = get_condition_set_string_rep(kbcs)
        assert kbcs_string_rep == ["X,test1-A", "Y,test1-B"]

    def test_no_dupes_when_requiredKbMatches_not_sorted(self):
        input_fields = [
            {
                "variant": "A",
                "kbVariantId": "test1",
                "requiredKbMatches": ["test1", "test2"],
            },
            {
                "variant": "B",
                "kbVariantId": "test2",
                "requiredKbMatches": ["test1", "test2"],
            },
            {
                "variant": "A",
                "kbVariantId": "test1",
                "requiredKbMatches": ["test2", "test1"],
            },
            {
                "variant": "B",
                "kbVariantId": "test2",
                "requiredKbMatches": ["test2", "test1"],
            },
        ]
        for item in input_fields:  # we don't care about these for this test
            item["variantType"] = "test"
            item["kbVariant"] = "test"
            item["kbStatementId"] = "X"
        gkb_matches = create_gkb_matches(input_fields)
        stmts = get_kb_matched_statements(gkb_matches)
        kbcs = get_kb_statement_matched_conditions(gkb_matches)
        assert len(stmts) == 1
        assert len(kbcs) == 1

    def test_partial_matches_omitted(self):
        """check statements that are only partially supported
        are omitted when allow_partial_matches=False"""
        input_fields = [
            {
                "variant": "A",
                "kbVariantId": "test1",
                "kbStatementId": "X",
                "requiredKbMatches": ["test1", "test2"],
            },
            {
                "variant": "B",
                "kbVariantId": "test2",
                "kbStatementId": "X",
                "requiredKbMatches": ["test1", "test2"],
            },
            {
                "variant": "A",
                "kbVariantId": "test1",
                "kbStatementId": "Y",
                "requiredKbMatches": ["test1", "test3"],
            },
        ]
        for item in input_fields:  # we don't care about these for this test
            item["variantType"] = "test"
            item["kbVariant"] = "test"
        gkb_matches = create_gkb_matches(input_fields)
        stmts = get_kb_matched_statements(gkb_matches)
        kbcs = get_kb_statement_matched_conditions(gkb_matches)
        assert len(stmts) == 2
        assert len(kbcs) == 1  # X only
        assert kbcs[0]["kbStatementId"] == "X"

    def test_partial_matches_omitted_even_when_var_used_elsewhere(self):
        """edge case -
        checks that vars that satisfy other conditions, but aren't explicitly used
        to satisfy conditions for some statement in the input,
        are not used in the satisfying condition sets for the statement
        so that it shows up in the results when it otherwise wouldn't.
        Eg here for statement Y, requirement test1 is satisfied
        but requirement test3 is not considered satisfied, even though it is
        satisfied for statement Z"""
        input_fields = [
            {
                "variant": "A",
                "kbVariantId": "test1",
                "kbStatementId": "X",
                "requiredKbMatches": ["test1", "test2"],
            },
            {
                "variant": "B",
                "kbVariantId": "test2",
                "kbStatementId": "X",
                "requiredKbMatches": ["test1", "test2"],
            },
            {
                "variant": "A",
                "kbVariantId": "test1",
                "kbStatementId": "Y",
                "requiredKbMatches": ["test1", "test3"],
            },
            {
                "variant": "C",
                "kbVariantId": "test3",
                "kbStatementId": "Z",
                "requiredKbMatches": ["test3"],
            },
        ]
        for item in input_fields:  # we don't care about these for this test
            item["variantType"] = "test"
            item["kbVariant"] = "test"
        gkb_matches = create_gkb_matches(input_fields)
        stmts = get_kb_matched_statements(gkb_matches)
        kbcs = get_kb_statement_matched_conditions(gkb_matches)
        assert len(stmts) == 3
        assert len(kbcs) == 2  # X and Z but not Y
        assert "Y" not in [item["kbStatementId"] for item in kbcs]

    def test_partial_matches_included(self):
        """check statements that are only partially supported
        are included when allow_partial_matches=True"""
        input_fields = [
            {
                "variant": "A",
                "kbVariantId": "test1",
                "kbStatementId": "X",
                "requiredKbMatches": ["test1", "test2"],
            },
            {
                "variant": "B",
                "kbVariantId": "test2",
                "kbStatementId": "X",
                "requiredKbMatches": ["test1", "test2"],
            },
            {
                "variant": "A",
                "kbVariantId": "test1",
                "kbStatementId": "Y",
                "requiredKbMatches": ["test1", "test3"],
            },
        ]
        for item in input_fields:  # we don't care about these for this test
            item["variantType"] = "test"
            item["kbVariant"] = "test"
        gkb_matches = create_gkb_matches(input_fields)
        stmts = get_kb_matched_statements(gkb_matches)
        kbcs = get_kb_statement_matched_conditions(gkb_matches, allow_partial_matches=True)
        assert len(stmts) == 2  # X and Y
        assert len(kbcs) == 2
