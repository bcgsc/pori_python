import json
import os
import pandas as pd
import pytest
import sys
import uuid
from typing import Generator
from unittest.mock import patch

from pori_python.ipr.connection import IprConnection
from pori_python.ipr.main import command_interface
from pori_python.types import IprGene

from .constants import EXCLUDE_INTEGRATION_TESTS

EXCLUDE_BCGSC_TESTS = os.environ.get("EXCLUDE_BCGSC_TESTS") == "1"
EXCLUDE_ONCOKB_TESTS = os.environ.get("EXCLUDE_ONCOKB_TESTS") == "1"
INCLUDE_UPLOAD_TESTS = os.environ.get("INCLUDE_UPLOAD_TESTS", 0) == "1"


def get_test_spec():
    ipr_spec = {"components": {"schemas": {"genesCreate": {"properties": {}}}}}
    ipr_gene_keys = IprGene.__required_keys__ | IprGene.__optional_keys__
    for key in ipr_gene_keys:
        ipr_spec["components"]["schemas"]["genesCreate"]["properties"][key] = ""
    return ipr_spec


def get_test_file(name: str) -> str:
    return os.path.join(os.path.dirname(__file__), "test_data", name)


@pytest.fixture(scope="module")
def loaded_reports(tmp_path_factory) -> Generator:
    json_file = tmp_path_factory.mktemp("inputs") / "content.json"
    async_json_file = tmp_path_factory.mktemp("inputs") / "async_content.json"
    patient_id = f"TEST_{str(uuid.uuid4())}"
    async_patient_id = f"TEST_ASYNC_{str(uuid.uuid4())}"
    json_contents = {
        "comparators": [
            {"analysisRole": "expression (disease)", "name": "1"},
            {"analysisRole": "expression (primary site)", "name": "2"},
            {"analysisRole": "expression (biopsy site)", "name": "3"},
            {
                "analysisRole": "expression (internal pancancer cohort)",
                "name": "4",
            },
        ],
        "patientId": patient_id,
        "project": "TEST",
        "sampleInfo": [
            {
                "sample": "Constitutional",
                "biopsySite": "Normal tissue",
                "sampleName": "SAMPLE1-PB",
                "primarySite": "Blood-Peripheral",
                "collectionDate": "11-11-11",
            },
            {
                "sample": "Tumour",
                "pathoTc": "90%",
                "biopsySite": "hepatic",
                "sampleName": "SAMPLE2-FF-1",
                "primarySite": "Vena Cava-Hepatic",
                "collectionDate": "12-12-12",
            },
        ],
        "expressionVariants": json.loads(
            pd.read_csv(get_test_file("expression.short.tab"), sep="\t").to_json(orient="records")
        ),
        "smallMutations": json.loads(
            pd.read_csv(get_test_file("small_mutations.short.tab"), sep="\t").to_json(
                orient="records"
            )
        ),
        "copyVariants": json.loads(
            pd.read_csv(get_test_file("copy_variants.short.tab"), sep="\t").to_json(
                orient="records"
            )
        ),
        "structuralVariants": json.loads(
            pd.read_csv(get_test_file("fusions.tab"), sep="\t").to_json(orient="records")
        ),
        "kbDiseaseMatch": "colorectal cancer",
        "cosmicSignatures": json.loads(
            pd.read_csv(get_test_file("cosmic_variants.tab"), sep="\t").to_json(orient="records")
        ),
        "hlaTypes": json.loads(
            pd.read_csv(get_test_file("hla_variants.tab"), sep="\t").to_json(orient="records")
        ),
        "images": [
            {
                "key": "cnvLoh.circos",
                "path": "test/testData/images/cnvLoh.png",
                "caption": "Test adding a caption to an image",
            }
        ],
    }
    json_file.write_text(
        json.dumps(
            json_contents,
            allow_nan=False,
        )
    )

    json_contents["patientId"] = async_patient_id
    async_json_file.write_text(
        json.dumps(
            json_contents,
            allow_nan=False,
        )
    )

    argslist = [
        "ipr",
        "--username",
        os.environ.get("IPR_USER", os.environ["USER"]),
        "--password",
        os.environ["IPR_PASS"],
        "--graphkb_username",
        os.environ.get("GRAPHKB_USER", os.environ.get("IPR_USER", os.environ["USER"])),
        "--graphkb_password",
        os.environ.get("GRAPHKB_PASS", os.environ["IPR_PASS"]),
        "--ipr_url",
        os.environ["IPR_TEST_URL"],
        "--graphkb_url",
        os.environ.get("GRAPHKB_URL", False),
        "--therapeutics",
        "--allow_partial_matches",
    ]

    sync_argslist = argslist.copy()
    sync_argslist.extend(["--content", str(json_file)])
    with patch.object(sys, "argv", sync_argslist):
        with patch.object(IprConnection, "get_spec", return_value=get_test_spec()):
            command_interface()

    async_argslist = argslist.copy()
    async_argslist.extend(["--content", str(async_json_file), "--async_upload"])
    with patch.object(sys, "argv", async_argslist):
        with patch.object(IprConnection, "get_spec", return_value=get_test_spec()):
            command_interface()

    ipr_conn = IprConnection(
        username=os.environ.get("IPR_USER", os.environ["USER"]),
        password=os.environ["IPR_PASS"],
        url=os.environ["IPR_TEST_URL"],
    )
    loaded_report = ipr_conn.get(uri=f"reports?searchText={patient_id}")
    async_loaded_report = ipr_conn.get(uri=f"reports?searchText={async_patient_id}")

    loaded_reports_result = {
        "sync": (patient_id, loaded_report),
        "async": (async_patient_id, async_loaded_report),
    }
    yield loaded_reports_result
    return
    ipr_conn.delete(uri=f"reports/{loaded_report['reports'][0]['ident']}")
    ipr_conn.delete(uri=f"reports/{async_loaded_report['reports'][0]['ident']}")


def get_section(loaded_report, section_name):
    ident = loaded_report[1]["reports"][0]["ident"]
    ipr_conn = IprConnection(
        username=os.environ.get("IPR_USER", os.environ["USER"]),
        password=os.environ["IPR_PASS"],
        url=os.environ["IPR_TEST_URL"],
    )
    return ipr_conn.get(uri=f"reports/{ident}/{section_name}")


def stringify_sorted(obj):
    """
    stringifies a (json) object
    in such a way that it can be compared for equality
    with another json object"""
    if isinstance(obj, list):
        obj = [stringify_sorted(item) for item in obj]
        obj.sort()
        return str(obj)
    elif isinstance(obj, dict):
        for key in ("ident", "updatedAt", "createdAt", "deletedAt"):
            obj.pop(key, None)
        keys = obj.keys()
        for key in keys:
            if isinstance(obj[key], list):
                obj[key] = stringify_sorted(obj[key])
            elif isinstance(obj[key], dict):
                obj[key] = stringify_sorted(obj[key])
        return str(obj)
    elif isinstance(obj, str):
        return obj
    else:
        return str(obj)


@pytest.mark.skipif(
    not INCLUDE_UPLOAD_TESTS, reason="excluding tests of upload to live ipr instance"
)
@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
class TestCreateReport:
    def test_patient_id_loaded_once(self, loaded_reports) -> None:
        sync_patient_id = loaded_reports["sync"][0]
        assert loaded_reports["sync"][1]["total"] == 1
        assert loaded_reports["sync"][1]["reports"][0]["patientId"] == sync_patient_id
        async_patient_id = loaded_reports["async"][0]
        assert loaded_reports["async"][1]["total"] == 1
        assert loaded_reports["async"][1]["reports"][0]["patientId"] == async_patient_id

    def test_expression_variants_loaded(self, loaded_reports) -> None:
        section = get_section(loaded_reports["sync"], "expression-variants")
        kbmatched = [item for item in section if item["kbMatches"]]
        assert "PTP4A3" in [item["gene"]["name"] for item in kbmatched]
        async_section = get_section(loaded_reports["async"], "expression-variants")
        async_equals_sync = stringify_sorted(section) == stringify_sorted(async_section)
        assert async_equals_sync

    def test_structural_variants_loaded(self, loaded_reports) -> None:
        section = get_section(loaded_reports["sync"], "structural-variants")
        kbmatched = [item for item in section if item["kbMatches"]]
        assert "(EWSR1,FLI1):fusion(e.7,e.4)" in [item["displayName"] for item in kbmatched]
        async_section = get_section(loaded_reports["async"], "structural-variants")
        async_equals_sync = stringify_sorted(section) == stringify_sorted(async_section)
        assert async_equals_sync

    def test_small_mutations_loaded(self, loaded_reports) -> None:
        section = get_section(loaded_reports["sync"], "small-mutations")
        kbmatched = [item for item in section if item["kbMatches"]]
        assert "FGFR2:p.R421C" in [item["displayName"] for item in kbmatched]
        assert "CDKN2A:p.T18M" in [item["displayName"] for item in kbmatched]
        async_section = get_section(loaded_reports["async"], "small-mutations")
        async_equals_sync = stringify_sorted(section) == stringify_sorted(async_section)
        assert async_equals_sync

    def test_copy_variants_loaded(self, loaded_reports) -> None:
        section = get_section(loaded_reports["sync"], "copy-variants")
        kbmatched = [item for item in section if item["kbMatches"]]
        assert ("ERBB2", "amplification") in [
            (item["gene"]["name"], item["displayName"]) for item in kbmatched
        ]
        async_section = get_section(loaded_reports["async"], "copy-variants")
        async_equals_sync = stringify_sorted(section) == stringify_sorted(async_section)
        assert async_equals_sync

    # # Uncomment when signatureVariants are supported in pori_ipr_api
    # def test_signature_variants_loaded(self, loaded_reports) -> None:
    #     section = get_section(loaded_reports["sync"], "signature-variants")
    #     kbmatched = [item for item in section if item["kbMatches"]]
    #     assert ("SBS2", "high signature") in [
    #         (item["signatureName"], item["variantTypeName"]) for item in kbmatched
    #     ]
    #     async_section = get_section(loaded_reports["async"], "signature-variants")
    #     assert compare_sections(section, async_section)

    def test_kb_matches_loaded(self, loaded_reports) -> None:
        section = get_section(loaded_reports["sync"], "kb-matches")
        observed_and_matched = set(
            [(item["kbVariant"], item["variant"]["displayName"]) for item in section]
        )
        for pair in [
            ("ERBB2 amplification", "amplification"),
            ("FGFR2 mutation", "FGFR2:p.R421C"),
            ("PTP4A3 overexpression", "increased expression"),
            ("EWSR1 and FLI1 fusion", "(EWSR1,FLI1):fusion(e.7,e.4)"),
            ("CDKN2A mutation", "CDKN2A:p.T18M"),
        ]:
            assert pair in observed_and_matched
        async_section = get_section(loaded_reports["async"], "kb-matches")
        async_equals_sync = stringify_sorted(section) == stringify_sorted(async_section)
        assert async_equals_sync

    def test_therapeutic_targets_loaded(self, loaded_reports) -> None:
        section = get_section(loaded_reports["sync"], "therapeutic-targets")
        therapeutic_target_genes = set([item["gene"] for item in section])
        for gene in ["CDKN2A", "ERBB2", "FGFR2", "PTP4A3"]:
            assert gene in therapeutic_target_genes
        async_section = get_section(loaded_reports["async"], "therapeutic-targets")
        async_equals_sync = stringify_sorted(section) == stringify_sorted(async_section)
        assert async_equals_sync

    def test_genomic_alterations_identified_loaded(self, loaded_reports) -> None:
        section = get_section(loaded_reports["sync"], "summary/genomic-alterations-identified")
        variants = set([item["geneVariant"] for item in section])
        for variant in [
            "FGFR2:p.R421C",
            "PTP4A3 (high_percentile)",
            "ERBB2 (Amplification)",
            "(EWSR1,FLI1):fusion(e.7,e.4)",
            "CDKN2A:p.T18M",
        ]:
            assert variant in variants
        async_section = get_section(
            loaded_reports["async"], "summary/genomic-alterations-identified"
        )
        async_equals_sync = stringify_sorted(section) == stringify_sorted(async_section)
        assert async_equals_sync

    def test_analyst_comments_loaded(self, loaded_reports) -> None:
        sync_section = get_section(loaded_reports["sync"], "summary/analyst-comments")
        assert sync_section["comments"]
        async_section = get_section(loaded_reports["async"], "summary/analyst-comments")
        assert async_section["comments"]
        assert sync_section["comments"] == async_section["comments"]

    def test_sample_info_loaded(self, loaded_reports) -> None:
        sync_section = get_section(loaded_reports["sync"], "sample-info")
        async_section = get_section(loaded_reports["async"], "sample-info")
        async_equals_sync = stringify_sorted(sync_section) == stringify_sorted(async_section)
        assert async_equals_sync

    def test_multivariant_multiconditionset_statements_loaded(self, loaded_reports) -> None:
        """
        Checks that multivariant statements and multiple condition sets prepared correctly
        by this package are handled as expected by the api.

        This test depends on the presence of a record for pmid:27302369 in the graphkb.
        This statement has three required variants. Two are present in the test data
        for test_upload.py. allow_partial_matches is passed as an arg so
        that this statement gets matched even though only 2/3 requirements
        are met.
        This is also a test of multiple condition sets since there are two variants
        in the test data that satisfy one of the conditions (the APC mutation)."""
        section = get_section(loaded_reports["sync"], "kb-matches/kb-matched-statements")
        multivariant_stmts = [item for item in section if item["reference"] == "pmid:27302369"]

        # if this statement is entered more than once there may be multiple sets of records to
        # check, so to make sure the count checks work, go stmt_id by stmt_id:
        stmt_ids = list(set([item["kbStatementId"] for item in multivariant_stmts]))
        for stmt_id in stmt_ids:
            stmts = [item for item in multivariant_stmts if item["kbStatementId"] == stmt_id]

            # we expect two stmts, one for each condition set
            assert len(stmts) == 2

            # we expect each condition set to have two kb variants in it
            # we expect the two kb variants to be the same in each stmt
            assert len(stmts[0]["kbMatches"]) == 2
            assert len(stmts[1]["kbMatches"]) == 2
            kbmatches1 = [item["kbVariant"] for item in stmts[0]["kbMatches"]]
            kbmatches2 = [item["kbVariant"] for item in stmts[1]["kbMatches"]]
            kbmatches1.sort()
            kbmatches2.sort()
            assert kbmatches1 == kbmatches2 == ["APC mutation", "KRAS mutation"]

            # we expect the two stmts to have different observed variant sets
            observedVariants1 = [item["variant"]["ident"] for item in stmts[0]["kbMatches"]]
            observedVariants2 = [item["variant"]["ident"] for item in stmts[1]["kbMatches"]]
            observedVariants1.sort()
            observedVariants2.sort()
            assert observedVariants1 != observedVariants2

            # we expect the two observed variant sets to have one element in common
            # (the kras mutation)
            assert len(observedVariants1) == len(observedVariants2) == 2
            assert len((set(observedVariants1 + observedVariants2))) == 3
