import json
import os
import pandas as pd
import pytest
import sys
import uuid
from typing import Dict, Tuple, List
from unittest.mock import MagicMock, patch

from pori_python.ipr.connection import IprConnection
from pori_python.ipr.main import command_interface
from pori_python.ipr.types import IprGene

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
def loaded_report(tmp_path_factory) -> Dict:
    mock = MagicMock()
    json_file = tmp_path_factory.mktemp("inputs") / "content.json"
    patient_id = f"TEST_{str(uuid.uuid4())}"
    json_file.write_text(
        json.dumps(
            {
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
                "expressionVariants": json.loads(
                    pd.read_csv(get_test_file("expression.short.tab"), sep="\t").to_json(
                        orient="records"
                    )
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
            },
            allow_nan=False,
        )
    )
    with patch.object(
        sys,
        "argv",
        [
            "ipr",
            "--username",
            os.environ.get("IPR_USER", os.environ["USER"]),
            "--password",
            os.environ["IPR_PASS"],
            "--ipr_url",
            os.environ["IPR_TEST_URL"],
            "--graphkb_url",
            os.environ.get("GRAPHKB_URL", False),
            "--content",
            str(json_file),
            "--therapeutics",
        ],
    ):
        with patch.object(IprConnection, "get_spec", return_value=get_test_spec()):
            command_interface()

    ipr_conn = IprConnection(
        username=os.environ.get("IPR_USER", os.environ["USER"]),
        password=os.environ["IPR_PASS"],
        url=os.environ["IPR_TEST_URL"],
    )
    loaded_report = ipr_conn.get(uri=f"reports?searchText={patient_id}")
    return (patient_id, loaded_report)


def get_section(loaded_report, section_name):
    ident = loaded_report[1]["reports"][0]["ident"]
    ipr_conn = IprConnection(
        username=os.environ.get("IPR_USER", os.environ["USER"]),
        password=os.environ["IPR_PASS"],
        url=os.environ["IPR_TEST_URL"],
    )
    return ipr_conn.get(uri=f"reports/{ident}/{section_name}")


@pytest.mark.skipif(
    not INCLUDE_UPLOAD_TESTS, reason="excluding tests of upload to live ipr instance"
)
@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
class TestCreateReport:
    def test_patient_id_loaded_once(self, loaded_report: Tuple) -> None:
        patient_id = loaded_report[0]
        assert loaded_report[1]["total"] == 1
        assert loaded_report[1]["reports"][0]["patientId"] == patient_id

    def test_analyst_comments_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, "summary/analyst-comments")
        assert section["comments"]
