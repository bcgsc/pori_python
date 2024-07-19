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
def loaded_reports(tmp_path_factory) -> Dict:
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
        "--ipr_url",
        os.environ["IPR_TEST_URL"],
        "--graphkb_url",
        os.environ.get("GRAPHKB_URL", False),
        "--therapeutics",
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

    return {
        "sync": (patient_id, loaded_report),
        "async": (async_patient_id, async_loaded_report),
    }


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
    def test_patient_id_loaded_once(self, loaded_reports) -> None:
        sync_patient_id = loaded_reports["sync"][0]
        assert loaded_reports["sync"][1]["total"] == 1
        assert loaded_reports["sync"][1]["reports"][0]["patientId"] == sync_patient_id
        async_patient_id = loaded_reports["async"][0]
        assert loaded_reports["async"][1]["total"] == 1
        assert loaded_reports["async"][1]["reports"][0]["patientId"] == async_patient_id

    def test_analyst_comments_loaded(self, loaded_reports) -> None:
        sync_section = get_section(loaded_reports["sync"], "summary/analyst-comments")
        assert sync_section["comments"]
        async_section = get_section(loaded_reports["async"], "summary/analyst-comments")
        assert async_section["comments"]
