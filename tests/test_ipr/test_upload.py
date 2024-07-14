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
    ipr_spec = {'components': {'schemas': {'genesCreate': {'properties': {}}}}}
    ipr_gene_keys = IprGene.__required_keys__ | IprGene.__optional_keys__
    for key in ipr_gene_keys:
        ipr_spec['components']['schemas']['genesCreate']['properties'][key] = ""
    return ipr_spec


def get_test_file(name: str) -> str:
    return os.path.join(os.path.dirname(__file__), 'test_data', name)


@pytest.fixture(scope='module')
def loaded_report(tmp_path_factory) -> Dict:
    mock = MagicMock()
    json_file = tmp_path_factory.mktemp('inputs') / 'content.json'
    patient_id = f'TEST_{str(uuid.uuid4())}'
    json_file.write_text(
        json.dumps(
            {
                'comparators': [
                    {'analysisRole': 'expression (disease)', 'name': '1'},
                    {'analysisRole': 'expression (primary site)', 'name': '2'},
                    {'analysisRole': 'expression (biopsy site)', 'name': '3'},
                    {'analysisRole': 'expression (internal pancancer cohort)', 'name': '4'},
                ],
                'patientId': patient_id,
                'project': 'TEST',
                'expressionVariants': json.loads(
                    pd.read_csv(get_test_file('expression.short.tab'), sep='\t').to_json(
                        orient='records'
                    )
                ),
                'smallMutations': json.loads(
                    pd.read_csv(get_test_file('small_mutations.short.tab'), sep='\t').to_json(
                        orient='records'
                    )
                ),
                'copyVariants': json.loads(
                    pd.read_csv(get_test_file('copy_variants.short.tab'), sep='\t').to_json(
                        orient='records'
                    )
                ),
                'structuralVariants': json.loads(
                    pd.read_csv(get_test_file('fusions.tab'), sep='\t').to_json(orient='records')
                ),
                'kbDiseaseMatch': 'colorectal cancer',
            },
            allow_nan=False,
        )
    )
    with patch.object(
        sys,
        'argv',
        [
            'ipr',
            '--username',
            os.environ.get('IPR_USER', os.environ['USER']),
            '--password',
            os.environ['IPR_PASS'],
            '--ipr_url',
            os.environ['IPR_TEST_URL'],
            '--graphkb_username',
            os.environ.get('GRAPHKB_USER', os.environ['USER']),
            '--graphkb_password',
            os.environ.get('GRAPHKB_PASS', os.environ['IPR_PASS']),
            '--graphkb_url',
            os.environ.get('GRAPHKB_URL', False),
            '--content',
            str(json_file),
            '--therapeutics',
        ],
    ):
        with patch.object(IprConnection, 'get_spec', return_value=get_test_spec()):
            command_interface()

    ipr_conn = IprConnection(
        username=os.environ.get('IPR_USER', os.environ['USER']),
        password=os.environ['IPR_PASS'],
        url=os.environ['IPR_TEST_URL'],
    )
    loaded_report = ipr_conn.get(uri=f"reports?searchText={patient_id}")
    yield (patient_id, loaded_report)

    report_ident = loaded_report['reports'][0]['ident']
    ipr_conn.delete(uri=f'reports/{report_ident}')


def get_section(loaded_report, section_name):
    ident = loaded_report[1]['reports'][0]['ident']
    ipr_conn = IprConnection(
        username=os.environ.get('IPR_USER', os.environ['USER']),
        password=os.environ['IPR_PASS'],
        url=os.environ['IPR_TEST_URL'],
    )
    return ipr_conn.get(uri=f"reports/{ident}/{section_name}")


@pytest.mark.skipif(
    not INCLUDE_UPLOAD_TESTS, reason="excluding tests of upload to live ipr instance"
)
@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
class TestCreateReport:
    def test_patient_id_loaded_once(self, loaded_report: Tuple) -> None:
        patient_id = loaded_report[0]
        assert loaded_report[1]['total'] == 1
        assert loaded_report[1]['reports'][0]['patientId'] == patient_id

    # TODO; add main section checks
    def test_main_sections_present(self, loaded_report: Tuple) -> None:
        return

    def test_expression_variants_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, 'expression-variants')
        kbmatched = [item for item in section if item['kbMatches']]
        assert 'PTP4A3' in [item['gene']['name'] for item in kbmatched]

    def test_structural_variants_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, 'structural-variants')
        kbmatched = [item for item in section if item['kbMatches']]
        assert '(EWSR1,FLI1):fusion(e.7,e.4)' in [item['displayName'] for item in kbmatched]

    def test_small_mutations_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, 'small-mutations')
        kbmatched = [item for item in section if item['kbMatches']]
        assert 'FGFR2:p.R421C' in [item['displayName'] for item in kbmatched]
        assert 'CDKN2A:p.T18M' in [item['displayName'] for item in kbmatched]

    def test_copy_variants_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, 'copy-variants')
        kbmatched = [item for item in section if item['kbMatches']]
        assert ('ERBB2', 'amplification') in [
            (item['gene']['name'], item['displayName']) for item in kbmatched
        ]

    def test_kb_matches_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, 'kb-matches')
        observed_and_matched = set(
            [(item['kbVariant'], item['variant']['displayName']) for item in section]
        )
        for pair in [
            ('ERBB2 amplification', 'amplification'),
            ('FGFR2 mutation', 'FGFR2:p.R421C'),
            ('PTP4A3 overexpression', 'increased expression'),
            ('EWSR1 and FLI1 fusion', '(EWSR1,FLI1):fusion(e.7,e.4)'),
            ('CDKN2A mutation', 'CDKN2A:p.T18M'),
        ]:
            assert pair in observed_and_matched

    def test_therapeutic_targets_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, 'therapeutic-targets')
        therapeutic_target_genes = set([item['gene'] for item in section])
        for gene in ['CDKN2A', 'ERBB2', 'FGFR2', 'PTP4A3']:
            assert gene in therapeutic_target_genes

    def test_genomic_alterations_identified_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, 'summary/genomic-alterations-identified')
        variants = set([item['geneVariant'] for item in section])
        for variant in [
            'FGFR2:p.R421C',
            'PTP4A3 (high_percentile)',
            'ERBB2 (Amplification)',
            '(EWSR1,FLI1):fusion(e.7,e.4)',
            'CDKN2A:p.T18M',
        ]:
            assert variant in variants

    def test_analyst_comments_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, 'summary/analyst-comments')
        assert section['comments']

    @pytest.mark.skip('todo pending ipr-api report-sample update')
    def test_sample_info_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, 'reports-sample-info')

    @pytest.mark.skipif(EXCLUDE_BCGSC_TESTS, reason="excluding tests requiring BCGSC loaders")
    def test_pharmacogenomic_variants_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, 'kb-matches?category=pharmacogenomic')
        assert section

    @pytest.mark.skipif(EXCLUDE_BCGSC_TESTS, reason="excluding tests requiring BCGSC loaders")
    def test_cancer_predisposition_variants_loaded(self, loaded_report: Tuple) -> None:
        section = get_section(loaded_report, 'kb-matches?category=cancer%20predisposition')
        assert section
