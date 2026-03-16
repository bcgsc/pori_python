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

EXCLUDE_BCGSC_TESTS = os.environ.get('EXCLUDE_BCGSC_TESTS') == '1'
EXCLUDE_ONCOKB_TESTS = os.environ.get('EXCLUDE_ONCOKB_TESTS') == '1'
INCLUDE_UPLOAD_TESTS = os.environ.get('INCLUDE_UPLOAD_TESTS', '0') == '1'
DELETE_UPLOAD_TEST_REPORTS = os.environ.get('DELETE_UPLOAD_TEST_REPORTS', '1') == '1'


def get_test_spec():
    ipr_spec = {'components': {'schemas': {'genesCreate': {'properties': {}}}}}
    ipr_gene_keys = IprGene.__required_keys__ | IprGene.__optional_keys__
    for key in ipr_gene_keys:
        ipr_spec['components']['schemas']['genesCreate']['properties'][key] = ''
    return ipr_spec


def get_test_file(name: str) -> str:
    return os.path.join(os.path.dirname(__file__), 'test_data', name)


@pytest.fixture(scope='module')
def loaded_reports(tmp_path_factory) -> Generator:
    """
    Load test data with selective flagging to enable end-to-end testing of:
    1. Flags from variant input TSVs (expressionVariants, smallMutations, etc.)
    2. Flags from transcript annotation TSV (transcript_flags)
    3. Variants without flags should NOT have observedVariantAnnotation entries

    This fixture:
    - Only flags a subset of variants per type (e.g. first expression variant, first mutation)
    - Creates a transcript flags TSV with 2 flags for APC gene variants
    - Verifies compatibility with pori_ipr_api's observedVariantAnnotation model
    """
    json_file = tmp_path_factory.mktemp('inputs') / 'content.json'
    async_json_file = tmp_path_factory.mktemp('inputs') / 'async_content.json'
    transcript_flags_file = tmp_path_factory.mktemp('inputs') / 'transcript_flags.tsv'
    patient_id = f'TEST_WITH_FLAGS{str(uuid.uuid4())}'
    async_patient_id = f'TEST_WITH_FLAGS_ASYNC_{str(uuid.uuid4())}'

    # Load data - only flag SOME variants to test that unflagged ones don't get annotations
    expvars = pd.read_csv(get_test_file('expression.short.tab'), sep='\t')
    # Flag only the first expression variant
    expvars['flags'] = ''
    expvars_variant_locs = expvars[~pd.isnull(expvars.kbCategory)].index[0:2].tolist()
    expvars.loc[expvars_variant_locs[0], 'flags'] = 'expression_flag_1'
    expvars.loc[expvars_variant_locs[1], 'flags'] = 'expression_flag_1,expression_flag_2'  # test multiple flags in one string
    expvars_json = expvars.to_json(orient='records')

    smallmuts = pd.read_csv(get_test_file('small_mutations.short.tab'), sep='\t')
    # Flag only the first small mutation
    smallmuts['flags'] = ''
    smallmuts.loc[0, 'flags'] = 'mutation_flag_1'

    # Find the first small mutation that is not on APC gene to avoid overlap with transcript flags test
    non_apc_indices = smallmuts[smallmuts['gene'] != 'APC'].index
    multi_flag_index = non_apc_indices[0]
    smallmuts.loc[multi_flag_index, 'flags'] = 'mutation_flag_2,mutation_flag_1'  # test multiple flags in one string

    # get transcript for this mutation to match in transcript flags file
    smallmut_gene = smallmuts.loc[multi_flag_index, 'gene']
    smallmut_transcript = smallmuts.loc[multi_flag_index, 'transcript']
    smallmuts_json = smallmuts.to_json(orient='records')

    copyvars = pd.read_csv(get_test_file('copy_variants.short.tab'), sep='\t')
    # Flag only the first copy variant
    copyvars['flags'] = ''
    copyvars.loc[0, 'flags'] = 'cnv_flag_1'
    copyvars.loc[1, 'flags'] = 'cnv_flag_1,cnv_flag_2'  # test multiple flags in one string
    copyvars_json = copyvars.to_json(orient='records')

    svs = pd.read_csv(get_test_file('fusions.tab'), sep='\t')
    # Flag only the first SV
    svs['flags'] = ''
    svs.loc[0, 'flags'] = 'sv_flag_1'
    svs.loc[1, 'flags'] = 'sv_flag_1,sv_flag_2'  # test multiple flags in one string
    svs_json = svs.to_json(orient='records')

    hla = pd.read_csv(get_test_file('hla_variants.tab'), sep='\t')
    hla_json = hla.to_json(orient='records')

    # Create a transcript flags file with flags for specific transcripts
    # Match transcripts from small_mutations.short.tab
    transcript_flags_df = pd.DataFrame({
        'gene': ['APC', 'APC', smallmut_gene, 'svgene1', 'svgene2', 'svgene3', 'svgene4'],
        'transcript': ['ENST00000457016', 'ENST00000257430', smallmut_transcript,'ENST00000358273', 'ENST00000397938', 'ENST00000373930', 'ENST00000457710'],
        'flags': ['transcript_flag_1', 'transcript_flag_2', 'additional_transcript_flag', 'sv_transcript_flag_1', 'sv_transcript_flag_2', 'sv_transcript_flag_3', 'sv_transcript_flag_4'],
    })
    transcript_flags_df.to_csv(transcript_flags_file, sep='\t', index=False)

    json_contents = {
        'comparators': [
            {'analysisRole': 'expression (disease)', 'name': '1'},
            {'analysisRole': 'expression (primary site)', 'name': '2'},
            {'analysisRole': 'expression (biopsy site)', 'name': '3'},
            {
                'analysisRole': 'expression (internal pancancer cohort)',
                'name': '4',
            },
        ],
        'patientId': patient_id,
        'project': 'TEST',
        'sampleInfo': [
            {
                'sample': 'Constitutional',
                'biopsySite': 'Normal tissue',
                'sampleName': 'SAMPLE1-PB',
                'primarySite': 'Blood-Peripheral',
                'collectionDate': '11-11-11',
            },
            {
                'sample': 'Tumour',
                'pathoTc': '90%',
                'biopsySite': 'hepatic',
                'sampleName': 'SAMPLE2-FF-1',
                'primarySite': 'Vena Cava-Hepatic',
                'collectionDate': '12-12-12',
            },
        ],
        'kbDiseaseMatch': 'colorectal cancer',
        'msi': [
            {
                'score': 1000.0,
                'kbCategory': 'microsatellite instability',
            }
        ],
        'hrd': {
            'score': 9999.0,
            'kbCategory': 'homologous recombination deficiency strong signature',
        },
        'expressionVariants': json.loads(expvars_json),
        'smallMutations': json.loads(smallmuts_json),
        'copyVariants': json.loads(copyvars_json),
        'structuralVariants': json.loads(svs_json),
        'cosmicSignatures': pd.read_csv(
            get_test_file('cosmic_variants.tab'), sep='\t'
        ).signature.tolist(),
        'hlaTypes': json.loads(hla_json),
    }

    json_contents['patientId'] = async_patient_id
    async_json_file.write_text(
        json.dumps(
            json_contents,
            allow_nan=False,
        )
    )

    argslist = [
        'ipr',
        '--username',
        os.environ.get('IPR_USER', os.environ['USER']),
        '--password',
        os.environ['IPR_PASS'],
        '--graphkb_username',
        os.environ.get('GRAPHKB_USER', os.environ.get('IPR_USER', os.environ['USER'])),
        '--graphkb_password',
        os.environ.get('GRAPHKB_PASS', os.environ['IPR_PASS']),
        '--ipr_url',
        os.environ['IPR_TEST_URL'],
        '--graphkb_url',
        os.environ.get('GRAPHKB_URL', False),
        '--therapeutics',
        '--allow_partial_matches',
        '-o upload_with_flags.json',
        '--transcript_flags',
        str(transcript_flags_file),
    ]

    async_argslist = argslist.copy()
    async_argslist.extend(['--content', str(async_json_file), '--async_upload'])
    with patch.object(sys, 'argv', async_argslist):
        with patch.object(IprConnection, 'get_spec', return_value=get_test_spec()):
            command_interface()

    ipr_conn = IprConnection(
        username=os.environ.get('IPR_USER', os.environ['USER']),
        password=os.environ['IPR_PASS'],
        url=os.environ['IPR_TEST_URL'],
    )
    async_loaded_report = ipr_conn.get(uri=f'reports?searchText={async_patient_id}')

    # Collect expected flagged genes for each variant type
    expected_flagged = {
        'expression': expvars[expvars['flags'] != '']['gene'].tolist(),
        'small_mutations': smallmuts[smallmuts['flags'] != '']['gene'].tolist(),
        'copy_variants': copyvars[copyvars['flags'] != '']['gene'].tolist(),
        'structural_variants_cterm': svs[svs['flags'] != '']['gene2'].tolist(),
        'structural_variants_nterm': svs[svs['flags'] != '']['gene1'].tolist()
    }

    loaded_reports_result = {
        'async': (async_patient_id, async_loaded_report),
        'expected_flagged': expected_flagged,
    }
    yield loaded_reports_result
    if DELETE_UPLOAD_TEST_REPORTS:
        ipr_conn.delete(uri=f'reports/{async_loaded_report["reports"][0]["ident"]}')


def get_section(loaded_report, section_name):
    ident = loaded_report[1]['reports'][0]['ident']
    ipr_conn = IprConnection(
        username=os.environ.get('IPR_USER', os.environ['USER']),
        password=os.environ['IPR_PASS'],
        url=os.environ['IPR_TEST_URL'],
    )
    return ipr_conn.get(uri=f'reports/{ident}/{section_name}')


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
        for key in ('ident', 'updatedAt', 'createdAt', 'deletedAt', 'reportId', 'variantId', 'id'):
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
    not INCLUDE_UPLOAD_TESTS, reason='excluding tests of upload to live ipr instance'
)
@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason='excluding long running integration tests')
class TestCreateReport:
    def test_patient_id_loaded_once(self, loaded_reports) -> None:
        async_patient_id = loaded_reports['async'][0]
        assert loaded_reports['async'][1]['total'] == 1
        assert loaded_reports['async'][1]['reports'][0]['patientId'] == async_patient_id

    def test_observed_variant_annotations_loaded(self, loaded_reports) -> None:
        """Test that flagged variants have observedVariantAnnotation with correct flags."""
        variants_section = get_section(loaded_reports['async'], 'variants')
        expected_flagged = loaded_reports['expected_flagged']

        # Check that expression variant with input flag has annotation
        exp_vars_with_annot = [v for v in variants_section if v['variantType'] == 'exp' and 'observedVariantAnnotation' in v and v['observedVariantAnnotation'] is not None]
        assert len(exp_vars_with_annot) > 0, "Should have at least one expression variant with annotation"

        # Find the flagged expression variants
        flagged_exp_genes = expected_flagged['expression']
        for gene in flagged_exp_genes:
            flagged_exp = [v for v in exp_vars_with_annot if v['gene']['name'] == gene]
            assert len(flagged_exp) >= 1, f"{gene} should be flagged with input flag"
            for var in flagged_exp:
                assert any(['expression_flag' in str(var['observedVariantAnnotation'].get('flags', []))])

        # Check that the variant with multiple flags has both flags correctly split
        multi_flag_exp = [
            v for v in exp_vars_with_annot if len(v['observedVariantAnnotation']['flags'])>1
        ]
        assert len(multi_flag_exp) > 0, "Should have at least one expression variant with multiple flags"
        for var in multi_flag_exp:
            flags = var['observedVariantAnnotation'].get('flags', [])
            if len(flags) > 1:
                assert 'expression_flag_1' in flags and 'expression_flag_2' in flags, \
                    f"Variant with multiple flags should have both expression_flag_1 and expression_flag_2, got {flags}"

        # Check that unflagged expression variants don't have observedVariantAnnotation
        unflagged_exp = [v for v in variants_section if v['variantType'] == 'exp' and v['gene']['name'] not in flagged_exp_genes]
        for var in unflagged_exp:
            assert 'observedVariantAnnotation' not in var or var['observedVariantAnnotation'] is None, \
                f"Unflagged expression variant {var['gene']['name']} should not have annotations"

    def test_variant_transcript_annotations_loaded(self, loaded_reports) -> None:
        """Test that variants with transcript flags have observedVariantAnnotation from transcript file."""
        variants_section = get_section(loaded_reports['async'], 'variants')

        # Find small mutation variants
        mut_vars = [v for v in variants_section if v['variantType'] == 'mut']
        assert len(mut_vars) > 0, "Should have small mutations"

        # each APC with ENST00000457016 transcript should have transcript_flag_1
        apc_mut_enst1 = [v for v in mut_vars if v['gene']['name'] == 'APC' and v.get('transcript') == 'ENST00000457016']
        assert len(apc_mut_enst1) > 0, "Should find at least one APC mutation with ENST00000457016 transcript"

        has_transcript_flag = True
        for var in apc_mut_enst1:
            if 'observedVariantAnnotation' in var and var['observedVariantAnnotation'] is not None:
                if 'transcript_flag_1' not in var['observedVariantAnnotation'].get('flags', []):
                    has_transcript_flag = False
                    break

        assert has_transcript_flag, "All mutations with transcript ENST00000457016 should have transcript_flag_1 from transcript file"

        # APC with ENST00000257430 transcript should have transcript_flag_2
        apc_mut_enst2 = [v for v in mut_vars if v['gene']['name'] == 'APC' and v.get('transcript') == 'ENST00000257430']
        assert len(apc_mut_enst2) > 0, "Should find APC mutation with ENST00000257430 transcript"

        has_second_flag = True
        for var in apc_mut_enst2:
            if 'observedVariantAnnotation' in var and var['observedVariantAnnotation'] is not None:
                if 'transcript_flag_2' not in var['observedVariantAnnotation'].get('flags', []):
                    has_second_flag = False
                    break

        assert has_second_flag, "All APC mutations with transcript ENST00000257430 should have transcript_flag_2 from transcript file"

        # Check that the variant with a transcript flag and multiple input flags
        # has all flags correctly represented in observedVariantAnnotation
        annotated_mut_vars = [v for v in mut_vars if 'observedVariantAnnotation' in v and v['observedVariantAnnotation'] is not None]
        multi_flag_mut = [v for v in annotated_mut_vars if len(v.get('observedVariantAnnotation').get('flags', [])) > 2]
        assert len(multi_flag_mut) > 0, "Should have at least one small mutation variant with multiple flags"
        for var in multi_flag_mut:
            flags = var['observedVariantAnnotation'].get('flags', [])
            if len(flags) > 2:
                assert 'mutation_flag_1' in flags and 'mutation_flag_2' in flags, \
                    f"Variant with multiple flags should have both mutation_flag_1 and mutation_flag_2, got {flags}"

    def test_fusion_variants_have_multiple_transcript_annotations_loaded(self, loaded_reports) -> None:
        """Test that variants with transcript flags have observedVariantAnnotation from transcript file."""
        variants_section = get_section(loaded_reports['async'], 'variants')

        # Find small mutation variants
        svs = [v for v in variants_section if v['variantType'] == 'sv']
        annotated_svs = [v for v in svs if 'observedVariantAnnotation' in v and v['observedVariantAnnotation'] is not None]

        assert len(annotated_svs) > 0, "Should have annotated svs"
