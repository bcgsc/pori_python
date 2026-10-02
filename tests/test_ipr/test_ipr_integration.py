import os
import json
import jsonschema
import pytest

from pori_python.ipr.connection import IprConnection

from .constants import EXCLUDE_INTEGRATION_TESTS


mock_ipr_results = [
    [
        {
            'text': '<p>no cancerType</p>',
            'variantName': 'ERBB2 amplification',
            'cancerType': [],
            'template': {'name': 'test3'},
            'projects': [{'name': 'test2'}],
        },
        {
            'text': '<p>normal</p>',
            'variantName': 'ERBB2 amplification',
            'cancerType': ['test1', 'test'],
            'template': {'name': 'test3'},
            'projects': [{'name': 'test2'}],
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
            'projects': [{'name': 'test2'}],
        },
    ],
    [
        {
            'text': '<p>normal, second variant</p>',
            'variantName': 'second variant',
            'cancerType': ['test1', 'test'],
            'template': {'name': 'test3'},
            'projects': [{'name': 'test2'}],
        },
    ],
]


def validate_mock_ipr_results_against_schema(schema, results):
    # Since the IPR generated schema has required fields, but the purpose of the test is to prevent
    # schema drift and not to validate POST payloads to variant-text, we create a custom validator
    # that ignores the required fields and only validates the properties that are defined in the schema.
    custom_validators = dict(jsonschema.Draft7Validator.VALIDATORS)
    custom_validators.pop('required', None)
    non_required_validator = jsonschema.validators.create(
        meta_schema=jsonschema.Draft7Validator.META_SCHEMA,
        validators=custom_validators,
        version='draft7-lazy',
    )
    validator = non_required_validator(schema)
    declared_properties = set(schema.get('properties', {}))

    for group_index, group in enumerate(results):
        for record_index, record in enumerate(group):
            location = f'group {group_index} record {record_index}'
            unknown = sorted(set(record) - declared_properties)
            assert not unknown, (
                f'{location}: properties {unknown} are not defined in the IPR variant-text schema'
            )
            errors = sorted(validator.iter_errors(record), key=lambda err: list(err.path))
            assert not errors, f'{location}: ' + '; '.join(
                f'{list(err.path)}: {err.message}' for err in errors
            )


def make_ipr_connection():
    return IprConnection(
        username=os.environ.get('IPR_USER', os.environ['USER']),
        password=os.environ['IPR_PASS'],
        url=os.environ['IPR_INTEGRATION_TEST_URL'],
    )


# DEVSU-3011 adding IPR integration test to validate variant text schema from the IPR API side to prevent schema drift
@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason='excluding long running integration tests')
class TestVariantTextSchema:
    def test_mock_ipr_results_match_variant_text_schema(self):
        ipr_conn = make_ipr_connection()
        schema = ipr_conn.request(
            'variant-text/schema',
            method='GET',
            headers=json.dumps(
                {
                    'Content-Length': '0',
                    'Accept': 'application/json',
                }
            ),
        )
        validate_mock_ipr_results_against_schema(schema, mock_ipr_results)
