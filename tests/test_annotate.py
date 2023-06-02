import os
import pytest
from graphkb import GraphKBConnection

from ipr.annotate import (
    annotate_positional_variants,
    get_gene_information,
    get_therapeutic_associated_genes,
)
from ipr.types import IprSmallMutationVariant

# TP53 examples from https://www.bcgsc.ca/jira/browse/SDEV-3122
# Mutations are actually identical but on alternate transcripts.

TP53_MUT_DICT = {
    'pref': IprSmallMutationVariant(  # type: ignore
        {
            'key': 'SDEV-3122_preferred',
            'gene': 'TP53',
            'hgvsGenomic': 'chr17:g.7674252C>T',
            'hgvsCds': 'ENST00000269305:c.711G>A',
            'hgvsProtein': 'TP53:p.M237I',
        }
    ),
    'intersect': IprSmallMutationVariant(  # type: ignore
        {
            'key': 'SDEV-3122_alt',
            'gene': 'TP53',
            'hgvsGenomic': 'chr17:g.7674252C>T',
            'hgvsCds': 'ENST00000610292:c.594G>A',
            'hgvsProtein': 'TP53:p.M198I',
        }
    ),
    'prot_only': IprSmallMutationVariant(  # type: ignore
        {
            'key': 'prot_only',
            'gene': 'TP53',
            'hgvsProtein': 'TP53:p.M237I',
        }
    ),
    'cds_only': IprSmallMutationVariant(  # type: ignore
        {
            'key': 'cds_only',
            'gene': 'TP53',
            'hgvsCds': 'ENST00000269305:c.711G>A',
        }
    ),
    'genome_only': IprSmallMutationVariant(  # type: ignore
        {
            'key': 'genome_only',
            'gene': 'TP53',
            'hgvsGenomic': 'chr17:g.7674252C>T',
        }
    ),
}


@pytest.fixture(scope='module')
def graphkb_conn():
    username = os.environ['IPR_USER']
    password = os.environ['IPR_PASS']
    graphkb_conn = GraphKBConnection()
    graphkb_conn.login(username, password)
    return graphkb_conn


def test_annotate_nonsense_vs_missense(graphkb_conn):
    """Verify missense (point mutation) is not mistaken for a nonsense (stop codon) mutation."""
    disease = 'cancer'
    for key in ('prot_only', 'cds_only', 'genome_only', 'pref'):
        matched = annotate_positional_variants(graphkb_conn, [TP53_MUT_DICT[key]], disease)
        # nonsense - stop codon - should not match.  This is missense not nonsense (#164:933).
        nonsense = [a for a in matched if a['kbVariant'] == 'TP53 nonsense']
        assert not nonsense, f"nonsense matched to {key}: {TP53_MUT_DICT[key]}"
        assert matched, f"should have matched in {key}: {TP53_MUT_DICT[key]}"


def test_annotate_nonsense_vs_missense_protein(graphkb_conn):
    """Verify missense (point mutation) is not mistaken for a nonsense (stop codon) mutation."""
    disease = 'cancer'
    for key in ('prot_only', 'pref'):
        matched = annotate_positional_variants(graphkb_conn, [TP53_MUT_DICT[key]], disease)
        # nonsense - stop codon - should not match.  This is missense not nonsense (#164:933).
        nonsense = [a for a in matched if 'nonsense' in a['kbVariant']]
        assert not nonsense, f"nonsense matched to {key}: {TP53_MUT_DICT[key]}"
        assert matched, f"should have matched in {key}: {TP53_MUT_DICT[key]}"


def test_annotate_structural_variants_tp53(graphkb_conn):
    """Verify alternate TP53 variants match."""
    disease = 'cancer'
    ref_key = 'prot_only'
    pref = annotate_positional_variants(graphkb_conn, [TP53_MUT_DICT[ref_key]], disease)
    # GERO-299 - nonsense - stop codon - should not match.  This is missense not nonsense (#164:933).
    nonsense = [a for a in pref if a['kbVariant'] == 'TP53 nonsense']
    assert not nonsense
    pref_vars = set([m['kbVariant'] for m in pref])
    assert pref_vars, f"No matches to {TP53_MUT_DICT[pref]}"
    print(pref_vars)
    for key, alt_rep in TP53_MUT_DICT.items():
        if key == ref_key:
            continue
        alt = annotate_positional_variants(graphkb_conn, [alt_rep], disease)
        alt_vars = set([m['kbVariant'] for m in alt])
        diff = pref_vars.symmetric_difference(alt_vars)
        missing = pref_vars.difference(alt_vars)

        known_issues = set()
        if key == 'genome_only':
            # genome_only matched to more precise type 'TP53 deleterious mutation' but not 'TP53 mutation'
            known_issues.add('TP53 mutation')

        missing = pref_vars.difference(alt_vars).difference(known_issues)
        print(alt_vars)
        assert not missing, f"{key} missing{missing}: {diff}"


def test_get_therapeutic_associated_genes(graphkb_conn):
    gene_list = get_therapeutic_associated_genes(graphkb_conn=graphkb_conn)
    assert gene_list, 'No get_therapeutic_associated_genes found'
    assert (
        len(gene_list) > 500
    ), f'Expected over 500 get_therapeutic_associated_genes but found {len(gene_list)}'


def test_get_gene_information(graphkb_conn):
    GENES = {
        'not_a_gene': set([]),
        'AKT3': set(['cancerRelated', 'knownSmallMutation', 'oncogene', 'therapeuticAssociated']),
        'TP53': set(
            [
                'knownFusionPartner',
                'cancerRelated',
                'therapeuticAssociated',
                'knownSmallMutation',
                'tumourSuppressor',
            ]
        ),
        'BRAF': set(
            [
                'knownFusionPartner',
                'cancerRelated',
                'therapeuticAssociated',
                'knownSmallMutation',
                'oncogene',
            ]
        ),
        'ERBB2': set(
            [
                'knownFusionPartner',
                'cancerRelated',
                'therapeuticAssociated',
                'knownSmallMutation',
                'oncogene',
            ]
        ),
        'ABCB1': set(['cancerRelated', 'therapeuticAssociated', 'knownSmallMutation']),
    }
    info = get_gene_information(graphkb_conn, GENES.keys())

    assert info, f"get_gene_information failed for {GENES.keys()}"
    fails = []

    for gene, flags in GENES.items():
        gene_info = [i for i in info if i['name'] == gene]
        if not flags:
            assert not gene_info, f"{gene} assumed to have all False (blank) got {gene_info}"
        else:
            gene_info_flags = set([key for key in gene_info[0].keys() if key != 'name'])
            if flags != gene_info_flags:
                fails.append(f"Expected {gene}: {sorted(flags)} - got {sorted(gene_info_flags)}")
    assert not fails, '\n'.join(fails)
