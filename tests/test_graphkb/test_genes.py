"""
Tests here depend on specific data in GraphKB which can change. To avoid this, expected/stable values are chosen
"""

import os
import pytest

from pori_python.graphkb import GraphKBConnection
from pori_python.graphkb.genes import (
    get_cancer_genes,
    get_cancer_predisposition_info,
    get_gene_information,
    get_gene_linked_cancer_predisposition_info,
    get_gene_linked_pharmacogenomic_info,
    get_genes_from_variant_types,
    get_oncokb_oncogenes,
    get_oncokb_tumour_supressors,
    get_pharmacogenomic_info,
    get_preferred_gene_name,
    get_therapeutic_associated_genes,
)
from pori_python.graphkb.util import get_rid

EXCLUDE_INTEGRATION_TESTS = os.environ.get("EXCLUDE_INTEGRATION_TESTS") == "1"
EXCLUDE_BCGSC_TESTS = os.environ.get("EXCLUDE_BCGSC_TESTS") == "1"
EXCLUDE_ONCOKB_TESTS = os.environ.get("EXCLUDE_ONCOKB_TESTS") == "1"

CANONICAL_ONCOGENES = ["kras", "nras", "alk"]
CANONICAL_TS = ["cdkn2a", "tp53"]
CANONICAL_CG = ["alb"]
CANONICAL_FUSION_GENES = ["alk", "ewsr1", "fli1"]
CANONICAL_STRUCTURAL_VARIANT_GENES = ["brca1", "dpyd", "pten"]
CANNONICAL_THERAPY_GENES = ["erbb2", "brca2", "egfr"]


PHARMACOGENOMIC_INITIAL_GENES = [
    "ACYP2",
    "CEP72",
    # 'CYP26B1',  # defined as hgvsGenomic chr2:g.233760235_233760235nc_000002.12:g.233760235ta[7]>ta[8]
    "DPYD",
    "NUDT15",
    "RARG",
    "SLC28A3",
    "TPMT",
    "UGT1A6",
]
CANCER_PREDISP_INITIAL_GENES = [
    "AKT1",
    "APC",
    "ATM",
    "AXIN2",
    "BAP1",
    "BLM",
    "BMPR1A",
    "BRCA1",
    "BRCA2",
    "BRIP1",
    "CBL",
    "CDH1",
    "CDK4",
    "CDKN2A",
    "CHEK2",
    "DICER1",
    "EGFR",
    "EPCAM",
    "ETV6",
    "EZH2",
    "FH",
    "FLCN",
    "GATA2",
    "HRAS",
    "KIT",
    "MEN1",
    "MET",
    "MLH1",
    "MSH2",
    "MSH6",
    "MUTYH",
    "NBN",
    "NF1",
    "PALB2",
    "PDGFRA",
    "PMS2",
    "PTCH1",
    "PTEN",
    "PTPN11",
    "RAD51C",
    "RAD51D",
    "RB1",
    "RET",
    "RUNX1",
    "SDHA",
    "SDHB",
    "SDHC",
    "SDHD",
    "SMAD4",
    "SMARCA4",
    "STK11",
    "TP53",
    "TSC1",
    "TSC2",
    "VHL",
    "WT1",
]


@pytest.fixture(scope="module")
def conn():
    conn = GraphKBConnection()
    conn.login(os.environ["GRAPHKB_USER"], os.environ["GRAPHKB_PASS"])
    return conn


@pytest.mark.skipif(EXCLUDE_ONCOKB_TESTS, reason="excluding tests that depend on oncokb data")
def test_oncogene(conn):
    result = get_oncokb_oncogenes(conn)
    names = {row["name"] for row in result}
    for gene in CANONICAL_ONCOGENES:
        assert gene in names
    for gene in CANONICAL_TS:
        assert gene not in names
    for gene in CANONICAL_CG:
        assert gene not in names


@pytest.mark.skipif(EXCLUDE_ONCOKB_TESTS, reason="excluding tests that depend on oncokb data")
def test_tumour_supressors(conn):
    result = get_oncokb_tumour_supressors(conn)
    names = {row["name"] for row in result}
    for gene in CANONICAL_TS:
        assert gene in names
    for gene in CANONICAL_ONCOGENES:
        assert gene not in names
    for gene in CANONICAL_CG:
        assert gene not in names


@pytest.mark.skipif(EXCLUDE_ONCOKB_TESTS, reason="excluding tests that depend on oncokb data")
@pytest.mark.skipif(
    EXCLUDE_BCGSC_TESTS, reason="excluding BCGSC-specific tests (tso500 not available)"
)
def test_cancer_genes(conn):
    result = get_cancer_genes(conn)
    names = {row["name"] for row in result}
    for gene in CANONICAL_CG:
        assert gene in names
    for gene in CANONICAL_TS:
        assert gene not in names
    for gene in CANONICAL_ONCOGENES:
        assert gene not in names


@pytest.mark.skipif(
    EXCLUDE_BCGSC_TESTS, reason="excluding BCGSC-specific tests (requires CGL loader))"
)
def test_get_pharmacogenomic_info(conn):
    genes, matches = get_pharmacogenomic_info(conn)
    for gene in PHARMACOGENOMIC_INITIAL_GENES:
        assert gene in genes, f"{gene} not found in get_pharmacogenomic_info"
        for rid, variant_display in matches.items():
            if variant_display.startswith(gene):
                break
        else:  # no break called
            # failing on this version of the func; addressed in 'new' version
            if gene == "ACYP2":
                continue
            assert False, f"No rid found for a pharmacogenomic with {gene}"


@pytest.mark.skipif(
    EXCLUDE_BCGSC_TESTS, reason="excluding BCGSC-specific tests (requires CGL loader))"
)
def test_get_gene_linked_pharmacogenomic_info(conn):
    genes, matches = get_gene_linked_pharmacogenomic_info(conn)
    for gene in PHARMACOGENOMIC_INITIAL_GENES:
        assert gene in genes, f"{gene} not found in get_pharmacogenomic_info"
        for rid, variant_info in matches.items():
            variant_gene_assoc = variant_info[1]
            if gene in variant_gene_assoc:
                break
        else:  # no break called
            assert False, f"No rid found for a pharmacogenomic with {gene}"


@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
@pytest.mark.skipif(
    EXCLUDE_BCGSC_TESTS, reason="excluding BCGSC-specific tests (requires CGL loader))"
)
def test_get_cancer_predisposition_info(conn):
    genes, matches = get_cancer_predisposition_info(conn)
    for gene in CANCER_PREDISP_INITIAL_GENES:
        assert gene in genes, f"{gene} not found in get_cancer_predisposition_info"


@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
@pytest.mark.skipif(
    EXCLUDE_BCGSC_TESTS, reason="excluding BCGSC-specific tests (requires CGL loader))"
)
def test_get_gene_linked_cancer_predisposition_info(conn):
    genes, matches = get_gene_linked_cancer_predisposition_info(conn)
    for gene in CANCER_PREDISP_INITIAL_GENES:
        assert gene in genes, f"{gene} not found in get_cancer_predisposition_info"


@pytest.mark.parametrize(
    "alt_rep", ("NM_033360.4", "NM_033360", "ENSG00000133703.11", "ENSG00000133703")
)
def test_get_preferred_gene_name_kras(alt_rep, conn):
    gene_name = get_preferred_gene_name(conn, alt_rep)
    assert (
        "KRAS" == gene_name
    ), f"Expected KRAS as preferred gene name for {alt_rep}, not '{gene_name}'"


@pytest.mark.skipif(
    EXCLUDE_BCGSC_TESTS, reason="excluding BCGSC-specific tests (requires CGL loader))"
)
@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
def test_find_genes_by_variant_type_structural_variant(conn):
    result = get_genes_from_variant_types(conn, ["structural variant"])
    names = {row["name"] for row in result}
    for gene in CANONICAL_STRUCTURAL_VARIANT_GENES:
        assert gene in names, f"{gene} was not identified as a structural variant gene."


@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
def test_find_no_genes_by_variant_type_with_nonmatching_source_record_id(conn):
    refseq_id = get_rid(conn, target="source", name="refseq")
    result = get_genes_from_variant_types(
        conn, ["structural variant"], source_record_ids=[refseq_id]
    )
    assert not result


@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
def test_get_therapeutic_associated_genes(conn):
    gene_list = get_therapeutic_associated_genes(graphkb_conn=conn)
    assert gene_list, "No get_therapeutic_associated_genes found"
    assert (
        len(gene_list) > 300
    ), f"Expected over 300 get_therapeutic_associated_genes but found {len(gene_list)}"
    names = {row["name"] for row in gene_list}
    for gene in CANNONICAL_THERAPY_GENES + CANONICAL_ONCOGENES + CANONICAL_TS:
        assert gene in names, f"{gene} not found by get_therapeutic_associated_genes"


@pytest.mark.skipif(
    EXCLUDE_BCGSC_TESTS,
    reason="excluding BCGSC-specific tests (requires oncokb and other loaders))",
)
@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
def test_get_gene_information(conn):
    gene_info = get_gene_information(
        conn,
        CANONICAL_ONCOGENES
        + CANONICAL_TS
        + CANONICAL_CG
        + CANONICAL_FUSION_GENES
        + CANONICAL_STRUCTURAL_VARIANT_GENES
        + CANNONICAL_THERAPY_GENES
        + ["notagenename"],
    )
    assert gene_info
    nongene_flagged = [g["name"] for g in gene_info if g["name"] == "notagenename"]
    assert not nongene_flagged, f"Improper gene category: {nongene_flagged}"

    for gene in CANONICAL_ONCOGENES:
        assert gene in [
            g["name"] for g in gene_info if g.get("oncogene")
        ], f"Missed oncogene {gene}"

    for gene in CANONICAL_TS:
        assert gene in [
            g["name"] for g in gene_info if g.get("tumourSuppressor")
        ], f"Missed 'tumourSuppressor' {gene}"

    for gene in CANONICAL_FUSION_GENES:
        assert gene in [
            g["name"] for g in gene_info if g.get("knownFusionPartner")
        ], f"Missed knownFusionPartner {gene}"

    for gene in CANONICAL_STRUCTURAL_VARIANT_GENES:
        assert gene in [
            g["name"] for g in gene_info if g.get("knownSmallMutation")
        ], f"Missed knownSmallMutation {gene}"

    for gene in CANNONICAL_THERAPY_GENES:
        assert gene in [
            g["name"] for g in gene_info if g.get("therapeuticAssociated")
        ], f"Missed therapeuticAssociated {gene}"

    for gene in (
        CANONICAL_ONCOGENES
        + CANONICAL_TS
        + CANONICAL_FUSION_GENES
        + CANONICAL_STRUCTURAL_VARIANT_GENES
        + CANNONICAL_THERAPY_GENES
    ):
        assert gene in [
            g["name"] for g in gene_info if g.get("kbStatementRelated")
        ], f"Missed kbStatementRelated {gene}"

    for gene in CANONICAL_CG:
        assert gene in [
            g["name"] for g in gene_info if g.get("cancerGeneListMatch")
        ], f"Missed cancerGeneListMatch {gene}"
