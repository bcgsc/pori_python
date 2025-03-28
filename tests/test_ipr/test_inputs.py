import json
import numpy as np
import os
import pandas as pd
import pytest
from unittest import mock

from pori_python.graphkb.match import INPUT_COPY_CATEGORIES
from pori_python.ipr.inputs import (
    COPY_OPTIONAL,
    check_comparators,
    check_variant_links,
    create_graphkb_sv_notation,
    preprocess_copy_variants,
    preprocess_cosmic,
    preprocess_expression_variants,
    preprocess_hla,
    preprocess_signature_variants,
    preprocess_small_mutations,
    preprocess_structural_variants,
    validate_report_content,
)
from pori_python.ipr.util import logger
from pori_python.types import IprFusionVariant, IprGeneVariant

DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data")
NON_EMPTY_STRING_NULLS = ["", None, np.nan, pd.NA]
EXPECTED_COSMIC = {"DBS9", "DBS11", "ID2", "ID7", "ID10", "SBS2", "SBS5", "DMMR"}
EXPECTED_HLA = {
    "HLA-A*02:01",
    "HLA-A*02",
    "HLA-A*30:01",
    "HLA-A*30",
    "HLA-B*27:01",
    "HLA-B*27",
    "HLA-B*15:01",
    "HLA-B*15",
    "HLA-C*03:03",
    "HLA-C*03",
    "HLA-C*06:02",
    "HLA-C*06",
}


def read_data_file(filename):
    pass


class TestPreProcessSmallMutations:
    def test_load_test_file(self) -> None:
        records = preprocess_small_mutations(
            pd.read_csv(os.path.join(DATA_DIR, "small_mutations.tab"), sep="\t").to_dict("records")
        )
        assert records
        assert len(records) == 2614

    def test_maintains_optional_fields(self):
        original = {
            "gene": "A1BG",
            "proteinChange": "p.V460M",
            "zygosity": "het",
            "tumourAltCount": 42,
            "tumourRefCount": 48,
            "hgvsProtein": "",
            "transcript": "ENST1000",
            "hgvsCds": "",
            "hgvsGenomic": "",
            "key": "02fe85a3477784b5ac0f8ecffb300d10",
            "variant": "blargh",
            "chromosome": "2",
            "startPosition": 1234,
        }
        records = preprocess_small_mutations([original])
        record = records[0]
        assert record["variantType"] == "mut"
        for col in original:
            assert col in record
        assert record["variant"] == "A1BG:p.V460M"
        assert "endPosition" in record
        assert record["endPosition"] == record["startPosition"]
        assert "tumourDepth" in record
        assert record["tumourDepth"] == 90

    def test_null(self):
        original = {
            "gene": "A1BG",
            "proteinChange": "p.V460M",
            "tumourAltCount": 42,
            "tumourRefCount": 48,
            "startPosition": 1234,
        }
        # Make sure TEST_KEYS are appropriate.
        # For some fields, like 'ref' and 'alt', NA is _not_ equivalent to a null string.
        TEST_KEYS = ["startPosition", "endPosition", "tumourAltCount", "tumourRefCount"]
        for key in TEST_KEYS:
            for null in NON_EMPTY_STRING_NULLS:
                small_mut = original.copy()
                small_mut[key] = null
                records = preprocess_small_mutations([small_mut])
                record = records[0]
                assert record["variantType"] == "mut"
                for col in original:
                    assert col in record
                assert record["variant"] == "A1BG:p.V460M"
                assert "endPosition" in record

    def test_load_small_mutations_probe(self) -> None:
        records = preprocess_small_mutations(
            pd.read_csv(os.path.join(DATA_DIR, "small_mutations_probe.tab"), sep="\t").to_dict(
                "records"
            )
        )
        assert records
        assert len(records) == 4
        assert records[0]["variantType"] == "mut"
        assert "variant" in records[0]


class TestPreProcessCopyVariants:
    def test_load_copy_variants(self) -> None:
        records = preprocess_copy_variants(
            pd.read_csv(os.path.join(DATA_DIR, "copy_variants.tab"), sep="\t").to_dict("records")
        )

        assert records
        assert len(records) == 4603
        assert records[0]["variantType"] == "cnv"
        assert "variant" in records[0]

    def test_add_chr_to_chrband(self) -> None:
        df1 = pd.read_csv(os.path.join(DATA_DIR, "copy_variants.tab"), sep="\t")
        df1 = df1.to_dict("records")
        records = preprocess_copy_variants(df1)
        assert records
        assert len(records) == 4603
        assert records[0]["chromosomeBand"] == "1q22.1"
        assert "chromosome" not in records[0]

    def test_add_int_chr_to_chrband(self) -> None:
        df1 = pd.read_csv(os.path.join(DATA_DIR, "copy_variants.tab"), sep="\t")
        df1["chromosome"] = df1["chromosome"].apply(lambda x: x.split("chr")[1])
        df1 = df1.to_dict("records")
        records = preprocess_copy_variants(df1)
        assert records
        assert len(records) == 4603
        assert records[0]["chromosomeBand"] == "1q22.1"
        assert "chromosome" not in records[0]

    def test_add_chr_to_chrband_if_chromosome_not_present(self) -> None:
        df1 = pd.read_csv(os.path.join(DATA_DIR, "copy_variants.tab"), sep="\t")
        df1["chr"] = df1["chromosome"].copy()
        df1.drop("chromosome", axis=1, inplace=True)
        df1 = df1.to_dict("records")
        records = preprocess_copy_variants(df1)
        assert records
        assert len(records) == 4603
        assert records[0]["chromosomeBand"] == "1q22.1"
        assert "chr" not in records[0]
        assert "chromosome" not in records[0]

    def test_do_not_add_chr_if_chr_already_in_chrband(self) -> None:
        df1 = pd.read_csv(os.path.join(DATA_DIR, "copy_variants.tab"), sep="\t")
        df2 = df1.copy()
        df1["chromosomeBand"] = df1["chromosomeBand"].apply(lambda x: "chr99" + x)
        df1 = df1.to_dict("records")
        df2["chromosomeBand"] = df2["chromosomeBand"].apply(lambda x: "99" + x)
        df2 = df2.to_dict("records")
        records = preprocess_copy_variants(df1)
        assert records
        assert len(records) == 4603
        assert records[0]["chromosomeBand"] == "chr99q22.1"
        assert "chr" not in records[0]  # make sure these cols are still getting removed
        assert "chromosome" not in records[0]
        records2 = preprocess_copy_variants(df2)
        assert records2[0]["chromosomeBand"] == "99q22.1"

    def test_no_error_if_chr_column_not_present(self) -> None:
        df1 = pd.read_csv(os.path.join(DATA_DIR, "copy_variants.tab"), sep="\t")
        df1.drop("chromosome", axis=1, inplace=True)
        df1 = df1.to_dict("records")
        records = preprocess_copy_variants(df1)
        assert records
        assert len(records) == 4603
        assert records[0]["chromosomeBand"] == "q22.1"

    def test_null(self):
        for kb_cat in list(INPUT_COPY_CATEGORIES.values()) + NON_EMPTY_STRING_NULLS:
            original = {"gene": "ERBB2", "kbCategory": kb_cat}
            for key in COPY_OPTIONAL:
                for null in NON_EMPTY_STRING_NULLS:
                    copy_var = original.copy()
                    copy_var[key] = null
                    records = preprocess_copy_variants([copy_var])
                    record = records[0]
                    assert record["variantType"] == "cnv"


class TestPreProcessSignatureVariants:
    def test_preprocess_cosmic(self) -> None:
        records = preprocess_cosmic(
            pd.read_csv(os.path.join(DATA_DIR, "cosmic_variants.tab"), sep="\t").to_dict("records")
        )
        assert records
        assert len(records) == len(EXPECTED_COSMIC)
        assert "variantTypeName" in records[0]
        assert "displayName" in records[0]

        signatureNames = {r.get("signatureName", "") for r in records}
        assert len(EXPECTED_COSMIC.symmetric_difference(signatureNames)) == 0

    def test_preprocess_hla(self) -> None:
        records = preprocess_hla(
            pd.read_csv(os.path.join(DATA_DIR, "hla_variants.tab"), sep="\t").to_dict("records")
        )
        assert records
        assert len(records) == len(EXPECTED_HLA)
        assert "variantTypeName" in records[0]
        assert "displayName" in records[0]

        signatureNames = {r.get("signatureName", "") for r in records}
        assert len(EXPECTED_HLA.symmetric_difference(signatureNames)) == 0

    def test_preprocess_signature_variants(self) -> None:
        records = preprocess_signature_variants(
            [
                *preprocess_cosmic(
                    pd.read_csv(os.path.join(DATA_DIR, "cosmic_variants.tab"), sep="\t").to_dict(
                        "records"
                    )
                ),
                *preprocess_hla(
                    pd.read_csv(os.path.join(DATA_DIR, "hla_variants.tab"), sep="\t").to_dict(
                        "records"
                    )
                ),
            ]
        )
        assert records
        assert len(records) == len(EXPECTED_COSMIC) + len(EXPECTED_HLA)
        assert "key" in records[0]


def test_load_structural_variants() -> None:
    records = preprocess_structural_variants(
        pd.read_csv(os.path.join(DATA_DIR, "fusions.tab"), sep="\t").to_dict("records")
    )
    assert records
    assert len(records) == 7
    assert records[0]["variantType"] == "sv"
    assert "variant" in records[0]


def test_load_expression_variants() -> None:
    records = preprocess_expression_variants(
        pd.read_csv(os.path.join(DATA_DIR, "expression.tab"), sep="\t").to_dict("records")
    )
    assert records
    assert len(records) == 4603
    assert records[0]["variantType"] == "exp"
    assert "variant" in records[0]


class TestCheckVariantLinks:
    def test_sm_missing_copy_empty_ok(self) -> None:
        genes = check_variant_links(
            small_mutations=[IprGeneVariant({"gene": "KRAS"})],  # type: ignore
            copy_variants=[],
            expression_variants=[IprGeneVariant({"gene": "KRAS", "variant": ""})],  # type: ignore
            structural_variants=[],
        )
        assert genes == {"KRAS"}

    def test_sm_missing_exp_empty_ok(self) -> None:
        genes = check_variant_links(
            small_mutations=[IprGeneVariant({"gene": "KRAS"})],  # type: ignore
            copy_variants=[IprGeneVariant({"gene": "KRAS", "variant": ""})],  # type: ignore
            expression_variants=[],
            structural_variants=[],
        )
        assert genes == {"KRAS"}

    def test_sm_missing_copy(self) -> None:
        with mock.patch.object(logger, "debug") as mock_debug:
            check_variant_links(
                small_mutations=[IprGeneVariant({"gene": "KRAS"})],  # type: ignore
                copy_variants=[IprGeneVariant({"gene": "CDK", "variant": ""})],  # type: ignore
                expression_variants=[IprGeneVariant({"gene": "KRAS", "variant": ""})],  # type: ignore
                structural_variants=[],
            )
            assert mock_debug.called

    def test_sm_missing_exp(self) -> None:
        with mock.patch.object(logger, "debug") as mock_debug:
            check_variant_links(
                small_mutations=[IprGeneVariant({"gene": "KRAS"})],  # type: ignore
                copy_variants=[IprGeneVariant({"gene": "KRAS", "variant": ""})],  # type: ignore
                expression_variants=[IprGeneVariant({"gene": "CDK", "variant": ""})],  # type: ignore
                structural_variants=[],
            )
            assert mock_debug.called

    def test_with_valid_inputs(self) -> None:
        genes = check_variant_links(
            small_mutations=[IprGeneVariant({"gene": "KRAS"})],  # type: ignore
            copy_variants=[
                IprGeneVariant({"gene": "KRAS", "variant": ""}),  # type: ignore
                IprGeneVariant({"gene": "CDK", "variant": ""}),  # type: ignore
            ],
            expression_variants=[IprGeneVariant({"gene": "KRAS", "variant": ""})],  # type: ignore
            structural_variants=[],
        )
        assert genes == {"KRAS"}

    def test_copy_missing_exp(self) -> None:
        with mock.patch.object(logger, "debug") as mock_debug:
            check_variant_links(
                small_mutations=[],
                copy_variants=[
                    IprGeneVariant({"gene": "BRAF", "variant": "copy gain"}),  # type: ignore
                    IprGeneVariant({"gene": "KRAS", "variant": ""}),  # type: ignore
                ],
                expression_variants=[IprGeneVariant({"gene": "KRAS", "variant": ""})],  # type: ignore
                structural_variants=[],
            )
            assert mock_debug.called

    def test_exp_missing_copy(self) -> None:
        with mock.patch.object(logger, "debug") as mock_debug:
            check_variant_links(
                small_mutations=[],
                copy_variants=[IprGeneVariant({"gene": "KRAS", "variant": ""})],  # type: ignore
                expression_variants=[
                    IprGeneVariant({"gene": "BRAF", "variant": "increased expression"})  # type: ignore
                ],
                structural_variants=[],
            )
            assert mock_debug.called


class TestCreateGraphkbSvNotation:
    def test_both_genes_and_exons(self) -> None:
        notation = create_graphkb_sv_notation(
            IprFusionVariant({"gene1": "A", "gene2": "B", "exon1": 1, "exon2": 2})  # type: ignore
        )
        assert notation == "(A,B):fusion(e.1,e.2)"

    def test_one_exon_missing(self) -> None:
        notation = create_graphkb_sv_notation(
            IprFusionVariant({"gene1": "A", "gene2": "B", "exon1": "", "exon2": 2})  # type: ignore
        )
        assert notation == "(A,B):fusion(e.?,e.2)"

    def test_one_gene_missing(self) -> None:
        notation = create_graphkb_sv_notation(
            IprFusionVariant({"gene1": "A", "gene2": "", "exon1": 1, "exon2": 2})  # type: ignore
        )
        assert notation == "(A,?):fusion(e.1,e.2)"

    def test_first_gene_missing(self) -> None:
        notation = create_graphkb_sv_notation(
            IprFusionVariant({"gene1": "", "gene2": "B", "exon1": 1, "exon2": 2})  # type: ignore
        )
        assert notation == "(B,?):fusion(e.2,e.1)"

    def test_no_genes_error(self) -> None:
        with pytest.raises(ValueError):
            create_graphkb_sv_notation(
                IprFusionVariant({"gene1": "", "gene2": "", "exon1": 1, "exon2": 2, "key": "x"})  # type: ignore
            )


class TestCheckComparators:
    def test_missing_disease_expression_error(self):
        content = {"comparators": [{"analysisRole": "expression (primary site)"}]}
        variants = [{}]

        with pytest.raises(ValueError):
            check_comparators(content, variants)

    def test_missing_primary_expression_error(self):
        content = {"comparators": [{"analysisRole": "expression (disease)"}]}
        variants = [{"primarySiteFoldChange": 1}]

        with pytest.raises(ValueError):
            check_comparators(content, variants)

    def test_missing_biopsy_expression_error(self):
        content = {"comparators": [{"analysisRole": "expression (disease)"}]}
        variants = [{"biopsySitePercentile": 1}]

        with pytest.raises(ValueError):
            check_comparators(content, variants)

    def test_expression_not_required_without_variants(self):
        content = {"comparators": []}
        variants = []

        assert check_comparators(content, variants) is None

    def test_missing_mutation_burden(self):
        content = {
            "comparators": [{"analysisRole": "mutation burden (secondary)"}],
            "images": [{"key": "mutationBurden.density_snv.primary"}],
        }
        variants = []

        with pytest.raises(ValueError):
            check_comparators(content, variants)


@pytest.mark.parametrize("example_name", ["no_variants", "sm_and_exp", "sm_only"])
def test_valid_json_inputs(example_name: str):
    with open(os.path.join(DATA_DIR, "json_examples", f"{example_name}.json"), "r") as fh:
        content = json.load(fh)
    validate_report_content(content)
