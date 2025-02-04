"""
handles annotating variants with annotation information from graphkb
"""

from __future__ import annotations

from requests.exceptions import HTTPError

from pandas import isnull
from tqdm import tqdm
from typing import Dict, List, Sequence, cast

from pori_python.graphkb import GraphKBConnection
from pori_python.graphkb import match as gkb_match
from pori_python.graphkb.match import INPUT_COPY_CATEGORIES
from pori_python.graphkb.statement import get_statements_from_variants
from pori_python.graphkb.util import FeatureNotFoundError
from pori_python.types import (
    Hashabledict,
    IprCopyVariant,
    IprExprVariant,
    IprStructuralVariant,
    KbMatch,
    Statement,
    Variant,
)

from .constants import TMB_HIGH_CATEGORY
from .ipr import convert_statements_to_alterations
from .util import convert_to_rid_set, logger

REPORTED_COPY_VARIANTS = (INPUT_COPY_CATEGORIES.AMP, INPUT_COPY_CATEGORIES.DEEP)


def get_second_pass_variants(
    graphkb_conn: GraphKBConnection, statements: List[Statement]
) -> List[Variant]:
    """Given a list of statements that have been matched, convert these to
    new category variants to be used in a second-pass matching.
    """
    # second-pass matching
    all_inferred_matches: Dict[str, Variant] = {}
    inferred_variants = {
        (s["subject"]["@rid"], s["relevance"]["name"])
        for s in statements
        if s["subject"] and s["subject"]["@class"] in ("Feature", "Signature")
    }

    for reference1, variant_type in inferred_variants:
        variants = gkb_match.match_category_variant(
            graphkb_conn, reference1, variant_type
        )

        for variant in variants:
            all_inferred_matches[variant["@rid"]] = variant
    inferred_matches: List[Variant] = list(all_inferred_matches.values())
    return inferred_matches


def get_ipr_statements_from_variants(
    graphkb_conn: GraphKBConnection, matches: List[Variant], disease_name: str
) -> List[KbMatch]:
    """IPR upload formatted GraphKB statements from the list of variants.

    Matches to GraphKB statements from the list of input variants. From these results matches
    again with the inferred variants. Then returns the results formatted for upload to IPR
    """
    if not matches:
        return []
    rows = []

    statements = get_statements_from_variants(graphkb_conn, matches)
    existing_statements = {s["@rid"] for s in statements}

    for ipr_row in convert_statements_to_alterations(
        graphkb_conn, statements, disease_name, convert_to_rid_set(matches)
    ):
        rows.append(ipr_row)

    # second-pass matching
    inferred_matches = get_second_pass_variants(graphkb_conn, statements)

    inferred_statements = [
        s
        for s in get_statements_from_variants(graphkb_conn, inferred_matches)
        if s["@rid"]
        not in existing_statements  # do not duplicate if non-inferred match
    ]

    for ipr_row in convert_statements_to_alterations(
        graphkb_conn,
        inferred_statements,
        disease_name,
        convert_to_rid_set(inferred_matches),
    ):
        ipr_row["kbData"]["inferred"] = True
        rows.append(ipr_row)

    return rows


def annotate_expression_variants(
    graphkb_conn: GraphKBConnection,
    variants: List[IprExprVariant],
    disease_name: str,
    show_progress: bool = False,
) -> List[KbMatch]:
    """Annotate expression variants with GraphKB in the IPR alterations format.

    Args:
        graphkb_conn: the graphkb api connection object
        variants: list of variants

    Returns:
        list of kbMatches records for IPR
    """
    skipped = 0
    alterations = []
    problem_genes = set()
    logger.info(f"Starting annotation of {len(variants)} expression category_variants")
    iterfunc = tqdm if show_progress else iter
    for row in iterfunc(variants):
        gene = row["gene"]
        variant = row["variant"]

        if not variant:
            skipped += 1
            logger.debug(f"Skipping malformed Expression {gene}: {row}")
            continue
        try:
            matches = gkb_match.match_expression_variant(graphkb_conn, gene, variant)
            for ipr_row in get_ipr_statements_from_variants(
                graphkb_conn, matches, disease_name
            ):
                ipr_row["variant"] = row["key"]
                ipr_row["variantType"] = row.get("variantType", "exp")
                alterations.append(ipr_row)
        except FeatureNotFoundError as err:
            problem_genes.add(gene)
            logger.debug(f"Unrecognized gene ({gene} {variant}): {err}")
        except ValueError as err:
            logger.error(f"failed to match variants ({gene} {variant}): {err}")

    if skipped:
        logger.info(f"skipped matching {skipped} expression information rows")
    if problem_genes:
        logger.error(f"gene finding failures for expression {sorted(problem_genes)}")
        logger.error(f"gene finding falure for {len(problem_genes)} expression genes")
    logger.info(
        f"matched {len(variants)} expression variants to {len(alterations)} graphkb annotations"
    )
    return alterations


def annotate_copy_variants(
    graphkb_conn: GraphKBConnection,
    variants: List[IprCopyVariant],
    disease_name: str,
    show_progress: bool = False,
) -> List[KbMatch]:
    """Annotate allowed copy variants with GraphKB in the IPR alterations format.

    Args:
        graphkb_conn: the graphkb api connection object
        variants: list of variants

    Returns:
        list of kbMatches records for IPR
    """
    skipped = 0
    alterations = []
    problem_genes = set()

    logger.info(f"Starting annotation of {len(variants)} copy category_variants")
    iterfunc = tqdm if show_progress else iter
    for row in iterfunc(variants):
        gene = row["gene"]
        variant = row["variant"]

        if variant not in REPORTED_COPY_VARIANTS:
            # https://www.bcgsc.ca/jira/browse/GERO-77
            skipped += 1
            logger.debug(
                f"Dropping {gene} copy change '{variant}' - not in REPORTED_COPY_VARIANTS"
            )
            continue
        try:
            matches = gkb_match.match_copy_variant(graphkb_conn, gene, variant)
            for ipr_row in get_ipr_statements_from_variants(
                graphkb_conn, matches, disease_name
            ):
                ipr_row["variant"] = row["key"]
                ipr_row["variantType"] = row.get("variantType", "cnv")
                alterations.append(ipr_row)
        except FeatureNotFoundError as err:
            problem_genes.add(gene)
            logger.debug(f"Unrecognized gene ({gene} {variant}): {err}")
        except ValueError as err:
            logger.error(f"failed to match variants ({gene} {variant}): {err}")

    if skipped:
        logger.info(
            f"skipped matching {skipped} copy number variants not in {REPORTED_COPY_VARIANTS}"
        )
    if problem_genes:
        logger.error(f"gene finding failures for copy variants {sorted(problem_genes)}")
        logger.error(
            f"gene finding failure for {len(problem_genes)} copy variant genes"
        )
    logger.info(
        f"matched {len(variants)} copy category variants to {len(alterations)} graphkb annotations"
    )
    return alterations


def annotate_positional_variants(
    graphkb_conn: GraphKBConnection,
    variants: Sequence[IprStructuralVariant] | Sequence[Hashabledict],
    disease_name: str,
    show_progress: bool = False,
) -> List[Hashabledict]:
    """Annotate SNP, INDEL or fusion variant calls with GraphKB and return in IPR match format.

    Hashable type is required to turn lists into sets.
    Args:
        graphkb_conn (GraphKBConnection): the graphkb api connection object
        variants (list.<dict>): list of variants. Defaults to [].
        disease_name (str): GraphKB disease name for statement matching.  'cancer' is most general
        show_progress (bool): Progressbar displayed for long runs.

    Returns:
        Hashable list of kbMatches records for IPR
    """
    VARIANT_KEYS = ("variant", "hgvsProtein", "hgvsCds", "hgvsGenomic")
    errors = 0
    alterations: List[Hashabledict] = []
    problem_genes = set()

    iterfunc = tqdm if show_progress else iter
    for row in iterfunc(variants):
        if not row.get("gene") and (not row.get("gene1") or not row.get("gene2")):
            # https://www.bcgsc.ca/jira/browse/GERO-56?focusedCommentId=1234791&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-1234791
            # should not match single gene SVs
            continue

        for var_key in VARIANT_KEYS:
            variant = row.get(var_key)
            matches = []
            if not variant or isnull(variant):
                continue
            try:
                try:
                    matches = gkb_match.match_positional_variant(graphkb_conn, variant)
                except HTTPError as parse_err:
                    # DEVSU-1885 - fix malformed single deletion described as substitution of blank
                    # eg. deletion described as substitution with nothing: 'chr1:g.150951027T>'
                    if (
                        variant[-1] == ">"
                        and "g." in variant
                        and variant[-2].isalpha()
                        and variant[-3].isnumeric()
                    ):
                        logger.warning(
                            f"Assuming malformed deletion variant {variant} is {variant[:-2] + 'del'}"
                        )
                        variant = variant[:-2] + "del"
                        matches = gkb_match.match_positional_variant(
                            graphkb_conn, variant
                        )
                    else:
                        raise parse_err

                for ipr_row in get_ipr_statements_from_variants(
                    graphkb_conn, matches, disease_name
                ):
                    ipr_row["variant"] = row["key"]
                    ipr_row["variantType"] = row.get(
                        "variantType", "mut" if row.get("gene") else "sv"
                    )
                    alterations.append(Hashabledict(ipr_row))

            except FeatureNotFoundError as err:
                logger.debug(f"failed to match positional variants ({variant}): {err}")
                errors += 1
                if "gene" in row:
                    problem_genes.add(row["gene"])
                elif "gene1" in row and f"({row['gene1']})" in str(err):
                    problem_genes.add(row["gene1"])
                elif "gene2" in row and f"({row['gene2']})" in str(err):
                    problem_genes.add(row["gene2"])
                elif "gene1" in row and "gene2" in row:
                    problem_genes.add(row["gene1"])
                    problem_genes.add(row["gene2"])
                else:
                    raise err
            except HTTPError as err:
                errors += 1
                logger.error(f"failed to match positional variants ({variant}): {err}")

    if problem_genes:
        logger.error(f"gene finding failures for {sorted(problem_genes)}")
        logger.error(
            f"{len(problem_genes)} gene finding failures for positional variants"
        )
    if errors:
        logger.error(f"skipped {errors} positional variants due to errors")

    # drop duplicates
    alterations = list(set(alterations))

    variant_types = ", ".join(sorted(set([alt["variantType"] for alt in alterations])))
    logger.info(
        f"matched {len(variants)} {variant_types} positional variants to {len(alterations)} graphkb annotations"
    )

    return alterations


def annotate_msi(
    graphkb_conn: GraphKBConnection,
    disease_name: str = "cancer",
    msi_category: str = "microsatellite instability",
) -> List[KbMatch]:
    """Annotate microsatellite instablity from GraphKB in the IPR alterations format.

    Match to GraphKb Category variants with similar names
    Args:
        graphkb_conn: the graphkb api connection object
        msi_category: such as 'microsatellite instability'

    Returns:
        list of kbMatches records for IPR
    """
    gkb_matches = []
    msi_categories = graphkb_conn.query(
        {
            "target": {
                "target": "CategoryVariant",
                "filters": {
                    "reference1": {
                        "target": "Signature",
                        "filters": {"name": msi_category},
                    }
                },
            },
            "queryType": "similarTo",
            "returnProperties": ["@rid", "displayName"],
        }
    )
    if msi_categories:
        msi_variants = [cast(Variant, var) for var in msi_categories]
        for ipr_row in get_ipr_statements_from_variants(
            graphkb_conn, msi_variants, disease_name
        ):
            ipr_row["variant"] = msi_category
            ipr_row["variantType"] = "msi"
            gkb_matches.append(ipr_row)
    return gkb_matches


def annotate_tmb(
    graphkb_conn: GraphKBConnection,
    disease_name: str = "cancer",
    category: str = TMB_HIGH_CATEGORY,
) -> List[KbMatch]:
    """Annotate Tumour Mutation Burden (tmb) categories from GraphKB in the IPR alterations format.

    Match to GraphKb Category variants with similar names
    Args:
        graphkb_conn: the graphkb api connection object
        disease_name: oncotree disease name for graphkb matching.
        category: such as 'high mutation burden'

    Returns:
        list of kbMatches records for IPR
    """
    gkb_matches = []
    categories = graphkb_conn.query(
        {
            "target": {
                "target": "CategoryVariant",
                "filters": {
                    "reference1": {
                        "target": "Signature",
                        "filters": {
                            "OR": [{"name": category}, {"displayName": category}]
                        },
                    }
                },
            },
            "queryType": "similarTo",
            "returnProperties": ["@rid", "displayName"],
        }
    )
    if categories:
        cat_variants = [cast(Variant, var) for var in categories]
        for ipr_row in get_ipr_statements_from_variants(
            graphkb_conn, cat_variants, disease_name
        ):
            ipr_row["variant"] = category
            ipr_row["variantType"] = "tmb"
            gkb_matches.append(ipr_row)
    return gkb_matches
