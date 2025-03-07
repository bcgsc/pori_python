"""
Contains functions specific to formatting reports for IPR that are unlikely to be used
by other reporting systems
"""

from __future__ import annotations
from itertools import product
from copy import copy
from typing import Dict, Iterable, List, Sequence, Set, Tuple, cast
import uuid
from pori_python.graphkb import GraphKBConnection
from pori_python.graphkb import statement as gkb_statement
from pori_python.graphkb import vocab as gkb_vocab
from pori_python.types import (
    Hashabledict,
    ImageDefinition,
    IprFusionVariant,
    IprGene,
    IprVariant,
    KbMatch,
    Statement,
    Variant,
    KbVariantMatch,
    KbMatchedStatement,
    KbMatchedStatementConditionSet,
    KbMatchSections,
)

from .constants import GERMLINE_BASE_TERMS, VARIANT_CLASSES
from .util import find_variant, logger


def display_evidence_levels(statement: Statement) -> str:
    result = []
    for evidence_level in statement.get("evidenceLevel", []) or []:
        if isinstance(evidence_level, str):
            result.append(evidence_level)
        elif "displayName" in evidence_level:
            result.append(evidence_level["displayName"])
    return ";".join(sorted(result))


def filter_structural_variants(
    structural_variants: List[IprFusionVariant],
    kb_matches: List[KbMatch] | List[Hashabledict],
    gene_annotations: List[IprGene],
) -> List[IprFusionVariant]:
    """
    Filter structural variants to remove non-high quality events unless they are matched/annotated or
    they involve a gene that is a known fusion partner
    """
    matched_svs = {match["variant"] for match in kb_matches if match["variantType"] == "sv"}
    fusion_genes = {
        gene["name"] for gene in gene_annotations if gene.get("knownFusionPartner", False)
    }

    result = []

    for structural_variant in structural_variants:
        if any(
            [
                structural_variant["highQuality"],
                structural_variant["key"] in matched_svs,
                structural_variant["gene1"] in fusion_genes,
                structural_variant["gene2"] in fusion_genes,
            ]
        ):
            result.append(structural_variant)
    return result


def get_evidencelevel_mapping(graphkb_conn: GraphKBConnection) -> Dict[str, str]:
    """IPR evidence level equivalents of GraphKB evidence returned as a dictionary.

    Args:
        graphkb_conn (GraphKBConnection): the graphkb api connection object

    Returns:
        dictionary mapping all EvidenceLevel RIDs to corresponding IPR EvidenceLevel displayName
    """
    # Get all EvidenceLevel from GraphKB
    # Note: not specifying any returnProperties allows for retreiving in/out_CrossReferenceOf
    evidence_levels = graphkb_conn.query({"target": "EvidenceLevel"})

    # Map EvidenceLevel RIDs to list of incoming CrossReferenceOf
    evidence_levels_mapping = dict(
        map(lambda d: (d["@rid"], d.get("in_CrossReferenceOf", [])), evidence_levels)
    )

    # Filter IPR EvidenceLevel and map each outgoing CrossReferenceOf to displayName
    ipr_source_rid = graphkb_conn.get_source("ipr")["@rid"]
    ipr_evidence_levels = filter(lambda d: d.get("source") == ipr_source_rid, evidence_levels)
    cross_references_mapping: Dict[str, str] = dict()
    ipr_rids_to_displayname: Dict[str, str] = dict()
    for level in ipr_evidence_levels:
        d = map(lambda i: (i, level["displayName"]), level.get("out_CrossReferenceOf", []))  # type: ignore
        cross_references_mapping.update(d)
        ipr_rids_to_displayname[level["@rid"]] = level["displayName"]  # type: ignore

    # Update EvidenceLevel mapping to corresponding IPR EvidenceLevel displayName
    def link_refs(refs) -> Tuple[str, str]:
        for rid in refs[1]:
            if cross_references_mapping.get(rid):
                return (refs[0], cross_references_mapping[rid])
        if refs[0] in ipr_rids_to_displayname:  # self-referencing IPR levels
            return (refs[0], ipr_rids_to_displayname[refs[0]])
        return (refs[0], "")

    evidence_levels_mapping = dict(map(link_refs, evidence_levels_mapping.items()))
    evidence_levels_mapping[""] = ""

    return evidence_levels_mapping  # type: ignore


# TODO for DEVSU-2550
def convert_statements_to_alterations(
    graphkb_conn: GraphKBConnection,
    statements: List[Statement],
    disease_matches: List[str],
    variant_matches: Iterable[str],
) -> List[KbMatch]:
    """Convert statements matched from graphkb into IPR equivalent representations.

    Args:
        graphkb_conn: the graphkb connection object
        statements: list of statement records from graphkb
        disease_matches: GraphKB disease RIDs
        variant_matches: the list of RIDs the variant matched for these statements

    Raises:
        ValueError: could not find the disease type in GraphKB

    Returns:
        IPR graphkb row representations

    Notes:
        - only report disease matched prognostic markers https://www.bcgsc.ca/jira/browse/GERO-72 and GERO-196
    """
    rows = []
    ev_map = get_evidencelevel_mapping(graphkb_conn)
    # GERO-318 - add all IPR-A evidence equivalents to the approvedTherapy flag
    approved = set([ev for (ev, ipr) in ev_map.items() if ipr == "IPR-A"])

    # get the recruitment status for any trial associated with a statement
    clinical_trials = [
        s["subject"]["@rid"] for s in statements if s["subject"]["@class"] == "ClinicalTrial"
    ]
    recruitment_statuses = {}
    if clinical_trials:
        clinical_trials = list(set(clinical_trials))
        for rid in clinical_trials:
            query_result = graphkb_conn.query(
                {
                    "target": {"target": "ClinicalTrial", "filters": {"@rid": rid}},
                    "returnProperties": ["@rid", "recruitmentStatus"],
                }
            )
            if query_result:
                recruitment_statuses[rid] = query_result[0]["recruitmentStatus"]  # type: ignore

    for statement in statements:
        variants = [
            cast(Variant, c) for c in statement["conditions"] if c["@class"] in VARIANT_CLASSES
        ]
        diseases = [c for c in statement["conditions"] if c["@class"] == "Disease"]
        disease_match = len(diseases) == 1 and diseases[0]["@rid"] in disease_matches
        pmid = ";".join([e["displayName"] for e in statement["evidence"]])

        ipr_section = gkb_statement.categorize_relevance(
            graphkb_conn, statement["relevance"]["@rid"]
        )
        approved_therapy = False
        if ipr_section == "therapeutic":
            for level in statement["evidenceLevel"] or []:
                if level["@rid"] in approved:
                    approved_therapy = True
                    break

        if ipr_section == "prognostic" and not disease_match:
            continue  # GERO-72 / GERO-196

        evidence_level_str = display_evidence_levels(statement)
        evidence_levels = statement.get("evidenceLevel") or []
        ipr_evidence_levels = [ev_map[el.get("@rid", "")] for el in evidence_levels if el]
        ipr_evidence_levels_str = ";".join(sorted(set([el for el in ipr_evidence_levels])))

        for variant in variants:
            if variant["@rid"] not in variant_matches:
                continue
            row = KbMatch(
                {
                    "approvedTherapy": approved_therapy or False,
                    "category": ipr_section or "unknown",
                    "context": (
                        statement["subject"]["displayName"] if statement["subject"] else ""
                    ),
                    "kbContextId": (statement["subject"]["@rid"] if statement["subject"] else ""),
                    "disease": ";".join(sorted(d.get("displayName", "") for d in diseases)),
                    "evidenceLevel": evidence_level_str or "",
                    "iprEvidenceLevel": ipr_evidence_levels_str or "",
                    "kbStatementId": statement["@rid"],
                    "kbVariant": str(variant.get("displayName", "")) or "",
                    "variant": str(variant.get("displayName", "")) or "",
                    "variantType": "",
                    "kbVariantId": variant["@rid"],
                    "matchedCancer": disease_match,
                    "reference": pmid,
                    "relevance": statement["relevance"]["displayName"],
                    "kbRelevanceId": statement["relevance"]["@rid"],
                    "externalSource": (
                        str(statement["source"].get("displayName", ""))
                        if statement["source"]
                        else ""
                    ),
                    "requiredKbMatches": [item["@rid"] for item in variants],
                    "externalStatementId": statement.get("sourceId", "") or "",
                    "reviewStatus": statement.get("reviewStatus", "") or "",
                    "kbData": {},
                }
            )
            if statement["relevance"]["name"] == "eligibility":
                row["kbData"]["recruitment_status"] = recruitment_statuses.get(
                    row["kbContextId"], "not found"
                )
            rows.append(row)
    return rows


def select_expression_plots(
    kb_matches: List[KbMatch] | List[Hashabledict], all_variants: Sequence[IprVariant]
) -> List[ImageDefinition]:
    """
    Given the list of expression variants, determine which expression
    historgram plots should be included in the IPR upload. This filters them
    based on the graphkb annotations to avoid loading more images than are required

    Args:
        kb_matches: the IPR graphkb annoations for all variants
        expression_variants: the list of expression variants loaded

    Returns:
        list of expression images to be loaded by IPR
    """

    selected_variants = {
        (match["variantType"], match["variant"])
        for match in kb_matches
        if match["category"] == "therapeutic"
    }
    images_by_gene: Dict[str, ImageDefinition] = {}
    selected_genes = set()
    for variant in all_variants:
        if (variant["variantType"], variant["key"]) in selected_variants:
            for key in ["gene", "gene1", "gene2"]:
                gene = variant.get(key)
                if gene:
                    selected_genes.add(str(gene))
        gene = str(variant.get("gene", ""))
        hist = str(variant.get("histogramImage", ""))
        if hist:
            images_by_gene[gene] = ImageDefinition({"key": f"expDensity.{gene}", "path": hist})
    return [images_by_gene[gene] for gene in selected_genes if gene in images_by_gene]


def create_key_alterations(
    kb_matches: List[Hashabledict], all_variants: Sequence[IprVariant]
) -> Tuple[List[Dict], Dict]:
    """Create the list of significant variants matched by the KB.

    This list of matches is also used to create the variant counts.
    """
    alterations = []
    type_mapping = {
        "mut": "smallMutations",
        "cnv": "CNVs",
        "sv": "SVs",
        "exp": "expressionOutliers",
    }
    counts: Dict[str, Set] = {v: set() for v in type_mapping.values()}
    skipped_variant_types = []
    for kb_match in kb_matches:
        variant_type = kb_match["variantType"]
        variant_key = kb_match["variant"]
        if kb_match["category"] == "unknown":
            continue

        if variant_type not in type_mapping.keys():
            if variant_type not in skipped_variant_types:
                skipped_variant_types.append(variant_type)
                logger.warning(
                    f"No summary key alterations for {variant_type}.  Skipping {variant_key}"
                )
            continue
        try:
            variant = find_variant(all_variants, variant_type, variant_key)
        except KeyError as err:
            logger.error(err)
            logger.error(f"No variant match found for {variant_key}")
            continue

        counts[type_mapping[variant_type]].add(variant_key)

        if variant_type == "exp":
            alterations.append(f'{variant.get("gene","")} ({variant.get("expressionState")})')
        elif variant_type == "cnv":
            alterations.append(f'{variant.get("gene","")} ({variant.get("cnvState")})')
        # only show germline if relevant
        elif kb_match["category"] in GERMLINE_BASE_TERMS and variant.get("germline"):
            alterations.append(f"germline {variant['variant']}")
        else:
            alterations.append(variant["variant"])

    counted_variants = set.union(*counts.values())
    counts["variantsUnknown"] = set()

    # count the un-matched variants
    for variant in all_variants:
        if variant["variant"] and variant["key"] not in counted_variants:
            counts["variantsUnknown"].add(variant["key"])

    return (
        [{"geneVariant": alt} for alt in set(alterations)],
        {k: len(v) for k, v in counts.items()},
    )


def germline_kb_matches(
    kb_matches: List[Hashabledict],
    all_variants: Sequence[IprVariant],
    assume_somatic: bool = True,
) -> List[Hashabledict]:
    """Filter kb_matches for matching to germline or somatic events using the 'germline' optional property.

    Statements related to pharmacogenomic toxicity or cancer predisposition are only relevant if
    the variant is present in the germline of the patient.
    Other statements, such as diagnostic or recurrent oncogenic mutations, are only relevant as
    somatic events in cancer.  Germline variants are excluded from these matches.

    Params:
        kb_matches: KbMatch statements to be filtered.  'variant' properties must match 'key' in all_variants.
        all_variants: IprVariants, with a 'germline' property, that were used for kb_matches creation.
        assume_somatic: Whether to assume somatic or germline when no 'germline' property exists in the variant.
    Returns:
        filtered list of kb_matches
    """
    ret_list = []
    germ_alts = [alt for alt in kb_matches if alt["category"] in GERMLINE_BASE_TERMS]
    somatic_alts = [alt for alt in kb_matches if alt not in germ_alts]
    if germ_alts:
        logger.info(f"checking germline status of {GERMLINE_BASE_TERMS}")
        for alt in germ_alts:
            var_list = [v for v in all_variants if v["key"] == alt["variant"]]
            germline_var_list = [v for v in var_list if v.get("germline")]
            unknown_var_list = [v for v in var_list if "germline" not in v]
            if germline_var_list:
                logger.debug(
                    f"germline kbStatementId:{alt['kbStatementId']}: {alt['kbVariant']} {alt['category']}"
                )
                ret_list.append(alt)
            elif unknown_var_list:
                logger.warning(
                    f"germline no data fail for: {alt['kbStatementId']}: {alt['kbVariant']} {alt['category']}"
                )
                if not assume_somatic:
                    logger.debug(
                        f"Keeping unverified match to germline kbStatementId:{alt['kbStatementId']}: {alt['kbVariant']} {alt['category']}"
                    )
                    ret_list.append(alt)
                else:
                    logger.debug(
                        f"Dropping unverified match to germline kbStatementId:{alt['kbStatementId']}: {alt['kbVariant']} {alt['category']}"
                    )
            else:
                logger.debug(
                    f"Dropping somatic match to germline kbStatementId:{alt['kbStatementId']}: {alt['kbVariant']} {alt['category']}"
                )
    if somatic_alts:
        # Remove any matches to germline events
        for alt in somatic_alts:
            var_list = [v for v in all_variants if v["key"] == alt["variant"]]
            somatic_var_list = [v for v in var_list if not v.get("germline", not assume_somatic)]
            if var_list and not somatic_var_list:
                logger.debug(
                    f"Dropping germline match to somatic statement kbStatementId:{alt['kbStatementId']}: {alt['kbVariant']} {alt['category']}"
                )
            elif somatic_var_list:
                ret_list.append(alt)  # match to somatic variant
            else:
                ret_list.append(alt)  # alteration not in any specific keys matches to check.

    return ret_list


def multi_variant_filtering(
    graphkb_conn: GraphKBConnection,
    gkb_matches: List[KbMatch],
    excludedTypes: List[str] = ["wildtype"],
) -> List[KbMatch]:
    """Filters out GraphKB matches that doesn't match to all required variants on multi-variant statements

    DEVSU-2477
    GKB Statements can be conditional to more than one variant, with implicit 'AND' operator. Since variants
    are matched only one at a time, any multi-variant statement get matched if one of their conditional
    variants is matching the observed ones, making de facto an 'OR' operator between conditions. The current
    function is filtering out these incomplete matches.

    Note: Wildtype variants are not taken into account at the moment.

    Params:
        graphkb_conn: the graphkb connection object
        gkb_matches: KbMatch statements to be filtered
        excludedTypes: List of variant type terms to exclude from filtering. Default to Wildtype
    Returns:
        filtered list of KbMatch statements
    """
    # All matching statements & variants (GKB RIDs)
    matching_statement_rids = {match["kbStatementId"] for match in gkb_matches}
    matching_variant_rids = {match["kbVariantId"] for match in gkb_matches}

    # Get conditions detail on all matching statements
    res = graphkb_conn.post(
        uri="query",
        data={
            "target": "Statement",
            "filters": {
                "@rid": list(matching_statement_rids),
                "operator": "IN",
            },
            "history": True,
            "returnProperties": [
                "@rid",
                "conditions.@rid",
                "conditions.@class",
                "conditions.type",
            ],
        },
    )
    statements = res["result"]

    # Get set of excluded Vocabulary RIDs for variant types
    excluded = {}
    if len(excludedTypes) != 0 and excludedTypes[0] != "":
        excluded = gkb_vocab.get_terms_set(graphkb_conn, excludedTypes)

    # Mapping statements to their conditional variants
    # (discarding non-variant conditions & variant conditions from excluded types)
    statement_to_variants = {}
    for statement in statements:
        statement_to_variants[statement["@rid"]] = {
            el["@rid"]
            for el in statement["conditions"]
            if (el["@class"] in VARIANT_CLASSES and el.get("type", "") not in excluded)
        }

    # Set of statements with complete matching
    complete_matching_statements = {
        statementRid
        for statementRid, variantRids in statement_to_variants.items()
        if variantRids.issubset(matching_variant_rids)
    }

    # Filtering out incompleted matches of gkb_matches
    return [
        match for match in gkb_matches if match["kbStatementId"] in complete_matching_statements
    ]


def get_kb_variants(
    gkb_matches: List[KbMatch] | List[Hashabledict],
) -> List[KbVariantMatch]:
    """Extracts the set of distinct kb variant records from the input
    list of gkb_matches records, which combine statement and variant matches.

    Params:
        gkb_matches: KbMatch statements to be processed
    Returns:
        set of distinct kbVariant records
    """
    kbVariants = {}
    for item in gkb_matches:
        kbv = KbVariantMatch(
            {
                "kbVariant": item["kbVariant"],
                "variant": item["variant"],
                "variantType": item["variantType"],
                "kbVariantId": item["kbVariantId"],
            }
        )
        kbVariants[str(kbv)] = kbv

    return [item for item in kbVariants.values()]


def get_kb_matched_statements(
    gkb_matches: List[KbMatch] | List[Hashabledict],
) -> List[KbMatchedStatement]:
    """Extracts the set of distinct kb statement records from the input
    list of gkb_matches records, which combine statement and variant matches.

    Params:
        gkb_matches: KbMatch statements to be processed
    Returns:
        set of distinct kbMatchedStatement records
    """
    kbMatchedStatements = {}
    kbs_keys = KbMatchedStatement.__annotations__.keys()
    for item in gkb_matches:
        stmt = copy(item)
        stmt["requiredKbMatches"].sort()
        kbs = KbMatchedStatement({key: val for (key, val) in stmt.items() if key in kbs_keys})
        dict_key = str(kbs)
        kbMatchedStatements[dict_key] = kbs
    return [*kbMatchedStatements.values()]


def get_kb_statement_matched_conditions(
    gkb_matches: List[KbMatch] | List[Hashabledict],
    allow_partial_matches: bool = False,
) -> List[KbMatchedStatementConditionSet]:
    """
    Prepares the kbMatchedStatementConditions section, with expected format
    kbStatementId: #999:999
    matchedConditions: [{'observedVariantKey': 'test1', 'kbVariantId': '#111:111'}]

    where the kbStatementId is a gkb statement rid
    and each of the observed variant keys is a reference to
    a kbMatch (ie, an observed variant/kb variant pair).

    Each record in the output from this function should represent
    one set of observed variants that satisfies the gkb variants in the
    conditions of the statement.

    If more than one set of observed variants satisfies the gkb variant conditions of the
    statement, the output from this function should include one record for each possible set.

    Eg if the stmt requires gkb variants A and B, and the observed variants include
    X which matches A, and Y and Z which both match B,
    then we expect two records for that statement in the output of this function,
    one with matchedConditions = [X, Y] and one with [X, Z].

    Expected format of one kbMatchedStatementCondition element:
    {
        "kbStatementId": "#multivariantstmt_singleconditionset",
        "matchedConditions": [
            {
                "observedVariantKey": "test1",
                "kbVariantId": "#111:111"
            },
            {
                "observedVariantKey": "tmb",
                "kbVariantId": "#333:333"
            }
        ]
    }

    Params:
        gkb_matches: KbMatch statements to be processed
        allow_partial_matches: include statements where not all requirements are satisfied
    Returns:
        list of KbStatementMatchedConditionSet records

    """

    kbMatchedStatements = get_kb_matched_statements(gkb_matches)
    kbMatchedStatementConditions = {}

    for kbStmt in kbMatchedStatements:
        stmts = [item for item in gkb_matches if item["kbStatementId"] == kbStmt["kbStatementId"]]
        requirements = {}
        for requirement in stmts[0]["requiredKbMatches"]:
            if not requirements.get(requirement, False):
                # only use explicit variant/statement links
                reqlist = [
                    {
                        "kbVariantId": requirement,
                        "observedVariantKey": item["variant"],
                    }
                    for item in gkb_matches
                    if (
                        item["kbVariantId"] == requirement
                        and item["kbStatementId"] == kbStmt["kbStatementId"]
                    )
                ]
                requirements[requirement] = reqlist

        # remove empty sets from requirements if allowing partial matches
        if allow_partial_matches:
            requirements = {key: val for (key, val) in requirements.items() if len(val) > 0}

        variantConditionSets = list(product(*requirements.values()))
        conditionSets = [
            {"kbStatementId": kbStmt["kbStatementId"], "matchedConditions": item}
            for item in variantConditionSets
        ]
        for conditionSet in conditionSets:
            matchedConditions = sorted(
                conditionSet["matchedConditions"],
                key=lambda x: (x["kbVariantId"], x["observedVariantKey"]),
            )
            kbmc = KbMatchedStatementConditionSet(
                {
                    "kbStatementId": conditionSet["kbStatementId"],
                    "matchedConditions": matchedConditions,
                }
            )
            key = str(
                uuid.uuid5(uuid.NAMESPACE_DNS, str(kbmc))
            )  # to make it more readable when debugging
            kbMatchedStatementConditions[key] = kbmc
    return [*kbMatchedStatementConditions.values()]


def get_kb_matches_sections(
    gkb_matches: List[KbMatch] | List[Hashabledict],
    allow_partial_matches=False,
) -> KbMatchSections:
    kb_variants = get_kb_variants(gkb_matches)
    kb_matched_statements = get_kb_matched_statements(gkb_matches)
    kb_statement_matched_conditions = get_kb_statement_matched_conditions(
        gkb_matches, allow_partial_matches
    )
    return {
        "kbMatches": kb_variants,
        "kbMatchedStatements": kb_matched_statements,
        "kbStatementMatchedConditions": kb_statement_matched_conditions,
    }


def get_kb_disease_matches(
    graphkb_conn: GraphKBConnection, kb_disease_match: str = None, verbose: bool = True
) -> list[str]:

    if not kb_disease_match:
        kb_disease_match = 'cancer'
        if verbose:
            logger.warning(f"No disease provided; will use '{kb_disease_match}'")

    if verbose:
        logger.info(f"Matching disease ({kb_disease_match}) to graphkb")

    disease_matches = {
        r["@rid"]
        for r in gkb_vocab.get_term_tree(
            graphkb_conn,
            kb_disease_match,
            ontology_class="Disease",
        )
    }

    if not disease_matches:
        msg = f"failed to match disease ({kb_disease_match}) to graphkb"
        if verbose:
            logger.error(msg)
        raise ValueError(msg)

    return disease_matches
