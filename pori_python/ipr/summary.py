from __future__ import annotations

import base64
import json
from typing import Callable, Dict, List, Sequence, Set, Tuple, cast
from urllib.parse import urlencode

from pori_python.graphkb import GraphKBConnection
from pori_python.graphkb.constants import RELEVANCE_BASE_TERMS
from pori_python.graphkb.statement import categorize_relevance
from pori_python.graphkb.util import convert_to_rid_list
from pori_python.ipr.inputs import create_graphkb_sv_notation
from pori_python.ipr.connection import IprConnection
from pori_python.types import (
    Hashabledict,
    IprVariant,
    KbMatch,
    Ontology,
    Record,
    Statement,
)

from .util import (
    convert_to_rid_set,
    generate_ontology_preference_key,
    get_preferred_drug_representation,
    get_preferred_gene_name,
    logger,
)

OTHER_DISEASES = "other disease types"
ENTREZ_GENE_URL = "https://www.ncbi.nlm.nih.gov/gene"
# TODO: https://www.bcgsc.ca/jira/browse/DEVSU-1181
GRAPHKB_GUI = "https://graphkb.bcgsc.ca"


def filter_by_record_class(
    record_list: Sequence[Record], *record_classes, exclude: bool = False
) -> List[Record]:
    """Given a list of records, return the subset matching a class or list of classes."""

    def check(name: str) -> bool:
        if exclude:
            return name not in record_classes
        else:
            return name in record_classes

    return [rec for rec in record_list if check(rec["@class"])]


def natural_join(word_list: List[str]) -> str:
    if len(word_list) > 1:
        return ", ".join(word_list[:-1]) + ", and " + word_list[-1]
    return "".join(word_list)


def get_displayname(rec: Record) -> str:
    ret_val = rec.get("displayName", rec["@rid"])
    return str(ret_val)


def natural_join_records(
    records: Sequence[Record], covert_to_word: Callable[[Record], str] = get_displayname
) -> str:
    word_list = sorted(list({covert_to_word(rec) for rec in records}))
    return natural_join(word_list)


def create_graphkb_link(record_ids: List[str], record_class: str = "Statement") -> str:
    """
    Create a link for a set of statements to the GraphKB client
    """
    record_ids = sorted(list(set(record_ids)))
    if len(record_ids) == 1:
        return f'{GRAPHKB_GUI}/view/{record_class}/{record_ids[0].replace("#", "")}'
    complex_param = base64.b64encode(json.dumps({"target": record_ids}).encode("utf-8"))
    search_params = {"complex": complex_param, "@class": record_class}
    return f"{GRAPHKB_GUI}/data/table?{urlencode(search_params)}"


def merge_diseases(
    diseases: List[Ontology] | List[Record], disease_matches: Set[str] = set()
) -> str:
    if len(convert_to_rid_set(diseases) - disease_matches) >= 2 and all(
        [d["@class"] == "Disease" for d in diseases]
    ):
        words = sorted(
            list(set([get_displayname(s) for s in diseases if s["@rid"] in disease_matches]))
        )
        words.append(OTHER_DISEASES)
        return natural_join(words)
    else:
        return natural_join_records(diseases)


def substitute_sentence_template(
    template: str,
    conditions: List[Ontology],
    subjects: List[Ontology],
    relevance: Ontology,
    evidence: List[Ontology],
    statement_rids: List[str] = [],
    disease_matches: Set[str] = set(),
) -> str:
    """Create the filled-in sentence template for a given template and list of substitutions
    which may be the result of the aggregation of 1 or more statements.
    """
    disease_conditions = filter_by_record_class(conditions, "Disease")
    variant_conditions = filter_by_record_class(
        conditions, "CategoryVariant", "CatalogueVariant", "PositionalVariant"
    )
    other_conditions = filter_by_record_class(
        conditions,
        "CategoryVariant",
        "CatalogueVariant",
        "PositionalVariant",
        "Disease",
        exclude=True,
    )
    result = template.replace(r"{relevance}", relevance["displayName"])

    if r"{subject}" in template:
        # remove subject from the conditions replacements
        subjects_ids = convert_to_rid_set(subjects)
        disease_conditions = [
            cast(Ontology, d) for d in disease_conditions if d["@rid"] not in subjects_ids
        ]
        variant_conditions = [
            cast(Ontology, d) for d in variant_conditions if d["@rid"] not in subjects_ids
        ]
        other_conditions = [d for d in other_conditions if d["@rid"] not in subjects_ids]

        result = result.replace(r"{subject}", merge_diseases(subjects, disease_matches))

    if r"{conditions:disease}" in template:
        result = result.replace(
            r"{conditions:disease}", merge_diseases(disease_conditions, disease_matches)
        )
    else:
        other_conditions.extend(disease_conditions)

    if r"{conditions:variant}" in template:
        result = result.replace(r"{conditions:variant}", natural_join_records(variant_conditions))
    else:
        other_conditions.extend(variant_conditions)

    result = result.replace(r"{conditions}", natural_join_records(other_conditions))

    link_url = create_graphkb_link(statement_rids) if statement_rids else ""

    if r"{evidence}" in template:
        evidence_str = ", ".join(sorted(list({e["displayName"] for e in evidence})))
        if link_url:
            evidence_str = f'<a href="{link_url}" target="_blank" rel="noopener">{evidence_str}</a>'
        result = result.replace(r"{evidence}", evidence_str)

    return result


def aggregate_statements(
    graphkb_conn: GraphKBConnection,
    template: str,
    statements: List[Statement],
    disease_matches: Set[str],
) -> Dict[str, str]:
    """
    Group Statements that only differ in disease conditions and evidence
    """
    hash_other: Dict[Tuple, List[Statement]] = {}

    def generate_key(statement: Statement) -> Tuple:
        result = [
            cond.get("displayName", cond["@rid"])
            for cond in filter_by_record_class(statement["conditions"], "Disease", exclude=True)
            if cond["@rid"] != statement["subject"]["@rid"]
        ]
        if statement.get("subject", {}).get("@class", "Disease") != "Disease":
            subject = statement["subject"]
            if subject["@class"] == "Therapy":
                alt = get_preferred_drug_representation(graphkb_conn, subject["@rid"])
                statement["subject"] = alt
            result.append(statement["subject"]["displayName"])
        result.append(statement["relevance"]["displayName"])
        result.append(statement["displayNameTemplate"])
        return tuple(sorted(set(result)))

    for statement in statements:
        key = generate_key(statement)
        hash_other.setdefault(key, []).append(statement)

    result = {}
    for key, group in hash_other.items():
        conditions = []
        subjects = []
        evidence = []
        relevance = group[0]["relevance"]
        template = group[0]["displayNameTemplate"]
        for statement in group:
            conditions.extend(statement["conditions"])
            evidence.extend(statement["evidence"])
            subjects.append(statement["subject"])

        sentence = substitute_sentence_template(
            template,
            conditions,
            subjects,
            relevance,
            evidence,
            statement_rids=convert_to_rid_list(group),
            disease_matches=disease_matches,
        )

        for statement in group:
            result[statement["@rid"]] = sentence
    return result


def display_variant(variant: IprVariant) -> str:
    """Short, human readable variant description string."""
    gene = variant.get("gene", "")
    if not gene and "gene1" in variant and "gene2" in variant:
        gene = f'({variant.get("gene1", "")},{variant.get("gene2", "")})'

    if variant.get("kbCategory"):
        return f'{variant.get("kbCategory")} of {gene}'

    # Special display of IprFusionVariant with exons
    if variant.get("exon1") or variant.get("exon2"):
        return create_graphkb_sv_notation(variant)  # type: ignore

    # Use chosen legacy 'proteinChange' or an hgvs description of lowest detail.
    hgvs = variant.get(
        "proteinChange",
        variant.get("hgvsProtein", variant.get("hgvsCds", variant.get("hgvsGenomic", ""))),
    )

    if gene and hgvs:
        return f"{gene}:{hgvs}"
    elif variant.get("variant"):
        return str(variant.get("variant"))

    raise ValueError(f"Unable to form display_variant of {variant}")


def display_variants(gene_name: str, variants: List[IprVariant]) -> str:
    result = sorted(list({v for v in [display_variant(e) for e in variants] if gene_name in v}))
    variants_text = natural_join(result)
    if len(result) > 1:
        return (
            f"Multiple variants of the gene {gene_name} were observed in this case: {variants_text}"
        )
    elif result:
        return f"{variants_text[0].upper()}{variants_text[1:]} was observed in this case."
    return ""


def create_section_html(
    graphkb_conn: GraphKBConnection,
    gene_name: str,
    sentences_by_statement_id: Dict[str, str],
    statements: Dict[str, Statement],
    exp_variants: List[IprVariant],
) -> str:
    """
    Generate HTML for a gene section of the comments
    """
    output = [f"<h2>{gene_name}</h2>"]

    sentence_categories: Dict[str, str] = {}

    for statement_id, sentence in sentences_by_statement_id.items():
        relevance = statements[statement_id]["relevance"]["@rid"]
        category = categorize_relevance(
            graphkb_conn,
            relevance,
            RELEVANCE_BASE_TERMS + [("resistance", ["no sensitivity"])],
        )
        sentence_categories[sentence] = category

    # get the entrez gene descriptive hugo name
    genes = graphkb_conn.query(
        {
            "target": "Feature",
            "filters": {
                "AND": [
                    {
                        "source": {
                            "target": "Source",
                            "filters": {"name": "entrez gene"},
                        }
                    },
                    {"name": gene_name},
                    {"biotype": "gene"},
                ]
            },
        }
    )
    genes = sorted(genes, key=generate_ontology_preference_key)  # type: ignore

    variants_text = display_variants(gene_name, exp_variants)
    if not variants_text:
        # exclude sections where they are not linked to an experimental variant. this can occur when there are co-occurent statements collected
        return ""
    if genes and genes[0].get("description", ""):
        description = ". ".join(genes[0]["description"].split(". ")[:2])  # type: ignore
        sourceId = genes[0].get("sourceId", "")

        output.append(
            f"""
<blockquote class="entrez_description" cite="{ENTREZ_GENE_URL}/{sourceId}">
    {description}.
</blockquote>
<p>
    {variants_text}
</p>
"""
        )

    sentences_used: Set[str] = set()

    for section in [
        {s for (s, v) in sentence_categories.items() if v == "diagnostic"},
        {s for (s, v) in sentence_categories.items() if v == "biological"},
        {s for (s, v) in sentence_categories.items() if v in ["therapeutic", "prognostic"]},
        {
            s
            for (s, v) in sentence_categories.items()
            if v
            not in [
                "diagnostic",
                "biological",
                "therapeutic",
                "prognostic",
                "resistance",
            ]
        },
        {s for (s, v) in sentence_categories.items() if v == "resistance"},
    ]:
        content = ". ".join(sorted(list(section - sentences_used)))
        sentences_used.update(section)
        output.append(f"<p>{content}</p>")
    return "\n".join(output)


def section_statements_by_genes(
    graphkb_conn: GraphKBConnection, statements: Sequence[Statement]
) -> Dict[str, Set[str]]:
    """Create Dict of statement @rid sets indexed by preferred gene names in conditions."""
    genes: Dict[str, Set[str]] = {}

    for statement in statements:
        for condition in statement["conditions"]:
            if condition.get("biotype", "") == "gene":
                gene = get_preferred_gene_name(graphkb_conn, condition["@rid"])
                genes.setdefault(gene, set()).add(statement["@rid"])
            else:
                for cond_ref_key in ("reference1", "reference2"):
                    cond_ref_gene = condition.get(cond_ref_key)
                    if cond_ref_gene:
                        gene = get_preferred_gene_name(graphkb_conn, str(cond_ref_gene))
                        genes.setdefault(gene, set()).add(statement["@rid"])

    return genes


def prep_single_ipr_variant_comment(variant_text):
    """Formats single item of custom variant text for inclusion in the analyst comments.

    Params:
        variant_text:

    Returns:
        section: html-formatted string
    """
    cancer_type = ",".join(variant_text["cancerType"])
    if not cancer_type:
        cancer_type = "no specific cancer types"
    cancer_type = f" ({cancer_type})"
    section = [f"<h2>{variant_text['variantName']}{cancer_type}</h2>"]
    section.append(f"<p>{variant_text['text']}</p>")
    return section


def get_ipr_analyst_comments(
    ipr_conn: IprConnection,
    matches: Sequence[KbMatch] | Sequence[Hashabledict],
    disease_name: str,
    project_name: str,
    report_type: str,
    include_nonspecific_disease: bool = False,
    include_nonspecific_project: bool = False,
    include_nonspecific_template: bool = False,
) -> str:
    """
    Given a list of kbmatches, checks the variant_texts table in IPR-API to get any
    pre-prepared text for this variant for inclusion in the analyst comments.
    Matches on template, project and variant_name. Matches on project, disease and template
    if possible. If no match is found and the related include_nonspecific arg is True,
    uses a result with no specified value for that field if a result is found (eg
    a result with no cancer type specified, if it exists).

    Params:
        ipr_conn: connection to the ipr db
        matches: list of kbmatches which will be included in the report
        disease_name: str, eg 'colorectal cancer'
        project_name: str, eg TEST or pog
        report_type: str, eg genomic or rapid
        include_nonspecific_disease: bool - true if variant texts that don't explicitly
            name a cancer type should be included
        include_nonspecific_project: bool - true if variant texts that don't explicitly
            name a project should be included
        include_nonspecific_template: bool - true if variant texts that don't explicitly
            name a project should be included
    Returns:
        html-formatted string
    """
    output_header = "<h3>The comments below were automatically drawn from curated text stored in IPR for variant matches in this report, and have not been manually reviewed</h3>"
    no_comments_found_output = "No comments found in IPR for variants in this report"
    output = []
    # get the list of variants to check for custom text for
    match_set = list(set([item["kbVariant"] for item in matches]))

    for variant in match_set:
        data = {
            "variantName": variant,
        }
        itemlist: list[dict] = []
        itemlist = ipr_conn.get("variant-text", data=data)  # type: ignore
        if itemlist:
            project_matches = [
                item
                for item in itemlist
                if 'project' in item.keys() and item['project']['name'] == project_name
            ]
            if project_matches:
                itemlist = project_matches
            elif include_nonspecific_project:
                itemlist = [item for item in itemlist if 'project' not in item.keys()]
            else:
                itemlist = []

            template_matches = [
                item
                for item in itemlist
                if 'template' in item.keys() and item['template']['name'] == report_type
            ]
            if template_matches:
                itemlist = template_matches
            elif include_nonspecific_template:
                itemlist = [item for item in itemlist if 'template' not in item.keys()]
            else:
                itemlist = []

            disease_matches = [item for item in itemlist if disease_name in item['cancerType']]
            if disease_matches:
                itemlist = disease_matches
            elif include_nonspecific_disease:
                itemlist = [item for item in itemlist if not item['cancerType']]
            else:
                itemlist = []

            for item in itemlist:
                section = prep_single_ipr_variant_comment(item)
                output.extend(section)

    if not output:
        return no_comments_found_output
    output.insert(0, output_header)
    return "\n".join(output)


def auto_analyst_comments(
    graphkb_conn: GraphKBConnection,
    matches: Sequence[KbMatch] | Sequence[Hashabledict],
    disease_matches: set[str],
    variants: Sequence[IprVariant],
) -> str:
    """Given a list of GraphKB matches, generate a text summary to add to the report."""
    templates: Dict[str, List[Statement]] = {}
    statements: Dict[str, Statement] = {}
    variants_by_keys = {v["key"]: v for v in variants}
    variant_keys_by_statement_ids: Dict[str, Set[str]] = {}

    for match in matches:
        rid = match["kbStatementId"]
        exp_variant = match["variant"]
        variant_keys_by_statement_ids.setdefault(rid, set()).add(exp_variant)

    exp_variants_by_statements: Dict[str, List[IprVariant]] = {}
    for rid, keys in variant_keys_by_statement_ids.items():
        try:
            exp_variants_by_statements[rid] = [variants_by_keys[key] for key in keys]
        except KeyError as err:
            logger.warning(f"No specific variant matched for {rid}:{keys} - {err}")
            exp_variants_by_statements[rid] = []

    # get details for statements
    for match in matches:
        rid = match["kbStatementId"].replace("#", "")
        result = graphkb_conn.request(f"/statements/{rid}?neighbors=1")["result"]

        templates.setdefault(result["displayNameTemplate"], []).append(result)
        statements[result["@rid"]] = result

    # aggregate similar sentences
    sentences = {}
    for template, group in templates.items():
        sentences.update(aggregate_statements(graphkb_conn, template, group, disease_matches))

    # section statements by genes
    statements_by_genes = section_statements_by_genes(graphkb_conn, list(statements.values()))

    output: List[str] = [
        "<h3>The comments below were automatically generated from matches to GraphKB and have not been manually reviewed</h3>"
    ]

    for section, statement_rids in sorted(
        statements_by_genes.items(), key=lambda x: len(x[1]), reverse=True
    ):
        exp_variants = {}
        for variant_list in [exp_variants_by_statements[r] for r in statement_rids]:
            for variant in variant_list:
                exp_variants[variant["key"]] = variant

        output.append(
            create_section_html(
                graphkb_conn,
                section,
                {r: sentences[r] for r in statement_rids},
                {r: statements[r] for r in statement_rids},
                list(exp_variants.values()),
            )
        )

    return "\n".join(output)
