"""
upload variant and report information to IPR
"""

from __future__ import annotations

import pandas
from typing import Any, Dict, List, Sequence

from pori_python.graphkb import GraphKBConnection
from pori_python.types import Hashabledict, IprVariant, KbMatch

from .util import (
    create_variant_name_tuple,
    find_variant,
    get_preferred_drug_representation,
    get_terms_set,
)


def create_therapeutic_options(
    graphkb_conn: GraphKBConnection,
    kb_matches: List[KbMatch] | List[Hashabledict],
    variants: Sequence[IprVariant],
) -> List[Dict]:
    """
    Generate therapeutic options summary from the list of kb-matches
    """
    options: List[Dict[str, Any]] = []
    resistance_markers = get_terms_set(graphkb_conn, ["no sensitivity"])

    for match in kb_matches:
        row_type = "therapeutic"
        for stmt in match["kbMatchedStatements"]:
            if stmt["category"] != "therapeutic" or stmt["relevance"] == "eligibility":
                continue

            if stmt["kbRelevanceId"] in resistance_markers:
                row_type = "chemoresistance"

            # do this inside the loop because we expect that most matches will only have
            # one statement - so in general will only be run once per outer loop anyway -
            # and that running it in the outer loop will result in more executions anyway
            # because the category/relevance check will not have been done yet
            variant = find_variant(variants, match["variantType"], match["variant"])

            # could possibly extract all drugs and run this outside loop to avoid
            # adding to runtime but not sure how long it takes or how many dupes
            # there are likely to be
            drug = get_preferred_drug_representation(graphkb_conn, stmt["kbContextId"])

            gene, variant_string = create_variant_name_tuple(variant)

            # TODO this may need updating when we update the ptt table
            options.append(
                {
                    "gene": gene,
                    "type": row_type,
                    "therapy": drug["displayName"],
                    "therapyGraphkbId": drug["@rid"],
                    "context": stmt["relevance"],
                    "contextGraphkbId": stmt["kbRelevanceId"],
                    "variantGraphkbId": match["kbVariantId"],
                    "variant": variant_string,
                    "evidenceLevel": stmt["evidenceLevel"],
                    "kbStatementIds": stmt["kbStatementId"],
                    "notes": "",
                }
            )
    if not options:
        return options
    options_df = pandas.DataFrame.from_records(options)

    def delimited_list(inputs: List, delimiter: str = " / ") -> str:
        return delimiter.join(sorted(list({i for i in inputs if i})))

    options_df = options_df.groupby(["gene", "type", "therapy", "variant"]).agg(
        {
            "evidenceLevel": delimited_list,
            "context": delimited_list,
            "notes": lambda x: delimited_list(x, " "),
        }
    )
    options_df = options_df.reset_index()
    options = options_df.to_dict("records")  # type: ignore
    therapeutic_rank = 0
    chemoresistance_rank = 0
    for option in options:
        if option["type"] == "therapeutic":
            option["rank"] = therapeutic_rank
            therapeutic_rank += 1
        else:
            option["rank"] = chemoresistance_rank
            chemoresistance_rank += 1
    return options
