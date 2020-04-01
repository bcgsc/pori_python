import argparse
import json
import logging
import os
import datetime
from typing import Dict, Optional

from argparse_env import ArgumentParser, Action
from graphkb import GraphKBConnection

from .inputs import (
    load_copy_variants,
    load_small_mutations,
    load_expression_variants,
    load_structural_variants,
    check_variant_links,
)
from .annotate import annotate_category_variants, annotate_positional_variants, get_gene_information
from .util import logger, LOG_LEVELS, trim_empty_values
from . import ipr


def file_path(path: str) -> str:
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(f'{repr(path)} is not a valid filename. does not exist')
    return path


def timestamp() -> str:
    return datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')


def command_interface() -> None:
    parser = ArgumentParser()
    parser.add_argument(
        '--username',
        env=True,
        action=Action,
        required=True,
        help='username to use connecting to graphkb/ipr',
    )
    parser.add_argument(
        '--password',
        env=True,
        action=Action,
        required=True,
        sensitive=True,
        help='password to use connecting to graphkb/ipr',
    )
    parser.add_argument('-c', '--copy_variants', required=False, type=file_path)
    parser.add_argument('-m', '--small_mutations', required=False, type=file_path)
    parser.add_argument('-s', '--structural_variants', required=False, type=file_path)
    parser.add_argument('-e', '--expression_variants', required=False, type=file_path)
    parser.add_argument(
        '-d',
        '--kb_disease_match',
        required=True,
        help='Disease name to be used in matching to GraphKB',
    )
    parser.add_argument('--ipr_url', default=ipr.DEFAULT_URL)
    parser.add_argument('--log_level', default='info', choices=LOG_LEVELS.keys())

    # TODO: upload JSON to IPR instead of writing output
    parser.add_argument(
        '-o',
        '--output_json',
        default='ipr_input.json',
        help='file path to write the output json content to',
    )

    args = parser.parse_args()

    main(args)


def clean_variant_rows(variants):
    # IPR cannot take these as input, may add variant later but will always drop type
    for v in variants:
        del v['variant']
        del v['variantType']
    return variants


def main(args, optional_content=None):
    """
    Run the matching and create the report JSON for upload to IPR

    Args:
        args (argparse.Namespace): Namespace of arguments with file names for variant inputs etc
        optional_content (dict): Pass-through content to include in the JSON upload
    """
    # set the default logging configuration
    logging.basicConfig(
        level=LOG_LEVELS[args.log_level],
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
        datefmt='%m-%d-%y %H:%M:%S',
    )
    ipr_conn = ipr.IprConnection(args.username, args.password, args.ipr_url)
    graphkb_conn = GraphKBConnection()
    graphkb_conn.login(args.username, args.password)

    copy_variants = load_copy_variants(args.copy_variants) if args.copy_variants else []
    logger.info(f'loaded {len(copy_variants)} copy variants')

    small_mutations = load_small_mutations(args.small_mutations) if args.small_mutations else []
    logger.info(f'loaded {len(small_mutations)} small mutations')

    expression_variants = (
        load_expression_variants(args.expression_variants) if args.expression_variants else []
    )
    logger.info(f'loaded {len(expression_variants)} expression variants')

    structural_variants = (
        load_structural_variants(args.structural_variants) if args.structural_variants else []
    )
    logger.info(f'loaded {len(structural_variants)} structural variants')

    genes_with_variants = check_variant_links(
        small_mutations, expression_variants, copy_variants, structural_variants
    )

    # filter excess variants not required for extra gene information
    logger.verbose('annotating small mutations')
    alterations = annotate_positional_variants(graphkb_conn, small_mutations, args.kb_disease_match)

    logger.verbose('annotating structural variants')
    alterations.extend(
        annotate_positional_variants(graphkb_conn, structural_variants, args.kb_disease_match)
    )

    logger.verbose('annotating copy variants')
    alterations.extend(
        annotate_category_variants(graphkb_conn, copy_variants, args.kb_disease_match)
    )

    logger.verbose('annotating expression variants')
    alterations.extend(
        annotate_category_variants(graphkb_conn, expression_variants, args.kb_disease_match, False)
    )
    logger.verbose('fetching gene annotations')
    gene_information = get_gene_information(graphkb_conn, genes_with_variants)

    logger.info(f'writing: {args.output_json}')
    output = optional_content or dict()

    key_alterations, variant_counts = ipr.create_key_alterations(
        alterations, expression_variants, copy_variants, structural_variants, small_mutations
    )

    with open(args.output_json, 'w') as fh:

        output.update(
            {
                'kbMatches': [trim_empty_values(a) for a in alterations],
                'copyVariants': [
                    trim_empty_values(c) for c in copy_variants if c['gene'] in genes_with_variants
                ],
                'smallMutations': [trim_empty_values(s) for s in small_mutations],
                'expressionVariants': [
                    trim_empty_values(e)
                    for e in expression_variants
                    if e['gene'] in genes_with_variants
                ],
                'kbDiseaseMatch': args.kb_disease_match,
                'kbUrl': graphkb_conn.url,
                'kbVersion': timestamp(),
                'structuralVariants': [trim_empty_values(s) for s in structural_variants],
                'genes': gene_information,
                'genomicAlterationsIdentified': key_alterations,
                'variantCounts': variant_counts,
            }
        )
        for section in output:
            logger.info(f'section {section} has {len(output[section])} rows')
        fh.write(json.dumps(output, indent='  ', sort_keys=True))
    logger.info(f'made {graphkb_conn.request_count} requests to graphkb')
    ipr_conn.upload_report(output)
