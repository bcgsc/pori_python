import argparse
import json
import logging
import os
import datetime

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
from .util import logger, LOG_LEVELS
from . import ipr


def file_path(path):
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(f'{repr(path)} is not a valid filename. does not exist')
    return path


def timestamp():
    return datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')


def command_interface():
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


def main(args, optional_content=None):
    """
    Run the matching and create the report JSON for upload (TODO) to IPR

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
    alterations = []
    logger.verbose('annotating small mutations')
    annotate_positional_variants(graphkb_conn, small_mutations, args.kb_disease_match)

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
    # TODO: Append gene level information to each variant type (until IPR does this itself?)

    logger.info(f'writing: {args.output_json}')
    output = optional_content or dict()

    with open(args.output_json, 'w') as fh:

        output.update(
            {
                'kbMatches': alterations,
                'copyVariants': [c for c in copy_variants if c['gene'] in genes_with_variants],
                'smallMutations': small_mutations,
                'expressionVariants': [
                    e for e in expression_variants if e['gene'] in genes_with_variants
                ],
                'kbDiseaseMatch': args.kb_disease_match,
                'kbUrl': graphkb_conn.url,
                'kbVersion': timestamp(),
                'structuralVariants': structural_variants,
                'genes': gene_information,
            }
        )
        for section in output:
            logger.info(f'section {section} has {len(output[section])} rows')
        fh.write(json.dumps(output, indent='  ', sort_keys=True))
    logger.info(f'made {graphkb_conn.request_count} requests to graphkb')
    ipr_conn.upload_report(output)
