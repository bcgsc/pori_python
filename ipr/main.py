import argparse
import datetime
import json
import logging
import os
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from graphkb import GraphKBConnection
from typing import Dict, List

from .annotate import annotate_category_variants, annotate_positional_variants, get_gene_information
from .connection import IprConnection
from .constants import DEFAULT_URL, GERMLINE_BASE_TERMS
from .inputs import (
    check_comparators,
    check_variant_links,
    preprocess_copy_variants,
    preprocess_expression_variants,
    preprocess_small_mutations,
    preprocess_structural_variants,
    validate_report_content,
)
from .ipr import create_key_alterations, filter_structural_variants, select_expression_plots
from .summary import summarize
from .therapeutic_options import create_therapeutic_options
from .types import KbMatch
from .util import LOG_LEVELS, logger, trim_empty_values

CACHE_GENE_MINIMUM = 5000


def file_path(path: str) -> str:
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(f'{repr(path)} is not a valid filename. does not exist')
    return path


def timestamp() -> str:
    return datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')


def command_interface() -> None:
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    req = parser.add_argument_group('required arguments')
    (req if not os.environ.get('USER') else parser).add_argument(
        '--username',
        required=not os.environ.get('USER'),
        default=os.environ.get('USER'),
        help='username to use connecting to graphkb/ipr',
    )
    req.add_argument('--password', required=True, help='password to use connecting to graphkb/ipr')
    req.add_argument(
        '-c', '--content', required=True, type=file_path, help="Report Content as JSON"
    )
    parser.add_argument('--ipr_url', default=DEFAULT_URL)
    parser.add_argument('--graphkb_url', default=None)
    parser.add_argument('--log_level', default='info', choices=LOG_LEVELS.keys())
    parser.add_argument(
        '--therapeutics', default=False, help='Generate therapeutic options', action='store_true'
    )
    parser.add_argument(
        '--skip_comments',
        default=False,
        action='store_true',
        help='Turn off generating the analyst comments section of the report',
    )
    parser.add_argument(
        '-o', '--output_json_path', help='path to a JSON to output the report upload body'
    )
    parser.add_argument(
        '-w',
        '--always_write_output_json',
        action="store_true",
        help='Write to output_json_path on successful IPR uploads instead of just when the upload fails',
    )

    args = parser.parse_args()

    with open(args.content, 'r') as fh:
        content = json.load(fh)

    create_report(
        username=args.username,
        password=args.password,
        content=content,
        ipr_url=args.ipr_url,
        graphkb_url=args.graphkb_url,
        log_level=args.log_level,
        output_json_path=args.output_json_path,
        always_write_output_json=args.always_write_output_json,
        generate_therapeutics=args.therapeutics,
        generate_comments=not args.skip_comments,
    )


def clean_unsupported_content(upload_content: Dict) -> Dict:
    """
    Remove unsupported content. This content is either added to facilitate creation
    or to support upcoming and soon to be supported content that we would like
    to implement but is not yet supported by the upload
    """
    drop_columns = ['variant', 'variantType', 'histogramImage']
    for variant_section in [
        'expressionVariants',
        'smallMutations',
        'copyVariants',
        'structuralVariants',
    ]:
        for variant in upload_content[variant_section]:
            for col in drop_columns:
                if col in variant:
                    del variant[col]

    for row in upload_content['kbMatches']:
        del row['kbContextId']
        del row['kbRelevanceId']
    return upload_content


def create_report(
    username: str,
    password: str,
    content: Dict,
    ipr_url: str = DEFAULT_URL,
    log_level: str = 'info',
    output_json_path: str = None,
    always_write_output_json: bool = False,
    ipr_upload: bool = True,
    interactive: bool = False,
    graphkb_url: str = '',
    generate_therapeutics: bool = False,
    generate_comments: bool = True,
    match_germline: bool = True,
    kb_match_filters: List[Dict] = [],
) -> Dict:
    """
    Run the matching and create the report JSON for upload to IPR

    Args:
        username: the username for connecting to GraphKB and IPR
        password: the password for connecting to GraphKB and IPR
        ipr_url: base URL to use in connecting to IPR
        log_level: the logging level
        content: report content
        output_json_path: path to a JSON file to output the report upload body.
        always_write_output_json: with successful IPR upload
        ipr_upload: upload report to ipr
        interactive: progressbars for interactive users
        cache_gene_minimum: minimum number of genes required for gene name caching optimization
        generate_therapeutics: create therapeutic options for upload with the report
        generate_comments: create the analyst comments section for upload with the report
        match_germline: check for germline status for matching
        kb_match_filters: filtered out properties [{alteration_field:(bool=require/exclude, values)}]
                          eg. Exclude cancer predisposition matches if 'externalSource' not 'CGL'
                          [{'category': (True, ['cancer predisposition']), 'externalSource': (False, ['CGL'])}]

    Returns:
        ipr_conn.upload_report return dictionary
    """
    # set the default logging configuration
    logging.basicConfig(
        level=LOG_LEVELS[log_level],
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
        datefmt='%m-%d-%y %H:%M:%S',
    )
    # validate the JSON content follows the specification
    validate_report_content(content)
    kb_disease_match = content['kbDiseaseMatch']
    # validate the input variants
    small_mutations = preprocess_small_mutations(content.get('smallMutations', []))
    copy_variants = preprocess_copy_variants(content.get('copyVariants', []))
    structural_variants = preprocess_structural_variants(content.get('structuralVariants', []))
    expression_variants = preprocess_expression_variants(content.get('expressionVariants', []))
    check_comparators(content, expression_variants)

    ipr_conn = IprConnection(username, password, ipr_url)
    if graphkb_url:
        logger.info(f'connecting to graphkb: {graphkb_url}')
        graphkb_conn = GraphKBConnection(graphkb_url)
    else:
        graphkb_conn = GraphKBConnection()
    graphkb_conn.login(username, password)

    genes_with_variants = check_variant_links(
        small_mutations, expression_variants, copy_variants, structural_variants
    )

    # filter excess variants not required for extra gene information
    logger.info(f'annotating small mutations')
    alterations: List[KbMatch] = annotate_positional_variants(
        graphkb_conn, small_mutations, kb_disease_match, show_progress=interactive
    )

    logger.info(f'annotating structural variants')
    alterations.extend(
        annotate_positional_variants(
            graphkb_conn, structural_variants, kb_disease_match, show_progress=interactive
        )
    )

    logger.info(f'annotating copy variants')
    alterations.extend(
        annotate_category_variants(
            graphkb_conn, copy_variants, kb_disease_match, show_progress=interactive
        )
    )

    logger.info(f'annotating expression variants')
    alterations.extend(
        annotate_category_variants(
            graphkb_conn,
            expression_variants,
            kb_disease_match,
            copy_variant=False,
            show_progress=interactive,
        )
    )

    all_variants = expression_variants + copy_variants + structural_variants + small_mutations

    if match_germline:  # verify germline kb statements matched germline observed variants
        germ_alts = [alt for alt in alterations if alt['category'] in GERMLINE_BASE_TERMS]
        if germ_alts:
            logger.info(f"checking germline status of {GERMLINE_BASE_TERMS}")
            alterations = [alt for alt in alterations if alt not in germ_alts]
            for alt in germ_alts:
                var_list = [v for v in all_variants if v['key'] == alt['variant']]
                germline_var_list = [v for v in var_list if 'germline' in v and v['germline']]
                if germline_var_list:
                    logger.info(
                        f"germline kbStatementId:{alt['kbStatementId']}: {alt['kbVariant']} {alt['category']}"
                    )
                    alterations.append(alt)
                elif var_list:
                    logger.info(
                        f"Dropping somatic match to kbStatementId:{alt['kbStatementId']}: {alt['kbVariant']} {alt['category']}"
                    )
                else:
                    logger.error(
                        f"germline check fail for: {alt['kbStatementId']}: {alt['kbVariant']} {alt['category']}"
                    )
                    alterations.append(alt)

    if kb_match_filters:
        logger.info('filtering kbmatched alterations.')
        filtered_alts = []
        for alt in alterations:
            for filter_dict in kb_match_filters:
                filtered_pos = all(
                    [alt[field] in values for field, (comp, values) in filter_dict.items() if comp]
                )
                filtered_exclude = all(
                    [
                        alt[field] not in values
                        for field, (comp, values) in filter_dict.items()
                        if not comp
                    ]
                )
                if filtered_pos and filtered_exclude:
                    logger.info(
                        f"Dropping kbStatementId:{alt['kbStatementId']}: {alt['kbVariant']} - {filter_dict}"
                    )
                    filtered_alts.append(alt)
        alterations = [alt for alt in alterations if alt not in filtered_alts]
    key_alterations, variant_counts = create_key_alterations(alterations, all_variants)

    logger.info('fetching gene annotations')
    gene_information = get_gene_information(graphkb_conn, genes_with_variants)

    if generate_therapeutics:
        logger.info('generating therapeutic options')
        targets = create_therapeutic_options(graphkb_conn, alterations, all_variants)
    else:
        targets = []

    logger.info('generating analyst comments')
    if generate_comments:
        comments = {
            'comments': summarize(
                graphkb_conn, alterations, disease_name=kb_disease_match, variants=all_variants
            )
        }
    else:
        comments = {'comments': ''}

    # thread safe deep-copy the original content
    output = json.loads(json.dumps(content))
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
            'kbDiseaseMatch': kb_disease_match,
            'kbUrl': graphkb_conn.url,
            'kbVersion': timestamp(),
            'structuralVariants': [
                trim_empty_values(s)
                for s in filter_structural_variants(
                    structural_variants, alterations, gene_information
                )
            ],
            'genes': gene_information,
            'genomicAlterationsIdentified': key_alterations,
            'variantCounts': variant_counts,
            'analystComments': comments,
            'therapeuticTarget': targets,
        }
    )
    output.setdefault('images', []).extend(
        select_expression_plots(
            alterations, expression_variants + copy_variants + structural_variants + small_mutations
        )
    )

    output = clean_unsupported_content(output)
    ipr_result = None
    upload_error = None

    if ipr_upload:
        try:
            logger.info(f'Uploading to IPR {ipr_conn.url}')
            ipr_result = ipr_conn.upload_report(output)
            logger.info(ipr_result)
            output.update(ipr_result)
        except Exception as err:
            upload_error = err
            logger.error(f"ipr_conn.upload_report failed: {err}", exc_info=True)
    if output_json_path:
        if always_write_output_json or not ipr_result:
            logger.info(f'Writing IPR upload json to: {output_json_path}')
            with open(output_json_path, 'w') as fh:
                fh.write(json.dumps(output))
    logger.info(f'made {graphkb_conn.request_count} requests to graphkb')
    logger.info(f'average load {int(graphkb_conn.load or 0)} req/s')
    if upload_error:
        raise upload_error
    return output
