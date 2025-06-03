DEFAULT_URL = "https://iprstaging-api.bcgsc.ca/api"
GERMLINE_BASE_TERMS = ("pharmacogenomic", "cancer predisposition")  # based on graphkb.constants
VARIANT_CLASSES = {"Variant", "CategoryVariant", "PositionalVariant", "CatalogueVariant"}

# all possible values for review status are: ['pending', 'not required', 'passed', 'failed', 'initial']
FAILED_REVIEW_STATUS = "failed"

# Signatures
COSMIC_SIGNATURE_VARIANT_TYPE = "high signature"
HLA_SIGNATURE_VARIANT_TYPE = "signature present"
TMB_SIGNATURE = "mutation burden"
TMB_SIGNATURE_HIGH_THRESHOLD = (
    10.0  # genomic mutations per mb - https://www.bcgsc.ca/jira/browse/GERO-296
)
TMB_SIGNATURE_VARIANT_TYPE = "high signature"
# Mapping micro-satellite from pipeline terms to GraphKB terms
MSI_MAPPING = {
    'microsatellite instability': {  # MSI
        'displayName': 'microsatellite instability high signature',
        'signatureName': 'microsatellite instability',
        'variantTypeName': 'high signature',
    },
    'microsatellite stable': {  # MSS
        'displayName': 'microsatellite stable signature present',
        'signatureName': 'microsatellite stable',
        'variantTypeName': 'signature present',
    },
}
