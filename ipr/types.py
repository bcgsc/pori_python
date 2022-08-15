from typing import List, Optional, Union

try:
    from typing import TypedDict  # type: ignore
except ImportError:
    from typing_extensions import TypedDict

from graphkb.types import Ontology, Record


class KbMatch(TypedDict):
    variant: str
    variantType: str
    approvedTherapy: bool
    category: str
    context: str
    kbContextId: str
    disease: str
    evidenceLevel: str
    kbStatementId: str
    kbVariant: str
    kbVariantId: str
    matchedCancer: bool
    reference: str
    relevance: str
    kbRelevanceId: str
    externalSource: str
    externalStatementId: str
    reviewStatus: str


class IprGene(TypedDict):
    name: str
    cancerRelated: Optional[bool]
    knownFusionPartner: Optional[bool]
    knownSmallMutation: Optional[bool]
    tumourSuppressor: Optional[bool]
    oncogene: Optional[bool]
    therapeuticAssociated: Optional[bool]


class IprVariantBase(TypedDict):
    """Required properties of all variants for IPR"""

    key: str
    variantType: str
    variant: str


class IprGeneVariant(IprVariantBase):
    gene: str


class IprCopyVariant(IprGeneVariant):
    # variantType == 'cnv'
    kbCategory: str
    cnvState: str


class IprExprVariant(IprGeneVariant):
    # variantType == 'exp'
    kbCategory: str
    expressionState: str
    histogramImage: Optional[str]


class IprStructVarBase(IprVariantBase):
    proteinChange: str  # TODO - make proteinChange optional and use hgvsNotations
    hgvsProtein: Optional[str]
    hgvsCds: Optional[str]
    hgvsGenomic: Optional[str]


class IprSmallMutationVariant(IprStructVarBase):
    """SNPs and small INDELs"""

    gene: str  # equivalent of gene1 in IprFusionVariant
    germline: Optional[bool]
    startPosition: Optional[int]
    endPosition: Optional[int]  # Must equal startPosition for SNPs
    normalAltCount: Optional[int]
    normalDepth: Optional[int]
    normalRefCount: Optional[int]
    rnaAltCount: Optional[int]
    rnaDepth: Optional[int]
    rnaRefCount: Optional[int]
    tumourAltCount: Optional[int]
    tumourDepth: Optional[int]
    tumourRefCount: Optional[int]


class IprFusionVariant(IprStructVarBase):
    gene1: str
    gene2: str
    exon1: int
    exon2: int
    highQuality: Optional[bool]  # high quality event found by multiple tools.
    svg: Optional[str]  # path to svg image of fusion


class ImageDefinition(TypedDict):
    key: str
    path: str


class GkbStatement(Record):
    """No 'links' handled."""

    relevance: Ontology
    subject: Ontology
    conditions: List[Ontology]
    evidence: List[Ontology]
    evidenceLevel: List[Ontology]
    source: Record
    sourceId: str
    reviewStatus: Optional[str]


IprStructuralVariant = Union[IprSmallMutationVariant, IprFusionVariant]
IprVariant = Union[IprCopyVariant, IprExprVariant, IprStructuralVariant]
