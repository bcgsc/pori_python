from typing import Dict, List, Optional, Tuple, TypedDict, Union

# TODO: Can constants in inputs.py like COPY_REQ, SMALL_MUT_REQ, just be replaced by types?

CategoryBaseTermMapping = List[Tuple[str, List[str]]]
Record = TypedDict("Record", {"@rid": str, "@class": str})
EmbeddedRecord = TypedDict("EmbeddedRecord", {"@class": str})


class Ontology(Record):
    sourceId: str
    sourceIdVersion: Optional[str]
    name: str
    source: Record
    displayName: str


class BasicPosition(EmbeddedRecord):
    pos: int


class CytobandPosition(EmbeddedRecord):
    arm: str
    majorBand: str
    minorBand: str


Position = Union[BasicPosition, CytobandPosition]


class Variant(Record):
    reference1: Ontology
    reference2: Optional[Ontology]
    type: Ontology
    zygosity: str
    germline: bool
    displayName: str


class PositionalVariant(Variant):
    break1Start: Union[Position, CytobandPosition]
    break1End: Optional[Union[Position, CytobandPosition]]
    break2Start: Optional[Union[Position, CytobandPosition]]
    break2End: Optional[Union[Position, CytobandPosition]]
    refSeq: Optional[str]
    untemplatedSeq: Optional[str]
    untemplatedSeqSize: Optional[int]


class ParsedVariant(TypedDict):
    reference1: str
    reference2: Optional[str]
    type: str
    zygosity: str
    germline: bool
    break1Start: Union[Position, CytobandPosition]
    break1End: Optional[Union[Position, CytobandPosition]]
    break2Start: Optional[Union[Position, CytobandPosition]]
    break2End: Optional[Union[Position, CytobandPosition]]
    refSeq: Optional[str]
    untemplatedSeq: Optional[str]
    untemplatedSeqSize: Optional[int]


class Statement(Record):
    relevance: Ontology
    subject: Ontology
    conditions: List[Ontology]
    evidence: List[Ontology]
    evidenceLevel: List[Ontology]
    source: Record
    sourceId: str


class KbMatch(TypedDict):
    variant: str
    variantType: str
    approvedTherapy: bool
    category: str
    context: str
    kbContextId: str
    disease: str
    evidenceLevel: str
    iprEvidenceLevel: Optional[str]
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
    kbData: Optional[Dict]


class IprGene(TypedDict):
    name: str
    kbStatementRelated: Optional[bool]
    knownFusionPartner: Optional[bool]
    knownSmallMutation: Optional[bool]
    tumourSuppressor: Optional[bool]
    oncogene: Optional[bool]
    therapeuticAssociated: Optional[bool]
    cancerGeneListMatch: Optional[bool]


class IprVariantBase(TypedDict):
    """Required properties of all variants for IPR."""

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
    """One of the hgvs notations or proteinChange is required."""

    hgvsProtein: Optional[str]
    hgvsCds: Optional[str]
    hgvsGenomic: Optional[str]
    proteinChange: Optional[str]  # Older - being deprecated


class IprSmallMutationVariant(IprStructVarBase):
    """SNPs and small INDELs"""

    # variantType == 'mut'
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
    # variantType = 'sv
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
    displayNameTemplate: str


IprStructuralVariant = Union[IprSmallMutationVariant, IprFusionVariant]
IprVariant = Union[IprCopyVariant, IprExprVariant, IprStructuralVariant]
