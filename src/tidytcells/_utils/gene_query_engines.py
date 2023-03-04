"""
Gene query engines
"""


from abc import ABC, abstractmethod
from .._resources import *
from typing import FrozenSet
from warnings import warn


class GeneQueryEngine(ABC):
    """
    Abstract base class for gene query engines.
    """

    @staticmethod
    @abstractmethod
    def query(precision: str) -> FrozenSet[str]:
        """
        List all unique genes/alleles of a particular type and species to the
        level of precision specified.
        """


class TCRQueryEngine(GeneQueryEngine):
    ref_dict = dict()

    @classmethod
    def query(cls, precision: str, functionality: str) -> FrozenSet[str]:
        tcrs = []

        for gene in cls.ref_dict:
            if precision == "gene":
                if (
                    functionality == "any"
                    or (
                        functionality in ("F", "P", "ORF")
                        and functionality in cls.ref_dict[gene].values()
                    )
                    or (
                        functionality == "NF"
                        and {"P", "ORF"}.intersection(cls.ref_dict[gene].values)
                    )
                ):
                    tcrs.append(gene)

                continue

            for d in cls.ref_dict[gene]:
                if (
                    functionality == "any"
                    or (
                        functionality in ("F", "P", "ORF")
                        and cls.ref_dict[gene][d] == functionality
                    )
                    or (functionality == "NF" and cls.ref_dict[gene][d] in ("P", "ORF"))
                ):
                    tcrs.append(gene + "*" + d)

        return frozenset(tcrs)


class HomoSapiensTCRQueryEngine(TCRQueryEngine):
    ref_dict = HOMOSAPIENS_TCR


class MusMusculusTCRQueryEngine(TCRQueryEngine):
    ref_dict = MUSMUSCULUS_TCR


class HLAQueryEngine(GeneQueryEngine):
    @staticmethod
    def query(precision: str, functionality: str) -> FrozenSet[str]:
        if precision == "allele":
            warn(
                "tidytcells is not fully aware of all HLA alleles, and the "
                "highest resolution it can provide is up to the level of the "
                "protein (two allele designations)."
            )

        hlas = []

        for gene in HOMOSAPIENS_MHC:
            if precision == "gene":
                hlas.append(gene)
                continue

            for d1 in HOMOSAPIENS_MHC[gene]:
                for d2 in HOMOSAPIENS_MHC[gene][d1]:
                    if not d2.isdigit():  # ignore G/P groups
                        continue

                    hlas.append(gene + "*" + d1 + ":" + d2)

        return frozenset(hlas)


class MusMusculusMHCQueryEngine(GeneQueryEngine):
    @staticmethod
    def query(precision: str, functionality: str) -> FrozenSet[str]:
        if precision == "allele":
            warn(
                "tidytcells is not aware of Mus musculus MHC alleles at all, "
                "and can only provide up to the level of the gene."
            )

        return frozenset(MUSMUSCULUS_MHC)
