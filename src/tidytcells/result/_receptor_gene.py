from tidytcells._resources import IG_SYMBOLS_BY_SP_DB, TR_SYMBOLS_BY_SP_DB
from tidytcells._utils.alignment import get_compatible_symbols


def _get_allele_tree(species: str, receptor_type: str, germline_database_source: str) -> dict:
    registry = TR_SYMBOLS_BY_SP_DB if receptor_type == "TR" else IG_SYMBOLS_BY_SP_DB
    return registry[species][germline_database_source]["symbols"]


def _get_allele_data(allele_tree: dict, allele_symbol: str) -> dict:
    gene_name, allele_designation = allele_symbol.split("*")
    return allele_tree[gene_name][allele_designation]


class ReceptorGene:
    """
    A wrapper object for the receptor gene.

    If standardization was successful, this object provides access to the standardized allele/gene/subgroup and other properties.
    When failed, the error message(s) and attempted partially standardized gene symbol can be retrieved.
    """

    def __init__(
        self,
        original_input: str,
        error: str | None,
        gene_name: str | None = None,
        allele_designation: str | None = None,
        subgroup_name: str | None = None,
        species: str | None = None,
        germline_database: dict | None = None,
    ) -> None:
        self._original_input = original_input
        self._error = error
        self._gene_name = gene_name
        self._allele_designation = allele_designation
        self._subgroup_name = subgroup_name
        self._species = species
        self._germline_database = germline_database

        self._highest_precision_symbol = None

        if self._gene_name is not None and self._allele_designation is not None:
            self._highest_precision_symbol = (
                f"{self._gene_name}*{self._allele_designation}"
            )
        elif self._gene_name is not None:
            self._highest_precision_symbol = self._gene_name
        elif self._subgroup_name is not None:
            self._highest_precision_symbol = self._subgroup_name

    def __str__(self) -> str:
        str_repr = self.symbol

        if str_repr is not None:
            return str_repr
        else:
            return ""

    @property
    def original_input(self) -> str | None:
        """The original input symbol."""
        return self._original_input

    @property
    def error(self) -> str | None:
        """The error message, only if standardization failed, otherwise ``None``."""
        return self._error

    @property
    def is_standardized(self) -> bool:
        """``True`` if the standardization was successful, ``False`` otherwise."""
        return self.error is None

    @property
    def attempted_fix(self) -> str | None:
        """The best attempt at fixing the input symbol, only of standardization failed, if the standardization was a success this returns ``None``."""
        if not self.is_standardized:
            return self._highest_precision_symbol

    @property
    def symbol(self) -> str | None:
        """The allele, gene or subgroup (whichever is most precise) if standardization was successful, otherwise ``None``."""
        if self.is_standardized:
            return self._highest_precision_symbol

    @property
    def allele(self) -> str | None:
        """The allele name, if standardization was successful and allele-level information is available, otherwise ``None``."""
        if (
            self.is_standardized
            and self._allele_designation is not None
            and self._gene_name is not None
        ):
            return f"{self._gene_name}*{self._allele_designation}"

    @property
    def gene(self) -> str | None:
        """The gene name, if standardization was successful and gene-level information is available, otherwise ``None``."""
        if self.is_standardized:
            return self._gene_name

    @property
    def subgroup(self) -> str | None:
        """The subgroup name, if standardization was successful, otherwise ``None``."""
        if self.is_standardized:
            return self._subgroup_name

    @property
    def locus(self) -> str | None:
        """
        The locus of the gene.
        This is typically the three-letter code ('TRA', 'TRB', 'TRG', 'TRD', 'IGH', 'IGL', 'IGK'),
        but for TRAV/DV genes, 'TRA/D' is returned.
        """
        if self.is_standardized:
            locus = self.symbol[0:3]
            if "/D" in self.symbol:
                locus += "/D"
            return locus

    @property
    def receptor_type(self) -> str | None:
        """'TR' for T cell receptor genes, or 'IG' for antibody genes if standardization was successful, otherwise ``None``."""
        if self.is_standardized:
            return self.symbol[0:2]

    @property
    def gene_type(self) -> str | None:
        """The gene type ('V', 'D' or 'J'), if standardization was successful, otherwise ``None``."""
        if self.is_standardized:
            return self.symbol[3]

    @property
    def species(self) -> str | None:
        """The species used to validate the gene name."""
        return self._species

    @property
    def germline_database(self) -> dict | None:
        """
        The germline database that was used to standardize this receptor gene.
        The return value will contain a dictionary specifying the source (``"IMGT"`` / ``"OGRDB"``), resource and version number.
        """
        return self._germline_database

    def get_all_alleles(self, enforce_functional: bool = True) -> list[str] | None:
        """
        Get all alleles related to the standardized symbol.

        Parameters
        ----------
        enforce_functional
            If ``True``, only functional alleles are returned.

        Returns
        -------
        list[str] | None
            A list of allele names.
        """
        if self.is_standardized:
            allele_tree = _get_allele_tree(
                self.species, self.receptor_type, self.germline_database["source"]
            )

            return get_compatible_symbols(
                self.symbol, allele_tree, self.gene_type, self.locus, enforce_functional
            )

    def get_aa_sequences(
        self, sequence_type: str = "ALL", enforce_functional: bool = True
    ) -> dict[str, dict[str, str]] | dict[str, str | None] | None:
        """
        Get amino acid sequence information related to the alleles of the standardized symbol.

        Parameters
        ----------
        sequence_type
            Which sequence to return. This can be:
            - For V genes: 'FR1', 'FR2', 'FR3', 'CDR1', 'CDR2', 'V-REGION'
            - For D genes: 'D-REGION'
            - For J genes: 'J-REGION', 'J-MOTIF'
            - Or 'ALL' to return all available sequences
        enforce_functional
            If ``True``, only information for functional alleles is returned.

        Returns
        -------
        dict[str, dict[str, str]] | dict[str, str | None] | None
            A dictionary with allele names as keys and sequences as values.
            When sequence_type is 'ALL', the result is a nested dictionary with allele names
            as outer keys, sequence types as inner keys, and sequences as inner values.
        """
        sequence_type = sequence_type.upper()
        sequence_type = (
            sequence_type + "-IMGT"
            if sequence_type in {"FR1", "FR2", "FR3", "CDR1", "CDR2"}
            else sequence_type
        )

        if self.is_standardized:
            allele_tree = _get_allele_tree(
                self.species, self.receptor_type, self.germline_database["source"]
            )

            alleles_of_interest = self.get_all_alleles(enforce_functional)

            if sequence_type == "ALL":
                return {
                    allele: _get_allele_data(allele_tree, allele)
                    for allele in alleles_of_interest
                }
            else:
                return {
                    allele: _get_allele_data(allele_tree, allele).get(sequence_type)
                    for allele in alleles_of_interest
                }
