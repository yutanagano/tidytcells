class MhGene:
    """
    A wrapper object for the MH gene.

    If standardization was successful, this object provides access to the
    standardized allele/gene and other properties. When failed, the error
    message(s) and attempted partially standardized gene symbol can be
    retrieved.
    """

    def __init__(
        self,
        original_input: str,
        error: str | None,
        gene_name: str | None = None,
        allele_designation: list[str] | None = None,
        species: str | None = None,
        mutation: str | None = None,
        mh_database: dict | None = None,
    ) -> None:
        self._original_input = original_input
        self._error = error
        self._gene_name = gene_name
        self._allele_designation = (
            allele_designation
            if allele_designation is not None and len(allele_designation) > 0
            else None
        )
        self._species = species
        self._mutation = mutation
        self._mh_database = mh_database

        self._highest_precision_symbol = self._gene_name

        if self._gene_name is not None and self._allele_designation is not None:
            self._highest_precision_symbol = (
                f"{self._gene_name}*{':'.join(self._allele_designation)}"
            )

        if self._highest_precision_symbol is not None and self._mutation is not None:
            self._highest_precision_symbol += f" {self._mutation} mutant"

    def __str__(self) -> str:
        str_repr = self.symbol

        if str_repr is not None:
            return str_repr
        else:
            return ""

    @property
    def original_input(self) -> str | None:
        """
        The original input symbol.
        """
        return self._original_input

    @property
    def error(self) -> str | None:
        """
        The error message, only if standardization failed, otherwise ``None``.
        """
        return self._error

    @property
    def is_standardized(self) -> bool:
        """
        ``True`` if the standardization was successful, ``False`` otherwise.
        """
        return self.error is None

    @property
    def attempted_fix(self) -> str | None:
        """
        The best attempt at fixing the input symbol.

        This is only set if standardization failed. Otherwise this is ``None``.
        """
        if not self.is_standardized:
            return self._highest_precision_symbol

    @property
    def symbol(self) -> str | None:
        """
        The allele or gene (whichever is most precise), with the mutation
        description appended if the input was a mutant (e.g. ``'HLA-A*02:01
        K66A mutant'``).

        Set to ``None`` if standardization was not successful.
        """
        if self.is_standardized:
            return self._highest_precision_symbol

    @property
    def allele(self) -> str | None:
        """
        The allele name, if standardization was successful and allele-level
        information is available, otherwise ``None``.
        """
        if (
            self.is_standardized
            and self._allele_designation is not None
            and self._gene_name is not None
        ):
            return f"{self._gene_name}*{':'.join(self._allele_designation)}"

    @property
    def gene(self) -> str | None:
        """
        The gene name, if standardization was successful, otherwise ``None``.
        """
        if self.is_standardized and self._gene_name is not None:
            return self._gene_name

    @property
    def protein(self) -> str | None:
        """
        The protein name, if standardization was successful and protein-level
        information is available, otherwise ``None``.
        """
        if (
            self.is_standardized
            and self._allele_designation is not None
            and self._gene_name is not None
        ):
            return f"{self._gene_name}*{':'.join(self._allele_designation[:2])}"

    @property
    def mutation(self) -> str | None:
        """
        The mutation description, if standardization was successful and the
        input was a mutant, otherwise ``None``.

        For example, ``'K66A'`` for ``'HLA-A*02:01 K66A mutant'``.
        """
        if self.is_standardized:
            return self._mutation

    @property
    def species(self) -> str | None:
        """
        The species used to validate the gene name.
        """
        return self._species

    @property
    def mh_database(self) -> dict | None:
        """
        The MH database that was used to standardize this receptor gene.
        The return value will contain a dictionary specifying the source (``"MRO"`` / ``"IPD"``), resource and version number.
        """
        return self._mh_database
