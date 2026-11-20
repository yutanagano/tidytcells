import itertools
import re
from functools import partial
from collections.abc import Callable

from tidytcells import _utils
from tidytcells._resources import MH_SYMBOLS_BY_SP_DB
from tidytcells._standardized_gene_symbol.mutation import interpret_remainder
from tidytcells.result import MhGene

ALLELE_SYMBOL_REGEX = re.compile(
    r"^(?P<gene>[A-Z0-9\-\.\:\/]+)(\*(?P<designation>[A-Z\d:]+))?"
)
HLA_SYMBOL_REGEX = re.compile(
    r"^(?P<gene>[A-Z0-9\-\.\:\/]+)(\*(?P<designation>[\d:]+G?P?)[LSCAQN]?)?"
)
HLA_CLASS_II_REGEX = re.compile(
    r"^(?P<gene>(HLA-)?(D[PQ][AB]|DRB|TAP)\d)(\*?(?P<designation>[\d:]+G?P?)[LSCAQN]?)?"
)


def standardize_allele(
    symbol: str,
    *,
    species: str,
    mh_database_info: dict,
    allow_zero_padded_three_digit_fields: bool,
) -> MhGene:
    """
    Standardize an MH symbol with numbered alleles (e.g. HLA-A*02:01,
    Mamu-A1*001:01). Input is matched case-insensitively; results are
    returned in the case used by the database.
    """
    valid_symbols, allele_synonyms, gene_synonyms = _mh_data(species, mh_database_info["source"])
    canonical_gene_names = {gene.upper(): gene for gene in valid_symbols}
    prefixes = sorted(
        {gene.split("-")[0].upper() for gene in valid_symbols if "-" in gene}
    )

    if species == "homosapiens":
        symbol_regex = HLA_SYMBOL_REGEX
        class_ii_regex = HLA_CLASS_II_REGEX
    else:
        symbol_regex = ALLELE_SYMBOL_REGEX
        class_ii_regex = None

    cleaned_symbol, gene_name, allele_designation, remainder = parse_allele_symbol(
        symbol, symbol_regex, class_ii_regex
    )

    validate = partial(
        _validation_error,
        species=species,
        canonical_gene_names=canonical_gene_names,
        valid_symbols=valid_symbols,
    )

    if validate(gene_name, allele_designation) is not None:
        # the symbol is not valid as parsed; try to repair it

        if cleaned_symbol in allele_synonyms:
            _, gene_name, allele_designation, _ = parse_allele_symbol(
                allele_synonyms[cleaned_symbol], symbol_regex, class_ii_regex
            )

        if (
            validate(gene_name, allele_designation) is not None
            and gene_name in gene_synonyms
        ):
            gene_name = gene_synonyms[gene_name].upper()

        if validate(gene_name, allele_designation) is not None:
            gene_name = _add_missing_prefix(gene_name, prefixes, canonical_gene_names)
            if species == "homosapiens":
                gene_name = gene_name.replace("CW", "C")

        if validate(gene_name, allele_designation) is not None:
            gene_name, allele_designation = _handle_forgotten_asterisk(
                gene_name, allele_designation, canonical_gene_names
            )

        if validate(gene_name, allele_designation) is not None:
            is_valid = partial(validate, gene_name)
            allele_designation = _handle_forgotten_colon(
                allele_designation, allow_zero_padded_three_digit_fields, is_valid
            )

        if validate(gene_name, allele_designation) is not None:
            is_valid = partial(validate, gene_name)
            allele_designation = _try_different_amounts_of_leading_zeros(
                allele_designation, allow_zero_padded_three_digit_fields, is_valid
            )

        if validate(gene_name, allele_designation) is not None:
            candidate = gene_name[:-1] if gene_name.endswith("1") else gene_name + "1"
            if validate(candidate, allele_designation) is None:
                gene_name = candidate

    error = validate(gene_name, allele_designation)

    # only if the symbol part resolves to a standard symbol is the
    # trailing text examined as a mutation specification
    mutation = None
    if error is None:
        mutation, mutation_error = interpret_remainder(remainder)
        if mutation_error is not None:
            error = mutation_error
        elif mutation is None and remainder:
            error = f"Unexpected trailing text: {remainder}"

    return MhGene(
        original_input=symbol,
        error=error,
        gene_name=canonical_gene_names.get(gene_name, gene_name),
        allele_designation=_canonical_allele_designation(
            gene_name, allele_designation, canonical_gene_names, valid_symbols
        ),
        species=species,
        mutation=mutation,
        mh_database=mh_database_info,
    )


def parse_allele_symbol(
    symbol: str,
    symbol_regex: re.Pattern,
    class_ii_regex: re.Pattern | None = None,
) -> tuple[str, str, list[str], str]:
    """
    Split an MH allele-style symbol (e.g. HLA-A*02:01, Mamu-A1*001:01) into
    its parts.

    Returns ``(cleaned_symbol, gene_name, allele_designation, remainder)``.
    Internal whitespace is kept while matching, so the symbol regex's match
    end (not the first space) defines where the symbol ends; anything after
    it is returned as the remainder.
    """
    symbol = _utils.clean_preserving_whitespace(symbol).upper()
    symbol = _replace_periods_between_digits_with_colon(symbol)

    for regex in filter(None, (class_ii_regex, symbol_regex)):
        parse_attempt = regex.match(symbol)
        if parse_attempt:
            symbol_part = symbol[: parse_attempt.end()]
            remainder = symbol[parse_attempt.end() :].strip()
            return (
                _utils.clean_and_uppercase(symbol_part),
                parse_attempt.group("gene"),
                _listify_allele_designation(parse_attempt.group("designation")),
                remainder,
            )

    return _utils.clean_and_uppercase(symbol), symbol, [], ""


def _mh_data(
    species: str, mh_database_source: str
) -> tuple[dict[str, dict], dict[str, str], dict[str, str]]:
    """
    The valid-symbol, synonym, and HGNC gene-synonym data for a species and
    MH database.
    """
    data = MH_SYMBOLS_BY_SP_DB[species][mh_database_source]
    return (
        data["symbols"],
        data.get("synonyms_allele", {}),
        data.get("synonyms_gene", {}),
    )


def _replace_periods_between_digits_with_colon(string: str) -> str:
    return re.sub(r"(?<=\d)\.(?=\d)", ":", string)


def _listify_allele_designation(allele_designation: str | None) -> list[str]:
    if allele_designation is None:
        return []

    allele_designation = re.sub(r"(?<=\d)[LSCAQN]$", "", allele_designation)

    return [d.zfill(2) if d.isdigit() else d for d in allele_designation.split(":")]


def _is_group(allele_designation: list[str]) -> bool:
    if not allele_designation:
        return False

    return allele_designation[-1].endswith(("G", "P"))


def _validation_error(
    gene_name: str,
    allele_designation: list[str],
    *,
    species: str,
    canonical_gene_names: dict[str, str],
    valid_symbols: dict[str, dict],
) -> str | None:
    if species == "homosapiens" and gene_name == "B2M" and not allele_designation:
        return None

    canonical_gene_name = canonical_gene_names.get(gene_name)
    if canonical_gene_name is None:
        return "Unrecognized gene name"

    # Verify allele designators up to the level of the protein (or G/P)
    designation = list(allele_designation)
    if not _is_group(allele_designation):
        designation = designation[:2]
    current_root = valid_symbols[canonical_gene_name]

    while designation:
        # case-insensitive: a few designator keys in the data are not
        # uppercase (e.g. "Mamu-B*074:03Sp")
        field = designation.pop(0).upper()
        current_root = next(
            (value for key, value in current_root.items() if key.upper() == field),
            None,
        )
        if current_root is None:
            return "Nonexistent allele for recognized gene"

    # If there are designator fields past the protein level, just make sure
    # they look like legitimate designator field values
    if not _is_group(allele_designation) and len(allele_designation) > 2:
        further_designators = allele_designation[2:]

        if len(further_designators) > 2:
            return "Too many allele designators"

        for field in further_designators:
            if not field.isdigit():
                return "Non-numerical allele designators"

            if len(field) < 2:
                return "Non-2-digit allele designators"

    return None


def _add_missing_prefix(
    gene_name: str, prefixes: list[str], canonical_gene_names: dict[str, str]
) -> str:
    if any(gene_name.startswith(f"{prefix}-") for prefix in prefixes):
        return gene_name

    for prefix in prefixes:
        candidate = f"{prefix}-{gene_name}"
        if candidate in canonical_gene_names:
            return candidate

    return f"{prefixes[0]}-{gene_name}" if len(prefixes) == 1 else gene_name


def _handle_forgotten_asterisk(
    gene_name: str,
    allele_designation: list[str],
    canonical_gene_names: dict[str, str],
) -> tuple[str, list[str]]:
    if allele_designation:
        return gene_name, allele_designation

    for split_idx in range(len(gene_name) - 1, 0, -1):
        gene, designation = gene_name[:split_idx], gene_name[split_idx:]

        if gene in canonical_gene_names and re.fullmatch(r"[\d:]+", designation):
            return gene, designation.split(":")

    return gene_name, allele_designation


def _split_into_designator_fields(
    digits: str, allow_zero_padded_three_digit_fields: bool, max_fields: int = 4
) -> list[list[str]]:
    candidates = [
        fields
        for fields in _all_designator_field_splits(
            digits, allow_zero_padded_three_digit_fields
        )
        if 2 <= len(fields) <= max_fields
    ]
    return sorted(candidates, key=lambda fields: [len(field) for field in fields])


def _all_designator_field_splits(
    digits: str, allow_zero_padded_three_digit_fields: bool
) -> list[list[str]]:
    if not digits:
        return [[]]

    splits = []
    for n in (2, 3):
        field = digits[:n]
        zero_padded = n == 3 and field.startswith("0")

        if len(field) == n and (
            allow_zero_padded_three_digit_fields or not zero_padded
        ):
            for tail in _all_designator_field_splits(
                digits[n:], allow_zero_padded_three_digit_fields
            ):
                splits.append([field] + tail)

    return splits


def _handle_forgotten_colon(
    allele_designation: list[str],
    allow_zero_padded_three_digit_fields: bool,
    is_valid: Callable[[list[str]], str | None],
) -> list[str]:
    if not allele_designation or not allele_designation[0].isdigit():
        return allele_designation

    original = allele_designation

    for fields in _split_into_designator_fields(
        original[0], allow_zero_padded_three_digit_fields
    ):
        candidate = fields + original[1:]
        if is_valid(candidate) is None:
            return candidate

    return original


def _try_different_amounts_of_leading_zeros(
    allele_designation: list[str],
    allow_zero_padded_three_digit_fields: bool,
    is_valid: Callable[[list[str]], str | None],
) -> list[str]:
    original = allele_designation
    reformatted_allele_designators = [
        [f"{int(ad):02}", f"{int(ad):03}"] if ad.isdigit() else [ad]
        for ad in original[:2]
    ]

    for new_designators in itertools.product(*reformatted_allele_designators):
        candidate = list(new_designators) + original[2:]
        if is_valid(candidate) is None:
            return candidate

    return original


def _canonical_allele_designation(
    gene_name: str,
    allele_designation: list[str],
    canonical_gene_names: dict[str, str],
    valid_symbols: dict[str, dict],
) -> list[str]:
    canonical_gene_name = canonical_gene_names.get(gene_name, gene_name)
    node = valid_symbols.get(canonical_gene_name, {})
    canonical = []

    for field in allele_designation:
        canonical_field = next((key for key in node if key.upper() == field), field)
        canonical.append(canonical_field)
        node = node.get(canonical_field, {})

    return canonical
