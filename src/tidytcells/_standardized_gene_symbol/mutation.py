import re

SUBSTITUTION = r"(?:[A-Z]:)?[ACDEFGHIKLMNPQRSTVWY]\d+[ACDEFGHIKLMNPQRSTVWY]"
MUTATION_REGEX = re.compile(
    rf"{SUBSTITUTION}(?:(?:\s*,\s*|\s+){SUBSTITUTION})*", re.IGNORECASE
)
MUTANT_KEYWORD_REGEX = re.compile(r"^(.+?)\s+(?:mutant|mutation)$", re.IGNORECASE)


def _normalise(mutation: str) -> str:
    return re.sub(r"\s*,\s*|\s+", ", ", mutation).upper()


def validate_mutation(mutation: str) -> str | None:
    """
    Return an error message if the mutation is not a valid substitution
    description, otherwise ``None``.
    """
    if not MUTATION_REGEX.fullmatch(mutation):
        return f"Invalid mutation: {mutation}"

    return None


def interpret_remainder(remainder: str) -> tuple[str | None, str | None]:
    """
    Interpret the text that follows a gene symbol as a potential mutation
    specification.

    Returns a ``(mutation, error)`` tuple:

    - ``(mutation, None)`` if the remainder is a valid mutation specification,
      with or without a trailing "mutant"/"mutation" keyword.
    - ``(None, "Invalid mutation: X")`` if a keyword is present but the body
      is not a valid substitution list.
    - ``(None, None)`` if the remainder does not look like a mutation
      specification; the caller should treat it as trailing text.
    """
    remainder = remainder.strip()
    if not remainder:
        return None, None

    keyword_match = MUTANT_KEYWORD_REGEX.match(remainder)
    if keyword_match:
        mutation = _normalise(keyword_match.group(1))
        error = validate_mutation(mutation)
        if error is None:
            return mutation, None
        return None, error

    if MUTATION_REGEX.fullmatch(remainder):
        return _normalise(remainder), None

    return None, None
