"""
Global configuration for tidytcells.

The currently selected germline databases are stored as module-level variables
and can be changed with the `set_*` functions (and restored to their defaults
with the `reset_*` functions), e.g.::

    import tidytcells
    tidytcells.config.set_receptor_database("OGRDB")
    tidytcells.config.reset_receptor_database()
"""

VALID_RECEPTOR_DATABASES = {"IMGT", "OGRDB"}

VALID_MH_DATABASES = {"IPD", "MRO"}

DEFAULT_RECEPTOR_DATABASE = "IMGT"
DEFAULT_MH_DATABASE = "MRO"

_receptor_database = DEFAULT_RECEPTOR_DATABASE
_mh_database = DEFAULT_MH_DATABASE


def get_receptor_database() -> str:
    return _receptor_database


def get_mh_database() -> str:
    return _mh_database


def set_receptor_database(name: str) -> None:
    global _receptor_database
    if name not in VALID_RECEPTOR_DATABASES:
        raise ValueError(f'"receptor_database" must be one of {sorted(VALID_RECEPTOR_DATABASES)}, got {name}.')
    _receptor_database = name


def set_mh_database(name: str) -> None:
    global _mh_database
    if name not in VALID_MH_DATABASES:
        raise ValueError(f'"mh_database" must be one of {sorted(VALID_MH_DATABASES)}, got {name}.')
    _mh_database = name


def reset_receptor_database() -> None:
    global _receptor_database
    _receptor_database = DEFAULT_RECEPTOR_DATABASE


def reset_mh_database() -> None:
    global _mh_database
    _mh_database = DEFAULT_MH_DATABASE


def reset_all_defaults() -> None:
    reset_receptor_database()
    reset_mh_database()
