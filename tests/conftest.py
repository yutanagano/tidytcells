import pytest

from tidytcells import config


@pytest.fixture(autouse=True)
def reset_receptor_database_config():
    """
    Guard against config-state leaks: if a test fails after
    set_receptor_database(...) but before its manual reset, the OGRDB state
    must not pollute subsequent tests.
    """
    config.reset_receptor_database()
    yield
    config.reset_receptor_database()
