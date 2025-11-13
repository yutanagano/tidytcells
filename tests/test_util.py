import pytest
from tidytcells._utils.alignment import get_compatible_symbols


class TestUtil:

    pass # todo add tests for alignment / get_compatible_symbols

    #
    # @pytest.mark.parametrize(
    #     ("seq", "extension", "expected"),
    #     (
    #         ("TRAJ47", "TRAJ47*02", True),
    #         ("TRBJ2-2", "TRBJ2-2*01", True),
    #         ("TRBJ2-2", "TRBJ2-2P", False),
    #     ))
    # def test_is_valid_extension(self, seq, extension, expected):
    #     result = get_compatible_symbols(seq, extension)
    #
    #     assert result == expected
    #
    #
