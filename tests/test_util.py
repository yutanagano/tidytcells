import pytest
from tidytcells import aa
from tidytcells._utils.conserved_aa_lookup import is_valid_extension


class TestUtil:

    @pytest.mark.parametrize(
        ("seq", "extension", "expected"),
        (
            ("TRAJ47", "TRAJ47*02", True),
            ("TRBJ2-2", "TRBJ2-2*01", True),
            ("TRBJ2-2", "TRBJ2-2P", False),
        ))
    def test_is_valid_extension(self, seq, extension, expected):
        result = is_valid_extension(seq, extension)

        assert result == expected

