import pytest
from tidytcells import aa
import warnings


class TestStandardise:
    @pytest.mark.parametrize("seq", ("KLGGALQAK", "LLQTGIHVRVSQPSL", "SQLLNAKYL"))
    def test_alreadY_correct(self, seq):
        result = aa.standardize(seq=seq)

        assert result == seq

    @pytest.mark.parametrize("seq", ("123456", "ASDFGHJKL", "A?AAAA", "AAAXAA"))
    def test_various_rejections(self, seq):
        with pytest.warns(UserWarning, match="not a valid amino acid"):
            result = aa.standardize(seq=seq)

        assert result == None

    @pytest.mark.parametrize(("seq", "expected"), (("klgak", "KLGAK"),))
    def test_various_corrections(self, seq, expected):
        result = aa.standardize(seq=seq)

        assert result == expected

    @pytest.mark.parametrize("seq", (None, 1, True, 5.4))
    def test_bad_input_type(self, seq):
        with pytest.raises(TypeError):
            aa.standardize(seq=seq)

    def test_standardise(self):
        result = aa.standardise(seq="KLGAK")

        assert result == "KLGAK"

    def test_suppress_warnings(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            aa.standardize(seq="123456", suppress_warnings=True)

    def test_on_fail(self):
        with pytest.warns(UserWarning):
            result = aa.standardize("foobarbaz", on_fail="keep")

        assert result == "foobarbaz"
