import pytest
from tidytcells import junction
import warnings


class TestStandardise:
    @pytest.mark.parametrize(
        "seq",
        (
            "CASSPGGADRRIDGYTF",
            "CASSLMPGQGSYEQYF",
            "CSARDPGDDKPQHF",
            "CAALRATGGNNKLTF",
            "CVVAVSNFGNEKLTF",
            "CAASASFGDNSKLIW",
        ),
    )
    def test_already_correct(self, seq):
        result = junction.standardise(seq=seq)

        assert result == seq

    @pytest.mark.parametrize("seq", ("123456", "ASDFGHJKL", "A?AAAA", "AAAXAA"))
    def test_various_rejections(self, seq):
        with pytest.warns(UserWarning, match="not a valid amino acid sequence"):
            result = junction.standardise(seq=seq)

        assert result is None

    @pytest.mark.parametrize(
        ("seq", "expected"),
        (
            ("casqyf", "CASQYF"),
            ("ASQY", "CASQYF"),
            ("CASQY", "CCASQYF"),
            ("ASQYF", "CASQYFF"),
        ),
    )
    def test_various_corrections(self, seq, expected):
        result = junction.standardise(seq=seq)

        assert result == expected

    @pytest.mark.parametrize("seq", (None, 1, True, 5.4))
    def test_bad_input_type(self, seq):
        with pytest.raises(TypeError):
            junction.standardise(seq=seq)

    @pytest.mark.parametrize("seq", ("ASQY", "CASQY", "ASQYF", "ASQYW"))
    def test_strict(self, seq):
        with pytest.warns(UserWarning, match="not a valid junction"):
            result = junction.standardise(seq=seq, strict=True)

        assert result is None

    def test_standardize(self):
        result = junction.standardize(seq="CASSPGGADRRIDGYTF")

        assert result == "CASSPGGADRRIDGYTF"

    def test_suppress_warnings(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            junction.standardise(seq="123456", suppress_warnings=True)

    def test_on_fail(self):
        with pytest.warns(UserWarning):
            result = junction.standardize("foobarbaz", on_fail="keep")

        assert result == "foobarbaz"
