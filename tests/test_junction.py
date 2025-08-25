import pytest
from tidytcells import junction


class Teststandardize:
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
        result = junction.standardize(seq=seq)

        assert result == seq

    @pytest.mark.parametrize("seq", ("123456", "ASDFGHJKL", "A?AAAA", "AAAXAA"))
    def test_various_rejections(self, seq, caplog):
        result = junction.standardize(seq=seq)
        assert "not a valid amino acid sequence" in caplog.text
        assert result is None

    @pytest.mark.parametrize(
        ("seq", "expected"),
        (
            ("casqyf", "CASQYF"),
            ("ASQY", "CASQYF"),
            ("CASQY", "CASQYF"),
            ("ASQYF", "CASQYF"),
        ),
    )
    def test_various_corrections(self, seq, expected):
        result = junction.standardize(seq=seq)

        assert result == expected

    @pytest.mark.parametrize("seq", (None, 1, True, 5.4))
    def test_bad_input_type(self, seq):
        with pytest.raises(TypeError):
            junction.standardize(seq=seq)

    @pytest.mark.parametrize("seq", ("ASQY", "CASQY", "ASQYF", "ASQYW"))
    def test_strict(self, seq, caplog):
        result = junction.standardize(seq=seq, strict=True)
        assert "not a valid junction" in caplog.text
        assert result is None

    def test_standardize(self):
        result = junction.standardize(seq="CASSPGGADRRIDGYTF")

        assert result == "CASSPGGADRRIDGYTF"

    def test_log_failures(self, caplog):
        junction.standardize(seq="123456", log_failures=False)
        assert len(caplog.records) == 0

    def test_on_fail(self, caplog):
        result = junction.standardize("foobarbaz", on_fail="keep")
        assert "Failed to standardize foobarbaz" in caplog.text
        assert result == "foobarbaz"

    @pytest.mark.parametrize(
        ("seq", "j_symbol", "species", "expected"),
        (
                ("AELNAGNNRKLI", "TRAJ38*01", "homosapiens", "CAELNAGNNRKLIW"),
                ("AELNAGNNRKLI", "TRAJ38*01", "musmusculus", "CAELNAGNNRKLIW"),
                ("AELNAGNNRKLI", None, "homosapiens", "CAELNAGNNRKLIF"),
                ("AELNAGNNRKLI", None, "musmusculus", "CAELNAGNNRKLIF"),
                ("AELNAGNNRKLI", None, "musmusculus", "CAELNAGNNRKLIF"),
                ("AAAAWF", "IGHJ5*01", "homosapiens", "CAAAAWFW"),
                ("AAAAWF", None, "homosapiens", "CAAAAWF"),
                ("AAAAAA", "IGHJ5", "homosapiens", "CAAAAAAW"),
                ("AAAAAA", "IGHJ", "homosapiens", "CAAAAAAW"), # all of IGH have W, all of IGL have F
                ("AAAAAA", "IGLJ", "homosapiens", "CAAAAAAF"),
                ("AAAAAA", "TRAJ", "homosapiens", None), # but TRA is ambiguous
        ),
    )
    def test_j_symbol(self, seq, j_symbol, species, expected):
        result = junction.standardize(seq=seq, j_symbol=j_symbol, species=species)

        assert result == expected

    def test_j_symbol_fail(self, caplog):
        result = junction.standardize("AAAAAA", j_symbol="TRAJ", on_fail="keep")
        assert "conserved amino acid was ambiguous" in caplog.text
        assert result == "AAAAAA"

    def test_j_symbol_fail_species(self):
        with pytest.raises(ValueError, match="not supported for IG genes for species musmusculus"):
            junction.standardize("AAAAAA", j_symbol="IGH", species="musmusculus", on_fail="keep")

