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
            ("CASQY", "CCASQYF"),
            ("ASQYF", "CASQYFF"),
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
                ("AAAAWF", None, "homosapiens", "CAAAAWFF"),
                ("AAAAAA", "IGHJ5", "homosapiens", "CAAAAAAW"),
                ("AAAAAA", "IGHJ", "homosapiens", "CAAAAAAW"), # all of IGH have W, all of IGL have F
                ("AAAAAA", "IGLJ", "homosapiens", "CAAAAAAF"),
                ("AAAAAA", "TRAJ", "homosapiens", "CAAAAAAF"), # but TRA is ambiguous
                ("AAAAAA", "TRAJ35", "homosapiens", "CAAAAAAC"), # non-canonical C ending for TRAJ35*01 human
                ("AAAAAA", "TRAJ35", "musmusculus", "CAAAAAAF"), # ...but not for mouse
        ),
    )
    def test_j_symbol(self, seq, j_symbol, species, expected):
        result = junction.standardize(seq=seq, j_symbol=j_symbol, species=species)

        assert result == expected

    def test_j_symbol_ambiguous_default(self, caplog):
        result = junction.standardize("AAAAAA", j_symbol="TRAJ", on_fail="keep", j_strict=False)
        assert result == "CAAAAAAF"

    def test_j_symbol_fail(self, caplog):
        result = junction.standardize("AAAAAA", j_symbol="TRAJ", on_fail="keep", j_strict=True)
        assert "conserved amino acid was ambiguous" in caplog.text
        assert result == "AAAAAA"

    def test_j_symbol_fail_species(self):
        with pytest.raises(ValueError, match="not supported for IG genes for species musmusculus"):
            junction.standardize("AAAAAA", j_symbol="IGH", species="musmusculus", on_fail="keep")

    @pytest.mark.parametrize(
        ("seq", "trimming", "species", "expected"),
        (
                ("AELNAGNNRKLI", False, "homosapiens", "CAELNAGNNRKLIF"),
                ("AELNAGNNRKLI", True, "homosapiens", "CAELNAGNNRKLIF"),
                ("CAELNAGNNRKLI", True, "homosapiens", "CCAELNAGNNRKLIF"),
                ("FCAVVFNMDSNYQLIW", True, "homosapiens", "CAVVFNMDSNYQLIW"), # trimming at start
                ("YFCAVVFNMDSNYQLIW", True, "homosapiens", "CAVVFNMDSNYQLIW"), # trimming at start
                ("YFCAVVFNMDSNYQLIWGAGTKL", True, "homosapiens", "CAVVFNMDSNYQLIW"), # should trim until the outermost conserved aa (W)
                ("DSSIYLCSVEATRADTQYFGPGFTRLTVL", True, "homosapiens", "CSVEATRADTQYF"), # trimming in middle
                ("CCEFTGGGNKLTF", True, "homosapiens", "CCEFTGGGNKLTF"), # stays the same
                ("CVGWHEQY", True, "homosapiens", "CCVGWHEQYF"), # stays the same ('CVGW' is too short)
                ("CASSPGVFGANVLTFG", True, "homosapiens", "CASSPGVFGANVLTF"), #make sure to trim at the last F
        ),
    )
    def test_trimming(self, seq, trimming, species, expected):
        result = junction.standardize(seq=seq, trimming=trimming, species=species)

        assert result == expected

    @pytest.mark.parametrize(
        ("seq", "trimming", "j_symbol", "species", "expected"),
        (
                ("AELNAGNNRKLI", True, "TRAJ38*01", "homosapiens", "CAELNAGNNRKLIW"), # trimming shouldnt matter here
                ("AELNAGNNRKLI", False, "TRAJ38*01", "homosapiens", "CAELNAGNNRKLIW"),
                ("YFCAVVFNMDSNYQLIWFGAGTKL", True, "TRA", "homosapiens", "CAVVFNMDSNYQLIWF"), # default: look for F
                ("YFCAVVFNMDSNYQLIWGQGGTKL", True, "TRAJ38*01", "homosapiens", "CAVVFNMDSNYQLIW"), # J gene with conserved W
                ("YFCAVVFNMDSNYQLIWFGAGTKL", True, "IGLJ", "homosapiens", "CAVVFNMDSNYQLIWF"), # J gene with conserved F
                ("YFCAELNAGNNRKLI", True, "TRAJ35", "homosapiens", "CYFCAELNAGNNRKLIC"), # J gene with conserved C

        ),
    )
    def test_trimming_and_j_symbol(self, seq, trimming, j_symbol, species, expected):
        result = junction.standardize(seq=seq, trimming=trimming, j_symbol=j_symbol, species=species)

        assert result == expected