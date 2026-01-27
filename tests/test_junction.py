import pytest
from tidytcells import junction


class Teststandardize:
    @pytest.mark.parametrize(
        ("seq", "locus"),
        (
                ("CASSPGGADRRIDGYTF", "TR"),
                ("CASSLMPGQGSYEQYF", "TR"),
                ("CSARDPGDDKPQHF", "TR"),
                ("CAALRATGGNNKLTF", "TR"),
                ("CVVAVSNFGNEKLTF", "TR"),
                ("CAASASFGDNSKLIW", "IG"),
        ),
    )
    def test_already_correct(self, seq, locus):
        result = junction.standardize(seq=seq, locus=locus)

        assert result.junction == seq

    @pytest.mark.parametrize("seq", ("123456", "ASDFGHJKL", "A?AAAA", "AAAXAA"))
    def test_various_rejections(self, seq, caplog):
        result = junction.standardize(seq=seq, locus="TR")
        assert "Not a valid amino acid sequence" in caplog.text
        assert "Not a valid amino acid sequence" in result.error
        assert result.junction is None
        assert not result.is_standardized

    @pytest.mark.parametrize(
        ("seq", "expected"),
        (
            ("casqyf", "CASQYF"),
            ("ASQY", "CASQYF"),
            ("ASQGELF", "CASQGELFF"),
            ("CASQY", "CASQYF"),
            ("ASQYF", "CASQYF"),
        ),
    )
    def test_various_corrections(self, seq, expected):
        result = junction.standardize(seq=seq, locus="TR")

        assert result.junction == expected

    @pytest.mark.parametrize("seq", (None, 1, True, 5.4))
    def test_bad_input_type(self, seq):
        with pytest.raises(TypeError):
            junction.standardize(seq=seq, locus="TR")

    def test_standardize(self):
        result = junction.standardize(seq="CASSPGGADRRIDGYTF", locus="TR")

        assert result.junction == "CASSPGGADRRIDGYTF"
        assert result.is_standardized
        assert result.error is None

    def test_log_failures(self, caplog):
        junction.standardize(seq="123456", log_failures=False, locus="TR")
        assert len(caplog.records) == 0

    def test_failed_input(self, caplog):
        result = junction.standardize("foobarbaz", locus="TR")
        assert "Failed to standardize foobarbaz" in caplog.text
        assert "Not a valid amino acid sequence" in result.error
        assert result.junction == None
        assert result.cdr3 == None
        assert result.original_input == "foobarbaz"

    def test_bad_locus(self, caplog):
        with pytest.raises(ValueError):
            junction.standardize("AELNAGNNRKLI", locus="TRB", j_symbol="TRA")

    def test_j_symbol_fail(self, caplog):
        result = junction.standardize(
            "AAAAAA", j_symbol="TRAJ", locus="TR"
        )

        assert "J side reconstruction unsuccessful" in caplog.text
        assert "J side reconstruction unsuccessful" in result.error
        assert result.original_input == "AAAAAA"
        assert result.junction == None
        assert not result.is_standardized

    def test_j_symbol_fail_locus(self):
        with pytest.raises(
            ValueError, match='"j_symbol" IGH is not a valid J gene for "locus" TR'
        ):
            junction.standardize(
                "AAAAAA", j_symbol="IGH", species="musmusculus", locus="TR"
            )

    def test_j_symbol_fail_species(self,caplog):
        result = junction.standardize(
            "AAAAAA", j_symbol="IGHJ", species="nonexistent", locus="IG"
        )
        assert "Unsupported" in caplog.text
        assert "Unsupported" in result.error
        assert result.junction is None
        assert not result.is_standardized


    @pytest.mark.parametrize(
        ("seq", "j_symbol", "v_symbol", "locus", "species", "expected"),
        (
                ("CASSPGVFGANVLTFG", None, None, "TR", "homosapiens", "CASSPGVFGANVLTF"),
                ("AELNAGNNR", "TRAJ38*01", None, "TR", "homosapiens", "CAELNAGNNRKLIW"),
                ("CASSPGVFGANVLTL", None, None, "TR", "homosapiens", "CASSPGVFGANVLTF"), # fw mismatch
                ("WASSPGVFGANVLTFG", None, None, "TR", "homosapiens", "CASSPGVFGANVLTF"), # C mismatch
                ("RASSPGVFGANVLTFG", None, None, "TR", "homosapiens", "CASSPGVFGANVLTF"),
                ("SASSPGVFGANVLTFG", None, "TRBV20-1", "TR", "homosapiens", "CSASSPGVFGANVLTF"),
                ("SASSPGVFGANVLTFG", None, None, "TR", "homosapiens", "CASSPGVFGANVLTF"),
                ("DASSPGVFGANVLTFG", None, None, "TR", "homosapiens", None), # do not correct D to C (only W/S/R/G/Y/F/A)
                ("MRESENMDSNYQYVF", None, None, "TRA", "homosapiens", "CAMRESENMDSNYQYVF"),
                ("AELNAGNNRKLI", "TRAJ38*01", None, "TR", "homosapiens", "CAELNAGNNRKLIW"), # for human this is CAE-
                ("AELNAGNNRKLI", "TRAJ38*01", None, "TR", "musmusculus", "CAAELNAGNNRKLIW"), # functional mouse gene best match contributes CAAE-
                ("AELNAGNNRKLI", "TRAJ38*01", "TRAV5-2*01", "TR", "musmusculus", "CAELNAGNNRKLIW"), # specifying a specific allele which contributes CAE- but is pseudogene overrides 'enforce_functional'
                ("YFCAVVFNMDSNYQLIWGLGTSL", "TRAJ38*01", "TRAV", "TR", "homosapiens", "CAVVFNMDSNYQLIW"),
                ("YFCAVVFNMDSNYQLIWGLGTSL", None, None, "TRA", "homosapiens", "CAVVFNMDSNYQLIW"),
                ("YFCAVVFNMDSNYQLIWLLLL", "TRAJ38*01", None, "TR", "homosapiens", None), # J anchor invalid, no output
                ("VYYCARFNMDSNYQYFQHWGQGTLVTVSS", "IGHJ", "IGHV", "IG", "homosapiens", "CARFNMDSNYQYFQHW"),
                ("YFCAELNAGNNVLH", "TRAJ35", None, "TR", "homosapiens", "CAELNAGNNVLHC"),
                ("AGGYQNFYFGKGTMLLVSP", None, None, "TR", "homosapiens", "CAGGYQNFYF"),
                ("VVNRGTGGFKTIFGAG", None, None, "TR", "homosapiens", "CVVNRGTGGFKTIF"),
                ("DSSIYLCSVEATRADTQYFGPGTRLTVL", None, None, "TR", "homosapiens", "CSVEATRADTQYF"),
                ("ASSTRSSGEL", None, None, "TR", "homosapiens", "CASSTRSSGELFF"),
                ("YSTDSSGDIWV", None, None, "IG", "homosapiens", "CYSTDSSGDIWVF"),
                ("QQYGSSPLT", "IGKJ4", "IGKV3", "IG", "homosapiens", "CQQYGSSPLTF"),
                ("CASTGSYGYTFGSGTRLTV", None, None, "TR", "homosapiens", "CASTGSYGYTF"),
                ("ASRPVAG", None, None, "TR", "homosapiens", None),
                ("csadaf", None, None, "TR", "homosapiens", "CSADAF"),
                ("yicsadafg", None, None, "TR", "homosapiens", "CSADAF"),
                ("csadkli", "TRAJ38", None, "TR", "homosapiens", "CSADKLIW"),
                ("csadkli", "TRAJ37", None, "TR", "homosapiens", "CSADKLIF"),
                ("CSYAYVF", None, "IGLV2-11", "IG", "homosapiens", "CCSYAYVF"),
                ("SYAYVF", None, "IGLV2-11", "IG", "homosapiens", "CCSYAYVF"),
                ("CSYAYVF", None, "IGLV2-14", "IG", "homosapiens", "CSYAYVF"),
                ("MRESENMD", None, None, "TRA", "homosapiens", "CAMRESENMDSSYKLIF"),
                ("CAGGSYGLTFGAG", None, None, "IG", "musmusculus", "CAGGSYGLTF"),
                ("ASGDWQGQNT", None, "TRBV13-2*01", "TRB", "musmusculus", "CASGDWQGQNTLYF"),
                ("AIWDYEFHVMDY", "IGHJ4*01", "IGHV1S34*01", "IG", "musmusculus", "CAIWDYEFHVMDYW"),
        )
    )
    def test_various_examples(self, seq, j_symbol, v_symbol, locus, species, expected):
        result = junction.standardize(seq=seq, j_symbol=j_symbol, v_symbol=v_symbol, species=species, locus=locus, allow_c_correction=True, allow_fw_correction=True, allow_v_reconstruction=True, allow_j_reconstruction=True, enforce_functional_j=True)

        assert result.junction == expected
        assert result.species == species

        if expected is not None:
            assert result.cdr3 == expected[1:-1]
            assert result.is_standardized
            assert result.error is None
            assert str(result) == expected
        else:
            assert not result.is_standardized
            assert result.error is not None
            assert len(result.error) > 0
            assert str(result) == ""

