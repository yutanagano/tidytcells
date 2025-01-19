import pytest
from tidytcells import mh
from tidytcells._resources import *


class TestStandardize:
    @pytest.mark.parametrize("species", ("foobar", "yoinkdoink", ""))
    def test_unsupported_species(self, species, caplog):
        result = mh.standardize(symbol="HLA-A*01:01:01:01", species=species)
        assert "Unsupported" in caplog.text
        assert result == "HLA-A*01:01:01:01"

    @pytest.mark.parametrize("symbol", (1234, None))
    def test_bad_type(self, symbol):
        with pytest.raises(TypeError):
            mh.standardize(symbol=symbol)

    def test_default_homosapiens(self):
        result = mh.standardize("HLA-B*07")

        assert result == "HLA-B*07"

    @pytest.mark.parametrize(
        ("symbol", "expected", "precision"),
        (
            ("HLA-DRB3*01:01:02:01", "HLA-DRB3*01:01:02:01", "allele"),
            ("HLA-DRB3*01:01:02:01", "HLA-DRB3*01:01", "protein"),
            ("HLA-DRB3*01:01:02:01", "HLA-DRB3", "gene"),
        ),
    )
    def test_precision(self, symbol, expected, precision):
        result = mh.standardize(
            symbol=symbol, species="homosapiens", precision=precision
        )

        assert result == expected

    def test_standardise(self):
        result = mh.standardise("HLA-B*07")

        assert result == "HLA-B*07"

    def test_log_failures(self, caplog):
        mh.standardize("foobarbaz", log_failures=False)
        assert len(caplog.records) == 0

    def test_on_fail(self, caplog):
        result = mh.standardize("foobarbaz", on_fail="keep")
        assert "Failed to standardize" in caplog.text
        assert result == "foobarbaz"


class TestStandardizeHomoSapiens:
    @pytest.mark.parametrize("symbol", [*VALID_HOMOSAPIENS_MH, "B2M"])
    def test_already_correctly_formatted(self, symbol):
        result = mh.standardize(symbol=symbol, species="homosapiens")

        assert result == symbol

    @pytest.mark.parametrize(
        "symbol", ("foobar", "yoinkdoink", "HLA-FOOBAR123456", "=======")
    )
    def test_invalid_mh(self, symbol, caplog):
        result = mh.standardize(symbol=symbol, species="homosapiens")
        assert "Failed to standardize" in caplog.text
        assert result == None

    @pytest.mark.parametrize("symbol", ("HLA-A*01:01:1:1:1:1:1:1",))
    def test_bad_allele_designation(self, symbol, caplog):
        result = mh.standardize(symbol=symbol, species="homosapiens")
        assert "Failed to standardize" in caplog.text
        assert result == None

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("D6S204", "HLA-C"),
            ("HLA-DQA*01:01", "HLA-DQA1*01:01"),
            ("HLA-DQB*05:01", "HLA-DQB1*05:01"),
            ("HLA-DRA1*01:01", "HLA-DRA*01:01"),
        ),
    )
    def test_fix_deprecated_names(self, symbol, expected):
        result = mh.standardize(symbol=symbol, species="homosapiens")

        assert result == expected

    @pytest.mark.parametrize(
        "symbol",
        (
            "HLA-A*01:01:01:01N",
            "HLA-A*01:01:01:01L",
            "HLA-A*01:01:01:01S",
            "HLA-A*01:01:01:01C",
            "HLA-A*01:01:01:01A",
            "HLA-A*01:01:01:01Q",
        ),
    )
    def test_remove_expression_qualifier(self, symbol):
        result = mh.standardize(symbol=symbol, species="homosapiens")

        assert result == "HLA-A*01:01:01:01"

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("HLA-B8", "HLA-B*08"),
            ("A*01:01", "HLA-A*01:01"),
            ("A1", "HLA-A*01"),
            ("HLA-B*5701", "HLA-B*57:01"),
            ("HLA-DQA1*0501", "HLA-DQA1*05:01"),
            ("B35.3", "HLA-B*35:03"),
            ("HLA-DQB103:01", "HLA-DQB1*03:01"),
            ("HLA-A*01:01:1:1", "HLA-A*01:01:01:01"),
            ("HLA-DRB*07:01", "HLA-DRB1*07:01"),
        ),
    )
    def test_various_typos(self, symbol, expected):
        result = mh.standardize(symbol=symbol, species="homosapiens")

        assert result == expected


class TestStandardizeMusMusculus:
    @pytest.mark.parametrize("symbol", VALID_MUSMUSCULUS_MH)
    def test_already_correctly_formatted(self, symbol):
        result = mh.standardize(symbol=symbol, species="musmusculus")

        assert result == symbol

    @pytest.mark.parametrize("symbol", ("foobar", "yoinkdoink", "MH1-ABC", "======="))
    def test_invalid_mh(self, symbol, caplog):
        result = mh.standardize(symbol=symbol, species="musmusculus")
        assert "Failed to standardize" in caplog.text
        assert result == None

    @pytest.mark.parametrize(
        ("symbol", "expected"), (("H-2Eb1", "MH2-EB1"), ("H-2Aa", "MH2-AA"))
    )
    def test_fix_deprecated_names(self, symbol, expected):
        result = mh.standardize(symbol=symbol, species="musmusculus")

        assert result == expected


class TestQuery:
    @pytest.mark.filterwarnings("ignore:tidytcells is not.+aware")
    @pytest.mark.parametrize(
        ("species", "precision", "expected_len", "expected_in", "expected_not_in"),
        (
            ("homosapiens", "allele", 25194, "HLA-DRB3*03:04", "HLA-DRB3*03:04P"),
            ("homosapiens", "gene", 46, "HLA-B", "HLA-FOO"),
            ("musmusculus", "allele", 70, "MH1-M10-1", "HLA-A"),
            ("musmusculus", "gene", 70, "MH1-Q8", "H2-Aa"),
        ),
    )
    def test_query_all(
        self, species, precision, expected_len, expected_in, expected_not_in
    ):
        result = mh.query(species=species, precision=precision)

        assert type(result) == frozenset
        assert len(result) == expected_len
        assert expected_in in result
        assert not expected_not_in in result

    @pytest.mark.parametrize(
        (
            "species",
            "precision",
            "contains",
            "expected_len",
            "expected_in",
            "expected_not_in",
        ),
        (
            ("homosapiens", "gene", "DR", 10, "HLA-DRA", "HLA-A"),
            ("musmusculus", "gene", "T", 24, "MH1-T10", "MH1-Q10"),
        ),
    )
    def test_query_contains(
        self, species, precision, contains, expected_len, expected_in, expected_not_in
    ):
        result = mh.query(
            species=species, precision=precision, contains_pattern=contains
        )

        assert len(result) == expected_len
        assert expected_in in result
        assert not expected_not_in in result

    def test_query_default_species(self):
        result = mh.query(precision="gene", contains_pattern="DR")

        assert len(result) == 10
        assert "HLA-DRA" in result


class TestGetChain:
    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("HLA-A", "alpha"),
            ("HLA-B", "alpha"),
            ("HLA-C", "alpha"),
            ("HLA-DPA1", "alpha"),
            ("HLA-DQA2", "alpha"),
            ("HLA-DRA", "alpha"),
            ("HLA-E", "alpha"),
            ("HLA-F", "alpha"),
            ("HLA-G", "alpha"),
            ("HLA-DPB2", "beta"),
            ("HLA-DQB1", "beta"),
            ("HLA-DRB3", "beta"),
            ("B2M", "beta"),
        ),
    )
    def test_get_chain(self, symbol, expected):
        result = mh.get_chain(symbol=symbol)

        assert result == expected

    @pytest.mark.parametrize("symbol", ("foo", "HLA", "0"))
    def test_unrecognised_gene_names(self, symbol, caplog):
        result = mh.get_chain(symbol=symbol)
        assert "Unrecognized gene" in caplog.text
        assert result == None

    @pytest.mark.parametrize("symbol", (1234, None))
    def test_bad_type(self, symbol):
        with pytest.raises(TypeError):
            mh.get_chain(symbol)

    def test_log_failures(self, caplog):
        mh.get_chain("foobarbaz", log_failures=False)
        assert len(caplog.records) == 0


class TestGetClass:
    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("HLA-A", 1),
            ("HLA-B", 1),
            ("HLA-C", 1),
            ("HLA-DPA1", 2),
            ("HLA-DQA2", 2),
            ("HLA-DRA", 2),
            ("HLA-DPB2", 2),
            ("HLA-DQB1", 2),
            ("HLA-DRB3", 2),
            ("HLA-E", 1),
            ("HLA-F", 1),
            ("HLA-G", 1),
            ("B2M", 1),
        ),
    )
    def test_get_class(self, symbol, expected):
        result = mh.get_class(symbol=symbol)

        assert result == expected

    @pytest.mark.parametrize("symbol", ("foo", "HLA", "0"))
    def test_unrecognised_gene_names(self, symbol, caplog):
        result = mh.get_class(symbol=symbol)
        assert "Unrecognized gene" in caplog.text
        assert result == None

    @pytest.mark.parametrize("symbol", (1234, None))
    def test_bad_type(self, symbol):
        with pytest.raises(TypeError):
            mh.get_class(symbol)

    def test_log_failures(self, caplog):
        mh.get_class("foobarbaz", log_failures=False)
        assert len(caplog.records) == 0
