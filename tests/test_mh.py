import pytest
from tidytcells import mh
from tidytcells._resources import *


class TestStandardize:
    @pytest.mark.parametrize("species", ("foobar", "yoinkdoink", ""))
    def test_unsupported_species(self, species, caplog):
        result = mh.standardize(symbol="HLA-A*01:01:01:01", species=species)
        assert "Unsupported" in caplog.text
        assert result.original_input == "HLA-A*01:01:01:01"
        assert not result.is_standardized
        assert result.symbol is None
        assert str(result) == ""

    @pytest.mark.parametrize("symbol", (1234, None))
    def test_bad_type(self, symbol):
        with pytest.raises(TypeError):
            mh.standardize(symbol=symbol)

    def test_default_homosapiens(self):
        result = mh.standardize("HLA-B*07")
        assert result.symbol == "HLA-B*07"
        assert result.allele == "HLA-B*07"
        assert result.protein == "HLA-B*07"
        assert result.gene == "HLA-B"
        assert result.is_standardized

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("HLA-B8", "HLA-B*08"),
            ("A1", "HLA-A*01"),
            ("H-2Eb1", "MH2-EB1"),
            ("H-2Aa", "MH2-AA"),
        ),
    )
    def test_any_species(self, symbol, expected):
        result = mh.standardize(symbol, species="any", database="IMGT")
        assert result.symbol == expected
        assert str(result) == expected
        assert result.is_standardized

    @pytest.mark.parametrize(
        ("symbol", "expected", "precision_level"),
        (
            ("HLA-DRB3*01:01:02:01", "HLA-DRB3*01:01:02:01", "allele"),
            ("HLA-DRB3*01:01:02:01", "HLA-DRB3*01:01", "protein"),
            ("HLA-DRB3*01:01:02:01", "HLA-DRB3", "gene"),
        ),
    )
    def test_precision(self, symbol, expected, precision_level):
        result = mh.standardize(
            symbol=symbol, species="homosapiens"
        )

        assert result.__getattribute__(precision_level) == expected

    def test_standardise(self):
        result = mh.standardise("HLA-B*07")

        assert result.symbol == "HLA-B*07"
        assert result.is_standardized
        assert result.error is None

    def test_log_failures(self, caplog):
        mh.standardize("foobarbaz", log_failures=False)
        assert len(caplog.records) == 0

    @pytest.mark.parametrize("database", ("foobar", "imgt", ""))
    def test_unsupported_database(self, database):
        with pytest.raises(ValueError):
            mh.standardize("HLA-A*02:01", database=database)

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("HLA-A*0201", "HLA-A*02:01"),
            ("HLA-Cw*0301", "HLA-C*03:04"),
            ("H-2Kb", "H2-Kb"),
            ("I-Ed", "H2-IEd"),
        ),
    )
    def test_any_species_mro(self, symbol, expected):
        result = mh.standardize(symbol, species="any", database="MRO")
        assert result.symbol == expected
        assert result.is_standardized


class TestStandardizeHomoSapiens:
    @pytest.mark.parametrize("symbol", [*VALID_HOMOSAPIENS_MH, "B2M"])
    def test_already_correctly_formatted(self, symbol):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="IMGT")

        assert result.symbol == symbol

    @pytest.mark.parametrize(
        "symbol", ("foobar", "yoinkdoink", "HLA-FOOBAR123456", "=======")
    )
    def test_invalid_mh(self, symbol, caplog):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="IMGT")
        assert "Failed to standardize" in caplog.text
        assert result.symbol is None
        assert result.error is not None
        assert not result.is_standardized

    @pytest.mark.parametrize("symbol", ("HLA-A*01:01:1:1:1:1:1:1",))
    def test_bad_allele_designation(self, symbol, caplog):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="IMGT")
        assert "Failed to standardize" in caplog.text
        assert result.symbol is None
        assert result.error is not None
        assert not result.is_standardized


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
        result = mh.standardize(symbol=symbol, species="homosapiens", database="IMGT")

        assert result.is_standardized
        assert result.error is None
        assert result.symbol == expected

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
        result = mh.standardize(symbol=symbol, species="homosapiens", database="IMGT")

        assert result.is_standardized
        assert result.error is None
        assert result.symbol == "HLA-A*01:01:01:01"

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
        result = mh.standardize(symbol=symbol, species="homosapiens", database="IMGT")

        assert result.symbol == expected


class TestStandardizeHomoSapiensMro:
    @pytest.mark.parametrize("symbol", VALID_HOMOSAPIENS_MH_MRO)
    def test_already_correctly_formatted(self, symbol):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="MRO")

        assert result.symbol == symbol
        assert result.is_standardized

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("HLA-A*0201", "HLA-A*02:01"),
            ("HLA-A0201", "HLA-A*02:01"),
            ("HLA-A02:01", "HLA-A*02:01"),
            ("HLA-Cw*0301", "HLA-C*03:04"),
            ("HLA-Cw*0601", "HLA-C*06:02"),
            ("hla-cw*0701", "HLA-C*07:01"),
        ),
    )
    def test_allele_synonyms(self, symbol, expected):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="MRO")

        assert result.is_standardized
        assert result.symbol == expected

    @pytest.mark.parametrize("symbol", ("HLA-A*01:01:1:1:1:1:1:1",))
    def test_bad_allele_designation(self, symbol, caplog):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="MRO")
        assert "Failed to standardize" in caplog.text
        assert result.symbol is None
        assert result.error is not None
        assert not result.is_standardized

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("HLA-DQA*01:01", "HLA-DQA1*01:01"),
            ("HLA-DQB*05:01", "HLA-DQB1*05:01"),
            ("HLA-DRA1*01:01", "HLA-DRA*01:01"),
        ),
    )
    def test_fix_deprecated_names(self, symbol, expected):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="MRO")

        assert result.is_standardized
        assert result.error is None
        assert result.symbol == expected

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
        result = mh.standardize(symbol=symbol, species="homosapiens", database="MRO")

        assert result.is_standardized
        assert result.error is None
        assert result.symbol == "HLA-A*01:01:01:01"

    def test_allele_synonyms_not_used_for_imgt(self, caplog):
        result = mh.standardize(symbol="HLA-Cw*0301", species="homosapiens", database="IMGT")

        assert not result.is_standardized

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("HLA-B8", "HLA-B*08"),
            ("A*01:01", "HLA-A*01:01"),
            ("HLA-B*5701", "HLA-B*57:01"),
            ("HLA-DQA1*0501", "HLA-DQA1*05:01"),
            ("HLA-DQB103:01", "HLA-DQB1*03:01"),
            ("HLA-DRB*07:01", "HLA-DRB1*07:01"),
            ("HLA-DQA*01:01", "HLA-DQA1*01:01"),
            ("HLA-A*3351", "HLA-A*33:51"),
        ),
    )
    def test_various_typos(self, symbol, expected):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="MRO")

        assert result.symbol == expected

    @pytest.mark.parametrize("symbol", ("HLA-DMA", "HLA-MICA*001", "HLA-TAP1*01:01", "D6S204", "foobar"))
    def test_not_in_mro(self, symbol, caplog):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="MRO")

        assert "Failed to standardize" in caplog.text
        assert result.symbol is None
        assert not result.is_standardized

    @pytest.mark.parametrize(
        ("symbol", "expected", "expected_mutation"),
        (
            ("HLA-Cw*0301 W167A mutant", "HLA-C*03:04 W167A mutant", "W167A"),
            ("HLA-A*02:01 K66A MUTANT", "HLA-A*02:01 K66A mutant", "K66A"),
        ),
    )
    def test_mutant_with_synonym(self, symbol, expected, expected_mutation):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="MRO")

        assert result.is_standardized
        assert result.symbol == expected
        assert result.mutation == expected_mutation


class TestStandardizeMusMusculus:
    @pytest.mark.parametrize("symbol", VALID_MUSMUSCULUS_MH)
    def test_already_correctly_formatted(self, symbol):
        result = mh.standardize(symbol=symbol, species="musmusculus", database="IMGT")

        assert result.symbol == symbol
        assert result.is_standardized
        assert result.error is None

    @pytest.mark.parametrize("symbol", ("foobar", "yoinkdoink", "MH1-ABC", "======="))
    def test_invalid_mh(self, symbol, caplog):
        result = mh.standardize(symbol=symbol, species="musmusculus", database="IMGT")
        assert "Failed to standardize" in caplog.text
        assert result.symbol is None
        assert not result.is_standardized
        assert result.error is not None

    @pytest.mark.parametrize(
        ("symbol", "expected"), (("H-2Eb1", "MH2-EB1"), ("H-2Aa", "MH2-AA"))
    )
    def test_fix_deprecated_names(self, symbol, expected):
        result = mh.standardize(symbol=symbol, species="musmusculus", database="IMGT")

        assert result.symbol == expected


class TestStandardizeMusMusculusMro:
    @pytest.mark.parametrize("symbol", VALID_MUSMUSCULUS_MH_MRO)
    def test_already_correctly_formatted(self, symbol):
        result = mh.standardize(symbol=symbol, species="musmusculus", database="MRO")

        assert result.symbol == symbol
        assert result.gene == symbol
        assert result.allele is None
        assert result.is_standardized
        assert result.error is None
        assert result.species == "musmusculus"

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("H-2Kb", "H2-Kb"),
            ("H-2-Kb", "H2-Kb"),
            ("H2Kb", "H2-Kb"),
            ("h2-kb", "H2-Kb"),
            ("Kb", "H2-Kb"),
            ("H-2Db", "H2-Db"),
            ("I-Ab", "H2-IAb"),
            ("IAg7", "H2-IAg7"),
            ("H-2 IEk", "H2-IEk"),
            ("H-2Lw16", "H2-Lq"),
            ("H2-Lw16", "H2-Lq"),
            ("Qa1b", "H2-Qa-1b"),
        ),
    )
    def test_various_typos_and_synonyms(self, symbol, expected):
        result = mh.standardize(symbol=symbol, species="musmusculus", database="MRO")

        assert result.is_standardized
        assert result.symbol == expected

    @pytest.mark.parametrize("symbol", ("foobar", "DR", "KB", "H2-Kz", "MH1-M5", "======="))
    def test_invalid_mh(self, symbol, caplog):
        result = mh.standardize(symbol=symbol, species="musmusculus", database="MRO")

        assert "Failed to standardize" in caplog.text
        assert result.symbol is None
        assert not result.is_standardized
        assert result.error is not None

    def test_imgt_names_not_used_for_mro(self, caplog):
        result = mh.standardize(symbol="H2-Kb", species="musmusculus", database="IMGT")

        assert not result.is_standardized


class TestQuery:
    @pytest.mark.filterwarnings("ignore:tidytcells is not.+aware")
    @pytest.mark.parametrize(
        ("species", "precision", "expected_len", "expected_in", "expected_not_in"),
        (
            ("homosapiens", "allele", 26961, "HLA-DRB3*03:04", "HLA-DRB3*03:04P"),
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

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("HLA-A*3351", "HLA-A*33:51"),
            ("HLA-A3058", "HLA-A*30:58"),
            ("HLA-B*8101", "HLA-B*81:01"),
            ("HLA-A*02101", "HLA-A*02:101"),
            ("HLA-B*390101", "HLA-B*39:01:01"),
            ("HLA-A*03010101", "HLA-A*03:01:01:01"),
        ),
    )
    def test_missing_colons_ambiguous_split(self, symbol, expected):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="IMGT")

        assert result.is_standardized
        assert result.symbol == expected

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("HLA-A*3351", "HLA-A*33:51"),
            ("HLA-A3058", "HLA-A*30:58"),
            ("HLA-B*8101", "HLA-B*81:01"),
            ("HLA-A*02101", "HLA-A*02:101"),
            ("HLA-B*390101", "HLA-B*39:01"),
            ("HLA-A*03010101", "HLA-A*03:01"),
        ),
    )
    def test_missing_colons_ambiguous_split_mro(self, symbol, expected):
        result = mh.standardize(symbol=symbol, species="homosapiens", database="MRO")

        assert result.is_standardized
        assert result.symbol == expected

    @pytest.mark.parametrize("database", ("IMGT", "MRO"))
    @pytest.mark.parametrize(
        ("symbol", "expected", "expected_allele", "expected_mutation"),
        (
                ("HLA-A*02:01 K66A mutant", "HLA-A*02:01 K66A mutant", "HLA-A*02:01", "K66A"),
                ("HLA-A0201 K66A, E63Q mutant", "HLA-A*02:01 K66A, E63Q mutant", "HLA-A*02:01", "K66A, E63Q"),
                ("HLA-B*08:01 B:I66A Mutant", "HLA-B*08:01 B:I66A mutant", "HLA-B*08:01", "B:I66A"),
                ("HLA-DRB1*01:01 G86Y mutant", "HLA-DRB1*01:01 G86Y mutant", "HLA-DRB1*01:01", "G86Y"),
        ),
    )
    def test_mutant(self, symbol, expected, expected_allele, expected_mutation, database):
        result = mh.standardize(symbol=symbol, species="homosapiens", database=database)

        assert result.is_standardized
        assert result.symbol == expected
        assert str(result) == expected
        assert result.allele == expected_allele
        assert result.protein == expected_allele
        assert result.gene == expected_allele.split("*")[0]
        assert result.mutation == expected_mutation

    def test_mutant_invalid_allele(self):
        result = mh.standardize(symbol="HLA-FOO*01:01 K66A mutant", species="homosapiens")

        assert not result.is_standardized
        assert result.symbol is None
        assert result.mutation is None
        assert result.attempted_fix == "HLA-FOO*01:01 K66A mutant"

    def test_no_mutation(self):
        result = mh.standardize(symbol="HLA-A*02:01", species="homosapiens")

        assert result.mutation is None