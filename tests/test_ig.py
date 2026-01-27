import pytest
from tidytcells import ig
from tidytcells._resources import VALID_HOMOSAPIENS_IG, VALID_MUSMUSCULUS_IG


class TestStandardize:
    @pytest.mark.parametrize("species", ("foobar", "yoinkdoink", ""))
    def test_unsupported_species(self, species, caplog):
        result = ig.standardize(symbol="foobarbaz", species=species)
        assert "Unsupported" in caplog.text
        assert "Unsupported" in result.error
        assert result.original_input == "foobarbaz"
        assert result.symbol is None
        assert not result.is_standardized
        assert result.species is None

    @pytest.mark.parametrize("symbol", (1234, None))
    def test_bad_type(self, symbol):
        with pytest.raises(TypeError):
            ig.standardize(symbol=symbol)

    def test_default_homosapiens(self):
        result = ig.standardize("IGHV1/OR15-1*01")
        assert result.is_standardized
        assert result.error is None
        assert result.symbol == "IGHV1/OR15-1*01"
        assert result.species == "homosapiens"

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("IGHV1/OR15-1*01", "IGHV1/OR15-1*01"),
            ("IGHV6-1", "IGHV6-1"),
            ("LV2-18", "IGLV2-18"),
        )
    )
    def test_any_species(self, symbol, expected):
        result = ig.standardize(symbol, species="any")
        assert result.is_standardized
        assert result.error is None
        assert result.symbol == expected

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (("IGLJ1*01&nbsp;", "IGLJ1*01"), (" IGLJ1 * 01 ", "IGLJ1*01")),
    )
    def test_remove_pollutants(self, symbol, expected):
        result = ig.standardize(symbol=symbol, species="homosapiens")

        assert result.is_standardized
        assert result.error is None
        assert result.symbol == expected
        assert result.allele == expected

    @pytest.mark.filterwarnings("ignore:Failed to standardize")
    @pytest.mark.parametrize(
        ("symbol", "expected", "enforce_functional"),
        (
            ("IGHV6-1*01", "IGHV6-1*01", True),
            ("IGHV6-1*03", None, True),
            ("IGHV6-1*03", "IGHV6-1*03", False),
            ("IGLV6-57", "IGLV6-57", True),
            ("IGLV7-35", None, True),
            ("IGLV7-35", "IGLV7-35", False),
        ),
    )
    def test_enforce_functional(self, symbol, expected, enforce_functional):
        result = ig.standardize(
            symbol=symbol, species="homosapiens", enforce_functional=enforce_functional
        )

        assert result.symbol == expected
        assert result.allele == expected or result.gene == expected
        assert result.species == "homosapiens"
        if result.is_standardized:
            assert result.gene_type == "V"


    @pytest.mark.filterwarnings("ignore:Failed to standardize")
    @pytest.mark.parametrize(
        ("symbol", "expected", "allow_subgroup"),
        (
            ("IGLV6-57", "IGLV6-57", True),
            ("IGLV6-57", "IGLV6-57", False),
            ("IGLV6", "IGLV6", True),
            ("IGLV6", None, False),
        ),
    )
    def test_allow_subgroup(self, symbol, expected, allow_subgroup):
        result = ig.standardize(
            symbol=symbol, species="homosapiens", allow_subgroup=allow_subgroup,
        )

        assert result.symbol == expected
        assert result.allele is None
        assert result.gene == expected or result.subgroup == expected


    @pytest.mark.filterwarnings("ignore:Failed to standardize")
    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("IGHV5-10-1*01", "IGHV5-10-1*01"),
            ("IGHV5-10-1", "IGHV5-10-1"),
            ("IGHJ1-1", "IGHJ1"),
            ("IGHJ1-1*01", "IGHJ1*01"),
            ("IGHJ4-1*02", "IGHJ4*02"),
            ("IGLV(VI)-22-1-1", "IGLV(VI)-22-1"),
        ),
    )
    def test_remove_illegal_dash1(self, symbol, expected, ):
        result = ig.standardize(
            symbol=symbol, species="homosapiens",
        )

        assert result.is_standardized
        assert result.symbol == expected

    @pytest.mark.parametrize(
        ("symbol", "allow_subgroup", "expected_allele", "expected_gene", "expected_subgroup", "expected_symbol"),
        (
            ("IGLV7-43*01", True, "IGLV7-43*01", "IGLV7-43", "IGLV7", "IGLV7-43*01"),
            ("IGLV7-43*01", False, "IGLV7-43*01", "IGLV7-43", "IGLV7", "IGLV7-43*01"),
            ("IGLV8-61", True, None, "IGLV8-61", "IGLV8", "IGLV8-61"),
            ("IGLV8-61", False, None, "IGLV8-61", "IGLV8", "IGLV8-61"),
            ("IGLV8", True, None, None, "IGLV8", "IGLV8"),
            ("IGLV8", False, None, None, None, None),
        ),
    )
    def test_precision(self, symbol, allow_subgroup, expected_allele, expected_gene, expected_subgroup, expected_symbol):
        result = ig.standardize(
            symbol=symbol, species="homosapiens", allow_subgroup=allow_subgroup,
        )

        assert result.allele == expected_allele
        assert result.gene == expected_gene
        assert result.subgroup == expected_subgroup
        assert result.symbol == expected_symbol
        if result.is_standardized:
            assert str(result) == expected_symbol
        else:
            assert str(result) == ""

    def test_standardise(self):
        result = ig.standardise("IGLV8/OR8-1*02")

        assert result.symbol == "IGLV8/OR8-1*02"

    def test_log_failures(self, caplog):
        ig.standardize("foobarbaz", log_failures=False)
        assert len(caplog.records) == 0


class TestStandardizeHomoSapiens:
    @pytest.mark.parametrize("symbol", VALID_HOMOSAPIENS_IG)
    def test_already_correctly_formatted(self, symbol):
        result = ig.standardize(symbol=symbol, species="homosapiens")

        assert result.symbol == symbol

    @pytest.mark.parametrize("symbol", ("foobar", "IGHD4-3*01"))
    def test_invalid_ig(self, symbol, caplog):
        result = ig.standardize(symbol=symbol, species="homosapiens")
        assert "Failed to standardize" in caplog.text
        assert result.symbol == None

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            (
                "A10",
                "IGKV6D-21",
            ),
            ("IGKV1OR2108", "IGKV1/OR2-108"),
            ("IGO1", "IGKV1/OR2-108"),
        ),
    )
    def test_resolve_alternate_ig_names(self, symbol, expected):
        result = ig.standardize(symbol=symbol, species="homosapiens")

        assert result.symbol == expected

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            (" IGKJ2*01", "IGKJ2*01"),
            ("IGLJ2-1", "IGLJ2"),
            ("IGHD01-01*01", "IGHD1-1*01"),
            ("IGKV19", "IGKV1-9"),
            ("IGHV1OR15-1*01", "IGHV1/OR15-1*01"),
            ("LV2-18", "IGLV2-18"),
            ("IGHV1/OR15-1*01", "IGHV1/OR15-1*01"),
            ("IGLV5.37", "IGLV5-37"),
        ),
    )
    def test_various_typos(self, symbol, expected):
        result = ig.standardize(symbol=symbol, species="homosapiens")

        assert result.symbol == expected
        assert result.species == "homosapiens"

class TestStandardizeMusMusculus:

    @pytest.mark.parametrize("symbol", VALID_MUSMUSCULUS_IG)
    def test_already_correctly_formatted(self, symbol):
        result = ig.standardize(symbol=symbol, species="musmusculus")

        assert result.symbol == symbol

    @pytest.mark.parametrize("symbol", ("foobar", "noice"))
    def test_invalid_ig(self, symbol, caplog):
        result = ig.standardize(symbol=symbol, species="musmusculus")
        assert "Failed to standardize" in caplog.text
        assert not result.is_standardized
        assert result.symbol is None
        assert result.species == "musmusculus"
        assert str(result) == ""

class TestQuery:
    @pytest.mark.parametrize(
        ("species", "precision", "expected_len", "expected_in", "expected_not_in"),
        (
            ("homosapiens", "allele", 1323, "IGKV3-15*02", "IGKV3-15"),
            ("homosapiens", "gene", 496, "IGKV3-15", "IGKV3-15*02"),
        ),
    )
    def test_query_all(
        self, species, precision, expected_len, expected_in, expected_not_in
    ):
        result = ig.query(species=species, precision=precision)

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
            ("homosapiens", "gene", "HJ", 9, "IGHJ1", "IGLJ1-6"),
            ("homosapiens", "gene", "LC", 10, "IGLC1", "IGHV1"),
        ),
    )
    def test_query_contains(
        self, species, precision, contains, expected_len, expected_in, expected_not_in
    ):
        result = ig.query(
            species=species, precision=precision, contains_pattern=contains
        )

        assert len(result) == expected_len
        assert expected_in in result
        assert not expected_not_in in result

    @pytest.mark.parametrize(
        (
            "species",
            "precision",
            "functionality",
            "expected_len",
            "expected_in",
            "expected_not_in",
        ),
        (
            ("homosapiens", "gene", "F", 190, "IGHD6-19", "IGHV1/OR16-21*01"),
            ("homosapiens", "allele", "NF", 630, "IGHV1/OR16-21*01", "IGL35*01"),
            ("homosapiens", "gene", "P", 283, "IGKV3D-31", "IGHV3-20*02"),
            ("homosapiens", "allele", "ORF", 101, "IGHV3-20*02", "IGKV3D-31"),
        ),
    )
    def test_query_functionality(
        self,
        species,
        precision,
        functionality,
        expected_len,
        expected_in,
        expected_not_in,
    ):
        result = ig.query(
            species=species, precision=precision, functionality=functionality
        )

        assert len(result) == expected_len
        assert expected_in in result
        assert not expected_not_in in result

    def test_query_default_species(self):
        result = ig.query(precision="gene", contains_pattern="HJ")
        assert len(result) == 9
        assert "IGHJ3" in result


class TestGetAaSequence:
    @pytest.mark.parametrize(
        ("symbol", "species", "expected"),
        (
            (
                "IGHV1-18*02",
                "homosapiens",
                {
                    "FR1-IMGT": "QVQLVQSGAEVKKPGASVKVSCKAS",
                    "CDR1-IMGT": "GYTFTSYG",
                    "FR2-IMGT": "ISWVRQAPGQGLEWMGW",
                    "CDR2-IMGT": "ISAYNGNT",
                    "FR3-IMGT": "NYAQKLQGRVTMTTDTSTSTAYMELRSLRSDDTA",
                    "V-REGION": "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYGISWVRQAPGQGLEWMGWISAYNGNTNYAQKLQGRVTMTTDTSTSTAYMELRSLRSDDTA",
                    'functionality': 'F',
                },
            ),
            (
                "IGHD1-1*01",
                "homosapiens",
                {
                    "D-REGION": "GTTGT",
                    'functionality': 'F',
                },
            ),
            (
                "IGHJ1*01",
                "homosapiens",
                {
                    'J-MOTIF': 'WGQG',
                    "J-REGION": "AEYFQHWGQGTLVTVSS",
                    'functionality': 'F',
                },
            )
        ),
    )
    def test_get_aa_sequence(self, symbol, species, expected):
        result = ig.get_aa_sequence(symbol=symbol, species=species)

        assert result == expected
