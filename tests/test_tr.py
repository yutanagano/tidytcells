import pytest
from tidytcells import tr
from tidytcells._resources import VALID_HOMOSAPIENS_TR, VALID_MUSMUSCULUS_TR


class TestStandardize:
    @pytest.mark.parametrize("species", ("foobar", "yoinkdoink", ""))
    def test_unsupported_species(self, species, caplog):
        result = tr.standardize(symbol="foobarbaz", species=species)
        assert "Unsupported" in caplog.text
        assert "Unsupported" in result.error
        assert result.original_input == "foobarbaz"
        assert result.symbol is None
        assert not result.is_standardized
        assert result.species is None

    @pytest.mark.parametrize("symbol", (1234, None))
    def test_bad_type(self, symbol):
        with pytest.raises(TypeError):
            tr.standardize(symbol=symbol)

    def test_default_homosapiens(self):
        result = tr.standardize("TRBV20/OR9-2*01")
        assert result.is_standardized
        assert result.error is None
        assert result.symbol == "TRBV20/OR9-2*01"
        assert result.species == "homosapiens"
        assert result.gene_type == "V"


    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("TRBV20/OR9-2*01", "TRBV20/OR9-2*01"),
            ("TCRAV14/4", "TRAV14/DV4"),
            ("TRAV15-1-DV6-1", "TRAV15-1/DV6-1"),
        )
    )
    def test_any_species(self, symbol, expected):
        result = tr.standardize(symbol, species="any")
        assert result.is_standardized
        assert result.error is None
        assert result.symbol == expected

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (("TRAV3*01&nbsp;", "TRAV3*01"), (" TRAV3 * 01 ", "TRAV3*01")),
    )
    def test_remove_pollutants(self, symbol, expected):
        result = tr.standardize(symbol=symbol, species="homosapiens")

        assert result.is_standardized
        assert result.error is None
        assert result.symbol == expected
        assert result.allele == expected

    @pytest.mark.filterwarnings("ignore:Failed to standardize")
    @pytest.mark.parametrize(
        ("symbol", "expected", "enforce_functional"),
        (
            ("TRAV35*01", "TRAV35*01", True),
            ("TRAV35*03", None, True),
            ("TRAV35*03", "TRAV35*03", False),
            ("TRAV35", "TRAV35", True),
            ("TRAV8-7", None, True),
            ("TRAV8-7", "TRAV8-7", False),
        ),
    )
    def test_enforce_functional(self, symbol, expected, enforce_functional):
        result = tr.standardize(
            symbol=symbol, species="homosapiens", enforce_functional=enforce_functional,
        )

        assert result.symbol == expected
        assert result.allele == expected or result.gene == expected
        if result.is_standardized:
            assert result.gene_type == "V"

    @pytest.mark.filterwarnings("ignore:Failed to standardize")
    @pytest.mark.parametrize(
        ("symbol", "expected", "allow_subgroup"),
        (
            ("TRBJ2-7", "TRBJ2-7", True),
            ("TRBJ2-7", "TRBJ2-7", False),
            ("TRBJ2", "TRBJ2", True),
            ("TRBJ2", None, False),
        ),
    )
    def test_allow_subgroup(self, symbol, expected, allow_subgroup):
        result = tr.standardize(
            symbol=symbol, species="homosapiens", allow_subgroup=allow_subgroup,
        )

        assert result.symbol == expected
        assert result.allele is None
        assert result.gene == expected or result.subgroup == expected

    @pytest.mark.filterwarnings("ignore:Failed to standardize")
    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("TRAV14/DV4*01", "TRAV14/DV4*01"),
            ("TRAV14/DV4-1*01", "TRAV14/DV4*01"),
            ("TRAV14-1/DV4-1*01", "TRAV14/DV4*01"),
            ("TRAV14-1/DV4*01", "TRAV14/DV4*01"),
            ("TRAV14-1/DV4*01", "TRAV14/DV4*01"),
            ("TRAV1-1*02", "TRAV1-1*02"),
            ("TRAV10-1", "TRAV10"),
            ("TRAV10-1*01", "TRAV10*01"),
            ("TCRAV29-01", "TRAV29/DV5"),
            ("TCRAV36-01*01", "TRAV36/DV7*01"),
        ),
    )
    def test_remove_illegal_dash1(self, symbol, expected, ):
        result = tr.standardize(
            symbol=symbol, species="homosapiens",
        )

        assert result.is_standardized
        assert result.symbol == expected
        assert result.species == "homosapiens"


    # ("IGLV7-43*01", True, "IGLV7-43*01", "IGLV7-43", "IGLV7", "IGLV7-43*01"),
    # ("IGLV7-43*01", False, "IGLV7-43*01", "IGLV7-43", "IGLV7", "IGLV7-43*01"),
    # ("IGLV8-61", True, None, "IGLV8-61", "IGLV8", "IGLV8-61"),
    # ("IGLV8-61", False, None, "IGLV8-61", "IGLV8", "IGLV8-61"),
    # ("IGLV8", True, None, None, "IGLV8", "IGLV8"),
    # ("IGLV8", False, None, None, None, None),

    @pytest.mark.parametrize(
        ("symbol", "allow_subgroup", "expected_allele", "expected_gene", "expected_subgroup", "expected_symbol"),
        (
            ("TRBV24/OR9-2*01", True, "TRBV24/OR9-2*01", "TRBV24/OR9-2", "TRBV24/OR9", "TRBV24/OR9-2*01"),
            ("TRBV24/OR9-2*01", False, "TRBV24/OR9-2*01", "TRBV24/OR9-2", "TRBV24/OR9", "TRBV24/OR9-2*01"),
            ("TRAV16", True, None, "TRAV16", "TRAV16", "TRAV16"),
            ("TRAV16", False, None, "TRAV16", "TRAV16", "TRAV16"),
            ("TRBV24/OR9", True, None, None, "TRBV24/OR9", "TRBV24/OR9"),
            ("TRBV24/OR9", False, None, None, None, None),
        ),
    )
    def test_precision(self, symbol, allow_subgroup, expected_allele, expected_gene, expected_subgroup, expected_symbol):
        result = tr.standardize(
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
        result = tr.standardise("TRBV20/OR9-2*01")

        assert result.symbol == "TRBV20/OR9-2*01"

    def test_log_failures(self, caplog):
        tr.standardize("foobarbaz", log_failures=False)
        assert len(caplog.records) == 0


class TestStandardizeHomoSapiens:
    @pytest.mark.parametrize("symbol", VALID_HOMOSAPIENS_TR)
    def test_already_correctly_formatted(self, symbol):
        result = tr.standardize(symbol=symbol, species="homosapiens")

        assert result.symbol == symbol

    @pytest.mark.parametrize("symbol", ("foobar", "TRAV3D-3*01"))
    def test_invalid_tr(self, symbol, caplog):
        result = tr.standardize(symbol=symbol, species="homosapiens")
        assert "Failed to standardize" in caplog.text
        assert result.symbol is None
        assert str(result) == ""

    @pytest.mark.parametrize(
        ("symbol", "expected"),
        (
            ("TCRAV32S1", "TRAV25"),
            ("TCRAV14S2", "TRAV38-1"),
            ("TCRBV21S1*01", "TRBV11-1*01"),
        ),
    )
    def test_resolve_alternate_tr_names(self, symbol, expected):
        result = tr.standardize(symbol=symbol, species="homosapiens")

        assert result.symbol == expected

    @pytest.mark.parametrize(
        ("symbol", "expected", "species"),
        (
            ("TRAV14DV4", "TRAV14/DV4", "homosapiens"),
            ("TRBV20OR9-2", "TRBV20/OR9-2", "homosapiens"),
            ("TRBV01", "TRBV1", "homosapiens"),
            ("TCRBV1", "TRBV1", "homosapiens"),
            ("TRAV14", "TRAV14/DV4", "homosapiens"),
            ("TRDV4", "TRAV14/DV4", "homosapiens"),
            ("TCRAV13S2", "TRAV13-2", "homosapiens"),
            ("TCRAV38S2", "TRAV38-2/DV8", "homosapiens"),
            ("TCRAV30-1", "TRAV30", "homosapiens"),
            ("TCRDV01-01*01", "TRDV1*01", "homosapiens"),
            ("TCRAV14/4", "TRAV14/DV4", "homosapiens"),
            ("TCRAV36*01", "TRAV36/DV7*01", "homosapiens"),
            ("29/DV5*01", "TRAV29/DV5*01", "homosapiens"),
            ("TCRBJ2.7", "TRBJ2-7", "homosapiens"),
            ("TRAV15-1", "TRAV15", "homosapiens"),
            ("TRAV15-1", "TRAV15-1/DV6-1", "musmusculus"),
            ("TCRBJ2S1", "TRBJ2-1", "musmusculus"),

        ),
    )
    def test_various_typos(self, symbol, expected, species):
        result = tr.standardize(symbol=symbol, species=species)

        assert result.symbol == expected
        assert result.species == species


class TestStandardizeMusMusculus:
    @pytest.mark.parametrize("symbol", VALID_MUSMUSCULUS_TR)
    def test_already_correctly_formatted(self, symbol):
        result = tr.standardize(symbol=symbol, species="musmusculus")

        assert result.symbol == symbol
        assert result.species == "musmusculus"

    @pytest.mark.parametrize("symbol", ("foobar", "noice"))
    def test_invalid_tr(self, symbol, caplog):
        result = tr.standardize(symbol=symbol, species="musmusculus")
        assert "Failed to standardize" in caplog.text
        assert not result.is_standardized
        assert result.symbol is None
        assert result.species == "musmusculus"


class TestQuery:
    @pytest.mark.parametrize(
        ("species", "precision", "expected_len", "expected_in", "expected_not_in"),
        (
            ("homosapiens", "allele", 521, "TRAJ8*02", "TRAJ8"),
            ("homosapiens", "gene", 250, "TRAJ8", "TRAJ8*02"),
            ("musmusculus", "allele", 556, "TRAJ4*01", "TRAJ4"),
            ("musmusculus", "gene", 273, "TRAJ4", "TRAJ4*01"),
        ),
    )
    def test_query_all(
        self, species, precision, expected_len, expected_in, expected_not_in
    ):
        result = tr.query(species=species, precision=precision)

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
            ("homosapiens", "gene", "AJ", 61, "TRAJ11", "TRBJ1-6"),
            ("musmusculus", "gene", "GC", 4, "TRGC1", "TRGV1"),
        ),
    )
    def test_query_contains(
        self, species, precision, contains, expected_len, expected_in, expected_not_in
    ):
        result = tr.query(
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
            ("homosapiens", "gene", "F", 187, "TRBJ2-7", "TRBV12-1"),
            ("homosapiens", "allele", "NF", 133, "TRAV35*03", "TRAV35*01"),
            ("homosapiens", "gene", "NF", 76, "TRAV35", "TRAJ30"),
            ("musmusculus", "gene", "P", 59, "TRGC3", "TRDV5"),
            ("musmusculus", "allele", "ORF", 24, "TRBV24*03", "TRBV24*01"),
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
        result = tr.query(
            species=species, precision=precision, functionality=functionality
        )

        assert len(result) == expected_len
        assert expected_in in result
        assert not expected_not_in in result

    def test_query_default_species(self):
        result = tr.query(precision="gene", contains_pattern="AJ")
        assert len(result) == 61
        assert "TRAJ11" in result


class TestGetAaSequence:
    @pytest.mark.parametrize(
        ("symbol", "species", "expected"),
        (
            (
                "TRAV10*02",
                "homosapiens",
                {
                    "CDR1-IMGT": "VSPFSN",
                    "CDR2-IMGT": "MTFSENT",
                    "FR1-IMGT": "KNQVEQSPQSLIILEGKNCTLQCNYT",
                    "FR2-IMGT": "LRWYKQDTGRGPVSLTI",
                    "FR3-IMGT": "KSNGRYTATLDADTKQSSLHITASQLSDSASYIC",
                    "V-REGION": "KNQVEQSPQSLIILEGKNCTLQCNYTVSPFSNLRWYKQDTGRGPVSLTIMTFSENTKSNGRYTATLDADTKQSSLHITASQLSDSASYICVVS",
                    'functionality': 'F'
                },
            ),
            (
                "TRBD1*01",
                "homosapiens",
                {
                    "D-REGION": "GTGG",
                    'functionality': 'F'
                },
            ),
            (
                "TRAJ47*02",
                "homosapiens",
                {
                    "J-REGION": "EYGNKLVFGAGTILRVKS",
                    'J-MOTIF': 'FGAG',
                    'functionality': 'F'
                },
            ),
            (
                "TRAV1*01",
                "musmusculus",
                {
                    "CDR1-IMGT": "TSGFNG",
                    "CDR2-IMGT": "VVLDGL",
                    "FR1-IMGT": "GQGVEQPDNLMSVEGTFARVNCTYS",
                    "FR2-IMGT": "LSWYQQREGHAPVFLSY",
                    "FR3-IMGT": "KDSGHFSTFLSRSNGYSYLLLTELQIKDSASYLC",
                    "V-REGION": "GQGVEQPDNLMSVEGTFARVNCTYSTSGFNGLSWYQQREGHAPVFLSYVVLDGLKDSGHFSTFLSRSNGYSYLLLTELQIKDSASYLCAVR",
                    'functionality': 'F'
                },
            ),
        ),
    )
    def test_get_aa_sequence(self, symbol, species, expected):
        result = tr.get_aa_sequence(symbol=symbol, species=species)

        assert result == expected
