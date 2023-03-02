import pytest
from tidytcells import tcr
from tidytcells._resources import (
    HOMOSAPIENS_TCR,
    MUSMUSCULUS_TCR
)


class TestStandardise:
    @pytest.mark.parametrize(
        'species',
        (
            'foobar',
            'yoinkdoink',
            ''
        )
    )
    def test_unsupported_species(self, species):
        with pytest.warns(UserWarning, match='Unsupported'):
            result = tcr.standardise(
                gene='foobarbaz',
                species=species
            )

        assert result == 'foobarbaz'
    

    @pytest.mark.parametrize(
        'gene',
        (
            1234, 
            None
        )
    )
    def test_bad_type(self, gene):
        with pytest.raises(TypeError):
            tcr.standardise(gene=gene)


    def test_default_homosapiens(self):
        result = tcr.standardise('TRBV20/OR9-2*01')

        assert result == 'TRBV20/OR9-2*01'


    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('TRAV3*01&nbsp;', 'TRAV3*01'),
            (' TRAV3 * 01 ', 'TRAV3*01')
        )
    )
    def test_remove_pollutants(self, gene, expected):
        result = tcr.standardise(
            gene=gene,
            species='homosapiens'
        )

        assert result == expected


    @pytest.mark.filterwarnings('ignore:Failed to standardise')
    @pytest.mark.parametrize(
        ('gene', 'expected', 'enforce_functional'),
        (
            ('TRAV35*01', 'TRAV35*01', True),
            ('TRAV35*03', None, True),
            ('TRAV35*03', 'TRAV35*03', False),
            ('TRAV35', 'TRAV35', True),
            ('TRAV8-7', None, True),
            ('TRAV8-7', 'TRAV8-7', False)
        )
    )
    def test_enforce_functional(self, gene, expected, enforce_functional):
        result = tcr.standardise(
            gene=gene,
            species='homosapiens',
            enforce_functional=enforce_functional
        )

        assert result == expected


    @pytest.mark.parametrize(
        ('gene', 'expected', 'precision'),
        (
            ('TRBV24/OR9-2*01', 'TRBV24/OR9-2*01', 'allele'),
            ('TRAV16*01', 'TRAV16', 'gene')
        )
    )
    def test_precision(self, gene, expected, precision):
        result = tcr.standardise(
            gene=gene,
            species='homosapiens',
            precision=precision
        )

        assert result == expected
    

    def test_standardize(self):
        result = tcr.standardize('TRBV20/OR9-2*01')

        assert result == 'TRBV20/OR9-2*01'


    def test_gene_name(self):
        result = tcr.standardise(gene_name='TRBV20/OR9-2*01')

        assert result == 'TRBV20/OR9-2*01'


class TestStandardiseHomoSapiens:
    @pytest.mark.parametrize('gene', HOMOSAPIENS_TCR)
    def test_already_correctly_formatted(self, gene):
        result = tcr.standardise(
            gene=gene,
            species='homosapiens'
        )

        assert result == gene


    @pytest.mark.parametrize(
        'gene',
        (
            'foobar', 
            'TRAV3D-3*01'
        )
    )
    def test_invalid_tcr(self, gene):
        with pytest.warns(UserWarning, match='Failed to standardise'):
            result = tcr.standardise(
                gene=gene,
                species='homosapiens'
            )

        assert result == None


    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('TCRAV32S1', 'TRAV25'),
            ('TCRAV14S2', 'TRAV38-1'),
            ('TCRBV21S1*01', 'TRBV11-1*01')
        )
    )
    def test_resolve_alternate_tcr_names(self, gene, expected):
        result = tcr.standardise(
            gene=gene,
            species='homosapiens'
        )

        assert result == expected


    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('TRAV14DV4', 'TRAV14/DV4'),
            ('TRBV20OR9-2', 'TRBV20/OR9-2'),
            ('TRBV01', 'TRBV1'),
            ('TCRBV1', 'TRBV1'),
            ('TRAV14', 'TRAV14/DV4'),
            ('TRDV4', 'TRAV14/DV4')
        )
    )
    def test_various_typos(self, gene, expected):
        result = tcr.standardise(
            gene=gene,
            species='homosapiens'
        )

        assert result == expected


class TestStandardiseMusMusculus:
    @pytest.mark.parametrize('gene', MUSMUSCULUS_TCR)
    def test_already_correctly_formatted(self, gene):
        result = tcr.standardise(
            gene=gene,
            species='musmusculus'
        )

        assert result == gene


    @pytest.mark.parametrize(
        'gene',
        (
            'foobar', 
            'TRAV3*01'
        )
    )
    def test_inivalid_tcr(self, gene):
        with pytest.warns(UserWarning, match='Failed to standardise'):
            result = tcr.standardise(
                gene=gene,
                species='musmusculus'
            )

        assert result == None