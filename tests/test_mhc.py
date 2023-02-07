import pytest
from tidytcells import mhc


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
            result = mhc.standardise(
                gene_name='HLA-A*01:01:01:01',
                species=species
            )

        assert result == 'HLA-A*01:01:01:01'


    @pytest.mark.parametrize(
        'gene',
        (
            1234, 
            None
        )
    )
    def test_bad_type(self, gene):
        with pytest.raises(TypeError):
            mhc.standardise(gene_name=gene)


    def test_default_homosapiens(self):
        result = mhc.standardise('HLA-B*07')

        assert result == 'HLA-B*07'



class TestStandardiseHomoSapiens:
    @pytest.mark.parametrize(
        'gene',
        (
            ('HLA-A'),
            ('HLA-B*07'),
            ('HLA-C*01:02'),
            ('HLA-DRB3*01:01:02:01'),
            ('HLA-A*01:01:01G'),
            ('HLA-A*01:01P'),
            ('B2M')
        )
    )
    def test_already_correctly_formatted(self, gene):
        result = mhc.standardise(
            gene_name=gene,
            species='HomoSapiens'
        )

        assert result == gene


    @pytest.mark.parametrize(
        'gene',
        (
            'foobar',
            'yoinkdoink',
            'HLA-FOOBAR123456'
        )
    )
    def test_invalid_mhc(self, gene):
        with pytest.warns(UserWarning, match='Unrecognised'):
            result = mhc.standardise(
                gene_name=gene,
                species='HomoSapiens'
            )
        
        assert result == None


    @pytest.mark.parametrize(
        'gene',
        (
            'HLA-A*01:01:1:1:1:1:1:1',
        )
    )
    def test_bad_allele_designation(self, gene):
        with pytest.warns(UserWarning, match='Unrecognised'):
            result = mhc.standardise(
                gene_name=gene,
                species='HomoSapiens'
            )
        
        assert result == None


    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('D6S204', 'HLA-C'),
            ('HLA-DQA*01:01', 'HLA-DQA1*01:01'),
            ('HLA-DQB*05:01', 'HLA-DQB1*05:01'),
            ('HLA-DRA1*01:01', 'HLA-DRA*01:01')
        )
    )
    def test_fix_deprecated_names(self, gene, expected):
        result = mhc.standardise(
            gene_name=gene,
            species='HomoSapiens'
        )

        assert result == expected


    @pytest.mark.parametrize(
        'gene',
        (
            'HLA-A*01:01:01:01N',
            'HLA-A*01:01:01:01L',
            'HLA-A*01:01:01:01S',
            'HLA-A*01:01:01:01C',
            'HLA-A*01:01:01:01A',
            'HLA-A*01:01:01:01Q',
        )
    )
    def test_remove_expression_qualifier(self, gene):
        result = mhc.standardise(
            gene_name=gene,
            species='HomoSapiens'
        )

        assert result == 'HLA-A*01:01:01:01'