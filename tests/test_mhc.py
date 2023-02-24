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
                gene='HLA-A*01:01:01:01',
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
            mhc.standardise(gene=gene)


    def test_default_homosapiens(self):
        result = mhc.standardise('HLA-B*07')

        assert result == 'HLA-B*07'


    @pytest.mark.parametrize(
        ('gene', 'expected', 'precision'),
        (
            ('HLA-DRB3*01:01:02:01', 'HLA-DRB3*01:01:02:01', 'allele'),
            ('HLA-DRB3*01:01:02:01', ('HLA-DRB3*01:01', ':02:01'), 'protein'),
            ('HLA-DRB3*01:01:02:01', 'HLA-DRB3', 'gene')
        )
    )
    def test_precision(self, gene, expected, precision):
        result = mhc.standardise(
            gene=gene,
            species='homosapiens',
            precision=precision
        )

        assert result == expected



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
            gene=gene,
            species='homosapiens'
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
        with pytest.warns(UserWarning, match='Failed to standardise'):
            result = mhc.standardise(
                gene=gene,
                species='homosapiens'
            )
        
        assert result == None


    @pytest.mark.parametrize(
        'gene',
        (
            'HLA-A*01:01:1:1:1:1:1:1',
        )
    )
    def test_bad_allele_designation(self, gene):
        with pytest.warns(UserWarning, match='Failed to standardise'):
            result = mhc.standardise(
                gene=gene,
                species='homosapiens'
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
            gene=gene,
            species='homosapiens'
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
            gene=gene,
            species='homosapiens'
        )

        assert result == 'HLA-A*01:01:01:01'


    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('HLA-B8', 'HLA-B*08'),
            ('A*01:01', 'HLA-A*01:01'),
            ('A1', 'HLA-A*01')
        )
    )
    def test_various_typos(self, gene, expected):
        result = mhc.standardise(
            gene=gene,
            species='homosapiens'
        )

        assert result == expected


class TestGetChain:
    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('HLA-A', 'alpha'),
            ('HLA-B', 'alpha'),
            ('HLA-C', 'alpha'),
            ('HLA-DPA1', 'alpha'),
            ('HLA-DQA2', 'alpha'),
            ('HLA-DRA', 'alpha'),
            ('HLA-E', 'alpha'),
            ('HLA-F', 'alpha'),
            ('HLA-G', 'alpha'),
            ('HLA-DPB2', 'beta'),
            ('HLA-DQB1', 'beta'),
            ('HLA-DRB3', 'beta'),
            ('B2M', 'beta')
        )
    )
    def test_get_chain(self, gene, expected):
        result = mhc.get_chain(gene_name=gene)

        assert result == expected

    
    @pytest.mark.parametrize(
        'gene', ('foo', 'HLA', '0')
    )
    def test_unrecognised_gene_names(self, gene):
        with pytest.warns(UserWarning):
            result = mhc.get_chain(gene_name=gene)

        assert result == None

    
    @pytest.mark.parametrize(
        'gene', (1234, None)
    )
    def test_bad_type(self, gene):
        with pytest.raises(TypeError):
            mhc.get_chain(gene)


class TestGetClass:
    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('HLA-A', 1),
            ('HLA-B', 1),
            ('HLA-C', 1),
            ('HLA-DPA1', 2),
            ('HLA-DQA2', 2),
            ('HLA-DRA', 2),
            ('HLA-DPB2', 2),
            ('HLA-DQB1', 2),
            ('HLA-DRB3', 2),
            ('HLA-E', 1),
            ('HLA-F', 1),
            ('HLA-G', 1),
            ('B2M', 1)
        )
    )
    def test_get_class(self, gene, expected):
        result = mhc.get_class(gene)

        assert result == expected
    

    @pytest.mark.parametrize(
        'gene', ('foo', 'HLA', '0')
    )
    def test_unrecognised_gene_names(self, gene):
        with pytest.warns(UserWarning):
            result = mhc.get_chain(gene_name=gene)

        assert result == None

    
    @pytest.mark.parametrize(
        'gene', (1234, None)
    )
    def test_bad_type(self, gene):
        with pytest.raises(TypeError):
            mhc.get_chain(gene)