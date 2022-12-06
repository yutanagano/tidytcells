import pytest
from tidytcells import mhc


class TestDecomposedHla:
    @pytest.mark.parametrize(
        ('decomp', 'expected'),
        (
            (('A', ['01', '01', '01', '01'], False, False), ('HLA-A*01:01', ':01:01')),
            (('B', ['01', '01'], False, False), ('HLA-B*01:01', None)),
            (('A', None, False, False), ('HLA-A', None)),
            (('A', ['01', '01', '01G'], True, False), ('HLA-A*01:01:01G', None))
        )
    )
    def test_compile(self, decomp, expected):
        decomp_hla = mhc._DecomposedHla(*decomp)

        assert decomp_hla.compile() == expected


    @pytest.mark.parametrize(
        ('gene_name', 'expected'),
        (
            ('Cw', 'C'),
            ('DPB', 'DPB1'),
            ('DRA1', 'DRA')
        )
    )
    def test_resolves_gene_name(self, gene_name, expected):
        decomp_hla = mhc._DecomposedHla(
            gene=gene_name,
            spec_fields=None,
            g_group=False,
            p_group=False
        )

        decomp_hla.resolve()

        assert decomp_hla.gene == expected
    

    def test_resolves_likely_missing_field_separators(self):
        decomp_hla = mhc._DecomposedHla(
            gene='A',
            spec_fields=['0101', '01'],
            g_group=False,
            p_group=False
        )

        decomp_hla.resolve()

        assert decomp_hla.spec_fields == ['01', '01', '01']
    

    def test_resolves_leading_zeros(self):
        decomp_hla = mhc._DecomposedHla(
            gene='A',
            spec_fields=['1', '1'],
            g_group=False,
            p_group=False
        )

        decomp_hla.resolve()

        assert decomp_hla.spec_fields == ['01', '01']


    @pytest.mark.parametrize(
        'decomp',
        (
            ('A', None, False, False),
            ('DRB1', ['14', '126', '02'], False, False),
            ('DQA2', ['01', '04', '12345', '12345'], False, False)
        )
    )
    def test_valid(self, decomp):
        decomp_mhc = mhc._DecomposedHla(*decomp)

        assert decomp_mhc.valid


    @pytest.mark.parametrize(
        'decomp',
        (
            ('FOOBAR', None, False, False),
            ('DRB1', ['123456789', '01'], False, False),
            ('DQA2', ['01', '99', '01', '03'], False, False),
            ('A', ['01', '01', '1', '1'], False, False),
            ('A', ['01', '01', 'foo', 'bar'], False, False)
        )
    )
    def test_invalid(self, decomp):
        decomp_mhc = mhc._DecomposedHla(*decomp)

        assert not decomp_mhc.valid


    @pytest.mark.parametrize(
        ('decomp', 'expected'),
        (
            (('A', ['1', '1'], False, False), True),
            (('A', ['01', '999'], False, False), False),
            (('A', ['1', '1', 'foo', 'bar'], False, False), False)
        )
    )
    def test_returns_resolve_success(self, decomp, expected):
        decomp_hla = mhc._DecomposedHla(*decomp)

        assert decomp_hla.resolve() == expected


class TestStandardise:
    @pytest.mark.parametrize(
        ('gene', 'species'),
        (
            ('HLA-A', 'HomoSapiens'),
            ('HLA-B*07', 'HomoSapiens'),
            ('HLA-C*01:02', 'HomoSapiens'),
            ('HLA-A*01:01:01G', 'HomoSapiens'),
            ('HLA-A*01:01P', 'HomoSapiens'),
            ('B2M', 'HomoSapiens')
        )
    )
    def test_already_correctly_formatted(self, gene, species):
        prot, allele = mhc.standardise(
            gene_name=gene,
            species=species
        )

        assert prot == gene
        assert allele == None

    @pytest.mark.parametrize(
        ('gene', 'species', 'ex_p', 'ex_a'),
        (
            ('HLA-DRA*01:01:01', 'HomoSapiens', 'HLA-DRA*01:01', ':01'),
            ('HLA-DQA1*01:01:01:01', 'HomoSapiens', 'HLA-DQA1*01:01', ':01:01')
        )
    )
    def test_split_allele_fields(self, gene, species, ex_p, ex_a):
        prot, allele = mhc.standardise(
            gene_name=gene,
            species=species
        )

        assert prot == ex_p
        assert allele == ex_a


    @pytest.mark.parametrize(
        ('gene', 'species'),
        (
            ('foobar', 'HomoSapiens'),
            ('yoinkdoink', 'HomoSapiens'),
            ('', 'HomoSapiens'),
            ('HLA-FOOBAR123456', 'HomoSapiens'),
            (1234, 'HomoSapiens'),
            (None, 'HomoSapiens')
        )
    )
    def test_bad_mhc(self, gene, species):
        with pytest.warns(UserWarning, match='Unrecognised'):
            prot, allele = mhc.standardise(
                gene_name=gene,
                species=species
            )
        
        assert prot == None
        assert allele == None


    @pytest.mark.parametrize(
        ('gene', 'species'),
        (
            ('HLA-A*01:01:1:1:1:1:1:1', 'HomoSapiens'),
        )
    )
    def test_bad_allele_designation(self, gene, species):
        with pytest.warns(UserWarning, match='Unrecognised'):
            prot, allele = mhc.standardise(
                gene_name=gene,
                species=species
            )
        
        assert prot == None
        assert allele == None


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
            prot, allele = mhc.standardise(
                gene_name='HLA-A*01:01:01:01',
                species=species
            )

        assert prot == 'HLA-A*01:01:01:01'
        assert allele == None


    @pytest.mark.parametrize(
        ('gene', 'species', 'ex_p'),
        (
            ('A01', 'HomoSapiens', 'HLA-A*01'),
            ('A1', 'HomoSapiens', 'HLA-A*01'),
            ('B*07', 'HomoSapiens', 'HLA-B*07'),
            ('DRB1*15:01', 'HomoSapiens', 'HLA-DRB1*15:01')
        )
    )
    def test_convert_from_shorthand(self, gene, species, ex_p):
        prot, allele = mhc.standardise(
            gene_name=gene,
            species=species
        )

        assert prot == ex_p
        assert allele == None


    @pytest.mark.parametrize(
        ('gene', 'species', 'ex_p'),
        (
            ('HLA-Cw*06:02', 'HomoSapiens', 'HLA-C*06:02'),
            ('HLA-DQA*01:01', 'HomoSapiens', 'HLA-DQA1*01:01'),
            ('HLA-DQB*05:01', 'HomoSapiens', 'HLA-DQB1*05:01'),
            ('HLA-DRA1*01:01', 'HomoSapiens', 'HLA-DRA*01:01'),
            ('HLA-BEL/COQ/DEL', 'HomoSapiens', 'HLA-Y')
        )
    )
    def test_fix_deprecated_names(self, gene, species, ex_p):
        prot, allele = mhc.standardise(
            gene_name=gene,
            species=species
        )

        assert prot == ex_p
        assert allele == None


    @pytest.mark.parametrize(
        ('gene', 'species', 'ex_p'),
        (
            ('HLA-DPB*01:01', 'HomoSapiens', 'HLA-DPB1*01:01'),
            ('HLA-DRB*01:01', 'HomoSapiens', 'HLA-DRB1*01:01')
        )
    )
    def test_desperate_gene_name_resolution(self, gene, species, ex_p):
        prot, allele = mhc.standardise(
            gene_name=gene,
            species=species
        )

        assert prot == ex_p
        assert allele == None


    @pytest.mark.parametrize(
        ('gene', 'species', 'ex_p', 'ex_a'),
        (
            ('HLA-A*01:01:01:01N', 'HomoSapiens', 'HLA-A*01:01', ':01:01'),
            ('HLA-A*01:01:01:01L', 'HomoSapiens', 'HLA-A*01:01', ':01:01'),
            ('HLA-A*01:01:01:01S', 'HomoSapiens', 'HLA-A*01:01', ':01:01'),
            ('HLA-A*01:01:01:01C', 'HomoSapiens', 'HLA-A*01:01', ':01:01'),
            ('HLA-A*01:01:01:01A', 'HomoSapiens', 'HLA-A*01:01', ':01:01'),
            ('HLA-A*01:01:01:01Q', 'HomoSapiens', 'HLA-A*01:01', ':01:01')
        )
    )
    def test_remove_expression_qualifier(self, gene, species, ex_p, ex_a):
        prot, allele = mhc.standardise(
            gene_name=gene,
            species=species
        )

        assert prot == ex_p
        assert allele == ex_a


    @pytest.mark.parametrize(
        ('gene', 'species', 'ex_p'),
        (
            ('HLA-B*5701', 'HomoSapiens', 'HLA-B*57:01'),
        )
    )
    def test_add_missing_colon(self, gene, species, ex_p):
        prot, allele = mhc.standardise(
            gene_name=gene,
            species=species
        )

        assert prot == ex_p
        assert allele == None


    @pytest.mark.parametrize(
        ('gene', 'species', 'ex_p', 'ex_a'),
        (
            ('HLA-B*7', 'HomoSapiens', 'HLA-B*07', None),
            ('HLA-B*07:2', 'HomoSapiens', 'HLA-B*07:02', None),
            ('HLA-B*07:02:1', 'HomoSapiens', 'HLA-B*07:02', ':01'),
            ('HLA-B*07:02:01:1', 'HomoSapiens', 'HLA-B*07:02', ':01:01')
        )
    )
    def test_add_leading_zero(self, gene, species, ex_p, ex_a):
        prot, allele = mhc.standardise(
            gene_name=gene,
            species=species
        )

        assert prot == ex_p
        assert allele == ex_a


    def test_default_homosapiens(self):
        prot, allele = mhc.standardise('HLA-B*07')

        assert prot == 'HLA-B*07'
        assert allele == None


class TestGetChain:
    @pytest.mark.parametrize(
        ('gene_name', 'expected'),
        (
            ('HLA-A#', 'alpha'),
            ('HLA-B#', 'alpha'),
            ('HLA-C#', 'alpha'),
            ('HLA-DPA#', 'alpha'),
            ('HLA-DQA#', 'alpha'),
            ('HLA-DRA#', 'alpha'),
            ('HLA-E#', 'alpha'),
            ('HLA-F#', 'alpha'),
            ('HLA-G#', 'alpha'),
            ('HLA-DPB#', 'beta'),
            ('HLA-DQB#', 'beta'),
            ('HLA-DRB#', 'beta'),
            ('B2M', 'beta')
        )
    )
    def test_get_chain(self, gene_name, expected):
        result = mhc.get_chain(gene_name=gene_name)

        assert result == expected

    
    @pytest.mark.parametrize(
        'gene_name', ('foo', 'HLA', '0', 1234, None)
    )
    def test_unrecognised_gene_names(self, gene_name):
        with pytest.warns(UserWarning):
            mhc.get_chain(gene_name=gene_name)


class TestClassify:
    @pytest.mark.parametrize(
        ('gene_name', 'expected'),
        (
            ('HLA-A#', 1),
            ('HLA-B#', 1),
            ('HLA-C#', 1),
            ('HLA-DPA#', 2),
            ('HLA-DQA#', 2),
            ('HLA-DRA#', 2),
            ('HLA-DPB#', 2),
            ('HLA-DQB#', 2),
            ('HLA-DRB#', 2),
            ('HLA-E#', 1),
            ('HLA-F#', 1),
            ('HLA-G#', 1),
            ('B2M', 1)
        )
    )
    def test_classify(self, gene_name, expected):
        result = mhc.classify(gene_name)

        assert result == expected