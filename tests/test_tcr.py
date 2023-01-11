import pytest
from tidytcells import tcr


class TestDecomposedHomoSapiensTcr:
    def test_compile(self):
        decomp_tcr = tcr._DecomposedHomoSapiensTcr(
            base = 'AV',
            num1 = '1',
            num2 = '1',
            p = False,
            or92 = False,
            d_designation = None,
            allele_num = '01'
        )

        assert decomp_tcr.compile() == 'TRAV1-1*01'


    def test_compile_with_dv(self):
        decomp_tcr = tcr._DecomposedHomoSapiensTcr(
            base='BV',
            num1='20',
            num2=None,
            p=False,
            or92=True,
            d_designation=None,
            allele_num='01'
        )

        assert decomp_tcr.compile() == 'TRBV20/OR9-2*01'
    

    def test_compile_with_or(self):
        decomp_tcr = tcr._DecomposedHomoSapiensTcr(
            base='AV',
            num1='14',
            num2=None,
            p=False,
            or92=False,
            d_designation='4',
            allele_num='01'
        )

        assert decomp_tcr.compile() == 'TRAV14/DV4*01'


    def test_compile_without_allele(self):
        decomp_tcr = tcr._DecomposedHomoSapiensTcr(
            base = 'AV',
            num1 = '1',
            num2 = '1',
            p = False,
            or92 = False,
            d_designation = None,
            allele_num = '01'
        )

        assert decomp_tcr.compile(allele=False) == 'TRAV1-1'


    def test_resolves_leading_zeros(self):
        decomp_tcr = tcr._DecomposedHomoSapiensTcr(
            base = 'AV',
            num1 = '000002',
            num2 = None,
            p = False,
            or92 = False,
            d_designation = None,
            allele_num = '01'
        )

        decomp_tcr.resolve()

        assert decomp_tcr.num1 == '2'
        assert decomp_tcr.num2 == None


    def test_resolves_d_designation(self):
        decomp_tcr = tcr._DecomposedHomoSapiensTcr(
            base = 'AV',
            num1 = '38',
            num2 = '2',
            p = False,
            or92 = False,
            d_designation = None,
            allele_num = '01'
        )

        decomp_tcr.resolve()

        assert decomp_tcr.d_designation == '8'


    @pytest.mark.parametrize(
        ('initial', 'expected'),
        (
            (None, None),
            ('1', '01'),
            ('01', '01'),
            ('10', '10'),
            ('foobar', None)
        )
    )
    def test_resolves_allele_number(self, initial, expected):
        decomp_tcr = tcr._DecomposedHomoSapiensTcr(
            base = 'AV',
            num1 = '1',
            num2 = '1',
            p = False,
            or92 = False,
            d_designation = None,
            allele_num = initial
        )

        decomp_tcr.resolve()

        assert decomp_tcr.allele_num == expected


    @pytest.mark.parametrize(
        ('decomp', 'expected'),
        (
            (('AV', '1', '1', False, False, None, '01'), '1'),
            (('AV', '2', '1', False, False, None, '01'), None)
        )
    )
    def test_resolves_num2_if_appropriate(self, decomp, expected):
        decomp_tcr = tcr._DecomposedHomoSapiensTcr(*decomp)

        decomp_tcr.resolve()

        assert decomp_tcr.num2 == expected


    @pytest.mark.parametrize(
        'decomp',
        (
            ('AV', '1', '1', False, False, None, '01'),
            ('BV', '20', None, False, True, None, '01'),
            ('AJ', '58', None, False, False, None, '01')
        )
    )
    def test_valid(self, decomp):
        decomp_tcr = tcr._DecomposedHomoSapiensTcr(*decomp)

        assert decomp_tcr.valid


    @pytest.mark.parametrize(
        'decomp',
        (
            ('AV', '99', None, False, False, None, None),
            ('BV', '20', '2', False, False, None, None),
            ('AJ', '58', None, False, False, None, '02')
        )
    )
    def test_invalid(self, decomp):
        decomp_tcr = tcr._DecomposedHomoSapiensTcr(*decomp)

        assert not decomp_tcr.valid
    

    @pytest.mark.parametrize(
        ('decomp', 'expected'),
        (
            (('AV', '1', '1', False, False, None, '01'), True),
            (('AV', '01', '00001', False, False, None, '1'), True)
        )
    )
    def test_returns_resolve_success(self, decomp, expected):
        decomp_tcr = tcr._DecomposedHomoSapiensTcr(*decomp)

        assert decomp_tcr.resolve() == expected


class TestStandardise:
    @pytest.mark.parametrize(
        'species',
        (
            'foobar',
            'yoinkdoink',
            ''
        )
    )
    def test_warning_unsupported_species(self, species):
        with pytest.warns(UserWarning, match='Unsupported'):
            result = tcr.standardise(
                gene_name='TRAV1-1*01',
                species=species
            )

        assert result == 'TRAV1-1*01'


class TestStandardiseHomoSapiens:
    @pytest.mark.parametrize(
        'gene',
        (
            'TRAV3*01',
            'TRAJ15*02',
            'TRBV5-1*02',
            'TRBJ2-1*01',
            'TRAJ15',
            'TRAV36/DV7',
            'TRBV29/OR9-2'
        )
    )
    def test_already_correctly_formatted(self, gene):
        result = tcr.standardise(
            gene_name=gene,
            species='HomoSapiens'
        )

        assert result == gene


    @pytest.mark.parametrize(
        'gene',
        (
            'foobar', 
            'TRAV3D-3*01', 
            '', 
            1234, 
            None
        )
    )
    def test_bad_tcr(self, gene):
        with pytest.warns(UserWarning, match='Unrecognised'):
            result = tcr.standardise(
                gene_name=gene,
                species='HomoSapiens'
            )

        assert result == None


    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('TRAV3*01&nbsp;', 'TRAV3*01'),
            ('TRBV07-01*01', 'TRBV7-1*01')
        )
    )
    def test_remove_pollutants(self, gene, expected):
        result = tcr.standardise(
            gene_name=gene,
            species='HomoSapiens'
        )

        assert result == expected


    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (  
            ('TRAV38-2DV8*01', 'TRAV38-2/DV8*01'),
            ('TRAV14DV4', 'TRAV14/DV4'),
            ('TRAV29DV5', 'TRAV29/DV5'),
            ('TRBV20OR9-2', 'TRBV20/OR9-2')
        )
    )
    def test_adds_slashes(self, gene, expected):
        result = tcr.standardise(
            gene_name=gene,
            species='HomoSapiens'
        )

        assert result == expected


    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('TCRAV32S1', 'TRAV25'),
            ('TCRAV14-2', 'TRAV38-1'),
            ('TCRBV21S1*01', 'TRBV11-1*01'),
            ('TCRBV12S3', 'TRBV12-3')
        )
    )
    def test_resolve_alternate_tcr_names(self, gene, expected):
        result = tcr.standardise(
            gene_name=gene,
            species='HomoSapiens'
        )

        assert result == expected


    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('TCRAV14/4', 'TRAV14/DV4'),
            ('36DV7*02', 'TRAV36/DV7*02'),
            ('TCRAV9-2*01/02/03/04', 'TRAV9-2'),
            ('TRAJ42.1', 'TRAJ42*01')
        )
    )
    def test_various_typos(self, gene, expected):
        result = tcr.standardise(
            gene_name=gene,
            species='HomoSapiens'
        )

        assert result == expected


    @pytest.mark.parametrize(
        'gene', ('TRAV36/DV1', 'TRAV1-1/DV1')
    )
    def test_fails_if_detected_av_dv_designation_mismatch(self, gene):
        with pytest.warns(UserWarning, match='Unrecognised'):
            result = tcr.standardise(
                gene_name=gene,
                species='HomoSapiens'
            )

        assert result == None


    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('TRAV2-1', 'TRAV2'),
            ('TRAV14-1/DV4', 'TRAV14/DV4')
        )
    )
    def test_removes_unnecessary_specifiers(self, gene, expected):
        result = tcr.standardise(
            gene_name=gene,
            species='HomoSapiens'
        )

        assert result == expected


    def test_default_homosapiens(self):
        result = tcr.standardise('TRAV3*01')

        assert result == 'TRAV3*01'