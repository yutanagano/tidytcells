import pytest
from tidytcells import tcr


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
                gene_name='foobarbaz',
                species=species
            )

        assert result == 'foobarbaz'


    def test_default_homosapiens(self):
        result = tcr.standardise('TRBV20/OR9-2*01')

        assert result == 'TRBV20/OR9-2*01'


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
            ''
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
        'gene',
        (
            1234, 
            None
        )
    )
    def test_bad_type(self, gene):
        with pytest.raises(TypeError):
            result = tcr.standardise(
                gene_name=gene,
                species='HomoSapiens'
            )


    @pytest.mark.parametrize(
        ('gene', 'expected'),
        (
            ('TRAV3*01&nbsp;', 'TRAV3*01'),
            (' TRAV3 * 01 ', 'TRAV3*01')
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
            ('TCRAV32S1', 'TRAV25'),
            ('TCRAV14S2', 'TRAV38-1'),
            ('TCRBV21S1*01', 'TRBV11-1*01')
        )
    )
    def test_resolve_alternate_tcr_names(self, gene, expected):
        result = tcr.standardise(
            gene_name=gene,
            species='HomoSapiens'
        )

        assert result == expected


# class TestStandardiseMusMusculus:
#     @pytest.mark.parametrize(
#         'gene',
#         (
#             'TRAV3-1*01',
#             'TRAJ15*01',
#             'TRBV5*02',
#             'TRBJ2-1*01',
#             'TRAJ15',
#             'TRAV21/DV12',
#             'TRAV3D-3',
#             'TRAV15D-1/DV6D-1'
#         )
#     )
#     def test_already_correctly_formatted(self, gene):
#         result = tcr.standardise(
#             gene_name=gene,
#             species='MusMusculus'
#         )

#         assert result == gene


#     @pytest.mark.parametrize(
#         'gene',
#         (
#             'foobar', 
#             'TRAV3*01', 
#             '', 
#             1234, 
#             None
#         )
#     )
#     def test_bad_tcr(self, gene):
#         with pytest.warns(UserWarning, match='Unrecognised'):
#             result = tcr.standardise(
#                 gene_name=gene,
#                 species='MusMusculus'
#             )

#         assert result == None


#     @pytest.mark.parametrize(
#         ('gene', 'expected'),
#         (
#             ('TRAV3D-3*01&nbsp;', 'TRAV3D-3*01'),
#             ('TRAV09D-01*01', 'TRAV9D-1*01')
#         )
#     )
#     def test_remove_pollutants(self, gene, expected):
#         result = tcr.standardise(
#             gene_name=gene,
#             species='MusMusculus'
#         )

#         assert result == expected


#     @pytest.mark.parametrize(
#         ('gene', 'expected'),
#         (  
#             ('TRAV4-4DV10*01', 'TRAV4-4/DV10*01'),
#             ('TRAV16DDV11', 'TRAV16D/DV11'),
#             ('TRAV15D-2DV6D-2', 'TRAV15D-2/DV6D-2')
#         )
#     )
#     def test_adds_slashes(self, gene, expected):
#         result = tcr.standardise(
#             gene_name=gene,
#             species='MusMusculus'
#         )

#         assert result == expected


#     @pytest.mark.parametrize(
#         ('gene', 'expected'),
#         (
#             ('TRAV14D-3/8', 'TRAV14D-3/DV8'),
#             ('21DV12*02', 'TRAV21/DV12*02'),
#             ('TRAV9D-2*01/02/03/04', 'TRAV9D-2'),
#             ('TRAJ42.1', 'TRAJ42*01')
#         )
#     )
#     def test_various_typos(self, gene, expected):
#         result = tcr.standardise(
#             gene_name=gene,
#             species='MusMusculus'
#         )

#         assert result == expected


#     def test_resolves_d_designation(self):
#         result = tcr.standardise(
#             gene_name='TRAV21',
#             species='MusMusculus'
#         )

#         assert result == 'TRAV21/DV12'


#     @pytest.mark.parametrize(
#         'gene', ('TRAV13-4/DV1', 'TRAV4-4/DV1')
#     )
#     def test_fails_if_detected_av_dv_designation_mismatch(self, gene):
#         with pytest.warns(UserWarning, match='Unrecognised'):
#             result = tcr.standardise(
#                 gene_name=gene,
#                 species='MusMusculus'
#             )

#         assert result == None


#     @pytest.mark.parametrize(
#         ('gene', 'expected'),
#         (
#             ('TRAV19-1', 'TRAV19'),
#             ('TRAV21-1/DV12', 'TRAV21/DV12')
#         )
#     )
#     def test_removes_unnecessary_specifiers(self, gene, expected):
#         result = tcr.standardise(
#             gene_name=gene,
#             species='MusMusculus'
#         )

#         assert result == expected