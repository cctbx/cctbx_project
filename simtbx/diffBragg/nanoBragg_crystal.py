"""
organizer for setting the nanoBragg crystal properties
"""

from simtbx.nanoBragg import shapetype
from scitbx.matrix import sqr, col


class nanoBragg_crystal(object):

    def __init__(self):

        self.dxtbx_crystal = nanoBragg_crystal.dxtbx_crystal_from_ucell_and_symbol()
        self.miller_array = nanoBragg_crystal.dummie_Fhkl()
        self.xtal_shape = shapetype.Gauss
        self.Ncells_abc = 10, 10, 10
        self.mos_spread_deg = 0
        self.n_mos_domains = 1
        self.thick_mm = 0.1
        self.missetting_matrix = sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

    @property
    def Omatrix(self):
        """
        Change of basis operator
        """
        sgi = self.dxtbx_crystal.get_space_group().info()
        to_p1 = sgi.change_of_basis_op_to_primitive_setting()
        return sqr(to_p1.c_inv().r().transpose().as_double())

    @property
    def dxtbx_crystal(self):
        return self._dxtbx_crystal

    @dxtbx_crystal.setter
    def dxtbx_crystal(self, val):
        #sginfo = val.get_space_group().info()
        #self._dxtbx_crystal = val.change_basis(sginfo.change_of_basis_op_to_primitive_setting())
        self._dxtbx_crystal = val

    @property
    def missetting_matrix(self):
        return self._missetting_matrix

    @missetting_matrix.setter
    def missetting_matrix(self, val):
        self._missetting_matrix = val

    @property
    def miller_array(self):
        return self._miller_array

    @miller_array.setter
    def miller_array(self, val):
        if str(val.observation_type) == "xray.intensity":
            val = val.as_amplitude_array()
        val = val.expand_to_p1()
        val = val.generate_bijvoet_mates()
        self._miller_array = val

    @property
    def Ncells_abc(self):
        return self._Ncells_abc

    @Ncells_abc.setter
    def Ncells_abc(self, val):
        self._Ncells_abc = val

    @property
    def mos_spread_deg(self):
        return self._mos_spread_deg

    @mos_spread_deg.setter
    def mos_spread_deg(self, val):
        self._mos_spread_deg = val

    @property
    def n_mos_domains(self):
        return self._n_mos_domains

    @n_mos_domains.setter
    def n_mos_domains(self, val):
        self._n_mos_domains = val

    @property
    def a_b_c_realspace_misset(self):
        # NOTE: this can be done in two lines, but this is explicit and better to read
        A = sqr(self.dxtbx_crystal.get_A()).inverse()
        a, b, c = A.as_list_of_lists()
        arot = self.missetting_matrix * col(a)
        brot = self.missetting_matrix * col(b)
        crot = self.missetting_matrix * col(c)
        return list(arot), list(brot), list(crot)

    @property
    def astar_bstar_cstar_misset(self):
        arot, brot, crot = self.a_b_c_realspace_misset
        self.Amat_star_row_vecs = sqr(arot + brot + crot).inverse()
        astar_rot, bstar_rot, cstar_rot = self.Amat_star_row_vecs.transpose().as_list_of_lists()  # trick to get the rows as lists
        return astar_rot, bstar_rot, cstar_rot

    @property
    def xtal_shape(self):
        return self._xtal_shape

    @xtal_shape.setter
    def xtal_shape(self, val):
        self._xtal_shape = val

    @property
    def thick_mm(self):
        return self._thick_mm

    @thick_mm.setter
    def thick_mm(self, val):
        self._thick_mm = val

    @staticmethod
    def dxtbx_crystal_from_ucell_and_symbol(ucell_tuple_Adeg=(65, 40, 75, 90, 110, 90), symbol="P1211"):
        """
        :param ucell_tuple_Adeg:  unit cell tuple a,b,c al, be, ga in Angstom and degrees
        :param symbol: lookup symbol for space group, e.g. 'P1'
        :return:a default crystal in conventional orientation, a along x-axis
        """
        from cctbx import crystal
        from dxtbx.model.crystal import CrystalFactory
        symm = crystal.symmetry("%f,%f,%f,%f,%f,%f" % ucell_tuple_Adeg, symbol)

        ucell = symm.unit_cell()
        O = ucell.orthogonalization_matrix()
        real_space_a = O[0], O[3], O[6]
        real_space_b = O[1], O[4], O[7]
        real_space_c = O[2], O[5], O[8]

        hall_symbol = symm.space_group_info().type().hall_symbol()

        return CrystalFactory.from_dict(
            {'__id__': 'crystal',
                       'real_space_a': real_space_a,
                       'real_space_b': real_space_b,
                       'real_space_c': real_space_c,
                       'space_group_hall_symbol': hall_symbol})

    @staticmethod
    def dummie_Fhkl():
        from simtbx.nanoBragg.tst_nanoBragg_basic import fcalc_from_pdb
        return fcalc_from_pdb(resolution=2, algorithm="fft", wavelength=1)

    def dxtbx_crystal_with_missetting(self):
        from copy import deepcopy
        C = deepcopy(self.dxtbx_crystal)
        astar, bstar, cstar = self.astar_bstar_cstar_misset
        C.set_A(tuple(astar) + tuple(bstar) + tuple(cstar))
        return C

    @staticmethod
    def abcstar_from_abc(a, b, c):
        Amat_star_row_vecs = sqr(tuple(a) + tuple(b) + tuple(c)).inverse()
        astar_rot, bstar_rot, cstar_rot = Amat_star_row_vecs.transpose().as_list_of_lists()
        return astar_rot, bstar_rot, cstar_rot
