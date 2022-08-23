from __future__ import division

from simtbx.nanoBragg import sim_data
from dxtbx.model.crystal import CrystalFactory
from simtbx.nanoBragg import nanoBragg_crystal
from simtbx.diffBragg import utils
from cctbx import miller, crystal


def test_S(S):
    print("instantiate")
    S.instantiate_diffBragg()
    print("get map")
    Amap = S.D.get_ASUid_map()

    print("verify")
    inds, _ = S.D.Fhkl_tuple
#   NOTE: are symm equivalent indices the same in different bases ?
    print("SYMBOL:", S.crystal.symbol, "UCELL:", S.crystal.dxtbx_crystal.get_unit_cell())
    sym = crystal.symmetry(S.crystal.dxtbx_crystal.get_unit_cell(), S.crystal.symbol)
    mset = miller.set(sym, inds, True)
    mset_asu = mset.map_to_asu()
    print("count")
    n1 = len(set(mset_asu.indices()))
    n2 = len(Amap)
    assert n1==n2, "%d, %d" % (n1, n2)

    print("compare indices")
    A = []
    for item in Amap.keys():
        h,k,l = map(int, item.split(","))
        A.append ((h,k,l))
    A = set(A)
    B = set(mset_asu.indices())
    assert A==B

    print("OK")


if __name__=="__main__":
    S = sim_data.SimData(True)
    test_S(S)

    S2 = sim_data.SimData()
    UCELL_A=40.3
    UCELL_B=180.3
    UCELL_C=142.6
    HALL=' C 2c 2'

    CRYSTAL_DICT={
        '__id__': 'crystal',
        'real_space_a': (UCELL_A, 0.0, 0.0),
        'real_space_b': (0.0, UCELL_B, 0.0),
        'real_space_c': (0.0, 0.0, UCELL_C),
        'space_group_hall_symbol': HALL}
    C = CrystalFactory.from_dict(CRYSTAL_DICT)
    nbC = nanoBragg_crystal.NBcrystal()
    nbC.dxtbx_crystal = C
    nbC.symbol = C.get_space_group().info().type().lookup_symbol()
    nbC.miller_array = utils.make_miller_array(nbC.symbol, C.get_unit_cell().parameters())
    S2.crystal = nbC
    test_S(S2)

    p65_cryst = {'__id__': 'crystal',
                 'real_space_a': (43.32309880004587, 25.5289818883498, 60.49634260901813),
                 'real_space_b': (34.201635357808115, -38.82573591182249, -59.255697149884924),
                 'real_space_c': (41.42476391176581, 229.70849483520402, -126.60059788183489),
                 'space_group_hall_symbol': ' P 65 2 (x,y,z+1/12)',
                 'ML_half_mosaicity_deg': 0.06671930026192037,
                 'ML_domain_size_ang': 6349.223840307989}
    p65_C = CrystalFactory.from_dict(p65_cryst)
    S3 = sim_data.SimData()

    nbC = nanoBragg_crystal.NBcrystal()
    nbC.dxtbx_crystal = p65_C
    nbC.symbol = p65_C.get_space_group().info().type().lookup_symbol()
    nbC.miller_array = utils.make_miller_array(nbC.symbol, p65_C.get_unit_cell().parameters())
    S3.crystal = nbC
    test_S(S3)
