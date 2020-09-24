

from simtbx.diffBragg.refiners.crystal_systems import OrthorhombicManager, TetragonalManager, MonoclinicManager, HexagonalManager
import numpy as np


def from_crystal(crystal):
    """

    :param crystal:  dxtbx crystal model
    :return:
    """

    a, b, c, al, be, ga = crystal.get_unit_cell().parameters()
    if a == b and a != c and np.allclose([al, be, ga], [90]*3):
        manager = TetragonalManager(a=a, c=c)

    elif a != b and b != c and a != c and np.allclose([al, be, ga], [90]*3):
        manager = OrthorhombicManager(a=a, b=b, c=c)

    elif a == b and a != c and np.allclose([al, ga], [90]*2) and np.allclose([be], [120]):
        manager = HexagonalManager(a=a, c=c)

    elif a != b and b != c and a != c and np.allclose([al, ga], [90]*2) and not np.allclose([be], [120]):
        manager = MonoclinicManager(a=a, b=b, c=c, beta=be)

    else:
        raise NotImplementedError("Not yet implemented for crystal model")

    return manager