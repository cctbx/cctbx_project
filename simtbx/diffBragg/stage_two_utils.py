from __future__ import division

import numpy as np
from simtbx.diffBragg import utils
from scitbx.matrix import sqr
from simtbx.diffBragg.refiners.parameters import RangedParameter


def PAR_from_params(params, experiment, best=None):
    ParameterType = RangedParameter
    PAR = StageTwoParams()

    # each modeler has an experiment self.E, with a crystal attached
    if best is not None:
        experiment.crystal.set_A(best.Amats.values[0])
    PAR.Umatrix = sqr(experiment.crystal.get_U())
    PAR.Bmatrix = sqr(experiment.crystal.get_B())

    # per shot Scale factor
    p = ParameterType()
    p.sigma = params.sigmas.G
    if best is not None:
        p.init = best.spot_scales.values[0]
    else:
        p.init = params.init.G
    p.init = np.sqrt(p.init)
    p.minval = params.mins.G
    p.maxval = params.maxs.G
    p.fix = params.fix.G
    PAR.Scale = p

    PAR.Nabc = []
    PAR.Ndef = []
    for i_N in range(3):
        p = ParameterType()
        p.sigma = params.sigmas.Nabc[i_N]
        if best is not None:
            p.init = best.ncells.values[0][i_N]
        else:
            p.init = params.init.Nabc[i_N]
        p.minval = params.mins.Nabc[i_N]
        p.maxval = params.maxs.Nabc[i_N]
        p.fix = params.fix.Nabc
        PAR.Nabc.append(p)

        p = ParameterType()
        p.sigma = params.sigmas.Ndef[i_N]
        if best is not None:
            p.init = best.ncells_def.values[0][i_N]
        else:
            p.init = params.init.Ndef[i_N]
        p.minval = params.mins.Ndef[i_N]
        p.maxval = params.maxs.Ndef[i_N]
        p.fix = params.fix.Ndef
        PAR.Ndef.append(p)

    PAR.RotXYZ_params = []
    for i_rot in range(3):
        p = ParameterType()
        p.sigma = params.sigmas.RotXYZ[i_rot]
        if best is not None:
            p.init = best.values[0][i_rot]
        else:
            p.init = 0
        p.minval =params.mins.RotXYZ[i_rot]
        p.maxval = params.maxs.RotXYZ[i_rot]
        p.fix = params.fix.RotXYZ
        PAR.RotXYZ_params.append(p)

    # unit cell parameters
    ucell_man = utils.manager_from_crystal(experiment.crystal)  # Note ucell man contains the best parameters (if best is not None)
    ucell_vary_perc = params.ucell_edge_perc / 100.
    PAR.ucell = []
    for i_uc, (name, val) in enumerate(zip(ucell_man.variable_names, ucell_man.variables)):
        if "Ang" in name:
            minval = val - ucell_vary_perc * val
            maxval = val + ucell_vary_perc * val
        else:
            val_in_deg = val * 180 / np.pi
            minval = (val_in_deg - params.ucell_ang_abs) * np.pi / 180.
            maxval = (val_in_deg + params.ucell_ang_abs) * np.pi / 180.
        p = ParameterType()
        p.sigma = params.sigmas.ucell[i_uc]
        p.init = val
        p.minval = minval
        p.maxval = maxval
        p.fix = params.fix.ucell
        if not params.quiet: print(
            "Unit cell variable %s (currently=%f) is bounded by %f and %f" % (name, val, minval, maxval))
        PAR.ucell.append(p)
    PAR.ucell_man = ucell_man

    # detector distance param:
    p = ParameterType()
    if best is not None:
        p.init = best.detz_shift_mm.values[0]*1e-3
    else:
        p.init = params.init.detz_shift *1e-3
    p.sigma = params.sigmas.detz_shift
    p.minval = params.mins.detz_shift * 1e-3
    p.maxval = params.maxs.detz_shift * 1e-3
    p.fix = params.fix.detz_shift
    PAR.detz_shift = p

    p = ParameterType()
    p.init = params.init.B
    p.sigma = params.sigmas.B
    p.minval = params.mins.B
    p.maxval = params.maxs.B
    p.fix = params.fix.B
    PAR.B = p

    return PAR


class StageTwoParams:
    def __init__(self):
        self.ucell = None
        self.RotXYZ = None
        self.Scale = None
        self.Nabc = None
        self.Ndef = None
        self.eta = None
        self.spec_coef = None
        self.Fhkl = None
        self.detz_shift = None
        self.paneRot = None
        self.PanXYZ = None
