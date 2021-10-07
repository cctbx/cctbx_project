from __future__ import division

import numpy as np
from simtbx.diffBragg import utils
from scitbx.matrix import sqr
from simtbx.diffBragg.refiners.parameters import RangedParameter


def PAR_from_params(params, experiment, best=None):
    """

    :param params:  diffBragg params phil
    :param experiment: dxtbx Experiment list
    :param best: optional row of the stage 1 pandas dataframe corresponding to the experiment list
                This would contain all of the model parameters from stage 1
    :return: PAR object for storing refinement parameters
    """
    ParameterType = RangedParameter
    PAR = StageTwoParams()

    # each modeler has an experiment self.E, with a crystal attached
    if best is not None:
        experiment.crystal.set_A(best.Amats.values[0])
    PAR.Umatrix = sqr(experiment.crystal.get_U())
    PAR.Bmatrix = sqr(experiment.crystal.get_B())

    # per shot Scale factor
    initG = params.init.G if best is None else best.spot_scales.values[0]
    #TODO reconcile the need to sqrt the gain for stage_two_refiner
    PAR.Scale = ParameterType(init=initG, minval=params.mins.G, maxval=params.maxs.G, fix=params.fix.G, sigma=params.sigmas.G,
                              center=params.centers.G, beta=params.betas.G)

    PAR.Nabc = [None]*3
    PAR.Ndef = [None]*3
    PAR.RotXYZ_params = [None]*3
    for i in range(3):
        initN = params.init.Nabc[i] if best is None else best.ncells.values[0][i]
        PAR.Nabc[i] = ParameterType(init=initN, minval=params.mins.Nabc[i],
                                    maxval=params.maxs.Nabc[i], fix=params.fix.Nabc, sigma=params.sigmas.Nabc[i],
                                    center=params.centers.Nabc[i], beta=params.betas.Nabc[i])

        PAR.RotXYZ_params[i] = ParameterType(init=0, minval=params.mins.RotXYZ[i],
                                            maxval=params.maxs.RotXYZ[i], fix=params.fix.RotXYZ,
                                            sigma=params.sigmas.RotXYZ[i],
                                            center=0, beta=params.betas.RotXYZ)
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
        p = ParameterType(init=val, minval=minval, maxval=maxval, fix=params.fix.ucell, sigma=params.sigmas.ucell[i_uc],
                          center=params.centers.ucell[i_uc], beta=params.betas.ucell[i_uc])
        PAR.ucell.append(p)
    PAR.ucell_man = ucell_man

    # detector distance param:
    init_shiftZ = params.init.detz_shift *1e-3  if best is None else best.detz_shift_mm.values[0]*1e-3
    PAR.detz_shift = ParameterType(init=init_shiftZ, sigma=params.sigmas.detz_shift, minval=params.mins.detz_shift*1e-3,
                                   maxval=params.maxs.detz_shift*1e-3, fix=params.fix.detz_shift,
                                   center=params.centers.detz_shift, beta=params.betas.detz_shift)

    PAR.B = ParameterType(init=params.init.B, sigma=params.sigmas.B, minval=params.mins.B, maxval=params.maxs.B, fix=True,
                          center=0, beta=1e8)

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
