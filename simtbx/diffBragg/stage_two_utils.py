from __future__ import division

import numpy as np
from simtbx.diffBragg import utils
from scitbx.matrix import sqr
from simtbx.diffBragg.refiners.parameters import RangedParameter
from simtbx.nanoBragg.utils import H5AttributeGeomWriter
from simtbx.diffBragg import hopper_utils


def create_gain_map(gain_map_file, expt=None, outname=None, convert_to_photons=True):
    """
    :param gain_map_file: output from stage_two_refiner, contains info for constructing the optimized gain map
    :param expt: optional, apply the gain correction to an image in the experiment imageset, and writes an
    image file to disk containint before/after images, as well as the correction itself
    :param outname: the output file name, if expt is provided
    :param convert_to_photons: bool, usually expt's contain data that are in ADU units,
    however the optimized gain map was determined against data in photon units. Only set this to False if the expt points to
    data that are in photon units
    :return: gain map as a numpy array, same shape as tye detector
    """
    gain_data = np.load(gain_map_file)
    region_shape = tuple(gain_data["region_shape"])
    det_shape = tuple(gain_data["det_shape"])
    gain_per_region = gain_data["gain_per_region"]
    adu_per_photon = gain_data["adu_per_photon"]
    gain_map = regionize_detector(det_shape, region_shape, gain_per_region)
    region_marker = gain_data["num_times_region_was_modeled"][()]
    ntimes = gain_data["num_times_pixel_was_modeled"][()]
    mask = ntimes > 0
    region_marker = regionize_detector(det_shape, region_shape, region_marker) #.astype(np.int))
    region_marker = region_marker > 0
    if expt is not None:
        assert outname is not None
        iset = expt.imageset
        data = np.array([a.as_numpy_array() for a in iset.get_raw_data(0)])
        if convert_to_photons:
            data /= adu_per_photon
        with H5AttributeGeomWriter(outname, det_shape, 5, detector=expt.detector, beam=expt.beam) as writer:
            writer.add_image(data)
            writer.add_image(data*gain_map)
            writer.add_image(gain_map)
            writer.add_image(gain_data["regions"])
            writer.add_image(ntimes)
        maskname = outname+".mask"
        utils.save_numpy_mask_as_flex(region_marker, maskname)
        maskname2 = outname+".mask2"
        utils.save_numpy_mask_as_flex(mask, maskname2)
        print("Wrote %s and %s and %s, open with dials.image_viewer" %(outname, maskname, maskname2))

    vals = gain_map[region_marker]
    print(np.mean(vals), np.std(vals), np.min(vals), np.max(vals))
    return gain_map


def regionize_detector(det_shape, region_shape, gain_map_values=None):
    """
    Create a 3-d numpy array (shaped after a multi panel detector)
    where each pixel corresponds to a unique region of shape `region_shape`
    :param det_shape: 3 tuple, shape of the detecor (num_panels, slow_dim, fast_dim)
    :param region_shape: 2 tuple, the desired shape of each region in pixels (nslow, nfast)
    :param gain_map_values: array of gains for each region (advanced usage, one would never use this directly)
    :return: 3-d numpy array of regions
    """
    Y,X = region_shape
    regions = np.zeros(det_shape)
    if gain_map_values is None:
        regions = regions.astype(np.int32)
    numPan, slowDim, fastDim = det_shape
    nx = np.array_split(np.arange(fastDim),int(fastDim/X)+1)
    ny = np.array_split(np.arange(slowDim),int(slowDim/Y)+1)
    i_region = 0
    for pid in range(numPan):
        for j in range(len(ny)):
            for i in range(len(nx)):
                istart = nx[i][0]
                istop = nx[i][-1]+1
                jstart = ny[j][0]
                jstop = ny[j][-1]+1
                if gain_map_values is not None:
                    regions[pid,jstart:jstop, istart:istop] = gain_map_values[i_region]
                else:
                    regions[pid, jstart:jstop, istart:istop] = i_region
                i_region += 1
    return regions


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
    PAR.eta = [None]*3
    PAR.RotXYZ_params = [None]*3
    PAR.diffuse_sigma = [None]*3
    PAR.diffuse_gamma = [None]*3

    eta_min = params.mins.eta_abc
    init_eta = params.init.eta_abc if best is None else best.eta_abc.values[0]
    if tuple(init_eta) == (0,0,0):
        eta_min=-1e-10,-1e-10,-1e-10
        if not params.simulator.crystal.num_mosaicity_samples == 1:
            raise ValueError("if all eta_abc are 0,0,0, num_mosaicity_samples should be 1")

    # TODO allow setting diffuse gamma/sigma from stage 1 (e.g. from the `best` dataframe)
    init_diffuse_sigma = params.init.diffuse_sigma
    init_diffuse_gamma = params.init.diffuse_gamma

    for i in range(3):
        initN = params.init.Nabc[i] if best is None else best.ncells.values[0][i]
        PAR.Nabc[i] = ParameterType(init=initN, minval=params.mins.Nabc[i],
                                    maxval=params.maxs.Nabc[i], fix=params.fix.Nabc, sigma=params.sigmas.Nabc[i],
                                    center=params.centers.Nabc[i] if params.centers.Nabc is not None else None,
                                    beta=params.betas.Nabc[i] if params.betas.Nabc is not None else None)

        initN = params.init.Ndef[i] if best is None else best.ncells_def.values[0][i]
        PAR.Ndef[i] = ParameterType(init=initN, minval=params.mins.Ndef[i],
                                    maxval=params.maxs.Ndef[i], fix=params.fix.Ndef, sigma=params.sigmas.Ndef[i],
                                    center=params.centers.Ndef[i] if params.centers.Ndef is not None else None,
                                    beta=params.betas.Ndef[i] if params.betas.Ndef is not None else None)

        PAR.RotXYZ_params[i] = ParameterType(init=0, minval=params.mins.RotXYZ[i],
                                            maxval=params.maxs.RotXYZ[i], fix=params.fix.RotXYZ,
                                            sigma=params.sigmas.RotXYZ[i],
                                            center=0 if params.betas.RotXYZ is not None else None, beta=params.betas.RotXYZ)

        PAR.eta[i] = ParameterType(init=init_eta[i], minval=eta_min[i],
                                       maxval=params.maxs.eta_abc[i], fix=params.fix.eta_abc,
                                       sigma=params.sigmas.eta_abc[i],
                                       center=params.betas.eta_abc[i] if params.centers.eta_abc is not None else None,
                                   beta=params.betas.eta_abc[i] if params.betas.eta_abc is not None else None)

        # TODO: diffuse scattering terms
        PAR.diffuse_sigma[i] = ParameterType(init=init_diffuse_sigma[i], minval=params.mins.diffuse_sigma[i],
                                        maxval=params.maxs.diffuse_sigma[i], fix=params.fix.diffuse_sigma,
                                        sigma=params.sigmas.diffuse_sigma[i],
                                        center=params.centers.diffuse_sigma[i] if params.centers.diffuse_sigma is not None else None,
                                             beta=params.betas.diffuse_sigma[i] if params.betas.diffuse_sigma is not None else None)

        PAR.diffuse_gamma[i] = ParameterType(init=init_diffuse_gamma[i], minval=params.mins.diffuse_gamma[i],
                                             maxval=params.maxs.diffuse_gamma[i], fix=params.fix.diffuse_gamma,
                                             sigma=params.sigmas.diffuse_gamma[i],
                                             center=params.centers.diffuse_gamma[i] if params.centers.diffuse_gamma is not None else None,
                                             beta=params.betas.diffuse_gamma[i] if params.betas.diffuse_gamma is not None else None)

    # unit cell parameters
    ucell_man = utils.manager_from_crystal(experiment.crystal)  # Note ucell man contains the best parameters (if best is not None)
    ucell_vary_perc = params.ucell_edge_perc / 100.
    PAR.ucell = []
    for i_uc, (name, val) in enumerate(zip(ucell_man.variable_names, ucell_man.variables)):
        if "Ang" in name:
            minval = val - ucell_vary_perc * val
            maxval = val + ucell_vary_perc * val
            if name == 'a_Ang':
                cent = params.centers.ucell_a
                beta = params.betas.ucell_a
            elif name == 'b_Ang':
                cent = params.centers.ucell_b
                beta = params.betas.ucell_b
            else:
                cent = params.centers.ucell_c
                beta = params.betas.ucell_c
        else:
            val_in_deg = val * 180 / np.pi
            minval = (val_in_deg - params.ucell_ang_abs) * np.pi / 180.
            maxval = (val_in_deg + params.ucell_ang_abs) * np.pi / 180.
            if name == 'alpha_rad':
                cent = params.centers.ucell_alpha
                beta = params.betas.ucell_alpha
            elif name == 'beta_rad':
                cent = params.centers.ucell_beta
                beta = params.betas.ucell_beta
            else:
                cent = params.centers.ucell_gamma
                beta = params.betas.ucell_gamma
            assert cent is not None
            assert beta is not None
            cent = cent * np.pi / 180.

        p = ParameterType(init=val, minval=minval, maxval=maxval, fix=params.fix.ucell, sigma=params.sigmas.ucell[i_uc],
                          center=cent, beta=beta)
        PAR.ucell.append(p)
    PAR.ucell_man = ucell_man

    # detector distance param:
    init_shiftZ = params.init.detz_shift *1e-3  if best is None else best.detz_shift_mm.values[0]*1e-3
    PAR.detz_shift = ParameterType(init=init_shiftZ, sigma=params.sigmas.detz_shift, minval=params.mins.detz_shift*1e-3,
                                   maxval=params.maxs.detz_shift*1e-3, fix=params.fix.detz_shift,
                                   center=params.centers.detz_shift, beta=params.betas.detz_shift)

    PAR.B = ParameterType(init=params.init.B, sigma=params.sigmas.B, minval=params.mins.B, maxval=params.maxs.B, fix=True,
                          center=params.centers.B, beta=params.betas.B)

    lam0, lam1 = params.init.spec
    if best is not None:
        lam0, lam1 = hopper_utils.get_lam0_lam1_from_pandas(best)
    PAR.spec_coef = []
    for i_p, init_val in enumerate((lam0, lam1)):
        p = ParameterType(init=init_val, sigma=params.sigmas.spec[i_p],
                          center=params.centers.spec[i_p] if params.centers.spec is not None else None,
                          beta=params.betas.spec[i_p] if params.betas.spec is not None else None,
                          fix=params.fix.spec,
                          minval=params.mins.spec[i_p], maxval=params.maxs.spec[i_p])
        PAR.spec_coef.append(p)

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
        self.diffuse_sigma = None
        self.diffuse_gamma = None
