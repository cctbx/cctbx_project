from __future__ import absolute_import, division, print_function
from dials.algorithms.shoebox import MaskCode
from copy import deepcopy
from dials.model.data import Shoebox
import numpy as np

from scipy.optimize import dual_annealing, basinhopping
from collections import Counter
from scitbx.matrix import sqr, col
from dxtbx.model.experiment_list import ExperimentListFactory
from simtbx.nanoBragg.utils import downsample_spectrum
from dials.array_family import flex
from simtbx.nanoBragg.anisotropic_mosaicity import AnisoUmats
from simtbx.diffBragg import utils
from simtbx.diffBragg.refiners.parameters import NormalParameter, RangedParameter

ROTX_ID = 0
ROTY_ID = 1
ROTZ_ID = 2
ROTXYZ_IDS = ROTX_ID, ROTY_ID, ROTZ_ID
NCELLS_ID = 9
UCELL_ID_OFFSET = 3
DETZ_ID = 10
FHKL_ID = 11


class DataModeler:

    def __init__(self, params):
        """ params is a simtbx.diffBragg.hopper phil"""
        self.params = params
        self.SIM = None
        self.E = None
        self.pan_fast_slow =None
        self.all_background =None
        self.roi_id =None
        self.u_id = None
        self.all_data =None
        self.all_sigmas =None
        self.all_trusted =None
        self.npix_total =None
        self.all_fast =None
        self.all_slow =None
        self.rois=None
        self.pids=None
        self.tilt_abc=None
        self.selection_flags=None
        self.background=None
        self.tilt_cov = None
        self.simple_weights = None
        self.refls_idx = None
        self.refls = None

    def clean_up(self):
        if self.SIM is not None:
            self.SIM.D.free_all()
            self.SIM.D.gpu_free()
            self.SIM.D.free_Fhkl2()

    def set_experiment(self, exp, load_imageset=True):
        if isinstance(exp, str):
            self.E = ExperimentListFactory.from_json_file(exp, load_imageset)[0]
        else:
            self.E = exp
        if self.params.opt_det is not None:
            opt_det_E = ExperimentListFactory.from_json_file(self.params.opt_det, False)[0]
            self.E.detector = opt_det_E.detector

    def load_refls(self, ref):
        if isinstance(ref, str):
            refls = flex.reflection_table.from_file(ref)
        else:
            # assert is a reflection table. ..
            refls = ref
        return refls

    def is_duplicate_hkl(self, refls):
        nref = len(refls)
        is_duplicate = np.zeros(nref, bool)
        if len(set(refls['miller_index'])) < nref:
            hkls = refls['miller_index']
            dupe_hkl = {h for h, count in Counter(hkls).items() if count > 1}
            for i_ref in range(nref):
                hh = refls[i_ref]['miller_index']
                is_duplicate[i_ref] = hh in dupe_hkl

        return is_duplicate

    def GatherFromReflectionTable(self, exp, ref):
        self.set_experiment(exp, load_imageset=False)
        self.refls = self.load_refls(ref)
        nref = len(self.refls)
        if nref ==0:
            return False
        self.refls_idx = list(range(nref))
        self.rois = [(x1, x2, y1, y2) for x1,x2,y1,y2,_,_ in self.refls["shoebox"].bounding_boxes()]
        self.pids = list(self.refls["panel"])

        npan = len(self.E.detector)
        nfast, nslow = self.E.detector[0].get_image_size()  # NOTE assumes all panels same shape
        img_data = np.zeros((npan, nslow, nfast))
        background = np.zeros_like(img_data)
        is_trusted = np.zeros((npan, nslow, nfast), bool)
        for i_ref in range(nref):
            ref = self.refls[i_ref]
            pid = ref['panel']
            x1,x2,y1,y2 = self.rois[i_ref]
            sb = ref['shoebox']
            img_data[pid,   y1:y2, x1:x2] = sb.data.as_numpy_array()[0]
            background[pid, y1:y2, x1:x2] = sb.background.as_numpy_array()[0]
            is_trusted[pid, y1:y2, x1:x2] = sb.mask.as_numpy_array()[0] == MaskCode.Valid  # I believe this is 1

        # can be used for Bfactor modeling
        self.Q = np.linalg.norm(self.refls["rlp"], axis=1)
        self.sigma_rdout = self.params.refiner.sigma_r / self.params.refiner.adu_per_photon

        self.data_to_one_dim(img_data, is_trusted, background)
        return True

    def GatherFromExperiment(self, exp, ref, remove_duplicate_hkl=True):
        self.set_experiment(exp, load_imageset=True)

        refls = self.load_refls(ref)

        img_data = utils.image_data_from_expt(self.E)
        img_data /= self.params.refiner.adu_per_photon
        is_trusted = utils.load_mask(self.params.roi.hotpixel_mask)
        hotpix_mask = None
        if is_trusted is not None:
            hotpix_mask = ~is_trusted
        self.sigma_rdout = self.params.refiner.sigma_r / self.params.refiner.adu_per_photon

        roi_packet = utils.get_roi_background_and_selection_flags(
            refls, img_data, shoebox_sz=self.params.roi.shoebox_size,
            reject_edge_reflections=self.params.roi.reject_edge_reflections,
            reject_roi_with_hotpix=self.params.roi.reject_roi_with_hotpix,
            background_mask=None, hotpix_mask=hotpix_mask,
            bg_thresh=self.params.roi.background_threshold,
            use_robust_estimation=not self.params.roi.fit_tilt,
            set_negative_bg_to_zero=self.params.roi.force_negative_background_to_zero,
            pad_for_background_estimation=self.params.roi.pad_shoebox_for_background_estimation,
            sigma_rdout=self.sigma_rdout, deltaQ=self.params.roi.deltaQ, experiment=self.E,
            weighted_fit=self.params.roi.fit_tilt_using_weights,
            tilt_relative_to_corner=self.params.relative_tilt, ret_cov=True)

        if roi_packet is None:
            return False

        self.rois, self.pids, self.tilt_abc, self.selection_flags, self.background, self.tilt_cov = roi_packet

        if remove_duplicate_hkl:
            is_not_a_duplicate = ~self.is_duplicate_hkl(refls)
            #from IPython import embed;embed();exit()
            self.selection_flags = np.logical_and( self.selection_flags, is_not_a_duplicate)

        if sum(self.selection_flags) == 0:
            if not self.params.quiet: print("No pixels slected, continuing")
            return False
        self.refls = refls
        self.refls_idx = [i_roi for i_roi in range(len(refls)) if self.selection_flags[i_roi]]

        self.rois = [roi for i_roi, roi in enumerate(self.rois) if self.selection_flags[i_roi]]
        self.tilt_abc = [abc for i_roi, abc in enumerate(self.tilt_abc) if self.selection_flags[i_roi]]
        self.pids = [pid for i_roi, pid in enumerate(self.pids) if self.selection_flags[i_roi]]
        self.tilt_cov = [cov for i_roi, cov in enumerate(self.tilt_cov) if self.selection_flags[i_roi]]
        self.Q = [np.linalg.norm(refls[i_roi]["rlp"]) for i_roi in range(len(refls)) if self.selection_flags[i_roi]]

        self.data_to_one_dim(img_data, is_trusted, self.background)
        return True

    def data_to_one_dim(self, img_data, is_trusted, background):
        all_data = []
        all_pid = []
        all_fast = []
        all_slow = []
        all_fast_relative = []
        all_slow_relative = []
        all_trusted = []
        all_sigmas = []
        all_background = []
        roi_id = []
        all_a, all_b, all_c = [], [], []
        all_q_perpix = []
        for i_roi in range(len(self.rois)):
            pid = self.pids[i_roi]
            x1, x2, y1, y2 = self.rois[i_roi]
            Y, X = np.indices((y2 - y1, x2 - x1))
            data = img_data[pid, y1:y2, x1:x2].copy()

            data = data.ravel()
            all_background += list(background[pid, y1:y2, x1:x2].ravel())
            trusted = is_trusted[pid, y1:y2, x1:x2].ravel()

            # TODO implement per-shot masking here
            #lower_cut = np.percentile(data, 20)
            #trusted[data < lower_cut] = False

            #d_strong_order = np.argsort(data)
            #trusted[d_strong_order[-1:]] = False
            all_trusted += list(trusted)
            #TODO ignore invalid value warning (handled below), or else mitigate it!
            all_sigmas += list(np.sqrt(data + self.sigma_rdout ** 2))
            all_fast += list(X.ravel() + x1)
            all_fast_relative += list(X.ravel())
            all_slow += list(Y.ravel() + y1)
            all_slow_relative += list(Y.ravel())
            all_data += list(data)
            npix = len(data)  # np.sum(trusted)
            all_pid += [pid] * npix
            roi_id += [i_roi] * npix
            #if self.tilt_abc is not None:
            #    a, b, c = self.tilt_abc[i_roi]
            #    all_a += [a] * npix
            #    all_b += [b] * npix
            #    all_c += [c] * npix
            all_q_perpix += [self.Q[i_roi]]*npix

        self.all_q_perpix = np.array(all_q_perpix)
        pan_fast_slow = np.ascontiguousarray((np.vstack([all_pid, all_fast, all_slow]).T).ravel())
        self.pan_fast_slow = flex.size_t(pan_fast_slow)
        self.all_background = np.array(all_background)
        self.roi_id = np.array(roi_id)
        self.all_data = np.array(all_data)
        self.all_sigmas = np.array(all_sigmas)
        # note rare chance for sigmas to be nan if the args of sqrt is below 0
        self.all_trusted = np.logical_and(np.array(all_trusted), ~np.isnan(all_sigmas))
        self.npix_total = len(all_data)
        self.all_fast = np.array(all_fast)
        self.all_slow = np.array(all_slow)
        self.simple_weights = 1/self.all_sigmas**2
        self.u_id = set(self.roi_id)

    def dump_gathered_to_refl(self, output_name, do_xyobs_sanity_check=False):
        """after running GatherFromExperiment, dump the gathered results
        (data, background etc) to a new reflection file which can then be used to run
        diffBragg without the raw data in the experiment (this exists mainly for portability, and
        unit tests)"""
        shoeboxes = []
        R = flex.reflection_table()
        for i_roi, i_ref in enumerate(self.refls_idx):
            roi_sel = self.roi_id==i_roi
            x1, x2, y1, y2 = self.rois[i_roi]
            roi_shape = y2-y1, x2-x1
            roi_img = self.all_data[roi_sel].reshape(roi_shape).astype(np.float32)  #NOTE this has already been converted to photon units
            roi_bg = self.all_background[roi_sel].reshape(roi_shape).astype(np.float32)

            sb = Shoebox((x1, x2, y1, y2, 0, 1))
            sb.allocate()
            sb.data = flex.float(np.ascontiguousarray(roi_img[None]))
            sb.background = flex.float(np.ascontiguousarray(roi_bg[None]))

            dials_mask = np.zeros(roi_img.shape).astype(np.int32)
            mask = self.all_trusted[roi_sel].reshape(roi_shape)
            dials_mask[mask] = dials_mask[mask] + MaskCode.Valid
            sb.mask = flex.int(np.ascontiguousarray(dials_mask[None]))

            # quick sanity test
            if do_xyobs_sanity_check:
                ref = self.refls[i_ref]
                x,y,_ = ref['xyzobs.px.value']
                assert x1 <= x <= x2, "exp %s; refl %d, %f %f %f" % (output_name, i_ref, x1,x,x2)
                assert y1 <= y <= y2, "exp %s; refl %d, %f %f %f" % (output_name, i_ref, y1,y,y2)

            R.extend(self.refls[i_ref: i_ref+1])
            shoeboxes.append(sb)

        R['shoebox'] = flex.shoebox(shoeboxes)
        R.as_file(output_name)

    def SimulatorFromExperiment(self, best=None):
        """optional best parameter is a single row of a pandas datafame containing the starting
        models, presumably optimized from a previous minimzation using this program"""

        ParameterType = RangedParameter if self.params.rescale_params else NormalParameter

        if best is not None:
            # set the crystal Umat (rotational displacement) and Bmat (unit cell)
            # Umatrix
            # NOTE: just set the best Amatrix here
            if self.params.apply_best_crystal_model:
                xax = col((-1, 0, 0))
                yax = col((0, -1, 0))
                zax = col((0, 0, -1))
                rotX,rotY,rotZ = best[["rotX", "rotY", "rotZ"]].values[0]
                RX = xax.axis_and_angle_as_r3_rotation_matrix(rotX, deg=False)
                RY = yax.axis_and_angle_as_r3_rotation_matrix(rotY, deg=False)
                RZ = zax.axis_and_angle_as_r3_rotation_matrix(rotZ, deg=False)
                M = RX * RY * RZ
                U = M * sqr(self.E.crystal.get_U())
                self.E.crystal.set_U(U)

                # Bmatrix:
                ucparam = best[["a","b","c","al","be","ga"]].values[0]
                ucman = utils.manager_from_params(ucparam)
                self.E.crystal.set_B(ucman.B_recipspace)

            ## TODO , currently need this anyway
            ucparam = best[["a","b","c","al","be","ga"]].values[0]
            ucman = utils.manager_from_params(ucparam)
            self.E.crystal.set_B(ucman.B_recipspace)

            # mosaic block
            self.params.init.Nabc = tuple(best.ncells.values[0])
            # scale factor
            self.params.init.G = best.spot_scales.values[0]

            if "detz_shift_mm" in list(best):
                self.params.init.detz_shift = best.detz_shift_mm.values[0]

        self.SIM = utils.simulator_from_expt_and_params(self.E, self.params)

        if self.params.spectrum_from_imageset:
            downsamp_spec(self.SIM, self.params, self.E)

        self.SIM.D.no_Nabc_scale = self.params.no_Nabc_scale
        self.SIM.num_xtals = self.params.number_of_xtals
        if self.params.eta_refine:
            self.SIM.umat_maker = AnisoUmats(num_random_samples=self.params.num_mosaic_blocks)
        self.SIM.Nabc_params = []
        self.SIM.RotXYZ_params = []
        self.SIM.Scale_params = []
        for i_xtal in range(self.SIM.num_xtals):
            for ii in range(3):
                p = ParameterType()
                p.sigma = self.params.sigmas.Nabc[ii]
                p.init = self.params.init.Nabc[ii]
                # set the mosaic block size
                p.minval = self.params.mins.Nabc[ii]
                p.maxval = self.params.maxs.Nabc[ii]
                self.SIM.Nabc_params.append(p)

                p = ParameterType()
                p.sigma = self.params.sigmas.RotXYZ[ii]
                p.init = 0
                p.minval = self.params.mins.RotXYZ[ii] * np.pi / 180.
                p.maxval = self.params.maxs.RotXYZ[ii] * np.pi / 180.
                self.SIM.RotXYZ_params.append(p)

            p = ParameterType()
            p.sigma = self.params.sigmas.G
            p.init = self.params.init.G
            p.minval = self.params.mins.G
            p.maxval = self.params.maxs.G
            self.SIM.Scale_params.append(p)

        ucell_man = utils.manager_from_crystal(self.E.crystal)
        ucell_vary_perc = self.params.ucell_edge_perc / 100.
        self.SIM.ucell_params = []
        for i_uc, (name, val) in enumerate(zip(ucell_man.variable_names, ucell_man.variables)):
            if "Ang" in name:
                minval = val - ucell_vary_perc * val
                maxval = val + ucell_vary_perc * val
            else:
                val_in_deg = val * 180 / np.pi
                minval = (val_in_deg - self.params.ucell_ang_abs) * np.pi / 180.
                maxval = (val_in_deg + self.params.ucell_ang_abs) * np.pi / 180.
            p = ParameterType()
            p.sigma = self.params.sigmas.ucell[i_uc]
            p.init = val
            p.minval = minval
            p.maxval = maxval
            if not self.params.quiet: print(
                "Unit cell variable %s (currently=%f) is bounded by %f and %f" % (name, val, minval, maxval))
            self.SIM.ucell_params.append(p)
        self.SIM.ucell_man = ucell_man

        p = ParameterType()
        p.init = self.params.init.detz_shift *1e-3
        p.sigma = self.params.sigmas.detz_shift
        p.minval = self.params.mins.detz_shift * 1e-3
        p.maxval = self.params.maxs.detz_shift * 1e-3
        self.SIM.DetZ_param = p

        # eta_max = self.params.maxs.eta
        # P.add("eta_a", value=0, min=0, max=eta_max * rad, vary=self.params.eta_refine)
        # P.add("eta_b", value=0, min=0, max=eta_max * rad, vary=self.params.eta_refine)
        # P.add("eta_c", value=0, min=0, max=eta_max * rad, vary=self.params.eta_refine)


    def Minimize(self, x0):
        target = TargetFunc(SIM=self.SIM, niter_per_J=self.params.niter_per_J)

        if self.params.method is None:
            method = "Nelder-Mead"
        else:
            method = self.params.method

        maxfev = self.params.nelder_mead_maxfev * self.npix_total

        at_min = target.at_minimum
        if self.params.quiet:
            at_min = target.at_minimum_quiet

        if method in ["L-BFGS-B", "BFGS", "CG", "dogleg", "SLSQP", "Newton-CG", "trust-ncg", "trust-krylov", "trust-exact", "trust-ncg"]:
            self.SIM.D.refine(ROTX_ID)
            self.SIM.D.refine(ROTY_ID)
            self.SIM.D.refine(ROTZ_ID)
            self.SIM.D.refine(NCELLS_ID)
            for i_ucell in range(len(self.SIM.ucell_man.variables)):
                self.SIM.D.refine(UCELL_ID_OFFSET + i_ucell)
            self.SIM.D.refine(DETZ_ID)

            args = (self.SIM, self.pan_fast_slow, self.all_data,
                    self.all_sigmas, self.all_trusted, self.all_background, not self.params.quiet, self.params, True)
            min_kwargs = {'args': args, "method": method, "jac": target.jac,
                          'hess': self.params.hess}
        else:
            args = (self.SIM, self.pan_fast_slow, self.all_data,
                    self.all_sigmas, self.all_trusted, self.all_background, not self.params.quiet, self.params, False)
            min_kwargs = {'args': args, "method": method, 'options':{'maxfev': maxfev}}

        if self.params.global_method=="basinhopping":
            out = basinhopping(target, x0,
                               niter=self.params.niter,
                               minimizer_kwargs=min_kwargs,
                               T=self.params.temp,
                               callback=at_min,
                               disp=not self.params.quiet,
                               stepsize=self.params.stepsize)
        else:
            bounds = [(-100,100)] * len(x0)  # TODO decide about bounds, usually x remains close to 1 during refinement
            print("Beginning the annealing process")
            args = min_kwargs.pop("args")
            if self.params.dual.no_local_search:
                compute_grads = args[-1]
                if compute_grads:
                    print("Warning, parameters setup to compute gradients, swicthing off because no_local_search=True")
                args = list(args)
                args[-1] = False  # switch off grad
                args = tuple(args)
            out = dual_annealing(target, bounds=bounds, args=args,
                                 no_local_search=self.params.dual.no_local_search,
                                 x0=x0,
                                 accept=self.params.dual.accept,
                                 visit=self.params.dual.visit,
                                 maxiter=self.params.niter,
                                 local_search_options=min_kwargs,
                                 callback=at_min)


        if not self.params.rescale_params:
            X = np.array(target.all_x)
            sig = 1 / np.std(X, 0)
            sig2 = sig/ sig.sum()
            print("G", sig[0], sig2[0])
            print("rotX", sig[1], sig2[1])
            print("rotY", sig[2], sig2[2])
            print("rotZ", sig[3], sig2[3])
            print("Na", sig[4], sig2[4])
            print("Nb", sig[5], sig2[5])
            print("Nc", sig[6], sig2[6])
            for i_uc, name in enumerate(self.SIM.ucell_man.variable_names):
                print(name, sig[7+i_uc], sig2[7+i_uc])
            n = 7+len(self.SIM.ucell_man.variables)
            print("DetZ", sig[n], sig2[n])

        P = out.x
        return P



def model(x, SIM, pfs, verbose=True, compute_grad=True):

    verbose = False
    num_per_xtal_params = SIM.num_xtals * 7
    n_ucell_param = len(SIM.ucell_man.variables)
    n_detector_param = 1

    assert n_ucell_param+num_per_xtal_params+n_detector_param == len(x)
    params_per_xtal = np.array_split(x[:num_per_xtal_params], SIM.num_xtals)

    # get the unit cell variables
    unitcell_var_reparam = x[num_per_xtal_params:num_per_xtal_params+n_ucell_param]
    unitcell_variables = [SIM.ucell_params[i].get_val(xval) for i, xval in enumerate(unitcell_var_reparam)]
    SIM.ucell_man.variables = unitcell_variables

    Bmatrix = SIM.ucell_man.B_recipspace
    SIM.D.Bmatrix = Bmatrix
    if compute_grad:
        for i_ucell in range(len(unitcell_variables)):
            SIM.D.set_ucell_derivative_matrix(
                i_ucell + UCELL_ID_OFFSET,
                SIM.ucell_man.derivative_matrices[i_ucell])

#   detector parameters
    x_shiftZ = x[num_per_xtal_params + n_ucell_param]
    shiftZ = SIM.DetZ_param.get_val(x_shiftZ)
    SIM.D.shift_origin_z(SIM.detector, shiftZ)

    npix = int(len(pfs) / 3)
    nparam = len(x)
    J = np.zeros((nparam, npix))  # note: order is: scale, rotX, rotY, rotZ, Na, Nb, Nc, ... (for each xtal), then ucell0, ucell1 , ucell2, .. detshift,
    model_pix = None
    for i_xtal in range(SIM.num_xtals):
        #SIM.D.raw_pixels_roi *= 0 #todo do i matter?
        scale_reparam, rotX_reparam, rotY_reparam, rotZ_reparam, \
        Na_reparam, Nb_reparam, Nc_reparam = params_per_xtal[i_xtal]

        rotX = SIM.RotXYZ_params[i_xtal * 3].get_val(rotX_reparam)
        rotY = SIM.RotXYZ_params[i_xtal * 3 + 1].get_val(rotY_reparam)
        rotZ = SIM.RotXYZ_params[i_xtal * 3 + 2].get_val(rotZ_reparam)

        ## update parameters:

        SIM.D.set_value(ROTX_ID, rotX)
        SIM.D.set_value(ROTY_ID, rotY)
        SIM.D.set_value(ROTZ_ID, rotZ)

        scale = SIM.Scale_params[i_xtal].get_val(scale_reparam)

        Na = SIM.Nabc_params[i_xtal * 3].get_val(Na_reparam)
        Nb = SIM.Nabc_params[i_xtal * 3 + 1].get_val(Nb_reparam)
        Nc = SIM.Nabc_params[i_xtal * 3 + 2].get_val(Nc_reparam)
        SIM.D.set_ncells_values(tuple([Na, Nb, Nc]))

        # SIM.D.verbose = 1
        # SIM.D.printout_pixel_fastslow = pfs[1],pfs[2]
        if verbose: print("\tXtal %d:" % i_xtal)
        if verbose: print("\tNcells=%f %f %f" % (Na, Nb, Nc))
        if verbose: print("\tspot scale=%f" % (scale))
        angles = tuple([x * 180 / np.pi for x in [rotX, rotY, rotZ]])
        if verbose: print("\trotXYZ= %f %f %f (degrees)" % angles)
        SIM.D.add_diffBragg_spots(pfs)

        pix = SIM.D.raw_pixels_roi[:npix]
        pix = pix.as_numpy_array()

        if model_pix is None:
            model_pix = scale*pix #SIM.D.raw_pixels_roi.as_numpy_array()[:npix]
        else:
            model_pix += scale*pix #SIM.D.raw_pixels_roi.as_numpy_array()[:npix]

        if compute_grad:
            scale_grad = model_pix / scale
            scale_grad = SIM.Scale_params[i_xtal].get_deriv(scale_reparam, scale_grad)
            J[7*i_xtal] += scale_grad

            rotX_grad = scale*SIM.D.get_derivative_pixels(ROTX_ID).as_numpy_array()[:npix]
            rotY_grad = scale*SIM.D.get_derivative_pixels(ROTY_ID).as_numpy_array()[:npix]
            rotZ_grad = scale*SIM.D.get_derivative_pixels(ROTZ_ID).as_numpy_array()[:npix]
            rotX_grad = SIM.RotXYZ_params[i_xtal*3].get_deriv(rotX_reparam, rotX_grad)
            rotY_grad = SIM.RotXYZ_params[i_xtal*3+1].get_deriv(rotY_reparam, rotY_grad)
            rotZ_grad = SIM.RotXYZ_params[i_xtal*3+2].get_deriv(rotZ_reparam, rotZ_grad)
            J[7*i_xtal + 1] += rotX_grad
            J[7*i_xtal + 2] += rotY_grad
            J[7*i_xtal + 3] += rotZ_grad

            Nabc_grad = SIM.D.get_ncells_derivative_pixels()
            #Na_grad = scale*SIM.D.get_Na_derivative_pixels()[:npix]
            Na_grad = scale*(Nabc_grad[0][:npix].as_numpy_array())
            Nb_grad = scale*(Nabc_grad[1][:npix].as_numpy_array())
            Nc_grad = scale*(Nabc_grad[2][:npix].as_numpy_array())

            #Na_grad, Nb_grad, Nc_grad = [scale*d.as_numpy_array()[:npix] for d in SIM.D.get_ncells_derivative_pixels()]
            Na_grad = SIM.Nabc_params[i_xtal * 3].get_deriv(Na_reparam, Na_grad)
            Nb_grad = SIM.Nabc_params[i_xtal * 3 + 1].get_deriv(Nb_reparam, Nb_grad)
            Nc_grad = SIM.Nabc_params[i_xtal * 3 + 2].get_deriv(Nc_reparam, Nc_grad)
            J[7*i_xtal + 4] += Na_grad
            J[7*i_xtal + 5] += Nb_grad
            J[7*i_xtal + 6] += Nc_grad

            for i_ucell in range(n_ucell_param):
                d = scale*SIM.D.get_derivative_pixels(UCELL_ID_OFFSET+i_ucell).as_numpy_array()[:npix]
                d = SIM.ucell_params[i_ucell].get_deriv(unitcell_var_reparam[i_ucell], d)
                J[7*SIM.num_xtals + i_ucell] += d

            d = SIM.D.get_derivative_pixels(DETZ_ID).as_numpy_array()[:npix]
            d = SIM.DetZ_param.get_deriv(x_shiftZ, d)
            J[7*SIM.num_xtals + n_ucell_param] += d


    #if verbose: print("\tunitcell= %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f" % SIM.ucell_man.unit_cell_parameters)
    return model_pix, J


def look_at_x(x, SIM):
    num_per_xtal_params = SIM.num_xtals * (7)
    n_ucell_param = len(SIM.ucell_man.variables)
    params_per_xtal = np.array_split(x[:num_per_xtal_params], SIM.num_xtals)

    for i_xtal in range(SIM.num_xtals):
        scale_reparam, rotX_reparam, rotY_reparam, rotZ_reparam, \
        Na_reparam, Nb_reparam, Nc_reparam = params_per_xtal[i_xtal]

        rotX = SIM.RotXYZ_params[i_xtal * 3].get_val(rotX_reparam)
        rotY = SIM.RotXYZ_params[i_xtal * 3 + 1].get_val(rotY_reparam)
        rotZ = SIM.RotXYZ_params[i_xtal * 3 + 2].get_val(rotZ_reparam)

        scale = SIM.Scale_params[i_xtal].get_val(scale_reparam)

        Na = SIM.Nabc_params[i_xtal * 3].get_val(Na_reparam)
        Nb = SIM.Nabc_params[i_xtal * 3 + 1].get_val(Nb_reparam)
        Nc = SIM.Nabc_params[i_xtal * 3 + 2].get_val(Nc_reparam)

        print("\tXtal %d:" % i_xtal)
        print("\tNcells=%f %f %f" % (Na, Nb, Nc))
        print("\tspot scale=%f" % (scale))
        angles = tuple([x * 180 / np.pi for x in [rotX, rotY, rotZ]])
        print("\trotXYZ= %f %f %f (degrees)" % angles)
    print("\tunitcell= %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f" % SIM.ucell_man.unit_cell_parameters)

    shiftZ = SIM.DetZ_param.get_val(x[num_per_xtal_params + n_ucell_param])
    print("\tshiftZ = %3.5f" % shiftZ)



def get_param_from_x(x, SIM):
    num_per_xtal_params = SIM.num_xtals * (7)
    n_ucell_param = len(SIM.ucell_man.variables)
    params_per_xtal = np.array_split(x[:num_per_xtal_params], SIM.num_xtals)
    unitcell_var_reparam = x[num_per_xtal_params:num_per_xtal_params+n_ucell_param]
    unitcell_variables = [SIM.ucell_params[i].get_val(xval) for i, xval in enumerate(unitcell_var_reparam)]
    SIM.ucell_man.variables = unitcell_variables
    a,b,c,al,be,ga = SIM.ucell_man.unit_cell_parameters

    detz_reparam = x[num_per_xtal_params + n_ucell_param]
    detz = SIM.DetZ_param.get_val(detz_reparam)


    #TODO generalize for n xtals
    i_xtal = 0

    scale_reparam, rotX_reparam, rotY_reparam, rotZ_reparam, \
        Na_reparam, Nb_reparam, Nc_reparam = params_per_xtal[i_xtal]

    rotX = SIM.RotXYZ_params[i_xtal * 3].get_val(rotX_reparam)
    rotY = SIM.RotXYZ_params[i_xtal * 3 + 1].get_val(rotY_reparam)
    rotZ = SIM.RotXYZ_params[i_xtal * 3 + 2].get_val(rotZ_reparam)

    scale = SIM.Scale_params[i_xtal].get_val(scale_reparam)

    Na = SIM.Nabc_params[i_xtal * 3].get_val(Na_reparam)
    Nb = SIM.Nabc_params[i_xtal * 3 + 1].get_val(Nb_reparam)
    Nc = SIM.Nabc_params[i_xtal * 3 + 2].get_val(Nc_reparam)

    return scale, rotX, rotY, rotZ, Na, Nb, Nc,a,b,c,al,be,ga, detz



class TargetFunc:
    def __init__(self, SIM, niter_per_J=1):
        self.niter_per_J = niter_per_J
        self.global_x = []
        self.all_x = []
        self.old_J = None
        self.old_model = None
        self.delta_x = None
        self.iteration = 0
        self.minima = []
        self.SIM = SIM

    def at_minimum_quiet(self, x, f, accept):
        self.iteration = 0
        self.all_x = []
        self.minima.append((f,x,accept))

    def at_minimum(self, x, f, accept):
        self.iteration = 0
        self.all_x = []
        look_at_x(x,self.SIM)
        self.minima.append((f,x,accept))

    def jac(self, x, *args):
        return self.g

    def __call__(self, x, *args, **kwargs):
        if self.all_x:
            self.delta_x = x - self.all_x[-1]
        update_terms = None
        if not self.iteration % (self.niter_per_J) == 0:
            update_terms = (self.delta_x, self.old_J, self.old_model)
        self.all_x.append(x)
        f, g, model, J = target_func(x, update_terms, *args, **kwargs)
        self.old_model = model
        self.old_J = J
        self.iteration += 1
        self.g = g
        return f
        #if g is not None:
        #    return f, g
        #else:
        #    return f

def target_func(x, udpate_terms, SIM, pfs, data, sigmas, trusted, background, verbose=True, params=None, compute_grad=True):

    #for i_x, xval in enumerate(x):
    #    all_x[xidx[i_x]] = xval

    if udpate_terms is not None:
        # if approximating the gradients, then fix the parameter refinment managers in diffBragg
        # so we dont waste time computing them
        _compute_grad = False
        SIM.D.fix(NCELLS_ID)
        SIM.D.fix(ROTX_ID)
        SIM.D.fix(ROTY_ID)
        SIM.D.fix(ROTZ_ID)
        for i_ucell in range(len(SIM.ucell_man.variables)):
            SIM.D.fix(UCELL_ID_OFFSET + i_ucell)
        SIM.D.fix(DETZ_ID)
    elif compute_grad:
        # actually compute the gradients
        _compute_grad = True
        SIM.D.let_loose(NCELLS_ID)
        SIM.D.let_loose(ROTX_ID)
        SIM.D.let_loose(ROTY_ID)
        SIM.D.let_loose(ROTZ_ID)
        for i_ucell in range(len(SIM.ucell_man.variables)):
            SIM.D.let_loose(UCELL_ID_OFFSET + i_ucell)
        SIM.D.let_loose(DETZ_ID)
    else:
        _compute_grad = False
    model_bragg, Jac = model(x, SIM, pfs, verbose=verbose, compute_grad=_compute_grad)

    if udpate_terms is not None:
        # try a Broyden update ?
        # https://people.duke.edu/~hpgavin/ce281/lm.pdf  equation 19
        delta_x, prev_J, prev_model_bragg = udpate_terms
        if prev_J is not None:
            delta_y = model_bragg - prev_model_bragg

            delta_J = (delta_y - np.dot(prev_J.T, delta_x))
            delta_J /= np.dot(delta_x,delta_x)
            Jac = prev_J + delta_J
    # Jac has shape of num_param x num_pix

    model_pix = model_bragg + background

    LL = params.use_likelihood_target

    #if not LL:
    W = 1/sigmas**2
    if LL:
        resid = (data - model_pix)  #minor technicality , to accommodate hand-written notes
    else:
        resid = (model_pix - data)

    G, rotX,rotY, rotZ, Na,Nb,Nc,a,b,c,al,be,ga,detz_shift = get_param_from_x(x, SIM)

    #TODO vectorize  / generalized framework for restraints
    ucvar = SIM.ucell_man.variables
    n_uc_param = len(ucvar)

    del_detz = detz_shift - params.centers.detz_shift

    G0 = params.centers.G
    delG = (G0-G)


    deg = 180 / np.pi
    rotX = deg*rotX
    rotY = deg*rotY
    rotZ = deg*rotZ
    rotX0,rotY0,rotZ0 = params.centers.RotXYZ
    Na0,Nb0,Nc0 = params.centers.Nabc
    del_rX = rotX0-rotX
    del_rY = rotY0-rotY
    del_rZ = rotZ0-rotZ

    del_Na = Na0 - Na
    del_Nb = Nb0 - Nb
    del_Nc = Nc0 - Nc

    if LL:
        sigma_rdout = params.refiner.sigma_r / params.refiner.adu_per_photon
        V = model_pix + sigma_rdout**2
        resid_square = resid**2
        fchi = (.5*(np.log(2*np.pi*V) + resid_square / V))[trusted].sum()   # negative log Likelihood target
        zscore = np.std((data - model_pix) / np.sqrt(V))
        # TODo make this a method the __call__ method of a class, and cache these terms
        Na_V = params.betas.Nabc[0]
        Nb_V = params.betas.Nabc[1]
        Nc_V = params.betas.Nabc[2]
        rx_V = params.betas.RotXYZ
        ry_V = params.betas.RotXYZ
        rz_V = params.betas.RotXYZ
        fN = .5*(np.log(2*np.pi*Na_V) + del_Na**2  / Na_V)
        fN += .5*(np.log(2*np.pi*Nb_V) + del_Nb**2  / Nb_V)
        fN += .5*(np.log(2*np.pi*Nc_V) + del_Nc**2  / Nc_V)

        frot = .5*(np.log(2*np.pi*rx_V) + del_rX**2  / rx_V)
        frot += .5*(np.log(2*np.pi*ry_V) + del_rY**2  / ry_V)
        frot += .5*(np.log(2*np.pi*rz_V) + del_rZ**2  / rz_V)

        G_V = params.betas.G
        fG = .5*(np.log(2*np.pi*G_V) + delG**2/G_V)


        detz_V = params.betas.detz_shift
        fz = .5*(np.log(2*np.pi*detz_V) + del_detz**2/detz_V)

        fucell = [0]*n_uc_param
        for i_ucell in range(n_uc_param):
            beta = params.betas.ucell[i_ucell]
            cent = params.centers.ucell[i_ucell]
            fucell[i_ucell] = .5*(np.log(2*np.pi*beta) + (cent-ucvar[i_ucell])**2/beta)


    else:
        fchi = (resid[trusted] ** 2 * W[trusted]).sum()   # weighted least squares target
        zscore = 0
        fN = params.betas.Nabc[0]*(del_Na )**2 + \
             params.betas.Nabc[1]*(del_Nb )**2 + \
             params.betas.Nabc[2]*(del_Nc )**2
        frot = params.betas.RotXYZ*((del_rX)**2+ (del_rY)**2 + (del_rZ )**2)
        fG = params.betas.G*delG**2
        fucell = [0]*n_uc_param
        for i_ucell in range(n_uc_param):
            beta = params.betas.ucell[i_ucell]
            if beta == 0:
                continue
            cent = params.centers.ucell[i_ucell]
            fucell[i_ucell] += beta * (cent - ucvar[i_ucell]) ** 2

        fz = 0

    fucell = sum(fucell)  # TODO distinguish betweem edge terms and angle terms
    f = fchi + frot + fN + fG + fucell + fz
    chi = fchi / f *100
    rot = frot / f*100
    uc = fucell / f*100
    n = fN / f*100
    gg = fG / f *100
    zz = fz / f * 100.
    g = None
    gnorm = -1
    if compute_grad:
        if LL:
            grad_term = (0.5 /V * (1-2*resid - resid_square / V))[trusted]
        else:
            grad_term = (2*resid*W)[trusted]
        Jac_t = Jac[:,trusted]
        g = np.array([np.sum(grad_term*Jac_t[param_idx]) for param_idx in range(Jac_t.shape[0])])
        if LL:
            g[0] += SIM.Scale_params[0].get_deriv(x[0], -delG / G_V)
            g[1] += SIM.RotXYZ_params[0].get_deriv(x[1], -del_rX / rx_V)
            g[2] += SIM.RotXYZ_params[1].get_deriv(x[2], -del_rY / ry_V)
            g[3] += SIM.RotXYZ_params[2].get_deriv(x[3], -del_rZ / rz_V)
            g[4] += SIM.Nabc_params[0].get_deriv(x[4], -del_Na / Na_V)
            g[5] += SIM.Nabc_params[1].get_deriv(x[5], del_Nb / Nb_V)
            g[6] += SIM.Nabc_params[2].get_deriv(x[6], -del_Nc / Nc_V)
            for i_uc in range(n_uc_param):
                beta = params.betas.ucell[i_uc]
                del_uc = params.centers.ucell[i_uc] - ucvar[i_uc]
                g[7+i_uc] += SIM.ucell_params[i_uc].get_deriv(x[7+i_uc], -del_uc / beta)
            g[7+n_uc_param] += SIM.DetZ_param.get_deriv(x[7+n_uc_param], -del_detz/detz_V)

        else:
            # TODO apply change of variable correction, as done for Likelihood restraint gradients above
            ber = params.betas.RotXYZ
            g[0] += -2*params.betas.G*delG
            g[1] += -ber*2*deg*del_rX
            g[2] += -ber*2*deg*del_rY
            g[3] += -ber*2*deg*del_rZ
            g[4] += -params.betas.Nabc[0]*2*del_Na
            g[5] += -params.betas.Nabc[1]*2*del_Nb
            g[6] += -params.betas.Nabc[2]*2*del_Nc
            for i_uc in range(n_uc_param):
                beta = params.betas.ucell[i_uc]
                if beta == 0:
                    continue
                del_uc = params.centers.ucell[i_uc] - ucvar[i_uc]
                g[7+i_uc] += -2*beta*del_uc
            #TODO detz gradient update for detz restraint,
        gnorm = np.linalg.norm(g)

    if verbose:
        print("F=%10.7g Z=%10.7g (chi: %.1f%%, rot: %.1f%% N: %.1f%%, G: %.1f%%, uc: %.1f%%, detz: %.1f%%), |g|=%10.7g" \
              % (f, zscore, chi, rot, n, gg, uc,zz,gnorm))

    return f, g, model_bragg, Jac


def refine(exp, ref, params, spec=None, gpu_device=None):
    if gpu_device is None:
        gpu_device = 0
    params.simulator.spectrum.filename = spec
    Modeler = DataModeler(params)
    if params.load_data_from_refls:
        Modeler.GatherFromReflectionTable(exp, ref)
    else:
        assert Modeler.GatherFromExperiment(exp, ref)

    Modeler.SimulatorFromExperiment()

    Modeler.SIM.D.device_Id = gpu_device

    # initial parameters (all set to 1, 7 parameters (scale, rotXYZ, Ncells_abc) per crystal (sausage) and then the unit cell parameters
    n_per_xtal = 7  # number of parameters that change per "shot" crystal, if multiple crystals per shot
    nparam = n_per_xtal * Modeler.SIM.num_xtals + len(Modeler.SIM.ucell_man.variables) + 1
    if params.rescale_params:
        x0 = [1] * nparam
    else:
        x0 = [np.nan] * nparam
        for i_xtal in range(Modeler.SIM.num_xtals):
            x0[n_per_xtal * i_xtal] = Modeler.SIM.Scale_params[i_xtal].init
            x0[n_per_xtal * i_xtal + 1] = Modeler.SIM.RotXYZ_params[3 * i_xtal].init
            x0[n_per_xtal * i_xtal + 2] = Modeler.SIM.RotXYZ_params[3 * i_xtal + 1].init
            x0[n_per_xtal * i_xtal + 3] = Modeler.SIM.RotXYZ_params[3 * i_xtal + 2].init
            x0[n_per_xtal * i_xtal + 4] = Modeler.SIM.Nabc_params[3 * i_xtal].init
            x0[n_per_xtal * i_xtal + 5] = Modeler.SIM.Nabc_params[3 * i_xtal + 1].init
            x0[n_per_xtal * i_xtal + 6] = Modeler.SIM.Nabc_params[3 * i_xtal + 2].init

        nucell = len(Modeler.SIM.ucell_man.variables)
        for i_ucell in range(nucell):
            x0[n_per_xtal * Modeler.SIM.num_xtals + i_ucell] = Modeler.SIM.ucell_params[i_ucell].init
        x0[n_per_xtal * Modeler.SIM.num_xtals + nucell] = Modeler.SIM.DetZ_param.init

    x = Modeler.Minimize(x0)
    best_model, _ = model(x, Modeler.SIM, Modeler.pan_fast_slow, compute_grad=False)

    new_crystal = update_crystal_from_x(Modeler.SIM, x)
    new_exp = deepcopy(Modeler.E)
    new_exp.crystal = new_crystal

    try:
        new_exp.beam.set_wavelength(Modeler.SIM.dxtbx_spec.get_weighted_wavelength())
    except:pass
    # if we strip the thickness from the detector, then update it here:
    #new_exp.detector. shift Z mm
    new_det = update_detector_from_x(Modeler.SIM, x)
    new_exp.detector = new_det

    new_refl = get_new_xycalcs(Modeler, best_model, new_exp)

    Modeler.clean_up()

    return new_exp, new_refl


def update_detector_from_x(SIM, x):
    scale, rotX, rotY, rotZ, Na, Nb, Nc, a, b, c, al, be, ga, detz_shift = get_param_from_x(x, SIM)
    detz_shift_mm = detz_shift*1e3
    det = SIM.detector
    det = utils.shift_panelZ(det, detz_shift_mm)
    return det


def update_crystal_from_x(SIM, x):
    scale, rotX, rotY, rotZ, Na, Nb, Nc, a, b, c, al, be, ga, detz_shift = get_param_from_x(x, SIM)

    xax = col((-1, 0, 0))
    yax = col((0, -1, 0))
    zax = col((0, 0, -1))
    ## update parameters:
    RX = xax.axis_and_angle_as_r3_rotation_matrix(rotX, deg=False)
    RY = yax.axis_and_angle_as_r3_rotation_matrix(rotY, deg=False)
    RZ = zax.axis_and_angle_as_r3_rotation_matrix(rotZ, deg=False)
    M = RX * RY * RZ
    U = M * sqr(SIM.crystal.dxtbx_crystal.get_U())
    new_C = deepcopy(SIM.crystal.dxtbx_crystal)
    new_C.set_U(U)

    ucparam = a, b, c, al, be, ga
    ucman = utils.manager_from_params(ucparam)
    new_C.set_B(ucman.B_recipspace)

    return new_C


def get_new_xycalcs(Modeler, best_model, new_exp):
    _,_,_, bragg_subimg = get_data_model_pairs(Modeler.rois, Modeler.pids, Modeler.roi_id, best_model, Modeler.all_data, background=Modeler.all_background)
    new_refls = deepcopy(Modeler.refls)
    new_refls['dials.xyzcal.px'] = deepcopy(new_refls['xyzcal.px'])
    new_refls['dials.xyzcal.mm'] = deepcopy(new_refls['xyzcal.mm'])
    new_refls['dials.xyzobs.mm.value'] = deepcopy(new_refls['xyzobs.mm.value'])
    new_xycalcs = flex.vec3_double(len(Modeler.refls), (np.nan, np.nan, np.nan))
    new_xycalcs_mm = flex.vec3_double(len(Modeler.refls), (np.nan, np.nan, np.nan))
    new_xyobs_mm = flex.vec3_double(len(Modeler.refls), (np.nan, np.nan, np.nan))
    for i_roi in range(len(bragg_subimg)):

        ref_idx = Modeler.refls_idx[i_roi]

        if np.any(bragg_subimg[i_roi] > 0):
            I = bragg_subimg[i_roi]
            Y, X = np.indices(bragg_subimg[i_roi].shape)
            x1, _, y1, _ = Modeler.rois[i_roi]
            X += x1
            Y += y1
            Isum = I.sum()
            xcom = (X * I).sum() / Isum + .5
            ycom = (Y * I).sum() / Isum + .5
            com = xcom, ycom, 0

            pid = Modeler.pids[i_roi]
            assert pid == new_refls[ref_idx]['panel']
            panel = new_exp.detector[pid]
            xmm, ymm = panel.pixel_to_millimeter((xcom, ycom))
            com_mm = xmm, ymm, 0
            xobs, yobs, _ = new_refls[ref_idx]["xyzobs.px.value"]
            xobs_mm, yobs_mm = panel.pixel_to_millimeter((xobs, yobs))
            obs_com_mm = xobs_mm, yobs_mm, 0

            new_xycalcs[ref_idx] = com
            new_xycalcs_mm[ref_idx] = com_mm
            new_xyobs_mm[ref_idx] = obs_com_mm

    new_refls["xyzcal.px"] = new_xycalcs
    new_refls["xyzcal.mm"] = new_xycalcs_mm
    new_refls["xyzobs.mm.value"] = new_xyobs_mm

    if Modeler.params.filter_unpredicted_refls_in_output:
        sel = [not np.isnan(x) for x,_,_ in new_refls['xyzcal.px']]
        new_refls = new_refls.select(flex.bool(sel))

    return new_refls


def get_data_model_pairs(rois, pids, roi_id, best_model, all_data, strong_flags=None, background=None):
    all_dat_img, all_mod_img = [], []
    all_strong = []
    all_bragg = []
    for i_roi in range(len(rois)):
        x1, x2, y1, y2 = rois[i_roi]
        mod = best_model[roi_id == i_roi].reshape((y2 - y1, x2 - x1))
        if strong_flags is not None:
            strong = strong_flags[roi_id == i_roi].reshape((y2 - y1, x2 - x1))
            all_strong.append(strong)
        else:
            all_strong.append(None)

        # dat = img_data[pid, y1:y2, x1:x2]
        dat = all_data[roi_id == i_roi].reshape((y2 - y1, x2 - x1))
        all_dat_img.append(dat)
        if background is not None:
            bg = background[roi_id==i_roi].reshape((y2-y1, x2-x1))
            # assume mod does not contain background
            all_bragg.append(mod)
            all_mod_img.append(mod+bg)
        else:  # assume mod contains background
            all_mod_img.append(mod)
            all_bragg.append(None)
        # print("Roi %d, max in data=%f, max in model=%f" %(i_roi, dat.max(), mod.max()))
    return all_dat_img, all_mod_img, all_strong, all_bragg


# set the X-ray spectra for this shot
def downsamp_spec(SIM, params, expt):
    SIM.dxtbx_spec = expt.imageset.get_spectrum(0)
    spec_en = SIM.dxtbx_spec.get_energies_eV()
    spec_wt = SIM.dxtbx_spec.get_weights()
    # ---- downsample the spectrum
    method2_param = {"filt_freq": params.downsamp_spec.filt_freq,
                     "filt_order": params.downsamp_spec.filt_order,
                     "tail": params.downsamp_spec.tail,
                     "delta_en": params.downsamp_spec.delta_en}
    downsamp_en, downsamp_wt = downsample_spectrum(spec_en.as_numpy_array(),
                                                   spec_wt.as_numpy_array(),
                                                   method=2, method2_param=method2_param)

    stride = params.simulator.spectrum.stride
    if stride > len(downsamp_en) or stride == 0:
        raise ValueError("Incorrect value for pinkstride")
    downsamp_en = downsamp_en[::stride]
    downsamp_wt = downsamp_wt[::stride]
    downsamp_wt = downsamp_wt / sum(downsamp_wt) * params.simulator.total_flux

    downsamp_wave = utils.ENERGY_CONV / downsamp_en
    SIM.beam.spectrum = list(zip(downsamp_wave, downsamp_wt))
    # the nanoBragg beam has an xray_beams property that is used internally in diffBragg
    SIM.D.xray_beams = SIM.beam.xray_beams