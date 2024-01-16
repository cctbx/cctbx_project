from __future__ import division
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--kokkos", action="store_true")
parser.add_argument("--curvatures", action='store_true')
parser.add_argument("--readout", type=float, default=0)
parser.add_argument("--perturb", choices=["G", "Nabc", "detz_shift", "crystal", "eta", "spec"], type=str, nargs="+", default=["crystal"] )
parser.add_argument("--eta", type=float, nargs=3, default=[.2, .1, .4])
parser.add_argument("--plot", action="store_true", help="shows the ground truth image")
parser.add_argument("--typeG", choices=["ranged", "positive"], default="ranged", type=str,  help="shows the ground truth image")
parser.add_argument("--typeNabc", choices=["ranged", "positive"], default="ranged", type=str,  help="shows the ground truth image")
parser.add_argument("--cmdlineHopper", action="store_true", help="test the command line program simtbx/command_line/hopper.py")
args = parser.parse_args()
name = "hopper_refine_%s" % "-".join(args.perturb)
import os

if args.kokkos:
    os.environ["DIFFBRAGG_USE_KOKKOS"]="1"
from simtbx.diffBragg.utils import find_diffBragg_instances
from simtbx.diffBragg.device import DeviceWrapper
with DeviceWrapper(0) as _:

    from dxtbx.model.crystal import Crystal
    from cctbx import uctbx
    from scitbx.matrix import sqr, rec, col
    import numpy as np
    from scipy.spatial.transform import Rotation
    from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
    from simtbx.nanoBragg.sim_data import SimData
    from simtbx.diffBragg import utils
    from scipy.signal import windows
    from dxtbx.model import Experiment
    from simtbx.nanoBragg import make_imageset
    from simtbx.diffBragg.phil import hopper_phil, philz
    from libtbx.phil import parse

    phil_scope = parse(hopper_phil+philz)

    ucell = (55, 65, 75, 90, 95, 90)
    ucell2 = (55.1, 65.2, 74.9, 90, 94.9, 90)
    symbol = "P121"

    # generate a random raotation
    rotation = Rotation.random(num=1, random_state=100)[0]
    Q = rec(rotation.as_quat(), n=(4, 1))
    rot_ang, rot_axis = Q.unit_quaternion_as_axis_and_angle()

    # generate a small perturbation rotation
    np.random.seed(1)
    perturb_rot_axis = np.random.random(3)
    perturb_rot_axis /= np.linalg.norm(perturb_rot_axis)
    perturb_rot_ang = 0.15  # degree random perturbtation

    # make the ground truth crystal:
    a_real, b_real, c_real = sqr(uctbx.unit_cell(ucell).orthogonalization_matrix()).transpose().as_list_of_lists()
    C = Crystal(a_real, b_real, c_real, symbol)
    C.rotate_around_origin(rot_axis, rot_ang)

    a2_real, b2_real, c2_real = sqr(uctbx.unit_cell(ucell2).orthogonalization_matrix()).transpose().as_list_of_lists()
    C2 = Crystal(a2_real, b2_real, c2_real, symbol)
    C2.rotate_around_origin(rot_axis, rot_ang)
    assert np.allclose(C2.get_U(), C.get_U())
    C2.rotate_around_origin(col(perturb_rot_axis), perturb_rot_ang)

    # Setup the simulation and create a realistic image
    # with background and noise
    # <><><><><><><><><><><><><><><><><><><><><><><><><>
    nbcryst = NBcrystal()
    nbcryst.dxtbx_crystal = C   # simulate ground truth
    nbcryst.thick_mm = 0.1
    nbcryst.isotropic_ncells = False
    if "eta" in args.perturb:
        nbcryst.n_mos_domains = 1000
        ETA_ABC_GT = args.eta
        nbcryst.anisotropic_mos_spread_deg = ETA_ABC_GT
        NCELLS_GT = 12,12,11
    else:
        NCELLS_GT = 12,12,11
    nbcryst.Ncells_abc = NCELLS_GT

    SIM = SimData(use_default_crystal=True)
    if "spec" in args.perturb:
        # initialize the simulator
        spec = SIM.beam.spectrum
        total_flux = spec[0][1]
        wave = spec[0][0]

        en = utils.ENERGY_CONV / wave
        delta_en = 1.5
        ens_gt = np.arange(en - 5, en + 6, delta_en)
        waves_gt = utils.ENERGY_CONV / ens_gt
        num_energies = len(ens_gt)
        fluxes_gt = np.ones(num_energies) * total_flux / num_energies
        fluxes_gt = fluxes_gt*windows.hann(num_energies)
        fluxes_gt /= fluxes_gt.sum()
        fluxes_gt *= total_flux

        spectrum_GT = list(zip(waves_gt, fluxes_gt))
        gt_lambda0 = waves_gt[0]
        gt_lambda1 = waves_gt[1] - waves_gt[0]
        spec_idx = np.arange(num_energies)
        assert np.allclose(waves_gt, gt_lambda0 + spec_idx*gt_lambda1)

        lam0 = np.random.normal(gt_lambda0, gt_lambda0 * 0.002)
        lam1 = np.random.normal(gt_lambda1, abs(gt_lambda1) * 0.002)
        waves_perturbed = lam0 + spec_idx * lam1
        print("ENERGY TRUTH=%.4f" % (utils.ENERGY_CONV / gt_lambda0))
        print("ENERGY PERTURBED=%.4f" % (utils.ENERGY_CONV / lam0))
        perturbed_spec = list(zip(waves_perturbed, fluxes_gt))
        SIM.beam.spectrum = spectrum_GT

    #SIM.detector = SimData.simple_detector(150, 0.1, (513, 512))
    if "eta" in args.perturb:
        shape = 513*3, 512*3
        #detdist = 70
    else:
        shape = 513, 512
    detdist = 150
    SIM.detector = SimData.simple_detector(detdist, 0.1, shape)
    SIM.crystal = nbcryst
    SIM.instantiate_diffBragg(oversample=0, auto_set_spotscale=True)
    SIM.D.default_F = 0
    SIM.D.F000 = 0
    SIM.D.progress_meter = False
    SIM.water_path_mm = 0.005
    SIM.air_path_mm = 0.1
    SIM.add_air = True
    SIM.add_Water = True
    SIM.include_noise = True
    SIM.D.verbose = 2
    SIM.D.add_diffBragg_spots()
    SIM.D.verbose = 0
    spots = SIM.D.raw_pixels.as_numpy_array()
    SIM._add_background()
    SIM.D.readout_noise_adu=args.readout
    SIM._add_noise()

    # This is the ground truth image:
    img = SIM.D.raw_pixels.as_numpy_array()
    if args.plot:
        import pylab as plt
        plt.imshow(img, vmax=100)
        plt.title("Ground truth image")
        plt.show()

    cbf_name = name + ".cbf"
    if args.cmdlineHopper:
        from IPython import embed
        embed()
        SIM.D.to_cbf(cbf_name)

    SIM.D.raw_pixels *= 0

    P = phil_scope.extract()
    E = Experiment()

    GT_spot_scale = SIM.D.spot_scale
    if "G" in args.perturb:
        P.init.G = GT_spot_scale*10
    else:
        P.init.G = GT_spot_scale

    P.types.G = args.typeG
    P.types.Nabc = args.typeNabc

    if "crystal" in args.perturb:
        E.crystal = C2
    else:
        E.crystal = C

    if "Nabc" in args.perturb:
        P.init.Nabc = 20,20,20
    else:
        P.init.Nabc = SIM.crystal.Ncells_abc

    if "detz_shift" in args.perturb:
        P.init.detz_shift = 1
    else:
        P.init.detz_shift = 0

    if "eta" in args.perturb:
        P.init.eta_abc = [0.12, 0.13, 0.14]
        P.simulator.crystal.num_mosaicity_samples = 250  # in practive, the number of mosaic domains we model should be smaller than whats in the crystal .. .
        P.simulator.crystal.has_isotropic_mosaicity = False
        P.fix.eta_abc = False

    if "spec" in args.perturb:
        P.fix.spec = False
        P.init.spec = [0,1]
        P.fix.Nabc=True
        P.fix.G=True
        P.fix.RotXYZ=True
        P.fix.ucell = True
        P.fix.detz_shift = True
        P.ftol=1e-15

    if args.perturb == ["detz_shift"]:
        P.fix.detz_shift = False
        P.fix.ucell=True
        P.fix.Nabc=True
        P.fix.G=True
        P.fix.RotXYZ=True

    E.detector = SIM.detector
    E.beam = SIM.D.beam
    E.imageset = make_imageset([img], E.beam, E.detector)
    #refls = utils.refls_from_sims([img], E.detector, E.beam, thresh=18)
    refls = utils.refls_from_sims([spots], E.detector, E.beam, thresh=18)
    print("%d REFLS" % len(refls))
    utils.refls_to_q(refls, E.detector, E.beam, update_table=True)
    utils.refls_to_hkl(refls, E.detector, E.beam, E.crystal, update_table=True)

    P.roi.shoebox_size = 20
    P.relative_tilt = False
    P.roi.fit_tilt = False
    P.roi.pad_shoebox_for_background_estimation=10
    P.roi.reject_edge_reflections = False
    P.refiner.sigma_r = SIM.D.readout_noise_adu
    P.refiner.adu_per_photon = SIM.D.quantum_gain
    P.simulator.init_scale = 1 #SIM.D.spot_scale
    P.simulator.beam.size_mm = SIM.beam.size_mm
    P.simulator.total_flux = SIM.D.flux
    P.use_restraints = False
    mtz_name = name +".mtz"
    SIM.crystal.miller_array.as_mtz_dataset(column_root_label="F").mtz_object().write(mtz_name)
    P.simulator.structure_factors.mtz_name = mtz_name
    P.simulator.structure_factors.mtz_column = "F(+),F(-)"
    P.niter = 0
    P.sigmas.RotXYZ = [1,1,1]
    P.logging.parameters=True
    P.niter_per_J = 1
    P.method="L-BFGS-B"
    P.ftol = 1e-10
    if "eta" in args.perturb:
        P.ftol=1e-8
    #P.method="Nelder-Mead"
    #P.fix.G = True
    #P.fix.Nabc =True
    #P.fix.detz_shift=True

    import logging
    import sys
    h = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.DEBUG, handlers=[h])

    from simtbx.diffBragg import hopper_utils
    spec = None
    if "spec" in args.perturb:
        spec = "tst_hopper_refine_spec.lam"
        wave, wt = map(np.array, zip(*perturbed_spec))
        utils.save_spectra_file(spec, wave, wt)

    if args.cmdlineHopper:
        from dxtbx.model import ExperimentList
        el_name = "%s.expt" % name
        import_cmd = "dials.import %s output.experiments=%s" % (cbf_name, el_name)
        os.system(import_cmd)
        # add the crystal to the imported expt
        El = ExperimentList.from_file(el_name)
        El[0].crystal = E.crystal
        refl_name = "%s.refl" % name
        El.as_file(el_name)
        refls.as_file(refl_name)
        hopper_input_lst = "%s.lst" % name
        with open(hopper_input_lst, "w") as o:
            o.write("%s %s\n" % (os.path.abspath(el_name), os.path.abspath(refl_name)))

        # save the above modified phil params to a file to be read in by hopper
        phil_file = name+ ".phil"
        modified_phil = phil_scope.format(python_object=P)
        with open(phil_file, 'w') as o:
            modified_phil.show(o)
        outdir = name + ".outdir"
        cmd = "hopper %s exp_ref_spec_file=%s outdir=%s logging.rank0_level=high" % (phil_file, hopper_input_lst, outdir)
        os.system(cmd)
        # TODO open the pandas output file and optimized expt in outdir and verify the optimized parameters are similar to ground
        exit()

    P.record_device_timings = True
    Eopt,_, Mod, SIM_used_by_hopper, x = hopper_utils.refine(E, refls, P, spec=spec, return_modeler=True)
    if SIM_used_by_hopper.D.record_timings:
        SIM_used_by_hopper.D.show_timings(MPI_RANK=0)

    G, rotX,rotY, rotZ, Na,Nb,Nc,_,_,_,_,_,_,_,_,_,a,b,c,al,be,ga,detz_shift = hopper_utils.get_param_from_x(x, Mod)
    eta_abc_opt = hopper_utils.get_mosaicity_from_x(x, Mod, SIM_used_by_hopper)

    print("Na, Nb, Nc= %f %f %f" % (Na, Nb, Nc))
    print("eta_abc optimized:", eta_abc_opt)

    # check crystal
    Copt = Eopt.crystal
    misset, misset_init = utils.compare_with_ground_truth(*C.get_real_space_vectors(), dxcryst_models=[Copt, E.crystal], symbol=symbol)
    print(misset_init, "init misset with ground truth")
    print(misset, "misset with ground truth")
    if "detz_shift" in args.perturb or "spec" in args.perturb:
        assert misset < 0.007, misset
    else:
        assert misset < 0.005, misset

    # check mosaic domain
    assert all(np.subtract(NCELLS_GT, [Na,Nb,Nc]) < 0.2), "%d, %d, %d" % (Na,Nb,Nb)

    # check spot scale
    perc_diff_G = abs(GT_spot_scale - G)/ GT_spot_scale * 100
    print("spot scale gt: %f; spot scale opt: %f; percent diff: %f %%" % (GT_spot_scale, G, perc_diff_G))
    max_Gperc = 1
    if "eta" in args.perturb:
        max_Gperc = 2
    assert perc_diff_G < max_Gperc, perc_diff_G

    # check detz
    print("detdist shift %f (should be 0)" % detz_shift)
    assert detz_shift < 0.2, detz_shift

    ucell_diff_init = np.abs(np.subtract(ucell , ucell2))
    ucell_diff = np.abs(np.subtract(ucell , Copt.get_unit_cell().parameters()))

    init_dev, init_dev_ang = ucell_diff_init[:3].sum(), ucell_diff_init[-3:].sum()
    dev, dev_ang = ucell_diff[:3].sum(), ucell_diff[-3:].sum()
    print("initial ucell dev: %f Angstrom; %f degree" % (init_dev, init_dev_ang))
    print("optimized ucell dev: %f Angstrom; %f degree" % (dev, dev_ang))
    assert dev_ang < init_dev_ang and dev_ang < 0.025, "init: %f curr: %f" % (init_dev_ang, dev_ang)
    if "detz_shift" not in args.perturb:
        assert dev < init_dev and dev < 0.025, "init: %f  curr: %f" % (init_dev, dev)

    if "eta" in args.perturb:
        print("eta_abc GT:", ETA_ABC_GT)
        u = np.array(eta_abc_opt)
        v = np.array(ETA_ABC_GT)
        perc_diff = np.abs(u-v) / v * 100.
        assert np.all(perc_diff < 22)  # this is acceptable for now, as we simulated with 5000 blocks, yet modeled with 600
    print("OK")

    if "spec" in args.perturb:
        p0 = Mod.P["lambda_offset"]
        p1 = Mod.P["lambda_scale"]
        coef = p0.get_val(x[p0.xpos]), p1.get_val(x[p1.xpos])
        waves_refined = coef[0] + coef[1] * waves_perturbed
        fluxsum = sum(fluxes_gt)
        en_ref_com = utils.ENERGY_CONV / (sum(fluxes_gt * waves_refined) / fluxsum)
        en_com = utils.ENERGY_CONV / (sum(fluxes_gt * waves_gt) / fluxsum)
        en_init_com = utils.ENERGY_CONV / (sum(fluxes_gt * waves_perturbed) / fluxsum)

        print("Before refinement: COM energy=%f" % en_init_com)
        print("AFTER refinement: COM energy=%f" % en_ref_com)
        print("Ground truth COM energy = %f" % en_com)
        assert abs(en_ref_com - en_com) < 1

    del SIM_used_by_hopper.D
    for name in find_diffBragg_instances(globals()): del globals()[name]
