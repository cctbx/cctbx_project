from __future__ import division
import glob
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--kokkos", action="store_true")
parser.add_argument("--nolog", action="store_true")
parser.add_argument("--readout", type=float, default=3)
parser.add_argument("--scale", type=float, default=1)
parser.add_argument("--perturb", choices=["G", "Nabc"], type=str, nargs="+", default=None)
parser.add_argument("--plot", action="store_true")
parser.add_argument("--beta", default=None, type=float)
parser.add_argument("--sigmaFhkl", default=1, type=float)
parser.add_argument("--sigmaG", default=1, type=float)
parser.add_argument("--maxiter", default=None, type=int)
parser.add_argument("--geo", action="store_true")
args = parser.parse_args()
import os

if args.kokkos:
    os.environ["DIFFBRAGG_USE_KOKKOS"]="1"
from simtbx.diffBragg.utils import find_diffBragg_instances
from simtbx.diffBragg.device import DeviceWrapper
with DeviceWrapper(0) as _:

    import logging
    import sys
    import pandas
    import numpy as np

    from cctbx import miller
    from dials.array_family import flex
    from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
    from simtbx.nanoBragg.sim_data import SimData
    from simtbx.diffBragg import hopper_utils
    from simtbx.diffBragg import utils
    from simtbx.diffBragg.hopper_ensemble_utils import load_inputs
    from dxtbx.model import Experiment
    from simtbx.nanoBragg import make_imageset
    from simtbx.diffBragg.phil import hopper_phil, philz
    from libtbx.phil import parse

    phil_scope = parse(hopper_phil+philz)


    p65_cryst = {'__id__': 'crystal',
                 'real_space_a': (43.32309880004587, 25.5289818883498, 60.49634260901813),
                 'real_space_b': (34.201635357808115, -38.82573591182249, -59.255697149884924),
                 'real_space_c': (41.42476391176581, 229.70849483520402, -126.60059788183489),
                 'space_group_hall_symbol': ' P 65 2 (x,y,z+1/12)',
                 'ML_half_mosaicity_deg': 0.06671930026192037,
                 'ML_domain_size_ang': 6349.223840307989}
    from dxtbx.model.crystal import CrystalFactory
    p65_C = CrystalFactory.from_dict(p65_cryst)
    ucell = p65_C.get_unit_cell().parameters()
    symbol = p65_C.get_space_group().info().type().lookup_symbol()

    # Setup the simulation and create a realistic image
    # with background and noise
    # <><><><><><><><><><><><><><><><><><><><><><><><><>
    nbcryst = NBcrystal()
    nbcryst.dxtbx_crystal = p65_C
    nbcryst.thick_mm = 0.005
    nbcryst.isotropic_ncells = False
    NCELLS_GT = 12,12,11
    nbcryst.Ncells_abc = NCELLS_GT
    nbcryst.space_group = "P6522"
    ma = utils.make_miller_array(symbol, ucell, d_min=1.5)
    np.random.seed(0)
    new_data = ma.d_spacings().data()*10

    ma = miller.array(ma.set(), new_data).set_observation_type_xray_amplitude()
    ma_map = {h:v for h,v in zip(ma.indices(), ma.data())}
    nbcryst.miller_array = ma
    assert ma.is_xray_amplitude_array()

    SIM = SimData(use_default_crystal=False)
    shape = 1000, 1001
    detdist = 140
    SIM.detector = SimData.simple_detector(detdist, 0.1, shape)
    SIM.crystal = nbcryst
    SIM.instantiate_diffBragg(oversample=1, auto_set_spotscale=True, default_F=0)

    # test the code for computing the acerage structure factor intensity with resolution
    # (this is why we set the structure factor data to be the same as the resolution (x10)
    num_dspace_bins = 10
    SIM.set_dspace_binning(num_dspace_bins, verbose=True)
    dspace_bins = SIM.D.dspace_bins
    ave_I_cell = SIM.D.ave_I_cell()[0]
    assert len(ave_I_cell) == num_dspace_bins
    assert len(dspace_bins) == num_dspace_bins + 1
    aves = []
    dspaces = []
    for i, (d1,d2) in enumerate(zip(dspace_bins, dspace_bins[1:])):
        ave_val = np.sqrt(ave_I_cell[i]) / 10.
        if not args.geo: assert d1 < ave_val < d2
        aves.append(ave_val)
        dspaces.append(.5*(d1+d2))

    print("1 0 7: ", ma.value_at_index((1,0,7)))
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
        plt.figure()
        plt.plot(dspaces, aves)
        plt.xlabel("Angstrom")
        plt.show()
    SIM.D.raw_pixels *= 0
    #pfs = 0,270,175
    #utils.show_diffBragg_state(SIM.D, pfs)
    SIM.D.raw_pixels *= 0

    P = phil_scope.extract()
    P.debug_mode=True
    E = Experiment()

    P.init.G = SIM.D.spot_scale
    E.crystal = p65_C

    P.init.Nabc = SIM.crystal.Ncells_abc
    P.init.detz_shift = 0

    E.detector = SIM.detector
    E.beam = SIM.D.beam
    E.imageset = make_imageset([img], E.beam, E.detector)
    refls = utils.refls_from_sims([spots], E.detector, E.beam, thresh=18)
    print("%d REFLS" % len(refls))
    utils.refls_to_q(refls, E.detector, E.beam, update_table=True)
    utils.refls_to_hkl(refls, E.detector, E.beam, E.crystal, update_table=True)

    P.roi.shoebox_size = 10
    P.roi.allow_overlapping_spots = True
    P.relative_tilt = False
    P.roi.fit_tilt = False
    P.roi.pad_shoebox_for_background_estimation=10
    P.roi.reject_edge_reflections = False
    P.refiner.sigma_r = SIM.D.readout_noise_adu
    P.refiner.adu_per_photon = SIM.D.quantum_gain
    P.simulator.init_scale = 1
    P.simulator.beam.size_mm = SIM.beam.size_mm
    P.simulator.oversample = SIM.D.oversample
    P.simulator.total_flux = SIM.D.flux
    P.use_restraints = False


    mset = ma.set()
    ma_map_keys, ma_map_values = list(ma_map.keys()), np.array(list(ma_map.values()))
    ma_map_values2 = np.random.normal(ma_map_values, scale=args.scale*ma_map_values)
    bad_map = {h:v for h,v in zip(ma_map_keys, ma_map_values2)}
    if args.scale ==0:
        assert np.allclose(ma_map_values, ma_map_values2)

    new_amps = flex.double()
    for h in mset.indices():
        amp = bad_map[h]
        new_amps.append(amp)

    ma2 = miller.array(mset, new_amps).set_observation_type_xray_amplitude()

    ma2_map = {h:v for h,v in zip(ma2.indices(), ma2.data())}
    name = "hopper_refine_Fhkl.mtz"
    print("1 0 7: ", ma2.value_at_index((1,0,7)))
    assert ma2.is_xray_amplitude_array()
    ma2.as_mtz_dataset(column_root_label="F").mtz_object().write(name)
    P.simulator.structure_factors.mtz_name = name
    P.simulator.structure_factors.mtz_column = "F(+),F(-)"
    P.logging.parameters=False
    P.method="L-BFGS-B"
    P.ftol = 1e-10
    P.space_group = symbol
    P.fix.Fhkl = False
    P.betas.Fhkl = args.beta
    P.fix.G = True
    P.types.G = "positive"
    P.centers.G = SIM.D.spot_scale*2
    P.betas.G=1e8
    P.use_restraints = args.beta is not None
    P.sigmas.G = args.sigmaG
    P.sigmas.Fhkl = args.sigmaFhkl
    if args.perturb is not None and "G" in args.perturb:
        P.fix.G = False
        P.init.G = SIM.D.spot_scale*10
        #P.maxs.G = SIM.D.spot_scale*100
    P.use_geometric_mean_Fhkl = args.geo
    P.fix.ucell=True
    P.fix.RotXYZ=True
    P.fix.Nabc=True
    if args.perturb is not None and "Nabc" in args.perturb:
        P.fix.Nabc = False
        P.init.Nabc = 20,20,18
    P.fix.detz_shift=True

    if not args.nolog:
        h = logging.StreamHandler(sys.stdout)
        logging.basicConfig(level=logging.DEBUG, handlers=[h])
    #del SIM.D

    P.outdir="_temp_fhkl_refine"
    if args.maxiter is not None:
        P.lbfgs_maxiter = args.maxiter
    P.record_device_timings = True
    Eopt,_, Mod,SIM_from_hopper, x = hopper_utils.refine(E, refls, P, return_modeler=True, free_mem=False)
    SIM_from_hopper.D.show_timings(0)

    logging.disable()
    print("\nResults\n<><><><><><>")

    Mod.exper_name = "dummie.expt"
    Mod.refl_name = "dummie.refl"
    Mod.save_up(x, SIM_from_hopper)

    # we can track the dominant hkls in each shoebox occuring within the diffBragg model
    #count_stats = utils.track_fhkl(Mod)
    #
    #
    #main_hkls = []
    #for i_roi in count_stats:
    #    stats = count_stats[i_roi]
    #    stats = sorted( list(stats.items()), key=lambda x: x[1])
    #    main_hkl, frac = stats[-1]
    #    print(main_hkl, frac)
    #    main_hkls.append(main_hkl)

    # this should agree with what we put into diffBragg in the reflection tables
    main_hkls_from_refls = utils.map_hkl_list(list(Mod.refls["miller_index"]), symbol=P.space_group)
    #assert len(set(main_hkls)) == len(set(main_hkls_from_refls))
    # good.

    # Now, this should also agree with the refined fhkl values, stored in the data table
    # the modeler save_up method creates an output file containing the refined fhkl values

    fnames = glob.glob("%s/Fhkl_scale/rank*/*.npz" % P.outdir)
    assert len(fnames)==1
    fhkl_f= fnames[0]
    fhkl_dat = np.load(fhkl_f)
    asu = list(map(tuple,fhkl_dat['asu_hkl']))
    asu_corrections = fhkl_dat['scale_fac']
    asu_corrections_var = fhkl_dat['scale_var']
    asu_is_nominal = fhkl_dat["is_nominal_hkl"]
    scale = {h:s for h,s in zip(asu, asu_corrections)}
    scale_var = {h:s for h,s in zip(asu, asu_corrections_var)}
    is_nominal = {h:is_nom for h,is_nom in zip(asu,  asu_is_nominal)}

    num_not_nominal = 0
    # TODO: figure out why some main_hkls are missing from is_nominal
    for hkl in main_hkls_from_refls:
        if hkl not in is_nominal:
            continue
        if not is_nominal[hkl]:
            num_not_nominal += 1

    nominal_hkl_corrections = {h:s for h,s in zip(asu, asu_corrections) if is_nominal[h]}
    not_nominal_hkl_corrections = {h:s for h,s in zip(asu, asu_corrections) if not is_nominal[h]}

    nominal_hkl_init = {h:1 for h in asu if is_nominal[h]}
    not_nominal_hkl_init = {h:1 for h in asu if not is_nominal[h]}


    def compute_r_factor_with_gt(corrections):
        gt_data = flex.double()
        opt_data = flex.double()
        flx_hkls = flex.miller_index()
        for hkl, scale in corrections.items():
            gt_amp = ma_map[hkl]
            gt_data.append(gt_amp)

            opt_amp = np.sqrt(scale) * ma2_map[hkl]
            opt_data.append(opt_amp)

            h,k,l = map(int, hkl)
            flx_hkls.append((h,k,l))
        mset = ma.miller_set(flx_hkls, ma.anomalous_flag())
        gt_arr = miller.array(mset, gt_data).set_observation_type_xray_amplitude()
        opt_arr = miller.array(mset, opt_data).set_observation_type_xray_amplitude()
        return gt_arr.r1_factor(opt_arr)


    r1_nominal_init = compute_r_factor_with_gt(nominal_hkl_init)
    r1_nominal = compute_r_factor_with_gt(nominal_hkl_corrections)

    r1_not_nominal_init = compute_r_factor_with_gt(not_nominal_hkl_init)
    r1_not_nominal = compute_r_factor_with_gt(not_nominal_hkl_corrections)

    print("\nResults\n<><><><><><>")
    print("For the dominant HKLs within each modeled shoebox (e.g. those with indexed reflections)")
    print("initial R1 factor=%.2f%%" % (r1_nominal_init*100))
    print("optimized R1 factor=%.2f%%" % (r1_nominal*100))

    diffs = []
    all_opts = []
    all_gts = []
    dsp_map = {d:val for d,val in zip(ma.d_spacings().indices(), ma.d_spacings().data())}
    ds = []
    for hkl in nominal_hkl_corrections:
        dsp = dsp_map[hkl]
        gt_val = ma_map[hkl]
        opt_val = np.sqrt(nominal_hkl_corrections[hkl]) * ma2_map[hkl]
        d = abs(gt_val - opt_val) / gt_val * 100
        diffs.append(d)

        all_opts.append(opt_val)
        all_gts.append(gt_val)
        ds.append(dsp)

    print("mean percent diff", np.mean(diffs))

    assert r1_nominal < r1_nominal_init, "r1_nom=%f, r1_not_nom=%f" %(r1_nominal, r1_nominal_init)
    assert r1_nominal < 0.04

    # test hopper_ensemble_refiner using this one shot
    # dump the refinement data to the reflection table format (e.g. the pixel data and background estimates)
    input_refl = os.path.join(P.outdir, "input_data.refl")
    Mod.dump_gathered_to_refl(input_refl)
    df = pandas.read_pickle("%s/pandas/rank0/stage1_dummie_0.pkl"% P.outdir)
    refl_col = "input_refls"
    df[refl_col] = [input_refl]
    P.refiner.load_data_from_refl = True
    P.refiner.check_expt_format = False

    #from simtbx.diffBragg import mpi_logger
    #P.logging.rank0_level="high"
    #mpi_logger.setup_logging_from_params(P)
    modelers = load_inputs(df, P, exper_key="opt_exp_name", refls_key=refl_col)
    modelers.outdir=P.outdir
    modelers.prep_for_refinement()
    print("Minimizing using hopper_ensemble_utils...")
    modelers.Minimize(save=True)
    if modelers.SIM.D.record_timings:
        modelers.SIM.D.show_timings(MPI_RANK=0)
    print("Done!")

    from iotbx.reflection_file_reader import any_reflection_file
    opt_F = any_reflection_file("_temp_fhkl_refine/optimized_channel0.mtz").as_miller_arrays()[0]
    opt_map = {h:v for h,v in zip(opt_F.indices(), opt_F.data())}
    hcommon = set(ma_map).intersection(opt_map)

    mset_common = ma.miller_set(flex.miller_index(list(hcommon)), ma.anomalous_flag())
    ma_vals = flex.double([ma_map[h] for h in hcommon])
    opt_vals = flex.double([opt_map[h] for h in hcommon])

    ma_common = miller.array(mset_common, ma_vals).set_observation_type_xray_amplitude()
    opt_common = miller.array(mset_common, opt_vals).set_observation_type_xray_amplitude()
    r1 = ma_common.r1_factor(opt_common)
    assert r1 < 0.04

    print("OK")
    del modelers.SIM.D
    del SIM_from_hopper.D
    for name in find_diffBragg_instances(globals()): del globals()[name]
