from __future__ import division
from cctbx import omz
import cctbx.omz.dev
from cctbx.array_family import flex
import libtbx.phil.command_line
from libtbx.test_utils import approx_equal
from libtbx import easy_run
from libtbx import easy_pickle
from libtbx.utils import date_and_time, user_plus_sys_time
import libtbx.load_env
from libtbx import Auto
from cStringIO import StringIO
import traceback
import sys, os
op = os.path

def get_master_phil(
      max_atoms=99,
      f_calc_options_algorithm="*direct fft",
      bulk_solvent_correction=False):
  return omz.dev.get_master_phil(
    iteration_limit=100,
    show_distances_threshold=0.5,
    bulk_solvent_correction=bulk_solvent_correction,
    grads_mean_sq_threshold=1e-6,
    f_calc_options_algorithm=f_calc_options_algorithm,
    additional_phil_string="""\
      max_atoms = %(max_atoms)s
        .type = int
      f_obs_f_calc_fan_outliers = *remove keep
        .type = choice
        .optional = False
      use_f_calc_as_f_obs = False
        .type = bool
      reset_u_iso = 0.05
        .type = float
      sites_mod_short = True
        .type = bool
      optimizers = *dev ls_simple ls_lm shelxl_fm shelxl_cg shelx76
        .type = choice(multi=True)
      ls_simple_iterations = 12
        .type = int
      shelxl_wght = None
        .type = str
        .help = '''
          SHELX-97 Manual 7-31:
            Refinement against F2 requires different weights to refinement
            against F; in particular, making all the weights equal ('unit
            weights'), although useful in the initial stages of refinement
            against F, is NEVER a sensible option for F2.'''
      shelxl_reset_sigmas = None
        .type = float
      shelxl_fm_iterations = 12
        .type = int
      shelxl_cg_iterations = 12
        .type = int
      shelx76_iterations = 12
        .type = int
      apply_iteration_limit_to_all = False
        .type = bool
      keep_tmp_files = False
        .type = bool
      export_refined = False
        .type = bool
      pickle_refined_dir = None
        .type = str
      wdir_root = None
        .type = str
      sorting_of_pickle_files = *down up
        .type = choice
        .optional = True
      random_subset {
        size = None
          .type = int
        seed = 0
          .type = int
      }
      tardy_samples {
        iq = None
          .type = int
        qmin = -180
          .type = float
        qmax = 180
          .type = float
        qstep = 3
          .type = float
      }
""" % vars())

def shelxl_weights_a_b(fo_sq, sigmas, fc_sq, osf_sq, a, b):
  assert sigmas.size() == fo_sq.size()
  assert fc_sq.size() == fo_sq.size()
  from cctbx.xray.targets.tst_shelxl_wght_ls import calc_w
  assert sigmas.all_ge(0.01)
  return calc_w(
    wa=a, wb=b, i_obs=fo_sq, i_sig=sigmas, i_calc=fc_sq, k=osf_sq**0.5)

def shelxl_weights(fo_sq, sigmas, fc_sq, osf_sq, shelxl_wght):
  if (shelxl_wght is None):
    shelxl_wght = ""
  vals = [float(s) for s in shelxl_wght.split()]
  assert len(vals) <= 6
  a, b, c, d, e, f = vals + [0.1, 0, 0, 0, 0, 0.33333][len(vals):]
  assert c == 0
  assert d == 0
  assert e == 0
  assert f == 0.33333
  return shelxl_weights_a_b(fo_sq, sigmas, fc_sq, osf_sq, a, b)

def show_cc_r1(
      params,
      label,
      f_obs,
      xray_structure=None,
      fc_abs=None,
      scale_factor=Auto):
  assert [xray_structure, fc_abs].count(None) == 1
  if (fc_abs is None):
    p = params.f_calc_options
    fc_abs = f_obs.structure_factors_from_scatterers(
      xray_structure=xray_structure,
      algorithm=p.algorithm,
      cos_sin_table=p.cos_sin_table).f_calc().amplitudes()
  corr = flex.linear_correlation(x=f_obs.data(), y=fc_abs.data())
  assert corr.is_well_defined()
  cc = corr.coefficient()
  r1 = f_obs.r1_factor(
    other=fc_abs, scale_factor=scale_factor, assume_index_matching=True)
  print "%-12s cc, r1: %.4f %.4f" % (label, cc, r1)
  sys.stdout.flush()
  return fc_abs, cc, r1

def run_smtbx_ls(mode, cod_id, i_obs, f_obs, xray_structure, params):
  import smtbx.refinement
  fo_sq = i_obs
  assert fo_sq.sigmas() is not None
  sel = (fo_sq.data() == 0) & (fo_sq.sigmas() == 0)
  fo_sq = fo_sq.select(~sel)
  fo_sq.select(fo_sq.sigmas() <= 0).show_array()
  assert fo_sq.sigmas().all_gt(0)
  if (1): # work around bug currently in smtbx weighting scheme implementation
    fo_sq = fo_sq.customized_copy(sigmas=flex.double(fo_sq.data().size(), 1))
  xobs = fo_sq.as_xray_observations()
  tm = user_plus_sys_time()
  rm = smtbx.refinement.model(
    fo_sq=xobs,
    xray_structure=xray_structure,
    constraints=[],
    restraints_manager=smtbx.refinement.restraints.manager(),
    weighting_scheme=smtbx.refinement.least_squares.unit_weighting())
  ls = rm.least_squares()
  if (mode == "simple"):
    for i_cycle in xrange(params.ls_simple_iterations):
      ls.build_up()
      try:
        ls.solve_and_step_forward()
      except RuntimeError, e:
        if (str(e).find("cholesky.failure") <= 0): raise
        print 'Aborting run_smtbx_ls("simple"): cholesky.failure: %s' \
          % cod_id
        break
      for sc in xray_structure.scatterers():
        if (sc.u_iso <= 0 or sc.u_iso > 1):
          sc.u_iso = 0.05
      show_cc_r1(params, "ls%02d" % (i_cycle+1), f_obs, xray_structure)
    tm.show_elapsed(prefix="time smtbx_ls_simple_iterations: ")
  elif (mode == "lm"):
    from scitbx.lstbx import normal_eqns_solving
    thresh = 1e-6
    try:
      cycles = normal_eqns_solving.levenberg_marquardt_iterations(
        ls,
        gradient_threshold=thresh,
        step_threshold=thresh,
        tau=1e-7)
    except RuntimeError, e:
      if (not str(e).startswith(
            "cctbx::adptbx::debye_waller_factor_exp: max_arg exceeded")):
        raise
      print 'Aborting run_smtbx_ls("lm"):' \
        ' debye_waller_factor_exp failure: %s' % cod_id
    show_cc_r1(params, "smtbx_lm", f_obs, xray_structure)
    tm.show_elapsed(prefix="time levenberg_marquardt_iterations: ")
  else:
    raise RuntimeError('Unknown run_smtbx_ls(mode="%s")' % mode)

def remove_tmp_files(file_names):
  for fn in file_names:
    if (op.isfile(fn)):
      os.remove(fn)
    assert not op.exists(fn)

def run_shelxl(
      mode,
      cod_id,
      i_obs,
      f_obs,
      xray_structure,
      params,
      reference_structure,
      expected_n_refinable_parameters):
  if (mode == "fm"):
    if (params.apply_iteration_limit_to_all):
      fm_cycles = params.iteration_limit
    else:
      fm_cycles = params.shelxl_fm_iterations
    cg_cycles = None
  elif (mode == "cg"):
    fm_cycles = None
    if (params.apply_iteration_limit_to_all):
      cg_cycles = params.iteration_limit
    else:
      cg_cycles = params.shelxl_cg_iterations
  else:
    raise RuntimeError("Unknown mode: " + mode)
  cwd_orig = os.getcwd()
  wdir = "wdir_%s_shelxl_%s_%s" % (cod_id, mode, os.getpid())
  if (params.wdir_root is not None):
    wdir = op.join(params.wdir_root, wdir)
  wdir_is_new = False
  if (not op.isdir(wdir)):
    os.mkdir(wdir)
    wdir_is_new = True
  remove_wdir = False
  try:
    os.chdir(wdir)
    tmp_file_names = ["tmp.ins", "tmp.hkl", "tmp.res", "tmp.lst"]
    remove_tmp_files(tmp_file_names)
    fo_sq = i_obs
    if (params.shelxl_reset_sigmas):
      fo_sq = fo_sq.customized_copy(
        sigmas=flex.double(fo_sq.indices().size(), params.shelxl_reset_sigmas))
    import iotbx.shelx
    open("tmp.ins", "w").writelines(iotbx.shelx.writer.generator(
      xray_structure=xray_structure,
      data_are_intensities=True,
      title="cod_id=%s mode=%s" % (cod_id, mode),
      wavelength=fo_sq.minimum_wavelength_based_on_d_min(),
      full_matrix_least_squares_cycles=fm_cycles,
      conjugate_gradient_least_squares_cycles=cg_cycles,
      weighting_scheme_params=params.shelxl_wght,
      sort_scatterers=False))
    fo_sq.export_as_shelx_hklf(file_object=open("tmp.hkl", "w"))
    import iotbx.shelx.hklf
    fo_sq = iotbx.shelx.hklf.reader(file_name="tmp.hkl") \
      .as_miller_arrays(crystal_symmetry=fo_sq)[0]
    buffers = easy_run.fully_buffered("shelxl tmp")
    buffers.raise_if_errors()
    refinement_unstable = False
    for line in buffers.stdout_lines:
      if (line.find("** REFINEMENT UNSTABLE **") >= 0):
        refinement_unstable = True
        print "Aborted: shelxl %s refinement unstable: %s" % (mode, cod_id)
        break
    res = open("tmp.res").read()
    try:
      refined = xray_structure.from_shelx(
        file=StringIO(res),
        min_distance_sym_equiv=0,
        strictly_shelxl=False)
    except iotbx.shelx.error, e:
      if (str(e).find("scatterer parameter") < 0):
        raise
      print "Aborted: shelxl %s refinement apparently unstable: %s" % (
        mode, cod_id)
      refined = None
    if (refined is not None):
      assert refined.crystal_symmetry().is_similar_symmetry(
        xray_structure)
      for sc,rsc in zip(xray_structure.scatterers(), refined.scatterers()):
        assert rsc.label == sc.label
        assert approx_equal(rsc.occupancy, sc.weight(), 1e-4)
        rsc.occupancy = sc.occupancy # XXX bug in res file reader
      def check_special_positions():
        result = True
        uc = xray_structure.unit_cell()
        sstab = xray_structure.site_symmetry_table()
        for i_sc in xray_structure.special_position_indices():
          sc = refined.scatterers()[i_sc]
          site_symmetry = sstab.get(i_sc)
          assert not site_symmetry.is_point_group_1()
          site_special = site_symmetry.special_op() * sc.site
          d = uc.mod_short_distance(sc.site, site_special)
          if (d > 1e-3):
            print "site moved off special position:"
            print "  %s" % sc.label
            print "    shelxl res: %11.6f %11.6f %11.6f" % sc.site
            print "    special_op: %11.6f %11.6f %11.6f" % site_special
            print "    distance moved: %.3f" % d
            result = False
        return result
      assert check_special_positions()
      xray_structure.replace_scatterers(refined.scatterers())
      res_osf = None
      res_hkl_count = None
      res_r1 = None
      res_n_parameters = None
      res_n_restraints = None
      for line in res.splitlines():
        if (line.startswith("FVAR ")):
          flds = line.split()
          assert len(flds) == 2
          res_osf = float(flds[1])
          continue
        if (not line.startswith("REM ")): continue
        assert not refinement_unstable
        if (line.startswith("REM R1 =")):
          flds = line.split()
          assert len(flds) == 15
          res_hkl_count = int(flds[13])
          res_r1 = float(flds[10])
        elif (line.find(" parameters refined ") >= 0):
          assert line.endswith(" restraints")
          flds = line.split()
          assert len(flds) == 7
          res_n_parameters = int(flds[1])
          res_n_restraints = int(flds[-2])
      if (not refinement_unstable):
        assert res_osf is not None
        assert res_hkl_count is not None
        assert res_r1 is not None
        assert res_n_parameters is not None
        assert res_n_restraints is not None
        #
        assert res_hkl_count == fo_sq.indices().size()
        def raise_unexpected_restraints(n_expected):
          raise RuntimeError(
            "Unexpected number of SHELXL restraints: %d (vs. %d expected)" % (
              res_n_restraints, n_expected))
        if (mode == "fm"):
          n_caos = fo_sq.space_group_info() \
            .number_of_continuous_allowed_origin_shifts()
          if (res_n_restraints != n_caos):
            sg_symbol = str(fo_sq.space_group_info())
            if (sg_symbol in ["P 63 m c", "P 63 c m"]):
              assert n_caos == 1
              assert res_n_restraints == 0
              print "INFO: SHELXL restraint count incorrect? code_code:", \
                cod_id
            else:
              raise_unexpected_restraints(n_caos)
        elif (mode == "cg"):
          if (res_n_restraints != 0):
            raise_unexpected_restraints(0)
        else:
          raise RuntimeError("Unknown mode: " + mode)
        assert res_n_parameters == expected_n_refinable_parameters + 1
        fc_abs, _, r1_fvar = show_cc_r1(
          params, "fvar_"+mode, f_obs, xray_structure, scale_factor=res_osf)
        r1_diff = r1_fvar - res_r1
        print "R1 recomputed - shelxl_%s.res: %.4f - %.4f = %.4f %s" % (
          mode, r1_fvar, res_r1, r1_diff, cod_id)
        if (abs(r1_diff) > 0.01):
          raise RuntimeError("R1 MISMATCH %s" % cod_id)
        _, _, r1_auto = show_cc_r1(
          params, "shelxl_"+mode, f_obs, fc_abs=fc_abs)
        print "R1 FVAR-Auto %s: %.4f" % (cod_id, r1_fvar - r1_auto)
        #
        lst_r1 = None
        lst_wr2 = None
        for line in open("tmp.lst").read().splitlines():
          l = line.strip()
          if (l.startswith("R1 = ")):
            lst_r1 = float(l.split()[9])
          elif (l.startswith("wR2 = ") and l.endswith(" for all data")):
            lst_wr2 = float(l.replace(","," ").split()[2])
        assert lst_r1 is not None
        assert lst_wr2 is not None
        assert lst_r1 == res_r1
        #
        fc_sq = fc_abs.f_as_f_sq()
        weights = shelxl_weights(
          fo_sq=fo_sq.data(),
          sigmas=fo_sq.sigmas(),
          fc_sq=fc_sq.data(),
          osf_sq=res_osf**2,
          shelxl_wght=params.shelxl_wght)
        num = flex.sum(
          weights * flex.pow2(fo_sq.data() / res_osf**2 - fc_sq.data()))
        den = flex.sum(
          weights * flex.pow2(fo_sq.data() / res_osf**2))
        assert den != 0
        wr2 = (num / den)**0.5
        wr2_diff = wr2 - lst_wr2
        if (abs(wr2_diff) > 0.01):
          info = " significantly different"
        else:
          info = ""
        print "wR2 recomputed - shelxl_%s.lst: %.4f - %.4f = %.4f %s%s" % (
          mode, wr2, lst_wr2, wr2_diff, cod_id, info)
        if (abs(wr2_diff) / max(lst_wr2, wr2) > 0.2):
          raise RuntimeError("wR2 MISMATCH %s" % cod_id)
    if (not params.keep_tmp_files):
      remove_tmp_files(tmp_file_names)
      remove_wdir = wdir_is_new
  finally:
    os.chdir(cwd_orig)
    if (remove_wdir):
      try: os.rmdir(wdir)
      except Exception: pass

def run_shelx76(
      cod_id,
      f_obs,
      xray_structure,
      fvars,
      encoded_sites,
      params,
      reference_structure):
  if (params.apply_iteration_limit_to_all):
    ls_cycles = params.iteration_limit
  else:
    ls_cycles = params.shelx76_iterations
  cwd_orig = os.getcwd()
  wdir = "wdir_%s_shelx76_%s" % (cod_id, os.getpid())
  if (params.wdir_root is not None):
    wdir = op.join(params.wdir_root, wdir)
  wdir_is_new = False
  if (not op.isdir(wdir)):
    os.mkdir(wdir)
    wdir_is_new = True
  remove_wdir = False
  try:
    os.chdir(wdir)
    tmp_file_names = [
      "tmp.ins", "tmp.lst", "fort.2", "fort.3", "fort.4", "fort.7"]
    remove_tmp_files(tmp_file_names)
    assert not op.exists("tmp.ins")
    from cctbx.development import run_shelx76
    run_shelx76.write_shelx76_ls(
      f_obs=f_obs,
      xray_structure=xray_structure,
      fvars=fvars,
      encoded_sites=encoded_sites,
      l_s_parameters=str(ls_cycles))
    assert op.exists("tmp.ins")
    buffers = easy_run.fully_buffered("shelx76 < tmp.ins > tmp.lst")
    buffers.raise_if_errors_or_output()
    lst = open("tmp.lst").read().splitlines()
    r_from_lst = None
    for line in lst:
      l = line.lstrip()
      if (l.startswith("R = ")):
        print l
        flds = l.split()
        assert len(flds) == 12
        if (flds[2].lower() == "nan"):
          print "Aborted: shelx76 refinement apparently unstable: %s" % (
            cod_id)
          r_from_lst = "nan"
          break
        r_from_lst = float(flds[2])
    assert r_from_lst is not None
    if (r_from_lst != "nan"):
      print "%-12s cc, r1: None %.4f" % ("shelx76", r_from_lst)
      if (not params.keep_tmp_files):
        remove_tmp_files(tmp_file_names)
        remove_wdir = wdir_is_new
  finally:
    os.chdir(cwd_orig)
    if (remove_wdir):
      try: os.rmdir(wdir)
      except Exception: pass

def process(params, pickle_file_name):
  cod_id = op.basename(pickle_file_name).split(".",1)[0]
  print "cod_id:", cod_id
  c_obs, structure_prep, edge_list = easy_pickle.load(
    file_name=pickle_file_name)
  changes = structure_prep.make_scatterer_labels_shelx_compatible_in_place()
  if (params.sites_mod_short):
    structure_prep = structure_prep.sites_mod_short()
  from iotbx.shelx import fvar_encoding
  structure_prep = \
    fvar_encoding.move_sites_if_necessary_for_shelx_fvar_encoding(
      xray_structure=structure_prep)
  structure_prep.show_summary().show_scatterers()
  if (len(changes) != 0):
    from libtbx.utils import plural_s
    print "INFO: %d atom name%s changed for compatibility with SHELXL:" \
      % plural_s(len(changes))
    for change in changes:
      print '  changed: "%s" -> "%s"' % change
  structure_prep.scattering_type_registry(table="it1992").show()
  fvar_encoding.dev_build_shelx76_fvars(structure_prep) # only an exercise
  print "."*79
  #
  if (len(params.optimizers) == 0):
    return
  #
  assert c_obs.is_xray_intensity_array() or c_obs.is_xray_amplitude_array()
  if (c_obs.is_xray_intensity_array()):
    i_obs = c_obs
    f_obs = c_obs.f_sq_as_f(algorithm="xtal_3_7")
  else:
    f_obs = c_obs
    i_obs = c_obs.f_as_f_sq(algorithm="shelxl")
  process_continue(
    params=params,
    cod_id=cod_id,
    c_obs=c_obs, i_obs=i_obs, f_obs=f_obs,
    structure_prep=structure_prep)

def process_continue(params, cod_id, c_obs, i_obs, f_obs, structure_prep):
  p = params.f_calc_options
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure=structure_prep,
    algorithm=p.algorithm,
    cos_sin_table=p.cos_sin_table).f_calc()
  sel = f_obs.f_obs_f_calc_fan_outlier_selection(f_calc=f_calc)
  assert sel is not None
  n_outliers = sel.count(True)
  if (n_outliers != 0):
    action = params.f_obs_f_calc_fan_outliers
    print "INFO: f_obs_f_calc_fan_outliers = %s: %d" % (action, n_outliers)
    if (action == "remove"):
      i_obs = i_obs.select(~sel)
      f_obs = f_obs.select(~sel)
  if (f_obs.anomalous_flag()):
    print "INFO: converting anomalous i+f_obs to non-anomalous."
    i_obs = i_obs.average_bijvoet_mates()
    f_obs = f_obs.average_bijvoet_mates()
  sel = ((i_obs.data() == 0) & (i_obs.sigmas() == 0)) \
      | ((f_obs.data() == 0) & (f_obs.sigmas() == 0))
  n_zero_d_and_s = sel.count(True)
  if (n_zero_d_and_s != 0):
    print "INFO: removing reflections with i+f_obs=0 and sigma=0:", \
      n_zero_d_and_s
    i_obs = i_obs.select(~sel)
    f_obs = f_obs.select(~sel)
  p = params.f_calc_options
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure=structure_prep,
    algorithm=p.algorithm,
    cos_sin_table=p.cos_sin_table).f_calc()
  if (params.use_f_calc_as_f_obs):
    print "INFO: using f_calc as i+f_obs"
    i_obs = f_calc.intensities().customized_copy(
      sigmas=flex.double(f_calc.indices().size(), 0.01))
    f_obs = f_calc.amplitudes().customized_copy(
      sigmas=flex.double(f_calc.indices().size(), 0.01))
  else:
    # scaling applied so that the data written in shelx hklf format
    # have sufficient significant digits, and FVAR is 1 (shelx76 seems
    # to be especially sensitive to FVAR >> 1)
    k = f_obs.scale_factor(f_calc=f_calc)
    assert k != 0
    s = 1/k**2
    print "INFO: scaling i_obs to f_calc by multiplying i_obs with: %.6g" % s
    i_obs = i_obs.apply_scaling(factor=s)
    s = 1/k
    print "INFO: scaling f_obs to f_calc by multiplying f_obs with: %.6g" % s
    f_obs = f_obs.apply_scaling(factor=s)
  def show(obs):
    obs.show_comprehensive_summary()
    from cod_select_and_pickle import \
      report_fraction_of_negative_observations_if_any as _
    _(cod_id, obs)
  if (c_obs.is_xray_intensity_array()):
    show(i_obs)
  else:
    show(f_obs)
  print "."*79
  #
  structure_work = structure_prep.deep_copy_scatterers()
  sel = structure_work.hd_selection()
  print "Removing hydrogen atoms:", sel.count(True)
  structure_work = structure_work.select(selection=~sel)
  sdt = params.show_distances_threshold
  if (sdt > 0):
    print "Distances smaller than %.6g A:" % sdt
    structure_work.show_distances(distance_cutoff=sdt)
    print "."*79
  #
  if (params.tardy_samples.iq is not None):
    from cctbx.omz import tardy_adaptor
    print
    tardy_adaptor.sample_e_pot(
      id_code=cod_id,
      f_obs=f_obs,
      xray_structure=structure_prep,
      edge_list=edge_list,
      params=params.tardy_samples)
    print
    return
  #
  from iotbx.shelx import fvar_encoding
  fvars, encoded_sites = fvar_encoding.dev_build_shelx76_fvars(structure_work)
  print "Number of FVARs for special position constraints:", len(fvars)-1
  print "."*79
  #
  show_cc_r1(params, "prep", f_obs, structure_prep)
  def cc_r1(label):
    show_cc_r1(params, label, f_obs, structure_work)
  cc_r1("no_h")
  structure_work.convert_to_isotropic()
  cc_r1("iso")
  structure_iso = structure_work.deep_copy_scatterers()
  #
  if (params.reset_u_iso is not None):
    structure_work.set_u_iso(value=params.reset_u_iso)
    cc_r1("setu")
  if (params.shake_sites_rmsd is not None):
    mt = flex.mersenne_twister(seed=0)
    structure_work.shift_sites_in_place(
      shift_length=params.shake_sites_rmsd,
      mersenne_twister=mt)
    print "rms difference after shift_sites_in_place: %.3f" \
      % structure_iso.rms_difference(structure_work)
    cc_r1("shift_xyz")
  #
  if (params.max_atoms is not None):
    n = structure_work.scatterers().size()
    if (n > params.max_atoms):
      print "Skipping refinement of large model: %d atoms COD %s" % (
        n, cod_id)
      return
  #
  structure_work.scatterers().flags_set_grads(state=False)
  for sc in structure_work.scatterers():
    sc.flags.set_grad_site(True)
    assert sc.flags.use_u_iso_only()
    sc.flags.set_grad_u_iso(True)
  n_refinable_parameters = structure_work.n_parameters(
    considering_site_symmetry_constraints=True)
  print "Number of refinable parameters:", n_refinable_parameters
  #
  if (params.iteration_limit < 1):
    return
  #
  if ("dev" not in params.optimizers):
    structure_dev = None
  else:
    structure_dev = structure_work.deep_copy_scatterers()
    omz.dev.refinement(
      i_obs=i_obs,
      f_obs=f_obs,
      xray_structure=structure_dev,
      params=params,
      reference_structure=structure_iso,
      expected_n_refinable_parameters=n_refinable_parameters,
      plot_samples_id=cod_id)
    show_cc_r1(params, "dev", f_obs, structure_dev)
    if (params.export_refined):
      file_name = "dev_%s_%s_%s.pdb" % (
        params.target_type, params.target_obs_type.lower(), cod_id)
      open(file_name, "w").write(structure_dev.as_pdb_file(
        remarks=[file_name]))
    if (params.pickle_refined_dir is not None):
      easy_pickle.dump(
        file_name=op.join(params.pickle_refined_dir, cod_id+".pickle"),
        obj=(c_obs, structure_dev, None))
      print >> open("%s/qi_%s" % (params.pickle_refined_dir, cod_id), "w"), (
        structure_dev.scatterers().size(),
        c_obs.space_group().order_p(),
        c_obs.indices().size(),
        c_obs.d_min())
  #
  def use_smtbx_ls(mode):
    if ("ls_"+mode not in params.optimizers):
      return None
    if (not libtbx.env.has_module(name="smtbx")):
      print "INFO: smtbx not available: refinement skipped."
      return None
    result = structure_work.deep_copy_scatterers()
    run_smtbx_ls(
      mode=mode,
      cod_id=cod_id,
      i_obs=i_obs,
      f_obs=f_obs,
      xray_structure=result,
      params=params)
    show_cc_r1(params, "ls_"+mode, f_obs, result)
    return result
  structure_ls_simple = use_smtbx_ls("simple")
  structure_ls_lm = use_smtbx_ls("lm")
  #
  def use_shelxl(mode):
    if ("shelxl_"+mode not in params.optimizers):
      return None
    result = structure_work.deep_copy_scatterers()
    run_shelxl(
      mode=mode,
      cod_id=cod_id,
      i_obs=i_obs,
      f_obs=f_obs,
      xray_structure=result,
      params=params,
      reference_structure=structure_iso,
      expected_n_refinable_parameters=n_refinable_parameters)
    if (params.export_refined):
      file_name = "shelxl_%s_%s.pdb" % (mode, cod_id)
      open(file_name, "w").write(result.as_pdb_file(
        remarks=[file_name]))
    return result
  structure_shelxl_fm = use_shelxl("fm")
  structure_shelxl_cg = use_shelxl("cg")
  #
  if ("shelx76" not in params.optimizers):
    structure_shelx76 = None
  else:
    structure_shelx76 = structure_work.deep_copy_scatterers()
    run_shelx76(
      cod_id=cod_id,
      f_obs=f_obs,
      xray_structure=structure_shelx76,
      fvars=fvars,
      encoded_sites=encoded_sites,
      params=params,
      reference_structure=structure_iso)
    if (params.export_refined):
      file_name = "shelx76_%s.pdb" % cod_id
      open(file_name, "w").write(structure_shelx76.as_pdb_file(
        remarks=[file_name]))

def run(args):
  from iotbx.option_parser import option_parser as iotbx_option_parser
  import libtbx.utils
  show_times = libtbx.utils.show_times(time_start="now")
  command_call = ["iotbx.python", __file__]
  command_line = (iotbx_option_parser(
    usage=" ".join(command_call) + " [options] directory|file...")
    .enable_chunk(easy_all=True)
    .enable_multiprocessing()
  ).process(args=args, min_nargs=1)
  if (command_line.run_multiprocessing_chunks_if_applicable(
        command_call=command_call)):
    show_times()
    return
  co = command_line.options
  #
  print "TIME BEGIN cod_refine:", date_and_time()
  print
  #
  master_phil = get_master_phil()
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil)
  phil_objects = []
  remaining_args = []
  for arg in command_line.args:
    if (arg.find("=") >= 0):
      phil_objects.append(argument_interpreter.process(arg=arg))
    else:
      remaining_args.append(arg)
  work_phil = master_phil.fetch(sources=phil_objects)
  work_phil.show()
  print
  params = work_phil.extract()
  #
  qi_dict = {}
  all_pickles = []
  for arg in remaining_args:
    if (op.isdir(arg)):
      for node in sorted(os.listdir(arg)):
        if (node.endswith(".pickle")):
          all_pickles.append(op.join(arg, node))
        elif (node.startswith("qi_") and len(node) == 10):
          qi = open(op.join(arg, node)).read().splitlines()
          if (len(qi) == 1):
            cod_id = node[3:]
            quick_info = eval(qi[0])
            assert cod_id not in qi_dict
            qi_dict[cod_id] = quick_info
    elif (op.isfile(arg)):
      all_pickles.append(arg)
    else:
      raise RuntimeError("Not a file or directory: %s" % arg)
  print "Number of pickle files:", len(all_pickles)
  print "Number of quick_infos:", len(qi_dict)
  sort_choice = params.sorting_of_pickle_files
  if (len(qi_dict) != 0 and sort_choice is not None):
    print "Sorting pickle files by n_atoms * n_refl:", sort_choice
    assert sort_choice in ["down", "up"]
    def sort_pickle_files():
      if (sort_choice == "down"): i_sign = -1
      else:                       i_sign = 1
      buffer = []
      for i,path in enumerate(all_pickles):
        cod_id = op.basename(path).split(".",1)[0]
        qi = qi_dict.get(cod_id)
        if (qi is None): nn = 2**31
        else:            nn = qi[0] * qi[1] * qi[2]
        buffer.append((nn, i_sign*i, path))
      buffer.sort()
      if (i_sign < 0):
        buffer.reverse()
      result = []
      for elem in buffer:
        result.append(elem[-1])
      return result
    all_pickles = sort_pickle_files()
  print
  #
  rss = params.random_subset.size
  if (rss is not None and rss > 0):
    seed = params.random_subset.seed
    print "Selecting subset of %d pickle files using random seed %d" % (
      rss, seed)
    mt = flex.mersenne_twister(seed=seed)
    perm = mt.random_permutation(size=len(all_pickles))[:rss]
    flags = flex.bool(len(all_pickles), False).set_selected(perm, True)
    all_pickles = flex.select(all_pickles, permutation=flags.iselection())
    print
  #
  from libtbx.path import makedirs_race
  if (params.wdir_root is not None):
    makedirs_race(path=params.wdir_root)
  if (params.pickle_refined_dir is not None):
    makedirs_race(path=params.pickle_refined_dir)
  #
  n_caught = 0
  for i_pickle,pickle_file_name in enumerate(all_pickles):
    if (i_pickle % command_line.chunk.n != command_line.chunk.i): continue
    tm = user_plus_sys_time()
    try:
      process(params, pickle_file_name)
    except KeyboardInterrupt:
      print >> sys.stderr, "CAUGHT EXCEPTION: KeyboardInterrupt"
      traceback.print_exc()
      print >> sys.stderr
      sys.stderr.flush()
      return
    except Exception:
      sys.stdout.flush()
      print >> sys.stderr, "CAUGHT EXCEPTION: %s" % pickle_file_name
      traceback.print_exc()
      print >> sys.stderr
      sys.stderr.flush()
      n_caught += 1
    else:
      print "done_with: %s (%.2f seconds)" % (pickle_file_name, tm.elapsed())
      print
      sys.stdout.flush()
  print
  print "Number of exceptions caught:", n_caught
  #
  show_times()
  print
  print "TIME END cod_refine:", date_and_time()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
