from __future__ import division
from cctbx import omz
import cctbx.omz.dev
from cctbx.array_family import flex
import libtbx.phil.command_line
from libtbx.test_utils import approx_equal
from libtbx import easy_run
from libtbx import easy_pickle
from libtbx.utils import user_plus_sys_time
from libtbx import Auto
from cStringIO import StringIO
import traceback
import sys, os
op = os.path

def get_master_phil():
  return omz.dev.get_master_phil(
    iteration_limit=100,
    grads_mean_sq_threshold=1e-6,
    additional_phil_string="""\
      max_atoms = 99
        .type = int
      f_obs_f_calc_fan_outliers = *remove keep
        .type = choice
        .optional = False
      reset_u_iso = 0.05
        .type = float
      sites_mod_positive = True
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
      keep_tmp_files = False
        .type = bool
      export_refined = False
        .type = bool
      try_special_op_simplifier = False
        .type = bool
""")

def shelxl_weights(fo_sq, sigmas, fc_sq, scale_factor, a=0.1, b=0):
  assert sigmas.size() == fo_sq.size()
  assert fc_sq.size() == fo_sq.size()
  result = flex.double()
  for o,s,c in zip(fo_sq, sigmas, fc_sq):
    o /= scale_factor
    s /= scale_factor
    p = (max(o, 0) + 2 * c) / 3
    den = s**2 + (a*p)**2 + b*p
    assert den > 1e-8
    w = 1 / den
    result.append(w)
  return result

def show_cc_r1(
      label, f_obs, xray_structure=None, fc_abs=None, scale_factor=Auto):
  assert [xray_structure, fc_abs].count(None) == 1
  if (fc_abs is None):
    fc_abs = f_obs.structure_factors_from_scatterers(
      xray_structure=xray_structure,
      algorithm="direct",
      cos_sin_table=False).f_calc().amplitudes()
  corr = flex.linear_correlation(x=f_obs.data(), y=fc_abs.data())
  assert corr.is_well_defined()
  cc = corr.coefficient()
  r1 = f_obs.r1_factor(
    other=fc_abs, scale_factor=scale_factor, assume_index_matching=True)
  print "%-12s cc, r1: %.3f %.3f" % (label, cc, r1)
  sys.stdout.flush()
  return fc_abs, cc, r1

def run_smtbx_ls(mode, cod_code, f_obs, xray_structure, params):
  import smtbx.refinement
  fo_sq = f_obs.f_as_f_sq(algorithm="shelxl")
  assert fo_sq.sigmas() is not None
  sel = (fo_sq.data() == 0) & (fo_sq.sigmas() == 0)
  fo_sq = fo_sq.select(~sel)
  fo_sq.select(fo_sq.sigmas() <= 0).show_array()
  assert fo_sq.sigmas().all_gt(0)
  if (1): # work around bug currently in smtbx weighting scheme implementation
    fo_sq = fo_sq.customized_copy(sigmas=flex.double(fo_sq.data().size(), 1))
  tm = user_plus_sys_time()
  rm = smtbx.refinement.model(
    fo_sq=fo_sq,
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
          % cod_code
        break
      for sc in xray_structure.scatterers():
        if (sc.u_iso <= 0 or sc.u_iso > 1):
          sc.u_iso = 0.05
      show_cc_r1("ls%02d" % (i_cycle+1), f_obs, xray_structure)
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
        ' debye_waller_factor_exp failure: %s' % cod_code
    show_cc_r1("smtbx_lm", f_obs, xray_structure)
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
      cod_code,
      f_obs,
      xray_structure,
      params,
      reference_structure,
      expected_n_refinable_parameters):
  if (mode == "fm"):
    fm_cycles = params.shelxl_fm_iterations
    cg_cycles = None
  elif (mode == "cg"):
    fm_cycles = None
    cg_cycles = params.shelxl_cg_iterations
  else:
    raise RuntimeError("Unknown mode: " + mode)
  cwd_orig = os.getcwd()
  wdir = "wdir_%s" % cod_code
  wdir_is_new = False
  if (not op.isdir(wdir)):
    os.mkdir(wdir)
    wdir_is_new = True
  remove_wdir = False
  try:
    os.chdir(wdir)
    tmp_file_names = ["tmp.ins", "tmp.hkl", "tmp.res", "tmp.lst"]
    remove_tmp_files(tmp_file_names)
    if (params.shelxl_reset_sigmas):
      f_obs = f_obs.customized_copy(
        sigmas=flex.double(f_obs.indices().size(), params.shelxl_reset_sigmas))
    import iotbx.shelx
    open("tmp.ins", "w").writelines(iotbx.shelx.writer.generator(
      xray_structure=xray_structure,
      data_are_intensities=False,
      title="cod_code=%s mode=%s" % (cod_code, mode),
      wavelength=f_obs.minimum_wavelength_based_on_d_min(),
      full_matrix_least_squares_cycles=fm_cycles,
      conjugate_gradient_least_squares_cycles=cg_cycles,
      weighting_scheme_params=params.shelxl_wght,
      sort_scatterers=False))
    f_obs.export_as_shelx_hklf(file_object=open("tmp.hkl", "w"))
    buffers = easy_run.fully_buffered("shelxl tmp")
    buffers.raise_if_errors()
    refinement_unstable = False
    for line in buffers.stdout_lines:
      if (line.find("** REFINEMENT UNSTABLE **") >= 0):
        refinement_unstable = True
        print "Aborted: shelxl %s refinement unstable: %s" % (mode, cod_code)
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
        mode, cod_code)
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
      res_fvar = None
      res_hkl_count = None
      res_r1 = None
      res_n_parameters = None
      res_n_restraints = None
      for line in res.splitlines():
        if (line.startswith("FVAR ")):
          flds = line.split()
          assert len(flds) == 2
          res_fvar = float(flds[1])
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
        assert res_fvar is not None
        assert res_hkl_count is not None
        assert res_r1 is not None
        assert res_n_parameters is not None
        assert res_n_restraints is not None
        #
        assert res_hkl_count == f_obs.indices().size()
        def raise_unexpected_restraints(n_expected):
          raise RuntimeError(
            "Unexpected number of SHELXL restraints: %d (vs. %d expected)" % (
              res_n_restraints, n_expected))
        if (mode == "fm"):
          n_caos = f_obs.space_group_info() \
            .number_of_continuous_allowed_origin_shifts()
          if (res_n_restraints != n_caos):
            sg_symbol = str(f_obs.space_group_info())
            if (sg_symbol in ["P 63 m c", "P 63 c m"]):
              assert n_caos == 1
              assert res_n_restraints == 0
              print "INFO: SHELXL restraint count incorrect? code_code:", \
                cod_code
            else:
              raise_unexpected_restraints(n_caos)
        elif (mode == "cg"):
          if (res_n_restraints != 0):
            raise_unexpected_restraints(0)
        else:
          raise RuntimeError("Unknown mode: " + mode)
        assert res_n_parameters == expected_n_refinable_parameters + 1
        fc_abs, _, r1_fvar = show_cc_r1(
          "fvar_"+mode, f_obs, xray_structure, scale_factor=res_fvar)
        r1_diff = r1_fvar - res_r1
        print "R1 recomputed - shelxl_%s.res: %.4f - %.4f = %.4f %s" % (
          mode, r1_fvar, res_r1, r1_diff, cod_code)
        if (abs(r1_diff) > 0.01):
          raise RuntimeError("R1 MISMATCH %s" % cod_code)
        _, _, r1_auto = show_cc_r1("shelxl_"+mode, f_obs, fc_abs=fc_abs)
        print "R1 FVAR-Auto %s: %.4f" % (cod_code, r1_fvar - r1_auto)
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
        import iotbx.shelx.hklf
        fo_wr = iotbx.shelx.hklf.reader(file_name="tmp.hkl") \
          .as_miller_arrays(crystal_symmetry=f_obs)[0]
        fo_sq = fo_wr.f_as_f_sq(algorithm="shelxl")
        fc_sq = fc_abs.f_as_f_sq()
        weights = shelxl_weights(
          fo_sq=fo_sq.data(),
          sigmas=fo_sq.sigmas(),
          fc_sq=fc_sq.data(),
          scale_factor=res_fvar)
        num = flex.sum(
          weights * flex.pow2(fo_sq.data() / res_fvar - fc_sq.data()))
        den = flex.sum(
          weights * flex.pow2(fo_sq.data() / res_fvar))
        assert den != 0
        wr2 = (num / den)**0.5
        wr2_diff = wr2 - lst_wr2
        if (abs(wr2_diff) > 0.01):
          info = " significantly different"
        else:
          info = ""
        print "wR2 recomputed - shelxl_%s.lst: %.4f - %.4f = %.4f %s%s" % (
          mode, lst_wr2, wr2, wr2_diff, cod_code, info)
        if (abs(wr2_diff) / max(lst_wr2, wr2) > 0.2):
          raise RuntimeError("wR2 MISMATCH %s" % cod_code)
    if (not params.keep_tmp_files):
      remove_tmp_files(tmp_file_names)
      remove_wdir = wdir_is_new
  finally:
    os.chdir(cwd_orig)
    if (remove_wdir):
      os.rmdir(wdir)

def run_shelx76(
      cod_code,
      f_obs,
      xray_structure,
      params,
      reference_structure,
      expected_n_refinable_parameters):
  cwd_orig = os.getcwd()
  wdir = "wdir_%s" % cod_code
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
      xray_structure=xray_structure,
      f_obs=f_obs,
      l_s_parameters=str(params.shelx76_iterations))
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
            cod_code)
          r_from_lst = "nan"
          break
        r_from_lst = float(flds[2])
    assert r_from_lst is not None
    if (r_from_lst != "nan"):
      print "%-12s cc, r1: None %.3f" % ("shelx76", r_from_lst)
      if (not params.keep_tmp_files):
        remove_tmp_files(tmp_file_names)
        remove_wdir = wdir_is_new
  finally:
    os.chdir(cwd_orig)
    if (remove_wdir):
      os.rmdir(wdir)

class row_expr(object):
  __slots__ = ["indices", "multipliers", "constant"]
  def __init__(O, indices, multipliers, constant):
    O.indices = indices
    O.multipliers = multipliers
    O.constant = constant
  def __str__(O):
    s = ""
    for i,m in zip(O.indices, O.multipliers):
      assert m != 0
      if (m < 0):
        s += "-"
        m = -m
      elif (len(s) != 0):
        s += "+"
      if (m != 1):
        s += str(m) + "*"
      s += "xyz"[i]
    c = O.constant
    if (c != 0 or len(s) == 0):
      if (c < 0):
        s += "-"
        c = -c
      elif (len(s) != 0):
        s += "+"
      s += str(c)
    return s

def vector_multiplier(a, b):
  result = None
  for i,j in zip(a,b):
    if (i == 0):
      if (j != 0): return None
    else:
      if (j == 0): return None
      m = j / i
      if (result is None):
        result = m
      elif (result != m):
        return None
  return result

def special_op_simplifier(special_op):
  rt = special_op.as_rational()
  r = rt.r
  t = rt.t
  rows = [r[:3], r[3:6], r[6:]]
  result = [None, None, None]
  import boost.rational
  r0 = boost.rational.int(0)
  r1 = boost.rational.int(1)
  n_done = 0
  for i_row,row in enumerate(rows):
    if (row == (0,0,0)):
      result[i_row] = row_expr([], [], t[i_row])
      n_done += 1
  if (n_done == 3):
    return result
  if (n_done == 0):
    m = []
    v = []
    for i in xrange(3):
      m.append([r[i+0], r[i+3]])
      v.append(r[i+6])
    from scitbx.matrix import row_echelon
    free_vars = row_echelon.form_rational(m, v)
    if (len(free_vars) == 0):
      sol = row_echelon.back_substitution_rational(m, v, free_vars, [None]*2)
      if (sol is not None and sol.count(0) == 0):
        for i_row in [0,1]:
          result[i_row] = row_expr([i_row], [r1], r0)
        result[2] = row_expr([0,1], sol, t[2] - sol[0]*t[0] - sol[1]*t[1])
        return result
  for i_row in xrange(3):
    if (result[i_row] is not None): continue
    result[i_row] = row_expr([i_row], [r1], r0)
    for j_row in xrange(i_row+1,3):
      if (result[j_row] is not None): continue
      m = vector_multiplier(rows[i_row], rows[j_row])
      if (m is None): continue
      result[j_row] = row_expr([i_row], [m], t[j_row] - m*t[i_row])
  return result

def sx76ss(xs):
  scs = xs.scatterers()
  sstab = xs.site_symmetry_table()
  for i_sc in xs.special_position_indices():
    ss = sstab.get(i_sc)
    con = ss.site_constraints()
    sc = scs[i_sc]
    site_indep = con.independent_params(all_params=sc.site)
    print "%-10s" % sc.label, numstr(sc.site)
    print ss.special_op()
    sos = special_op_simplifier(ss.special_op())
    expr = ",".join([str(e) for e in sos])
    print "LOOK", expr
    ns = dict(zip("xyz", sc.site))
    expr_site = eval(expr, ns, {})
    print "expr_site:", numstr(expr_site)
    assert approx_equal(expr_site, sc.site, 1e-4)
    print

def process(params, pickle_file_name):
  cod_code = op.basename(pickle_file_name).split(".",1)[0]
  print "cod_code:", cod_code
  f_obs, structure_cod = easy_pickle.load(file_name=pickle_file_name)
  changes = structure_cod.make_scatterer_labels_shelx_compatible_in_place()
  structure_cod.show_summary().show_scatterers()
  if (params.try_special_op_simplifier):
    sx76ss(structure_cod)
    sx76ss(structure_cod.niggli_cell())
    return
  if (len(changes) != 0):
    from libtbx.utils import plural_s
    print "INFO: %d atom name%s changed for compatibility with SHELXL:" \
      % plural_s(len(changes))
    for change in changes:
      print '  changed: "%s" -> "%s"' % change
  structure_cod.scattering_type_registry(table="it1992").show()
  print "."*79
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure=structure_cod,
    algorithm="direct",
    cos_sin_table=False).f_calc()
  sel = f_obs.f_obs_f_calc_fan_outlier_selection(f_calc=f_calc)
  assert sel is not None
  n_outliers = sel.count(True)
  if (n_outliers != 0):
    action = params.f_obs_f_calc_fan_outliers
    print "INFO: f_obs_f_calc_fan_outliers = %s: %d" % (action, n_outliers)
    if (action == "remove"):
      f_obs = f_obs.select(~sel)
  if (f_obs.anomalous_flag()):
    print "INFO: converting anomalous f_obs to non-anomalous."
    f_obs = f_obs.average_bijvoet_mates()
  sel = (f_obs.data() == 0) & (f_obs.sigmas() == 0)
  n_zero_f_and_s = sel.count(True)
  if (n_zero_f_and_s != 0):
    print "INFO: removing reflections with f_obs=0 and sigma=0:", \
      n_zero_f_and_s
    f_obs = f_obs.select(~sel)
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure=structure_cod,
    algorithm="direct",
    cos_sin_table=False).f_calc()
  k = f_obs.scale_factor(f_calc=f_calc)
  assert k != 0
  s = 1/k
  print "INFO: scaling f_obs to f_calc by multiplying f_obs with: %.6g" % s
    # scaling applied so that the data written in shelx hklf format
    # have sufficient significant digits, and FVAR is 1 (shelx76 seems
    # to be especially sensitive to FVAR >> 1)
  f_obs = f_obs.customized_copy(data=f_obs.data()*s, sigmas=f_obs.sigmas()*s)
  f_obs.show_comprehensive_summary()
  print "."*79
  #
  if (params.sites_mod_positive):
    structure_work = structure_cod.sites_mod_positive()
  else:
    structure_work = structure_cod.deep_copy_scatterers()
  structure_work.scattering_type_registry(table="it1992")
  def cc_r1(label):
    show_cc_r1(label, f_obs, structure_work)
  #
  cc_r1("cod")
  #
  sel = structure_work.hd_selection()
  print "Removing hydrogen atoms:", sel.count(True)
  structure_work = structure_work.select(selection=~sel)
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
        n, cod_code)
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
  if ("dev" not in params.optimizers):
    structure_dev = None
  else:
    structure_dev = structure_work.deep_copy_scatterers()
    omz.dev.refinement(
      f_obs=f_obs,
      xray_structure=structure_dev,
      params=params,
      reference_structure=structure_iso,
      expected_n_refinable_parameters=n_refinable_parameters)
    show_cc_r1("dev", f_obs, structure_dev)
    if (params.export_refined):
      file_name = "dev_%s_%s_%s.pdb" % (
        params.target_type, params.target_obs_type.lower(), cod_code)
      open(file_name, "w").write(structure_dev.as_pdb_file(
        remarks=[file_name]))
  #
  def use_smtbx_ls(mode):
    if ("ls_"+mode not in params.optimizers):
      return None
    result = structure_work.deep_copy_scatterers()
    run_smtbx_ls(
      mode=mode,
      cod_code=cod_code,
      f_obs=f_obs,
      xray_structure=result,
      params=params)
    show_cc_r1("ls_"+mode, f_obs, result)
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
      cod_code=cod_code,
      f_obs=f_obs,
      xray_structure=result,
      params=params,
      reference_structure=structure_iso,
      expected_n_refinable_parameters=n_refinable_parameters)
    if (params.export_refined):
      file_name = "shelxl_%s_%s.pdb" % (mode, cod_code)
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
      cod_code=cod_code,
      f_obs=f_obs,
      xray_structure=structure_shelx76,
      params=params,
      reference_structure=structure_iso,
      expected_n_refinable_parameters=n_refinable_parameters)
    if (params.export_refined):
      file_name = "shelx76_%s.pdb" % cod_code
      open(file_name, "w").write(structure_shelx76.as_pdb_file(
        remarks=[file_name]))

def run(args):
  from iotbx.option_parser import option_parser as iotbx_option_parser
  from libtbx import easy_pickle
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
  all_pickles = []
  for arg in remaining_args:
    if (op.isdir(arg)):
      for node in sorted(os.listdir(arg)):
        if (not node.endswith(".pickle")): continue
        all_pickles.append(op.join(arg, node))
    elif (op.isfile(arg)):
      all_pickles.append(arg)
    else:
      raise RuntimeError("Not a file or directory: %s" % arg)
  print "Number of pickle files:", len(all_pickles)
  print
  #
  n_caught = 0
  for i_pickle,pickle_file_name in enumerate(all_pickles):
    if (i_pickle % command_line.chunk.n != command_line.chunk.i): continue
    try:
      process(params, pickle_file_name)
    except KeyboardInterrupt:
      print "CAUGHT EXCEPTION: KeyboardInterrupt"
      return
    except Exception:
      sys.stdout.flush()
      print >> sys.stderr, "CAUGHT EXCEPTION: %s" % pickle_file_name
      traceback.print_exc()
      print >> sys.stderr
      sys.stderr.flush()
      n_caught += 1
    else:
      print "done_with: %s" % pickle_file_name
      print
      sys.stdout.flush()
  print
  print "Number of exceptions caught:", n_caught
  #
  show_times()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
