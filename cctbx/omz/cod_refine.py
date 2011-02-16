from cctbx import omz
import cctbx.omz.dev
import libtbx.phil.command_line
from libtbx import easy_pickle
from libtbx.utils import user_plus_sys_time
from libtbx import Auto
import traceback
import sys, os
op = os.path

def get_master_phil():
  return omz.dev.get_master_phil(
    iteration_limit=100,
    grads_mean_sq_threshold=1e-6,
    additional_phil_string="""\
      max_atoms = 100
        .type = int
      reset_u_iso = 0.05
        .type = float
      optimizers = *dev smtbx_lsq
        .type = choice(multi=True)
      smtbx_lsq_iterations = 12
        .type = int
""")

def show_cc_r1(label, f_obs, xray_structure):
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure=xray_structure).f_calc().amplitudes()
  from cctbx.array_family import flex
  corr = flex.linear_correlation(x=f_obs.data(), y=f_calc.data())
  assert corr.is_well_defined()
  cc = corr.coefficient()
  r1 = f_obs.r1_factor(
    other=f_calc, scale_factor=Auto, assume_index_matching=True)
  print "%-12s cc, r1: %.3f %.3f" % (label, cc, r1)
  sys.stdout.flush()

def smtbx_lsq_refine(cod_code, f_obs, xray_structure, n_cycles):
  import smtbx.refinement
  xray_structure.scatterers().flags_set_grads(state=False)
  for sc in xray_structure.scatterers():
    sc.flags.set_grad_site(True)
    assert sc.flags.use_u_iso_only()
    sc.flags.set_grad_u_iso(True)
  fo_sq = f_obs.f_as_f_sq()
  assert fo_sq.sigmas() is not None
  sel = (fo_sq.data() == 0) & (fo_sq.sigmas() == 0)
  fo_sq = fo_sq.select(~sel)
  fo_sq.select(fo_sq.sigmas() <= 0).show_array()
  assert fo_sq.sigmas().all_gt(0)
  if (1): # work around bug currently in smtbx weighting scheme implementation
    from cctbx.array_family import flex
    fo_sq = fo_sq.customized_copy(sigmas=flex.double(fo_sq.data().size(), 1))
  tm = user_plus_sys_time()
  rm = smtbx.refinement.model(
    fo_sq=fo_sq,
    xray_structure=xray_structure,
    constraints=[],
    restraints_manager=smtbx.refinement.restraints.manager(),
    weighting_scheme=smtbx.refinement.least_squares.unit_weighting())
  ls = rm.least_squares()
  for i_cycle in xrange(n_cycles):
    ls.build_up()
    try:
      ls.solve_and_step_forward()
    except RuntimeError, e:
      if (str(e).find("cholesky.failure") <= 0): raise
      print "Aborting smtbx_lsq_refine: cholesky.failure: %s" % cod_code
      break
    for sc in xray_structure.scatterers():
      if (sc.u_iso <= 0 or sc.u_iso > 1):
        sc.u_iso = 0.05
    show_cc_r1("ls%02d" % (i_cycle+1), f_obs, xray_structure)
  tm.show_elapsed(prefix="time smtbx lsq: ")

def process(params, pickle_file_name):
  cod_code = op.basename(pickle_file_name).split(".",1)[0]
  print "cod_code:", cod_code
  f_obs, structure_cod = easy_pickle.load(file_name=pickle_file_name)
  structure_cod.show_summary().show_scatterers()
  print "."*79
  f_obs.show_comprehensive_summary()
  print "."*79
  #
  structure_work = structure_cod.deep_copy_scatterers()
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
    from scitbx.array_family import flex
    mt = flex.mersenne_twister(seed=0)
    structure_work.shake_sites_in_place(
      rms_difference=params.shake_sites_rmsd,
      allow_all_fixed=True,
      random_double=mt.random_double)
    cc_r1("shake_xyz")
  #
  if (params.max_atoms is not None):
    n = structure_work.scatterers().size()
    if (n > params.max_atoms):
      print "Skipping refinement of large model: %d atoms COD %s" % (
        n, cod_code)
      return
  #
  if ("smtbx_lsq" in params.optimizers):
    structure_smtbx_lsq = structure_work.deep_copy_scatterers()
    smtbx_lsq_refine(
      cod_code=cod_code,
      f_obs=f_obs,
      xray_structure=structure_smtbx_lsq,
      n_cycles=params.smtbx_lsq_iterations)
  #
  if ("dev" in params.optimizers):
    structure_dev = structure_work.deep_copy_scatterers()
    omz.dev.ls_refinement(
      f_obs=f_obs,
      xray_structure=structure_dev,
      params=params,
      reference_structure=structure_iso)
    show_cc_r1("dev", f_obs, structure_dev)

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
