import traceback
import sys, os
op = os.path

class cod_data(object):

  def __init__(O, cod_code, hkl_cif_pair):
    O.cod_code = cod_code
    refl_file, model_file = hkl_cif_pair
    import iotbx.cif.builders
    print "refl_file:", refl_file
    print "model_file:", model_file
    refl_cif = iotbx.cif.reader(file_path=refl_file)
    model_cif = iotbx.cif.reader(file_path=model_file)
    from_coordinate_files = []
    from_reflection_files = []
    def get_cs(cif, buffer):
      for cif_block in cif.model().values():
        cs = iotbx.cif.builders.crystal_symmetry_builder(
          cif_block=cif_block).crystal_symmetry
        buffer.append(cs)
    get_cs(refl_cif, from_reflection_files)
    get_cs(model_cif, from_coordinate_files)
    import cctbx.crystal
    combined_cs = cctbx.crystal.select_crystal_symmetry(
      from_coordinate_files=from_coordinate_files,
      from_reflection_files=from_reflection_files)
    assert combined_cs.unit_cell() is not None
    assert combined_cs.space_group_info() is not None
    #
    miller_arrays = refl_cif.as_miller_arrays()
    meas_a = []
    meas_i = []
    for ma in miller_arrays:
      s = str(ma.info())
      if (s.find("_meas") >= 0):
        if (ma.is_xray_amplitude_array()):
          meas_a.append(ma)
        elif (ma.is_xray_intensity_array()):
          meas_i.append(ma)
    if (len(meas_a) != 0):
      O.f_obs = meas_a[0]
    elif (len(meas_i) != 0):
      O.f_obs = meas_i[0].f_sq_as_f()
    else:
      raise RuntimeError("Missing f_obs array.")
    O.f_obs = O.f_obs.customized_copy(
      crystal_symmetry=combined_cs)
    print "."*79
    O.f_obs.show_comprehensive_summary()
    print "."*79
    #
    O.xray_structure = cctbx.xray.structure.from_cif(file_path=model_file)
    O.xray_structure = O.xray_structure.customized_copy(
      crystal_symmetry=combined_cs)
    O.xray_structure.show_summary().show_scatterers()
    print "."*79

  def have_zero_occupancies(O):
    return not O.xray_structure.scatterers().extract_occupancies().all_ne(0)

  def have_shelxl_compatible_scattering_types(O):
    from cctbx.eltbx.xray_scattering import \
      shelxl_97_2_980324_tabulated_chemical_elements as known
    for sc in O.xray_structure.scatterers():
      if (sc.scattering_type.strip().upper() not in known):
        return False
    return True

  def have_sys_absent(O):
    return (O.f_obs.sys_absent_flags().data().count(True) != 0)

  def have_redundant_data(O):
    return not O.f_obs.is_unique_set_under_symmetry()

  def have_bad_sigmas(O):
    if (O.f_obs.sigmas() is None):
      print "Missing sigmas:", O.cod_code
      return True
    f_sq = O.f_obs.f_as_f_sq()
    sel = (f_sq.data() == 0) & (f_sq.sigmas() == 0)
    result = not f_sq.select(~sel).sigmas().all_gt(0)
    if (result):
      print "Zero or negative sigmas:", O.cod_code
    return result

  def f_obs_f_calc_correlation(O):
    f_calc = O.f_obs.structure_factors_from_scatterers(
      xray_structure=O.xray_structure).f_calc()
    from cctbx.array_family import flex
    lc = flex.linear_correlation(
      x=O.f_obs.data(),
      y=f_calc.amplitudes().data())
    assert lc.is_well_defined()
    result = lc.coefficient()
    print "f_obs_f_calc_correlation: %.3f %s" % (
      lc.coefficient(), O.cod_code)
    return result

def build_hkl_cif(cod_codes):
  cod_svn = os.environ.get("COD_SVN_WORKING_COPY")
  assert cod_svn is not None
  cif_dir = op.join(cod_svn, "cif")
  hkl_dir = op.join(cod_svn, "hkl")
  hkl_cif = []
  if (len(cod_codes) == 0):
    hkl_only = []
    for sub_dir in sorted(os.listdir(hkl_dir)):
      if (sub_dir.startswith(".")): continue
      hkl_sub_dir = op.join(hkl_dir, sub_dir)
      for node in sorted(os.listdir(hkl_sub_dir)):
        if (node.startswith(".")): continue
        if (not node.endswith(".hkl")): continue
        cod_code = node[:-4]
        hkl_path = op.join(hkl_sub_dir, node)
        cif_path = op.join(cif_dir, sub_dir, cod_code+".cif")
        if (not op.isfile(cif_path)):
          hkl_only.append(hkl_path)
        else:
          hkl_cif.append((hkl_path, cif_path))
    print "Number of hkl without cif:", len(hkl_only)
  else:
    n_missing_all = 0
    for cod_code in cod_codes:
      hkl_path = op.join(hkl_dir, cod_code[0], cod_code+".hkl")
      cif_path = op.join(cif_dir, cod_code[0], cod_code+".cif")
      n_missing = 0
      if (not op.isfile(cif_path)):
        print "Missing COD cif file:", cif_path
        n_missing += 1
      if (not op.isfile(hkl_path)):
        print "Missing COD hkl file:", hkl_path
        n_missing += 1
      if (n_missing == 0):
        hkl_cif.append((hkl_path, cif_path))
      else:
        n_missing_all += n_missing
    if (n_missing_all != 0):
      raise RuntimeError("Number of missing COD files: %d" % n_missing_all)
  print "Number of hkl+cif:", len(hkl_cif)
  return hkl_cif

def run(args):
  from iotbx.option_parser import option_parser as iotbx_option_parser
  from libtbx import easy_pickle
  import libtbx.utils
  show_times = libtbx.utils.show_times(time_start="now")
  command_call = ["iotbx.python", __file__]
  command_line = (iotbx_option_parser(
    usage=" ".join(command_call) + " [options] [cod_code...]")
    .enable_chunk(easy_all=True)
    .enable_multiprocessing()
    .option(None, "--min_f_obs_f_calc_correlation",
      type="float",
      default=0.99,
      metavar="FLOAT")
  ).process(args=args)
  if (command_line.run_multiprocessing_chunks_if_applicable(
        command_call=command_call)):
    show_times()
    return
  co = command_line.options
  #
  hkl_cif = build_hkl_cif(cod_codes=command_line.args)
  #
  if (not op.isdir("cod_ma_xs")):
    os.makedirs("cod_ma_xs")
  n_caught = 0
  for i_pair,pair in enumerate(hkl_cif):
    cod_code = op.basename(pair[0])[:-4]
    if (i_pair % command_line.chunk.n != command_line.chunk.i): continue
    try:
      cd = cod_data(cod_code=cod_code, hkl_cif_pair=pair)
    except KeyboardInterrupt:
      print "CAUGHT EXCEPTION: KeyboardInterrupt"
      return
    except Exception:
      sys.stdout.flush()
      print >> sys.stderr, "CAUGHT EXCEPTION: cod.py: %s" % cod_code
      traceback.print_exc()
      print >> sys.stderr
      sys.stderr.flush()
      n_caught += 1
    else:
      if (    not cd.have_zero_occupancies()
          and cd.have_shelxl_compatible_scattering_types()
          and not cd.have_sys_absent()
          and not cd.have_redundant_data()
          and not cd.have_bad_sigmas()
          and cd.f_obs_f_calc_correlation() >= co.min_f_obs_f_calc_correlation):
        easy_pickle.dump(
          file_name="cod_ma_xs/%s.pickle" % cod_code,
          obj=(cd.f_obs, cd.xray_structure))
      else:
        print "filtering out:", cod_code
      print "done_with:", cod_code
      print
  print
  print "Number of exceptions caught:", n_caught
  #
  show_times()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
