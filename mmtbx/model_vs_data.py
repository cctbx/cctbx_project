import sys, os, time, random
from cctbx.array_family import flex
from iotbx import pdb
from cctbx import adptbx
from iotbx.option_parser import iotbx_option_parser
from libtbx.utils import Sorry
from iotbx import reflection_file_utils
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.f_model
from iotbx.pdb import crystal_symmetry_from_pdb
from iotbx import crystal_symmetry_from_any
from iotbx import reflection_file_utils
from iotbx import reflection_file_reader
from libtbx.str_utils import format_value
import iotbx
from mmtbx import utils
from iotbx import pdb
from libtbx import easy_run
from cStringIO import StringIO
from mmtbx import model_statistics
from iotbx.pdb import extract_rfactors_resolutions_sigma
import iotbx.pdb.remark_3_interpretation
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from iotbx.pdb import combine_unique_pdb_files
from libtbx import group_args
from mmtbx import masks
from mmtbx import maps
from mmtbx import real_space_correlation

if (1):
  random.seed(0)
  flex.set_random_seed(0)

class mvd(object):

  def __init__(self):
    self.crystal       = None
    self.models        = None
    self.data          = None
    self.model_vs_data = None
    self.pdb_header    = None
    self.misc          = None
    self.pdb_file      = None
    self.fmodel        = None

  def collect(self,
              crystal       = None,
              models        = None,
              data          = None,
              model_vs_data = None,
              pdb_header    = None,
              misc          = None):
    if(crystal       is not None): self.crystal       = crystal
    if(models        is not None): self.models        = models
    if(data          is not None): self.data          = data
    if(model_vs_data is not None): self.model_vs_data = model_vs_data
    if(pdb_header    is not None): self.pdb_header    = pdb_header
    if(misc          is not None): self.misc          = misc

  def show(self, log = None):
    if(log is None): log = sys.stdout
    # crystal
    print >> log, "  Unit cell:       ", self.crystal.uc
    print >> log, "  Space group:     ", self.crystal.sg, \
                  "number of symmetry operations:", self.crystal.n_sym_op
    # models
    print >> log, "  Number of models:", len(self.models)
    for i_seq, i_model in enumerate(self.models):
      print >> log, "  Model #%s:"%str("%d"%(i_seq+1)).strip()
      print >> log, "    Number of residues in alternative conformations:", \
        i_model.overall_counts_i_seq.n_alt_conf_pure + \
        i_model.overall_counts_i_seq.n_alt_conf_proper + \
        i_model.overall_counts_i_seq.n_alt_conf_improper
      rc = i_model.overall_counts_i_seq.resname_classes
      print >> log, "    Residue content:"
      len_max = 0
      for k in rc.keys():
        k=k.replace("common_","")
        if(len(k)>len_max): len_max = len(k)
      fmt = "      %-"+str(len_max)+"s : %d"
      for k,v in zip(rc.keys(), rc.values()):
        print >> log, fmt%(k.replace("common_",""), v)
      x = i_model.xray_structure_stat
      print >> log, "    Atoms:"
      print >> log, "      atom_number_(type:count:occ_sum) : %s (%s)"%(x.n_atoms,x.atom_counts_str)
      print >> log, "      ADP_(min,max,mean)               : %s %s %s"%(x.b_min,x.b_max,x.b_mean)
      print >> log, "      occupancies_(min,max,mean)       : %s %s %s"%(x.o_min,x.o_max,x.o_mean)
      print >> log, "      number_of_anisotropic            : "+format_value("%-7s",x.n_aniso)
      print >> log, "      number_of_non_positive_definite  : %s"%x.n_npd
      g = i_model.model_statistics_geometry
      if(g is not None):
        print >> log, "    Stereochemistry statistics (mean, max, count):"
        print >> log, "      bonds            : %8.4f %8.4f %d" % (g.b_mean, g.b_max, g.b_number)
        print >> log, "      angles           : %8.4f %8.4f %d" % (g.a_mean, g.a_max, g.a_number)
        print >> log, "      dihedrals        : %8.4f %8.4f %d" % (g.d_mean, g.d_max, g.d_number)
        print >> log, "      chirality        : %8.4f %8.4f %d" % (g.c_mean, g.c_max, g.c_number)
        print >> log, "      planarity        : %8.4f %8.4f %d" % (g.p_mean, g.p_max, g.p_number)
        print >> log, "      non-bonded (min) : %8.4f" % (g.n_min)
      if(i_model.ramalyze is not None):
        outl = i_model.ramalyze.get_outliers_count_and_fraction()
        allo = i_model.ramalyze.get_allowed_count_and_fraction()
        favo = i_model.ramalyze.get_favored_count_and_fraction()
        print >> log, "      Ramachandran plot, number of:"
        print >> log, "        outliers : %-5d (%-5.2f %s)"%(outl[0],outl[1]*100.,"%")
        print >> log, "        allowed  : %-5d (%-5.2f %s)"%(allo[0],allo[1]*100.,"%")
        print >> log, "        favored  : %-5d (%-5.2f %s)"%(favo[0],favo[1]*100.,"%")
      if(i_model.rotalyze is not None):
        print >> log, "      Rotamer outliers        : %d (%s %s) goal: %s" %(
          i_model.rotalyze.get_outliers_count_and_fraction()[0],
          str("%6.2f"%(i_model.rotalyze.get_outliers_count_and_fraction()[1]*100.)).strip(),
          "%", i_model.rotalyze.get_outliers_goal())
      if(i_model.cbetadev is not None):
        print >> log, "      Cbeta deviations >0.25A : %d"%i_model.cbetadev.get_outlier_count()
      if(i_model.clashscore is not None):
        print >> log, "      All-atom clashscore     : %.2f (steric overlaps >0.4A per 1000 atoms)"% \
          i_model.clashscore.get_clashscore()
    #
    print >> log, "  Data:"
    result = " \n    ".join([
      "data_label              : %s"%                    self.data.data_label,
      "high_resolution         : "+format_value("%-5.2f",self.data.high_resolution),
      "low_resolution          : "+format_value("%-6.2f",self.data.low_resolution),
      "completeness_in_range   : "+format_value("%-6.2f",self.data.completeness_in_range),
      "completeness(d_min-inf) : "+format_value("%-6.2f",self.data.completeness_d_min_inf),
      "completeness(6A-inf)    : "+format_value("%-6.2f",self.data.completeness_6A_inf),
      "wilson_b                : "+format_value("%-6.1f",self.data.wilson_b),
      "number_of_reflections   : "+format_value("%-8d",  self.data.number_of_reflections),
      "test_set_size           : "+format_value("%-8.4f",self.data.test_set_size),
      "test_flag_value         : "+format_value("%-d",   self.data.test_flag_value),
      "number_of_Fobs_outliers : "+format_value("%-8d",  self.data.number_of_Fobs_outliers),
      "twinned                 : "+format_value("%-s",   self.data.twinned),
      "anomalous_flag          : "+format_value("%-6s",  self.data.anomalous_flag)
      ])
    print >> log, "   ", result
    #
    print >> log, "  Model_vs_Data:"
    b_cart = " ".join([("%8.2f"%v).strip() for v in self.model_vs_data.b_cart])
    result = [
      "r_work(re-computed)                : %s"%format_value("%-6.4f",self.model_vs_data.r_work).strip(),
      "r_free(re-computed)                : %s"%format_value("%-6.4f",self.model_vs_data.r_free).strip(),
      "bulk_solvent_(k_sol,b_sol)         : %s %s"%(format_value("%-5.2f",self.model_vs_data.k_sol),
                                                    format_value("%-7.2f",self.model_vs_data.b_sol)),
      "overall_anisotropic_scale_(b_cart) : %-s"%b_cart]
    sc = self.model_vs_data.solvent_content_via_mask
    if (sc is not None): sc *= 100
    result.append("solvent_content_estimated_via_mask : %-s %%"
      % format_value("%.1f", sc))
    result = " \n    ".join(result)
    print >> log, "   ", result
    #
    if(self.pdb_header is not None):
      print >> log, "  Information extracted from PDB file header:"
      print >> log, "    program_name    : %-s"%format_value("%s",self.pdb_header.program_name)
      print >> log, "    year            : %-s"%format_value("%s",self.pdb_header.year)
      print >> log, "    r_work          : %-s"%format_value("%s",self.pdb_header.r_work)
      print >> log, "    r_free          : %-s"%format_value("%s",self.pdb_header.r_free)
      print >> log, "    high_resolution : %-s"%format_value("%s",self.pdb_header.high_resolution)
      print >> log, "    low_resolution  : %-s"%format_value("%s",self.pdb_header.low_resolution)
      print >> log, "    sigma_cutoff    : %-s"%format_value("%s",self.pdb_header.sigma_cutoff)
      print >> log, "    matthews_coeff  : %-s"%format_value("%s",self.pdb_header.matthews_coeff)
      if(self.pdb_header.solvent_cont is not None):
        print >> log, "    solvent_cont    : %-s %%"%format_value("%s",self.pdb_header.solvent_cont)
      else:
        print >> log, "    solvent_cont    : %-s"%format_value("%s",self.pdb_header.solvent_cont)
      if(self.pdb_header.tls is not None):
        print >> log, "    TLS             : %-s"%format_value("%s",
          " ".join([str(self.pdb_header.tls.pdb_inp_tls.tls_present),
          "(number of groups: %s)"%
          str(len(self.pdb_header.tls.tls_selections))]))
      else:
        print >> log, "    TLS             : None"
    #
    print >> log, "  After applying resolution and sigma cutoffs:"
    print >> log, "    n_refl_cutoff : %-s"%format_value("%d",self.misc.n_refl_cutoff).strip()
    print >> log, "    r_work_cutoff : %-s"%format_value("%6.4f",self.misc.r_work_cutoff).strip()
    print >> log, "    r_free_cutoff : %-s"%format_value("%6.4f",self.misc.r_free_cutoff).strip()

def get_program_name(file_lines):
  result = None
  for line in file_lines:
    line = line.strip()
    result = iotbx.pdb.remark_3_interpretation.get_program(st = line)
    if(result is not None): return result
  if(result is not None):
    result = "_".join(result.split())
  return result

def get_solvent_content(file_lines):
  mc = []
  for remark in file_lines:
    remark = remark.upper()
    if(remark.count("SOLVENT")==1 and
       remark.count("CONTENT")==1 and
       remark.count("REMARK 280")==1):
      try:
        mc.append(remark.split()[6])
      except:
        try:
          mc.append(remark[remark.index(":")+1:])
        except:
          mc.append(remark)
  result = None
  if(len(mc) == 1):
    try: result = float(mc[0])
    except IndexError: pass
    except ValueError: pass
  return result

def get_matthews_coeff(file_lines):
  mc = []
  for remark in file_lines:
    remark = remark.upper()
    if(remark.count("MATTHEWS")==1 and
       remark.count("COEFFICIENT")==1 and
       remark.count("REMARK 280")==1):
      try:
        mc.append(remark.split()[6])
      except:
        try:
          mc.append(remark[remark.index(":")+1:])
        except:
          mc.append(remark)
  result = None
  if(len(mc) == 1):
    try: result = float(mc[0])
    except IndexError: pass
    except ValueError: pass
  return result

def show_geometry(processed_pdb_file, scattering_table, hierarchy,
                  model_selections, show_geometry_statistics, mvd_obj):
  xray_structures = processed_pdb_file.xray_structure()
  hd_sel_all = xray_structures.hd_selection()
  if(show_geometry_statistics):
    sctr_keys = \
      xray_structures.scattering_type_registry().type_count_dict().keys()
    has_hd = "H" in sctr_keys or "D" in sctr_keys
    geometry = processed_pdb_file.geometry_restraints_manager(
      show_energies                = False,
      plain_pairs_radius           = 5.0,
      assume_hydrogens_all_missing = not has_hd)
    restraints_manager_all = mmtbx.restraints.manager(
      geometry      = geometry,
      normalization = True)
  models = hierarchy.models()
  geometry_statistics = []
  for i_seq, model_selection in enumerate(model_selections):
    hierarchy_i_seq = pdb.hierarchy.root()
    hierarchy_i_seq.append_model(models[i_seq].detached_copy())
    overall_counts_i_seq = hierarchy_i_seq.overall_counts()
    xray_structure = xray_structures.select(model_selection)
    hd_sel = xray_structure.hd_selection()
    if(hd_sel.count(True) > 0 and scattering_table != "neutron"):
      xray_structure_stat = show_xray_structure_statistics(
        xray_structure = xray_structure,
        hd_sel = ~hd_sel)
    else:
      xray_structure_stat = show_xray_structure_statistics(
        xray_structure = xray_structure)
    model_statistics_geometry = None
    ramalyze_obj = None
    rotalyze_obj = None
    cbetadev_obj = None
    clashscore_obj = None
    if(show_geometry_statistics):
      # exclude hydrogens
      if(hd_sel.count(True) > 0 and scattering_table != "neutron"):
        xray_structure = xray_structure.select(~hd_sel)
        model_selection = model_selection.select(~hd_sel)
        geometry = restraints_manager_all.geometry.select(selection = ~hd_sel_all)
      model_selection_as_bool = flex.bool(xray_structures.scatterers().size(),
        model_selection)
      geometry = restraints_manager_all.geometry.select(selection =
        model_selection_as_bool)
      restraints_manager = mmtbx.restraints.manager(
        geometry      = geometry,
        normalization = True)
      restraints_manager.geometry.pair_proxies(sites_cart =
        xray_structure.sites_cart())
      model_statistics_geometry = model_statistics.geometry(
        sites_cart         = xray_structure.sites_cart(),
        hd_selection       = hd_sel,
        ignore_hd          = False,
        restraints_manager = restraints_manager)
      #
      from mmtbx.validation.ramalyze import ramalyze
      from mmtbx.validation.rotalyze import rotalyze
      from mmtbx.validation.cbetadev import cbetadev
      from mmtbx.validation.clashscore import clashscore
      need_ramachandran = False
      ramalyze_obj = None
      rotalyze_obj = None
      cbetadev_obj = None
      clashscore_obj = None
      rc = overall_counts_i_seq.resname_classes
      n_residues = 0
      for k in rc.keys():
        if(k.count('amino_acid')):
          need_ramachandran = True
          n_residues = int(rc[k])
          break
      if(need_ramachandran):
        ramalyze_obj = ramalyze()
        output, output_list = ramalyze_obj.analyze_pdb(hierarchy =
          hierarchy_i_seq, outliers_only = False)
        rotalyze_obj = rotalyze()
        output_rotalyze, output_list_rotalyze = rotalyze_obj.analyze_pdb(
          hierarchy = hierarchy_i_seq, outliers_only = False)
        cbetadev_obj = cbetadev()
        output_cbetadev, sum_cbtadev, output_list_cbetadev = \
          cbetadev_obj.analyze_pdb(hierarchy = hierarchy_i_seq,
          outliers_only = False)
        clashscore_obj = clashscore()
        output_clashscore, bad_clashes = clashscore_obj.analyze_clashes(
          hierarchy = hierarchy_i_seq)
    geometry_statistics.append(group_args(
      overall_counts_i_seq      = overall_counts_i_seq,
      xray_structure_stat       = xray_structure_stat,
      model_statistics_geometry = model_statistics_geometry,
      ramalyze                  = ramalyze_obj,
      rotalyze                  = rotalyze_obj,
      cbetadev                  = cbetadev_obj,
      clashscore                = clashscore_obj))
  mvd_obj.collect(models = geometry_statistics)
  return geometry_statistics

def show_xray_structure_statistics(xray_structure, hd_sel = None):
  atom_counts = xray_structure.scattering_types_counts_and_occupancy_sums()
  atom_counts_strs = []
  for ac in atom_counts:
    atom_counts_strs.append("%s:%s:%s"%(ac.scattering_type,str(ac.count),
      str("%10.2f"%ac.occupancy_sum).strip()))
  atom_counts_str = " ".join(atom_counts_strs)
  if(hd_sel is not None):
    xray_structure = xray_structure.select(~hd_sel)
  b_isos = xray_structure.extract_u_iso_or_u_equiv()
  n_aniso = xray_structure.use_u_aniso().count(True)
  n_not_positive_definite = xray_structure.is_positive_definite_u().count(False)
  b_mean = format_value("%-6.1f",adptbx.u_as_b(flex.mean(b_isos))).strip()
  b_min = format_value("%-6.1f",adptbx.u_as_b(flex.min(b_isos))).strip()
  b_max = format_value("%-6.1f",adptbx.u_as_b(flex.max(b_isos))).strip()
  n_atoms = format_value("%-8d",xray_structure.scatterers().size()).strip()
  n_npd = format_value("%-8s",n_not_positive_definite).strip()
  occ = xray_structure.scatterers().extract_occupancies()
  o_mean = format_value("%-6.2f",flex.mean(occ)).strip()
  o_min = format_value("%-6.2f",flex.min(occ)).strip()
  o_max = format_value("%-6.2f",flex.max(occ)).strip()
  return group_args(n_atoms         = n_atoms,
                    atom_counts_str = atom_counts_str,
                    b_min           = b_min,
                    b_max           = b_max,
                    b_mean          = b_mean,
                    o_min           = o_min,
                    o_max           = o_max,
                    o_mean          = o_mean,
                    n_aniso         = n_aniso,
                    n_npd           = n_npd)

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def show_data(fmodel, n_outl, test_flag_value, f_obs_labels, fmodel_cut):
  info = fmodel.info()
  flags_pc = \
   fmodel.r_free_flags.data().count(True)*1./fmodel.r_free_flags.data().size()
  twinned = str(fmodel_cut.twin)
  if(fmodel_cut.twin != fmodel.twin):
    twinned = "May be, %s or %s"%(str(fmodel_cut.twin), str(fmodel.twin))
  return group_args(data_label              = f_obs_labels,
                    high_resolution         = info.d_min,
                    low_resolution          = info.d_max,
                    completeness_in_range   = info.completeness_in_range,
                    completeness_d_min_inf  = info.completeness_d_min_inf,
                    completeness_6A_inf     = info.completeness_6_inf,
                    wilson_b                = fmodel.wilson_b(),
                    number_of_reflections   = info.number_of_reflections,
                    test_set_size           = flags_pc,
                    test_flag_value         = test_flag_value,
                    number_of_Fobs_outliers = n_outl,
                    twinned                 = twinned,
                    anomalous_flag          = fmodel.f_obs.anomalous_flag())

def show_model_vs_data(fmodel):
  d_max, d_min = fmodel.f_obs.d_max_min()
  flags_pc = fmodel.r_free_flags.data().count(True)*100./\
    fmodel.r_free_flags.data().size()
  if(flags_pc == 0): r_free = None
  else: r_free = fmodel.r_free()
  sc = None
  mm = getattr(fmodel, "mask_manager", None)
  if (mm is not None):
    sc = mm.solvent_content_via_mask
  return group_args(
    r_work = fmodel.r_work(),
    r_free = r_free,
    k_sol  = fmodel.k_sol(),
    b_sol  = fmodel.b_sol(),
    b_cart = fmodel.b_cart(),
    solvent_content_via_mask=sc)

def run(args,
        command_name             = "mmtbx.model_vs_data",
        show_geometry_statistics = True,
        model_size_max_atoms     = 80000,
        data_size_max_reflections= 1000000,
        unit_cell_max_dimension  = 800.,
        return_fmodel_and_pdb    = False):
  if(len(args) == 0): args = ["--help"]
  command_line = (iotbx_option_parser(
    usage="%s reflection_file pdb_file [options]" % command_name,
    description='Example: %s data.mtz model.pdb'%command_name)
    .option(None, "--f_obs_label",
      action="store",
      default=None,
      type="string",
      help="Label for F-obs (or I-obs).")
    .option(None, "--r_free_flags_label",
      action="store",
      default=None,
      type="string",
      help="Label for free R flags.")
    .option(None, "--scattering_table",
      action="store",
      default="n_gaussian",
      type="string",
      help="Choice for scattering table: n_gaussian (default) or wk1995 or it1992 or neutron.")
    .option(None, "--comprehensive",
      action="store",
      default=True,
      type="bool",
      help="Show detailed statistics per residue (map CC, etc).")
    .option(None, "--r_factor_per_reflection",
      action="store",
      default=False,
      type="bool",
      help="Output R-factor for each reflection individually.")
    .option(None, "--map",
      action="store",
      default="None",
      type="string",
      help="Map type string: [p][m]Fo+[q][D]Fc[kick][filled]. Examples: 2mFo-DFc, 3.2Fo-2.3Fc, Fc, anom, fo-fc_kick.")
    .option("--ignore_giant_models_and_datasets",
      action="store_true",
      help="Ignore too big models and data files to avoid potential memory problems.")
    ).process(args=args)
  if(command_line.options.scattering_table not in ["n_gaussian","wk1995",
     "it1992","neutron"]):
    raise Sorry("Incorrect scattering_table.")
  #
  mvd_obj = mvd()
  #
  map_type_obj = None
  if(command_line.options.map != str(None)):
    map_type_obj = mmtbx.map_names(map_name_string = command_line.options.map)
  #
  processed_args = utils.process_command_line_args(args = args,
    log = sys.stdout)
  reflection_files = processed_args.reflection_files
  if(len(reflection_files) == 0):
    raise Sorry("No reflection file found.")
  crystal_symmetry = processed_args.crystal_symmetry
  if(crystal_symmetry is None):
    raise Sorry("No crystal symmetry found.")
  if(len(processed_args.pdb_file_names) == 0):
    raise Sorry("No PDB file found.")
  pdb_file_names = processed_args.pdb_file_names
  #
  rfs = reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    reflection_files = reflection_files)
  parameters = utils.data_and_flags_master_params().extract()
  if(command_line.options.f_obs_label is not None):
    parameters.labels = command_line.options.f_obs_label
  if(command_line.options.r_free_flags_label is not None):
    parameters.r_free_flags.label = command_line.options.r_free_flags_label
  determine_data_and_flags_result = utils.determine_data_and_flags(
    reflection_file_server  = rfs,
    parameters              = parameters,
    data_parameter_scope    = "refinement.input.xray_data",
    flags_parameter_scope   = "refinement.input.xray_data.r_free_flags",
    data_description        = "X-ray data",
    keep_going              = True,
    log                     = StringIO())
  f_obs = determine_data_and_flags_result.f_obs
  number_of_reflections = f_obs.indices().size()
  if(command_line.options.ignore_giant_models_and_datasets and
     number_of_reflections > data_size_max_reflections):
    raise Sorry("Too many reflections: %d"%number_of_reflections)
  #
  max_unit_cell_dimension = max(f_obs.unit_cell().parameters()[:3])
  if(command_line.options.ignore_giant_models_and_datasets and
     max_unit_cell_dimension > unit_cell_max_dimension):
    raise Sorry("Too large unit cell (max dimension): %s"%
      str(max_unit_cell_dimension))
  #
  r_free_flags = determine_data_and_flags_result.r_free_flags
  test_flag_value = determine_data_and_flags_result.test_flag_value
  if(r_free_flags is None):
    r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
    test_flag_value=None
  #
  mmtbx_pdb_file = mmtbx.utils.pdb_file(
    pdb_file_names   = pdb_file_names,
    cif_objects      = processed_args.cif_objects,
    crystal_symmetry = crystal_symmetry,
    use_elbow        = show_geometry_statistics,
    log              = sys.stdout)
  mmtbx_pdb_file.set_ppf(stop_if_duplicate_labels = False)
  processed_pdb_file = mmtbx_pdb_file.processed_pdb_file
  pdb_raw_records = mmtbx_pdb_file.pdb_raw_records
  pdb_inp = mmtbx_pdb_file.pdb_inp
  #
  xsfppf = mmtbx.utils.xray_structures_from_processed_pdb_file(
    processed_pdb_file = processed_pdb_file,
    scattering_table   = command_line.options.scattering_table,
    d_min              = f_obs.d_min())
  xray_structures = xsfppf.xray_structures
  if(0): #XXX normalize occupancies if all models have occ=1 so the total=1
    n_models = len(xray_structures)
    for xrs in xray_structures:
      occ = xrs.scatterers().extract_occupancies()
      occ = occ/n_models
      xrs.set_occupancies(occ)
  model_selections = xsfppf.model_selections
  mvd_obj.collect(crystal = group_args(
    uc       = f_obs.unit_cell(),
    sg       = f_obs.crystal_symmetry().space_group_info().symbol_and_number(),
    n_sym_op = f_obs.crystal_symmetry().space_group_info().type().group().order_z(),
    uc_vol   = f_obs.unit_cell().volume()))
  #
  hierarchy = pdb_inp.construct_hierarchy()
  pdb_atoms = hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  geometry_statistics = show_geometry(
    processed_pdb_file       = processed_pdb_file,
    scattering_table         = command_line.options.scattering_table,
    hierarchy                = hierarchy,
    model_selections         = model_selections,
    show_geometry_statistics = show_geometry_statistics,
    mvd_obj                  = mvd_obj)
  #
  # Extract TLS
  pdb_tls = None
  pdb_inp_tls = mmtbx.tls.tools.tls_from_pdb_inp(
    remark_3_records = pdb_inp.extract_remark_iii_records(3),
    pdb_hierarchy    = hierarchy)
  pdb_tls = group_args(pdb_inp_tls           = pdb_inp_tls,
                       tls_selections        = [],
                       tls_selection_strings = [])
  # XXX no TLS + multiple models
  if(pdb_inp_tls.tls_present and pdb_inp_tls.error_string is None and
     len(xray_structures)==1):
    pdb_tls = mmtbx.tls.tools.extract_tls_from_pdb(
      pdb_inp_tls       = pdb_inp_tls,
      all_chain_proxies = mmtbx_pdb_file.processed_pdb_file.all_chain_proxies,
      xray_structure    = xsfppf.xray_structure_all)
    if(len(pdb_tls.tls_selections)==len(pdb_inp_tls.tls_params) and
       len(pdb_inp_tls.tls_params) > 0):
      xray_structures = [utils.extract_tls_and_u_total_from_pdb(
        f_obs          = f_obs,
        r_free_flags   = r_free_flags,
        xray_structure = xray_structures[0], # XXX no TLS + multiple models
        tls_selections = pdb_tls.tls_selections,
        tls_groups     = pdb_inp_tls.tls_params)]
  #
  fmodel = utils.fmodel_simple(xray_structures = xray_structures,
                               f_obs           = f_obs,
                               r_free_flags    = r_free_flags)
  n_outl = f_obs.data().size() - fmodel.f_obs.data().size()
  mvd_obj.collect(model_vs_data = show_model_vs_data(fmodel))
  #
  # Extract information from PDB file header and output (if any)
  pub_r_work       = None
  pub_r_free       = None
  pub_high         = None
  pub_low          = None
  pub_sigma        = None
  pub_program_name = None
  pub_solv_cont    = None
  pub_matthews     = None
  published_results = extract_rfactors_resolutions_sigma.extract(
    file_name = pdb_file_names[0])
  if(published_results is not None):
    pub_r_work = published_results.r_work
    pub_r_free = published_results.r_free
    pub_high   = published_results.high
    pub_low    = published_results.low
    pub_sigma  = published_results.sigma
  pub_program_name = get_program_name(file_lines = pdb_raw_records)
  pub_solv_cont    = get_solvent_content(file_lines = pdb_raw_records)
  pub_matthews     = get_matthews_coeff(file_lines = pdb_raw_records)
  mvd_obj.collect(pdb_header = group_args(
    program_name    = pub_program_name,
    year            = pdb_inp.extract_header_year(),
    r_work          = pub_r_work,
    r_free          = pub_r_free,
    high_resolution = pub_high,
    low_resolution  = pub_low,
    sigma_cutoff    = pub_sigma,
    matthews_coeff  = pub_matthews,
    solvent_cont    = pub_solv_cont,
    tls             = pdb_tls))
  #
  # Recompute R-factors using published cutoffs
  fmodel_cut = fmodel
  tmp_sel = flex.bool(fmodel.f_obs.data().size(), True)
  if(pub_sigma is not None and fmodel.f_obs.sigmas() is not None):
    tmp_sel &= fmodel.f_obs.data() > fmodel.f_obs.sigmas()*pub_sigma
  if(pub_high is not None and abs(pub_high-fmodel.f_obs.d_min()) > 0.03):
    tmp_sel &= fmodel.f_obs.d_spacings().data() > pub_high
  if(pub_low is not None and abs(pub_low-fmodel.f_obs.d_max_min()[0]) > 0.03):
    tmp_sel &= fmodel.f_obs.d_spacings().data() < pub_low
  if(tmp_sel.count(True) != tmp_sel.size() and tmp_sel.count(True) > 0):
    fmodel_cut = utils.fmodel_simple(
      xray_structures = xray_structures,
      f_obs           = fmodel.f_obs.select(tmp_sel),
      r_free_flags    = fmodel.r_free_flags.select(tmp_sel))
  mvd_obj.collect(misc = group_args(
    r_work_cutoff = fmodel_cut.r_work(),
    r_free_cutoff = fmodel_cut.r_free(),
    n_refl_cutoff = fmodel_cut.f_obs.data().size()))
  mvd_obj.collect(data =
    show_data(fmodel          = fmodel,
              n_outl          = n_outl,
              test_flag_value = test_flag_value,
              f_obs_labels    = f_obs.info().label_string(),
              fmodel_cut      = fmodel_cut))
  mvd_obj.show()
  if return_fmodel_and_pdb :
    mvd_obj.pdb_file = processed_pdb_file
    mvd_obj.fmodel = fmodel
  if(map_type_obj is not None):
    map_params = maps.map_and_map_coeff_master_params().fetch(
      maps.cast_map_coeff_params(map_type_obj)).extract()
    maps_obj = maps.compute_map_coefficients(fmodel = fmodel_cut, params =
      map_params.map_coefficients)
    fn = os.path.basename(processed_args.reflection_file_names[0])
    if(fn.count(".")):
      prefix = fn[:fn.index(".")]
    else: prefix= fn
    file_name = prefix+"_%s_map_coeffs.mtz"%map_type_obj.format()
    maps_obj.write_mtz_file(file_name = file_name)
  # report map cc
  if(command_line.options.comprehensive and not fmodel_cut.twin and
     fmodel_cut.xray_structure is not None):
    show_hydrogens = False
    if(fmodel_cut.f_calc().d_min() <= 1.0 or
       command_line.options.scattering_table == "neutron"): show_hydrogens=True
    real_space_correlation.simple(
      fmodel                = fmodel_cut,
      pdb_hierarchy         = hierarchy,
      map_1_name            = "Fc",
      map_2_name            = "2mFo-DFc",
      details_level         = "automatic",
      atom_radius           = None,
      number_of_grid_points = 100,
      show                  = True,
      log                   = None,
      show_hydrogens        = show_hydrogens,
      selection             = None,
      set_cc_to_zero_if_n_grid_points_less_than = 50,
      poor_cc_threshold                         = 0.7,
      poor_map_value_threshold                  = 1.0)
  #
  if(command_line.options.r_factor_per_reflection):
    fmodel_cut.r_work_per_reflection()
  return mvd_obj

def read_mvd_output(file_lines, name):
  unit_cell        = None
  space_group      = None
  unit_cell_volume = None
  number_of_models = None
  # XXX multiple models are not supported
  n_altloc = None
  amino_acid     = None
  rna_dna        = None
  water          = None
  small_molecule = None
  element        = None
  other          = None
  #
  n_atoms     = None
  atom_counts = None
  adp_min     = None
  adp_max     = None
  adp_mean    = None
  occ_min     = None
  occ_max     = None
  occ_mean    = None
  n_aniso     = None
  n_npd       = None
  #
  bonds_rmsd     = None
  bonds_max      = None
  bonds_cnt      = None
  angles_rmsd    = None
  angles_max     = None
  angles_cnt     = None
  dihedrals_rmsd = None
  dihedrals_max  = None
  dihedrals_cnt  = None
  chirality_rmsd = None
  chirality_max  = None
  chirality_cnt  = None
  planarity_rmsd = None
  planarity_max  = None
  planarity_cnt  = None
  non_bonded_min = None
  #
  rama_outliers = None
  rama_allowed  = None
  rama_favored  = None
  #
  rota_outl  = None
  cbeta_dev  = None
  clashscore = None
  #
  data_label      = None
  d_min           = None
  d_max           = None
  cmpl_in_range   = None
  cmpl_d_min_inf  = None
  cmpl_6A_inf     = None
  wilson_b        = None
  n_refl          = None
  test_set_size   = None
  test_flag_value = None
  n_fobs_outl     = None
  twinned         = None
  anom_flag       = None
  #
  r_work_re_computed = None
  r_free_re_computed = None
  k_sol              = None
  b_sol              = None
  b_cart             = None
  #
  program_name   = None
  year           = None
  r_work_pdb     = None
  r_free_pdb     = None
  d_min_pdb      = None
  d_max_pdb      = None
  sigma_cutoff   = None
  matthews_coeff = None
  solvent_cont   = None
  tls            = None
  n_tls_groups   = None
  #
  n_refl_cutoff = None
  r_work_cutoff = None
  r_free_cutoff = None
  ###
  def helper(x):
    if(x.lower()=="none"): return None
    if(x.isdigit()): return int(x)
    if(x.lower()=="true"): return True
    if(x.lower()=="false"): return False
    return float(x)
  for line in file_lines:
    x = line.strip()
    xs = x.split()
    if(x.startswith("Unit cell:        ")):
      vtmp = " ".join(xs[2:]).replace(","," ").replace("("," ").replace(")"," ")
      unit_cell        = vtmp
    if(x.startswith("Space group:      ")): space_group      = xs[2:]
    if(x.startswith("Unit cell volume: ")): unit_cell_volume = float(xs[3])
    if(x.startswith("Number of models: ")): number_of_models = int(xs[3])
    if(number_of_models == 1): # XXX multiple models are not supported
      if(x.startswith("Number of residues in alternative conformations:")): 2
      if(x.startswith("amino_acid")):     amino_acid     = int(xs[2])
      if(x.startswith("rna_dna")):        rna_dna        = int(xs[2])
      if(x.startswith("water")):          water          = int(xs[2])
      if(x.startswith("small_molecule")): small_molecule = int(xs[2])
      if(x.startswith("element")):        element        = int(xs[2])
      if(x.startswith("other")):          other          = int(xs[2])
      if(x.startswith("atom_number_(type:count:occ_sum) :")):
        n_atoms     = int(xs[2])
        atom_counts = xs[3:]
      if(x.startswith("ADP_(min,max,mean)               :")):
        adp_min     = float(xs[2])
        adp_max     = float(xs[3])
        adp_mean    = float(xs[4])
      if(x.startswith("occupancies_(min,max,mean)       :")):
        occ_min     = float(xs[2])
        occ_max     = float(xs[3])
        occ_mean    = float(xs[4])
      if(x.startswith("number_of_anisotropic            :")): n_aniso = int(xs[2])
      if(x.startswith("number_of_non_positive_definite  :")): n_npd   = int(xs[2])
      if(x.startswith("bonds            :")):
        bonds_rmsd     = float(xs[2])
        bonds_max      = float(xs[3])
        bonds_cnt      = int(xs[4])
      if(x.startswith("angles           :")):
        angles_rmsd    = float(xs[2])
        angles_max     = float(xs[3])
        angles_cnt     = int(xs[4])
      if(x.startswith("dihedrals        :")):
        dihedrals_rmsd = float(xs[2])
        dihedrals_max  = float(xs[3])
        dihedrals_cnt  = int(xs[4])
      if(x.startswith("chirality        :")):
        chirality_rmsd = float(xs[2])
        chirality_max  = float(xs[3])
        chirality_cnt  = int(xs[4])
      if(x.startswith("planarity        :")):
        planarity_rmsd = float(xs[2])
        planarity_max  = float(xs[3])
        planarity_cnt  = int(xs[4])
      if(x.startswith("non-bonded (min) :")): non_bonded_min = float(xs[3])
      if(x.startswith("outliers :")): rama_outliers = float(xs[3].replace("(",""))
      if(x.startswith("allowed  :")): rama_allowed  = float(xs[3].replace("(",""))
      if(x.startswith("favored  :")): rama_favored  = float(xs[3].replace("(",""))
      if(x.startswith("Rotamer outliers :")):        rota_outl = float(xs[3].replace("(",""))
      if(x.startswith("Cbeta deviations >0.25A :")): cbeta_dev = float(xs[4].replace("(",""))
      if(x.startswith("All-atom clashscore     :")): cbeta_dev = float(xs[3].replace("(",""))
      if(x.startswith("data_label              :")): data_label      = xs[2]
      if(x.startswith("high_resolution         :")): d_min           = float(xs[2])
      if(x.startswith("low_resolution          :")): d_max           = float(xs[2])
      if(x.startswith("completeness_in_range   :")): cmpl_in_range   = float(xs[2])
      if(x.startswith("completeness(d_min-inf) :")): cmpl_d_min_inf  = float(xs[2])
      if(x.startswith("completeness(6A-inf)    :")): cmpl_6A_inf     = float(xs[2])
      if(x.startswith("wilson_b                :")): wilson_b        = helper(xs[2])
      if(x.startswith("number_of_reflections   :")): n_refl          = float(xs[2])
      if(x.startswith("test_set_size           :")): test_set_size   = float(xs[2])
      if(x.startswith("test_flag_value         :")): test_flag_value = helper(xs[2])
      if(x.startswith("number_of_Fobs_outliers :")): n_fobs_outl     = float(xs[2])
      if(x.startswith("twinned                 :")):
        if(xs[2].lower()=="false"): twinned = False
        else:                       twinned = xs[2]
      if(x.startswith("anomalous_flag          :")):
        if(xs[2].lower()=="false"): anom_flag = False
        else:                       anom_flag = True
      if(x.startswith("r_work(re-computed)                :")): r_work_re_computed = float(xs[2])
      if(x.startswith("r_free(re-computed)                :")): r_free_re_computed = helper(xs[2])
      if(x.startswith("bulk_solvent_(k_sol,b_sol)         :")):
        k_sol, b_sol = float(xs[2]), float(xs[3])
      if(x.startswith("overall_anisotropic_scale_(b_cart) :")): b_cart = flex.double([float(i) for i in xs[2:]])
      if(x.startswith("program_name    :")):
        if(xs[2].lower()=="none"): program_name = None
        else:                      program_name = "_".join(xs[2:])
      if(x.startswith("year            :")): year           = helper(xs[2])
      if(x.startswith("r_work          :")): r_work_pdb     = helper(xs[2])
      if(x.startswith("r_free          :")): r_free_pdb     = helper(xs[2])
      if(x.startswith("high_resolution :")): d_min_pdb      = helper(xs[2])
      if(x.startswith("low_resolution  :")): d_max_pdb      = helper(xs[2])
      if(x.startswith("sigma_cutoff    :")): sigma_cutoff   = helper(xs[2])
      if(x.startswith("matthews_coeff  :")): matthews_coeff = helper(xs[2])
      if(x.startswith("solvent_cont    :")): solvent_cont   = helper(xs[2])
      if(x.startswith("TLS             :")):
        tls          = helper(xs[2])
        n_tls_groups = int(xs[6].replace(")",""))
      if(x.startswith("n_refl_cutoff :")): n_refl_cutoff = helper(xs[2])
      if(x.startswith("r_work_cutoff :")): r_work_cutoff = helper(xs[2])
      if(x.startswith("r_free_cutoff :")): r_free_cutoff = helper(xs[2])

  ###
  return group_args(
    name             = name            ,
    unit_cell        = unit_cell       ,
    space_group      = space_group     ,
    unit_cell_volume = unit_cell_volume,
    number_of_models = number_of_models,
    n_altloc = n_altloc,
    amino_acid     = amino_acid    ,
    rna_dna        = rna_dna       ,
    water          = water         ,
    small_molecule = small_molecule,
    element        = element       ,
    other          = other         ,
    n_atoms     = n_atoms    ,
    atom_counts = atom_counts,
    adp_min     = adp_min    ,
    adp_max     = adp_max    ,
    adp_mean    = adp_mean   ,
    occ_min     = occ_min    ,
    occ_max     = occ_max    ,
    occ_mean    = occ_mean   ,
    n_aniso     = n_aniso    ,
    n_npd       = n_npd      ,
    bonds_rmsd     = bonds_rmsd    ,
    bonds_max      = bonds_max     ,
    angles_rmsd    = angles_rmsd   ,
    angles_max     = angles_max    ,
    dihedrals_rmsd = dihedrals_rmsd,
    dihedrals_max  = dihedrals_max ,
    chirality_rmsd = chirality_rmsd,
    chirality_max  = chirality_max ,
    planarity_rmsd = planarity_rmsd,
    planarity_max  = planarity_max ,
    non_bonded_min = non_bonded_min,
    rama_outliers = rama_outliers,
    rama_allowed  = rama_allowed ,
    rama_favored  = rama_favored ,
    data_label      = data_label     ,
    d_min           = d_min          ,
    d_max           = d_max          ,
    cmpl_in_range   = cmpl_in_range  ,
    cmpl_d_min_inf  = cmpl_d_min_inf ,
    cmpl_6A_inf     = cmpl_6A_inf    ,
    wilson_b        = wilson_b       ,
    n_refl          = n_refl         ,
    test_set_size   = test_set_size  ,
    test_flag_value = test_flag_value,
    n_fobs_outl     = n_fobs_outl    ,
    twinned         = twinned        ,
    anom_flag       = anom_flag      ,
    r_work_re_computed = r_work_re_computed,
    r_free_re_computed = r_free_re_computed,
    k_sol              = k_sol             ,
    b_sol              = b_sol             ,
    b_cart             = b_cart            ,
    program_name   = program_name  ,
    year           = year          ,
    r_work_pdb     = r_work_pdb    ,
    r_free_pdb     = r_free_pdb    ,
    d_min_pdb      = d_min_pdb     ,
    d_max_pdb      = d_max_pdb     ,
    sigma_cutoff   = sigma_cutoff  ,
    matthews_coeff = matthews_coeff,
    solvent_cont   = solvent_cont  ,
    tls            = tls           ,
    n_tls_groups   = n_tls_groups  ,
    n_refl_cutoff = n_refl_cutoff,
    r_work_cutoff = r_work_cutoff,
    r_free_cutoff = r_free_cutoff)
