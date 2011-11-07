import sys, os, random, re
from cctbx.array_family import flex
from iotbx import pdb
from cctbx import adptbx, sgtbx
from libtbx.utils import Sorry
from iotbx import reflection_file_utils
from mmtbx import real_space_correlation
from iotbx import reflection_file_utils
from libtbx.str_utils import format_value
import iotbx
from mmtbx import utils
from iotbx import pdb
from libtbx import easy_pickle
from cStringIO import StringIO
from mmtbx import model_statistics
from iotbx.pdb import extract_rfactors_resolutions_sigma
import iotbx.pdb.remark_3_interpretation
from libtbx import group_args
import mmtbx.restraints
import mmtbx.find_peaks
import mmtbx.maps
import mmtbx.masks

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
    self.maps          = None

  def collect(self,
              crystal       = None,
              models        = None,
              data          = None,
              model_vs_data = None,
              pdb_header    = None,
              misc          = None,
              maps          = None):
    if(crystal       is not None): self.crystal       = crystal
    if(models        is not None): self.models        = models
    if(data          is not None): self.data          = data
    if(model_vs_data is not None): self.model_vs_data = model_vs_data
    if(pdb_header    is not None): self.pdb_header    = pdb_header
    if(misc          is not None): self.misc          = misc
    if(maps          is not None): self.maps          = maps

  def get_summary (self) :
    return summarize_results(self)

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
        i_model.n_residues_in_altlocs
      print >> log, "    Residue content:"
      len_max = 0
      for k in i_model.resname_classes:
        k = k.split()
        if(len(k[0])>len_max): len_max = len(k[0])
      fmt = "      %-"+str(len_max)+"s : %s"
      for rc in i_model.resname_classes:
        k,v = rc.split()
        print >> log, fmt%(k, v)
      x = i_model.xray_structure_stat
      print >> log, "    Atoms:"
      print >> log, "      content: %s"%x.all.n_atoms
      max_field_for_element_count = 1
      for element_occupancy_sum in x.all.atom_counts_str.split():
        individual_element_count_len = len(element_occupancy_sum.split(":")[1])
        if(individual_element_count_len>max_field_for_element_count):
          max_field_for_element_count = individual_element_count_len
      for element_occupancy_sum in x.all.atom_counts_str.split():
        iecl = "%" + str(max_field_for_element_count)+"s"
        fmt = "        %s "+"count: %s "%iecl+"occupancy sum: %s"
        print >> log, fmt%tuple(element_occupancy_sum.split(":"))
      print >> log, "      ADP (min,max,mean):"
      print >> log, "        all           (%s atoms): %s %s %s"%(x.all.n_atoms, x.all.b_min,x.all.b_max,x.all.b_mean)

      if(x.sidechain     is not None):print >> log, "        side chains   (%s atoms): %s %s %s"%(x.sidechain.n_atoms, x.sidechain.b_min,x.sidechain.b_max,x.sidechain.b_mean)
      if(x.backbone      is not None):print >> log, "        main chains   (%s atoms): %s %s %s"%(x.backbone.n_atoms, x.backbone.b_min,x.backbone.b_max,x.backbone.b_mean)
      if(x.macromolecule is not None):print >> log, "        macromolecule (%s atoms): %s %s %s"%(x.macromolecule.n_atoms, x.macromolecule.b_min,x.macromolecule.b_max,x.macromolecule.b_mean)
      if(x.ligand        is not None):print >> log, "        ligands       (%s atoms): %s %s %s"%(x.ligand.n_atoms, x.ligand.b_min,x.ligand.b_max,x.ligand.b_mean)
      if(x.solvent       is not None):print >> log, "        solvent       (%s atoms): %s %s %s"%(x.solvent.n_atoms, x.solvent.b_min,x.solvent.b_max,x.solvent.b_mean)
      if(i_model.rms_b_iso_or_b_equiv_bonded is not None):
        print >> log, "      mean bonded (Bi-Bj) : %s"%format_value("%8.2f",i_model.rms_b_iso_or_b_equiv_bonded).strip()
      print >> log, "      occupancies (min,max,mean)       : %s %s %s"%(x.all.o_min,x.all.o_max,x.all.o_mean)
      print >> log, "      number_of_anisotropic            : "+format_value("%-7s",x.all.n_aniso)
      print >> log, "      number_of_non_positive_definite  : %s"%x.all.n_npd
      g = i_model.geometry_all
      if(g is not None):
        print >> log, "    Stereochemistry statistics (mean, max, count) - overall:"
        print >> log, "      bonds            : %8.4f %8.4f %d" % (g.b_mean, g.b_max, g.b_number)
        print >> log, "      angles           : %8.4f %8.4f %d" % (g.a_mean, g.a_max, g.a_number)
        print >> log, "      dihedrals        : %8.4f %8.4f %d" % (g.d_mean, g.d_max, g.d_number)
        print >> log, "      chirality        : %8.4f %8.4f %d" % (g.c_mean, g.c_max, g.c_number)
        print >> log, "      planarity        : %8.4f %8.4f %d" % (g.p_mean, g.p_max, g.p_number)
        print >> log, "      non-bonded (min) : %8.4f" % (g.n_min)
      if([i_model.geometry_solvent,i_model.geometry_ligand].count(None)==0):
        g = i_model.geometry_macromolecule
        if(g is not None):
          print >> log, "    Stereochemistry statistics (mean, max, count) - macromolecule:"
          print >> log, "      bonds            : %8.4f %8.4f %d" % (g.b_mean, g.b_max, g.b_number)
          print >> log, "      angles           : %8.4f %8.4f %d" % (g.a_mean, g.a_max, g.a_number)
          print >> log, "      dihedrals        : %8.4f %8.4f %d" % (g.d_mean, g.d_max, g.d_number)
          print >> log, "      chirality        : %8.4f %8.4f %d" % (g.c_mean, g.c_max, g.c_number)
          print >> log, "      planarity        : %8.4f %8.4f %d" % (g.p_mean, g.p_max, g.p_number)
          print >> log, "      non-bonded (min) : %8.4f" % (g.n_min)
        g = i_model.geometry_ligand
        if(g is not None):
          print >> log, "    Stereochemistry statistics (mean, max, count) - ligands:"
          print >> log, "      bonds            : %8.4f %8.4f %d" % (g.b_mean, g.b_max, g.b_number)
          print >> log, "      angles           : %8.4f %8.4f %d" % (g.a_mean, g.a_max, g.a_number)
          print >> log, "      dihedrals        : %8.4f %8.4f %d" % (g.d_mean, g.d_max, g.d_number)
          print >> log, "      chirality        : %8.4f %8.4f %d" % (g.c_mean, g.c_max, g.c_number)
          print >> log, "      planarity        : %8.4f %8.4f %d" % (g.p_mean, g.p_max, g.p_number)
          print >> log, "      non-bonded (min) : %8.4f" % (g.n_min)
        g = i_model.geometry_solvent
        if(g is not None):
          print >> log, "    Stereochemistry statistics - solvent:"
          print >> log, "      non-bonded (min) : %8.4f" % (g.n_min)
      if(i_model.molprobity is not None):
        outl = i_model.molprobity.ramalyze_outliers
        allo = i_model.molprobity.ramalyze_allowed
        favo = i_model.molprobity.ramalyze_favored
        print >> log, "    Molprobity statistics:"
        print >> log, "      Ramachandran plot, number of:"
        print >> log, "        outliers : %-5d (%-5.2f %s)"%(outl[0],outl[1]*100.,"%")
        print >> log, "        allowed  : %-5d (%-5.2f %s)"%(allo[0],allo[1]*100.,"%")
        print >> log, "        favored  : %-5d (%-5.2f %s)"%(favo[0],favo[1]*100.,"%")
      if(i_model.molprobity is not None):
        print >> log, "      Rotamer outliers        : %d (%s %s)" %(
          i_model.molprobity.rotalyze[0],
          str("%6.2f"%(i_model.molprobity.rotalyze[1]*100.)).strip(),"%")
      if(i_model.molprobity is not None):
        print >> log, "      Cbeta deviations >0.25A : %d"%i_model.molprobity.cbetadev
      if(i_model.molprobity is not None):
        print >> log, "      All-atom clashscore     : %.2f (steric overlaps >0.4A per 1000 atoms)"% \
          i_model.molprobity.clashscore
    #
    print >> log, "  Data:"
    result = " \n    ".join([
      "data_label                           : %s"%                    self.data.data_label,
      "high_resolution                      : "+format_value("%-5.2f",self.data.high_resolution),
      "low_resolution                       : "+format_value("%-6.2f",self.data.low_resolution),
      "completeness_in_range                : "+format_value("%-6.2f",self.data.completeness_in_range),
      "completeness(d_min-inf)              : "+format_value("%-6.2f",self.data.completeness_d_min_inf),
      "completeness(6A-inf)                 : "+format_value("%-6.2f",self.data.completeness_6A_inf),
      "wilson_b                             : "+format_value("%-6.1f",self.data.wilson_b),
      "number_of_reflections                : "+format_value("%-8d",  self.data.number_of_reflections),
      "number_of_reflections(non-anomalous) : "+format_value("%-8d",  self.data.number_of_reflections_merged),
      "test_set_size                        : "+format_value("%-8.4f",self.data.test_set_size),
      "test_flag_value                      : "+format_value("%-s",   self.data.test_flag_value),
      "number_of_Fobs_outliers              : "+format_value("%-8d",  self.data.number_of_Fobs_outliers),
      "twinned                              : "+format_value("%-s",   self.data.twinned),
      "anomalous_flag                       : "+format_value("%-6s",  self.data.anomalous_flag)
      ])
    print >> log, "   ", result
    #
    print >> log, "  Model_vs_Data:"
    b_cart = " ".join([("%8.2f"%v).strip() for v in self.model_vs_data.b_cart])
    result = [
      "r_work(re-computed)                : %s"%format_value("%-6.4f",self.model_vs_data.r_work).strip(),
      "r_free(re-computed)                : %s"%format_value("%-6.4f",self.model_vs_data.r_free).strip(),
      "bulk_solvent_(k_sol,b_sol)         : %s %s"%(" ".join(["%-5.2f"%i for i in self.model_vs_data.k_sol]).strip(),
                                                    format_value("%-7.2f",self.model_vs_data.b_sol)),
      "overall_anisotropic_scale_(b_cart) : %-s"%b_cart]
    sc = self.model_vs_data.solvent_content_via_mask
    if (sc is not None): sc *= 100
    result.append("solvent_content_estimated_via_mask : %-s %%"
      % format_value("%.1f", sc))
    result = " \n    ".join(result)
    print >> log, "   ", result
    if(self.maps is not None):
      print >> log, "    mFo-DFc map: positive and negative peak numbers:"
      print >> log, "      >  3 sigma: ", self.maps.peaks_plus_3
      print >> log, "      >  6 sigma: ", self.maps.peaks_plus_6
      print >> log, "      >  9 sigma: ", self.maps.peaks_plus_9
      print >> log, "      < -3 sigma: ", self.maps.peaks_minus_3
      print >> log, "      < -6 sigma: ", self.maps.peaks_minus_6
      print >> log, "      < -9 sigma: ", self.maps.peaks_minus_9
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
    #
    if(self.model_vs_data.sigmaa_plot is not None):
      print >> log, "  Statistics in resolution bins:"
      print >> log, "    SIGMAA vs Resolution:"
      print >> log, "     resolution(A)  sigmaa"
      for resolution, sigmaa in zip(self.model_vs_data.sigmaa_plot.resolution,
                                     self.model_vs_data.sigmaa_plot.sigmaa):
        print >> log, "        %10.3f%8.3f"%(float(resolution), float(sigmaa))

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
      except Exception:
        try:
          mc.append(remark[remark.index(":")+1:])
        except Exception:
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
      except Exception:
        try:
          mc.append(remark[remark.index(":")+1:])
        except Exception:
          mc.append(remark)
  result = None
  if(len(mc) == 1):
    try: result = float(mc[0])
    except IndexError: pass
    except ValueError: pass
  return result

def molprobity_stats(model_statistics_geometry, resname_classes):
  result = None
  need_ramachandran = False
  ramalyze_obj = None
  rotalyze_obj = None
  cbetadev_obj = None
  clashscore_obj = None
  rc = resname_classes
  n_residues = 0
  for k in rc.keys():
    if(k.count('amino_acid')):
      need_ramachandran = True
      n_residues = int(rc[k])
      break
  if(need_ramachandran):
    msg = model_statistics_geometry
    return group_args(
      ramalyze_outliers = msg.ramalyze_obj.get_outliers_count_and_fraction(),
      ramalyze_allowed  = msg.ramalyze_obj.get_allowed_count_and_fraction(),
      ramalyze_favored  = msg.ramalyze_obj.get_favored_count_and_fraction(),
      rotalyze          = msg.rotalyze_obj.get_outliers_count_and_fraction(),
      cbetadev          = msg.cbetadev_obj.get_outlier_count(),
      clashscore        = msg.clashscore_obj.get_clashscore())
  else: return None

def show_geometry(xray_structures, processed_pdb_file, scattering_table, hierarchy,
                  model_selections, show_geometry_statistics, mvd_obj,
                  atom_selections):
  if(len(xray_structures)>1):
    tmp = xray_structures[0]
    for xi in xray_structures[1:]:
      tmp = tmp.concatenate(xi)
    xray_structures = tmp
  else:
    xray_structures = xray_structures[0]
  ##
  utils.assert_xray_structures_equal(
    x1 = xray_structures,
    x2 = processed_pdb_file.xray_structure(),
    sites = True,
    adp = False,
    occupancies = True)
  ##
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
  n_residues_in_altlocs = None
  for i_seq, model_selection in enumerate(model_selections):
    hierarchy_i_seq = pdb.hierarchy.root()
    hierarchy_i_seq.append_model(models[i_seq].detached_copy())
    #
    overall_counts_i_seq = hierarchy_i_seq.overall_counts()
    n_residues_in_altlocs = \
      overall_counts_i_seq.n_alt_conf_pure + \
      overall_counts_i_seq.n_alt_conf_proper + \
      overall_counts_i_seq.n_alt_conf_improper
    #
    resname_classes = []
    for k,v in zip(overall_counts_i_seq.resname_classes.keys(),
                   overall_counts_i_seq.resname_classes.values()):
      resname_classes.append(" ".join([k.replace("common_",""), str(v)]))
    #
    xray_structure = xray_structures.select(model_selection)
    hd_sel = xray_structure.hd_selection()
    def select_atom_selections(selection       = model_selection,
                               atom_selections = atom_selections):
      result = group_args()
      result.all           = atom_selections.all          .select(selection)
      result.macromolecule = atom_selections.macromolecule.select(selection)
      result.solvent       = atom_selections.solvent      .select(selection)
      result.ligand        = atom_selections.ligand       .select(selection)
      result.backbone      = atom_selections.backbone     .select(selection)
      result.sidechain     = atom_selections.sidechain    .select(selection)
      return result
    atom_selections_i_model = select_atom_selections(selection = model_selection,
      atom_selections = atom_selections)
    if(hd_sel.count(True) > 0 and scattering_table != "neutron"):
      xray_structure_stat = show_xray_structure_statistics(
        xray_structure  = xray_structure,
        hd_sel          = hd_sel,
        atom_selections = atom_selections_i_model)
    else:
      xray_structure_stat = show_xray_structure_statistics(
        xray_structure = xray_structure,
        atom_selections=atom_selections_i_model)
    model_statistics_geometry_macromolecule = None
    model_statistics_geometry_solvent = None
    model_statistics_geometry_ligand = None
    model_statistics_geometry_all = None
    molprobity_stats_i_seq = None
    rms_b_iso_or_b_equiv_bonded = None
    if(show_geometry_statistics):
      # exclude hydrogens
      if(hd_sel.count(True) > 0 and scattering_table != "neutron"):
        xray_structure = xray_structure.select(~hd_sel)
        model_selection = model_selection.select(~hd_sel)
        geometry = restraints_manager_all.geometry.select(selection=~hd_sel_all)
        atom_selections_i_model = select_atom_selections(selection =~hd_sel_all,
           atom_selections = atom_selections_i_model)
      model_selection_as_bool = flex.bool(xray_structures.scatterers().size(),
        model_selection)
      geometry = restraints_manager_all.geometry.select(selection =
        model_selection_as_bool)
      restraints_manager = mmtbx.restraints.manager(
        geometry      = geometry,
        normalization = True)
      restraints_manager.geometry.pair_proxies(sites_cart =
        xray_structure.sites_cart())
      ###
      model_statistics_geometry_all = model_statistics.geometry(
        sites_cart         = xray_structure.sites_cart(),
        pdb_hierarchy      = hierarchy,
        hd_selection       = hd_sel,
        ignore_hd          = False,
        molprobity_scores  = True,
        restraints_manager = restraints_manager)
      #
      if(atom_selections.macromolecule.count(True)>0):
        mac_sel = atom_selections_i_model.macromolecule
        model_statistics_geometry_macromolecule = model_statistics.geometry(
          sites_cart         = xray_structure.sites_cart().select(mac_sel),
          pdb_hierarchy      = hierarchy,
          hd_selection       = xray_structure.select(mac_sel).hd_selection(),
          ignore_hd          = False,
          molprobity_scores  = True,
          restraints_manager = restraints_manager.select(mac_sel))
      #
      if(atom_selections.solvent.count(True)>0):
        sol_sel = atom_selections_i_model.solvent
        model_statistics_geometry_solvent = model_statistics.geometry(
          sites_cart         = xray_structure.sites_cart().select(sol_sel),
          pdb_hierarchy      = hierarchy,
          hd_selection       = xray_structure.select(sol_sel).hd_selection(),
          ignore_hd          = False,
          molprobity_scores  = True,
          restraints_manager = restraints_manager.select(sol_sel))
      #
      if(atom_selections.ligand.count(True)>0):
        lig_sel = atom_selections_i_model.ligand
        model_statistics_geometry_ligand = model_statistics.geometry(
          sites_cart         = xray_structure.sites_cart().select(lig_sel),
          pdb_hierarchy      = hierarchy,
          hd_selection       = xray_structure.select(lig_sel).hd_selection(),
          ignore_hd          = False,
          molprobity_scores  = True,
          restraints_manager = restraints_manager.select(lig_sel))
      ###
      rms_b_iso_or_b_equiv_bonded = utils.rms_b_iso_or_b_equiv_bonded(
        restraints_manager = restraints_manager,
        xray_structure     = xray_structure)
      #
      molprobity_stats_i_seq = molprobity_stats(
        model_statistics_geometry = model_statistics_geometry_all,
        resname_classes = overall_counts_i_seq.resname_classes)
    geometry_statistics.append(group_args(
      n_residues_in_altlocs       = n_residues_in_altlocs,
      resname_classes             = resname_classes,
      xray_structure_stat         = xray_structure_stat,
      rms_b_iso_or_b_equiv_bonded = rms_b_iso_or_b_equiv_bonded,
      geometry_all                = model_statistics_geometry_all,
      geometry_macromolecule      = model_statistics_geometry_macromolecule,
      geometry_solvent            = model_statistics_geometry_solvent,
      geometry_ligand             = model_statistics_geometry_ligand,
      molprobity                  = molprobity_stats_i_seq))
  mvd_obj.collect(models = geometry_statistics)
  return geometry_statistics

def show_xray_structure_statistics(xray_structure, atom_selections, hd_sel = None):
  result = group_args(
    all           = None,
    macromolecule = None,
    sidechain     = None,
    solvent       = None,
    ligand        = None,
    backbone      = None)
  if(hd_sel is not None):
    xray_structure = xray_structure.select(~hd_sel)
  for key in atom_selections.__dict__.keys():
    value = atom_selections.__dict__[key]
    if(value.count(True) > 0):
      if(hd_sel is not None):
        value = value.select(~hd_sel)
      xrs = xray_structure.select(value)
      atom_counts = xrs.scattering_types_counts_and_occupancy_sums()
      atom_counts_strs = []
      for ac in atom_counts:
        atom_counts_strs.append("%s:%s:%s"%(ac.scattering_type,str(ac.count),
          str("%10.2f"%ac.occupancy_sum).strip()))
      atom_counts_str = " ".join(atom_counts_strs)
      b_isos = xrs.extract_u_iso_or_u_equiv()
      n_aniso = xrs.use_u_aniso().count(True)
      n_not_positive_definite = xrs.is_positive_definite_u().count(False)
      b_mean = format_value("%-6.1f",adptbx.u_as_b(flex.mean(b_isos)))
      b_min = format_value("%-6.1f",adptbx.u_as_b(flex.min(b_isos)))
      b_max = format_value("%-6.1f",adptbx.u_as_b(flex.max(b_isos)))
      n_atoms = format_value("%-8d",xrs.scatterers().size()).strip()
      n_npd = format_value("%-8s",n_not_positive_definite).strip()
      occ = xrs.scatterers().extract_occupancies()
      o_mean = format_value("%-6.2f",flex.mean(occ)).strip()
      o_min = format_value("%-6.2f",flex.min(occ)).strip()
      o_max = format_value("%-6.2f",flex.max(occ)).strip()
      tmp_result = group_args(
        n_atoms         = n_atoms,
        atom_counts_str = atom_counts_str,
        b_min           = b_min,
        b_max           = b_max,
        b_mean          = b_mean,
        o_min           = o_min,
        o_max           = o_max,
        o_mean          = o_mean,
        n_aniso         = n_aniso,
        n_npd           = n_npd)
      setattr(result,key,tmp_result)
  return result

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def show_data(fmodel, n_outl, test_flag_value, f_obs_labels, fmodel_cut):
  info = fmodel.info()
  flags_pc = \
   fmodel.r_free_flags().data().count(True)*1./fmodel.r_free_flags().data().size()
  twinned = str(fmodel_cut.twin)
  if(fmodel_cut.twin != fmodel.twin):
    twinned = "May be, %s or %s"%(str(fmodel_cut.twin), str(fmodel.twin))
  return group_args(
    data_label                   = f_obs_labels,
    high_resolution              = info.d_min,
    low_resolution               = info.d_max,
    completeness_in_range        = info.completeness_in_range,
    completeness_d_min_inf       = info.completeness_d_min_inf,
    completeness_6A_inf          = info.completeness_6_inf,
    wilson_b                     = fmodel.wilson_b(),
    number_of_reflections        = info.number_of_reflections,
    number_of_reflections_merged = info.number_of_reflections_merged,
    test_set_size                = flags_pc,
    test_flag_value              = test_flag_value,
    number_of_Fobs_outliers      = n_outl,
    twinned                      = twinned,
    anomalous_flag               = fmodel.f_obs().anomalous_flag())

def show_model_vs_data(fmodel):
  d_max, d_min = fmodel.f_obs().d_max_min()
  flags_pc = fmodel.r_free_flags().data().count(True)*100./\
    fmodel.r_free_flags().data().size()
  if(flags_pc == 0): r_free = None
  else: r_free = fmodel.r_free()
  sc = None
  mm = getattr(fmodel, "mask_manager", None)
  if (mm is not None):
    sc = mm.solvent_content_via_mask
  sigmaa_plot = None
  if(not fmodel.twin):
    sigmaa_plot = fmodel.sigmaa().show_short(silent=True)
  return group_args(
    r_work                   = fmodel.r_work(),
    r_free                   = r_free,
    k_sol                    = fmodel.k_sols(),
    b_sol                    = fmodel.b_sol(),
    b_cart                   = fmodel.b_cart(),
    solvent_content_via_mask = sc,
    sigmaa_plot              = sigmaa_plot)

def maps(fmodel, mvd_obj, map_cutoff = 3.0, map_type = "mFo-DFc"):
  result = group_args(
    peaks_plus_3  = 0,
    peaks_plus_6  = 0,
    peaks_plus_9  = 0,
    peaks_minus_3 = 0,
    peaks_minus_6 = 0,
    peaks_minus_9 = 0)
  params = mmtbx.find_peaks.master_params.extract()
  params.peak_search.min_cross_distance = 1.2
  if(fmodel.twin and fmodel.r_free_flags().data().count(True)==0): # XXX
    rff = fmodel.r_free_flags().generate_r_free_flags()
    fmodel.update(r_free_flags = rff) # XXX
  peaks_plus = mmtbx.find_peaks.manager(
    fmodel     = fmodel,
    map_type   = map_type,
    map_cutoff = map_cutoff,
    params     = params,
    silent     = True).peaks()
  if(peaks_plus is not None): peaks_plus = peaks_plus.heights
  peaks_minus = mmtbx.find_peaks.manager(
    fmodel     = fmodel,
    map_type   = map_type,
    map_cutoff = -1.*map_cutoff,
    params     = params,
    silent     = True).peaks()
  if(peaks_minus is not None): peaks_minus = peaks_minus.heights
  if([peaks_minus,peaks_plus].count(None) == 0):
    result = group_args(
      peaks_plus_3  = (peaks_plus > 3).count(True),
      peaks_plus_6  = (peaks_plus > 6).count(True),
      peaks_plus_9  = (peaks_plus > 9).count(True),
      peaks_minus_3 = (peaks_minus < -3).count(True),
      peaks_minus_6 = (peaks_minus < -6).count(True),
      peaks_minus_9 = (peaks_minus < -9).count(True))
  return result

msg="""\

phenix.model_vs_data: compute model, data and model-to-data fit statistics.

Reference:
  phenix.model_vs_data: a high-level tool for the calculation of
  crystallographic model and data statistics. P.V. Afonine, R.W.
  Grosse-Kunstleve, V.B. Chen, J.J. Headd, N.W. Moriarty, J.S. Richardson,
  D.C. Richardson, A. Urzhumtsev, P.H. Zwart, P.D. Adams
  J. Appl. Cryst. 43, 677-685 (2010).

Inputs:
  - File with reflection data (Fobs or Iobs), and R-free flags (optionally);
  - label(s) selecting which reflection data arrays should be used (in case
    there are multiple choices in input file, there is no need to provide labels
    otherwise);
  - PDB file with input model;
  - some other optional parameters.

Usage examples:
  1. phenix.model_vs_data model.pdb data.hkl
  2. phenix.model_vs_data model.pdb data.hkl f_obs_label="F" r_free_flags_label="FREE"
  3. phenix.model_vs_data model.pdb data.hkl scattering_table=neutron
  4. phenix.model_vs_data model.pdb data.hkl map="2mFo-DFc"
  5. phenix.model_vs_data model.pdb data.hkl map="3Fo-2Fc" map="mFo-DFc"

  Note: Map type string: [p][m]Fo+[q][D]Fc[kick][filled]. Examples: 2mFo-DFc,
  3.2Fo-2.3Fc, Fc, anom, fo-fc_kick, etc.
"""

master_params_str="""\
f_obs_label = None
  .type = str
r_free_flags_label = None
  .type = str
scattering_table = wk1995  it1992  *n_gaussian  neutron
  .type = choice
map = None
  .type = str
  .multiple = True
high_resolution = None
  .type = float
comprehensive = False
  .type = bool
dump_result_object_as_pickle = False
  .type = bool
ignore_giant_models_and_datasets = True
  .type = bool
"""

def defaults(log, silent):
  if(not silent): print >> log, "Default params::\n"
  parsed = iotbx.phil.parse(master_params_str)
  if(not silent): parsed.show(prefix="  ", out=log)
  if(not silent): print >> log
  return parsed

def run(args,
        command_name             = "mmtbx.model_vs_data",
        show_geometry_statistics = True,
        model_size_max_atoms     = 80000,
        data_size_max_reflections= 1000000,
        unit_cell_max_dimension  = 800.,
        return_fmodel_and_pdb    = False,
        out                      = None,
        log                      = sys.stdout):
  import mmtbx.f_model_info
  if(len(args)==0):
    print >> log, msg
    defaults(log=log, silent=False)
    return
  parsed = defaults(log=log, silent=True)
  #
  mvd_obj = mvd()
  #
  processed_args = utils.process_command_line_args(args = args,
    log = sys.stdout, master_params = parsed)
  params = processed_args.params.extract()
  #
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
  if(params.f_obs_label is not None):
    parameters.labels = params.f_obs_label
  if(params.r_free_flags_label is not None):
    parameters.r_free_flags.label = params.r_free_flags_label
  if (params.high_resolution is not None) :
    parameters.high_resolution = params.high_resolution
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
  if(params.ignore_giant_models_and_datasets and
     number_of_reflections > data_size_max_reflections):
    raise Sorry("Too many reflections: %d"%number_of_reflections)
  #
  max_unit_cell_dimension = max(f_obs.unit_cell().parameters()[:3])
  if(params.ignore_giant_models_and_datasets and
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
  # just to avoid going any further with bad PDB file....
  pdb_inp.xray_structures_simple()
  #
  acp = processed_pdb_file.all_chain_proxies
  atom_selections = group_args(
    all           = acp.selection(string = "all"),
    macromolecule = acp.selection(string = "protein or dna or rna"),
    solvent       = acp.selection(string = "water"), # XXX single_atom_residue
    ligand        = acp.selection(string = "not (protein or dna or rna or water)"),
    backbone      = acp.selection(string = "backbone"),
    sidechain     = acp.selection(string = "sidechain"))
  #
  xsfppf = mmtbx.utils.xray_structures_from_processed_pdb_file(
    processed_pdb_file = processed_pdb_file,
    scattering_table   = params.scattering_table,
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
  ###########################
  geometry_statistics = show_geometry(
    xray_structures          = xray_structures,
    processed_pdb_file       = processed_pdb_file,
    scattering_table         = params.scattering_table,
    hierarchy                = hierarchy,
    model_selections         = model_selections,
    show_geometry_statistics = show_geometry_statistics,
    mvd_obj                  = mvd_obj,
    atom_selections          = atom_selections)
  ###########################
  mp = mmtbx.masks.mask_master_params.extract()
  fmodel = utils.fmodel_simple(
    xray_structures  = xray_structures,
    scattering_table = params.scattering_table,
    mask_params      = mp,
    f_obs            = f_obs,
    r_free_flags     = r_free_flags)
  #
  n_outl = f_obs.data().size() - fmodel.f_obs().data().size()
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
  tmp_sel = flex.bool(fmodel.f_obs().data().size(), True)
  if(pub_sigma is not None and fmodel.f_obs().sigmas() is not None):
    tmp_sel &= fmodel.f_obs().data() > fmodel.f_obs().sigmas()*pub_sigma
  if(pub_high is not None and abs(pub_high-fmodel.f_obs().d_min()) > 0.03):
    tmp_sel &= fmodel.f_obs().d_spacings().data() > pub_high
  if(pub_low is not None and abs(pub_low-fmodel.f_obs().d_max_min()[0]) > 0.03):
    tmp_sel &= fmodel.f_obs().d_spacings().data() < pub_low
  if(tmp_sel.count(True) != tmp_sel.size() and tmp_sel.count(True) > 0):
    fmodel_cut = utils.fmodel_simple(
      xray_structures  = xray_structures,
      scattering_table = params.scattering_table,
      f_obs            = fmodel.f_obs().select(tmp_sel),
      r_free_flags     = fmodel.r_free_flags().select(tmp_sel))
  mvd_obj.collect(misc = group_args(
    r_work_cutoff = fmodel_cut.r_work(),
    r_free_cutoff = fmodel_cut.r_free(),
    n_refl_cutoff = fmodel_cut.f_obs().data().size()))
  mvd_obj.collect(data =
    show_data(fmodel          = fmodel,
              n_outl          = n_outl,
              test_flag_value = test_flag_value,
              f_obs_labels    = f_obs.info().label_string(),
              fmodel_cut      = fmodel_cut))
  # map statistics
  if(len(xray_structures)==1): # XXX no multi-model support yet
    mvd_obj.collect(maps = maps(fmodel = fmodel, mvd_obj = mvd_obj))
  #
  mvd_obj.show(log=out)
  if return_fmodel_and_pdb :
    mvd_obj.pdb_file = processed_pdb_file
    mvd_obj.fmodel = fmodel
  if(len(params.map) > 0):
    for map_name_string in params.map:
      map_type_obj = mmtbx.map_names(map_name_string = map_name_string)
      map_params = mmtbx.maps.map_and_map_coeff_master_params().fetch(
        mmtbx.maps.cast_map_coeff_params(map_type_obj)).extract()
      maps_obj = mmtbx.maps.compute_map_coefficients(fmodel = fmodel_cut, params =
        map_params.map_coefficients)
      fn = os.path.basename(processed_args.reflection_file_names[0])
      if(fn.count(".")):
        prefix = fn[:fn.index(".")]
      else: prefix= fn
      file_name = prefix+"_%s_map_coeffs.mtz"%map_type_obj.format()
      maps_obj.write_mtz_file(file_name = file_name)
  # statistics in bins
  if(not fmodel.twin):
    mmtbx.f_model_info.r_work_and_completeness_in_resolution_bins(fmodel = fmodel,
      out = log)
  # report map cc
  if(params.comprehensive and not fmodel_cut.twin and
     fmodel_cut.xray_structure is not None):
    show_hydrogens = False
    if(fmodel_cut.f_calc().d_min() <= 1.0 or
       params.scattering_table == "neutron"): show_hydrogens=True
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
      poor_map_1_value_threshold                = 1.0,
      poor_map_2_value_threshold                = 1.0)
  #
  if(params.dump_result_object_as_pickle):
    output_prefixes = []
    for op in processed_args.pdb_file_names+processed_args.reflection_file_names:
      op = os.path.basename(op)
      try: op = op[:op.index(".")]
      except Exception: pass
      if(not op in output_prefixes): output_prefixes.append(op)
    output_prefix = "_".join(output_prefixes)
    easy_pickle.dump("%s.pickle"%output_prefix, mvd_obj)
  return mvd_obj

def summarize_results (mvd_obj) :
  space_group_str = getattr(mvd_obj.crystal, "sg", None)
  if (space_group_str is None) :
    space_group = None
  else :
    space_group = sgtbx.space_group_info(re.sub("\ \(.*", "", space_group_str))
  model_stats = mvd_obj.models[0]
  geometry_stats = getattr(model_stats, "geometry_all", None)
  molprobity_stats = getattr(model_stats, "molprobity", None)
  xs_stats = getattr(getattr(model_stats, "xray_structure_stat", None),
                     "all", None)
  mm_stats = getattr(getattr(model_stats, "xray_structure_stat", None),
      "macromolecule", None)
  bb_stats = getattr(getattr(model_stats, "xray_structure_stat", None),
    "backbone", None)
  if (molprobity_stats is not None) :
    c_beta_deviations = molprobity_stats.cbetadev
    clashscore = molprobity_stats.clashscore
    rama_allowed = molprobity_stats.ramalyze_allowed[1] * 100.0
    rama_favored = molprobity_stats.ramalyze_favored[1] * 100.0
    rama_outliers = molprobity_stats.ramalyze_outliers[1] * 100.0
    rotamer_outliers = molprobity_stats.rotalyze[1] * 100.0
  else:
    c_beta_deviations = clashscore = rama_allowed = rama_favored = \
      rama_outliers = rotamer_outliers =None
  def convert_float (value) :
    if (value is None) : return None
    try :
      return float(value)
    except ValueError :
      return None
  return group_args(
    space_group=space_group,
    unit_cell=getattr(mvd_obj.crystal, "uc", None),
    r_work=mvd_obj.model_vs_data.r_work,
    r_free=mvd_obj.model_vs_data.r_free,
    pdb_header_r_work=getattr(mvd_obj.pdb_header, "r_work", None),
    pdb_header_r_free=getattr(mvd_obj.pdb_header, "r_free", None),
    r_work_cutoffs=getattr(mvd_obj.misc, "r_work_cutoff", None),
    r_free_cutoffs=getattr(mvd_obj.misc, "r_free_cutoff", None),
    program_name=getattr(mvd_obj.pdb_header, "program_name", None),
    d_max_pdb=getattr(mvd_obj.pdb_header, "low_resolution", None),
    d_min_pdb=getattr(mvd_obj.pdb_header, "high_resolution", None),
    d_min=mvd_obj.data.high_resolution,
    d_max=mvd_obj.data.low_resolution,
    completeness_6A_inf=mvd_obj.data.completeness_6A_inf,
    completeness_d_min_inf=mvd_obj.data.completeness_d_min_inf,
    completeness_in_range=mvd_obj.data.completeness_in_range,
    adp_max_all=convert_float(getattr(xs_stats, "b_max",None)),
    adp_mean_all=convert_float(getattr(xs_stats, "b_mean",None)),
    adp_min_all=convert_float(getattr(xs_stats, "b_min", None)),
    adp_max_mm=convert_float(getattr(mm_stats, "b_max",None)),
    adp_mean_mm=convert_float(getattr(mm_stats, "b_mean",None)),
    adp_min_mm=convert_float(getattr(mm_stats, "b_min", None)),
    adp_max_bb=convert_float(getattr(bb_stats, "b_max",None)),
    adp_mean_bb=convert_float(getattr(bb_stats, "b_mean",None)),
    adp_min_bb=convert_float(getattr(bb_stats, "b_min", None)),
    n_aniso=getattr(xs_stats, "n_aniso", None),
    wilson_b=mvd_obj.data.wilson_b,
    b_sol=mvd_obj.model_vs_data.b_sol,
    k_sol=mvd_obj.model_vs_data.k_sol,
    solvent_content_via_mask=mvd_obj.model_vs_data.solvent_content_via_mask,
    bond_rmsd=getattr(geometry_stats, "b_mean", None),
    bond_max_deviation=getattr(geometry_stats, "b_max", None),
    angle_rmsd=getattr(geometry_stats, "a_mean", None),
    angle_max_deviation=getattr(geometry_stats, "a_max", None),
    dihedral_rmsd=getattr(geometry_stats, "d_mean", None),
    dihedral_max_deviation=getattr(geometry_stats, "d_max", None),
    chirality_rmsd=getattr(geometry_stats, "c_mean", None),
    chirality_max_deviation=getattr(geometry_stats, "c_max", None),
    planarity_rmsd=getattr(geometry_stats, "p_mean", None),
    planarity_max_deviation=getattr(geometry_stats, "p_max", None),
    rama_favored=rama_favored,
    rama_allowed=rama_allowed,
    rama_outliers=rama_outliers,
    rotamer_outliers=rotamer_outliers,
    c_beta_deviations=c_beta_deviations,
    clashscore=clashscore,
    twin_law=mvd_obj.data.twinned,
    tls=getattr(getattr(getattr(mvd_obj.pdb_header, "tls", None), "pdb_inp_tls",
      None), "tls_present"))
