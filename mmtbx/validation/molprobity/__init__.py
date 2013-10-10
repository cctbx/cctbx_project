
# TODO combine with some parts of mmtbx.kinemage.validation
# TODO incorporate unmerged data
# TODO merge this with Table 1 output, PDB deposition, etc.

from __future__ import division
from mmtbx.validation import validation, residue
from mmtbx.validation import model_properties
from mmtbx.validation import experimental
from mmtbx.validation import clashscore
from mmtbx.validation import restraints
from mmtbx.validation import ramalyze
from mmtbx.validation import rotalyze
from mmtbx.validation import cbetadev
from mmtbx.validation import waters
from libtbx.str_utils import make_header, make_sub_header, format_value
from libtbx import slots_getstate_setstate
from libtbx.utils import null_out
import libtbx.load_env
import os.path
import sys

class molprobity (slots_getstate_setstate) :
  """
  Comprehensive validation.  At a minimum this performs the standard MolProbity
  analyses (ramalyze, rotalyze, cbetadev, clashscore).  If a geometry
  restraints manager is available, the deviations from standard covalent
  geometry will also be displayed.  Passing an fmodel object enables the
  re-calculation of R-factors and real-space correlation.
  """

  __slots__ = [
    "ramalyze",
    "rotalyze",
    "cbetadev",
    "clashscore",
    "restraints",
    "missing_atoms",
    "data_stats",
    "real_space",
    #"real_space_atoms",
    "model_stats",
    "waters",
    "header_info",
    "merging",
    "_multi_criterion",
  ]

  def __init__ (self,
      pdb_hierarchy,
      xray_structure=None,
      fmodel=None,
      geometry_restraints_manager=None,
      header_info=None,
      unmerged_data=None,
      keep_hydrogens=True,
      nuclear=False,
      save_probe_unformatted_file=None,
      show_hydrogen_outliers=False,
      min_cc_two_fofc=0.8,
      outliers_only=True) :
    for name in self.__slots__ :
      setattr(self, name, None)
    self.header_info = header_info
    if pdb_hierarchy.contains_protein() :
      self.ramalyze = ramalyze.ramalyze(
        pdb_hierarchy=pdb_hierarchy,
        outliers_only=False,
        out=null_out(),
        quiet=True)
      self.rotalyze = rotalyze.rotalyze(
        pdb_hierarchy=pdb_hierarchy,
        outliers_only=False,
        out=null_out(),
        quiet=True)
      self.cbetadev = cbetadev.cbetadev(
        pdb_hierarchy=pdb_hierarchy,
        outliers_only=True,
        out=null_out(),
        quiet=True)
    self.clashscore = clashscore.clashscore(
      pdb_hierarchy=pdb_hierarchy,
      save_probe_unformatted_file=save_probe_unformatted_file,
      nuclear=nuclear,
      keep_hydrogens=keep_hydrogens,
      out=null_out(),
      verbose=False)
    self.model_stats = model_properties.model_statistics(
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      ignore_hd=(not nuclear))
    if (geometry_restraints_manager is not None) :
      assert (xray_structure is not None)
      self.restraints = restraints.combined(
        pdb_hierarchy=pdb_hierarchy,
        xray_structure=xray_structure,
        geometry_restraints_manager=geometry_restraints_manager,
        ignore_hd=(not nuclear))
    if (fmodel is not None) :
      self.data_stats = experimental.data_statistics(fmodel)
      self.waters = waters.waters(
        pdb_hierarchy=pdb_hierarchy,
        xray_structure=xray_structure,
        fmodel=fmodel,
        collect_all=False)
      self.real_space = experimental.real_space(
        fmodel=fmodel,
        pdb_hierarchy=pdb_hierarchy,
        cc_min=min_cc_two_fofc)
    self._multi_criterion = multi_criterion_view(pdb_hierarchy)

  def molprobity_score (self) :
    if (self.ramalyze is None) : return None
    from mmtbx.validation import utils
    return utils.molprobity_score(
      clashscore=self.clashscore.get_clashscore(),
      rota_out=self.rotalyze.out_percent,
      rama_fav=self.ramalyze.fav_percent)

  def show (self, out=sys.stdout, outliers_only=True) :
    """
    Comprehensive output with individual outlier lists, plus summary.
    """
    if (self.data_stats is not None) :
      make_header("Experimental data", out=out)
      self.data_stats.show(out=out, prefix="  ")
      if (self.real_space is not None) :
        make_sub_header("Residues with poor real-space CC", out=out)
        self.real_space.show(out=out, prefix="  ")
      if (self.waters is not None) :
        make_sub_header("Suspicious water molecules", out=out)
        self.waters.show(out=out, prefix="  ")
    if (self.restraints is not None) :
      make_header("Geometry restraints", out=out)
      self.restraints.show(out=out, prefix="  ")
    make_header("Molprobity validation", out=out)
    if (self.ramalyze is not None) :
      make_sub_header("Ramachandran angles", out=out)
      self.ramalyze.show(out=out, prefix="  ", outliers_only=outliers_only)
    if (self.rotalyze is not None) :
      make_sub_header("Sidechain rotamers", out=out)
      self.rotalyze.show(out=out, prefix="  ", outliers_only=outliers_only)
    if (self.cbetadev is not None) :
      make_sub_header("C-beta deviations", out=out)
      self.cbetadev.show(out=out, prefix="  ", outliers_only=outliers_only)
    make_sub_header("Bad clashes", out=out)
    self.clashscore.show(out=out, prefix="  ")
    make_header("Summary", out=out)
    self.show_summary(out=out, prefix="  ")

  def show_summary (self, out=sys.stdout, prefix="") :
    """
    Summarize outliers or scores for each analysis.
    """
    def fs (format, value) :
      return format_value(format, value, replace_none_with=("(none)"))
    if (self.ramalyze is not None) :
      print >> out, "%sRamachandran outliers = %s %%" % (prefix,
        fs("%6.2f", self.ramalyze.out_percent))
      print >> out, "%s             favored  = %s %%" % (prefix,
        fs("%6.2f", self.ramalyze.fav_percent))
    if (self.rotalyze is not None) :
      print >> out, "%sRotamer outliers      = %s %%" % (prefix,
        fs("%6.2f", self.rotalyze.out_percent))
    if (self.cbetadev is not None) :
      print >> out, "%sC-beta deviations     = %s" % (prefix,
        fs("%6d", self.cbetadev.n_outliers))
    print >> out, "%sClashscore            = %6.2f" % (prefix,
      self.clashscore.get_clashscore())
    mpscore = self.molprobity_score()
    if (mpscore is not None) :
      print >> out, "%sMolprobity score      = %6.2f" % (prefix, mpscore)
    if (self.restraints is not None) :
      rms_bonds, rms_angles = self.restraints.get_bonds_angles_rmsds()
      print >> out, "%sRMS(bonds)            = %8.4f" % (prefix, rms_bonds)
      print >> out, "%sRMS(angles)           = %6.2f" % (prefix, rms_angles)
    if (self.data_stats is not None) :
      self.data_stats.show_summary(prefix=prefix, out=out)
    if (self.restraints is None) or (self.data_stats is None) :
      if (self.header_info is not None) :
        self.header_info.show(out=out,
          prefix=prefix,
          include_r_factors=(self.data_stats is None),
          include_rms_geom=(self.restraints is None))

  def as_mmcif_records (self) : # TODO
    raise NotImplementedError()

  def as_table1 (self) : # TODO
    raise NotImplementedError()

  def as_multi_criterion_view (self) :
    if (not self._multi_criterion.is_populated) :
      if (self.real_space is not None) :
        self._multi_criterion.process_outliers(self.real_space.results)
      if (self.waters is not None) :
        self._multi_criterion.process_outliers(self.waters.results)
      if (self.ramalyze is not None) :
        self._multi_criterion.process_outliers(self.ramalyze.results)
      if (self.rotalyze is not None) :
        self._multi_criterion.process_outliers(self.rotalyze.results)
      if (self.cbetadev is not None) :
        self._multi_criterion.process_outliers(self.cbetadev.results)
      if (self.clashscore is not None) :
        self._multi_criterion.process_outliers(self.clashscore.results)
    return self._multi_criterion

  def display_wx_plots (self) :
    if (self.ramalyze is not None) :
      self.ramalyze.display_wx_plots()
    if (self.rotalyze is not None) :
      self.rotalyze.display_wx_plots()
    mc = self.as_multi_criterion_view()
    mc.display_wx_plots()

  def write_coot_script (self, file_name) :
    coot_script = libtbx.env.find_in_repositories(
      relative_path="cctbx_project/cootbx/validation_lists.py",
      test=os.path.isfile)
    if (coot_script is None) :
      raise Sorry("Can't find template Python script for Coot.")
    f = open(file_name, "w")
    f.write("# script auto-generated by phenix.molprobity\n")
    f.write("\n")
    f.write(open(coot_script).read())
    f.write("\n")
    f.write("data = {}\n")
    if (self.ramalyze is not None) :
      f.write("data['rama'] = %s\n" % self.ramalyze.as_coot_data())
    if (self.rotalyze is not None) :
      f.write("data['rota'] = %s\n" % self.rotalyze.as_coot_data())
    if (self.cbetadev is not None) :
      f.write("data['cbeta'] = %s\n" % self.cbetadev.as_coot_data())
    if (self.clashscore is not None) :
      f.write("data['probe'] = %s\n" % self.clashscore.as_coot_data())
      if (self.clashscore.probe_file is not None) :
        f.write("handle_read_draw_probe_dots_unformatted(\"%s\", 0, 0)\n" %
          self.clashscore.probe_file)
        f.write("show_probe_dots(True, True)\n")
    f.write("gui = coot_molprobity_todo_list_gui(data=data)\n")
    f.close()

class pdb_header_info (slots_getstate_setstate) :
  """
  Container for information extracted from the PDB header (if available).
  """
  __slots__ = ["d_min", "d_max", "r_work", "r_free", "rms_bonds", "rms_angles",
    "refinement_program"]
  def __init__ (self, pdb_file) :
    for name in self.__slots__ :
      setattr(self, name, None)
    if (pdb_file is not None) :
      import iotbx.pdb.hierarchy
      from iotbx.pdb import extract_rfactors_resolutions_sigma
      pdb_in = iotbx.pdb.hierarchy.input(file_name=pdb_file)
      published_results = extract_rfactors_resolutions_sigma.extract(
        file_lines=pdb_in.input.remark_section(), file_name=None)
      if (published_results is not None) :
        self.r_work = published_results.r_work
        self.r_free = published_results.r_free
        self.d_min = published_results.high
        self.d_max = published_results.low
      self.refinement_program = pdb_in.input.get_program_name()
      # XXX phenix.refine hack, won't work for other programs
      lines = open(pdb_file).readlines()
      for line in lines :
        if (line.startswith("REMARK Final:")) :
          fields = line.strip().split()
          self.rms_bonds = float(fields[-4])
          self.rms_angles = float(fields[-1])
          break

  def show (self, out=sys.stdout, prefix="", include_r_factors=True,
      include_rms_geom=True) :
    if (self.refinement_program is not None) :
      print >> out, "%sRefinement program    = %s" % (prefix,
        self.refinement_program)
    if (include_r_factors) :
      if (self.d_min is not None) :
        print >> out, "%sHigh resolution       = %6.2f" % (prefix, self.d_min)
      if (self.r_work is not None) :
        print >> out, "%sR-work                = %8.4f" % (prefix, self.r_work)
      if (self.r_free is not None) :
        print >> out, "%sR-free                = %8.4f" % (prefix, self.r_free)
    if (include_rms_geom) :
      if (self.rms_bonds is not None) :
        print >> out, "%sRMS(bonds)            = %8.4f" % (prefix,
          self.rms_bonds)
      if (self.rms_angles is not None) :
        print >> out, "%sRMS(angles)           = %6.2f" % (prefix,
          self.rms_angles)

class residue_multi_criterion (residue) :
  """
  Container for multiple outliers associated with a single residue.  If data
  are used, this may include real-space statistics regardless of whether the
  residue is technically an outlier or not.
  """
  __slots__ = residue.__slots__ + ["outliers", "n_confs", "i_seq"]
  def __init__ (self, **kwds) :
    residue.__init__(self, **kwds)
    self.outliers = []

  def add_outlier (self, outlier) :
    if isinstance(outlier, residue) :
      assert self.is_same_residue_group(outlier)
    self.outliers.append(outlier)

  def _find_outlier_type (self, outlier_type=None, outlier_types=(),
      retrieve_all=False) :
    assert (outlier_type is not None) or (len(outlier_types) > 0)
    for outlier in self.outliers :
      if (not outlier.is_outlier()) and (not retrieve_all) :
        continue
      otype = type(outlier).__name__
      if (otype == outlier_type) or (otype in outlier_types) :
        return True
    return False

  def is_ramachandran_outlier (self) :
    return self._find_outlier_type("ramachandran")

  def is_rotamer_outlier (self) :
    return self._find_outlier_type("rotamer")

  def is_cbeta_outlier (self) :
    return self._find_outlier_type("cbeta")

  def is_clash_outlier (self) :
    return self._find_outlier_type("clash")

  def is_geometry_outlier (self) :
    return self._find_outlier_type(
      outlier_types=["bond","angle","dihedral","chirality","planarity"])

  def __str__ (self) :
    outliers = []
    if self.is_ramachandran_outlier() : outliers.append("rama")
    if self.is_rotamer_outlier() : outliers.append("rota")
    if self.is_cbeta_outlier() : outliers.append("cb")
    if self.is_clash_outlier() : outliers.append("clash")
    if self.is_geometry_outlier() : outliers.append("geo")
    if (len(outliers) == 0) : outliers = ["---"]
    return "%s  %s" % (self.id_str(), ",".join(outliers))

  def __hash__ (self) :
    return self.residue_group_id_str().__hash__()

  def __cmp__ (self, other) :
    return cmp(self.i_seq, other.i_seq)

  def get_real_space_plot_values (self) :
    import numpy
    for outlier in self.outliers :
      if (type(outlier).__name__ == 'residue_real_space') :
        values = [ outlier.b_iso, outlier.cc, outlier.two_fofc, outlier.fmodel ]
        return values
    return [ numpy.NaN ] * 4

  def get_outlier_plot_values (self) :
    import numpy
    y = []
    if self.is_ramachandran_outlier() : y.append(1)
    else : y.append(numpy.NaN)
    if self.is_rotamer_outlier() : y.append(1)
    else : y.append(numpy.NaN)
    if self.is_cbeta_outlier() : y.append(1)
    else : y.append(numpy.NaN)
    if self.is_clash_outlier() : y.append(1)
    else : y.append(numpy.NaN)
    return y

class multi_criterion_view (slots_getstate_setstate) :
  """
  Container for generating multi-criterion plots and tables from separate lists
  of outliers.
  """
  __slots__ = ["residues", "is_populated"]
  def __init__ (self, pdb_hierarchy, include_all=False) :
    self.is_populated = False
    self.residues = {}
    i_seq = 0
    for chain in pdb_hierarchy.only_model().chains() :
      if (not include_all) :
        if (not chain.is_protein()) and (not chain.is_na()) :
          continue
      for residue_group in chain.residue_groups() :
        resname = residue_group.atom_groups()[0].resname
        if (resname == "HOH") : continue
        combined = residue_multi_criterion(
          chain_id=chain.id,
          resseq=residue_group.resseq,
          icode=residue_group.icode,
          resname=residue_group.atom_groups()[0].resname,
          altloc="",
          i_seq=i_seq,
          n_confs=len(residue_group.atom_groups()))
        id_str = combined.residue_group_id_str()
        self.residues[id_str] = combined
        i_seq += 1

  def process_outliers (self, outliers, log=sys.stderr) :
    self.is_populated = True
    for outlier in outliers :
      if outlier.is_single_residue_object() :
        if (outlier.resname == "HOH") : continue
        id_str = outlier.residue_group_id_str()
        if (id_str in self.residues) :
          self.residues[id_str].add_outlier(outlier)
        else :
          print >> log, "missing residue group '%s'" % id_str
      else :
        have_ids = set([])
        for atom in outlier.atoms_info :
          id_str = atom.residue_group_id_str()
          if (atom.resname == "HOH") or (id_str in have_ids) : continue
          if (id_str in self.residues) :
            self.residues[id_str].add_outlier(outlier)
            have_ids.add(id_str)
          else :
            print >> log, "missing residue group '%s'" % id_str

  def data (self) :
    return sorted(self.residues.values())

  def binned_data (self) :
    from mmtbx.validation import graphics
    return graphics.residue_binner(self.data())

  def get_y_limits (self) :
    import numpy
    values = []
    for outlier in self.data() :
      values.append(outlier.get_real_space_plot_values())
    values = numpy.array(values).transpose()
    rho_min = min(min(values[2]), min(values[3]))
    rho_max = max(max(values[2]), max(values[3]))
    return {
      "rho" : (rho_min, rho_max),
      "b" : (min(values[0]), max(values[0])),
      "cc" : (min(values[1]), max(values[1])),
    }

  def display_wx_plots (self) :
    import wxtbx.plots.molprobity
    frame = wxtbx.plots.molprobity.multi_criterion_frame(
      parent=None,
      title="MolProbity multi-criterion plot",
      validation=self)
    frame.Show()
