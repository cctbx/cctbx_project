from __future__ import division
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import sys, math, mmtbx
from cctbx import geometry_restraints
from mmtbx.tls import tools
from libtbx.str_utils import line_breaker
from libtbx.str_utils import format_value
from itertools import count
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.validation.cbetadev import cbetadev
from mmtbx.validation.clashscore import clashscore
from mmtbx.validation.utils import molprobity_score
from mmtbx.validation import omegalyze
from mmtbx.validation import cablam
import iotbx.cif.model
from libtbx import group_args

class statistics(object):
  def __init__(self,
               pdb_hierarchy,
               geometry_restraints_manager=None):
    self.pdb_hierarchy = pdb_hierarchy
    self.geometry_restraints_manager = geometry_restraints_manager
    self.from_restraints = None
    self.update()

  def update(self, pdb_hierarchy=None):
    if(pdb_hierarchy is not None):
      self.pdb_hierarchy = pdb_hierarchy
    if(self.geometry_restraints_manager is not None):
      sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
      self.from_restraints = \
        self.geometry_restraints_manager.energies_sites(
          sites_cart        = sites_cart,
          compute_gradients = False)
      assert approx_equal(
        self.from_restraints.target,
        self.from_restraints.angle_residual_sum+
        self.from_restraints.bond_residual_sum+
        self.from_restraints.chirality_residual_sum+
        self.from_restraints.dihedral_residual_sum+
        self.from_restraints.nonbonded_residual_sum+
        self.from_restraints.planarity_residual_sum+
        self.from_restraints.parallelity_residual_sum+
        self.from_restraints.reference_coordinate_residual_sum+
        self.from_restraints.reference_dihedral_residual_sum+
        self.from_restraints.ncs_dihedral_residual_sum+
        self.from_restraints.den_residual_sum+
        self.from_restraints.ramachandran_residual_sum)

  def angle(self):
    mi,ma,me,n = 0,0,0,0
    if(self.from_restraints is not None):
      mi,ma,me = self.from_restraints.angle_deviations()
      n = self.from_restraints.get_filtered_n_angle_proxies()
    return group_args(min = mi, max = ma, mean = me, n = n)

  def bond(self):
    mi,ma,me,n = 0,0,0,0
    if(self.from_restraints is not None):
      mi,ma,me = self.from_restraints.bond_deviations()
      n = self.from_restraints.get_filtered_n_bond_proxies()
    return group_args(min = mi, max = ma, mean = me, n = n)

  def chirality(self):
    mi,ma,me,n = 0,0,0,0
    if(self.from_restraints is not None):
      mi,ma,me = self.from_restraints.chirality_deviations()
      n = self.from_restraints.n_chirality_proxies
    return group_args(min = mi, max = ma, mean = me, n = n)

  def dihedral(self):
    mi,ma,me,n = 0,0,0,0
    if(self.from_restraints is not None):
      mi,ma,me = self.from_restraints.dihedral_deviations()
      n = self.from_restraints.n_dihedral_proxies
    return group_args(min = mi, max = ma, mean = me, n = n)

  def planarity(self):
    mi,ma,me,n = 0,0,0,0
    if(self.from_restraints is not None):
      mi,ma,me = self.from_restraints.planarity_deviations()
      n = self.from_restraints.n_planarity_proxies
    return group_args(min = mi, max = ma, mean = me, n = n)

  def parallelity(self):
    mi,ma,me,n = 0,0,0,0
    if(self.from_restraints is not None):
      mi,ma,me = self.from_restraints.parallelity_deviations()
      n = self.from_restraints.n_parallelity_proxies
    return group_args(min = mi, max = ma, mean = me, n = n)

  def nonbonded(self):
    mi,ma,me,n = 0,0,0,0
    if(self.from_restraints is not None):
      mi,ma,me = self.from_restraints.nonbonded_deviations()
      n = self.from_restraints.n_nonbonded_proxies
    return group_args(min = mi, max = ma, mean = me, n = n)

  def ramachandran(self):
    result = ramalyze(pdb_hierarchy = self.pdb_hierarchy, outliers_only = False)
    return group_args(
      outliers = result.percent_outliers,
      allowed  = result.percent_allowed,
      favored  = result.percent_favored)

  def rotamer(self):
    result = rotalyze(pdb_hierarchy = self.pdb_hierarchy, outliers_only = False)
    return group_args(outliers = result.percent_outliers)

  def c_beta(self):
    result = cbetadev(pdb_hierarchy = self.pdb_hierarchy,
       outliers_only = True, out = null_out()) # XXX Why it is different from others?
    return group_args(outliers = result.get_weighted_outlier_percent())

  def clash(self):
    result = clashscore(pdb_hierarchy = self.pdb_hierarchy)
    return group_args(
      score = result.get_clashscore())

  def cablam(self):
    result = cablam.cablamalyze(self.pdb_hierarchy, outliers_only=False,
      out=null_out(), quiet=True) # XXX Why it is different from others?
    return group_args(
      outliers    = result.percent_outliers(),
      disfavored  = result.percent_disfavored(),
      ca_outliers = result.percent_ca_outliers())

  def omega(self):
    result = omegalyze.omegalyze(pdb_hierarchy=self.pdb_hierarchy, quiet=True) # XXX
    # XXX Move this to omegalyze function.
    n_proline         = result.n_proline()
    n_general         = result.n_general()
    n_cis_proline     = result.n_cis_proline()
    n_cis_general     = result.n_cis_general()
    n_twisted_proline = result.n_twisted_proline()
    n_twisted_general = result.n_twisted_general()
    cis_general     = 0
    twisted_general = 0
    cis_proline     = 0
    twisted_proline = 0
    if(n_proline != 0):
      cis_proline     = n_cis_proline    *100./n_proline
      twisted_proline = n_twisted_proline*100./n_proline
    if(n_general != 0):
      cis_general     = n_cis_general    *100./n_general
      twisted_general = n_twisted_general*100./n_general
    return group_args(
      cis_proline     = cis_proline,
      cis_general     = cis_general,
      twisted_general = twisted_general,
      twisted_proline = twisted_proline)

  def show(self, log=None, prefix="", lowercase=False):
    if(log is None): log = sys.stdout
    def fmt(f1,f2,d1):
      fmt_str= "%6.3f %7.3f %6d"
      if f1 is None  : return '   -       -       -  '
      return fmt_str%(f1,f2,d1)
    def fmt2(f1):
      if f1 is None: return '  -   '
      return "%-6.3f"%(f1)
    a,b,c,d,p,n = self.angle(), self.bond(), self.chirality(), self.dihedral(), \
      self.planarity(), self.nonbonded()
    result = "%s" % prefix
    result += """
%sDEVIATIONS FROM IDEAL VALUES.
%s  BOND      : %s
%s  ANGLE     : %s
%s  CHIRALITY : %s
%s  PLANARITY : %s
%s  DIHEDRAL  : %s
%s  MIN NONBONDED DISTANCE : %s
%s"""%(prefix,
       prefix, fmt(b.mean, b.max, b.n),
       prefix, fmt(a.mean, a.max, a.n),
       prefix, fmt(c.mean, c.max, c.n),
       prefix, fmt(p.mean, p.max, p.n),
       prefix, fmt(d.mean, d.max, d.n),
       prefix, fmt2(n.min),
       prefix)
    result += "%s" % prefix
    result += """
%sMOLPROBITY STATISTICS.
%s  ALL-ATOM CLASHSCORE : %s
%s  RAMACHANDRAN PLOT:
%s    OUTLIERS : %-5.2f %s
%s    ALLOWED  : %-5.2f %s
%s    FAVORED  : %-5.2f %s
%s  ROTAMER OUTLIERS : %s %s
%s  CBETA DEVIATIONS : %-d
%s  PEPTIDE PLANE:
%s    CIS-PROLINE     : %s
%s    CIS-GENERAL     : %s
%s    TWISTED PROLINE : %s
%s    TWISTED GENERAL : %s"""%(
        prefix,
        prefix, format_value("%-6.2f", self.clash().score).strip(),
        prefix,
        prefix, self.ramachandran().outliers, "%",
        prefix, self.ramachandran().allowed, "%",
        prefix, self.ramachandran().favored, "%",
        prefix, str("%6.2f"%(self.rotamer().outliers)).strip(),"%",
        prefix, self.c_beta().outliers,
        prefix,
        prefix, str(self.omega().cis_proline),
        prefix, str(self.omega().cis_general),
        prefix, str(self.omega().twisted_general),
        prefix, str(self.omega().twisted_proline))
    if(lowercase):
      result = result.swapcase()
    print >> log, result

  def as_cif_block(self, cif_block=None):
    if cif_block is None:
      cif_block = iotbx.cif.model.block()
    cif_block["_refine.pdbx_stereochemistry_target_values"] = "GeoStd + Monomer Library"
    loop = iotbx.cif.model.loop(header=(
      "_refine_ls_restr.type",
      "_refine_ls_restr.number",
      "_refine_ls_restr.dev_ideal",
      #"_refine_ls_restr.dev_ideal_target",
      "_refine_ls_restr.weight",
      #"_refine_ls_restr.pdbx_refine_id",
      "_refine_ls_restr.pdbx_restraint_function",
    ))
    a,b,c,d,p,n = self.angle(), self.bond(), self.chirality(), self.dihedral(), \
      self.planarity(), self.nonbonded()
    loop.add_row(("f_bond_d",           b.n, b.mean, "?", "?"))
    loop.add_row(("f_angle_d",          a.n, a.mean, "?", "?"))
    loop.add_row(("f_chiral_restr",     c.n, c.mean, "?", "?"))
    loop.add_row(("f_plane_restr",      p.n, p.mean, "?", "?"))
    loop.add_row(("f_dihedral_angle_d", d.n, d.mean, "?", "?"))
    cif_block.add_loop(loop)
    return cif_block

# XXX TO GO NEXT
class geometry(statistics):
  def __init__(
        self,
        pdb_hierarchy,
        restraints_manager,
        molprobity_scores=False,
        n_histogram_slots=10,
        cdl_restraints=False,
        ignore_hydrogens=False,  # XXX only used by amber
        automatically_use_amber=True,
        ):

    super(geometry, self).__init__(
        pdb_hierarchy=pdb_hierarchy,
        geometry_restraints_manager=restraints_manager)

  def as_cif_block(self, cif_block=None):
    if cif_block is None:
      cif_block = iotbx.cif.model.block()
    cif_block["_refine.pdbx_stereochemistry_target_values"] = "GeoStd + Monomer Library"
    if self.cdl_restraints:
      cif_block["_refine.pdbx_stereochemistry_target_values"] += " + CDL v1.2"
    loop = iotbx.cif.model.loop(header=(
      "_refine_ls_restr.type",
      "_refine_ls_restr.number",
      "_refine_ls_restr.dev_ideal",
      #"_refine_ls_restr.dev_ideal_target",
      "_refine_ls_restr.weight",
      #"_refine_ls_restr.pdbx_refine_id",
      "_refine_ls_restr.pdbx_restraint_function",
    ))
    a,b,c,d,p,n = self.angle(), self.bond(), self.chirality(), self.dihedral(), \
      self.planarity(), self.nonbonded()
    loop.add_row(("f_bond_d",           self.b.n, self.b.mean, "?", "?"))
    loop.add_row(("f_angle_d",          self.a.n, self.a.mean, "?", "?"))
    loop.add_row(("f_chiral_restr",     self.c.n, self.c.mean, "?", "?"))
    loop.add_row(("f_plane_restr",      self.p.n, self.p.mean, "?", "?"))
    loop.add_row(("f_dihedral_angle_d", self.d.n, self.d.mean, "?", "?"))
    cif_block.add_loop(loop)
    return cif_block

class model_content(object):
  def __init__(self, model):
    self.atoms_count = model.xray_structure.scatterers().size()
    self.atoms_occupancy_sum = \
      flex.sum(model.xray_structure.scatterers().extract_occupancies())
    self.scattering_types_counts_and_occupancy_sums = \
      model.xray_structure.scattering_types_counts_and_occupancy_sums()

  def show(self, out = None, prefix = "", pdb_deposition = False):
    if(out is None): out = sys.stdout
    if(pdb_deposition):
      prefix = "REMARK   3  "
    fmt = "   %5s               %10d        %8.2f"
    print >> out, prefix+"MODEL CONTENT."
    print >> out, prefix+" ELEMENT        ATOM RECORD COUNT   OCCUPANCY SUM"
    for item in self.scattering_types_counts_and_occupancy_sums:
      print >> out, prefix+fmt % (item.scattering_type, item.count,
        item.occupancy_sum)
    print >> out,prefix+fmt%("TOTAL",self.atoms_count,self.atoms_occupancy_sum)

class adp(object):
  def __init__(self, model,  n_histogram_slots = 10, file_name=None,
      selection=None):
    self.wilson_b = model.wilson_b
    self.file_name = file_name
    self.selection = selection
    self.rms_b_iso_or_b_equiv_bonded = model.rms_b_iso_or_b_equiv_bonded()
    eps = math.pi**2*8
    solvent_selection = model.solvent_selection()
    hd_selection = model.xray_structure.hd_selection()
    m_noH_sel = ((~solvent_selection) & (~hd_selection))
    s_noH_sel = ((solvent_selection) & (~hd_selection))
    #
    xs_a     = model.xray_structure
    xs_a_noH = model.xray_structure.select(~hd_selection)
    xs_s_noH = model.xray_structure.select(s_noH_sel)
    xs_m_noH = model.xray_structure.select(m_noH_sel)
    xs_h     = model.xray_structure.select(hd_selection)
    #
    u_a     = xs_a    .extract_u_iso_or_u_equiv()
    u_a_noH = xs_a_noH.extract_u_iso_or_u_equiv()
    u_s_noH = xs_s_noH.extract_u_iso_or_u_equiv()
    u_m_noH = xs_m_noH.extract_u_iso_or_u_equiv()
    u_h     = xs_h    .extract_u_iso_or_u_equiv()
    self.b_min_a,    self.b_max_a,    self.b_mean_a    = self.mmmd(u_a,    eps)
    self.b_min_a_noH,self.b_max_a_noH,self.b_mean_a_noH= self.mmmd(u_a_noH,eps)
    self.b_min_s_noH,self.b_max_s_noH,self.b_mean_s_noH= self.mmmd(u_s_noH,eps)
    self.b_min_m_noH,self.b_max_m_noH,self.b_mean_m_noH= self.mmmd(u_m_noH,eps)
    self.b_min_h,    self.b_max_h,    self.b_mean_h    = self.mmmd(u_h,    eps)
    #
    uc = model.xray_structure.unit_cell()
    a_a     = xs_a    .scatterers().anisotropy(unit_cell =uc).select(xs_a    .use_u_aniso())
    a_a_noH = xs_a_noH.scatterers().anisotropy(unit_cell =uc).select(xs_a_noH.use_u_aniso())
    a_s_noH = xs_s_noH.scatterers().anisotropy(unit_cell =uc).select(xs_s_noH.use_u_aniso())
    a_m_noH = xs_m_noH.scatterers().anisotropy(unit_cell =uc).select(xs_m_noH.use_u_aniso())
    a_h     = xs_h    .scatterers().anisotropy(unit_cell =uc).select(xs_h    .use_u_aniso())
    #
    self.n_aniso_a     = xs_a    .use_u_aniso().count(True)
    self.n_aniso_a_noH = xs_a_noH.use_u_aniso().count(True)
    self.n_aniso_s_noH = xs_s_noH.use_u_aniso().count(True)
    self.n_aniso_m_noH = xs_m_noH.use_u_aniso().count(True)
    self.n_aniso_h     = xs_h    .use_u_aniso().count(True)
    self.n_iso_a       = xs_a    .use_u_iso().count(True)
    self.n_iso_a_noH   = xs_a_noH.use_u_iso().count(True)
    self.n_iso_s_noH   = xs_s_noH.use_u_iso().count(True)
    self.n_iso_m_noH   = xs_m_noH.use_u_iso().count(True)
    self.n_iso_h       = xs_h    .use_u_iso().count(True)
    #
    self.a_min_a,    self.a_max_a,    self.a_mean_a    = self.mmmd(a_a)
    self.a_min_a_noH,self.a_max_a_noH,self.a_mean_a_noH= self.mmmd(a_a_noH)
    self.a_min_s_noH,self.a_max_s_noH,self.a_mean_s_noH= self.mmmd(a_s_noH)
    self.a_min_m_noH,self.a_max_m_noH,self.a_mean_m_noH= self.mmmd(a_m_noH)
    self.a_min_h,    self.a_max_h,    self.a_mean_h    = self.mmmd(a_h)
    #
    self.b_a_noH_histogram = flex.histogram(data = u_a_noH * eps,
      n_slots = n_histogram_slots)
    self.b_a_noH = u_a_noH * eps # need this for phenix gui
    self.a_a_noH_histogram = flex.histogram(data = a_a_noH,
      n_slots = n_histogram_slots)
    #
    self._show_anisotropy = (xs_a.use_u_aniso()).count(True)

  def mmmd(self, values, eps = None):
    if(eps is not None): values = values * eps
    return flex.min_default(values, None), \
           flex.max_default(values, None), \
           flex.mean_default(values, None)

  def _fmtl(self,a,b,c,d,e,f,g,h, pad):
    n = []
    for item in [a,b,c,d,e,f,g,h]:
      if(item is None): item = str(item)
      elif(str(item).count(".")): item = str("%8.2f"%item).strip()
      else: item = str("%6d"%item).strip()
      n.append(item)
    fmt = "%-6s %-6s %-8s%-8s%-8s %-6s%-7s%-7s"+pad+"|"
    return fmt%(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7])

  def _histogram(self, histogram, out, pad_l, pad_r):
    result = []
    # n_slots must be even
    low_cutoff_1 = histogram.data_min()
    for (i_1,n_1) in enumerate(histogram.slots()):
      high_cutoff_1 = histogram.data_min() + histogram.slot_width()*(i_1+1)
      result.append([i_1,low_cutoff_1,high_cutoff_1,n_1])
      low_cutoff_1 = high_cutoff_1
    result1, result2 = result[:len(result)//2], result[len(result)//2:]
    fmt = "|"+pad_l+"   %d:%10.3f -%8.3f:%5d   |   %d:%10.3f -%8.3f:%5d    "+pad_r+"|"
    for r1, r2 in zip(result1, result2):
      print >> out, fmt%(r1[0],r1[1],r1[2],r1[3], r2[0],r2[1],r2[2],r2[3])
    print >> out, "|"+pad_l+\
     "                            =>continue=>                              "+\
     pad_r+"|"

  def show(self, out=None, prefix="", padded=None, pdb_deposition=False):
    if(pdb_deposition):
      self.show_2(out = out, prefix="REMARK   3  ")
    else:
      self.show_1(out = out, padded = padded, prefix = prefix)

  # XXX for GUI - see wxtbx.adp_statistics
  def format_tables (self) :
    def fs (value) :
      if (isinstance(value, int)) :
        return format_value("%d", value, replace_none_with="---")
      elif (isinstance(value, float)) :
        return format_value("%.2f", value, replace_none_with="---")
      else :
        return "---"
    table1 = []
    row1 = ["all"]
    for value in [self.n_iso_a,self.n_aniso_a,self.b_min_a,self.b_max_a,
                  self.b_mean_a,self.a_min_a,self.a_max_a,self.a_mean_a] :
      row1.append(fs(value))
    row2 = ["all(noH)"]
    for value in [self.n_iso_a_noH,self.n_aniso_a_noH,self.b_min_a_noH,
                  self.b_max_a_noH,self.b_mean_a_noH,self.a_min_a_noH,
                  self.a_max_a_noH,self.a_mean_a_noH] :
      row2.append(fs(value))
    row3 = ["solvent"]
    for value in [self.n_iso_s_noH,self.n_aniso_s_noH,self.b_min_s_noH,
                  self.b_max_s_noH,self.b_mean_s_noH,self.a_min_s_noH,
                  self.a_max_s_noH,self.a_mean_s_noH] :
      row3.append(fs(value))
    row4 = ["macro."]
    for value in [self.n_iso_m_noH,self.n_aniso_m_noH,self.b_min_m_noH,
                  self.b_max_m_noH,self.b_mean_m_noH,self.a_min_m_noH,
                  self.a_max_m_noH,self.a_mean_m_noH] :
      row4.append(fs(value))
    row5 = ["hydrogen"]
    for value in [self.n_iso_h,self.n_aniso_h,self.b_min_h,self.b_max_h,
                  self.b_mean_h,self.a_min_h,self.a_max_h,self.a_mean_h] :
      row5.append(fs(value))
    table1 = [row1,row2,row3,row4,row5]
    def get_bins (histogram) :
      result = []
      low_cutoff_1 = histogram.data_min()
      for (i_1,n_1) in enumerate(histogram.slots()):
        high_cutoff_1 = histogram.data_min() + histogram.slot_width()*(i_1+1)
        result.append([i_1, "%.3f - %.3f" % (low_cutoff_1,high_cutoff_1), n_1])
        low_cutoff_1 = high_cutoff_1
      return result
    table2 = get_bins(self.b_a_noH_histogram)
    table3 = get_bins(self.a_a_noH_histogram)
    return (table1,table2,table3)

  def format_plots (self) :
    def get_values (histogram) :
      y = histogram.slots()
      y_range = histogram.data_min(), histogram.data_max()
      return y, y_range
    y1 = get_values(self.b_a_noH_histogram)
    y2 = get_values(self.a_a_noH_histogram)
    return y1, y2

  def show_1(self, out = None, padded = None, prefix = ""):
    if(out is None): out = sys.stdout
    pad_l, pad_r = "", ""
    if(padded):
       pad_l = " "*3
       pad_r = " "*4
    #
    if(self.wilson_b is not None):
       bw = str("%9.3f"%self.wilson_b).strip()
       s1 = "|-ADP statistics (Wilson B = "+bw+")"
    else:
       s1 = "|-ADP statistics"
    if(padded): s1 = s1 + "-"*( 79-len(s1+"|") ) + "|"
    else: s1 = s1 + "-"*( 72-len(s1+"|") ) + "|"
    print >> out, prefix+s1
    #
    pp = prefix+"|"+pad_l
    print >> out, pp+" Atom    | Number of   | Isotropic or equivalent| Anisotropy lmin/max "+pad_r+"|"
    print >> out, pp+" type    |iso    aniso | min     max     mean   | min   max    mean   "+pad_r+"|"
    print >> out, pp+" - - - - |- - - - - - -| - - - - - - - - - - - -| - - - - - - - - - - "+pad_r+"|"
    #
    print >> out, pp+" all     : "+self._fmtl(self.n_iso_a,self.n_aniso_a,self.b_min_a,self.b_max_a,self.b_mean_a,self.a_min_a,self.a_max_a,self.a_mean_a,pad_r)
    print >> out, pp+" all(noH): "+self._fmtl(self.n_iso_a_noH,self.n_aniso_a_noH,self.b_min_a_noH,self.b_max_a_noH,self.b_mean_a_noH,self.a_min_a_noH,self.a_max_a_noH,self.a_mean_a_noH,pad_r)
    print >> out, pp+" Sol.    : "+self._fmtl(self.n_iso_s_noH,self.n_aniso_s_noH,self.b_min_s_noH,self.b_max_s_noH,self.b_mean_s_noH,self.a_min_s_noH,self.a_max_s_noH,self.a_mean_s_noH,pad_r)
    print >> out, pp+" Mac.    : "+self._fmtl(self.n_iso_m_noH,self.n_aniso_m_noH,self.b_min_m_noH,self.b_max_m_noH,self.b_mean_m_noH,self.a_min_m_noH,self.a_max_m_noH,self.a_mean_m_noH,pad_r)
    print >> out, pp+" Hyd.    : "+self._fmtl(self.n_iso_h,self.n_aniso_h,self.b_min_h,self.b_max_h,self.b_mean_h,self.a_min_h,self.a_max_h,self.a_mean_h,pad_r)
    print >> out, pp+" "+"- "*34+pad_r+" |"
    #
    name = "|"+pad_l+"    Distribution of isotropic (or equivalent) ADP for non-H atoms:    "+pad_r+"|"
    print >> out, prefix+name
    print >> out, prefix+"|"+pad_l+" Bin#      value range     #atoms | Bin#      value range     #atoms  "+pad_r+"|"
    self._histogram(histogram = self.b_a_noH_histogram, out = out, pad_l=pad_l,pad_r=pad_r)
    #
    if(self._show_anisotropy):
       print >> out, prefix+"|"+pad_l+" "+"- "*34+pad_r+" |"
       name = "|"+pad_l+"                     Distribution of anisotropy:                      "+pad_r+"|"
       print >> out, prefix+name
       print >> out, prefix+"|"+pad_l+" Bin#      value range     #atoms | Bin#      value range     #atoms  "+pad_r+"|"
       self._histogram(histogram = self.a_a_noH_histogram, out = out, pad_l=pad_l,pad_r=pad_r)
    #
    if(padded): print >> out, prefix+"|"+"-"*77+"|"
    else: print >> out, prefix+"|"+"-"*70+"|"
    out.flush()

  def _fmtl2(self,a,b):
    n = []
    for item in [a,b]:
      if(item is None): item = str(item)
      elif(str(item).count(".")): item = str("%8.2f"%item).strip()
      else: item = str("%7d"%item).strip()
      n.append(item)
    fmt = " %6s %7s"
    return fmt%(n[0],n[1])

  def show_2(self, out, prefix):
    if(self.wilson_b is None):
      wilson_b = str(self.wilson_b)
    elif(str(self.wilson_b).count(".")):
      wilson_b = str("%8.2f"%self.wilson_b).strip()
    print >> out, prefix+"ATOMIC DISPLACEMENT PARAMETERS."
    print >> out, prefix+" WILSON B : %-s"%wilson_b
    print >> out, prefix+" RMS(B_ISO_OR_EQUIVALENT_BONDED) : %-s"%format_value(
      "%8.2f", self.rms_b_iso_or_b_equiv_bonded).strip()
    print >> out, prefix+" ATOMS          NUMBER OF ATOMS"
    print >> out, prefix+"                  ISO.  ANISO. "
    print >> out, prefix+"  ALL         :"+self._fmtl2(self.n_iso_a    ,self.n_aniso_a    )
    print >> out, prefix+"  ALL (NO H)  :"+self._fmtl2(self.n_iso_a_noH,self.n_aniso_a_noH)
    print >> out, prefix+"  SOLVENT     :"+self._fmtl2(self.n_iso_s_noH,self.n_aniso_s_noH)
    print >> out, prefix+"  NON-SOLVENT :"+self._fmtl2(self.n_iso_m_noH,self.n_aniso_m_noH)
    print >> out, prefix+"  HYDROGENS   :"+self._fmtl2(self.n_iso_h    ,self.n_aniso_h    )
    out.flush()

  def as_cif_block(self, cif_block=None):
    if cif_block is None:
      cif_block = iotbx.cif.model.block()
    cif_block["_reflns.B_iso_Wilson_estimate"] = self.wilson_b
    cif_block["_refine.B_iso_mean"] = self.b_mean_a
    #_refine.aniso_B[1][1]
    #_refine.aniso_B[2][2]
    #_refine.aniso_B[3][3]
    #_refine.aniso_B[1][2]
    #_refine.aniso_B[1][3]
    #_refine.aniso_B[2][3]
    return cif_block

class model(object):
  def __init__(self,
               model,
               ignore_hd,
               use_molprobity=True,
               ncs_manager=None,
               cdl_restraints=False,
               general_selection=None,
               ):
    self.geometry = model.geometry_statistics(
      ignore_hd = ignore_hd,
      molprobity_scores=use_molprobity,
      cdl_restraints=cdl_restraints,
      general_selection=general_selection,
      )
    self.content = model_content(model)
    self.adp = adp(model)
    self.tls_groups = model.tls_groups
    self.anomalous_scatterer_groups = model.anomalous_scatterer_groups
    self.ncs_groups = model.extract_ncs_groups()
    self.ncs_manager = ncs_manager
    self.pdb_hierarchy = model.pdb_hierarchy(sync_with_xray_structure=True)
    self.cdl_restraints = cdl_restraints

  def show(self, out=None, prefix="", padded=None, pdb_deposition=False):
    if(out is None): out = sys.stdout
    if(pdb_deposition): prefix="REMARK   3  "
    if(self.geometry is not None):
      self.geometry.show(log=out, prefix=prefix)
    print >> out, prefix
    self.adp.show(out=out, prefix=prefix, padded=padded,
      pdb_deposition=pdb_deposition)
    if(self.tls_groups is not None):
      print >> out, prefix
      tools.remark_3_tls(tlsos             = self.tls_groups.tlsos,
                         selection_strings = self.tls_groups.selection_strings,
                         out               = out)
    if(self.anomalous_scatterer_groups is not None):
      print >> out, prefix
      self.show_anomalous_scatterer_groups(out = out)
    if(self.ncs_groups is not None):
      print >> out, prefix
      self.show_ncs_groups(out = out)
    if(self.ncs_manager is not None):
      print >> out, prefix
      self.show_torsion_ncs_groups(out = out)

  def show_ncs_groups(self, out = None):
    if(out is None): out = sys.stdout
    pr = "REMARK   3  "
    print >>out,pr+"NCS DETAILS."
    print >>out,pr+" NUMBER OF NCS GROUPS : %-6d"%len(self.ncs_groups)
    for i_group, ncs_group in enumerate(self.ncs_groups):
      print >>out,pr+" NCS GROUP : %-6d"%(i_group+1)
      selection_strings = ncs_group.group.selection_strings
      for i_op,pair,mx,rms in zip(
          count(1),
          ncs_group.group.selection_pairs,
          ncs_group.matrices,
          ncs_group.rms):
        print >> out,pr+"  NCS OPERATOR : %-d" % i_op
        lines = line_breaker(selection_strings[0], width=34)
        for i_line, line in enumerate(lines):
          if(i_line == 0):
            print >> out, pr+"   REFERENCE SELECTION: %s"%line
          else:
            print >> out, pr+"                      : %s"%line
        lines = line_breaker(selection_strings[i_op], width=34)
        for i_line, line in enumerate(lines):
          if(i_line == 0):
            print >> out, pr+"   SELECTION          : %s"%line
          else:
            print >> out, pr+"                      : %s"%line
        print >> out,pr+"   ATOM PAIRS NUMBER  : %-d" % len(pair[0])
        print >> out,pr+"   RMSD               : %-10.3f" % rms

  def ncs_as_cif_block(self, cif_block=None):
    if cif_block is None:
      cif_block = iotbx.cif.model.block()

    ncs_ens_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_ens.id",
      "_struct_ncs_ens.details"))
    ncs_dom_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_dom.id",
      "_struct_ncs_dom.pdbx_ens_id",
      "_struct_ncs_dom.details"))
    ncs_dom_lim_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_dom_lim.pdbx_ens_id",
      "_struct_ncs_dom_lim.dom_id",
      #"_struct_ncs_dom_lim.pdbx_component_id",
      #"_struct_ncs_dom_lim.pdbx_refine_code",
      "_struct_ncs_dom_lim.beg_auth_asym_id",
      "_struct_ncs_dom_lim.beg_auth_seq_id",
      "_struct_ncs_dom_lim.end_auth_asym_id",
      "_struct_ncs_dom_lim.end_auth_seq_id",
      "_struct_ncs_dom_lim.selection_details"))

    ncs_oper_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_oper.id",
      "_struct_ncs_oper.code",
      "_struct_ncs_oper.matrix[1][1]",
      "_struct_ncs_oper.matrix[1][2]",
      "_struct_ncs_oper.matrix[1][3]",
      "_struct_ncs_oper.matrix[2][1]",
      "_struct_ncs_oper.matrix[2][2]",
      "_struct_ncs_oper.matrix[2][3]",
      "_struct_ncs_oper.matrix[3][1]",
      "_struct_ncs_oper.matrix[3][2]",
      "_struct_ncs_oper.matrix[3][3]",
      "_struct_ncs_oper.vector[1]",
      "_struct_ncs_oper.vector[2]",
      "_struct_ncs_oper.vector[3]",
      "_struct_ncs_oper.details"))

    ncs_ens_gen_loop = iotbx.cif.model.loop(header=(
      "_struct_ncs_ens_gen.dom_id_1",
      "_struct_ncs_ens_gen.dom_id_2",
      "_struct_ncs_ens_gen.ens_id",
      "_struct_ncs_ens_gen.oper_id"))

    oper_id = 0
    if self.ncs_groups is not None:
      for i_group, ncs_group in enumerate(self.ncs_groups):
        ncs_ens_loop.add_row((i_group+1, "?"))
        selection_strings = ncs_group.group.selection_strings
        matrices = ncs_group.matrices
        rms = ncs_group.rms
        pair_count = len(ncs_group.group.selection_pairs[0])
        for i_domain, domain_selection in enumerate(selection_strings):
          ncs_dom_loop.add_row((i_domain+1, i_group+1, "?"))
          # XXX TODO: export individual sequence ranges from selection
          ncs_dom_lim_loop.add_row(
            (i_group+1, i_domain+1, "?", "?", "?", "?", domain_selection))
          if i_domain > 0:
            rt_mx = ncs_group.matrices[i_domain-1]
            oper_id += 1
            row = [oper_id, "given"]
            row.extend(rt_mx.r)
            row.extend(rt_mx.t)
            row.append("?")
            ncs_oper_loop.add_row(row)
            ncs_ens_gen_loop.add_row((1, i_domain+1, i_group+1, oper_id))
    elif self.ncs_manager is not None:
      for i_group, ncs_group in enumerate(self.ncs_manager.ncs_groups):
        ncs_ens_loop.add_row((i_group+1, "?"))
        for i_domain, ncs_domain in enumerate(ncs_group):
          ncs_dom_loop.add_row((i_domain+1, i_group+1, "?"))
          segments = self.ncs_manager.master_ranges[ncs_domain]
          asym_id = ncs_domain.split("and")[0].split()[1].strip()
          for i_segment, segment_range in enumerate(segments):
            ncs_dom_lim_loop.add_row(
              (i_group+1, i_domain+1, asym_id, segment_range[0],
               asym_id, segment_range[1], ncs_domain))
    cif_block.add_loop(ncs_ens_loop)
    cif_block.add_loop(ncs_dom_loop)
    cif_block.add_loop(ncs_dom_lim_loop)
    if self.ncs_groups is not None:
      cif_block.add_loop(ncs_oper_loop)
      cif_block.add_loop(ncs_ens_gen_loop)
    return cif_block

  def show_torsion_ncs_groups(self, out = None):
    if(out is None): out = sys.stdout
    restraint_groups = self.ncs_manager.ncs_groups
    torsion_counts=self.ncs_manager.get_number_of_restraints_per_group(
      pdb_hierarchy=self.pdb_hierarchy)
    sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    self.ncs_manager.get_torsion_rmsd(sites_cart=sites_cart)
    pr = "REMARK   3  "
    print >>out,pr+"TORSION NCS DETAILS."
    print >>out,pr+" NUMBER OF NCS GROUPS : %-6d"%len(restraint_groups)
    for i_group, ncs_group in enumerate(restraint_groups):
      count = 0
      print >>out,pr+" NCS GROUP : %-6d"%(i_group+1)
      selection_strings = ncs_group
      for selection in selection_strings:
        lines = line_breaker(selection, width=34)
        for i_line, line in enumerate(lines):
          if (i_line == 0):
            print >> out, pr+"   SELECTION          : %s"%line
          else:
            print >> out, pr+"                      : %s"%line
        count += torsion_counts[selection]
      print >> out,pr+"   RESTRAINED TORSIONS: %-d" % count
      if self.ncs_manager.torsion_rmsd is not None:
        print >> out,pr+"   BELOW LIMIT RMSD   : %-10.3f" % \
          self.ncs_manager.torsion_rmsd
      if self.ncs_manager.all_torsion_rmsd is not None:
        print >> out,pr+"   ALL RESTRAINT RMSD : %-10.3f" % \
          self.ncs_manager.all_torsion_rmsd
    if self.ncs_manager.histogram_under_limit is not None:
      print >> out, pr + "  Histogram of differences under limit:"
      self.ncs_manager.histogram_under_limit.show(
        f=out,
        prefix=pr+"  ",
        format_cutoffs="%8.3f")
    if self.ncs_manager.histogram_over_limit is not None:
      print >> out, pr + "  Histogram of differences over limit:"
      self.ncs_manager.histogram_over_limit.show(
        f=out,
        prefix=pr+"  ",
        format_cutoffs="%8.3f")

  def show_anomalous_scatterer_groups(self, out = None):
    if(out is None): out = sys.stdout
    pr = "REMARK   3  "
    print >>out,pr+"ANOMALOUS SCATTERER GROUPS DETAILS."
    print >>out,pr+" NUMBER OF ANOMALOUS SCATTERER GROUPS : %-6d"%\
      len(self.anomalous_scatterer_groups)
    counter = 0
    for group in self.anomalous_scatterer_groups:
      counter += 1
      print >>out,pr+" ANOMALOUS SCATTERER GROUP : %-6d"%counter
      lines = line_breaker(group.selection_string, width=45)
      for i_line, line in enumerate(lines):
        if(i_line == 0):
          print >> out, pr+"  SELECTION: %s"%line
        else:
          print >> out, pr+"           : %s"%line
      print >>out,pr+"  fp  : %-15.4f"%group.f_prime
      print >>out,pr+"  fdp : %-15.4f"%group.f_double_prime
      out.flush()

  def as_cif_block(self, cif_block=None):
    import iotbx.pdb.mmcif
    if self.geometry is not None:
      cif_block = self.geometry.as_cif_block(cif_block=cif_block)
    cif_block = self.adp.as_cif_block(cif_block=cif_block)
    if self.tls_groups is not None:
      cif_block = iotbx.pdb.mmcif.tls_as_cif_block(
        tlsos=self.tls_groups.tlsos,
        selection_strings=self.tls_groups.selection_strings,
        cif_block=cif_block)
    if self.anomalous_scatterer_groups is not None:
      pass
      #self.show_anomalous_scatterer_groups(out = out)
    if self.ncs_groups is not None or self.ncs_manager is not None:
      self.ncs_as_cif_block(cif_block=cif_block)
    return cif_block

class info(object):
  def __init__(self, model,
                     fmodel_x          = None,
                     fmodel_n          = None,
                     refinement_params = None,
                     ignore_hd         = True,
                     general_selection = None,
                     use_molprobity    = True):
    ref_par = refinement_params
    ncs_manager = None
    if model != None:
      if model.restraints_manager != None:
        if model.restraints_manager.geometry != None:
          ncs_manager = model.restraints_manager.geometry.ncs_dihedral_manager
    self.model = mmtbx.model_statistics.model(
      model = model,
      ignore_hd = ignore_hd,
      use_molprobity = use_molprobity,
      ncs_manager = ncs_manager,
      cdl_restraints = ref_par.pdb_interpretation.restraints_library.cdl,
      general_selection = general_selection,
      )
    self.data_x, self.data_n = None, None
    if(fmodel_x is not None):
      self.data_x = fmodel_x.info(
        free_reflections_per_bin = ref_par.alpha_beta.free_reflections_per_bin,
        max_number_of_bins       = ref_par.main.max_number_of_resolution_bins)
    if(fmodel_n is not None):
      self.data_n = fmodel_n.info(
        free_reflections_per_bin = ref_par.alpha_beta.free_reflections_per_bin,
        max_number_of_bins       = ref_par.main.max_number_of_resolution_bins)

  def show_remark_3(self, out = None):
    prefix = "REMARK   3  "
    if(out is None): out = sys.stdout
    if(self.data_n is not None):
      print >> out, prefix+"X-RAY DATA."
      print >> out, prefix
    self.data_x.show_remark_3(out = out)
    print >> out, prefix
    if(self.data_n is not None):
      print >> out, prefix+"NEUTRON DATA."
      print >> out, prefix
      self.data_n.show_remark_3(out = out)
      print >> out, prefix
    self.model.show(out = out, pdb_deposition = True)

  def as_cif_block(self, cif_block=None):
    cif_block = self.data_x.as_cif_block(cif_block=cif_block)
    # XXX Neutron data?
    cif_block = self.model.as_cif_block(cif_block=cif_block)
    return cif_block
