from __future__ import division
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import sys, math, mmtbx
from cctbx import geometry_restraints
from mmtbx.tls import tools
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

class adp(object):
  def __init__(self, model, wilson_b = None, n_histogram_slots = 10, file_name=None,
      selection=None):
    self.wilson_b = wilson_b
    self.file_name = file_name
    self.selection = selection
    self.rms_b_iso_or_b_equiv_bonded = model.rms_b_iso_or_b_equiv_bonded()
    eps = math.pi**2*8
    solvent_selection = model.solvent_selection()
    hd_selection = model.get_hd_selection()
    m_noH_sel = ((~solvent_selection) & (~hd_selection))
    s_noH_sel = ((solvent_selection) & (~hd_selection))
    #
    xs_a     = model.get_xray_structure()
    xs_a_noH = model.get_xray_structure().select(~hd_selection)
    xs_s_noH = model.get_xray_structure().select(s_noH_sel)
    xs_m_noH = model.get_xray_structure().select(m_noH_sel)
    xs_h     = model.get_xray_structure().select(hd_selection)
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
    uc = model.get_xray_structure().unit_cell()
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
    print >> out, pp+" all     : "+self._fmtl(self.n_iso_a, self.n_aniso_a, self.b_min_a,self.b_max_a,self.b_mean_a,self.a_min_a,self.a_max_a,self.a_mean_a,pad_r)
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
  """
  This class does not make any sence and should be removed.
  Please do not add anything here.
  Most likely this class should be merged into class info() below.
  """
  def __init__(self,
               model,
               wilson_b = None,
               use_molprobity=True,
               general_selection=None,
               ):
    self.model_manager = model
    self.geometry = model.geometry_statistics(
      general_selection=general_selection)
    self.adp = adp(model, wilson_b=wilson_b)
    self.tls_groups = model.tls_groups
    self.anomalous_scatterer_groups = model.anomalous_scatterer_groups

  def show(self, out=None, prefix="", padded=None, pdb_deposition=False):
    if(out is None): out = sys.stdout
    if(pdb_deposition): prefix="REMARK   3  "
    if(self.geometry is not None):
      self.geometry.show(log=out, prefix=prefix)
    print >> out, prefix
    self.adp.show(out=out, prefix=prefix, padded=padded,
      pdb_deposition=pdb_deposition)

    for info_pdb_str in [self.model_manager.tls_groups_as_pdb(),
        self.model_manager.anomalous_scatterer_groups_as_pdb(),
        self.model_manager.cartesian_NCS_as_pdb(),
        self.model_manager.torsion_NCS_as_pdb()]:
      if len(info_pdb_str) > 0:
        print >> out, prefix
        print >> out, info_pdb_str

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
    # adding NCS information.
    # It is not clear why we dump cartesian NCS first, and if it is absent,
    # Torsion NCS next. What about NCS constraints?
    if self.model_manager.cartesian_NCS_present():
      self.model_manager.cartesian_NCS_as_cif_block(cif_block=cif_block)
    elif self.model_manager.torsion_NCS_present():
      self.model_manager.torsion_NCS_as_cif_block(cif_block=cif_block)
    return cif_block

class info(object):
  def __init__(self, model,
                     fmodel_x          = None,
                     fmodel_n          = None,
                     refinement_params = None,
                     general_selection = None,
                     use_molprobity    = True):
    ref_par = refinement_params
    wilson_b = None
    if fmodel_x is not None:
      wilson_b = fmodel_x.wilson_b()
    elif fmodel_n is not None:
      wilson_b = fmodel_n.wilson_b()
    self.model = mmtbx.model_statistics.model(
      model = model,
      wilson_b = wilson_b,
      use_molprobity = use_molprobity,
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
