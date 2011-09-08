from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import sys, math, mmtbx
from cctbx import geometry_restraints
from mmtbx.tls import tools
from libtbx.str_utils import line_breaker
from libtbx.str_utils import format_value
from itertools import count

class geometry(object):
  def __init__(self,
               sites_cart,
               pdb_hierarchy,
               restraints_manager,
               hd_selection,
               ignore_hd = True,
               main_chain_selection = None,
               ignore_side_chain = False,
               n_histogram_slots = 10,
               molprobity_scores=False):
    zero = [0, 0, 0, 0, 0]
    self.a_target, self.a_mean, self.a_max, self.a_min, self.a_number = zero
    self.b_target, self.b_mean, self.b_max, self.b_min, self.b_number = zero
    self.c_target, self.c_mean, self.c_max, self.c_min, self.c_number = zero
    self.d_target, self.d_mean, self.d_max, self.d_min, self.d_number = zero
    self.p_target, self.p_mean, self.p_max, self.p_min, self.p_number = zero
    self.n_target, self.n_mean, self.n_max, self.n_min, self.n_number = zero
    self.dr_target, self.dr_number = 0, 0
    self.target = 0.0
    self.number_of_restraints = 0.0
    if(ignore_side_chain): # and main_chain_selection.count(True) > 0):
      restraints_manager = restraints_manager.select(
        selection = main_chain_selection)
      sites_cart = sites_cart.select(main_chain_selection)
    elif(ignore_hd and hd_selection.count(True) > 0):
      restraints_manager = restraints_manager.select(selection = ~hd_selection)
      sites_cart = sites_cart.select(~hd_selection)
    energies_sites = \
      restraints_manager.energies_sites(sites_cart        = sites_cart,
                                        compute_gradients = True)
    esg = energies_sites.geometry
    b_deviations = esg.bond_deviations()
    n_deviations = esg.nonbonded_deviations()
    a_deviations = esg.angle_deviations()
    d_deviations = esg.dihedral_deviations()
    c_deviations = esg.chirality_deviations()
    p_deviations = esg.planarity_deviations()
    self.a_min, self.a_max, self.a_mean = a_deviations
    self.b_min, self.b_max, self.b_mean = b_deviations
    self.c_min, self.c_max, self.c_mean = c_deviations
    self.d_min, self.d_max, self.d_mean = d_deviations
    self.p_min, self.p_max, self.p_mean = p_deviations
    self.n_min, self.n_max, self.n_mean = n_deviations
    self.a_target,self.a_number=esg.angle_residual_sum,    esg.n_angle_proxies
    self.b_target,self.b_number=esg.bond_residual_sum,     esg.n_bond_proxies
    self.c_target,self.c_number=esg.chirality_residual_sum,esg.n_chirality_proxies
    self.d_target,self.d_number=esg.dihedral_residual_sum, esg.n_dihedral_proxies
    self.p_target,self.p_number=esg.planarity_residual_sum,esg.n_planarity_proxies
    self.n_target,self.n_number=esg.nonbonded_residual_sum,esg.n_bond_proxies
    self.rd_target, self.rd_number= \
      esg.reference_dihedral_residual_sum, esg.n_reference_dihedral_proxies
    self.nd_target, self.nd_number= \
      esg.ncs_dihedral_residual_sum, esg.n_ncs_dihedral_proxies
    self.ra_target, self.ra_number = \
      esg.generic_restraint_residual_sum, esg.n_generic_proxies
    self.target = esg.target # XXX normalization ?
    assert approx_equal(self.target, self.a_target+self.b_target+self.c_target+
      self.d_target+self.p_target+self.n_target+self.rd_target+
      self.nd_target+self.ra_target)
    self.norm_of_gradients = esg.gradients.norm()
    self.number_of_restraints = esg.number_of_restraints # XXX normalization ?
    #
    if(self.number_of_restraints > 0):
      self.target = self.target / self.number_of_restraints
    if(self.a_number > 0): self.a_target = self.a_target/self.a_number
    if(self.b_number > 0): self.b_target = self.b_target/self.b_number
    if(self.c_number > 0): self.c_target = self.c_target/self.c_number
    if(self.d_number > 0): self.d_target = self.d_target/self.d_number
    if(self.p_number > 0): self.p_target = self.p_target/self.p_number
    if(self.n_number > 0): self.n_target = self.n_target/self.n_number
    #
    rmg = restraints_manager.geometry
    bond_deltas = geometry_restraints.bond_deltas(
                          sites_cart         = sites_cart,
                          sorted_asu_proxies = rmg.pair_proxies().bond_proxies)
    angle_deltas = geometry_restraints.angle_deltas(
                                                sites_cart = sites_cart,
                                                proxies    = rmg.angle_proxies)
    nonbonded_distances = esg.nonbonded_distances()
    self.bond_deltas_histogram = \
      flex.histogram(data = flex.abs(bond_deltas), n_slots = n_histogram_slots)
    self.angle_deltas_histogram = \
      flex.histogram(data = flex.abs(angle_deltas), n_slots =n_histogram_slots)
    self.nonbonded_distances_histogram = flex.histogram(
      data = flex.abs(nonbonded_distances), n_slots = n_histogram_slots)
    # molprobity scores
    if(molprobity_scores):
      from mmtbx.validation.ramalyze import ramalyze
      from mmtbx.validation.rotalyze import rotalyze
      from mmtbx.validation.cbetadev import cbetadev
      from mmtbx.validation.clashscore import clashscore
      self.ramalyze_obj = ramalyze()
      self.ramalyze_obj.analyze_pdb(hierarchy = pdb_hierarchy,
        outliers_only = False)
      self.rotalyze_obj = rotalyze()
      self.rotalyze_obj.analyze_pdb(hierarchy = pdb_hierarchy,
        outliers_only = False)[1]
      self.cbetadev_obj = cbetadev()
      self.cbetadev_obj.analyze_pdb(hierarchy = pdb_hierarchy,
        outliers_only = False)
      self.clashscore_obj = clashscore()
      self.clashscore_obj.analyze_clashes(hierarchy = pdb_hierarchy)
      self.clashscore = self.clashscore_obj.get_clashscore()
      self.ramachandran_outliers = self.ramalyze_obj.get_outliers_count_and_fraction()[1]*100.
      self.ramachandran_allowed  = self.ramalyze_obj.get_allowed_count_and_fraction()[1]*100.
      self.ramachandran_favored  = self.ramalyze_obj.get_favored_count_and_fraction()[1]*100.
      self.rotamer_outliers = self.rotalyze_obj.get_outliers_count_and_fraction()[1]*100.
      self.c_beta_dev = self.cbetadev_obj.get_outlier_count()

  def show(self, out=None, prefix="", pdb_deposition=False, message = ""):
    if(out is None): out = sys.stdout
    if(pdb_deposition):
      prefix = "REMARK   3  "
      self.show_overall_2(out = out, prefix = prefix)
    else:
      self.show_bond_angle_nonbonded_histogram_1(out = out, message = message)
      self.show_overall_1(out = out, message = message)

  def show_bond_angle_nonbonded_histogram_1(self, out = None, message = ""):
    if(out is None): out = sys.stdout
    if(len(message) > 0): message_ = "|-Geometry statistics: "+message
    else: message_ = "|-Geometry statistics"
    line_len = len(message_+"|")
    fill_len = 80-line_len-1
    print >> out, message_+"-"*(fill_len)+"|"
    print >> out, "|            Histogram of deviations from ideal values for                    |"
    print >> out, "| Bonds                | Angles                   |  Nonbonded contacts       |"
    fmt = "| %5.3f - %5.3f: %5d | %7.3f - %7.3f: %5d | %6.3f - %6.3f: %8d |"
    self._bond_angle_nonbonded_histogram(out=out, prefix="", format=fmt)
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def show_overall_1(self, out = None, message = ""):
    if(out is None): out = sys.stdout
    fmt = "| %s | %s | %s %s %s | %s | %s |"
    if(len(message) > 0): message_ = "|-Geometry statistics: "+message
    else: message_ = "|-Geometry statistics"
    line_len = len(message_+"|")
    fill_len = 80-line_len-1
    s0 = message_+"-"*(fill_len)+"|"
    s1 = "|      Type |    Count |    Deviation from ideal |     Targets | Target (sum) |"
    s2 = "|           |          |    rmsd      max    min |             |              |"
    s3 = fmt%(format_value("%9s","bond"),
              format_value("%8d",self.b_number),
              format_value("%7.3f",self.b_mean),
              format_value("%8.3f",self.b_max),
              format_value("%6.3f",self.b_min),
              format_value("%11.3f",self.b_target),
              format_value("%12s"," "))
    s4 = fmt%(format_value("%9s","angle"),
              format_value("%8d",self.a_number),
              format_value("%7.3f",self.a_mean),
              format_value("%8.3f",self.a_max),
              format_value("%6.3f",self.a_min),
              format_value("%11.3f",self.a_target),
              format_value("%12s"," "))
    s5 = fmt%(format_value("%9s","chirality"),
              format_value("%8d",self.c_number),
              format_value("%7.3f",self.c_mean),
              format_value("%8.3f",self.c_max),
              format_value("%6.3f",self.c_min),
              format_value("%11.3f",self.c_target),
              format_value("%12.3f",self.target))
    s6 = fmt%(format_value("%9s","planarity"),
              format_value("%8d",self.p_number),
              format_value("%7.3f",self.p_mean),
              format_value("%8.3f",self.p_max),
              format_value("%6.3f",self.p_min),
              format_value("%11.3f",self.p_target),
              format_value("%12s"," "))
    s7 = fmt%(format_value("%9s","dihedral"),
              format_value("%8d",self.d_number),
              format_value("%7.3f",self.d_mean),
              format_value("%8.3f",self.d_max),
              format_value("%6.3f",self.d_min),
              format_value("%11.3f",self.d_target),
              format_value("%12s"," "))
    s8 = fmt%(format_value("%9s","nonbonded"),
              format_value("%8d",self.n_number),
              format_value("%7.3f",self.n_mean),
              format_value("%8.3f",self.n_max),
              format_value("%6.3f",self.n_min),
              format_value("%11.3f",self.n_target),
              format_value("%12s"," "))
    s9 = "|"+"-"*77+"|"
    for line in [s0,s1,s2,s3,s4,s5,s6,s7,s8,s9]:
      print >> out, line
    out.flush()

  def show_bond_angle_nonbonded_histogram_2(self, out = None, prefix = ""):
    if(out is None): out = sys.stdout
    print >> out, prefix+"DEVIATIONS FROM IDEAL VALUES (HISTOGRAM).       "
    print >> out, prefix+" BONDS            ANGLES               NONBONDED"
    fmt = " %-4.2f-%-4.2f %5d  %-6.2f-%-6.2f %5d  %-4.2f-%-4.2f %5d"
    self._bond_angle_nonbonded_histogram(out=out, prefix=prefix, format=fmt)
    out.flush()

  def show_overall_2(self, out = None, prefix = ""):
    if(out is None): out = sys.stdout
    fmt = "%6.3f %7.3f %6d"
    pr = prefix
    print >> out, pr+"DEVIATIONS FROM IDEAL VALUES."
    print >> out, pr+"               RMSD     MAX  COUNT"
    print >> out, pr+" BOND      : "+fmt%(self.b_mean,self.b_max,self.b_number)
    print >> out, pr+" ANGLE     : "+fmt%(self.a_mean,self.a_max,self.a_number)
    print >> out, pr+" CHIRALITY : "+fmt%(self.c_mean,self.c_max,self.c_number)
    print >> out, pr+" PLANARITY : "+fmt%(self.p_mean,self.p_max,self.p_number)
    print >> out, pr+" DIHEDRAL  : "+fmt%(self.d_mean,self.d_max,self.d_number)
    print >> out, pr+" MIN NONBONDED DISTANCE : %-6.3f"%self.n_min
    print >> out, pr
    self.show_molprobity_scores(out = out, prefix = prefix)
    out.flush()

  def show_molprobity_scores(self, out = None, prefix = ""):
    pr = prefix
    print >> out, pr+"MOLPROBITY STATISTICS."
    print >> out, pr+" ALL-ATOM CLASHSCORE : %-6.2f"%self.clashscore
    print >> out, pr+" RAMACHANDRAN PLOT:"
    print >> out, pr+"   OUTLIERS : %-5.2f %s"%(self.ramachandran_outliers,"%")
    print >> out, pr+"   ALLOWED  : %-5.2f %s"%(self.ramachandran_allowed,"%")
    print >> out, pr+"   FAVORED  : %-5.2f %s"%(self.ramachandran_favored,"%")
    print >> out, pr+" ROTAMER OUTLIERS : %s %s"%(
      str("%6.2f"%(self.rotamer_outliers)).strip(),"%")
    print >> out, pr+" CBETA DEVIATIONS : %-d"%self.c_beta_dev
    out.flush()

  def _bond_angle_nonbonded_histogram(self, out, prefix, format):
    h_1 = self.bond_deltas_histogram
    h_2 = self.angle_deltas_histogram
    h_3 = self.nonbonded_distances_histogram
    lc_1 = h_1.data_min()
    lc_2 = h_2.data_min()
    lc_3 = h_3.data_min()
    s_1 = enumerate(h_1.slots())
    s_2 = enumerate(h_2.slots())
    s_3 = enumerate(h_3.slots())
    for (i_1,n_1),(i_2,n_2),(i_3,n_3) in zip(s_1, s_2, s_3):
      hc_1 = h_1.data_min() + h_1.slot_width() * (i_1+1)
      hc_2 = h_2.data_min() + h_2.slot_width() * (i_2+1)
      hc_3 = h_3.data_min() + h_3.slot_width() * (i_3+1)
      print >> out, prefix+format%(lc_1,hc_1,n_1,lc_2,hc_2,n_2,lc_3,hc_3,n_3)
      lc_1 = hc_1
      lc_2 = hc_2
      lc_3 = hc_3
    out.flush()

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
  def __init__(self, model,  n_histogram_slots = 10):
    self.wilson_b = model.wilson_b
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

class model(object):
  def __init__(self, model, ignore_hd):
    self.geometry = model.geometry_statistics(ignore_hd = ignore_hd,
      molprobity_scores=True)
    self.content = model_content(model)
    self.adp = adp(model)
    self.tls_groups = model.tls_groups
    self.anomalous_scatterer_groups = model.anomalous_scatterer_groups
    self.ncs_groups = model.extract_ncs_groups()

  def show(self, out=None, prefix="", padded=None, pdb_deposition=False):
    if(out is None): out = sys.stdout
    if(pdb_deposition): prefix="REMARK   3  "
    if(self.geometry is not None):
      self.geometry.show(out=out, prefix=prefix, pdb_deposition=pdb_deposition)
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

class info(object):
  def __init__(self, model,
                     fmodel_x          = None,
                     fmodel_n          = None,
                     refinement_params = None,
                     ignore_hd         = True):
    ref_par = refinement_params
    self.model = mmtbx.model_statistics.model(model = model,
      ignore_hd = ignore_hd)
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
