from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import sys, math, mmtbx
from cctbx import geometry_restraints
from mmtbx import bulk_solvent
from libtbx import adopt_init_args
from mmtbx.tls import tools
from libtbx.itertbx import count
from libtbx.str_utils import line_breaker


class geometry(object):
  def __init__(self,
               sites_cart,
               restraints_manager,
               n_histogram_slots = 10):
    zero = [0, 0, 0, 0, 0]
    self.a_target, self.a_mean, self.a_max, self.a_min, self.a_number = zero
    self.b_target, self.b_mean, self.b_max, self.b_min, self.b_number = zero
    self.c_target, self.c_mean, self.c_max, self.c_min, self.c_number = zero
    self.d_target, self.d_mean, self.d_max, self.d_min, self.d_number = zero
    self.p_target, self.p_mean, self.p_max, self.p_min, self.p_number = zero
    self.n_target, self.n_mean, self.n_max, self.n_min, self.n_number = zero
    self.target = 0.0
    self.number_of_restraints = 0.0
    energies_sites = \
      restraints_manager.energies_sites(sites_cart        = sites_cart,
                                        compute_gradients = False)
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
    self.target = esg.target # XXX normalization ?
    assert approx_equal(self.target, self.a_target+self.b_target+self.c_target+
      self.d_target+self.p_target+self.n_target)
    self.number_of_restraints = esg.number_of_restraints # XXX normalization ?
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

  def show(self, out=None, prefix="", pdb_deposition=False):
    if(out is None): out = sys.stdout
    if(pdb_deposition):
      prefix = "REMARK   3  "
      self.show_overall_2(out = out, prefix = prefix)
      print >> out, prefix
      self.show_bond_angle_nonbonded_histogram_2(out = out, prefix = prefix)
    else:
      self.show_bond_angle_nonbonded_histogram_1(out = out, prefix = prefix)
      self.show_overall_1(out = out, prefix = prefix)

  def show_bond_angle_nonbonded_histogram_1(self, out = None, prefix = ""):
    if(out is None): out = sys.stdout
    print >> out, \
        prefix+"|------------------------------------------------------------|"
    print >> out, \
        prefix+"|        Histogram of deviations from ideal values for       |"
    print >> out, \
        prefix+"|Bonds             |Angles                |Nonbonded contacts|"
    fmt = "|%5.3f-%5.3f:%6d|%7.3f-%7.3f:%6d|%5.3f-%5.3f:%6d|"
    self._bond_angle_nonbonded_histogram(out=out, prefix=prefix, format=fmt)
    print >> out, \
        prefix+"|------------------------------------------------------------|"
    out.flush()

  def show_overall_1(self, out = None, prefix = ""):
    if(out is None): out = sys.stdout
    fmt1 = "|%7.3f%8.3f%7.3f|%12.3f|              |"
    fmt2 = "|%7.3f%8.3f%7.3f|%12.3f|%12.3f  |"
    print >> out, prefix+"|-Geometry statistics----------------------------------------|"
    print >> out, prefix+"|Type     | Deviation from ideal |   Targets  |Target (sum)  |"
    print >> out, prefix+"|         |  rmsd     max    min |            |              |"
    print >> out, prefix+"|bond     "+fmt1%(self.b_mean,self.b_max,self.b_min,self.b_target)
    print >> out, prefix+"|angle    "+fmt1%(self.a_mean,self.a_max,self.a_min,self.a_target)
    print >> out, prefix+"|chirality"+fmt2%(self.c_mean,self.c_max,self.c_min,self.c_target,self.target)
    print >> out, prefix+"|planarity"+fmt1%(self.p_mean,self.p_max,self.p_min,self.p_target)
    print >> out, prefix+"|dihedral "+fmt1%(self.d_mean,self.d_max,self.d_min,self.d_target)
    print >> out, prefix+"|nonbonded"+fmt1%(self.n_mean,self.n_max,self.n_min,self.n_target)
    print >> out, prefix+"|------------------------------------------------------------|"
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
    self.element_types_and_counts = self._get_element_types_and_counts(
      xray_structure = model.xray_structure)

  def _get_element_types_and_counts(self, xray_structure):
    reg = xray_structure.scattering_type_registry()
    result = []
    scattering_types = xray_structure.scatterers().extract_scattering_types()
    occupancies = xray_structure.scatterers().extract_occupancies()
    for scattering_type in reg.type_index_pairs_as_dict().keys():
      selection = scattering_types == scattering_type
      count = selection.count(True)
      occupancy_sum = flex.sum(occupancies.select(selection))
      result_ = (scattering_type, count, occupancy_sum)
      result.append(result_)
    return result

  def show(self, out = None, prefix = "", pdb_deposition = False):
    if(out is None): out = sys.stdout
    if(pdb_deposition):
      prefix = "REMARK   3  "
    fmt = "   %5s               %10d        %8.2f"
    print >> out, prefix+"MODEL CONTENT."
    print >> out, prefix+" ELEMENT        ATOM RECORD COUNT   OCCUPANCY SUM"
    for item in self.element_types_and_counts:
      print >> out, prefix+fmt % (item[0], item[1], item[2])
    print >> out,prefix+fmt%("TOTAL",self.atoms_count,self.atoms_occupancy_sum)

class adp(object):
  def __init__(self, model,  n_histogram_slots = 10):
    self.wilson_b = model.wilson_b
    eps = math.pi**2*8
    solvent_selection = model.solvent_selection
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
    u_a     = xs_a    .extract_u_iso_or_u_equiv().select(xs_a    .use_u_iso())
    u_a_noH = xs_a_noH.extract_u_iso_or_u_equiv().select(xs_a_noH.use_u_iso())
    u_s_noH = xs_s_noH.extract_u_iso_or_u_equiv().select(xs_s_noH.use_u_iso())
    u_m_noH = xs_m_noH.extract_u_iso_or_u_equiv().select(xs_m_noH.use_u_iso())
    u_h     = xs_h    .extract_u_iso_or_u_equiv().select(xs_h    .use_u_iso())
    self.b_min_a,    self.b_max_a,    self.b_mean_a    = self.mmmd(u_a,    eps)
    self.b_min_a_noH,self.b_max_a_noH,self.b_mean_a_noH= self.mmmd(u_a_noH,eps)
    self.b_min_s_noH,self.b_max_s_noH,self.b_mean_s_noH= self.mmmd(u_s_noH,eps)
    self.b_min_m_noH,self.b_max_m_noH,self.b_mean_m_noH= self.mmmd(u_m_noH,eps)
    self.b_min_h,    self.b_max_h,    self.b_mean_h    = self.mmmd(u_h,    eps)
    #
    uc = model.xray_structure.unit_cell()
    a_a     = xs_a    .scatterers().anisotropy(
                                  unit_cell =uc).select(xs_a    .use_u_aniso())
    a_a_noH = xs_a_noH.scatterers().anisotropy(
                                  unit_cell =uc).select(xs_a_noH.use_u_aniso())
    a_s_noH = xs_s_noH.scatterers().anisotropy(
                                  unit_cell =uc).select(xs_s_noH.use_u_aniso())
    a_m_noH = xs_m_noH.scatterers().anisotropy(
                                  unit_cell =uc).select(xs_m_noH.use_u_aniso())
    a_h     = xs_h    .scatterers().anisotropy(
                                  unit_cell =uc).select(xs_h    .use_u_aniso())
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
    result1, result2 = result[:len(result)/2], result[len(result)/2:]
    fmt = "|"+pad_l+"   %d:%10.3f -%8.3f:%5d   |   %d:%10.3f -%8.3f:%5d    "+pad_r+"|"
    for r1, r2 in zip(result1, result2):
      print >> out, self.prefix+fmt%(r1[0],r1[1],r1[2],r1[3], r2[0],r2[1],r2[2],r2[3])
    print >> out, self.prefix+"|"+pad_l+\
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
    self._histogram(values = self.b_a_noH_histogram, out = out, pad_l=pad_l,pad_r=pad_r)
    #
    if(self._show_anisotropy):
       print >> out, prefix+"|"+pad_l+" "+"- "*34+pad_r+" |"
       name = "|"+pad_l+"                     Distribution of anisotropy:                      "+pad_r+"|"
       print >> out, self.prefix+name
       print >> out, self.prefix+"|"+pad_l+" Bin#      value range     #atoms | Bin#      value range     #atoms  "+pad_r+"|"
       self._histogram(values = self.a_a_noH_histogram, out = out, pad_l=pad_l,pad_r=pad_r)
    #
    if(padded): print >> out, prefix+"|"+"-"*77+"|"
    else: print >> out, prefix+"|"+"-"*70+"|"
    out.flush()

  def _fmtl2(self,a,b,c,d,e):
    n = []
    for item in [a,b,c,d,e]:
      if(item is None): item = str(item)
      elif(str(item).count(".")): item = str("%8.2f"%item).strip()
      else: item = str("%7d"%item).strip()
      n.append(item)
    fmt = " %6s %7s  %7s %7s %7s"
    return fmt%(n[0],n[1],n[2],n[3],n[4])

  def show_2(self, out, prefix):
    if(self.wilson_b is None):
      wilson_b = str(self.wilson_b)
    elif(str(self.wilson_b).count(".")):
      wilson_b = str("%8.2f"%self.wilson_b).strip()
    print >> out, prefix+"ATOMIC DISPLACEMENT PARAMETERS."
    print >> out, prefix+" WILSON B : %-s"%wilson_b
    print >> out, prefix+" ATOMS          NUMBER OF ATOMS  ISOTROPIC OR EQUIVALENT"
    print >> out, prefix+"                  ISO.  ANISO.      MIN     MAX    MEAN"
    print >> out, prefix+"  ALL         :"+self._fmtl2(self.n_iso_a,self.n_aniso_a,self.b_min_a,self.b_max_a,self.b_mean_a)
    print >> out, prefix+"  ALL (NO H)  :"+self._fmtl2(self.n_iso_a_noH,self.n_aniso_a_noH,self.b_min_a_noH,self.b_max_a_noH,self.b_mean_a_noH)
    print >> out, prefix+"  SOLVENT     :"+self._fmtl2(self.n_iso_s_noH,self.n_aniso_s_noH,self.b_min_s_noH,self.b_max_s_noH,self.b_mean_s_noH)
    print >> out, prefix+"  NON-SOLVENT :"+self._fmtl2(self.n_iso_m_noH,self.n_aniso_m_noH,self.b_min_m_noH,self.b_max_m_noH,self.b_mean_m_noH)
    print >> out, prefix+"  HYDROGENS   :"+self._fmtl2(self.n_iso_h,self.n_aniso_h,self.b_min_h,self.b_max_h,self.b_mean_h)
    print >> out, prefix
    print >> out, prefix+"ATOMIC DISPLACEMENT PARAMETERS (HISTOGRAM, NON-H)."
    if(self._show_anisotropy):
      print >> out, prefix+" ISOTROPIC OR EQUIVALENT             ANISOTROPY"
    else:
      print >> out, prefix+" ISOTROPIC OR EQUIVALENT"
    h_1 = self.b_a_noH_histogram
    h_2 = self.a_a_noH_histogram
    lc_1 = h_1.data_min()
    lc_2 = h_2.data_min()
    s_1 = enumerate(h_1.slots())
    s_2 = enumerate(h_2.slots())
    for (i_1,n_1),(i_2,n_2) in zip(s_1, s_2):
      hc_1 = h_1.data_min() + h_1.slot_width() * (i_1+1)
      hc_2 = h_2.data_min() + h_2.slot_width() * (i_2+1)
      if(self._show_anisotropy):
        fmt = "    %6.2f-%-6.2f %6d     %5.2f-%-5.2f %6d"
        print >> out, prefix+fmt%(lc_1,hc_1,n_1, lc_2,hc_2,n_2)
      else:
        fmt = "    %6.2f-%-6.2f %6d"
        print >> out, prefix+fmt%(lc_1,hc_1,n_1)
      lc_1 = hc_1
      lc_2 = hc_2
    out.flush()

class model(object):
  def __init__(self, model):
    self.geometry = geometry(sites_cart = model.xray_structure.sites_cart(),
                             restraints_manager = model.restraints_manager)
    self.content = model_content(model)
    self.adp = adp(model)
    self.tls_groups = model.tls_groups
    self.anomalous_scatterer_groups = model.anomalous_scatterer_groups
    self.ncs_groups = model.extract_ncs_groups()

  def show(self, out=None, prefix="", padded=None, pdb_deposition=False):
    if(out is None): out = sys.stdout
    if(pdb_deposition): prefix="REMARK   3  "
    self.content.show(out=out, prefix=prefix, pdb_deposition=pdb_deposition)
    print >> out, prefix
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
        #print >> out,pr+"   REFERENCE SELECTION: %-s" % selection_strings[0]
        #print >> out,pr+"   SELECTION          : %-s" % selection_strings[i_op]
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
      #print >>out,pr
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

class resolution_bin(object):
  def __init__(self,
               i_bin        = None,
               d_range      = None,
               completeness = None,
               alpha_work   = None,
               beta_work    = None,
               r_work       = None,
               r_free       = None,
               target_work  = None,
               target_free  = None,
               n_work       = None,
               n_free       = None,
               mean_f_obs   = None,
               fom          = None,
               pher_work    = None,
               pher_free    = None):
    adopt_init_args(self, locals())

class data(object):
  def __init__(self, fmodel, free_reflections_per_bin, max_number_of_bins,
               sigma_fobs_rejection_criterion, sigma_iobs_rejection_criterion):
    mp = fmodel.mask_params
    self.twin_fraction = None # XXX
    self.twin_law = None # XXX
    self.r_work = fmodel.r_work()
    self.r_free = fmodel.r_free()
    self.r_all = fmodel.r_all()
    self.overall_scale_k1 = fmodel.scale_k1()
    self.number_of_test_reflections = fmodel.f_obs_t.data().size()
    self.number_of_work_reflections = fmodel.f_obs_w.data().size()
    self.number_of_reflections = fmodel.f_obs.data().size()
    self.k_sol = fmodel.k_sol()
    self.b_sol = fmodel.b_sol()
    self.b_cart = fmodel.b_cart()
    self.mask_solvent_radius = mp.solvent_radius
    self.mask_shrink_radius = mp.shrink_truncation_radius
    self.mask_grid_step_factor = mp.grid_step_factor
    self.ml_phase_error = flex.mean(fmodel.phase_errors())
    self.ml_coordinate_error = fmodel.model_error_ml()
    self.d_max, self.d_min = fmodel.f_obs.resolution_range()
    self.completeness_in_range = fmodel.f_obs.completeness(d_max = self.d_max)
    self.sigma_fobs_rejection_criterion = sigma_fobs_rejection_criterion
    self.sigma_iobs_rejection_criterion = sigma_iobs_rejection_criterion
    self.target_name = fmodel.target_name
    self.sf_algorithm = fmodel.sfg_params.algorithm
    self.bins = self.statistics_in_resolution_bins(
      fmodel = fmodel,
      free_reflections_per_bin = free_reflections_per_bin,
      max_number_of_bins = max_number_of_bins)

  def statistics_in_resolution_bins(self, fmodel, free_reflections_per_bin,
                                    max_number_of_bins):
    result = []
    fo_t = fmodel.f_obs_t
    fc_t = fmodel.f_model_t()
    fo_w = fmodel.f_obs_w
    fc_w = fmodel.f_model_w()
    alpha_w, beta_w = fmodel.alpha_beta_w()
    alpha_t, beta_t = fmodel.alpha_beta_t()
    pher_w = fmodel.phase_errors_work()
    pher_t = fmodel.phase_errors_test()
    fom = fmodel.figures_of_merit()
    fmodel.f_obs.setup_binner(n_bins=fmodel.determine_n_bins(
      free_reflections_per_bin=free_reflections_per_bin,
      max_n_bins=max_number_of_bins))
    fo_t.use_binning_of(fmodel.f_obs)
    fc_t.use_binning_of(fo_t)
    fo_w.use_binning_of(fo_t)
    fc_w.use_binning_of(fo_t)
    alpha_w.use_binning_of(fo_t)
    alpha_t.use_binning_of(fo_t)
    beta_w.use_binning_of(fo_t)
    beta_t.use_binning_of(fo_t)
    target_functor = fmodel.target_functor()
    target_result = target_functor(compute_gradients=False)
    tpr = target_result.target_per_reflection()
    if(tpr.size() != 0):
      tpr_w = tpr.select(fmodel.work)
      tpr_t = tpr.select(fmodel.test)
    for i_bin in fo_t.binner().range_used():
      sel_t = fo_t.binner().selection(i_bin)
      sel_w = fo_w.binner().selection(i_bin)
      sel_all = fmodel.f_obs.binner().selection(i_bin)
      sel_fo_all = fmodel.f_obs.select(sel_all)
      sel_fo_t = fo_t.select(sel_t)
      sel_fc_t = fc_t.select(sel_t)
      sel_fo_w = fo_w.select(sel_w)
      sel_fc_w = fc_w.select(sel_w)
      if (tpr.size() == 0):
        sel_tpr_w = None
        sel_tpr_t = None
      else:
        denom_w = flex.sum(sel_fo_w.data()*sel_fo_w.data())
        denom_t = flex.sum(sel_fo_t.data()*sel_fo_t.data())
        if(denom_w != 0):
           sel_tpr_w = flex.sum(tpr_w.select(sel_w))/denom_w
        else:
           sel_tpr_w = flex.sum(tpr_w.select(sel_w))
        if(denom_t != 0):
           sel_tpr_t = flex.sum(tpr_t.select(sel_t))/denom_t
        else:
           sel_tpr_t = flex.sum(tpr_t.select(sel_t))
      d_max_,d_min_ = sel_fo_all.d_max_min()
      completeness = fmodel.f_obs.resolution_filter(
        d_min= d_min_,d_max= d_max_).completeness(d_max = d_max_)
      d_range = fo_t.binner().bin_legend(
        i_bin=i_bin, show_bin_number=False, show_counts=False)
      bin = resolution_bin(
        i_bin        = i_bin,
        d_range      = d_range,
        completeness = completeness,
        alpha_work   = flex.mean(alpha_w.select(sel_w).data()),
        beta_work    = flex.mean(beta_w.select(sel_w).data()),
        r_work       = bulk_solvent.r_factor(sel_fo_w.data(), sel_fc_w.data()),
        r_free       = bulk_solvent.r_factor(sel_fo_t.data(), sel_fc_t.data()),
        target_work  = sel_tpr_w,
        target_free  = sel_tpr_t,
        n_work       = sel_fo_w.data().size(),
        n_free       = sel_fo_t.data().size(),
        mean_f_obs   = flex.mean(sel_fo_all.data()),
        fom          = flex.mean(fom.select(sel_all)),
        pher_work    = flex.mean(pher_w.select(sel_w)),
        pher_free    = flex.mean(pher_t.select(sel_t)))
      result.append(bin)
    return result

  def show_remark_3(self, out = None):
    if(out is None): out = sys.stdout
    pr = "REMARK   3  "
    print >> out,pr+"REFINEMENT TARGET : %-s"%self.target_name.upper()
    print >> out,pr
    print >> out,pr+"DATA USED IN REFINEMENT."
    print >> out,pr+" RESOLUTION RANGE HIGH (ANGSTROMS) : %-8.3f"%self.d_min
    print >> out,pr+" RESOLUTION RANGE LOW  (ANGSTROMS) : %-8.3f"%self.d_max
    if(self.sigma_fobs_rejection_criterion is None and
       self.sigma_iobs_rejection_criterion is None):
      print >> out,pr+" DATA CUTOFF            (SIGMA(F)) : None"
    elif(self.sigma_fobs_rejection_criterion is not None):
      print >> out,pr+" DATA CUTOFF            (SIGMA(F)) : %-6.2f"%\
        self.sigma_fobs_rejection_criterion
    elif(self.sigma_iobs_rejection_criterion is not None):
      print >> out,pr+" DATA CUTOFF            (SIGMA(I)) : %-6.2f"%\
        self.sigma_iobs_rejection_criterion
    print >> out,pr+" COMPLETENESS FOR RANGE        (%s) : %-6.2f"%\
      ("%", self.completeness_in_range*100.0)
    print >> out,pr+" NUMBER OF REFLECTIONS             : %-10d"%\
      self.number_of_reflections
    print >> out,pr
    print >> out,pr+"FIT TO DATA USED IN REFINEMENT."
    print >> out,pr+" R VALUE     (WORKING + TEST SET) : %-6.4f"%self.r_all
    print >> out,pr+" R VALUE            (WORKING SET) : %-6.4f"%self.r_work
    print >> out,pr+" FREE R VALUE                     : %-6.4f"%self.r_free
    print >> out,pr+" FREE R VALUE TEST SET SIZE   (%s) : %-6.2f"%("%",
      float(self.number_of_test_reflections)/self.number_of_reflections*100.)
    print >> out,pr+" FREE R VALUE TEST SET COUNT      : %-10d"%\
      self.number_of_test_reflections
    print >> out,pr
    print >> out,pr+"FIT TO DATA USED IN REFINEMENT (IN BINS)."
    print >> out,pr+" BIN  RESOLUTION RANGE  COMPL.   NWORK NFREE  RWORK RFREE"
    fmt = " %3d %-17s %6.2f %8d %5d  %5.2f %5.2f"
    for bin in self.bins:
      print >> out,pr+fmt%(bin.i_bin, bin.d_range, bin.completeness*100.,
        bin.n_work, bin.n_free, bin.r_work*100., bin.r_free*100.)
    print >> out,pr
    print >> out,pr+"BULK SOLVENT MODELLING."
    print >> out,pr+" METHOD USED        : FLAT BULK SOLVENT MODEL"
    print >> out,pr+" SOLVENT RADIUS     : %-8.2f"%self.mask_solvent_radius
    print >> out,pr+" SHRINKAGE RADIUS   : %-8.2f"%self.mask_shrink_radius
    print >> out,pr+" GRID STEP FACTOR   : %-8.2f"%self.mask_grid_step_factor
    print >> out,pr+" K_SOL              : %-8.3f"%self.k_sol
    print >> out,pr+" B_SOL              : %-8.3f"%self.b_sol
    print >> out,pr
    print >> out,pr+"ERROR ESTIMATES."
    print >> out,pr+\
      " COORDINATE ERROR (MAXIMUM-LIKELIHOOD BASED)     : %-8.2f"%\
      self.ml_coordinate_error
    print >> out,pr+\
      " PHASE ERROR (DEGREES, MAXIMUM-LIKELIHOOD BASED) : %-8.2f"%\
      self.ml_phase_error
    print >> out,pr
    print >> out,pr+"OVERALL SCALE FACTORS."
    print >> out,pr+\
      " SCALE = SUM(|F_OBS|*|F_MODEL|)/SUM(|F_MODEL|**2) : %-12.4f"%\
      self.overall_scale_k1
    print >> out,pr+" ANISOTROPIC SCALE MATRIX ELEMENTS (IN CARTESIAN BASIS)."
    print >> out,pr+"  B11 : %-15.4f"%self.b_cart[0]
    print >> out,pr+"  B22 : %-15.4f"%self.b_cart[1]
    print >> out,pr+"  B33 : %-15.4f"%self.b_cart[2]
    print >> out,pr+"  B12 : %-15.4f"%self.b_cart[3]
    print >> out,pr+"  B13 : %-15.4f"%self.b_cart[4]
    print >> out,pr+"  B23 : %-15.4f"%self.b_cart[5]
    print >> out,pr
    print >> out,pr+"R FACTOR FORMULA."
    print >> out,pr+" R = SUM(||F_OBS|-SCALE*|F_MODEL||)/SUM(|F_OBS|)"
    print >> out,pr
    print >> out,pr+"TOTAL MODEL STRUCTURE FACTOR (F_MODEL)."
    print >> out,pr+" F_MODEL = FB_CART * (F_CALC_ATOMS + F_BULK)"
    print >> out,pr+"  F_BULK = K_SOL * EXP(-B_SOL * S**2 / 4) * F_MASK"
    print >> out,pr+"  F_CALC_ATOMS = ATOMIC MODEL STRUCTURE FACTORS"
    print >> out,pr+"  FB_CART = EXP(-H(t) * A(-1) * B * A(-1t) * H)"
    print >> out,pr+"   A = orthogonalization matrix, H = MILLER INDEX"
    print >> out,pr+"   (t) = TRANSPOSE, (-1) = INVERSE"
    print >> out,pr
    print >> out,pr+"STRUCTURE FACTORS CALCULATION ALGORITHM : %-s"%\
      self.sf_algorithm.upper()

class info(object):
  def __init__(self, model,
                     fmodel_x          = None,
                     fmodel_n          = None,
                     refinement_params = None):
    ref_par = refinement_params
    self.model = mmtbx.model_statistics.model(model = model)
    self.data_x, self.data_n = None, None
    if(fmodel_x is not None):
      self.data_x = data(fmodel = fmodel_x,
        free_reflections_per_bin = ref_par.alpha_beta.free_reflections_per_bin,
        max_number_of_bins = ref_par.main.max_number_of_resolution_bins,
        sigma_fobs_rejection_criterion =
          ref_par.main.sigma_fobs_rejection_criterion,
        sigma_iobs_rejection_criterion =
          ref_par.main.sigma_iobs_rejection_criterion)
    if(fmodel_n is not None):
      self.data_n = data(fmodel = fmodel_n,
        free_reflections_per_bin = ref_par.alpha_beta.free_reflections_per_bin,
        max_number_of_bins = ref_par.main.max_number_of_resolution_bins,
        sigma_fobs_rejection_criterion =
          ref_par.main.sigma_fobs_rejection_criterion,
        sigma_iobs_rejection_criterion =
          ref_par.main.sigma_iobs_rejection_criterion)

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
