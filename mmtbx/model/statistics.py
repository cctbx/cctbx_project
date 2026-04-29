from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import sys, math
from libtbx.str_utils import format_value, round_2_for_cif, round_4_for_cif
from itertools import count
import iotbx.cif.model
from libtbx.test_utils import approx_equal
from libtbx import group_args
from libtbx.utils import null_out
from mmtbx.validation import rama_z
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.validation.rotalyze import rotalyze
from mmtbx.validation.cbetadev import cbetadev
from mmtbx.validation.clashscore import clashscore
from mmtbx.validation.utils import molprobity_score
from mmtbx.validation import omegalyze
from mmtbx.validation import cablam
from cctbx import adptbx
import six
from libtbx import Auto

class geometry(object):
  def __init__(self,
               model,
               fast_clash=False,
               condensed_probe=False,
               use_hydrogens=True):
    self.model = model
    self.pdb_hierarchy = model.get_hierarchy()
    self.fast_clash = fast_clash
    self.condensed_probe = condensed_probe
    self.restraints_source = None
    self.from_restraints = None
    self.use_hydrogens = use_hydrogens
    self.cached_result = None
    self.cached_clash = None
    self.cached_rama = None
    self.cached_rota = None
    self._init(self.pdb_hierarchy, model.restraints_manager.geometry)

  def _init(self, pdb_hierarchy=None, geometry_restraints_manager=None):
    # XXX Really, this should be part of constructor (to avoid confusion)!
    if(pdb_hierarchy is not None):
      self.pdb_hierarchy = pdb_hierarchy
    if(geometry_restraints_manager is not None):
      self.restraints_source = geometry_restraints_manager.get_source()
      sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
      working_geometry_restraints_manager = geometry_restraints_manager
      if(not self.use_hydrogens):
        asc = self.pdb_hierarchy.atom_selection_cache()
        sel = asc.selection("not (element H or element D)")
        sites_cart = sites_cart.select(sel)
        working_geometry_restraints_manager = geometry_restraints_manager.select(sel)
      self.from_restraints = \
        working_geometry_restraints_manager.energies_sites(
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

  def angle(self, origin_id=Auto, return_rmsZ=False):
    mi,ma,me,n = 0,0,0,0
    outliers = 0
    if(self.from_restraints is not None):
      if return_rmsZ:
        mi,ma,me,n = self.from_restraints.angle_deviations_z(origin_id=origin_id)
      else:
        mi,ma,me,n = self.from_restraints.angle_deviations(origin_id=origin_id)
      # n = self.from_restraints.get_filtered_n_angle_proxies(origin_id=origin_id)
      outliers = self.from_restraints.get_angle_outliers(
        sites_cart = self.pdb_hierarchy.atoms().extract_xyz(),
        sigma_threshold=4,
        origin_id=origin_id,
        )
    return group_args(min = mi, max = ma, mean = me, n = n, outliers = outliers)

  def bond(self, origin_id=Auto, return_rmsZ=False):
    mi,ma,me,n = 0,0,0,0
    outliers = 0
    if(self.from_restraints is not None):
      if return_rmsZ:
        mi,ma,me,n = self.from_restraints.bond_deviations_z(origin_id=origin_id)
      else:
        mi,ma,me,n = self.from_restraints.bond_deviations(origin_id=origin_id)
      # n = self.from_restraints.get_filtered_n_bond_proxies(origin_id=origin_id)
      outliers = self.from_restraints.get_bond_outliers(
        sites_cart = self.pdb_hierarchy.atoms().extract_xyz(),
        sigma_threshold=4,
        origin_id=origin_id,
        )
    return group_args(min = mi, max = ma, mean = me, n = n, outliers = outliers)

  def chirality(self):
    mi,ma,me,n = 0,0,0,0
    if(self.from_restraints is not None):
      mi,ma,me = self.from_restraints.chirality_deviations()
      n = self.from_restraints.n_chirality_proxies
      outliers = self.from_restraints.get_chirality_outliers(
        sites_cart = self.pdb_hierarchy.atoms().extract_xyz(),
        sigma_threshold=4)
    return group_args(min = mi, max = ma, mean = me, n = n, outliers = outliers)

  def dihedral(self, return_rmsZ=False):
    mi,ma,me,n = 0,0,0,0
    outliers = 0
    if(self.from_restraints is not None):
      if return_rmsZ:
        mi,ma,me,n = self.from_restraints.dihedral_deviations_z()
      else:
        mi,ma,me = self.from_restraints.dihedral_deviations()
        n = self.from_restraints.get_filtered_n_dihedral_proxies()
      outliers = self.from_restraints.get_dihedral_outliers(
        sites_cart = self.pdb_hierarchy.atoms().extract_xyz(),
        sigma_threshold=4)
    return group_args(min = mi, max = ma, mean = me, n = n, outliers = outliers)

  def planarity(self):
    mi,ma,me,n = 0,0,0,0
    if(self.from_restraints is not None):
      mi,ma,me = self.from_restraints.planarity_deviations()
      n = self.from_restraints.get_filtered_n_planarity_proxies()
    pp = protein_planes(hierarchy = self.pdb_hierarchy)
    asc = self.pdb_hierarchy.atom_selection_cache()
    sel = asc.selection("element H or element D")
    pp_noH = None
    if sel.count(True) > 0:
      pp_noH = protein_planes(hierarchy = self.pdb_hierarchy.select(~sel))
    return group_args(min = mi, max = ma, mean = me, n = n,
      protein_planes_max_dev = pp,
      protein_planes_max_dev_noH = pp_noH)

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
    if self.cached_rama is None:
      self.cached_rama = ramalyze(
          pdb_hierarchy = self.pdb_hierarchy,
          outliers_only = False)
    return group_args(
      outliers = self.cached_rama.percent_outliers,
      allowed  = self.cached_rama.percent_allowed,
      favored  = self.cached_rama.percent_favored,
      ramalyze = self.cached_rama #XXX Bulky object -- REMOVE!
      )

  def rotamer(self):
    if self.cached_rota is None:
      self.cached_rota = rotalyze(
          pdb_hierarchy = self.pdb_hierarchy,
          outliers_only = False)
    return group_args(
      outliers = self.cached_rota.percent_outliers,
      favored  = self.cached_rota.percent_favored,
      allowed  = self.cached_rota.percent_allowed,
      rotalyze = self.cached_rota #XXX Bulky object -- REMOVE!
      )

  def c_beta(self):
    result = cbetadev(pdb_hierarchy = self.pdb_hierarchy,
      outliers_only = True, out = null_out()) # XXX Why it is different from others?
    return group_args(
      outliers = result.get_outlier_percent(),
      cbetadev = result #XXX Bulky object -- REMOVE!
      )

  def clash(self):
    if self.cached_clash is None:
      self.cached_clash = clashscore(
        pdb_hierarchy   = self.pdb_hierarchy,
        fast            = self.fast_clash,
        condensed_probe = self.condensed_probe,
        keep_hydrogens  = self.use_hydrogens,
        nuclear         = self.model.is_neutron())
    return group_args(
      score   = self.cached_clash.get_clashscore(),
      clashes = self.cached_clash #XXX Bulky object -- REMOVE! - Not kidding,
          # it contains 1 probe output, which can be GigaBytes!
      )

  def cablam(self):
    result = cablam.cablamalyze(self.pdb_hierarchy, outliers_only=False,
      out=null_out(), quiet=True) # XXX Why it is different from others?
    gui_table = group_args(column_labels = result.gui_list_headers,
                           column_formats = result.gui_formats,
                           data = result.as_gui_table_data(include_zoom=True),
                           column_widths = result.wx_column_widths)
    return group_args(
      outliers    = result.percent_outliers(),
      disfavored  = result.percent_disfavored(),
      ca_outliers = result.percent_ca_outliers(),
      gui_table   = gui_table)

  def rama_z_score(self):
    return rama_z.rama_z(models = [self.model], log = null_out()).get_result()

  def omega(self):
    result = omegalyze.omegalyze(pdb_hierarchy=self.pdb_hierarchy, quiet=True)
    # XXX Move this to omegalyze function.
    n_proline         = result.n_proline()
    n_general         = result.n_general()
    n_cis_proline     = result.n_cis_proline()
    n_cis_general     = result.n_cis_general()
    n_twisted_proline = result.n_twisted_proline()
    n_twisted_general = result.n_twisted_general()
    cis_general       = 0
    twisted_general   = 0
    cis_proline       = 0
    twisted_proline   = 0
    if(n_proline != 0):
      cis_proline     = n_cis_proline    *100./n_proline
      twisted_proline = n_twisted_proline*100./n_proline
    if(n_general != 0):
      cis_general     = n_cis_general    *100./n_general
      twisted_general = n_twisted_general*100./n_general
    return group_args(
      cis_proline       = cis_proline,
      cis_general       = cis_general,
      twisted_general   = twisted_general,
      twisted_proline   = twisted_proline,
      n_twisted_general = n_twisted_general,
      omegalyze         = result #XXX Bulky object -- REMOVE!
      )

  def mp_score(self):
    return molprobity_score(
      clashscore = self.clash().score,
      rota_out   = self.rotamer().outliers,
      rama_fav   = self.ramachandran().favored)

  def result(self, slim=False):
    from libtbx import Auto
    if(self.cached_result is None):
      self.cached_result = group_args(
         angle            = self.angle(),
         bond             = self.bond(),
         angle_z          = self.angle(return_rmsZ=True),
         bond_z           = self.bond(return_rmsZ=True),
         chirality        = self.chirality(),
         dihedral         = self.dihedral(),
         planarity        = self.planarity(),
         parallelity      = self.parallelity(),
         nonbonded        = self.nonbonded(),
         ramachandran     = self.ramachandran(),
         rotamer          = self.rotamer(),
         c_beta           = self.c_beta(),
         clash            = self.clash(),
         molprobity_score = self.mp_score(),
         cablam           = self.cablam(),
         omega            = self.omega(),
         rama_z           = self.rama_z_score())
    if(slim):
      delattr(self.cached_result.ramachandran, "ramalyze")
      delattr(self.cached_result.clash,        "clashes")
      delattr(self.cached_result.rotamer,      "rotalyze")
      delattr(self.cached_result.cablam,       "gui_table")
      delattr(self.cached_result.omega,        "omegalyze")
      delattr(self.cached_result.c_beta,       "cbetadev")
      delattr(self.cached_result.angle,        "outliers")
      delattr(self.cached_result.bond,         "outliers")
      delattr(self.cached_result.dihedral,     "outliers")
    return self.cached_result

  def show_bond_and_angle(self):
    r = self.result()
    f="bond: %6.3f angle: %6.2f"
    return f%(r.bond.mean, r.angle.mean)

  def show_bond_and_angle_and_dihedral(self, assert_zero=False):
    r = self.result()
    f="bond: %6.3f angle: %6.2f dihedral: %6.2f"
    if assert_zero:
      assert abs(r.bond.mean)<1e-3, 'bond rmsd is not zero : %0.5f' % r.bond.mean
      assert abs(r.angle.mean)<1e-3, 'angle rmsd is not zero : %0.5f' % r.angle.mean
      assert abs(r.dihedral.mean)<1e-3, 'dihedral rmsd is not zero : %0.5f' % r.dihedral.mean
    return f%(r.bond.mean, r.angle.mean, r.dihedral.mean)

  def show_planarity_details(self):
    r = self.result()
    f='planarity (n=%d) mean: %6.3f min:  %6.3f max:  %6.3f'
    return f%(r.planarity.n, r.planarity.mean, r.planarity.min, r.planarity.max)

  def show_short(self):
    r = self.result()
    f="%s %s %s %s %s %s %s %s"
    return f%(
      format_value("bond: %6.3f",r.bond.mean),
      format_value("angle: %6.2f", r.angle.mean),
      format_value("clash: %5.1f", r.clash.score),
      format_value("rota: %5.2f", r.rotamer.outliers),
      format_value("rama: f: %6.2f", r.ramachandran.favored),
      format_value("o: %6.2f", r.ramachandran.outliers),
      format_value("Z: %6.2f", r.rama_z.whole.value),
      format_value("cb: %6.2f", r.c_beta.outliers) )

  def number_of_bonds(self):
    bonds = self.bond()
    return bonds.n

  def show(self, log=None, prefix="",
           exclude_protein_only_stats=False,
           include_rmsd_details=True,
           uppercase=True):
    if(log is None): log = sys.stdout
    def fmt(f1,f2,d1,z1=None):
      if f1 is None  : return '   -       -       -  '
      if z1 is None:
        fmt_str= "%6.3f %7.3f %6d"
        return fmt_str%(f1,f2,d1)
      else:
        fmt_str= "%6.3f %7.3f %6d  Z=%6.3f"
        return fmt_str%(f1,f2,d1,z1)
    def fmt2(f1):
      if f1 is None: return '  -   '
      return "%-6.3f"%(f1)
    res = self.result()
    a,b,c,d,p,n = res.angle, res.bond, res.chirality, res.dihedral, \
      res.planarity, res.nonbonded
    az, bz = res.angle_z, res.bond_z
    assert b.n>=bz.n, 'rmsd.n != rmsZ.n %s %s' % (b.n, bz.n)
    assert a.n==az.n
    result = """
%sGeometry Restraints Library: %s
%sDeviations from Ideal Values - rmsd, rmsZ for bonds and angles.
%s  Bond      : %s
%s  Angle     : %s
%s  Chirality : %s
%s  Planarity : %s
%s  Dihedral  : %s
%s  Min Nonbonded Distance : %s
%s\n"""%(prefix,
       # prefix,
       self.restraints_source,
       prefix,
       prefix, fmt(b.mean, b.max, b.n, bz.mean),
       prefix, fmt(a.mean, a.max, a.n, az.mean),
       prefix, fmt(c.mean, c.max, c.n),
       prefix, fmt(p.mean, p.max, p.n),
       prefix, fmt(d.mean, d.max, d.n),
       prefix, fmt2(n.min).strip(),
       prefix)
    if not exclude_protein_only_stats:
      result += """%sMolprobity Statistics.
%s  All-atom Clashscore : %s
%s  Ramachandran Plot:
%s    Outliers : %5.2f %%
%s    Allowed  : %5.2f %%
%s    Favored  : %5.2f %%
%s  Rotamer:
%s    Outliers : %5.2f %%
%s    Allowed  : %5.2f %%
%s    Favored  : %5.2f %%
%s  Cbeta Deviations : %5.2f %%
%s  Peptide Plane:
%s    Cis-proline     : %s %%
%s    Cis-general     : %s %%
%s    Twisted Proline : %s %%
%s    Twisted General : %s %%"""%(
        prefix,
        # prefix,
        prefix, format_value("%5.2f", res.clash.score).strip(),
        prefix,
        prefix, res.ramachandran.outliers,
        prefix, res.ramachandran.allowed,
        prefix, res.ramachandran.favored,
        prefix,
        prefix, res.rotamer.outliers,
        prefix, res.rotamer.allowed,
        prefix, res.rotamer.favored,
        prefix, res.c_beta.outliers,
        prefix,
        prefix, format_value("%5.2f", res.omega.cis_proline).strip(),
        prefix, format_value("%5.2f", res.omega.cis_general).strip(),
        prefix, format_value("%5.2f", res.omega.twisted_proline).strip(),
        prefix, format_value("%5.2f", res.omega.twisted_general).strip())
      result += "\n%s"%prefix
      result += res.rama_z.as_string(prefix=prefix)
      result += "\n%s"%prefix
    #
      # result+="\n%s"%prefix
      if p.protein_planes_max_dev:
        result += "\n%sMax deviation from planes:"%prefix
        result += "\n%s   Type  MaxDev  MeanDev LineInFile"%prefix
        for pp in p.protein_planes_max_dev:
          result += "\n%s %s %s %s  %s"%(
            prefix,
            format_value("%s", pp.resname),
            format_value("%7.3f", pp.max_dev),
            format_value("%7.3f", pp.mean_dev),
            pp.id)
        result += "\n%s\n"%prefix
      if p.protein_planes_max_dev_noH is not None:
        result += "\n\n%sMax deviation from planes (no H):"%prefix
        result += "\n%s   Type  MaxDev  MeanDev LineInFile"%prefix
        for pp in p.protein_planes_max_dev_noH:
          result += "\n%s %s %s %s  %s"%(
            prefix,
            format_value("%s", pp.resname),
            format_value("%7.3f", pp.max_dev),
            format_value("%7.3f", pp.mean_dev),
            pp.id)
        result += "\n%s\n"%prefix
    if include_rmsd_details:
      from cctbx.geometry_restraints.linking_class import linking_class
      origin_ids = linking_class()
      result += '%sDetails of bonding type rmsd' % (prefix)
      for key, i in origin_ids.items():
        bond_rc=self.bond(origin_id=-i)
        angle_rc=self.angle(origin_id=-i)
        if bond_rc.n:
          result += '\n%s  %-20s : bond   %12.5f (%5d)' % (prefix,
                                                           key,
                                                           bond_rc.mean,
                                                           bond_rc.n)
        if angle_rc.n:
          result += '\n%s  %-20s : angle  %12.5f (%5d)' % (prefix,
                                                           key,
                                                           angle_rc.mean,
                                                           angle_rc.n)
    #
    if( uppercase ):
      result = result.upper()
    print(result, file=log)

  def as_cif_block(self, cif_block=None, pdbx_refine_id=''):
    if cif_block is None:
      cif_block = iotbx.cif.model.block()
    cif_block["_refine.pdbx_stereochemistry_target_values"] = \
        self.restraints_source if self.restraints_source is not None else '?'
    loop = iotbx.cif.model.loop(header=(
      "_refine_ls_restr.pdbx_refine_id",
      "_refine_ls_restr.type",
      "_refine_ls_restr.number",
      "_refine_ls_restr.dev_ideal",
      #"_refine_ls_restr.dev_ideal_target",
      "_refine_ls_restr.weight",
      #"_refine_ls_restr.pdbx_refine_id",
      "_refine_ls_restr.pdbx_restraint_function",
    ))
    res = self.result()
    a,b,c,d,p,n = res.angle, res.bond, res.chirality, res.dihedral, \
      res.planarity, res.nonbonded
    loop.add_row((pdbx_refine_id, "f_bond_d",           b.n, round_4_for_cif(b.mean), "?", "?"))
    loop.add_row((pdbx_refine_id, "f_angle_d",          a.n, round_4_for_cif(a.mean), "?", "?"))
    loop.add_row((pdbx_refine_id, "f_chiral_restr",     c.n, round_4_for_cif(c.mean), "?", "?"))
    loop.add_row((pdbx_refine_id, "f_plane_restr",      p.n, round_4_for_cif(p.mean), "?", "?"))
    loop.add_row((pdbx_refine_id, "f_dihedral_angle_d", d.n, round_4_for_cif(d.mean), "?", "?"))
    cif_block.add_loop(loop)
    return cif_block

class composition(object):
  def __init__(self, pdb_hierarchy):
    self._result = pdb_hierarchy.composition()

  def result(self):
    return self._result

  def show(self, log, prefix=""):
    r = self.result()
    ligs=(",".join([("%s:%s"%k).strip() for k in r.other_cnts.items()])).strip()
    if(len(ligs)==0): ligs=None
    print(prefix, "Number of:", file=log)
    print(prefix, "  all atoms      :", r.n_atoms, file=log)
    print(prefix, "  H or D atoms   :", r.n_hd, file=log)
    print(prefix, "  chains         :", r.n_chains, file=log)
    print(prefix, "  a.a. residues  :", r.n_protein, file=log)
    print(prefix, "  nucleotides    :", r.n_nucleotide, file=log)
    print(prefix, "  water          :", r.n_water, file=log)
    print(prefix, "  other (ligands):", r.n_other, file=log)
    print(prefix, "Ligands:", ligs, file=log)

class occupancy(object):
  def __init__(self, hierarchy):
    self._result = hierarchy.occupancy_counts()

  def result(self):
    return self._result

  def show(self, log, prefix=""):
    r = self.result()
    def p(m): print(prefix, m, file=log)
    p("mean     : %4.2f"%r.mean)
    p("negative : %d"%r.negative)
    p("zero     : %d (%-6.2f%s)"%(r.zero_count,r.zero_fraction,"%"))
    p("occ>1    : %d (%-6.2f%s)"%(r.greater_than_1_count,r.greater_than_1_fraction,"%"))
    p("altlocs  : %-6.2f(%s)"%(r.alt_conf_frac,"%"))

def rms_b_iso_or_b_equiv_bonded(
      geometry_restraints_manager,
      sites_cart,
      b_isos):
  bond_proxies_simple, asu = geometry_restraints_manager.geometry.\
    get_covalent_bond_proxies(sites_cart=sites_cart)
  values = flex.double()
  result = None
  for proxy in bond_proxies_simple:
    i_seq, j_seq = proxy.i_seqs
    b_i = b_isos[i_seq]
    b_j = b_isos[j_seq]
    abs_diff_sq = abs(b_i-b_j)**2
    values.append(abs_diff_sq)
  if(values.size() > 0):
    result = math.sqrt(flex.sum(values) / values.size())
  return result

class adp(object):
  def __init__(self,
               pdb_hierarchy,
               xray_structure,
               use_hydrogens=False,
               geometry_restraints_manager=None,
               use_z_scores=False):
    if(not use_hydrogens) or use_z_scores:
      not_hd_sel = ~xray_structure.hd_selection()
      pdb_hierarchy  = pdb_hierarchy.select(not_hd_sel)
      xray_structure = xray_structure.select(not_hd_sel)
      if(geometry_restraints_manager is not None):
        geometry_restraints_manager = \
          geometry_restraints_manager.select(not_hd_sel)
    b_isos = xray_structure.extract_u_iso_or_u_equiv() * adptbx.u_as_b(1.)
    if use_z_scores: mean = b_isos.as_z_scores()
    sites_cart = xray_structure.sites_cart()
    asc = pdb_hierarchy.atom_selection_cache()
    def get_stats(sel_str, rms_bonded=False):
      sel = asc.selection(sel_str)
      xrs = xray_structure.select(sel)
      n_iso   = xrs.use_u_iso().count(True)
      n_aniso = xrs.use_u_aniso().count(True)
      anisotropy = xrs.scatterers().anisotropy(unit_cell = xrs.unit_cell())
      if(sel.count(True)==0): return None
      b_isos_selected = b_isos.select(sel)
      sites_cart_selected = sites_cart.select(sel)
      mi,ma,me = b_isos_selected.min_max_mean().as_tuple()
      rms_b_iso_bonded = None
      if(rms_bonded and geometry_restraints_manager is not None):
        grm = geometry_restraints_manager.select(sel)
        rms_b_iso_bonded = rms_b_iso_or_b_equiv_bonded(
          geometry_restraints_manager = geometry_restraints_manager.select(sel),
          sites_cart                  = sites_cart_selected,
          b_isos                      = b_isos_selected)
      return group_args(
        min              = mi,
        max              = ma,
        mean             = me,
        n_iso            = n_iso,
        n_aniso          = n_aniso,
        n_zero           = (b_isos_selected < 0.01).count(True),
        rms_b_iso_bonded = rms_b_iso_bonded)
    overall    = get_stats(sel_str="all",        rms_bonded=True)
    protein    = get_stats(sel_str="protein",    rms_bonded=True)
    nucleotide = get_stats(sel_str="nucleotide", rms_bonded=True)
    hd         = get_stats(sel_str="element H or element D")
    water      = get_stats(sel_str="water")
    other      = get_stats(sel_str="not (water or nucleotide or protein)")
    chains = {}
    for chain in pdb_hierarchy.chains():
      chains[chain.id] = get_stats(sel_str="chain '%s'"%chain.id)
    histogram = flex.histogram(data = b_isos, n_slots = 10)
    self._result = group_args(
      overall    = overall,
      protein    = protein,
      nucleotide = nucleotide,
      hd         = hd,
      water      = water,
      other      = other,
      chains     = chains,
      histogram  = histogram)

  def result(self):
    return self._result

  def show(self, log, prefix=""):
    def format_str(v):
      rms = "   N/A"
      if(v.rms_b_iso_bonded is not None): rms = "%6.2f"%v.rms_b_iso_bonded
      return "%6.2f %6.2f %6.2f %s %5d %5d"%(
        v.min, v.max, v.mean, rms, v.n_iso, v.n_aniso)
    r = self.result()
    pad=" "*14
    print(prefix, "Individual atomic B", file=log)
    print(prefix, pad, "min    max   mean <Bi,j>   iso aniso", file=log)
    if(r.overall is not None):
      print(prefix, "  Overall: ", format_str(r.overall), file=log)
    if(r.protein is not None):
      print(prefix, "  Protein: ", format_str(r.protein), file=log)
    if(r.nucleotide is not None):
      print(prefix, "  RNA/DNA: ", format_str(r.nucleotide), file=log)
    if(r.hd is not None):
      print(prefix, "  H and D: ", format_str(r.hd), file=log)
    if(r.water is not None):
      print(prefix, "  Water:   ", format_str(r.water), file=log)
    if(r.other is not None):
      print(prefix, "  Other:   ", format_str(r.other), file=log)
    for k,v in six.iteritems( r.chains):
      print(prefix, "  Chain %2s:"%k, format_str(v), file=log)
    def show_hist(h):
      print(prefix, "      Values      Number of atoms", file=log)
      lc_1 = h.data_min()
      s_1 = enumerate(h.slots())
      for (i_1,n_1) in s_1:
        hc_1 = h.data_min() + h.slot_width() * (i_1+1)
        print(prefix, "  %6.2f - %-6.2f %8d" % (lc_1, hc_1, n_1), file=log)
        lc_1 = hc_1
    if(r.overall is not None and abs(r.overall.min-r.overall.max)>1.e-3):
      print(prefix, "  Histogram:", file=log)
      show_hist(r.histogram)

  def as_cif_block(self, cif_block=None):
    if cif_block is None:
      cif_block = iotbx.cif.model.block()
    if self._result.overall is not None:
      cif_block["_refine.B_iso_mean"] = round_2_for_cif(self._result.overall.mean)
    else:
      cif_block["_refine.B_iso_mean"] = '?'
    return cif_block

class info(object):
  def __init__(self, model,
                     fmodel_x          = None,
                     fmodel_n          = None,
                     refinement_params = None):
    ref_par = refinement_params
    self.wilson_b = None
    self.model = model
    if fmodel_x is not None:
      self.wilson_b = fmodel_x.wilson_b()
    elif fmodel_n is not None:
      self.wilson_b = fmodel_n.wilson_b()
    self.geometry = model.geometry_statistics()
    self.adp = model.adp_statistics()
    self.data_x, self.data_n = None, None
    if(fmodel_x is not None):
      self.data_x = fmodel_x.info(
        free_reflections_per_bin = ref_par.alpha_beta.free_reflections_per_bin,
        max_number_of_bins       = ref_par.main.max_number_of_resolution_bins)
    if(fmodel_n is not None):
      self.data_n = fmodel_n.info(
        free_reflections_per_bin = ref_par.alpha_beta.free_reflections_per_bin,
        max_number_of_bins       = ref_par.main.max_number_of_resolution_bins)

    self._pdbx_refine_id = ''
    if self.data_x is not None:
      self._pdbx_refine_id = 'X-RAY DIFFRACTION'
    if self.data_n is not None:
      # !!! Warning: "X-ray+Neutron" is not compliant with mmCIF dictionary:
      # http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_exptl.method.html
      self._pdbx_refine_id = 'NEUTRON DIFFRACTION' if self.data_x is None else 'X-ray+Neutron'
    if self._pdbx_refine_id == '':
      # most likely electron microscopy, but checking scattering table anyway
      if self.model.get_xray_structure().scattering_type_registry().last_table() == "electron":
        self._pdbx_refine_id = 'ELECTRON MICROSCOPY'


  def show_remark_3(self, out = None):
    prefix = "REMARK   3  "
    if(out is None): out = sys.stdout
    if(self.data_x is not None):
      print(prefix+"X-RAY DATA.", file=out)
      print(prefix, file=out)
      self.data_x.show_remark_3(out = out)
    if(self.data_n is not None):
      print(prefix+"NEUTRON DATA.", file=out)
      print(prefix, file=out)
      self.data_n.show_remark_3(out = out)
      print(prefix, file=out, end='')
    if self.adp is not None:
      self.adp.show(log=out, prefix=prefix)
      print(prefix, file=out, end='')
    if(self.geometry is not None):
      self.geometry.show(log=out, prefix=prefix)
    for info_pdb_str in [self.model.tls_groups_as_pdb(),
        self.model.anomalous_scatterer_groups_as_pdb(),
        self.model.cartesian_NCS_as_pdb(),
        self.model.torsion_NCS_as_pdb()]:
      if len(info_pdb_str) > 0:
        print(prefix, file=out)
        print(info_pdb_str, end='', file=out)

  def get_pdbx_refine_id(self):
    return self._pdbx_refine_id

  def as_cif_block(self, cif_block=None):
    if cif_block is None:
      cif_block = iotbx.cif.model.block()
    if self.data_x is not None:
      cif_block = self.data_x.as_cif_block(cif_block=cif_block, scattering_type=self._pdbx_refine_id)
    if self.wilson_b is not None:
      cif_block["_reflns.B_iso_Wilson_estimate"] = round_2_for_cif(self.wilson_b)
    # XXX Neutron data?
    if self.geometry is not None:
      cif_block = self.geometry.as_cif_block(cif_block=cif_block, pdbx_refine_id=self._pdbx_refine_id)
    if self.adp is not None:
      cif_block = self.adp.as_cif_block(cif_block=cif_block)
    cif_block = self.model.tls_groups_as_cif_block(cif_block=cif_block)
    # What about anomalous_scatterer_groups here?
    return cif_block

def protein_planes(hierarchy):
  from cctbx import geometry_restraints
  # XXX This needs to be generalized !!! XXX
  planes = {
  "ARG":     [' NE ', ' CZ ', ' NH1', ' NH2', ' HE ', 'HH11', 'HH12', 'HH21', 'HH22'],
  "TYR":     ["CB","CG", "CD1","CD2","CE1","CE2","CZ", "OH", "HD1","HD2","HE1","HE2"],
  "PHE":     ["CB","CG", "CD1","CD2","CE1","CE2","CZ", "OH", "HD1","HD2","HE1","HE2"],
  "TRP":     ["CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2","HD1","HE1","HE3","HZ2","HZ3","HH2"],
  "HIS":     [ "CG","ND1","CE1","NE2","CD2","HCD2","HNE2","HND1"]
  }
  asc = hierarchy.atom_selection_cache()
  result = []
  for resname in planes.keys():
    MAX = -1
    ID = None
    RESNAME = None
    deltas = flex.double()
    names = planes[resname]
    ss = ss = "resname %s and (%s)"%(resname," or ".join(["name %s"%it for it in names]))
    sel = asc.selection(string=ss)
    h = hierarchy.select(sel)
    for residue_group in h.residue_groups():
      if len(residue_group.atom_groups())>1: continue
      atom_group = residue_group.atom_groups()[0]
      assert atom_group.resname==resname
      atoms = atom_group.atoms()
      size = atoms.size()
      sites_cart = atoms.extract_xyz()
      pproxy = geometry_restraints.planarity_proxy(
        i_seqs  = list(range(size)),
        weights = flex.double([1]*size))
      planarity = geometry_restraints.planarity(
        sites_cart = sites_cart,
        proxy      = pproxy)
      abs_deltas = flex.abs(planarity.deltas())
      deltas.extend(abs_deltas)
      max_delta  = flex.max(abs_deltas)
      if max_delta > MAX:
        MAX = max_delta
        ID = atom_group.id_str()
        RESNAME = resname
    if(MAX > 0 and deltas.size() > 0):
      assert approx_equal(MAX, flex.max(deltas), 1.e-4)
    r = group_args(
      resname  = RESNAME,
      max_dev  = MAX,
      mean_dev = flex.mean(deltas) if deltas.size()>0 else None,
      id       = ID)
    if MAX>=0:
      result.append(r)
  return result
