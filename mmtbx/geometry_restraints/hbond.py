
from libtbx import easy_mp
from libtbx.str_utils import make_header
from libtbx.str_utils import pad_string as ps
import libtbx.phil
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from libtbx import adopt_init_args, group_args
from math import sqrt
import sys

master_phil = libtbx.phil.parse("""
  restraint_type = *Auto simple lennard_jones implicit
    .type = choice
    .short_caption = Hydrogen bond restraint type
    .caption = Automatic Simple_(H-O) Simple_(N-O) Angle-dependent_(N-O)
  include_side_chains = True
    .type = bool
  optimize_hbonds = False
    .type = bool
  optimize_hbonds_thorough = False
    .type = bool
  optimize_mode = *first last every_macro_cycle
    .type = choice
    .expert_level = 3
  restraints_weight = 1.0
    .type = float
  falloff_distance = 0.05
    .type = float
  exclude_nonbonded = True
    .type = bool
  distance_ideal_h_o = 1.975
    .type = float
  distance_cut_h_o = 2.5
    .type = float
  distance_ideal_n_o = 2.9
    .type = float
  distance_cut_n_o = 3.5
    .type = float
  implicit
    .short_caption = Implicit hydrogens
    .help = Based on H-bond potential for CNS by Chapman lab
  {
    theta_high = 155
      .type = float
    theta_low = 115
      .type = float
    theta_cut = 90
      .type = float
  }
  explicit
    .short_caption = Explicit hydrogens
    .help = Similar to Rosetta H-bond energy (Kortemme & Baker)
  {
    theta_ideal = 180
      .type = float
    theta_sigma = 5
      .type = float
    psi_ideal = 155
      .type = float
    psi_sigma = 5
      .type = float
    relative_weights = 1.0 1.0 1.0
      .type = floats(size=3)
  }
  lennard_jones {
    potential = *4_6 6_12
      .type = choice
  }
  simple
    .short_caption = Simple distance-based potentials
    .help = Pseudo-bond restraints
  {
    sigma = 0.05
      .type = float
    slack = 0.0
      .type = float
  }
""")

# XXX this is gross, but I'd prefer to delay importing yet another shared
# library until the last minute, in order to reduce startup time
class core (object) :
  def __init__ (self) :
    self.proxies = []
    self.exclude_nb_list = []

  def add_nonbonded_exclusion (self, i_seq, j_seq) :
    self.exclude_nb_list.append((i_seq, j_seq))

class build_proxies (core) :
  proxy_array_type = None #"shared_h_bond_simple_proxy"
  proxy_type = None
  def __init__ (self) :
    core.__init__(self)
    import cctbx.geometry_restraints # import dependency
    import boost.python
    self.ext = boost.python.import_ext("mmtbx_hbond_restraints_ext")
    self.proxies = getattr(self.ext, self.proxy_array_type)()

  def add_proxy (self, **kwds) :
    proxy_class = getattr(self.ext, self.proxy_type)
    proxy = proxy_class(**kwds)
    self.proxies.append(proxy)

class build_simple_hbond_proxies (build_proxies) :
  proxy_array_type = "shared_h_bond_simple_proxy"
  proxy_type = "h_bond_simple_proxy"

class build_lennard_jones_proxies (build_proxies) :
  proxy_array_type = "shared_h_bond_lennard_jones_proxy"
  proxy_type = "h_bond_lj_proxy"

# Fabiola et al. (2002) Protein Sci. 11:1415-23
# http://www.ncbi.nlm.nih.gov/pubmed/12021440
class build_implicit_hbond_proxies (build_proxies) :
  proxy_array_type = "shared_h_bond_implicit_proxy"
  proxy_type = "h_bond_implicit_proxy"

class explicit_proxy (object) :
  def __init__ (self,
                i_seqs, # donor, H, acceptor, acceptor base
                distance_ideal,
                distance_cut,
                theta_ideal,
                psi_ideal,
                weight=1.0,
                relative_weights=(1.0,1.0,1.0)) :
    assert (len(relative_weights) == 3)
    assert (len(i_seqs) == 4)
    assert (distance_cut is None) or (distance_cut > distance_ideal)
    adopt_init_args(self, locals())

class build_explicit_hbond_proxies (core) :
  def add_proxy (self, **kwds) :
    proxy = explicit_proxy(**kwds)
    self.proxies.append(proxy)

# This is not used for actual restraints, but for storing data to be output
# in other formats (e.g. PyMOL, REFMAC, kinemage, etc.)
class distance_proxy (group_args) :
  pass

class build_distance_proxies (core) :
  def add_proxy (self, **kwds) :
    proxy = distance_proxy(**kwds)
    self.proxies.append(proxy)

def target_and_gradients (proxies,
                          sites_cart,
                          gradient_array=None,
                          falloff_distance=0.05,
                          use_finite_differences=True,
                          lennard_jones_potential="4_6") :
  from scitbx.array_family import flex
  import boost.python
  ext = boost.python.import_ext("mmtbx_hbond_restraints_ext")
  if (gradient_array is None) :
    gradient_array = flex.vec3_double(sites_cart.size(), (0,0,0))
  sum = 0.0
  if (type(proxies).__name__ == "shared_h_bond_simple_proxy") :
    sum = ext.h_bond_simple_residual_sum(
      sites_cart=sites_cart,
      proxies=proxies,
      gradient_array=gradient_array,
      falloff_distance=falloff_distance)
  elif (type(proxies).__name__ == "shared_h_bond_lennard_jones_proxy") :
    if (lennard_jones_potential == "4_6") :
      weight_scale = 400
      a = 6
      b = 4
      sigma_base = 0.81649658092772603
    elif (lennard_jones_potential == "6_12") :
      weight_scale = 250
      a = 12
      b = 6
      sigma_base = 0.89089871814033927
    sum = ext.h_bond_lennard_jones_residual_sum(
      sites_cart=sites_cart,
      proxies=proxies,
      gradient_array=gradient_array,
      falloff_distance=falloff_distance,
      a=a,
      b=b,
      scale=weight_scale,
      sigma_base=sigma_base,
      use_finite_differences=use_finite_differences)
  elif (type(proxies).__name__ == "shared_h_bond_implicit_proxy") :
    sum = ext.h_bond_implicit_residual_sum(
      sites_cart=sites_cart,
      proxies=proxies,
      gradient_array=gradient_array,
      falloff_distance=falloff_distance)
  else :
    assert 0
  return sum

def get_simple_bonds (proxies) :
  from scitbx.array_family import shared # import dependency
  import boost.python
  ext = boost.python.import_ext("mmtbx_hbond_restraints_ext")
  if (type(proxies).__name__ == "shared_h_bond_simple_proxy") :
    return ext.simple_hbonds_as_simple_bonds(proxies)
  elif (type(proxies).__name__ == "shared_h_bond_lennard_jones_proxy") :
    return ext.lj_hbonds_as_simple_bonds(proxies)
  elif (type(proxies).__name__ == "shared_h_bond_implicit_proxy") :
    return ext.implicit_hbonds_as_simple_bonds(proxies)
  elif isinstance(proxies, list) :
    return shared.stl_set_unsigned(get_simple_bond_equivalents(proxies))
  else :
    assert 0

def get_simple_bond_equivalents (proxies) :
  """
Given a list of proxies, extract the pair of atom indices representing the
actual "bond", e.g. H-O or N-O.  These are used to exclude the atom pair from
the nonbonded interaction restraints function.
"""
  bond_pairs = []
  for proxy in proxies :
    bond_pairs.append(_get_simple_bond(proxy))
  return bond_pairs

def _get_simple_bond (proxy) :
  import scitbx.array_family # import dependency
  if isinstance(proxy, distance_proxy) :
    pair = proxy.i_seqs[0], proxy.i_seqs[1]
  else :
    assert 0
  return pair

def _distance (sites_cart, i_seq, j_seq) :
  (x1, y1, z1) = sites_cart[i_seq]
  (x2, y2, z2) = sites_cart[j_seq]
  dist = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
  return dist

def filter_excessive_distances (proxies, sites_cart) :
  filtered_proxies = []
  for proxy in proxies :
    i_seq, j_seq = _get_simple_bond(proxy)
    distance = _distance(sites_cart, i_seq, j_seq)
    if (proxy.distance_cut > 0) and (distance > proxy.distance_cut) :
      continue
    filtered_proxies.append(proxy)
  return filtered_proxies

def as_pymol_dashes (proxies, pdb_hierarchy, filter=True, out=sys.stdout) :
  pdb_atoms = pdb_hierarchy.atoms()
  sites_cart = pdb_atoms.extract_xyz()
  if (filter) :
    proxies = filter_excessive_distances(proxies, sites_cart)
  for proxy in proxies :
    i_seq, j_seq = _get_simple_bond(proxy)
    atom1 = pdb_atoms[i_seq].fetch_labels()
    atom2 = pdb_atoms[j_seq].fetch_labels()
    base_sele = """chain "%s" and resi %s and name %s"""
    sele1 = base_sele % (atom1.chain_id, atom1.resseq, atom1.name)
    sele2 = base_sele % (atom2.chain_id, atom2.resseq, atom2.name)
    print >>out, "dist %s, %s" % (sele1, sele2)

def as_refmac_restraints (proxies, pdb_hierarchy, filter=True, out=sys.stdout,
    sigma=0.05) :
  pdb_atoms = pdb_hierarchy.atoms()
  sites_cart = pdb_atoms.extract_xyz()
  if (filter) :
    proxies = filter_excessive_distances(proxies, sites_cart)
  for proxy in proxies :
    i_seq, j_seq = _get_simple_bond(proxy)
    atom1 = pdb_atoms[i_seq].fetch_labels()
    atom2 = pdb_atoms[j_seq].fetch_labels()
    cmd = (("exte dist first chain %s residue %s atom %s " +
            "second chain %s residue %s atom %s value %.3f sigma %.2f") %
      (atom1.chain_id, atom1.resseq, atom1.name, atom2.chain_id,
       atom2.resseq, atom2.name, proxy.distance_ideal, 0.05))
    print >> out, cmd

def as_kinemage (proxies, pdb_hierarchy, filter=True, out=sys.stdout) :
  pdb_atoms = pdb_hierarchy.atoms()
  sites_cart = pdb_atoms.extract_xyz()
  print >> out, """\
@group {PHENIX H-bonds}
@subgroup {H-bond dots} dominant"""
  if (filter) :
    proxies = filter_excessive_distances(proxies, sites_cart)
  for proxy in proxies :
    i_seq, j_seq = _get_simple_bond(proxy)
    a = pdb_atoms[i_seq].xyz
    b = pdb_atoms[j_seq].xyz
    ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    print >> out, """@dotlist {Drawn dots} color= green"""
    for x in range(1, 12) :
      fac = float(x) / 12
      vec = (a[0] + (ab[0]*fac), a[1] + (ab[1]*fac), a[2] + (ab[2]*fac))
      if (x == 1) :
        print >> out, "{drawn} %.4f %.4f %.4f" % vec
      else :
        print >> out, "{''} %.4f %.4f %.4f" % vec

def find_implicit_hydrogen_bonds (pdb_hierarchy,
                                  xray_structure,
                                  params,
                                  log=None) :
  if (log is None) :
    log = null_out()
  donor_selection = "name N or (resname HOH and name O)"
  if (params.include_side_chains) :
    donor_selection += " or (resname ASN and name ND2) or "+ \
        "(resname GLN and name NE2) or (resname TRP and name NE1) or "+ \
        "(resname HIS and name NE2) or (resname LYS and name NZ) or "+ \
        "(resname ARG and name N*) or (resname TYR and name OH) or " +\
        "(resname SER and name OG) or (resname THR and name OG1) or "+\
        "(resname HOH and name O)"
  # What about waters?  No acceptor base, but could distances alone be used?
  acceptor_selection = "name O"
  if (params.include_side_chains) :
    acceptor_selection += " or (resname SER and name OG) or "+ \
        "(resname THR and name OG1) or (resname TYR and name OH) or "+ \
        "(resname ASN and name OD1) or (resname GLN and name OE1)"
  acceptor_base_selection = "name C"
  if (params.include_side_chains) :
    acceptor_base_selection += " or (resname SER and name CB) or "+ \
        "(resname THR and name CB) or (resname TYR and name CZ) or "+ \
        "(resname ASN and name CG) or (resname GLN and name CD)"
  selection_cache = pdb_hierarchy.atom_selection_cache()
  donors = selection_cache.selection(donor_selection)
  acceptors = selection_cache.selection(acceptor_selection)
  acceptor_bases = selection_cache.selection(acceptor_base_selection)
  restraint_type = params.restraint_type
  if (restraint_type == "simple") :
    build_proxies = build_simple_hbond_proxies()
  else :
    assert (restraint_type in ["Auto", "implicit"])
    build_proxies = build_implicit_hbond_proxies()
  if (len(donors) == 0) or (len(acceptors) == 0) or (len(acceptor_bases)==0) :
    return build_proxies # None
  pair_asu_table = xray_structure.pair_asu_table(
    distance_cutoff=params.distance_cut_n_o)
  pair_sym_table = pair_asu_table.extract_pair_sym_table(
    skip_j_seq_less_than_i_seq=False)
  atoms = pdb_hierarchy.atoms()
  unit_cell = xray_structure.unit_cell()
  sites_frac = xray_structure.sites_frac()
  sites_cart = xray_structure.sites_cart()
  for j_seq, pair_sym_dict in enumerate(pair_sym_table) :
    if (not acceptors[j_seq]) :
      continue
    site_j = sites_frac[j_seq]
    atom_j = atoms[j_seq]
    atom_labels_j = atom_j.fetch_labels()
    resid_j = atom_labels_j.resid()
    chain_j = atom_labels_j.chain_id
    acceptor_i_seq = j_seq
    donor_i_seq = acceptor_base_i_seq = None
    for k_seq, sym_ops in pair_sym_dict.items() :
      if (not donors[k_seq]) and (not acceptor_bases[k_seq]) :
        continue
      atom_k = atoms[k_seq]
      site_k = sites_frac[k_seq]
      atom_labels_k = atom_k.fetch_labels()
      resid_k = atom_labels_k.resid()
      chain_k = atom_labels_k.chain_id
      # skip donor-acceptor pairs in the same residue - for proteins this is
      # a safe assumption
      if (donors[k_seq]) and (resid_j == resid_k) and (chain_j == chain_k) :
        continue
      # FIXME need to include H-bonds between symmetry mates
      for sym_op in sym_ops :
        if (sym_op.is_unit_mx()) :
          if (acceptor_bases[k_seq]) :
            distance = unit_cell.distance(site_j, site_k)
            if (acceptor_base_i_seq is not None) :
              pass
              #print >> log, "    already have acceptor base for %s" % \
              #  atom_j.id_str()
              #break
            elif (distance < 2.0) :
              #print >> log, "  bonded to %s:" % atom_labels_j.id_str()
              #print >> log, "    %s" % atom_labels_k.id_str()
              acceptor_base_i_seq = k_seq
          else :
            donor_i_seq = k_seq
    if (not None in [donor_i_seq, acceptor_base_i_seq]) :
      i_seqs = [donor_i_seq, acceptor_i_seq, acceptor_base_i_seq]
      angle = unit_cell.angle(
        sites_frac[i_seqs[0]], sites_frac[i_seqs[1]], sites_frac[i_seqs[2]])
      if (angle >= params.implicit.theta_cut) :
        if (restraint_type == "simple") :
          build_proxies.add_proxy(
            i_seqs=[donor_i_seq, acceptor_i_seq],
            distance_ideal=params.distance_ideal_n_o,
            distance_cut=params.distance_cut_n_o,
            weight=params.restraints_weight/(params.simple.sigma**2),
            slack=params.simple.slack)
        else :
          build_proxies.add_proxy(
            i_seqs=[donor_i_seq, acceptor_i_seq, acceptor_base_i_seq],
            distance_ideal=params.distance_ideal_n_o,
            distance_cut=params.distance_cut_n_o,
            theta_low=params.implicit.theta_low,
            theta_high=params.implicit.theta_high,
            weight=params.restraints_weight)
        build_proxies.add_nonbonded_exclusion(donor_i_seq, acceptor_i_seq)
  print >> log, ""
  print >> log, "  %d hydrogen bond restraints generated." % \
    len(build_proxies.proxies)
  return build_proxies

# See Fabiola et al. (2002) Protein Sci. 11:1415-23 for motivation
class optimize_hbond_restraints (object) :
  def __init__ (self, model, fmodels, monitors, target_weights, params, nproc,
      log) :
    adopt_init_args(self, locals())
    assert (params.optimize_hbonds or params.optimize_hbonds_thorough)
    grm = model.restraints_manager.geometry
    assert ((grm is not None) and (grm.generic_restraints_manager is not None))
    self.rm = grm.generic_restraints_manager
    from mmtbx import utils
    utils.assert_xray_structures_equal(
      x1 = self.fmodels.fmodel_xray().xray_structure,
      x2 = self.model.xray_structure)
    self.save_scatterers = self.fmodels.fmodel_xray().xray_structure.\
        deep_copy_scatterers().scatterers()
    make_header("Optimizing H-bond refinement", out=self.log)
    fmodels.create_target_functors()
    assert approx_equal(self.fmodels.fmodel_xray().target_w(),
      self.fmodels.target_functor_result_xray(
        compute_gradients=False).target_work())
    trial_settings = []
    weights = [0.25, 0.5, 1.0, 1.5, 2.0, 4.0]
    distances = [2.8, 2.85, 2.9, 2.95, 3.0]
    thetas_high = [150, 155, 160, 165, 170]
    thetas_low = [95, 100, 105, 110, 115]
    cutoffs = [3.4, 3.5, 3.6, 3.7]
    # first grid search: distance_ideal_n_o vs. restraints_weight
    # XXX the paper does this a little differently:
    #  1. optimize weight (separately for MC and SC bonds)
    #  2. optimize theta_high and theta_low
    #  3. optimize ideal distance, cutoff, and angle cutoff (not used here)
    #  4. optimize the weight a final time with parameters from (2) and (3)
    for w in weights :
      for rij in distances :
        trial_settings.append(group_args(
          distance_ideal_n_o=rij,
          restraints_weight=w,
          distance_cut_n_o=self.params.distance_cut_n_o,
          theta_high=self.params.implicit.theta_high,
          theta_low=self.params.implicit.theta_low))
    stdout_and_results = easy_mp.pool_map(
      fixed_func=self._trial_minimization,
      args=trial_settings,
      buffer_stdout_stderr=True)
    results = [ r for (so, r) in stdout_and_results ]
    print >> self.log, ps("-" * 70)
    print >> self.log, ps(" Weight vs. N-O distance:")
    self.process_results(results, (not params.optimize_hbonds_thorough))
    if (params.optimize_hbonds_thorough) :
      # second grid search: theta_high vs. theta_low
      self.reset_sites()
      print >> self.log, ps("")
      trial_settings = []
      for theta1 in thetas_high :
        for theta2 in thetas_low :
          for distance in cutoffs :
            trial_settings.append(group_args(
              distance_ideal_n_o=self.params.distance_ideal_n_o,
              restraints_weight=self.params.restraints_weight,
              distance_cut_n_o=distance,
              theta_high=theta1,
              theta_low=theta2))
      stdout_and_results = easy_mp.pool_map(
        fixed_func=self._trial_minimization,
        args=trial_settings,
        buffer_stdout_stderr=True)
      results = [ r for (so, r) in stdout_and_results ]
      print >> self.log, ps(" Theta high vs. theta low vs. N-O cutoff:")
      self.process_results(results)
    print >> self.log, ps("-" * 70)

  def process_results (self, results, print_final_params=True) :
    results_final = []
    best_result = None
    best_r_gap = sys.maxint
    for result in results :
      if (result is not None) :
        results_final.append(result)
        if (result.r_gap < best_r_gap) :
          best_r_gap = result.r_gap
          best_result = result
    if (len(results_final) != 0) :
      print >> self.log, ps("      wt  dist   cut theta1 theta2 r_work r_free r_gap")
      print >> self.log, ps("   " + "-" * 53)
      for result in results_final :
        if (result is best_result) :
          result.show(out=self.log, mark="<<<")
        else :
          result.show(out=self.log)
      self.copy_settings(best_result.settings)
      self.fmodels.fmodel_xray().xray_structure.set_sites_frac(
        results[-1].sites_frac)
      self.fmodels.update_xray_structure(
        xray_structure = self.fmodels.fmodel_xray().xray_structure,
        update_f_calc  = True,
        update_f_mask  = True)
      self.model.xray_structure = self.fmodels.fmodel_xray().xray_structure
      print >> self.log, ps("   " + "-" * 53)
      if (print_final_params) :
        print >> self.log, ps("    Optimized H-bond parameters:")
        print >> self.log, ps("            Weight (restraints_weight) = %4.2f"%
          self.params.restraints_weight)
        print >> self.log, ps("     N-O distance (distance_ideal_n_o) = %5.3f"%
          self.params.distance_ideal_n_o)
        print >> self.log, ps("             Cutoff (distance_cut_n_o) = %5.3f"%
          self.params.distance_cut_n_o)
        print >> self.log, ps("               High angle (theta_high) = %3.0d"%
          self.params.implicit.theta_high)
        print >> self.log, ps("                 Low angle (theta_low) = %3.0d"%
          self.params.implicit.theta_low)
    else :
      print "  No hydrogen bonds found (works for proteins only at present)."

  def copy_settings (self, settings) :
    self.params.distance_ideal_n_o = settings.distance_ideal_n_o
    self.params.restraints_weight = settings.restraints_weight
    self.params.distance_cut_n_o = settings.distance_cut_n_o
    self.params.implicit.theta_high = settings.theta_high
    self.params.implicit.theta_low = settings.theta_low

  def reset_sites (self) :
    self.fmodels.fmodel_xray().xray_structure.replace_scatterers(
      self.save_scatterers.deep_copy())
    self.fmodels.update_xray_structure(
      xray_structure = self.fmodels.fmodel_xray().xray_structure,
      update_f_calc  = True)

  def _trial_minimization (self, settings) :
    self.copy_settings(settings)
    self.rm.update_hydrogen_bonds(
      pdb_hierarchy=self.model.pdb_hierarchy(),
      xray_structure=self.model.xray_structure,
      params=self.params,
      log=null_out())
    n_hbonds = self.rm.get_n_hbonds()
    if (n_hbonds == 0) :
      return None
    self.reset_sites()
    new_sites = self.minimize()
    return trial_result(
      r_work=self.fmodels.fmodel_xray().r_work(),
      r_free=self.fmodels.fmodel_xray().r_free(),
      settings=settings,
      sites_frac=new_sites)

  def minimize(self):
    import mmtbx.refinement.minimization
    from mmtbx import utils
    import scitbx.lbfgs
    utils.assert_xray_structures_equal(
      x1 = self.fmodels.fmodel_xray().xray_structure,
      x2 = self.model.xray_structure)
    self.model.set_refine_individual_sites()
    minimized = mmtbx.refinement.minimization.lbfgs(
      restraints_manager       = self.model.restraints_manager,
      fmodels                  = self.fmodels,
      model                    = self.model,
      refine_xyz               = True,
      target_weights           = self.target_weights)
    self.model.xray_structure = self.fmodels.fmodel_xray().xray_structure
    assert minimized.xray_structure is self.model.xray_structure
    utils.assert_xray_structures_equal(
      x1 = minimized.xray_structure,
      x2 = self.model.xray_structure)
    return minimized.xray_structure.sites_frac()

class trial_result (object) :
  def __init__ (self, r_work, r_free, settings, sites_frac) :
    adopt_init_args(self, locals())
    self.r_gap = r_free - r_work

  def show (self, out=None, mark="") :
    if (out is None) : out = null_out()
    fs = "    %4.2f %4.3f %4.3f %6.1f %6.1f %6.4f %6.4f %6.4f %3s"
    print >> out, ps(fs % (self.settings.restraints_weight,
      self.settings.distance_ideal_n_o, self.settings.distance_cut_n_o,
      self.settings.theta_high, self.settings.theta_low, self.r_work,
      self.r_free, self.r_gap, mark))
