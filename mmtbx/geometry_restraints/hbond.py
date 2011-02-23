
from __future__ import division
import libtbx.phil
from libtbx import adopt_init_args, group_args
from math import sqrt
import sys

master_phil = libtbx.phil.parse("""
  restraints_weight = 1.0
    .type = float
  falloff_distance = 0.05
    .type = float
  exclude_nonbonded = True
    .type = bool
  implicit
    .short_caption = Implicit hydrogens
    .help = Based on H-bond potential for CNS by Chapman lab
  {
    distance_ideal = 2.9
      .type = float
    distance_cut = 3.5
      .type = float
    theta_high = 155
      .type = float
    theta_low = 115
      .type = float
  }
  explicit
    .short_caption = Explicit hydrogens
    .help = Similar to Rosetta H-bond energy (Kortemme & Baker)
  {
    distance_ideal = 1.975
      .type = float
    distance_cut = None
      .type = float
    distance_sigma = 0.05
      .type = float
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
  simple
    .short_caption = Simple distance-based potential
    .help = Pseudo-bond restraints
  {
    distance_ideal_h_o = 1.975
      .type = float
    distance_cut_h_o = 2.5
      .type = float
    distance_ideal_n_o = 2.9
      .type = float
    distance_cut_n_o = 3.5
      .type = float
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
                          use_finite_differences=True) :
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
  elif (type(proxies).__name__ == "shared_h_bond_implicit_proxy") :
    sum = ext.h_bond_implicit_residual_sum(
      sites_cart=sites_cart,
      proxies=proxies,
      gradient_array=gradient_array,
      falloff_distance=falloff_distance,
      use_finite_differences=use_finite_differences)
  else :
    assert 0
  return sum

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
  if (type(proxy).__name__ == "h_bond_simple_proxy") :
    pair = proxy.i_seqs
  elif isinstance(proxy, implicit_proxy) :
    pair = proxy.i_seqs[0], proxy.i_seqs[1]
  elif isinstance(proxy, explicit_proxy) :
    pair = proxy.i_seqs[1], proxy.i_seqs[2]
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
    if (distance > proxy.distance_cut) :
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
    atom1 = atoms[i_seq].fetch_labels()
    atom2 = atoms[i_seq].fetch_labels()
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
    atom2 = pdb_atoms[i_seq].fetch_labels()
    cmd = (("exte dist first chain %s residue %s atom %s " +
            "second chain %s residue %s atom %s value %.3f sigma %.2f") %
      (atom1.chain_id, atom1.resseq, atom1.name, atom2.chain_id,
       atom2.resseq, atom2.name, bond.distance_ideal, 0.05))
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
