
import libtbx.phil
from libtbx.utils import Sorry
import sys

master_phil = libtbx.phil.parse("""
max_rmsd = None
  .type = float
sigma = 5.0
  .type = float
slack = 0
  .type = float
""")


rotamer_db = None
def load_sidechain_properties () :
  global rotamer_db
  if (rotamer_db is None) :
    from mmtbx.rotamer import sidechain_angles
    rotamer_db = sidechain_angles.SidechainAngles(False)

# TODO unify this with Ramachandran proxy generation?
def extract_proxies (pdb_hierarchy,
                     atom_selection=None,
                     log=None) :
  if (log is None) : log = sys.stdout
  from mmtbx.rotamer import sidechain_angles
  import mmtbx.rotamer
  import boost.python
  ext = boost.python.import_ext("mmtbx_rotamer_restraints_ext")
  residues = sidechain_angles.collect_sidechain_chi_angles(pdb_hierarchy)
  proxies = ext.shared_rotamer_proxy()
  for residue_info in residues :
    residue_name = residue_info.residue_name
    if (residue_name == "MSE") :
      residue_name = "MET"
    n_angles = 0
    i_seqs = [0] * 16
    last_chi = 0
    for chi in residue_info.chis :
      if (chi.chi_id > (last_chi + 1)) :
        break
      for k in range(4) :
        i_seqs[(n_angles*4) + k] = chi.i_seqs[k]
      n_angles += 1
      last_chi = chi.chi_id
    proxy = ext.rotamer_proxy(
      residue_name=residue_name,
      n_angles=n_angles,
      chi1=i_seqs[0:4],
      chi2=i_seqs[4:8],
      chi3=i_seqs[8:12],
      chi4=i_seqs[12:16])
    proxies.append(proxy)
  print >> log, ""
  print >> log, "  %d rotamer restraints generated." % proxies.size()
  return proxies

class manager (object) :
  def __init__ (self, pdb_hierarchy, params=None, log=None) :
    if (params is None) :
      params = master_phil.fetch().extract()
    if (params.sigma is None) or (params.slack is None) :
      raise Sorry("Rotamer restraint sigma and slack must be decimal numbers.")
    self.params = params
    self.pdb_hierarchy = pdb_hierarchy
    self.rotamer_proxies = extract_proxies(pdb_hierarchy, log=log)
    self._proxy_indices = None

  def cross_reference_proxies (self, dihedral_proxies, log=None) :
    if (log is None) : log = sys.stdout
    self._proxy_indices = []
    self._invert_signs = []
    have_proxies = set([])
    atoms = self.pdb_hierarchy.atoms()
    for proxy in self.rotamer_proxies :
      if (proxy.n_angles == 0) :
        self._proxy_indices.append(None)
        continue
      chi_proxies = [None] * 4
      invert = [False] * 4
      for j_seq, dihedral_proxy in enumerate(dihedral_proxies) :
        if (j_seq in have_proxies) :
          continue
        chi_id = proxy.find_dihedral_proxy(dihedral_proxy)
        if (chi_id != 0) :
          chi_proxies[abs(chi_id)-1] = j_seq
          have_proxies.add(j_seq)
          if (chi_id < 0) :
            invert[abs(chi_id)-1] = True
          if (not None in chi_proxies) :
            break
      n_missing = chi_proxies.count(None) - (4 - proxy.n_angles)
      if (n_missing > 0) :
        first_i_seq = proxy.first_i_seq()
        print >> log, "    warning: can't find %d dihedral proxies for %s" % (
          n_missing, atoms[first_i_seq].id_str()[9:-1])
      self._proxy_indices.append(chi_proxies)
      self._invert_signs.append(invert)

  def update_dihedral_proxies (self,
                               sites_cart,
                               dihedral_proxies,
                               verbose=True,
                               log=None) :
    n_updates = 0
    if (log is None) :
      log = sys.stdout
    params = self.params
    load_sidechain_properties()
    pdb_atoms = self.pdb_hierarchy.atoms()
    if (self._proxy_indices is None) :
      self.cross_reference_proxies(dihedral_proxies)
    assert (len(self.rotamer_proxies) == len(self._proxy_indices))
    print >> log, " Updating dihedral proxies from rotamer library"
    weight = 1 / params.sigma**2
    for k, proxy in enumerate(self.rotamer_proxies) :
      if (proxy.n_angles == 0) or (self._proxy_indices[k].count(None) == 4) :
        continue
      rotamers = rotamer_db.get_rotamers(proxy.residue_name)
      if (rotamers is None) :
        print >> log, "  skipping residue %s" % proxy.residue_name
        continue
      best_rmsd = sys.maxint
      best_angles = None
      best_rotamer = None
      for rotamer, angles in rotamers.iteritems() :
        angles_ = angles
        if (len(angles) < 4) :
          angles_.extend([0] * (4 - len(angles)))
        rmsd = proxy.get_rotamer_rmsd(
          angles=angles_,
          sites_cart=sites_cart)
        if (rmsd < best_rmsd) :
          best_rmsd = rmsd
          best_angles = angles
          best_rotamer = rotamer
      if (best_angles is None) :
        print >> log, "  can't find angles for residue %s" % resname
      elif (params.max_rmsd is not None) and (best_rmsd > params.max_rmsd) :
        print >> log, "  lowest RMSD exceeds max"
      else :
        print_rotamer = verbose
        for chi_id, j_seq in enumerate(self._proxy_indices[k]) :
          if (j_seq is None) :
            continue
          dihedral_proxy = dihedral_proxies[j_seq]
          if (print_rotamer) :
            first_i_seq = dihedral_proxy.i_seqs[0]
            atom = pdb_atoms[first_i_seq]
            print >> log, "  %s : %s (rmsd=%.3f)" % (atom.id_str()[9:-1],
              best_rotamer, best_rmsd)
            print_rotamer = False
          if (self._invert_signs[k][chi_id]) :
            dihedral_proxy.angle_ideal = - angles[chi_id]
          else :
            dihedral_proxy.angle_ideal = angles[chi_id]
          dihedral_proxy.weight = weight
          n_updates += 1
    print >> log, " %d proxies modified." % n_updates
    return n_updates
