
from __future__ import division
from __future__ import print_function
from mmtbx.building import generate_sidechain_clusters
from libtbx.str_utils import format_value
from libtbx.math_utils import ifloor, iceil
from libtbx import adopt_init_args

torsion_search_params = """
chi_increment_degrees = 10 #5
  .type = int
  .help = Number of degrees per rotation of each chi angle
backrub = True # ignore & use direct input arguments for debugging...
  .type = bool
  .help = Sample backrub motion
backrub_range = 20
  .type = float
  .help = This specifies the number of degrees to sample on either direction \
    from the current position
backrub_increment = 5
  .type = float
  .help = Number of degrees per backrub motion
shear = True # ignore & use direct input arguments for debugging...
  .type = bool
  .help = Sample shear motion
shear_range = 8
  .type = float
  .help = This specifies the number of degrees to sample on either direction \
    from the current position
shear_increment = 2
  .type = float
  .help = Number of degrees per shear motion
"""

pseudosymmetric_sidechains = ["PHE","TYR","HIS","ASP","ASN","GLU","GLN"]

class conformation(object):
  """
  Stores information about a specific conformation of a residue (and optionally
  the adjacent backrub atoms) or of a pair of residues (and optionally the
  adjacent shear and/or backrub atoms), and evalutes fit to density.
  """
  def __init__(self,
      residue,
      prev_residue,
      next_residue,
      next_next_residue,
      backrub,
      shear,
      chis,
      rotamer,
      rmsd,
      sites_cart,
      sites_selection,
      sidechain_selection):
    adopt_init_args(self, locals())
    self.mean_fofc = self.mean_2fofc = self.translation = None

  def show_summary(self, out=None):
    if (out is None) : out = sys.stdout
    chi_angles = None
    if (self.chis is not None):
      chi_angles = ",".join([ "%.1f" % x for x in self.chis ])
    density_str = ""
    if (self.mean_fofc is not None):
      density_str = "%6.2f  %6.2f  " % (self.mean_2fofc, self.mean_fofc)
    translation_str = ""
    if (self.translation is not None):
      translation_str = "%4.2f,%4.2f,%4.2f  " % tuple(self.translation)
    print("""%12s  %6s  %6s  %6s  %6.3f  %s%s%s""" % (
      self.residue.id_str(), format_value("%6.1f", self.backrub),
      format_value("%6.1f", self.shear),
      self.rotamer, self.rmsd, translation_str, density_str, chi_angles), file=out)

  def set_translation(self, translation):
    self.translation = translation
    return self

  def translate_sites(self, sites_new, translation):
    sites_cart_new = self.sites_cart.deep_copy()
    if (translation is not None):
      sites_cart_new.set_selected(self.sites_selection, sites_new)
    return conformation(
      residue=self.residue,
      prev_residue=self.prev_residue,
      next_residue=self.next_residue,
      next_next_residue=self.next_next_residue,
      backrub=self.backrub,
      shear=self.shear,
      chis=self.chis,
      rotamer=self.rotamer,
      rmsd=self.rmsd,
      sites_cart=sites_cart_new,
      sites_selection=self.sites_selection,
      sidechain_selection=self.sidechain_selection).set_translation(translation)

  def sites_selected(self):
    return self.sites_cart.select(self.sites_selection)

  def set_atom_sites(self, pdb_atoms):
    pdb_atoms.set_xyz(self.sites_cart)

  def as_pdb_fragment(self, pdb_atoms):
    conf_atoms = pdb_atoms.select(self.sites_selection)
    return fragment(conf_atoms)

  # FIXME I really need to handle this better
  def update_rmsd(self, sites_start):
    sites_sc = self.sites_cart.select(self.sidechain_selection)
    sites_sc_start = sites_start.select(self.sidechain_selection)
    self.rmsd = sites_sc.rms_difference(sites_sc_start)
    return self

class fragment(object):
  """
  Pseudo-residue object for an arbitrary collection of atoms.
  """
  def __init__(self, atoms):
    self._atoms = atoms
    self._i_seqs = atoms.extract_i_seq()
    assert (not self._i_seqs.all_eq(0))

  def atoms(self):
    return self._atoms

def generate_single_residue_confs(
      atom_group,
      sites_cart,
      mon_lib_srv,
      params,
      prev_residue=None,
      next_residue=None,
      next_next_residue=None,
      backrub=False,
      shear=False):
  confs = []
  confgen = residue_conformation_generator(
    residue=atom_group,
    sites_cart=sites_cart,
    mon_lib_srv=mon_lib_srv,
    params=params,
    prev_residue=prev_residue,
    next_residue=next_residue,
    next_next_residue=next_next_residue)
  if (confgen.has_motions()):
    for conf in confgen(backrub=backrub, shear=shear):
      yield conf

def torsion_search_nested(
      clusters,
      sites_cart,
      last_chi_symmetric=None,
      increment_degrees=10):
  """
  Iterate over all possible sidechain Chi angle combinations.
  """
  from scitbx.array_family import flex
  from scitbx.matrix import rotate_point_around_axis
  n_angles = len(clusters)
  assert (n_angles >= 1)
  angle_range = 180
  r1 = [ifloor(-angle_range/increment_degrees)] * n_angles
  r2 = [iceil(angle_range/increment_degrees)] * n_angles
  if (last_chi_symmetric):
    r1[-1] = ifloor(-90/increment_degrees)
    r2[-1] = iceil(90/increment_degrees)
  nested_loop = flex.nested_loop(begin=r1, end=r2, open_range=False)
  selection = clusters[0].atoms_to_rotate
  for angles in nested_loop:
    xyz_moved = sites_cart.deep_copy()
    for i, angle_fraction in enumerate(angles):
      cl = clusters[i]
      for atom in cl.atoms_to_rotate:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = xyz_moved[cl.axis[0]],
          axis_point_2 = xyz_moved[cl.axis[1]],
          point        = xyz_moved[atom],
          angle        = angle_fraction*increment_degrees,
          deg=True)
        xyz_moved[atom] = new_xyz
    yield xyz_moved

def backrub_rotate(
    sites,
    i_seqs_primary,
    primary_axis,
    primary_angle,
    secondary_axis1=None,
    secondary_axis2=None,
    secondary_angle1=None,
    secondary_angle2=None):
  """
  Performs a "backrub" move of magnitude primary_angle.
  Rotates the atoms in sites indexed by i_seqs_primary around primary_axis.
  Directly modifies sites -- does not return anything.
  """
  from scitbx.array_family import flex
  from scitbx import matrix
  assert isinstance(primary_angle, float)
  assert (len(primary_axis) == 2)
  sites_new = flex.vec3_double()
  for i_seq in i_seqs_primary :
    xyz = matrix.rotate_point_around_axis(
      axis_point_1=primary_axis[0],
      axis_point_2=primary_axis[1],
      point=sites[i_seq],
      angle=primary_angle,
      deg=True)
    sites_new.append(xyz)
  sites.set_selected(i_seqs_primary, sites_new)

# TODO
# shear_rotate that takes just 1 residue instead of all the i_seqs stuff?
# be smarter about which primary2 rot'ns to try (should correlate to primary1)
def shear_rotate(
    sites,
    i_seqs_primary1,
    i_seqs_primary2,
    i_seqs_middle,
    i_seqs_ca234,
    primary_axis1,
    primary_angle,
    secondary_axis1=None,
    secondary_axis2=None,
    secondary_axis3=None):
  """
  Performs a "shear" move.
  First, rotates the atoms in sites indexed by i_seqs_primary1 within the
  Calpha 1-2-3 plane (around primary_axis1).
  Then, rotates the atoms in sites indexed by i_seqs_primary2 within the new
  Calpha 2'-3'-4 plane.
  Finally, moves the atoms in sites indexed by i_seqs_middle to plug the middle
  peptide into the Calpha 2'-3' gap.
  Directly modifies sites -- does not return anything.
  """
  from scitbx.array_family import flex
  from scitbx import matrix
  from sys import float_info
  assert isinstance(primary_angle, float)
  assert (len(primary_axis1) == 2)
  if ((i_seqs_primary1 is None) or
      (i_seqs_primary2 is None) or
      (i_seqs_middle is None)):
    print("TODO set shear_rotate i_seqs automatically if not provided")
    return
  # prep for second primary rotation axis
  sites_start = sites.deep_copy()
  orig_dist = abs(matrix.col(sites_start[i_seqs_ca234[1]]) - \
                  matrix.col(sites_start[i_seqs_ca234[0]]))
  # this must be happening b/c previous shears didn't actually move atoms past
  # ca2 for some reason...
  assert (orig_dist > 3.7 and orig_dist < 3.9), \
    "original ca2-ca3 distance for shear looks fishy: %.3f" % orig_dist
  # do first primary rotation
  sites_new = flex.vec3_double()
  for i_seq in i_seqs_primary1 :
    xyz = matrix.rotate_point_around_axis(
      axis_point_1=primary_axis1[0],
      axis_point_2=primary_axis1[1],
      point=sites[i_seq],
      angle=primary_angle,
      deg=True)
    sites_new.append(xyz)
  sites.set_selected(i_seqs_primary1, sites_new)
  # prep for second primary rotation axis (cont'd)
  ca23 = matrix.col(sites[i_seqs_ca234[1]]) - matrix.col(sites[i_seqs_ca234[0]])
  ca34 = matrix.col(sites[i_seqs_ca234[2]]) - matrix.col(sites[i_seqs_ca234[1]])
  normal2 = ca23.cross(ca34)
  ca4_normal2 = matrix.col(sites_start[i_seqs_ca234[2]]) + matrix.col(normal2)
  primary_axis2 = [sites_start[i_seqs_ca234[2]], ca4_normal2]
  # do second primary rotation:
  # find magnitude that best preserves original ca2-ca3 distance
  best_diff_orig_dist = float_info.max
  best_sites_new = None
  best_theta = None
  theta = -1.5 * abs(primary_angle)
  theta_increment = 0.1
  while (theta < 1.5 * abs(primary_angle)):
    ca3_new = matrix.rotate_point_around_axis(
      axis_point_1=primary_axis2[0],
      axis_point_2=primary_axis2[1],
      point=sites[i_seqs_ca234[1]],
      angle=theta,
      deg=True)
    dist = abs(matrix.col(ca3_new) - matrix.col(sites[i_seqs_ca234[0]]))
    diff_orig_dist = abs(dist - orig_dist)
    if (diff_orig_dist < best_diff_orig_dist):
      best_theta = theta
      best_diff_orig_dist = diff_orig_dist
      # this rotation is best so far,
      # so do it to the *entire* third peptide and save the result
      sites_new = flex.vec3_double()
      for i_seq in i_seqs_primary2 :
        xyz = matrix.rotate_point_around_axis(
          axis_point_1=primary_axis2[0],
          axis_point_2=primary_axis2[1],
          point=sites[i_seq],
          angle=theta,
          deg=True)
        sites_new.append(xyz)
      best_sites_new = sites_new
    theta += theta_increment
  sites.set_selected(i_seqs_primary2, best_sites_new)
  final_dist = abs(matrix.col(sites[i_seqs_ca234[1]]) - \
                   matrix.col(sites[i_seqs_ca234[0]]))
  assert (abs(final_dist - orig_dist) < 0.1), "failed to close ca2-ca3 gap!"
  # plug in middle peptide by averaging the positions that result from
  # applying the two primary rotations to the middle peptide
  sites_new = flex.vec3_double()
  for i_seq in i_seqs_middle :
    xyz1 = matrix.rotate_point_around_axis(
      axis_point_1=primary_axis1[0],
      axis_point_2=primary_axis1[1],
      point=sites[i_seq],
      angle=primary_angle,
      deg=True)
    xyz2 = matrix.rotate_point_around_axis(
      axis_point_1=primary_axis2[0],
      axis_point_2=primary_axis2[1],
      point=sites[i_seq],
      angle=best_theta,
      deg=True)
    xyz_avg = (matrix.col(xyz1) + matrix.col(xyz2)) / 2
    sites_new.append(xyz_avg)
  sites.set_selected(i_seqs_middle, sites_new)

class residue_conformation_generator(object):
  """
  Class for sampling all rotameric conformations of a residue, including the
  backrub and shear motions if possible.
  """
  def __init__(self,
      residue,
      sites_cart,
      mon_lib_srv,
      params,
      prev_residue=None,
      next_residue=None,
      next_next_residue=None,
      evaluate_backbone_callback=None):
    adopt_init_args(self, locals())
    from mmtbx.rotamer import rotamer_eval
    from mmtbx.rotamer import ramachandran_eval
    import iotbx.pdb
    from scitbx.array_family import flex
    get_class = iotbx.pdb.common_residue_names_get_class
    assert get_class(residue.resname) == "common_amino_acid", residue.resname
    self.rotamer_scorer = rotamer_eval.RotamerEval(data_version="8000")
    self.ramachandran_scorer = ramachandran_eval.RamachandranEval()
    self.sidechain_clusters = generate_sidechain_clusters(
      residue=residue,
      mon_lib_srv=mon_lib_srv)
    self.sites_start = sites_cart.deep_copy()
    self.i_seqs_residue = residue.atoms().extract_i_seq()
    self.i_seqs_sidechain = flex.size_t()
    for atom in self.residue.atoms():
      if (not atom.name.strip() in ["C","N","H","CA","CB","self"]):
        self.i_seqs_sidechain.append(atom.i_seq)
    self.i_seqs_primary = flex.size_t()
    if (not None in [prev_residue, next_residue]):
      self.set_up_backrub()
      if (next_next_residue is not None):
        self.set_up_shear()
        for i_seq in self.shear_i_seqs_primary1 :
          self.i_seqs_primary.append(i_seq)
        for i_seq in self.shear_i_seqs_primary2 :
          self.i_seqs_primary.append(i_seq)
        for i_seq in self.shear_i_seqs_middle :
          self.i_seqs_primary.append(i_seq)
      else:
        for i_seq in self.backrub_i_seqs : self.i_seqs_primary.append(i_seq)
    else :
      self.i_seqs_primary = self.i_seqs_residue

  def has_rotamers(self):
    return (len(self.sidechain_clusters) > 0)

  def has_backrub(self):
    return ((self.params.backrub) and
      (not None in [self.prev_residue, self.next_residue]))

  def has_shear(self):
    return ((self.params.shear) and
      (not None in [self.prev_residue, self.next_residue, self.next_next_residue]))

  def has_motions(self):
    return ((self.has_rotamers()) or (self.has_backrub()) or (self.has_shear()))

  def set_up_backrub(self):
    from scitbx.array_family import flex
    self.backrub_primary_axis = []
    self.backrub_secondary_axis1 = []
    self.backrub_secondary_axis2 = []
    self.backrub_i_seqs = flex.size_t()
    self.backrub_i_seqs_secondary1 = flex.size_t()
    self.backrub_i_seqs_secondary2 = flex.size_t()
    for atom in self.prev_residue.atoms():
      if (atom.name == " CA "):
        self.backrub_primary_axis.append(atom.xyz)
      elif (atom.name.strip() in ["C","self"]):
        self.backrub_i_seqs.append(atom.i_seq)
    for atom in self.residue.atoms():
      self.backrub_i_seqs.append(atom.i_seq)
      if (not atom.name.strip() in ["C","N","H","CA","CB","self"]):
        if (not atom.i_seq in self.i_seqs_sidechain):
          self.i_seqs_sidechain.append(atom.i_seq)
    for atom in self.next_residue.atoms():
      if (atom.name == " CA "):
        self.backrub_primary_axis.append(atom.xyz)
      elif (atom.name.strip() in ["N", "H"]):
        self.backrub_i_seqs.append(atom.i_seq)
    assert (len(self.backrub_primary_axis) == 2)

  def set_up_shear(self):
    from scitbx.array_family import flex
    from scitbx import matrix
    self.shear_primary_axis1 = []
    self.shear_secondary_axis1 = []
    self.shear_secondary_axis2 = []
    self.shear_secondary_axis3 = []
    self.shear_i_seqs_primary1 = flex.size_t()
    self.shear_i_seqs_primary2 = flex.size_t()
    self.shear_i_seqs = flex.size_t() # union of primary1, primary2, middle
    self.shear_i_seqs_middle = flex.size_t()
    self.shear_i_seqs_secondary1 = flex.size_t()
    self.shear_i_seqs_secondary2 = flex.size_t()
    self.shear_i_seqs_secondary3 = flex.size_t()
    self.i_seqs_ca234 = flex.size_t()
    ca123_xyz = []
    for atom in self.prev_residue.atoms():
      if (atom.name == " CA "):
        self.shear_primary_axis1.append(atom.xyz)
        ca123_xyz.append(atom.xyz)
      if (atom.name.strip() in ["C","self"]):
        self.shear_i_seqs_primary1.append(atom.i_seq)
    for atom in self.residue.atoms():
      if (atom.name == " CA "):
        ca123_xyz.append(atom.xyz)
      if (not atom.name.strip() in ["C","self"]):
        self.shear_i_seqs_primary1.append(atom.i_seq)
      if (not atom.name.strip() in ["C","N","H","CA","CB","self"]):
        if (not atom.i_seq in self.i_seqs_sidechain):
          self.i_seqs_sidechain.append(atom.i_seq)
      if (atom.name.strip() in ["C","self"]):
        self.shear_i_seqs_middle.append(atom.i_seq)
    for atom in self.next_residue.atoms():
      if (atom.name == " CA "):
        ca123_xyz.append(atom.xyz)
        ca12 = matrix.col(ca123_xyz[1]) - matrix.col(ca123_xyz[0])
        ca23 = matrix.col(ca123_xyz[2]) - matrix.col(ca123_xyz[1])
        normal1 = ca12.cross(ca23)
        ca1_normal1 = matrix.col(ca123_xyz[0]) + normal1
        self.shear_primary_axis1.append(ca1_normal1)
      if (atom.name.strip() in ["N","H"]):
        self.shear_i_seqs_middle.append(atom.i_seq)
      else :
        self.shear_i_seqs_primary2.append(atom.i_seq)
    for atom in self.next_next_residue.atoms():
      if (atom.name.strip() in ["N","H"]):
        self.shear_i_seqs_primary2.append(atom.i_seq)
    for i_seq in self.shear_i_seqs_primary1:
      self.shear_i_seqs.append(i_seq)
    for i_seq in self.shear_i_seqs_primary2:
      self.shear_i_seqs.append(i_seq)
    for i_seq in self.shear_i_seqs_middle:
      self.shear_i_seqs.append(i_seq)
    # prep for keeping track of ca2 and ca3
    for atom in self.residue.atoms():
      if (atom.name == " CA "):
        self.i_seqs_ca234.append(atom.i_seq)
    for atom in self.next_residue.atoms():
      if (atom.name == " CA "):
        self.i_seqs_ca234.append(atom.i_seq)
    for atom in self.next_next_residue.atoms():
      if (atom.name == " CA "):
        self.i_seqs_ca234.append(atom.i_seq)
    assert (len(self.shear_primary_axis1) == 2)

  def do_backrub_rotate(self, angle):
    backrub_rotate(
      sites=self.sites_cart,
      i_seqs_primary=self.backrub_i_seqs,
      primary_axis=self.backrub_primary_axis,
      primary_angle=angle)

  def do_shear_rotate(self, angle):
    shear_rotate(
      sites=self.sites_cart,
      i_seqs_primary1=self.shear_i_seqs_primary1,
      i_seqs_primary2=self.shear_i_seqs_primary2,
      i_seqs_middle=self.shear_i_seqs_middle,
      i_seqs_ca234=self.i_seqs_ca234,
      primary_axis1=self.shear_primary_axis1,
      primary_angle=angle)

  def sites_selected(self):
    return self.sites_cart.select(self.i_seqs_primary)

  def get_rotamer_info(self):
    """
    Retrieve the rotamer ID (or selfUTLIER) and list of chi angles.
    """
    chi_angles = rotamer_flag = None
    if (self.has_rotamers()):
       rotamer_flag = self.rotamer_scorer.evaluate_residue(self.residue)
       chi_angles = self.rotamer_scorer.chi_angles(self.residue)
    return rotamer_flag, chi_angles

  def all_valid_ramachandran(self):
    """
    Check whether phi/psi falls into the "favored" or "allowed" Ramachandran
    region for ALL of the involved residues (one to four, depending).
    """
    #rama_flag = self.ramachandran_scorer.evaluate(self.residue)
    #if (rama_flag == "selfUTLIER"):
    #  return False
    #if (self.prev_residue):
    #  rama_flag = self.ramachandran_scorer.evaluate(self.prev_residue)
    #  if (rama_flag == "selfUTLIER"):
    #    return False
    #if (self.next_residue):
    #  rama_flag = self.ramachandran_scorer.evaluate(self.next_residue)
    #  if (rama_flag == "selfUTLIER"):
    #    return False
    #if (self.next_next_residue is not None):
    #  rama_flag = self.ramachandran_scorer.evaluate(self.next_next_residue)
    #  if (rama_flag == "selfUTLIER"):
    #    return False
    return True

  # TODO:
  # return confs with 2nd sidechain's rotamers also enumerated for shears!
  # try different permutations of shears & backrubs?
  def __call__(self, backrub=False, shear=False):
    """
    Generator for rotameric conformations.
    Includes backrub and shear enumeration if requested.
    """
    from scitbx.array_family import flex
    id_str = self.residue.id_str()
    id_str_prev = id_str_next = id_str_next_next = None
    if (not None in [self.prev_residue, self.next_residue]):
      id_str_prev = self.prev_residue.id_str()
      id_str_next = self.next_residue.id_str()
      if (not self.next_next_residue is None):
        id_str_next_next = self.next_next_residue.id_str()
    ag_atoms = self.residue.atoms()
    ag_i_seqs = ag_atoms.extract_i_seq()
    assert (not ag_i_seqs.all_eq(0))
    ag_sites = self.sites_cart.select(ag_i_seqs).deep_copy()
    sites_start = self.sites_cart.select(self.i_seqs_residue).deep_copy()
    if (shear and backrub):
      if (not (self.has_shear() and self.has_backrub())):
        return
      theta_shear = - self.params.shear_range
      self.do_shear_rotate(theta_shear)
      while (theta_shear <= self.params.shear_range):
        if (not self.all_valid_ramachandran()):
          continue
        sites_preshear = self.sites_cart.select(self.i_seqs_primary).deep_copy()
        ag_sites_preshear = self.sites_cart.select(self.i_seqs_residue).deep_copy()
        ag_atoms.set_xyz(ag_sites_preshear)
        theta_backrub = - self.params.backrub_range
        self.do_backrub_rotate(theta_backrub)
        while (theta_backrub <= self.params.backrub_range):
          if (not self.all_valid_ramachandran()):
            continue
          sites_prebackrub = self.sites_cart.select(self.backrub_i_seqs).deep_copy()
          ag_sites_prebackrub = \
            self.sites_cart.select(self.i_seqs_residue).deep_copy()
          ag_atoms.set_xyz(ag_sites_prebackrub)
          for sites_new in self.iter_sidechain_confs():
            self.sites_cart.set_selected(self.i_seqs_residue, sites_new)
            ag_atoms.set_xyz(sites_new)
            rotamer_flag, chi_angles = self.get_rotamer_info()
            if (rotamer_flag == "selfUTLIER"):
              continue
            conf = conformation(
              residue=self.residue,
              prev_residue=self.prev_residue,
              next_residue=self.next_residue,
              next_next_residue=self.next_next_residue,
              shear=theta_shear,
              backrub=theta_backrub,
              chis=chi_angles,
              rmsd=self.sidechain_rmsd(),
              rotamer=rotamer_flag,
              sites_cart=self.sites_cart.deep_copy(),
              sites_selection=self.shear_i_seqs,
              sidechain_selection=self.i_seqs_sidechain)
            yield conf
          ag_atoms.set_xyz(ag_sites_prebackrub)
          self.sites_cart.set_selected(self.backrub_i_seqs, sites_prebackrub)
          self.do_backrub_rotate(self.params.backrub_increment)
          theta_backrub += self.params.backrub_increment
        ag_atoms.set_xyz(ag_sites_preshear)
        self.sites_cart.set_selected(self.i_seqs_primary, sites_preshear)
        self.do_shear_rotate(self.params.shear_increment)
        theta_shear += self.params.shear_increment
    elif (shear):
      if (not self.has_shear()):
        return
      theta_shear = - self.params.shear_range
      self.do_shear_rotate(theta_shear)
      while (theta_shear <= self.params.shear_range):
        if (not self.all_valid_ramachandran()):
          continue
        sites_preshear = self.sites_cart.select(self.i_seqs_primary).deep_copy()
        ag_sites_preshear = self.sites_cart.select(self.i_seqs_residue).deep_copy()
        ag_atoms.set_xyz(ag_sites_preshear)
        for sites_new in self.iter_sidechain_confs():
          self.sites_cart.set_selected(self.i_seqs_residue, sites_new)
          ag_atoms.set_xyz(sites_new)
          rotamer_flag, chi_angles = self.get_rotamer_info()
          if (rotamer_flag == "selfUTLIER"):
            continue
          conf = conformation(
            residue=self.residue,
            prev_residue=self.prev_residue,
            next_residue=self.next_residue,
            next_next_residue=self.next_next_residue,
            shear=theta_shear,
            backrub=None,
            chis=chi_angles,
            rmsd=self.sidechain_rmsd(),
            rotamer=rotamer_flag,
            sites_cart=self.sites_cart.deep_copy(),
            sites_selection=self.shear_i_seqs,
            sidechain_selection=self.i_seqs_sidechain)
          yield conf
        ag_atoms.set_xyz(ag_sites_preshear)
        self.sites_cart.set_selected(self.i_seqs_primary, sites_preshear)
        self.do_shear_rotate(self.params.shear_increment)
        theta_shear += self.params.shear_increment
    elif (backrub):
      if (not self.has_backrub()):
        return
      theta_backrub = - self.params.backrub_range
      self.do_backrub_rotate(theta_backrub)
      while (theta_backrub <= self.params.backrub_range):
        if (not self.all_valid_ramachandran()):
          continue
        sites_prebackrub = self.sites_cart.select(self.backrub_i_seqs).deep_copy()
        ag_sites_prebackrub = self.sites_cart.select(self.i_seqs_residue).deep_copy()
        ag_atoms.set_xyz(ag_sites_prebackrub)
        for sites_new in self.iter_sidechain_confs():
          self.sites_cart.set_selected(self.i_seqs_residue, sites_new)
          ag_atoms.set_xyz(sites_new)
          rotamer_flag, chi_angles = self.get_rotamer_info()
          if (rotamer_flag == "selfUTLIER"):
            continue
          conf = conformation(
            residue=self.residue,
            prev_residue=self.prev_residue,
            next_residue=self.next_residue,
            next_next_residue=None,
            shear=None,
            backrub=theta_backrub,
            chis=chi_angles,
            rmsd=self.sidechain_rmsd(),
            rotamer=rotamer_flag,
            sites_cart=self.sites_cart.deep_copy(),
            sites_selection=self.backrub_i_seqs,
            sidechain_selection=self.i_seqs_sidechain)
          yield conf
        ag_atoms.set_xyz(ag_sites_prebackrub)
        self.sites_cart.set_selected(self.backrub_i_seqs, sites_prebackrub)
        self.do_backrub_rotate(self.params.backrub_increment)
        theta_backrub += self.params.backrub_increment
    else :
      for sites_new in self.iter_sidechain_confs():
        self.sites_cart.set_selected(self.i_seqs_residue, sites_new)
        ag_atoms.set_xyz(sites_new)
        rotamer_flag, chi_angles = self.get_rotamer_info()
        if (rotamer_flag == "selfUTLIER"):
          continue
        conf = conformation(
          residue=self.residue,
          prev_residue=None,
          next_residue=None,
          next_next_residue=None,
          backrub=None,
          shear=None,
          chis=chi_angles,
          rmsd=self.sidechain_rmsd(),
          rotamer=rotamer_flag,
          sites_cart=self.sites_cart.deep_copy(),
          sites_selection=ag_i_seqs,
          sidechain_selection=self.i_seqs_sidechain)
        yield conf

  def iter_sidechain_confs(self):
    """
    Generate new conformations for this residue starting from its current
    coordinates, but do not change any variables in this class.
    """
    ag_sites = self.residue.atoms().extract_xyz()
    if (len(self.sidechain_clusters) == 0):
      yield ag_sites
    else :
      last_chi_symmetric = False
      if (self.residue.resname in pseudosymmetric_sidechains):
        last_chi_symmetric = True
      for sites_new in torsion_search_nested(
          clusters=self.sidechain_clusters,
          sites_cart=ag_sites,
          last_chi_symmetric=last_chi_symmetric,
          increment_degrees=self.params.chi_increment_degrees):
        yield sites_new

  def sidechain_rmsd(self):
    sites_start = self.sites_start.select(self.i_seqs_sidechain)
    sites_current = self.sites_cart.select(self.i_seqs_sidechain)
    return sites_start.rms_difference(sites_current)

  def evaluate_backbone_conformation(self):
    if (self.evaluate_backbone_callback is not None):
      return self.evaluate_backbone_callback(self.sites_cart.select(self.i_seqs_residue))
    return True
