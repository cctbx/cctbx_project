from __future__ import division
from scitbx.array_family import flex
from scitbx.math import superpose
import iotbx.pdb
import scitbx.matrix

class asu_ncs_converter(object):
  """
  Simple method of converting multiple NCS copy ASU content into one NCS copy
  and NCS transformations. Limitations:
    - Assumes each chain to be individual NCS copy.
    - Takes first NCS copy as the master copy.
    - Assume all copies are identical (number and order of atoms)
    - Anisotropic ADP are not handlesd.
  Further improvements: use atom selections to identify NCS copies.
  """

  def __init__(self, pdb_hierarchy, eps = 0.01, add_identity=True):
    self.pdb_hierarchy = pdb_hierarchy
    n_atoms_per_chain = flex.int()
    sites_cart_chain_0 = None
    self.rotation_matrices = []
    self.translation_vectors = []
    self.back_rotation_matrices = []
    self.back_translation_vectors = []
    self.ph_first_chain = None
    #
    for i_chain, chain in enumerate(pdb_hierarchy.chains()):
      n_atoms_per_chain.append(chain.atoms_size())
    #
    outlier_found = False
    if(n_atoms_per_chain.all_eq(n_atoms_per_chain[0])):
      for i_chain, chain in enumerate(pdb_hierarchy.chains()):
        if(chain.is_na() or chain.is_protein()):
          n_atoms_per_chain.append(chain.atoms_size())
          if(sites_cart_chain_0 is None and i_chain==0):
            sites_cart_chain_0 = chain.atoms().extract_xyz()
            sel = flex.size_t(xrange(sites_cart_chain_0.size()))
            self.ph_first_chain = pdb_hierarchy.select(sel)
            if(add_identity):
              um = scitbx.matrix.sqr((
                1,0,0,
                0,1,0,
                0,0,1))
              zv = scitbx.matrix.col((0, 0, 0))
              self.rotation_matrices.append(um)
              self.translation_vectors.append(zv)
              self.back_rotation_matrices.append(um)
              self.back_translation_vectors.append(zv)
          if(i_chain > 0):
            # first copy onto others
            lsq_fit_obj = superpose.least_squares_fit(
              reference_sites = sites_cart_chain_0,
              other_sites     = chain.atoms().extract_xyz())
            self.rotation_matrices.append(lsq_fit_obj.r.transpose())
            self.translation_vectors.append(lsq_fit_obj.t)
            d =  flex.sqrt((sites_cart_chain_0-
              lsq_fit_obj.other_sites_best_fit()).dot()).min_max_mean().as_tuple()
            if(d[1]>2):
              outlier_found=True
            # others onto first copy
            lsq_fit_obj = superpose.least_squares_fit(
              reference_sites = chain.atoms().extract_xyz(),
              other_sites     = sites_cart_chain_0)
            self.back_rotation_matrices.append(lsq_fit_obj.r)
            self.back_translation_vectors.append(lsq_fit_obj.t)
    if(outlier_found): self._init()

  def _init(self):
    self.rotation_matrices = []
    self.translation_vectors = []
    self.back_rotation_matrices = []
    self.back_translation_vectors = []
    self.ph_first_chain = None

  def is_ncs_present(self):
    return len(self.rotation_matrices)>1

  def update_sites_cart(self, sites_cart_master_ncs_copy):
    cm = flex.double(sites_cart_master_ncs_copy.mean()) - \
         flex.double(self.ph_first_chain.atoms().extract_xyz().mean())
    cm = flex.vec3_double(sites_cart_master_ncs_copy.size(),list(cm))
    for i_chain, chain in enumerate(self.pdb_hierarchy.chains()):
      if(i_chain==0):
        self.ph_first_chain.atoms().set_xyz(sites_cart_master_ncs_copy)
      chain.atoms().set_xyz(
        self.back_rotation_matrices[i_chain].elems*(sites_cart_master_ncs_copy-cm) +
          self.back_translation_vectors[i_chain]+cm)

  def write_pdb_file(self, file_name, crystal_symmetry, mode):
    assert mode in ["asu", "ncs"]
    of = open(file_name, "w")
    if(mode=="ncs"):
      print >> of, iotbx.pdb.format_MTRIX_pdb_string(
        rotation_matrices   = self.back_rotation_matrices,
        translation_vectors = self.back_translation_vectors)
      print >> of, self.ph_first_chain.as_pdb_string(
        crystal_symmetry = crystal_symmetry)
    else:
      print >> of, self.pdb_hierarchy.as_pdb_string(
        crystal_symmetry = crystal_symmetry)
    of.close()
