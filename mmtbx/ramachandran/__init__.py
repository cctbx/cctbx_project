
from __future__ import division
import libtbx.load_env
from libtbx.math_utils import ifloor, iceil
import sys
import os

# TODO convert to C++
class lookup_table (object) :
  def __init__ (self, lines) :
    from scitbx.array_family import flex
    self._plot = flex.double()
    for line in lines :
      val, phi, psi = line.split()
      assert ((int(phi) % 2 == 1) and (int(psi) % 2 == 1))
      self._plot.append(float(val))
    assert (self._plot.size() == 32400)
    grid = flex.grid((180,180))
    self._plot.reshape(grid)

  def get_score (self, phi, psi) :
    return self.__getitem__((phi, psi))

  def __getitem__ (self, phi_psi) :
    # http://en.wikipedia.org/wiki/Bilinear_interpolation
    (phi, psi) = phi_psi
    phi_1 = ifloor(phi)
    phi_2 = iceil(phi)
    psi_1 = ifloor(psi)
    psi_2 = iceil(psi)
    if ((phi_1 % 2) == 0) :
      phi_1 -= 1
    elif ((phi_2 % 2) == 0) :
      phi_2 += 1
    if ((psi_1 % 2) == 0) :
      psi_1 -= 1
    elif ((psi_2 % 2) == 0) :
      psi_2 += 1
    r11 = self._getitem(phi_1, psi_1)
    r12 = self._getitem(phi_1, psi_2)
    r21 = self._getitem(phi_2, psi_1)
    r22 = self._getitem(phi_2, psi_2)
    d_phi_d_psi = float((phi_2 - phi_1) * (psi_2 - psi_1))
    r_phi_psi = ((r11 / d_phi_d_psi) * (phi_2 - phi) * (psi_2 - psi)) + \
                ((r21 / d_phi_d_psi) * (phi - phi_1) * (psi_2 - psi)) + \
                ((r12 / d_phi_d_psi) * (phi_2 - phi) * (psi - psi_1)) + \
                ((r22 / d_phi_d_psi) * (phi - phi_1) * (psi - psi_1))
    return r_phi_psi

  def _getitem (self, phi, psi) :
    assert (phi < 180 and phi > -180 and psi < 180 and psi > -180)
    assert ((phi % 2) == 1 and (psi % 2) == 1)
    i = int((phi + 179) / 2)
    j = int((psi + 179) / 2)
    return self._plot[(i,j)]

  def compute_angular_gradients (self, phi, psi, delta=2) :
    r_start = self.get_score(phi, psi)
    r_phi_1 = self.get_score(phi - delta, psi)
    r_phi_2 = self.get_score(phi + delta, psi)
    d_r_d_phi = (r_phi_2 - r_phi_1) / (delta * 2)
    r_psi_1 = self.get_score(phi, psi - delta)
    r_psi_2 = self.get_score(phi, psi + delta)
    d_r_d_psi = (r_psi_2 - r_psi_1) / (delta * 2)
    return (d_r_d_phi, d_r_d_psi)

tables = {}

def load_tables () :
  for residue_type in ["ala", "gly", "prepro", "pro"] :
    file_name = libtbx.env.find_in_repositories(
      relative_path="chem_data/rotarama_data/%s.rama.combined.data" %
        residue_type,
      test=os.path.isfile)
    f = open(file_name, "r")
    t = lookup_table(f.readlines())
    tables[residue_type] = t

if (len(tables) == 0) :
  load_tables()

class proxy (object) :
  def __init__ (self, i_seqs, residue_type) :
    assert (len(i_seqs) == 5)
    self.i_seqs = i_seqs
    self.residue_type = residue_type

  def _extract_sites (self, sites_cart) :
    from scitbx.array_family import flex
    sites = flex.vec3_double()
    for i_seq in self.i_seqs :
      sites.append(sites_cart[i_seq])
    return sites

  def compute_gradients_finite_differences (self,
                                            sites_cart,
                                            gradients,
                                            rama_table,
                                            delta=0.1,
                                            weight=1.0) :
    assert ((delta > 0) and (delta < 1.0))
    from cctbx import geometry_restraints
    def dihedral (sites) :
      return geometry_restraints.dihedral(
        sites=sites,
        angle_ideal=0,
        weight=1).angle_model
    #from scitbx.math import dihedral
    sites = self._extract_sites(sites_cart)
    #gradients = flex.vec3_double()
    for k, i_seq in enumerate(self.i_seqs) :
      k_grad = list(gradients[k])
      xyz = list(sites[k])
      for u in [0,1,2] :
        xyz[u] -= delta
        sites[k] = xyz
        phi1 = dihedral(list(sites[0:4]))
        psi1 = dihedral(list(sites[1:]))
        xyz[u] += (delta * 2)
        sites[k] = xyz
        phi2 = dihedral(list(sites[0:4]))
        psi2 = dihedral(list(sites[1:]))
        xyz[u] -= delta
        sites[k] = xyz
        phi0 = dihedral(list(sites[0:4]))
        psi0 = dihedral(list(sites[1:]))
        r1 = rama_table.get_score(phi1, psi1)
        r2 = rama_table.get_score(phi2, psi2)
        w = weight * rama_table.get_score(phi0, psi0)
        k_grad[u] = w * ((r2 - r1) / delta)
      gradients[k] = k_grad
    return gradients

class restraints_manager (object) :
  def __init__ (self, pdb_hierarchy, log=sys.stdout) :
    from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter
    from cctbx import geometry_restraints
    self._phi_psi = []
    assert (len(pdb_hierarchy.models()) == 1)
    for chain in pdb_hierarchy.models()[0].chains() :
      for conformer in chain.conformers() :
        residues = conformer.residues()
        for i, residue in enumerate(residues) :
          if (not residue.resname in one_letter_given_three_letter) :
            continue
          next_res, prev_res = None, None
          resseq2 = residue.resseq_as_int()
          resseq1, resseq3 = None, None
          if (i > 0):
            resseq1 = residues[i-1].resseq_as_int()
            if (resseq2 != (resseq1 + 1)) :
              continue
            prev_res = residues[i-1]
          if (i < (len(residues) - 1)) :
            resseq3 = residues[i+1].resseq_as_int()
            if (resseq2 != (resseq3 - 1)) :
              continue
            next_res = residues[i+1]
          if (next_res is not None) and (prev_res is not None) :
            c1, n2, ca2, c2, n3 = None, None, None, None, None
            for atom in prev_res.atoms() :
              if (atom.name == " C  ") :
                c1 = atom
                break
            for atom in residue.atoms() :
              if (atom.name == " N  ") :
                n2 = atom
              elif (atom.name == " CA ") :
                ca2 = atom
              elif (atom.name == " C  ") :
                c2 = atom
            for atom in next_res.atoms() :
              if (atom.name == " N  ") :
                n3 = atom
            if (None in [c1, n2, ca2, c2, n3]) :
              #print >> log, "  incomplete backbone for %s %d-%d, skipping." % \
              #  (chain.id, resseq1, resseq3)
              continue
            pep1 = geometry_restraints.bond(
              sites=[c1.xyz,n2.xyz],
              distance_ideal=1,
              weight=1)
            pep2 = geometry_restraints.bond(
              sites=[c2.xyz,n3.xyz],
              distance_ideal=1,
              weight=1)
            if (pep1.distance_model > 4) or (pep2.distance_model > 4) :
              continue
            i_seqs = [c1.i_seq,n2.i_seq,ca2.i_seq,c2.i_seq,n3.i_seq]
            if (residue.resname == "PRO") :
              residue_type = "pro"
            elif (residue.resname == "GLY") :
              residue_type = "gly"
            elif (residues[i+1].resname == "PRO") :
              residue_type = "prepro"
            else :
              residue_type = "ala"
            phi_psi = proxy(i_seqs, residue_type)
            self._phi_psi.append(phi_psi)
    print >> log, "%d Ramachandran restraints generated." % len(self._phi_psi)

  def compute_gradients_finite_differences (self, sites_cart, delta=0.1) :
    from scitbx.array_family import flex
    gradients = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    for phi_psi in self._phi_psi :
      table = tables[phi_psi.residue_type]
      phi_psi.compute_gradients_finite_differences(
        sites_cart=sites_cart,
        gradients=gradients,
        rama_table=tables[phi_psi.residue_type],
        delta=0.1)

def exercise() :
  from iotbx import file_reader
  from libtbx.test_utils import approx_equal
  (phi, psi) = (60.12, -57.9)
  assert approx_equal(tables["ala"][(phi, psi)], -3.558, eps=0.001)
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  pdb_in = file_reader.any_file(pdb_file).file_object
  pdb_hierarchy = pdb_in.construct_hierarchy()
  m = restraints_manager(pdb_hierarchy)
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  m.compute_gradients_finite_differences(sites_cart)

if __name__ == "__main__" :
  exercise()
  print "OK"
