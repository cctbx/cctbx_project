
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
    print self._plot.size()
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
    if (phi_1 % 2) == 0 :
      phi_1 -= 1
    elif (phi_2 % 2) == 0 :
      phi_2 += 1
    if (psi_1 % 2) == 0 :
      psi_1 -= 1
    elif (psi_2 % 2) == 0 :
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
  def __init__ (self, i_seqs) :
    assert (len(i_seqs) == 5)
    self.i_seqs = i_seqs

  def extract_sites (self, sites_cart) :
    from scitbx.array_family import flex
    sites = flex.double()
    for i_seq in i_seqs :
      sites.append(sites_cart[i_seq])
    return sites

class restraints_manager (object) :
  def __init__ (self, pdb_hierarchy, log=sys.stdout) :
    self._phi_psi = []
    assert (len(pdb_hierarchy.models()) == 1)
    for chain in pdb_hierarchy.models()[0].chains() :
      for conformer in chain.conformers() :
        residues = conformer.residues()
        for i, residue in enumerate(residues) :
          next_res, prev_res = None, None
          if (i > 0):
            if (residue.resseq_as_int() != residues[i-1].resseq_as_int()) :
              continue
            prev_res = residues[i-1]
          if (i < (len(residues) - 1)) :
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
            phi_psi = proxy([c1.i_seq,n2.i_seq,ca2.i_seq,c2.i_seq,n3.i_seq])
            self._phi_psi.append(phi_psi)
    print >> log, "%d Ramachandran restraints generated." % len(self._phi_psi)

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

if __name__ == "__main__" :
  exercise()
  print "OK"
