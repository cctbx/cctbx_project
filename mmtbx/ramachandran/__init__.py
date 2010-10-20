
from __future__ import division
import libtbx.load_env
from libtbx.math_utils import ifloor, iceil
import sys
import os

import boost.python
ext = boost.python.import_ext("mmtbx_ramachandran_ext")
from mmtbx_ramachandran_ext import *

tables = {}
def load_tables () :
  from scitbx.array_family import flex
  for residue_type in ["ala", "gly", "prepro", "pro"] :
    file_name = libtbx.env.find_in_repositories(
      relative_path="chem_data/rotarama_data/%s.rama.combined.data" %
        residue_type,
      test=os.path.isfile)
    f = open(file_name, "r")
    data = flex.double()
    for line in f.readlines() :
      val, phi, psi = line.split()
      assert ((int(phi) % 2 == 1) and (int(psi) % 2 == 1))
      data.append(float(val))
    t = lookup_table(data, 180)
    tables[residue_type] = t

if (len(tables) == 0) :
  load_tables()

class generic_proxy (object) :
  restraint_type = None

class proxy (generic_proxy) :
  restraint_type = "ramachandran"
  def __init__ (self, i_seqs, residue_type) :
    assert (len(i_seqs) == 5)
    self.i_seqs = i_seqs
    self.residue_type = residue_type

  def extract_sites (self, sites_cart) :
    from scitbx.array_family import flex
    sites = flex.vec3_double()
    for i_seq in self.i_seqs :
      sites.append(sites_cart[i_seq])
    return sites

def phi_psi_restraints_residual_sum (proxies,
                                     sites_cart,
                                     gradient_array=None,
                                     delta=0.1) :
  from scitbx.array_family import flex
  sum = 0
  for phi_psi in proxies :
    sum += _compute_gradients(
      proxy=phi_psi,
      sites_cart=sites_cart,
      gradients=gradient_array,
      delta=0.1)
  return sum

def dihedral (sites) :
  from cctbx import geometry_restraints
  return geometry_restraints.dihedral(
    sites=sites,
    angle_ideal=0,
    weight=1).angle_model

def _compute_gradients (proxy,
                        sites_cart,
                        gradients,
                        delta=0.1,
                        weight=1.0,
                        use_finite_differences=False) :
  rama_table = tables[proxy.residue_type]
  assert ((delta > 0) and (delta < 1.0))
  from scitbx.array_family import flex
  residual = 0
  new_gradients = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
  #from scitbx.math import dihedral
  sites = proxy.extract_sites(sites_cart)
  phi0 = dihedral(list(sites[0:4]))
  psi0 = dihedral(list(sites[1:]))
  residual = rama_table.get_energy(phi0, psi0)
  #print phi0, psi0, residual
  if use_finite_differences :
    _compute_gradients_finite_differences(
      proxy=proxy,
      sites=sites,
      new_gradients=new_gradients,
      delta=delta,
      weight=weight)
  else :
    #assert 0
    _compute_gradients_analytical(
      proxy=proxy,
      phi_start=phi0,
      psi_start=psi0,
      sites=sites,
      new_gradients=new_gradients,
      delta_cart=delta,
      delta_ang=5.0,
      weight=weight)
  gradients += new_gradients
  return residual

def _compute_gradients_finite_differences (proxy,
                                           sites,
                                           new_gradients,
                                           delta,
                                           weight) :
  rama_table = tables[proxy.residue_type]
  for k, i_seq in enumerate(proxy.i_seqs) :
    grad = list(new_gradients[i_seq])
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
      r1 = rama_table.get_energy(phi1, psi1)
      r2 = rama_table.get_energy(phi2, psi2)
      w = weight # * (0.1 * abs(residual))
      grad[u] += w * ((r2 - r1) / (delta * 2))
    new_gradients[i_seq] = grad

def _compute_gradients_analytical (proxy,
                                   phi_start,
                                   psi_start,
                                   sites,
                                   new_gradients,
                                   delta_ang,
                                   delta_cart,
                                   weight) :
  rama_table = tables[proxy.residue_type]
  r_phi_1 = rama_table.get_energy(phi_start - delta_ang, psi_start)
  r_phi_2 = rama_table.get_energy(phi_start + delta_ang, psi_start)
  d_r_d_phi = (r_phi_2 - r_phi_1) / (delta_ang * 2)
  r_psi_1 = rama_table.get_energy(phi_start, psi_start - delta_ang)
  r_psi_2 = rama_table.get_energy(phi_start, psi_start + delta_ang)
  d_r_d_psi = (r_psi_2 - r_psi_1) / (delta_ang * 2)
  for k, i_seq in enumerate(proxy.i_seqs) :
    grad = list(new_gradients[i_seq])
    xyz = list(sites[k])
    for u in [0,1,2] :
      xyz[u] -= delta_cart
      sites[k] = xyz
      phi_1 = dihedral(list(sites[0:4]))
      psi_1 = dihedral(list(sites[1:]))
      xyz[u] += (delta_cart * 2)
      sites[k] = xyz
      phi_2 = dihedral(list(sites[0:4]))
      psi_2 = dihedral(list(sites[1:]))
      xyz[u] -= delta_cart
      sites[k] = xyz
      d_phi_d_u = (phi_2 - phi_1) / (delta_cart * 2)
      d_psi_d_u = (psi_2 - psi_1) / (delta_cart * 2)
      d_r_d_u = (d_r_d_phi * d_phi_d_u) + (d_r_d_psi * d_psi_d_u)
      grad[u] = weight * d_r_d_u
    new_gradients[i_seq] = grad
  #return (d_r_d_phi, d_r_d_psi)

def extract_proxies (pdb_hierarchy, log=sys.stdout) :
  from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter
  from cctbx import geometry_restraints
  proxies = []
  for model in pdb_hierarchy.models() :
    for chain in model.chains() :
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
            proxies.append(phi_psi)
  print >> log, "%d Ramachandran restraints generated." % len(proxies)
  return proxies

class generic_restraints_helper (object) :
  def __init__ (self) :
    pass

  def restraints_residual_sum (self,
                               sites_cart,
                               proxies,
                               gradient_array=None,
                               unit_cell=None) :
    ramachandran_proxies = []
    for proxy in proxies :
      if (proxy.restraint_type == "ramachandran") :
        ramachandran_proxies.append(proxy)
    return phi_psi_restraints_residual_sum(
      sites_cart=sites_cart,
      proxies=ramachandran_proxies,
      gradient_array=gradient_array)
