from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
import math
import iotbx.pdb
from scitbx.array_family import flex
from scitbx import matrix
from cctbx import adptbx
from mmtbx.utils import rotatable_bonds

def get_axes_and_atoms_i_seqs(pdb_hierarchy, mon_lib_srv):
  get_class = iotbx.pdb.common_residue_names_get_class
  axes_and_atoms_i_seqs = []
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for conformer in residue_group.conformers():
          for residue in conformer.residues():
            if(get_class(residue.resname) == "common_amino_acid"):
              aaa = rotatable_bonds.axes_and_atoms_aa_specific(
                residue = residue, mon_lib_srv = mon_lib_srv)
              if(aaa is not None):
                for aaa_ in aaa:
                  tmp_ = []
                  axis = flex.size_t()
                  moving = flex.size_t()
                  axis_ = aaa_[0]
                  atoms_ = aaa_[1]
                  for i_seq, atom in enumerate(residue.atoms()):
                    if(i_seq in axis_): axis.append(atom.i_seq)
                    elif(i_seq in atoms_): moving.append(atom.i_seq)
                  assert len(axis) == 2
                  assert len(moving) > 0
                  tmp_.append([axis, moving])
                  axes_and_atoms_i_seqs.append(tmp_)
  return axes_and_atoms_i_seqs

def check_u_cart(axis, u_cart):
  ev = adptbx.eigenvalues(u_cart)
  assert ev[0] > 0
  assert approx_equal(ev[1:],[0,0])
  es = adptbx.eigensystem(u_cart)
  evs = es.vectors(0),es.vectors(1),es.vectors(2)
  a1 = math.acos((axis[0]*evs[0][0]+axis[1]*evs[0][1]+axis[2]*evs[0][2])/ \
    (evs[0][0]**2+evs[0][1]**2+evs[0][2]**2)) *180/math.pi
  a2 = math.acos((axis[0]*evs[1][0]+axis[1]*evs[1][1]+axis[2]*evs[1][2])/ \
    (evs[1][0]**2+evs[1][1]**2+evs[1][2]**2)) *180/math.pi
  a3 = math.acos((axis[0]*evs[2][0]+axis[1]*evs[2][1]+axis[2]*evs[2][2])/ \
    (evs[2][0]**2+evs[2][1]**2+evs[2][2]**2)) *180/math.pi
  assert approx_equal(a1, 90.)

def set_ladp(xray_structure, axes_and_atoms_i_seqs, value, depth,
             enable_recursion=True):
  sc = (math.pi/180)
  sites_cart = xray_structure.sites_cart()
  scatterers = xray_structure.scatterers()
  all_selections = flex.size_t()
  u_carts = flex.sym_mat3_double(sites_cart.size(), [0,0,0,0,0,0])
  for i_seq, aaa_ in enumerate(axes_and_atoms_i_seqs):
    if(enable_recursion): query = i_seq >= depth
    else: query = i_seq == depth
    if(query): all_selections.extend(aaa_[0][1])
  for i_seq, r in enumerate(axes_and_atoms_i_seqs):
    if(enable_recursion): query = i_seq >= depth
    else: query = i_seq == depth
    if(query):
      for aaai in r:
        G1 = flex.double(sites_cart[aaai[0][0]])
        G2 = flex.double(sites_cart[aaai[0][1]])
        g = G2-G1
        dg = math.sqrt(g[0]**2+g[1]**2+g[2]**2)
        lx,ly,lz = g/dg
        l = [lx,ly,lz]
        L = matrix.sqr((lx**2,lx*ly,lx*lz, lx*ly,ly**2,ly*lz, lx*lz,ly*lz,lz**2))
        for i_seq_moving in aaai[1]:
          site_cart = sites_cart[i_seq_moving]
          delta = flex.double(site_cart) - G1
          A = matrix.sqr(
            (0,delta[2],-delta[1], -delta[2],0,delta[0], delta[1],-delta[0],0))
          u_cart = (value * A * L * A.transpose() * sc).as_sym_mat3()
          check_u_cart(axis = l, u_cart = u_cart)
          scatterers[i_seq_moving].flags.set_use_u_aniso(True)
          u_carts[i_seq_moving] = list(flex.double(u_carts[i_seq_moving]) +
            flex.double(u_cart))
  xray_structure.set_u_cart(u_cart = u_carts, selection = all_selections)
  return xray_structure
