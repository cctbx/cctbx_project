from __future__ import absolute_import, division, print_function
import sys
import time

from mmtbx.conformation_dependent_library.hpdl_database import get_hpdl_database
from six.moves import range

def get_histidine_protonation(ag):
  lookup = {"HD1" : 0, # ND1
            "HE2" : 0, # NE2
            }
  for atom in ag.atoms():
    if atom.name.strip() in lookup:
      lookup[atom.name.strip()]+=1
  if lookup["HD1"] and lookup["HE2"]: return "ND1 and NE2 protonated"
  elif lookup["HD1"]: return "Only ND1 protonated"
  elif lookup["HE2"]: return "Only NE2 protonated"
  return None

def setup_restraints():
  return None

def update_restraints(hierarchy,
                      geometry, # restraints_manager,
                      current_geometry=None, # xray_structure!!
                      sites_cart=None,
                      esd_factor=1.,
                      log=None,
                      verbose=False,
                      ):
  #
  def _set_or_reset_bond_restraints(geometry,
                                    lookup,
                                    ignore_esd=True,
                                    log=None,
                                    verbose=False,
                                    ):
    count = 0
    for i_seqs, values in lookup.items():
      if len(i_seqs)!=2: continue
      bond=geometry.bond_params_table.lookup(*list(i_seqs))
      assert bond
      if verbose:
        from math import sqrt
        key = list(i_seqs)
        key.append("names")
        key = tuple(key)
        names = lookup.get(key, "")
        old_bond_ideal=bond.distance_ideal
        old_bond_weight=bond.weight
        print(" i_seqs %-15s %s initial %8.3f %8.3f %8.3f final %8.3f %8.3f %8.3f" % (
          i_seqs,
          names,
          bond.distance_ideal,
          bond.weight,
          1/sqrt(bond.weight),
          values[0],
          1/values[1]**2,
          values[1],
        ), file=log)
        if bond.distance_ideal!=old_bond_ideal or bond.weight!=old_bond_weight:
          print('RESETTING', file=log)
      bond.distance_ideal=values[0]
      if not ignore_esd:
        bond.weight = 1/values[1]**2
      count+=1
    return count
  #
  def _set_or_reset_angle_restraints(geometry,
                                     lookup,
                                     ignore_esd=True,
                                     log=None,
                                     verbose=False,
                                     ):
    count = 0
    for angle_proxy in geometry.angle_proxies:
      if angle_proxy.i_seqs in lookup:
        if verbose:
          old_angle_ideal=angle_proxy.angle_ideal
          old_angle_weight=angle_proxy.weight
          print(" i_seqs %-15s initial %12.3f %12.3f" % (
            angle_proxy.i_seqs,
            angle_proxy.angle_ideal,
            angle_proxy.weight,
            ), end=' ', file=log)
        assert angle_proxy.angle_ideal<181
        angle_proxy.angle_ideal = lookup[angle_proxy.i_seqs][0]
        if not ignore_esd:
          angle_proxy.weight = esd_factor/lookup[angle_proxy.i_seqs][1]**2
        if verbose:
          print("final   %12.3f %12.3f" % (
            angle_proxy.angle_ideal,
            angle_proxy.weight,
          ), file=log)
          if old_angle_ideal!=angle_proxy.angle_ideal or old_angle_weight!=angle_proxy.weight:
            print("RESETTING", file=log)
        count += 1
    return count
  #
  t0=time.time()
  sites_cart = None
  if current_geometry:
    sites_cart = current_geometry.sites_cart()
  i_seqs_restraints = {}
  #
  def _alt_loc_atom_generator(residue_group, atom_group):
    atoms = []
    for ag in residue_group.atom_groups():
      if ag.altloc.strip()=="" or ag.altloc.strip()==atom_group.altloc.strip():
        for atom in ag.atoms(): yield atom
  #
  hpdl_database = get_hpdl_database() #include_hydrogens=False)
  count=0
  counts = {}
  for model in hierarchy.models():
    #if verbose: print 'model: "%s"' % model.id
    for chain in model.chains():
      #if verbose: print 'chain: "%s"' % chain.id
      for residue_group in chain.residue_groups():
        for atom_group in residue_group.atom_groups():
          if atom_group.resname!="HIS": continue
          protonation = get_histidine_protonation(atom_group)
          if verbose or 1:
            rc = predict_protonation(atom_group)
            if rc is None:
              s = "%s" % rc
            else:
              s = "%0.1f, %0.1f" % tuple(rc)
            print('%satom group "%s" has %-22s (%s)' % (
              ' '*6,
              atom_group.id_str(),
              protonation,
              s,
            ), file=log)
            #interpret_his1_his2(*tuple(rc))
          counts.setdefault(protonation, 0)
          if protonation is None: continue
          counts[protonation]+=1
          count+=1
          restraints = hpdl_database[protonation]
          for names, values in restraints.items():
            i_seqs = []
            atoms = []
            for name in names:
              # need to test this...
              for atom in _alt_loc_atom_generator(residue_group, atom_group):
                if name.strip()==atom.name.strip():
                  i_seqs.append(atom.i_seq)
                  atoms.append(atom.name)
                  break
            if len(i_seqs)!=len(names): continue
            i_seqs_restraints[tuple(i_seqs)] = values
            if len(i_seqs)!=2:
              i_seqs.reverse()
              i_seqs_restraints[tuple(i_seqs)] = values
            if verbose:
              i_seqs.append("names")
              i_seqs_restraints[tuple(i_seqs)]=atoms

  count_b = _set_or_reset_bond_restraints(geometry,
                                          i_seqs_restraints,
                                          log=log,
                                          #verbose=verbose,
                                          )
  count_a = _set_or_reset_angle_restraints(geometry,
                                           i_seqs_restraints,
                                           log=log,
                                           #verbose=verbose,
                                           )
  #
  print("    Number of bonds, angles adjusted : %d, %d in %s HIS" % (
    count_b,
    count_a,
    count,
    ), file=log)
  #return rdl_proxies

def predict_protonation(ag, verbose=False):
  from cctbx import geometry
  from mmtbx.monomer_library.linking_utils import get_distance2
  from math import sqrt
#  from cctbx_geometry_restraints_ext import *
#  from cctbx.array_family import flex
  """
his1 = -37.35*X1 + 15.57*X2 - 0.64*X3 + 0.76*X4 + 17.30

his2 = -2.16*X1 - 6.08*X2 + 0.56*X3 + 0.42*X4 - 94.46

X1 = ND1-CE1
X2 = CE1-NE2
X3 = -ND1-
X4 = -NE2-

his1<0 ~> ND1 protonated
his1>0
  his2<0 ~> NE2 protonated
  his2>0 ~> doubly protonated
  """
  bonds = {
    ("ND1", "CE1") : None,
    ("NE2", "CE1") : None,
    }
  if verbose:
    for atom in ag.atoms():
      print(atom.quote())
  for i, tmp_atom in enumerate(bonds):
    atoms = []
    for j in range(2):
      for atom in ag.atoms():
        if atom.name.strip()==tmp_atom[j]:
          atoms.append(atom)
          break
    if len(atoms)==2:
      d2=get_distance2(*atoms)
      bonds[tmp_atom]=sqrt(d2)
  angles = {
    ("CG", "ND1", "CE1") : None,
    ("CE1","NE2", "CD2") : None,
    }
  for i, tmp_atom in enumerate(angles):
    atoms = []
    for j in range(3):
      for atom in ag.atoms():
        if atom.name.strip()==tmp_atom[j]:
          atoms.append(atom.xyz)
          break
    if len(atoms)==3:
      angle=geometry.angle((atoms)).angle_model
      angles[tmp_atom]=angle
  if None in bonds.values() or None in angles.values(): return None
  his1 =  17.30 - 0.64*angles[("CG", "ND1", "CE1")]
  his1 +=         0.76*angles[("CE1","NE2", "CD2")]
  his1 -=        37.35*bonds[("ND1", "CE1")]
  his1 +=        15.57*bonds[("NE2", "CE1")]
  his2 = -94.46 + 0.56*angles[("CG", "ND1", "CE1")]
  his2 +=         0.42*angles[("CE1","NE2", "CD2")]
  his2 -=         2.61*bonds[("ND1", "CE1")]
  his2 -=         6.08*bonds[("NE2", "CE1")]
  if verbose:
    print('his1',his1)
    print('his2',his2)
  return (his1, his2)

def interpret_his1_his2(his1, his2):
  def _interpret_his1_his2(his1, his2, limit=0):
    if his1<-limit:
      return "Only ND1 protonated"
    elif his1>0:
      if his2<-limit:
        return "Only NE2 protonated"
      elif his2>0:
        return "ND1 and NE2 protonated"
    else:
      return None
  limit = 1
  print(limit, _interpret_his1_his2(his1, his2, limit))
  limit = 0
  print(limit, _interpret_his1_his2(his1, his2, limit))

def run(pdb_filename):
  from iotbx import pdb
  print(pdb_filename)
  pdb_inp = pdb.input(pdb_filename)
  pdb_hieratchy = pdb_inp.construct_hierarchy()
  for residue_group in pdb_hieratchy.residue_groups():
    for atom_group in residue_groups.atom_groups():
      print(get_histidine_protonation(atom_group))

if __name__=="__main__":
  run(sys.argv[1])
