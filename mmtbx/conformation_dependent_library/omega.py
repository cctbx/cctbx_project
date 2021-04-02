from __future__ import absolute_import, division, print_function
import copy

from mmtbx.conformation_dependent_library.bond_angle_registry import \
  bond_angle_registry
from mmtbx.conformation_dependent_library import generate_protein_threes
from mmtbx.conformation_dependent_library.multi_residue_class import \
  RestraintsRegistry

from mmtbx.conformation_dependent_library.omega_database import omega_database

from mmtbx.conformation_dependent_library.cdl_utils import get_res_type_group
from six.moves import range

columns = [
  "",
  "",
  "mACNA", # Ca(0) - C(0) - N(1) - Ca(1)
  "sACNA",
  ]

def get_restraint_values(threes, verbose=False):

  res_type_group = get_res_type_group(
    threes[-2].resname,
    threes[-1].resname,
  )
  if verbose: print(threes[-2].resname,threes[-1].resname,res_type_group)
  if res_type_group is None: return None
  restraint_values = []
  key = threes.get_cdl_key(force_plus_one=True,
                           omega_cdl=True,
                          )
  if verbose: print('key',key)
  if key is None: return None
  previous_key = None
  if len(key)==4:
    previous_key = key[:2]
  key = key[-2:]
  if previous_key:
    rv = omega_database[res_type_group][previous_key]
    restraint_values.append(rv)
  else:
    restraint_values.append(None)
  rv = omega_database[res_type_group][key]
  restraint_values.append(rv)
  return restraint_values

def setup_restraints(geometry,
                     verbose=False,
                     ):
  ba_registry = bond_angle_registry()
  for angle in geometry.dihedral_proxies:
    ba_registry[angle.i_seqs]=angle
  return ba_registry

def apply_updates(self,
                  restraint_value_pairs,
                  cdl_proxies,
                  ideal=True,
                  esd=True,
                  esd_factor=1.,
                  average=True,
                  verbose=False,
                  ):
    if not average:
      if restraint_value_pairs[0]=="I":
        print(restraint_values)
        assert 0
        return
    atoms = self.get_i_seqs()
    for j, restraint_values in enumerate(restraint_value_pairs):
      if restraint_values is None: continue
      for i, value in enumerate(restraint_values):
        if verbose:
          print('i,value',i,value)
        if i<2: continue
        if columns[i][0]=="s": continue
        code = columns[i][1:]
        names = []
        if code=="ACNA":
          if j:
            names = ["CA_i", "C_i", "N_plus_1", "CA_plus_1"]
          else:
            names = ["CA_minus_1", "C_minus_1", "N_i", "CA_i"]
        for j in range(len(names)):
          names[j] = atoms[names[j]].i_seq
        if len(names)==4:
          rnames = copy.deepcopy(names)
          rnames.reverse()
          angle_proxy = cdl_proxies.get(tuple(names), None)
          if angle_proxy is None:
            angle_proxy = cdl_proxies.get(tuple(rnames), None)
          if angle_proxy is None:
            outl=""
            for key in atoms:
              outl += "\n    %-10s %s" % ( key, atoms[key].quote())
            raise Sorry("""CDL angle to be changed not set in model.
  Possible problems:
    Residue on special positions.

  Check:%s""" % outl)
          if verbose:
            print(" i_seqs %-15s initial %12.3f %12.3f final %12.3f %12.3f" % (
              angle_proxy.i_seqs,
              angle_proxy.angle_ideal,
              angle_proxy.weight,
              restraint_values[i],
              1/restraint_values[i+1]**2,
              ))
          names.sort()
          self.registry[tuple(names)] = restraint_values[i]
          if ideal: angle_proxy.angle_ideal = restraint_values[i]
          if esd: angle_proxy.weight = esd_factor * 1/restraint_values[i+1]**2
        else:
          assert 0

def update_restraints(hierarchy,
                      geometry, # restraints_manager,
                      current_geometry=None, # xray_structure!!
                      sites_cart=None,
                      cdl_proxies=None,
                      ideal=True,
                      esd=True,
                      esd_factor=1.,
                      log=None,
                      verbose=False,
                      ):
  global registry
  registry = RestraintsRegistry()
  if current_geometry:
    assert not sites_cart
    sites_cart = current_geometry.sites_cart()
  if sites_cart:
    pdb_atoms = hierarchy.atoms()
    # XXX PDB_TRANSITION VERY SLOW
    for j_seq, atom in enumerate(pdb_atoms):
      atom.xyz = sites_cart[j_seq]

  threes = None
  average_updates = 0
  total_updates = 0
  for threes in generate_protein_threes(hierarchy,
                                        geometry,
                                        cdl_class=True,
                                        omega_cdl=True,
                                        #verbose=verbose,
                                        ):
    threes.apply_updates = apply_updates
    if threes.cis_group(omega_cdl=True):
      if verbose and 0:
        print('cis '*20)
        print(threes)
      continue

    restraint_values = get_restraint_values(threes, verbose=verbose)
    if restraint_values is None: continue
    if restraint_values[0]=="I":
      average_updates += 1
    else:
      total_updates += 1
    threes.apply_updates(threes,
                         restraint_values,
                         cdl_proxies,
                         ideal=ideal,
                         esd=esd,
                         esd_factor=esd_factor,
                         verbose=verbose,
                         )
  if registry.n:
    threes.apply_average_updates(registry)
    assert 0
  geometry.reset_internals()
  if verbose and threes and threes.errors:
    if log:
      log.write("  Residues not completely updated with CDL restraints\n\n")
    for line in threes.errors:
      if log:
        log.write("%s\n" % line)
      else:
        print(line)
  return geometry

def run(filename):
  if False:
    for i in range(-188,188):
      print(i,round_to_ten(i),abs(i-round_to_ten(i)))
    assert 0

  pdb_inp = iotbx.pdb.input(filename)
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.atoms().reset_serial()
  assert 0, "broken run method"
  update_restraints(hierarchy,
                    #verbose=True,
                    )

if __name__=="__main__":
  if 0:
    psi = -180
    lookup = "Gly_nonxpro"
    print(lookup)
    for phi in range(170,181):
      key = (round_to_ten(psi),round_to_ten(phi))
      print('key',psi,phi,round_to_ten(psi),round_to_ten(phi),key, end=' ')
      print(cdl_database[lookup][key][:4])

  run(sys.argv[1])
