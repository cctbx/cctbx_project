import os, sys

def update_restraints(hierarchy,
                      geometry, # restraints_manager,
                      # current_geometry=None, # xray_structure!!
                      # sites_cart=None,
                      # rdl_proxies=None,
                      use_phenix_esd=True,
                      log=None,
                      verbose=False,
                      ):
  from restraintlib import launcher
  from restraintlib.printer import TuplePrinter
  from restraintlib.restraints import analyze_pdb_hierarhy
  restraint_groups = launcher.load_restraints_lib()
  printer = TuplePrinter(override_sigma=use_phenix_esd) # use_phenix_esd not really needed here
  rc = analyze_pdb_hierarhy(hierarchy, restraint_groups, restraint_groups, printer)
  bond_restraints = []
  angle_restraints = {}
  atoms = hierarchy.atoms()
  for i_seqs, ideal, esd in rc:
    if verbose:
      if len(i_seqs)==2:
        print('%-20s %s %s %s %s' % (i_seqs,
              atoms[i_seqs[0]].quote(),
              atoms[i_seqs[1]].quote(),
              ideal,
              esd,
              ))
      else:
        print('%-20s %s %s %s %s %s' % (i_seqs,
              atoms[i_seqs[0]].quote(),
              atoms[i_seqs[1]].quote(),
              atoms[i_seqs[2]].quote(),
              ideal,
              esd,
              ))
    if len(i_seqs)==2:
      bond_restraints.append([i_seqs, ideal, esd])
    elif len(i_seqs)==3:
      angle_restraints[i_seqs]=[ideal, esd]
    else:
      assert 0
  remove=[]
  n_bonds=0
  for i, (i_seqs, ideal, esd) in enumerate(bond_restraints):
    bond=geometry.bond_params_table.lookup(*list(i_seqs))
    remove.append(i)
    bond.distance_ideal=ideal
    if not use_phenix_esd:
      bond.weight = 1/esd**2
    n_bonds+=1
  remove.reverse()
  for r in remove:
    del bond_restraints[r]
  n_angles=0
  for angle_proxy in geometry.angle_proxies:
    i_seqs = list(angle_proxy.i_seqs)
    ap = None
    if tuple(i_seqs) in angle_restraints:
      ap = angle_proxy
      i_seqs = tuple(i_seqs)
    if ap is None:
      i_seqs.reverse()
      i_seqs = tuple(i_seqs)
      if tuple(i_seqs) in angle_restraints:
        ap = angle_proxy
    if ap:
      if verbose:
        old_angle_ideal=angle_proxy.angle_ideal
        old_angle_weight=angle_proxy.weight
        print(" i_seqs %-15s initial %12.3f %12.3f" % (
          angle_proxy.i_seqs,
          angle_proxy.angle_ideal,
          angle_proxy.weight,
          ), end=' ', file=log)
      assert angle_proxy.angle_ideal<181
      angle_proxy.angle_ideal = angle_restraints[i_seqs][0]
      if not use_phenix_esd:
        angle_proxy.weight = 1/angle_restraints[i_seqs][1]**2
      del angle_restraints[i_seqs]
      n_angles+=1
      if verbose:
        print(" i_seqs %-15s final %12.3f %12.3f\n" % (
          angle_proxy.i_seqs,
          angle_proxy.angle_ideal,
          angle_proxy.weight,
          ), end=' ', file=log)
  assert not bond_restraints
  for i_seqs in angle_restraints:
    print(i_seqs,
          atoms[i_seqs[0]].quote(),
          atoms[i_seqs[1]].quote(),
          atoms[i_seqs[2]].quote(),
          ideal,
          esd,
          )
  assert not angle_restraints, 'not finished angle_restraints: %s' % angle_restraints
  if n_bonds or n_angles:
    print('''  CDL for nucleotides adjusted restraints counts
    bonds  : %-5d
    angles : %-5d''' % (n_bonds, n_angles),
      file=log)
  return
