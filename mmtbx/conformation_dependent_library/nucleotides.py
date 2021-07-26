import os, sys

def update_restraints(hierarchy,
                      geometry, # restraints_manager,
                      current_geometry=None, # xray_structure!!
                      sites_cart=None,
                      rdl_proxies=None,
                      log=None,
                      verbose=False,
                      ):
  from restraintlib import launcher
  from restraintlib.printer import TuplePrinter
  from restraintlib.restraints import analyze_pdb_hierarhy
  # rc = launcher.cdl_nucleotides(hierarchy)
  override_sigma=False
  restraint_groups = launcher.load_restraints_lib()
  printer = TuplePrinter(override_sigma=override_sigma)
  rc = analyze_pdb_hierarhy(hierarchy, restraint_groups, restraint_groups, printer)
  bond_restraints = []
  angle_restraints = {}
  atoms = hierarchy.atoms()
  for i_seqs, ideal, esd in rc:
    if 0:
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
  for i, (i_seqs, ideal, esd) in enumerate(bond_restraints):
    bond=geometry.bond_params_table.lookup(*list(i_seqs))
    remove.append(i)
    bond.distance_ideal=ideal
    bond.weight = 1/esd**2
  remove.reverse()
  for r in remove:
    del bond_restraints[r]
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
      if 1:
        old_angle_ideal=angle_proxy.angle_ideal
        old_angle_weight=angle_proxy.weight
        print(" i_seqs %-15s initial %12.3f %12.3f" % (
          angle_proxy.i_seqs,
          angle_proxy.angle_ideal,
          angle_proxy.weight,
          ), end=' ', file=log)
      assert angle_proxy.angle_ideal<181
      angle_proxy.angle_ideal = angle_restraints[i_seqs][0]
      # if not ignore_esd:
      angle_proxy.weight = 1/angle_restraints[i_seqs][1]**2
      del angle_restraints[i_seqs]
      if 1:
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
    assert 0
  assert not angle_restraints, 'not finished angle_restraints: %s' % angle_restraints
  return
