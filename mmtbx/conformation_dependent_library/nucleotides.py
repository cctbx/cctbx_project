from __future__ import division, print_function

restraintlib_installed=True
try:
  from restraintlib import launcher
  from restraintlib.printer import TuplePrinter
  from restraintlib.restraints import analyze_pdb_hierarhy
except ImportError as e:
  restraintlib_installed = False

cofsky_esd = {
  'AD' : {
    ("C1'", 'N1') : 0.02,
    ("C3'", "O3'"): 0.02,
    ("O4'", "C1'", 'N9') : 0.3,
    ("O4'", "C1'", 'N1') : 0.3,
  },
  'U' : {
    ('C5', 'C6', 'N1') : 0.5,
    ('C6', 'N1', 'C2') : 0.5,
  },
  'C' : {
    ('C6', 'N1', 'C2') : 0.5,
  },
  'A' : {
    ('N7', 'C8', 'N9') : 0.5,
  },
}
cofsky_esd['TD']=cofsky_esd['AD']
cofsky_esd['CD']=cofsky_esd['AD']
cofsky_esd['GD']=cofsky_esd['AD']

def update_restraints(hierarchy,
                      geometry, # restraints_manager,
                      use_phenix_esd=True,
                      csd_factor=1.,
                      use_cofsky_esd=False,
                      log=None,
                      verbose=False,
                      ):
  if not restraintlib_installed:
    print('  RestraintLib not installed\n', file=log)
    return False
  if use_cofsky_esd:
    assert 0, 'use_cofsky_esd not implemented'
  restraint_groups = launcher.load_restraints_lib()
  printer = TuplePrinter(override_sigma=use_phenix_esd) # use_phenix_esd not really needed here
  rc = analyze_pdb_hierarhy(hierarchy, restraint_groups, restraint_groups, printer)
  bond_restraints = []
  angle_restraints = {}
  if verbose:
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
  c_bonds=0
  non_standard_issues=False
  for i, (i_seqs, ideal, esd) in enumerate(bond_restraints):
    if verbose:
      print(i, i_seqs, atoms[i_seqs[0]].quote(), atoms[i_seqs[1]].quote(), ideal, esd)
    bond=geometry.bond_params_table.lookup(*list(i_seqs))
    remove.append(i)
    if bond is None:
      atoms = hierarchy.atoms()
      print('bond between %s and %s not found in restraints' % (
        atoms[i_seqs[0]].quote(),
        atoms[i_seqs[1]].quote()
        ))
      non_standard_issues=True
      continue
    if bond.distance_ideal!=ideal:
      c_bonds+=1
    bond.distance_ideal=ideal
    if not use_phenix_esd:
      if csd_factor!=1: esd*=csd_factor
      bond.weight = 1/esd**2
    n_bonds+=1
  remove.reverse()
  for r in remove:
    del bond_restraints[r]
  n_angles=0
  c_angles=0
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
      if angle_proxy.angle_ideal!=angle_restraints[i_seqs][0]:
        c_angles+=1
      angle_proxy.angle_ideal = angle_restraints[i_seqs][0]
      if not use_phenix_esd:
        esd = angle_restraints[i_seqs][1]
        if csd_factor!=1: esd*=csd_factor
        angle_proxy.weight = 1/esd**2
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
  if not non_standard_issues:
    assert not angle_restraints, 'not finished angle_restraints: %s' % angle_restraints
  if n_bonds or n_angles:
    print('''  CDL for nucleotides adjusted restraints counts
    coord. : total (unchanged)
    bonds  : %5d (%5d)
    angles : %5d (%5d)''' % (n_bonds, c_bonds, n_angles, c_angles),
      file=log)
  return (n_bonds, c_bonds, n_angles, c_angles)
