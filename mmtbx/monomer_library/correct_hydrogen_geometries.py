from __future__ import absolute_import, division, print_function
from six.moves import zip

def get_bad_hydrogen_i_seqs(hierarchy,
                            restraints_manager=None,
                            xray_structure=None,
                            sites_cart=None,
                            verbose=False,
                            ):
  from cctbx import geometry

  def get_bonded(atom_group, atom, exclude=[]): # slow
    min_d2=1000
    min_atom=None
    for atom1 in atom_group.atoms():
      if atom1.name.strip()==atom.name.strip(): continue
      process=False
      for atom2 in exclude:
        if atom1.name.strip()==atom2.name.strip():
          break
      else:
        process=True
      if not process: continue
      if(atom.element.strip() in ["H", "D", "T"] and
         atom1.element.strip() in ["H", "D", "T"]
         ): continue
      d2 = (atom.xyz[0]-atom1.xyz[0])**2
      d2 += (atom.xyz[1]-atom1.xyz[1])**2
      d2 += (atom.xyz[2]-atom1.xyz[2])**2
      if d2<min_d2:
        min_d2=d2
        min_atom=atom1
    if min_atom:
      return [min_atom]
    else:
      return []

  def get_angles(atom_group, atom, verbose=False):
    bad_hydrogen_count = 0
    corrected_hydrogen_count = []
    bonded = get_bonded(atom_group, atom)
    if not bonded:
      bad_hydrogen_count+=1
      if verbose: print('not bonded: %s' % atom.format_atom_record())
      return bad_hydrogen_count, corrected_hydrogen_count
    for ba in bonded:
      angled = get_bonded(atom_group, ba, exclude=[atom])
    if not angled:
      bad_hydrogen_count+=1
      if verbose: print('not angled: %s' % atom.format_atom_record())
      return bad_hydrogen_count, corrected_hydrogen_count
    if angled[0].element.strip() in ["H", "D", "T"]:
      return bad_hydrogen_count, corrected_hydrogen_count
    try: angle = geometry.angle((atom.xyz,ba.xyz,angled[0].xyz)).angle_model
    except Exception:
      print('  Bad angle "%s"' % (atom.format_atom_record()[:26]))
      bad_hydrogen_count +=1
      return bad_hydrogen_count, corrected_hydrogen_count
    if angle<85.:
      xyz=[0,0,0]
      xyz[0] = ba.xyz[0]*2-atom.xyz[0]
      xyz[1] = ba.xyz[1]*2-atom.xyz[1]
      xyz[2] = ba.xyz[2]*2-atom.xyz[2]
      angle = geometry.angle((xyz,ba.xyz,angled[0].xyz)).angle_model
      if angle>95.:
        atom.xyz = tuple(xyz)
        if verbose: print('  Inverted "%s" ' % (atom.format_atom_record()[:26]))
        corrected_hydrogen_count.append(atom.i_seq)
    return bad_hydrogen_count, corrected_hydrogen_count

  def get_i_seqs(hierarchy):
    bad_hydrogen_count=0
    corrected_hydrogen_count=[]
    for model in hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          #if len(residue_group.atom_groups())>1: continue
          for atom_group_i, atom_group in enumerate(residue_group.atom_groups()):
            for i, atom in enumerate(atom_group.atoms()):
              if atom.element.strip() in ["H", "D", "T"]:
                rc = get_angles(atom_group, atom, verbose=verbose)
                bad_hydrogen_count += rc[0]
                corrected_hydrogen_count += rc[1]
    return corrected_hydrogen_count, None

  def get_i_seqs_from_restraints_manager(
      hierarchy,
      restraints_manager,
      xray_structure = None,
      sites_cart = None,
      ):
    if sites_cart is None and xray_structure:
      sites_cart = xray_structure.sites_cart()
    if sites_cart is None:
      xray_structure = hierarchy.extract_xray_structure()
      sites_cart = xray_structure.sites_cart()

    i_seqs=[]
    xyzs=[]
    angle_proxies_simple = restraints_manager.geometry.angle_proxies
    atoms = hierarchy.atoms()
    for proxy in angle_proxies_simple:
      i_seq, j_seq, k_seq = proxy.i_seqs
      if(atoms[i_seq].element.strip() in ["H", "D"] or
         atoms[k_seq].element.strip() in ["H", "D"]
         ):
        if(atoms[i_seq].element.strip() in ["H", "D"] and
           atoms[k_seq].element.strip() in ["H", "D"]
           ):
          continue
        if(atoms[i_seq].element.strip() in ["H", "D"]):
          i_h = i_seq
          site_i = sites_cart[i_seq]
          site_k = sites_cart[k_seq]
        else:
          i_h = k_seq
          site_i = sites_cart[k_seq]
          site_k = sites_cart[i_seq]
        site_j = sites_cart[j_seq]

        if i_h in i_seqs: continue

        angle = geometry.angle((site_i, site_j, site_k)).angle_model
        if angle<85.:
          xyz=[0,0,0]
          xyz[0] = site_j[0]*2-site_i[0]
          xyz[1] = site_j[1]*2-site_i[1]
          xyz[2] = site_j[2]*2-site_i[2]
          angle = geometry.angle((xyz, site_j, site_k)).angle_model
          if angle>95.:
            xyzs.append(tuple(xyz))
            i_seqs.append(i_h)
    return i_seqs, xyzs

  if restraints_manager:
    i_seqs, xyzs = get_i_seqs_from_restraints_manager(
      hierarchy,
      restraints_manager,
      xray_structure=xray_structure,
      sites_cart=sites_cart,
      )
  else:
    i_seqs, xyzs = get_i_seqs(hierarchy)
  return i_seqs, xyzs

def correct_hydrogen_geometries(hierarchy,
                                restraints_manager=None,
                                xray_structure=None,
                                sites_cart=None,
                                verbose=False,
                                ):
  assert xray_structure or sites_cart
  bad_hydrogen_count=0
  corrected_hydrogen_count=[]
  if len(hierarchy.models())>1:
    print("  \nModel files with more than one model are ignored\n")
    return bad_hydrogen_count, corrected_hydrogen_count

  i_seqs, xyzs = get_bad_hydrogen_i_seqs(hierarchy,
                                         restraints_manager=restraints_manager,
                                         xray_structure=xray_structure,
                                         sites_cart=sites_cart,
                                         )
  if xyzs is None:
    return 0, i_seqs
  assert len(i_seqs)==len(xyzs)
  if xray_structure:
    sites_cart = xray_structure.sites_cart()
  for i_seq, xyz in zip(i_seqs, xyzs):
    sites_cart[i_seq] = xyz
  if xray_structure:
    xray_structure.set_sites_cart(sites_cart)
  else:
    for i, atom in enumerate(hierarchy.atoms()):
      atom.xyz = sites_cart[i]
  corrected_hydrogen_count = i_seqs
  atoms = hierarchy.atoms()
  for i, i_seq in enumerate(corrected_hydrogen_count):
    corrected_hydrogen_count[i] = atoms[i_seq].id_str()
  return bad_hydrogen_count, corrected_hydrogen_count
