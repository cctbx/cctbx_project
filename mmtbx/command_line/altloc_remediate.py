from __future__ import absolute_import, division, print_function
import os, sys
from libtbx import phil
import libtbx.phil.command_line
import iotbx.pdb
from libtbx.utils import Sorry
from six.moves import range

master_phil_string = """

altloc_remediate
  .caption = None
{
  input
  {
    pdb_file_name = None
      .type = path
      .short_caption = model
      .help = PDB filename
      .style = bold file_type:pdb
    residue_selection = None
      .type = str
      .multiple = True
      .short_caption = Use this selection to
      .help = The
  }
  control
  {
    correct_alt_loc = True
      .type = bool
      .help = General improvement
    spread_alt_loc = False
      .type = bool
      .help = expand the alt. loc.
    use_geometry_as_spread_criteria = True
      .type = bool
    block_alt_loc = False
      .type = bool
  }
  output
  {
    file_name = None
      .type = path
      .short_caption = Output file
      .help = Defaults to current directory
      .style = bold new_file file_type:pdb
    display = False
      .type = bool
      .short_caption = Display the residues only, no remediation
  }
}
"""

master_params = master_phil_string # need for auto documentation
if False:
  print('-'*80)
  print(master_phil_string)
  print('-'*80)
master_phil = phil.parse(master_phil_string)

mainchain = set(["CA", "C", "N", "O"])
partial1 = set(["C", "O"])
partial2 = set(["N", "H", "HA"])

def check_for_exchange_hd(residue_group, remove):
  if len(remove)!=2: return remove
  ag1 = residue_group.atom_groups()[remove[0]]
  ag2 = residue_group.atom_groups()[remove[1]]
  if len(ag1.atoms())!=len(ag2.atoms()): return remove
  isotopes1 = []
  for atom in ag1.atoms():
    if atom.element.strip() not in isotopes1:
      isotopes1.append(atom.element.strip())
  if len(isotopes1)!=1: return remove
  isotopes2 = []
  for atom in ag2.atoms():
    if atom.element.strip() not in isotopes2:
      isotopes2.append(atom.element.strip())
  if len(isotopes2)!=1: return remove
  if isotopes1[0]=="H" and isotopes2[0]=="D": return []
  if isotopes1[0]=="D" and isotopes2[0]=="H": return []

def adjust_hydrogen_alt_locs(hierarchy):
  for residue_group in hierarchy.residue_groups():
    remove = []
    for atom_group_i, atom_group in enumerate(residue_group.atom_groups()):
      if not atom_group.altloc.strip(): continue
      hydrogens = []
      for atom in atom_group.atoms():
        if atom.element.strip() not in ["H", "D"]: continue
        hydrogens.append(atom)
      if len(hydrogens)==len(atom_group.atoms()):
        remove.append(atom_group_i)
    if remove:
      # check for H/D exchangeable
      remove = check_for_exchange_hd(residue_group, remove)
    if remove:
      for atom_group_i, atom_group in enumerate(residue_group.atom_groups()):
        if not atom_group.altloc.strip(): blank = atom_group
        if atom_group_i in remove:
          blank_atom_names = []
          for atom in blank.atoms(): blank_atom_names.append(atom.name.strip())
          for atom in atom_group.atoms():
            if not atom.name.strip() in blank_atom_names:
              blank.append_atom(atom.detached_copy())
            atom_group.remove_atom(atom)
          residue_group.remove_atom_group(atom_group)
  return hierarchy

def all_mainchain_alt_loc(hierarchy):
  for residue_group in hierarchy.residue_groups():
    blank = None
    atoms = set()
    move_blank_to_alt_locs = False
    move_CO_to_alt_locs = False
    move_NH_to_alt_locs = False
    for atom_group_i, atom_group in enumerate(residue_group.atom_groups()):
      if not atom_group.altloc.strip():
        blank = atom_group
        continue
      if blank is None: break
      for atom in atom_group.atoms():
        atoms.add(atom.name.strip())
      if atoms.intersection(mainchain):
        if atoms.intersection(partial1) and not atoms.intersection(partial2):
          move_CO_to_alt_locs = True
        elif atoms.intersection(partial2) and not atoms.intersection(partial1):
          move_NH_to_alt_locs = True
        else:
          move_blank_to_alt_locs = True
        break
    if 0:
      print('move_blank_to_alt_locs',move_blank_to_alt_locs)
      print('move_CO_to_alt_locs   ',move_CO_to_alt_locs)
      print('move_NH_to_alt_locs   ',move_NH_to_alt_locs)
    moving_atoms = []
    moving_group = None
    if move_blank_to_alt_locs:
      for atom_group_i, atom_group in enumerate(residue_group.atom_groups()):
        if not atom_group.altloc.strip():
          moving_group = atom_group
          moving_atoms = atom_group.atoms()
          continue
        if moving_group is None: break # blank group must be first?
    else:
      if move_CO_to_alt_locs:
        for atom_group_i, atom_group in enumerate(residue_group.atom_groups()):
          if atom_group.altloc.strip(): continue
          for atom in atom_group.atoms():
            if atom.name.strip() in partial1:
              moving_atoms.append(atom)
              moving_group = atom_group
      if move_NH_to_alt_locs:
        for atom_group_i, atom_group in enumerate(residue_group.atom_groups()):
          if atom_group.altloc.strip(): continue
          for atom in atom_group.atoms():
            if atom.name.strip() in partial2:
              moving_atoms.append(atom)
              moving_group = atom_group

    if moving_atoms:
      for atom in moving_atoms:
        for atom_group_i, atom_group in enumerate(residue_group.atom_groups()):
          if not atom_group.altloc.strip(): continue
          atom_group.append_atom(atom.detached_copy())
        moving_group.remove_atom(atom)
      if len(moving_group.atoms())==0:
        residue_group.remove_atom_group(moving_group)
  return hierarchy

def general_corrections(hierarchy):
  hierarchy = all_mainchain_alt_loc(hierarchy)
  hierarchy = adjust_hydrogen_alt_locs(hierarchy)
  return hierarchy

def get_alt_loc_type(residue_group, verbose=False):
  altloc_list = set()
  altloc_list_h = set()
  atom_names = []
  atom_names_h = []
  for atom in residue_group.atoms():
    if atom.name.strip() not in atom_names_h:
        atom_names_h.append(atom.name.strip())
    if atom.parent().altloc.strip():
      altloc_list_h.add(atom.name.strip())
    if atom.element.strip() in ["H", "D"]: continue
    if atom.name.strip() not in atom_names: atom_names.append(atom.name.strip())
    if atom.parent().altloc.strip():
      altloc_list.add(atom.name.strip())
  for atom_group in residue_group.atom_groups(): break
  if not altloc_list:
    return None
  if len(altloc_list)==len(atom_names): # this should be all atoms
    return "all"
  if altloc_list.intersection(mainchain)==mainchain:
    if altloc_list.difference(mainchain):
      print("mainchain only",mainchain)
      print('intersection',altloc_list.intersection(mainchain))
      print('difference',altloc_list.difference(mainchain))
    else:
      assert 0
  elif altloc_list.intersection(partial1)==partial1:
    return 'partial'
  elif altloc_list.intersection(partial2)==partial2:
    return 'partial'
  else:
    return "sidechain only"
  assert 0

def get_alt_locs(hierarchy, verbose=False):
  altlocs = {}
  for model in hierarchy.models():
    if verbose: print('model: "%s"' % model.id)
    for chain in model.chains():
      if verbose: print('chain: "%s"' % chain.id)
      altlocs.setdefault(chain.id, {})
      three = []
      for residue_group in chain.residue_groups():
        if verbose: print('  residue_group: resseq="%s" icode="%s"' % (
          residue_group.resseq, residue_group.icode))
        if len(residue_group.atom_groups())>1:
          altlocs[chain.id].setdefault((residue_group.resseq, residue_group.icode),
                                       len(residue_group.atom_groups()),
                                       )
  return altlocs

def generate_threes(hierarchy, verbose=False):
  altlocs = get_alt_locs(hierarchy)
  for model in hierarchy.models():
    if verbose: print('model: "%s"' % model.id)
    for chain in model.chains():
      if verbose: print('chain: "%s"' % chain.id)
      three = []
      for residue_group in hierarchy.residue_groups():
        if verbose: print('  residue_group: resseq="%s" icode="%s"' % (
          residue_group.resseq, residue_group.icode))
        if chain.id != residue_group.parent().id: break
        three.append(residue_group)
        if len(three)>3: del three[0]
        if len(three)<3: continue
        chain_dict = altlocs.get(chain.id, {})
        key = (three[1].resseq,
               three[1].icode,
               )
        is_altloc = chain_dict.get(key, None)
        if is_altloc not in [2, 3]: continue
        yield three

def get_distance(x, y):
  from math import sqrt
  d = 0
  for i in range(3):
    d += (x[i]-y[i])**2
  return sqrt(d)

def get_geometry_flags(three):
  rc = {}
  cbs = []
  for atom in three[1].atoms():
    if atom.name.strip()=="CB": cbs.append(atom)
  if len(cbs)==2:
    d = get_distance(cbs[0].xyz, cbs[1].xyz)
    rc["c_beta_dist"] = d
  else:
    rc["c_beta_dist"] = None
  return rc

def expand_altloc(hierarchy, verbose=False):
  for three in generate_threes(hierarchy):
    types = []
    for residue_group in three:
      types.append(get_alt_loc_type(residue_group))
    if types[1] in ["sidechain only"]: continue
    spread_altlocs = []
    for atom_group in three[1].atom_groups():
      spread_altlocs.append(atom_group.altloc)
    for j in range(0,3,2):
      if len(three[j].atom_groups())==1:
        clone = three[j].atom_groups()[0].detached_copy()
        three[j].append_atom_group(clone)
        for i, atom_group in enumerate(three[j].atom_groups()):
          atom_group.altloc = spread_altlocs[i]
  hierarchy.atoms().reset_serial()
  if verbose: hierarchy.show()
  return hierarchy

def get_altloc_data(hierarchy):
  rc = []
  for residue_group in hierarchy.residue_groups():
    rc.append([])
    for atom_group in residue_group.atom_groups():
      rc[-1].append(atom_group.altloc)
  return rc

def print_residue_group(residue_group):
  print('residue_group', residue_group.resseq)
  atom_names = []
  for atom_group in residue_group.atom_groups():
    print('  altloc "%s"' % atom_group.altloc)
    for i, atom in enumerate(atom_group.atoms()):
      print("    %2d %s" % (i, atom.format_atom_record()))
      if atom.name.strip() not in atom_names: atom_names.append(atom.name.strip())
  #print 'atoms',len(atom_names)
  #print 'end'

def display_model(hierarchy):
  for residue_group in hierarchy.residue_groups():
    print_residue_group(residue_group)

def atom_name_in_atom_group(atom_group, atom):
  for tmp in atom_group.atoms():
    if tmp.name.strip()==atom.name.strip(): return True
  return False

def merge_altloc_into_space(residue_group):
  space_atom_group = None
  for atom_group in residue_group.atom_groups():
    if not atom_group.altloc.strip():
      space_atom_group = atom_group
      continue
    assert space_atom_group
    assert space_atom_group.resname == atom_group.resname
    for atom in atom_group.atoms():
      if not atom_name_in_atom_group(space_atom_group, atom):
        tmp = atom.detached_copy()
        space_atom_group.append_atom(tmp)
      atom_group.remove_atom(atom)
    residue_group.remove_atom_group(atom_group)

def spread_to_c_alpha(residue_group,
                      spread_altlocs,
                      pre_peptide=True,
                      verbose=1,
  ):
  def _check_atom_groups_ok(residue_group, atoms):
    for atom_group in residue_group.atom_groups():
      if not atom_group.altloc.strip(): continue
      for atom in atom_group.atoms():
        if atom.name.strip() not in atoms: return False
    return True
  atoms = partial2
  if pre_peptide:
    atoms = partial1
  if not _check_atom_groups_ok(residue_group, atoms):
    print('not spreading to')
    assert 0
    return
  if len(residue_group.atom_groups())!=1:
    # merge
    merge_altloc_into_space(residue_group)
  assert len(residue_group.atom_groups())==1
  resname = residue_group.atom_groups()[0].resname
  spread_atoms = {}
  for sa in spread_altlocs:
    spread_atoms[sa] = []
    ag = iotbx.pdb.hierarchy.atom_group()
    ag.altloc = sa
    ag.resname = resname
    residue_group.append_atom_group(ag)
  if verbose:
    print('duplicating',atoms)
  for atom_group in residue_group.atom_groups():
    for atom in atom_group.atoms():
      if atom.name.strip() in atoms:
        for sa in spread_atoms:
          spread_atoms[sa].append(atom.detached_copy())
        atom_group.remove_atom(atom)
  for sa in sorted(spread_atoms):
    if not sa.strip(): continue
    for atom_group in residue_group.atom_groups():
      if atom_group.altloc==sa: break
    else: assert 0
    for new_atom in spread_atoms[sa]:
      if verbose: print('adding to "%s" <- %s' % (sa, atom.quote()))
      atom_group.append_atom(new_atom)

def spread_to_residue(residue_group):
  blank_atom_group = residue_group.atom_groups()[0]
  altlocs = []
  for atom_group in residue_group.atom_groups():
    if atom_group == blank_atom_group: continue
    altlocs.append(atom_group.altloc)
  for atom in blank_atom_group.atoms():
    for atom_group in residue_group.atom_groups():
      if atom_group == blank_atom_group: continue
      d_atom = atom.detached_copy()
      atom_group.append_atom(d_atom)
    blank_atom_group.remove_atom(atom)
  residue_group.remove_atom_group(blank_atom_group)

def correct_altloc(hierarchy, max_c_beta_deviation=0.2, verbose=False):
  for three in generate_threes(hierarchy):
    types = []
    for residue_group in three:
      types.append(get_alt_loc_type(residue_group))
    spread_altlocs = []
    for atom_group in three[1].atom_groups():
      spread_altlocs.append(atom_group.altloc)
    geometry_flags = get_geometry_flags(three)
    if types[:2] == [None, "all"]:
      spread_to_c_alpha(three[0], spread_altlocs)
    elif types[:2] == [None, "sidechain only"]:
      if geometry_flags["c_beta_dist"]>max_c_beta_deviation:
        spread_to_residue(three[1])
        spread_to_c_alpha(three[0], spread_altlocs)

    if types[1:] == ["all", None]:
      spread_to_c_alpha(three[2], spread_altlocs, pre_peptide=False)
    elif types[1:] == ["sidechain only", None]:
      if (geometry_flags["c_beta_dist"] and
          geometry_flags["c_beta_dist"] > max_c_beta_deviation):
        spread_to_c_alpha(three[2], spread_altlocs, pre_peptide=False)

  hierarchy.atoms().reset_serial()
  if verbose: hierarchy.show()
  return hierarchy

def block_alt_loc(hierarchy):
  for residue_group in hierarchy.residue_groups():
    moving_group = None
    if len(residue_group.atom_groups())==1: continue
    for atom_group in residue_group.atom_groups():
      if not atom_group.altloc.strip():
        moving_group = atom_group
        break
    if moving_group:
      for atom in moving_group.atoms():
        for atom_group_i, atom_group in enumerate(residue_group.atom_groups()):
          if not atom_group.altloc.strip(): continue
          atom_group.append_atom(atom.detached_copy())
        moving_group.remove_atom(atom)
      if len(moving_group.atoms())==0:
        residue_group.remove_atom_group(moving_group)

  hierarchy.atoms().reset_serial()
  return hierarchy

def correct_occupancies(hierarchy):
  for atom_group in hierarchy.atom_groups():
    if not atom_group.altloc.strip(): continue
    occ=None
    for atom in atom_group.atoms():
      if atom.occ!=1:
        occ = atom.occ
    for atom in atom_group.atoms():
      if atom.occ==1:
        if occ:
          atom.occ=occ
        else:
          atom.occ=.5

def run(rargs):
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="altloc_remediate")
  #
  phils = []
  phil_args = []
  pdbs = []
  for arg in args:
    if os.path.isfile(arg):
      if iotbx.pdb.is_pdb_file(arg):
        pdbs.append(arg)
        continue
      try :
        file_phil = phil.parse(file_name=arg)
      except RuntimeError :
        pass
      else :
        phils.append(file_phil)
    else :
      phil_args.append(arg)
      phils.append(argument_interpreter.process(arg))
  working_phil = master_phil.fetch(sources=phils)
  #working_phil.show()
  working_params = working_phil.extract()

  if not getattr(working_params, "altloc_remediate", False):
    raise Sorry('Must have a "altloc_remediate" scope in phil file')
  in_scope = working_params.altloc_remediate.input
  control = working_params.altloc_remediate.control
  output = working_params.altloc_remediate.output
  if control.block_alt_loc and control.correct_alt_loc:
    raise Sorry('''Block residue altlocs are not consider "correct".
  To block residues set control.correct_alt_loc=False''')
  #
  for i, pdb in enumerate(pdbs):
    if i==0 and not in_scope.pdb_file_name:
      in_scope.pdb_file_name = pdbs[i]
  #
  if in_scope.pdb_file_name is None:
    raise Sorry("Must supply a protein PDB file")
  #
  if not output.file_name:
    output.file_name = in_scope.pdb_file_name
    d = os.path.dirname(output.file_name)
    output.file_name = os.path.basename(output.file_name)
    output.file_name = output.file_name.split(".")[0]
    if control.spread_alt_loc:
      output.file_name += "_spread"
    if control.correct_alt_loc:
      output.file_name += "_correct"
    if control.block_alt_loc:
      output.file_name += "_block"
    output.file_name += ".pdb"
    output.file_name = os.path.join(d, output.file_name)
  #
  preamble = output.file_name.split(".")[0]
  print("\n  Writing effective parameters to %s.eff\n" % preamble)
  print("#phil __ON__")
  working_phil.format(python_object=working_params).show()
  print("#phil __OFF__\n")
  f=open("%s.eff" % preamble, "w")
  f.write(working_phil.format(python_object=working_params).as_str())
  f.close()

  pdb_inp = iotbx.pdb.input(in_scope.pdb_file_name)
  hierarchy = pdb_inp.construct_hierarchy()
  #hierarchy.show()
  if output.display:
    display_model(hierarchy)
    return
  hierarchy = general_corrections(hierarchy)

  if control.spread_alt_loc:
    hierarchy = expand_altloc(hierarchy)

  if control.correct_alt_loc:
    hierarchy = correct_altloc(hierarchy)

  if control.block_alt_loc:
    hierarchy = block_alt_loc(hierarchy)

  correct_occupancies(hierarchy)

  print("  Writing output to %s" % output.file_name)
  f=open(output.file_name, "w")
  f.write(hierarchy.as_pdb_string(
    crystal_symmetry=pdb_inp.crystal_symmetry()),
          )
  f.close()

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(args)
