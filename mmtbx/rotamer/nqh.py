from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
import libtbx.load_env
import iotbx.pdb
import sys, re
from mmtbx.refinement.flip_peptide_side_chain import flip_residue
from iotbx.pdb.utils import all_chain_ids
from mmtbx.utils import run_reduce_with_timeout
from mmtbx.validation.clashscore import check_and_report_reduce_failure

def seg_id_to_chain_id(pdb_hierarchy):
  import string
  two_character_chain_ids = []
  segid_list = []
  seg_dict = {}
  for atom in pdb_hierarchy.atoms():
    if atom.segid not in segid_list:
      segid_list.append(atom.segid)
  lower_letters = string.ascii_lowercase
  upper_letters = string.ascii_uppercase
  two_character_chain_ids = all_chain_ids()
  for id in segid_list:
    chainID = two_character_chain_ids[0]
    seg_dict[id] = chainID
    two_character_chain_ids.remove(chainID)
  return seg_dict

def find_bare_chains_with_segids(pdb_hierarchy):
  bare_chains = False
  for chain in pdb_hierarchy.chains():
    if chain.id in ['', ' ', '  ']:
      segid = None
      for atom in chain.atoms():
        if segid == None:
          segid = atom.segid
        elif segid != None and segid != atom.segid:
          #require that each chain have a unique segid for this logic
          return False
      if segid != None and segid not in ['', ' ', '  ', '   ', '    ']:
        bare_chains = True
  return bare_chains

def assign_chain_ids(pdb_hierarchy, seg_dict):
  rename_txt = ""
  for chain in pdb_hierarchy.chains():
    if chain.id in ['', ' ', '  ']:
      segid = None
      for atom in chain.atoms():
        if segid == None:
          segid = atom.segid
        elif segid != atom.segid:
          print(segid, atom.segid)
          raise Sorry("multiple segid values defined for chain")
      new_id = seg_dict[segid]
      chain.id = new_id
      rename_txt = rename_txt + \
      "segID %s renamed chain %s for Reduce N/Q/H analysis\n" % (segid, new_id)
  return rename_txt

def force_unique_chain_ids(pdb_hierarchy):
  used_chain_ids = []
  two_char = all_chain_ids()
  #filter all used chains
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      cur_id = chain.id
      if cur_id in two_char:
        two_char.remove(cur_id)
  #force unique chain ids
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      cur_id = chain.id
      if cur_id not in used_chain_ids:
        used_chain_ids.append(cur_id)
      else:
        new_id = two_char[0]
        chain.id = new_id
        two_char.remove(new_id)

def get_nqh_flips(pdb_hierarchy):
  pdb_string = pdb_hierarchy.as_pdb_string()
  output = run_reduce(pdb_string=pdb_string, remove_hydrogens=True)
  user_mods = []
  score_dict = {}
  atom_notes = []
  for line in output.stdout_lines:
    if re.match(r'USER  MOD', line):
      flip = re.search(r':(.{15}):FLIP(.*):sc=.*[F|C]\(o=(.*),f=(.*)\)',line)
      if(flip):
        if flip.group(3)[-1:] == '!':
          score_o = float(flip.group(3)[:-1])
        else:
          score_o = float(flip.group(3))
        if flip.group(4)[-1:] == '!':
          score_f = float(flip.group(4)[:-1])
        else:
          score_f = float(flip.group(4))
        score_diff = score_f - score_o
        if score_f <= -2.0:
          user_mods.append('!'+flip.group(1).strip())
        else:
          user_mods.append(flip.group(1).strip())
        score_dict[flip.group(1).strip()] = score_diff
        if re.search(r'HIS',line):
          atom_notes.append([flip.group(1).strip(),flip.group(2).strip()])
  return user_mods, atom_notes, score_dict

# return true if H atoms are inconsistent with Reduce protonation state
# if no D or E H's are present, always return False to allow flip
def check_for_his_h(residue, atom_notes, key):
  atom_list =[]
  for atom in residue.atoms():
    atom_list.append(atom.name)

  #check for no H atoms
  if not (" HD1" in atom_list or " HD2" in atom_list or " HE1" in atom_list or
          " HE2" in atom_list):
    return False
  for entry in atom_notes:
    if entry[0] == key:
      if entry[1].strip() == "no HD1":
        if not ( (" HD2" in atom_list and " HE1" in atom_list and " HE2" in atom_list) and
                 not " HD1" in atom_list ):
          return True
      elif entry[1].strip() == "no HE2":
        if not ( (" HD1" in atom_list and " HD2" in atom_list and " HE1" in atom_list) and
                 not " HE2" in atom_list):
          return True
      elif entry[1].strip() == "+bothHN":
        if not (" HD1" in atom_list and " HD2" in atom_list and " HE1" in atom_list and
                " HE2" in atom_list):
          return True
  return False

def run_reduce(pdb_string, remove_hydrogens=True):
  assert (libtbx.env.has_module(name="reduce"))
  trim = " -quiet -trim -"
  build = " -quiet -build -allalt -"
  input_str = ""
  if(remove_hydrogens):
    clean = run_reduce_with_timeout(parameters=trim,
                                    stdin_lines=pdb_string)
    check_and_report_reduce_failure(clean, pdb_string, "reduce_failure.pdb")
    input_str = "\n".join(clean.stdout_lines)
  else:
    input_str = pdb_string
  output = run_reduce_with_timeout(parameters=build,
                                   stdin_lines=input_str)
  check_and_report_reduce_failure(output, input_str, "reduce_failure.pdb")
  return output

def flip_selected(pdb_hierarchy, # changed in-place
                  mon_lib_srv,
                  flip_list,
                  atom_notes,
                  score_dict,
                  log):
  sites_cart_start = pdb_hierarchy.atoms().extract_xyz()
  total_flipped = 0
  #
  # XXX TEMPORARY DISABLED. THIS DOES NOT COVER ALL. USE GRM INSTEAD OF LINK STUFF. XXX
  #
  ##find excluded residues that are covalently linked
  exclude_sidechain = []
  #for arg in apply_cif_links:
  #  for res in arg.pdbres_pair:
  #    resname = res[8:11]
  #    if resname.upper() in ['ASN', 'GLN', 'HIS']:
  #      chainID = res[11:13]
  #      resnum = res[13:-1]
  #      key = chainID.strip()+resnum+resname
  #      #exclude_sidechain.append(res[8:-1].strip())
  #      exclude_sidechain.append(key)
  #print exclude_sidechain
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        conformers = residue_group.conformers()
        for conformer in residue_group.conformers():
          residue = conformer.only_residue()
          key = chain.id+residue.resid()+residue.resname+"    "+conformer.altloc
          key = key.strip()
          if key[0:9] in exclude_sidechain:
            print(key[0:9]+' is covalently modified...skipping', file=log)
          elif '!'+key in flip_list:
            print('** '+key+' **'+' both conformations clash, **PLEASE CHECK MANUALLY**', file=log)
          elif key in flip_list:
            if residue.resname=="HIS":
              if check_for_his_h(residue, atom_notes, key)==True:
                continue
            print(key, file=log)
            # here we are going to do flips in hierarchy and then simply
            # extract xray_structure with updated coordinates
            flip_residue(residue,mon_lib_srv)
            total_flipped += 1
  print("\nTotal number of N/Q/H flips: %d\n" % total_flipped, file=log)

def flip(pdb_hierarchy, log=None, mon_lib_srv=None):
  if mon_lib_srv is None:
    mon_lib_srv = mmtbx.monomer_library.server.server()
  if(log is None): log = sys.stdout
  pdb_hierarchy.atoms().reset_i_seq()
  print("Analyzing N/Q/H residues for possible flip corrections...", file=log)
  tmp_pdb_hierarchy = pdb_hierarchy.deep_copy()
  tmp_pdb_hierarchy.atoms().reset_i_seq()
  #analyze chain/segid relationship
  bare_chains = find_bare_chains_with_segids(pdb_hierarchy=pdb_hierarchy)
  if bare_chains:
    seg_dict = seg_id_to_chain_id(pdb_hierarchy=tmp_pdb_hierarchy)
    rename_txt = assign_chain_ids(pdb_hierarchy=tmp_pdb_hierarchy,
                                        seg_dict=seg_dict)
    print(rename_txt, file=log)
  if tmp_pdb_hierarchy.overall_counts().n_duplicate_chain_ids > 0:
    force_unique_chain_ids(
      pdb_hierarchy=tmp_pdb_hierarchy)
  try:
    flip_list, atom_notes, score_dict = get_nqh_flips(
      pdb_hierarchy = tmp_pdb_hierarchy)
  except OSError as e :
    if (e.errno == 32) :
      print("WARNING: system error attempting to run Reduce", file=log)
      print("         this may indicate a crash or OS conflict", file=log)
      print("         phenix.refine will continue running, but you", file=log)
      print("         should check Asn/Gln/His residues visually.", file=log)
      print("         contact help@phenix-online.org if this problem", file=log)
      print("         occurs frequently.", file=log)
      return pdb_hierarchy
    else :
      raise e
  if len(flip_list) == 0:
    print("\nNo N/Q/H corrections needed this macrocycle", file=log)
  else:
    print("\nFlipped N/Q/H residues before XYZ refinement:", file=log)
    flip_selected(
      pdb_hierarchy = pdb_hierarchy, # changed in-place
      mon_lib_srv   = mon_lib_srv,
      flip_list     = flip_list,
      atom_notes    = atom_notes,
      score_dict    = score_dict,
      log           = log)
  return pdb_hierarchy
