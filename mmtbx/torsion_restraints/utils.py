from libtbx.utils import format_exception, Sorry
from libtbx import Auto
from iotbx.pdb import common_residue_names_get_class
from mmtbx.validation.cbetadev import cbetadev
from iotbx.pdb import amino_acid_codes
from libtbx import group_args
import ccp4io_adaptbx
import math
import mmtbx.alignment

def selection(string, cache):
  return cache.selection(
    string=string)

def iselection(string, cache=None):
  return selection(string=string, cache=cache).iselection()

def phil_atom_selection_multiple(
      cache,
      string_list,
      allow_none=False,
      allow_auto=False,
      raise_if_empty_selection=True):
  result = []
  for string in string_list:
    if (string is None):
      if (allow_none): return None
      raise Sorry('Atom selection cannot be None:\n  =None')
    elif (string is Auto):
      if (allow_auto): return Auto
      raise Sorry('Atom selection cannot be Auto:\n  %s=Auto')
    try:
        result.append(selection(string=string, cache=cache).iselection())
    except KeyboardInterrupt: raise
    except Exception, e: # keep e alive to avoid traceback
      fe = format_exception()
      raise Sorry('Invalid atom selection:\n  %s=%s\n  (%s)' % (
        'reference_group', string, fe))
    if (raise_if_empty_selection and result.count(True) == 0):
      raise Sorry('Empty atom selection:\n  %s=%s' % (
        'reference_group', string))
  return result

def phil_atom_selections_as_i_seqs_multiple(cache,
                                            string_list):
  result = []
  iselection = phil_atom_selection_multiple(
        cache=cache,
        string_list=string_list,
        raise_if_empty_selection=False)
  for i in iselection:
    if (i.size() == 0):
      raise Sorry("No atom selected")
    for atom in i:
      result.append(atom)
  return result

def is_residue_in_selection(i_seqs, selection):
  for i_seq in i_seqs:
    if i_seq not in selection:
      return False
  return True

def get_i_seqs(atoms):
  i_seqs = []
  for atom in atoms:
    i_seqs.append(atom.i_seq)
  return i_seqs

def modernize_rna_resname(resname):
  if common_residue_names_get_class(resname,
       consider_ccp4_mon_lib_rna_dna=True) == "common_rna_dna" or \
     common_residue_names_get_class(resname,
       consider_ccp4_mon_lib_rna_dna=True) == "ccp4_mon_lib_rna_dna":
    tmp_resname = resname.strip()
    if len(tmp_resname) == 1:
      return "  "+tmp_resname
    elif len(tmp_resname) == 2:
      if tmp_resname[0:1].upper() == 'D':
        return " "+tmp_resname.upper()
      elif tmp_resname[1:].upper() == 'D':
        return " D"+tmp_resname[0:1].upper()
      elif tmp_resname[1:].upper() == 'R':
        return "  "+tmp_resname[0:1].upper()
    elif tmp_resname in ["ADE", "CYT", "GUA", "URI"]:
      return "  "+tmp_resname[0:1].upper()
  return resname

def modernize_rna_atom_name(atom):
   new_atom = atom.replace('*',"'")
   if new_atom == " O1P":
     new_atom = " OP1"
   elif new_atom == " O2P":
     new_atom = " OP2"
   return new_atom

def build_name_hash(pdb_hierarchy):
  i_seq_name_hash = dict()
  for atom in pdb_hierarchy.atoms():
    atom_name = atom.pdb_label_columns()[0:4]
    resname = atom.pdb_label_columns()[5:8]
    if resname.upper() == "MSE":
      resname = "MET"
      if atom_name == " SE ":
        atom_name = " SD "
    updated_resname = modernize_rna_resname(resname)
    if common_residue_names_get_class(updated_resname) == "common_rna_dna":
      updated_atom = modernize_rna_atom_name(atom=atom_name)
    else:
      updated_atom = atom_name
    key = updated_atom+atom.pdb_label_columns()[4:5]+\
          updated_resname+atom.pdb_label_columns()[8:]
    i_seq_name_hash[atom.i_seq]=key
  return i_seq_name_hash

def build_i_seq_hash(pdb_hierarchy):
  name_i_seq_hash = dict()
  for atom in pdb_hierarchy.atoms():
    atom_name = atom.pdb_label_columns()[0:4]
    resname = atom.pdb_label_columns()[5:8]
    if resname.upper() == "MSE":
      resname = "MET"
      if atom_name == " SE ":
        atom_name == " SD "
    updated_resname = modernize_rna_resname(resname)
    if common_residue_names_get_class(updated_resname) == "common_rna_dna":
      updated_atom = modernize_rna_atom_name(atom=atom_name)
    else:
      updated_atom = atom_name
    key = updated_atom+atom.pdb_label_columns()[4:5]+\
          updated_resname+atom.pdb_label_columns()[8:]
    name_i_seq_hash[key]=atom.i_seq
  return name_i_seq_hash

def build_xyz_hash(pdb_hierarchy):
  name_xyz_hash = dict()
  for atom in pdb_hierarchy.atoms():
    name_xyz_hash[atom.pdb_label_columns()]=atom.xyz
  return name_xyz_hash

def build_resid_hash(pdb_hierarchy):
  resid_hash = dict()
  for rg in pdb_hierarchy.residue_groups():
    resid = rg.resseq_as_int()
    for atom in rg.atoms():
      resid_hash[atom.i_seq]=resid
  return resid_hash

def build_i_seq_xyz_hash(pdb_hierarchy):
  i_seq_xyz_hash = dict()
  for atom in pdb_hierarchy.atoms():
    i_seq_xyz_hash[atom.i_seq] = atom.xyz
  return i_seq_xyz_hash

def build_element_hash(pdb_hierarchy):
  i_seq_element_hash = dict()
  for atom in pdb_hierarchy.atoms():
    i_seq_element_hash[atom.i_seq]=atom.element
  return i_seq_element_hash

def build_cbetadev_hash(pdb_hierarchy):
  cb = cbetadev()
  cbetadev_hash = dict()
  cbeta_out = cb.analyze_pdb(hierarchy=pdb_hierarchy)
  for line in cbeta_out[0].splitlines():
    temp = line.split(':')
    dev = temp[5]
    if dev == "dev":
      continue
    key = temp[1].upper()+temp[2].upper()+temp[3]+temp[4].rstrip()
    cbetadev_hash[key] = dev
  return cbetadev_hash

def build_chain_hash(pdb_hierarchy):
  chain_hash = dict()
  for chain in pdb_hierarchy.chains():
    for atom in chain.atoms():
      chain_hash[atom.i_seq] = chain.id
  return chain_hash

def build_segid_hash(pdb_hierarchy):
  segid_hash = dict()
  for atom in pdb_hierarchy.atoms():
    segid_hash[atom.i_seq] = atom.segid
  return segid_hash

def build_sym_atom_hash(pdb_hierarchy):
  sym_atom_hash = dict()
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        for residue in conformer.residues():
          if residue.resname.upper() in ['ASP', 'GLU', 'PHE', 'TYR']:
            if residue.resname.upper() == 'ASP':
              atom1 = ' OD1'
              atom2 = ' OD2'
            elif residue.resname.upper() == 'GLU':
              atom1 = ' OE1'
              atom2 = ' OE2'
            elif residue.resname.upper() in ['PHE', 'TYR']:
              atom1 = ' CD1'
              atom2 = ' CD2'
            atom1_i_seq = None
            atom2_i_seq = None
            for atom in residue.atoms():
              if atom.name == atom1:
                atom1_i_seq = atom.i_seq
              elif atom.name == atom2:
                atom2_i_seq = atom.i_seq
            if atom1_i_seq != None and atom2_i_seq != None:
              sym_atom_hash[atom1_i_seq] = atom2_i_seq
              sym_atom_hash[atom2_i_seq] = atom1_i_seq
  return sym_atom_hash

def angle_distance(angle1, angle2):
  distance = math.fabs(angle1 - angle2)
  if distance > 180.0:
    distance -= 360.0
  return math.fabs(distance)

def get_angle_average(angles):
  n_angles = len(angles)
  sum = 0.0
  a1 = angles[0]
  if a1 > 180.0:
    a1 -= 360.0
  elif a1 < -180.0:
    a1 += 360.0
  sum += a1
  for angle in angles[1:]:
    a2 = angle
    if (a1 - a2) > 180.0:
      a2 += 360.0
    elif (a1 - a2) < -180.0:
      a2 -= 360.0
    sum += a2
  average = sum / n_angles
  return average

def _ssm_align(reference_chain,
               moving_chain):
  ssm = ccp4io_adaptbx.SecondaryStructureMatching(
          reference=reference_chain,
          moving=moving_chain)
  ssm_alignment = ccp4io_adaptbx.SSMAlignment.residue_groups(match=ssm)
  return ssm, ssm_alignment

def _alignment(pdb_hierarchy,
               params,
               sequences,
               padded_sequences,
               structures,
               log=None):
  if(log is None): log = sys.stdout
  res_match_hash = {}
  model_mseq_res_hash = {}
  ref_mseq_res_hash = {}
  model_seq = sequences[0]
  model_seq_padded = padded_sequences[0]
  model_structures = structures[0]
  ref_seq = sequences[1]
  ref_seq_padded = padded_sequences[1]
  ref_structures = structures[1]
  for struct in model_structures:
    model_mseq_res_hash[struct.i_seq] = \
      struct.rg.atoms()[0].pdb_label_columns()[4:]
  for struct in ref_structures:
    ref_mseq_res_hash[struct.i_seq] = \
      struct.rg.atoms()[0].pdb_label_columns()[4:]
  if model_seq == ref_seq:
    pg = mmtbx.alignment.pairwise_global(
           model_seq,
           ref_seq)
  else:
    pg = mmtbx.alignment.pairwise_global(
           model_seq_padded,
           ref_seq_padded)
  offset_i = 0
  offset_j = 0
  i = 0
  j = 0
  seq_j = pg.result2[j]
  for seq_i in pg.result1:
    seq_j = pg.result2[j]
    if seq_i == seq_j and seq_i != 'X' and seq_j != 'X':
      res_match_hash[model_mseq_res_hash[i-offset_i]] = \
        ref_mseq_res_hash[j-offset_j]
      i += 1
      j += 1
    else:
      if seq_i == 'X' and seq_j == 'X':
        i += 1
        j += 1
        offset_i += 1
        offset_j += 1
      elif (seq_i == 'X' and seq_j == '-') or \
           (seq_i == '-' and seq_j == 'X'):
        i += 1
        j += 1
        offset_i += 1
        offset_j += 1
      elif seq_i == 'X':
        i += 1
        j += 1
        offset_i += 1
      elif seq_j == 'X':
        i += 1
        j += 1
        offset_j += 1
      elif seq_i == '-':
        i += 1
        j += 1
        offset_i += 1
      elif seq_j == '-':
        i += 1
        j += 1
        offset_j += 1
      else:
        i += 1
        j += 1
  return res_match_hash

def chain_from_selection(chain, selection):
  from iotbx.pdb.hierarchy import new_hierarchy_from_chain
  new_hierarchy = new_hierarchy_from_chain(chain=chain).select(selection)
  print dir(new_hierarchy)

def hierarchy_from_selection(pdb_hierarchy, selection, log):
  import iotbx.pdb.hierarchy
  temp_hierarchy = pdb_hierarchy.select(selection)
  altloc = None
  hierarchy = iotbx.pdb.hierarchy.root()
  model = iotbx.pdb.hierarchy.model()
  for chain in temp_hierarchy.chains():
    for conformer in chain.conformers():
      if not conformer.is_protein() and not conformer.is_na():
        continue
      elif altloc is None or conformer.altloc == altloc:
        model.append_chain(chain.detached_copy())
        altloc = conformer.altloc
      else:
        print >> log, \
        "* Multiple alternate conformations found, using altid %s *" \
        % altloc
        continue
  if len(model.chains()) != 1:
    raise Sorry("more than one chain in selection")
  hierarchy.append_model(model)
  return hierarchy

def extract_sequence_and_sites(pdb_hierarchy, selection):
  seq = []
  result = []
  padded_seq = []
  last_resseq = 0
  counter = 0
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      is_na = False
      for conformer in chain.conformers():
        if conformer.is_na():
          is_na = True
      for rg in chain.residue_groups():
        if(len(rg.unique_resnames())==1):
          resname = rg.unique_resnames()[0]
          if is_na:
            olc = get_nucleic_acid_one_letter_code(resname)
          else:
            olc= \
            amino_acid_codes.one_letter_given_three_letter.get(resname,"X")
          atoms = rg.atoms()
          i_seqs = get_i_seqs(atoms)
          if(olc!="X") and is_residue_in_selection(i_seqs, selection):
            seq.append(olc)
            resseq = rg.resseq_as_int()
            if (resseq > (last_resseq + 1)) :
              for x in range(resseq - last_resseq - 1) :
                padded_seq.append('X')
            last_resseq = resseq
            result.append(group_args(i_seq = counter, rg = rg))
            padded_seq.append(olc)
            counter += 1
  return "".join(seq), "".join(padded_seq), result

def get_nucleic_acid_one_letter_code(resname):
  olc=amino_acid_codes.one_letter_given_three_letter.get(resname,"X")
  if olc != "X":
    return "X"
  if resname[0:2] == "  ":
    return resname[2]
  elif resname[0] == " " and (resname[1] == "D" or resname[1] == "d"):
    return resname[2]
  else:
    return resname[0]

def get_unique_segid(chain):
  segid = None
  for atom in chain.atoms():
    if segid is None:
      segid = atom.segid
    elif segid != atom.segid:
      return None
  return segid
