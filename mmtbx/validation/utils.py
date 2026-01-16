
from __future__ import absolute_import, division, print_function
from libtbx import easy_pickle
from libtbx import group_args
from libtbx.utils import Sorry
from collections import defaultdict
from iotbx.pdb import common_residue_names_get_class
import os.path
import math
import sys
from six.moves import range

def export_ramachandran_distribution(n_dim_table, scale_factor=0.25):
  """
  Convert a MolProbity Ramachandran distribution to a format suitable for
  display using matplotlib (see wxtbx.plots).
  """
  import numpy
  z = n_dim_table.lookupTable
  z_grid = [ [ z[i + (j * 180)] for j in range(180) ]
                          for i in range(180) ]
  npz = numpy.array(z_grid)
  return npz ** scale_factor

def export_rotamer_distribution(n_dim_table, scale_factor=0.5):
  """
  Convert a MolProbity rotamer distribution to a format suitable for
  display using matplotlib (see wxtbx.plots).  Will reduce dimensionality to
  2 if necessary.
  """
  import numpy
  z = n_dim_table.lookupTable
  n_dim = n_dim_table.nDim
  assert n_dim >= 2
  x_offset = 1
  for nbins in n_dim_table.nBins[1:] :
    x_offset *= nbins
  y_width = 1
  if n_dim > 2 :
    for nbins in n_dim_table.nBins[2:] :
      y_width *= nbins
  z_grid = [ [] for i in range(n_dim_table.nBins[1]) ]
  for i in range(n_dim_table.nBins[0]):
    for j in range(n_dim_table.nBins[1]):
      z_total = 0
      for k in range(y_width):
        z_total += z[(i * x_offset) + (j * y_width) + k]
      z_grid[j].append(z_total)
  npz = numpy.array(z_grid)
  return npz ** scale_factor

def get_rotarama_data(residue_type=None, pos_type=None, db="rama",
    convert_to_numpy_array=False):
  from mmtbx.rotamer import ramachandran_eval
  from mmtbx.rotamer.rotamer_eval import find_rotarama_data_dir
  # backwards compatibility
  if (pos_type == "proline") : pos_type = "trans-proline"
  if (pos_type == "prepro") : pos_type = "pre-proline"
  assert (pos_type in ["general", "cis-proline", "trans-proline", "glycine",
    "isoleucine or valine", "pre-proline",None])
  assert (db in ["rama", "rota"])
  assert (residue_type is not None) or (pos_type is not None)
  if pos_type is not None :
    residue_type = ramachandran_eval.aminoAcids_8000[pos_type]
  if residue_type.lower() in ["phe", "tyr"] :
    residue_type = "phetyr"
  assert (residue_type is not None)
  rama_data_dir = find_rotarama_data_dir()
  if (db == "rama"):
    pkl_file = "%s8000-%s.pickle" % (db, residue_type)
  else :
    pkl_file = "%s8000-%s.pickle" % (db, residue_type.lower())
  ndt = easy_pickle.load(os.path.join(rama_data_dir, pkl_file))
  if convert_to_numpy_array :
    if (db == "rama"):
      return export_ramachandran_distribution(ndt)
    else :
      return export_rotamer_distribution(ndt)
  else :
    return ndt

def decode_atom_str(atom_id):
  chain_id = atom_id[8:10].strip()
  if (chain_id == ""):
    chain_id = " "
  return group_args(
    name = atom_id[0:4],
    altloc = atom_id[4],
    resname = atom_id[5:8],
    chain_id = chain_id,
    resid = atom_id[10:],
    resseq = atom_id[10:-1].strip())

def find_sequence_mismatches(pdb_hierarchy,
                              sequences,
                              assume_same_order=True,
                              expected_sequence_identity=0.8,
                              log=sys.stdout):
  chains = pdb_hierarchy.models()[0].chains()
  chain_ids = []
  actual_seqs = []
  expected_seqs = []
  if (len(chains) != len(sequences)) or (not assume_same_order):
    print("Can't determine sequence->chain mapping autoamtically", file=log)
    print("Running sequence alignments. . .", file=log)
    from mmtbx.alignment import pairwise_global_wrapper
    for chain in chains :
      chain_seq = chain.as_padded_sequence()
      actual_seqs.append(chain_seq)
      chain_ids.append(chain.id)
      best_identity = 0
      best_sequence = None
      for sequence in sequences :
        pg = pairwise_global_wrapper(chain_seq, sequence)
        identity = pg.calculate_sequence_identity()
        if (identity >= expected_sequence_identity):
          if (identity >= best_identity):
            best_identity = identity
            best_sequence = sequence
      expected_seqs.append(best_sequence)
  mismatches = []

def molprobity_score(clashscore, rota_out, rama_fav):
  """
  Calculate the overall Molprobity score, as described here:
    http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2877634/?tool=pubmed
    http://kinemage.biochem.duke.edu/suppinfo/CASP8/methods.html
  """
  if clashscore is not None and rota_out is not None and rama_fav is not None \
       and (clashscore >= 0) and (rota_out >= 0) and (rama_fav >= 0):
    rama_iffy = 100. - rama_fav
    mpscore = (( 0.426 * math.log(1 + clashscore) ) +
             ( 0.33 * math.log(1 + max(0, rota_out - 1)) ) +
             ( 0.25 * math.log(1 + max(0, rama_iffy - 2)) )) + 0.5
  else :
    return -1 # FIXME prevents crashing on RNA and None in inputs
  return mpscore

from typing import Dict, Any, Optional

def calculate_overall_residue_quality_score(
    residue_data: Dict[str, Any],
    weights: Optional[Dict[str, float]] = None
) -> Optional[float]:
    """
    Calculates an aggregated quality score for a single residue based on various validation metrics.

    Args:
        residue_data (Dict[str, Any]): A dictionary containing all the extracted
                                         validation metrics for a single residue.
                                         Expected keys should match those in ResidueResult,
                                         e.g., 'ramalyze_category', 'is_cbeta_outlier', etc.
        weights (Optional[Dict[str, float]]): A dictionary of weights for each validation category.
                                              If None, default weights will be used.

    Returns:
        Optional[float]: The aggregated quality score (0-1), or None if no applicable metrics.
    """

    # --- Default Weights (can be overridden by the 'weights' argument) ---
    DEFAULT_WEIGHTS = {
        'RAMA': 0.10,
        'ROTA': 0.10,
        'CBETA': 0.10,
        'CABLAM': 0.20,
        'OMEGA': 0.20,
        'RNA_SUITE': 0.35,
        'RNA_PUCKER': 0.35,
        'BONDS': 0.05,
        'ANGLES': 0.05,
        'CHIRAL': 0.05,
        'CLASH': 0.15,
    }

    current_weights = weights if weights is not None else DEFAULT_WEIGHTS

    total_quality_points = 0.0
    total_possible_points = 0.0

    # --- Helper to safely get data ---
    def get_data(key, default=None):
        return residue_data.get(key, default)

    # --- 1. Ramachandran Contribution ---
    s_rama = 1.0
    if get_data('ramalyze_category') == 'outlier':
        s_rama = 0.0
    elif get_data('ramalyze_category') == 'allowed':
        s_rama = 0.8
    if get_data('ramalyze_type') not in ['not_applicable', 'not_evaluated']:
        total_quality_points += s_rama * current_weights['RAMA']
        total_possible_points += current_weights['RAMA']

    # --- 2. Rotamer Contribution ---
    s_rota = 1.0
    if get_data('rotalyze_category') == 'outlier':
        s_rota = 0.0
    elif get_data('rotalyze_category') == 'allowed':
        s_rota = 0.8
    if get_data('rotalyze_category') != 'not_evaluated':
        total_quality_points += s_rota * current_weights['ROTA']
        total_possible_points += current_weights['ROTA']

    # --- 3. C-beta Deviation Contribution ---
    s_cbeta = 1.0
    if get_data('is_cbeta_outlier'):
        s_cbeta = 0.0
    if not get_data('is_glycine', False): # Glycine doesn't have a C-beta
        total_quality_points += s_cbeta * current_weights['CBETA']
        total_possible_points += current_weights['CBETA']

    # --- 4. CaBLAM Contribution ---
    s_cablam = 1.0
    if get_data('cablam_outlier_type') == 'outlier':
        s_cablam = 0.0
    elif get_data('cablam_outlier_type') == 'disfavored':
        s_cablam = 0.5
    if get_data('cablam_outlier_type') != 'not_evaluated':
        total_quality_points += s_cablam * current_weights['CABLAM']
        total_possible_points += current_weights['CABLAM']

    # --- 5. Omega Angle Contribution ---
    s_omega = 1.0
    if get_data('omega_type') == 'twisted':
        s_omega = 0.0
    elif get_data('omega_type') == 'cis' and not get_data('is_proline'):
        s_omega = 0.0
    if get_data('omega_res_type') != 'not_applicable':
        total_quality_points += s_omega * current_weights['OMEGA']
        total_possible_points += current_weights['OMEGA']

    # --- 6. RNA Suite Contribution ---
    s_rna_suite = 1.0
    if get_data('is_rna_residue', False):
        if get_data('is_rna_suite_outlier', False):
            s_rna_suite = 0.1
        total_quality_points += s_rna_suite * current_weights['RNA_SUITE']
        total_possible_points += current_weights['RNA_SUITE']

    # --- 7. RNA Pucker Contribution ---
    s_rna_pucker = 1.0
    if get_data('is_rna_residue', False):
        if get_data('is_rna_pucker_outlier', False):
            s_rna_pucker = 0.1
        total_quality_points += s_rna_pucker * current_weights['RNA_PUCKER']
        total_possible_points += current_weights['RNA_PUCKER']

    # --- 8. Bond Lengths Contribution ---
    s_bonds = 1.0
    if get_data('num_bond_outliers_res', 0) > 0:
        s_bonds = 0.0 # Simple penalty for any bond outlier
    total_quality_points += s_bonds * current_weights['BONDS']
    total_possible_points += current_weights['BONDS']

    # --- 9. Bond Angles Contribution ---
    s_angles = 1.0
    if get_data('num_angle_outliers_res', 0) > 0:
        s_angles = 0.0 # Simple penalty for any angle outlier
    total_quality_points += s_angles * current_weights['ANGLES']
    total_possible_points += current_weights['ANGLES']

    # --- 10. Chirality Contribution ---
    s_chiral = 1.0
    if get_data('num_chiral_outliers_res', 0) > 0:
        s_chiral = 0.0
    total_quality_points += s_chiral * current_weights['CHIRAL']
    total_possible_points += current_weights['CHIRAL']

    # --- 11. Steric Clash Contribution ---
    s_clash = 1.0
    num_bad_clashes = get_data('num_bad_clashes_res', 0) # Assumes you are still populating this simple count during parsing
    if num_bad_clashes > 0:
        s_clash = max(0.0, 1.0 - (0.5 * num_bad_clashes))
        if num_bad_clashes > 2: s_clash = 0.0 # Critical failure threshold
    total_quality_points += s_clash * current_weights['CLASH']
    total_possible_points += current_weights['CLASH']


    # --- Final Normalization ---
    if total_possible_points > 0:
        final_score = total_quality_points / total_possible_points
        return round(final_score, 4) # Round for consistency
    else:
        # If no metrics were applicable (e.g., a pure ligand without any standard validations)
        return None # Or 1.0 if you prefer a perfect score for non-evaluated items

def use_segids_in_place_of_chainids(hierarchy, strict=False):
  use_segids = False
  for model in hierarchy.models():
    for chain in model.chains():
      if chain.id in [' ', '  ']:
        cur_segid = None
        for atom in chain.atoms():
          # new as of 20150203
          if atom.segid not in ['    ', '']:
            return True
          # It makes no sense to require indentical segID for
          # Chains with blank chainID. This was commented out by BJH on 20150203
          #if cur_segid is None:
          #  cur_segid = atom.segid
          #if atom.segid not in ['    ', '']:
          #  if atom.segid != cur_segid:
          #    if strict:
          #      raise Sorry("Chains with blank chainID may not have multiple"+
          #                  " segID values")
          #    else:
          #      return False
        #if len(cur_segid.strip()) > 0:
        #  use_segids = True
        #else:
        #  return False
  return use_segids

#this function assumes that use_segids_in_place_of_chainids() is True
def get_segid_as_chainid(chain):
  for atom in chain.atoms():
    return atom.segid

def get_rna_backbone_dihedrals(processed_pdb_file,
      geometry=None, pdb_hierarchy=None):
  # at present, this will only return measurements for angles arising from
  # atoms with altloc ' ' or altloc 'A'
  # TO-DO: extend to more alternates JJH 140108
  from cctbx import geometry_restraints
  bb_dihedrals = defaultdict(dict)
  formatted_out = []
  alt_tracker = {}
  if (processed_pdb_file is not None):
    sites_cart = processed_pdb_file.all_chain_proxies.sites_cart
    geometry = processed_pdb_file.geometry_restraints_manager()
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  else :
    assert (not None in [geometry, pdb_hierarchy])
    sites_cart = pdb_hierarchy.atoms().extract_xyz()
  dihedral_proxies = geometry.dihedral_proxies
  i_seq_name_hash = build_name_hash(pdb_hierarchy=pdb_hierarchy)

  def is_blank_or_alt_a(proxy):
    for i in proxy.i_seqs:
       alt = i_seq_name_hash[i][4:5]
       if alt not in [' ', 'A']:
         return False
    return True

  for dp in dihedral_proxies:
    atoms = []
    debug_key = ""
    invert_sign = False
    dp.sort_i_seqs()
    for i in dp.i_seqs:
      atoms.append(i_seq_name_hash[i][0:4].strip())
      debug_key = debug_key+i_seq_name_hash[i]
    if len(atoms) != 4:
      continue
    name = match_dihedral_to_name(atoms=atoms)
    #handle dihedral equivalences
    if name == None:
      inverted_atoms = get_inverted_atoms(atoms=atoms, improper=False)
      name = match_dihedral_to_name(atoms=inverted_atoms)
      if name == None:
        inverted_atoms = get_inverted_atoms(atoms=atoms, improper=True)
        name = match_dihedral_to_name(atoms=inverted_atoms)
        if name is not None:
          invert_sign = True
    if (name is not None) and (is_blank_or_alt_a(dp)):
      restraint = geometry_restraints.dihedral(
                                               sites_cart=sites_cart,
                                               proxy=dp)
      key = i_seq_name_hash[dp.i_seqs[1]][4:]
      if alt_tracker.get(key[1:]) is None:
        alt_tracker[key[1:]] = []
      if key[0:1] not in alt_tracker[key[1:]]:
        alt_tracker[key[1:]].append(key[0:1])
      bb_dihedrals[key][name] = restraint.angle_model
      if invert_sign:
        bb_dihedrals[key][name] = bb_dihedrals[key][name] * -1.0
  for key in list(bb_dihedrals.keys()):
    altloc = key[0:1]
    resname = key[1:4]
    chainID = key[4:6]
    resnum = key[6:10]
    i_code = key[10:]
    if 'A' in alt_tracker[key[1:]]:
      if altloc != 'A':
        continue
    if bb_dihedrals[key].get('alpha') is not None:
      alpha = "%.3f" % bb_dihedrals[key]['alpha']
    # FIXME will the lookup below ever actually work?
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('alpha') is not None:
      alpha = "%.3f" % bb_dihedrals[' '+key[1:]]['alpha']
    else:
      alpha = '__?__'
    if bb_dihedrals[key].get('beta') is not None:
      beta = "%.3f" % bb_dihedrals[key]['beta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('beta') is not None:
      beta = "%.3f" % bb_dihedrals[' '+key[1:]]['beta']
    else:
      beta = '__?__'
    if bb_dihedrals[key].get('gamma') is not None:
      gamma = "%.3f" % bb_dihedrals[key]['gamma']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('gamma') is not None:
      gamma = "%.3f" % bb_dihedrals[' '+key[1:]]['gamma']
    else:
      gamma = '__?__'
    if bb_dihedrals[key].get('delta'):
      delta = "%.3f" % bb_dihedrals[key]['delta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('delta') is not None:
      delta = "%.3f" % bb_dihedrals[' '+key[1:]]['delta']
    else:
      delta = '__?__'
    if bb_dihedrals[key].get('epsilon'):
      epsilon = "%.3f" % bb_dihedrals[key]['epsilon']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('epsilon') is not None:
      epsilon = "%.3f" % bb_dihedrals[' '+key[1:]]['epsilon']
    else:
      epsilon = '__?__'
    if bb_dihedrals[key].get('zeta'):
      zeta = "%.3f" % bb_dihedrals[key]['zeta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('zeta') is not None:
      zeta = "%.3f" % bb_dihedrals[' '+key[1:]]['zeta']
    else:
      zeta = '__?__'
    eval = "%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s" \
           % (" ",
              "1",
              chainID,
              resnum,
              i_code,
              altloc,
              resname,
              alpha,
              beta,
              gamma,
              delta,
              epsilon,
              zeta)
    formatted_out.append(eval)
  formatted_out.sort()
  backbone_dihedrals = ""
  for line in formatted_out:
    backbone_dihedrals += line+'\n'
  return backbone_dihedrals

def get_inverted_atoms(atoms, improper=False):
  temp = []
  if not improper:
    temp.append(atoms[3])
    temp.append(atoms[2])
    temp.append(atoms[1])
    temp.append(atoms[0])
  else:
    temp.append(atoms[3])
    temp.append(atoms[1])
    temp.append(atoms[2])
    temp.append(atoms[0])
  return temp

def match_dihedral_to_name(atoms):
  name = None
  alpha = ["O3'","P","O5'","C5'"]
  beta = ["P","O5'","C5'","C4'"]
  gamma = ["O5'","C5'","C4'","C3'"]
  delta = ["C5'","C4'","C3'","O3'"]
  epsilon = ["C4'","C3'","O3'","P"]
  zeta = ["C3'","O3'","P","O5'"]
  if atoms == alpha:
    name = "alpha"
  elif atoms == beta:
    name = "beta"
  elif atoms == gamma:
    name = "gamma"
  elif atoms == delta:
    name = "delta"
  elif atoms == epsilon:
    name = "epsilon"
  elif atoms == zeta:
    name = "zeta"
  return name

def build_name_hash(pdb_hierarchy):
  i_seq_name_hash = dict()
  for atom in pdb_hierarchy.atoms():
    i_seq_name_hash[atom.i_seq]=atom.pdb_label_columns()
  return i_seq_name_hash

def get_mmtype_from_resname(resname):
  class_string = common_residue_names_get_class(resname)
  if "amino_acid" in class_string:
    return "PROTEIN"
  elif "rna_dna" in class_string:
    return "NA"
  else:
    return "OTHER"

def exercise():
  from libtbx.test_utils import approx_equal
  try :
    import numpy
  except ImportError :
    test_numpy = False
    print("Numpy not installed, will skip array conversion.")
  else :
    test_numpy = True
  # ramachandran
  z_data = get_rotarama_data(pos_type="general",
    convert_to_numpy_array=test_numpy)
  z_data = get_rotarama_data(pos_type="pre-proline",
    convert_to_numpy_array=test_numpy)
  # rotamer
  z_data = get_rotarama_data(residue_type="arg",
    db="rota",
    convert_to_numpy_array=test_numpy)
  z_data = get_rotarama_data(residue_type="phe",
    db="rota",
    convert_to_numpy_array=test_numpy)
  atom_info = decode_atom_str(" OD2 ASP A  14L")
  assert (atom_info.name == " OD2") and (atom_info.resname == "ASP")
  assert (atom_info.altloc == " ") and (atom_info.chain_id == "A")
  assert (atom_info.resid == "  14L") and (atom_info.resseq == "14")
  mpscore = molprobity_score(48.0, 9.95, 86.44) # 2hr0
  assert approx_equal(mpscore, 3.55, eps=0.01)
  mpscore = molprobity_score(215.8, 17.99, 52.18) # 3mku
  assert approx_equal(mpscore, 4.71, eps=0.01)

if __name__ == "__main__" :
  exercise()
  print("OK")
