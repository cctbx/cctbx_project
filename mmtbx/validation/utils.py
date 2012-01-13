
from libtbx import easy_pickle
from libtbx import group_args
import os.path
import math
import sys

def export_ramachandran_distribution (n_dim_table, scale_factor=0.25) :
  import numpy
  z = n_dim_table.lookupTable
  z_grid = [ [ z[i + (j * 180)] for j in range(180) ]
                          for i in range(180) ]
  npz = numpy.array(z_grid)
  return npz ** scale_factor

def export_rotamer_distribution (n_dim_table, scale_factor=0.5) :
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
  for i in range(n_dim_table.nBins[0]) :
    for j in range(n_dim_table.nBins[1]) :
      z_total = 0
      for k in range(y_width) :
        z_total += z[(i * x_offset) + (j * y_width) + k]
      z_grid[j].append(z_total)
  npz = numpy.array(z_grid)
  return npz ** scale_factor

def get_rotarama_data (residue_type=None, pos_type=None, db="rama",
    convert_to_numpy_array=False) :
  from mmtbx.rotamer import ramachandran_eval
  from mmtbx.rotamer.rotamer_eval import find_rotarama_data_dir
  assert (pos_type in ["general", "proline", "glycine", "prepro", None])
  assert (db in ["rama", "rota"])
  assert (residue_type is not None) or (pos_type is not None)
  if pos_type is not None :
    residue_type = ramachandran_eval.aminoAcids[pos_type]
  if residue_type.lower() in ["phe", "tyr"] :
    residue_type = "phetyr"
  assert (residue_type is not None)
  rama_data_dir = find_rotarama_data_dir()
  pkl_file = "%s500-%s.pickle" % (db, residue_type.lower())
  ndt = easy_pickle.load(os.path.join(rama_data_dir, pkl_file))
  if convert_to_numpy_array :
    if (db == "rama") :
      return export_ramachandran_distribution(ndt)
    else :
      return export_rotamer_distribution(ndt)
  else :
    return ndt

def decode_atom_str (atom_id) :
  chain_id = atom_id[8:10].strip()
  if (chain_id == "") :
    chain_id = " "
  return group_args(
    name = atom_id[0:4],
    altloc = atom_id[4],
    resname = atom_id[5:8],
    chain_id = chain_id,
    resid = atom_id[10:],
    resseq = atom_id[10:-1].strip())

def find_sequence_mismatches (pdb_hierarchy,
                              sequences,
                              assume_same_order=True,
                              expected_sequence_identity=0.8,
                              log=sys.stdout) :
  chains = pdb_hierarchy.models()[0].chains()
  chain_ids = []
  actual_seqs = []
  expected_seqs = []
  if (len(chains) != len(sequences)) or (not assume_same_order) :
    print >> log, "Can't determine sequence->chain mapping autoamtically"
    print >> log, "Running sequence alignements. . ."
    from mmtbx.alignment import pairwise_global_wrapper
    for chain in chains :
      main_conf = chain.conformers()[0]
      chain_seq = main_conf.as_padded_sequence()
      actual_seqs.append(chain_seq)
      chain_ids.append(chain.id)
      best_identity = 0
      best_sequence = None
      for sequence in sequences :
        pg = pairwise_global_wrapper(chain_seq, sequence)
        identity = pg.calculate_sequence_identity()
        if (identity >= expected_sequence_identity) :
          if (identity >= best_identity) :
            best_identity = identity
            best_sequence = sequence
      expected_seqs.append(best_sequence)
  mismatches = []

def molprobity_score (clashscore, rota_out, rama_fav) :
  """
  Calculate the overall Molprobity score, as described here:
    http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2877634/?tool=pubmed
    http://kinemage.biochem.duke.edu/suppinfo/CASP8/methods.html
  """
  if (clashscore >= 0) and (rota_out >= 0) and (rama_fav >= 0) :
    rama_iffy = 100. - rama_fav
    mpscore = (( 0.426 * math.log(1 + clashscore) ) +
             ( 0.33 * math.log(1 + max(0, rota_out - 1)) ) +
             ( 0.25 * math.log(1 + max(0, rama_iffy - 2)) )) + 0.5
  else :
    return -1 # FIXME prevents crashing on RNA
  return mpscore

def exercise () :
  from libtbx.test_utils import approx_equal
  try :
    import numpy
  except ImportError :
    test_numpy = False
    print "Numpy not installed, will skip array conversion."
  else :
    test_numpy = True
  # ramachandran
  z_data = get_rotarama_data(pos_type="general",
    convert_to_numpy_array=test_numpy)
  z_data = get_rotarama_data(pos_type="prepro",
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
  print "OK"
