from scitbx import matrix
from scitbx.math import superpose
import iotbx.pdb
from cctbx.array_family import flex
from mmtbx.monomer_library import idealized_aa      
from libtbx.utils import Sorry


alpha_pdb_str = """\
ATOM      1  N   ALA     2       1.643  -2.366  -1.408  1.00
ATOM      3  CA  ALA     2       1.280  -3.608  -2.069  1.00
ATOM      6  CB  ALA     2       1.361  -4.762  -1.068  1.00
ATOM     10  C   ALA     2      -0.114  -3.466  -2.684  1.00
ATOM     11  O   ALA     2      -0.327  -3.827  -3.840  1.00
ATOM     12  N   ALA     3      -1.028  -2.938  -1.882  1.00
ATOM     14  CA  ALA     3      -2.395  -2.743  -2.333  1.00
ATOM     17  CB  ALA     3      -3.228  -2.150  -1.194  1.00
ATOM     21  C   ALA     3      -2.396  -1.855  -3.579  1.00
ATOM     22  O   ALA     3      -3.059  -2.167  -4.567  1.00
"""

alpha310_pdb_str = """\
ATOM      1  N   ALA    1       -1.204  -0.514   0.643  1.0
ATOM      1  CA  ALA    1        0.000   0.000   0.000  1.0
ATOM      1  CB  ALA    1        0.870   0.757   1.006  1.0
ATOM      1  C   ALA    1        0.804  -1.124  -0.644  1.0
ATOM      1  O   ALA    1        1.628  -0.884  -1.526  1.0
ATOM      1  N   ALA    2        0.559  -2.352  -0.197  1.0
ATOM      1  CA  ALA    2        1.260  -3.515  -0.728  1.0
ATOM      1  CB  ALA    2        0.743  -4.801  -0.079  1.0
ATOM      1  C   ALA    2        1.116  -3.602  -2.244  1.0
ATOM      1  O   ALA    2        1.905  -4.266  -2.915  1.0
"""

beta_pdb_str = """\
ATOM      1  N   ALA    1      -1.204  -0.514   0.643  1.0
ATOM      1  CA  ALA    1       0.000   0.000   0.000  1.0
ATOM      1  CB  ALA    1      -0.000   1.530   0.000  1.0
ATOM      1  C   ALA    1       1.258  -0.522   0.685  1.0
ATOM      1  O   ALA    1       1.316  -0.616   1.911  1.0
ATOM      1  N   ALA    2       2.265  -0.860  -0.115  1.0
ATOM      1  CA  ALA    2       3.524  -1.373   0.412  1.0
ATOM      1  CB  ALA    2       3.770  -2.802  -0.076  1.0
ATOM      1  C   ALA    2       4.694  -0.478   0.018  1.0
ATOM      1  O   ALA    2       4.758   0.018  -1.107  1.0
"""


def get_r_t_matrices_from_structure(pdb_str):
  """ Return rotation and translation matrices for the ideal structure.        

  The function determines r and t matrices from alingment of 1st and 2nd 
  residues of the structure passed in pdb_str. 
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  conformer = pdb_hierarchy.models()[0].chains()[0].conformers()[0]
  residues = conformer.residues()
  fixed_sites = flex.vec3_double()
  moving_sites = flex.vec3_double()
  main_chain_atoms = ["N","CA","C","O"]
  if len(residues)>=2:
    for (r, arr) in [(residues[0], fixed_sites), (residues[1], moving_sites)]:
      for a in r.atoms():
        if a.name.strip() in main_chain_atoms:
          arr.append(a.xyz)
  else:
    raise Sorry('pdb_str should contain at least 2 residues')
  lsq_fit_obj = superpose.least_squares_fit(reference_sites = fixed_sites, 
                                            other_sites = moving_sites)
  return lsq_fit_obj.r, lsq_fit_obj.t

def make_ss_structure_from_seq(pdb_str, seq):    
  """ Return pdb.hierarchy with secondary structure according to sequence.        
  
  pdb_str - "ideal" structure at least 2 residues long.
  seq - string with sequence (one-letter codes)
  
  worth putting asserts on pdb_str not to be empty and len(seq)>1 
  """

  if len(seq)<1:
    raise Sorry('seq should contain at least one residue.')
  r, t = get_r_t_matrices_from_structure(pdb_str)

  aac = iotbx.pdb.amino_acid_codes.three_letter_given_one_letter    
  ideal_res_dict = idealized_aa.residue_dict()                      
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)        
  pdb_hierarchy = pdb_inp.construct_hierarchy()                     
  chain = pdb_hierarchy.models()[0].chains()[0]
  current_ala_ag = chain.residue_groups()[0].atom_groups()[0]
  new_chain = iotbx.pdb.hierarchy.chain(id=" ")
  new_chain.pre_allocate_residue_groups(number_of_additional_residue_groups=\
                                                                      len(seq))

  i=1
  for c in seq:
    # put ALA
    rg = iotbx.pdb.hierarchy.residue_group(resseq="%i" % i, icode="")
    new_chain.append_residue_group(residue_group=rg)
    ag_to_place = current_ala_ag.detached_copy()
    rg.append_atom_group(atom_group=ag_to_place)
    current_ala_ag.atoms().set_xyz(
                          r.elems*current_ala_ag.atoms().extract_xyz()+t.elems)
    i += 1
    if c == 'A':
      continue
    ag_to_place.resname = aac[c]
    if c == 'G':
      for a in ag_to_place.atoms():
        if a.name.strip() == "CB":
          ag_to_place.remove_atom(atom=a)
          break
      continue
    # align residue from ideal_res_dict to just placed ALA (ag_to_place)
    fixed_sites = flex.vec3_double()
    moving_sites = flex.vec3_double()
    reper_atoms = ["CB","CA", "N"]
    for (ag, arr) in [(ag_to_place, fixed_sites), 
                      (ideal_res_dict[aac[c].lower()], moving_sites)]:
      for a in ag.atoms():
        if a.name.strip() in reper_atoms:
          arr.append(a.xyz)
    lsq_fit_obj = superpose.least_squares_fit(reference_sites = fixed_sites, 
                                              other_sites = moving_sites)
    ideal_correct_ag = ideal_res_dict[aac[c].lower()].\
                      models()[0].\
                      chains()[0].\
                      residue_groups()[0].\
                      atom_groups()[0].detached_copy()
    ideal_correct_ag.atoms().set_xyz(
      lsq_fit_obj.r.elems*ideal_correct_ag.atoms().extract_xyz()+\
      lsq_fit_obj.t.elems)
    ag_to_place.pre_allocate_atoms(number_of_additional_atoms=\
                                                len(ideal_correct_ag.atoms())-5)
    for a in ideal_correct_ag.atoms():
      if a.name.strip() not in ["N","CA","C","O", "CB"]:
        ag_to_place.append_atom(atom=a.detached_copy())
		
  new_pdb_h = iotbx.pdb.hierarchy.new_hierarchy_from_chain(new_chain)
  new_pdb_h.atoms().reset_i_seq()
  new_pdb_h.atoms().reset_serial()
  return new_pdb_h


def beta():
  pdb_hierarchy = make_ss_structure_from_seq(beta_pdb_str, "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_beta_seq.pdb")


def alpha_310():
  pdb_hierarchy = make_ss_structure_from_seq(alpha310_pdb_str, "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_helix310_seq.pdb")
	

def alpha():
  pdb_hierarchy = make_ss_structure_from_seq(alpha_pdb_str, "ACEDGFIHKMLNQPSRTWVY")
  pdb_hierarchy.write_pdb_file(file_name = "o_helix_seq.pdb")


def run():
  beta()
  alpha_310()
  alpha()
  
  
if (__name__ == "__main__"):
  run()
