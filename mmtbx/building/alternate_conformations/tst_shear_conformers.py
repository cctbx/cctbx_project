
from __future__ import absolute_import, division, print_function
import libtbx.phil
from six.moves import cStringIO as StringIO

master_params = """
set_partial_occupancy = None
  .type = float
rmsd_min = 1.0
  .type = float
rigid_body_refine = False
  .type = bool
translation = 0.5
  .type = float
  .help = Specifies the radius for the optional translation search, in \
    Angstroms
translation_sampling = 10
  .type = int
  .optional = False
  .help = Number of distances within the translation search to sample.
torsion_search {
  include scope mmtbx.building.alternate_conformations.conformer_generation.torsion_search_params
}"""

def master_phil():
  return libtbx.phil.parse(master_params, process_includes=True)

def exercise():
  from mmtbx.building.alternate_conformations import conformer_generation
  from mmtbx.monomer_library import server
  import iotbx.pdb.hierarchy
  generate_inputs()
  params = master_phil().extract()
  params.torsion_search.chi_increment_degrees=90 # Does not make sense but lets
                                                 # test to run quick. With
                                                 # default value test never
                                                 # finishes.
  mon_lib_srv = server.server()
  pdb_in = iotbx.pdb.hierarchy.input(file_name="shear_frag_single.pdb")
  hierarchy = pdb_in.hierarchy
  pdb_atoms = hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  sites_cart = pdb_atoms.extract_xyz()
  xrs = pdb_in.input.xray_structure_simple()
  models = []
  prev_res = next_res = next_next_res = None
  for chain in hierarchy.only_model().chains():
    residue_groups = chain.residue_groups()
    n_rg = len(residue_groups) # should be 4
    for i_res, residue_group in enumerate(residue_groups):
      sites_orig = sites_cart.deep_copy()
      next_res = next_next_res = None
      if (i_res < (n_rg - 1)):
        next_res = residue_groups[i_res+1].atom_groups()[0]
      if (i_res < (n_rg - 2)):
        next_next_res = residue_groups[i_res+2].atom_groups()[0]
      atom_groups = residue_group.atom_groups()
      primary_conf = atom_groups[0]
      out = StringIO()
      confs = []
      tmp = list(conformer_generation.generate_single_residue_confs(
        atom_group=primary_conf,
        sites_cart=sites_cart.deep_copy(),
        mon_lib_srv=mon_lib_srv,
        params=params.torsion_search,
        prev_residue=prev_res,
        next_residue=next_res,
        next_next_residue=next_next_res,
        backrub=False,
        shear=True))
      for conf in tmp :
          conf.show_summary(out=out)
          confs.append(conf)
      prev_res = primary_conf
      if (confs is None):
        continue
      if (i_res == 1):
        assert ("""  A ILE   7     None     4.0      mt""")
      for conf in confs :
        sites_new = sites_cart.set_selected(conf.sites_selection,
          conf.sites_selected())
        pdb_atoms.set_xyz(sites_new)
        models.append(hierarchy.only_model().detached_copy())
  new_hierarchy = iotbx.pdb.hierarchy.root()
  for i_model, conf in enumerate(models):
    conf.id = str(i_model + 1)
    new_hierarchy.append_model(conf)
  open("shear_frag_naive_ensemble.pdb", "w").write(
    new_hierarchy.as_pdb_string())

def generate_inputs():
  import iotbx.pdb.hierarchy
  pdb_in = iotbx.pdb.hierarchy.input(source_info=None, pdb_string="""\
REMARK derived from shear region in crambin (1ejg)
ATOM     89  N   SER A   6       8.926   6.705  14.756  1.00  1.84           N
ATOM     90  CA  SER A   6       8.755   5.347  15.231  1.00  1.70           C
ATOM     91  C   SER A   6       8.832   4.429  14.008  1.00  1.84           C
ATOM     92  O   SER A   6       8.757   4.867  12.860  1.00  2.14           O
ATOM     93  CB  SER A   6       7.409   5.161  15.930  1.00  1.78           C
ATOM     94  OG  SER A   6       6.371   5.265  14.953  1.00  2.15           O
ATOM     95  H   SER A   6       8.700   6.916  13.796  1.00  2.20           H
ATOM     96  HA  SER A   6       9.567   5.098  15.928  1.00  2.03           H
ATOM     97  HB2 SER A   6       7.371   4.176  16.415  1.00  2.13           H
ATOM     98  HB3 SER A   6       7.277   5.932  16.702  1.00  2.13           H
ATOM     99  HG  SER A   6       5.560   4.870  15.300  1.00  3.21           H
ATOM    100  N   ILE A   7       8.936   3.129  14.302  1.00  2.20           N
ATOM    101  CA  ILE A   7       8.829   2.039  13.300  1.00  2.58           C
ATOM    103  C   ILE A   7       7.559   2.141  12.460  1.00  2.53           C
ATOM    105  O   ILE A   7       7.573   2.124  11.205  1.00  2.78           O
ATOM    107  CB  ILE A   7       8.928   0.644  13.932  1.00  3.24           C
ATOM    109  CG1 ILE A   7      10.309   0.385  14.573  1.00  5.58           C
ATOM    111  CG2 ILE A   7       8.510  -0.464  12.908  1.00  4.24           C
ATOM    113  CD1 ILE A   7      11.430   0.295  13.473  1.00  7.60           C
ATOM    115  H   ILE A   7       9.078   2.878  15.269  1.00  2.65           H
ATOM    116  HA  ILE A   7       9.714   2.123  12.623  1.00  3.11           H
ATOM    118  HB  ILE A   7       8.200   0.598  14.759  1.00  3.92           H
ATOM    120 HG12 ILE A   7      10.550   1.170  15.240  1.00  6.88           H
ATOM    122 HG13 ILE A   7      10.274  -0.579  15.106  1.00  6.88           H
ATOM    124 HG21 ILE A   7       7.863   0.004  12.146  1.00  6.50           H
ATOM    126 HG22 ILE A   7       7.955  -1.251  13.426  1.00  6.50           H
ATOM    128 HG23 ILE A   7       9.402  -0.873  12.434  1.00  6.50           H
ATOM    130 HD11 ILE A   7      12.389   0.055  13.960  1.00 11.57           H
ATOM    132 HD12 ILE A   7      11.499   1.217  12.922  1.00 11.57           H
ATOM    134 HD13 ILE A   7      11.196  -0.544  12.761  1.00 11.57           H
ATOM    136  N   VAL A   8       6.382   2.222  13.070  1.00  1.92           N
ATOM    138  CA  VAL A   8       5.099   2.259  12.380  1.00  2.32           C
ATOM    140  C   VAL A   8       5.208   3.386  11.373  1.00  2.63           C
ATOM    141  O   VAL A   8       4.712   3.302  10.238  1.00  2.67           O
ATOM    142  CB  VAL A   8       3.944   2.394  13.375  1.00  3.36           C
ATOM    144  CG1 VAL A   8       2.629   2.981  12.701  1.00  4.22           C
ATOM    146  CG2 VAL A   8       3.635   1.036  13.955  1.00  4.82           C
ATOM    148  H   VAL A   8       6.374   2.231  14.084  1.00  2.27           H
ATOM    150  HA  VAL A   8       4.948   1.364  11.845  1.00  2.65           H
ATOM    152  HB  VAL A   8       4.204   3.077  14.197  1.00  4.16           H
ATOM    154 HG11 VAL A   8       2.862   3.979  12.286  1.00  6.22           H
ATOM    156 HG12 VAL A   8       1.845   3.103  13.478  1.00  6.22           H
ATOM    158 HG13 VAL A   8       2.286   2.318  11.926  1.00  6.22           H
ATOM    160 HG21 VAL A   8       4.536   0.671  14.467  1.00  7.38           H
ATOM    162 HG22 VAL A   8       3.367   0.343  13.145  1.00  7.38           H
ATOM    164 HG23 VAL A   8       2.805   1.098  14.673  1.00  7.38           H
ATOM    166  N   ALA A   9       5.668   4.556  11.848  1.00  2.19           N
ATOM    167  CA  ALA A   9       5.637   5.732  10.997  1.00  1.92           C
ATOM    168  C   ALA A   9       6.496   5.510   9.742  1.00  1.64           C
ATOM    169  O   ALA A   9       6.087   5.887   8.640  1.00  1.98           O
ATOM    170  CB  ALA A   9       6.088   6.968  11.772  1.00  2.30           C
ATOM    171  H   ALA A   9       6.031   4.617  12.787  1.00  2.64           H
ATOM    172  HA  ALA A   9       4.599   5.892  10.675  1.00  2.30           H
ATOM    173  HB1 ALA A   9       5.450   7.098  12.658  1.00  3.45           H
ATOM    174  HB2 ALA A   9       6.005   7.855  11.128  1.00  3.45           H
ATOM    175  HB3 ALA A   9       7.132   6.841  12.089  1.00  3.45           H
END
""")
  xrs = pdb_in.input.xray_structure_simple(cryst1_substitution_buffer_layer=5)
  f = open("shear_frag.pdb", "w")
  f.write(pdb_in.hierarchy.as_pdb_string(crystal_symmetry=xrs))
  f.close()
  for chain in pdb_in.hierarchy.only_model().chains():
    for residue_group in chain.residue_groups():
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) > 1):
        residue_group.remove_atom_group(atom_groups[1])
        atom_groups[0].altloc = " "
        for atom in residue_group.atoms():
          atom.occ = 1.0
  f = open("shear_frag_single.pdb", "w")
  f.write(pdb_in.hierarchy.as_pdb_string(crystal_symmetry=xrs))
  f.close()

if (__name__ == "__main__"):
  exercise()
  print("OK")
