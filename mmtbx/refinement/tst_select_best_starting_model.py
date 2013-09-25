
from __future__ import division
from libtbx.utils import null_out
from libtbx import easy_run

def exercise () :
  from mmtbx.refinement import select_best_starting_model
  from iotbx import file_reader
  from cctbx import uctbx
  from cctbx import sgtbx
  from scitbx.array_family import flex
  import random
  unit_cell = (24.937, 8.866, 25.477, 90.00, 107.08, 90.00)
  space_group = "P21"
  pdb_base = """\
CRYST1   24.937    8.866   25.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.040101  0.000000  0.012321        0.00000
SCALE2      0.000000  0.112790  0.000000        0.00000
SCALE3      0.000000  0.000000  0.041062        0.00000
ATOM      1  N   GLY A   1       8.992   0.474  -6.096  1.00 16.23           N
ATOM      2  CA  GLY A   1       9.033   0.047  -4.707  1.00 16.20           C
ATOM      3  C   GLY A   1       7.998  -1.029  -4.448  1.00 15.91           C
ATOM      4  O   GLY A   1       7.548  -1.689  -5.385  1.00 16.11           O
ATOM      5  N   ASN A   2       7.625  -1.218  -3.185  1.00 15.02           N
ATOM      6  CA  ASN A   2       6.523  -2.113  -2.848  1.00 13.92           C
ATOM      7  C   ASN A   2       5.220  -1.618  -3.428  1.00 12.24           C
ATOM      8  O   ASN A   2       4.955  -0.418  -3.432  1.00 11.42           O
ATOM      9  CB  ASN A   2       6.376  -2.261  -1.340  1.00 14.42           C
ATOM     10  CG  ASN A   2       7.620  -2.786  -0.697  1.00 13.92           C
ATOM     11  OD1 ASN A   2       8.042  -3.915  -0.978  1.00 14.39           O
ATOM     12  ND2 ASN A   2       8.232  -1.975   0.168  1.00 12.78           N
ATOM     13  N   ASN A   3       4.406  -2.553  -3.904  1.00 12.20           N
ATOM     14  CA  ASN A   3       3.164  -2.226  -4.594  1.00 11.81           C
ATOM     15  C   ASN A   3       1.925  -2.790  -3.910  1.00 10.59           C
ATOM     16  O   ASN A   3       1.838  -3.991  -3.653  1.00 10.32           O
ATOM     17  CB  ASN A   3       3.231  -2.727  -6.046  1.00 12.51           C
ATOM     18  CG  ASN A   3       1.973  -2.405  -6.848  1.00 12.59           C
ATOM     19  OD1 ASN A   3       1.662  -1.239  -7.106  1.00 13.64           O
ATOM     20  ND2 ASN A   3       1.260  -3.443  -7.268  1.00 12.39           N
ATOM     21  N   GLN A   4       0.973  -1.913  -3.608  1.00 10.34           N
ATOM     22  CA  GLN A   4      -0.366  -2.335  -3.208  1.00 10.00           C
ATOM     23  C   GLN A   4      -1.402  -1.637  -4.085  1.00 10.21           C
ATOM     24  O   GLN A   4      -1.514  -0.414  -4.070  1.00  8.99           O
ATOM     25  CB  GLN A   4      -0.656  -2.027  -1.736  1.00 10.00           C
ATOM     26  CG  GLN A   4      -1.927  -2.705  -1.229  1.00 10.50           C
ATOM     27  CD  GLN A   4      -2.482  -2.102   0.060  1.00 11.36           C
ATOM     28  OE1 GLN A   4      -2.744  -0.900   0.151  1.00 12.29           O
ATOM     29  NE2 GLN A   4      -2.684  -2.951   1.055  1.00 10.43           N
ATOM     30  N   GLN A   5      -2.154  -2.406  -4.857  1.00 10.48           N
ATOM     31  CA  GLN A   5      -3.247  -1.829  -5.630  1.00 11.24           C
ATOM     32  C   GLN A   5      -4.591  -2.382  -5.178  1.00 11.40           C
ATOM     33  O   GLN A   5      -4.789  -3.599  -5.092  1.00 11.94           O
ATOM     34  CB  GLN A   5      -3.024  -2.023  -7.129  1.00 11.14           C
ATOM     35  CG  GLN A   5      -1.852  -1.222  -7.653  1.00 10.65           C
ATOM     36  CD  GLN A   5      -1.338  -1.748  -8.965  1.00 10.73           C
ATOM     37  OE1 GLN A   5      -0.794  -2.845  -9.028  1.00 10.14           O
ATOM     38  NE2 GLN A   5      -1.511  -0.968 -10.027  1.00 11.31           N
ATOM     39  N   ASN A   6      -5.504  -1.471  -4.872  1.00 11.56           N
ATOM     40  CA  ASN A   6      -6.809  -1.838  -4.359  1.00 12.07           C
ATOM     41  C   ASN A   6      -7.856  -1.407  -5.353  1.00 13.18           C
ATOM     42  O   ASN A   6      -8.257  -0.251  -5.362  1.00 13.64           O
ATOM     43  CB  ASN A   6      -7.053  -1.149  -3.017  1.00 12.12           C
ATOM     44  CG  ASN A   6      -5.966  -1.446  -1.998  1.00 12.31           C
ATOM     45  OD1 ASN A   6      -5.833  -2.579  -1.517  1.00 13.43           O
ATOM     46  ND2 ASN A   6      -5.198  -0.423  -1.645  1.00 11.88           N
ATOM     47  N   TYR A   7      -8.298  -2.332  -6.193  1.00 14.34           N
ATOM     48  CA  TYR A   7      -9.162  -1.980  -7.317  1.00 15.00           C
ATOM     49  C   TYR A   7     -10.603  -1.792  -6.893  1.00 15.64           C
ATOM     50  O   TYR A   7     -11.013  -2.278  -5.838  1.00 15.68           O
ATOM     51  CB  TYR A   7      -9.064  -3.041  -8.412  1.00 15.31           C
ATOM     52  CG  TYR A   7      -7.657  -3.197  -8.931  1.00 15.06           C
ATOM     53  CD1 TYR A   7      -6.785  -4.118  -8.368  1.00 15.24           C
ATOM     54  CD2 TYR A   7      -7.193  -2.400  -9.960  1.00 14.96           C
ATOM     55  CE1 TYR A   7      -5.489  -4.253  -8.830  1.00 14.94           C
ATOM     56  CE2 TYR A   7      -5.905  -2.526 -10.429  1.00 15.13           C
ATOM     57  CZ  TYR A   7      -5.055  -3.451  -9.861  1.00 14.97           C
ATOM     58  OH  TYR A   7      -3.768  -3.572 -10.335  1.00 14.93           O
ATOM     59  OXT TYR A   7     -11.378  -1.149  -7.601  1.00 15.89           O
TER
"""
  pdb_base_water = """\
HETATM   64  O   HOH S   1     -10.466  -2.347  -3.168  1.00 17.57           O
HETATM   65  O   HOH S   2       6.469   1.081  -7.070  1.00 21.27           O
HETATM   66  O   HOH S   3     -11.809   0.108  -9.956  1.00 27.52           O
HETATM   67  O   HOH S   4       1.580  -3.455 -11.035  1.00 44.76           O
END
"""
  open("tst_start_model_base.pdb", "w").write(pdb_base+pdb_base_water)
  params = """
    high_resolution = 1.75
    add_sigmas = True
    pdb_file = tst_start_model_base.pdb
    output {
      label = F
      type = *real complex
      file_name = tst_start_model_base.mtz
    }
  """
  open("tst_start_model_fmodel.eff", "w").write(params)
  assert (easy_run.fully_buffered(
    "phenix.fmodel tst_start_model_fmodel.eff"
  ).raise_if_errors().return_code == 0)
  mtz_in = file_reader.any_file("tst_start_model_base.mtz")
  f_obs = mtz_in.file_server.miller_arrays[0]
  symm = f_obs.crystal_symmetry().customized_copy(
    space_group_info=sgtbx.space_group_info("P2"))
  f_obs = f_obs.customized_copy(crystal_symmetry=symm)
  random.seed(12345) # XXX makes results more predictable
  flags = f_obs.generate_r_free_flags(fraction=0.1)
  mtz_data = f_obs.as_mtz_dataset(
    column_root_label="F")
  mtz_data.add_miller_array(flags,
    column_root_label="FreeR_flag")
  mtz_data.mtz_object().write("tst_start_model.mtz")
  pdb_in = file_reader.any_file("tst_start_model_base.pdb")
  hierarchy_in = pdb_in.file_object.construct_hierarchy()
  xrs_in = pdb_in.file_object.xray_structure_simple()
  selection = hierarchy_in.atom_selection_cache().selection
  # Model 1: very few changes, but shifted by (1,0,0.5)
  symm2 = xrs_in.crystal_symmetry().customized_copy(
    unit_cell=uctbx.unit_cell((24.932, 8.841, 25.501, 90.00, 107.5, 90.00)))
  #u_iso = xrs_in.extract_u_iso_or_u_equiv()
  #xrs_out = xrs_in.deep_copy_scatterers().set_u_iso(
  #  selection=flex.bool(u_iso.size(), True), values=u_iso*1.1)
  xrs_out = xrs_in.deep_copy_scatterers()
  sites_cart = xrs_out.sites_cart()
  sites_cart += flex.vec3_double(sites_cart.size(), (1.0, 0.0, 0.5))
  xrs_out.set_sites_cart(sites_cart)
  hierarchy_out = hierarchy_in.deep_copy()
  hierarchy_out.adopt_xray_structure(xrs_out)
  open("tst_start_model_1.pdb", "w").write(
    hierarchy_out.as_pdb_string(crystal_symmetry=symm2))
  # Model 2: no sidechains
  mc_sele = selection("(name N or name C or name O or name CA or name CB)")
  hierarchy_out = hierarchy_in.select(mc_sele)
  xrs_out = xrs_in.select(mc_sele)
  open("tst_start_model_2.pdb", "w").write(
    hierarchy_out.as_pdb_string(crystal_symmetry=xrs_out))
  # Model 3: P1 symmetry
  symm3 = xrs_in.crystal_symmetry().customized_copy(
    space_group_info=sgtbx.space_group_info("P1"))
  open("tst_start_model_3.pdb", "w").write(
    hierarchy_out.as_pdb_string(crystal_symmetry=symm3))
  # Model 4: shaken coordinates and ADPs
  def random_double (size, factor=1) :
    d = flex.double()
    for x in range(size) :
      d.append(random.random() * factor)
    return d
  xrs_out = xrs_in.customized_copy()
  xrs_out.shake_sites_in_place(0.2, random_double=random_double)
  xrs_out.shake_adp()
  hierarchy_out = hierarchy_in.deep_copy()
  hierarchy_out.adopt_xray_structure(xrs_out)
  open("tst_start_model_4.pdb", "w").write(
    hierarchy_out.as_pdb_string(crystal_symmetry=xrs_out))
  # Model 5: perfect, but missing CRYST1
  open("tst_start_model_5.pdb", "w").write(hierarchy_in.as_pdb_string())
  # run method
  params = select_best_starting_model.master_phil.extract()
  params.rigid_body_refine = False
  result = select_best_starting_model.select_model(
    model_file_names=[
      "tst_start_model_1.pdb",
      "tst_start_model_2.pdb",
      "tst_start_model_3.pdb",
      "tst_start_model_4.pdb",
      "tst_start_model_5.pdb"],
    f_obs=f_obs,
    r_free_flags=flags,
    params=params,
    skip_twin_detection=True,
    log=null_out())
  #result.show(verbose=True)
  assert (result.best_model_file == "tst_start_model_4.pdb")
  params.rigid_body_refine = True
  result = select_best_starting_model.select_model(
    model_file_names=[
      "tst_start_model_1.pdb",
      "tst_start_model_2.pdb",
      "tst_start_model_3.pdb",
      "tst_start_model_4.pdb",
      "tst_start_model_5.pdb"],
    f_obs=f_obs,
    r_free_flags=flags,
    params=params,
    skip_twin_detection=True,
    log=null_out())
  #result.show(verbose=True)
  assert (result.best_model_file == "tst_start_model_1.pdb")

if (__name__ == "__main__") :
  exercise()
  print "OK"
