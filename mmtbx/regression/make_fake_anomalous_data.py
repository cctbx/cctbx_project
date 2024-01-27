
from __future__ import absolute_import, division, print_function
from mmtbx.ions.utils import anonymize_ions
from libtbx.utils import null_out
from libtbx import group_args
from libtbx import easy_run
import iotbx.pdb
import os

def generate_calcium_inputs(file_base="ca_frag", anonymize=True):
  """
  Generate both a PDB file and an MTZ file for the calcium-bound structure,
  with the calcium optionally  replaced by solvent after F(model) was
  calculated. Calcium is simulated as an anomalous scatterer at 1.025 A.

  (This method is also a suitable template for generating any other model/data
  files for specific structures.)

  Parameters
  ----------
  file_base : str, optional
  anonymize : bool, optional
      Replace all ions in the returned pdb file with waters.

  Returns
  -------
  mtz_path : str
  pdb_path : str
  """
  pdb_file = write_pdb_input_calcium_binding(file_base=file_base)
  mtz_file = generate_mtz_file(
    file_base=file_base,
    d_min=1.5,
    anomalous_scatterers=[
      group_args(selection="element CA", fp=0.25, fdp=0.5)
      ])
  assert (os.path.isfile(pdb_file)) and (os.path.isfile(mtz_file))
  if anonymize:
    pdb_in = iotbx.pdb.input(pdb_file)
    hierarchy = pdb_in.construct_hierarchy()
    hierarchy, n = anonymize_ions(hierarchy, log=null_out())
    pdb_file = file_base + "_hoh.pdb"
    hierarchy.write_pdb_file(
      file_name=pdb_file,
      crystal_symmetry=pdb_in.crystal_symmetry())
    assert os.path.isfile(pdb_file)
  return os.path.abspath(mtz_file), os.path.abspath(pdb_file)

def generate_cd_cl_inputs(file_base="cd_cl_frag"):
  """
  Creates a fake model and reflection data for a structure containing cadmium
  and chloride ions.

  Parameters
  ----------
  file_base : str, optional

  Returns
  -------
  mtz_path : str
  pdb_path : str
  """
  pdb_file = write_pdb_input_cd_cl(file_base=file_base)
  mtz_file = generate_mtz_file(
    file_base=file_base,
    d_min=1.3,
    anomalous_scatterers=[
      group_args(selection="element CD", fp=-0.29, fdp=2.676),
      group_args(selection="element CL", fp=0.256, fdp=0.5),
    ])
  assert (os.path.isfile(pdb_file)) and (os.path.isfile(mtz_file))
  return mtz_file, pdb_file

def generate_zinc_inputs(file_base="zn_frag", anonymize=True,
    wavelength=None):
  """
  Generate both a PDB file and an MTZ file for the zinc-bound structure,
  with the zinc optionally replaced by solvent after F(model) was
  calculated. Zinc is simulated as an anomalous scatterer at 1.54 A.

  Parameters
  ----------
  file_base : str, optional
  anonymize : bool, optional
      Replace all ions in the returned pdb file with waters.

  Returns
  -------
  mtz_path : str
  pdb_path : str
  """
  anom = None
  if (wavelength is None):
    fp = -1.3
    fdp = 0.47
    anom = [ group_args(selection="element ZN", fp=fp, fdp=fdp), ]
  pdb_file = write_pdb_input_zinc_binding(file_base=file_base)
  mtz_file = generate_mtz_file(
    file_base=file_base,
    d_min=1.9,
    anomalous_scatterers=anom,
    wavelength=wavelength)
  assert os.path.isfile(pdb_file) and os.path.isfile(mtz_file)
  if anonymize:
    pdb_in = iotbx.pdb.input(pdb_file)
    hierarchy = pdb_in.construct_hierarchy()
    hierarchy, n = anonymize_ions(hierarchy, log=null_out())
    pdb_file = file_base + "_hoh.pdb"
    hierarchy.write_pdb_file(
      file_name=pdb_file,
      crystal_symmetry=pdb_in.crystal_symmetry())
    assert os.path.isfile(pdb_file)
  return os.path.abspath(mtz_file), os.path.abspath(pdb_file)

def generate_magnessium_inputs(file_base="mg_frag", anonymize=True):
  """
  Creates a fake model and reflection data for a structure containing magnesium
  ions.

  Parameters
  ----------
  file_base : str, optional
  anonymize : bool, optional
      Replace all ions in the returned pdb file with waters.

  Returns
  -------
  mtz_path : str
  pdb_path : str
  """
  pdb_file = write_pdb_input_magnessium_binding (file_base=file_base)
  mtz_file = generate_mtz_file(
    file_base=file_base,
    d_min=1.5)
  assert os.path.isfile(pdb_file) and os.path.isfile(mtz_file)
  if anonymize:
    pdb_in = iotbx.pdb.input(pdb_file)
    hierarchy = pdb_in.construct_hierarchy()
    hierarchy, n = anonymize_ions(hierarchy, log=null_out())
    pdb_file = file_base + "_hoh.pdb"
    hierarchy.write_pdb_file(
      file_name=pdb_file,
      crystal_symmetry=pdb_in.crystal_symmetry())
    assert os.path.isfile(pdb_file)
  return os.path.abspath(mtz_file), os.path.abspath(pdb_file)

def generate_mtz_file(file_base, d_min, anomalous_scatterers=None,
    wavelength=None):
  """
  Create an MTZ file containing amplitudes (and R-free-flags) calculated from
  the PDB file, at the specified resolution and with enumerated anomalous
  scatterer groups.

  Parameters
  ----------
  file_base : str
  d_min : float
  anomalous_scatterers : group_args

  Returns
  -------
  mtz_path : str

  Examples
  --------
  >>> from libtbx import group_args
  >>> from mmtbx.regression import make_fake_anomalous_data as mfad
  >>> file_base = "tst_gen_file"
  >>> pdb_file = mfad.write_pdb_input_calcium_binding(file_base=file_base)
  >>> scatterers = group_args(selection="element CA", fp=0.25, fdp=0.5)
  >>> mtz_file = mfad.generate_mtz_file(
  ...   file_base, 1.5, anomalous_scatterers=[scatterers]
  ...   )
  >>> assert os.path.isfile(pdb_file) and os.path.isfile(mtz_file)
  """
  mtz_file = file_base + ".mtz"
  if (os.path.isfile(mtz_file)):
    os.remove(mtz_file)
  params = """
    high_resolution = %g
    r_free_flags_fraction = 0.1
    add_sigmas = True
    wavelength = %s
    output {
      label = F
      type = *real complex
      file_name = %s.mtz
    }\n""" % (d_min, wavelength, file_base)
  if anomalous_scatterers is not None:
    params += "anomalous_scatterers {\n"
    for group in anomalous_scatterers :
      params += ("  group {\n"
                 "    selection = %s\n"
                 "    f_prime = %g\n"
                 "    f_double_prime = %g\n"
                 "  }\n") % (group.selection, group.fp, group.fdp)
    params += "}\n"
  eff_file = file_base + "_fmodel.eff"
  with open(eff_file, "w") as f:
    f.write(params)
  assert (easy_run.fully_buffered("phenix.fmodel %s.pdb %s" % (file_base, eff_file),
    ).raise_if_errors().return_code == 0)
  return mtz_file

def write_pdb_input_calcium_binding(file_base="ca_frag", write_files=True):
  """
  Outputs a selection of atoms from a structure of a calcium-binding fold.
  Original data: d_min=1.4, wavelength=1.116

  Parameters
  ----------
  file_base : str, optional
  write_files : bool, optional
      If false, just return the pdb hierarchy and x-ray crystal structure,
      instead of writing it to a file.

  Returns
  -------
  str or tuple of iotbx.pdb.hierarchy.root, cctbx.xray.structure.structure
  """
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N   ASP A  37      10.710  14.456   9.568  1.00 15.78           N
ATOM      2  CA  ASP A  37       9.318  14.587   9.999  1.00 18.38           C
ATOM      3  C   ASP A  37       8.402  13.523   9.395  1.00 15.46           C
ATOM      4  O   ASP A  37       7.295  13.309   9.892  1.00 16.65           O
ATOM      5  CB  ASP A  37       8.771  16.006   9.784  1.00 24.16           C
ATOM      6  CG  ASP A  37       8.900  16.478   8.357  1.00 28.30           C
ATOM      7  OD1 ASP A  37       9.145  15.643   7.457  1.00 30.98           O
ATOM      8  OD2 ASP A  37       8.747  17.700   8.132  1.00 36.64           O-1
ATOM      9  N   LEU A  38       8.862  12.840   8.346  1.00 16.14           N
ATOM     10  CA  LEU A  38       8.092  11.733   7.778  1.00 16.45           C
ATOM     11  C   LEU A  38       7.990  10.538   8.735  1.00 15.91           C
ATOM     12  O   LEU A  38       7.160   9.650   8.546  1.00 16.71           O
ATOM     13  CB  LEU A  38       8.675  11.305   6.425  1.00 17.74           C
ATOM     14  CG  LEU A  38       8.715  12.386   5.345  1.00 18.16           C
ATOM     15  CD1 LEU A  38       9.091  11.771   4.003  1.00 20.50           C
ATOM     16  CD2 LEU A  38       7.391  13.166   5.250  1.00 19.36           C
ATOM     17  N   ALA A  39       8.828  10.526   9.770  1.00 13.52           N
ATOM     18  CA  ALA A  39       8.825   9.454  10.756  1.00 15.05           C
ATOM     19  C   ALA A  39       7.758   9.662  11.829  1.00 16.22           C
ATOM     20  O   ALA A  39       7.489   8.754  12.609  1.00 15.36           O
ATOM     21  CB  ALA A  39      10.201   9.341  11.405  1.00 13.76           C
ATOM     22  N   ILE A  40       7.164  10.855  11.883  1.00 14.32           N
ATOM     23  CA  ILE A  40       6.164  11.156  12.915  1.00 15.53           C
ATOM     24  C   ILE A  40       4.885  11.819  12.386  1.00 18.72           C
ATOM     25  O   ILE A  40       4.178  12.497  13.136  1.00 21.36           O
ATOM     26  CB  ILE A  40       6.752  12.039  14.047  1.00 17.76           C
ATOM     27  CG1 ILE A  40       7.253  13.375  13.493  1.00 19.23           C
ATOM     28  CG2 ILE A  40       7.871  11.314  14.773  1.00 19.67           C
ATOM     29  CD1 ILE A  40       7.538  14.403  14.583  1.00 22.28           C
ATOM     30  N   ASP A  41       4.569  11.592  11.115  1.00 17.33           N
ATOM     31  CA  ASP A  41       3.461  12.305  10.474  1.00 20.42           C
ATOM     32  C   ASP A  41       2.132  11.539  10.466  1.00 23.63           C
ATOM     33  O   ASP A  41       1.117  12.047   9.976  1.00 24.40           O
ATOM     34  CB  ASP A  41       3.852  12.742   9.051  1.00 19.02           C
ATOM     35  CG  ASP A  41       4.135  11.566   8.114  1.00 20.05           C
ATOM     36  OD1 ASP A  41       4.047  10.386   8.529  1.00 17.79           O
ATOM     37  OD2 ASP A  41       4.465  11.830   6.938  1.00 19.57           O-1
ATOM     38  N   GLY A  42       2.132  10.328  11.012  1.00 19.52           N
ATOM     39  CA  GLY A  42       0.917   9.531  11.073  1.00 19.91           C
ATOM     40  C   GLY A  42       0.622   8.804   9.775  1.00 25.32           C
ATOM     41  O   GLY A  42      -0.390   8.104   9.659  1.00 25.73           O
ATOM     42  N   ASN A  43       1.505   8.967   8.794  1.00 20.57           N
ATOM     43  CA  ASN A  43       1.326   8.337   7.488  1.00 20.84           C
ATOM     44  C   ASN A  43       2.269   7.149   7.310  1.00 20.54           C
ATOM     45  O   ASN A  43       3.484   7.322   7.215  1.00 18.86           O
ATOM     46  CB  ASN A  43       1.563   9.367   6.382  1.00 20.63           C
ATOM     47  CG  ASN A  43       1.087   8.896   5.020  1.00 19.17           C
ATOM     48  ND2 ASN A  43       1.027   9.825   4.069  1.00 24.27           N
ATOM     49  OD1 ASN A  43       0.793   7.718   4.813  1.00 21.70           O
ATOM     50  N   PRO A  44       1.711   5.933   7.258  1.00 20.33           N
ATOM     51  CA  PRO A  44       2.542   4.737   7.081  1.00 20.08           C
ATOM     52  C   PRO A  44       3.221   4.664   5.716  1.00 20.25           C
ATOM     53  O   PRO A  44       4.166   3.891   5.562  1.00 26.14           O
ATOM     54  CB  PRO A  44       1.538   3.589   7.222  1.00 20.20           C
ATOM     55  CG  PRO A  44       0.211   4.201   6.901  1.00 25.57           C
ATOM     56  CD  PRO A  44       0.288   5.596   7.442  1.00 19.95           C
ATOM     57  N   ALA A  45       2.750   5.441   4.743  1.00 20.05           N
ATOM     58  CA  ALA A  45       3.316   5.406   3.393  1.00 23.75           C
ATOM     59  C   ALA A  45       4.525   6.323   3.231  1.00 25.35           C
ATOM     60  O   ALA A  45       5.211   6.281   2.205  1.00 26.63           O
ATOM     61  CB  ALA A  45       2.246   5.741   2.353  1.00 24.89           C
ATOM     62  N   THR A  46       4.783   7.153   4.239  1.00 19.37           N
ATOM     63  CA  THR A  46       5.952   8.027   4.226  1.00 19.58           C
ATOM     64  C   THR A  46       6.923   7.556   5.292  1.00 19.22           C
ATOM     65  O   THR A  46       6.516   6.968   6.298  1.00 16.37           O
ATOM     66  CB  THR A  46       5.588   9.477   4.547  1.00 17.37           C
ATOM     67  CG2 THR A  46       4.599  10.035   3.514  1.00 21.06           C
ATOM     68  OG1 THR A  46       4.999   9.537   5.849  1.00 18.86           O
ATOM     69  N   SER A  47       8.206   7.816   5.068  1.00 17.65           N
ATOM     70  CA  SER A  47       9.224   7.383   6.009  1.00 16.72           C
ATOM     71  C   SER A  47      10.503   8.190   5.885  1.00 17.50           C
ATOM     72  O   SER A  47      10.765   8.834   4.859  1.00 15.22           O
ATOM     73  CB  SER A  47       9.541   5.916   5.774  1.00 17.26           C
ATOM     74  OG  SER A  47      10.208   5.738   4.530  1.00 18.61           O
ATOM     75  N   TRP A  48      11.299   8.139   6.947  1.00 14.32           N
ATOM     76  CA  TRP A  48      12.634   8.723   6.955  1.00 13.15           C
ATOM     77  C   TRP A  48      13.652   7.627   6.683  1.00 13.64           C
ATOM     78  O   TRP A  48      13.756   6.646   7.433  1.00 12.78           O
ATOM     79  CB  TRP A  48      12.912   9.395   8.300  1.00 13.01           C
ATOM     80  CG  TRP A  48      14.327   9.861   8.482  1.00 13.43           C
ATOM     81  CD1 TRP A  48      14.924  10.924   7.876  1.00 14.36           C
ATOM     82  CD2 TRP A  48      15.313   9.291   9.357  1.00 13.44           C
ATOM     83  CE2 TRP A  48      16.485  10.060   9.221  1.00 12.26           C
ATOM     84  CE3 TRP A  48      15.315   8.202  10.239  1.00 12.49           C
ATOM     85  NE1 TRP A  48      16.228  11.048   8.307  1.00 13.42           N
ATOM     86  CZ2 TRP A  48      17.653   9.778   9.940  1.00 13.96           C
ATOM     87  CZ3 TRP A  48      16.476   7.922  10.948  1.00 11.99           C
ATOM     88  CH2 TRP A  48      17.628   8.705  10.790  1.00 13.64           C
TER      89      TRP A  48
ATOM     90  N   ILE A 153      12.422   4.326   8.671  1.00 13.05           N
ATOM     91  CA  ILE A 153      11.440   4.432   9.754  1.00 14.48           C
ATOM     92  C   ILE A 153      10.194   5.191   9.295  1.00 12.58           C
ATOM     93  O   ILE A 153      10.264   6.388   8.994  1.00 13.74           O
ATOM     94  CB  ILE A 153      12.019   5.189  10.975  1.00 15.38           C
ATOM     95  CG1 ILE A 153      13.320   4.550  11.466  1.00 15.63           C
ATOM     96  CG2 ILE A 153      10.986   5.280  12.090  1.00 15.20           C
ATOM     97  CD1 ILE A 153      13.950   5.283  12.647  1.00 15.77           C
ATOM     98  N   SER A 154       9.051   4.511   9.278  1.00 14.97           N
ATOM     99  CA  SER A 154       7.820   5.135   8.801  1.00 14.50           C
ATOM    100  C   SER A 154       7.001   5.838   9.874  1.00 15.23           C
ATOM    101  O   SER A 154       6.353   6.843   9.589  1.00 15.50           O
ATOM    102  CB  SER A 154       6.944   4.134   8.052  1.00 16.67           C
ATOM    103  OG  SER A 154       6.522   3.100   8.920  1.00 16.82           O
ATOM    104  N   GLU A 155       7.018   5.314  11.097  1.00 15.34           N
ATOM    105  CA  GLU A 155       6.218   5.885  12.176  1.00 15.01           C
ATOM    106  C   GLU A 155       6.825   5.575  13.533  1.00 12.95           C
ATOM    107  O   GLU A 155       7.234   4.443  13.789  1.00 15.35           O
ATOM    108  CB  GLU A 155       4.780   5.346  12.139  1.00 17.13           C
ATOM    109  CG  GLU A 155       3.862   6.037  11.124  1.00 22.41           C
ATOM    110  CD  GLU A 155       3.866   7.552  11.272  1.00 24.24           C
ATOM    111  OE1 GLU A 155       3.618   8.050  12.389  1.00 26.14           O
ATOM    112  OE2 GLU A 155       4.132   8.250  10.271  1.00 23.46           O-1
ATOM    113  N   ILE A 156       6.882   6.596  14.382  1.00 13.75           N
ATOM    114  CA  ILE A 156       7.220   6.440  15.788  1.00 12.29           C
ATOM    115  C   ILE A 156       6.066   7.007  16.602  1.00 11.05           C
ATOM    116  O   ILE A 156       5.584   8.113  16.320  1.00 13.14           O
ATOM    117  CB  ILE A 156       8.506   7.196  16.147  1.00 11.71           C
ATOM    118  CG1 ILE A 156       9.714   6.573  15.433  1.00 12.75           C
ATOM    119  CG2 ILE A 156       8.742   7.158  17.657  1.00 12.85           C
ATOM    120  CD1 ILE A 156      10.932   7.490  15.398  1.00 13.97           C
TER     121      ILE A 156
HETATM  122 CA   CA  S   1       5.334   8.357   8.032  1.00  5.89          CA+2
HETATM  123  O   HOH S   2       5.396  15.243  10.734  1.00 22.95           O
HETATM  124  O   HOH S   3       3.629   8.994  14.414  1.00 25.28           O
END""")
  xrs = pdb_in.xray_structure_simple(cryst1_substitution_buffer_layer=5)
  if (write_files):
    pdb_file = file_base + ".pdb"
    if (os.path.exists(pdb_file)):
      os.remove(pdb_file)
    f = open(pdb_file, "w")
    f.write(pdb_in.construct_hierarchy().as_pdb_string(crystal_symmetry=xrs))
    f.close()
    return pdb_file
  else :
    return pdb_in, xrs

def write_pdb_input_cd_cl(file_base="cd_cl_frag"):
  """
  Outputs a selection of atoms from a structure crystallized in cadmium cloride.
  Original data: d_min=1.37, wavelength=1.116

  Parameters
  ----------
  file_base : str, optional

  Returns
  -------
  str
  """
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N   GLU A  18     -13.959  12.159  -6.598  1.00 10.08           N
ATOM      2  CA  GLU A  18     -13.297  13.465  -6.628  1.00  9.83           C
ATOM      3  C   GLU A  18     -11.946  13.282  -7.309  1.00  9.18           C
ATOM      4  CB  GLU A  18     -13.128  14.035  -5.210  1.00 11.96           C
ATOM      5  CG  GLU A  18     -14.455  14.401  -4.522  1.00 13.56           C
ATOM      6  CD  GLU A  18     -14.291  15.239  -3.242  1.00 14.89           C
ATOM      7  OE1 GLU A  18     -14.172  14.646  -2.143  1.00 14.24           O
ATOM      8  OE2 GLU A  18     -14.309  16.498  -3.306  1.00 14.37           O1-
ATOM      9  N   PRO A  19     -11.418  14.350  -7.928  1.00  9.57           N
ATOM     10  CD  PRO A  19     -12.043  15.673  -8.082  1.00 11.94           C
ATOM     11  CB  SER A  20      -7.426  14.613  -4.832  1.00 12.48           C
ATOM     12  CG  ASP A  21     -12.019   9.806  -4.306  1.00 11.79           C
ATOM     13  OD1 ASP A  21     -11.866   9.680  -3.077  1.00 12.53           O
ATOM     14  OD2 ASP A  21     -13.123  10.103  -4.806  1.00 14.00           O1-
TER      15      ASP A  21
ATOM     16  C   GLU C  16     -17.717  21.556   1.133  1.00  9.93           C
ATOM     17  O   GLU C  16     -17.663  20.929   0.078  1.00 13.55           O
ATOM     18  N   VAL C  17     -16.747  21.509   2.043  1.00  9.02           N
ATOM     19  CA  VAL C  17     -15.504  20.772   1.813  1.00 10.30           C
ATOM     20  C   VAL C  17     -15.061  20.062   3.080  1.00 10.37           C
ATOM     21  O   VAL C  17     -15.568  20.347   4.166  1.00 12.61           O
ATOM     22  CB  VAL C  17     -14.357  21.714   1.368  1.00 11.35           C
ATOM     23  CG1 VAL C  17     -14.653  22.284  -0.019  1.00 14.25           C
ATOM     24  CG2 VAL C  17     -14.155  22.832   2.385  1.00 11.46           C
ATOM     25  N   GLU C  18     -14.120  19.134   2.943  1.00 10.25           N
ATOM     26  CA  GLU C  18     -13.458  18.521   4.087  1.00  9.78           C
ATOM     27  C   GLU C  18     -12.108  19.210   4.277  1.00  9.04           C
ATOM     28  O   GLU C  18     -11.528  19.727   3.312  1.00 10.29           O
ATOM     29  CB  GLU C  18     -13.263  17.016   3.870  1.00 11.44           C
ATOM     30  CG  GLU C  18     -14.568  16.219   3.805  1.00 14.22           C
ATOM     31  CD  GLU C  18     -14.339  14.700   3.837  1.00 18.19           C
ATOM     32  OE1 GLU C  18     -14.090  14.093   2.769  1.00 17.19           O
ATOM     33  OE2 GLU C  18     -14.408  14.090   4.930  1.00 18.29           O1-
ATOM     34  N   PRO C  19     -11.589  19.217   5.512  1.00  8.64           N
ATOM     35  N   SER C  20      -9.117  18.227   4.541  1.00  8.92           N
ATOM     36  CA  SER C  20      -8.002  17.726   3.740  1.00  9.11           C
ATOM     37  C   SER C  20      -8.323  17.664   2.251  1.00 10.30           C
ATOM     38  O   SER C  20      -7.579  17.059   1.490  1.00 12.39           O
ATOM     39  N   ASP C  21      -9.424  18.283   1.826  1.00  9.10           N
ATOM     40  CA  ASP C  21      -9.687  18.394   0.390  1.00  8.68           C
ATOM     41  C   ASP C  21      -8.518  19.104  -0.299  1.00  8.06           C
ATOM     42  O   ASP C  21      -7.956  20.072   0.240  1.00  9.29           O
ATOM     43  CB  ASP C  21     -10.987  19.165   0.130  1.00  9.00           C
ATOM     44  CG  ASP C  21     -12.175  18.250  -0.151  1.00 12.21           C
ATOM     45  OD1 ASP C  21     -11.996  17.207  -0.825  1.00 13.40           O
ATOM     46  OD2 ASP C  21     -13.298  18.562   0.293  1.00 14.53           O1-
ATOM     47  N   THR C  22      -8.134  18.619  -1.473  1.00  7.93           N
ATOM     48  CA  THR C  22      -7.123  19.287  -2.268  1.00  8.33           C
ATOM     49  C   THR C  22      -7.719  20.485  -2.992  1.00  8.59           C
ATOM     50  O   THR C  22      -8.935  20.600  -3.185  1.00  7.86           O
ATOM     51  OG1 THR C  22      -7.612  17.883  -4.158  1.00  8.83           O
TER      52      THR C  22
ATOM     53  N   ASN C  25      -9.688  19.311  -5.824  1.00  7.73           N
ATOM     54  CA  ASN C  25     -11.065  18.970  -5.481  1.00  8.22           C
ATOM     55  C   ASN C  25     -11.961  20.210  -5.409  1.00  7.89           C
ATOM     56  O   ASN C  25     -13.107  20.192  -5.859  1.00  8.41           O
ATOM     57  CB  ASN C  25     -11.124  18.234  -4.145  1.00  8.10           C
ATOM     58  CG  ASN C  25     -10.667  16.786  -4.233  1.00  8.60           C
ATOM     59  ND2 ASN C  25     -10.637  16.128  -3.077  1.00 11.20           N
ATOM     60  OD1 ASN C  25     -10.352  16.264  -5.306  1.00  9.78           O
ATOM     61  N   VAL C  26     -11.433  21.291  -4.843  1.00  7.87           N
ATOM     62  CA  VAL C  26     -12.169  22.548  -4.780  1.00  8.48           C
ATOM     63  CG2 VAL C  26     -11.277  23.116  -2.518  1.00 13.29           C
TER      64      VAL C  26
ATOM     65  CB  LYS C  29     -16.653  22.060  -5.667  1.00  9.82           C
ATOM     66  CG  LYS C  29     -16.565  20.717  -4.972  1.00 12.48           C
ATOM     67  CD  LYS C  29     -16.305  20.909  -3.485  1.00 16.22           C
ATOM     68  CE  LYS C  29     -16.187  19.579  -2.753  1.00 22.76           C
ATOM     69  NZ  LYS C  29     -17.525  18.994  -2.465  1.00 26.67           N1+
TER      70      LYS C  29
ATOM     71  CD1 LEU C  56     -10.613  22.962   1.772  1.00 12.91           C
TER      72      LEU C  56
HETATM   73  O   HOH S   4     -10.355  10.447   0.977  1.00 18.33           O
HETATM   74  O   HOH S  12     -14.439  17.636  -5.803  1.00 14.31           O
HETATM   75  O   HOH S  19     -13.395  12.091  -1.881  1.00 18.37           O
HETATM   76  O   HOH S  28     -13.593  15.089   0.608  1.00 19.93           O
HETATM   77  O   HOH S  41     -10.275  12.953  -2.172  1.00 18.72           O
HETATM   78  O   HOH S 121     -18.485  11.552  -2.969  1.00 23.16           O
HETATM   79  O   HOH S 124     -16.053  16.459  -7.803  1.00 32.54           O
HETATM   80  O   HOH S 134     -10.338  14.782   1.447  1.00 18.96           O
HETATM   81  O   HOH S 162     -17.702  13.108  -4.810  1.00 41.20           O
HETATM   82  O   HOH S 204     -17.671  15.077  -3.159  1.00 42.64           O
HETATM   84 CD   CD  X   2     -14.095  16.409  -0.841  1.00 15.00      ION CD2+
HETATM   86 CL   CL  Y   3     -16.374  16.899  -0.285  1.00 10.22      ION CL1-
END""")
  xrs = pdb_in.xray_structure_simple(cryst1_substitution_buffer_layer=5)
  pdb_file = file_base + ".pdb"
  if (os.path.exists(pdb_file)):
    os.remove(pdb_file)
  f = open(pdb_file, "w")
  f.write("""REMARK 200  WAVELENGTH OR RANGE        (A) : 1.116\n""")
  f.write(pdb_in.construct_hierarchy().as_pdb_string(crystal_symmetry=xrs))
  f.close()
  return pdb_file

def write_pdb_input_zinc_binding(file_base="zn_frag", write_files=True):
  """
  Outputs a selection of atoms from a structure crystallized with zinc.

  Parameters
  ----------
  file_base : str, optional
  write_files : bool, optional
      If false, just return the pdb hierarchy and x-ray crystal structure,
      instead of writing it to a file.

  Returns
  -------
  str or tuple of iotbx.pdb.hierarchy.root, cctbx.xray.structure.structure
  """
  # XXX extracted from 2whz (thermolysin, wavelength = 1.54A), but with an
  # incomplete coordination shell
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  OD1 ASN A 112      34.902  38.381  -2.709  1.00 12.10           O
ATOM      2  C   ALA A 113      31.533  40.581  -3.237  1.00 12.06           C
ATOM      3  O   ALA A 113      32.560  40.184  -3.842  1.00 12.75           O
ATOM      4  N   PHE A 114      30.551  41.251  -3.854  1.00 11.88           N
ATOM      5  CA  PHE A 114      30.707  41.621  -5.277  1.00 12.52           C
ATOM      6  C   PHE A 114      29.633  42.613  -5.705  1.00 12.63           C
ATOM      7  CB  PHE A 114      30.665  40.371  -6.202  1.00 12.72           C
ATOM      8  N   TRP A 115      29.933  43.311  -6.794  1.00 13.30           N
ATOM      9  CA  TRP A 115      28.968  44.084  -7.575  1.00 13.94           C
ATOM     10  O   TRP A 115      29.375  42.549  -9.378  1.00 15.06           O
ATOM     11  CB  TRP A 115      29.646  45.329  -8.137  1.00 14.57           C
ATOM     13  O   VAL A 139      33.773  47.006  -1.063  1.00 11.73           O
ATOM     15  N   HIS A 142      34.994  49.588  -2.455  1.00 13.21           N
ATOM     16  CA  HIS A 142      35.593  48.855  -3.582  1.00 11.38           C
ATOM     17  C   HIS A 142      34.455  48.358  -4.484  1.00 12.11           C
ATOM     18  O   HIS A 142      34.468  48.596  -5.693  1.00 11.57           O
ATOM     19  CB  HIS A 142      36.442  47.680  -3.082  1.00 11.93           C
ATOM     20  CG  HIS A 142      36.945  46.756  -4.150  1.00  9.54           C
ATOM     21  CD2 HIS A 142      36.303  45.839  -4.912  1.00  9.50           C
ATOM     22  ND1 HIS A 142      38.289  46.628  -4.440  1.00 11.23           N
ATOM     23  CE1 HIS A 142      38.446  45.706  -5.381  1.00 14.08           C
ATOM     24  NE2 HIS A 142      37.258  45.206  -5.676  1.00 14.47           N
ATOM     25  N   GLU A 143      33.448  47.732  -3.879  1.00 11.17           N
ATOM     26  CA  GLU A 143      32.363  47.129  -4.683  1.00 12.04           C
ATOM     27  C   GLU A 143      31.546  48.215  -5.395  1.00 12.79           C
ATOM     28  O   GLU A 143      31.237  48.078  -6.584  1.00 12.97           O
ATOM     29  CB  GLU A 143      31.461  46.221  -3.834  1.00 11.99           C
ATOM     30  CG  GLU A 143      32.199  45.096  -3.030  1.00 12.94           C
ATOM     31  CD  GLU A 143      32.974  44.075  -3.867  1.00 15.46           C
ATOM     32  OE1 GLU A 143      33.029  44.153  -5.112  1.00 16.48           O
ATOM     33  OE2 GLU A 143      33.563  43.175  -3.251  1.00 16.80           O1-
ATOM     35  C   THR A 145      33.305  50.876  -8.608  1.00 12.68           C
ATOM     36  N   HIS A 146      33.071  49.560  -8.540  1.00 12.36           N
ATOM     37  CA  HIS A 146      32.869  48.794  -9.756  1.00 13.33           C
ATOM     38  CB  HIS A 146      32.643  47.318  -9.470  1.00 12.80           C
ATOM     39  CG  HIS A 146      33.919  46.567  -9.214  1.00 14.00           C
ATOM     40  CD2 HIS A 146      34.290  45.747  -8.206  1.00 12.81           C
ATOM     41  ND1 HIS A 146      34.992  46.623 -10.084  1.00 12.00           N
ATOM     42  CE1 HIS A 146      35.970  45.863  -9.617  1.00 15.40           C
ATOM     43  NE2 HIS A 146      35.552  45.289  -8.506  1.00 11.49           N
ATOM     45  CG  TYR A 157      36.255  42.005 -14.146  1.00 19.06           C
ATOM     46  CD1 TYR A 157      36.607  43.119 -13.385  1.00 20.01           C
ATOM     47  CD2 TYR A 157      36.279  40.754 -13.542  1.00 19.69           C
ATOM     48  CE1 TYR A 157      36.970  42.993 -12.047  1.00 18.89           C
ATOM     49  CE2 TYR A 157      36.645  40.612 -12.193  1.00 20.96           C
ATOM     50  CZ  TYR A 157      36.971  41.732 -11.462  1.00 20.87           C
ATOM     51  OH  TYR A 157      37.304  41.604 -10.142  1.00 21.10           O
ATOM     53  O   GLY A 162      40.788  44.717 -13.313  1.00 13.40           O
ATOM     55  CA  ASN A 165      39.280  47.961 -13.045  1.00 11.01           C
ATOM     56  C   ASN A 165      39.422  47.647 -11.555  1.00 11.21           C
ATOM     57  O   ASN A 165      38.940  48.402 -10.709  1.00 10.60           O
ATOM     58  CB  ASN A 165      38.414  46.889 -13.730  1.00 12.11           C
ATOM     59  CG  ASN A 165      36.977  46.958 -13.330  1.00 13.51           C
ATOM     60  ND2 ASN A 165      36.082  46.882 -14.307  1.00 12.19           N
ATOM     61  OD1 ASN A 165      36.661  47.056 -12.155  1.00 13.28           O
ATOM     62  N   GLU A 166      40.119  46.539 -11.250  1.00 11.43           N
ATOM     63  CA  GLU A 166      40.413  46.155  -9.851  1.00 10.59           C
ATOM     64  C   GLU A 166      41.166  47.218  -9.055  1.00 10.24           C
ATOM     65  O   GLU A 166      40.744  47.581  -7.940  1.00  9.90           O
ATOM     66  CB  GLU A 166      41.140  44.803  -9.783  1.00 10.41           C
ATOM     67  CG  GLU A 166      40.169  43.654  -9.741  1.00 11.01           C
ATOM     68  CD  GLU A 166      39.267  43.743  -8.537  1.00 15.78           C
ATOM     69  OE1 GLU A 166      39.745  43.921  -7.380  1.00 14.70           O
ATOM     70  OE2 GLU A 166      38.063  43.637  -8.765  1.00 12.05           O1-
ATOM     71  N   ALA A 167      42.243  47.759  -9.628  1.00 11.39           N
ATOM     72  C   ILE A 168      39.558  51.420  -8.190  1.00 10.98           C
ATOM     73  N   SER A 169      39.076  50.192  -8.040  1.00 10.62           N
ATOM     74  CA  SER A 169      38.198  49.860  -6.904  1.00 11.23           C
ATOM     75  C   SER A 169      38.998  49.850  -5.589  1.00 11.38           C
ATOM     76  O   SER A 169      38.484  50.325  -4.547  1.00 11.01           O
ATOM     77  CB  SER A 169      37.462  48.529  -7.144  1.00 12.04           C
ATOM     78  OG  SER A 169      36.259  48.790  -7.899  1.00 11.91           O
ATOM     79  N   ASP A 170      40.233  49.353  -5.645  1.00 11.74           N
ATOM     80  CA  ASP A 170      41.142  49.411  -4.436  1.00 11.68           C
ATOM     81  CB  ASP A 170      42.424  48.581  -4.612  1.00 11.14           C
ATOM     82  CG  ASP A 170      42.166  47.082  -4.475  1.00 10.59           C
ATOM     83  OD1 ASP A 170      41.076  46.659  -3.990  1.00 11.83           O
ATOM     84  OD2 ASP A 170      43.010  46.278  -4.909  1.00 12.10           O1-
ATOM     86  NE  ARG A 203      42.877  43.716  -3.830  1.00 12.26           N
ATOM     87  CZ  ARG A 203      41.692  43.111  -3.748  1.00 11.93           C
ATOM     88  NH1 ARG A 203      41.631  41.782  -3.754  1.00 12.22           N1+
ATOM     89  NH2 ARG A 203      40.579  43.829  -3.675  1.00 11.99           N
ATOM     91  C   VAL A 230      43.021  41.431 -10.925  1.00 14.80           C
ATOM     92  O   VAL A 230      43.311  42.599 -10.716  1.00 14.05           O
ATOM     93  CB  VAL A 230      41.294  41.138 -12.775  1.00 14.48           C
ATOM     94  CG2 VAL A 230      40.404  40.176 -11.957  1.00 14.82           C
ATOM     95  N   HIS A 231      42.923  40.503  -9.959  1.00 15.61           N
ATOM     96  CA  HIS A 231      43.233  40.790  -8.542  1.00 15.70           C
ATOM     97  CB  HIS A 231      42.434  39.898  -7.566  1.00 15.19           C
ATOM     98  CG  HIS A 231      40.955  40.000  -7.734  1.00 16.56           C
ATOM     99  CD2 HIS A 231      40.056  40.844  -7.174  1.00 17.34           C
ATOM    100  ND1 HIS A 231      40.241  39.168  -8.573  1.00 16.47           N
ATOM    101  CE1 HIS A 231      38.960  39.505  -8.529  1.00 15.36           C
ATOM    102  NE2 HIS A 231      38.825  40.518  -7.687  1.00 15.29           N
HETATM  104  N   ILE A1001      35.271  41.187  -4.048  1.00 16.78           N
HETATM  105  CA  ILE A1001      36.473  41.816  -3.412  1.00 17.77           C
HETATM  106  C   ILE A1001      37.764  41.046  -3.756  1.00 18.36           C
HETATM  107  O   ILE A1001      38.832  41.640  -3.838  1.00 16.99           O
HETATM  108  CB  ILE A1001      36.323  42.040  -1.875  1.00 17.98           C
HETATM  109  CG1 ILE A1001      37.464  42.899  -1.320  1.00 16.76           C
HETATM  110  CG2 ILE A1001      36.360  40.753  -1.098  1.00 20.77           C
HETATM  111  CD1 ILE A1001      37.576  44.234  -1.913  1.00 17.79           C
HETATM  112  N   TYR A1002      37.653  39.739  -3.981  1.00 18.85           N
HETATM  113  CA  TYR A1002      38.834  38.984  -4.422  1.00 20.77           C
HETATM  114  C   TYR A1002      38.430  37.858  -5.360  1.00 21.31           C
HETATM  115  O   TYR A1002      39.302  37.219  -5.983  1.00 22.07           O
HETATM  116  CB  TYR A1002      39.632  38.473  -3.225  1.00 21.20           C
HETATM  117  OXT TYR A1002      37.226  37.594  -5.543  1.00 22.10           O
HETATM  118  O   HOH A2080      32.747  43.340 -10.935  1.00 44.97           O
HETATM  119  O   HOH A2097      33.150  44.796 -12.822  1.00 35.29           O
HETATM  120  O   HOH A2190      32.628  42.780  -7.708  1.00 17.31           O
HETATM  121  O   HOH A2206      36.319  49.625 -10.585  1.00 16.17           O
HETATM  122  O   HOH A2210      33.342  47.349 -12.959  1.00 24.18           O
HETATM  123  O   HOH A2284      42.614  44.224  -6.731  1.00 15.85           O
HETATM  125 ZN   ZN  A1317      36.762  44.115  -7.277  1.00 15.32          ZN+2
END
""")
# XXX this is the missing water
#HETATM  124  O   HOH A2351      36.007  42.038  -6.782  1.00 10.54           O
  xrs = pdb_in.xray_structure_simple(cryst1_substitution_buffer_layer=5)
  if write_files:
    pdb_file = file_base + ".pdb"
    if (os.path.exists(pdb_file)):
      os.remove(pdb_file)
    with open(pdb_file, "w") as f:
      f.write(pdb_in.construct_hierarchy().as_pdb_string(crystal_symmetry=xrs))
    return pdb_file
  else:
    return pdb_in, xrs

def write_pdb_input_magnessium_binding(file_base="mg_frag", write_files=True):
  """
  Outputs a selection of atoms from a structure crystallized with magnesium.

  Parameters
  ----------
  file_base : str, optional
  write_files : bool, optional
      If false, just return the pdb hierarchy and x-ray crystal structure,
      instead of writing it to a file.

  Returns
  -------
  str or tuple of iotbx.pdb.hierarchy.root, cctbx.xray.structure.structure
  """
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
REMARK  fragment derived from 3ip0
ATOM      1  N   ASP A  95     -12.581   4.250  10.804  1.00  5.08           N
ATOM      2  CA  ASP A  95     -11.950   3.273   9.927  1.00  5.04           C
ATOM      3  C   ASP A  95     -11.505   2.089  10.784  1.00  4.97           C
ATOM      4  O   ASP A  95     -10.680   2.251  11.678  1.00  5.43           O
ATOM      5  CB  ASP A  95     -10.737   3.916   9.223  1.00  5.26           C
ATOM      6  CG  ASP A  95     -10.103   3.032   8.160  1.00  5.03           C
ATOM      7  OD1 ASP A  95     -10.879   2.437   7.372  1.00  5.23           O
ATOM      8  OD2 ASP A  95      -8.836   2.949   8.149  1.00  5.44           O1-
TER       9      ASP A  95
ATOM     10  N   ASP A  97      -9.724  -1.416  10.219  1.00  4.47           N
ATOM     11  CA  ASP A  97      -8.951  -2.179   9.250  1.00  4.64           C
ATOM     12  C   ASP A  97      -8.480  -3.469   9.910  1.00  4.67           C
ATOM     13  O   ASP A  97      -7.947  -3.449  11.042  1.00  5.63           O
ATOM     14  CB  ASP A  97      -7.746  -1.380   8.744  1.00  5.20           C
ATOM     15  CG  ASP A  97      -8.132  -0.223   7.828  1.00  5.06           C
ATOM     16  OD1 ASP A  97      -9.317  -0.118   7.445  1.00  5.31           O
ATOM     17  OD2 ASP A  97      -7.197   0.558   7.508  1.00  5.74           O1-
TER      18      ASP A  97
HETATM   19  PA  APC A 171     -11.494   1.693   3.101  1.00  5.93           P
HETATM   20  PB  APC A 171      -8.904   2.838   3.989  1.00  5.40           P
HETATM   21  PG  APC A 171      -6.189   1.883   3.824  1.00  5.44           P
HETATM   22  C5' APC A 171     -11.917  -0.062   1.209  1.00  6.23           C
HETATM   23  O5' APC A 171     -11.907   0.228   2.621  1.00  5.98           O
HETATM   24  C4' APC A 171     -12.005  -1.551   1.015  1.00  6.31           C
HETATM   25  O4' APC A 171     -10.811  -2.160   1.564  1.00  6.21           O
HETATM   26  C3' APC A 171     -13.142  -2.283   1.727  1.00  6.94           C
HETATM   27  O3' APC A 171     -14.374  -2.203   1.097  1.00  8.42           O
HETATM   28  C2' APC A 171     -12.584  -3.692   1.886  1.00  6.56           C
HETATM   29  O2' APC A 171     -13.250  -4.576   0.968  1.00  8.16           O
HETATM   30  C1' APC A 171     -11.120  -3.552   1.530  1.00  6.15           C
HETATM   31  N1  APC A 171      -9.299  -4.552   6.318  1.00  5.61           N
HETATM   32  O1A APC A 171     -12.241   2.714   2.317  1.00  6.59           O
HETATM   33  O1B APC A 171      -9.519   4.202   3.986  1.00  6.49           O
HETATM   34  O1G APC A 171      -5.931   1.987   5.301  1.00  5.83           O
HETATM   35  C2  APC A 171     -10.353  -3.827   5.903  1.00  5.64           C
HETATM   36  O2A APC A 171     -11.698   1.666   4.588  1.00  5.88           O
HETATM   37  O2B APC A 171      -8.855   2.079   5.313  1.00  5.21           O
HETATM   38  O2G APC A 171      -5.035   2.418   2.993  1.00  6.01           O
HETATM   39  C3A APC A 171      -9.734   1.834   2.758  1.00  5.70           C
HETATM   40  N3  APC A 171     -10.798  -3.610   4.683  1.00  5.70           N
HETATM   41  O3B APC A 171      -7.396   2.924   3.462  1.00  5.24           O
HETATM   42  O3G APC A 171      -6.675   0.503   3.400  1.00  6.98           O
HETATM   43  C4  APC A 171     -10.034  -4.246   3.763  1.00  5.51           C
HETATM   44  C5  APC A 171      -8.964  -5.056   4.033  1.00  5.40           C
HETATM   45  C6  APC A 171      -8.569  -5.207   5.373  1.00  5.37           C
HETATM   46  N6  APC A 171      -7.500  -5.912   5.734  1.00  5.71           N
HETATM   47  N7  APC A 171      -8.414  -5.587   2.862  1.00  5.73           N
HETATM   48  C8  APC A 171      -9.167  -5.072   1.920  1.00  6.02           C
HETATM   49  N9  APC A 171     -10.168  -4.254   2.383  1.00  5.71           N
HETATM   50  N1 AHHR A 181      -3.249   8.723   8.653  0.45  3.59           N
HETATM   51  C2 AHHR A 181      -3.170   8.473   9.936  0.45  3.27           C
HETATM   52  N2 AHHR A 181      -2.312   9.188  10.698  0.45  4.00           N
HETATM   53  N3 AHHR A 181      -3.822   7.446  10.541  0.45  3.91           N
HETATM   54  C4 AHHR A 181      -4.648   6.566   9.912  0.45  3.74           C
HETATM   55  O4 AHHR A 181      -5.221   5.635  10.511  0.45  4.44           O
HETATM   56  N5 AHHR A 181      -5.483   5.915   7.737  0.45  3.23           N
HETATM   57  C6 AHHR A 181      -5.582   6.218   6.448  0.45  3.92           C
HETATM   58  C6AAHHR A 181      -6.336   5.270   5.533  0.45  4.40           C
HETATM   59  O6AAHHR A 181      -7.145   4.378   6.216  0.45  4.54           O
HETATM   60  C7 AHHR A 181      -4.896   7.322   5.892  0.45  5.79           C
HETATM   61  N8 AHHR A 181      -4.134   8.150   6.611  0.45  5.53           N
HETATM   62  C9 AHHR A 181      -4.046   7.908   7.922  0.45  3.58           C
HETATM   63  C10AHHR A 181      -4.762   6.791   8.479  0.45  3.12           C
HETATM   64  N1 BHHS A 182      -3.604   8.719   8.748  0.55  6.49           N
HETATM   65  C2 BHHS A 182      -3.478   8.455  10.036  0.55  6.53           C
HETATM   66  N2 BHHS A 182      -2.673   9.243  10.809  0.55  6.01           N
HETATM   67  N3 BHHS A 182      -4.109   7.429  10.645  0.55  5.05           N
HETATM   68  C4 BHHS A 182      -4.896   6.564   9.963  0.55  5.93           C
HETATM   69  O4 BHHS A 182      -5.536   5.719  10.560  0.55  5.29           O
HETATM   70  N5 BHHS A 182      -5.756   5.984   7.764  0.55  7.30           N
HETATM   71  C6 BHHS A 182      -6.698   5.438   5.588  0.55  6.59           C
HETATM   72  C6ABHHS A 182      -5.887   6.282   6.483  0.55  5.89           C
HETATM   73  O6ABHHS A 182      -7.411   4.466   6.183  0.55  6.94           O
HETATM   74  O6BBHHS A 182      -6.769   5.796   4.398  0.55  7.23           O
HETATM   75  C7 BHHS A 182      -5.274   7.421   5.955  0.55  5.00           C
HETATM   76  N8 BHHS A 182      -4.530   8.228   6.695  0.55  4.94           N
HETATM   77  C9 BHHS A 182      -4.383   7.927   7.998  0.55  5.58           C
HETATM   78  C10BHHS A 182      -5.029   6.797   8.531  0.55  6.96           C
HETATM   79  O   HOH A 201     -12.056  -0.140   6.826  1.00  5.82           O
HETATM   80  O   HOH A 202     -10.078  -0.655   4.762  1.00  6.16           O
HETATM   81  O   HOH A 203      -5.851   3.057   8.110  1.00  6.43           O
HETATM   82  O   HOH A 205      -7.904   1.644  10.650  1.00  6.56           O
HETATM   83  O   HOH A 216     -13.376  -2.280   5.233  1.00  7.21           O
HETATM   84  O   HOH A 219      -5.353   2.683  10.917  1.00  8.33           O
HETATM   85  O   HOH A 243      -8.078   5.439   9.608  1.00  9.50           O
HETATM   86  O   HOH A 264     -14.344   0.622   4.874  1.00 14.43           O
HETATM   87  O  AHOH A 499      -7.735  -1.646   4.892  0.50  5.41           O
HETATM   88  O  BHOH A 499      -7.601  -1.751   4.259  0.50  8.56           O
HETATM   89  O  BHOH A 501     -13.472   3.920   6.947  0.56  7.54           O
HETATM   90  O  BHOH A 504      -9.634   5.594   6.224  0.25  5.98           O
HETATM   91  O  CHOH A 504      -9.119   6.099   7.366  0.25  7.90           O
HETATM   92 MG   MG  A 161     -10.440   0.924   6.059  1.00  5.29      MG  MG
HETATM   93 MG   MG  A 162      -7.326   2.422   6.749  1.00  5.40      MG  MG
""")
  xrs = pdb_in.xray_structure_simple(cryst1_substitution_buffer_layer=5)
  if write_files:
    pdb_file = file_base + ".pdb"
    if (os.path.exists(pdb_file)):
      os.remove(pdb_file)
    with open(pdb_file, "w") as f:
      f.write(pdb_in.construct_hierarchy().as_pdb_string(crystal_symmetry=xrs))
    return pdb_file
  else:
    return pdb_in, xrs
