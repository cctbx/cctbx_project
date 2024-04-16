
from __future__ import absolute_import, division, print_function
import os
from six.moves import cStringIO as StringIO
from libtbx.utils import null_out
from libtbx import group_args, Auto
import iotbx.pdb
import mmtbx.ions.identify
from mmtbx.ions.svm.dump_sites import master_phil
from iotbx.data_manager import DataManager


def exercise():
  from mmtbx.regression import make_fake_anomalous_data
  import mmtbx.ions.utils
  pdb_in = """\
CRYST1   51.491   51.491   35.389  90.00  90.00 120.00 P 31 2 1
SCALE1      0.019421  0.011213  0.000000        0.00000
SCALE2      0.000000  0.022425  0.000000        0.00000
SCALE3      0.000000  0.000000  0.028257        0.00000
HETATM   32  CA  CGU A  17       7.453  25.360  36.702  1.00 25.21           C
HETATM   33  C   CGU A  17       6.252  24.666  36.060  1.00 24.08           C
HETATM   34  O   CGU A  17       6.408  23.698  35.327  1.00 22.85           O
HETATM   35  CB  CGU A  17       7.547  24.924  38.163  1.00 28.34           C
HETATM   36  CG  CGU A  17       8.807  24.090  38.525  1.00 29.46           C
HETATM   37  CD1 CGU A  17       9.396  23.286  37.336  1.00 28.04           C
HETATM   38  CD2 CGU A  17       8.411  23.255  39.740  1.00 32.29           C
HETATM   39 OE11 CGU A  17      10.339  23.775  36.690  1.00 31.46           O
HETATM   40 OE12 CGU A  17       8.917  22.160  37.075  1.00 26.97           O
HETATM   41 OE21 CGU A  17       7.958  23.926  40.668  1.00 35.00           O
HETATM   42 OE22 CGU A  17       8.527  22.036  39.780  1.00 33.69           O
ATOM     43  N   PRO A  18       5.029  25.135  36.349  1.00 23.16           N
ATOM     62  CA  ARG A  20       7.902  23.943  32.052  1.00 22.37           C
ATOM     63  C   ARG A  20       7.515  22.468  32.019  1.00 24.90           C
ATOM     64  O   ARG A  20       7.956  21.738  31.130  1.00 24.00           O
ATOM     65  CB  ARG A  20       9.024  24.136  33.067  1.00 26.75           C
ATOM     67  CD  ARG A  20      10.812  25.597  34.000  1.00 36.42           C
HETATM   72  N   CGU A  21       6.701  22.022  32.980  1.00 24.22           N
HETATM   73  CA  CGU A  21       6.293  20.612  33.012  1.00 23.24           C
HETATM   74  C   CGU A  21       5.432  20.293  31.805  1.00 23.70           C
HETATM   75  O   CGU A  21       5.561  19.221  31.216  1.00 20.30           O
HETATM   76  CB  CGU A  21       5.506  20.267  34.289  1.00 24.58           C
HETATM   77  CG  CGU A  21       6.392  20.445  35.528  1.00 26.52           C
HETATM   78  CD1 CGU A  21       7.353  19.249  35.754  1.00 27.96           C
HETATM   79  CD2 CGU A  21       5.507  20.718  36.738  1.00 29.78           C
HETATM   80 OE11 CGU A  21       8.366  19.406  36.482  1.00 27.23           O
HETATM   81 OE12 CGU A  21       7.056  18.159  35.217  1.00 25.25           O
HETATM   82 OE21 CGU A  21       4.695  21.625  36.586  1.00 36.91           O
HETATM   83 OE22 CGU A  21       5.664  20.139  37.797  1.00 32.02           O
ATOM     93  C   CYS A  23       7.212  20.248  27.692  1.00 25.63           C
ATOM     94  O   CYS A  23       7.306  19.599  26.656  1.00 22.02           O
HETATM   97  N   CGU A  24       7.761  19.852  28.842  1.00 26.69           N
HETATM   98  CA  CGU A  24       8.527  18.607  28.931  1.00 29.70           C
HETATM   99  C   CGU A  24       7.665  17.456  28.476  1.00 31.08           C
HETATM  100  O   CGU A  24       8.143  16.541  27.812  1.00 32.94           O
HETATM  101  CB  CGU A  24       8.981  18.304  30.367  1.00 26.05           C
HETATM  102  CG  CGU A  24       9.966  19.357  30.876  1.00 26.18           C
HETATM  103  CD1 CGU A  24      11.275  19.290  30.093  1.00 24.75           C
HETATM  104  CD2 CGU A  24      10.148  19.172  32.390  1.00 27.43           C
HETATM  105 OE11 CGU A  24      12.023  18.293  30.233  1.00 29.79           O
HETATM  106 OE12 CGU A  24      11.537  20.244  29.348  1.00 24.99           O
HETATM  107 OE21 CGU A  24       9.100  19.190  33.043  1.00 28.87           O
HETATM  108 OE22 CGU A  24      11.260  19.084  32.908  1.00 24.87           O
ATOM    143  O   CYS A  29      10.353  21.841  23.789  1.00 30.74           O
ATOM    146  N   ASP A  30       9.604  19.770  24.234  1.00 32.83           N
ATOM    147  CA  ASP A  30      10.776  19.402  25.014  1.00 34.15           C
ATOM    148  C   ASP A  30      12.026  19.580  24.177  1.00 36.29           C
ATOM    149  O   ASP A  30      12.937  20.322  24.544  1.00 34.50           O
ATOM    150  CB  ASP A  30      10.685  17.949  25.464  1.00 33.18           C
ATOM    151  CG  ASP A  30      11.714  17.607  26.523  1.00 32.22           C
ATOM    152  OD1 ASP A  30      12.621  18.428  26.752  1.00 32.53           O
ATOM    153  OD2 ASP A  30      11.608  16.524  27.125  1.00 31.78           O
ATOM    154  N   GLU A  31      12.056  18.885  23.045  1.00 39.34           N
ATOM    155  CA  GLU A  31      13.186  18.954  22.135  1.00 40.16           C
ATOM    172  CA  ALA A  33      13.225  23.877  24.346  1.00 39.26           C
ATOM    173  C   ALA A  33      14.746  23.914  24.481  1.00 38.24           C
ATOM    175  CB  ALA A  33      12.600  23.326  25.630  1.00 37.33           C
ATOM    176  N   ASP A  34      15.400  22.799  24.170  1.00 39.56           N
ATOM    177  CA  ASP A  34      16.857  22.723  24.258  1.00 40.96           C
ATOM    180  CB  ASP A  34      17.352  21.300  23.976  1.00 40.20           C
ATOM    181  CG  ASP A  34      17.006  20.327  25.083  1.00 38.93           C
ATOM    182  OD1 ASP A  34      16.981  20.742  26.262  1.00 41.79           O
ATOM    183  OD2 ASP A  34      16.777  19.140  24.778  1.00 37.45           O
TER
HETATM  316 CA    CA A  71      13.077  17.433  32.271  1.00 22.23          CA
HETATM  317 CA    CA A  72      13.835  18.867  28.887  1.00 30.50          CA
HETATM  318 CA    CA A  73      10.897  18.813  35.385  1.00 50.79          CA
HETATM  320  O   HOH A  75      13.387  22.461  33.530  1.00 24.93           O
HETATM  323  O   HOH A  78      10.578  15.304  29.567  1.00 23.15           O
HETATM  324  O   HOH A  79       5.020  20.563  40.636  1.00 44.02           O
HETATM  325  O   HOH A  80       2.823  22.144  38.546  1.00 36.74           O
HETATM  326  O   HOH A  81      10.434  22.631  29.604  1.00 25.89           O
HETATM  327  O   HOH A  82       6.522  15.691  36.473  1.00 27.82           O
HETATM  332  O   HOH A  87      11.624  15.358  31.822  1.00 24.92           O
HETATM  333  O   HOH A  88      13.763  16.798  28.667  1.00 29.47           O
HETATM  334  O   HOH A  89       6.350  16.973  32.340  1.00 37.83           O
HETATM  338  O   HOH A  93      10.474  21.054  34.739  1.00 25.48           O
HETATM  342  O   HOH A  97      16.203  18.688  27.720  1.00 28.10           O
HETATM  343  O   HOH A  98       8.186  14.327  30.477  1.00 49.44           O
HETATM  344  O   HOH A  99       8.625  16.477  33.868  1.00 48.13           O
HETATM  347  O   HOH A 102      15.462  16.714  24.789  1.00 42.90           O
HETATM  356  O   HOH A 111       4.757  17.423  38.879  1.00 34.26           O
HETATM  358  O   HOH A 113      10.313  14.495  25.452  1.00 40.66           O
HETATM  359  O   HOH A 114       1.979  18.616  37.760  1.00 34.25           O
HETATM  363  O   HOH A 118      13.926  20.627  27.271  1.00 29.62           O
HETATM  365  O   HOH A 120      16.240  23.471  27.700  1.00 48.79           O
HETATM  370  O   HOH A 125       2.747  18.823  35.170  1.00 50.30           O
HETATM  372  O   HOH A 127       5.228  23.553  41.559  1.00 42.02           O
HETATM  373  O   HOH A 128       5.298  21.833  43.473  1.00 41.96           O
HETATM  377  O   HOH A 132      13.181  22.613  29.210  1.00 35.43           O
TER
"""
  file_base = "tst_symmetry_axis"
  with open(file_base + ".pdb", "w") as f:
    f.write(pdb_in)
  pdb_inp = iotbx.pdb.input(file_base + ".pdb")
  hierarchy = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  hierarchy, n = mmtbx.ions.utils.anonymize_ions(hierarchy, log=null_out())
  assert (n == 3)
  with open(file_base + "_in.pdb", "w") as f:
    f.write(hierarchy.as_pdb_string(crystal_symmetry=xrs))
  mtz_file = make_fake_anomalous_data.generate_mtz_file(
    file_base=file_base,
    d_min=1.5,
    anomalous_scatterers=[group_args(
      selection="element CA",
      fp=0.0,
      fdp=0.4)])

  dm = DataManager()
  m = dm.get_model(file_base + "_in.pdb")
  dm.process_miller_array_file(mtz_file)
  fmo = dm.get_fmodel(scattering_table="n_gaussian")
  params = master_phil().extract()
  params.use_phaser=False
  #params.skip_twin_test=True
  params.elements='CA'
  params.input.wavelength = 0.9792
  out = StringIO()
  results = get_analyze_waters_result(m,fmo,params,out)
  assert "Valence sum:  1.916" in out.getvalue()
  assert out.getvalue().count("Probable cation: CA+2") >= 1
  os.remove(file_base + ".pdb")
  os.remove(file_base + "_in.pdb")
  os.remove(file_base + ".mtz")
  os.remove(file_base + "_fmodel.eff")


def get_analyze_waters_result(m,fmo,params,out,manager_class = None):
  m.process(make_restraints=True)
  grm = m.get_restraints_manager()
  manager = mmtbx.ions.identify.create_manager(
    pdb_hierarchy = m.get_hierarchy(),
    fmodel = fmo,
    geometry_restraints_manager = grm.geometry,
    wavelength = params.input.wavelength,
    params = params,
    nproc = params.nproc,
    log = out,
    manager_class = manager_class)
  candidates = Auto
  if (params.elements is not Auto) and (params.elements is not None):
    from cctbx.eltbx import chemical_elements
    lu = chemical_elements.proper_upper_list()
    elements = params.elements.replace(",", " ")
    candidates = elements.split()
    for elem in candidates :
      if (elem.upper() not in lu):
        raise Sorry("Unrecognized element '%s'" % elem)

  results = manager.analyze_waters(
    out = out,
    candidates = candidates)

  return results

if (__name__ == "__main__"):
  exercise()
  print("OK")
