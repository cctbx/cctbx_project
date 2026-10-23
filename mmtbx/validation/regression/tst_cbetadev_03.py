from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry
import iotbx.pdb
from mmtbx.validation import cbetadev
from mmtbx.validation import cbd_utils
from mmtbx.validation import cbd_mahalanobis

# cyclosporin A, chain D of 1m63: a D-amino acid (DAL 1), N-methyl residues
# (MLE, MVA, BMT), ABA, SAR and VAL 9 between two modified residues
pdb_str = """\
HETATM 5544  N   DAL D   1      86.641  15.755  29.240  1.00 30.97           N
HETATM 5545  CA  DAL D   1      87.439  16.415  28.216  1.00 31.47           C
HETATM 5546  C   DAL D   1      87.053  17.885  28.034  1.00 34.80           C
HETATM 5547  O   DAL D   1      87.007  18.643  29.015  1.00 37.80           O
HETATM 5548  CB  DAL D   1      88.907  16.300  28.567  1.00 29.17           C
HETATM 5549  N   MLE D   2      86.751  18.338  26.774  1.00 35.05           N
HETATM 5550  CA  MLE D   2      86.377  19.774  26.558  1.00 28.06           C
HETATM 5551  C   MLE D   2      84.983  19.904  25.891  1.00 26.97           C
HETATM 5552  O   MLE D   2      84.873  20.261  24.704  1.00 24.31           O
HETATM 5553  CB  MLE D   2      87.494  20.340  25.662  1.00 20.50           C
HETATM 5554  CG  MLE D   2      88.738  21.050  26.209  1.00 19.47           C
HETATM 5555  CD1 MLE D   2      89.424  20.285  27.285  1.00 24.18           C
HETATM 5556  CD2 MLE D   2      89.705  21.262  25.098  1.00 17.37           C
HETATM 5557  CN  MLE D   2      86.712  17.365  25.610  1.00 33.56           C
HETATM 5558  N   MLE D   3      83.886  19.688  26.655  1.00 24.64           N
HETATM 5559  CA  MLE D   3      82.566  19.713  26.039  1.00 20.55           C
HETATM 5560  C   MLE D   3      81.669  20.857  26.541  1.00 21.04           C
HETATM 5561  O   MLE D   3      81.949  21.479  27.576  1.00 19.41           O
HETATM 5562  CB  MLE D   3      81.963  18.337  26.365  1.00 19.71           C
HETATM 5563  CG  MLE D   3      82.344  17.065  25.605  1.00 20.65           C
HETATM 5564  CD1 MLE D   3      81.134  16.154  25.692  1.00 18.50           C
HETATM 5565  CD2 MLE D   3      82.722  17.321  24.145  1.00 17.87           C
HETATM 5566  CN  MLE D   3      84.018  19.533  28.161  1.00 28.13           C
HETATM 5567  N   MVA D   4      80.601  21.225  25.781  1.00 19.02           N
HETATM 5568  CA  MVA D   4      79.722  22.269  26.244  1.00 12.98           C
HETATM 5569  C   MVA D   4      78.297  21.731  26.401  1.00 12.68           C
HETATM 5570  O   MVA D   4      77.535  21.733  25.432  1.00 14.24           O
HETATM 5571  CB  MVA D   4      79.741  23.495  25.240  1.00 12.86           C
HETATM 5572  CG1 MVA D   4      78.895  24.630  25.773  1.00 17.04           C
HETATM 5573  CG2 MVA D   4      81.149  24.022  24.999  1.00  6.38           C
HETATM 5574  CN  MVA D   4      80.364  20.594  24.420  1.00 20.69           C
HETATM 5575  N   BMT D   5      77.921  21.164  27.565  1.00 11.78           N
HETATM 5576  CA  BMT D   5      76.551  20.663  27.696  1.00 12.89           C
HETATM 5577  C   BMT D   5      75.742  21.427  28.709  1.00 12.54           C
HETATM 5578  O   BMT D   5      76.271  22.198  29.496  1.00 13.68           O
HETATM 5579  CB  BMT D   5      76.474  19.103  27.945  1.00 16.62           C
HETATM 5580  OG1 BMT D   5      76.958  18.761  29.246  1.00 24.36           O
HETATM 5581  CG2 BMT D   5      77.180  18.235  26.828  1.00 13.25           C
HETATM 5582  CD1 BMT D   5      76.754  18.608  25.408  1.00 17.00           C
HETATM 5583  CD2 BMT D   5      76.991  16.744  27.094  1.00 11.35           C
HETATM 5584  CE  BMT D   5      75.525  16.382  27.250  1.00 17.06           C
HETATM 5585  CZ  BMT D   5      75.163  15.184  27.909  1.00 25.87           C
HETATM 5586  CH  BMT D   5      73.708  14.909  28.241  1.00 26.92           C
HETATM 5587  CN  BMT D   5      78.858  21.149  28.748  1.00  8.19           C
HETATM 5588  N   ABA D   6      74.442  21.169  28.694  1.00 16.92           N
HETATM 5589  CA  ABA D   6      73.496  21.855  29.552  1.00 20.29           C
HETATM 5590  C   ABA D   6      72.794  20.942  30.561  1.00 21.80           C
HETATM 5591  O   ABA D   6      71.615  20.603  30.381  1.00 24.35           O
HETATM 5592  CB  ABA D   6      72.442  22.561  28.676  1.00 26.62           C
HETATM 5593  CG  ABA D   6      72.986  23.451  27.568  1.00 22.26           C
HETATM 5594  N   SAR D   7      73.464  20.544  31.666  1.00 23.09           N
HETATM 5595  CA  SAR D   7      72.854  19.658  32.661  1.00 25.98           C
HETATM 5596  C   SAR D   7      73.804  18.553  33.103  1.00 26.04           C
HETATM 5597  O   SAR D   7      74.041  18.408  34.311  1.00 30.20           O
HETATM 5598  CN  SAR D   7      74.676  21.254  32.192  1.00 26.44           C
HETATM 5599  N   MLE D   8      74.259  17.657  32.181  1.00 24.57           N
HETATM 5600  CA  MLE D   8      75.247  16.647  32.580  1.00 20.62           C
HETATM 5601  C   MLE D   8      76.543  16.870  31.839  1.00 16.15           C
HETATM 5602  O   MLE D   8      76.648  16.595  30.650  1.00 12.75           O
HETATM 5603  CB  MLE D   8      74.601  15.270  32.287  1.00 17.53           C
HETATM 5604  CG  MLE D   8      73.424  14.778  33.148  1.00 15.59           C
HETATM 5605  CD1 MLE D   8      73.057  13.407  32.722  1.00 17.71           C
HETATM 5606  CD2 MLE D   8      73.719  14.808  34.624  1.00 10.51           C
HETATM 5607  CN  MLE D   8      73.778  17.716  30.745  1.00 23.20           C
ATOM   5608  N   VAL D   9      77.456  17.537  32.542  1.00 17.88           N
ATOM   5609  CA  VAL D   9      78.794  17.904  32.071  1.00 16.81           C
ATOM   5610  C   VAL D   9      79.626  16.627  31.883  1.00 16.38           C
ATOM   5611  O   VAL D   9      79.727  15.827  32.809  1.00 16.21           O
ATOM   5612  CB  VAL D   9      79.517  18.851  33.111  1.00 15.13           C
ATOM   5613  CG1 VAL D   9      80.783  19.458  32.515  1.00 20.21           C
ATOM   5614  CG2 VAL D   9      78.582  19.962  33.588  1.00 13.54           C
HETATM 5615  N   MLE D  10      80.162  16.357  30.670  1.00 14.56           N
HETATM 5616  CA  MLE D  10      80.971  15.154  30.508  1.00 16.29           C
HETATM 5617  C   MLE D  10      82.469  15.465  30.489  1.00 21.41           C
HETATM 5618  O   MLE D  10      82.856  16.602  30.732  1.00 25.90           O
HETATM 5619  CB  MLE D  10      80.435  14.483  29.240  1.00 16.62           C
HETATM 5620  CG  MLE D  10      79.006  13.929  29.309  1.00 13.42           C
HETATM 5621  CD1 MLE D  10      78.656  13.295  28.001  1.00 17.45           C
HETATM 5622  CD2 MLE D  10      78.870  12.956  30.446  1.00 16.36           C
HETATM 5623  CN  MLE D  10      80.129  17.398  29.564  1.00 13.53           C
ATOM   5624  N   ALA D  11      83.310  14.458  30.265  1.00 27.79           N
ATOM   5625  CA  ALA D  11      84.771  14.632  30.264  1.00 29.35           C
ATOM   5626  C   ALA D  11      85.398  15.312  29.048  1.00 31.45           C
ATOM   5627  O   ALA D  11      84.771  15.450  27.992  1.00 33.15           O
ATOM   5628  CB  ALA D  11      85.463  13.293  30.527  1.00 27.78           C
TER
END
"""

def get_hierarchy():
  return iotbx.pdb.input(lines=pdb_str.splitlines(),
                         source_info=None).construct_hierarchy()

def get_deviations(**kwds):
  validation = cbetadev.cbetadev(pdb_hierarchy=get_hierarchy(),
                                 out=null_out(),
                                 **kwds)
  return dict(((r.resseq.strip(), r.resname), r.deviation)
              for r in validation.results)

def tst_01():
  """
  Table residue names, with and without use_parent_for_modified_residues.
  """
  for resname in ['ALA', 'PRO', 'VAL', 'MSE']:
    assert cbd_utils.get_table_resname(resname) == resname
    assert cbd_utils.get_table_resname(resname,
      use_parent_for_modified_residues=True) == resname
  for resname, parent in [('HYP', 'PRO'),
                          ('TPO', 'THR'),
                          ('SEP', 'SER'),
                          ('MLY', 'LYS'),
                          ('MVA', 'VAL'),
                          ('MLE', 'LEU'),
                          ('BMT', 'THR'),
                          ]:
    assert cbd_utils.get_table_resname(resname) is None
    assert cbd_utils.get_table_resname(resname,
      use_parent_for_modified_residues=True) == parent, resname
  # no single parent, including D-amino acids
  for resname in ['DAL', 'DPR', 'MLU', 'UNK']:
    assert cbd_utils.get_table_resname(resname,
      use_parent_for_modified_residues=True) is None, resname
  # no table is not scored against the others table
  assert cbd_mahalanobis.get_bin_residual(None, (-60., -40.), 0.1, 0.) is None

def tst_02():
  """
  Modified and D-amino acids get phi/psi, so their neighbours do too.
  """
  phi_psi = cbetadev.get_phi_psi_dict(get_hierarchy())
  resseqs = set([key.split()[2] for key in phi_psi])
  # every residue except the two chain ends
  assert resseqs == set([str(i) for i in range(2, 11)]), resseqs
  # VAL 9 is between MLE 8 and MLE 10
  phi, psi = phi_psi['|pdbres="VAL D   9 ": |']
  assert approx_equal(phi, -67.1, eps=0.1)
  assert approx_equal(psi, 122.3, eps=0.1)
  # MLE 2 follows DAL 1
  phi, psi = phi_psi['|pdbres="MLE D   2 ": |']
  assert approx_equal(phi, -122.8, eps=0.1)
  assert approx_equal(psi, 75.5, eps=0.1)

def tst_03():
  """
  Corrections with and without use_parent_for_modified_residues.
  """
  default = get_deviations()
  off = get_deviations(apply_phi_psi_correction=True)
  on = get_deviations(apply_phi_psi_correction=True,
                      use_parent_for_modified_residues=True)
  positive_phi = get_deviations(apply_positive_phi_correction=True,
                                use_parent_for_modified_residues=True)
  assert len(default) == len(off) == len(on) == len(positive_phi) == 10
  # the standard residue next to modified residues is corrected either way
  assert approx_equal(default[('9', 'VAL')], 0.0619, eps=1e-4)
  assert approx_equal(off[('9', 'VAL')], 0.0562, eps=1e-4)
  assert approx_equal(on[('9', 'VAL')], 0.0562, eps=1e-4)
  # modified residues are only corrected with the parent table
  expected = {
    ('2', 'MLE')  : (0.1567, 0.1325),
    ('3', 'MLE')  : (0.1443, 0.1374),
    ('4', 'MVA')  : (0.0555, 0.0533),
    ('5', 'BMT')  : (0.1716, 0.1808),
    ('6', 'ABA')  : (0.0224, 0.0020),
    ('8', 'MLE')  : (0.1255, 0.1034),
    ('10', 'MLE') : (0.1558, 0.1379),
  }
  for key, (uncorrected, corrected) in expected.items():
    assert approx_equal(default[key], uncorrected, eps=1e-4), key
    assert approx_equal(off[key], uncorrected, eps=1e-4), key
    assert approx_equal(on[key], corrected, eps=1e-4), key
  # D-amino acids are never corrected; ALA 11 is the chain end
  for key in [('1', 'DAL'), ('11', 'ALA')]:
    assert approx_equal(on[key], default[key]), key
  # no residue here has 0 < phi < 150
  for key in default:
    assert approx_equal(positive_phi[key], default[key]), key

def tst_04():
  """
  use_parent_for_modified_residues needs a correction or Mahalanobis.
  """
  from iotbx.cli_parser import run_program
  from mmtbx.programs import cbetadev as cbetadev_program
  file_name = 'tst_cbetadev_03.pdb'
  with open(file_name, 'w') as f:
    f.write(pdb_str)
  try:
    run_program(program_class=cbetadev_program.Program,
                args=[file_name, 'use_parent_for_modified_residues=True'],
                logger=null_out())
  except Sorry as e:
    assert str(e).find('use_parent_for_modified_residues requires') > -1
  else:
    raise AssertionError('Sorry not raised')
  run_program(program_class=cbetadev_program.Program,
              args=[file_name,
                    'apply_phi_psi_correction=True',
                    'use_parent_for_modified_residues=True'],
              logger=null_out())

if (__name__ == "__main__"):
  tst_01()
  tst_02()
  tst_03()
  tst_04()
  print("OK")
