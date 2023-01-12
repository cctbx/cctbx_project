from __future__ import absolute_import, division, print_function
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.test_utils import assert_lines_in_text

def tst_1():
  pdb_str = """\
ATOM   3705  N   PRO A 839     -72.061   7.475  91.699  1.00 19.22           N
ATOM   3706  CA  PRO A 839     -71.162   6.864  90.749  1.00 15.72           C
ATOM   3707  C   PRO A 839     -70.566   7.946  89.845  1.00 19.35           C
ATOM   3708  O   PRO A 839     -71.127   9.016  89.760  1.00 16.63           O
ATOM   3709  CB  PRO A 839     -72.065   5.933  89.906  1.00 21.17           C
ATOM   3710  CG  PRO A 839     -73.211   5.597  90.858  1.00 24.24           C
ATOM   3711  CD  PRO A 839     -73.429   6.917  91.585  1.00 16.35           C
ATOM   3719  N  AMET A 840     -69.536   7.589  89.072  0.22 22.15           N
ATOM   3720  CA AMET A 840     -68.907   8.564  88.135  0.22 20.67           C
ATOM   3721  C  AMET A 840     -69.805   8.754  86.915  0.22 19.08           C
ATOM   3722  O  AMET A 840     -70.492   7.790  86.550  0.22 20.65           O
ATOM   3723  CB AMET A 840     -67.538   8.082  87.621  0.22 21.34           C
ATOM   3724  CG AMET A 840     -66.422   8.081  88.677  0.22 26.22           C
ATOM   3725  SD AMET A 840     -64.849   7.459  88.020  0.22 35.67           S
ATOM   3726  CE AMET A 840     -64.456   8.716  86.804  0.22 27.23           C
ATOM   3736  N  BMET A 840     -69.531   7.591  89.077  0.42 22.17           N
ATOM   3737  CA BMET A 840     -68.906   8.565  88.135  0.42 20.62           C
ATOM   3738  C  BMET A 840     -69.810   8.756  86.920  0.42 18.97           C
ATOM   3739  O  BMET A 840     -70.488   7.791  86.549  0.42 20.61           O
ATOM   3740  CB BMET A 840     -67.543   8.081  87.610  0.42 21.27           C
ATOM   3741  CG BMET A 840     -66.418   8.078  88.653  0.42 26.27           C
ATOM   3742  SD BMET A 840     -64.866   7.435  87.968  0.42 35.93           S
ATOM   3743  CE BMET A 840     -65.243   5.686  87.845  0.42 17.98           C
ATOM   3753  N  CMET A 840     -69.533   7.590  89.075  0.36 22.17           N
ATOM   3754  CA CMET A 840     -68.906   8.564  88.135  0.36 20.64           C
ATOM   3755  C  CMET A 840     -69.808   8.756  86.919  0.36 19.02           C
ATOM   3756  O  CMET A 840     -70.489   7.791  86.549  0.36 20.62           O
ATOM   3757  CB CMET A 840     -67.542   8.079  87.611  0.36 21.30           C
ATOM   3758  CG CMET A 840     -66.426   8.032  88.663  0.36 26.24           C
ATOM   3759  SD CMET A 840     -64.891   7.366  87.964  0.36 35.95           S
ATOM   3760  CE CMET A 840     -63.880   7.218  89.436  0.36 19.93           C
ATOM   3770  N   ASP A 841     -69.924   9.986  86.413  1.00 20.11           N
ATOM   3771  CA  ASP A 841     -70.719  10.250  85.211  1.00 23.77           C
ATOM   3772  C   ASP A 841     -72.158   9.760  85.347  1.00 26.73           C
ATOM   3773  O   ASP A 841     -72.719   9.126  84.444  1.00 26.93           O
ATOM   3774  CB  ASP A 841     -70.043   9.623  83.985  1.00 20.96           C
ATOM   3775  CG  ASP A 841     -68.601  10.057  83.840  1.00 38.39           C
ATOM   3776  OD1 ASP A 841     -68.357  11.284  83.832  1.00 51.85           O
ATOM   3777  OD2 ASP A 841     -67.711   9.179  83.749  1.00 53.14           O
"""
  h = iotbx.pdb.input(lines=pdb_str, source_info=None).construct_hierarchy()
  h_atoms = h.atoms()
  ogs = h.occupancy_groups_simple()
  # set occupancy values that will be rounded wrong
  test_occ = [0.225, 0.225, 0.55]
  for i, g in enumerate(ogs[0]):
    h_atoms.select(flex.size_t(g)).set_occ(flex.double([test_occ[i]]*len(g)))
  # observe wrong rounding
  wrong_pdb_str = h.as_pdb_string()
  assert_lines_in_text(wrong_pdb_str, """\
ATOM      8  N  AMET A 840     -69.536   7.589  89.072  0.23 22.15           N
ATOM      9  CA AMET A 840     -68.907   8.564  88.135  0.23 20.67           C
ATOM     10  C  AMET A 840     -69.805   8.754  86.915  0.23 19.08           C
ATOM     11  O  AMET A 840     -70.492   7.790  86.550  0.23 20.65           O
ATOM     12  CB AMET A 840     -67.538   8.082  87.621  0.23 21.34           C
ATOM     13  CG AMET A 840     -66.422   8.081  88.677  0.23 26.22           C
ATOM     14  SD AMET A 840     -64.849   7.459  88.020  0.23 35.67           S
ATOM     15  CE AMET A 840     -64.456   8.716  86.804  0.23 27.23           C
ATOM     16  N  BMET A 840     -69.531   7.591  89.077  0.23 22.17           N
ATOM     17  CA BMET A 840     -68.906   8.565  88.135  0.23 20.62           C
ATOM     18  C  BMET A 840     -69.810   8.756  86.920  0.23 18.97           C
ATOM     19  O  BMET A 840     -70.488   7.791  86.549  0.23 20.61           O
ATOM     20  CB BMET A 840     -67.543   8.081  87.610  0.23 21.27           C
ATOM     21  CG BMET A 840     -66.418   8.078  88.653  0.23 26.27           C
ATOM     22  SD BMET A 840     -64.866   7.435  87.968  0.23 35.93           S
ATOM     23  CE BMET A 840     -65.243   5.686  87.845  0.23 17.98           C
ATOM     24  N  CMET A 840     -69.533   7.590  89.075  0.55 22.17           N
ATOM     25  CA CMET A 840     -68.906   8.564  88.135  0.55 20.64           C
ATOM     26  C  CMET A 840     -69.808   8.756  86.919  0.55 19.02           C
ATOM     27  O  CMET A 840     -70.489   7.791  86.549  0.55 20.62           O
ATOM     28  CB CMET A 840     -67.542   8.079  87.611  0.55 21.30           C
ATOM     29  CG CMET A 840     -66.426   8.032  88.663  0.55 26.24           C
ATOM     30  SD CMET A 840     -64.891   7.366  87.964  0.55 35.95           S
ATOM     31  CE CMET A 840     -63.880   7.218  89.436  0.55 19.93           C
    """)
  h.round_occupancies_in_place(2)
  correct_pdb_str = h.as_pdb_string()
  assert_lines_in_text(correct_pdb_str, """\
ATOM      8  N  AMET A 840     -69.536   7.589  89.072  0.23 22.15           N
ATOM      9  CA AMET A 840     -68.907   8.564  88.135  0.23 20.67           C
ATOM     10  C  AMET A 840     -69.805   8.754  86.915  0.23 19.08           C
ATOM     11  O  AMET A 840     -70.492   7.790  86.550  0.23 20.65           O
ATOM     12  CB AMET A 840     -67.538   8.082  87.621  0.23 21.34           C
ATOM     13  CG AMET A 840     -66.422   8.081  88.677  0.23 26.22           C
ATOM     14  SD AMET A 840     -64.849   7.459  88.020  0.23 35.67           S
ATOM     15  CE AMET A 840     -64.456   8.716  86.804  0.23 27.23           C
ATOM     16  N  BMET A 840     -69.531   7.591  89.077  0.22 22.17           N
ATOM     17  CA BMET A 840     -68.906   8.565  88.135  0.22 20.62           C
ATOM     18  C  BMET A 840     -69.810   8.756  86.920  0.22 18.97           C
ATOM     19  O  BMET A 840     -70.488   7.791  86.549  0.22 20.61           O
ATOM     20  CB BMET A 840     -67.543   8.081  87.610  0.22 21.27           C
ATOM     21  CG BMET A 840     -66.418   8.078  88.653  0.22 26.27           C
ATOM     22  SD BMET A 840     -64.866   7.435  87.968  0.22 35.93           S
ATOM     23  CE BMET A 840     -65.243   5.686  87.845  0.22 17.98           C
ATOM     24  N  CMET A 840     -69.533   7.590  89.075  0.55 22.17           N
ATOM     25  CA CMET A 840     -68.906   8.564  88.135  0.55 20.64           C
ATOM     26  C  CMET A 840     -69.808   8.756  86.919  0.55 19.02           C
ATOM     27  O  CMET A 840     -70.489   7.791  86.549  0.55 20.62           O
ATOM     28  CB CMET A 840     -67.542   8.079  87.611  0.55 21.30           C
ATOM     29  CG CMET A 840     -66.426   8.032  88.663  0.55 26.24           C
ATOM     30  SD CMET A 840     -64.891   7.366  87.964  0.55 35.95           S
ATOM     31  CE CMET A 840     -63.880   7.218  89.436  0.55 19.93           C
    """)

if (__name__ == "__main__"):
  tst_1()
  print('OK')
