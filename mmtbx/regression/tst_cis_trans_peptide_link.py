from __future__ import absolute_import, division, print_function
from six.moves import cStringIO as StringIO
from libtbx import easy_run
from six.moves import range
pdbs = ["""
ATOM   3239  N   ASN A 411       8.430  37.928 107.306  1.00 14.13           N
ATOM   3240  CA  ASN A 411       9.390  38.201 106.233  1.00 15.21           C
ATOM   3241  C   ASN A 411      10.787  38.396 106.824  1.00 16.78           C
ATOM   3242  O   ASN A 411      11.038  39.457 107.406  1.00 16.47           O
ATOM   3243  CB  ASN A 411       8.929  39.444 105.442  1.00 14.52           C
ATOM   3244  CG  ASN A 411       9.920  39.867 104.360  1.00 14.85           C
ATOM   3245  OD1 ASN A 411      10.958  39.246 104.178  1.00 16.53           O
ATOM   3246  ND2 ASN A 411       9.605  40.952 103.654  1.00 14.53           N
ATOM   3247  N   SER A 412      11.326  37.241 107.164  1.00 18.80           N
ATOM   3248  CA  SER A 412      11.922  36.184 106.321  1.00 17.93           C
ATOM   3249  C   SER A 412      13.005  35.451 107.046  1.00 19.48           C
ATOM   3250  O   SER A 412      13.620  35.982 107.956  1.00 20.17           O
ATOM   3251  CB  SER A 412      12.279  36.520 104.861  1.00 18.82           C
ATOM   3252  OG  SER A 412      13.554  37.143 104.791  1.00 18.04           O
ATOM   3253  N   ALA A 413      13.170  34.187 106.696  1.00 19.01           N
ATOM   3254  CA  ALA A 413      14.295  33.427 107.221  1.00 20.12           C
ATOM   3255  C   ALA A 413      15.125  32.912 106.082  1.00 20.71           C
ATOM   3256  O   ALA A 413      14.595  32.411 105.098  1.00 20.46           O
ATOM   3257  CB  ALA A 413      13.795  32.257 108.046  1.00 19.57           C
""",
  """
ATOM     71  N   SER A 411      39.662  65.271  37.503  1.00 61.40           N
ATOM     72  CA  SER A 411      38.219  65.068  37.450  1.00 59.89           C
ATOM     73  C   SER A 411      37.919  63.580  37.301  1.00 57.15           C
ATOM     74  O   SER A 411      38.565  62.895  36.507  1.00 53.39           O
ATOM     75  CB  SER A 411      37.599  65.840  36.289  1.00 55.80           C
ATOM     76  N   PRO A 412      37.102  63.122  37.964  1.00 59.58           N
ATOM     77  CA  PRO A 412      36.829  61.665  37.929  1.00 63.69           C
ATOM     78  C   PRO A 412      35.870  61.273  36.808  1.00 63.21           C
ATOM     79  O   PRO A 412      34.706  60.913  37.019  1.00 68.01           O
ATOM     80  CB  PRO A 412      36.236  61.415  39.313  1.00 64.16           C
ATOM     81  CG  PRO A 412      35.516  62.685  39.620  1.00 68.82           C
ATOM     82  CD  PRO A 412      36.321  63.804  39.007  1.00 64.77           C
ATOM     83  N   THR A 413      36.378  61.338  35.577  1.00 65.04           N
ATOM     84  CA  THR A 413      35.605  61.101  34.371  1.00 72.28           C
ATOM     85  C   THR A 413      36.415  60.173  33.463  1.00 76.81           C
ATOM     86  O   THR A 413      37.017  60.604  32.483  1.00 75.12           O
ATOM     87  CB  THR A 413      35.271  62.418  33.666  1.00 75.72           C
ATOM     88  OG1 THR A 413      36.491  63.024  33.216  1.00 83.26           O
ATOM     89  CG2 THR A 413      34.574  63.397  34.618  1.00 68.25           C
""",
]

params = {"" : [-180],
          """
pdb_interpretation {
 apply_cis_trans_specification {
   cis_trans_mod = cis
   residue_selection = chain a and resseq 412
  }
}
""" : [0],
          """
pdb_interpretation {
 apply_cis_trans_specification {
   cis_trans_mod = cis
   residue_selection = chain a and resseq 412 and name CA
  }
}
""" : [0],
  }

cmd_args = {
  "" : [0],
  "peptide_link.apply_all_trans=False" : [0],
  "peptide_link.apply_all_trans=True" : [-180],
  }

geo_specs = ["""dihedral pdb=" CA  ASN A 411 "
         pdb=" C   ASN A 411 "
         pdb=" N   SER A 412 "
         pdb=" CA  SER A 412 "
    ideal   model   delta  harmonic     sigma   weight residual
  %4s.00""",
  """dihedral pdb=" CA  SER A 411 "
         pdb=" C   SER A 411 "
         pdb=" N   PRO A 412 "
         pdb=" CA  PRO A 412 "
    ideal   model   delta  harmonic     sigma   weight residual
  %4s.00"""
  ]

def cis_trans_specification():
  for i, pdb in enumerate(pdbs):
    preamble = "bad_cis_peptide_%02d" % i
    f=open("%s.pdb" % preamble, "w")
    f.write(pdb)
    f.close()
    for param, results in params.items():
      f=open("%s.params" % preamble, "w")
      f.write(param)
      f.close()
      cmd = "phenix.geometry_minimization %(preamble)s.pdb %(preamble)s.params" % locals()
      print(cmd)
      ero = easy_run.fully_buffered(command=cmd)
      out = StringIO()
      ero.show_stdout(out=out)

      with open("%(preamble)s_minimized.geo" % locals(), "r") as f:
        lines = f.read()
      geo_spec=geo_specs[i]
      print(geo_spec % results[0])
      if lines.find(geo_spec % results[0])==1:
        if lines.find(geo_spec % abs(results[0]))==1:
          assert 0, ".geo specification not found"

def trans_only_specification():
  # must be run after cis_trans_specification
  for i in range(len(pdbs)):
    preamble = "bad_cis_peptide_%02d" % i
    for arg, results in cmd_args.items():
      cmd = "phenix.geometry_minimization %(preamble)s_minimized.pdb %(arg)s" % locals()
      print(cmd)
      ero = easy_run.fully_buffered(command=cmd)
      out = StringIO()
      ero.show_stdout(out=out)

      with open("%(preamble)s_minimized_minimized.geo" % locals(), "r") as f:
        lines = f.read()
      geo_spec=geo_specs[i]
      print(geo_spec % results[0])
      if lines.find(geo_spec % results[0])==1:
        if lines.find(geo_spec % abs(results[0]))==1:
          assert 0, ".geo specification not found"

def run():
  cis_trans_specification()
  trans_only_specification()

if __name__=="__main__":
  run()#sys.argv[1])
