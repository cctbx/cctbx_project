from __future__ import division
import StringIO
from libtbx import easy_run
pdb = """
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
"""

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

geo_spec = """dihedral pdb=" CA  ASN A 411 "
         pdb=" C   ASN A 411 "
         pdb=" N   SER A 412 "
         pdb=" CA  SER A 412 "
    ideal   model   delta  harmonic     sigma   weight residual
  %4s.00"""

def cis_trans_specification():
  preamble = "bad_cis_peptide"
  f=file("%s.pdb" % preamble, "wb")
  f.write(pdb)
  f.close()
  for param, results in params.items():
    f=file("%s.params" % preamble, "wb")
    f.write(param)
    f.close()
    cmd = "phenix.geometry_minimization %(preamble)s.pdb %(preamble)s.params" % locals()
    print cmd
    ero = easy_run.fully_buffered(command=cmd)
    out = StringIO.StringIO()
    ero.show_stdout(out=out)

    lines = file("%(preamble)s_minimized.geo" % locals(), "rb").read()
    print geo_spec % results[0]
    if lines.find(geo_spec % results[0])==1:
      if lines.find(geo_spec % abs(results[0]))==1
        assert 0, ".geo specification not found"

def trans_only_specification():
  # must be run after cis_trans_specification
  preamble = "bad_cis_peptide"
  for arg, results in cmd_args.items():
    cmd = "phenix.geometry_minimization %(preamble)s_minimized.pdb %(arg)s" % locals()
    print cmd
    ero = easy_run.fully_buffered(command=cmd)
    out = StringIO.StringIO()
    ero.show_stdout(out=out)

    lines = file("%(preamble)s_minimized_minimized.geo" % locals(), "rb").read()
    print geo_spec % results[0]
    if lines.find(geo_spec % results[0])==1:
      if lines.find(geo_spec % abs(results[0]))==1
        assert 0, ".geo specification not found"

def run():
  cis_trans_specification()
  trans_only_specification()

if __name__=="__main__":
  run()#sys.argv[1])
