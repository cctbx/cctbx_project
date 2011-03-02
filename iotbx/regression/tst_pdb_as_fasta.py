
def exercise () :
  from iotbx.command_line import pdb_as_fasta
  pdb_str = """\
ATOM      2  CA  GLY A   3      -9.052   4.207   4.651  1.00 16.57           C
ATOM      6  CA  ASN A   4      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  ASN A   5      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN A   6       0.384   1.888   3.199  1.00 10.53           C
ATOM     22  CA  ALA A   6A      0.384   1.888   3.199  1.00 10.53           C
ATOM     22  CA  GLY A   6B      0.384   1.888   3.199  1.00 10.53           C
ATOM     31  CA  GLN A   7       3.270   2.361   5.640  1.00 11.39           C
ATOM     40  CA  ASN A   8       6.831   2.310   4.318  1.00 12.30           C
ATOM     48  CA  PTR A   9       9.159   2.144   7.299  1.00 15.18           C
TER
ATOM     50  O   HOH     1       0.000   0.000   0.000  1.0  20.0            O
END"""
  pdb_str_2 = """\
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  MSE A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN A   4       0.384   1.888   3.199  1.00 10.53           C
TER
ATOM     22  P   G  BB  83       0.000   0.000   0.000  1.00 10.53           P
ATOM     22  P   A  BB  84       0.000   0.000   0.000  1.00 10.53           P
ATOM     22  P   C  BB  87       0.000   0.000   0.000  1.00 10.53           P
ATOM     22  P   U  BB  88       0.000   0.000   0.000  1.00 10.53           P
TER
ATOM     90  O   HOH X   1       0.000   0.000   0.000  1.00 20.00           0
END
"""
  open("tmp1.pdb", "w").write(pdb_str)
  open("tmp2.pdb", "w").write(pdb_str_2)
  params = pdb_as_fasta.master_phil.fetch().extract()
  params.pdb_as_fasta.file_name.extend(["tmp1.pdb", "tmp2.pdb"])
  of = pdb_as_fasta.run(params=params)
  seq_in = open(of).read()
  assert (seq_in == """\
tmp1 chain ' A'
XXGNNQAGQNX
tmp2 chain ' A'
GNMQ
tmp2 chain 'BB'
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXGAXXCU""")
  params.pdb_as_fasta.pad_missing_residues=False
  of = pdb_as_fasta.run(params=params)
  seq_in = open(of).read()
  assert (seq_in == """\
tmp1 chain ' A'
GNNQAGQNX
tmp2 chain ' A'
GNMQ
tmp2 chain 'BB'
GACU""")
  params.pdb_as_fasta.include_insertion_residues = False
  of = pdb_as_fasta.run(params=params)
  seq_in = open(of).read()
  assert (seq_in == """\
tmp1 chain ' A'
GNNQQNX
tmp2 chain ' A'
GNMQ
tmp2 chain 'BB'
GACU""")
  params.pdb_as_fasta.pad_missing_residues = True
  params.pdb_as_fasta.include_insertion_residues = True
  params.pdb_as_fasta.ignore_missing_residues_at_start = True
  of = pdb_as_fasta.run(params=params)
  seq_in = open(of).read()
  assert (seq_in == """\
tmp1 chain ' A'
GNNQAGQNX
tmp2 chain ' A'
GNMQ
tmp2 chain 'BB'
GAXXCU""")

if (__name__ == "__main__") :
  exercise()
  print "OK"
