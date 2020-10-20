from __future__ import absolute_import, division, print_function
import iotbx.pdb

pdb_str = """
SSBOND   1 CYS     30    CYS    115
CRYST1   28.174   52.857   68.929  90.00  90.00  90.00 P 21 21 21
ATOM    233  N   CYS    30      20.478 -18.192  14.641  1.00  7.59           N
ATOM    234  CA  CYS    30      20.362 -19.361  15.504  1.00  7.59           C
ATOM    235  C   CYS    30      20.010 -18.961  16.933  1.00  7.59           C
ATOM    236  O   CYS    30      19.363 -19.717  17.658  1.00  7.59           O
ATOM    237  CB  CYS    30      21.661 -20.169  15.489  1.00  7.59           C
ATOM    238  SG  CYS    30      23.138 -19.221  15.922  1.00  7.59           S
ATOM    884  N   CYS   115      24.507 -23.447  17.578  1.00  7.59           N
ATOM    885  CA  CYS   115      25.208 -22.644  16.583  1.00  7.59           C
ATOM    886  C   CYS   115      26.689 -23.004  16.526  1.00  7.59           C
ATOM    887  O   CYS   115      27.240 -23.233  15.449  1.00  7.59           O
ATOM    888  CB  CYS   115      25.040 -21.153  16.883  1.00  7.59           C
ATOM    889  SG  CYS   115      25.798 -20.055  15.662  1.00  7.59           S
ATOM    910  N   ASP   119      28.272 -25.962  10.279  1.00  7.59           N
ATOM    911  CA  ASP   119      28.924 -25.690   9.004  1.00  7.59           C
ATOM    912  C   ASP   119      27.987 -24.944   8.060  1.00  7.59           C
ATOM    913  O   ASP   119      26.943 -25.464   7.667  1.00  7.59           O
ATOM    914  CB  ASP   119      29.397 -26.993   8.355  1.00  7.59           C
ATOM    915  CG  ASP   119      30.145 -26.760   7.057  1.00  7.59           C
ATOM    916  OD1 ASP   119      30.431 -25.589   6.731  1.00  7.59           O
ATOM    917  OD2 ASP   119      30.453 -27.752   6.362  1.00  7.59           O
TER
END
"""

def run():
  ph = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
  asc = ph.atom_selection_cache()
  sel = asc.selection("name SG and (resseq 30 or resseq 115)")
  h = ph.select(~sel)
  h.flip_symmetric_amino_acids()

if (__name__ == "__main__"):
  run()
