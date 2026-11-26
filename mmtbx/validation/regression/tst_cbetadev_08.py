from __future__ import absolute_import, division, print_function
import io
from libtbx.utils import null_out
import iotbx.pdb
from mmtbx.validation import cbetadev

# Outliers removed / new outliers are counted per conformer, like the number
# of outliers, and each conformer is listed with its altloc.

# 3noq chain A 10-14 (as in tst_cbetadev_07) with GLU 12 split into A and B
# conformers that share the backbone. The CB of each is placed along the
# GLU 12 correction vector, so with strength 1 both cross the 0.25 A cutoff
# the same way: removed (0.28 -> 0.20, 0.30 -> 0.22) or new (0.20 -> 0.28,
# 0.22 -> 0.30).
pdb_template = """\
ATOM     78  N   PHE A  10       1.018  25.773  25.906  1.00  7.05           N
ATOM     79  CA  PHE A  10       1.266  26.317  27.227  1.00  6.77           C
ATOM     80  C   PHE A  10       2.741  26.229  27.573  1.00  5.97           C
ATOM     81  O   PHE A  10       3.483  25.403  27.044  1.00  6.02           O
ATOM     82  CB  PHE A  10       0.449  25.585  28.302  1.00  6.81           C
ATOM     89  N   PRO A  11       3.187  27.102  28.485  1.00  6.59           N
ATOM     90  CA  PRO A  11       4.592  27.045  28.938  1.00  6.67           C
ATOM     91  C   PRO A  11       4.963  25.654  29.424  1.00  5.60           C
ATOM     92  O   PRO A  11       4.219  25.024  30.171  1.00  6.09           O
ATOM     93  CB  PRO A  11       4.637  28.085  30.056  1.00  8.19           C
ATOM     96  N  AGLU A  12       6.144  25.180  29.020  0.50  5.53           N
ATOM     96  N  BGLU A  12       6.144  25.180  29.020  0.50  5.53           N
ATOM     97  CA AGLU A  12       6.663  23.872  29.432  0.50  5.74           C
ATOM     97  CA BGLU A  12       6.663  23.872  29.432  0.50  5.74           C
ATOM     98  C  AGLU A  12       5.862  22.715  28.829  0.50  5.27           C
ATOM     98  C  BGLU A  12       5.862  22.715  28.829  0.50  5.27           C
ATOM     99  O  AGLU A  12       5.905  21.597  29.338  0.50  5.87           O
ATOM     99  O  BGLU A  12       5.905  21.597  29.338  0.50  5.87           O
%s
ATOM    105  N   VAL A  13       5.218  22.955  27.680  1.00  5.26           N
ATOM    106  CA  VAL A  13       4.665  21.861  26.891  1.00  5.07           C
ATOM    107  C   VAL A  13       5.795  20.934  26.440  1.00  4.61           C
ATOM    108  O   VAL A  13       6.918  21.352  26.172  1.00  5.20           O
ATOM    109  CB  VAL A  13       3.897  22.443  25.681  1.00  5.43           C
ATOM    112  N   GLN A  14       5.451  19.653  26.286  1.00  4.78           N
ATOM    113  CA  GLN A  14       6.303  18.683  25.617  1.00  5.03           C
ATOM    114  C   GLN A  14       6.040  18.876  24.131  1.00  4.88           C
ATOM    115  O   GLN A  14       4.990  18.476  23.636  1.00  5.66           O
ATOM    116  CB  GLN A  14       5.945  17.271  26.043  1.00  5.89           C
END
"""

cb_removed = """\
ATOM    100  CB AGLU A  12       6.946  23.744  30.956  0.50  7.38           C
ATOM    100  CB BGLU A  12       6.966  23.743  30.956  0.50  7.38           C"""

cb_new = """\
ATOM    100  CB AGLU A  12       6.466  23.773  30.959  0.50  7.38           C
ATOM    100  CB BGLU A  12       6.447  23.775  30.959  0.50  7.38           C"""

def run_cbetadev(cb_lines, out=None):
  if out is None: out = null_out()
  h = iotbx.pdb.input(source_info=None,
    lines=(pdb_template % cb_lines).splitlines()).construct_hierarchy()
  return cbetadev.cbetadev(h, out=out, quiet=False,
                           apply_phi_psi_correction=True,
                           phi_psi_correction_strength=1.)

def glu12(c):
  return dict((r.altloc, r.deviation) for r in c.results
              if r.resseq.strip() == '12')

def exercise_both_conformers_removed():
  out = io.StringIO()
  c = run_cbetadev(cb_removed, out=out)
  devs = glu12(c)
  assert sorted(devs) == ['A', 'B']
  assert max(devs.values()) < 0.25
  assert [x[0].strip() for x in c.outliers_removed] == ['A  12 A', 'A  12 B']
  assert c.new_outliers == []
  text = out.getvalue()
  assert 'Outliers removed :     2   A  12 A, A  12 B' in text, text
  assert 'New outliers     :     0' in text, text

def exercise_both_conformers_new():
  out = io.StringIO()
  c = run_cbetadev(cb_new, out=out)
  devs = glu12(c)
  assert min(devs.values()) >= 0.25
  assert [x[0].strip() for x in c.new_outliers] == ['A  12 A', 'A  12 B']
  assert c.outliers_removed == []
  # counted per conformer, like the number of outliers
  assert c.n_outliers == len(c.new_outliers) == 2
  text = out.getvalue()
  assert 'New outliers     :     2   A  12 A, A  12 B' in text, text

def run():
  exercise_both_conformers_removed()
  exercise_both_conformers_new()
  print('OK')

if __name__ == '__main__':
  run()
