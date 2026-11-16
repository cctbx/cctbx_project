"""
The dictionary lookup in add_missing_H_atoms_at_bogus_position() is cached per
residue name and atom names. Placement with and without the cache must give the
same H and the same restraint objects.
Fragments: P, Q: the same ALA twice, one lookup. R, S: A from 1rrr with O2'
(RNA) and without it (DNA): the atom names choose the dictionary. E, F: BGM from
2e4i, in the CCD only: each residue passes its own auto-generated CIF on.
"""
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen

pdb_str = """\
CRYST1  120.000  120.000  120.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA P   1      11.204   2.315  10.312  1.00 20.00           N
ATOM      2  CA  ALA P   1      12.555   2.835  10.528  1.00 20.00           C
ATOM      3  C   ALA P   1      13.588   1.726  10.381  1.00 20.00           C
ATOM      4  O   ALA P   1      13.325   0.567  10.721  1.00 20.00           O
ATOM      5  CB  ALA P   1      12.866   3.965   9.556  1.00 20.00           C
ATOM      6  N   ALA Q   1      11.204  22.315  10.312  1.00 20.00           N
ATOM      7  CA  ALA Q   1      12.555  22.835  10.528  1.00 20.00           C
ATOM      8  C   ALA Q   1      13.588  21.726  10.381  1.00 20.00           C
ATOM      9  O   ALA Q   1      13.325  20.567  10.721  1.00 20.00           O
ATOM     10  CB  ALA Q   1      12.866  23.965   9.556  1.00 20.00           C
ATOM     11  P     A R   1      40.342  -5.036   5.198  1.00 20.00           P
ATOM     12  OP1   A R   1      41.732  -5.326   5.616  1.00 20.00           O
ATOM     13  OP2   A R   1      39.245  -5.144   6.185  1.00 20.00           O
ATOM     14  O5'   A R   1      40.301  -3.545   4.582  1.00 20.00           O
ATOM     15  C5'   A R   1      41.387  -2.992   3.859  1.00 20.00           C
ATOM     16  C4'   A R   1      41.046  -1.560   3.437  1.00 20.00           C
ATOM     17  O4'   A R   1      39.807  -1.514   2.753  1.00 20.00           O
ATOM     18  C3'   A R   1      40.900  -0.617   4.629  1.00 20.00           C
ATOM     19  O3'   A R   1      42.167  -0.168   5.083  1.00 20.00           O
ATOM     20  C2'   A R   1      40.048   0.493   4.017  1.00 20.00           C
ATOM     21  O2'   A R   1      40.862   1.442   3.353  1.00 20.00           O
ATOM     22  C1'   A R   1      39.217  -0.242   2.962  1.00 20.00           C
ATOM     23  N9    A R   1      37.794  -0.369   3.344  1.00 20.00           N
ATOM     24  C8    A R   1      37.080  -1.486   3.707  1.00 20.00           C
ATOM     25  N7    A R   1      35.798  -1.279   3.832  1.00 20.00           N
ATOM     26  C5    A R   1      35.646   0.067   3.520  1.00 20.00           C
ATOM     27  C6    A R   1      34.521   0.913   3.401  1.00 20.00           C
ATOM     28  N6    A R   1      33.267   0.490   3.594  1.00 20.00           N
ATOM     29  N1    A R   1      34.723   2.204   3.071  1.00 20.00           N
ATOM     30  C2    A R   1      35.965   2.631   2.857  1.00 20.00           C
ATOM     31  N3    A R   1      37.096   1.935   2.908  1.00 20.00           N
ATOM     32  C4    A R   1      36.863   0.637   3.246  1.00 20.00           C
ATOM     33  P     A S   1      40.342  24.964   5.198  1.00 20.00           P
ATOM     34  OP1   A S   1      41.732  24.674   5.616  1.00 20.00           O
ATOM     35  OP2   A S   1      39.245  24.856   6.185  1.00 20.00           O
ATOM     36  O5'   A S   1      40.301  26.455   4.582  1.00 20.00           O
ATOM     37  C5'   A S   1      41.387  27.008   3.859  1.00 20.00           C
ATOM     38  C4'   A S   1      41.046  28.440   3.437  1.00 20.00           C
ATOM     39  O4'   A S   1      39.807  28.486   2.753  1.00 20.00           O
ATOM     40  C3'   A S   1      40.900  29.383   4.629  1.00 20.00           C
ATOM     41  O3'   A S   1      42.167  29.832   5.083  1.00 20.00           O
ATOM     42  C2'   A S   1      40.048  30.493   4.017  1.00 20.00           C
ATOM     43  C1'   A S   1      39.217  29.758   2.962  1.00 20.00           C
ATOM     44  N9    A S   1      37.794  29.631   3.344  1.00 20.00           N
ATOM     45  C8    A S   1      37.080  28.514   3.707  1.00 20.00           C
ATOM     46  N7    A S   1      35.798  28.721   3.832  1.00 20.00           N
ATOM     47  C5    A S   1      35.646  30.067   3.520  1.00 20.00           C
ATOM     48  C6    A S   1      34.521  30.913   3.401  1.00 20.00           C
ATOM     49  N6    A S   1      33.267  30.490   3.594  1.00 20.00           N
ATOM     50  N1    A S   1      34.723  32.204   3.071  1.00 20.00           N
ATOM     51  C2    A S   1      35.965  32.631   2.857  1.00 20.00           C
ATOM     52  N3    A S   1      37.096  31.935   2.908  1.00 20.00           N
ATOM     53  C4    A S   1      36.863  30.637   3.246  1.00 20.00           C
HETATM   54  P   BGM E   1      68.359  -6.888   0.692  1.00 20.00           P
HETATM   55  OP1 BGM E   1      67.274  -7.004  -0.316  1.00 20.00           O
HETATM   56  OP2 BGM E   1      68.051  -7.102   2.128  1.00 20.00           O
HETATM   57  O5' BGM E   1      69.530  -7.884   0.276  1.00 20.00           O
HETATM   58  C5' BGM E   1      70.684  -8.044   1.100  1.00 20.00           C
HETATM   59  C4' BGM E   1      71.707  -8.912   0.405  1.00 20.00           C
HETATM   60  O4' BGM E   1      72.390  -8.119  -0.598  1.00 20.00           O
HETATM   61  C1' BGM E   1      73.786  -8.121  -0.343  1.00 20.00           C
HETATM   62  N9  BGM E   1      74.314  -6.811  -0.712  1.00 20.00           N
HETATM   63  C8  BGM E   1      74.952  -6.468  -1.881  1.00 20.00           C
HETATM   64  N7  BGM E   1      75.311  -5.214  -1.921  1.00 20.00           N
HETATM   65  C5  BGM E   1      74.886  -4.697  -0.705  1.00 20.00           C
HETATM   66  C4  BGM E   1      74.269  -5.670   0.053  1.00 20.00           C
HETATM   67  N3  BGM E   1      73.726  -5.562   1.283  1.00 20.00           N
HETATM   68  C2  BGM E   1      73.835  -4.336   1.766  1.00 20.00           C
HETATM   69  N2  BGM E   1      73.347  -4.052   2.982  1.00 20.00           N
HETATM   70  N1  BGM E   1      74.430  -3.296   1.093  1.00 20.00           N
HETATM   71  C6  BGM E   1      74.997  -3.385  -0.175  1.00 20.00           C
HETATM   72  O6  BGM E   1      75.510  -2.385  -0.693  1.00 20.00           O
HETATM   73  C2' BGM E   1      73.936  -8.472   1.123  1.00 20.00           C
HETATM   74  C3' BGM E   1      72.806  -9.468   1.312  1.00 20.00           C
HETATM   75  O3' BGM E   1      73.183 -10.779   0.874  1.00 20.00           O
HETATM   76 BR   BGM E   1      75.281  -7.696  -3.314  1.00 20.00          BR
HETATM   77  P   BGM F   1      68.359  23.112   0.692  1.00 20.00           P
HETATM   78  OP1 BGM F   1      67.274  22.996  -0.316  1.00 20.00           O
HETATM   79  OP2 BGM F   1      68.051  22.898   2.128  1.00 20.00           O
HETATM   80  O5' BGM F   1      69.530  22.116   0.276  1.00 20.00           O
HETATM   81  C5' BGM F   1      70.684  21.956   1.100  1.00 20.00           C
HETATM   82  C4' BGM F   1      71.707  21.088   0.405  1.00 20.00           C
HETATM   83  O4' BGM F   1      72.390  21.881  -0.598  1.00 20.00           O
HETATM   84  C1' BGM F   1      73.786  21.879  -0.343  1.00 20.00           C
HETATM   85  N9  BGM F   1      74.314  23.189  -0.712  1.00 20.00           N
HETATM   86  C8  BGM F   1      74.952  23.532  -1.881  1.00 20.00           C
HETATM   87  N7  BGM F   1      75.311  24.786  -1.921  1.00 20.00           N
HETATM   88  C5  BGM F   1      74.886  25.303  -0.705  1.00 20.00           C
HETATM   89  C4  BGM F   1      74.269  24.330   0.053  1.00 20.00           C
HETATM   90  N3  BGM F   1      73.726  24.438   1.283  1.00 20.00           N
HETATM   91  C2  BGM F   1      73.835  25.664   1.766  1.00 20.00           C
HETATM   92  N2  BGM F   1      73.347  25.948   2.982  1.00 20.00           N
HETATM   93  N1  BGM F   1      74.430  26.704   1.093  1.00 20.00           N
HETATM   94  C6  BGM F   1      74.997  26.615  -0.175  1.00 20.00           C
HETATM   95  O6  BGM F   1      75.510  27.615  -0.693  1.00 20.00           O
HETATM   96  C2' BGM F   1      73.936  21.528   1.123  1.00 20.00           C
HETATM   97  C3' BGM F   1      72.806  20.532   1.312  1.00 20.00           C
HETATM   98  O3' BGM F   1      73.183  19.221   0.874  1.00 20.00           O
HETATM   99 BR   BGM F   1      75.281  22.304  -3.314  1.00 20.00          BR
END
"""

def placed():
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  obj = reduce_hydrogen.place_hydrogens(model=model)
  obj.run()
  hs = []
  for a in obj.get_model().get_hierarchy().atoms():
    if a.element.strip() == "H":
      hs.append((a.parent().parent().parent().id, a.name.strip(), a.xyz))
  return sorted(hs), [name for name, cif in obj.restraint_objects]

def h_names(hs, chain):
  return sorted(h[1] for h in hs if h[0] == chain)

def uncached(residue, names, mon_lib_srv, cache):
  return reduce_hydrogen.mon_lib_query(residue=residue,
    mon_lib_srv=mon_lib_srv, raise_sorry=False)

def test_same_with_and_without_cache():
  saved = reduce_hydrogen._mon_lib_query_cached
  hs, objects = placed()
  reduce_hydrogen._mon_lib_query_cached = uncached
  try:
    hs_ref, objects_ref = placed()
  finally:
    reduce_hydrogen._mon_lib_query_cached = saved
  assert hs == hs_ref
  assert len(hs) == 60, len(hs)
  assert objects == objects_ref == ["auto_BGM", "auto_BGM"], objects
  # the reference itself: RNA and DNA differ at C2', the BGM are the same
  assert "HO2'" in h_names(hs_ref, "R") and "H2''" not in h_names(hs_ref, "R")
  assert "H2''" in h_names(hs_ref, "S") and "HO2'" not in h_names(hs_ref, "S")
  assert h_names(hs_ref, "P") == h_names(hs_ref, "Q")
  assert h_names(hs_ref, "E") == h_names(hs_ref, "F")

def test_repeated_residue_is_one_lookup():
  saved = reduce_hydrogen.mon_lib_query
  calls = []
  def counting(residue, **kwds):
    calls.append(residue.parent().parent().id)
    return saved(residue, **kwds)
  reduce_hydrogen.mon_lib_query = counting
  try:
    placed()
  finally:
    reduce_hydrogen.mon_lib_query = saved
  # add_missing_H first; place_n_terminal_propeller() has its own lookups after
  assert calls[:5] == ["P", "R", "S", "E", "F"], calls

if (__name__ == "__main__"):
  test_same_with_and_without_cache()
  test_repeated_residue_is_one_lookup()
  print("OK")
