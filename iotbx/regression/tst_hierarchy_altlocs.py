from __future__ import absolute_import, division, print_function
import iotbx.pdb
from six.moves import cStringIO as StringIO
from libtbx.test_utils import show_diff

pdb_str1 = """
ATOM      1  N  AVAL A   1      -4.898   0.072  13.387  0.49  7.34           N
ATOM      2  CA AVAL A   1      -4.626   0.703  12.080  0.49  7.71           C
ATOM      3  C  AVAL A   1      -3.475   1.680  12.227  0.49  7.52           C
ATOM      4  O  AVAL A   1      -3.125   2.100  13.335  0.49  7.59           O
ATOM      5  CB AVAL A   1      -5.882   1.390  11.495  0.49  7.79           C
ATOM      6  CG1AVAL A   1      -6.979   0.367  11.258  0.49  7.20           C
ATOM      7  CG2AVAL A   1      -6.359   2.496  12.408  0.49  7.81           C
ATOM      8  H1 AVAL A   1      -4.794  -0.809  13.318  0.49  8.81           H
ATOM      9  H2 AVAL A   1      -4.331   0.391  13.995  0.49  8.81           H
ATOM     10  H3 AVAL A   1      -5.733   0.254  13.636  0.49  8.81           H
ATOM     11  HA AVAL A   1      -4.371   0.017  11.444  0.49  9.26           H
ATOM     12  HB AVAL A   1      -5.655   1.790  10.641  0.49  9.35           H
HETATM 2211  O   HOH S 216       9.105  -6.647  -4.343  0.25 56.37           O
HETATM 2219  O  BHOH S 224       6.977   3.045   9.044  0.31 55.60           O
"""

pdb_str2 = """
HETATM 2174  O   HOH S 179      -5.781   7.569   9.276  0.24 14.70           O
HETATM 2174  O  AHOH S 179      -5.781   7.569   9.276  0.24 14.70           O
HETATM 2174  O   HOH S 179      -4.781   7.569   9.276  0.24 14.70           H
HETATM 2174  O  CHOH S 179      -6.781   7.569   9.276  0.24 14.70           H
END
"""

pdb_str3 = """
HETATM 2174  O   HOH S 179      -5.781   7.569   9.276  0.24 14.70           O
HETATM 2174  O   HOH S 179      -5.781   7.569   9.276  0.24 14.70           O
HETATM 2174  O   HOH S 179      -4.781   7.569   9.276  0.24 14.70           H
HETATM 2174  O   HOH S 179      -6.781   7.569   9.276  0.24 14.70           H
END
"""

pdb_str4 = """
ATOM      1  N   ALA A   1      -4.898   0.072  13.387  0.49  7.34           N
ATOM      2  CA  ALA A   1      -4.626   0.703  12.080  0.49  7.71           C
ATOM      3  C   ALA A   1      -3.475   1.680  12.227  0.49  7.52           C
ATOM      4  O   ALA A   1      -3.125   2.100  13.335  0.49  7.59           O
ATOM      5  CB  ALA A   1      -5.882   1.390  11.495  0.49  7.79           C
ATOM      5  CB AALA A   1      -5.882   1.390  11.495  0.49  7.79           C
"""

pdb_str5 = """
ATOM      1  N   ALA A   1      -4.898   0.072  13.387  0.49  7.34           N
ATOM      2  CA  ALA A   1      -4.626   0.703  12.080  0.49  7.71           C
ATOM      3  C   ALA A   1      -3.475   1.680  12.227  0.49  7.52           C
ATOM      4  O   ALA A   1      -3.125   2.100  13.335  0.49  7.59           O
ATOM      5  CB AALA A   1      -5.882   1.390  11.495  0.49  7.79           C
ATOM      5  CB BALA A   1      -5.882   1.390  11.495  0.49  7.79           C
"""

pdb_str6 = """
HETATM 2174  O   HOH S 179      -5.781   7.569   9.276  0.24 14.70           O
"""
def get_h_oc(s):
  h = iotbx.pdb.input(source_info=None, lines=s).construct_hierarchy()
  oc = h.overall_counts()
  return h, oc

#
# This is current state of affairs. This was also the case in
# Phenix 1.20 and 1.18.
#

def tst1():
  """Everything is good here, empty altloc as expected.
  """
  h, oc = get_h_oc(pdb_str1)
  altlocs = [ag.altloc for ag in h.only_model().atom_groups()]
  assert altlocs == ['A', '', 'B'], altlocs
  assert oc.n_alt_conf_improper == 0, oc.n_alt_conf_improper
  assert oc.duplicate_atom_labels == [], list(oc.duplicate_atom_labels)
  #

def tst2():
  """There's duplicated atom - HOH without altloc. It gets whitespace: ' ',
  but overall_counts has several warnings.
  """
  h, oc = get_h_oc(pdb_str2)
  altlocs = [ag.altloc for ag in h.only_model().atom_groups()]
  assert altlocs == [' ', 'A', 'C'], altlocs
  # for a in h.atoms():
  #   print("%s '%s'" % (a.id_str(), a.parent().altloc))

  # Things that let one know about improper alt conf:
  assert oc.n_alt_conf_improper == 1, oc.n_alt_conf_improper
  assert len(oc.duplicate_atom_labels) == 1, list(oc.duplicate_atom_labels)

  duplicate_output = StringIO()
  oc.show_duplicate_atom_labels(out=duplicate_output)
  do_value = duplicate_output.getvalue()
  assert not show_diff(do_value, """\
number of groups of duplicate atom labels: 1
  total number of affected atoms:          2
  group "HETA    .*.  O   HOH S 179 .*.     O  "
        "HETA    .*.  O   HOH S 179 .*.     H  "
""")

def tst3():
  """ This is a simple case of duplicated atom labels:
  no n_alt_conf_improper
  """
  h, oc = get_h_oc(pdb_str3)
  altlocs = [ag.altloc for ag in h.only_model().atom_groups()]
  assert altlocs == [''], altlocs
  assert oc.n_alt_conf_improper == 0, oc.n_alt_conf_improper
  assert len(oc.duplicate_atom_labels) == 1, list(oc.duplicate_atom_labels)

  duplicate_output = StringIO()
  oc.show_duplicate_atom_labels(out=duplicate_output)
  do_value = duplicate_output.getvalue()
  assert not show_diff(do_value, """\
number of groups of duplicate atom labels: 1
  total number of affected atoms:          4
  group "HETA    .*.  O   HOH S 179 .*.     O  "
        "HETA    .*.  O   HOH S 179 .*.     O  "
        "HETA    .*.  O   HOH S 179 .*.     H  "
        "HETA    .*.  O   HOH S 179 .*.     H  "
""")

def tst4():
  """This is pure case where only n_alt_conf_improper shows problem,
  but no duplicated atom labels.
  """
  # Here we have a residue with one atom duplicated out of several.
  # CB atoms get ' ' and 'A' altlocs, while the rest get ''.
  h, oc = get_h_oc(pdb_str4)
  altlocs = [ag.altloc for ag in h.only_model().atom_groups()]
  # for a in h.atoms():
  #   print("%s '%s'" % (a.id_str(), a.parent().altloc))
  assert altlocs == ['', ' ', 'A'], altlocs
  assert oc.n_alt_conf_improper == 1, oc.n_alt_conf_improper
  assert len(oc.duplicate_atom_labels) == 0, list(oc.duplicate_atom_labels)

  duplicate_output = StringIO()
  oc.show_duplicate_atom_labels(out=duplicate_output)
  do_value = duplicate_output.getvalue()
  assert not show_diff(do_value, "")

def tst5():
  """putting whitespace into already created hierarchy results
  in exactly same behavior as having incorrect one:
  here we transform good hierarchy into one from tst4 by changing
  altloc A to ' ' and get the same n_alt_conf_improper == 1."""
  h, oc = get_h_oc(pdb_str5)
  altlocs = [ag.altloc for ag in h.only_model().atom_groups()]
  assert altlocs == ['', 'A', 'B'], altlocs
  assert oc.n_alt_conf_improper == 0, oc.n_alt_conf_improper
  assert len(oc.duplicate_atom_labels) == 0, list(oc.duplicate_atom_labels)

  # Now we change altloc A to ' '
  h.atoms()[-2].parent().altloc = ' '
  # h.atoms()[-1].parent().altloc = ' '

  altlocs = [ag.altloc for ag in h.only_model().atom_groups()]
  assert altlocs == ['', ' ', 'B'], altlocs
  # print(h.as_pdb_string())
  oc2 = h.overall_counts()
  assert oc2.n_alt_conf_improper == 1, oc2.n_alt_conf_improper
  assert len(oc2.duplicate_atom_labels) == 0, list(oc2.duplicate_atom_labels)

def tst6():
  """Putting whitespace into altloc seem to always result
  in n_alt_conf_improper.
  """
  h, oc = get_h_oc(pdb_str6)
  altlocs = [ag.altloc for ag in h.only_model().atom_groups()]
  assert altlocs == [''], altlocs
  assert oc.n_alt_conf_improper == 0, oc.n_alt_conf_improper
  assert len(oc.duplicate_atom_labels) == 0, list(oc.duplicate_atom_labels)

  # Now we change altloc A to ' '
  h.atoms()[0].parent().altloc = ' '
  altlocs = [ag.altloc for ag in h.only_model().atom_groups()]
  assert altlocs == [' '], altlocs
  oc2 = h.overall_counts()
  assert oc2.n_alt_conf_improper == 1, oc2.n_alt_conf_improper
  assert len(oc2.duplicate_atom_labels) == 0, list(oc2.duplicate_atom_labels)


if (__name__ == "__main__"):
  tst1()
  tst2()
  tst3()
  tst4()
  tst5()
  tst6()
  print("OK")
