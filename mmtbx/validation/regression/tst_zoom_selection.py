"""
Validation results carry an atom selection for the residue or atoms they
describe.  The Phenix GUI sends it to the graphics programs when a row in an
outlier table is clicked; ChimeraX can only center on atoms, not on the
coordinates that Coot uses.
"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from libtbx.utils import format_cpu_times
import mmtbx.validation
from mmtbx.validation import cablam, ramalyze
from mmtbx.regression import model_1yjp

# i_seq 2-4 are residue 2, 5-6 the insertion 2A, 9-10 hybrid-36 A000 (10000),
# 11-12 a blank chain ID, 13-14 primed RNA atom names, 15-16 a negative resseq,
# 17-18 O5' next to an old-style O5*
pdb_str = """\
ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  GLY A   1       1.500   0.000   0.000  1.00 20.00           C
ATOM      3  N   ASN A   2       3.000   0.000   0.000  1.00 20.00           N
ATOM      4  CA  ASN A   2       4.500   0.000   0.000  1.00 20.00           C
ATOM      5  CB  ASN A   2       6.000   0.000   0.000  1.00 20.00           C
ATOM      6  N   ASN A   2A      7.500   0.000   0.000  1.00 20.00           N
ATOM      7  CA  ASN A   2A      9.000   0.000   0.000  1.00 20.00           C
ATOM      8  N   ALA A   3      10.500   0.000   0.000  1.00 20.00           N
ATOM      9  CA  ALA A   3      12.000   0.000   0.000  1.00 20.00           C
ATOM     10  N   GLY BA000      13.500   0.000   0.000  1.00 20.00           N
ATOM     11  CA  GLY BA000      15.000   0.000   0.000  1.00 20.00           C
ATOM     12  N   GLY     5      16.500   0.000   0.000  1.00 20.00           N
ATOM     13  CA  GLY     5      18.000   0.000   0.000  1.00 20.00           C
ATOM     14  O5'   U C   1      19.500   0.000   0.000  1.00 20.00           O
ATOM     15  C5'   U C   1      21.000   0.000   0.000  1.00 20.00           C
ATOM     16  N   GLY D  -3      22.500   0.000   0.000  1.00 20.00           N
ATOM     17  CA  GLY D  -3      24.000   0.000   0.000  1.00 20.00           C
ATOM     18  O5'   U C   2      25.500   0.000   0.000  1.00 20.00           O
ATOM     19  O5*   U C   2      27.000   0.000   0.000  1.00 20.00           O
"""

def _selected(cache, selection_string):
  assert selection_string is not None
  return list(cache.selection(selection_string).iselection())

def _cache():
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str)
  return pdb_in.construct_hierarchy().atom_selection_cache()

def exercise_residue():
  """A residue result selects that residue only: residue 2 not the
  insertion 2A, hybrid-36 and negative numbers, a blank chain ID."""
  cache = _cache()
  for chain_id, resseq, icode, expected in [
      ("A", "   2", " ", [2, 3, 4]),
      ("A", "   2", "A", [5, 6]),
      ("B", "A000", " ", [9, 10]),
      (" ", "   5", " ", [11, 12]),
      ("D", "  -3", " ", [15, 16]),
      ]:
    r = mmtbx.validation.residue(chain_id=chain_id, resseq=resseq,
      icode=icode, resname="GLY", altloc="", segid=None)
    assert _selected(cache, r.as_selection_string()) == expected, \
      (chain_id, resseq, icode, r.as_selection_string())

def exercise_atom():
  """A single-atom result selects that atom, including primed names; a *
  in an old-style name is not a wildcard."""
  cache = _cache()
  for chain_id, resseq, name, expected in [
      ("A", "   2", " CB ", [4]),
      ("C", "   1", " O5'", [13]),
      ("C", "   2", " O5*", [18]),
      ]:
    a = mmtbx.validation.atom(chain_id=chain_id, resseq=resseq, icode=" ",
      resname="ASN", altloc="", name=name)
    assert _selected(cache, a.as_selection_string()) == expected, \
      (name, a.as_selection_string())

def exercise_atoms():
  """A multi-atom result (bond, clash) selects every atom it involves."""
  cache = _cache()
  atoms_info = [
    mmtbx.validation.atom_info(chain_id="A", resseq="   2", icode=" ",
      resname="ASN", altloc="", name=" CA "),
    mmtbx.validation.atom_info(chain_id="A", resseq="   2", icode="A",
      resname="ASN", altloc="", name=" N  "),
  ]
  r = mmtbx.validation.atoms(atoms_info=atoms_info)
  assert _selected(cache, r.as_selection_string()) == [3, 5], \
    r.as_selection_string()

def exercise_restored_results():
  """Results restored from old pickles: py2-era slot objects keep bytes
  (fix_py2_pickle skips them), and slots added after the pickle was
  written are unset.  A residue without a number has nothing to select."""
  cache = _cache()
  r = mmtbx.validation.residue(chain_id=b"A", resseq=b"   2", icode=b"A",
    resname=b"ASN", altloc=b"", segid=None)
  assert _selected(cache, r.as_selection_string()) == [5, 6]
  r = mmtbx.validation.residue.__new__(mmtbx.validation.residue)
  r.__setstate__({"chain_id": "A", "resseq": "   3"})
  assert _selected(cache, r.as_selection_string()) == [7, 8]
  r = mmtbx.validation.residue(chain_id="A", resseq=None)
  assert r.as_selection_string() is None

def exercise_other_atom_records():
  """clashscore2 stores resseq as an int; elbow's carbohydrate results keep
  their own atom objects in atoms_info.  Neither may make zoom_info raise."""
  cache = _cache()
  a = mmtbx.validation.atom_info(chain_id="A", resseq=3, icode="",
    resname="ALA", altloc="", name="CA")
  r = mmtbx.validation.atoms(atoms_info=[a])
  assert _selected(cache, r.zoom_info()[0]) == [8]
  # numbers past 9999 are hybrid-36 in a selection, numpy integers too
  import numpy
  for chain_id, resseq, expected in [("B", 10000, [9, 10]),
                                     ("B", numpy.int64(10000), [9, 10]),
                                     ("D", -3, [15, 16])]:
    a = mmtbx.validation.atom_info(chain_id=chain_id, resseq=resseq,
      icode="", resname="GLY", altloc="", name="CA")
    r = mmtbx.validation.atoms(atoms_info=[a])
    assert _selected(cache, r.zoom_info()[0]) == expected[1:], (resseq,
      r.zoom_info()[0])
  class carbo_atom(object):
    chain_id, resseq, icode, name, resName = "A", "   3", " ", " N  ", "ALA"
  r = mmtbx.validation.atoms(atoms_info=[carbo_atom()])
  assert _selected(cache, r.zoom_info()[0]) == [7]

def exercise_table_rows():
  """Every row of a GUI outlier table carries its residue's selection."""
  pdb_in = iotbx.pdb.input(source_info=None, lines=model_1yjp)
  hierarchy = pdb_in.construct_hierarchy()
  cache = hierarchy.atom_selection_cache()
  atoms = hierarchy.atoms()
  rama = ramalyze.ramalyze(pdb_hierarchy=hierarchy, outliers_only=False)
  rows = rama.as_gui_table_data(outliers_only=False, include_zoom=True)
  resseqs = []
  for row in rows:
    isel = _selected(cache, row[-2])
    labels = set((atoms[i].fetch_labels().chain_id,
                  atoms[i].fetch_labels().resseq.strip()) for i in isel)
    assert len(labels) == 1, (row[-2], labels)
    resseqs.append(labels.pop())
  assert resseqs == [("A", "2"), ("A", "3"), ("A", "4"), ("A", "5"),
                     ("A", "6")], resseqs

def exercise_old_cablam_rows():
  """CaBLAM tables pickled by older versions have no selection, only the
  cells of as_table_row_phenix(): chain ' A', residue 'SER    3 '."""
  cache = _cache()
  for chain_cell, residue_cell, expected in [
      (" A", "ALA    3 ", [7, 8]),
      (" A", "ASN    2A", [5, 6]),
      ("B", "GLY A000 ", [9, 10]),
      ]:
    row = [chain_cell, residue_cell, "CaBLAM Outlier", 0.005, 0.009, "",
           0.0, 0.0, 0.0, None, (1.0, 2.0, 3.0)]
    sel = cablam.selection_from_table_row(row)
    assert _selected(cache, sel) == expected, (residue_cell, sel)

if (__name__ == "__main__"):
  exercise_residue()
  exercise_atom()
  exercise_atoms()
  exercise_restored_results()
  exercise_other_atom_records()
  exercise_table_rows()
  exercise_old_cablam_rows()
  print(format_cpu_times())
  print("OK")
