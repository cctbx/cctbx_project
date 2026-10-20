"""place_hydrogens on a model returned by placement or reduce2 gives the same
result as on that model read from a file. Residues that get H but have no
restraints are listed."""
from __future__ import absolute_import, division, print_function
import iotbx.cif
import iotbx.pdb
import mmtbx.model
import mmtbx.monomer_library.server
from libtbx.utils import null_out
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from mmtbx.hydrogens import reduce_hydrogen

cryst1 = "CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1\n"

# 1yjp residues 1-3
pdb_protein = """\
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
"""

# BGM: in the CCD only (no geostd entry), 13 H
pdb_bgm = """\
HETATM   25  P   BGM C   1      -9.873  19.915  28.793  1.00 20.00           P
HETATM   26  OP1 BGM C   1     -10.967  19.590  29.736  1.00 20.00           O
HETATM   27  OP2 BGM C   1      -9.665  21.510  28.742  1.00 20.00           O
HETATM   28  O5' BGM C   1      -8.512  19.211  29.289  1.00 20.00           O
HETATM   29  C5' BGM C   1      -7.969  19.403  30.597  1.00 20.00           C
HETATM   30  C4' BGM C   1      -6.677  18.594  30.735  1.00 20.00           C
HETATM   31  O4' BGM C   1      -5.653  19.118  29.863  1.00 20.00           O
HETATM   32  C1' BGM C   1      -4.389  18.682  30.389  1.00 20.00           C
HETATM   33  N9  BGM C   1      -3.370  19.708  30.154  1.00 20.00           N
HETATM   34  C8  BGM C   1      -3.522  21.053  30.334  1.00 20.00           C
HETATM   35  N7  BGM C   1      -2.414  21.666  30.030  1.00 20.00           N
HETATM   36  C5  BGM C   1      -1.490  20.754  29.638  1.00 20.00           C
HETATM   37  C4  BGM C   1      -2.090  19.493  29.718  1.00 20.00           C
HETATM   38  N3  BGM C   1      -1.388  18.403  29.389  1.00 20.00           N
HETATM   39  C2  BGM C   1      -0.138  18.494  28.992  1.00 20.00           C
HETATM   40  N2  BGM C   1       0.542  17.348  28.664  1.00 20.00           N
HETATM   41  N1  BGM C   1       0.501  19.694  28.896  1.00 20.00           N
HETATM   42  C6  BGM C   1      -0.147  20.839  29.205  1.00 20.00           C
HETATM   43  O6  BGM C   1       0.417  21.917  29.124  1.00 20.00           O
HETATM   44  C2' BGM C   1      -4.574  18.465  31.904  1.00 20.00           C
HETATM   45  C3' BGM C   1      -6.071  18.772  32.151  1.00 20.00           C
HETATM   46  O3' BGM C   1      -6.639  17.837  33.070  1.00 20.00           O
HETATM   47  OP3 BGM C   1     -10.252  19.375  27.325  1.00 20.00           O
HETATM   48  BR  BGM C   1      -5.103  21.898  30.935  1.00 20.00          BR
"""

# ethanol under a code that is not in the CCD
pdb_04 = """\
HETATM   49  C1   04 B   1      15.000  10.000  15.000  1.00 20.00           C
HETATM   50  C2   04 B   1      16.520  10.000  15.000  1.00 20.00           C
HETATM   51  O1   04 B   1      17.220  11.210  15.000  1.00 20.00           O
"""

# its dictionary: atoms only, no bonds
cif_04_no_bonds = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 04  04  'test ligand'  ligand  9  3  .

data_comp_04
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
 04  C1   C  CH3   0
 04  C2   C  CH2   0
 04  O1   O  OH1   0
 04  H11  H  HCH3  0
 04  H12  H  HCH3  0
 04  H13  H  HCH3  0
 04  H21  H  HCH2  0
 04  H22  H  HCH2  0
 04  HO1  H  HOH1  0
"""

cases = {
  "protein_bgm": cryst1 + pdb_protein + pdb_bgm,
  "bgm":         cryst1 + pdb_bgm,
}

def from_text(text, restraint_objects=None):
  return mmtbx.model.manager(
    model_input=iotbx.pdb.input(source_info=None, lines=text),
    restraint_objects=restraint_objects, stop_for_unknowns=False,
    log=null_out())

def place(model, **kwargs):
  obj = reduce_hydrogen.place_hydrogens(model=model, **kwargs)
  obj.run()
  return obj

def n_h(model, resname):
  return sum(1 for a in model.get_hierarchy().atoms()
    if a.parent().resname.strip() == resname and a.element.strip() == "H")

def atoms(model):
  return [(a.id_str(), "%.3f %.3f %.3f" % a.xyz)
          for a in model.get_hierarchy().atoms()]

def check_reuse(model, keep, label):
  from_file = place(from_text(model.model_as_pdb()),
    keep_existing_H=keep).get_model()
  assert n_h(from_file, "BGM") == 13, label
  for m in (model, model.select(flex.bool(model.size(), True)),
            model.deep_copy()):
    placed = place(m, keep_existing_H=keep).get_model()
    assert atoms(placed) == atoms(from_file), (label, n_h(placed, "BGM"))

def test_placement_twice():
  for name, text in cases.items():
    for keep in (False, True):
      returned = place(from_text(text), keep_existing_H=keep).get_model()
      check_reuse(returned, keep, (name, keep))

def test_after_reduce2():
  from iotbx.cli_parser import CCTBXParser
  from mmtbx.programs import reduce2
  for name, text in cases.items():
    prefix = "tst_add_hydrogen_19_" + name
    with open(prefix + ".pdb", "w") as f:
      f.write(text)
    parser = CCTBXParser(program_class=reduce2.Program, logger=null_out())
    parser.parse_args([prefix + ".pdb", "output.filename=%s_out.pdb" % prefix,
      "output.description_file_name=%s.txt" % prefix, "output.overwrite=True"])
    prog = reduce2.Program(parser.data_manager, parser.working_phil.extract(),
      logger=StringIO())
    prog.validate()
    prog.run()
    for keep in (False, True):
      check_reuse(prog.get_results().model, keep, (name, keep))

def test_no_restraints_listed():
  ro = [("04.cif", iotbx.cif.reader(input_string=cif_04_no_bonds).model())]
  obj = place(from_text(cryst1 + pdb_protein + pdb_04, ro))
  assert n_h(obj.get_model(), "04") == 0
  assert obj.get_counts().residues_no_restraints == ["04 B 1"]
  assert "04" not in obj.no_H_placed_mlq
  log = StringIO()
  obj.show(log)
  assert ("No H atoms were placed on the following residues because the "
    "processed model has no restraints for them:\n04 B 1\n") in log.getvalue()
  # truncated residue: not listed
  text = "".join(l for l in (cryst1 + pdb_protein).splitlines(True)
                 if " ND2 ASN A   3" not in l)
  assert place(from_text(text)).get_counts().residues_no_restraints == []

class counting_server(mmtbx.monomer_library.server.server):
  n_queries = 0
  def get_comp_comp_id_and_atom_name_interpretation(self, **kwargs):
    self.n_queries += 1
    return mmtbx.monomer_library.server.server.\
      get_comp_comp_id_and_atom_name_interpretation(self, **kwargs)

def test_given_server():
  srv = counting_server()
  place(from_text(cases["protein_bgm"]), mon_lib_srv=srv)
  assert srv.n_queries > 0

if (__name__ == "__main__"):
  test_placement_twice()
  test_after_reduce2()
  test_no_restraints_listed()
  test_given_server()
  print("OK")
