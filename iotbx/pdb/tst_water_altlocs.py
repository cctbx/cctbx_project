"""iotbx.pdb.water_altlocs: altloc labels for unlabelled partial-occupancy
alternates around waters; labels only, nothing else changes."""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from iotbx.pdb import water_altlocs

CRYST1 = "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1\n"

def water(serial, resseq, x, y, z, occ, altloc=" "):
  return ("HETATM%5d  O  %sHOH A%4d    %8.3f%8.3f%8.3f%6.2f 20.00           O"
          % (serial, altloc, resseq, x, y, z, occ))

def hierarchy(lines, cryst1=CRYST1):
  inp = iotbx.pdb.input(source_info=None,
    lines=(cryst1 + "\n".join(lines) + "\nEND\n").splitlines())
  return inp.construct_hierarchy(), inp.crystal_symmetry()

def altlocs(h):
  return dict(((a.parent().parent().resseq.strip(), a.name.strip()),
               a.parent().altloc.strip()) for a in h.atoms())

def test_water_pair():
  h, cs = hierarchy([water(1, 1, 5, 5, 5, 0.5), water(2, 2, 5.8, 5, 5, 0.5),
                     water(3, 3, 12, 12, 12, 0.5), water(4, 4, 15, 15, 15, 1.0)])
  before = [(a.xyz, a.occ, a.b) for a in h.atoms()]
  r = water_altlocs.run(h, cs)
  assert r.n_labelled == 2 and r.manual == [], (r.n_labelled, r.manual)
  assert r.labelled == [("HOH A 1", "A"), ("HOH A 2", "B")], r.labelled
  a = altlocs(h)
  assert a[("1", "O")] == "A" and a[("2", "O")] == "B", a
  # an isolated partial water and a full one stay blank
  assert a[("3", "O")] == "" and a[("4", "O")] == "", a
  assert [(x.xyz, x.occ, x.b) for x in h.atoms()] == before
  # a second run changes nothing
  assert water_altlocs.run(h, cs).n_labelled == 0

def test_three_waters():
  h, cs = hierarchy([water(1, 1, 5, 5, 5, 0.3), water(2, 2, 5.8, 5, 5, 0.4),
                     water(3, 3, 5.4, 5.7, 5, 0.3)])
  r = water_altlocs.run(h, cs)
  assert r.n_labelled == 3, r.labelled
  a = altlocs(h)
  # the highest occupancy is A; three mutual overlaps need three labels
  assert a[("2", "O")] == "A", a
  assert sorted(a.values()) == ["A", "B", "C"], a

def test_full_occupancy_untouched():
  h, cs = hierarchy([water(1, 1, 5, 5, 5, 1.0), water(2, 2, 5.8, 5, 5, 1.0),
                     water(3, 3, 8, 8, 8, 1.0), water(4, 4, 8.5, 8, 8, 0.5)])
  r = water_altlocs.run(h, cs)
  assert r.n_labelled == 0 and r.manual == [], (r.labelled, r.manual)

def test_symmetry_mate():
  # 0.6 A apart through the lattice translation along x
  h, cs = hierarchy([water(1, 1, 0.3, 5, 5, 0.5), water(2, 2, 19.7, 5, 5, 0.5)])
  r = water_altlocs.run(h, cs)
  assert r.n_labelled == 2, r.labelled
  # and without symmetry (a box) they are 19.4 A apart: not alternates
  h, cs = hierarchy([water(1, 1, 0.3, 5, 5, 0.5), water(2, 2, 19.7, 5, 5, 0.5)],
                    cryst1="")
  assert water_altlocs.run(h, None).n_labelled == 0

ser = """\
ATOM      1  N   SER A   5       5.000   5.000   5.000  1.00 10.00           N
ATOM      2  CA  SER A   5       6.450   5.000   5.000  1.00 10.00           C
ATOM      3  C   SER A   5       7.000   6.400   5.000  1.00 10.00           C
ATOM      4  O   SER A   5       6.300   7.400   5.000  1.00 10.00           O
ATOM      5  CB ASER A   5       7.000   4.200   6.200  0.60 10.00           C
ATOM      6  OG ASER A   5       6.500   2.900   6.300  0.60 10.00           O
ATOM      7  CB BSER A   5       7.000   4.200   3.800  0.40 10.00           C
ATOM      8  OG BSER A   5       6.500   2.900   3.700  0.40 10.00           O"""

def test_water_at_split_residue():
  # a water on the site of OG B: it belongs to conformer A
  h, cs = hierarchy(ser.splitlines() + [water(9, 10, 6.6, 2.5, 3.6, 0.6)])
  r = water_altlocs.run(h, cs)
  assert r.labelled == [("HOH A 10", "A")], r.labelled
  assert altlocs(h)[("10", "O")] == "A"
  # three conformers at the atom: left to the user
  three = ser.splitlines() + [
    "ATOM      9  CB CSER A   5       7.000   4.200   5.000  0.20 10.00           C",
    "ATOM     10  OG CSER A   5       6.500   2.900   5.000  0.20 10.00           O",
    water(11, 10, 6.6, 2.5, 3.6, 0.6)]
  h, cs = hierarchy(three)
  r = water_altlocs.run(h, cs)
  assert r.n_labelled == 0, r.labelled
  # one line per overlapping atom (CB B and OG B are both within 2.2 A)
  assert len(r.manual) == 2, r.manual
  assert "HOH A 10 overlaps OG SER A 5 (residue altlocs ABC)" in r.manual[1], \
    r.manual

def test_water_and_ligand():
  # a chloride at 0.6 and a water at 0.4 on overlapping sites: the ligand
  # takes A, the water B; a water on a blank partial protein atom is
  # reported, not labelled
  lines = ser.splitlines() + [
    "HETATM   20 CL    CL A  20      12.000  12.000  12.000  0.60 10.00          CL",
    water(21, 21, 12.9, 12.0, 12.0, 0.4),
    "ATOM     22  N   GLY A   6      15.000  15.000  15.000  0.50 10.00           N",
    water(23, 22, 15.7, 15.0, 15.0, 0.5)]
  h, cs = hierarchy(lines)
  r = water_altlocs.run(h, cs)
  assert r.labelled == [("CL A 20", "A"), ("HOH A 21", "B")], r.labelled
  assert len(r.manual) == 1 and "blank partial N GLY A 6" in r.manual[0], \
    r.manual
  assert altlocs(h)[("22", "O")] == ""

if __name__ == "__main__":
  test_water_pair()
  test_three_waters()
  test_full_occupancy_untouched()
  test_symmetry_mate()
  test_water_at_split_residue()
  test_water_and_ligand()
  print("OK")
