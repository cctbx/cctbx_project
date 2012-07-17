from iotbx import pdb
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
from libtbx.str_utils import show_string
from libtbx.utils import hashlib_md5, Sorry, format_cpu_times
import libtbx.load_env
from libtbx import Auto
from cStringIO import StringIO
import math
import cPickle
import pickle
import random
import sys, os

def exercise_atom():
  a = pdb.hierarchy.atom()
  assert a.name == ""
  a.name = "abcd"
  assert a.name == "abcd"
  try: a.name = "xyzhkl"
  except (ValueError, RuntimeError), e:
    assert str(e) == "string is too long for target variable " \
      "(maximum length is 4 characters, 6 given)."
  else: raise Exception_expected
  assert a.segid == ""
  a.segid = "stuv"
  assert a.segid == "stuv"
  assert a.element == ""
  a.element = "ca"
  assert a.element == "ca"
  assert a.charge == ""
  a.charge = "2+"
  assert a.charge == "2+"
  assert a.serial == ""
  a.serial = "A0000"
  assert a.serial == "A0000"
  assert a.serial_as_int() == 100000
  assert a.xyz == (0,0,0)
  a.xyz = (1,-2,3)
  assert a.xyz == (1,-2,3)
  assert a.sigxyz == (0,0,0)
  a.sigxyz = (-2,3,1)
  assert a.sigxyz == (-2,3,1)
  assert a.occ == 0
  a.occ = 0.5
  assert a.occ == 0.5
  assert a.sigocc == 0
  a.sigocc = 0.7
  assert a.sigocc == 0.7
  assert a.b == 0
  a.b = 5
  assert a.b == 5
  assert a.sigb == 0
  a.sigb = 7
  assert a.sigb == 7
  assert a.uij == (-1,-1,-1,-1,-1,-1)
  assert not a.uij_is_defined()
  a.uij = (1,-2,3,4,-5,6)
  assert a.uij == (1,-2,3,4,-5,6)
  assert a.uij_is_defined()
  a.uij_erase()
  assert not a.uij_is_defined()
  assert not a.siguij_is_defined()
  if (pdb.hierarchy.atom.has_siguij()):
    assert a.siguij == (-1,-1,-1,-1,-1,-1)
    a.siguij = (-2,3,4,-5,6,1)
    assert a.siguij == (-2,3,4,-5,6,1)
    assert a.siguij_is_defined()
  a.siguij_erase()
  assert not a.siguij_is_defined()
  assert not a.hetero
  a.hetero = True
  assert a.hetero
  assert a.i_seq == 0
  try: a.i_seq = 1
  except AttributeError: pass
  else: raise Exception_expected
  assert a.tmp == 0
  a.tmp = 3
  assert a.tmp == 3
  #
  a = (pdb.hierarchy.atom()
    .set_name(new_name="NaMe")
    .set_segid(new_segid="sEgI")
    .set_element(new_element="El")
    .set_charge(new_charge="cH")
    .set_serial(new_serial="B1234")
    .set_xyz(new_xyz=(1.3,2.1,3.2))
    .set_sigxyz(new_sigxyz=(.1,.2,.3))
    .set_occ(new_occ=0.4)
    .set_sigocc(new_sigocc=0.1)
    .set_b(new_b=4.8)
    .set_sigb(new_sigb=0.7)
    .set_uij(new_uij=(1.3,2.1,3.2,4.3,2.7,9.3))
    .set_hetero(new_hetero=True))
  if (pdb.hierarchy.atom.has_siguij()):
    assert a.set_siguij(new_siguij=(.1,.2,.3,.6,.1,.9)) is a
  assert a.name == "NaMe"
  assert a.segid == "sEgI"
  assert a.element == "El"
  assert a.charge == "cH"
  assert a.serial == "B1234"
  assert approx_equal(a.xyz, (1.3,2.1,3.2))
  assert approx_equal(a.sigxyz, (.1,.2,.3))
  assert approx_equal(a.occ, 0.4)
  assert approx_equal(a.sigocc, 0.1)
  assert approx_equal(a.b, 4.8)
  assert approx_equal(a.sigb, 0.7)
  assert approx_equal(a.uij, (1.3,2.1,3.2,4.3,2.7,9.3))
  if (pdb.hierarchy.atom.has_siguij()):
    assert approx_equal(a.siguij, (.1,.2,.3,.6,.1,.9))
  assert a.hetero
  assert a.tmp == 0
  try: a.set_name(new_name="12345")
  except (ValueError, RuntimeError), e:
    assert str(e) == "string is too long for target variable " \
      "(maximum length is 4 characters, 5 given)."
  else: raise Exception_expected
  #
  a.tmp = 7
  ac = a.detached_copy()
  assert ac.tmp == 0
  assert ac.name == "1234"
  assert ac.segid == "sEgI"
  assert ac.element == "El"
  assert ac.charge == "cH"
  assert ac.serial == "B1234"
  assert approx_equal(ac.xyz, (1.3,2.1,3.2))
  assert approx_equal(ac.sigxyz, (.1,.2,.3))
  assert approx_equal(ac.occ, 0.4)
  assert approx_equal(ac.sigocc, 0.1)
  assert approx_equal(ac.b, 4.8)
  assert approx_equal(ac.sigb, 0.7)
  assert approx_equal(ac.uij, (1.3,2.1,3.2,4.3,2.7,9.3))
  if (pdb.hierarchy.atom.has_siguij()):
    assert approx_equal(ac.siguij, (.1,.2,.3,.6,.1,.9))
  assert ac.hetero
  #
  for e in ["H", "H ", " H", "D", "D ", " D"]:
    a.element = e
    assert a.element_is_hydrogen()
  for e in ["", "h", "h ", " h", "d", "d ", " d"]:
    a.element = e
    assert not a.element_is_hydrogen()
  #
  a.name = "1234"
  a.element = "El"
  assert a.determine_chemical_element_simple() is None
  a.name = "NA  "
  a.element = " N"
  assert a.determine_chemical_element_simple() == " N"
  a.element = "CU"
  assert a.determine_chemical_element_simple() == "CU"
  for e in ["", " ", "  "]:
    a.element = e
    a.name = "NA  "
    assert a.determine_chemical_element_simple() == "NA"
    a.name = " D"
    assert a.determine_chemical_element_simple() == " D"
    for d in "0123456789":
      a.name = d+"H"
      assert a.determine_chemical_element_simple() == " H"
  a.set_name(new_name=None)
  a.set_segid(new_segid=None)
  a.set_element(new_element=None)
  a.set_charge(new_charge=None)
  a.set_serial(new_serial=None)
  assert a.name == ""
  assert a.segid == ""
  assert a.element == ""
  assert a.charge == ""
  assert a.serial == ""
  #
  assert not a.set_chemical_element_simple_if_necessary()
  assert a.element == ""
  a.name = " N  "
  assert a.set_chemical_element_simple_if_necessary()
  assert a.element == " N"
  a.name = " X  "
  a.element = ""
  assert not a.set_chemical_element_simple_if_necessary()
  assert a.element == ""
  a.element = "Na"
  assert a.set_chemical_element_simple_if_necessary()
  assert a.element == "NA"
  a.element = "Na"
  assert not a.set_chemical_element_simple_if_necessary(tidy_existing=False)
  assert a.element == "Na"
  a.element = "N"
  assert a.set_chemical_element_simple_if_necessary()
  assert a.element == " N"
  a.element = "N"
  assert not a.set_chemical_element_simple_if_necessary(tidy_existing=False)
  assert a.element == "N"
  a.element = " N"
  assert not a.set_chemical_element_simple_if_necessary()
  assert a.element == " N"
  a.name = ""
  a.element = ""
  #
  assert a.charge == ""
  assert a.charge_tidy() == "  "
  assert a.charge_tidy(strip=True) == ""
  def check(inp, out):
    a.charge = inp
    assert a.charge_tidy() == out
  check("", "  ")
  check(" ", "  ")
  check("  ", "  ")
  check("0", "  ")
  check(" 0", "  ")
  check("0 ", "  ")
  check("00", "  ")
  check(" +", "1+")
  check("+ ", "1+")
  check("++", "2+")
  check(" -", "1-")
  check("- ", "1-")
  check("--", "2-")
  check("0+", "0+")
  check("+0", "0+")
  check("0-", "0-")
  check("-0", "0-")
  check("+-", None)
  check("-+", None)
  check("3+", "3+")
  check("+4", "4+")
  check("5-", "5-")
  check("-6", "6-")
  check("7 ", None)
  check(" 8", None)
  check("12", None)
  check("1a", None)
  check("a1", None)
  check("a ", None)
  check(" a", None)
  a.charge = ""
  #
  assert a.distance(other=a) == 0
  b = a.detached_copy()
  b.xyz = (3,5,2)
  assert approx_equal(a.distance(other=b), 3.56931365951)
  assert approx_equal(b.distance(other=a), 3.56931365951)
  assert approx_equal(a.distance(other_xyz=(3,5,2)), 3.56931365951)
  assert approx_equal(a.distance(other_xyz=(2,3,5)), 2.13072757527)
  assert a.angle(atom_1=a, atom_3=a) is None
  assert a.angle(atom_1=a, atom_3=a, deg=True) is None
  assert a.angle(atom_1_xyz=a.xyz, atom_3_xyz=a.xyz) is None
  assert a.angle(atom_1_xyz=a.xyz, atom_3_xyz=a.xyz, deg=True) is None
  assert approx_equal(a.angle(atom_1=b, atom_3=b), 0)
  assert approx_equal(a.angle(atom_1=b, atom_3=b, deg=True), 0)
  c = a.detached_copy()
  c.xyz = (5,3,1)
  assert approx_equal(a.angle(atom_1=b, atom_3=c, deg=True), 42.6776898341)
  assert approx_equal(math.degrees(a.angle(atom_1=c, atom_3=b, deg=False)),
    42.6776898341)
  #
  ag = pdb.hierarchy.atom_group()
  ac = pdb.hierarchy.atom(parent=ag, other=a)
  assert ac.memory_id() != a.memory_id()
  assert ac.parent().memory_id() == ag.memory_id()
  assert ac.name == a.name
  assert ac.segid == a.segid
  assert ac.element == a.element
  assert ac.charge == a.charge
  assert ac.serial == a.serial
  assert ac.xyz == a.xyz
  assert ac.sigxyz == a.sigxyz
  assert ac.occ == a.occ
  assert ac.sigocc == a.sigocc
  assert ac.b == a.b
  assert ac.sigb == a.sigb
  assert ac.uij == a.uij
  if (pdb.hierarchy.atom.has_siguij()):
    assert ac.siguij == a.siguij
  assert ac.hetero == a.hetero
  assert ac.tmp == 0
  #
  assert a.pdb_label_columns() == "               "
  #
  a = pdb.hierarchy.atom()
  assert a.set_serial(new_serial="ABCDE") is a
  assert a.serial == "ABCDE"
  a.serial = "CDEFG"
  assert a.serial == "CDEFG"
  assert a.set_serial(new_serial=200000) is a
  assert a.serial == "A255S"
  assert a.serial_as_int() == 200000
  a.serial = 100000
  assert a.serial == "A0000"
  try: a.set_serial(new_serial="ABCDEF")
  except (ValueError, RuntimeError), e:
    assert str(e) == "string is too long for target variable " \
      "(maximum length is 5 characters, 6 given)."
  else: raise Exception_expected
  try: a.set_serial(new_serial=-10000)
  except ValueError, e:
    assert str(e) == "value is less than -9999"
  else: raise Exception_expected
  try: a.set_serial(new_serial=87440032)
  except ValueError, e:
    assert str(e) == "value is greater than 87440031"
  else: raise Exception_expected
  try: a.set_serial(new_serial=sys)
  except TypeError, e:
    assert str(e) == "value must be a Python str or int."
  else: raise Exception_expected
  try: a.serial = sys
  except TypeError, e:
    assert str(e) == "value must be a Python str or int."
  else: raise Exception_expected
  #
  atoms = pdb.hierarchy.af_shared_atom()
  atoms.reset_serial()
  atoms.reset_i_seq()
  atoms.reset_tmp()
  atoms.append(pdb.hierarchy.atom())
  assert [atom.serial for atom in atoms] == [""]
  assert [atom.i_seq for atom in atoms] == [0]
  assert [atom.tmp for atom in atoms] == [0]
  atoms.reset_serial(first_value=2)
  atoms.reset_i_seq()
  sentinel = atoms.reset_tmp(first_value=2)
  assert [atom.serial for atom in atoms] == ["    2"]
  assert [atom.i_seq for atom in atoms] == [0]
  assert [atom.tmp for atom in atoms] == [2]
  atoms.append(pdb.hierarchy.atom())
  atoms.append(pdb.hierarchy.atom())
  assert [atom.serial for atom in atoms] == ["    2", "", ""]
  assert [atom.i_seq for atom in atoms] == [0,0,0]
  assert [atom.tmp for atom in atoms] == [2,0,0]
  del sentinel
  assert [atom.tmp for atom in atoms] == [0,0,0]
  atoms.reset_serial()
  atoms.reset_i_seq()
  sentinel = atoms.reset_tmp()
  assert [atom.serial for atom in atoms] == ["    1", "    2", "    3"]
  assert [atom.i_seq for atom in atoms] == [0,1,2]
  assert [atom.tmp for atom in atoms] == [0,1,2]
  del sentinel
  assert [atom.tmp for atom in atoms] == [0,0,0]
  sentinel = atoms.reset_tmp(first_value=0, increment=0)
  assert [atom.tmp for atom in atoms] == [0] * 3
  del sentinel
  sentinel = atoms.reset_tmp(first_value=5, increment=-3)
  assert [atom.tmp for atom in atoms] == [5,2,-1]
  del sentinel
  #
  sentinel = atoms.reset_tmp_for_occupancy_groups_simple()
  assert [atom.tmp for atom in atoms] == [0,1,2]
  del sentinel
  atoms[0].element = "D"
  atoms[2].element = "H"
  sentinel = atoms.reset_tmp_for_occupancy_groups_simple()
  assert [atom.tmp for atom in atoms] == [-1,1,-1]
  del sentinel
  assert [atom.tmp for atom in atoms] == [0,0,0]
  #
  sentinel = atoms.reset_tmp()
  try: atoms.reset_tmp()
  except RuntimeError, e:
    assert str(e) == \
      "Another associated atom_tmp_sentinel instance still exists."
  else: raise Exception_expected
  #
  for atom,n in zip(atoms, ["", "NA  ", " X  "]): atom.name = n
  for atom,e in zip(atoms, ["", "", "N"]): atom.element = e
  assert atoms.set_chemical_element_simple_if_necessary() == 2
  assert [atom.element for atom in atoms] == ["", "NA", " N"]
  for atom,e in zip(atoms, ["", "", "N"]): atom.element = e
  assert atoms.set_chemical_element_simple_if_necessary(
    tidy_existing=False) == 1
  assert [atom.element for atom in atoms] == ["", "NA", "N"]
  #
  assert sorted(atoms.build_dict().keys()) == ["", " X  ", "NA  "]
  assert sorted(atoms.build_dict(strip_names=True).keys()) == ["", "NA", "X"]
  atoms[0].name = "x   "
  assert sorted(atoms.build_dict(upper_names=True).keys()) \
      == [" X  ", "NA  ", "X   "]
  d = atoms.build_dict(
    strip_names=True,
    upper_names=True,
    throw_runtime_error_if_duplicate_keys=False)
  assert sorted([atom.name for atom in d.values()]) == ["NA  ", "x   "]
  try:
    atoms.build_dict(strip_names=True, upper_names=True)
  except RuntimeError, e:
    assert not show_diff(str(e), '''\
Duplicate keys in build_dict(strip_names=true, upper_names=true,\
 convert_stars_to_primes=false):
  pdb="x              "
  pdb=" X             "''')
  else: raise Exception_expected
  atoms[0].name = "x*  "
  d = atoms.build_dict(
    strip_names=True,
    upper_names=True,
    convert_stars_to_primes=True)
  assert sorted(d.keys()) == ["NA", "X", "X'"]
  #
  a = pdb.hierarchy.atom()
  def check(scattering_type, e, c):
    a.set_element_and_charge_from_scattering_type_if_necessary(
      scattering_type=scattering_type)
    assert not show_diff(a.element, e)
    assert not show_diff(a.charge, c)
  a.element = ""
  a.charge = ""
  check("", "  ", "  ")
  a.element = "Na"
  a.charge = ""
  check("", "Na", "")
  a.charge = "  "
  check("", "Na", "  ")
  a.element = ""
  a.charge = ""
  check("Na+", "NA", "1+")
  a.element = "Na"
  a.charge = "+"
  check("NA1+", "Na", "+")
  a.element = "Na"
  a.charge = "2+"
  check("NA1+", "NA", "1+")

def exercise_atom_group():
  ag = pdb.hierarchy.atom_group()
  assert ag.altloc == ""
  assert ag.resname == ""
  ag = pdb.hierarchy.atom_group(altloc=None, resname=None)
  assert ag.altloc == ""
  assert ag.resname == ""
  ag = pdb.hierarchy.atom_group(altloc="a", resname="xyz")
  assert ag.altloc == "a"
  assert ag.resname == "xyz"
  ag.altloc = None
  ag.resname = None
  assert ag.altloc == ""
  assert ag.resname == ""
  assert ag.confid() == "    "
  #
  ag.altloc = "l"
  ag.resname = "res"
  assert ag.confid() == "lres"
  ag.append_atom(atom=pdb.hierarchy.atom().set_name(new_name="n"))
  rg = pdb.hierarchy.residue_group()
  for i,agc in enumerate([
                 pdb.hierarchy.atom_group(parent=rg, other=ag),
                 ag.detached_copy()]):
    assert agc.memory_id() != ag.memory_id()
    assert ag.parent() is None
    if (i == 0):
      assert agc.parent().memory_id() == rg.memory_id()
    else:
      assert agc.parent() is None
    assert agc.altloc == "l"
    assert agc.resname == "res"
    assert agc.atoms_size() == 1
    assert agc.atoms()[0].memory_id() != ag.atoms()[0].memory_id()
    assert agc.atoms()[0].name == "n"
    ag.append_atom(atom=pdb.hierarchy.atom().set_name(new_name="o"))
    assert ag.atoms_size() == 2+i
    assert agc.atoms_size() == 1
  #
  ag = pdb.hierarchy.atom_group()
  assert ag.parent() is None
  rg1 = pdb.hierarchy.residue_group()
  rg2 = pdb.hierarchy.residue_group()
  assert rg1.memory_id() != rg2.memory_id()
  ag = pdb.hierarchy.atom_group(parent=rg1)
  assert ag.parent(optional=False).memory_id() == rg1.memory_id()
  del rg1
  assert ag.parent() is None
  try:
    ag.parent(optional=False)
  except RuntimeError, e:
    assert not show_diff(str(e), "atom_group has no parent residue_group")
  else: raise Exception_expected
  #
  rg1 = pdb.hierarchy.residue_group()
  ag = pdb.hierarchy.atom_group(altloc="a", resname="xyz")
  rg1.append_atom_group(atom_group=ag)
  assert ag.altloc == "a"
  assert ag.resname == "xyz"
  assert ag.parent().memory_id() == rg1.memory_id()
  del rg1
  assert ag.parent() is None
  #
  ag.pre_allocate_atoms(number_of_additional_atoms=2)
  assert ag.atoms_size() == 0
  assert ag.atoms().size() == 0
  ag.append_atom(atom=pdb.hierarchy.atom().set_name(new_name="ca"))
  assert ag.atoms_size() == 1
  assert ag.atoms().size() == 1
  ag.append_atom(atom=pdb.hierarchy.atom().set_name(new_name="n"))
  assert ag.atoms_size() == 2
  assert ag.atoms().size() == 2
  assert [atom.name for atom in ag.atoms()] == ["ca", "n"]
  for i in xrange(3):
    ag.append_atom(pdb.hierarchy.atom())
  assert ag.atoms_size() == 5
  assert ag.atoms().size() == 5
  for atom in ag.atoms():
    assert atom.parent(optional=False).memory_id() == ag.memory_id()
  assert [a.name for a in ag.atoms()] == ["ca", "n", "", "", ""]
  #
  ag.insert_atom(i=0, atom=pdb.hierarchy.atom().set_name(new_name="0"))
  assert [a.name for a in ag.atoms()] == ["0", "ca", "n", "", "", ""]
  ag.insert_atom(i=-1, atom=pdb.hierarchy.atom().set_name(new_name="x"))
  assert [a.name for a in ag.atoms()] == ["0", "ca", "n", "", "", "x", ""]
  a = ag.atoms()[-1]
  assert a.parent().memory_id() == ag.memory_id()
  ag.remove_atom(i=-1)
  assert a.parent() is None
  assert [a.name for a in ag.atoms()] == ["0", "ca", "n", "", "", "x"]
  ag.remove_atom(i=1)
  assert [a.name for a in ag.atoms()] == ["0", "n", "", "", "x"]
  a = ag.atoms()[-2]
  assert a.parent().memory_id() == ag.memory_id()
  assert ag.find_atom_index(atom=a, must_be_present=True) == 3
  ag.remove_atom(i=-2)
  assert a.parent() is None
  assert [a.name for a in ag.atoms()] == ["0", "n", "", "x"]
  a = pdb.hierarchy.atom().set_name(new_name="y")
  assert ag.find_atom_index(atom=a) == -1
  try: ag.find_atom_index(atom=a, must_be_present=True)
  except RuntimeError, e:
    assert str(e) == "atom not in atom_group."
  else: raise Exception_expected
  ag.insert_atom(i=4, atom=a)
  assert ag.find_atom_index(atom=a) == 4
  assert [a.name for a in ag.atoms()] == ["0", "n", "", "x", "y"]
  #
  try: pdb.hierarchy.atom_group(altloc="ab")
  except (ValueError, RuntimeError), e:
    assert str(e) == "string is too long for target variable " \
      "(maximum length is 1 character, 2 given)."
  else: raise Exception_expected
  #
  ag1 = pdb.hierarchy.atom_group()
  atom = pdb.hierarchy.atom()
  assert atom.parent() is None
  ag1.append_atom(atom=atom)
  assert atom.parent().memory_id() == ag1.memory_id()
  ag2 = pdb.hierarchy.atom_group()
  ag2.append_atom_with_other_parent(atom=atom)
  assert atom.parent().memory_id() == ag1.memory_id()
  del ag1
  assert atom.parent() is None
  try:
    atom.parent(optional=False)
  except RuntimeError, e:
    assert not show_diff(str(e), "atom has no parent atom_group")
  else: raise Exception_expected

def exercise_residue_group():
  rg = pdb.hierarchy.residue_group()
  assert rg.resseq == ""
  assert rg.icode == ""
  assert rg.link_to_previous
  rg.resseq = 10000
  assert rg.resseq == "A000"
  assert rg.resseq_as_int() == 10000
  rg = pdb.hierarchy.residue_group(
    resseq="   1", icode="i", link_to_previous=False)
  assert rg.resseq == "   1"
  rg.resseq = "   2"
  assert rg.resseq == "   2"
  assert rg.icode == "i"
  rg.icode = "j"
  assert rg.icode == "j"
  assert not rg.link_to_previous
  rg.link_to_previous = True
  assert rg.link_to_previous
  rg.link_to_previous = False
  #
  ag = pdb.hierarchy.atom_group(altloc="a")
  assert ag.parent() is None
  rg.append_atom_group(atom_group=ag)
  assert ag.parent().memory_id() == rg.memory_id()
  c = pdb.hierarchy.chain()
  for i,rgc in enumerate([
                 pdb.hierarchy.residue_group(parent=c, other=rg),
                 rg.detached_copy()]):
    assert rgc.memory_id() != rg.memory_id()
    assert rg.parent() is None
    if (i == 0):
      assert rgc.parent().memory_id() == c.memory_id()
    else:
      assert rgc.parent() is None
    assert rgc.resseq == "   2"
    assert rgc.icode == "j"
    assert not rgc.link_to_previous
    assert rgc.atom_groups_size() == 1
    assert rgc.atom_groups()[0].memory_id() != rg.atom_groups()[0].memory_id()
    assert rgc.atom_groups()[0].altloc == "a"
    rg.append_atom_group(atom_group=pdb.hierarchy.atom_group(altloc="%d"%i))
    assert rg.atom_groups_size() == 2+i
    assert rgc.atom_groups_size() == 1
    assert [ag.altloc for ag in rg.atom_groups()] == ["a", "0", "1"][:i+2]
  #
  c1 = pdb.hierarchy.chain(id="a")
  c2 = pdb.hierarchy.chain(id="b")
  assert c1.memory_id() != c2.memory_id()
  rg = pdb.hierarchy.residue_group()
  assert rg.parent() is None
  rg = pdb.hierarchy.residue_group(parent=c1)
  assert rg.parent(optional=False).memory_id() == c1.memory_id()
  del c1
  assert rg.parent() is None
  try:
    rg.parent(optional=False)
  except RuntimeError, e:
    assert not show_diff(str(e), "residue_group has no parent chain")
  else: raise Exception_expected
  #
  c1 = pdb.hierarchy.chain(id="p")
  rg13l = pdb.hierarchy.residue_group(resseq="13", icode="l")
  c1.append_residue_group(rg13l)
  assert rg13l.resseq == "13"
  assert rg13l.icode == "l"
  #
  c1 = pdb.hierarchy.chain(id="a")
  c1.pre_allocate_residue_groups(number_of_additional_residue_groups=2)
  assert c1.residue_groups_size() == 0
  assert len(c1.residue_groups()) == 0
  for i in xrange(2):
    c1.append_residue_group(residue_group=pdb.hierarchy.residue_group())
  assert c1.residue_groups_size() == 2
  assert len(c1.residue_groups()) == 2
  for residue_group in c1.residue_groups():
    assert residue_group.parent().memory_id() == c1.memory_id()
  assert c1.atoms_size() == 0
  assert c1.atoms().size() == 0
  #
  for altloc in ["w", "v", "u"]:
    rg.insert_atom_group(
      i=0, atom_group=pdb.hierarchy.atom_group(altloc=altloc))
  assert [ag.altloc for ag in rg.atom_groups()] == ["u", "v", "w"]
  rg.remove_atom_group(i=-1)
  assert [ag.altloc for ag in rg.atom_groups()] == ["u", "v"]
  ag = rg.atom_groups()[1]
  assert ag.parent().memory_id() == rg.memory_id()
  assert rg.find_atom_group_index(atom_group=ag) == 1
  rg.remove_atom_group(atom_group=ag)
  assert ag.parent() is None
  assert rg.find_atom_group_index(atom_group=ag) == -1
  try: rg.find_atom_group_index(atom_group=ag, must_be_present=True)
  except RuntimeError, e:
    assert str(e) == "atom_group not in residue_group."
  else: raise Exception_expected
  #
  ag1 = pdb.hierarchy.atom_group()
  ag2 = pdb.hierarchy.atom_group()
  a = pdb.hierarchy.atom()
  ag1.append_atom(atom=a)
  try: ag2.append_atom(atom=a)
  except RuntimeError, e:
    assert str(e) == "atom has another parent atom_group already."
  else: raise Exception_expected
  #
  rg = pdb.hierarchy.residue_group()
  assert rg.resid() == "     "
  rg = pdb.hierarchy.residue_group(resseq="1", icode="i")
  assert rg.resid() == "   1i"
  rg = pdb.hierarchy.residue_group(resseq=" 1 ", icode="j")
  assert rg.resid() == "  1 j"
  rg = pdb.hierarchy.residue_group(resseq="ABCD", icode="")
  assert rg.resid() == "ABCD "
  rg = pdb.hierarchy.residue_group(resseq="ABCD", icode="E")
  assert rg.resid() == "ABCDE"
  #
  rg = pdb.hierarchy.residue_group()
  ag = pdb.hierarchy.atom_group(altloc=" ")
  rg.append_atom_group(atom_group=ag)
  assert not rg.have_conformers()
  ag = pdb.hierarchy.atom_group(altloc="")
  rg.append_atom_group(atom_group=ag)
  assert not rg.have_conformers()
  ag = pdb.hierarchy.atom_group(altloc="a")
  rg.append_atom_group(atom_group=ag)
  assert rg.have_conformers()
  #
  rg = pdb.hierarchy.residue_group()
  assert rg.move_blank_altloc_atom_groups_to_front() == 0
  ag = pdb.hierarchy.atom_group(altloc="a")
  rg.append_atom_group(atom_group=ag)
  assert rg.move_blank_altloc_atom_groups_to_front() == 0
  ag = pdb.hierarchy.atom_group(altloc=" ")
  rg.append_atom_group(atom_group=ag)
  assert rg.move_blank_altloc_atom_groups_to_front() == 1
  #
  rg = pdb.hierarchy.residue_group()
  rg.resseq = "x"
  try: rg.resseq_as_int()
  except (ValueError, RuntimeError), e:
    assert not show_diff(str(e), 'invalid residue sequence number: "x"')
  else: raise Exception_expected
  #
  rg = pdb.hierarchy.residue_group()
  assert len(rg.unique_resnames()) == 0
  def rga(altloc, resname):
    rg.append_atom_group(atom_group=pdb.hierarchy.atom_group(
      altloc=altloc, resname=resname))
  rga("", "RN1")
  assert list(rg.unique_resnames()) == ["RN1"]
  rga("", "RN")
  assert list(rg.unique_resnames()) == ["RN", "RN1"]
  rga("A", "RN1")
  assert list(rg.unique_resnames()) == ["RN", "RN1"]
  rga("A", "RN2")
  assert list(rg.unique_resnames()) == ["RN", "RN1", "RN2"]

def exercise_chain():
  c = pdb.hierarchy.chain()
  assert c.id == ""
  c = pdb.hierarchy.chain(id="a")
  assert c.id == "a"
  c.id = "x"
  assert c.id == "x"
  #
  m1 = pdb.hierarchy.model(id="1")
  m2 = pdb.hierarchy.model(id="2")
  assert m1.memory_id() != m2.memory_id()
  c = pdb.hierarchy.chain()
  assert c.parent() is None
  c = pdb.hierarchy.chain(parent=m1)
  assert c.parent(optional=False).memory_id() == m1.memory_id()
  del m1
  assert c.parent() is None
  try:
    c.parent(optional=False)
  except RuntimeError, e:
    assert not show_diff(str(e), "chain has no parent model")
  else: raise Exception_expected
  #
  c = pdb.hierarchy.chain()
  #
  c = pdb.hierarchy.chain()
  c.pre_allocate_residue_groups(number_of_additional_residue_groups=2)
  assert c.residue_groups_size() == 0
  assert len(c.residue_groups()) == 0
  for i in xrange(2):
    c.append_residue_group(residue_group=pdb.hierarchy.residue_group())
  assert c.residue_groups_size() == 2
  assert len(c.residue_groups()) == 2
  for residue_group in c.residue_groups():
    assert residue_group.parent().memory_id() == c.memory_id()
  assert c.atoms_size() == 0
  assert c.atoms().size() == 0
  #
  c.residue_groups()[0].resseq = "ugh"
  c.id = "ci"
  m = pdb.hierarchy.model()
  for i,cc in enumerate([
                pdb.hierarchy.chain(parent=m, other=c),
                c.detached_copy()]):
    assert cc.memory_id() != c.memory_id()
    assert c.parent() is None
    if (i == 0):
      assert cc.parent().memory_id() == m.memory_id()
    else:
      assert cc.parent() is None
    assert cc.id == "ci"
    assert cc.residue_groups_size() == 2
    assert cc.residue_groups()[0].memory_id() \
         != c.residue_groups()[0].memory_id()
    assert cc.residue_groups()[0].resseq == "ugh"
    c.append_residue_group(
      residue_group=pdb.hierarchy.residue_group(resseq="%03d"%i))
    assert c.residue_groups_size() == 3+i
    assert cc.residue_groups_size() == 2
    assert [rg.resseq for rg in c.residue_groups()] \
        == ["ugh", "", "000", "001"][:i+3]
  #
  c.insert_residue_group(
    i=3, residue_group=pdb.hierarchy.residue_group(resseq="b012"))
  assert [rg.resseq for rg in c.residue_groups()] \
      == ["ugh", "", "000", "b012", "001"]
  c.remove_residue_group(i=1)
  assert [rg.resseq for rg in c.residue_groups()] \
      == ["ugh", "000", "b012", "001"]
  rg = c.residue_groups()[1]
  assert rg.parent().memory_id() == c.memory_id()
  assert c.find_residue_group_index(residue_group=rg) == 1
  c.remove_residue_group(residue_group=rg)
  assert rg.parent() is None
  assert c.find_residue_group_index(residue_group=rg) == -1
  try: c.find_residue_group_index(residue_group=rg, must_be_present=True)
  except RuntimeError, e:
    assert str(e) == "residue_group not in chain."
  else: raise Exception_expected
  #
  rg1 = pdb.hierarchy.residue_group()
  rg2 = pdb.hierarchy.residue_group()
  ag = pdb.hierarchy.atom_group()
  rg1.append_atom_group(atom_group=ag)
  try: rg2.append_atom_group(atom_group=ag)
  except RuntimeError, e:
    assert str(e) == "atom_group has another parent residue_group already."
  else: raise Exception_expected
  #
  r = pdb.hierarchy.root()
  m = pdb.hierarchy.model()
  r.append_model(m)
  c = pdb.hierarchy.chain(id="c")
  m.append_chain(c)
  assert r.as_pdb_string() == "TER\n"
  rg = pdb.hierarchy.residue_group(resseq="s", icode="j")
  c.append_residue_group(residue_group=rg)
  ag = pdb.hierarchy.atom_group(altloc="a", resname="r")
  rg.append_atom_group(atom_group=ag)
  ag.append_atom(pdb.hierarchy.atom().set_name("n"))
  assert ag.only_atom().pdb_label_columns() == "n   a  r c   sj"
  assert not show_diff(r.as_pdb_string(), """\
ATOM        n   a  r c   sj      0.000   0.000   0.000  0.00  0.00
TER
""")
  rg = pdb.hierarchy.residue_group(resseq="t", icode="k")
  c.append_residue_group(residue_group=rg)
  ag =  pdb.hierarchy.atom_group(altloc="b", resname="q")
  rg.append_atom_group(atom_group=ag)
  ag.append_atom(pdb.hierarchy.atom().set_name("m"))
  rg = pdb.hierarchy.residue_group(
    resseq="u", icode="l", link_to_previous=False)
  c.append_residue_group(residue_group=rg)
  ag = pdb.hierarchy.atom_group(altloc="d", resname="p")
  rg.append_atom_group(atom_group=ag)
  ag.append_atom(pdb.hierarchy.atom().set_name("o"))
  assert not show_diff(r.as_pdb_string(), """\
ATOM        n   a  r c   sj      0.000   0.000   0.000  0.00  0.00
ATOM        m   b  q c   tk      0.000   0.000   0.000  0.00  0.00
BREAK
ATOM        o   d  p c   ul      0.000   0.000   0.000  0.00  0.00
TER
""")
  #
  atoms = c.atoms()
  assert atoms.size() == 3
  atoms[0].set_sigxyz((1,2,3)).set_sigocc(4).set_sigb(5)
  atoms[1].set_uij((6,7,8,3,5,4))
  siguij_2_line = ""
  if (pdb.hierarchy.atom.has_siguij()):
    atoms[2].set_siguij((.6,.7,.8,.3,.5,.4))
    siguij_2_line = """\
SIGUIJ      o   d  p c   ul    6000   7000   8000   3000   5000   4000
"""
  assert not show_diff(r.as_pdb_string(), """\
ATOM        n   a  r c   sj      0.000   0.000   0.000  0.00  0.00
SIGATM      n   a  r c   sj      1.000   2.000   3.000  4.00  5.00
ATOM        m   b  q c   tk      0.000   0.000   0.000  0.00  0.00
ANISOU      m   b  q c   tk   60000  70000  80000  30000  50000  40000
BREAK
ATOM        o   d  p c   ul      0.000   0.000   0.000  0.00  0.00
%sTER
""" % siguij_2_line)
  atoms[0].set_uij((6,3,8,2,9,1))
  siguij_0_line = ""
  if (pdb.hierarchy.atom.has_siguij()):
    atoms[0].set_siguij((.6,.3,.8,.2,.9,.1))
    siguij_0_line = """\
SIGUIJ      n   a  r c   sj    6000   3000   8000   2000   9000   1000        Cg
"""
  atoms[0].set_charge("Cg")
  assert not show_diff(r.as_pdb_string(), """\
ATOM        n   a  r c   sj      0.000   0.000   0.000  0.00  0.00            Cg
SIGATM      n   a  r c   sj      1.000   2.000   3.000  4.00  5.00            Cg
ANISOU      n   a  r c   sj   60000  30000  80000  20000  90000  10000        Cg
%sATOM        m   b  q c   tk      0.000   0.000   0.000  0.00  0.00
ANISOU      m   b  q c   tk   60000  70000  80000  30000  50000  40000
BREAK
ATOM        o   d  p c   ul      0.000   0.000   0.000  0.00  0.00
%sTER
""" % (siguij_0_line, siguij_2_line))
  sio = StringIO()
  assert r.as_pdb_string(cstringio=sio, interleaved_conf=1) is sio
  assert r.as_pdb_string(cstringio=sio, interleaved_conf=2) is sio
  assert r.as_pdb_string(cstringio=sio, atom_hetatm=False) is sio
  assert r.as_pdb_string(cstringio=sio, sigatm=False) is sio
  assert r.as_pdb_string(cstringio=sio, anisou=False) is sio
  assert isinstance(
    r.as_pdb_string(cstringio=sio, siguij=False, return_cstringio=False),
    str)
  expected = """\
ATOM        n   a  r c   sj      0.000   0.000   0.000  0.00  0.00            Cg
SIGATM      n   a  r c   sj      1.000   2.000   3.000  4.00  5.00            Cg
ANISOU      n   a  r c   sj   60000  30000  80000  20000  90000  10000        Cg
SIGUIJ      n   a  r c   sj    6000   3000   8000   2000   9000   1000        Cg
ATOM        m   b  q c   tk      0.000   0.000   0.000  0.00  0.00
ANISOU      m   b  q c   tk   60000  70000  80000  30000  50000  40000
BREAK
ATOM        o   d  p c   ul      0.000   0.000   0.000  0.00  0.00
SIGUIJ      o   d  p c   ul    6000   7000   8000   3000   5000   4000
TER
ATOM        n   a  r c   sj      0.000   0.000   0.000  0.00  0.00            Cg
SIGATM      n   a  r c   sj      1.000   2.000   3.000  4.00  5.00            Cg
ANISOU      n   a  r c   sj   60000  30000  80000  20000  90000  10000        Cg
SIGUIJ      n   a  r c   sj    6000   3000   8000   2000   9000   1000        Cg
ATOM        m   b  q c   tk      0.000   0.000   0.000  0.00  0.00
ANISOU      m   b  q c   tk   60000  70000  80000  30000  50000  40000
BREAK
ATOM        o   d  p c   ul      0.000   0.000   0.000  0.00  0.00
SIGUIJ      o   d  p c   ul    6000   7000   8000   3000   5000   4000
TER
SIGATM      n   a  r c   sj      1.000   2.000   3.000  4.00  5.00            Cg
ANISOU      n   a  r c   sj   60000  30000  80000  20000  90000  10000        Cg
SIGUIJ      n   a  r c   sj    6000   3000   8000   2000   9000   1000        Cg
ANISOU      m   b  q c   tk   60000  70000  80000  30000  50000  40000
BREAK
SIGUIJ      o   d  p c   ul    6000   7000   8000   3000   5000   4000
TER
ATOM        n   a  r c   sj      0.000   0.000   0.000  0.00  0.00            Cg
ANISOU      n   a  r c   sj   60000  30000  80000  20000  90000  10000        Cg
SIGUIJ      n   a  r c   sj    6000   3000   8000   2000   9000   1000        Cg
ATOM        m   b  q c   tk      0.000   0.000   0.000  0.00  0.00
ANISOU      m   b  q c   tk   60000  70000  80000  30000  50000  40000
BREAK
ATOM        o   d  p c   ul      0.000   0.000   0.000  0.00  0.00
SIGUIJ      o   d  p c   ul    6000   7000   8000   3000   5000   4000
TER
ATOM        n   a  r c   sj      0.000   0.000   0.000  0.00  0.00            Cg
SIGATM      n   a  r c   sj      1.000   2.000   3.000  4.00  5.00            Cg
SIGUIJ      n   a  r c   sj    6000   3000   8000   2000   9000   1000        Cg
ATOM        m   b  q c   tk      0.000   0.000   0.000  0.00  0.00
BREAK
ATOM        o   d  p c   ul      0.000   0.000   0.000  0.00  0.00
SIGUIJ      o   d  p c   ul    6000   7000   8000   3000   5000   4000
TER
ATOM        n   a  r c   sj      0.000   0.000   0.000  0.00  0.00            Cg
SIGATM      n   a  r c   sj      1.000   2.000   3.000  4.00  5.00            Cg
ANISOU      n   a  r c   sj   60000  30000  80000  20000  90000  10000        Cg
ATOM        m   b  q c   tk      0.000   0.000   0.000  0.00  0.00
ANISOU      m   b  q c   tk   60000  70000  80000  30000  50000  40000
BREAK
ATOM        o   d  p c   ul      0.000   0.000   0.000  0.00  0.00
TER
"""
  if (pdb.hierarchy.atom.has_siguij()):
    assert not show_diff(sio.getvalue(), expected)
  else:
    def filter_expected():
      result = []
      for line in expected.splitlines():
        if (line.startswith("SIGUIJ")): continue
        result.append(line)
      return "\n".join(result)+"\n"
    assert not show_diff(sio.getvalue(), filter_expected())
  assert not show_diff(r.as_pdb_string(
    atom_hetatm=False, sigatm=False, anisou=False, siguij=False), """\
BREAK
TER
""")
  #
  a = pdb.hierarchy.atom()
  assert a.pdb_label_columns() == "               "
  a.set_name("n123")
  assert a.pdb_label_columns() == "n123           "
  ag = pdb.hierarchy.atom_group(altloc="a", resname="res")
  ag.append_atom(a)
  assert a.pdb_label_columns() == "n123ares       "
  rg = pdb.hierarchy.residue_group(resseq="a000", icode="i")
  rg.append_atom_group(ag)
  assert a.pdb_label_columns() == "n123ares  a000i"
  c = pdb.hierarchy.chain(id="ke")
  c.append_residue_group(rg)
  assert a.pdb_label_columns() == "n123areskea000i"
  #
  assert a.pdb_element_charge_columns() == "    "
  a.set_element("e")
  assert a.pdb_element_charge_columns() == " e  "
  a.set_charge("+")
  assert a.pdb_element_charge_columns() == " e+ "
  a.set_element("el")
  a.set_charge("2+")
  assert a.pdb_element_charge_columns() == "el2+"

def exercise_model():
  m = pdb.hierarchy.model()
  assert m.id == ""
  m = pdb.hierarchy.model(id="42")
  assert m.id == "42"
  m.id = "-23"
  assert m.id == "-23"
  #
  m = pdb.hierarchy.model(id="17")
  assert m.parent() is None
  m.pre_allocate_chains(number_of_additional_chains=2)
  assert m.chains_size() == 0
  assert len(m.chains()) == 0
  ch_a = pdb.hierarchy.chain(id="a")
  m.append_chain(chain=ch_a)
  assert ch_a.id == "a"
  assert ch_a.parent().memory_id() == m.memory_id()
  assert m.chains_size() == 1
  assert len(m.chains()) == 1
  ch_b = pdb.hierarchy.chain(id="b")
  assert ch_b.id == "b"
  assert ch_b.parent() is None
  m.append_chain(chain=ch_b)
  assert m.chains_size() == 2
  chains = m.chains()
  assert len(chains) == 2
  assert chains[0].memory_id() == ch_a.memory_id()
  assert chains[1].memory_id() == ch_b.memory_id()
  for i in xrange(3):
    m.append_chain(pdb.hierarchy.chain())
  assert m.chains_size() == 5
  assert len(m.chains()) == 5
  for chain in m.chains():
    assert chain.parent().memory_id() == m.memory_id()
  assert m.atoms_size() == 0
  assert m.atoms().size() == 0
  #
  r = pdb.hierarchy.root()
  for i,mc in enumerate([
                pdb.hierarchy.model(parent=r, other=m),
                m.detached_copy()]):
    assert mc.memory_id() != m.memory_id()
    assert m.parent() is None
    if (i == 0):
      assert mc.parent().memory_id() == r.memory_id()
    else:
      assert mc.parent() is None
    assert mc.id == "17"
    assert mc.chains_size() == 5
    assert mc.chains()[0].memory_id() != m.chains()[0].memory_id()
    assert mc.chains()[0].id == "a"
    m.append_chain(chain=pdb.hierarchy.chain(id="%d"%i))
    assert m.chains_size() == 6+i
    assert mc.chains_size() == 5
    assert [c.id for c in m.chains()] \
        == ["a", "b", "", "", "", "0", "1"][:i+6]
  #
  m.insert_chain(i=-3, chain=pdb.hierarchy.chain(id="3"))
  assert [c.id for c in m.chains()] \
      == ["a", "b", "", "", "3", "", "0", "1"]
  m.remove_chain(i=-2)
  assert [c.id for c in m.chains()] \
      == ["a", "b", "", "", "3", "", "1"]
  c = m.chains()[0]
  assert c.parent().memory_id() == m.memory_id()
  assert m.find_chain_index(chain=c) == 0
  m.remove_chain(chain=c)
  assert c.parent() is None
  assert m.find_chain_index(chain=c) == -1
  try: m.find_chain_index(chain=c, must_be_present=True)
  except RuntimeError, e:
    assert str(e) == "chain not in model."
  else: raise Exception_expected
  #
  m1 = pdb.hierarchy.model()
  m2 = pdb.hierarchy.model()
  c = pdb.hierarchy.chain()
  m1.append_chain(chain=c)
  try: m2.append_chain(chain=c)
  except RuntimeError, e:
    assert str(e) == "chain has another parent model already."
  else: raise Exception_expected

def exercise_root():
  r = pdb.hierarchy.root()
  m = pdb.hierarchy.model()
  assert m.parent() is None
  m = pdb.hierarchy.model(parent=r)
  assert m.parent().memory_id() == r.memory_id()
  assert m.id == ""
  m = pdb.hierarchy.model(parent=r, id="2")
  assert m.parent(optional=False).memory_id() == r.memory_id()
  assert m.id == "2"
  del r
  assert m.parent() is None
  try:
    m.parent(optional=False)
  except RuntimeError, e:
    assert not show_diff(str(e), "model has no parent root")
  else: raise Exception_expected
  #
  r = pdb.hierarchy.root()
  assert r.info.size() == 0
  r.info.append("abc")
  assert r.info.size() == 1
  r.info = flex.std_string(["a", "b"])
  assert r.info.size() == 2
  r.pre_allocate_models(number_of_additional_models=2)
  assert r.models_size() == 0
  assert len(r.models()) == 0
  m_a = pdb.hierarchy.model(id="3")
  r.append_model(model=m_a)
  assert m_a.id == "3"
  assert m_a.parent().memory_id() == r.memory_id()
  assert r.models_size() == 1
  assert len(r.models()) == 1
  m_b = pdb.hierarchy.model(id="5")
  assert m_b.parent() is None
  r.append_model(model=m_b)
  assert r.models_size() == 2
  models = r.models()
  assert len(models) == 2
  assert models[0].memory_id() == m_a.memory_id()
  assert models[1].memory_id() == m_b.memory_id()
  for i in xrange(3):
    r.append_model(model=pdb.hierarchy.model())
  assert r.models_size() == 5
  assert len(r.models()) == 5
  for model in r.models():
    assert model.parent().memory_id() == r.memory_id()
  assert r.atoms_size() == 0
  assert r.atoms().size() == 0
  #
  rc = r.deep_copy()
  assert rc.memory_id() != r.memory_id()
  assert list(rc.info) == ["a", "b"]
  assert rc.info.id() != r.info.id()
  assert rc.models_size() == 5
  assert rc.models()[0].memory_id() != r.models()[0].memory_id()
  assert rc.models()[0].id == "3"
  r.append_model(model=pdb.hierarchy.model(id="7"))
  assert r.models_size() == 6
  assert rc.models_size() == 5
  assert [m.id for m in r.models()] == ["3", "5", "", "", "", "7"]
  assert [m.id for m in rc.models()] == ["3", "5", "", "", ""]
  rc.append_model(model=pdb.hierarchy.model(id="8"))
  assert r.models_size() == 6
  assert rc.models_size() == 6
  assert [m.id for m in rc.models()] == ["3", "5", "", "", "", "8"]
  #
  r = rc.deep_copy()
  r.insert_model(i=4, model=pdb.hierarchy.model(id="M"))
  assert [m.id for m in r.models()] \
      == ["3", "5", "", "", "M", "", "8"]
  r.remove_model(i=1)
  assert [m.id for m in r.models()] \
      == ["3", "", "", "M", "", "8"]
  m = r.models()[-1]
  assert m.parent().memory_id() == r.memory_id()
  assert r.find_model_index(model=m) == 5
  r.remove_model(model=m)
  assert m.parent() is None
  assert r.find_model_index(model=m) == -1
  try: r.find_model_index(model=m, must_be_present=True)
  except RuntimeError, e:
    assert str(e) == "model not in root."
  else: raise Exception_expected
  #
  r1 = pdb.hierarchy.root()
  r2 = pdb.hierarchy.root()
  m = pdb.hierarchy.model()
  r1.append_model(model=m)
  try: r2.append_model(model=m)
  except RuntimeError, e:
    assert str(e) == "model has another parent root already."
  else: raise Exception_expected

def exercise_atom_id_str():
  a = pdb.hierarchy.atom()
  a.set_name(new_name="NaMe")
  a.set_serial(new_serial="B1234")
  assert a.id_str() == 'pdb="NaMe           "'
  assert a.id_str(pdbres=True) == 'pdbres="          "'
  ag = pdb.hierarchy.atom_group(altloc="A", resname="GLY")
  ag.append_atom(a)
  assert a.id_str() == 'pdb="NaMeAGLY       "'
  assert a.id_str(True) == 'pdbres="GLY       "'
  rg = pdb.hierarchy.residue_group(resseq="1234", icode="J")
  rg.append_atom_group(ag)
  assert a.id_str() == 'pdb="NaMeAGLY  1234J"'
  assert a.id_str(pdbres=True) == 'pdbres="GLY  1234J"'
  ch = pdb.hierarchy.chain(id="Ch")
  ch.append_residue_group(rg)
  assert a.id_str() == 'pdb="NaMeAGLYCh1234J"'
  assert a.id_str(pdbres=True) == 'pdbres="GLYCh1234J"'
  md = pdb.hierarchy.model()
  md.append_chain(ch)
  assert a.id_str() == 'pdb="NaMeAGLYCh1234J"'
  assert a.id_str(pdbres=True) == 'pdbres="GLYCh1234J"'
  md.id = ""
  assert a.id_str() == 'pdb="NaMeAGLYCh1234J"'
  assert a.id_str(pdbres=True) == 'pdbres="GLYCh1234J"'
  md.id = "1"
  assert a.id_str() == 'model="   1" pdb="NaMeAGLYCh1234J"'
  assert a.id_str(pdbres=True) == 'model="   1" pdbres="GLYCh1234J"'
  md.id = "12345678"
  assert a.id_str() == 'model="12345678" pdb="NaMeAGLYCh1234J"'
  assert a.id_str(pdbres=True) == 'model="12345678" pdbres="GLYCh1234J"'
  try: md.id = "123456789"
  except (ValueError, RuntimeError), e:
    assert str(e) == "string is too long for target variable " \
      "(maximum length is 8 characters, 9 given)."
  else: raise Exception_expected
  md.id = "12345678"
  a.segid = "1234"
  assert a.id_str(suppress_segid=False) \
      == 'model="12345678" pdb="NaMeAGLYCh1234J" segid="1234"'
  assert a.id_str(suppress_segid=True) \
      == 'model="12345678" pdb="NaMeAGLYCh1234J"'
  assert a.id_str(pdbres=True) \
      == 'model="12345678" pdbres="GLYCh1234J" segid="1234"'
  assert a.id_str(pdbres=True, suppress_segid=True) \
      == 'model="12345678" pdbres="GLYCh1234J"'
  md.id = ""
  assert a.id_str() == 'pdb="NaMeAGLYCh1234J" segid="1234"'
  assert a.id_str(pdbres=True) \
      == 'pdbres="GLYCh1234J" segid="1234"'
  assert a.id_str(pdbres=True, suppress_segid=True) \
      == 'pdbres="GLYCh1234J"'
  rt = pdb.hierarchy.root()
  rt.append_model(md)
  assert a.id_str() == 'pdb="NaMeAGLYCh1234J" segid="1234"'
  assert a.id_str(pdbres=True) == 'pdbres="GLYCh1234J" segid="1234"'
  md.id = "    "
  assert a.id_str() == 'model="    " pdb="NaMeAGLYCh1234J" segid="1234"'
  assert a.id_str(pdbres=True) \
      == 'model="    " pdbres="GLYCh1234J" segid="1234"'
  #
  cf = ch.only_conformer()
  rd = cf.only_residue()
  assert rd.id_str() == 'model="    " pdbres="GLYCh1234J" segid="1234"'
  assert rd.id_str(suppress_segid=1) == 'model="    " pdbres="GLYCh1234J"'
  md.id = "12345678"
  assert rd.id_str() == 'model="12345678" pdbres="GLYCh1234J" segid="1234"'
  assert rd.id_str(suppress_segid=1) == 'model="12345678" pdbres="GLYCh1234J"'
  del cf
  assert rd.id_str(suppress_segid=-1) == 'pdbres="GLY  1234J" segid="1234"'
  assert rd.id_str(suppress_segid=1) == 'pdbres="GLY  1234J"'
  #
  a2 = pdb.hierarchy.atom().set_segid(new_segid="abcd")
  ag.append_atom(atom=a2)
  cf = ch.only_conformer()
  rd = cf.only_residue()
  assert rd.id_str(suppress_segid=1) == 'model="12345678" pdbres="GLYCh1234J"'
  assert rd.id_str(suppress_segid=-1) == 'model="12345678" pdbres="GLYCh1234J"'
  try: rd.id_str()
  except ValueError, e:
    assert not show_diff(str(e), '''\
residue.id_str(suppress_segid=false): segid is not unique:
  model="12345678" pdbres="GLYCh1234J" segid="1234"''')
  else: raise Exception_expected

def exercise_format_atom_record():
  a = (pdb.hierarchy.atom()
    .set_name(new_name="NaMe")
    .set_serial(new_serial="B1234")
    .set_xyz(new_xyz=(1.3,2.1,3.2))
    .set_sigxyz(new_sigxyz=(.1,.2,.3))
    .set_occ(new_occ=0.4)
    .set_sigocc(new_sigocc=0.1)
    .set_b(new_b=4.8)
    .set_sigb(new_sigb=0.7)
    .set_uij(new_uij=(1.3,2.1,3.2,4.3,2.7,9.3)))
  if (pdb.hierarchy.atom.has_siguij()):
    assert a.set_siguij(new_siguij=(.1,.2,.3,.6,.1,.9)) is a
  for hetero,record_name in [(False, "ATOM  "), (True, "HETATM")]:
    a.set_hetero(new_hetero=hetero)
    for segid in ["", "    ", "s", "sEgI"]:
      a.set_segid(new_segid=segid)
      for element in ["", "  ", "e", "El"]:
        a.set_element(new_element=element)
        for charge in ["", "  ", "c", "cH"]:
          a.set_charge(new_charge=charge)
          segielch = "%-4s%2s%-2s" % (segid, element, charge)
          s = a.format_atom_record()
          assert not show_diff(s, ("""%s\
B1234 NaMe                 1.300   2.100   3.200  0.40  4.80      \
%s""" % (record_name, segielch)).rstrip())
          sc = a.format_atom_record(replace_floats_with=None)
          assert not show_diff(sc, s)
          sc = a.format_atom_record(replace_floats_with="")
          assert not show_diff(sc, "%-27.27s%-8.8s" % (s[:27], s[72:]))
          sc = a.format_atom_record(replace_floats_with=" ")
          assert not show_diff(sc, "%-27.27s %-8.8s" % (s[:27], s[72:]))
          assert not show_diff(a.format_sigatm_record(), ("""SIGATM\
B1234 NaMe                 0.100   0.200   0.300  0.10  0.70      \
%s""" % segielch).rstrip())
          assert not show_diff(a.format_anisou_record(), ("""ANISOU\
B1234 NaMe              13000  21000  32000  43000  27000  93000  \
%s""" % segielch).rstrip())
          if (pdb.hierarchy.atom.has_siguij()):
            assert not show_diff(a.format_siguij_record(), ("""SIGUIJ\
B1234 NaMe               1000   2000   3000   6000   1000   9000  \
%s""" % segielch).rstrip())
          else:
            assert not show_diff(a.format_siguij_record(), ("""SIGUIJ\
B1234 NaMe             -10000 -10000 -10000 -10000 -10000 -10000  \
%s""" % segielch).rstrip())
          ag = pdb.hierarchy.atom_group(altloc="x", resname="uvw")
          ag.append_atom(atom=a)
          s = a.format_atom_record()
          assert not show_diff(s, ("""%s\
B1234 NaMexuvw             1.300   2.100   3.200  0.40  4.80      \
%s""" % (record_name, segielch)).rstrip())
          sc = a.format_atom_record(replace_floats_with=".*.")
          assert not show_diff(sc, "%-27.27s.*.%-8.8s" % (s[:27], s[72:]))
          assert not show_diff(a.format_sigatm_record(), ("""SIGATM\
B1234 NaMexuvw             0.100   0.200   0.300  0.10  0.70      \
%s""" % segielch).rstrip())
          assert not show_diff(a.format_anisou_record(), ("""ANISOU\
B1234 NaMexuvw          13000  21000  32000  43000  27000  93000  \
%s""" % segielch).rstrip())
          if (pdb.hierarchy.atom.has_siguij()):
            assert not show_diff(a.format_siguij_record(), ("""SIGUIJ\
B1234 NaMexuvw           1000   2000   3000   6000   1000   9000  \
%s""" % segielch).rstrip())
          else:
            assert not show_diff(a.format_siguij_record(), ("""SIGUIJ\
B1234 NaMexuvw         -10000 -10000 -10000 -10000 -10000 -10000  \
%s""" % segielch).rstrip())
          rg = pdb.hierarchy.residue_group(resseq="pqrs", icode="t")
          rg.append_atom_group(atom_group=ag)
          s = a.format_atom_record()
          assert not show_diff(s, ("""%s\
B1234 NaMexuvw  pqrst      1.300   2.100   3.200  0.40  4.80      \
%s""" % (record_name, segielch)).rstrip())
          sc = a.format_atom_record(replace_floats_with=".*.")
          assert not show_diff(sc, "%-27.27s.*.%-8.8s" % (s[:27], s[72:]))
          assert not show_diff(a.format_sigatm_record(), ("""SIGATM\
B1234 NaMexuvw  pqrst      0.100   0.200   0.300  0.10  0.70      \
%s""" % segielch).rstrip())
          assert not show_diff(a.format_anisou_record(), ("""ANISOU\
B1234 NaMexuvw  pqrst   13000  21000  32000  43000  27000  93000  \
%s""" % segielch).rstrip())
          if (pdb.hierarchy.atom.has_siguij()):
            assert not show_diff(a.format_siguij_record(), ("""SIGUIJ\
B1234 NaMexuvw  pqrst    1000   2000   3000   6000   1000   9000  \
%s""" % segielch).rstrip())
          else:
            assert not show_diff(a.format_siguij_record(), ("""SIGUIJ\
B1234 NaMexuvw  pqrst  -10000 -10000 -10000 -10000 -10000 -10000  \
%s""" % segielch).rstrip())
          for chain_id in ["", "g", "hi"]:
            ch = pdb.hierarchy.chain(id=chain_id)
            ch.append_residue_group(residue_group=rg)
            s = a.format_atom_record()
            assert not show_diff(s, ("""%s\
B1234 NaMexuvw%2spqrst      1.300   2.100   3.200  0.40  4.80      \
%s""" % (record_name, chain_id, segielch)).rstrip())
            sc = a.format_atom_record(replace_floats_with=".*.")
            assert not show_diff(sc, "%-27.27s.*.%-8.8s" % (s[:27], s[72:]))
            assert not show_diff(a.format_sigatm_record(), ("""SIGATM\
B1234 NaMexuvw%2spqrst      0.100   0.200   0.300  0.10  0.70      \
%s""" % (chain_id, segielch)).rstrip())
            assert not show_diff(a.format_anisou_record(), ("""ANISOU\
B1234 NaMexuvw%2spqrst   13000  21000  32000  43000  27000  93000  \
%s""" % (chain_id, segielch)).rstrip())
            if (pdb.hierarchy.atom.has_siguij()):
              assert not show_diff(a.format_siguij_record(), ("""SIGUIJ\
B1234 NaMexuvw%2spqrst    1000   2000   3000   6000   1000   9000  \
%s""" % (chain_id, segielch)).rstrip())
            else:
              assert not show_diff(a.format_siguij_record(), ("""SIGUIJ\
B1234 NaMexuvw%2spqrst  -10000 -10000 -10000 -10000 -10000 -10000  \
%s""" % (chain_id, segielch)).rstrip())
          #
          assert not show_diff(
            a.quote(), '"'+a.format_atom_record(replace_floats_with=".*.")+'"')
          assert not show_diff(
            a.quote(full=True), '"'+a.format_atom_record()+'"')
          #
          del ag, rg, ch
          #
          rs = [
            a.format_atom_record(),
            a.format_sigatm_record(),
            a.format_anisou_record()]
          if (a.has_siguij()):
            rs.append(a.format_siguij_record())
          assert not show_diff(
            a.format_atom_record_group(), "\n".join(rs))
          assert not show_diff(
            a.format_atom_record_group(atom_hetatm=False), "\n".join(rs[1:]))
          assert not show_diff(
            a.format_atom_record_group(sigatm=False), "\n".join(rs[:1]+rs[2:]))
          assert not show_diff(
            a.format_atom_record_group(anisou=False), "\n".join(rs[:2]+rs[3:]))
          if (a.has_siguij()):
            assert not show_diff(
              a.format_atom_record_group(siguij=False), "\n".join(rs[:-1]))
          else:
            assert not show_diff(
              a.format_atom_record_group(siguij=False), "\n".join(rs))
  #
  atom = pdb.hierarchy.atom()
  atom.name = "NaMe"
  atom.xyz = (1.e6,-1.e5,1.e5)
  atom.occ = 9999.
  atom.b = 99999.
  assert not show_diff(atom.format_atom_record(), """\
ATOM        NaMe              1000000.-100000.100000.09999.099999.""")
  atom.xyz = (0,1.e7,0)
  try: atom.format_atom_record()
  except RuntimeError, e: assert not show_diff(str(e), """\
atom Y coordinate value does not fit into F8.3 format:
  "ATOM        NaMe           "
  value: 10000000.000""")
  else: raise Exception_expected
  atom.xyz = (0,0,0)
  atom.occ = 100000
  try: atom.format_atom_record()
  except RuntimeError, e: assert not show_diff(str(e), """\
atom occupancy factor does not fit into F6.2 format:
  "ATOM        NaMe           "
  occupancy factor: 100000.00""")
  else: raise Exception_expected
  atom.occ = 0
  atom.b = 200000
  try: atom.format_atom_record()
  except RuntimeError, e: assert not show_diff(str(e), """\
atom B-factor does not fit into F6.2 format:
  "ATOM        NaMe           "
  B-factor: 200000.00""")
  else: raise Exception_expected
  #
  atom.sigxyz = (1.e6,-1.e5,1.e5)
  atom.sigocc = 99999.
  atom.sigb = 9999.
  assert not show_diff(atom.format_sigatm_record(), """\
SIGATM      NaMe              1000000.-100000.100000.099999.9999.0""")
  atom.sigxyz = (0,0,2.e7)
  try: atom.format_sigatm_record()
  except RuntimeError, e: assert not show_diff(str(e), """\
atom sigma Z coordinate value does not fit into F8.3 format:
  "SIGATM      NaMe           "
  value: 20000000.000""")
  else: raise Exception_expected
  atom.sigxyz = (0,0,0)
  atom.sigocc = 200000
  try: atom.format_sigatm_record()
  except RuntimeError, e: assert not show_diff(str(e), """\
atom sigma occupancy factor does not fit into F6.2 format:
  "SIGATM      NaMe           "
  sigma occupancy factor: 200000.00""")
  else: raise Exception_expected
  atom.sigocc = 0
  atom.sigb = 300000
  try: atom.format_sigatm_record()
  except RuntimeError, e: assert not show_diff(str(e), """\
atom sigma B-factor does not fit into F6.2 format:
  "SIGATM      NaMe           "
  sigma B-factor: 300000.00""")
  else: raise Exception_expected
  #
  atom.uij = (0,1000,0,0,0,0)
  try: atom.format_anisou_record()
  except RuntimeError, e: assert not show_diff(str(e), """\
atom U22 value * 10000 does not fit into F7.0 format:
  "ANISOU      NaMe           "
  value * 10000: 10000000""")
  else: raise Exception_expected
  #
  atom.siguij = (0,0,0,0,3333,0)
  if (pdb.hierarchy.atom.has_siguij()):
    try: atom.format_siguij_record()
    except RuntimeError, e: assert not show_diff(str(e), """\
atom sigma U13 value * 10000 does not fit into F7.0 format:
  "SIGUIJ      NaMe           "
  value * 10000: 33330000""")
    else: raise Exception_expected
  else:
    assert not show_diff(atom.format_siguij_record(), """\
SIGUIJ      NaMe             -10000 -10000 -10000 -10000 -10000 -10000""")

def exercise_construct_hierarchy():
  def check(pdb_string,
        expected_root_as_str=None,
        expected_overall_counts_as_str=None,
        level_id=None,
        prefix=""):
    pdb_inp = pdb.input(source_info=None, lines=flex.split_lines(pdb_string))
    root = pdb_inp.construct_hierarchy()
    if (expected_root_as_str is not None):
      s = root.as_str(prefix=prefix, level_id=level_id)
      if (len(expected_root_as_str) == 0):
        sys.stdout.write(s)
      else:
        assert not show_diff(s, expected_root_as_str)
    if (expected_overall_counts_as_str is not None):
      s = root.overall_counts().as_str(
        prefix=prefix,
        residue_groups_max_show=3,
        duplicate_atom_labels_max_show=3)
      if (len(expected_overall_counts_as_str) == 0):
        sys.stdout.write(s)
      else:
        assert not show_diff(s, expected_overall_counts_as_str)
    return root
  #
  check("""\
MODEL        1
ATOM      1  N   MET A   1       6.215  22.789  24.067  1.00  0.00           N
ATOM      2  CA  MET A   1       6.963  22.789  22.822  1.00  0.00           C
HETATM    3  C   MET A   2       7.478  21.387  22.491  1.00  0.00           C
ATOM      4  O   MET A   2       8.406  20.895  23.132  1.00  0.00           O
ENDMDL
MODEL 3
HETATM    9 2H3  MPR B   5      16.388   0.289   6.613  1.00  0.08
SIGATM    9 2H3  MPR B   5       0.155   0.175   0.155  0.00  0.05
ANISOU    9 2H3  MPR B   5      848    848    848      0      0      0
SIGUIJ    9 2H3  MPR B   5      510    510    510      0      0      0
TER
ATOM     10  N   CYSCH   6      14.270   2.464   3.364  1.00  0.07
SIGATM   10  N   CYSCH   6       0.012   0.012   0.011  0.00  0.00
ANISOU   10  N   CYSCH   6      788    626    677   -344    621   -232
SIGUIJ   10  N   CYSCH   6        3     13      4     11      6     13
TER
ENDMDL
END
""", """\
model id="   1" #chains=1
  chain id="A" #residue_groups=2
    resid="   1 " #atom_groups=1
      altloc="" resname="MET" #atoms=2
        " N  "
        " CA "
    resid="   2 " #atom_groups=1
      altloc="" resname="MET" #atoms=2
        " C  "
        " O  "
model id="   3" #chains=2
  chain id="B" #residue_groups=1
    resid="   5 " #atom_groups=1
      altloc="" resname="MPR" #atoms=1
        "2H3 "
  chain id="CH" #residue_groups=1
    resid="   6 " #atom_groups=1
      altloc="" resname="CYS" #atoms=1
        " N  "
""", """\
total number of:
  models:     2
  chains:     3
  alt. conf.: 0
  residues:   4
  atoms:      6
  anisou:     2
number of atom element+charge types: 4
histogram of atom element+charge frequency:
  "    " 2
  " C  " 2
  " N  " 1
  " O  " 1
residue name classes:
  "common_amino_acid" 3
  "other"             1
number of chain ids: 3
histogram of chain id frequency:
  "A"  1
  "B"  1
  "CH" 1
number of alt. conf. ids: 0
number of residue names: 3
histogram of residue name frequency:
  "MET" 2
  "CYS" 1
  "MPR" 1    other
""")
  #
  check("""\
ATOM         N1 AR01
ATOM         N2 BR01
ATOM         N1 CR02
ATOM         N2  R02
""", """\
model id="" #chains=1
  chain id=" " #residue_groups=2
    resid="     " #atom_groups=2
      altloc="A" resname="R01" #atoms=1
        " N1 "
      altloc="B" resname="R01" #atoms=1
        " N2 "
    resid="     " #atom_groups=2  ### Info: same as previous resid ###
      altloc="" resname="R02" #atoms=1
        " N2 "
      altloc="C" resname="R02" #atoms=1
        " N1 "
""")
  #
  check("""\
ATOM         N1 BR01
ATOM         N1  R01
ATOM         N2  R01
ATOM         N3 BR01
ATOM         N3  R01
ATOM         N1  R02
ATOM         N1 BR02
ATOM         N2  R02
ATOM         N3 BR02
ATOM         N3  R02
ATOM         N1  R03
ATOM         N1 BR03
ATOM         N2  R03
ATOM         N3 BR03
ATOM         N3  R03
""", """\
model id="" #chains=1
  chain id=" " #residue_groups=3
    resid="     " #atom_groups=3
      altloc="" resname="R01" #atoms=1
        " N2 "
      altloc=" " resname="R01" #atoms=2
        " N1 "
        " N3 "
      altloc="B" resname="R01" #atoms=2
        " N1 "
        " N3 "
    resid="     " #atom_groups=3  ### Info: same as previous resid ###
      altloc="" resname="R02" #atoms=1
        " N2 "
      altloc=" " resname="R02" #atoms=2
        " N1 "
        " N3 "
      altloc="B" resname="R02" #atoms=2
        " N1 "
        " N3 "
    resid="     " #atom_groups=3  ### Info: same as previous resid ###
      altloc="" resname="R03" #atoms=1
        " N2 "
      altloc=" " resname="R03" #atoms=2
        " N1 "
        " N3 "
      altloc="B" resname="R03" #atoms=2
        " N1 "
        " N3 "
""")
  #
  root = check("""\
ATOM         N1 BR01
ATOM         N1  R01
ATOM         N2  R01
ATOM         N3 BR01
ATOM         N3  R01
ATOM         N1 AR02
ATOM         N1 BR02
ATOM         N2  R02
ATOM         N3 BR02
ATOM         N3 AR02
ATOM         N1  R03
ATOM         N1 BR03
ATOM         N2  R03
ATOM         N3 BR03
ATOM         N3  R03
""", """\
  model id="" #chains=1
    chain id=" " #residue_groups=3
      resid="     " #atom_groups=3
        altloc="" resname="R01" #atoms=1
          " N2 "
        altloc=" " resname="R01" #atoms=2
          " N1 "
          " N3 "
        altloc="B" resname="R01" #atoms=2
          " N1 "
          " N3 "
      resid="     " #atom_groups=3  ### Info: same as previous resid ###
        altloc="" resname="R02" #atoms=1
          " N2 "
        altloc="A" resname="R02" #atoms=2
          " N1 "
          " N3 "
        altloc="B" resname="R02" #atoms=2
          " N1 "
          " N3 "
      resid="     " #atom_groups=3  ### Info: same as previous resid ###
        altloc="" resname="R03" #atoms=1
          " N2 "
        altloc=" " resname="R03" #atoms=2
          " N1 "
          " N3 "
        altloc="B" resname="R03" #atoms=2
          " N1 "
          " N3 "
""", """\
  total number of:
    models:      1
    chains:      1
    alt. conf.:  3
    residues:    3
    atoms:      15
    anisou:      0
  number of atom element+charge types: 1
  histogram of atom element+charge frequency:
    "    " 15
  residue name classes:
    "other" 3
  number of chain ids: 1
  histogram of chain id frequency:
    " " 1
  number of alt. conf. ids: 3
  histogram of alt. conf. id frequency:
    " " 1
    "A" 1
    "B" 1
  residue alt. conf. situations:
    pure main conf.:     0
    pure alt. conf.:     0
    proper alt. conf.:   1
    ### ERROR: improper alt. conf. ###
    improper alt. conf.: 2
  chains with mix of proper and improper alt. conf.: 1
    residue with proper altloc
      "ATOM         N2  R02       .*.        "
      "ATOM         N1 AR02       .*.        "
      "ATOM         N3 AR02       .*.        "
      "ATOM         N1 BR02       .*.        "
      "ATOM         N3 BR02       .*.        "
    residue with improper altloc
      "ATOM         N2  R01       .*.        "
      "ATOM         N1  R01       .*.        "
      "ATOM         N3  R01       .*.        "
      "ATOM         N1 BR01       .*.        "
      "ATOM         N3 BR01       .*.        "
  number of residue names: 3
  histogram of residue name frequency:
    "R01" 1    other
    "R02" 1    other
    "R03" 1    other
  ### WARNING: consecutive residue_groups with same resid ###
  number of consecutive residue groups with same resid: 2
    residue group:
      "ATOM         N2  R01       .*.        "
      ... 3 atoms not shown
      "ATOM         N3 BR01       .*.        "
    next residue group:
      "ATOM         N2  R02       .*.        "
      ... 3 atoms not shown
      "ATOM         N3 BR02       .*.        "
    next residue group:
      "ATOM         N2  R03       .*.        "
      ... 3 atoms not shown
      "ATOM         N3 BR03       .*.        "
""", prefix="  ")
  oc = root.overall_counts()
  assert oc.warnings() == [
    '### WARNING: consecutive residue_groups with same resid ###']
  assert oc.errors() == ['### ERROR: improper alt. conf. ###']
  assert oc.errors_and_warnings() == [
    '### ERROR: improper alt. conf. ###',
    '### WARNING: consecutive residue_groups with same resid ###']
  #
  check("""\
ATOM         N1 BR01
ATOM         N1 AR01
ATOM         N2 CR01
ATOM         N3 BR01
ATOM         N3 AR01
ATOM         N1  R02
ATOM         N1 BR02
ATOM         N2  R02
ATOM         N3 BR02
ATOM         N3  R02
ATOM         N1 CR03
ATOM         N1 BR03
ATOM         N2 BR03
ATOM         N2 CR03
ATOM         N3  R03
""", """\
model id="" #chains=1
  chain id=" " #residue_groups=3
    resid="     " #atom_groups=3
      altloc="B" resname="R01" #atoms=2
        " N1 "
        " N3 "
      altloc="A" resname="R01" #atoms=2
        " N1 "
        " N3 "
      altloc="C" resname="R01" #atoms=1
        " N2 "
    resid="     " #atom_groups=3  ### Info: same as previous resid ###
      altloc="" resname="R02" #atoms=1
        " N2 "
      altloc=" " resname="R02" #atoms=2
        " N1 "
        " N3 "
      altloc="B" resname="R02" #atoms=2
        " N1 "
        " N3 "
    resid="     " #atom_groups=3  ### Info: same as previous resid ###
      altloc="" resname="R03" #atoms=1
        " N3 "
      altloc="C" resname="R03" #atoms=2
        " N1 "
        " N2 "
      altloc="B" resname="R03" #atoms=2
        " N1 "
        " N2 "
""", """\
total number of:
  models:      1
  chains:      1
  alt. conf.:  4
  residues:    3
  atoms:      15
  anisou:      0
number of atom element+charge types: 1
histogram of atom element+charge frequency:
  "    " 15
residue name classes:
  "other" 3
number of chain ids: 1
histogram of chain id frequency:
  " " 1
number of alt. conf. ids: 4
histogram of alt. conf. id frequency:
  " " 1
  "A" 1
  "B" 1
  "C" 1
residue alt. conf. situations:
  pure main conf.:     0
  pure alt. conf.:     1
  proper alt. conf.:   1
  ### ERROR: improper alt. conf. ###
  improper alt. conf.: 1
chains with mix of proper and improper alt. conf.: 1
  residue with proper altloc
    "ATOM         N3  R03       .*.        "
    "ATOM         N1 CR03       .*.        "
    "ATOM         N2 CR03       .*.        "
    "ATOM         N1 BR03       .*.        "
    "ATOM         N2 BR03       .*.        "
  residue with improper altloc
    "ATOM         N2  R02       .*.        "
    "ATOM         N1  R02       .*.        "
    "ATOM         N3  R02       .*.        "
    "ATOM         N1 BR02       .*.        "
    "ATOM         N3 BR02       .*.        "
number of residue names: 3
histogram of residue name frequency:
  "R01" 1    other
  "R02" 1    other
  "R03" 1    other
### WARNING: consecutive residue_groups with same resid ###
number of consecutive residue groups with same resid: 2
  residue group:
    "ATOM         N1 BR01       .*.        "
    ... 3 atoms not shown
    "ATOM         N2 CR01       .*.        "
  next residue group:
    "ATOM         N2  R02       .*.        "
    ... 3 atoms not shown
    "ATOM         N3 BR02       .*.        "
  next residue group:
    "ATOM         N3  R03       .*.        "
    ... 3 atoms not shown
    "ATOM         N2 BR03       .*.        "
""")
  #
  check("""\
REMARK    ANTIBIOTIC                              26-JUL-06   2IZQ
ATOM    220  N  ATRP A  11      20.498  12.832  34.558  0.50  6.03           N
ATOM    221  CA ATRP A  11      21.094  12.032  35.602  0.50  5.24           C
ATOM    222  C  ATRP A  11      22.601  12.088  35.532  0.50  6.49           C
ATOM    223  O  ATRP A  11      23.174  12.012  34.439  0.50  7.24           O
ATOM    234  H  ATRP A  11      20.540  12.567  33.741  0.50  7.24           H
ATOM    235  HA ATRP A  11      20.771  12.306  36.485  0.50  6.28           H
ATOM    244  N  CPHE A  11      20.226  13.044  34.556  0.15  6.35           N
ATOM    245  CA CPHE A  11      20.950  12.135  35.430  0.15  5.92           C
ATOM    246  C  CPHE A  11      22.448  12.425  35.436  0.15  6.32           C
ATOM    247  O  CPHE A  11      22.961  12.790  34.373  0.15  6.08           O
ATOM    255  N  BTYR A  11      20.553  12.751  34.549  0.35  5.21           N
ATOM    256  CA BTYR A  11      21.106  11.838  35.524  0.35  5.51           C
ATOM    257  C  BTYR A  11      22.625  11.920  35.572  0.35  5.42           C
ATOM    258  O  BTYR A  11      23.299  11.781  34.538  0.35  5.30           O
ATOM    262  HB2CPHE A  11      21.221  10.536  34.146  0.15  7.21           H
ATOM    263  CD2BTYR A  11      18.463  10.012  36.681  0.35  9.08           C
ATOM    264  HB3CPHE A  11      21.198  10.093  35.647  0.15  7.21           H
ATOM    265  CE1BTYR A  11      17.195   9.960  34.223  0.35 10.76           C
ATOM    266  HD1CPHE A  11      19.394   9.937  32.837  0.15 10.53           H
ATOM    267  CE2BTYR A  11      17.100   9.826  36.693  0.35 11.29           C
ATOM    268  HD2CPHE A  11      18.873  10.410  36.828  0.15  9.24           H
ATOM    269  CZ BTYR A  11      16.546   9.812  35.432  0.35 11.90           C
ATOM    270  HE1CPHE A  11      17.206   9.172  32.650  0.15 12.52           H
ATOM    271  OH BTYR A  11      15.178   9.650  35.313  0.35 19.29           O
ATOM    272  HE2CPHE A  11      16.661   9.708  36.588  0.15 11.13           H
ATOM    273  HZ CPHE A  11      15.908   9.110  34.509  0.15 13.18           H
ATOM    274  H  BTYR A  11      20.634  12.539  33.720  0.35  6.25           H
ATOM    275  HA BTYR A  11      20.773  12.116  36.402  0.35  6.61           H
ATOM    276  HB2BTYR A  11      20.949  10.064  34.437  0.35  6.78           H
""", """\
model id="" #chains=1
  chain id="A" #residue_groups=1
    resid="  11 " #atom_groups=3  ### Info: with mixed residue names ###
      altloc="A" resname="TRP" #atoms=6
      altloc="C" resname="PHE" #atoms=11
      altloc="B" resname="TYR" #atoms=12
""", """\
total number of:
  models:      1
  chains:      1
  alt. conf.:  3
  residues:    1 (1 with mixed residue names)
  atoms:      29
  anisou:      0
number of atom element+charge types: 4
histogram of atom element+charge frequency:
  " H  " 12
  " C  " 10
  " O  "  4
  " N  "  3
residue name classes:
  "common_amino_acid" 3
number of chain ids: 1
histogram of chain id frequency:
  "A" 1
number of alt. conf. ids: 3
histogram of alt. conf. id frequency:
  "A" 1
  "B" 1
  "C" 1
residue alt. conf. situations:
  pure main conf.:     0
  pure alt. conf.:     1
  proper alt. conf.:   0
  improper alt. conf.: 0
chains with mix of proper and improper alt. conf.: 0
number of residue names: 3
histogram of residue name frequency:
  "PHE" 1
  "TRP" 1
  "TYR" 1
""", level_id="atom_group")
  #
  root = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
BREAK
""")).construct_hierarchy()
  assert root.models_size() == 0
  root = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
BREAK
ATOM      1  CB  LYS   109
BREAK
TER
""")).construct_hierarchy()
  assert not root.only_residue_group().link_to_previous
  root = pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
BREAK
ATOM      1  CB  LYS   109
ATOM      2  CG  LYS   109
BREAK
TER
""")).construct_hierarchy()
  assert not root.only_residue_group().link_to_previous
  pdb_str = """\
ATOM      1  CB  LYS   109
ATOM      2  CG  LYS   109
ATOM      3  CA  LYS   110
ATOM      4  CB  LYS   110
BREAK
ATOM      5  CA  LYS   111
ATOM      6  CB  LYS   111
ATOM      7  CA  LYS   112
ATOM      8  CB  LYS   112
"""
  lines = flex.split_lines(pdb_str)
  for i_proc in [0,1]:
    root = pdb.input(source_info=None, lines=lines).construct_hierarchy()
    residue_groups = root.only_chain().residue_groups()
    assert len(residue_groups) == 4
    assert not residue_groups[0].link_to_previous
    assert residue_groups[1].link_to_previous
    assert not residue_groups[2].link_to_previous
    assert residue_groups[3].link_to_previous
    if (i_proc == 0):
      lines = lines.select(flex.size_t([0,2,4,5,7]))
  try: pdb.input(
    source_info=None,
    lines=flex.split_lines("""\
REMARK
ATOM      1  CB  LYS   109
BREAK
ATOM      2  CG  LYS   109
""")).construct_hierarchy()
  except RuntimeError, e:
    assert not show_diff(str(e), "Misplaced BREAK record (input line 3).")
  else: raise Exception_expected
  try: pdb.input(
    source_info="file abc",
    lines=flex.split_lines("""\
REMARK
ATOM      1  CA  LYS   109
ATOM      2  CB  LYS   109
BREAK
ATOM      3  CA  LYS   110
BREAK
ATOM      4  CB  LYS   110
""")).construct_hierarchy()
  except RuntimeError, e:
    assert not show_diff(str(e), "Misplaced BREAK record (file abc, line 6).")
  else: raise Exception_expected
  #
  check(pdb_str, """\
:=model id="" #chains=1
:=  chain id=" " #residue_groups=4
:=    resid=" 109 " #atom_groups=1
:=      altloc="" resname="LYS" #atoms=2
:=        " CB "
:=        " CG "
:=    resid=" 110 " #atom_groups=1
:=      altloc="" resname="LYS" #atoms=2
:=        " CA "
:=        " CB "
:=    ### chain break ###
:=    resid=" 111 " #atom_groups=1
:=      altloc="" resname="LYS" #atoms=2
:=        " CA "
:=        " CB "
:=    resid=" 112 " #atom_groups=1
:=      altloc="" resname="LYS" #atoms=2
:=        " CA "
:=        " CB "
""", """\
:=total number of:
:=  models:     1
:=  chains:     1 (1 explicit chain break)
:=  alt. conf.: 0
:=  residues:   4
:=  atoms:      8
:=  anisou:     0
:=number of atom element+charge types: 1
:=histogram of atom element+charge frequency:
:=  "    " 8
:=residue name classes:
:=  "common_amino_acid" 4
:=number of chain ids: 1
:=histogram of chain id frequency:
:=  " " 1
:=number of alt. conf. ids: 0
:=number of residue names: 1
:=histogram of residue name frequency:
:=  "LYS" 4
""", prefix=":=")
  #
  check(pdb_str, """\
model id="" #chains=1
  chain id=" " #residue_groups=4
    resid=" 109 " #atom_groups=1
      altloc="" resname="LYS" #atoms=2
    resid=" 110 " #atom_groups=1
      altloc="" resname="LYS" #atoms=2
    ### chain break ###
    resid=" 111 " #atom_groups=1
      altloc="" resname="LYS" #atoms=2
    resid=" 112 " #atom_groups=1
      altloc="" resname="LYS" #atoms=2
""", level_id="atom_group")
  #
  check(pdb_str, """\
model id="" #chains=1
  chain id=" " #residue_groups=4
    resid=" 109 " #atom_groups=1
    resid=" 110 " #atom_groups=1
    ### chain break ###
    resid=" 111 " #atom_groups=1
    resid=" 112 " #atom_groups=1
""", level_id="residue_group")
  #
  check(pdb_str, """\
model id="" #chains=1
  chain id=" " #residue_groups=4
""", level_id="chain")
  #
  check(pdb_str, """\
model id="" #chains=1
""", level_id="model")
  #
  check("""\
MODEL        1
ENDMDL
MODEL 1
ENDMDL
MODEL     1
ENDMDL
""", """\
model id="   1" #chains=0  ### ERROR: duplicate model id ###
  ### WARNING: empty model ###
model id="   1" #chains=0  ### ERROR: duplicate model id ###
  ### WARNING: empty model ###
model id="   1" #chains=0  ### ERROR: duplicate model id ###
  ### WARNING: empty model ###
""", """\
total number of:
  ### ERROR: duplicate model ids ###
  ### WARNING: empty model ###
  models:     3 (3 with duplicate model ids; 3 empty)
  chains:     0
  alt. conf.: 0
  residues:   0
  atoms:      0
  anisou:     0
number of atom element+charge types: 0
residue name classes: None
number of chain ids: 0
number of alt. conf. ids: 0
number of residue names: 0
""")
  #
  check("""\
MODEL        1
ATOM                 A
ATOM                 B
ATOM                 A
ENDMDL
MODEL 1
ATOM                 A   1
BREAK
ATOM                 A   2
ATOM                 B
ATOM                 A
ENDMDL
MODEL     2
ATOM                 A
BREAK
ATOM                 A    I
ATOM                 B
ENDMDL
""", """\
model id="   1" #chains=3  ### ERROR: duplicate model id ###
  chain id="A" #residue_groups=1  ### WARNING: duplicate chain id ###
    resid="     " #atom_groups=1
      altloc="" resname="   " #atoms=1
        "    "
  chain id="B" #residue_groups=1
    resid="     " #atom_groups=1
      altloc="" resname="   " #atoms=1
        "    "
  chain id="A" #residue_groups=1  ### WARNING: duplicate chain id ###
    resid="     " #atom_groups=1
      altloc="" resname="   " #atoms=1
        "    "
model id="   1" #chains=3  ### ERROR: duplicate model id ###
  chain id="A" #residue_groups=2  ### WARNING: duplicate chain id ###
    resid="   1 " #atom_groups=1
      altloc="" resname="   " #atoms=1
        "    "
    ### chain break ###
    resid="   2 " #atom_groups=1
      altloc="" resname="   " #atoms=1
        "    "
  chain id="B" #residue_groups=1
    resid="     " #atom_groups=1
      altloc="" resname="   " #atoms=1
        "    "
  chain id="A" #residue_groups=1  ### WARNING: duplicate chain id ###
    resid="     " #atom_groups=1
      altloc="" resname="   " #atoms=1
        "    "
model id="   2" #chains=2
  chain id="A" #residue_groups=2
    resid="     " #atom_groups=1
      altloc="" resname="   " #atoms=1
        "    "
    ### chain break ###
    resid="    I" #atom_groups=1
      altloc="" resname="   " #atoms=1
        "    "
  chain id="B" #residue_groups=1
    resid="     " #atom_groups=1
      altloc="" resname="   " #atoms=1
        "    "
""", """\
total number of:
  ### ERROR: duplicate model ids ###
  models:      3 (2 with duplicate model ids)
  ### WARNING: duplicate chain ids ###
  chains:      8 (4 with duplicate chain ids; 2 explicit chain breaks)
  alt. conf.:  0
  residues:   10
  ### ERROR: duplicate atom labels ###
  atoms:      10 (2 with duplicate labels)
  anisou:      0
number of atom element+charge types: 1
histogram of atom element+charge frequency:
  "    " 10
residue name classes:
  "other" 10
number of chain ids: 2
histogram of chain id frequency:
  "A" 5
  "B" 3
number of alt. conf. ids: 0
number of residue names: 1
histogram of residue name frequency:
  "   " 10    other
number of groups of duplicate atom labels: 1
  total number of affected atoms:          2
  group "ATOM                 A     .*.        "
        "ATOM                 A     .*.        "
""")
  #
  check("""\
ATOM     54  CA  GLY A   9
ATOM     55  CA  GLY A   9
ATOM     56  CA BGLY A   9
""", """\
model id="" #chains=1
  chain id="A" #residue_groups=1
    resid="   9 " #atom_groups=2
      altloc=" " resname="GLY" #atoms=2
        " CA "
        " CA "
      altloc="B" resname="GLY" #atoms=1
        " CA "
""", """\
total number of:
  models:     1
  chains:     1
  alt. conf.: 2
  residues:   1
  ### ERROR: duplicate atom labels ###
  atoms:      3 (2 with duplicate labels)
  anisou:     0
number of atom element+charge types: 1
histogram of atom element+charge frequency:
  "    " 3
residue name classes:
  "common_amino_acid" 1
number of chain ids: 1
histogram of chain id frequency:
  "A" 1
number of alt. conf. ids: 2
histogram of alt. conf. id frequency:
  " " 1
  "B" 1
residue alt. conf. situations:
  pure main conf.:     0
  pure alt. conf.:     0
  proper alt. conf.:   0
  ### ERROR: improper alt. conf. ###
  improper alt. conf.: 1
chains with mix of proper and improper alt. conf.: 0
residue with improper altloc
  "ATOM     54  CA  GLY A   9 .*.        "
  "ATOM     55  CA  GLY A   9 .*.        "
  "ATOM     56  CA BGLY A   9 .*.        "
number of residue names: 1
histogram of residue name frequency:
  "GLY" 1
number of groups of duplicate atom labels: 1
  total number of affected atoms:          2
  group "ATOM     54  CA  GLY A   9 .*.        "
        "ATOM     55  CA  GLY A   9 .*.        "
""")
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM     68  HD1 LEU B 441
ATOM     69  HD1 LEU B 441
ATOM     70  HD1 LEU B 441
ATOM     71  HD2 LEU B 441
ATOM     72  HD2 LEU B 441
ATOM     73  HD2 LEU B 441
"""))
  oc = pdb_inp.construct_hierarchy().overall_counts()
  assert oc.errors() == ['### ERROR: duplicate atom labels ###']
  assert len(oc.warnings()) == 0
  oc.raise_improper_alt_conf_if_necessary()
  oc.raise_chains_with_mix_of_proper_and_improper_alt_conf_if_necessary()
  try: oc.raise_duplicate_atom_labels_if_necessary(max_show=1)
  except Sorry, e:
    assert not show_diff(str(e), '''\
number of groups of duplicate atom labels: 2
  total number of affected atoms:          6
  group "ATOM     68  HD1 LEU B 441 .*.        "
        "ATOM     69  HD1 LEU B 441 .*.        "
        "ATOM     70  HD1 LEU B 441 .*.        "
  ... 1 remaining group not shown''')
  else: raise Exception_expected
  sio = StringIO()
  oc.show_duplicate_atom_labels(out=sio)
  assert not show_diff(sio.getvalue(), """\
number of groups of duplicate atom labels: 2
  total number of affected atoms:          6
  group "ATOM     68  HD1 LEU B 441 .*.        "
        "ATOM     69  HD1 LEU B 441 .*.        "
        "ATOM     70  HD1 LEU B 441 .*.        "
  group "ATOM     71  HD2 LEU B 441 .*.        "
        "ATOM     72  HD2 LEU B 441 .*.        "
        "ATOM     73  HD2 LEU B 441 .*.        "
""")
  #
  check("""\
ATOM     68  HD1 LEU B 441
ATOM     69  HD2 LEU B 441
ATOM     70  HD3 LEU B 441
ATOM     71  HD4 LEU B 441
ATOM     72  HD1 LEU B 441
ATOM     73  HD2 LEU B 441
ATOM     74  HD3 LEU B 441
ATOM     75  HD4 LEU B 441
""", None, """\
total number of:
  models:     1
  chains:     1
  alt. conf.: 0
  residues:   1
  ### ERROR: duplicate atom labels ###
  atoms:      8 (8 with duplicate labels)
  anisou:     0
number of atom element+charge types: 1
histogram of atom element+charge frequency:
  "    " 8
residue name classes:
  "common_amino_acid" 1
number of chain ids: 1
histogram of chain id frequency:
  "B" 1
number of alt. conf. ids: 0
number of residue names: 1
histogram of residue name frequency:
  "LEU" 1
number of groups of duplicate atom labels: 4
  total number of affected atoms:          8
  group "ATOM     68  HD1 LEU B 441 .*.        "
        "ATOM     72  HD1 LEU B 441 .*.        "
  group "ATOM     69  HD2 LEU B 441 .*.        "
        "ATOM     73  HD2 LEU B 441 .*.        "
  group "ATOM     70  HD3 LEU B 441 .*.        "
        "ATOM     74  HD3 LEU B 441 .*.        "
  ... 1 remaining group not shown
""")
  #
  for segid in ["", "SEGI"]:
    oc = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM     68  HD1 LEU   441                                             %s
ATOM     72  HD1 LEU   441
""" % segid)).construct_hierarchy().overall_counts()
    if (segid != ""):
      oc.raise_duplicate_atom_labels_if_necessary()
    else:
      try: oc.raise_duplicate_atom_labels_if_necessary()
      except Sorry, e:
        assert not show_diff(str(e), '''\
number of groups of duplicate atom labels: 1
  total number of affected atoms:          2
  group "ATOM     68  HD1 LEU   441 .*.        "
        "ATOM     72  HD1 LEU   441 .*.        "''')
      else: raise Exception_expected
  #
  check("""\
HEADER    HYDROLASE                               19-JUL-05   2BWX
ATOM   2038  N   CYS A 249      68.746  44.381  71.143  0.70 21.04           N
ATOM   2039  CA  CYS A 249      68.957  43.022  71.606  0.70 21.28           C
ATOM   2040  C   CYS A 249      70.359  42.507  71.362  0.70 19.80           C
ATOM   2041  O   CYS A 249      71.055  42.917  70.439  0.70 19.80           O
ATOM   2042  CB ACYS A 249      67.945  42.064  70.987  0.40 24.99           C
ATOM   2044  SG ACYS A 249      66.261  42.472  71.389  0.40 27.94           S
ATOM   2043  CB BCYS A 249      67.928  42.101  70.948  0.30 23.34           C
ATOM   2045  SG BCYS A 249      67.977  40.404  71.507  0.30 26.46           S
HETATM 2046  N  CCSO A 249      68.746  44.381  71.143  0.30 21.04           N
HETATM 2047  CA CCSO A 249      68.957  43.022  71.606  0.30 21.28           C
HETATM 2048  CB CCSO A 249      67.945  42.064  70.987  0.30 24.99           C
HETATM 2049  SG CCSO A 249      66.261  42.472  71.389  0.30 27.94           S
HETATM 2050  C  CCSO A 249      70.359  42.507  71.362  0.30 19.80           C
HETATM 2051  O  CCSO A 249      71.055  42.917  70.439  0.30 19.80           O
HETATM 2052  OD CCSO A 249      66.275  42.201  72.870  0.30 23.67           O
""", """\
model id="" #chains=1
  chain id="A" #residue_groups=2
    resid=" 249 " #atom_groups=3
      altloc="" resname="CYS" #atoms=4
      altloc="A" resname="CYS" #atoms=2
      altloc="B" resname="CYS" #atoms=2
    resid=" 249 " #atom_groups=1  ### Info: same as previous resid ###
      altloc="C" resname="CSO" #atoms=7
""", """\
total number of:
  models:      1
  chains:      1
  alt. conf.:  3
  residues:    2
  atoms:      15
  anisou:      0
number of atom element+charge types: 4
histogram of atom element+charge frequency:
  " C  " 7
  " O  " 3
  " S  " 3
  " N  " 2
residue name classes:
  "common_amino_acid" 1
  "other"             1
number of chain ids: 1
histogram of chain id frequency:
  "A" 1
number of alt. conf. ids: 3
histogram of alt. conf. id frequency:
  "A" 1
  "B" 1
  "C" 1
residue alt. conf. situations:
  pure main conf.:     0
  pure alt. conf.:     1
  proper alt. conf.:   1
  improper alt. conf.: 0
chains with mix of proper and improper alt. conf.: 0
number of residue names: 2
histogram of residue name frequency:
  "CSO" 1    other
  "CYS" 1
### WARNING: consecutive residue_groups with same resid ###
number of consecutive residue groups with same resid: 1
  residue group:
    "ATOM   2038  N   CYS A 249 .*.     N  "
    ... 6 atoms not shown
    "ATOM   2045  SG BCYS A 249 .*.     S  "
  next residue group:
    "HETATM 2046  N  CCSO A 249 .*.     N  "
    ... 5 atoms not shown
    "HETATM 2052  OD CCSO A 249 .*.     O  "
""", level_id="atom_group")
  #
  check("""\
HEADER    OXIDOREDUCTASE                          17-SEP-97   1OHJ
HETATM 1552  C3  COP   188      11.436  28.065  13.009  1.00  8.51           C
HETATM 1553  C1  COP   188      13.269  26.907  13.759  1.00  8.86           C
HETATM 1582  O24 COP   188      13.931  34.344  22.009  1.00 20.08           O
HETATM 1583  O25 COP   188      13.443  32.717  20.451  1.00 20.18           O
HETATM 1608  O   HOH   188      20.354  30.097  11.632  1.00 21.33           O
HETATM 1569  C28ACOP   188      14.231  36.006  18.087  0.50 25.20           C
HETATM 1571  C29ACOP   188      13.126  36.948  17.945  0.50 26.88           C
HETATM 1604  O40ACOP   188      15.720  40.117  14.909  0.50 31.54           O
HETATM 1606  O41ACOP   188      15.816  42.243  14.385  0.50 31.73           O
HETATM 1570  C28BCOP   188      14.190  36.055  18.102  0.50 24.97           C
HETATM 1572  C29BCOP   188      13.133  37.048  18.009  0.50 26.45           C
HETATM 1605  O40BCOP   188      10.794  41.093  18.747  0.50 30.51           O
HETATM 1607  O41BCOP   188      12.838  40.007  19.337  0.50 30.37           O
""", """\
model id="" #chains=1
  chain id=" " #residue_groups=3
    resid=" 188 " #atom_groups=1
      altloc="" resname="COP" #atoms=4
    resid=" 188 " #atom_groups=1  ### Info: same as previous resid ###
      altloc="" resname="HOH" #atoms=1
    resid=" 188 " #atom_groups=2  ### Info: same as previous resid ###
      altloc="A" resname="COP" #atoms=4
      altloc="B" resname="COP" #atoms=4
""", """\
total number of:
  models:      1
  chains:      1
  alt. conf.:  2
  residues:    3
  atoms:      13
  anisou:      0
number of atom element+charge types: 2
histogram of atom element+charge frequency:
  " O  " 7
  " C  " 6
residue name classes:
  "other"        2
  "common_water" 1
number of chain ids: 1
histogram of chain id frequency:
  " " 1
number of alt. conf. ids: 2
histogram of alt. conf. id frequency:
  "A" 1
  "B" 1
residue alt. conf. situations:
  pure main conf.:     2
  pure alt. conf.:     1
  proper alt. conf.:   0
  improper alt. conf.: 0
chains with mix of proper and improper alt. conf.: 0
number of residue names: 2
histogram of residue name frequency:
  "COP" 2    other
  "HOH" 1    common water
### WARNING: consecutive residue_groups with same resid ###
number of consecutive residue groups with same resid: 2
  residue group:
    "HETATM 1552  C3  COP   188 .*.     C  "
    ... 2 atoms not shown
    "HETATM 1583  O25 COP   188 .*.     O  "
  next residue group:
    "HETATM 1608  O   HOH   188 .*.     O  "
  next residue group:
    "HETATM 1569  C28ACOP   188 .*.     C  "
    ... 6 atoms not shown
    "HETATM 1607  O41BCOP   188 .*.     O  "
""", level_id="atom_group")
  #
  check("""\
ATOM      1  N   R01     1I
ATOM      2  N   R02     1I
ATOM      3  N   R03     1I
ATOM      4  N   R04     1I
ATOM      5  N   R05     1I
""", None, """\
total number of:
  models:     1
  chains:     1
  alt. conf.: 0
  residues:   5
  atoms:      5
  anisou:     0
number of atom element+charge types: 1
histogram of atom element+charge frequency:
  "    " 5
residue name classes:
  "other" 5
number of chain ids: 1
histogram of chain id frequency:
  " " 1
number of alt. conf. ids: 0
number of residue names: 5
histogram of residue name frequency:
  "R01" 1    other
  "R02" 1    other
  "R03" 1    other
  "R04" 1    other
  "R05" 1    other
### WARNING: consecutive residue_groups with same resid ###
number of consecutive residue groups with same resid: 4
  residue group:
    "ATOM      1  N   R01     1I.*.        "
  next residue group:
    "ATOM      2  N   R02     1I.*.        "
  next residue group:
    "ATOM      3  N   R03     1I.*.        "
  next residue group:
    "ATOM      4  N   R04     1I.*.        "
  ------------------------------------------
  ... 1 remaining instance not shown
""")
  #
  root = pdb.hierarchy.root()
  assert not show_diff(root.as_str(), """\
### WARNING: empty hierarchy ###
""")
  assert not show_diff(root.overall_counts().as_str(), """\
total number of:
  models:     0
  chains:     0
  alt. conf.: 0
  residues:   0
  atoms:      0
  anisou:     0
number of atom element+charge types: 0
residue name classes: None
number of chain ids: 0
number of alt. conf. ids: 0
number of residue names: 0
""")
  model = pdb.hierarchy.model()
  root.append_model(model=model)
  assert not show_diff(root.as_str(), """\
model id="" #chains=0
  ### WARNING: empty model ###
""")
  assert not show_diff(root.overall_counts().as_str(), """\
total number of:
  ### WARNING: empty model ###
  models:     1 (1 empty)
  chains:     0
  alt. conf.: 0
  residues:   0
  atoms:      0
  anisou:     0
number of atom element+charge types: 0
residue name classes: None
number of chain ids: 0
number of alt. conf. ids: 0
number of residue names: 0
""")
  chain = pdb.hierarchy.chain()
  model.append_chain(chain=chain)
  assert not show_diff(root.as_str(), """\
model id="" #chains=1
  chain id="" #residue_groups=0
    ### WARNING: empty chain ###
""")
  assert not show_diff(root.overall_counts().as_str(), """\
total number of:
  models:     1
  ### WARNING: empty chain ###
  chains:     1 (1 empty)
  alt. conf.: 0
  residues:   0
  atoms:      0
  anisou:     0
number of atom element+charge types: 0
residue name classes: None
number of chain ids: 1
histogram of chain id frequency:
  "" 1
number of alt. conf. ids: 0
number of residue names: 0
""")
  residue_group = pdb.hierarchy.residue_group()
  chain.append_residue_group(residue_group=residue_group)
  assert not show_diff(root.as_str(), """\
model id="" #chains=1
  chain id="" #residue_groups=1
    resid="     " #atom_groups=0
      ### WARNING: empty residue_group ###
""")
  assert not show_diff(root.overall_counts().as_str(), """\
total number of:
  models:     1
  chains:     1
  alt. conf.: 0
  residues:   1
  atoms:      0
  anisou:     0
  ### WARNING: empty residue_group ###
  empty residue_groups: 1
number of atom element+charge types: 0
residue name classes: None
number of chain ids: 1
histogram of chain id frequency:
  "" 1
number of alt. conf. ids: 0
number of residue names: 0
""")
  atom_group = pdb.hierarchy.atom_group()
  residue_group.append_atom_group(atom_group=atom_group)
  assert not show_diff(root.as_str(), """\
model id="" #chains=1
  chain id="" #residue_groups=1
    resid="     " #atom_groups=1
      altloc="" resname="" #atoms=0
        ### WARNING: empty atom_group ###
""")
  oc = root.overall_counts()
  assert not show_diff(oc.as_str(), """\
total number of:
  models:     1
  chains:     1
  alt. conf.: 0
  residues:   1
  atoms:      0
  anisou:     0
  ### WARNING: empty atom_group ###
  empty atom_groups: 1
number of atom element+charge types: 0
residue name classes:
  "other" 1
number of chain ids: 1
histogram of chain id frequency:
  "" 1
number of alt. conf. ids: 0
number of residue names: 1
histogram of residue name frequency:
  "" 1    other
""")
  assert len(oc.errors()) == 0
  assert oc.warnings() == ['### WARNING: empty atom_group ###']
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM         N1 BR01     1
ATOM         N1  R01     1
ATOM         N1 AR02     2
ATOM         N1 BR02     2
ATOM         N2  R02     2
"""))
  oc = pdb_inp.construct_hierarchy().overall_counts()
  assert len(oc.warnings()) == 0
  oc.raise_duplicate_atom_labels_if_necessary()
  try: oc.raise_improper_alt_conf_if_necessary()
  except Sorry, e:
    assert not show_diff(str(e), '''\
residue with proper altloc
  "ATOM         N2  R02     2 .*.        "
  "ATOM         N1 AR02     2 .*.        "
  "ATOM         N1 BR02     2 .*.        "
residue with improper altloc
  "ATOM         N1  R01     1 .*.        "
  "ATOM         N1 BR01     1 .*.        "''')
  else: raise Exception_expected
  try: oc.raise_chains_with_mix_of_proper_and_improper_alt_conf_if_necessary()
  except Sorry, e:
    assert not show_diff(str(e), '''\
chains with mix of proper and improper alt. conf.: 1
  residue with proper altloc
    "ATOM         N2  R02     2 .*.        "
    "ATOM         N1 AR02     2 .*.        "
    "ATOM         N1 BR02     2 .*.        "
  residue with improper altloc
    "ATOM         N1  R01     1 .*.        "
    "ATOM         N1 BR01     1 .*.        "''')
  else: raise Exception_expected
  #
  sio = StringIO()
  summary = pdb.hierarchy.show_summary(
    out=sio,
    residue_groups_max_show=None,
    duplicate_atom_labels_max_show=None,
    pdb_string="""\
HEADER    HYDROLASE (SERINE PROTEINASE)           24-APR-89   1P04
HETATM 1410  N   B2I P   1      14.927  32.740  15.704  1.00 12.51           N
HETATM 1411  CA  B2I P   1      14.664  32.247  14.329  1.00 13.81           C
HETATM 1412  CB  B2I P   1      15.498  33.033  13.284  1.00 15.18           C
HETATM 1413  CG1 B2I P   1      15.052  34.514  13.266  1.00 15.84           C
HETATM 1414  CG2 B2I P   1      17.012  33.028  13.465  1.00 13.83           C
HETATM 1415  CD1 B2I P   1      13.720  34.655  12.538  1.00 16.03           C
HETATM 1416  B   B2I P   1      14.742  30.705  14.151  1.00 13.29           B
HETATM 1417  O1  B2I P   1      14.883  30.181  12.718  1.00 12.46           O
HETATM 1418  O2  B2I P   1      13.824  30.019  15.080  1.00 12.11           O
TER    1419      B2I P   1
HETATM 1420  S   SO4     1      32.724  31.060  29.009  0.99 34.51           S
HETATM 1421  O1  SO4     1      33.557  31.197  30.235  0.95 37.71           O
HETATM 1422  O2  SO4     1      31.648  32.120  28.923  0.86 36.50           O
HETATM 1423  O3  SO4     1      33.681  31.215  27.852  0.96 31.36           O
HETATM 1424  O4  SO4     1      32.067  29.709  28.958  1.00 44.58           O
""")
  assert summary.input.atoms().size() == 14
  assert summary.hierarchy.atoms().size() == 14
  oc = summary.overall_counts
  assert len(oc.consecutive_residue_groups_with_same_resid) == 0
  assert len(sio.getvalue().splitlines()) == 27
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM         CA  ASN     1
ATOM         C   ASN     1
ATOM         O   HOH     2
ATOM         H1  HOH     2
ATOM         H2  HOH     2
ATOM         CA  GLU     3
ATOM         C   GLU     3
ATOM         O   H2O     4
ATOM         O   OH2     5
ATOM         O   DOD     6
ATOM         P   U       7
ATOM         O   D2O     8
ATOM         O   OD2     9
ATOM         CA  ALA    10
ATOM         C   ALA    10
ATOM         O   TIP    11
ATOM         O   TIP    12
"""))
  oc = pdb_inp.construct_hierarchy().overall_counts()
  assert oc.resname_classes == {
    'common_water': 8, 'common_rna_dna': 1, 'common_amino_acid': 3}
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM         N  AASN     1
ATOM         CA AASN     1
ATOM         N  AGLY     1
ATOM         CA AGLY     1
ATOM         CA AASN     2
ATOM         CA AGLY     2
"""))
  oc = pdb_inp.construct_hierarchy().overall_counts()
  assert not show_diff(oc.as_str(), """\
total number of:
  models:     1
  chains:     1
  alt. conf.: 1
  residues:   2 (2 with mixed residue names)
  atoms:      6
  anisou:     0
number of atom element+charge types: 1
histogram of atom element+charge frequency:
  "    " 6
residue name classes:
  "common_amino_acid" 4
number of chain ids: 1
histogram of chain id frequency:
  " " 1
number of alt. conf. ids: 1
histogram of alt. conf. id frequency:
  "A" 1
residue alt. conf. situations:
  pure main conf.:     0
  pure alt. conf.:     2
  proper alt. conf.:   0
  improper alt. conf.: 0
chains with mix of proper and improper alt. conf.: 0
number of residue names: 2
histogram of residue name frequency:
  "ASN" 2
  "GLY" 2
### ERROR: residue group with multiple resnames using same altloc ###
residue groups with multiple resnames using same altloc: 2
  residue group:
    "ATOM         N  AASN     1 .*.        "
    ... 2 atoms not shown
    "ATOM         CA AGLY     1 .*.        "
  residue group:
    "ATOM         CA AASN     2 .*.        "
    "ATOM         CA AGLY     2 .*.        "
""")
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM         N  AASN     1
ATOM         N  AGLY     1
ATOM         CA AGLY     1
ATOM         CA AASN     2
ATOM         CA AGLY     2
"""))
  oc = pdb_inp.construct_hierarchy().overall_counts()
  try: oc \
   .raise_residue_groups_with_multiple_resnames_using_same_altloc_if_necessary(
      max_show=1)
  except Sorry, e:
    assert not show_diff(str(e), """\
residue groups with multiple resnames using same altloc: 2
  residue group:
    "ATOM         N  AASN     1 .*.        "
    "ATOM         N  AGLY     1 .*.        "
    "ATOM         CA AGLY     1 .*.        "
  ... 1 remaining instance not shown""")
  else: raise Exception_expected

def exercise_convenience_generators():
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
MODEL        1
ATOM      1  N  AR11 A   1
ATOM      2  O  AR11 A   1
ATOM      3  N  BR21 A   1
ATOM      4  O  BR21 A   1
ATOM      5  N  AR12 A   2
ATOM      6  O  AR12 A   2
ATOM      7  N  BR22 A   2
ATOM      8  O  BR22 A   2
TER
ATOM      9  N  AR11 B   1
ATOM     10  O  AR11 B   1
ATOM     11  N  BR21 B   1
ATOM     12  O  BR21 B   1
ATOM     13  N  AR12 B   2
ATOM     14  O  AR12 B   2
ATOM     15  N  BR22 B   2
ATOM     16  O  BR22 B   2
TER
ENDMDL
MODEL        2
ATOM      1  N  AR11 A   1
ATOM      2  O  AR11 A   1
ATOM      3  N  BR21 A   1
ATOM      4  O  BR21 A   1
ATOM      5  N  AR12 A   2
ATOM      6  O  AR12 A   2
ATOM      7  N  BR22 A   2
ATOM      8  O  BR22 A   2
TER
ATOM      9  N  AR11 B   1
ATOM     10  O  AR11 B   1
ATOM     11  N  BR21 B   1
ATOM     12  O  BR21 B   1
ATOM     13  N  AR12 B   2
ATOM     14  O  AR12 B   2
ATOM     15  N  BR22 B   2
ATOM     16  O  BR22 B   2
TER
ENDMDL
"""))
  obj = pdb_inp.construct_hierarchy(residue_group_post_processing=False)
  assert [model.id for model in obj.models()] == ["   1", "   2"]
  assert [chain.id for chain in obj.chains()] == ["A", "B"] * 2
  assert [rg.resid() for rg in obj.residue_groups()] == ["   1 ", "   2 "] * 4
  assert [ag.confid() for ag in obj.atom_groups()] \
      == ["AR11", "BR21", "AR12", "BR22"] * 4
  assert obj.atoms_size() == 32
  assert obj.atoms().size() == 32
  assert [atom.name for atom in obj.atoms()] == [" N  ", " O  "] * 16
  obj = obj.models()[0]
  assert [chain.id for chain in obj.chains()] == ["A", "B"]
  assert [rg.resid() for rg in obj.residue_groups()] == ["   1 ", "   2 "] * 2
  assert [ag.confid() for ag in obj.atom_groups()] \
      == ["AR11", "BR21", "AR12", "BR22"] * 2
  assert obj.atoms_size() == 16
  assert obj.atoms().size() == 16
  assert [atom.name for atom in obj.atoms()] == [" N  ", " O  "] * 8
  obj = obj.chains()[0]
  assert [rg.resid() for rg in obj.residue_groups()] == ["   1 ", "   2 "]
  assert [ag.confid() for ag in obj.atom_groups()] \
      == ["AR11", "BR21", "AR12", "BR22"]
  assert obj.atoms_size() == 8
  assert obj.atoms().size() == 8
  assert [atom.name for atom in obj.atoms()] == [" N  ", " O  "] * 4
  obj = obj.residue_groups()[0]
  assert [ag.confid() for ag in obj.atom_groups()] \
      == ["AR11", "BR21"]
  assert obj.atoms_size() == 4
  assert obj.atoms().size() == 4
  assert [atom.name for atom in obj.atoms()] == [" N  ", " O  "] * 2
  obj = obj.atom_groups()[0]
  assert obj.atoms_size() == 2
  assert obj.atoms().size() == 2
  assert [atom.name for atom in obj.atoms()] == [" N  ", " O  "]

def exercise_only():
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
MODEL        1
ATOM      2  N  ARES C   3I
ENDMDL
"""))
  hierarchy = pdb_inp.construct_hierarchy(residue_group_post_processing=False)
  obj = hierarchy
  assert obj.only_model().id == "   1"
  assert obj.only_chain().id == "C"
  assert obj.only_residue_group().resid() == "   3I"
  assert obj.only_atom_group().altloc == "A"
  assert obj.only_atom().name == " N  "
  obj = obj.only_model()
  assert obj.only_chain().id == "C"
  assert obj.only_residue_group().resid() == "   3I"
  assert obj.only_atom_group().resname == "RES"
  assert obj.only_atom().name == " N  "
  obj = obj.only_chain()
  assert obj.only_residue_group().resid() == "   3I"
  assert obj.only_atom_group().altloc == "A"
  assert obj.only_atom().name == " N  "
  obj = obj.only_residue_group()
  assert obj.only_atom_group().resname == "RES"
  assert obj.only_atom().name == " N  "
  obj = obj.only_atom_group()
  assert obj.only_atom().name == " N  "
  #
  obj = hierarchy
  assert obj.only_conformer().altloc == "A"
  assert obj.only_residue().resname == "RES"
  obj = obj.only_model()
  assert obj.only_conformer().altloc == "A"
  assert obj.only_residue().resname == "RES"
  obj = obj.only_chain()
  assert obj.only_conformer().altloc == "A"
  assert obj.only_residue().resname == "RES"
  obj = obj.only_conformer()
  assert obj.only_residue().resname == "RES"
  assert obj.only_atom().name == " N  "
  obj = obj.only_residue()
  assert obj.only_atom().name == " N  "
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM      1      R01 A
ATOM      2      R01 A
ATOM      3      R02 A
ATOM      4      R02 A
ATOM      5      R01 B
"""))
  hierarchy = pdb_inp.construct_hierarchy(residue_group_post_processing=False)
  sio = StringIO()
  for chain in hierarchy.chains():
    for residue in chain.residues():
      for atom in residue.atoms():
        print >> sio, atom.quote()
  assert not show_diff(sio.getvalue(), """\
"ATOM      1      R01 A     .*.        "
"ATOM      2      R01 A     .*.        "
"ATOM      3      R02 A     .*.        "
"ATOM      4      R02 A     .*.        "
"ATOM      5      R01 B     .*.        "
""")

exercise_merge_pdb_inp = pdb.input(
  source_info=None, lines=flex.split_lines("""\
ATOM   1716  N  ALEU   190      28.628   4.549  20.230  0.70  3.78           N
ATOM   1717  CA ALEU   190      27.606   5.007  19.274  0.70  3.71           C
ATOM   1718  CB ALEU   190      26.715   3.852  18.800  0.70  4.15           C
ATOM   1719  CG ALEU   190      25.758   4.277  17.672  0.70  4.34           C
ATOM   1829  N  BLEU   190      28.428   4.746  20.343  0.30  5.13           N
ATOM   1830  CA BLEU   190      27.378   5.229  19.418  0.30  4.89           C
ATOM   1831  CB BLEU   190      26.539   4.062  18.892  0.30  4.88           C
ATOM   1832  CG BLEU   190      25.427   4.359  17.878  0.30  5.95           C
ATOM   1724  N  ATHR   191      27.350   7.274  20.124  0.70  3.35           N
ATOM   1725  CA ATHR   191      26.814   8.243  21.048  0.70  3.27           C
ATOM   1726  CB ATHR   191      27.925   9.229  21.468  0.70  3.73           C
ATOM   1727  OG1ATHR   191      28.519   9.718  20.259  0.70  5.22           O
ATOM   1728  CG2ATHR   191      28.924   8.567  22.345  0.70  4.21           C
ATOM   1729  C  ATHR   191      25.587   8.983  20.559  0.70  3.53           C
ATOM   1730  O  ATHR   191      24.872   9.566  21.383  0.70  3.93           O
ATOM   1833  CD1BLEU   190      26.014   4.711  16.521  0.30  6.21           C
ATOM   1835  C  BLEU   190      26.506   6.219  20.135  0.30  4.99           C
ATOM   1836  O  BLEU   190      25.418   5.939  20.669  0.30  5.91           O
ATOM   1721  CD2ALEU   190      24.674   3.225  17.536  0.70  5.31           C
ATOM   1722  C  ALEU   190      26.781   6.055  20.023  0.70  3.36           C
ATOM   1723  O  ALEU   190      25.693   5.796  20.563  0.70  3.68           O
ATOM   8722  C  DLEU   190      26.781   6.055  20.023  0.70  3.36           C
ATOM   8723  O  DLEU   190      25.693   5.796  20.563  0.70  3.68           O
ATOM   9722  C  CLEU   190      26.781   6.055  20.023  0.70  3.36           C
ATOM   9723  O  CLEU   190      25.693   5.796  20.563  0.70  3.68           O
"""))

def exercise_merge_atom_groups():
  lines = []
  root = exercise_merge_pdb_inp.construct_hierarchy(
    residue_group_post_processing=False)
  chain = root.models()[0].chains()[0]
  residue_groups = chain.residue_groups()
  assert len(residue_groups) == 3
  for i_ag in [0,1]:
    primary_atom_group = residue_groups[0].atom_groups()[i_ag]
    assert (primary_atom_group.altloc, primary_atom_group.resname) \
        == ("AB"[i_ag], "LEU")
    secondary_atom_group = residue_groups[2].atom_groups()[1-i_ag]
    try:
      residue_groups[0].merge_atom_groups(
        primary=secondary_atom_group,
        secondary=primary_atom_group)
    except RuntimeError, e:
      assert not show_diff(str(e), """\
"primary" atom_group has a different or no parent\
 (this residue_group must be the parent).""")
    else: raise Exception_expected
    try:
      residue_groups[0].merge_atom_groups(
        primary=primary_atom_group,
        secondary=primary_atom_group)
    except RuntimeError, e:
      assert not show_diff(str(e), """\
"primary" and "secondary" atom_groups are identical.""")
    else: raise Exception_expected
    try:
      residue_groups[0].merge_atom_groups(
        primary=primary_atom_group,
        secondary=residue_groups[2].atom_groups()[i_ag])
    except RuntimeError, e:
      assert str(e).find("secondary.data->altloc == primary.data->altloc") > 0
    else: raise Exception_expected
    assert primary_atom_group.atoms_size() == 4
    assert secondary_atom_group.atoms_size() == 3
    residue_groups[0].merge_atom_groups(
      primary=primary_atom_group,
      secondary=secondary_atom_group)
    assert primary_atom_group.atoms_size() == 7
    assert secondary_atom_group.atoms_size() == 0
    sio = StringIO()
    for atom in primary_atom_group.atoms():
      print >> sio, atom.format_atom_record()
    assert not show_diff(sio.getvalue(), ["""\
ATOM   1716  N  ALEU   190      28.628   4.549  20.230  0.70  3.78           N
ATOM   1717  CA ALEU   190      27.606   5.007  19.274  0.70  3.71           C
ATOM   1718  CB ALEU   190      26.715   3.852  18.800  0.70  4.15           C
ATOM   1719  CG ALEU   190      25.758   4.277  17.672  0.70  4.34           C
ATOM   1721  CD2ALEU   190      24.674   3.225  17.536  0.70  5.31           C
ATOM   1722  C  ALEU   190      26.781   6.055  20.023  0.70  3.36           C
ATOM   1723  O  ALEU   190      25.693   5.796  20.563  0.70  3.68           O
""", """\
ATOM   1829  N  BLEU   190      28.428   4.746  20.343  0.30  5.13           N
ATOM   1830  CA BLEU   190      27.378   5.229  19.418  0.30  4.89           C
ATOM   1831  CB BLEU   190      26.539   4.062  18.892  0.30  4.88           C
ATOM   1832  CG BLEU   190      25.427   4.359  17.878  0.30  5.95           C
ATOM   1833  CD1BLEU   190      26.014   4.711  16.521  0.30  6.21           C
ATOM   1835  C  BLEU   190      26.506   6.219  20.135  0.30  4.99           C
ATOM   1836  O  BLEU   190      25.418   5.939  20.669  0.30  5.91           O
"""][i_ag])

def exercise_merge_residue_groups():
  root = exercise_merge_pdb_inp.construct_hierarchy(
    residue_group_post_processing=False)
  chain = root.models()[0].chains()[0]
  residue_groups = chain.residue_groups()
  assert len(residue_groups) == 3
  try:
    chain.merge_residue_groups(
      primary=residue_groups[0],
      secondary=residue_groups[1])
  except RuntimeError, e:
    assert str(e).find("secondary.data->resseq == primary.data->resseq") > 0
  else: raise Exception_expected
  assert residue_groups[0].atom_groups_size() == 2
  assert residue_groups[2].atom_groups_size() == 4
  assert residue_groups[2].parent().memory_id() == chain.memory_id()
  assert chain.residue_groups_size() == 3
  chain.merge_residue_groups(
    primary=residue_groups[0],
    secondary=residue_groups[2])
  assert residue_groups[0].atom_groups_size() == 4
  assert residue_groups[2].atom_groups_size() == 0
  assert residue_groups[2].parent() is None
  assert chain.residue_groups_size() == 2
  sio = StringIO()
  for atom_group in residue_groups[0].atom_groups():
    for atom in atom_group.atoms():
      print >> sio, atom.format_atom_record()
  assert not show_diff(sio.getvalue(), """\
ATOM   1716  N  ALEU   190      28.628   4.549  20.230  0.70  3.78           N
ATOM   1717  CA ALEU   190      27.606   5.007  19.274  0.70  3.71           C
ATOM   1718  CB ALEU   190      26.715   3.852  18.800  0.70  4.15           C
ATOM   1719  CG ALEU   190      25.758   4.277  17.672  0.70  4.34           C
ATOM   1721  CD2ALEU   190      24.674   3.225  17.536  0.70  5.31           C
ATOM   1722  C  ALEU   190      26.781   6.055  20.023  0.70  3.36           C
ATOM   1723  O  ALEU   190      25.693   5.796  20.563  0.70  3.68           O
ATOM   1829  N  BLEU   190      28.428   4.746  20.343  0.30  5.13           N
ATOM   1830  CA BLEU   190      27.378   5.229  19.418  0.30  4.89           C
ATOM   1831  CB BLEU   190      26.539   4.062  18.892  0.30  4.88           C
ATOM   1832  CG BLEU   190      25.427   4.359  17.878  0.30  5.95           C
ATOM   1833  CD1BLEU   190      26.014   4.711  16.521  0.30  6.21           C
ATOM   1835  C  BLEU   190      26.506   6.219  20.135  0.30  4.99           C
ATOM   1836  O  BLEU   190      25.418   5.939  20.669  0.30  5.91           O
ATOM   8722  C  DLEU   190      26.781   6.055  20.023  0.70  3.36           C
ATOM   8723  O  DLEU   190      25.693   5.796  20.563  0.70  3.68           O
ATOM   9722  C  CLEU   190      26.781   6.055  20.023  0.70  3.36           C
ATOM   9723  O  CLEU   190      25.693   5.796  20.563  0.70  3.68           O
""")
  for i_rg,j_rg,s in [(2,2,"primary"),(0,2,"secondary")]:
    try:
      chain.merge_residue_groups(
        primary=residue_groups[i_rg],
        secondary=residue_groups[j_rg])
    except RuntimeError, e:
      assert not show_diff(str(e), """\
"%s" residue_group has a different or no parent\
 (this chain must be the parent).""" % s)
    else: raise Exception_expected

def exercise_chain_merge_residue_groups(n_trials=30):
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
HEADER    HYDROLASE                               22-NOV-07   2VHL
HETATM 6362  O   HOH B2048      47.616  10.724 150.212  1.00 46.48           O
HETATM 6363  O  AHOH B2049      46.408  16.672 146.066  0.50 12.81           O
HETATM 6364  O   HOH B2050      29.343  12.806 185.898  1.00 35.57           O
HETATM 6365  O  BHOH B2049      43.786  12.615 147.734  0.50 28.43           O
HETATM 6366  O   HOH B2052      35.068  19.167 155.349  1.00 15.97           O
"""))
  for rgpp in [False, True]:
    chain = pdb_inp.construct_hierarchy(
      residue_group_post_processing=rgpp).only_chain()
    if (not rgpp):
      assert chain.residue_groups_size() == 5
      indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
      assert list(indices) == [1]
    assert chain.residue_groups_size() == 4
    indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
    assert indices.size() == 0
    del chain
  lines = flex.split_lines("""\
HETATM 6363  O  AHOH B2049
HETATM 6364  O  ZHOH B2050
HETATM 6365  O  BHOH B2049
HETATM 6366  O  YHOH B2052
HETATM 9365  O  CHOH B2049
HETATM 9367  O  XHOH B2052
""")
  pdb_inp = pdb.input(source_info=None, lines=lines)
  chain = pdb_inp.construct_hierarchy(
    residue_group_post_processing=False).only_chain()
  assert chain.residue_groups_size() == 6
  indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
  assert list(indices) == [0, 2]
  assert chain.residue_groups_size() == 3
  indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
  assert indices.size() == 0
  for i_trial in xrange(n_trials):
    pdb_inp = pdb.input(
      source_info=None,
      lines=lines.select(flex.random_permutation(size=lines.size())))
    chain = pdb_inp.construct_hierarchy(
      residue_group_post_processing=False).only_chain()
    indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
    assert indices.size() <= 2
    indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
    assert indices.size() == 0
    del chain
    chain = pdb_inp.construct_hierarchy().only_chain()
    indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
    assert indices.size() == 0
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
HEADER    SERINE PROTEASE                         10-NOV-95   1RTF
HETATM 2397  P   PO4     1      -7.520  25.376  38.369  1.00 39.37           P
HETATM 2398  O1  PO4     1      -6.610  24.262  38.967  1.00 40.00           O
HETATM 2399  O2  PO4     1      -6.901  25.919  37.049  1.00 41.07           O
HETATM 2400  O3  PO4     1      -8.894  24.741  38.097  1.00 45.09           O
HETATM 2401  O4  PO4     1      -7.722  26.556  39.350  1.00 42.48           O
HETATM 2402  C1  BEN     1      -6.921  31.206  33.893  1.00 23.35           C
HETATM 2403  C2  BEN     1      -8.189  30.836  34.344  1.00 23.15           C
HETATM 2404  C3  BEN     1      -8.335  29.863  35.342  1.00 20.74           C
HETATM 2405  C4  BEN     1      -7.206  29.254  35.893  1.00 19.45           C
HETATM 2406  C5  BEN     1      -5.932  29.618  35.445  1.00 20.83           C
HETATM 2407  C6  BEN     1      -5.794  30.589  34.450  1.00 20.99           C
HETATM 2408  C   BEN     1      -6.767  32.249  32.859  1.00 24.30           C
HETATM 2409  N1  BEN     1      -5.570  32.641  32.497  1.00 24.56           N
HETATM 2410  N2  BEN     1      -7.824  32.785  32.299  1.00 24.58           N
HETATM 2415  O   HOH     1       4.020  20.521  19.336  1.00 38.74           O
HETATM 2418  O   WAT     2      14.154  16.852  21.753  1.00 49.41           O
"""))
  chain = pdb_inp.construct_hierarchy(
    residue_group_post_processing=False).only_chain()
  assert chain.residue_groups_size() == 4
  assert [residue_group.resid() for residue_group in chain.residue_groups()] \
      == ["   1 ", "   1 ", "   1 ", "   2 "]
  for residue_group in chain.residue_groups():
    assert residue_group.atom_groups_size() == 1
    assert residue_group.atom_groups()[0].parent().memory_id() \
        == residue_group.memory_id()
  assert [residue_group.atom_groups()[0].resname
           for residue_group in chain.residue_groups()] \
      == ["PO4", "BEN", "HOH", "WAT"]
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
HETATM 2418  O   WAT     2
HETATM 2397  P   PO4     1
HETATM 2398  O1  PO4     1
HETATM 2402  C1  BEN     1
HETATM 2403  C2  BEN     1
HETATM 2404  C3  BEN     1
HETATM 2415  O   HOH     1
HETATM 9418  O   WAT     2
HETATM 9397  P   PO4     1
HETATM 9398  O1  PO4     1
HETATM 9402  C1  BEN     1
HETATM 9403  C2  BEN     1
HETATM 9404  C3  BEN     1
HETATM 9415  O   HOH     1
"""))
  chain = pdb_inp.construct_hierarchy(
    residue_group_post_processing=False).only_chain()
  assert chain.residue_groups_size() == 8
  for residue_group in chain.residue_groups():
    assert residue_group.atom_groups_size() == 1
    assert residue_group.atom_groups()[0].parent().memory_id() \
        == residue_group.memory_id()
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
HETATM 6362  O  CHOH B   1
HETATM 6363  O  AHOH B   1
HETATM 6364  O   HOH B   2
HETATM 6365  O   HOH B   1
"""))
  chain = pdb_inp.construct_hierarchy(
    residue_group_post_processing=False).only_chain()
  assert chain.residue_groups_size() == 3
  assert chain.residue_groups()[2].only_atom_group().altloc == " "
  for rg in chain.residue_groups():
    rg.edit_blank_altloc()
  assert chain.residue_groups()[2].only_atom_group().altloc == ""
  indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
  assert indices.size() == 0
  assert chain.residue_groups_size() == 3

def exercise_edit_blank_altloc(n_trials=30):
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM         N1
ATOM         N2
"""))
  for rgpp in [False, True]:
    residue_group = pdb_inp.construct_hierarchy(
      residue_group_post_processing=rgpp).only_residue_group()
    for i_proc in [0,1]:
      assert residue_group.edit_blank_altloc() == (1,0)
    del residue_group
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM         N1 A
ATOM         N2 B
"""))
  for rgpp in [False, True]:
    residue_group = pdb_inp.construct_hierarchy(
      residue_group_post_processing=rgpp).only_residue_group()
    rgc = residue_group.detached_copy()
    assert rgc.move_blank_altloc_atom_groups_to_front() == 0
    for i_proc in [0,1]:
      assert residue_group.edit_blank_altloc() == (0,0)
    del residue_group
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM         N1
ATOM         N2 B
"""))
  for rgpp in [False, True]:
    residue_group = pdb_inp.construct_hierarchy(
      residue_group_post_processing=rgpp).only_residue_group()
    if (not rgpp):
      atom_groups = residue_group.atom_groups()
      assert len(atom_groups) == 2
      assert atom_groups[0].altloc == " "
      assert atom_groups[1].altloc == "B"
    for i_proc in [0,1]:
      if (not rgpp or i_proc != 0):
        assert residue_group.edit_blank_altloc() == (1,0)
      atom_groups = residue_group.atom_groups()
      assert len(atom_groups) == 2
      assert atom_groups[0].altloc == ""
      assert atom_groups[1].altloc == "B"
    del atom_groups
    del residue_group
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM         N1
ATOM         N1 B
"""))
  residue_group = pdb_inp.construct_hierarchy(
    residue_group_post_processing=False).only_residue_group()
  atom_groups = residue_group.atom_groups()
  assert len(atom_groups) == 2
  assert atom_groups[0].altloc == " "
  assert atom_groups[1].altloc == "B"
  rgc = residue_group.detached_copy()
  assert rgc.move_blank_altloc_atom_groups_to_front() == 1
  for i_proc in [0,1]:
    assert residue_group.edit_blank_altloc() == (0,1)
    atom_groups = residue_group.atom_groups()
    assert len(atom_groups) == 2
    assert atom_groups[0].altloc == " "
    assert atom_groups[1].altloc == "B"
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM         N1 B
ATOM         N1
"""))
  for edit_chain in [False, True]:
    chain = pdb_inp.construct_hierarchy(
      residue_group_post_processing=False).only_chain()
    residue_group = chain.only_residue_group()
    atom_groups = residue_group.atom_groups()
    assert len(atom_groups) == 2
    assert atom_groups[0].altloc == "B"
    assert atom_groups[1].altloc == " "
    rgc = residue_group.detached_copy()
    assert rgc.move_blank_altloc_atom_groups_to_front() == 1
    for i_proc in [0,1]:
      if (not edit_chain):
        assert residue_group.edit_blank_altloc() == (0,1)
      else:
        for rg in chain.residue_groups():
          rg.edit_blank_altloc()
      atom_groups = residue_group.atom_groups()
      assert len(atom_groups) == 2
      assert atom_groups[0].altloc == " "
      assert atom_groups[1].altloc == "B"
    del atom_groups
    del residue_group
    del chain
  #
  lines = flex.split_lines("""\
ATOM         N1 B
ATOM         N1
ATOM         N2
ATOM         N3 B
ATOM         N3
""")
  for i_trial in xrange(n_trials):
    pdb_inp = pdb.input(source_info=None, lines=lines)
    residue_group = pdb_inp.construct_hierarchy(
      residue_group_post_processing=False).only_residue_group()
    atom_groups = residue_group.atom_groups()
    assert len(atom_groups) == 2
    if (i_trial == 0):
      assert atom_groups[0].altloc == "B"
      assert atom_groups[1].altloc == " "
    else:
      assert sorted([atom_group.altloc for atom_group in atom_groups]) \
          == [" ", "B"]
    for i_proc in [0,1]:
      assert residue_group.edit_blank_altloc() == (1,1)
      atom_groups = residue_group.atom_groups()
      assert len(atom_groups) == 3
      assert atom_groups[0].altloc == ""
      assert atom_groups[1].altloc == " "
      assert atom_groups[2].altloc == "B"
      lines = lines.select(flex.random_permutation(size=lines.size()))

def exercise_find_pure_altloc_ranges():
  c = pdb.hierarchy.chain()
  assert c.find_pure_altloc_ranges().size() == 0
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM            A        1
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert c.find_pure_altloc_ranges().size() == 0
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM            A        1
ATOM            B        1
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert c.find_pure_altloc_ranges().size() == 0
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM            A        1
ATOM            B        2
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert list(c.find_pure_altloc_ranges()) == [(0,2)]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM            A        1
BREAK
ATOM            B        2
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert c.find_pure_altloc_ranges().size() == 0
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM            A        1
ATOM            B        1
ATOM                     2
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert c.find_pure_altloc_ranges().size() == 0
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM            A        1
ATOM            B        2
ATOM            C        3
ATOM                     4
ATOM            E        5
ATOM            F        6
ATOM            G        6
ATOM            H        6
ATOM                     7
ATOM                     8
ATOM            I        9
ATOM            J       10
BREAK
ATOM            L       11
ATOM            M       12
ATOM         N1 N       13
ATOM         N2         13
ATOM            O       14
ATOM                    14
ATOM            P       15
ATOM                    15
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert list(c.find_pure_altloc_ranges()) \
      == [(0,3),(4,6),(8,10),(10,12),(13,15)]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
HEADER    CELL CYCLE                              13-SEP-05   2B05
HETATM10989  O   HOH    29     -66.337 -28.299 -26.997  1.00 40.05           O
HETATM10990  O  AHOH    32     -57.432 -22.290 -45.876  0.50  2.46           O
HETATM10991  O  BHOH    32     -59.435 -22.422 -45.055  0.50 17.09           O
HETATM10992  O   HOH    36     -56.803 -18.433 -29.790  1.00 43.00           O
HETATM10993  O   HOH    37     -51.860 -26.755 -35.092  1.00 35.90           O
HETATM10994  O  AHOH    39     -68.867 -23.643 -49.077  0.50 12.37           O
HETATM10995  O  BHOH    39     -69.097 -21.979 -50.740  0.50 21.64           O
HETATM10996  O   HOH    40     -65.221 -13.774 -33.183  1.00 36.14           O
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert c.find_pure_altloc_ranges().size() == 0
  #
  caa = "common_amino_acid"
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM            AALA     1
ATOM            BGLY     1
ATOM            ATYR     2
ATOM            BTHR     2
ATOM            AHOH     3
ATOM            BHOH     3
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert list(c.find_pure_altloc_ranges()) == [(0,3)]
  assert list(c.find_pure_altloc_ranges(common_residue_name_class_only=caa)) \
    == [(0,2)]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM            AHOH     3
ATOM            BHOH     3
ATOM            AALA     1
ATOM            BGLY     1
ATOM            ATYR     2
ATOM            BTHR     2
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert list(c.find_pure_altloc_ranges()) == [(0,3)]
  assert list(c.find_pure_altloc_ranges(common_residue_name_class_only=caa)) \
    == [(1,3)]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM            AALA     1
ATOM            BGLY     1
ATOM            AHOH     3
ATOM            BHOH     3
ATOM            ATYR     2
ATOM            BTHR     2
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert list(c.find_pure_altloc_ranges()) == [(0,3)]
  assert c.find_pure_altloc_ranges(common_residue_name_class_only=caa).size() \
    == 0
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM            AALA     1
ATOM            BGLY     1
ATOM            AHOH     3
ATOM            BHOH     3
ATOM            ATYR     2
ATOM            BTHR     2
ATOM            ASER     4
ATOM            BSER     4
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert list(c.find_pure_altloc_ranges()) == [(0,4)]
  assert list(c.find_pure_altloc_ranges(common_residue_name_class_only=caa)) \
    == [(2,4)]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM            ASER     4
ATOM            BSER     4
ATOM            AALA     1
ATOM            BGLY     1
ATOM            AHOH     3
ATOM            BHOH     3
ATOM            ATYR     2
ATOM            BTHR     2
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert list(c.find_pure_altloc_ranges()) == [(0,4)]
  assert list(c.find_pure_altloc_ranges(common_residue_name_class_only=caa)) \
    == [(0,2)]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM             MET     1
ATOM            ALEU     2
ATOM            ASER     3
ATOM            BSER     3
ATOM            AALA     4
ATOM            BGLY     4
ATOM            AHOH     5
ATOM            BHOH     5
ATOM            AHOH     6
ATOM            ATYR     7
ATOM            BTHR     7
ATOM            AGLY     8
ATOM            BGLY     8
"""))
  c = pdb_inp.construct_hierarchy().only_chain()
  assert list(c.find_pure_altloc_ranges()) == [(1,8)]
  assert list(c.find_pure_altloc_ranges(common_residue_name_class_only=caa)) \
      == [(1,4), (6,8)]

pdb_1nym_60 = """\
HEADER    HYDROLASE                               12-FEB-03   1NYM
ATOM     60  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM     63  CB ALYS A  32       9.197   8.686  12.246  0.29 14.71           C
ATOM     64  CB BLYS A  32       9.193   8.732  12.170  0.71 12.23           C
ATOM     74  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM     77  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
ATOM     82  CA ALYS A  34      14.979   4.895  12.608  0.60 15.67           C
ATOM     83  CA BLYS A  34      14.977   5.207  12.331  0.40 16.38           C
ATOM     88  CB ALYS A  34      15.128   3.896  11.472  0.60 12.11           C
ATOM     89  CB BLYS A  34      15.132   4.867  10.839  0.40 13.86           C
ATOM    100  CA AASP A  35      15.328   8.688  12.044  0.60 16.75           C
ATOM    101  CA BASP A  35      15.474   8.937  12.096  0.40 17.43           C
ATOM    106  CB AASP A  35      14.367   9.683  11.373  0.60 16.80           C
ATOM    107  CB BASP A  35      14.491   9.903  11.431  0.40 18.66           C
ATOM    115  CA  ALA A  36      14.978   9.140  15.828  1.00 12.65           C
ATOM    118  CB  ALA A  36      13.768   8.688  16.639  1.00 13.00           C
ATOM    121  CA AGLU A  37      17.683   6.514  16.549  0.59 12.26           C
ATOM    122  CA BGLU A  37      17.999   6.949  16.048  0.41 12.47           C
ATOM    127  CB AGLU A  37      17.694   5.030  16.164  0.59 11.08           C
ATOM    128  CB BGLU A  37      18.148   5.560  15.440  0.41 12.53           C
ATOM    139  CA AASP A  38      19.923   8.463  14.202  0.59 17.31           C
ATOM    140  CA BASP A  38      19.789   9.284  13.597  0.41 19.32           C
ATOM    145  CB AASP A  38      19.615   8.739  12.727  0.59 24.06           C
ATOM    146  CB BASP A  38      19.279   9.626  12.201  0.41 26.28           C
ATOM    155  CA AGLN A  39      19.069  11.941  15.596  0.62 19.31           C
ATOM    156  CA BGLN A  39      18.919  12.283  15.753  0.38 20.06           C
ATOM    161  CB AGLN A  39      17.681  12.586  15.630  0.62 21.92           C
ATOM    162  CB BGLN A  39      17.560  12.987  15.681  0.38 21.79           C
ATOM    172  CA  LEU A  40      19.526  10.711  19.160  1.00 13.99           C
ATOM    175  CB  LEU A  40      18.478   9.858  19.880  1.00 13.56           C
"""

pdb_2izq_220 = """\
HEADER    ANTIBIOTIC                              26-JUL-06   2IZQ
ATOM    220  N  ATRP A  11      20.498  12.832  34.558  0.50  6.03           N
ATOM    221  CA ATRP A  11      21.094  12.032  35.602  0.50  5.24           C
ATOM    222  C  ATRP A  11      22.601  12.088  35.532  0.50  6.49           C
ATOM    223  O  ATRP A  11      23.174  12.012  34.439  0.50  7.24           O
ATOM    224  CB ATRP A  11      20.690  10.588  35.288  0.50  6.15           C
ATOM    225  CG ATRP A  11      19.252  10.269  35.140  0.50  5.91           C
ATOM    226  CD1ATRP A  11      18.524  10.178  33.986  0.50  7.01           C
ATOM    227  CD2ATRP A  11      18.371   9.973  36.236  0.50  5.97           C
ATOM    228  NE1ATRP A  11      17.252   9.820  34.321  0.50  9.83           N
ATOM    229  CE2ATRP A  11      17.132   9.708  35.665  0.50  7.37           C
ATOM    230  CE3ATRP A  11      18.543   9.924  37.615  0.50  6.38           C
ATOM    231  CZ2ATRP A  11      16.033   9.388  36.460  0.50  8.25           C
ATOM    232  CZ3ATRP A  11      17.448   9.586  38.402  0.50  8.04           C
ATOM    233  CH2ATRP A  11      16.240   9.320  37.784  0.50  8.66           C
ATOM    234  H  ATRP A  11      20.540  12.567  33.741  0.50  7.24           H
ATOM    235  HA ATRP A  11      20.771  12.306  36.485  0.50  6.28           H
ATOM    236  HB2ATRP A  11      21.135  10.330  34.466  0.50  7.38           H
ATOM    237  HB3ATRP A  11      21.045  10.023  35.993  0.50  7.38           H
ATOM    244  N  CPHE A  11      20.226  13.044  34.556  0.15  6.35           N
ATOM    245  CA CPHE A  11      20.950  12.135  35.430  0.15  5.92           C
ATOM    246  C  CPHE A  11      22.448  12.425  35.436  0.15  6.32           C
ATOM    247  O  CPHE A  11      22.961  12.790  34.373  0.15  6.08           O
ATOM    248  CB CPHE A  11      20.768  10.667  34.994  0.15  6.01           C
ATOM    249  CG CPHE A  11      19.330  10.235  34.845  0.15  7.05           C
ATOM    250  CD1CPHE A  11      18.847   9.877  33.587  0.15  8.78           C
ATOM    251  CD2CPHE A  11      18.533  10.174  35.995  0.15  7.70           C
ATOM    252  CE1CPHE A  11      17.551   9.436  33.473  0.15 10.43           C
ATOM    253  CE2CPHE A  11      17.230   9.752  35.854  0.15  9.27           C
ATOM    254  CZ CPHE A  11      16.789   9.396  34.594  0.15 10.98           C
ATOM    255  N  BTYR A  11      20.553  12.751  34.549  0.35  5.21           N
ATOM    256  CA BTYR A  11      21.106  11.838  35.524  0.35  5.51           C
ATOM    257  C  BTYR A  11      22.625  11.920  35.572  0.35  5.42           C
ATOM    258  O  BTYR A  11      23.299  11.781  34.538  0.35  5.30           O
ATOM    259  CB BTYR A  11      20.694  10.354  35.327  0.35  5.65           C
ATOM    260  CG BTYR A  11      19.188  10.175  35.507  0.35  7.68           C
ATOM    261  CD1BTYR A  11      18.548  10.134  34.268  0.35  9.45           C
ATOM    262  HB2CPHE A  11      21.221  10.536  34.146  0.15  7.21           H
ATOM    263  CD2BTYR A  11      18.463  10.012  36.681  0.35  9.08           C
ATOM    264  HB3CPHE A  11      21.198  10.093  35.647  0.15  7.21           H
ATOM    265  CE1BTYR A  11      17.195   9.960  34.223  0.35 10.76           C
ATOM    266  HD1CPHE A  11      19.394   9.937  32.837  0.15 10.53           H
ATOM    267  CE2BTYR A  11      17.100   9.826  36.693  0.35 11.29           C
ATOM    268  HD2CPHE A  11      18.873  10.410  36.828  0.15  9.24           H
ATOM    269  CZ BTYR A  11      16.546   9.812  35.432  0.35 11.90           C
ATOM    270  HE1CPHE A  11      17.206   9.172  32.650  0.15 12.52           H
ATOM    271  OH BTYR A  11      15.178   9.650  35.313  0.35 19.29           O
ATOM    272  HE2CPHE A  11      16.661   9.708  36.588  0.15 11.13           H
ATOM    273  HZ CPHE A  11      15.908   9.110  34.509  0.15 13.18           H
ATOM    274  H  BTYR A  11      20.634  12.539  33.720  0.35  6.25           H
ATOM    275  HA BTYR A  11      20.773  12.116  36.402  0.35  6.61           H
HETATM  283  N   DLE A  12      23.179  12.148  36.720  1.00  7.16           N
HETATM  284  CA  DLE A  12      24.625  12.084  36.893  1.00  8.29           C
HETATM  285  CB ADLE A  12      25.039  10.717  37.621  0.65  9.02           C
HETATM  286  CB BDLE A  12      25.209  10.741  37.032  0.35 12.70           C
HETATM  287  CG ADLE A  12      24.658   9.548  36.780  0.65 12.06           C
HETATM  288  CG BDLE A  12      25.429   9.378  36.572  0.35 15.20           C
HETATM  289  CD1ADLE A  12      25.656   9.433  35.596  0.65 16.84           C
HETATM  290  CD1BDLE A  12      26.192   8.543  37.585  0.35 16.77           C
HETATM  291  CD2ADLE A  12      24.682   8.288  37.613  0.65 15.34           C
HETATM  292  CD2BDLE A  12      24.065   8.724  36.277  0.35 16.96           C
HETATM  293  C   DLE A  12      25.029  13.153  37.899  1.00  8.11           C
HETATM  294  O   DLE A  12      24.343  13.330  38.907  1.00 11.62           O
HETATM  295  H  ADLE A  12      22.682  12.228  37.418  0.50  8.60           H
HETATM  296  HA ADLE A  12      25.095  12.196  36.041  0.50  9.94           H
HETATM  297  HB1ADLE A  12      25.997  10.708  37.775  0.65 10.83           H
HETATM  298  HB1BDLE A  12      26.135  11.000  37.162  0.35 15.23           H
HETATM  299  HB2ADLE A  12      24.595  10.659  38.481  0.65 10.83           H
HETATM  300  HB2BDLE A  12      24.897  10.541  37.929  0.35 15.23           H
HETATM  301  HG ADLE A  12      23.753   9.685  36.429  0.65 14.47           H
HETATM  302  HG BDLE A  12      25.946   9.409  35.740  0.35 18.24           H
"""

def exercise_occupancy_groups_simple():
  def atom_serials(atoms, list_of_occ_groups):
    result = []
    for groups in list_of_occ_groups:
      group_names = []
      for group in groups:
        group_names.append([int(atoms[i].serial) for i in group])
      result.append(group_names)
    return result
  #
  def grouped_serials(
        pdb_inp,
        common_residue_name_class_only="common_amino_acid"):
    hierarchy = pdb_inp.construct_hierarchy()
    atoms = hierarchy.atoms()
    sentinel = atoms.reset_tmp_for_occupancy_groups_simple()
    chain = hierarchy.only_chain()
    return atom_serials(atoms, chain.occupancy_groups_simple(
      common_residue_name_class_only=common_residue_name_class_only))
  #
  for altloc_o2_a in ["A", " "]:
    pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM      0  S   SO4
ATOM      1  O1  SO4
ATOM      2  O2 %sSO4
ATOM      3  O2 BSO4
ATOM      4  O3  SO4
ATOM      5  O4  SO4
""" % altloc_o2_a))
    assert grouped_serials(pdb_inp) == [[[2], [3]]]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM      6  S  ASO4     1       1.302   1.419   1.560  0.70 10.00           S
ATOM      7  O1 ASO4     1       1.497   1.295   0.118  0.70 10.00           O
ATOM      8  O2 ASO4     1       1.098   0.095   2.140  0.70 10.00           O
ATOM      9  O3 ASO4     1       2.481   2.037   2.159  0.70 10.00           O
ATOM     10  O4 ASO4     1       0.131   2.251   1.823  0.70 10.00           O
ATOM     11  S  BSO4     1       3.302   3.419   3.560  0.30 10.00           S
ATOM     12  O1 BSO4     1       3.497   3.295   2.118  0.30 10.00           O
ATOM     13  O2 BSO4     1       3.098   2.095   4.140  0.30 10.00           O
ATOM     14  O3 BSO4     1       4.481   4.037   4.159  0.30 10.00           O
ATOM     15  O4 BSO4     1       2.131   4.251   3.823  0.30 10.00           O
"""))
  assert grouped_serials(pdb_inp) == [[[6,7,8,9,10], [11,12,13,14,15]]]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM     16  O  AHOH     2       5.131   5.251   5.823  0.60 10.00           O
ATOM     17  O  BHOH     2       6.131   6.251   6.823  0.40 10.00           O
"""))
  assert grouped_serials(pdb_inp) == [[[16], [17]]]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM     18  O   HOH     3       1.132   5.963   7.065  1.00 15.00           O
ATOM     19  H1  HOH     3       1.160   5.211   6.437  1.00 15.00           H
ATOM     20  H2  HOH     3       1.122   5.579   7.967  1.00 15.00           H
"""))
  assert grouped_serials(pdb_inp) == []
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM     21  O   HOH     4       6.131   7.251   5.000  0.50 15.00           O
"""))
  assert grouped_serials(pdb_inp) == [[[21]]]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM     22  O   HOH     5       0.131   7.251   5.000  0.00 15.00           O
"""))
  assert grouped_serials(pdb_inp) == []
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM     23  S   SO4     6       6.302   6.419   1.560  0.50 10.00           S
ATOM     24  O1 ASO4     6       6.497   6.295   0.118  0.60 10.00           O
ATOM     25  O2 ASO4     6       6.098   5.095   2.140  0.60 10.00           O
ATOM     26  O3 ASO4     6       7.481   7.037   2.159  0.60 10.00           O
ATOM     27  O4 ASO4     6       5.131   7.251   1.823  0.60 10.00           O
ATOM     28  O1 BSO4     6       8.497   8.295   2.118  0.40 10.00           O
ATOM     29  O2 BSO4     6       8.098   7.095   4.140  0.40 10.00           O
ATOM     30  O3 BSO4     6       9.481   9.037   4.159  0.40 10.00           O
ATOM     31  O4 BSO4     6       7.131   9.251   3.823  0.40 10.00           O
"""))
  assert grouped_serials(pdb_inp) == [[[23]], [[24,25,26,27], [28,29,30,31]]]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM      1  O  AHOH     1                                                   O
ATOM      2  O  BHOH     1                                                   O
ATOM      3  H1 AHOH     1                                                   H
ATOM      4  H1 BHOH     1                                                   H
ATOM      5  H2 AHOH     1                                                   H
ATOM      6  H2 BHOH     1                                                   H
"""))
  assert grouped_serials(pdb_inp) == [[[1],[2]]]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM      1  O   HOH     1                              0.60                 O
ATOM      2  H1  HOH     1                              0.60                 H
ATOM      3  H2  HOH     1                              0.60                 H
"""))
  assert grouped_serials(pdb_inp) == [[[1]]]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines(pdb_1nym_60))
  assert grouped_serials(pdb_inp) == [
    [[63],[64]],
    [[82,88,100,106],[83,89,101,107]],
    [[121,127,139,145,155,161],[122,128,140,146,156,162]]]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM     82  CA AMOD A  34      14.979   4.895  12.608  0.60 15.67           C
ATOM     83  CA BMOD A  34      14.977   5.207  12.331  0.40 16.38           C
ATOM     88  CB AMOD A  34      15.128   3.896  11.472  0.60 12.11           C
ATOM     89  CB BMOD A  34      15.132   4.867  10.839  0.40 13.86           C
ATOM    100  CA AASP A  35      15.328   8.688  12.044  0.60 16.75           C
ATOM    101  CA BASP A  35      15.474   8.937  12.096  0.40 17.43           C
ATOM    106  CB AASP A  35      14.367   9.683  11.373  0.60 16.80           C
ATOM    107  CB BASP A  35      14.491   9.903  11.431  0.40 18.66           C
"""))
  assert grouped_serials(pdb_inp) == [
    [[82,88],[83,89]],
    [[100,106],[101,107]]]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines(pdb_2izq_220))
  assert grouped_serials(pdb_inp) == [
    [[220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233],
     [244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254],
     [255, 256, 257, 258, 259, 260, 261, 263, 265, 267, 269, 271]],
    [[285, 287, 289, 291], [286, 288, 290, 292]]]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM    221  CA ATRP A  11      21.094  12.032  35.602  0.50  5.24           C
ATOM    224  CB ATRP A  11      20.690  10.588  35.288  0.50  6.15           C
ATOM    245  CA CPHE A  11      20.950  12.135  35.430  0.15  5.92           C
ATOM    248  CB CPHE A  11      20.768  10.667  34.994  0.15  6.01           C
ATOM    256  CA BTYR A  11      21.106  11.838  35.524  0.35  5.51           C
ATOM    259  CB BTYR A  11      20.694  10.354  35.327  0.35  5.65           C
HETATM  285  CB ADLE A  12      25.039  10.717  37.621  0.65  9.02           C
HETATM  286  CB BDLE A  12      25.209  10.741  37.032  0.35 12.70           C
HETATM  287  CG ADLE A  12      24.658   9.548  36.780  0.65 12.06           C
HETATM  288  CG BDLE A  12      25.429   9.378  36.572  0.35 15.20           C
"""))
  assert grouped_serials(pdb_inp) == [
    [[221, 224], [245, 248], [256, 259]], [[285, 287], [286, 288]]]
  assert grouped_serials(pdb_inp, common_residue_name_class_only=None) == [
    [[221, 224, 285, 287], [245, 248], [256, 259, 286, 288]]]
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM      1  O  AHOH A   1                                                   O
ATOM      2  O  BHOH A   1                                                   O
ATOM      3  O  AHOH B   1                                                   O
ATOM      4  O  BHOH B   1                                                   O
"""))
  hierarchy = pdb_inp.construct_hierarchy()
  list_of_groups = hierarchy.occupancy_groups_simple(
    common_residue_name_class_only="common_amino_acid")
  assert list_of_groups == [[[0], [1]], [[2], [3]]]

def conformers_as_str(conformers):
  s = StringIO()
  for cf in conformers:
    print >> s, "conformer:", show_string(cf.altloc)
    for rd in cf.residues():
      assert rd.resseq_as_int() == pdb.hy36decode(width=4, s=rd.resseq)
      print >> s, "  residue:", \
        show_string(rd.resname), \
        show_string(rd.resseq), \
        show_string(rd.icode), \
        int(rd.link_to_previous), \
        int(rd.is_pure_main_conf)
      for atom in rd.atoms():
        print >> s, "    atom:", show_string(atom.name)
  return s.getvalue()

def exercise_conformers():
  assert len(pdb.hierarchy.chain().conformers()) == 0
  assert len(pdb.hierarchy.residue_group().conformers()) == 0
  #
  def check(pdb_string, expected):
    pdb_inp = pdb.input(source_info=None, lines=flex.split_lines(pdb_string))
    chain = pdb_inp.construct_hierarchy().only_chain()
    conformers = chain.conformers()
    s = conformers_as_str(conformers)
    if (len(expected) == 0):
      sys.stdout.write(s)
    else:
      assert not show_diff(s, expected)
    #
    for rg in chain.residue_groups():
      for cf in rg.conformers():
        assert cf.parent(optional=False).memory_id() == chain.memory_id()
        assert cf.residues_size() == 1
      rgc = rg.detached_copy()
      for cf in rgc.conformers():
        assert cf.parent() is None
        assert cf.residues_size() == 1
        try:
          cf.parent(optional=False)
        except RuntimeError, e:
          assert not show_diff(str(e), "conformer has no parent chain")
        else: raise Exception_expected
      #
      if (chain.residue_groups_size() == 1):
        conformers = rg.conformers()
        s = conformers_as_str(conformers)
        if (len(expected) != 0):
          assert not show_diff(s, expected)
  #
  check("""\
ATOM         N   RES     1I
""", """\
conformer: ""
  residue: "RES" "   1" "I" 0 1
    atom: " N  "
""")
  #
  check("""\
ATOM         N  ARES     1I
""", """\
conformer: "A"
  residue: "RES" "   1" "I" 0 0
    atom: " N  "
""")
  #
  check("""\
ATOM         N1  RES     1I
ATOM         N2 ARES     1I
""", """\
conformer: "A"
  residue: "RES" "   1" "I" 0 0
    atom: " N1 "
    atom: " N2 "
""")
  #
  for altloc_o2_a in ["A", " "]:
    check("""\
ATOM         S   SO4     1I
ATOM         O1 %sSO4     1I
ATOM         O1 BSO4     1I
""" % altloc_o2_a, """\
conformer: "%s"
  residue: "SO4" "   1" "I" 0 0
    atom: " S  "
    atom: " O1 "
conformer: "B"
  residue: "SO4" "   1" "I" 0 0
    atom: " S  "
    atom: " O1 "
""" % altloc_o2_a)
  #
  check("""\
ATOM         S  ASO4     1
ATOM         O1 ASO4     1
ATOM         O2 ASO4     1
ATOM         S  BSO4     1
ATOM         O1 BSO4     1
ATOM         O2 BSO4     1
""", """\
conformer: "A"
  residue: "SO4" "   1" " " 0 0
    atom: " S  "
    atom: " O1 "
    atom: " O2 "
conformer: "B"
  residue: "SO4" "   1" " " 0 0
    atom: " S  "
    atom: " O1 "
    atom: " O2 "
""")
  #
  check("""\
ATOM         S  ASO4     1
ATOM         O1 ASO4     1
ATOM         O2 ASO4     1
ATOM         S  BSO4     1
ATOM         O1 BSO4     1
ATOM         O2 BSO4     1
ATOM         O   HOH     2
""", """\
conformer: "A"
  residue: "SO4" "   1" " " 0 0
    atom: " S  "
    atom: " O1 "
    atom: " O2 "
  residue: "HOH" "   2" " " 1 1
    atom: " O  "
conformer: "B"
  residue: "SO4" "   1" " " 0 0
    atom: " S  "
    atom: " O1 "
    atom: " O2 "
  residue: "HOH" "   2" " " 1 1
    atom: " O  "
""")
  #
  check("""\
ATOM         S   SO4     6
ATOM         O1 ASO4     6
ATOM         O2 ASO4     6
ATOM         O3 BSO4     6
ATOM         O4 BSO4     6
""", """\
conformer: "A"
  residue: "SO4" "   6" " " 0 0
    atom: " S  "
    atom: " O1 "
    atom: " O2 "
conformer: "B"
  residue: "SO4" "   6" " " 0 0
    atom: " S  "
    atom: " O3 "
    atom: " O4 "
""")
  #
  check("""\
ATOM         N1  R01     1I
ATOM         N2  R01     1I
ATOM         N1  R02     1I
ATOM         N2  R02     1I
""", """\
conformer: ""
  residue: "R01" "   1" "I" 0 1
    atom: " N1 "
    atom: " N2 "
  residue: "R02" "   1" "I" 1 1
    atom: " N1 "
    atom: " N2 "
""")
  #
  check("""\
ATOM         N1 AR01     1I
ATOM         N2  R01     1I
ATOM         N1  R02     1I
ATOM         N2  R02     1I
""", """\
conformer: "A"
  residue: "R01" "   1" "I" 0 0
    atom: " N2 "
    atom: " N1 "
  residue: "R02" "   1" "I" 1 1
    atom: " N1 "
    atom: " N2 "
""")
  #
  check("""\
ATOM         N1 AR01     1I
ATOM         N2  R01     1I
ATOM         N1 AR02     1I
ATOM         N2  R02     1I
""", """\
conformer: "A"
  residue: "R01" "   1" "I" 0 0
    atom: " N2 "
    atom: " N1 "
  residue: "R02" "   1" "I" 1 0
    atom: " N2 "
    atom: " N1 "
""")
  #
  check("""\
ATOM         N1 AR01     1I
ATOM         N2  R01     1I
ATOM         N1 BR02     1I
ATOM         N2  R02     1I
""", """\
conformer: "A"
  residue: "R01" "   1" "I" 0 0
    atom: " N2 "
    atom: " N1 "
  residue: "R02" "   1" "I" 1 1
    atom: " N2 "
conformer: "B"
  residue: "R01" "   1" "I" 0 1
    atom: " N2 "
  residue: "R02" "   1" "I" 1 0
    atom: " N2 "
    atom: " N1 "
""")
  #
  check("""\
ATOM         N1 AR01     1I
ATOM         N2 BR01     1I
ATOM         N1  R02     1I
ATOM         N2  R02     1I
""", """\
conformer: "A"
  residue: "R01" "   1" "I" 0 0
    atom: " N1 "
  residue: "R02" "   1" "I" 1 1
    atom: " N1 "
    atom: " N2 "
conformer: "B"
  residue: "R01" "   1" "I" 0 0
    atom: " N2 "
  residue: "R02" "   1" "I" 1 1
    atom: " N1 "
    atom: " N2 "
""")
  #
  check("""\
ATOM         N1 AR01     1I
ATOM         N2 BR01     1I
ATOM         N1 CR02     1I
ATOM         N2  R02     1I
""", """\
conformer: "A"
  residue: "R01" "   1" "I" 0 0
    atom: " N1 "
  residue: "R02" "   1" "I" 1 1
    atom: " N2 "
conformer: "B"
  residue: "R01" "   1" "I" 0 0
    atom: " N2 "
  residue: "R02" "   1" "I" 1 1
    atom: " N2 "
conformer: "C"
  residue: "R02" "   1" "I" 1 0
    atom: " N2 "
    atom: " N1 "
""")
  #
  check("""\
ATOM         N1 AR01     1I
ATOM         N2 BR01     1I
ATOM         N1 CR02     1I
ATOM         N2 DR02     1I
""", """\
conformer: "A"
  residue: "R01" "   1" "I" 0 0
    atom: " N1 "
conformer: "B"
  residue: "R01" "   1" "I" 0 0
    atom: " N2 "
conformer: "C"
  residue: "R02" "   1" "I" 0 0
    atom: " N1 "
conformer: "D"
  residue: "R02" "   1" "I" 0 0
    atom: " N2 "
""")
  #
  check(pdb_1nym_60, """\
conformer: "A"
  residue: "LYS" "  32" " " 0 0
    atom: " CA "
    atom: " CB "
  residue: "VAL" "  33" " " 1 1
    atom: " CA "
    atom: " CB "
  residue: "LYS" "  34" " " 1 0
    atom: " CA "
    atom: " CB "
  residue: "ASP" "  35" " " 1 0
    atom: " CA "
    atom: " CB "
  residue: "ALA" "  36" " " 1 1
    atom: " CA "
    atom: " CB "
  residue: "GLU" "  37" " " 1 0
    atom: " CA "
    atom: " CB "
  residue: "ASP" "  38" " " 1 0
    atom: " CA "
    atom: " CB "
  residue: "GLN" "  39" " " 1 0
    atom: " CA "
    atom: " CB "
  residue: "LEU" "  40" " " 1 1
    atom: " CA "
    atom: " CB "
conformer: "B"
  residue: "LYS" "  32" " " 0 0
    atom: " CA "
    atom: " CB "
  residue: "VAL" "  33" " " 1 1
    atom: " CA "
    atom: " CB "
  residue: "LYS" "  34" " " 1 0
    atom: " CA "
    atom: " CB "
  residue: "ASP" "  35" " " 1 0
    atom: " CA "
    atom: " CB "
  residue: "ALA" "  36" " " 1 1
    atom: " CA "
    atom: " CB "
  residue: "GLU" "  37" " " 1 0
    atom: " CA "
    atom: " CB "
  residue: "ASP" "  38" " " 1 0
    atom: " CA "
    atom: " CB "
  residue: "GLN" "  39" " " 1 0
    atom: " CA "
    atom: " CB "
  residue: "LEU" "  40" " " 1 1
    atom: " CA "
    atom: " CB "
""")
  #
  check(pdb_2izq_220, """\
conformer: "A"
  residue: "TRP" "  11" " " 0 0
    atom: " N  "
    atom: " CA "
    atom: " C  "
    atom: " O  "
    atom: " CB "
    atom: " CG "
    atom: " CD1"
    atom: " CD2"
    atom: " NE1"
    atom: " CE2"
    atom: " CE3"
    atom: " CZ2"
    atom: " CZ3"
    atom: " CH2"
    atom: " H  "
    atom: " HA "
    atom: " HB2"
    atom: " HB3"
  residue: "DLE" "  12" " " 1 0
    atom: " N  "
    atom: " CA "
    atom: " C  "
    atom: " O  "
    atom: " CB "
    atom: " CG "
    atom: " CD1"
    atom: " CD2"
    atom: " H  "
    atom: " HA "
    atom: " HB1"
    atom: " HB2"
    atom: " HG "
conformer: "C"
  residue: "PHE" "  11" " " 0 0
    atom: " N  "
    atom: " CA "
    atom: " C  "
    atom: " O  "
    atom: " CB "
    atom: " CG "
    atom: " CD1"
    atom: " CD2"
    atom: " CE1"
    atom: " CE2"
    atom: " CZ "
    atom: " HB2"
    atom: " HB3"
    atom: " HD1"
    atom: " HD2"
    atom: " HE1"
    atom: " HE2"
    atom: " HZ "
  residue: "DLE" "  12" " " 1 1
    atom: " N  "
    atom: " CA "
    atom: " C  "
    atom: " O  "
conformer: "B"
  residue: "TYR" "  11" " " 0 0
    atom: " N  "
    atom: " CA "
    atom: " C  "
    atom: " O  "
    atom: " CB "
    atom: " CG "
    atom: " CD1"
    atom: " CD2"
    atom: " CE1"
    atom: " CE2"
    atom: " CZ "
    atom: " OH "
    atom: " H  "
    atom: " HA "
  residue: "DLE" "  12" " " 1 0
    atom: " N  "
    atom: " CA "
    atom: " C  "
    atom: " O  "
    atom: " CB "
    atom: " CG "
    atom: " CD1"
    atom: " CD2"
    atom: " HB1"
    atom: " HB2"
    atom: " HG "
""")
  #
  check("""\
HEADER    HORMONE                                 01-MAY-98   1ZEH
HETATM  878  C1 ACRS     5      12.880  14.021   1.197  0.50 33.23           C
HETATM  879  C1 BCRS     5      12.880  14.007   1.210  0.50 34.27           C
HETATM  880  C2 ACRS     5      12.755  14.853   0.093  0.50 33.88           C
HETATM  881  C2 BCRS     5      13.935  13.115   1.278  0.50 34.25           C
HETATM  882  C3 ACRS     5      13.668  14.754  -0.945  0.50 33.82           C
HETATM  883  C3 BCRS     5      14.848  13.014   0.238  0.50 34.30           C
HETATM  884  C4 ACRS     5      14.707  13.834  -0.888  0.50 33.46           C
HETATM  885  C4 BCRS     5      14.695  13.821  -0.884  0.50 34.40           C
HETATM  886  C5 ACRS     5      14.835  13.001   0.219  0.50 33.30           C
HETATM  887  C5 BCRS     5      13.635  14.719  -0.957  0.50 34.78           C
HETATM  888  C6 ACRS     5      13.916  13.105   1.252  0.50 33.26           C
HETATM  889  C6 BCRS     5      12.731  14.813   0.090  0.50 34.86           C
HETATM  890  C7 ACRS     5      13.552  15.660  -2.169  0.50 33.90           C
HETATM  891  C7 BCRS     5      16.001  12.014   0.353  0.50 34.77           C
HETATM  892  O1 ACRS     5      11.973  14.116   2.233  0.50 34.24           O
HETATM  893  O1 BCRS     5      11.973  14.107   2.248  0.50 35.28           O
HETATM  894  O   HOH     5      -0.924  19.122  -8.629  1.00 11.73           O
HETATM  895  O   HOH     6     -19.752  11.918   3.524  1.00 13.44           O
""", """\
conformer: "A"
  residue: "CRS" "   5" " " 0 0
    atom: " C1 "
    atom: " C2 "
    atom: " C3 "
    atom: " C4 "
    atom: " C5 "
    atom: " C6 "
    atom: " C7 "
    atom: " O1 "
  residue: "HOH" "   5" " " 1 1
    atom: " O  "
  residue: "HOH" "   6" " " 1 1
    atom: " O  "
conformer: "B"
  residue: "CRS" "   5" " " 0 0
    atom: " C1 "
    atom: " C2 "
    atom: " C3 "
    atom: " C4 "
    atom: " C5 "
    atom: " C6 "
    atom: " C7 "
    atom: " O1 "
  residue: "HOH" "   5" " " 1 1
    atom: " O  "
  residue: "HOH" "   6" " " 1 1
    atom: " O  "
""")
  check("""\
HEADER    HYDROLASE                               22-NOV-07   2VHL
HETATM 6362  O   HOH B2048      47.616  10.724 150.212  1.00 46.48           O
HETATM 6363  O  AHOH B2049      46.408  16.672 146.066  0.50 12.81           O
HETATM 6364  O   HOH B2050      29.343  12.806 185.898  1.00 35.57           O
HETATM 6365  O  BHOH B2049      43.786  12.615 147.734  0.50 28.43           O
HETATM 6366  O   HOH B2052      35.068  19.167 155.349  1.00 15.97           O
""", """\
conformer: "A"
  residue: "HOH" "2048" " " 0 1
    atom: " O  "
  residue: "HOH" "2049" " " 1 0
    atom: " O  "
  residue: "HOH" "2050" " " 1 1
    atom: " O  "
  residue: "HOH" "2052" " " 1 1
    atom: " O  "
conformer: "B"
  residue: "HOH" "2048" " " 0 1
    atom: " O  "
  residue: "HOH" "2049" " " 1 0
    atom: " O  "
  residue: "HOH" "2050" " " 1 1
    atom: " O  "
  residue: "HOH" "2052" " " 1 1
    atom: " O  "
""")

def exercise_residue():
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM      1  N   MET A   1
ATOM      2  CA  MET A   1
"""))
  h = pdb_inp.construct_hierarchy()
  residue = h.only_residue()
  assert residue.resid() == "   1 "
  assert residue.find_atom_by(name=None) is None
  assert residue.find_atom_by(name=" N  ").name == " N  "
  assert residue.find_atom_by(name="N   ") is None
  assert residue.find_atom_by(name=" CA ").name == " CA "
  try:
    residue.parent(optional=False)
  except RuntimeError, e:
    assert not show_diff(str(e), "residue has no parent conformer")
  else: raise Exception_expected

def exercise_is_identical_hierarchy():
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
MODEL        0
ATOM      1  N   MET A   1
ATOM      2  CA  MET A   1
ATOM      3  N   GLY A   2
ATOM      4  CA  GLY A   2
ENDMDL
MODEL        1
ATOM      1  N   MET A   1
ATOM      2  CA  MET A   1
ATOM      3  N   GLY A   2
ATOM      4  CA  GLY A   2
ENDMDL
MODEL        2
ATOM      1  N   MET A   1
HETATM    2  CA  MET A   1
ATOM      3  N   GLY A   2
ATOM      4  CA  GLY A   2
ENDMDL
MODEL        3
ATOM      1  N   MET A   1
ATOM     12  CA  MET A   1
ATOM      3  N   GLY A   2
ATOM      4  CA  GLY A   2
ENDMDL
MODEL        4
ATOM      1  N   MET A   1
ATOM      2  CA  MET A   1
ATOM      3  NX  GLY A   2
ATOM      4  CA  GLY A   2
ENDMDL
MODEL        5
ATOM      1  N   MET A   1                                                   E
ATOM      2  CA  MET A   1
ATOM      3  N   GLY A   2
ATOM      4  CA  GLY A   2
ENDMDL
MODEL        6
ATOM      1  N   MET A   1
ATOM      2  CA  MET A   1
ATOM      3  N   GLY A   2                                                    X
ATOM      4  CA  GLY A   2
ENDMDL
"""))
  models = pdb_inp.construct_hierarchy().models()
  assert models[0].is_identical_hierarchy(models[1])
  assert models[1].is_identical_hierarchy(models[0])
  assert models[0].only_chain().is_identical_hierarchy(
    other=models[1].only_chain())
  assert models[0].only_chain().residue_groups()[0].is_identical_hierarchy(
    other=models[1].only_chain().residue_groups()[0])
  for other in models[2:]:
    assert not models[0].is_identical_hierarchy(other=other)
    assert not other.is_identical_hierarchy(other=models[0])
    assert not models[0].only_chain().is_identical_hierarchy(
      other=other.only_chain())
    assert not models[0].only_chain().residue_groups()[0] \
      .is_identical_hierarchy(other=models[2].only_chain().residue_groups()[0])

def exercise_is_similar_hierarchy():
  s0 = """\
MODEL        0
ATOM      1  %s  MET A   1
ATOM      2  CA AMET A   1
ATOM      3  N   GLY A   2
ATOM      4  CA AGLY A   2
ENDMDL
MODEL        1
ATOM      1  N   MET A   1
ATOM      2  CA AMET A   1
ATOM      3  N   %3s A   2
ATOM      4  CA A%3s A   2
ENDMDL
"""
  i1 = pdb.input(source_info=None, lines=flex.split_lines(
    s0 % ("N ", "GLY", "GLY")))
  h1 = i1.construct_hierarchy()
  assert h1.is_similar_hierarchy(other=h1)
  assert h1.is_similar_hierarchy(other=h1.deep_copy())
  assert h1.models()[0].is_similar_hierarchy(
    other=h1.models()[1])
  assert h1.models()[0].only_chain().is_similar_hierarchy(
    other=h1.models()[1].only_chain())
  assert h1.models()[0].only_chain().residue_groups()[0].is_similar_hierarchy(
    other=h1.models()[1].only_chain().residue_groups()[0])
  for an,rn in [("C ", "GLY"), ("N ", "ALA")]:
    i2 = pdb.input(source_info=None, lines=flex.split_lines(
      s0 % (an, rn, rn)))
    h2 = i2.construct_hierarchy()
    assert not h1.is_similar_hierarchy(other=h2)
    assert not h2.is_similar_hierarchy(other=h1)
    assert not h1.models()[0].is_similar_hierarchy(
      other=h2.models()[0]) == (rn == "GLY")
    assert h1.models()[0].only_chain().is_similar_hierarchy(
      other=h2.models()[0].only_chain()) == (an == "N ")
    assert h1.models()[0].only_chain().residue_groups()[0] \
      .is_similar_hierarchy(
        other=h2.models()[0].only_chain().residue_groups()[0]) == (an == "N ")

def exercise_atoms():
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM      1  N   GLN A   3      35.299  11.075  19.070  1.00 36.89           N
ATOM      2  CA  GLN A   3      34.482   9.927  18.794  0.63 37.88           C
SIGATM    2  CA  GLN A   3       1.200   2.300   3.400  0.04  0.05           C
ANISOU    2  CA  GLN A   3     7794   3221   3376  -1227   1064   2601       C
ATOM      3  Q   GLN A   3      35.130   8.880  17.864  0.84 37.52           C
ANISOU    3  Q   GLN A   3     7875   3041   3340   -981    727   2663       C
SIGUIJ    3  Q   GLN A   3       75     41     40     -1      7     63       C
ATOM      4  O   GLN A   3      34.548   7.819  17.724  1.00 38.54      STUV
ATOM      5 1CB AGLN A   3      32.979  10.223  18.469  1.00 37.80          X
HETATM    6 CA  AION B   1      32.360  11.092  17.308  0.92 35.96          CA2+
HETATM    7 CA   ION B   2      30.822  10.665  17.190  1.00 36.87
"""))
  atoms = pdb_inp.atoms()
  assert list(atoms.extract_serial()) == [
    "    1", "    2", "    3", "    4", "    5", "    6", "    7"]
  assert list(atoms.extract_name()) == [
    " N  ", " CA ", " Q  ", " O  ", "1CB ", "CA  ", "CA  "]
  assert (list(atoms.extract_segid()) == [
    '    ', '    ', '    ', 'STUV', '    ', '    ', '    '])
  xyz = atoms.extract_xyz()
  assert approx_equal(xyz, [
    (35.299,11.075,19.070),
    (34.482,9.927,18.794),
    (35.130,8.880,17.864),
    (34.548,7.819,17.724),
    (32.979,10.223,18.469),
    (32.360,11.092,17.308),
    (30.822,10.665,17.190)])
  sigxyz = atoms.extract_sigxyz()
  assert approx_equal(sigxyz, [
    (0,0,0),
    (1.2,2.3,3.4),
    (0,0,0),
    (0,0,0),
    (0,0,0),
    (0,0,0),
    (0,0,0)])
  occ = atoms.extract_occ()
  assert approx_equal(occ,
    [1.00,0.63,0.84,1.00,1.00,0.92,1.00])
  sigocc = atoms.extract_sigocc()
  assert approx_equal(sigocc,
    [0,0.04,0,0,0,0,0])
  b = atoms.extract_b()
  assert approx_equal(b,
    [36.89,37.88,37.52,38.54,37.80,35.96,36.87])
  sigb = atoms.extract_sigb()
  assert approx_equal(sigb,
    [0,0.05,0,0,0,0,0])
  uij = atoms.extract_uij()
  assert approx_equal(uij, [
    (-1,-1,-1,-1,-1,-1),
    (0.7794, 0.3221, 0.3376, -0.1227, 0.1064, 0.2601),
    (0.7875, 0.3041, 0.3340, -0.0981, 0.0727, 0.2663),
    (-1,-1,-1,-1,-1,-1),
    (-1,-1,-1,-1,-1,-1),
    (-1,-1,-1,-1,-1,-1),
    (-1,-1,-1,-1,-1,-1)])
  expected_siguij = [
    (-1,-1,-1,-1,-1,-1),
    (-1,-1,-1,-1,-1,-1),
    (0.0075, 0.0041, 0.0040, -0.0001, 0.0007, 0.0063),
    (-1,-1,-1,-1,-1,-1),
    (-1,-1,-1,-1,-1,-1),
    (-1,-1,-1,-1,-1,-1),
    (-1,-1,-1,-1,-1,-1)]
  if (pdb.hierarchy.atom.has_siguij()):
    siguij = atoms.extract_siguij()
    assert approx_equal(siguij, expected_siguij)
  assert list(atoms.extract_hetero()) == [5,6]
  assert list(atoms.extract_element()) == [" N"," C"," C","  ","X ","CA","  "]
  assert atoms.extract_element(strip=False).all_eq(atoms.extract_element())
  assert list(atoms.extract_element(strip=True)) ==["N","C","C","","X","CA",""]
  atoms.reset_i_seq()
  assert list(atoms.extract_i_seq()) == range(7)
  assert list(
    pdb.hierarchy.af_shared_atom([atoms[3], atoms[6]]).extract_i_seq()) == \
      [3, 6]
  assert list(atoms.extract_tmp_as_size_t()) == [0]*7
  for i,atom in enumerate(atoms): atom.tmp = i+3
  assert list(atoms.extract_tmp_as_size_t()) == [3,4,5,6,7,8,9]
  atoms[3].tmp = -1
  try: atoms.extract_tmp_as_size_t()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "atom.tmp less than zero: cannot convert to unsigned value.")
  else: raise Exception_expected
  #
  assert atoms.set_xyz(new_xyz=xyz+(1,2,3)) is atoms
  assert approx_equal(atoms.extract_xyz(), [
    (36.299,13.075,22.070),
    (35.482,11.927,21.794),
    (36.130,10.880,20.864),
    (35.548,9.819,20.724),
    (33.979,12.223,21.469),
    (33.360,13.092,20.308),
    (31.822,12.665,20.190)])
  assert atoms.set_sigxyz(new_sigxyz=sigxyz+(1,2,3)) is atoms
  assert approx_equal(atoms.extract_sigxyz(), [
    (1,2,3),
    (2.2,4.3,6.4),
    (1,2,3),
    (1,2,3),
    (1,2,3),
    (1,2,3),
    (1,2,3)])
  assert atoms.set_occ(new_occ=occ+1.23) is atoms
  assert approx_equal(atoms.extract_occ(),
    [2.23,1.86,2.07,2.23,2.23,2.15,2.23])
  assert atoms.set_sigocc(new_sigocc=sigocc+3) is atoms
  assert approx_equal(atoms.extract_sigocc(),
    [3,3.04,3,3,3,3,3])
  assert atoms.set_b(new_b=b+10) is atoms
  assert approx_equal(atoms.extract_b(),
    [46.89,47.88,47.52,48.54,47.80,45.96,46.87])
  assert atoms.set_sigb(new_sigb=sigb+5) is atoms
  assert approx_equal(atoms.extract_sigb(),
    [5,5.05,5,5,5,5,5])
  assert atoms.set_uij(new_uij=flex.sym_mat3_double(expected_siguij)) is atoms
  assert approx_equal(atoms.extract_uij(), expected_siguij)
  if (pdb.hierarchy.atom.has_siguij()):
    assert atoms.set_siguij(new_siguij=uij) is atoms
    assert approx_equal(atoms.extract_siguij(), [
      (-1,-1,-1,-1,-1,-1),
      (0.7794, 0.3221, 0.3376, -0.1227, 0.1064, 0.2601),
      (0.7875, 0.3041, 0.3340, -0.0981, 0.0727, 0.2663),
      (-1,-1,-1,-1,-1,-1),
      (-1,-1,-1,-1,-1,-1),
      (-1,-1,-1,-1,-1,-1),
      (-1,-1,-1,-1,-1,-1)])
  #
  h = pdb_inp.construct_hierarchy()
  for i in xrange(2):
    s = h.as_pdb_string()
    d = hashlib_md5(s).hexdigest()
    if (pdb.hierarchy.atom.has_siguij()):
      assert d == "c4089359af431bb2962d6a8e457dd86f"
    else:
      assert d == "a1dd6605ed08b56862b9d7ae6b9a547b"
    h.write_pdb_file(file_name="tmp_tst_hierarchy.pdb")
    assert not show_diff(open("tmp_tst_hierarchy.pdb").read(), s)
    h = pdb.input(
      source_info=None, lines=flex.split_lines(s)).construct_hierarchy()
  #
  atoms = h.atoms().select(indices=flex.size_t([2,5,3,0]))
  assert [a.name for a in atoms] == [" Q  ", "CA  ", " O  ", " N  "]
  atoms = atoms.select(indices=flex.size_t([3,0,1,2]), reverse=True)
  assert [a.name for a in atoms] == ["CA  ", " O  ", " N  ", " Q  "]
  atoms = atoms.select(indices=flex.size_t([3,0,1,2]), reverse=False)
  assert [a.name for a in atoms] == [" Q  ", "CA  ", " O  ", " N  "]
  atoms = atoms.select(flex.bool([False,True,False,False]))
  assert [a.name for a in atoms] == ["CA  "]
  #
  assert [int(a.serial) for a in h.atoms_with_i_seq_mismatch()] \
      == [2, 3, 4, 5, 6, 7]
  h.atoms().reset_i_seq()
  assert h.atoms_with_i_seq_mismatch().size() == 0
  assert pdb.hierarchy.root().atoms_with_i_seq_mismatch().size() == 0

def check_wpf(hierarchy, kwargs={}, trailing=None, expected=None):
  if ("atoms_reset_serial_first_value" in kwargs):
    pdb_str = hierarchy.deep_copy().as_pdb_string(**kwargs)
  else:
    pdb_str = hierarchy.as_pdb_string(**kwargs)
  if (trailing is not None): pdb_str = pdb_str.replace(trailing, "")
  if (expected is None):
    sys.stdout.write(pdb_str)
  else:
    assert not show_diff(pdb_str, expected)
  hierarchy.write_pdb_file(file_name="tmp_tst_hierarchy.pdb", **kwargs)
  pdb_file = open("tmp_tst_hierarchy.pdb").read()
  if (trailing is not None): pdb_file = pdb_file.replace(trailing, "")
  assert not show_diff(pdb_file, pdb_str)
  #
  pdb_inp = pdb.input(file_name="tmp_tst_hierarchy.pdb")
  assert pdb_inp.atoms().size() == hierarchy.atoms_size()
  kwargs = dict(kwargs)
  for discard in ["atoms_reset_serial_first_value", "interleaved_conf"]:
    if (discard in kwargs): del kwargs[discard]
  pdb_inp.write_pdb_file(file_name="tmp2.pdb", **kwargs)
  pdb_str2 = pdb_inp.as_pdb_string(**kwargs)
  assert not show_diff(open("tmp2.pdb").read(), pdb_str2)
  pdb_inp2 = pdb.input(file_name="tmp2.pdb")
  assert pdb_inp2.atoms().size() == pdb_inp.atoms().size()
  assert pdb_inp.extract_cryst1_z_columns() \
      == pdb_inp2.extract_cryst1_z_columns()
  if ("cryst1_z" in kwargs):
    assert pdb_inp.extract_cryst1_z_columns() == ("%4s" % kwargs["cryst1_z"])
  #
  return pdb_str

def exercise_atoms_interleaved_conf():
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM      1  N  ATRP A   1I
ATOM      4  N  CPHE A   1I
ATOM      2  C  ATRP A   1I
ATOM      8  CA BTYR A   1I
ATOM      3  O  ATRP A   1I
ATOM      5  CA CPHE A   1I
ATOM      6  C  CPHE A   1I
ATOM      9  C  BTYR A   1I
ATOM      7  O  CPHE A   1I
ATOM     10  O  BTYR A   1I
"""))
  hierarchy = pdb_inp.construct_hierarchy()
  for obj in [hierarchy.only_residue_group(),
              hierarchy.only_chain(),
              hierarchy.only_model(),
              hierarchy]:
    for interleaved_conf in [-1,0,1]:
      rs = [atom.format_atom_record(replace_floats_with="")
        for atom in obj.atoms(interleaved_conf=interleaved_conf)]
      assert not show_diff("\n".join([r[:-8] for r in rs]), """\
ATOM      1  N  ATRP A   1I
ATOM      2  C  ATRP A   1I
ATOM      3  O  ATRP A   1I
ATOM      4  N  CPHE A   1I
ATOM      5  CA CPHE A   1I
ATOM      6  C  CPHE A   1I
ATOM      7  O  CPHE A   1I
ATOM      8  CA BTYR A   1I
ATOM      9  C  BTYR A   1I
ATOM     10  O  BTYR A   1I""")
  for obj in [hierarchy.only_residue_group(),
              hierarchy.only_chain(),
              hierarchy.only_model(),
              hierarchy]:
    rs = [atom.format_atom_record(replace_floats_with="")
      for atom in obj.atoms(interleaved_conf=2)]
    assert not show_diff("\n".join([r[:-8] for r in rs]), """\
ATOM      1  N  ATRP A   1I
ATOM      4  N  CPHE A   1I
ATOM      2  C  ATRP A   1I
ATOM      6  C  CPHE A   1I
ATOM      9  C  BTYR A   1I
ATOM      3  O  ATRP A   1I
ATOM      7  O  CPHE A   1I
ATOM     10  O  BTYR A   1I
ATOM      5  CA CPHE A   1I
ATOM      8  CA BTYR A   1I""")
  trailing = " A   1I      0.000   0.000   0.000  0.00  0.00"
  check_wpf(hierarchy, {"interleaved_conf":1}, trailing, """\
ATOM      1  N  ATRP
ATOM      2  C  ATRP
ATOM      3  O  ATRP
ATOM      4  N  CPHE
ATOM      5  CA CPHE
ATOM      6  C  CPHE
ATOM      7  O  CPHE
ATOM      8  CA BTYR
ATOM      9  C  BTYR
ATOM     10  O  BTYR
TER
""")
  check_wpf(hierarchy, {"interleaved_conf":2}, trailing, """\
ATOM      1  N  ATRP
ATOM      4  N  CPHE
ATOM      2  C  ATRP
ATOM      6  C  CPHE
ATOM      9  C  BTYR
ATOM      3  O  ATRP
ATOM      7  O  CPHE
ATOM     10  O  BTYR
ATOM      5  CA CPHE
ATOM      8  CA BTYR
TER
""")
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM      1  N  ATRP A   1I
ATOM      4  N  CPHE A   1I
ATOM      2  C  ATRP A   1I
ATOM      8  CA BTRP A   1I
ATOM      3  O  ATRP A   1I
ATOM      5  CA CPHE A   1I
ATOM      6  C  CPHE A   1I
ATOM      9  C  BTRP A   1I
ATOM      7  O  CPHE A   1I
ATOM     10  O  BTRP A   1I
"""))
  hierarchy = pdb_inp.construct_hierarchy()
  rs = [atom.format_atom_record(replace_floats_with="")
    for atom in hierarchy.only_residue_group().atoms(interleaved_conf=1)]
  assert not show_diff("\n".join([r[:-8] for r in rs]), """\
ATOM      1  N  ATRP A   1I
ATOM      2  C  ATRP A   1I
ATOM      9  C  BTRP A   1I
ATOM      3  O  ATRP A   1I
ATOM     10  O  BTRP A   1I
ATOM      8  CA BTRP A   1I
ATOM      4  N  CPHE A   1I
ATOM      5  CA CPHE A   1I
ATOM      6  C  CPHE A   1I
ATOM      7  O  CPHE A   1I""")
  rs = [atom.format_atom_record(replace_floats_with="")
    for atom in hierarchy.only_residue_group().atoms(interleaved_conf=2)]
  assert not show_diff("\n".join([r[:-8] for r in rs]), """\
ATOM      1  N  ATRP A   1I
ATOM      4  N  CPHE A   1I
ATOM      2  C  ATRP A   1I
ATOM      6  C  CPHE A   1I
ATOM      9  C  BTRP A   1I
ATOM      3  O  ATRP A   1I
ATOM      7  O  CPHE A   1I
ATOM     10  O  BTRP A   1I
ATOM      5  CA CPHE A   1I
ATOM      8  CA BTRP A   1I""")
  check_wpf(hierarchy, {"interleaved_conf":1}, trailing, """\
ATOM      1  N  ATRP
ATOM      2  C  ATRP
ATOM      9  C  BTRP
ATOM      3  O  ATRP
ATOM     10  O  BTRP
ATOM      8  CA BTRP
ATOM      4  N  CPHE
ATOM      5  CA CPHE
ATOM      6  C  CPHE
ATOM      7  O  CPHE
TER
""")
  for interleaved_conf in [2,3,4]:
    check_wpf(hierarchy, {"interleaved_conf":interleaved_conf}, trailing, """\
ATOM      1  N  ATRP
ATOM      4  N  CPHE
ATOM      2  C  ATRP
ATOM      6  C  CPHE
ATOM      9  C  BTRP
ATOM      3  O  ATRP
ATOM      7  O  CPHE
ATOM     10  O  BTRP
ATOM      5  CA CPHE
ATOM      8  CA BTRP
TER
""")
  check_wpf(hierarchy,
    {"interleaved_conf":2,
     "atoms_reset_serial_first_value": 13}, trailing, """\
ATOM     13  N  ATRP
ATOM     14  N  CPHE
ATOM     15  C  ATRP
ATOM     16  C  CPHE
ATOM     17  C  BTRP
ATOM     18  O  ATRP
ATOM     19  O  CPHE
ATOM     20  O  BTRP
ATOM     21  CA CPHE
ATOM     22  CA BTRP
TER
""")
  assert list(hierarchy.atoms().extract_serial()) == [
    "   13", "   15", "   18", "   14", "   21",
    "   16", "   19", "   22", "   17", "   20"]
  #
  hierarchy = pdb.input(source_info=None, lines=flex.split_lines("""\
MODEL        1
ATOM      1  N   MET A   1
ATOM      2  CA  MET A   1
ATOM      3  N   GLY A   2
ATOM      4  CA  GLY A   2
ENDMDL
MODEL        2
ATOM      1  N   MET A   1
ATOM      2  CA  MET A   1
ATOM      3  N   GLY A   2
ATOM      4  CA  GLY A   2
ENDMDL
""")).construct_hierarchy()
  hierarchy.atoms_reset_serial(interleaved_conf=0, first_value=10)
  assert list(hierarchy.atoms().extract_serial()) \
      == ["   10", "   11", "   12", "   13"] * 2
  hierarchy.atoms_reset_serial()
  assert list(hierarchy.atoms().extract_serial()) \
      == ["    1", "    2", "    3", "    4"] * 2

def exercise_as_pdb_string(pdb_file_names, comprehensive):
  pdb_string = """\
HETATM  145  C21 DA7  3014      18.627   3.558  25.202  0.50 29.50           C
ATOM    146  C8 ADA7  3015       9.021 -13.845  22.131  0.50 26.57           C
"""
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines(pdb_string))
  hierarchy = pdb_inp.construct_hierarchy()
  check_wpf(hierarchy, expected=pdb_string+"TER\n")
  rem = "REMARK EXERCISE"
  for obj in [pdb_inp, hierarchy]:
    print >> open("tmp_tst_hierarchy.pdb", "w"), rem
    obj.write_pdb_file(
      file_name="tmp_tst_hierarchy.pdb", open_append=True, append_end=True)
    assert not show_diff(
      open("tmp_tst_hierarchy.pdb").read(), rem+"\n"+pdb_string+"TER\nEND\n")
  check_wpf(
    hierarchy,
    kwargs={"crystal_symmetry": (2,3,4,80,90,100), "cryst1_z": 5},
    expected="""\
CRYST1    2.000    3.000    4.000  80.00  90.00 100.00 P 1           5
SCALE1      0.500000  0.088163 -0.015793        0.00000
SCALE2      0.000000  0.338476 -0.060632        0.00000
SCALE3      0.000000  0.000000  0.253979        0.00000
""" + pdb_string+"TER\n")
  check_wpf(
    hierarchy,
    kwargs={"cryst1_z": 7, "write_scale_records": False},
    expected="""\
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           7
""" + pdb_string+"TER\n")
  for obj in [pdb_inp, hierarchy]:
    print >> open("tmp_tst_hierarchy.pdb", "w"), rem
    obj.write_pdb_file(
      file_name="tmp_tst_hierarchy.pdb",
      cryst1_z="",
      write_scale_records=False,
      open_append=True)
    assert not show_diff(open("tmp_tst_hierarchy.pdb").read(), rem + """
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1
""" + pdb_string+"TER\n")
  #
  # only pdb entry (as of 2008 Mar 26) for which
  # s1 = h1.as_pdb_string(interleaved_conf=1)
  # s2 = h2.as_pdb_string(interleaved_conf=1)
  # leads to s1 != s2
  i1 = pdb.input(source_info=None, lines=flex.split_lines("""\
HEADER    TRANSPORT PROTEIN                       21-AUG-07   2V93
ATOM   1549  O  ACYS A 211      24.080  12.057  26.978  0.95  0.00           O
ATOM   1550  SG ACYS A 211      24.960  15.526  27.230  0.16  0.00           S
ATOM   1551  SG BCYS A 211      24.418  16.196  28.728  0.16  0.00           S
ATOM   1552  N  CCYS A 211      12.649  13.688  18.306  0.05  0.00           N
ATOM   1553  CA CCYS A 211      12.104  13.010  17.135  0.05  0.00           C
ATOM   1554  C  CCYS A 211      13.148  12.135  16.481  0.05  0.00           C
ATOM   1555  O  CCYS A 211      12.770  11.031  16.051  0.05  0.00           O
ATOM   1556  CB CCYS A 211      11.540  14.072  16.172  0.05  0.00           C
ATOM   1557  SG CCYS A 211      24.361  15.163  26.341  0.16  0.00           S
"""))
  h1 = i1.construct_hierarchy()
  s1 = h1.as_pdb_string(interleaved_conf=1)
  assert not show_diff(s1, """\
ATOM   1549  O  ACYS A 211      24.080  12.057  26.978  0.95  0.00           O
ATOM   1555  O  CCYS A 211      12.770  11.031  16.051  0.05  0.00           O
ATOM   1550  SG ACYS A 211      24.960  15.526  27.230  0.16  0.00           S
ATOM   1551  SG BCYS A 211      24.418  16.196  28.728  0.16  0.00           S
ATOM   1557  SG CCYS A 211      24.361  15.163  26.341  0.16  0.00           S
ATOM   1552  N  CCYS A 211      12.649  13.688  18.306  0.05  0.00           N
ATOM   1553  CA CCYS A 211      12.104  13.010  17.135  0.05  0.00           C
ATOM   1554  C  CCYS A 211      13.148  12.135  16.481  0.05  0.00           C
ATOM   1556  CB CCYS A 211      11.540  14.072  16.172  0.05  0.00           C
TER
""")
  i2 = pdb.input(source_info=None, lines=flex.split_lines(s1))
  h2 = i2.construct_hierarchy()
  s2 = h2.as_pdb_string(interleaved_conf=1)
  assert not show_diff(s2, """\
ATOM   1549  O  ACYS A 211      24.080  12.057  26.978  0.95  0.00           O
ATOM   1555  O  CCYS A 211      12.770  11.031  16.051  0.05  0.00           O
ATOM   1550  SG ACYS A 211      24.960  15.526  27.230  0.16  0.00           S
ATOM   1557  SG CCYS A 211      24.361  15.163  26.341  0.16  0.00           S
ATOM   1551  SG BCYS A 211      24.418  16.196  28.728  0.16  0.00           S
ATOM   1552  N  CCYS A 211      12.649  13.688  18.306  0.05  0.00           N
ATOM   1553  CA CCYS A 211      12.104  13.010  17.135  0.05  0.00           C
ATOM   1554  C  CCYS A 211      13.148  12.135  16.481  0.05  0.00           C
ATOM   1556  CB CCYS A 211      11.540  14.072  16.172  0.05  0.00           C
TER
""")
  assert h1.is_similar_hierarchy(other=h2)
  assert h2.is_similar_hierarchy(other=h1)
  #
  open("tmp_tst_hierarchy.pdb", "wb").write("""\
CRYST1    2.000    3.000    4.000  90.00  80.00  90.00 P 2           5
ATOM      0  S   SO4     0       3.302   8.419   8.560  1.00 10.00           S
ATOM      1  O1  SO4     0       3.497   8.295   7.118  1.00 10.00           O
ATOM      2  O2 ASO4     0       3.098   7.095   9.140  0.80 10.00           O
ATOM      3  O2 BSO4     0       3.498   7.495   9.440  0.20 10.00           O
ATOM      4  O3  SO4     0       4.481   9.037   9.159  1.00 10.00           O
ATOM      5  O4  SO4     0       2.131   9.251   8.823  1.00 10.00           O
""")
  pdb.rewrite_normalized(
    input_file_name="tmp_tst_hierarchy.pdb",
    output_file_name="tmp_norm.pdb")
  assert not show_diff(open("tmp_norm.pdb").read(), """\
CRYST1    2.000    3.000    4.000  90.00  80.00  90.00 P 1 2 1
SCALE1      0.500000  0.000000 -0.088163        0.00000
SCALE2      0.000000  0.333333  0.000000        0.00000
SCALE3      0.000000  0.000000  0.253857        0.00000
ATOM      1  S   SO4     0       3.302   8.419   8.560  1.00 10.00           S
ATOM      2  O1  SO4     0       3.497   8.295   7.118  1.00 10.00           O
ATOM      3  O3  SO4     0       4.481   9.037   9.159  1.00 10.00           O
ATOM      4  O4  SO4     0       2.131   9.251   8.823  1.00 10.00           O
ATOM      5  O2 ASO4     0       3.098   7.095   9.140  0.80 10.00           O
ATOM      6  O2 BSO4     0       3.498   7.495   9.440  0.20 10.00           O
TER
END
""")
  pdb.rewrite_normalized(
    input_file_name="tmp_tst_hierarchy.pdb",
    output_file_name="tmp_norm.pdb",
    keep_original_crystallographic_section=True,
    keep_original_atom_serial=True)
  assert not show_diff(open("tmp_norm.pdb").read(), """\
CRYST1    2.000    3.000    4.000  90.00  80.00  90.00 P 2           5
ATOM      0  S   SO4     0       3.302   8.419   8.560  1.00 10.00           S
ATOM      1  O1  SO4     0       3.497   8.295   7.118  1.00 10.00           O
ATOM      4  O3  SO4     0       4.481   9.037   9.159  1.00 10.00           O
ATOM      5  O4  SO4     0       2.131   9.251   8.823  1.00 10.00           O
ATOM      2  O2 ASO4     0       3.098   7.095   9.140  0.80 10.00           O
ATOM      3  O2 BSO4     0       3.498   7.495   9.440  0.20 10.00           O
TER
END
""")
  #
  if (pdb_file_names is None):
    print "Skipping exercise_as_pdb_string(): input files not available"
    return
  prev_file_name = None
  for file_name in pdb_file_names:
    if (not comprehensive and random.random() > 0.05):
      continue
    pdb_inp_1 = pdb.input(file_name=file_name)
    hierarchy_1 = pdb_inp_1.construct_hierarchy()
    pdb_str_1 = hierarchy_1.as_pdb_string(append_end=False)
    pdb_inp_2 = pdb.input(
      source_info=None, lines=flex.split_lines(pdb_str_1))
    hierarchy_2 = pdb_inp_2.construct_hierarchy()
    check_wpf(
      hierarchy=hierarchy_2,
      kwargs={"append_end": True},
      expected=pdb_str_1+"END\n")
    assert hierarchy_1.is_similar_hierarchy(other=hierarchy_2)
    assert hierarchy_2.is_similar_hierarchy(other=hierarchy_1)
    if (prev_file_name is not None):
      if (   hierarchy_1.is_similar_hierarchy(other=prev_hierarchy)
          or prev_hierarchy.is_similar_hierarchy(other=hierarchy_1)):
        # some files are known to be similar
        p = os.path.basename(prev_file_name)
        c = os.path.basename(     file_name)
        if (    p[:3] != c[:3]
            and (len(p) < 15 or len(c) < 15 or p[-12:] != c[-12:])):
          print "WARNING: similar hierarchies:"
          print " ", show_string(prev_file_name)
          print " ", show_string(file_name)
    prev_file_name = file_name
    prev_pdb_str = pdb_str_1
    prev_hierarchy = hierarchy_1
    #
    awls = []
    for obj in [pdb_inp_2, hierarchy_2]:
      sio = StringIO()
      for awl in obj.atoms_with_labels():
        print >> sio, awl.format_atom_record_group(), \
          int(awl.is_first_in_chain), \
          int(awl.is_first_after_break), \
          awl.model_id
      awls.append(sio.getvalue())
    assert not show_diff(*awls), file_name

def exercise_atom_with_labels():
  awl = pdb.hierarchy.atom_with_labels()
  assert not show_diff(awl.format_atom_record(), """\
ATOM                             0.000   0.000   0.000  0.00  0.00""")
  assert not show_diff(awl.format_atom_record(replace_floats_with="#"), """\
ATOM                       #        """)
  assert not show_diff(awl.format_atom_record_group(), """\
ATOM                             0.000   0.000   0.000  0.00  0.00""")
  assert awl.serial == ""
  assert awl.name == ""
  assert awl.model_id == ""
  assert awl.chain_id == ""
  assert awl.resseq == ""
  assert awl.icode == ""
  assert awl.altloc == ""
  assert awl.resname == ""
  assert not awl.is_first_in_chain
  assert not awl.is_first_after_break
  awlc = awl.detached_copy()
  awl.serial = "12345"
  assert awlc.serial == ""
  awl.name = "NaMe"
  awl.model_id = "MoDl"
  assert awl.model_id == "MoDl"
  assert awlc.model_id == ""
  awl.chain_id = "Ch"
  awl.resseq = "ABCD"
  assert awl.resseq_as_int() == 24701
  awl.icode = "I"
  assert awl.resid() == "ABCDI"
  awl.altloc = "l"
  awl.resname = "rNm"
  awl.is_first_in_chain = True
  assert awl.is_first_in_chain
  awl.is_first_after_break = True
  assert awl.is_first_after_break
  assert not show_diff(awl.format_atom_record_group(), """\
ATOM  12345 NaMelrNmChABCDI      0.000   0.000   0.000  0.00  0.00""")
  assert not show_diff(awl.format_atom_record(), """\
ATOM  12345 NaMelrNmChABCDI      0.000   0.000   0.000  0.00  0.00""")
  assert not show_diff(awl.quote(), '''\
"ATOM  12345 NaMelrNmChABCDI.*.        "''')
  assert not show_diff(awl.quote(full=True), '''\
"ATOM  12345 NaMelrNmChABCDI      0.000   0.000   0.000  0.00  0.00"''')
  awl.sigxyz = (0.1,0.2,0.3)
  awl.uij = (1,2,3,0.1,0.2,0.3)
  awl.siguij = (3,1,2,0.3,0.1,0.2)
  assert not show_diff(awl.format_sigatm_record(), """\
SIGATM12345 NaMelrNmChABCDI      0.100   0.200   0.300  0.00  0.00""")
  assert not show_diff(awl.format_anisou_record(), """\
ANISOU12345 NaMelrNmChABCDI   10000  20000  30000   1000   2000   3000""")
  if (pdb.hierarchy.atom.has_siguij()):
    assert not show_diff(awl.format_siguij_record(), """\
SIGUIJ12345 NaMelrNmChABCDI   30000  10000  20000   3000   1000   2000""")
    assert not show_diff(awl.format_atom_record_group(), """\
ATOM  12345 NaMelrNmChABCDI      0.000   0.000   0.000  0.00  0.00
SIGATM12345 NaMelrNmChABCDI      0.100   0.200   0.300  0.00  0.00
ANISOU12345 NaMelrNmChABCDI   10000  20000  30000   1000   2000   3000
SIGUIJ12345 NaMelrNmChABCDI   30000  10000  20000   3000   1000   2000""")
  else:
    assert not show_diff(awl.format_siguij_record(), """\
SIGUIJ12345 NaMelrNmChABCDI  -10000 -10000 -10000 -10000 -10000 -10000""")
    assert not show_diff(awl.format_atom_record_group(), """\
ATOM  12345 NaMelrNmChABCDI      0.000   0.000   0.000  0.00  0.00
SIGATM12345 NaMelrNmChABCDI      0.100   0.200   0.300  0.00  0.00
ANISOU12345 NaMelrNmChABCDI   10000  20000  30000   1000   2000   3000""")
  assert not show_diff(
    awl.format_atom_record_group(sigatm=False, anisou=False, siguij=False),
    """\
ATOM  12345 NaMelrNmChABCDI      0.000   0.000   0.000  0.00  0.00""")
  #
  assert not show_diff(awl.id_str(),
    'model="MoDl" pdb="NaMelrNmChABCDI"')
  assert not show_diff(awl.id_str(pdbres=True),
    'model="MoDl" pdbres="rNmChABCDI"')
  awl.segid = "sEgI"
  assert not show_diff(awl.id_str(),
    'model="MoDl" pdb="NaMelrNmChABCDI" segid="sEgI"')
  assert not show_diff(awl.id_str(suppress_segid=True),
    'model="MoDl" pdb="NaMelrNmChABCDI"')
  assert not show_diff(awl.id_str(pdbres=True),
    'model="MoDl" pdbres="rNmChABCDI" segid="sEgI"')
  awl.model_id = ""
  assert not show_diff(awl.id_str(suppress_segid=False),
    'pdb="NaMelrNmChABCDI" segid="sEgI"')
  assert not show_diff(awl.id_str(pdbres=True),
    'pdbres="rNmChABCDI" segid="sEgI"')
  assert not show_diff(awl.id_str(pdbres=True, suppress_segid=True),
    'pdbres="rNmChABCDI"')
  awl.segid = "    "
  assert not show_diff(awl.id_str(),
    'pdb="NaMelrNmChABCDI"')
  assert not show_diff(awl.id_str(pdbres=True),
    'pdbres="rNmChABCDI"')
  #
  a = pdb.hierarchy.atom()
  a.serial = "   1A"
  try: a.serial_as_int()
  except (RuntimeError, ValueError), e:
    assert not show_diff(str(e), 'invalid atom serial number: "   1A"')
  else: raise Exception_expected
  #
  awl.serial = "   1A"
  try: awl.serial_as_int()
  except (RuntimeError, ValueError), e:
    assert not show_diff(str(e), """\
invalid atom serial number:
  ATOM     1A NaMelrNmChABCDI      0.000   0.000   0.000  0.00  0.00
        ^^^^^""")
  else: raise Exception_expected
  #
  awl.resseq = " 18A"
  try: awl.resseq_as_int()
  except (RuntimeError, ValueError), e:
    assert not show_diff(str(e), """\
invalid residue sequence number:
  ATOM     1A NaMelrNmCh 18AI      0.000   0.000   0.000  0.00  0.00
                        ^^^^""")
  else: raise Exception_expected
  #
  h = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM                    1B
""")).construct_hierarchy()
  for r in [h.only_residue_group(), h.only_residue()]:
    try: r.resseq_as_int()
    except (RuntimeError, ValueError), e:
      assert not show_diff(str(e), """\
invalid residue sequence number:
  ATOM                    1B       0.000   0.000   0.000  0.00  0.00
                        ^^^^""")
    else: raise Exception_expected
  #
  ch = pdb.hierarchy.chain(id="C")
  ch.append_residue_group(pdb.hierarchy.residue_group(resseq="  1C"))
  ch.residue_groups()[0].append_atom_group(pdb.hierarchy.atom_group())
  for r in[ch.only_residue_group(), ch.only_conformer().only_residue()]:
    try: r.resseq_as_int()
    except (RuntimeError, ValueError), e:
      assert not show_diff(str(e), 'invalid residue sequence number: "  1C"')
    else: raise Exception_expected
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
MODEL        0
ATOM      1  C   MET A   1
ATOM      2  CA AMET A   1
BREAK
ATOM      3  N   GLY A   2
ATOM      4  CA AGLY A   2
ENDMDL
MODEL        1
ATOM      1  N   MET A   1
ATOM      2  CA AMET A   1
ATOM      3  N   ALA B   2
ATOM      4  CA AALA B   2
ENDMDL
"""))
  hierarchy = pdb_inp.construct_hierarchy()
  for obj in [hierarchy, pdb_inp]:
    for fmt in ["format_atom_record_group", "format_atom_record_group"]:
      sio = StringIO()
      for awl in obj.atoms_with_labels():
        assert awl.parent(optional=False) is not None
        awlc = awl.detached_copy()
        assert awlc.parent() is None
        try:
          awlc.parent(optional=False)
        except RuntimeError, e:
          assert not show_diff(str(e), "atom has no parent atom_group")
        else: raise Exception_expected
        print >> sio, getattr(awl, fmt)(), \
          int(awl.is_first_in_chain), \
          int(awl.is_first_after_break), \
          awl.model_id
      assert not show_diff(sio.getvalue(), """\
ATOM      1  C   MET A   1       0.000   0.000   0.000  0.00  0.00 1 0    0
ATOM      2  CA AMET A   1       0.000   0.000   0.000  0.00  0.00 0 0    0
ATOM      3  N   GLY A   2       0.000   0.000   0.000  0.00  0.00 0 1    0
ATOM      4  CA AGLY A   2       0.000   0.000   0.000  0.00  0.00 0 0    0
ATOM      1  N   MET A   1       0.000   0.000   0.000  0.00  0.00 1 0    1
ATOM      2  CA AMET A   1       0.000   0.000   0.000  0.00  0.00 0 0    1
ATOM      3  N   ALA B   2       0.000   0.000   0.000  0.00  0.00 1 0    1
ATOM      4  CA AALA B   2       0.000   0.000   0.000  0.00  0.00 0 0    1
""")
  #
  expected = """\
MODEL        0
ATOM      1  C   MET A   1       0.000   0.000   0.000  0.00  0.00
ATOM      2  CA AMET A   1       0.000   0.000   0.000  0.00  0.00
BREAK
ATOM      3  N   GLY A   2       0.000   0.000   0.000  0.00  0.00
ATOM      4  CA AGLY A   2       0.000   0.000   0.000  0.00  0.00
TER
ENDMDL
MODEL        1
ATOM      1  N   MET A   1       0.000   0.000   0.000  0.00  0.00
ATOM      2  CA AMET A   1       0.000   0.000   0.000  0.00  0.00
TER
ATOM      3  N   ALA B   2       0.000   0.000   0.000  0.00  0.00
ATOM      4  CA AALA B   2       0.000   0.000   0.000  0.00  0.00
TER
ENDMDL
"""
  for append_end in [False, True]:
    if (append_end): expected += "END\n"
    assert not show_diff(pdb_inp.as_pdb_string(append_end=append_end),expected)
    pdb_inp.write_pdb_file(file_name="tmp_tst_hierarchy.pdb", append_end=append_end)
    assert not show_diff(open("tmp_tst_hierarchy.pdb").read(), expected)
  #
  lines = flex.split_lines("""\
MODEL     SKDI
HETATMB1234 NaMeLResChUvwqI      1.300   2.100   3.200  0.40  4.80      sEgIElcH
ENDMDL
""")
  pdb_inp = pdb.input(source_info=None, lines=lines)
  hierarchy = pdb_inp.construct_hierarchy()
  awl = hierarchy.atoms()[0].fetch_labels()
  assert not show_diff(awl.format_atom_record(), lines[1])
  assert awl.model_id == "SKDI"
  assert not awl.is_first_in_chain
  assert not awl.is_first_after_break

def exercise_transfer_chains_from_other():
  atoms_x = """\
ATOM      1  X1      A
ATOM      2  X2      A
ATOM      3  X1      B
ATOM      4  X2      B"""
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
MODEL        1
%s
ENDMDL
MODEL        2
%s
ENDMDL
""" % (atoms_x, atoms_x.replace("X","Y"))))
  hierarchy = pdb_inp.construct_hierarchy()
  models = hierarchy.models()
  assert [md.chains_size() for md in models] == [2, 2]
  models[0].transfer_chains_from_other(other=models[1])
  assert [md.chains_size() for md in models] == [4, 0]
  hierarchy.remove_model(model=models[1])
  trailing = "           0.000   0.000   0.000  0.00  0.00"
  check_wpf(hierarchy, trailing=trailing, expected="""\
ATOM      1  X1      A
ATOM      2  X2      A
TER
ATOM      3  X1      B
ATOM      4  X2      B
TER
ATOM      1  Y1      A
ATOM      2  Y2      A
TER
ATOM      3  Y1      B
ATOM      4  Y2      B
TER
""")
  #
  for suffixes in [Auto, "FG"]:
    roots = [
      pdb.input(source_info=None, lines=flex.split_lines(lines))
        .construct_hierarchy()
          for lines in [atoms_x, atoms_x.replace("X","Y")]]
    joined = pdb.hierarchy.join_roots(
      roots=roots, chain_id_suffixes=suffixes)
    if (suffixes is Auto): f,g = "1", "2"
    else:                  f,g = "F", "G"
    trailing = "           0.000   0.000   0.000  0.00  0.00"
    check_wpf(joined, trailing=trailing, expected="""\
ATOM      1  X1     A%s
ATOM      2  X2     A%s
TER
ATOM      3  X1     B%s
ATOM      4  X2     B%s
TER
ATOM      1  Y1     A%s
ATOM      2  Y2     A%s
TER
ATOM      3  Y1     B%s
ATOM      4  Y2     B%s
TER
""" % (f,f,f,f,g,g,g,g))
  #
  roots = [pdb_inp.construct_hierarchy() for pdb_inp in [
    pdb.input(source_info=None, lines=flex.split_lines("""\
MODEL        1
%s
ENDMDL
MODEL        1
%s
ENDMDL
MODEL        1
%s
ENDMDL
""" % (atoms_x, atoms_x.replace("X","Y"), atoms_x.replace("X","Z")))),
    pdb.input(source_info=None, lines=flex.split_lines("""\
MODEL        3
%s
ENDMDL
MODEL        3
%s
ENDMDL
""" % (atoms_x.replace("X","P"), atoms_x.replace("X","Q"))))]]
  joined = pdb.hierarchy.join_roots(roots=roots)
  trailing = "           0.000   0.000   0.000  0.00  0.00"
  check_wpf(joined, trailing=trailing, expected="""\
MODEL        1
ATOM      1  X1     A1
ATOM      2  X2     A1
TER
ATOM      3  X1     B1
ATOM      4  X2     B1
TER
ATOM      1  P1     A2
ATOM      2  P2     A2
TER
ATOM      3  P1     B2
ATOM      4  P2     B2
TER
ENDMDL
MODEL        2
ATOM      1  Y1     A1
ATOM      2  Y2     A1
TER
ATOM      3  Y1     B1
ATOM      4  Y2     B1
TER
ATOM      1  Q1     A2
ATOM      2  Q2     A2
TER
ATOM      3  Q1     B2
ATOM      4  Q2     B2
TER
ENDMDL
MODEL        3
ATOM      1  Z1     A1
ATOM      2  Z2     A1
TER
ATOM      3  Z1     B1
ATOM      4  Z2     B1
TER
ENDMDL
""")

def exercise_root_select(n_trials=100):
  h_all = pdb.input(source_info=None, lines=flex.split_lines("""\
MODEL        1
ATOM      1
ATOM      2
ATOM      3
ATOM      4
ATOM      5
ATOM      6
ATOM      7
ATOM      8
ENDMDL
MODEL        2
ATOM     11
ATOM     12
ATOM     13
ATOM     14
ATOM     15
ATOM     16
ATOM     17
ATOM     18
ENDMDL
""")).construct_hierarchy()
  try: h_all.select(atom_selection=flex.bool())
  except (ValueError, RuntimeError), e:
    assert str(e) == "atom_selection array too short."
  else: raise Exception_expected
  try: h_all.select(atom_selection=flex.bool(17, False))
  except (ValueError, RuntimeError), e:
    assert str(e) == "atom_selection array too large."
  else: raise Exception_expected
  try: h_all.select(atom_selection=flex.size_t([1,0]))
  except (ValueError, RuntimeError), e:
    assert str(e) == "atom_selection indices not in strictly ascending order."
  else: raise Exception_expected
  try: h_all.select(atom_selection=flex.size_t([16]))
  except (ValueError, RuntimeError), e:
    assert str(e) \
        == "atom_selection indices greater than or equal to number of atoms."
  else: raise Exception_expected
  for atom_selections in [[flex.bool(16, False), flex.bool()],
                          [flex.size_t(), flex.size_t()]]:
    h_sel = h_all.select(atom_selection=atom_selections[0])
    assert h_sel.models_size() == 0
    assert h_sel.select(atom_selection=atom_selections[1]).models_size() == 0
  try: h_sel.select(atom_selection=flex.bool(1, False))
  except (ValueError, RuntimeError), e:
    assert str(e) == "atom_selection array too large."
  else: raise Exception_expected
  try: h_sel.select(atom_selection=flex.size_t([0]))
  except (ValueError, RuntimeError), e:
    assert str(e) \
        == "atom_selection indices greater than or equal to number of atoms."
  else: raise Exception_expected
  for atom_selection in [flex.bool(16, True),
                         flex.size_t_range(16)]:
    h_sel = h_all.select(atom_selection=atom_selection)
    assert h_sel.is_similar_hierarchy(other=h_all)
    assert h_all.is_similar_hierarchy(other=h_sel)
  for atom_selection in [flex.bool([True]*8+[False]*8),
                         flex.size_t_range(8)]:
    h_sel = h_all.select(atom_selection=atom_selection)
    assert h_sel.models_size() == 1
    assert h_sel.only_model().is_identical_hierarchy(
      other=h_all.models()[0])
    assert not h_sel.only_model().is_identical_hierarchy(
      other=h_all.models()[1])
  for atom_selection in [flex.bool([False]*8+[True]*8),
                         flex.size_t_range(8,16)]:
    h_sel = h_all.select(atom_selection=atom_selection)
    assert h_sel.models_size() == 1
    assert not h_sel.only_model().is_identical_hierarchy(
      other=h_all.models()[0])
    assert h_sel.only_model().is_identical_hierarchy(
      other=h_all.models()[1])
  for atom_selection in [flex.bool([False,True]*8),
                         flex.size_t_range(1,16,2)]:
    h_sel = h_all.select(atom_selection=atom_selection)
    assert [a.serial.lstrip() for a in h_sel.atoms()] \
        == ["2", "4", "6", "8", "12", "14", "16", "18"]
  for atom_selection in [flex.bool([True,False]*8),
                         flex.size_t_range(0,16,2)]:
    h_sel = h_all.select(atom_selection=atom_selection)
    assert [a.serial.lstrip() for a in h_sel.atoms()] \
        == ["1", "3", "5", "7", "11", "13", "15", "17"]
  a_all = h_all.atoms()
  a_all.reset_i_seq()
  sentinel = a_all.reset_tmp(first_value=1, increment=0)
  for copy_atoms in [False, True]:
    for i_trial in xrange(n_trials):
      sel = flex.random_bool(size=16, threshold=0.5)
      a_sel = a_all.select(sel)
      h_sel = h_all.select(sel, copy_atoms=copy_atoms)
      assert h_sel.atoms_size() == sel.count(True)
      def check(h_sel):
        for a,b in zip(a_sel, h_sel.atoms()):
          assert a.serial == b.serial
          if (copy_atoms):
            assert a.memory_id() != b.memory_id()
            assert a.tmp == 1
            assert b.tmp == 0
            assert b.i_seq == 0
          else:
            assert a.memory_id() == b.memory_id()
      check(h_sel)
      h_isel = h_all.select(sel.iselection(), copy_atoms=copy_atoms)
      check(h_isel)
      assert h_isel.is_similar_hierarchy(h_sel)
      assert h_sel.is_similar_hierarchy(h_isel)
      assert h_sel.as_pdb_string() == h_isel.as_pdb_string()

def exercise_root_altloc_indices(n_trials=10):
  lines = flex.split_lines("""\
ATOM      1  N
ATOM      2  CB A
ATOM      3  CA
ATOM      4  CB B
ATOM      5  CG
ATOM      6  CG A
ATOM      7  CD A
ATOM      8  C
ATOM      9  CD B
ATOM     10  O
""")
  for i_trial in xrange(n_trials):
    pdb_inp = pdb.input(
      source_info=None,
      lines=lines.select(flex.random_permutation(size=lines.size())))
    hierarchy = pdb_inp.construct_hierarchy()
    atoms = hierarchy.atoms()
    ai = dict([(k,sorted([int(atoms[i].serial) for i in v]))
      for k,v in hierarchy.altloc_indices().items()])
    assert ai == {
      "": [1, 3, 8, 10],
      " ": [5],
      "A": [2, 6, 7],
      "B": [4, 9]}

def exercise_root_pickling():
  pdb_inp = pdb.input(source_info=None, lines="""\
MODEL        1
ATOM      1  N   MET A   1       6.215  22.789  24.067  1.00  0.00           N
ATOM      2  CA  MET A   1       6.963  22.789  22.822  1.00  0.00           C
BREAK
HETATM    3  C   MET A   2       7.478  21.387  22.491  1.00  0.00           C
ATOM      4  O   MET A   2       8.406  20.895  23.132  1.00  0.00           O
ENDMDL
MODEL 3
HETATM    9 2H3  MPR B   5      16.388   0.289   6.613  1.00  0.08
SIGATM    9 2H3  MPR B   5       0.155   0.175   0.155  0.00  0.05
ANISOU    9 2H3  MPR B   5      848    848    848      0      0      0
SIGUIJ    9 2H3  MPR B   5      510    510    510      0      0      0
TER
ATOM     10  N   CYSCH   6      14.270   2.464   3.364  1.00  0.07
SIGATM   10  N   CYSCH   6       0.012   0.012   0.011  0.00  0.00
ANISOU   10  N   CYSCH   6      788    626    677   -344    621   -232
SIGUIJ   10  N   CYSCH   6        3     13      4     11      6     13
TER
ENDMDL
""")
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.info.append("a")
  hierarchy.info.append("b")
  for p in [pickle, cPickle]:
    s = p.dumps(hierarchy, 1)
    l = p.loads(s)
    assert not show_diff("\n".join(l.info), "\n".join(hierarchy.info))
    assert not show_diff(l.as_pdb_string(), hierarchy.as_pdb_string())

def exercise_residue_pickling():
  pdb_inp = pdb.input(source_info=None, lines="""\
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ANISOU    2  CA  GLY A   1      788    626    677   -344    621   -232       C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
SIGATM    9  N   ASN A   3       0.012   0.012   0.011  0.00  0.00           N
ATOM     10  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     11  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     12  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
""")
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  conformer = hierarchy.only_conformer()
  for residue in conformer.residues():
    eps = hierarchy.select(residue.atoms().extract_i_seq()).as_pdb_string()
    rsc = residue.standalone_copy()
    assert not show_diff(rsc.root().as_pdb_string(), eps)
    for p in [pickle, cPickle]:
      for rp in [residue, rsc]:
        for i_pass in xrange(2):
          s = p.dumps(rp, 1)
          l = p.loads(s)
          assert not show_diff(l.root().as_pdb_string(), eps)
          rp = l

def exercise_hierarchy_input():
  pdb_obj = pdb.hierarchy.input(pdb_string=pdb_2izq_220)
  i_atoms = pdb_obj.input.atoms()
  h_atoms = pdb_obj.hierarchy.atoms()
  assert i_atoms.size() == 70
  assert h_atoms.size() == 70
  h_i_perm = pdb_obj.hierarchy_to_input_atom_permutation()
  i_h_perm = pdb_obj.input_to_hierarchy_atom_permutation()
  def check(atoms, atoms_perm):
    for a,ap in zip(atoms, atoms_perm):
      assert ap.memory_id() == a.memory_id()
  check(h_atoms.select(h_i_perm), i_atoms)
  check(i_atoms.select(i_h_perm), h_atoms)

def exercise_other () :
  # XXX Nat's utility functions
  pdb_inp = pdb.input(source_info=None, lines="""\
CRYST1    2.000    3.000    4.000  90.00  80.00  90.00 P 2           5
ATOM      0  S   SO4     0       3.302   8.419   8.560  1.00 10.00           S
ANISOU    0  S   SO4     0     1000   2000   3000    400   -500    600
ATOM      1  O1  SO4     0       3.497   8.295   7.118  1.00 10.00           O
ATOM      4  O3  SO4     0       4.481   9.037   9.159  1.00 20.00           O
ATOM      5  O4  SO4     0       2.131   9.251   8.823  1.00 30.00           O
ATOM      2  O2 ASO4     0       3.098   7.095   9.140  0.80 40.00           O
ATOM      3  O2 BSO4     0       3.498   7.495   9.440  0.20 50.00           O
TER
END
""")
  hierarchy = pdb_inp.construct_hierarchy()
  xray_structure = hierarchy.extract_xray_structure()
  assert xray_structure.sites_cart().size() == hierarchy.atoms().size()
  xray_structure.scale_adps(2.0)
  atoms = hierarchy.atoms()
  atoms.set_adps_from_scatterers(xray_structure.scatterers(),
    xray_structure.unit_cell())
  assert approx_equal(atoms.extract_uij(),
    [(0.2, 0.4, 0.6, 0.08, -0.1, 0.12),
     (-1,-1,-1,-1,-1,-1), (-1,-1,-1,-1,-1,-1), (-1,-1,-1,-1,-1,-1),
     (-1,-1,-1,-1,-1,-1), (-1,-1,-1,-1,-1,-1)])
  assert approx_equal(atoms.extract_b(),
    [31.5827341, 20.0, 40.0, 60.0, 80.0, 100.0])
  xray_structure.convert_to_isotropic()
  atoms.set_adps_from_scatterers(xray_structure.scatterers(),
    xray_structure.unit_cell())
  assert approx_equal(atoms.extract_uij(), [(-1,-1,-1,-1,-1,-1)]*6)
  assert approx_equal(atoms.extract_b(),
    [31.5827341, 20.0, 40.0, 60.0, 80.0, 100.0])
  pdb_inp = pdb.input(source_info=None, lines="""\
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ANISOU    2  CA  GLY A   1      788    626    677   -344    621   -232       C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN B   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN B   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN B   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN B   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  N   ASN B   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     10  CA  ASN B   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     11  C   ASN B   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     12  O   ASN B   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM      0  S   SO4 C   4       3.302   8.419   8.560  1.00 10.00           S
ATOM      1  O1  SO4 C   4       3.497   8.295   7.118  1.00 10.00           O
ATOM      4  O3  SO4 C   4       4.481   9.037   9.159  1.00 10.00           O
ATOM      5  O4  SO4 C   4       2.131   9.251   8.823  1.00 10.00           O
""")
  hierarchy = pdb_inp.construct_hierarchy()
  pdb_inp_new = pdb.input(source_info=None, lines="""\
ATOM      5  N   ASN B   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN B   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN B   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN B   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  N  AASN B   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     10  CA AASN B   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     11  C  AASN B   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     12  O  AASN B   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM      9  N  BASN B   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     10  CA BASN B   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     11  C  BASN B   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     12  O  BASN B   3      -1.872   0.119   3.648  1.00 10.42           O
""")
  h2 = pdb_inp_new.construct_hierarchy()
  chain_b = h2.models()[0].chains()[0]
  partial_hierarchy = pdb.hierarchy.new_hierarchy_from_chain(chain_b)
  assert not show_diff(partial_hierarchy.as_pdb_string(), h2.as_pdb_string())
  pdb.hierarchy.find_and_replace_chains(hierarchy, partial_hierarchy)
  assert not show_diff(hierarchy.as_pdb_string(), """\
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ANISOU    2  CA  GLY A   1      788    626    677   -344    621   -232       C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
TER
ATOM      5  N   ASN B   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN B   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN B   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN B   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  N  AASN B   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     10  CA AASN B   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     11  C  AASN B   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     12  O  AASN B   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM      9  N  BASN B   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     10  CA BASN B   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     11  C  BASN B   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     12  O  BASN B   3      -1.872   0.119   3.648  1.00 10.42           O
TER
ATOM      0  S   SO4 C   4       3.302   8.419   8.560  1.00 10.00           S
ATOM      1  O1  SO4 C   4       3.497   8.295   7.118  1.00 10.00           O
ATOM      4  O3  SO4 C   4       4.481   9.037   9.159  1.00 10.00           O
ATOM      5  O4  SO4 C   4       2.131   9.251   8.823  1.00 10.00           O
TER
""")
  # distance-based connectivity
  pdb_hierarchy = pdb.input(source_info=None, lines="""\
ATOM      1  C01 LIG A   1      -2.986   0.015   1.643  1.00 20.00      A    C
ATOM      2  N02 LIG A   1      -1.545   0.015   1.643  1.00 20.00      A    N
ATOM      3  C03 LIG A   1      -0.733   0.015   2.801  1.00 20.00      A    C
ATOM      4  N04 LIG A   1       0.593   0.015   2.395  1.00 20.00      A    N
ATOM      5  C05 LIG A   1       0.618   0.015   1.034  1.00 20.00      A    C
ATOM      6  N06 LIG A   1       1.758   0.015   0.102  1.00 20.00      A    N
ATOM      7  C07 LIG A   1       3.092  -0.060   0.694  1.00 20.00      A    C
ATOM      8  C08 LIG A   1       1.525   0.015  -1.360  1.00 20.00      A    C
ATOM      9  O09 LIG A   1       2.489  -0.024  -2.139  1.00 20.00      A    O
ATOM     10  N10 LIG A   1       0.158   0.015  -1.888  1.00 20.00      A    N
ATOM     11  C11 LIG A   1      -0.025   0.024  -3.330  1.00 20.00      A    C
ATOM     12  C12 LIG A   1      -0.986   0.015  -0.959  1.00 20.00      A    C
ATOM     13  O13 LIG A   1      -2.155   0.008  -1.408  1.00 20.00      A    O
ATOM     14  C14 LIG A   1      -0.733   0.015   0.565  1.00 20.00      A    C
ATOM     15 H011 LIG A   1      -3.346   0.016   2.662  1.00 20.00      A    H
ATOM     16 H012 LIG A   1      -3.347   0.896   1.133  1.00 20.00      A    H
ATOM     17 H013 LIG A   1      -3.347  -0.868   1.136  1.00 20.00      A    H
ATOM     18 H031 LIG A   1      -1.083   0.020   3.822  1.00 20.00      A    H
ATOM     19 H071 LIG A   1       3.184  -0.975   1.260  1.00 20.00      A    H
ATOM     20 H072 LIG A   1       3.245   0.785   1.348  1.00 20.00      A    H
ATOM     21 H073 LIG A   1       3.835  -0.047  -0.090  1.00 20.00      A    H
ATOM     22 H111 LIG A   1       0.508   0.861  -3.756  1.00 20.00      A    H
ATOM     23 H112 LIG A   1      -1.076   0.113  -3.560  1.00 20.00      A    H
ATOM     24 H113 LIG A   1       0.358  -0.896  -3.748  1.00 20.00      A    H
""").construct_hierarchy()
  bonds = pdb_hierarchy.distance_based_simple_two_way_bond_sets()
  assert bonds.size() == pdb_hierarchy.atoms().size()
  #print list(bonds[0])
  assert list(bonds[0]) == [1, 14, 15, 16]
  pdb_hierarchy = pdb.input(source_info=None, lines="""\
ATOM      1  C1  EOH     1       3.108   0.653  -8.526  1.00  0.00           C
ATOM      2  C2  EOH     1       4.597   0.674  -8.132  1.00  0.00           C
ATOM      3 1H1  EOH     1       2.815  -0.349  -8.761  1.00  0.00           H
ATOM      4 2H1  EOH     1       2.517   1.015  -7.711  1.00  0.00           H
ATOM      5 3H1  EOH     1       2.956   1.278  -9.381  1.00  0.00           H
ATOM      6 1H2 AEOH     1       5.210   0.503  -9.017  1.00  0.00           H
ATOM      7 2H2 AEOH     1       4.790  -0.110  -7.400  1.00  0.00           H
ATOM      8  OH AEOH     1       4.922   1.945  -7.565  1.00  0.00           O
ATOM      9  HH AEOH     1       5.850   1.958  -7.320  1.00  0.00           H
ATOM     10 1H2 BEOH     1       5.198   0.305  -8.963  1.00  0.00           H
ATOM     11 2H2 BEOH     1       4.751   0.037  -7.261  1.00  0.00           H
ATOM     12  OH BEOH     1       4.988   2.012  -7.818  1.00  0.00           O
ATOM     13  HH BEOH     1       5.916   2.025  -7.573  1.00  0.00           H
""").construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  bonds = pdb_hierarchy.distance_based_simple_two_way_bond_sets()
  assert (list(bonds[7]) == [1,8])
  assert (list(bonds[10]) == [1])
  assert (list(bonds[11]) == [1,12])
  assert (list(bonds[12]) == [11])

  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/2hr0.pdb", test=os.path.isfile)
  if (pdb_file is not None) :
    (n_res,n_frag,n_wat)=pdb.hierarchy.get_residue_and_fragment_count(pdb_file)
    assert (n_res == 1548) and (n_frag == 7) and (n_wat == 708)
  #
  pdb_hierarchy = pdb.input(source_info=None, lines="""\
HETATM    0  PA  ANP A   1      27.977  39.209  67.441  1.00 26.15           P
HETATM    1  O2A ANP A   1      28.552  39.588  68.829  1.00 27.12           O
HETATM    2 MN    MN B   1      28.911  38.079  64.440  1.00 24.79          Mn
""").construct_hierarchy()
  bond_sets = pdb_hierarchy.distance_based_simple_two_way_bond_sets()
  assert [list(bond) for bond in bond_sets] == [[1], [0], []]
  # the following tests check order of xray scatterers vs. hierarchy atoms
  pdb_str = """\
ATOM      1  C1  EOH     1       3.108   0.653  -8.526  1.00  0.00           C
ATOM      2  C2  EOH     1       4.597   0.674  -8.132  1.00  0.00           C
ATOM      3 1H1  EOH     1       2.815  -0.349  -8.761  1.00  0.00           H
ATOM      4 2H1  EOH     1       2.517   1.015  -7.711  1.00  0.00           H
ATOM      5 3H1  EOH     1       2.956   1.278  -9.381  1.00  0.00           H
ATOM      6 1H2 AEOH     1       5.210   0.503  -9.017  1.00  0.00           H
ATOM      7 1H2 BEOH     1       5.198   0.305  -8.963  1.00  0.00           H
ATOM      8 2H2 AEOH     1       4.790  -0.110  -7.400  1.00  0.00           H
ATOM      9 2H2 BEOH     1       4.751   0.037  -7.261  1.00  0.00           H
ATOM     10  OH AEOH     1       4.922   1.945  -7.565  1.00  0.00           O
ATOM     11  OH BEOH     1       4.988   2.012  -7.818  1.00  0.00           O
ATOM     12  HH AEOH     1       5.850   1.958  -7.320  1.00  0.00           H
ATOM     13  HH BEOH     1       5.916   2.025  -7.573  1.00  0.00           H
"""
  pdb_in = pdb.hierarchy.input(pdb_string=pdb_str)
  xrs = pdb_in.xray_structure_simple()
  assert (xrs.sites_cart().size() == 13)
  assert (approx_equal(xrs.sites_cart()[-3][-1], -7.261, eps=0.0001))
  pdb_in = pdb.input(source_info=None, lines=pdb_str)
  hierarchy = pdb_in.construct_hierarchy()
  xrs = hierarchy.extract_xray_structure(
    crystal_symmetry=pdb_in.crystal_symmetry())
  assert (xrs.sites_cart().size() == 13)
  assert (approx_equal(xrs.sites_cart()[-3][-1], -7.261, eps=0.0001))
  # sequence extraction
  pdb_hierarchy = pdb.input(source_info=None, lines="""\
ATOM      2  CA  GLY A   3      -9.052   4.207   4.651  1.00 16.57           C
ATOM      6  CA  ASN A   4      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  ASN A   5      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN A   6       0.384   1.888   3.199  1.00 10.53           C
ATOM     31  CA  GLN A   7       3.270   2.361   5.640  1.00 11.39           C
ATOM     40  CA  ASN A   8       6.831   2.310   4.318  1.00 12.30           C
ATOM     48  CA  TYR A   9       9.159   2.144   7.299  1.00 15.18           C
""").construct_hierarchy()
  main_conf = pdb_hierarchy.models()[0].chains()[0].conformers()[0]
  assert (main_conf.as_padded_sequence() == "XXGNNQQNY")
  pdb_hierarchy = pdb.input(source_info=None, lines="""\
ATOM      2  CA  GLY A  -2      -9.052   4.207   4.651  1.00 16.57           C
ATOM      6  CA  ASN A  -1      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  ASN A   0      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN A   1       0.384   1.888   3.199  1.00 10.53           C
ATOM     31  CA  GLN A   4       3.270   2.361   5.640  1.00 11.39           C
ATOM     40  CA  ASN A   5       6.831   2.310   4.318  1.00 12.30           C
ATOM     48  CA  TYR A   6       9.159   2.144   7.299  1.00 15.18           C
""").construct_hierarchy()
  main_conf = pdb_hierarchy.models()[0].chains()[0].conformers()[0]
  assert (main_conf.as_padded_sequence() == "GNNQXXQNY")
  pdb_hierarchy = pdb.input(source_info=None, lines="""\
ATOM      2  CA  GLY A   3      -9.052   4.207   4.651  1.00 16.57           C
ATOM      6  CA  ASN A   4      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  ASN A   5      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN A   6       0.384   1.888   3.199  1.00 10.53           C
ATOM     22  CA  ALA A   6A      0.384   1.888   3.199  1.00 10.53           C
ATOM     22  CA  GLY A   6B      0.384   1.888   3.199  1.00 10.53           C
ATOM     31  CA  GLN A   7       3.270   2.361   5.640  1.00 11.39           C
ATOM     40  CA  ASN A   8       6.831   2.310   4.318  1.00 12.30           C
ATOM     48  CA  TYR A   9       9.159   2.144   7.299  1.00 15.18           C
""").construct_hierarchy()
  main_conf = pdb_hierarchy.models()[0].chains()[0].conformers()[0]
  assert (main_conf.as_padded_sequence() == "XXGNNQAGQNY")
  assert (main_conf.as_padded_sequence(skip_insertions=True) == "XXGNNQQNY")
  resids = main_conf.get_residue_ids()
  assert (len(resids) == 11)
  assert (resids[0] == resids[1] == None)
  assert (resids[-4].strip() == "6B")
  # sites_diff
  hierarchy_1 = pdb.input(source_info=None, lines="""\
ATOM      0  O   WAT B   1      17.523   2.521  10.381  1.10 16.78           O
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
TER
ATOM      5  O   HOH S   1      -7.523   2.521  10.381  0.10  6.78           O
""").construct_hierarchy()
  hierarchy_2 = pdb.input(source_info=None, lines="""\
ATOM      1  N   GLY A   1      -9.009   4.612   7.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   5.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   5.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -6.523   3.521   6.381  1.00 16.78           O
TER
ATOM      5  O   HOH S   1      -9.523   5.521  11.381  0.10  6.78           O
""").construct_hierarchy()
  hierarchy_new = pdb.hierarchy.sites_diff(hierarchy_1, hierarchy_2)
  assert approx_equal(hierarchy_new.atoms().extract_b(),
      [1.0, 1.0, 1.0, 1.7320508, -1.0])
  deltas = pdb.hierarchy.sites_diff(hierarchy_1, hierarchy_2,
      exclude_waters=False, return_hierarchy=False)
  assert approx_equal(deltas, [1.0, 1.0, 1.0, 1.7320508, 3.7416574])
  # show_file_summary
  out = StringIO()
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/2C30.pdb",
    test=os.path.isfile)
  if (pdb_file is None) :
    return
  pdb_in = pdb.input(file_name=pdb_file)
  info = pdb.show_file_summary(pdb_in, out=out)
  expected = """\
Number of atoms:          2691
Number of chains:         3
Chain IDs:                A, Z
Alternate conformations:  0
Amino acid residues:      289
Water molecules:          350
Elemental ions:           1 ( CL)
Other molecules:          2 (PO4, SEP)
Mean isotropic B-factor:  22.10 (range: 10.58 - 55.34)
Space group:              P 21 21 21
Unit cell:                59.781 66.674 96.999 90 90 90
"""
  assert (out.getvalue() == expected), "\n"+out.getvalue()
  # get_contiguous_ranges
  pdb_hierarchy = pdb.input(source_info=None, lines="""\
ATOM      2  CA  GLY A  -1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      6  CA  ASN A   0      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  ASN A   1      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN A   2       0.384   1.888   3.199  1.00 10.53           C
ATOM     22  CA  ALA A   2A      0.384   1.888   3.199  1.00 10.53           C
ATOM     22  CA  GLY A   2B      0.384   1.888   3.199  1.00 10.53           C
ATOM     31  CA  GLN A   4       3.270   2.361   5.640  1.00 11.39           C
ATOM     40  CA  ASN A   5       6.831   2.310   4.318  1.00 12.30           C
ATOM     48  CA  TYR A   6       9.159   2.144   7.299  1.00 15.18           C
TER
ATOM      5  O   HOH S   1      -9.523   5.521  11.381  0.10  6.78           O
""").construct_hierarchy()
  selections = pdb.hierarchy.get_contiguous_ranges(pdb_hierarchy)
  assert (selections == [
   "chain 'A' and ((resid   -1  through    2B) or (resid    4  through    6 ))",
   "chain 'S' and ((resid    1 ))"])

def exercise_equality_and_hashing():

  root = pdb.hierarchy.root()
  model = pdb.hierarchy.model()
  root.append_model( model )
  chain = pdb.hierarchy.chain()
  model.append_chain( chain )
  rg = pdb.hierarchy.residue_group()
  chain.append_residue_group( rg )
  ag = pdb.hierarchy.atom_group()
  rg.append_atom_group( ag )
  atom = pdb.hierarchy.atom()
  ag.append_atom( atom )
  collection = set( [ root, model, chain, rg, ag, atom ] )

  # Root
  pr = model.parent()
  assert pr is not root
  assert pr == root
  assert pr in collection

  # Model
  pm = chain.parent()
  assert pm is not model
  assert pm == model
  assert pm in collection

  # Chain
  pc = rg.parent()
  assert pc is not chain
  assert pc == chain
  assert pc in collection

  # Residue group
  prg = ag.parent()
  assert prg is not rg
  assert prg == rg
  assert prg in collection

  # Atom group
  pag = atom.parent()
  assert pag is not ag
  assert pag == ag
  assert pag in collection

  # Atom
  pa = ag.atoms()[0]
  assert pa is not atom
  assert pa == atom
  assert pa in collection

def exercise_atom_is_in_same_conformer_as():
  pdb_hierarchy = pdb.input(source_info=None, lines="""\
MODEL        1
ATOM      0  N   MET
ATOM      1  CA AMET
ATOM      2  CA BMET
ENDMDL
MODEL        2
ATOM      3  N   MET
ATOM      4  CA  MET
ATOM      5  CA BMET
ENDMDL
END
""").construct_hierarchy()
  atoms = pdb_hierarchy.atoms()
  for first in [0,3]:
    for i in xrange(3):
      assert atoms[first].is_in_same_conformer_as(atoms[first+i])
    assert not atoms[first+1].is_in_same_conformer_as(atoms[first+2])
  for i in xrange(3):
    for j in xrange(3,6):
      assert not atoms[i].is_in_same_conformer_as(atoms[j])

def get_phenix_regression_pdb_file_names():
  pdb_dir = libtbx.env.find_in_repositories("phenix_regression/pdb")
  if (pdb_dir is None): return None
  result = []
  for node in os.listdir(pdb_dir):
    if (not (node.endswith(".pdb") or node.endswith(".ent"))): continue
    result.append(os.path.join(pdb_dir, node))
  assert len(result) != 0
  return result

def exercise_adopt_xray_structure():
  from cctbx import adptbx
  pdb_inp = pdb.input(source_info=None, lines="""\
CRYST1   12.000   13.000   14.000  80.00  90.00 100.00 P 1
ATOM      0  O   WAT B   1      17.523   2.521  10.381  1.10 16.78           O
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ANISOU    2  CA  GLY A   1      788    626    677   -344    621   -232       C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
TER
ATOM      5  O   HOH S   1      -7.523   2.521  10.381  0.10  6.78           O
TER
ATOM      6  O   HOH     1      10.523   2.521   5.381  1.00 16.78           O
ANISOU    6  O   HOH     1      788    626    677   -344    621   -232       O
""")
  hierarchy = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  xrs_new1 = xrs.deep_copy_scatterers()
  xrs_new1.shake_adp()
  xrs_new1.shake_occupancies()
  xrs_new1.shake_sites_in_place(rms_difference=0.5)
  hierarchy.adopt_xray_structure(xray_structure = xrs_new1)
  xrs_new2 = hierarchy.extract_xray_structure(crystal_symmetry =
    xrs.crystal_symmetry())
  uc = xrs.unit_cell()
  orth = uc.orthogonalize
  for s1,s2 in zip(xrs_new1.scatterers(), xrs_new2.scatterers()):
    assert approx_equal(orth(s1.site), orth(s2.site), 1.e-3)
    assert approx_equal(adptbx.u_as_b(s1.u_iso), adptbx.u_as_b(s2.u_iso), 1.e-2)
    assert approx_equal(s1.occupancy, s2.occupancy, 1.e-2)
    assert approx_equal(s1.u_star, s2.u_star)
    assert s1.scattering_type == s2.scattering_type
  xrs_new3 = xrs.concatenate(other = xrs_new1)
  try: hierarchy.adopt_xray_structure(xray_structure = xrs_new3)
  except RuntimeError, e: pass
  assert str(e) == "Incompatible size of hierarchy and scatterers array."

def exercise(args):
  comprehensive = "--comprehensive" in args
  forever = "--forever" in args
  print "iotbx.pdb.hierarchy.atom.sizeof_data():", \
    pdb.hierarchy.atom.sizeof_data()
  offsets = pdb.hierarchy.atom.data_offsets()
  if (comprehensive):
    print "iotbx.pdb.hierarchy.atom.data_offsets():"
    prev = 0
    for key,value in sorted(offsets.items()):
      print "  %+3d %3d %s" % (key-prev, key, value)
      prev = key
  phenix_regression_pdb_file_names = get_phenix_regression_pdb_file_names()
  while True:
    exercise_adopt_xray_structure()
    exercise_atom()
    exercise_atom_group()
    exercise_residue_group()
    exercise_chain()
    exercise_model()
    exercise_root()
    exercise_atom_id_str()
    exercise_format_atom_record()
    exercise_construct_hierarchy()
    exercise_convenience_generators()
    exercise_only()
    exercise_merge_atom_groups()
    exercise_merge_residue_groups()
    exercise_chain_merge_residue_groups()
    exercise_edit_blank_altloc()
    exercise_find_pure_altloc_ranges()
    exercise_occupancy_groups_simple()
    exercise_conformers()
    exercise_residue()
    exercise_is_identical_hierarchy()
    exercise_is_similar_hierarchy()
    exercise_atoms()
    exercise_atoms_interleaved_conf()
    exercise_as_pdb_string(
      pdb_file_names=phenix_regression_pdb_file_names,
      comprehensive=comprehensive)
    exercise_atom_with_labels()
    exercise_transfer_chains_from_other()
    exercise_root_select()
    exercise_root_altloc_indices()
    exercise_root_pickling()
    exercise_residue_pickling()
    exercise_hierarchy_input()
    exercise_other()
    exercise_equality_and_hashing()
    exercise_atom_is_in_same_conformer_as()
    if (not forever): break
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
