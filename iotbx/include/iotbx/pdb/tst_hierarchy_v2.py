from iotbx import pdb
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
from libtbx.str_utils import show_string
from libtbx.utils import format_cpu_times
from cStringIO import StringIO
import sys

def exercise_atom():
  a = pdb.hierarchy_v2.atom()
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
  assert a.siguij == (-1,-1,-1,-1,-1,-1)
  assert not a.siguij_is_defined()
  a.siguij = (-2,3,4,-5,6,1)
  assert a.siguij == (-2,3,4,-5,6,1)
  assert a.siguij_is_defined()
  assert not a.hetero
  a.hetero = True
  assert a.hetero
  assert a.tmp == 0
  a.tmp = 3
  assert a.tmp == 3
  #
  a = (pdb.hierarchy_v2.atom()
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
    .set_siguij(new_siguij=(.1,.2,.3,.6,.1,.9))
    .set_hetero(new_hetero=True))
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
  assert approx_equal(ac.siguij, (.1,.2,.3,.6,.1,.9))
  assert ac.hetero
  #
  assert a.determine_chemical_element_simple() is None
  a.name = "NA  "
  a.element = " N"
  assert a.determine_chemical_element_simple() == " N"
  a.element = "CU"
  assert a.determine_chemical_element_simple() == "CU"
  a.element = "  "
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
  ag = pdb.hierarchy_v2.atom_group()
  ac = pdb.hierarchy_v2.atom(parent=ag, other=a)
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
  assert ac.siguij == a.siguij
  assert ac.hetero == a.hetero
  assert ac.tmp == 0

def exercise_atom_group():
  ag = pdb.hierarchy_v2.atom_group()
  assert ag.altloc == ""
  assert ag.resname == ""
  ag = pdb.hierarchy_v2.atom_group(altloc=None, resname=None)
  assert ag.altloc == ""
  assert ag.resname == ""
  ag = pdb.hierarchy_v2.atom_group(altloc="a", resname="xyz")
  assert ag.altloc == "a"
  assert ag.resname == "xyz"
  ag.altloc = None
  ag.resname = None
  assert ag.altloc == ""
  assert ag.resname == ""
  #
  ag.append_atom(atom=pdb.hierarchy_v2.atom().set_name(new_name="n"))
  rg = pdb.hierarchy_v2.residue_group()
  for i,agc in enumerate([
                 pdb.hierarchy_v2.atom_group(parent=rg, other=ag),
                 ag.detached_copy()]):
    assert agc.memory_id() != ag.memory_id()
    assert ag.parent() is None
    if (i == 0):
      assert agc.parent().memory_id() == rg.memory_id()
    else:
      assert agc.parent() is None
    assert agc.atoms_size() == 1
    assert agc.atoms()[0].memory_id() != ag.atoms()[0].memory_id()
    assert agc.atoms()[0].name == "n"
    ag.append_atom(atom=pdb.hierarchy_v2.atom().set_name(new_name="o"))
    assert ag.atoms_size() == 2+i
    assert agc.atoms_size() == 1
  #
  ag = pdb.hierarchy_v2.atom_group()
  assert ag.parent() is None
  rg1 = pdb.hierarchy_v2.residue_group()
  rg2 = pdb.hierarchy_v2.residue_group()
  assert rg1.memory_id() != rg2.memory_id()
  ag = pdb.hierarchy_v2.atom_group(parent=rg1)
  assert ag.parent().memory_id() == rg1.memory_id()
  del rg1
  assert ag.parent() is None
  #
  rg1 = pdb.hierarchy_v2.residue_group()
  ag = rg1.new_atom_group(altloc="a", resname="xyz")
  assert ag.altloc == "a"
  assert ag.resname == "xyz"
  assert ag.parent().memory_id() == rg1.memory_id()
  del rg1
  assert ag.parent() is None
  #
  ag.pre_allocate_atoms(number_of_additional_atoms=2)
  assert ag.atoms_size() == 0
  assert len(ag.atoms()) == 0
  ag.append_atom(atom=pdb.hierarchy_v2.atom().set_name(new_name="ca"))
  assert ag.atoms_size() == 1
  assert len(ag.atoms()) == 1
  ag.append_atom(atom=pdb.hierarchy_v2.atom().set_name(new_name="n"))
  assert ag.atoms_size() == 2
  assert len(ag.atoms()) == 2
  assert [atom.name for atom in ag.atoms()] == ["ca", "n"]
  ag.new_atoms(number_of_additional_atoms=3)
  assert ag.atoms_size() == 5
  assert len(ag.atoms()) == 5
  for atom in ag.atoms():
    assert atom.parent().memory_id() == ag.memory_id()
  for atom in ag.atoms():
    assert atom.tmp == 0
  assert ag.reset_atom_tmp(new_value=7) == 5
  for atom in ag.atoms():
    assert atom.tmp == 7
  assert [a.name for a in ag.atoms()] == ["ca", "n", "", "", ""]
  #
  ag.insert_atom(i=0, atom=pdb.hierarchy_v2.atom().set_name(new_name="0"))
  assert [a.name for a in ag.atoms()] == ["0", "ca", "n", "", "", ""]
  ag.insert_atom(i=-1, atom=pdb.hierarchy_v2.atom().set_name(new_name="x"))
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
  a = pdb.hierarchy_v2.atom().set_name(new_name="y")
  assert ag.find_atom_index(atom=a) == -1
  try: ag.find_atom_index(atom=a, must_be_present=True)
  except RuntimeError, e:
    assert str(e) == "atom not in atom_group."
  else: raise Exception_expected
  ag.insert_atom(i=4, atom=a)
  assert ag.find_atom_index(atom=a) == 4
  assert [a.name for a in ag.atoms()] == ["0", "n", "", "x", "y"]
  #
  try: pdb.hierarchy_v2.atom_group(altloc="ab")
  except (ValueError, RuntimeError), e:
    assert str(e) == "string is too long for target variable " \
      "(maximum length is 1 character, 2 given)."
  else: raise Exception_expected

def exercise_residue_group():
  rg = pdb.hierarchy_v2.residue_group()
  assert rg.resseq == ""
  assert rg.icode == ""
  assert rg.link_to_previous
  rg = pdb.hierarchy_v2.residue_group(
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
  #
  ag = pdb.hierarchy_v2.atom_group(altloc="a")
  assert ag.parent() is None
  rg.append_atom_group(atom_group=ag)
  assert ag.parent().memory_id() == rg.memory_id()
  c = pdb.hierarchy_v2.chain()
  for i,rgc in enumerate([
                 pdb.hierarchy_v2.residue_group(parent=c, other=rg),
                 rg.detached_copy()]):
    assert rgc.memory_id() != rg.memory_id()
    assert rg.parent() is None
    if (i == 0):
      assert rgc.parent().memory_id() == c.memory_id()
    else:
      assert rgc.parent() is None
    assert rgc.atom_groups_size() == 1
    assert rgc.atom_groups()[0].memory_id() != rg.atom_groups()[0].memory_id()
    assert rgc.atom_groups()[0].altloc == "a"
    rg.append_atom_group(atom_group=pdb.hierarchy_v2.atom_group(altloc="%d"%i))
    assert rg.atom_groups_size() == 2+i
    assert rgc.atom_groups_size() == 1
    assert [ag.altloc for ag in rg.atom_groups()] == ["a", "0", "1"][:i+2]
  #
  c1 = pdb.hierarchy_v2.chain(id="a")
  c2 = pdb.hierarchy_v2.chain(id="b")
  assert c1.memory_id() != c2.memory_id()
  rg = pdb.hierarchy_v2.residue_group()
  assert rg.parent() is None
  rg = pdb.hierarchy_v2.residue_group(parent=c1)
  assert rg.parent().memory_id() == c1.memory_id()
  del c1
  assert rg.parent() is None
  #
  c1 = pdb.hierarchy_v2.chain(id="p")
  rg13l = c1.new_residue_group(resseq="13", icode="l")
  assert rg13l.resseq == "13"
  assert rg13l.icode == "l"
  #
  c1 = pdb.hierarchy_v2.chain(id="a")
  c1.pre_allocate_residue_groups(number_of_additional_residue_groups=2)
  assert c1.residue_groups_size() == 0
  assert len(c1.residue_groups()) == 0
  c1.new_residue_groups(number_of_additional_residue_groups=2)
  assert c1.residue_groups_size() == 2
  assert len(c1.residue_groups()) == 2
  for residue_group in c1.residue_groups():
    assert residue_group.parent().memory_id() == c1.memory_id()
  assert c1.reset_atom_tmp(new_value=8) == 0
  #
  for altloc in ["w", "v", "u"]:
    rg.insert_atom_group(
      i=0, atom_group=pdb.hierarchy_v2.atom_group(altloc=altloc))
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
  ag1 = pdb.hierarchy_v2.atom_group()
  ag2 = pdb.hierarchy_v2.atom_group()
  a = pdb.hierarchy_v2.atom()
  ag1.append_atom(atom=a)
  try: ag2.append_atom(atom=a)
  except RuntimeError, e:
    assert str(e) == "atom has another parent atom_group already."
  else: raise Exception_expected
  #
  rg = pdb.hierarchy_v2.residue_group()
  assert rg.resid() == "     "
  rg = pdb.hierarchy_v2.residue_group(resseq="1", icode="i")
  assert rg.resid() == "   1i"
  rg = pdb.hierarchy_v2.residue_group(resseq=" 1 ", icode="j")
  assert rg.resid() == "  1 j"
  rg = pdb.hierarchy_v2.residue_group(resseq="ABCD", icode="")
  assert rg.resid() == "ABCD "
  rg = pdb.hierarchy_v2.residue_group(resseq="ABCD", icode="E")
  assert rg.resid() == "ABCDE"
  #
  rg = pdb.hierarchy_v2.residue_group()
  ag = pdb.hierarchy_v2.atom_group(altloc=" ")
  rg.append_atom_group(atom_group=ag)
  assert not rg.have_conformers()
  ag = pdb.hierarchy_v2.atom_group(altloc="a")
  rg.append_atom_group(atom_group=ag)
  assert rg.have_conformers()

def exercise_chain():
  c = pdb.hierarchy_v2.chain()
  assert c.id == ""
  c = pdb.hierarchy_v2.chain(id="a")
  assert c.id == "a"
  c.id = "x"
  assert c.id == "x"
  #
  m1 = pdb.hierarchy_v2.model(id="1")
  m2 = pdb.hierarchy_v2.model(id="2")
  assert m1.memory_id() != m2.memory_id()
  c = pdb.hierarchy_v2.chain()
  assert c.parent() is None
  c = pdb.hierarchy_v2.chain(parent=m1)
  assert c.parent().memory_id() == m1.memory_id()
  del m1
  assert c.parent() is None
  #
  c = pdb.hierarchy_v2.chain()
  #
  c = pdb.hierarchy_v2.chain()
  c.pre_allocate_residue_groups(number_of_additional_residue_groups=2)
  assert c.residue_groups_size() == 0
  assert len(c.residue_groups()) == 0
  c.new_residue_groups(number_of_additional_residue_groups=2)
  assert c.residue_groups_size() == 2
  assert len(c.residue_groups()) == 2
  for residue_group in c.residue_groups():
    assert residue_group.parent().memory_id() == c.memory_id()
  assert c.reset_atom_tmp(new_value=9) == 0
  #
  c.residue_groups()[0].resseq = "ugh"
  m = pdb.hierarchy_v2.model()
  for i,cc in enumerate([
                pdb.hierarchy_v2.chain(parent=m, other=c),
                c.detached_copy()]):
    assert cc.memory_id() != c.memory_id()
    assert c.parent() is None
    if (i == 0):
      assert cc.parent().memory_id() == m.memory_id()
    else:
      assert cc.parent() is None
    assert cc.residue_groups_size() == 2
    assert cc.residue_groups()[0].memory_id() \
         != c.residue_groups()[0].memory_id()
    assert cc.residue_groups()[0].resseq == "ugh"
    c.append_residue_group(
      residue_group=pdb.hierarchy_v2.residue_group(resseq="%03d"%i))
    assert c.residue_groups_size() == 3+i
    assert cc.residue_groups_size() == 2
    assert [rg.resseq for rg in c.residue_groups()] \
        == ["ugh", "", "000", "001"][:i+3]
  #
  c.insert_residue_group(
    i=3, residue_group=pdb.hierarchy_v2.residue_group(resseq="b012"))
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
  rg1 = pdb.hierarchy_v2.residue_group()
  rg2 = pdb.hierarchy_v2.residue_group()
  ag = pdb.hierarchy_v2.atom_group()
  rg1.append_atom_group(atom_group=ag)
  try: rg2.append_atom_group(atom_group=ag)
  except RuntimeError, e:
    assert str(e) == "atom_group has another parent residue_group already."
  else: raise Exception_expected

def exercise_model():
  m = pdb.hierarchy_v2.model()
  assert m.id == ""
  m = pdb.hierarchy_v2.model(id="42")
  assert m.id == "42"
  m.id = "-23"
  assert m.id == "-23"
  #
  m = pdb.hierarchy_v2.model(id="1")
  assert m.parent() is None
  m.pre_allocate_chains(number_of_additional_chains=2)
  assert m.chains_size() == 0
  assert len(m.chains()) == 0
  ch_a = m.new_chain(id="a")
  assert ch_a.id == "a"
  assert ch_a.parent().memory_id() == m.memory_id()
  assert m.chains_size() == 1
  assert len(m.chains()) == 1
  ch_b = pdb.hierarchy_v2.chain(id="b")
  assert ch_b.id == "b"
  assert ch_b.parent() is None
  m.append_chain(chain=ch_b)
  assert m.chains_size() == 2
  chains = m.chains()
  assert len(chains) == 2
  assert chains[0].memory_id() == ch_a.memory_id()
  assert chains[1].memory_id() == ch_b.memory_id()
  m.new_chains(number_of_additional_chains=3)
  assert m.chains_size() == 5
  assert len(m.chains()) == 5
  for chain in m.chains():
    assert chain.parent().memory_id() == m.memory_id()
  assert m.reset_atom_tmp(new_value=2) == 0
  #
  r = pdb.hierarchy_v2.root()
  for i,mc in enumerate([
                pdb.hierarchy_v2.model(parent=r, other=m),
                m.detached_copy()]):
    assert mc.memory_id() != m.memory_id()
    assert m.parent() is None
    if (i == 0):
      assert mc.parent().memory_id() == r.memory_id()
    else:
      assert mc.parent() is None
    assert mc.chains_size() == 5
    assert mc.chains()[0].memory_id() != m.chains()[0].memory_id()
    assert mc.chains()[0].id == "a"
    m.append_chain(chain=pdb.hierarchy_v2.chain(id="%d"%i))
    assert m.chains_size() == 6+i
    assert mc.chains_size() == 5
    assert [c.id for c in m.chains()] \
        == ["a", "b", "", "", "", "0", "1"][:i+6]
  #
  m.insert_chain(i=-3, chain=pdb.hierarchy_v2.chain(id="3"))
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
  m1 = pdb.hierarchy_v2.model()
  m2 = pdb.hierarchy_v2.model()
  c = pdb.hierarchy_v2.chain()
  m1.append_chain(chain=c)
  try: m2.append_chain(chain=c)
  except RuntimeError, e:
    assert str(e) == "chain has another parent model already."
  else: raise Exception_expected

def exercise_root():
  r = pdb.hierarchy_v2.root()
  m = pdb.hierarchy_v2.model()
  assert m.parent() is None
  m = pdb.hierarchy_v2.model(parent=r)
  assert m.parent().memory_id() == r.memory_id()
  assert m.id == ""
  m = pdb.hierarchy_v2.model(parent=r, id="2")
  assert m.parent().memory_id() == r.memory_id()
  assert m.id == "2"
  del r
  assert m.parent() is None
  #
  r = pdb.hierarchy_v2.root()
  assert r.info.size() == 0
  r.info.append("abc")
  assert r.info.size() == 1
  r.info = flex.std_string(["a", "b"])
  assert r.info.size() == 2
  r.pre_allocate_models(number_of_additional_models=2)
  assert r.models_size() == 0
  assert len(r.models()) == 0
  m_a = r.new_model(id="3")
  assert m_a.id == "3"
  assert m_a.parent().memory_id() == r.memory_id()
  assert r.models_size() == 1
  assert len(r.models()) == 1
  m_b = pdb.hierarchy_v2.model(id="5")
  assert m_b.parent() is None
  r.append_model(model=m_b)
  assert r.models_size() == 2
  models = r.models()
  assert len(models) == 2
  assert models[0].memory_id() == m_a.memory_id()
  assert models[1].memory_id() == m_b.memory_id()
  r.new_models(number_of_additional_models=3)
  assert r.models_size() == 5
  assert len(r.models()) == 5
  for model in r.models():
    assert model.parent().memory_id() == r.memory_id()
  r.reset_atom_tmp(new_value=1) == 0
  #
  rc = r.deep_copy()
  assert rc.memory_id() != r.memory_id()
  assert rc.models_size() == 5
  assert rc.models()[0].memory_id() != r.models()[0].memory_id()
  assert rc.models()[0].id == "3"
  r.append_model(model=pdb.hierarchy_v2.model(id="7"))
  assert r.models_size() == 6
  assert rc.models_size() == 5
  assert [m.id for m in r.models()] == ["3", "5", "", "", "", "7"]
  assert [m.id for m in rc.models()] == ["3", "5", "", "", ""]
  rc.append_model(model=pdb.hierarchy_v2.model(id="8"))
  assert r.models_size() == 6
  assert rc.models_size() == 6
  assert [m.id for m in rc.models()] == ["3", "5", "", "", "", "8"]
  #
  r = rc.deep_copy()
  r.insert_model(i=4, model=pdb.hierarchy_v2.model(id="M"))
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
  r1 = pdb.hierarchy_v2.root()
  r2 = pdb.hierarchy_v2.root()
  m = pdb.hierarchy_v2.model()
  r1.append_model(model=m)
  try: r2.append_model(model=m)
  except RuntimeError, e:
    assert str(e) == "model has another parent root already."
  else: raise Exception_expected

def exercise_format_atom_record_using_parents():
  for hetero,record_name in [(False, "ATOM  "), (True, "HETATM")]:
    a = (pdb.hierarchy_v2.atom()
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
      .set_siguij(new_siguij=(.1,.2,.3,.6,.1,.9))
      .set_hetero(new_hetero=hetero))
    s = a.format_atom_record_using_parents()
    assert not show_diff(s, """\
%sB1234 NaMe                 1.300   2.100   3.200  0.40  4.80      sEgIElcH"""
      % record_name)
    ag = pdb.hierarchy_v2.atom_group(altloc="x", resname="uvw")
    ag.append_atom(atom=a)
    s = a.format_atom_record_using_parents()
    assert not show_diff(s, """\
%sB1234 NaMexuvw             1.300   2.100   3.200  0.40  4.80      sEgIElcH"""
      % record_name)
    rg = pdb.hierarchy_v2.residue_group(resseq="pqrs", icode="t")
    rg.append_atom_group(atom_group=ag)
    s = a.format_atom_record_using_parents()
    assert not show_diff(s, """\
%sB1234 NaMexuvw  pqrst      1.300   2.100   3.200  0.40  4.80      sEgIElcH"""
      % record_name)
    for chain_id in ["", "g", "hi"]:
      ch = pdb.hierarchy_v2.chain(id=chain_id)
      ch.append_residue_group(residue_group=rg)
      s = a.format_atom_record_using_parents()
      assert not show_diff(s, """\
%sB1234 NaMexuvw%2spqrst      1.300   2.100   3.200  0.40  4.80      sEgIElcH"""
      % (record_name, chain_id))

def exercise_construct_hierarchy():
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
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
END
"""))
  assert [atom.name for atom in pdb_inp.atoms_v2()] \
      == [" N  ", " CA ", " C  ", " O  ", "2H3 ", " N  "]
  sio = StringIO()
  root = pdb_inp.construct_hierarchy_v2()
  for model in root.models():
    print >> sio, "m:", show_string(model.id)
    for chain in model.chains():
      print >> sio, "  c:", show_string(chain.id)
      for residue_group in chain.residue_groups():
        print >> sio, "    rg:", \
          show_string(residue_group.resseq + residue_group.icode)
        for atom_group in residue_group.atom_groups():
          print >> sio, "      ag:", \
            show_string(atom_group.altloc), show_string(atom_group.resname)
  assert not show_diff(sio.getvalue(), """\
m: "   1"
  c: "A"
    rg: "   1 "
      ag: " " "MET"
    rg: "   2 "
      ag: " " "MET"
m: "   3"
  c: "B"
    rg: "   5 "
      ag: " " "MPR"
  c: "CH"
    rg: "   6 "
      ag: " " "CYS"
""")
  #
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
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
"""))
  lines = []
  root = pdb_inp.construct_hierarchy_v2()
  for model in root.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        lines.append("@residue_group")
        for atom_group in residue_group.atom_groups():
          lines.append("@atom_group")
          for atom in atom_group.atoms():
            lines.append(atom.format_atom_record_using_parents())
  assert not show_diff("\n".join(lines)+"\n", """\
@residue_group
@atom_group
ATOM    220  N  ATRP A  11      20.498  12.832  34.558  0.50  6.03           N
ATOM    221  CA ATRP A  11      21.094  12.032  35.602  0.50  5.24           C
ATOM    222  C  ATRP A  11      22.601  12.088  35.532  0.50  6.49           C
ATOM    223  O  ATRP A  11      23.174  12.012  34.439  0.50  7.24           O
ATOM    234  H  ATRP A  11      20.540  12.567  33.741  0.50  7.24           H
ATOM    235  HA ATRP A  11      20.771  12.306  36.485  0.50  6.28           H
@atom_group
ATOM    244  N  CPHE A  11      20.226  13.044  34.556  0.15  6.35           N
ATOM    245  CA CPHE A  11      20.950  12.135  35.430  0.15  5.92           C
ATOM    246  C  CPHE A  11      22.448  12.425  35.436  0.15  6.32           C
ATOM    247  O  CPHE A  11      22.961  12.790  34.373  0.15  6.08           O
ATOM    262  HB2CPHE A  11      21.221  10.536  34.146  0.15  7.21           H
ATOM    264  HB3CPHE A  11      21.198  10.093  35.647  0.15  7.21           H
ATOM    266  HD1CPHE A  11      19.394   9.937  32.837  0.15 10.53           H
ATOM    268  HD2CPHE A  11      18.873  10.410  36.828  0.15  9.24           H
ATOM    270  HE1CPHE A  11      17.206   9.172  32.650  0.15 12.52           H
ATOM    272  HE2CPHE A  11      16.661   9.708  36.588  0.15 11.13           H
ATOM    273  HZ CPHE A  11      15.908   9.110  34.509  0.15 13.18           H
@atom_group
ATOM    255  N  BTYR A  11      20.553  12.751  34.549  0.35  5.21           N
ATOM    256  CA BTYR A  11      21.106  11.838  35.524  0.35  5.51           C
ATOM    257  C  BTYR A  11      22.625  11.920  35.572  0.35  5.42           C
ATOM    258  O  BTYR A  11      23.299  11.781  34.538  0.35  5.30           O
ATOM    263  CD2BTYR A  11      18.463  10.012  36.681  0.35  9.08           C
ATOM    265  CE1BTYR A  11      17.195   9.960  34.223  0.35 10.76           C
ATOM    267  CE2BTYR A  11      17.100   9.826  36.693  0.35 11.29           C
ATOM    269  CZ BTYR A  11      16.546   9.812  35.432  0.35 11.90           C
ATOM    271  OH BTYR A  11      15.178   9.650  35.313  0.35 19.29           O
ATOM    274  H  BTYR A  11      20.634  12.539  33.720  0.35  6.25           H
ATOM    275  HA BTYR A  11      20.773  12.116  36.402  0.35  6.61           H
ATOM    276  HB2BTYR A  11      20.949  10.064  34.437  0.35  6.78           H
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
  root = exercise_merge_pdb_inp.construct_hierarchy_v2()
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
      print >> sio, atom.format_atom_record_using_parents()
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
  lines = []
  root = exercise_merge_pdb_inp.construct_hierarchy_v2()
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
      print >> sio, atom.format_atom_record_using_parents()
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

def get_single_chain(root):
  assert root.models_size() == 1
  model = root.models()[0]
  assert model.chains_size() == 1
  return model.chains()[0]

def exercise_chain_merge_and_split_residue_groups():
  pdb_inp = pdb.input(source_info=None, lines=flex.split_lines("""\
HEADER    HYDROLASE                               22-NOV-07   2VHL
HETATM 6362  O   HOH B2048      47.616  10.724 150.212  1.00 46.48           O
HETATM 6363  O  AHOH B2049      46.408  16.672 146.066  0.50 12.81           O
HETATM 6364  O   HOH B2050      29.343  12.806 185.898  1.00 35.57           O
HETATM 6365  O  BHOH B2049      43.786  12.615 147.734  0.50 28.43           O
HETATM 6366  O   HOH B2052      35.068  19.167 155.349  1.00 15.97           O
"""))
  chain = get_single_chain(root=pdb_inp.construct_hierarchy_v2())
  assert chain.residue_groups_size() == 5
  indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
  assert list(indices) == [1]
  assert chain.residue_groups_size() == 4
  indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
  assert indices.size() == 0
  lines = flex.split_lines("""\
HETATM 6363  O  AHOH B2049
HETATM 6364  O  ZHOH B2050
HETATM 6365  O  BHOH B2049
HETATM 6366  O  YHOH B2052
HETATM 9365  O  CHOH B2049
HETATM 9367  O  XHOH B2052
""")
  pdb_inp = pdb.input(source_info=None, lines=lines)
  chain = get_single_chain(root=pdb_inp.construct_hierarchy_v2())
  assert chain.residue_groups_size() == 6
  indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
  assert list(indices) == [0, 2]
  assert chain.residue_groups_size() == 3
  indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
  assert indices.size() == 0
  for i_trial in xrange(100):
    pdb_inp = pdb.input(
      source_info=None,
      lines=lines.select(flex.random_permutation(size=lines.size())))
    chain = get_single_chain(root=pdb_inp.construct_hierarchy_v2())
    indices = chain.merge_disconnected_residue_groups_with_pure_altloc()
    assert indices.size() <= 2
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
  chain = get_single_chain(root=pdb_inp.construct_hierarchy_v2())
  assert chain.residue_groups_size() == 2
  assert list(
    chain.split_residue_groups_with_mixed_resnames_but_only_blank_altloc()) \
      == [(0,3)]
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
  assert list(
    chain.split_residue_groups_with_mixed_resnames_but_only_blank_altloc()) \
      == []
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
  chain = get_single_chain(root=pdb_inp.construct_hierarchy_v2())
  assert chain.residue_groups_size() == 4
  assert list(
    chain.split_residue_groups_with_mixed_resnames_but_only_blank_altloc()) \
      == [(1,3), (5,3)]
  assert list(
    chain.split_residue_groups_with_mixed_resnames_but_only_blank_altloc()) \
      == []
  for residue_group in chain.residue_groups():
    assert residue_group.atom_groups_size() == 1
    assert residue_group.atom_groups()[0].parent().memory_id() \
        == residue_group.memory_id()

def exercise(args):
  print "iotbx.pdb.hierarchy_v2.atom.sizeof_data():", \
    pdb.hierarchy_v2.atom.sizeof_data()
  forever = "--forever" in args
  while True:
    exercise_atom()
    exercise_atom_group()
    exercise_residue_group()
    exercise_chain()
    exercise_model()
    exercise_root()
    exercise_format_atom_record_using_parents()
    exercise_construct_hierarchy()
    exercise_merge_atom_groups()
    exercise_merge_residue_groups()
    exercise_chain_merge_and_split_residue_groups()
    if (not forever): break
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
