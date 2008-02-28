from iotbx import pdb
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, approx_equal
from libtbx.utils import format_cpu_times
import sys

def exercise_atom():
  a = pdb.hierarchy_v2.atom()
  assert a.name == ""
  a.name = "abcd"
  assert a.name == "abcd"
  try: a.name = "xyzhkl"
  except ValueError, e:
    assert str(e) == "string is too long for name attribute " \
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
  #
  a.tmp = 7
  ac = a.detached_copy()
  assert ac.tmp == 0
  assert ac.name == "NaMe"
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
  assert ch_a.parent().memory_id() == m.memory_id()
  assert m.chains_size() == 1
  assert len(m.chains()) == 1
  ch_b = pdb.hierarchy_v2.chain(id="b")
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

def exercise(args):
  assert len(args) == 0
  exercise_atom()
  exercise_atom_group()
  exercise_residue_group()
  exercise_chain()
  exercise_model()
  exercise_root()
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
