from __future__ import absolute_import, division, print_function
from cctbx import geometry_restraints
from cctbx.array_family import flex
from libtbx import phil
from libtbx.utils import format_cpu_times
from libtbx.test_utils import Exception_expected, show_diff
from six.moves import cStringIO as StringIO

def exercise_atom():
  a = geometry_restraints.motif_atom(name="a")
  assert a.name == "a"
  assert a.scattering_type == ""
  assert a.nonbonded_type == ""
  assert a.partial_charge == 0
  a = geometry_restraints.motif_atom(
    name="b", scattering_type="c", nonbonded_type="d", partial_charge=0.5)
  assert a.name == "b"
  assert a.scattering_type == "c"
  assert a.nonbonded_type == "d"
  assert a.partial_charge == 0.5
  a.name = "x"
  assert a.name == "x"
  a.scattering_type = "y"
  assert a.scattering_type == "y"
  a.nonbonded_type = "z"
  assert a.nonbonded_type == "z"
  a.partial_charge = 1.5
  assert a.partial_charge == 1.5

def exercise_bond():
  b = geometry_restraints.motif_bond(atom_names=["a", "b"])
  assert b.atom_names == ("a", "b")
  assert b.type == ""
  assert b.distance_ideal == 0
  assert b.weight == 0
  assert b.id == ""
  b = geometry_restraints.motif_bond(
    atom_names=["c", "d"], type="single", distance_ideal=1.5, weight=2.5,
    id="t")
  assert b.atom_names == ("c", "d")
  assert b.type == "single"
  assert b.distance_ideal == 1.5
  assert b.weight == 2.5
  assert b.id == "t"
  b.atom_names = ["x", "y"]
  assert b.atom_names == ("x", "y")
  b.type = "double"
  assert b.type == "double"
  b.distance_ideal = 3.5
  assert b.distance_ideal == 3.5
  b.weight = 6.5
  assert b.weight == 6.5
  b.id = "u"
  assert b.id == "u"

def exercise_angle():
  a = geometry_restraints.motif_angle(atom_names=["a", "b", "c"])
  assert a.atom_names == ("a", "b", "c")
  assert a.angle_ideal == 0
  assert a.weight == 0
  assert a.id == ""
  a = geometry_restraints.motif_angle(
    atom_names=["x", "y", "z"], angle_ideal=3.5, weight=4.5, id="t")
  assert a.atom_names == ("x", "y", "z")
  assert a.angle_ideal == 3.5
  assert a.weight == 4.5
  assert a.id == "t"
  a.atom_names = ["p", "q", "r"]
  assert a.atom_names == ("p", "q", "r")
  a.angle_ideal = 8.5
  assert a.angle_ideal == 8.5
  a.weight = 9.5
  assert a.weight == 9.5
  a.id = "u"
  assert a.id == "u"

def exercise_dihedral():
  d = geometry_restraints.motif_dihedral(atom_names=["a","b","c","d"])
  assert d.atom_names == ("a","b","c","d")
  assert d.angle_ideal == 0
  assert d.weight == 0
  assert d.periodicity == 0
  assert d.id == ""
  d = geometry_restraints.motif_dihedral(
    atom_names=["e","f","g","h"], angle_ideal=2.5, weight=5.5, periodicity=3,
    id="t")
  assert d.atom_names == ("e","f","g","h")
  assert d.angle_ideal == 2.5
  assert d.weight == 5.5
  assert d.periodicity == 3
  assert d.id == "t"
  d.atom_names = ["h","k","l","m"]
  assert d.atom_names == ("h","k","l","m")
  d.angle_ideal = 12.5
  assert d.angle_ideal == 12.5
  d.weight = 15.5
  assert d.weight == 15.5
  d.periodicity = 13
  assert d.periodicity == 13
  d.id = "u"
  assert d.id == "u"

def exercise_chirality():
  c = geometry_restraints.motif_chirality(atom_names=["a","b","c","d"])
  assert c.atom_names == ("a","b","c","d")
  assert c.volume_sign == ""
  assert not c.both_signs
  assert c.volume_ideal == 0
  assert c.weight == 0
  assert c.id == ""
  c = geometry_restraints.motif_chirality(
    atom_names=["e","f","g","h"], volume_sign="positive", both_signs=True,
    volume_ideal=2.5, weight=5.5, id="t")
  assert c.atom_names == ("e","f","g","h")
  assert c.volume_sign == "positive"
  assert c.both_signs
  assert c.volume_ideal == 2.5
  assert c.weight == 5.5
  assert c.id == "t"
  c.atom_names = ["h","k","l","m"]
  assert c.atom_names == ("h","k","l","m")
  c.volume_sign = "both"
  assert c.volume_sign == "both"
  c.both_signs = False
  assert not c.both_signs
  c.volume_ideal = 12.5
  assert c.volume_ideal == 12.5
  c.weight = 15.5
  assert c.weight == 15.5
  c.id = "u"
  assert c.id == "u"

def exercise_planarity():
  p = geometry_restraints.motif_planarity()
  assert p.atom_names.size() == 0
  assert p.weights.size() == 0
  assert p.id == ""
  p = geometry_restraints.motif_planarity(
    atom_names=flex.std_string(["a", "b"]),
    weights=flex.double([1.5, 2.5]),
    id="t")
  assert list(p.atom_names) == ["a", "b"]
  assert list(p.weights) == [1.5, 2.5]
  assert p.id == "t"
  p.atom_names = flex.std_string(["x", "y", "z"])
  assert list(p.atom_names) == ["x", "y", "z"]
  p.weights = flex.double([3.5, 6.5])
  assert list(p.weights) == [3.5, 6.5]
  p.id = "u"
  assert p.id == "u"

def exercise_motif():
  m = geometry_restraints.motif()
  assert m.id == ""
  assert m.description == ""
  assert m.info.size() == 0
  assert m.manipulation_ids.size() == 0
  assert len(m.atoms_as_list()) == 0
  assert len(m.bonds_as_list()) == 0
  assert len(m.angles_as_list()) == 0
  assert len(m.dihedrals_as_list()) == 0
  assert len(m.chiralities_as_list()) == 0
  assert len(m.planarities_as_list()) == 0
  m.id = "a"
  assert m.id == "a"
  m.description = "b"
  assert m.description == "b"
  m.info.append("c")
  m.info.append("d")
  assert m.info.size() == 2
  m.manipulation_ids.append("e")
  m.manipulation_ids.append("f")
  assert m.manipulation_ids.size() == 2
  l = [geometry_restraints.motif_atom(*a, **k) for a,k in [
    ((), {"name": "a", "scattering_type": "b", "nonbonded_type": "c",
     "partial_charge": -0.1}),
    (("d"),{})]]
  m.set_atoms(l)
  assert [atom.name for atom in m.atoms_as_list()] == ["a", "d"]
  l = [geometry_restraints.motif_bond(*a, **k) for a,k in [
    ((), {"atom_names": ("a", "b"), "type": "c", "distance_ideal": 2.3,
     "weight": 11.5, "id": "d"}),
    ((("e", "f"),),{})]]
  m.set_bonds(l)
  assert [bond.atom_names for bond in m.bonds_as_list()] \
      == [("a", "b"), ("e", "f")]
  l = [geometry_restraints.motif_angle(*a, **k) for a,k in [
    ((), {"atom_names": ("a", "b", "c"), "angle_ideal": 3.2,
     "weight": 15.1, "id": "d"}),
    ((("e", "f", "g"),),{})]]
  m.set_angles(l)
  assert [angle.atom_names for angle in m.angles_as_list()] \
      == [("a", "b", "c"), ("e", "f", "g")]
  l = [geometry_restraints.motif_dihedral(*a, **k) for a,k in [
    ((), {"atom_names": ("a", "b", "c", "d"), "angle_ideal": 4.3,
     "weight": 1.5, "periodicity": 3, "id": "e"}),
    ((("f", "g", "h", "i"),),{})]]
  m.set_dihedrals(l)
  assert [dihedral.atom_names for dihedral in m.dihedrals_as_list()] \
      == [("a", "b", "c", "d"), ("f", "g", "h", "i")]
  l = [geometry_restraints.motif_chirality(*a, **k) for a,k in [
    ((), {"atom_names": ("a", "b", "c", "d"), "volume_sign": "e",
     "both_signs": True, "volume_ideal": 10.3, "weight": 1.5, "id": "f"}),
    ((("g", "h", "i", "j"),),{})]]
  m.set_chiralities(l)
  assert [chirality.atom_names for chirality in m.chiralities_as_list()] \
      == [("a", "b", "c", "d"), ("g", "h", "i", "j")]
  l = [geometry_restraints.motif_planarity(
         atom_names=flex.std_string(["a", "b"]),
         weights=flex.double([1,2]),
         id="c"),
       geometry_restraints.motif_planarity(
         atom_names=flex.std_string(["d"]),
         weights=flex.double([3]))]
  m.set_planarities(l)
  assert [list(planarity.atom_names)
            for planarity in m.planarities_as_list()] \
      == [["a", "b"], ["d"]]
  out = StringIO()
  m.show(out=out, prefix="&*")
  assert not show_diff(out.getvalue(), """\
&*geometry_restraints.motif {
&*  id = "a"
&*  description = "b"
&*  info = "c"
&*  info = "d"
&*  manipulation_id = "e"
&*  manipulation_id = "f"
&*  atom = [name scattering_type nonbonded_type partial_charge]
&*  atom = "a" "b" "c" -0.1
&*  atom = "d" "" "" 0
&*  bond = [atom_name*2 type distance_ideal weight id]
&*  bond = "a" "b" "c" 2.3 11.5 "d"
&*  bond = "e" "f" "" 0 0 ""
&*  angle = [atom_name*3 angle_ideal weight id]
&*  angle = "a" "b" "c" 3.2 15.1 "d"
&*  angle = "e" "f" "g" 0 0 ""
&*  dihedral = [atom_name*4 angle_ideal weight periodicity id]
&*  dihedral = "a" "b" "c" "d" 4.3 1.5 3 "e"
&*  dihedral = "f" "g" "h" "i" 0 0 0 ""
&*  chirality = [atom_name*4 volume_sign both_signs volume_ideal weight id]
&*  chirality = "a" "b" "c" "d" "e" True 10.3 1.5 "f"
&*  chirality = "g" "h" "i" "j" "" False 0 0 ""
&*  planarity {
&*    id = "c"
&*    atom = [name weight]
&*    atom = "a" 1
&*    atom = "b" 2
&*  }
&*  planarity {
&*    id = ""
&*    atom = [name weight]
&*    atom = "d" 3
&*  }
&*}
""")
  out = StringIO()
  m.show(out=out)
  phil.parse(input_string=out.getvalue())

def exercise_alteration():
  a = geometry_restraints.motif_alteration()
  assert a.action == ""
  assert a.operand == ""
  a = geometry_restraints.motif_alteration(action="add")
  assert a.action == "add"
  assert a.operand == ""
  a = geometry_restraints.motif_alteration(action="change", operand="bond")
  assert a.action == "change"
  assert a.operand == "bond"
  assert a.motif_ids.size() == 0
  assert a.atom.name == ""
  assert a.motif_atom_name == ""
  assert a.bond.atom_names == ("", "")
  assert a.angle.atom_names == ("", "", "")
  assert a.dihedral.atom_names == ("", "", "", "")
  assert a.chirality.atom_names == ("", "", "", "")
  assert a.planarity.atom_names.size() == 0
  assert a.planarity_motif_id == ""
  assert a.planarity_atom_actions_as_list() == []
  assert not a.change_partial_charge()
  assert not a.change_distance_ideal()
  assert not a.change_weight()
  assert not a.change_angle_ideal()
  assert not a.change_periodicity()
  assert not a.change_both_signs()
  assert not a.change_volume_ideal()
  for action in ["", "add", "delete", "change"]:
    a.action = action
    assert a.action == action
  try: a.action = "chnage"
  except RuntimeError as e:
    assert str(e) == 'Unknown cctbx::geometry_restraints::motif::alteration' \
      '::action_type: "chnage"\n' \
      '  Possible action types are: "", "add", "delete", "change"'
  else: raise Exception_expected
  for operand in ["", "atom", "bond", "angle", "dihedral",
                  "chirality", "planarity"]:
    a.operand = operand
    assert a.operand == operand
  try: a.operand = "diehdral"
  except RuntimeError as e:
    assert str(e) == 'Unknown cctbx::geometry_restraints::motif::alteration' \
      '::operand_type: "diehdral"\n' \
      '  Possible operand types are: "", "atom", "bond", "angle",' \
      ' "dihedral", "chirality", "planarity"'
  else: raise Exception_expected
  a.motif_ids = flex.std_string(["-", "+"])
  assert list(a.motif_ids) == ["-", "+"]
  a.atom = geometry_restraints.motif_atom(name="a")
  assert a.atom.name == "a"
  a.motif_atom_name = "x"
  assert a.motif_atom_name == "x"
  a.bond = geometry_restraints.motif_bond(atom_names=["b", "c"])
  assert a.bond.atom_names == ("b", "c")
  a.angle = geometry_restraints.motif_angle(atom_names=["d", "e", "f"])
  assert a.angle.atom_names == ("d", "e", "f")
  a.dihedral = geometry_restraints.motif_dihedral(atom_names=["h","i","j","k"])
  assert a.dihedral.atom_names == ("h", "i", "j", "k")
  a.chirality = geometry_restraints.motif_chirality(
    atom_names=["l","m","n","o"])
  assert a.chirality.atom_names == ("l", "m", "n", "o")
  a.planarity = geometry_restraints.motif_planarity(
    atom_names=flex.std_string(["p", "q"]),
    weights=flex.double([1,2]))
  a.planarity_motif_id = "h"
  assert a.planarity_motif_id == "h"
  a.set_planarity_atom_actions(["add", "change", "delete"])
  assert a.planarity_atom_actions_as_list() == ["add", "change", "delete"]
  try: a.set_planarity_atom_actions(["add", "chnage", "delete"])
  except RuntimeError as e:
    assert str(e) == 'Unknown cctbx::geometry_restraints::motif::alteration' \
      '::action_type: "chnage"\n' \
      '  Possible action types are: "", "add", "delete", "change"'
  else: raise Exception_expected
  assert a.planarity_atom_actions_as_list() == []
  assert list(a.planarity.atom_names) == ["p", "q"]
  assert a.set_change_partial_charge(state=True) is a
  assert a.change_partial_charge()
  assert not a.change_distance_ideal()
  assert a.set_change_distance_ideal(state=True) is a
  assert a.change_distance_ideal()
  assert not a.change_weight()
  assert a.set_change_weight(state=True) is a
  assert a.change_weight()
  assert not a.change_angle_ideal()
  assert a.set_change_angle_ideal(state=True) is a
  assert a.change_angle_ideal()
  assert not a.change_periodicity()
  assert a.set_change_periodicity(state=True) is a
  assert a.change_periodicity()
  assert not a.change_both_signs()
  assert a.set_change_both_signs(state=True) is a
  assert a.change_both_signs()
  assert not a.change_volume_ideal()
  assert a.set_change_volume_ideal(state=True) is a
  assert a.change_volume_ideal()

def exercise_manipulation():
  m = geometry_restraints.motif_manipulation()
  assert m.id == ""
  assert m.description == ""
  assert m.info.size() == 0
  assert len(m.alterations_as_list()) == 0
  m.id = "a"
  assert m.id == "a"
  m.description = "b"
  assert m.description == "b"
  m.info.append("c")
  m.info.append("d")
  assert list(m.info) == ["c", "d"]
  alts = []
  a = geometry_restraints.motif_alteration()
  a.action = "add"
  a.operand = "bond"
  alts.append(a)
  a = geometry_restraints.motif_alteration()
  a.action = "delete"
  a.operand = "angle"
  alts.append(a)
  m.set_alterations(alts)
  l = m.alterations_as_list()
  assert len(l) == 2
  assert l[0].action == "add"
  assert l[0].operand == "bond"
  assert l[1].action == "delete"
  assert l[1].operand == "angle"
  #
  m = geometry_restraints.motif_manipulation()
  m.id = "a"
  m.description = "b"
  m.info.append("c")
  #
  alterations = []
  a = geometry_restraints.motif_alteration(action="add", operand="atom")
  a.motif_ids.append("1")
  a.atom.name = "x"
  a.atom.scattering_type = "y"
  a.atom.nonbonded_type = "z"
  a.atom.partial_charge = -0.2
  alterations.append(a)
  for change_partial_charge in [False, True]:
    a = geometry_restraints.motif_alteration(action="change", operand="atom")
    a.motif_ids.append("0")
    a.motif_atom_name = "p"
    a.atom.name = "q"
    a.atom.scattering_type = "r"
    a.atom.nonbonded_type = "s"
    a.atom.partial_charge = 0.3
    a.set_change_partial_charge(state=change_partial_charge)
    alterations.append(a)
  a = geometry_restraints.motif_alteration(action="delete", operand="atom")
  a.motif_ids.append("3")
  a.motif_atom_name = "d"
  alterations.append(a)
  #
  a = geometry_restraints.motif_alteration(action="add", operand="bond")
  a.motif_ids.append("1")
  a.motif_ids.append("2")
  a.bond.atom_names = ["j", "k"]
  a.bond.distance_ideal = 1.2
  a.bond.weight = 3.1
  a.bond.id = "s"
  alterations.append(a)
  for change_distance_ideal in [False, True]:
    for change_weight in [False, True]:
      a = geometry_restraints.motif_alteration(action="change", operand="bond")
      a.motif_ids.append("3")
      a.motif_ids.append("4")
      a.bond.atom_names = ["l", "m"]
      a.bond.distance_ideal = 4.5
      a.bond.weight = 2.8
      a.bond.id = "g"
      a.set_change_distance_ideal(state=change_distance_ideal)
      a.set_change_weight(state=change_weight)
      alterations.append(a)
  a = geometry_restraints.motif_alteration(action="delete", operand="bond")
  a.motif_ids.append("5")
  a.motif_ids.append("6")
  a.bond.atom_names = ["t", "n"]
  alterations.append(a)
  #
  a = geometry_restraints.motif_alteration(action="add", operand="angle")
  a.motif_ids.append("1")
  a.motif_ids.append("2")
  a.motif_ids.append("3")
  a.angle.atom_names = ["j", "k", "l"]
  a.angle.angle_ideal = 1.2
  a.angle.weight = 3.1
  a.angle.id = "s"
  alterations.append(a)
  for change_angle_ideal in [False, True]:
    for change_weight in [False, True]:
      a = geometry_restraints.motif_alteration(
        action="change", operand="angle")
      a.motif_ids.append("3")
      a.motif_ids.append("4")
      a.motif_ids.append("5")
      a.angle.atom_names = ["l", "m", "g"]
      a.angle.angle_ideal = 4.5
      a.angle.weight = 2.8
      a.angle.id = "g"
      a.set_change_angle_ideal(state=change_angle_ideal)
      a.set_change_weight(state=change_weight)
      alterations.append(a)
  a = geometry_restraints.motif_alteration(action="delete", operand="angle")
  a.motif_ids.append("6")
  a.motif_ids.append("7")
  a.motif_ids.append("8")
  a.angle.atom_names = ["t", "n", "d"]
  alterations.append(a)
  #
  a = geometry_restraints.motif_alteration(action="add", operand="dihedral")
  a.motif_ids.append("1")
  a.motif_ids.append("2")
  a.motif_ids.append("3")
  a.motif_ids.append("4")
  a.dihedral.atom_names = ["j", "k", "l", "r"]
  a.dihedral.angle_ideal = 1.2
  a.dihedral.weight = 3.1
  a.dihedral.periodicity = 20
  a.dihedral.id = "s"
  alterations.append(a)
  for change_angle_ideal in [False, True]:
    for change_weight in [False, True]:
      for change_periodicity in [False, True]:
        a = geometry_restraints.motif_alteration(
          action="change", operand="dihedral")
        a.motif_ids.append("5")
        a.motif_ids.append("6")
        a.motif_ids.append("7")
        a.motif_ids.append("8")
        a.dihedral.atom_names = ["l", "m", "g", "e"]
        a.dihedral.angle_ideal = 4.5
        a.dihedral.weight = 2.8
        a.dihedral.periodicity = 13
        a.dihedral.id = "g"
        a.set_change_angle_ideal(state=change_angle_ideal)
        a.set_change_weight(state=change_weight)
        a.set_change_periodicity(state=change_periodicity)
        alterations.append(a)
  a = geometry_restraints.motif_alteration(action="delete", operand="dihedral")
  a.motif_ids.append("9")
  a.motif_ids.append("10")
  a.motif_ids.append("11")
  a.motif_ids.append("12")
  a.dihedral.atom_names = ["t", "n", "d", "y"]
  alterations.append(a)
  #
  a = geometry_restraints.motif_alteration(action="add", operand="chirality")
  a.motif_ids.append("1")
  a.motif_ids.append("2")
  a.motif_ids.append("3")
  a.motif_ids.append("4")
  a.chirality.atom_names = ["j", "k", "l", "r"]
  a.chirality.volume_sign = "x"
  a.chirality.volume_ideal = 1.2
  a.chirality.weight = 3.1
  a.chirality.id = "w"
  alterations.append(a)
  for change_volume_ideal in [False, True]:
    for change_weight in [False, True]:
      a = geometry_restraints.motif_alteration(
        action="change", operand="chirality")
      a.motif_ids.append("5")
      a.motif_ids.append("6")
      a.motif_ids.append("7")
      a.motif_ids.append("8")
      a.chirality.atom_names = ["l", "m", "g", "e"]
      a.chirality.volume_sign = "u"
      a.chirality.volume_ideal = 4.5
      a.chirality.weight = 2.8
      a.chirality.id = "g"
      a.set_change_volume_ideal(state=change_volume_ideal)
      a.set_change_weight(state=change_weight)
      alterations.append(a)
  a = geometry_restraints.motif_alteration(
    action="delete", operand="chirality")
  a.motif_ids.append("9")
  a.motif_ids.append("10")
  a.motif_ids.append("11")
  a.motif_ids.append("12")
  a.chirality.atom_names = ["t", "n", "d", "y"]
  alterations.append(a)
  #
  a = geometry_restraints.motif_alteration(action="add", operand="planarity")
  a.planarity_motif_id = "7"
  a.planarity.id = "e"
  a.motif_ids.append("1")
  a.motif_ids.append("2")
  a.motif_ids.append("3")
  a.planarity.atom_names.append("s")
  a.planarity.atom_names.append("d")
  a.planarity.atom_names.append("f")
  a.planarity.weights.append(1.2)
  a.planarity.weights.append(2.3)
  a.planarity.weights.append(3.4)
  alterations.append(a)
  a = geometry_restraints.motif_alteration(
    action="change", operand="planarity")
  a.planarity_motif_id = "8"
  a.planarity.id = "f"
  a.set_planarity_atom_actions(["add", "change", "delete"])
  a.motif_ids.append("1")
  a.motif_ids.append("2")
  a.motif_ids.append("3")
  a.planarity.atom_names.append("s")
  a.planarity.atom_names.append("d")
  a.planarity.atom_names.append("f")
  a.planarity.weights.append(1.2)
  a.planarity.weights.append(2.3)
  a.planarity.weights.append(0)
  alterations.append(a)
  a = geometry_restraints.motif_alteration(
    action="delete", operand="planarity")
  a.planarity_motif_id = "4"
  a.planarity.id = "n"
  alterations.append(a)
  #
  m.set_alterations(alterations)
  out = StringIO()
  m.show(out=out, prefix="^")
  assert not show_diff(out.getvalue(), """\
^geometry_restraints.motif_manipulation {
^  id = "a"
^  description = "b"
^  info = "c"
^  atom = add [motif_id name scattering_type nonbonded_type partial_charge]
^  atom = add "1" "x" "y" "z" -0.2
^  atom = change [motif_id motif_atom_name \\
^                 name scattering_type nonbonded_type partial_charge]
^  atom = change "0" "p" \\
^                "q" "r" "s" None
^  atom = change "0" "p" \\
^                "q" "r" "s" 0.3
^  atom = delete [motif_id motif_atom_name]
^  atom = delete "3" "d"
^  bond = add [(motif_id atom_name)*2 type distance_ideal weight id]
^  bond = add "1" "j" "2" "k" 1.2 3.1 "s"
^  bond = change [(motif_id atom_name)*2 type distance_ideal weight id]
^  bond = change "3" "l" "4" "m" None None "g"
^  bond = change "3" "l" "4" "m" None 2.8 "g"
^  bond = change "3" "l" "4" "m" 4.5 None "g"
^  bond = change "3" "l" "4" "m" 4.5 2.8 "g"
^  bond = delete [(motif_id atom_name)*2]
^  bond = delete "5" "t" "6" "n"
^  angle = add [(motif_id atom_name)*3 type angle_ideal weight id]
^  angle = add "1" "j" "2" "k" "3" "l" 1.2 3.1 "s"
^  angle = change [(motif_id atom_name)*3 type angle_ideal weight id]
^  angle = change "3" "l" "4" "m" "5" "g" None None "g"
^  angle = change "3" "l" "4" "m" "5" "g" None 2.8 "g"
^  angle = change "3" "l" "4" "m" "5" "g" 4.5 None "g"
^  angle = change "3" "l" "4" "m" "5" "g" 4.5 2.8 "g"
^  angle = delete [(motif_id atom_name)*3]
^  angle = delete "6" "t" "7" "n" "8" "d"
^  dihedral = add [(motif_id atom_name)*4 angle_ideal weight periodicity id]
^  dihedral = add "1" "j" "2" "k" "3" "l" "4" "r" 1.2 3.1 20 "s"
^  dihedral = change [(motif_id atom_name)*4 angle_ideal weight periodicity id]
^  dihedral = change "5" "l" "6" "m" "7" "g" "8" "e" None None None "g"
^  dihedral = change "5" "l" "6" "m" "7" "g" "8" "e" None None 13 "g"
^  dihedral = change "5" "l" "6" "m" "7" "g" "8" "e" None 2.8 None "g"
^  dihedral = change "5" "l" "6" "m" "7" "g" "8" "e" None 2.8 13 "g"
^  dihedral = change "5" "l" "6" "m" "7" "g" "8" "e" 4.5 None None "g"
^  dihedral = change "5" "l" "6" "m" "7" "g" "8" "e" 4.5 None 13 "g"
^  dihedral = change "5" "l" "6" "m" "7" "g" "8" "e" 4.5 2.8 None "g"
^  dihedral = change "5" "l" "6" "m" "7" "g" "8" "e" 4.5 2.8 13 "g"
^  dihedral = delete [(motif_id atom_name)*4]
^  dihedral = delete "9" "t" "10" "n" "11" "d" "12" "y"
^  chirality = add [(motif_id atom_name)*4 \\
^                   volume_sign volume_ideal weight id]
^  chirality = add "1" "j" "2" "k" "3" "l" "4" "r" \\
^                  "x" 1.2 3.1 "w"
^  chirality = change [(motif_id atom_name)*4 \\
^                      volume_sign volume_ideal weight id]
^  chirality = change "5" "l" "6" "m" "7" "g" "8" "e" \\
^                     "u" None None "g"
^  chirality = change "5" "l" "6" "m" "7" "g" "8" "e" \\
^                     "u" None 2.8 "g"
^  chirality = change "5" "l" "6" "m" "7" "g" "8" "e" \\
^                     "u" 4.5 None "g"
^  chirality = change "5" "l" "6" "m" "7" "g" "8" "e" \\
^                     "u" 4.5 2.8 "g"
^  chirality = delete [(motif_id atom_name)*4]
^  chirality = delete "9" "t" "10" "n" "11" "d" "12" "y"
^  planarity {
^    action = add
^    motif_id = "7"
^    id = "e"
^    atom = [motif_id name weight]
^    atom = "1" "s" 1.2
^    atom = "2" "d" 2.3
^    atom = "3" "f" 3.4
^  }
^  planarity {
^    action = change
^    motif_id = "8"
^    id = "f"
^    atom = add [motif_id name weight]
^    atom = add "1" "s" 1.2
^    atom = change [motif_id name weight]
^    atom = change "2" "d" 2.3
^    atom = delete [motif_id name]
^    atom = delete "3" "f"
^  }
^  planarity {
^    action = delete
^    motif_id = "4"
^    id = "n"
^  }
^}
""")
  out = StringIO()
  m.show(out=out)
  phil.parse(input_string=out.getvalue())

def exercise():
  exercise_atom()
  exercise_bond()
  exercise_angle()
  exercise_dihedral()
  exercise_chirality()
  exercise_planarity()
  exercise_motif()
  exercise_alteration()
  exercise_manipulation()
  print(format_cpu_times())

if (__name__ == "__main__"):
  exercise()
