from iotbx.pdb import rna_dna_detection
from cctbx import geometry_restraints
from cctbx import uctbx
from cctbx.array_family import flex
from libtbx.utils import if_none
from libtbx import slots_getstate_setstate
import copy
import sys

peptide_comp_groups = ("L-peptide", "D-peptide")
dna_rna_comp_groups = ("DNA", "RNA")
non_chain_links = ("SS-bridge",)

class looped_data(object):
  "mix-in to extract keywords from self.__class__.__doc__"

  def __init__(self, **keyword_arguments):
    self._looped_ids = self.__class__.__doc__.split()
    self._cif_keywords = []
    for looped_id in self.looped_ids():
      key = looped_id.split(".")[1]
      typ = None
      key_typ = key.split(":")
      if (len(key_typ) == 2):
        key = key_typ[0]
        typ = key_typ[1]
        if (typ not in ("int", "float")):
          raise RuntimeError("Corrupt looped id: %s" % looped_id)
      val = keyword_arguments.get(key, keyword_arguments.get(
        looped_id.split(":")[0], None))
      if (val == "."): val = ""
      if (typ is not None):
        if (val == ""): val = None
        if (val is not None):
          def raise_value_error(expected_type):
            raise ValueError(
              "%s value expected in CIF file: %s %s" % (
                expected_type, looped_id.split(":",1)[0], val))
          if (typ == "int"):
            try:
              val = int(val)
            except ValueError:
              raise_value_error("Integer")
          else:
            try:
              val = float(val)
            except ValueError:
              raise_value_error("Floating-point")
      setattr(self, key, val)
      self._cif_keywords.append(key)

  def looped_ids(self):
    return self._looped_ids

  def cif_keywords(self):
    return self._cif_keywords

  def __copy__(self):
    args = {}
    for key in self.cif_keywords():
      args[key] = getattr(self, key)
    return self.__class__(**args)

  def show(self, f=None):
    if (f is None): f = sys.stdout
    for key in self.cif_keywords():
      print >> f, "_%s.%s: %s" % (
        self.__class__.__name__, key, getattr(self, key))

  def show_loop_header(self, f):
    print >> f, "loop_"
    for key in self.cif_keywords():
      print >> f, "_%s.%s" % (
        self.__class__.__name__, key)

  def show_loop_data(self, f):
    for key in self.cif_keywords():
      val = getattr(self, key)
      if (val in (None, "")): val = "."
      print >> f, val,
    print >> f

def show_loop(data_list, f):
  first = True
  for datum in data_list:
    if (first):
      datum.show_loop_header(f=f)
      first = False
    datum.show_loop_data(f=f)
  if (not first):
    print >> f

def apply_chem_mod_list(label, mod_list, result_list):
  for mod in mod_list:
    if (mod.function == "add"):
      result_list.append(mod.as_chem_comp())
    elif (mod.function == "delete"):
      new_list = []
      for orig in result_list:
        if (not mod.is_matching_mod_for(orig)):
          new_list.append(orig)
      result_list = new_list
    elif (mod.function == "change"):
      for orig in result_list:
        if (mod.is_matching_mod_for(orig)):
          mod.apply_change_in_place(orig)
    else:
      raise RuntimeError("Unknown _chem_mod_%s.function: %s"
        % (label, str(mod.function)))
  return result_list

class chem_plane(object):

  def __init__(self, plane_id):
    self.plane_id = plane_id
    self.plane_atoms = []

def get_bond(bond_list, atom_id_1, atom_id_2):
  for bond in bond_list:
    if (bond.atom_id_1 == atom_id_1 and bond.atom_id_2 == atom_id_2):
      return bond
    if (bond.atom_id_1 == atom_id_2 and bond.atom_id_2 == atom_id_1):
      return bond
  return None

def normalized_bond_list(bond_list):
  result = []
  for bond in bond_list:
    bond = copy.copy(bond)
    if (bond.atom_id_1 > bond.atom_id_2):
      bond.atom_id_1, bond.atom_id_2 = bond.atom_id_2, bond.atom_id_1
    result.append(bond)
  return result

def get_chir_volume_ideal(volume_sign, bonds, angles):
  if (None in bonds or None in angles):
    return None
  try:
    uc = uctbx.unit_cell(  [bond.value_dist for bond in bonds]
                         + [angle.value_angle for angle in angles])
  except RuntimeError:
    return None
  if (volume_sign.startswith("neg")):
    return -uc.volume()
  return uc.volume()

def group_planes(plane_atom_list):
  result = []
  result_dict = {}
  for plane_atom in plane_atom_list:
    plane = result_dict.get(plane_atom.plane_id, None)
    if (plane is None):
      plane = chem_plane(plane_atom.plane_id)
      result.append(plane)
      result_dict[plane_atom.plane_id] = plane
    plane.plane_atoms.append(plane_atom)
  return result

def esd_as_weight(esd):
  if (esd is None): return 0
  if (esd == 0): return 0
  return 1./(esd*esd)

class comp_comp_id(slots_getstate_setstate):

  __slots__ = [
    "source_info",
    "chem_comp",
    "atom_list",
    "bond_list",
    "angle_list",
    "tor_list",
    "chir_list",
    "plane_atom_list",
    "rotamer_info_phil_str_list",
    "__rotamer_info",
    "classification"]

  def __init__(self, source_info, chem_comp):
    self.source_info = source_info
    self.chem_comp = chem_comp
    self.atom_list = []
    self.bond_list = []
    self.angle_list = []
    self.tor_list = []
    self.chir_list = []
    self.plane_atom_list = []
    self.rotamer_info_phil_str_list = []
    self.__rotamer_info = None
    self.classification = None

  def normalize_atom_ids_in_place(self):
    atom_ids_mod = []
    atom_ids_mod_set = set()
    is_rna_dna = (self.get_classification() in ["RNA", "RNAv3", "DNA", "DNAv3"])
    for atom in self.atom_list:
      atom_id = atom.atom_id.replace("'", "*")
      if is_rna_dna:
        if atom_id == "OP1":
          atom_id = "O1P"
        elif atom_id == "OP2":
          atom_id = "O2P"
      if (atom_id in atom_ids_mod_set):
        return False # changing ids results in ambiguity
      atom_ids_mod.append(atom_id)
      atom_ids_mod_set.add(atom_id)
    del atom_ids_mod_set
    #
    atom_id_map = {}
    for atom,atom_id in zip(self.atom_list, atom_ids_mod):
      if (atom.atom_id != atom_id):
        atom_id_map[atom.atom_id] = atom_id
        atom.atom_id = atom_id
    if (len(atom_id_map) == 0):
      return False
    def replace(obj_list, attrs):
      for obj in obj_list:
        for attr in attrs:
          atom_id = atom_id_map.get(getattr(obj, attr))
          if (atom_id is not None):
            setattr(obj, attr, atom_id)
    replace(self.bond_list, (
      "atom_id_1", "atom_id_2"))
    replace(self.angle_list, (
      "atom_id_1", "atom_id_2", "atom_id_3"))
    replace(self.tor_list, (
      "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4"))
    replace(self.chir_list, (
      "atom_id_centre", "atom_id_1", "atom_id_2", "atom_id_3"))
    replace(self.plane_atom_list, (
      "atom_id",))
    return True

  def __copy__(self):
    result = comp_comp_id(source_info=None, chem_comp=self.chem_comp)
    result.atom_list = [copy.copy(e) for e in self.atom_list]
    result.bond_list = [copy.copy(e) for e in self.bond_list]
    result.angle_list = [copy.copy(e) for e in self.angle_list]
    result.tor_list = [copy.copy(e) for e in self.tor_list]
    result.chir_list = [copy.copy(e) for e in self.chir_list]
    result.plane_atom_list = [copy.copy(e) for e in self.plane_atom_list]
    result.rotamer_info_phil_str_list = [
      copy.copy(e) for e in self.rotamer_info_phil_str_list]
    result.__rotamer_info = copy.deepcopy(self.__rotamer_info)
    result.classification = None
    return result

  def atom_ids(self, sorted=False):
    result = []
    for atom in self.atom_list:
      result.append(atom.atom_id)
    if (sorted): result.sort()
    return result

  def atom_dict(self):
    result = {}
    for atom in self.atom_list:
      result[atom.atom_id] = atom
    return result

  def hydrogen_deuterium_aliases(self):
    result = {}
    for atom in self.atom_list:
      if (atom.type_symbol != "H"): continue
      atom_id = atom.atom_id
      if (atom_id[:1] != "H"): continue
      result["D"+atom_id[1:]] = atom_id
    return result

  def atom_by_id(self, atom_id):
    for atom in self.atom_list:
      if (atom.atom_id == atom_id):
        return atom
    return None

  def i_atom_by_id(self, atom_id):
    for i_atom in xrange(len(self.atom_list)):
      if (self.atom_list[i_atom].atom_id == atom_id):
        return i_atom
    return None

  def non_hydrogen_atoms(self):
    result = []
    for atom in self.atom_list:
      if (atom.type_symbol != "H"):
        result.append(atom)
    return result

  def get_bond(self, atom_id_1, atom_id_2):
    return get_bond(self.bond_list, atom_id_1, atom_id_2)

  def normalized_bond_list(self):
    return normalized_bond_list(bond_list=self.bond_list)

  def get_angle(self, atom_id_1, atom_id_2, atom_id_3):
    for angle in self.angle_list:
      if (angle.atom_id_2 == atom_id_2):
        if (angle.atom_id_1 == atom_id_1 and angle.atom_id_3 == atom_id_3):
          return angle
        if (angle.atom_id_1 == atom_id_3 and angle.atom_id_3 == atom_id_1):
          return angle
    return None

  def normalized_angle_list(self):
    result = []
    for angle in self.angle_list:
      angle = copy.copy(angle)
      if (angle.atom_id_1 > angle.atom_id_3):
        angle.atom_id_1, angle.atom_id_3 = angle.atom_id_3, angle.atom_id_1
      result.append(angle)
    return result

  def delete_atom_in_place(self, atom_id):
    i_atom = self.i_atom_by_id(atom_id)
    if (i_atom is None):
      raise RuntimeError("delete_atom_in_place: unknown atom_id: %s" % atom_id)
    del self.atom_list[i_atom]
    new_bond_list = []
    for bond in self.bond_list:
      if (atom_id not in (bond.atom_id_1, bond.atom_id_2)):
        new_bond_list.append(bond)
    self.bond_list = new_bond_list
    new_angle_list = []
    for angle in self.angle_list:
      if (atom_id not in (angle.atom_id_1, angle.atom_id_2, angle.atom_id_3)):
        new_angle_list.append(angle)
    self.angle_list = new_angle_list
    new_tor_list = []
    for tor in self.tor_list:
      if (atom_id not in (tor.atom_id_1, tor.atom_id_2,
                          tor.atom_id_3, tor.atom_id_4)):
        new_tor_list.append(tor)
    self.tor_list = new_tor_list
    new_chir_list = []
    for chir in self.chir_list:
      if (atom_id not in (chir.atom_id_centre, chir.atom_id_1,
                          chir.atom_id_2, chir.atom_id_3)):
        new_chir_list.append(chir)
    self.chir_list = new_chir_list
    new_plane_atom_list = []
    for plane_atom in self.plane_atom_list:
      if (atom_id != plane_atom.atom_id):
        new_plane_atom_list.append(plane_atom)
    self.plane_atom_list = new_plane_atom_list

  def change_atom_in_place(self, mod_atom):
    i_atom = self.i_atom_by_id(mod_atom.atom_id)
    if (i_atom is None):
      raise RuntimeError(
        "change_atom_in_place: unknown atom_id: %s" % mod_atom.atom_id)
    atom = copy.copy(self.atom_list[i_atom])
    old_atom_id = atom.atom_id
    for attr_name in atom.__dict__.keys():
      new_attr = getattr(mod_atom, "new_"+attr_name, None)
      if (new_attr not in (None, "")):
        setattr(atom, attr_name, new_attr)
    self.atom_list[i_atom] = atom
    if (atom.atom_id != old_atom_id):
      for bond in self.bond_list:
        if (bond.atom_id_1 == old_atom_id):
          bond.atom_id_1 = atom.atom_id
        if (bond.atom_id_2 == old_atom_id):
          bond.atom_id_2 = atom.atom_id
      for angle in self.angle_list:
        if (angle.atom_id_1 == old_atom_id):
          angle.atom_id_1 = atom.atom_id
        if (angle.atom_id_2 == old_atom_id):
          angle.atom_id_2 = atom.atom_id
        if (angle.atom_id_3 == old_atom_id):
          angle.atom_id_3 = atom.atom_id
      for tor in self.tor_list:
        if (tor.atom_id_1 == old_atom_id):
          tor.atom_id_1 = atom.atom_id
        if (tor.atom_id_2 == old_atom_id):
          tor.atom_id_2 = atom.atom_id
        if (tor.atom_id_3 == old_atom_id):
          tor.atom_id_3 = atom.atom_id
        if (tor.atom_id_4 == old_atom_id):
          tor.atom_id_4 = atom.atom_id
      for chir in self.chir_list:
        if (chir.atom_id_centre == old_atom_id):
          chir.atom_id_centre = atom.atom_id
        if (chir.atom_id_1 == old_atom_id):
          chir.atom_id_1 = atom.atom_id
        if (chir.atom_id_2 == old_atom_id):
          chir.atom_id_2 = atom.atom_id
        if (chir.atom_id_3 == old_atom_id):
          chir.atom_id_3 = atom.atom_id
      for plane_atom in self.plane_atom_list:
        if (plane_atom.atom_id == old_atom_id):
          plane_atom.atom_id = atom.atom_id

  def apply_mod(self, mod):
    result = copy.copy(self)
    for mod_atom in mod.atom_list:
      if (mod_atom.function == "add"):
        result.atom_list.append(mod_atom.as_chem_comp())
      elif (mod_atom.function == "delete"):
        result.delete_atom_in_place(mod_atom.atom_id)
      elif (mod_atom.function == "change"):
        result.change_atom_in_place(mod_atom)
      else:
        raise RuntimeError("Unknown _chem_mod_atom.function: "
          + str(mod_atom.function))
    result.bond_list = apply_chem_mod_list(
      "bond", mod.bond_list, result.bond_list)
    result.angle_list = apply_chem_mod_list(
      "angle", mod.angle_list, result.angle_list)
    result.tor_list = apply_chem_mod_list(
      "tor", mod.tor_list, result.tor_list)
    result.chir_list = apply_chem_mod_list(
      "chir", mod.chir_list, result.chir_list)
    result.plane_atom_list = apply_chem_mod_list(
      "plane_atom", mod.plane_atom_list, result.plane_atom_list)
    return result

  def test_for_peptide(self, atom_dict):
    for required_atom in ("N", "CA", "C", "O"):
      if (not required_atom in atom_dict): return None
    backbone_bonds = {}
    for bond in self.bond_list:
      pair = [bond.atom_id_1, bond.atom_id_2]
      pair.sort()
      pair = " ".join(pair)
      if (pair in ("CA N", "C CA")):
        backbone_bonds[pair] = 0
      elif (pair == "C O" and bond.type != "single"):
        backbone_bonds[pair] = 0
    if (len(backbone_bonds) != 3): return None
    return "peptide"

  def test_for_rna_dna(self, atom_dict):
    return rna_dna_detection.classification(
      atom_dict=atom_dict, bond_list=self.bond_list)

  def test_for_water(self, atom_dict):
    atom_list = atom_dict.keys()
    atom_list.sort()
    if (atom_list == ["H1", "H2", "O"]):
      return "water"
    return None

  def set_classification(self):
    atom_dict = self.atom_dict()
    pep = self.test_for_peptide(atom_dict)
    nuc = self.test_for_rna_dna(atom_dict)
    if (pep is not None and nuc is None):
      self.classification = pep
    elif (pep is None and nuc is not None):
      self.classification = nuc
    else:
      self.classification = self.test_for_water(atom_dict)
      if (self.classification is not None):
        return self
      self.classification = "undetermined"
    return self

  def get_classification(self):
    if (self.classification is None): self.set_classification()
    return self.classification

  def is_peptide(self):
    return self.get_classification() == "peptide"

  def is_rna(self):
    return self.get_classification() == "RNA"

  def is_dna(self):
    return self.get_classification() == "DNA"

  def is_rna_dna(self):
    return self.get_classification() in ("RNA", "DNA")

  def is_water(self):
    return self.get_classification() == "water"

  def get_chir_bonds_and_angles(self, chir):
    bonds = []
    for atom_id_i in [chir.atom_id_1, chir.atom_id_2, chir.atom_id_3]:
      bonds.append(self.get_bond(chir.atom_id_centre, atom_id_i))
    angles = []
    for atom_id_i,atom_id_j in [(chir.atom_id_2,chir.atom_id_3),
                                (chir.atom_id_1,chir.atom_id_3),
                                (chir.atom_id_1,chir.atom_id_2)]:
      angles.append(self.get_angle(atom_id_i,chir.atom_id_centre,atom_id_j))
    return bonds, angles

  def get_chir_volume_ideal(self, chir):
    bonds, angles = self.get_chir_bonds_and_angles(chir)
    return get_chir_volume_ideal(chir.volume_sign, bonds, angles)

  def get_planes(self):
    return group_planes(self.plane_atom_list)

  def show(self, f=None):
    if (f is None): f = sys.stdout
    show_loop(data_list=self.atom_list, f=f)
    show_loop(data_list=self.bond_list, f=f)
    show_loop(data_list=self.angle_list, f=f)
    show_loop(data_list=self.tor_list, f=f)
    show_loop(data_list=self.chir_list, f=f)
    show_loop(data_list=self.plane_atom_list, f=f)
    show_loop(data_list=self.rotamer_info_phil_str_list, f=f)
    return self

  def as_geometry_restraints_motif(self):
    result = geometry_restraints.motif()
    result.id = if_none(self.chem_comp.id, "")
    result.description = if_none(self.chem_comp.name, "").strip()
    if (self.source_info is not None):
      result.info.append(self.source_info)
    result.set_atoms([
      geometry_restraints.motif_atom(
        name=if_none(atom.atom_id, ""),
        scattering_type=if_none(atom.type_symbol, ""),
        nonbonded_type=if_none(atom.type_energy, ""),
        partial_charge=if_none(atom.partial_charge, 0))
          for atom in self.atom_list])
    result.set_bonds([
      geometry_restraints.motif_bond(
        atom_names=[
          if_none(bond.atom_id_1, ""),
          if_none(bond.atom_id_2, "")],
        type=if_none(bond.type, ""),
        distance_ideal=if_none(bond.value_dist, 0),
        weight=esd_as_weight(bond.value_dist_esd))
          for bond in self.bond_list])
    result.set_angles([
      geometry_restraints.motif_angle(
        atom_names=[
          if_none(angle.atom_id_1, ""),
          if_none(angle.atom_id_2, ""),
          if_none(angle.atom_id_3, "")],
        angle_ideal=if_none(angle.value_angle, 0),
        weight=esd_as_weight(angle.value_angle_esd))
          for angle in self.angle_list])
    result.set_dihedrals([
      geometry_restraints.motif_dihedral(
        atom_names=[
          if_none(tor.atom_id_1, ""),
          if_none(tor.atom_id_2, ""),
          if_none(tor.atom_id_3, ""),
          if_none(tor.atom_id_4, "")],
        angle_ideal=if_none(tor.value_angle, 0),
        weight=esd_as_weight(tor.value_angle_esd),
        periodicity=if_none(tor.period, 0),
        id=tor.id)
          for tor in self.tor_list])
    result.set_chiralities([
      geometry_restraints.motif_chirality(
        atom_names=[
          if_none(chir.atom_id_centre, ""),
          if_none(chir.atom_id_1, ""),
          if_none(chir.atom_id_2, ""),
          if_none(chir.atom_id_3, "")],
        volume_sign=chir.volume_sign,
        id=chir.id)
          for chir in self.chir_list])
    planarities = []
    for plane in self.get_planes():
      atom_names = flex.std_string([if_none(plane_atom.atom_id, "")
        for plane_atom in plane.plane_atoms])
      weights = flex.double([esd_as_weight(plane_atom.dist_esd)
        for plane_atom in plane.plane_atoms])
      planarities.append(
        geometry_restraints.motif_planarity(
          atom_names=atom_names,
          weights=weights,
          id=plane.plane_id))
    result.set_planarities(planarities)
    return result

  def rotamer_info(self):
    if (self.__rotamer_info is None):
      if (len(self.rotamer_info_phil_str_list) == 0):
        return None
      assert len(self.rotamer_info_phil_str_list) == 1
      from mmtbx.monomer_library.rotamer_utils import rotamer_info_master_phil
      import libtbx.phil
      self.__rotamer_info = rotamer_info_master_phil().fetch(
        source=libtbx.phil.parse(
          input_string=self.rotamer_info_phil_str_list[0].phil_str)).extract()
    return self.__rotamer_info

  def rotamer_iterator(self, atom_names, sites_cart, fine_sampling=False):
    from mmtbx.monomer_library.rotamer_utils import rotamer_iterator
    return rotamer_iterator(
      comp_comp_id=self,
      atom_names=atom_names,
      sites_cart=sites_cart, fine_sampling=fine_sampling)

class chem_comp(looped_data):
  """
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all:int
_chem_comp.number_atoms_nh:int
_chem_comp.desc_level
  """

class chem_comp_atom(looped_data):
  """
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge:float
  """

class chem_comp_tree(looped_data):
  """
_chem_comp_tree.atom_id
_chem_comp_tree.atom_back
_chem_comp_tree.atom_forward
_chem_comp_tree.connect_type
  """

class chem_comp_bond(looped_data):
  """
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist:float
_chem_comp_bond.value_dist_esd:float
  """

  def atom_ids(O):
    return (O.atom_id_1, O.atom_id_2)

class chem_comp_angle(looped_data):
  """
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle:float
_chem_comp_angle.value_angle_esd:float
  """

class chem_comp_tor(looped_data):
  """
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle:float
_chem_comp_tor.alt_value_angle
_chem_comp_tor.value_angle_esd:float
_chem_comp_tor.period:int
  """

  def atom_ids(O):
    return (O.atom_id_1, O.atom_id_2, O.atom_id_3, O.atom_id_4)

class chem_comp_chir(looped_data):
  """
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
  """

class chem_comp_plane_atom(looped_data):
  """
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd:float
  """

class chem_comp_rotamer_info(looped_data):
  """
_chem_comp_rotamer_info.phil_str
  """

class chem_comp_deriv(looped_data):
  """
_chem_comp_deriv.comp_id
_chem_comp_deriv.source_comp_id
_chem_comp_deriv.mod_id
_chem_comp_deriv.name
_chem_comp_deriv.group
  """

class chem_comp_synonym(looped_data):
  """
_chem_comp_synonym.comp_id
_chem_comp_synonym.comp_alternative_id
_chem_comp_synonym.mod_id
  """

class chem_comp_synonym_atom(looped_data):
  """
_chem_comp_synonym_atom.comp_id
_chem_comp_synonym_atom.comp_alternative_id
_chem_comp_synonym_atom.atom_id
_chem_comp_synonym_atom.atom_alternative_id
  """

class link_link_id:

  def __init__(self, source_info, chem_link):
    self.source_info = source_info
    self.chem_link = chem_link
    self.bond_list = []
    self.angle_list = []
    self.tor_list = []
    self.chir_list = []
    self.plane_list = []

  def get_bond(self, m_i, m_j,
                     atom_1_comp_id, atom_id_1,
                     atom_2_comp_id, atom_id_2):
    if (atom_2_comp_id == atom_1_comp_id):
      if (atom_1_comp_id == 1):
        return m_i.get_bond(atom_id_1, atom_id_2)
      if (atom_1_comp_id == 2):
        return m_j.get_bond(atom_id_1, atom_id_2)
      raise AssertionError
    for bond in self.bond_list:
      if (    bond.atom_1_comp_id == atom_1_comp_id
          and bond.atom_2_comp_id == atom_2_comp_id
          and bond.atom_id_1 == atom_id_1
          and bond.atom_id_2 == atom_id_2):
        return bond
      if (    bond.atom_1_comp_id == atom_2_comp_id
          and bond.atom_2_comp_id == atom_1_comp_id
          and bond.atom_id_1 == atom_id_2
          and bond.atom_id_2 == atom_id_1):
        return bond
    return None

  def get_angle(self, m_i, m_j,
                      atom_1_comp_id, atom_id_1,
                      atom_2_comp_id, atom_id_2,
                      atom_3_comp_id, atom_id_3):
    if (    atom_2_comp_id == atom_1_comp_id
        and atom_3_comp_id == atom_1_comp_id):
      if   (atom_1_comp_id == 1):
        return m_i.get_angle(atom_id_1, atom_id_2, atom_id_3)
      if (atom_1_comp_id == 2):
        return m_j.get_angle(atom_id_1, atom_id_2, atom_id_3)
      raise AssertionError
    for angle in self.angle_list:
      if (    angle.atom_2_comp_id == atom_2_comp_id
          and angle.atom_id_2 == atom_id_2):
        if (    angle.atom_1_comp_id == atom_1_comp_id
            and angle.atom_3_comp_id == atom_3_comp_id
            and angle.atom_id_1 == atom_id_1
            and angle.atom_id_3 == atom_id_3):
          return angle
        if (    angle.atom_1_comp_id == atom_3_comp_id
            and angle.atom_3_comp_id == atom_1_comp_id
            and angle.atom_id_1 == atom_id_3
            and angle.atom_id_3 == atom_id_1):
          return angle
    return None

  def get_chir_bonds_and_angles(self, m_i, m_j, chir):
    bonds = []
    for atom_i_comp_id,atom_id_i in [(chir.atom_1_comp_id, chir.atom_id_1),
                                     (chir.atom_2_comp_id, chir.atom_id_2),
                                     (chir.atom_3_comp_id, chir.atom_id_3)]:
      bonds.append(self.get_bond(
        m_i, m_j,
        chir.atom_centre_comp_id, chir.atom_id_centre,
        atom_i_comp_id, atom_id_i))
    angles = []
    for atom_i_comp_id,atom_id_i,atom_j_comp_id,atom_id_j in [
        (chir.atom_2_comp_id, chir.atom_id_2,
         chir.atom_3_comp_id, chir.atom_id_3),
        (chir.atom_1_comp_id, chir.atom_id_1,
         chir.atom_3_comp_id, chir.atom_id_3),
        (chir.atom_1_comp_id, chir.atom_id_1,
         chir.atom_2_comp_id, chir.atom_id_2)]:
      angles.append(self.get_angle(
        m_i, m_j,
        atom_i_comp_id, atom_id_i,
        chir.atom_centre_comp_id, chir.atom_id_centre,
        atom_j_comp_id, atom_id_j))
    return bonds, angles

  def get_chir_volume_ideal(self, m_i, m_j, chir):
    bonds, angles = self.get_chir_bonds_and_angles(m_i, m_j, chir)
    return get_chir_volume_ideal(chir.volume_sign, bonds, angles)

  def get_planes(self):
    return group_planes(self.plane_list)

  def as_geometry_restraints_motif_manipulation(self):
    result = geometry_restraints.motif_manipulation()
    result.id = if_none(self.chem_link.id, "")
    result.description = if_none(self.chem_link.name, "")
    if (self.source_info is not None):
      result.info.append(self.source_info)
    alts = []
    for bond in self.bond_list:
      a = geometry_restraints.motif_alteration(
        action="add", operand="bond")
      a.motif_ids.append(str(if_none(bond.atom_1_comp_id, "")))
      a.motif_ids.append(str(if_none(bond.atom_2_comp_id, "")))
      a.bond.atom_names = [
        if_none(bond.atom_id_1, ""),
        if_none(bond.atom_id_2, "")]
      a.bond.type = if_none(bond.type, "")
      a.bond.distance_ideal = if_none(bond.value_dist, 0)
      a.bond.weight = esd_as_weight(bond.value_dist_esd)
      alts.append(a)
    for angle in self.angle_list:
      a = geometry_restraints.motif_alteration(
        action="add", operand="angle")
      a.motif_ids.append(str(if_none(angle.atom_1_comp_id, "")))
      a.motif_ids.append(str(if_none(angle.atom_2_comp_id, "")))
      a.motif_ids.append(str(if_none(angle.atom_3_comp_id, "")))
      a.angle.atom_names = [
        if_none(angle.atom_id_1, ""),
        if_none(angle.atom_id_2, ""),
        if_none(angle.atom_id_3, "")]
      a.angle.angle_ideal = if_none(angle.value_angle, 0)
      a.angle.weight = esd_as_weight(angle.value_angle_esd)
      alts.append(a)
    for tor in self.tor_list:
      a = geometry_restraints.motif_alteration(
        action="add", operand="dihedral")
      a.motif_ids.append(str(if_none(tor.atom_1_comp_id, "")))
      a.motif_ids.append(str(if_none(tor.atom_2_comp_id, "")))
      a.motif_ids.append(str(if_none(tor.atom_3_comp_id, "")))
      a.motif_ids.append(str(if_none(tor.atom_4_comp_id, "")))
      a.dihedral.atom_names = [
        if_none(tor.atom_id_1, ""),
        if_none(tor.atom_id_2, ""),
        if_none(tor.atom_id_3, ""),
        if_none(tor.atom_id_4, "")]
      a.dihedral.angle_ideal = if_none(tor.value_angle, 0)
      a.dihedral.weight = esd_as_weight(tor.value_angle_esd)
      a.dihedral.periodicity = if_none(tor.period, 0)
      a.dihedral.id = if_none(tor.id, "")
      alts.append(a)
    for chir in self.chir_list:
      a = geometry_restraints.motif_alteration(
        action="add", operand="chirality")
      a.motif_ids.append(str(if_none(chir.atom_centre_comp_id, "")))
      a.motif_ids.append(str(if_none(chir.atom_1_comp_id, "")))
      a.motif_ids.append(str(if_none(chir.atom_2_comp_id, "")))
      a.motif_ids.append(str(if_none(chir.atom_3_comp_id, "")))
      a.chirality.atom_names = [
        if_none(chir.atom_id_centre, ""),
        if_none(chir.atom_id_1, ""),
        if_none(chir.atom_id_2, ""),
        if_none(chir.atom_id_3, "")]
      a.chirality.volume_sign = if_none(chir.volume_sign, 0)
      a.chirality.id = if_none(chir.id, "")
      alts.append(a)
    for plane in self.get_planes():
      a = geometry_restraints.motif_alteration(
        action="add", operand="planarity")
      for plane_atom in plane.plane_atoms:
        if (plane_atom.atom_comp_id is None):
          a.motif_ids.append("")
        else:
          a.motif_ids.append(str(plane_atom.atom_comp_id))
      a.planarity.atom_names = flex.std_string([if_none(plane_atom.atom_id, "")
        for plane_atom in plane.plane_atoms])
      a.planarity.weights = flex.double([esd_as_weight(plane_atom.dist_esd)
        for plane_atom in plane.plane_atoms])
      a.planarity.id = plane.plane_id
      alts.append(a)
    result.set_alterations(alts)
    return result

class chem_link(looped_data):
  """
_chem_link.id
_chem_link.comp_id_1
_chem_link.mod_id_1
_chem_link.group_comp_1
_chem_link.comp_id_2
_chem_link.mod_id_2
_chem_link.group_comp_2
_chem_link.name
  """

class chem_link_bond(looped_data):
  """
_chem_link_bond.atom_1_comp_id:int
_chem_link_bond.atom_id_1
_chem_link_bond.atom_2_comp_id:int
_chem_link_bond.atom_id_2
_chem_link_bond.type
_chem_link_bond.value_dist:float
_chem_link_bond.value_dist_esd:float
  """

class chem_link_angle(looped_data):
  """
_chem_link_angle.atom_1_comp_id:int
_chem_link_angle.atom_id_1
_chem_link_angle.atom_2_comp_id:int
_chem_link_angle.atom_id_2
_chem_link_angle.atom_3_comp_id:int
_chem_link_angle.atom_id_3
_chem_link_angle.value_angle:float
_chem_link_angle.value_angle_esd:float
  """

class chem_link_tor(looped_data):
  """
_chem_link_tor.id
_chem_link_tor.atom_1_comp_id:int
_chem_link_tor.atom_id_1
_chem_link_tor.atom_2_comp_id:int
_chem_link_tor.atom_id_2
_chem_link_tor.atom_3_comp_id:int
_chem_link_tor.atom_id_3
_chem_link_tor.atom_4_comp_id:int
_chem_link_tor.atom_id_4
_chem_link_tor.value_angle:float
_chem_link_tor.value_angle_esd:float
_chem_link_tor.period:int
  """

class chem_link_chir(looped_data):
  """
_chem_link_chir.id
_chem_link_chir.atom_centre_comp_id:int
_chem_link_chir.atom_id_centre
_chem_link_chir.atom_1_comp_id:int
_chem_link_chir.atom_id_1
_chem_link_chir.atom_2_comp_id:int
_chem_link_chir.atom_id_2
_chem_link_chir.atom_3_comp_id:int
_chem_link_chir.atom_id_3
_chem_link_chir.volume_sign
  """

class chem_link_plane(looped_data):
  """
_chem_link_plane.plane_id
_chem_link_plane.atom_comp_id:int
_chem_link_plane.atom_id
_chem_link_plane.dist_esd:float
  """

class mod_mod_id:

  def __init__(self, source_info, chem_mod):
    self.source_info = source_info
    self.chem_mod = chem_mod
    self.atom_list = []
    self.bond_list = []
    self.angle_list = []
    self.tor_list = []
    self.chir_list = []
    self.plane_atom_list = []

  def get_planes(self):
    return group_planes(self.plane_atom_list)

  def show(self, f=None):
    if (f is None): f = sys.stdout
    show_loop(data_list=self.atom_list, f=f)
    show_loop(data_list=self.bond_list, f=f)

  def as_geometry_restraints_motif_manipulation(self):
    result = geometry_restraints.motif_manipulation()
    result.id = if_none(self.chem_mod.id, "")
    result.description = if_none(self.chem_mod.name, "")
    if (self.source_info is not None):
      result.info.append(self.source_info)
    alts = []
    for atom in self.atom_list:
      a = geometry_restraints.motif_alteration(
        action=if_none(atom.function, ""), operand="atom")
      a.motif_ids.append("")
      if (a.action != "add"):
        a.motif_atom_name = if_none(atom.atom_id, "")
      a.atom.name = if_none(atom.new_atom_id, "")
      a.atom.scattering_type = if_none(atom.new_type_symbol, "")
      a.atom.nonbonded_type = if_none(atom.new_type_energy, "")
      a.atom.partial_charge = if_none(atom.new_partial_charge, 0)
      if (a.action == "change"):
        a.set_change_partial_charge(state=atom.new_partial_charge is not None)
      alts.append(a)
    for bond in self.bond_list:
      a = geometry_restraints.motif_alteration(
        action=if_none(bond.function, ""), operand="bond")
      a.motif_ids.append("")
      a.motif_ids.append("")
      a.bond.atom_names = [
        if_none(bond.atom_id_1, ""),
        if_none(bond.atom_id_2, "")]
      a.bond.type = if_none(bond.new_type, "")
      a.bond.distance_ideal = if_none(bond.new_value_dist, 0)
      a.bond.weight = esd_as_weight(bond.new_value_dist_esd)
      if (a.action == "change"):
        a.set_change_distance_ideal(state=bond.new_value_dist is not None)
        a.set_change_weight(state=bond.new_value_dist_esd is not None)
      alts.append(a)
    for angle in self.angle_list:
      a = geometry_restraints.motif_alteration(
        action=if_none(bond.function, ""), operand="angle")
      a.motif_ids.append("")
      a.motif_ids.append("")
      a.motif_ids.append("")
      a.angle.atom_names = [
        if_none(angle.atom_id_1, ""),
        if_none(angle.atom_id_2, ""),
        if_none(angle.atom_id_3, "")]
      a.angle.angle_ideal = if_none(angle.new_value_angle, 0)
      a.angle.weight = esd_as_weight(angle.new_value_angle_esd)
      if (a.action == "change"):
        a.set_change_angle_ideal(state=angle.new_value_angle is not None)
        a.set_change_weight(state=angle.new_value_angle_esd is not None)
      alts.append(a)
    for tor in self.tor_list:
      a = geometry_restraints.motif_alteration(
        action=if_none(tor.function, ""), operand="dihedral")
      a.motif_ids.append("")
      a.motif_ids.append("")
      a.motif_ids.append("")
      a.motif_ids.append("")
      a.dihedral.atom_names = [
        if_none(tor.atom_id_1, ""),
        if_none(tor.atom_id_2, ""),
        if_none(tor.atom_id_3, ""),
        if_none(tor.atom_id_4, "")]
      a.dihedral.angle_ideal = if_none(tor.new_value_angle, 0)
      a.dihedral.weight = esd_as_weight(tor.new_value_angle_esd)
      a.dihedral.periodicity = if_none(tor.new_period, 0)
      a.dihedral.id = if_none(tor.id, "")
      if (a.action == "change"):
        a.set_change_angle_ideal(state=tor.new_value_angle is not None)
        a.set_change_weight(state=tor.new_value_angle_esd is not None)
        a.set_change_periodicity(state=tor.new_period is not None)
      alts.append(a)
    for chir in self.chir_list:
      a = geometry_restraints.motif_alteration(
        action=if_none(chir.function, ""), operand="chirality")
      a.motif_ids.append("")
      a.motif_ids.append("")
      a.motif_ids.append("")
      a.motif_ids.append("")
      a.chirality.atom_names = [
        if_none(chir.atom_id_centre, ""),
        if_none(chir.atom_id_1, ""),
        if_none(chir.atom_id_2, ""),
        if_none(chir.atom_id_3, "")]
      a.chirality.volume_sign = if_none(chir.new_volume_sign, 0)
      a.chirality.id = if_none(chir.id, "")
      alts.append(a)
    for plane in self.get_planes():
      a = geometry_restraints.motif_alteration(
        action="change", operand="planarity")
      a.motif_ids.resize(len(plane.plane_atoms))
      a.set_planarity_atom_actions([if_none(plane_atom.function, "")
        for plane_atom in plane.plane_atoms])
      a.planarity.atom_names = flex.std_string([if_none(plane_atom.atom_id, "")
        for plane_atom in plane.plane_atoms])
      a.planarity.weights = flex.double([esd_as_weight(plane_atom.new_dist_esd)
        for plane_atom in plane.plane_atoms])
      a.planarity.id = plane.plane_id
      alts.append(a)
    result.set_alterations(alts)
    return result

class chem_mod(looped_data):
  """
_chem_mod.id
_chem_mod.name
_chem_mod.comp_id
_chem_mod.group_id
  """

class chem_mod_atom(looped_data):
  """
_chem_mod_atom.function
_chem_mod_atom.atom_id
_chem_mod_atom.new_atom_id
_chem_mod_atom.new_type_symbol
_chem_mod_atom.new_type_energy
_chem_mod_atom.new_partial_charge:float
  """

  def as_chem_comp(self):
    return chem_comp_atom(
      atom_id=self.new_atom_id,
      type_symbol=self.new_type_symbol,
      type_energy=self.new_type_energy,
      partial_charge=self.new_partial_charge)

def new_if_defined(old, new):
  if (new is None): return old
  if (new == ""): return old
  return new

class chem_mod_bond(looped_data):
  """
_chem_mod_bond.function
_chem_mod_bond.atom_id_1
_chem_mod_bond.atom_id_2
_chem_mod_bond.new_type
_chem_mod_bond.new_value_dist:float
_chem_mod_bond.new_value_dist_esd:float
  """

  def as_chem_comp(self):
    return chem_comp_bond(
      atom_id_1=self.atom_id_1,
      atom_id_2=self.atom_id_2,
      type=self.new_type,
      value_dist=self.new_value_dist,
      value_dist_esd=self.new_value_dist_esd)

  def is_matching_mod_for(self, bond):
    return (    self.atom_id_1 == bond.atom_id_1
            and self.atom_id_2 == bond.atom_id_2) \
        or (    self.atom_id_1 == bond.atom_id_2
            and self.atom_id_2 == bond.atom_id_1)

  def apply_change_in_place(self, bond):
    bond.type = new_if_defined(
      bond.type, self.new_type)
    bond.value_dist = new_if_defined(
      bond.value_dist, self.new_value_dist)
    bond.value_dist_esd = new_if_defined(
      bond.value_dist_esd, self.new_value_dist_esd)

class chem_mod_tree(looped_data):
  """
_chem_mod_tree.function
_chem_mod_tree.atom_id
_chem_mod_tree.atom_back
_chem_mod_tree.back_type
_chem_mod_tree.atom_forward
_chem_mod_tree.connect_type
  """

  def as_chem_comp(self):
    return chem_comp_tree(
      atom_id=self.atom_id,
      atom_back=self.atom_back,
      atom_forward=self.atom_forward,
      connect_type=self.connect_type)

  def is_matching_mod_for(self, tree_entry):
    return     self.atom_id == tree_entry.atom_id \
           and self.atom_back == tree_entry.atom_back \
           and self.atom_forward == tree_entry.atom_forward

  def apply_change_in_place(self, tree_entry):
    tree_entry.connect_type = new_if_defined(
      tree_entry.connect_type, self.connect_type)

class chem_mod_angle(looped_data):
  """
_chem_mod_angle.function
_chem_mod_angle.atom_id_1
_chem_mod_angle.atom_id_2
_chem_mod_angle.atom_id_3
_chem_mod_angle.new_value_angle:float
_chem_mod_angle.new_value_angle_esd:float
  """

  def as_chem_comp(self):
    return chem_comp_angle(
      atom_id_1=self.atom_id_1,
      atom_id_2=self.atom_id_2,
      atom_id_3=self.atom_id_3,
      value_angle=self.new_value_angle,
      value_angle_esd=self.new_value_angle_esd)

  def is_matching_mod_for(self, angle):
    if (self.atom_id_2 != angle.atom_id_2): return False
    return (    self.atom_id_1 == angle.atom_id_1
            and self.atom_id_3 == angle.atom_id_3) \
        or (    self.atom_id_1 == angle.atom_id_3
            and self.atom_id_3 == angle.atom_id_1)

  def apply_change_in_place(self, angle):
    angle.value_angle = new_if_defined(
      angle.value_angle, self.new_value_angle)
    angle.value_angle_esd = new_if_defined(
      angle.value_angle_esd, self.new_value_angle_esd)

class chem_mod_tor(looped_data):
  """
_chem_mod_tor.function
_chem_mod_tor.id
_chem_mod_tor.atom_id_1
_chem_mod_tor.atom_id_2
_chem_mod_tor.atom_id_3
_chem_mod_tor.atom_id_4
_chem_mod_tor.new_value_angle:float
_chem_mod_tor.new_value_angle_esd:float
_chem_mod_tor.new_period:int
  """

  def as_chem_comp(self):
    return chem_comp_tor(
      id=self.id,
      atom_id_1=self.atom_id_1,
      atom_id_2=self.atom_id_2,
      atom_id_3=self.atom_id_3,
      atom_id_4=self.atom_id_4,
      value_angle=self.new_value_angle,
      value_angle_esd=self.new_value_angle_esd,
      period=self.new_period)

  def is_matching_mod_for(self, tor):
    return     self.atom_id_1 == tor.atom_id_1 \
           and self.atom_id_2 == tor.atom_id_2 \
           and self.atom_id_3 == tor.atom_id_3 \
           and self.atom_id_4 == tor.atom_id_4

  def apply_change_in_place(self, tor):
    tor.value_angle = new_if_defined(
      tor.value_angle, self.new_value_angle)
    tor.value_angle_esd = new_if_defined(
      tor.value_angle_esd, self.new_value_angle_esd)
    tor.period = new_if_defined(
      tor.period, self.new_period)

class chem_mod_chir(looped_data):
  """
_chem_mod_chir.function
_chem_mod_chir.id
_chem_mod_chir.atom_id_centre
_chem_mod_chir.atom_id_1
_chem_mod_chir.atom_id_2
_chem_mod_chir.atom_id_3
_chem_mod_chir.new_volume_sign
  """

  def as_chem_comp(self):
    return chem_comp_chir(
      atom_id_centre=self.atom_id_centre,
      atom_id_1=self.atom_id_1,
      atom_id_2=self.atom_id_2,
      atom_id_3=self.atom_id_3,
      volume_sign=self.new_volume_sign)

  def is_matching_mod_for(self, chir):
    return     self.atom_id_centre == chir.atom_id_centre \
           and self.atom_id_1 == chir.atom_id_1 \
           and self.atom_id_2 == chir.atom_id_2 \
           and self.atom_id_3 == chir.atom_id_3

  def apply_change_in_place(self, chir):
    chir.volume_sign = new_if_defined(
      chir.volume_sign, self.new_volume_sign)

class chem_mod_plane_atom(looped_data):
  """
_chem_mod_plane_atom.function
_chem_mod_plane_atom.plane_id
_chem_mod_plane_atom.atom_id
_chem_mod_plane_atom.new_dist_esd:float
  """

  def as_chem_comp(self):
    return chem_comp_plane_atom(
      plane_id=self.plane_id,
      atom_id=self.atom_id,
      dist_esd=self.new_dist_esd)

  def is_matching_mod_for(self, plane_atom):
    return     self.plane_id == plane_atom.plane_id \
           and self.atom_id == plane_atom.atom_id

  def apply_change_in_place(self, plane_atom):
    plane_atom.dist_esd = new_if_defined(
      plane_atom.dist_esd, self.new_dist_esd)

class energy_lib_synonym(looped_data):
  """
_lib_synonym.atom_type
_lib_synonym.atom_alternative_type
  """

class energy_lib_atom(looped_data):
  """
_lib_atom.type
_lib_atom.weight:float
_lib_atom.hb_type
_lib_atom.vdw_radius:float
_lib_atom.vdwh_radius:float
_lib_atom.ion_radius:float
_lib_atom.element
_lib_atom.valency:int
_lib_atom.sp:int
  """

class energy_lib_vdw(looped_data):
  """
_lib_vdw.atom_type_1
_lib_vdw.atom_type_2
_lib_vdw.energy_min:float
_lib_vdw.radius_min:float
_lib_vdw.H_flag
  """
