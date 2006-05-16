from __future__ import generators
from mmtbx.monomer_library import cif_types
from mmtbx.monomer_library import mmCIF
from scitbx.python_utils import dicts
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, format_exception, windows_device_names
import libtbx.load_env
import libtbx.path
import copy
import os

try: import sets
except ImportError: pass # XXX Python 2.2 compatibility
else: windows_device_names = sets.Set(windows_device_names)

class MonomerLibraryServerError(RuntimeError): pass

mon_lib_env_vars = ["MMTBX_CCP4_MONOMER_LIB", "CLIBD_MON"]

def load_mon_lib_file(mon_lib_path, relative_path_components=[]):
  if (mon_lib_path is not None):
    cif_path = os.path.join(mon_lib_path, *relative_path_components)
    if (os.path.isfile(cif_path)):
      return cif_path
  return None

def find_mon_lib_file(env_vars=mon_lib_env_vars, relative_path_components=[]):
  result = load_mon_lib_file(
    mon_lib_path=os.environ.get(env_vars[0], None),
    relative_path_components=relative_path_components)
  if (result is not None): return result
  for relative_path in ["mon_lib", "ext_ref_files/mon_lib"]:
    result = load_mon_lib_file(
      mon_lib_path=libtbx.env.find_in_repositories(
        relative_path=relative_path),
      relative_path_components=relative_path_components)
    if (result is not None): return result
  for env_var in env_vars[1:]:
    result = load_mon_lib_file(
      mon_lib_path=os.environ.get(env_var, None),
      relative_path_components=relative_path_components)
    if (result is not None): return result
  return None

class mon_lib_cif_loader(object):

  def __init__(self, path=None, relative_path_components=[], strict=False):
    self.path = path
    if (self.path is None):
      self.path = find_mon_lib_file(
        relative_path_components=relative_path_components)
      if (self.path is None):
        raise MonomerLibraryServerError(
          "Cannot find CCP4 monomer library."
          " Please define one of these environment variables: "
          + ", ".join(mon_lib_env_vars))
    self.cif = mmCIF.mmCIFFile()
    self.cif.load_file(self.path, strict=strict)

def mon_lib_list_cif(path=None, strict=False):
  return mon_lib_cif_loader(
    path=path,
    relative_path_components=["list", "mon_lib_list.cif"],
    strict=strict)

def mon_lib_ener_lib_cif(path=None, strict=False):
  return mon_lib_cif_loader(
    path=path,
    relative_path_components=["ener_lib.cif"],
    strict=strict)

class trivial_html_tag_filter(object):

  def __init__(self, file_name):
    self.f = iter(open(file_name))

  def next(self):
    while 1:
      result = self.f.next()
      if (result[0] != "<"):
        return result

  def __iter__(self):
    return self

def read_cif(file_name):
  cif_object = mmCIF.mmCIFFile()
  cif_object.load_file(fil=trivial_html_tag_filter(file_name), strict=False)
  return cif_object

def convert_list_block(
      source_info,
      cif_object,
      list_name,
      list_item_name,
      data_prefix,
      cif_type_inner,
      cif_type_outer,
      outer_mappings):
  tabulated_items = {}
  list_block = cif_object.get(list_name)
  if (list_block is not None):
    list_item_block = list_block.get(list_item_name)
    if (list_item_block is not None):
      for item in list_item_block:
        obj_inner = cif_type_inner(**item)
        tabulated_items[obj_inner.id] = obj_inner
  for cif_data in cif_object:
    if (cif_data.name.startswith(data_prefix)):
      item_id = cif_data.name[len(data_prefix):]
      if (data_prefix + item_id == list_name): continue
      obj_inner = tabulated_items.get(item_id)
      if (obj_inner is None):
        obj_inner = cif_type_inner(id=item_id)
      obj_outer = None
      for loop_block,lst_name in outer_mappings:
        rows = cif_data.get(loop_block)
        if (rows is None): continue
        if (obj_outer is None):
          obj_outer = cif_type_outer(source_info, obj_inner)
        lst = getattr(obj_outer, lst_name)
        typ = getattr(cif_types, loop_block)
        for row in rows:
          lst.append(typ(**row))
      if (obj_outer is not None):
        yield obj_outer

def get_rows(cif_object, data_name, table_name):
  cif_data = cif_object.get(data_name)
  if (cif_data is None): return []
  return cif_data.get(table_name, [])

class process_cif_mixin(object):

  def process_cif_object(self, cif_object, file_name=None):
    if (file_name is None):
      source_info = None
    else:
      source_info = "file: "+file_name
    try: self.convert_all(source_info=source_info, cif_object=cif_object)
    except KeyboardInterrupt: raise
    except:
      if (file_name is None): file_name = "(file name not available)"
      raise Sorry(
        "Error processing CIF file:\n"
        "  %s\n"
        "  (%s)" % (show_string(file_name), format_exception()))

  def process_cif(self, file_name):
    try: cif_object = read_cif(file_name=file_name)
    except KeyboardInterrupt: raise
    except:
      raise Sorry(
        "Error reading CIF file:\n"
        "  %s\n"
        "  (%s)" % (show_string(file_name), format_exception()))
    self.process_cif_object(cif_object=cif_object, file_name=file_name)

class server(process_cif_mixin):

  def __init__(self, list_cif=None):
    if (list_cif is None):
      list_cif = mon_lib_list_cif()
    self.root_path = os.path.dirname(os.path.dirname(list_cif.path))
    self.deriv_list_dict = {}
    self.comp_synonym_list_dict = {}
    self.comp_synonym_atom_list_dict = dicts.with_default_factory(dict)
    self.comp_comp_id_dict = {}
    self.link_link_id_list = []
    self.link_link_id_dict = {}
    self.mod_mod_id_list = []
    self.mod_mod_id_dict = {}
    self.convert_all(
      source_info="file: "+list_cif.path,
      cif_object=list_cif.cif, skip_comp_list=True)
    self._create_rna_dna_placeholders()

  def convert_all(self, source_info, cif_object, skip_comp_list=False):
    self.convert_deriv_list_dict(cif_object=cif_object)
    self.convert_comp_synonym_list(cif_object=cif_object)
    self.convert_comp_synonym_atom_list(cif_object=cif_object)
    if (not skip_comp_list):
      self.convert_comp_list(
        source_info=source_info, cif_object=cif_object)
    self.convert_link_list(
      source_info=source_info, cif_object=cif_object)
    self.convert_mod_list(
      source_info=source_info, cif_object=cif_object)

  def convert_deriv_list_dict(self, cif_object):
    for row in get_rows(cif_object,
                 "deriv_list", "chem_comp_deriv"):
      deriv = cif_types.chem_comp_deriv(**row)
      self.deriv_list_dict[deriv.comp_id] = deriv

  def convert_comp_synonym_list(self, cif_object):
    for row in get_rows(cif_object,
                 "comp_synonym_list", "chem_comp_synonym"):
      self.comp_synonym_list_dict[row["comp_alternative_id"]] = row["comp_id"]

  def convert_comp_synonym_atom_list(self, cif_object):
    for row in get_rows(cif_object,
                 "comp_synonym_atom_list", "chem_comp_synonym_atom"):
      synonym = cif_types.chem_comp_synonym_atom(**row)
      d = self.comp_synonym_atom_list_dict[synonym.comp_id]
      d[synonym.atom_alternative_id] = synonym.atom_id
      if (synonym.comp_alternative_id != ""):
        d = self.comp_synonym_atom_list_dict[synonym.comp_alternative_id]
        d[synonym.atom_alternative_id] = synonym.atom_id

  def convert_comp_list(self, source_info, cif_object):
    for comp_comp_id in convert_list_block(
                          source_info=source_info,
                          cif_object=cif_object,
                          list_name="comp_list",
                          list_item_name="chem_comp",
                          data_prefix="comp_",
                          cif_type_inner=cif_types.chem_comp,
                          cif_type_outer=cif_types.comp_comp_id,
                          outer_mappings=[
                           ("chem_comp_atom","atom_list"),
                           ("chem_comp_tree","tree_list"),
                           ("chem_comp_bond","bond_list"),
                           ("chem_comp_angle","angle_list"),
                           ("chem_comp_tor","tor_list"),
                           ("chem_comp_chir","chir_list"),
                           ("chem_comp_plane_atom","plane_atom_list")]):
      chem_comp = comp_comp_id.chem_comp
      self.comp_comp_id_dict[chem_comp.id.strip().upper()] = comp_comp_id
      tlc = chem_comp.three_letter_code
      if (tlc is not None):
        tlc = tlc.strip()
        if (1 <= len(tlc) <= 3):
          self.comp_comp_id_dict[tlc.upper()] = comp_comp_id

  def convert_link_list(self, source_info, cif_object):
    for link_link_id in convert_list_block(
                          source_info=source_info,
                          cif_object=cif_object,
                          list_name="link_list",
                          list_item_name="chem_link",
                          data_prefix="link_",
                          cif_type_inner=cif_types.chem_link,
                          cif_type_outer=cif_types.link_link_id,
                          outer_mappings=[
                           ("chem_link_bond","bond_list"),
                           ("chem_link_angle","angle_list"),
                           ("chem_link_tor","tor_list"),
                           ("chem_link_chir","chir_list"),
                           ("chem_link_plane","plane_list")]):
      self.link_link_id_list.append(link_link_id)
      self.link_link_id_dict[link_link_id.chem_link.id] = link_link_id

  def convert_mod_list(self, source_info, cif_object):
    for mod_mod_id in convert_list_block(
                        source_info=source_info,
                        cif_object=cif_object,
                        list_name="mod_list",
                        list_item_name="chem_mod",
                        data_prefix="mod_",
                        cif_type_inner=cif_types.chem_mod,
                        cif_type_outer=cif_types.mod_mod_id,
                        outer_mappings=[
                          ("chem_mod_atom","atom_list"),
                          ("chem_mod_tree","tree_list"),
                          ("chem_mod_bond","bond_list"),
                          ("chem_mod_angle","angle_list"),
                          ("chem_mod_tor","tor_list"),
                          ("chem_mod_chir","chir_list"),
                          ("chem_mod_plane_atom","plane_atom_list")]):
      self.mod_mod_id_list.append(mod_mod_id)
      self.mod_mod_id_dict[mod_mod_id.chem_mod.id] = mod_mod_id

  def get_comp_comp_id(self, comp_id):
    comp_id = comp_id.strip().upper()
    if (len(comp_id) == 0): return None
    try: return self.comp_comp_id_dict[comp_id]
    except KeyError: pass
    std_comp_id = self.comp_synonym_list_dict.get(comp_id, "").strip().upper()
    def find_file():
      for i_pass in [0,1]:
        for trial_comp_id in [std_comp_id, comp_id]:
          if (len(trial_comp_id) == 0): continue
          dir_name = os.path.join(self.root_path, trial_comp_id[0].lower())
          if (os.path.isdir(dir_name)):
            if (trial_comp_id in windows_device_names):
              cif_name = "%s_%s.cif" % (trial_comp_id, trial_comp_id)
            else:
              cif_name = trial_comp_id + ".cif"
            if (i_pass == 0):
              file_name = os.path.join(dir_name, cif_name)
              if (os.path.isfile(file_name)):
                return file_name
            else:
              cif_name = cif_name.lower()
              for node in os.listdir(dir_name):
                if (node.lower() != cif_name): continue
                return os.path.join(dir_name, node)
      return None
    file_name = find_file()
    if (file_name is None): return None
    self.process_cif(file_name=file_name)
    if (len(std_comp_id) > 0):
      result = self.comp_comp_id_dict.get(std_comp_id)
    else:
      result = None
    if (result is not None):
      self.comp_comp_id_dict[comp_id] = result
    else:
      result = self.comp_comp_id_dict.get(comp_id)
      if (result is not None):
        self.comp_comp_id_dict[std_comp_id] = result
      else:
        or_std_comp_id = ""
        if (len(std_comp_id) > 0):
          or_std_comp_id = "or %s" % show_string(std_comp_id)
        raise Sorry(
          "Monomer library file %s does not define comp_id %s%s" % (
            show_string(file_name), show_string(comp_id), or_std_comp_id))
    return result

  def _create_rna_dna_placeholders(self):
    for base_code in ["A", "C", "G"]:
      rna = self.get_comp_comp_id(base_code+"r")
      dna = self.get_comp_comp_id(base_code+"d")
      chem_comp = cif_types.chem_comp(
        id=base_code+"?",
        three_letter_code=None,
        name=None,
        group="rna_dna_placeholder",
        number_atoms_all=None,
        number_atoms_nh=None,
        desc_level=None)
      comp_comp_id = cif_types.comp_comp_id(
        source_info=None, chem_comp=chem_comp)
      for atom in rna.atom_list:
        comp_comp_id.atom_list.append(copy.copy(atom))
      rna_atom_dict = rna.atom_dict()
      for atom in dna.atom_list:
        if (not rna.atom_dict().has_key(atom.atom_id)):
          comp_comp_id.atom_list.append(copy.copy(atom))
      self.comp_comp_id_dict[base_code] = comp_comp_id
      self.comp_comp_id_dict["+"+base_code] = comp_comp_id

class ener_lib(process_cif_mixin):

  def __init__(self, ener_lib_cif=None, source_info=None):
    if (ener_lib_cif is None):
      ener_lib_cif = mon_lib_ener_lib_cif()
    self.lib_synonym = {}
    self.lib_atom = {}
    self.lib_vdw = []
    self.convert_all(source_info=source_info, cif_object=ener_lib_cif.cif)
    self.source_infos = []

  def convert_all(self, source_info, cif_object):
    if (source_info is not None):
      self.source_infos.append(source_info)
    self.convert_lib_synonym(cif_object=cif_object)
    self.convert_lib_atom(cif_object=cif_object)
    self.convert_lib_vdw(cif_object=cif_object)

  def convert_lib_synonym(self, cif_object):
    for row in get_rows(cif_object, "energy", "lib_synonym"):
      syn = cif_types.energy_lib_synonym(**row)
      self.lib_synonym[syn.atom_alternative_type] = syn.atom_type

  def convert_lib_atom(self, cif_object):
    for row in get_rows(cif_object, "energy", "lib_atom"):
      entry = cif_types.energy_lib_atom(**row)
      self.lib_atom[entry.type] = entry

  def convert_lib_vdw(self, cif_object):
    for row in get_rows(cif_object, "energy", "lib_vdw"):
      vdw = cif_types.energy_lib_vdw(**row)
      self.lib_vdw.append(vdw)

  def vdw_lookup(self, atom_energy_types_pair):
    atom_energy_types_pair = [self.lib_synonym.get(t,t)
      for t in atom_energy_types_pair]
    for vdw in self.lib_vdw:
      if (   (    vdw.atom_type_1 == atom_energy_types_pair[0]
              and vdw.atom_type_2 == atom_energy_types_pair[1])
          or (    vdw.atom_type_1 == atom_energy_types_pair[1]
              and vdw.atom_type_2 == atom_energy_types_pair[0])):
        return vdw.radius_min
    entries = [self.lib_atom.get(t,None) for t in atom_energy_types_pair]
    if (None not in entries):
      radii = [entry.vdw_radius for entry in entries]
      if (None not in radii):
        return radii[0] + radii[1]
    return None
