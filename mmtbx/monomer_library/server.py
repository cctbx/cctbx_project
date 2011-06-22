from mmtbx.monomer_library import cif_types
import iotbx.cif
from iotbx.pdb import residue_name_plus_atom_names_interpreter
from scitbx.python_utils import dicts
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, format_exception, windows_device_names
import libtbx.load_env
import libtbx.path
import os

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
  for relative_path in [
        "chem_data/geostd",
        "chem_data/mon_lib",
        "mon_lib",
        "ext_ref_files/mon_lib"]:
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
    self.cif = read_cif(file_name=self.path)

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
  return iotbx.cif.reader(file_path=file_name, strict=False).model()

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
      for row in list_item_block.iterrows():
        obj_inner = cif_type_inner(**row)
        tabulated_items[obj_inner.id] = obj_inner
  for key, cif_data in cif_object.iteritems():
    if (key.startswith(data_prefix)):
      item_id = key[len(data_prefix):]
      if (data_prefix + item_id == list_name): continue
      obj_inner = tabulated_items.get(item_id)
      if (obj_inner is None):
        obj_inner = cif_type_inner(id=item_id)
      obj_outer = None
      for loop_block,lst_name in outer_mappings:
        rows = cif_data.get('_'+loop_block)
        if (rows is None): continue
        if (obj_outer is None):
          obj_outer = cif_type_outer(source_info, obj_inner)
        lst = getattr(obj_outer, lst_name)
        typ = getattr(cif_types, loop_block)
        for row in rows.iterrows():
          lst.append(typ(**row))
      if (obj_outer is not None):
        yield obj_outer

def convert_comp_list(source_info, cif_object):
  for comp_comp_id in convert_list_block(
                        source_info=source_info,
                        cif_object=cif_object,
                        list_name="comp_list",
                        list_item_name="_chem_comp",
                        data_prefix="comp_",
                        cif_type_inner=cif_types.chem_comp,
                        cif_type_outer=cif_types.comp_comp_id,
                        outer_mappings=[
                         ("chem_comp_atom","atom_list"),
                         ("chem_comp_bond","bond_list"),
                         ("chem_comp_angle","angle_list"),
                         ("chem_comp_tor","tor_list"),
                         ("chem_comp_chir","chir_list"),
                         ("chem_comp_plane_atom","plane_atom_list"),
                         ("chem_comp_rotamer_info",
                          "rotamer_info_phil_str_list")]):
    yield comp_comp_id

def convert_link_list(source_info, cif_object):
  for link_link_id in convert_list_block(
                        source_info=source_info,
                        cif_object=cif_object,
                        list_name="link_list",
                        list_item_name="_chem_link",
                        data_prefix="link_",
                        cif_type_inner=cif_types.chem_link,
                        cif_type_outer=cif_types.link_link_id,
                        outer_mappings=[
                         ("chem_link_bond","bond_list"),
                         ("chem_link_angle","angle_list"),
                         ("chem_link_tor","tor_list"),
                         ("chem_link_chir","chir_list"),
                         ("chem_link_plane","plane_list")]):
    yield link_link_id

def convert_mod_list(source_info, cif_object):
  for mod_mod_id in convert_list_block(
                      source_info=source_info,
                      cif_object=cif_object,
                      list_name="mod_list",
                      list_item_name="_chem_mod",
                      data_prefix="mod_",
                      cif_type_inner=cif_types.chem_mod,
                      cif_type_outer=cif_types.mod_mod_id,
                      outer_mappings=[
                        ("chem_mod_atom","atom_list"),
                        ("chem_mod_bond","bond_list"),
                        ("chem_mod_angle","angle_list"),
                        ("chem_mod_tor","tor_list"),
                        ("chem_mod_chir","chir_list"),
                        ("chem_mod_plane_atom","plane_atom_list")]):
    yield mod_mod_id

def get_rows(cif_object, data_name, table_name):
  cif_data = cif_object.get(data_name)
  if (cif_data is None): return []
  table = cif_data.get(table_name)
  if table is not None: return table.iterrows()
  else: return []

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
    self.geostd_path = os.path.join(os.path.dirname(self.root_path), "geostd")
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
    self.comp_comp_id_mod_dict = {}
    self.process_geostd_rna_dna()

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
                 "deriv_list", "_chem_comp_deriv"):
      deriv = cif_types.chem_comp_deriv(**row)
      self.deriv_list_dict[deriv.comp_id] = deriv

  def convert_comp_synonym_list(self, cif_object):
    for row in get_rows(cif_object,
                 "comp_synonym_list", "_chem_comp_synonym"):
      self.comp_synonym_list_dict[row["_chem_comp_synonym.comp_alternative_id"]] \
          = row["_chem_comp_synonym.comp_id"]

  def convert_comp_synonym_atom_list(self, cif_object):
    for row in get_rows(cif_object,
                 "comp_synonym_atom_list", "_chem_comp_synonym_atom"):
      synonym = cif_types.chem_comp_synonym_atom(**row)
      d = self.comp_synonym_atom_list_dict[synonym.comp_id]
      d[synonym.atom_alternative_id] = synonym.atom_id
      if (synonym.comp_alternative_id != ""):
        d = self.comp_synonym_atom_list_dict[synonym.comp_alternative_id]
        d[synonym.atom_alternative_id] = synonym.atom_id

  def convert_comp_list(self, source_info, cif_object):
    for comp_comp_id in convert_comp_list(
                          source_info=source_info, cif_object=cif_object):
      comp_comp_id.normalize_atom_ids_in_place()
      chem_comp = comp_comp_id.chem_comp
      self.comp_comp_id_dict[chem_comp.id.strip().upper()] = comp_comp_id

  def convert_link_list(self, source_info, cif_object):
    for link_link_id in convert_link_list(
                          source_info=source_info, cif_object=cif_object):
      self.link_link_id_list.append(link_link_id)
      self.link_link_id_dict[link_link_id.chem_link.id] = link_link_id

  def convert_mod_list(self, source_info, cif_object):
    for mod_mod_id in convert_mod_list(
                        source_info=source_info, cif_object=cif_object):
      self.mod_mod_id_list.append(mod_mod_id)
      self.mod_mod_id_dict[mod_mod_id.chem_mod.id] = mod_mod_id

  def process_geostd_rna_dna(self):
    if (not os.path.isdir(self.geostd_path)): return
    for file_name in [
          "chain_link_rna2p.cif",
          "chain_link_rna3p.cif",
          "mod_rna2p.cif",
          "mod_rna3p.cif"]:
      self.process_cif(
        file_name=os.path.join(self.geostd_path, "rna_dna", file_name))

  def get_comp_comp_id_direct(self, comp_id):
    comp_id = comp_id.strip().upper()
    if (len(comp_id) == 0): return None
    result = self.comp_comp_id_dict.get(comp_id)
    if (result is not None):
      return result
    std_comp_id = self.comp_synonym_list_dict.get(comp_id, "").strip().upper()
    def find_file():
      for i_pass in [0,1]:
        for trial_comp_id in [std_comp_id, comp_id]:
          if (len(trial_comp_id) == 0): continue
          dir_name = os.path.join(self.geostd_path, trial_comp_id[0].lower())
          # check the Geo Standard
          if (os.path.isdir(dir_name)):
            cif_name = "data_" + trial_comp_id + ".cif"
            if (i_pass == 0):
              file_name = os.path.join(dir_name, cif_name)
              if (os.path.isfile(file_name)):
                return file_name
            else:
              cif_name = cif_name.lower()
              for node in os.listdir(dir_name):
                if (node.lower() != cif_name): continue
                return os.path.join(dir_name, node)
          # check PHENIX Mon Lib
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
        if (len(std_comp_id) != 0):
          self.comp_comp_id_dict[std_comp_id] = result
      else:
        or_std_comp_id = ""
        if (len(std_comp_id) > 0):
          or_std_comp_id = "or %s" % show_string(std_comp_id)
        raise Sorry(
          "Monomer library file %s does not define comp_id %s%s" % (
            show_string(file_name), show_string(comp_id), or_std_comp_id))
    return result

  def get_comp_comp_id_and_atom_name_interpretation(self,
        residue_name,
        atom_names,
        translate_cns_dna_rna_residue_names=None):
    rnpani = residue_name_plus_atom_names_interpreter(
      residue_name=residue_name,
      atom_names=atom_names,
      translate_cns_dna_rna_residue_names=translate_cns_dna_rna_residue_names,
      return_mon_lib_dna_name=True)
    if (rnpani.work_residue_name is None): return (None, None)
    d_aa_rn = getattr(
      rnpani.atom_name_interpretation, "d_aa_residue_name", None)
    if (    d_aa_rn is not None
        and self.comp_comp_id_dict.get(d_aa_rn) is not None):
      return (self.get_comp_comp_id_direct(comp_id=d_aa_rn), None)
    return (
      self.get_comp_comp_id_direct(comp_id=rnpani.work_residue_name),
      rnpani.atom_name_interpretation)

  def rotamer_iterator(self,
        comp_id,
        atom_names,
        sites_cart,
        fine_sampling=False):
    comp_comp_id = self.get_comp_comp_id_direct(comp_id=comp_id)
    if (comp_comp_id is None) :
      return None
    return comp_comp_id.rotamer_iterator(
      atom_names=atom_names,
      sites_cart=sites_cart,
      fine_sampling=fine_sampling)

  def get_comp_comp_id_mod(self, comp_comp_id, mod_ids):
    key = "%".join((comp_comp_id.chem_comp.id,) + mod_ids)
    result = self.comp_comp_id_mod_dict.get(key)
    if (result is None):
      mod_comp_comp_id = comp_comp_id
      chem_mod_ids = []
      for mod_id in mod_ids:
        mod_mod_id = self.mod_mod_id_dict[mod_id]
        mod_comp_comp_id = mod_comp_comp_id.apply_mod(mod=mod_mod_id)
        chem_mod_ids.append(mod_mod_id.chem_mod.id)
      result = (mod_comp_comp_id, chem_mod_ids)
      self.comp_comp_id_mod_dict[key] = result
    return result

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
    for row in get_rows(cif_object, "energy", "_lib_synonym"):
      syn = cif_types.energy_lib_synonym(**row)
      self.lib_synonym[syn.atom_alternative_type] = syn.atom_type

  def convert_lib_atom(self, cif_object):
    for row in get_rows(cif_object, "energy", "_lib_atom"):
      entry = cif_types.energy_lib_atom(**row)
      self.lib_atom[entry.type] = entry

  def convert_lib_vdw(self, cif_object):
    for row in get_rows(cif_object, "energy", "_lib_vdw"):
      vdw = cif_types.energy_lib_vdw(**row)
      self.lib_vdw.append(vdw)
