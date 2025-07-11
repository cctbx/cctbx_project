from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import cif_types
import iotbx.cif
from iotbx.pdb import residue_name_plus_atom_names_interpreter
from scitbx.python_utils import dicts
from libtbx.str_utils import show_string
from libtbx.utils import Sorry, format_exception, windows_device_names
import libtbx.load_env
import libtbx.path
import os
import six

class MonomerLibraryServerError(RuntimeError): pass

# CLIBD_MON no longer supported because of the cascade to GeoStd
mon_lib_env_vars = ["MMTBX_CCP4_MONOMER_LIB", "CLIBD_MON"]

def load_mon_lib_file(mon_lib_path,
                      relative_path_components=[],
                     ):
  if (mon_lib_path is not None):
    cif_path = os.path.join(mon_lib_path, *relative_path_components)
    if (os.path.isfile(cif_path)):
      return cif_path
  return None

def find_mon_lib_file(env_vars=mon_lib_env_vars,
                      relative_path_components=[],
                      ):
  redirect_dir=os.environ.get(env_vars[0], None)
  result = load_mon_lib_file(
    mon_lib_path=redirect_dir,
    relative_path_components=relative_path_components)
  if (result is not None): return result
  relative_paths = [
    "chem_data/geostd",
    "chem_data/mon_lib",
    "mon_lib",
    'geostd',
    "ext_ref_files/mon_lib"]
  if ('mon_lib_list.cif' in relative_path_components
      # or 'ener_lib.cif' in relative_path_components
      ):
    relative_paths.reverse()
  if redirect_dir and 'geostd_list.cif' in relative_path_components: return None
  for relative_path in relative_paths:
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

  def __init__(self,
               path=None,
               relative_path_components=[],
               strict=False,
               ):
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

def geostd_list_cif(path=None, strict=False):
  try:
    return mon_lib_cif_loader(
      path=path,
      relative_path_components=["list", "geostd_list.cif"],
      strict=strict)
  except MonomerLibraryServerError as e:
    return None

def mon_lib_ener_lib_cif(path=None, strict=False):
  # actually does both
  mon_lib_ener = mon_lib_cif_loader(
    path=path,
    relative_path_components=["ener_lib.cif"],
    strict=strict)
  return mon_lib_ener
  try:
    geostd_ener = mon_lib_cif_loader(
      path=path,
      relative_path_components=["geostd_ener_lib.cif"],
      strict=strict)
  except MonomerLibraryServerError:
    geostd_ener = None
  if geostd_ener:
    mon_lib_ener = merge_and_overwrite_cifs(geostd_ener, mon_lib_ener)
  return mon_lib_ener

class trivial_html_tag_filter(object):

  def __init__(self, file_name):
    self.f = iter(open(file_name))

  def next(self):
    while 1:
      result = next(self.f)
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
        obj_inner = cif_type_inner(**dict(row))
        tabulated_items[obj_inner.id] = obj_inner
  for key, cif_data in six.iteritems(cif_object):
    if (key.startswith(data_prefix)):
      item_id = key[len(data_prefix):]
      if (data_prefix + item_id == list_name): continue
      obj_inner = tabulated_items.get(item_id)
      if (obj_inner is None):
        obj_inner = cif_type_inner(id=item_id)
      obj_outer = None
      for loop_block,lst_name in outer_mappings:
        rows = cif_data.get('_'+loop_block)
        if (rows is None):
          d = dict((k, v) for k, v in six.iteritems(cif_data)
                   if k.startswith("_"+loop_block))
          if len(d) > 0:
            lst = getattr(obj_outer, lst_name)
            typ = getattr(cif_types, loop_block)
            lst.append(typ(**dict(d)))
          continue
        if (obj_outer is None):
          obj_outer = cif_type_outer(source_info, obj_inner)
        lst = getattr(obj_outer, lst_name)
        typ = getattr(cif_types, loop_block)
        for row in rows.iterrows():
          lst.append(typ(**dict(row)))
      if (obj_outer is not None):
        obj_outer.cif_object = cif_data
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

def merge_and_overwrite_cifs(geostd_list_cif_obj, list_cif, verbose=False):
  mon_lib_cif = list_cif.cif
  if geostd_list_cif_obj is not None:
    geostd_cif = geostd_list_cif_obj.cif
    for geostd_block_key in geostd_cif:
      geostd_block = geostd_cif[geostd_block_key]
      mon_lib_block = mon_lib_cif.get(geostd_block_key)
      if verbose:
        print('+'*80)
        print(geostd_block_key)
        print('-'*80)
        print(mon_lib_block)
        print('+'*80)
      replace_block=False
      merge_block=False
      if geostd_block_key.endswith("_list"):
        pass
      elif geostd_block_key.endswith("energy"):
        merge_block=True
      else:
        replace_block=True
      if not mon_lib_block:
        replace_block=True
      if merge_block:
        def _row_match(row1, row2, key):
          match_keys = {
            "_lib_atom" : ["_lib_atom.type"],
            }
          assert key in match_keys
          k = match_keys[key]
          count=0
          for k1 in row1:
            if not k1 in k: continue
            if k1 not in row2: continue
            if row1[k1]==row2[k1]: count+=1
          if count==len(k): return True
          return False

        geostd_loop_keys = list(geostd_block.keys())
        if not geostd_loop_keys: continue
        done = []
        for geostd_loop_key in geostd_loop_keys:
          if geostd_loop_key.find(".")==-1: continue
          geostd_loop_key = geostd_loop_key.split(".")[0]
          if geostd_loop_key in done: continue
          done.append(geostd_loop_key)
          geostd_loop = geostd_block.get(geostd_loop_key)
          mon_lib_loop = mon_lib_block.get(geostd_loop_key)
          if not mon_lib_loop: continue
          if verbose: print('''
%s
  Merging rows to CIF
    block "%s"
      loop "%s
%s
''' % (
          "-"*80,
          geostd_block_key,
          geostd_loop_key,
          "-"*80,
          ))
          for row1 in geostd_loop.iterrows():
            if verbose:
              for key in row1:
                print("    %-20s : %s" % (
                  key.replace("%s." % geostd_loop_key, ""),
                  row1[key],
                  ))
              print()
            for i, row2 in enumerate(mon_lib_loop.iterrows()):
              if _row_match(row1, row2, geostd_loop_key):
                mon_lib_loop.delete_row(i)
                mon_lib_loop.add_row(list(row1.values()))
      elif not replace_block:
        # add to loops row-wise
        geostd_loop_keys = list(geostd_block.keys())
        if not geostd_loop_keys: continue
        done = []
        for geostd_loop_key in geostd_loop_keys:
          if geostd_loop_key.find(".")==-1: continue
          geostd_loop_key = geostd_loop_key.split(".")[0]
          if geostd_loop_key in done: continue
          done.append(geostd_loop_key)
          geostd_loop = geostd_block.get(geostd_loop_key)
          mon_lib_loop = mon_lib_block.get(geostd_loop_key)
          if not mon_lib_loop: continue
          if verbose: print('''
%s
  Adding rows to CIF
    block "%s"
      loop "%s
%s
''' % (
          "-"*80,
          geostd_block_key,
          geostd_loop_key,
          "-"*80,
          ))
          for row in geostd_loop.iterrows():
            if verbose:
              for key in row:
                print("    %-20s : %s" % (
                  key.replace("%s." % geostd_loop_key, ""),
                  row[key],
                  ))
              print()
            mon_lib_loop.add_row(list(row.values()))
      else:
        # add block
        if verbose:
          print('\n%s\n  Adding/replacing CIF block\n%s' % ("-"*80, "-"*80))
          print(geostd_block)
        mon_lib_cif[geostd_block_key] = geostd_block
      if verbose:
        print('+'*80)
        print(mon_lib_block)
        print('+'*80)
  return list_cif

def print_filtered_cif(cif):
  outl = ""
  for line in str(cif).split("\n"):
    if line.find("_")==-1:
      #outl += "%s\n" % line
      continue
    outl += "%s\n" % line
  print(outl)
  #assert 0

class process_cif_mixin(object):

  def process_cif_object(self, cif_object, file_name=None, cache=True, process_tor=False):
    if (file_name is None):
      source_info = None
    else:
      source_info = "file: "+file_name
    try:
      self.convert_all(source_info=source_info,
                       cif_object=cif_object,
                       cache=cache,
                       process_tor=process_tor,
        )
    except KeyboardInterrupt: raise
    except Exception:
      if (file_name is None): file_name = "(file name not available)"
      raise Sorry(
        "Error processing CIF file:\n"
        "  %s\n"
        "  (%s)" % (show_string(file_name), format_exception()))

  def process_cif(self, file_name, cache=True):
    try: cif_object = read_cif(file_name=file_name)
    except KeyboardInterrupt: raise
    except Exception:
      raise Sorry(
        "Error reading CIF file:\n"
        "  %s\n"
        "  (%s)" % (show_string(file_name), format_exception()))
    self.process_cif_object(cif_object=cif_object,
                            file_name=file_name,
                            cache=cache,
      )

class id_dict(dict):
  def __setitem__(self, key, item):
    print("DICT",key, item)
    if key in self:
      print('assert 0')
      assert 0
    dict.__setitem__(self, key, item)

class server(process_cif_mixin):

  def __init__(self, list_cif=None, another_list_cif=None, verbose=False):
    if (list_cif is None):
      list_cif = mon_lib_list_cif()
    if (another_list_cif is None):
      # get the geostd CIF links object and merge
      geostd_list_cif_obj = geostd_list_cif()
      list_cif = merge_and_overwrite_cifs(geostd_list_cif_obj,
                                          list_cif,
                                          )
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

  def convert_all(self,
                  source_info,
                  cif_object,
                  skip_comp_list=False,
                  cache=True,
                  process_tor=False,
                  ):
    if process_tor:
      self.convert_ccp4_tor_list(source_info=source_info, cif_object=cif_object)
    self.convert_deriv_list_dict(cif_object=cif_object)
    self.convert_comp_synonym_list(cif_object=cif_object)
    self.convert_comp_synonym_atom_list(cif_object=cif_object)
    if (not skip_comp_list):
      self.convert_comp_list(
        source_info=source_info, cif_object=cif_object, cache=cache)
    self.convert_link_list(
      source_info=source_info, cif_object=cif_object)
    self.convert_mod_list(
      source_info=source_info, cif_object=cif_object)

  def convert_deriv_list_dict(self, cif_object):
    for row in get_rows(cif_object,
                 "deriv_list", "_chem_comp_deriv"):
      deriv = cif_types.chem_comp_deriv(**dict(row))
      self.deriv_list_dict[deriv.comp_id] = deriv

  def convert_comp_synonym_list(self, cif_object):
    for row in get_rows(cif_object,
                 "comp_synonym_list", "_chem_comp_synonym"):
      self.comp_synonym_list_dict[row["_chem_comp_synonym.comp_alternative_id"]] \
          = row["_chem_comp_synonym.comp_id"]

  def convert_comp_synonym_atom_list(self, cif_object):
    for row in get_rows(cif_object,
                 "comp_synonym_atom_list", "_chem_comp_synonym_atom"):
      synonym = cif_types.chem_comp_synonym_atom(**dict(row))
      d = self.comp_synonym_atom_list_dict[synonym.comp_id]
      d[synonym.atom_alternative_id] = synonym.atom_id
      if (synonym.comp_alternative_id != ""):
        d = self.comp_synonym_atom_list_dict[synonym.comp_alternative_id]
        d[synonym.atom_alternative_id] = synonym.atom_id

  def convert_comp_list(self, source_info, cif_object, cache=True):
    for comp_comp_id in convert_comp_list(
                          source_info=source_info, cif_object=cif_object):
      comp_comp_id.normalize_atom_ids_in_place()
      chem_comp = comp_comp_id.chem_comp
      if cache:
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
          "mod_rna3p.cif",
          "mod_rna2p_pur.cif",
          "mod_rna3p_pur.cif",
          "mod_rna2p_pyr.cif",
          "mod_rna3p_pyr.cif"]:
      self.process_cif(
        file_name=os.path.join(self.geostd_path, "rna_dna", file_name))

  def convert_ccp4_tor_list(self, source_info, cif_object):
    def get_tor_key(rc):
      key = [rc['_chem_comp_tor.atom_id_1'],
             rc['_chem_comp_tor.atom_id_2'],
             rc['_chem_comp_tor.atom_id_3'],
             rc['_chem_comp_tor.atom_id_4'],
             ]
      key.sort()
      key=tuple(key)
      return key
    for key, block in cif_object.items():
      if key=='comp_list': continue
      tors = block.get_loop_or_row('_chem_comp_tor')
      if tors is None: continue
      if '_chem_comp_tor.alt_value_angle' in tors.keys(): continue
      alt_value_angle = {}
      remove=[]
      for i, rc in enumerate(tors.iterrows()):
        key = get_tor_key(rc)
        tmp = alt_value_angle.setdefault(key, [])
        if tmp:
          remove.append(i)
        if rc['_chem_comp_tor.value_angle'] not in tmp:
          tmp.append(rc['_chem_comp_tor.value_angle'])
      remove.reverse()
      for r in remove:
        tors.delete_row(r)
      column = []
      for row, key in zip(tors.iterrows(), alt_value_angle.keys()):
        assert key==get_tor_key(row)
        if len(alt_value_angle[key])==1:
          column.append('.')
        else:
          column.append(','.join(alt_value_angle[key][1:]))
      tors.add_column('_chem_comp_tor.alt_value_angle', column)

  def get_comp_comp_id_direct(self,
                              comp_id,
                              pH_range=None, # low *neutral high
                              specific_residue_restraints=None,
                              ad_hoc_single_atom_residues=False,
                              user_supplied_restraints_directory=None,
                              user_supplied_pre_post=None,
                              return_filename=False,
                             ):
    comp_id = comp_id.strip().upper()
    if (len(comp_id) == 0): return None
    #
    # user supplied file
    #
    result = None
    file_name = None
    if specific_residue_restraints:
      file_name = specific_residue_restraints
    elif pH_range is not None:
      if pH_range in ["neutral"]: rr = ""
      elif pH_range in ["low", "high"]:
        rr = "_pH_%s" % pH_range
      elif pH_range in ['neutron']:
        rr = "_%s" % pH_range
      else:
        assert 0, "unknown pH value : %s" % pH_range
      rr_cif_name = "data_%s%s.cif" % (comp_id.upper(), rr)
      rr_cif_name = os.path.join(self.geostd_path,
                                 comp_id[0].lower(),
                                 rr_cif_name)
      if os.path.exists(rr_cif_name):
        file_name = rr_cif_name
      else:
        return None
    else:
      result = self.comp_comp_id_dict.get(comp_id)
    if (result is not None):
      return result
    std_comp_id = self.comp_synonym_list_dict.get(comp_id, "").strip().upper()
    def find_user_file(user_path):
      for trial_comp_id in [std_comp_id, comp_id]:
        if not trial_comp_id: continue
        for dir_name in [user_path, os.path.join(user_path, trial_comp_id[0].lower())]:
          for cif_name in ["%s.cif" % (trial_comp_id),
                           "data_%s.cif" % (trial_comp_id)]:
            file_name=os.path.join(dir_name, cif_name)
            if (os.path.isfile(file_name)): return file_name
      return None
    def find_file():
      if user_supplied_restraints_directory and user_supplied_pre_post=='pre':
        file_name = find_user_file(user_supplied_restraints_directory)
        if file_name is not None: return file_name
      for i_pass in [0,1]:
        for trial_comp_id in [std_comp_id, comp_id]:
          if (len(trial_comp_id) == 0): continue
          dir_name = os.path.join(self.geostd_path, trial_comp_id[0].lower())
          # check the Geo Standard
          if (os.path.isdir(dir_name)):
            cif_name = "data_%s.cif" % (trial_comp_id)
            if (i_pass == 0):
              file_names = [os.path.join(dir_name, cif_name)]
              for file_name in file_names:
                if (os.path.isfile(file_name)):
                  return file_name
            else:
              cif_name = cif_name.lower()
              for node in os.listdir(dir_name):
                if (node.lower() != cif_name): continue
                print('node',node)
                print('i_pass'*10,i_pass)
                print(dir_name)
                print(os.listdir(dir_name))
                print('cif_name',cif_name)
                assert 0
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
      if user_supplied_restraints_directory:
        file_name = find_user_file(user_supplied_restraints_directory)
        if file_name is not None: return file_name
      return None
    if (file_name is None): file_name = find_file()
    if (file_name is None): return None
    cache=True
    if specific_residue_restraints:
      cache=False
    if return_filename: return file_name
    self.process_cif(file_name=file_name, cache=cache)
    if (len(std_comp_id) > 0):
      result = self.comp_comp_id_dict.get(std_comp_id)
    else:
      result = None
    if (result is not None):
      self.comp_comp_id_dict[comp_id] = result
    else:
      if specific_residue_restraints:
        cif_object = read_cif(file_name=specific_residue_restraints)
        for i_comp_comp, comp_comp_id in enumerate(convert_comp_list(
            source_info=specific_residue_restraints,
            cif_object=cif_object,
          )):
          comp_comp_id.normalize_atom_ids_in_place()
          result=comp_comp_id
        assert i_comp_comp==0
        return result
      else:
        result = self.comp_comp_id_dict.get(comp_id)
      if (result is not None):
        if (len(std_comp_id) != 0):
          self.comp_comp_id_dict[std_comp_id] = result
      elif ad_hoc_single_atom_residues:
        comp_id = comp_id.replace('_EL', '')
        result = self.comp_comp_id_dict.get(comp_id)
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
        translate_cns_dna_rna_residue_names=None,
        specific_residue_restraints=None,
        ad_hoc_single_atom_residues=False,
        user_supplied_restraints_directory=None,
        user_supplied_pre_post=None,
        ):
    # not sure this works with specific_residue_restraints
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
      self.get_comp_comp_id_direct(
        comp_id=rnpani.work_residue_name,
        specific_residue_restraints=specific_residue_restraints,
        ad_hoc_single_atom_residues=ad_hoc_single_atom_residues,
        user_supplied_restraints_directory=user_supplied_restraints_directory,
        user_supplied_pre_post=user_supplied_pre_post,
        ),
      rnpani.atom_name_interpretation,
      )

  def rotamer_iterator(self,
        comp_id,
        atom_names,
        sites_cart,
        fine_sampling=False):
    comp_comp_id = self.get_comp_comp_id_direct(comp_id=comp_id)
    if (comp_comp_id is None):
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

  def __init__(self,
               ener_lib_cif=None,
               source_info=None,
               use_neutron_distances=False,
               ):
    if (ener_lib_cif is None):
      ener_lib_cif = mon_lib_ener_lib_cif()
    self.lib_synonym = {}
    self.lib_atom = {}
    self.lib_vdw = []
    self.convert_all(source_info=source_info,
                     cif_object=ener_lib_cif.cif,
                     use_neutron_distances=use_neutron_distances,
                     )
    self.source_infos = []

  def convert_all(self,
                  source_info,
                  cif_object,
                  use_neutron_distances=False,
                  cache=True, # not used
                  process_tor=False, # not used
                  ):
    if (source_info is not None):
      self.source_infos.append(source_info)
    self.convert_lib_synonym(cif_object=cif_object)
    self.convert_lib_atom(cif_object=cif_object,
                          use_neutron_distances=use_neutron_distances,
                          )
    self.convert_lib_vdw(cif_object=cif_object)

  def convert_lib_synonym(self, cif_object):
    for row in get_rows(cif_object, "energy", "_lib_synonym"):
      syn = cif_types.energy_lib_synonym(**dict(row))
      self.lib_synonym[syn.atom_alternative_type] = syn.atom_type

  def convert_lib_atom(self, cif_object, use_neutron_distances=False):
    for row in get_rows(cif_object, "energy", "_lib_atom"):
      entry = cif_types.energy_lib_atom(**dict(row))
      if use_neutron_distances:
        if entry.vdw_radius_neutron is not None:
          entry.vdw_radius = entry.vdw_radius_neutron
      self.lib_atom[entry.type] = entry

  def convert_lib_vdw(self, cif_object):
    for row in get_rows(cif_object, "energy", "_lib_vdw"):
      vdw = cif_types.energy_lib_vdw(**dict(row))
      self.lib_vdw.append(vdw)
