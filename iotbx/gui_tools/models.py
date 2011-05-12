
# TODO: TESTS!, manage CIF files better, detect duplicates

import iotbx.gui_tools
from iotbx import file_reader
import cStringIO
import os

class cif_handler (iotbx.gui_tools.manager) :
  file_type = "cif"
  file_type_label = "CIF"

class model_handler (iotbx.gui_tools.manager) :
  file_type = "pdb"
  file_type_label = "PDB"
  def __init__ (self,
                allowed_param_names=None,
                allowed_multiple_params=None,
                debug=False,
                cif_param_names=None,
                multiple_cif_params=None,
                tmp_dir=None,
                construct_hierarchy=True) :
    iotbx.gui_tools.manager.__init__(self,
      allowed_param_names=allowed_param_names,
      allowed_multiple_params=allowed_multiple_params,
      debug=debug)
    self.cif_handler = cif_handler(
      allowed_param_names=cif_param_names,
      allowed_multiple_params=multiple_cif_params)
    self.tmp_dir = tmp_dir
    self.construct_hierarchy = construct_hierarchy
    self.add_callback = lambda file_name : True
    self.remove_callback = lambda file_name : True

  def clear_format_specific_cache (self) :
    if hasattr(self, "cif_handler") :
      self.cif_handler.clear_cache()
    self._cached_pdb_hierarchies = {}
    self._cached_bonds = {}
    self._cached_xray_structures = {}

  def set_callbacks (self, add_callback, remove_callback) :
    self.add_callback = add_callback
    self.remove_callback = remove_callback

  def add_file_callback (self, file_name) :
    self.add_callback(file_name)

  def remove_file_callback (self, file_name) :
    self.remove_callback(file_name)

  def save_other_file_data (self, input_file) :
    if self.construct_hierarchy :
      file_name = input_file.file_name
      pdb_hierarchy = input_file.file_object.construct_hierarchy()
      pdb_hierarchy.atoms().reset_i_seq()
      self._cached_pdb_hierarchies[file_name] = pdb_hierarchy
      return pdb_hierarchy

  def get_complete_model_file (self, file_param_name=None) :
    file_names = []
    if (file_param_name is not None) :
      file_names = self.get_param_files(file_param_name)
    else :
      file_names = self._cached_input_files.keys()
    if (len(file_names) == 1) :
      return file_names[0]
    elif (len(file_names) != 0) :
      pdb_hierarchies = []
      for file_name in file_names :
        pdb_hierarchies.append(self.get_pdb_hierarchy(file_name))
      tmp_dir = self.tmp_dir
      if tmp_dir is None :
        tmp_dir = "/var/tmp"
      assert os.path.isdir(tmp_dir)
      file_name = os.path.join(tmp_dir, "current_model.pdb")
      f = open(file_name, "w")
      for pdb_hierarchy in pdb_hierarchies :
        f.write(pdb_hierarchy.as_pdb_string())
      f.write("END")
      f.close()
      return file_name
    return None

  def combine_pdb_files (self, file_names) :
    symm = None
    pdb_str = cStringIO.StringIO()
    hierarchies = []
    for file_name in file_names :
      pdb_file = self._cached_input_files[file_name]
      file_symm = pdb_file.file_object.crystal_symmetry()
      if (file_symm is not None) and (symm is not None) :
        symm = file_symm
      hierarchy_str = self.get_pdb_hierarchy(file_name).as_pdb_string()
      hierarchies.append(hierarchy_str)
    if (symm is not None) :
      import iotbx.pdb
      cryst1 = iotbx.pdb.format_cryst1_record(symm)
      scale = iotbx.pdb.format_scale_records(symm.unit_cell())
      pdb_str.write(cryst1 + "\n")
      pdb_str.write(scale + "\n")
    pdb_str.write("\n".join(hierarchies))
    return pdb_str.getvalue()

  def get_combined_pdb_input (self, file_param_name=None) :
    if (file_param_name is not None) :
      file_names = self.get_param_files(file_param_name)
      if (len(file_names) == 1) :
        hierarchy = self.get_pdb_hierarchy(file_names[0])
        xray_structure = self.get_xray_structure(file_names[0])
        return (hierarchy, xray_structure)
    else :
      file_names = self._cached_input_files.keys()
    if (len(file_names) == 0) :
      raise RuntimeError("No PDB files loaded.")
    pdb_str = self.combine_pdb_files(file_names)
    import iotbx.pdb
    pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str)
    hierarchy = pdb_in.construct_hierarchy()
    hierarchy.atoms().reset_i_seq()
    xray_structure = pdb_in.xray_structure_simple()
    return (hierarchy, xray_structure)

  #--- CIF files
  def save_cif_file (self, *args, **kwds) :
    self.cif_handler.save_file(*args, **kwds)

  def set_param_cif_file (self, *args, **kwds) :
    self.cif_handler.set_param_file(*args, **kwds)

  def get_current_cif_file_names (self) :
    return self.cif_handler.get_current_file_names()

  def get_cifs (self) :
    files = []
    for file_name in self.cif_handler.get_current_file_names() :
      files.append(self.get_cif_file(file_name))
    return files

  def get_cif_objects (self) :
    return [ file.file_object for file in self.get_cifs() ]

  def get_cif_file (self, file_name) :
    return self.cif_handler.get_file(file_name)

  def remove_cif_file (self, file_name) :
    return self.cif_handler.remove_file(file_name)
  #---

  def get_pdb_hierarchy (self, file_name) :
    input_file = self.get_file(file_name)
    if (input_file is not None) :
      if (not file_name in self._cached_pdb_hierarchies) :
        pdb_hierarchy = input_file.file_object.construct_hierarchy()
        pdb_hierarchy.atoms().reset_i_seq()
        self._cached_pdb_hierarchies[file_name] = pdb_hierarchy
      return self._cached_pdb_hierarchies.get(file_name)
    return None

  def get_xray_structure (self, file_name) :
    pdb_file = self.get_file(file_name)
    if pdb_file is not None :
      xray_structure = self._cached_xray_structures.get(file_name, None)
      if xray_structure is None :
        xray_structure = pdb_file.file_object.xray_structure_simple()
        self._cached_xray_structures[file_name] = xray_structure
      return xray_structure
    return None

  def get_connectivity (self, file_name) :
    assert os.path.isfile(file_name)
    if (self.get_file(file_name) is None) :
      self.save_file(file_name=file_name)
    if ((file_name in self._cached_bonds) and
        (not self.file_is_modified(file_name))) :
      return self._cached_bonds[file_name]
    pdb_hierarchy = self.get_pdb_hierarchy(file_name)
    assert (pdb_hierarchy is not None)
    atomic_bonds = pdb_hierarchy.distance_based_simple_two_way_bond_sets()
    self._cached_bonds[file_name] = atomic_bonds
    return atomic_bonds

  def get_pdb_file_symmetry (self, file_name) :
    pdb_file = self.get_file(file_name)
    if (pdb_file is None) :
      pdb_file = file_reader.any_file(file_name)
      pdb_file.assert_file_type("pdb")
    return pdb_file.file_object.crystal_symmetry()

  def create_copy_with_fake_symmetry (self, file_name, tmp_dir=None) :
    import iotbx.pdb
    if (tmp_dir is None) :
      if (self.tmp_dir is None) :
        tmp_dir = "/var/tmp"
      else :
        tmp_dir = self.tmp_dir
    assert os.path.isdir(tmp_dir)
    pdb_hierarchy = self.get_pdb_hierarchy(file_name)
    if pdb_hierarchy is None :
      pdb_file = file_reader.any_file(file_name)
      pdb_file.assert_file_type("pdb")
      self.save_file(pdb_file)
      pdb_hierarchy = self._cached_pdb_hierarchies[file_name]
    xyz = pdb_hierarchy.atoms().extract_xyz()
    symm = get_fake_symmetry(xyz.min(), xyz.max())
    output_file = os.path.join(tmp_dir, os.path.basename(file_name))
    f = open(output_file, "w")
    f.write("%s\n" % iotbx.pdb.format_cryst1_and_scale_records(
                        crystal_symmetry=symm,
                        write_scale_records=True))
    f.write(pdb_hierarchy.as_pdb_string())
    f.close()
    return output_file

def get_fake_symmetry (xyz_min, xyz_max) :
  from iotbx import crystal_symmetry_from_any
  a = xyz_max[0] - xyz_min[0] + 10.0
  b = xyz_max[1] - xyz_min[1] + 10.0
  c = xyz_max[2] - xyz_min[2] + 10.0
  combined = "%.3f,%.3f,%.3f,90,90,90,P1" % (a, b, c)
  symm = crystal_symmetry_from_any.from_string(combined)
  return symm
