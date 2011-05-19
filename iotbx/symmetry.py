
from libtbx import group_args
from libtbx.utils import Sorry
import sys

class manager (object) :
  """
  Class for keeping track of symmetry information from multiple files.  This
  is particularly problematic in the phenix.refine GUI, where users may supply
  any number of PDB files as input, plus a PDB file representing a reference
  structure, and up to five reflection files.  Automatically taking the latest
  symmetry provided, without taking into account the symmetry information in
  other files and the current GUI state, may result in errors if files contain
  incompatible information.
  """
  def __init__ (self, prefer_pdb_space_group=True) :
    self.pdb_file_names = []
    self.reflection_file_names = []
    self.symmetry_by_file = {}
    self.current_space_group = None
    self.current_unit_cell = None
    self.prefer_pdb_space_group = prefer_pdb_space_group

  def set_current (self, space_group, unit_cell) :
    self.current_space_group = space_group
    self.current_unit_cell = unit_cell

  def get_current (self) :
    return (self.current_space_group, self.current_unit_cell)

  def as_symmetry_object (self) :
    if (self.current_space_group is None) or (self.current_unit_cell is None) :
      raise Sorry("Either the space group or the unit cell (or both) is "+
        "undefined.")
    import cctbx.crystal
    return cctbx.crystal.symmetry(
      space_group_info=self.current_space_group,
      unit_cell=self.current_unit_cell)

  def get_current_as_strings (self) :
    sg, uc = self.get_current()
    if (uc is None) :
      uc_str = ""
    else :
      uc_str = "%.3g %.3g %.3g %.3g %.3g %.3g" % uc.parameters()
    if (sg is None) :
      sg_str = ""
    else :
      sg_str = str(sg)
    return (sg_str, uc_str)

  def set_current_as_strings (self, space_group, unit_cell) :
    """Set symmetry from fields in the GUI."""
    if (space_group == "") or (unit_cell is None) :
      self.current_space_group = None
    else :
      from cctbx import sgtbx
      try :
        self.current_space_group = sgtbx.space_group_info(space_group)
      except RuntimeError, e :
        if ("symbol not recognized" in str(e)) :
          raise Sorry(("The current value for the space group parameter, "+
            "'%s', could not be recognized as a valid space group symbol.") %
            space_group)
        else :
          raise
    if (unit_cell == "") or (unit_cell is None) :
      self.current_unit_cell = None
    else :
      from cctbx import uctbx
      self.current_unit_cell = uctbx.unit_cell(unit_cell)

  def process_pdb_file (self, input_file) :
    """Extract symmetry info from iotbx.file_reader._any_file object"""
    symm = input_file.file_object.crystal_symmetry()
    if (symm is not None) :
      space_group = symm.space_group_info()
      unit_cell = symm.unit_cell()
    else :
      space_group, unit_cell = None, None
    file_name = input_file.file_name
    return self.add_pdb_file(file_name, space_group, unit_cell)

  def add_pdb_file (self, file_name, space_group, unit_cell) :
    self.pdb_file_names.append(file_name)
    self.symmetry_by_file[file_name] = (space_group, unit_cell)
    return self.check_consistency_and_set_symmetry(
      file_name=file_name,
      space_group=space_group,
      unit_cell=unit_cell,
      file_type="pdb")

  def process_reflections_file (self, input_file) :
    """Extract symmetry info from iotbx.file_reader._any_file object"""
    symm = input_file.file_server.miller_arrays[0].crystal_symmetry()
    if (symm is not None) :
      space_group = symm.space_group_info()
      unit_cell = symm.unit_cell()
    else :
      space_group, unit_cell = None, None
    file_name = input_file.file_name
    return self.add_reflections_file(file_name, space_group, unit_cell)

  def add_reflections_file (self, file_name, space_group, unit_cell) :
    self.reflection_file_names.append(file_name)
    self.symmetry_by_file[file_name] = (space_group, unit_cell)
    return self.check_consistency_and_set_symmetry(
      file_name=file_name,
      space_group=space_group,
      unit_cell=unit_cell,
      file_type="hkl")

  def check_cell_compatibility (self, program_name,
      raise_error_if_incomplete=False) :
    if (self.current_unit_cell is None) or (self.current_space_group is None) :
      if (raise_error_if_incomplete) :
        raise Sorry("Either the unit cell or the space group (or both) is "+
          "not set; these parameters are required to run %s." % program_name)
      return None
    else :
      from cctbx import crystal
      try :
        symm = crystal.symmetry(space_group=self.current_space_group.group(),
          unit_cell=self.current_unit_cell)
      except AssertionError, e :
        raise Sorry("Unit cell parameters are not consistent with the "+
          "currently set space group.  Please make sure that the symmetry "+
          "information is entered correctly.")
      else :
        return True

  def check_consistency_and_set_symmetry (self, file_name, space_group,
      unit_cell, file_type) :
    space_group_mismatch = False
    set_new_space_group = False
    unit_cell_mismatch = False
    incompatible_cell = False
    if (space_group is not None) :
      if (self.current_space_group is not None) :
        current_sgname = str(self.current_space_group)
        new_sgname = str(space_group)
        if (current_sgname != new_sgname) :
          group = self.current_space_group.group()
          derived_sg = group.build_derived_point_group()
          if (space_group.group().build_derived_point_group() != derived_sg) :
            space_group_mismatch = True
          elif (file_type == "pdb") and (self.prefer_pdb_space_group) :
            self.current_space_group = space_group
      else :
        self.current_space_group = space_group
    if (unit_cell is not None) :
      if (self.current_unit_cell is not None) :
        if (not self.current_unit_cell.is_similar_to(unit_cell)) :
          unit_cell_mismatch = True
      else :
        self.current_unit_cell = unit_cell
    return (space_group_mismatch, unit_cell_mismatch)

  def get_symmetry_choices (self) :
    sg_files = []
    uc_files = []
    all_file_names = self.pdb_file_names + self.reflection_file_names
    for file_name in all_file_names :
      space_group, unit_cell = self.symmetry_by_file[file_name]
      if (space_group is not None) :
        sg_files.append((file_name, str(space_group)))
      if (unit_cell is not None) :
        uc_files.append((file_name, str(unit_cell)))
    return group_args(
      current_space_group=str(self.current_space_group),
      current_unit_cell=str(self.current_unit_cell),
      space_group_files=sg_files,
      unit_cell_files=uc_files)

  def show (self, out=None) :
    if (out is None) :
      out = sys.stdout
    all_file_names = self.pdb_file_names + self.reflection_file_names
    for file_name in all_file_names :
      space_group, unit_cell = self.symmetry_by_file[file_name]
      print >> out, "%s: %s %s" % (os.path.basename(file_name), str(unit_cell),
        str(space_group))
