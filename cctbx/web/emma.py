from cctbx import euclidean_model_matching as emma
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.web import io_utils
from cctbx.web import cgi_utils
from itertools import count
import sys

def interpret_form_data(form):
  inp = cgi_utils.inp_from_form(form,
    (("ucparams_1", ""),
     ("sgsymbol_1", ""),
     ("convention_1", ""),
     ("format_1", None),
     ("coor_type_1", None),
     ("skip_columns_1", 0),
     ("ucparams_2", ""),
     ("sgsymbol_2", ""),
     ("convention_2", ""),
     ("format_2", None),
     ("coor_type_2", None),
     ("skip_columns_2", 0),
     ("tolerance", "3.0"),
     ("diffraction_index_equivalent", None)))
  inp.coordinates = []
  for suffix in ("1", "2"):
    inp.coordinates.append(cgi_utils.coordinates_from_form(form, suffix))
  return inp

def interpret_generic_coordinate_line(line, skip_columns):
  flds = line.replace(",", " ").split()
  try: site = [float(x) for x in flds[skip_columns: skip_columns+3]]
  except Exception: raise RuntimeError, "FormatError: " + line
  if (len(site) != 3):
    raise RuntimeError, "FormatError: " + line
  return " ".join(flds[:skip_columns]), site

def pdb_file_to_emma_model(crystal_symmetry, pdb_inp, other_symmetry):
  crystal_symmetry = pdb_inp.crystal_symmetry(
    crystal_symmetry=crystal_symmetry,
    weak_symmetry=False)
  if (other_symmetry is not None):
    crystal_symmetry = crystal_symmetry.join_symmetry(
      other_symmetry=other_symmetry,
      force=False)
  positions = []
  for atom in pdb_inp.atoms_with_labels():
    if (crystal_symmetry.unit_cell() is None):
      raise RuntimeError("Unit cell parameters unknown.")
    positions.append(emma.position(
      ":".join([str(len(positions)+1),
                atom.name, atom.resname, atom.chain_id]),
      crystal_symmetry.unit_cell().fractionalize(atom.xyz)))
  m = emma.model(
    crystal_symmetry.special_position_settings(),
    positions)
  m.label = "Other model"
  return m

def sdb_file_to_emma_model(crystal_symmetry, sdb_file):
  positions = []
  for i,site in zip(count(1),sdb_file.sites):
    if (crystal_symmetry.unit_cell() is None):
      raise RuntimeError("Unit cell parameters unknown.")
    positions.append(emma.position(
      ":".join((str(i), site.segid, site.type)),
      crystal_symmetry.unit_cell().fractionalize((site.x, site.y, site.z))))
  m = emma.model(
    crystal_symmetry.special_position_settings(),
    positions)
  m.label = sdb_file.file_name
  return m

class web_to_models(object):

  def __init__(self, ucparams, sgsymbol, convention,
               format, coor_type, skip_columns, coordinates,
               other_symmetry=None):
    self.ucparams = ucparams
    self.sgsymbol = sgsymbol
    self.convention = convention
    self.other_symmetry = other_symmetry
    self.coordinate_format = None
    if (format == "generic"):
      skip_columns = io_utils.interpret_skip_columns(skip_columns)
      self.positions = []
      for line in coordinates:
        label, site = interpret_generic_coordinate_line(line, skip_columns)
        if (label == ""): label = "Site" + str(len(self.positions)+1)
        if (coor_type != "Fractional"):
          site = self.get_unit_cell().fractionalize(site)
        self.positions.append(emma.position(label, site))
      self.coordinate_format = "generic"
    else:
      if (self.coordinate_format is None):
        try:
          import iotbx.pdb
          pdb_inp = iotbx.pdb.input(source_info=None, lines=coordinates)
        except KeyboardInterrupt: raise
        except Exception:
          pass
        else:
          self.pdb_model = pdb_file_to_emma_model(
            self.crystal_symmetry(), pdb_inp, other_symmetry)
          if (len(self.pdb_model.positions()) > 0):
            self.coordinate_format = "pdb"
      if (self.coordinate_format is None):
        try:
          from iotbx.cns import sdb_reader
          self.sdb_files = sdb_reader.multi_sdb_parser(coordinates)
        except KeyboardInterrupt: raise
        except Exception:
          pass
        else:
          self.coordinate_format = "sdb"
      if (self.coordinate_format is None):
        raise RuntimeError("Coordinate format unknown.")
    self.i_next_model = 0

  def get_unit_cell(self, other_unit_cell=None):
    if (self.ucparams): return uctbx.unit_cell(self.ucparams)
    assert other_unit_cell is not None, "Unit cell parameters unknown."
    return other_unit_cell

  def get_space_group_info(self, other_space_group_info=None):
    if (self.sgsymbol):
      return sgtbx.space_group_info(
        symbol=self.sgsymbol,
        table_id=self.convention)
    assert other_space_group_info is not None, "Space group symbol unknown."
    return other_space_group_info

  def crystal_symmetry(self):
    unit_cell = None
    space_group_symbol = None
    if (self.ucparams): unit_cell = uctbx.unit_cell(self.ucparams)
    if (self.sgsymbol): space_group_symbol = self.sgsymbol
    return crystal.symmetry(
      unit_cell=unit_cell,
      space_group_symbol=space_group_symbol)

  def get_next(self):
    if (self.coordinate_format == "generic"):
      if (self.i_next_model): return None
      crystal_symmetry = crystal.symmetry(
        unit_cell=self.get_unit_cell(),
        space_group_info=self.get_space_group_info())
      m = emma.model(
        crystal.special_position_settings(crystal_symmetry),
        self.positions)
      m.label = "Model 2"
      self.i_next_model += 1
      return m
    if (self.coordinate_format == "pdb"):
      if (self.i_next_model): return None
      self.i_next_model += 1
      return self.pdb_model
    if (self.coordinate_format == "sdb"):
      if (self.i_next_model >= len(self.sdb_files)): return None
      sdb = self.sdb_files[self.i_next_model]
      crystal_symmetry = crystal.symmetry(
        unit_cell=self.get_unit_cell(sdb.unit_cell),
        space_group_info=self.get_space_group_info(sdb.space_group_info))
      if (self.other_symmetry is not None):
        crystal_symmetry = crystal_symmetry.join_symmetry(
          other_symmetry=self.other_symmetry,
          force=False)
      m = sdb_file_to_emma_model(crystal_symmetry, sdb)
      self.i_next_model += 1
      return m
    raise RuntimeError("Internal error.")

def model_is_too_large(model, max_number_of_positions=200):
  if (len(model.positions()) <= max_number_of_positions):
    return False
  print "*"*79
  print "ERROR: too many sites in this model (given %d, limit is %d)" % (
    len(model.positions()), max_number_of_positions)
  print "Note: EMMA is designed to work with heavy-atom substructures,"
  print "      NOT complete macromolecular structures."
  print "Hint: install cctbx on your computer and use the command"
  print "    phenix.emma"
  print "which has no limit on the number of sites."
  print "cctbx downloads:"
  print "    http://cci.lbl.gov/cctbx_build/"
  print "*"*79
  print
  return True

def run_implementation(server_info, inp, status):
  if (inp.ucparams_2 == ""):
      inp.ucparams_2 = inp.ucparams_1
  if (inp.sgsymbol_2 == ""):
      inp.sgsymbol_2 = inp.sgsymbol_1
      inp.convention_2 = inp.convention_1

  tolerance = float(inp.tolerance)
  print "Tolerance:", tolerance
  if (tolerance <= 0.):
    raise ValueError, "Tolerance must be greater than zero."
  print

  diffraction_index_equivalent = int(inp.diffraction_index_equivalent)
  if (diffraction_index_equivalent):
    print "Models are diffraction index equivalent."
    print

  models1 = web_to_models(
    inp.ucparams_1,
    inp.sgsymbol_1, inp.convention_1,
    inp.format_1, inp.coor_type_1, inp.skip_columns_1,
    inp.coordinates[0])
  model1 = models1.get_next()
  assert model1, "Problems reading reference model."
  model1.show("Reference model")
  if (model_is_too_large(model=model1)):
    return
  assert not models1.get_next()
  if (model1.unit_cell() is None):
    raise RuntimeError("Unit cell parameters unknown (reference model).")
  if (model1.space_group_info() is None):
    raise RuntimeError("Space group unknown (reference model).")

  models2 = web_to_models(
    inp.ucparams_2,
    inp.sgsymbol_2, inp.convention_2,
    inp.format_2, inp.coor_type_2, inp.skip_columns_2,
    inp.coordinates[1],
    other_symmetry=model1)
  while 1:
    print "#" * 79
    print
    model2 = models2.get_next()
    if (not model2): break
    model2.show(model2.label)
    assert model2.unit_cell() is not None
    assert model2.space_group_info() is not None
    if (model_is_too_large(model=model2)):
      continue
    sys.stdout.flush()
    model_matches = emma.model_matches(
      model1=model1,
      model2=model2,
      tolerance=tolerance,
      models_are_diffraction_index_equivalent=diffraction_index_equivalent)
    if (model_matches.n_matches() == 0):
      print "No matches."
      print
    else:
      for match in model_matches.refined_matches:
        print "." * 79
        print
        match.show()

def run(server_info, inp, status):
  print "<pre>"
  run_implementation(server_info=server_info, inp=inp, status=status)
  print "</pre>"
