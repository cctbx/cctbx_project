from iotbx import crystal_symmetry_from_any
import iotbx.pdb
from iotbx.pdb import cryst1_interpretation
from iotbx.pdb import crystal_symmetry_from_pdb
from iotbx.cns import sdb_reader
from iotbx.shelx import from_ins
from iotbx.kriber import strudat
from iotbx.option_parser import iotbx_option_parser
from cctbx import euclidean_model_matching as emma
from cctbx import crystal
from libtbx.itertbx import count
import sys, os

class MultipleEntriesError(RuntimeError): pass

def get_emma_model_from_pdb(file_name=None,
                            pdb_records=None,
                            crystal_symmetry=None):
  assert [file_name, pdb_records].count(None) == 1
  if (pdb_records is None):
    pdb_records = iotbx.pdb.parser.collect_records(
      raw_records=open(file_name),
      ignore_master=0001)
  cryst1_symmetry = None
  for record in pdb_records:
    if (record.record_name.startswith("CRYST1")):
      try:
        cryst1_symmetry = cryst1_interpretation.crystal_symmetry(
          cryst1_record=record)
      except: pass
      break
  if (cryst1_symmetry is not None):
    crystal_symmetry = cryst1_symmetry.join_symmetry(
      other_symmetry=crystal_symmetry,
      force=00000)
  positions = []
  for record in pdb_records:
    if (not record.record_name in ("ATOM", "HETATM")): continue
    if (crystal_symmetry.unit_cell() is None):
      raise RuntimeError("Unit cell parameters unknown.")
    positions.append(emma.position(
      ":".join([str(len(positions)+1),
                record.name, record.resName, record.chainID]),
      crystal_symmetry.unit_cell().fractionalize(record.coordinates)))
  assert len(positions) > 0
  result = emma.model(
    crystal.special_position_settings(crystal_symmetry),
    positions)
  if (file_name is not None):
    result.label = file_name
  return result

def get_emma_model_from_sdb(file_name, crystal_symmetry):
  sdb_files = sdb_reader.multi_sdb_parser(open(file_name))
  if (len(sdb_files) > 1):
    raise MultipleEntriesError(
      "SDB file %s may contain only one structure." % file_name)
  assert len(sdb_files) == 1
  sdb_file = sdb_files[0]
  crystal_symmetry = crystal_symmetry.join_symmetry(
    other_symmetry=sdb_file.crystal_symmetry(),
    force=0001)
  positions = []
  for i,site in zip(count(1),sdb_file.sites):
    if (crystal_symmetry.unit_cell() is None):
      raise RuntimeError("Unit cell parameters unknown.")
    positions.append(emma.position(
      ":".join((str(i), site.segid, site.type)),
      crystal_symmetry.unit_cell().fractionalize((site.x, site.y, site.z))))
  assert len(positions) > 0
  result = emma.model(
    crystal.special_position_settings(crystal_symmetry),
    positions)
  result.label = sdb_file.file_name
  return result

def get_emma_model_from_solve(file_name, crystal_symmetry):
  positions = []
  for line in open(file_name):
    flds = line.split()
    if (len(flds) < 4 or flds[0].lower() != "xyz"): continue
    site = [float(x) for x in flds[1:4]]
    positions.append(emma.position("site"+str(len(positions)+1), site))
  assert len(positions) > 0
  result = emma.model(
    crystal.special_position_settings(crystal_symmetry),
    positions)
  result.label = file_name
  return result

def get_emma_model_from_ins(file_name):
  return from_ins.from_ins(file_name=file_name).as_emma_model()

def get_emma_model_from_strudat(file_name):
  strudat_entries = strudat.read_all_entries(open(file_name))
  if (len(strudat_entries.entries) > 1):
    raise MultipleEntriesError(
      "strudat file %s may contain only one structure." % file_name)
  assert len(strudat_entries.entries) == 1
  return strudat_entries.entries[0].as_xray_structure().as_emma_model()

def get_emma_model(file_name, crystal_symmetry):
  if (not os.path.isfile(file_name)):
    raise RuntimeError("File not found: %s" % file_name)
  try:
    return get_emma_model_from_pdb(
      file_name=file_name,
      crystal_symmetry=crystal_symmetry)
  except:
    pass
  try:
    return get_emma_model_from_sdb(
      file_name=file_name,
      crystal_symmetry=crystal_symmetry)
  except MultipleEntriesError:
    raise
  except:
    pass
  try:
    return get_emma_model_from_solve(
      file_name=file_name,
      crystal_symmetry=crystal_symmetry)
  except:
    pass
  try:
    return get_emma_model_from_ins(file_name=file_name)
  except:
    pass
  try:
    return get_emma_model_from_strudat(file_name=file_name)
  except MultipleEntriesError:
    raise
  except:
    pass
  raise RuntimeError("Coordinate file %s: unknown format." % file_name)

def run(args):
  command_line = (iotbx_option_parser(
    usage="iotbx.emma [options] reference_coordinates"
         +" reference_coordinates other_coordinates",
    description="Example: iotbx.emma model1.pdb model2.sdb")
    .enable_symmetry_comprehensive()
    .option(None, "--tolerance",
      action="store",
      type="float",
      dest="tolerance",
      default=3.,
      help="match tolerance",
      metavar="FLOAT")
    .option(None, "--diffraction_index_equivalent",
      action="store_true",
      dest="diffraction_index_equivalent",
      help="Use only if models are diffraction-index equivalent.")
  ).process(args=args, nargs=2)
  crystal_symmetry = command_line.symmetry
  if (   crystal_symmetry.unit_cell() is None
      or crystal_symmetry.space_group_info() is None):
    for file_name in command_line.args:
      crystal_symmetry = crystal_symmetry.join_symmetry(
        other_symmetry=crystal_symmetry_from_any.extract_from(
          file_name=file_name),
        force=00000)
  tolerance = command_line.options.tolerance
  print "Tolerance:", tolerance
  if (tolerance <= 0.):
    raise ValueError, "Tolerance must be greater than zero."
  print
  diffraction_index_equivalent = \
    command_line.options.diffraction_index_equivalent
  if (diffraction_index_equivalent):
    print "Models are diffraction index equivalent."
    print
  emma_models = []
  for file_name in command_line.args:
    emma_models.append(get_emma_model(
      file_name=file_name,
      crystal_symmetry=crystal_symmetry))
  emma_models[0].show("Reference model")
  emma_models[1].show("Other model")
  for model,label in zip(emma_models, ["reference", "other"]):
    if (model.unit_cell() is None):
      raise RuntimeError("Unit cell parameters unknown (%s model)." % label)
    if (model.space_group_info() is None):
      raise RuntimeError("Space group unknown (%s model)." % label)
  model_matches = emma.model_matches(
    model1=emma_models[0],
    model2=emma_models[1],
    tolerance=tolerance,
    models_are_diffraction_index_equivalent=diffraction_index_equivalent)
  if (model_matches.n_matches() == 0):
    print "No matches."
    print
  else:
    max_n_pairs = None
    for match in model_matches.refined_matches:
      if (max_n_pairs is None or len(match.pairs) > max_n_pairs*0.2):
        print "." * 79
        print
        match.show()
      if (max_n_pairs is None):
        max_n_pairs = len(match.pairs)

if (__name__ == "__main__"):
  run(sys.argv[1:])
