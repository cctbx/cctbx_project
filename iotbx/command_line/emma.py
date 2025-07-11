"""Match atoms in two models, allowing for origin shifts
and crystallographic relationships"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.emma
# LIBTBX_SET_DISPATCHER_NAME iotbx.emma

from iotbx import crystal_symmetry_from_any
import iotbx.pdb
from iotbx.cns import sdb_reader
from iotbx.kriber import strudat
from iotbx.option_parser import option_parser
from cctbx import euclidean_model_matching as emma
import sys, os
import cctbx.xray
from six.moves import zip

class MultipleEntriesError(RuntimeError): pass

def get_emma_model_from_pdb(file_name=None,
                            pdb_records=None,
                            crystal_symmetry=None):
  assert [file_name, pdb_records].count(None) == 1
  if (pdb_records is None):
    pdb_inp = iotbx.pdb.input(file_name=file_name)
  else:
    pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_records)
  crystal_symmetry = pdb_inp.crystal_symmetry(
    crystal_symmetry=crystal_symmetry,
    weak_symmetry=True)
  if (not crystal_symmetry or crystal_symmetry.unit_cell() is None):
    raise RuntimeError("Unit cell parameters unknown for model %s." %(
        file_name))
  if (crystal_symmetry.space_group_info() is None):
    raise RuntimeError("Space group unknown.")
  positions = []
  for atom in pdb_inp.atoms_with_labels():
    positions.append(emma.position(
      ":".join([str(len(positions)+1),
                atom.name, atom.resname, atom.chain_id]),
      crystal_symmetry.unit_cell().fractionalize(atom.xyz)))
  assert len(positions) > 0
  result = emma.model(
    crystal_symmetry.special_position_settings(),
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
    force=True)
  positions = []
  for i,site in enumerate(sdb_file.sites):
    if (crystal_symmetry.unit_cell() is None):
      raise RuntimeError("Unit cell parameters unknown.")
    positions.append(emma.position(
      ":".join((str(i+1), site.segid, site.type)),
      crystal_symmetry.unit_cell().fractionalize((site.x, site.y, site.z))))
  assert len(positions) > 0
  result = emma.model(
    crystal_symmetry.special_position_settings(),
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
    crystal_symmetry.special_position_settings(),
    positions)
  result.label = file_name
  return result

def get_emma_model_from_ins(file_name):
  return cctbx.xray.structure.from_shelx(file=open(file_name)).as_emma_model()

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
  except KeyboardInterrupt: raise
  except Exception:
    if (iotbx.pdb.is_pdb_file(file_name)): raise
  try:
    return get_emma_model_from_sdb(
      file_name=file_name,
      crystal_symmetry=crystal_symmetry)
  except MultipleEntriesError:
    raise
  except KeyboardInterrupt: raise
  except Exception:
    pass
  try:
    return get_emma_model_from_solve(
      file_name=file_name,
      crystal_symmetry=crystal_symmetry)
  except KeyboardInterrupt: raise
  except Exception:
    pass
  try:
    return get_emma_model_from_ins(file_name=file_name)
  except KeyboardInterrupt: raise
  except Exception:
    pass
  try:
    return get_emma_model_from_strudat(file_name=file_name)
  except MultipleEntriesError:
    raise
  except KeyboardInterrupt: raise
  except Exception:
    pass
  raise RuntimeError("Coordinate file %s: unknown format." % file_name)

def run(args, command_name="phenix.emma"):
  command_line = (option_parser(
    usage=command_name + " [options]"
         +" reference_coordinates other_coordinates",
    description="Example: %s model1.pdb model2.sdb" % command_name)
    .enable_symmetry_comprehensive()
    .option(None, "--output_pdb",
      action="store",
      type="str",
      default="",
      help="Output pdb: second model transformed to best match first model",
      metavar="STR")
    .option(None, "--tolerance",
      action="store",
      type="float",
      default=3.,
      help="match tolerance",
      metavar="FLOAT")
    .option(None, "--diffraction_index_equivalent",
      action="store_true",
      help="Use only if models are diffraction-index equivalent.")
  ).process(args=args, nargs=2)
  crystal_symmetry = command_line.symmetry
  if (   crystal_symmetry.unit_cell() is None
      or crystal_symmetry.space_group_info() is None):
    for file_name in command_line.args:
      crystal_symmetry = crystal_symmetry.join_symmetry(
        other_symmetry=crystal_symmetry_from_any.extract_from(
          file_name=file_name),
        force=False)
  output_pdb = command_line.options.output_pdb
  if output_pdb:
    print("Output pdb:",output_pdb)
  tolerance = command_line.options.tolerance
  print("Tolerance:", tolerance)
  if (tolerance <= 0.):
    raise ValueError("Tolerance must be greater than zero.")
  print()
  diffraction_index_equivalent = \
    command_line.options.diffraction_index_equivalent
  if (diffraction_index_equivalent):
    print("Models are diffraction index equivalent.")
    print()
  second_model_as_pdb_inp=None
  emma_models = []
  for file_name in command_line.args:
    emma_models.append(get_emma_model(
      file_name=file_name,
      crystal_symmetry=crystal_symmetry))
    if len(emma_models)==2 and os.path.isfile(file_name):
      try:
        second_model_as_pdb_inp=iotbx.pdb.input(
           file_name=file_name)
      except Exception as e:
        pass
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
    print("No matches.")
    print()
  else:
    max_n_pairs = None
    first=True
    for match in model_matches.refined_matches:
      if (max_n_pairs is None or len(match.pairs) > max_n_pairs*0.2):
        print("." * 79)
        print()
        match.show()
        if first and output_pdb: # 2013-01-25 tt
          if second_model_as_pdb_inp:
            match.get_transformed_model2(output_pdb=output_pdb,
              template_pdb_inp=second_model_as_pdb_inp,
              f=sys.stdout)
          else:
            print("No output model as input model was not PDB")
        first=False
      if (max_n_pairs is None):
        max_n_pairs = len(match.pairs)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
