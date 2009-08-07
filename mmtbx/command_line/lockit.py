from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import iotbx.pdb.atom_name_interpretation
import iotbx.mtz
import iotbx.phil, libtbx.phil
from cctbx import maptbx
import cctbx.maptbx.real_space_refinement_simple
from cctbx.array_family import flex
import scitbx.lbfgs
from libtbx.str_utils import show_string
from libtbx.utils import Sorry
import libtbx
import sys, os
op = os.path

def residue_is_suitable_for_scoring(residue):
  matched_atom_names = iotbx.pdb.atom_name_interpretation.interpreters[
    residue.resname].match_atom_names(
      atom_names=residue.atoms().extract_name())
  names = matched_atom_names.unexpected
  if (len(names) != 0):
    print residue.id_str(), "unexpected atoms:", " ".join(sorted(names))
    return False
  names = matched_atom_names.missing_atom_names(ignore_hydrogen=True)
  if (len(names) != 0):
    print residue.id_str(), "missing atoms:", " ".join(sorted(names))
    return False
  return True

def score_residue(
      density_map,
      pdb_hierarchy,
      geometry_restraints_manager,
      sites_cart,
      residue):
  rs_f = maptbx.real_space_target_simple(
    unit_cell=geometry_restraints_manager.crystal_symmetry.unit_cell(),
    density_map=density_map,
    sites_cart=sites_cart.select(indices=residue.atoms().extract_i_seq()))
  print residue.id_str(), "score: %.6g" % rs_f

def rotamer_scoring(
      density_map,
      pdb_hierarchy,
      geometry_restraints_manager,
      sites_cart):
  n_other_residues = 0
  n_amino_acids_ignored = 0
  n_amino_acids_scored = 0
  get_class = iotbx.pdb.common_residue_names_get_class
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue in chain.only_conformer().residues():
        if (get_class(residue.resname) != "common_amino_acid"):
          n_other_residues += 1
        elif (not residue_is_suitable_for_scoring(residue=residue)):
          n_amino_acids_ignored += 1
        else:
          score_residue(
            density_map=density_map,
            pdb_hierarchy=pdb_hierarchy,
            geometry_restraints_manager=geometry_restraints_manager,
            sites_cart=sites_cart,
            residue=residue)
          n_amino_acids_scored += 1
  if (n_amino_acids_scored != 0): print
  print "number of amino acid residues scored:", n_amino_acids_scored
  print "number of amino acid residues ignored:", n_amino_acids_ignored
  print "number of other residues:", n_other_residues
  print
  sys.stdout.flush()

class try_read_file(object):

  def __init__(O, file_name):
    def set(file_type, file_content):
      O.file_name = file_name
      O.file_type = file_type
      O.file_content = file_content
    lead = open(file_name, "rb").read(3)
    if (lead == "MTZ"):
      mtz_obj = iotbx.mtz.object(file_name=file_name)
      set(file_type="mtz", file_content=mtz_obj)
      return
    try:
      pdb_inp = iotbx.pdb.input(file_name=file_name)
    except KeyboardInterrupt: raise
    except:
      if (iotbx.pdb.is_pdb_file(file_name=file_name)):
        raise
      pdb_inp = None
    else:
      if (pdb_inp.atoms().size() != 0):
        set(file_type="pdb", file_content=pdb_inp)
        return
    try:
      cif_obj = mmtbx.monomer_library.server.read_cif(file_name=file_name)
    except KeyboardInterrupt: raise
    except: pass
    else:
      if (len(cif_obj) != 0):
        set(file_type="cif", file_content=cif_obj)
        return
    try:
      phil_obj = iotbx.phil.parse(file_name=file_name)
    except KeyboardInterrupt: raise
    except: pass
    else:
      set(file_type="phil", file_content=phil_obj)
      return
    if (pdb_inp is not None):
      if (pdb_inp.unknown_section().size() != 0):
        set(file_type=None, file_content=None)
        return
      if (pdb_inp.header_section().size() != 0):
        set(file_type="pdb", file_content=pdb_inp)
    set(file_type=None, file_content=None)
    return

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""\
symmetry_from_file = None
  .type = path
  .multiple = True
unit_cell = None
  .type=unit_cell
space_group = None
  .type=space_group

map_coeff_labels = 2FOFCWT PH2FOFCWT
  .type = strings
map_resolution_factor = 1/3
  .type = float

lbfgs_max_iterations = 0
  .type = int
real_space_target_weight = 1
  .type = float
real_space_gradients_delta_resolution_factor = 1/3
  .type = float

rotamer_scoring = False
  .type = bool

pdb_interpretation {
  include scope mmtbx.monomer_library.pdb_interpretation.master_params
}

geometry_restraints.edits {
  include scope \
    mmtbx.monomer_library.pdb_interpretation.geometry_restraints_edits_str
}

geometry_restraints.remove {
  include scope \
    mmtbx.monomer_library.pdb_interpretation.geometry_restraints_remove_str
}
""", process_includes=True)

def run(args):
  show_times = libtbx.utils.show_times(time_start="now")
  master_phil = get_master_phil()
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil)
  phil_objects = []
  file_objects = {
    "mtz": [],
    "pdb": [],
    "cif": []}
  for arg in args:
    if (len(arg) == 0): continue
    def try_as_file():
      if (not op.isfile(arg)): return False
      obj = try_read_file(file_name=arg)
      if (obj.file_type is None): return False
      if (obj.file_type == "phil"):
        phil_objects.append(obj.file_content)
      else:
        file_objects[obj.file_type].append(obj)
      return True
    def try_as_command_line_params():
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except:
        if (op.isfile(arg)):
          raise Sorry(
            "Error processing file: %s" % show_string(arg))
        raise Sorry(
          "Command-line argument not recognized: %s" % show_string(arg))
      phil_objects.append(command_line_params)
    if (not try_as_file()):
      try_as_command_line_params()
  work_phil = master_phil.fetch(sources=phil_objects)
  work_phil.show()
  print
  work_params = work_phil.extract()
  #
  assert len(work_params.symmetry_from_file) == 0 # TODO not implemented
  assert work_params.unit_cell is None # TODO not implemented
  assert work_params.space_group is None # TODO not implemented
  #
  assert len(file_objects["mtz"]) == 1
  miller_arrays = file_objects["mtz"][0].file_content.as_miller_arrays()
  map_coeffs = None
  for miller_array in miller_arrays:
    if (miller_array.info().labels == work_params.map_coeff_labels):
      map_coeffs = miller_array
      break
  assert map_coeffs is not None
  #
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  for file_obj in file_objects["cif"]:
    print "Processing CIF file: %s" % show_string(file_obj.file_name)
    for srv in [mon_lib_srv, ener_lib]:
      srv.process_cif_object(
        cif_object=file_obj.file_content,
        file_name=file_obj.file_name)
  #
  assert len(file_objects["pdb"]) == 1 # TODO not implemented
  file_obj = file_objects["pdb"][0]
  processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=work_params.pdb_interpretation,
    file_name=file_obj.file_name,
    pdb_inp=file_obj.file_content,
    strict_conflict_handling=True,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    log=sys.stdout)
  #
  grm = processed_pdb_file.geometry_restraints_manager(
    params_edits=work_params.geometry_restraints.edits,
    params_remove=work_params.geometry_restraints.remove)
  print
  sys.stdout.flush()
  #
  d_min = map_coeffs.d_min()
  fft_map = map_coeffs.fft_map(
    d_min=d_min,
    resolution_factor=work_params.map_resolution_factor)
  fft_map.apply_sigma_scaling()
  density_map = fft_map.real_map()
  real_space_gradients_delta = \
    d_min * work_params.real_space_gradients_delta_resolution_factor
  print "real_space_gradients_delta: %.6g" % real_space_gradients_delta
  print
  sys.stdout.flush()
  #
  sites_cart = processed_pdb_file.all_chain_proxies.sites_cart_exact()
  #
  if (work_params.lbfgs_max_iterations != 0):
    refined = maptbx.real_space_refinement_simple.lbfgs(
      sites_cart=sites_cart,
      density_map=density_map,
      geometry_restraints_manager=grm,
      real_space_target_weight=work_params.real_space_target_weight,
      real_space_gradients_delta=real_space_gradients_delta,
      lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=work_params.lbfgs_max_iterations))
    sites_cart = refined.sites_cart
    grm.energies_sites(sites_cart=sites_cart).show()
    print
    print "number_of_function_evaluations:", \
      refined.number_of_function_evaluations
    print "real+geo target start: %.6g" % refined.f_start
    print "real+geo target final: %.6g" % refined.f_final
    print
  if (work_params.rotamer_scoring):
    rotamer_scoring(
      density_map=density_map,
      pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
      geometry_restraints_manager=grm,
      sites_cart=sites_cart)
  show_times()
  sys.stdout.flush()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
