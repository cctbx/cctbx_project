from iotbx import pdb
import iotbx.phil
from iotbx.option_parser import option_parser
from cctbx import euclidean_model_matching
from cctbx import sgtbx
from scitbx import matrix
from libtbx.utils import Sorry
from libtbx.str_utils import show_string
from libtbx.test_utils import approx_equal
import sys, os

master_params = iotbx.phil.parse(input_string="""\
reference {
  file_name=None
    .type=path
  atom_selection=None
    .type=str
}
other {
  file_name=None
    .type=path
  atom_selection=None
    .type=str
}
output {
  file_name=None
    .type=path
  atom_selection=None
    .type=str
}
crystal_symmetry {
  unit_cell=None
    .type=unit_cell
  space_group=None
    .type=space_group
}
""")

def run(args, command_name="iotbx.pdb.superpose_centers_of_mass"):
  if (len(args) == 0): args = ["--help"]
  command_line = (option_parser(
    usage=
      "%s [options] [reference_file] [other_file] [parameter_file]" %
       command_name)
    .enable_show_defaults()
    .enable_symmetry_comprehensive()
  ).process(args=args)
  if (command_line.expert_level is not None):
    master_params.show(
      expert_level=command_line.expert_level,
      attributes_level=command_line.attributes_level)
    sys.exit(0)
  #
  # Loop over command-line arguments.
  #
  parameter_interpreter = master_params.command_line_argument_interpreter()
  parsed_params = []
  pdb_file_names = []
  command_line_params = []
  for arg in command_line.args:
    arg_is_processed = False
    if (os.path.isfile(arg)):
      params = None
      try: params = iotbx.phil.parse(file_name=arg)
      except KeyboardInterrupt: raise
      except RuntimeError: pass
      else:
        if (len(params.objects) == 0):
          params = None
      if (params is not None):
        parsed_params.append(params)
        arg_is_processed = True
      elif (pdb.is_pdb_file(file_name=arg)):
        pdb_file_names.append(arg)
        arg_is_processed = True
    if (not arg_is_processed):
      try:
        params = parameter_interpreter.process(arg=arg)
      except Sorry, e:
        if (not os.path.isfile(arg)): raise
        raise Sorry("Unknown file format: %s" % arg)
      else:
        command_line_params.append(params)
  #
  # Consolidation of inputs, resulting in effective phil_params.
  #
  phil_params = master_params.fetch(
    sources=parsed_params+command_line_params)
  params = phil_params.extract()
  for param_group in [params.reference, params.other, params.output]:
    if (param_group.file_name is None
        and len(pdb_file_names) > 0):
      param_group.file_name = pdb_file_names[0]
      pdb_file_names = pdb_file_names[1:]
  if (len(pdb_file_names) > 0):
    raise Sorry("Too many PDB file names: %s" % ", ".join([
      show_string(s) for s in pdb_file_names]))
  if (params.output.file_name is None
      and params.other.file_name is not None):
    name = os.path.basename(params.other.file_name)
    if (name.lower().endswith(".pdb")): name = name[:-4]
    name += "_superposed.pdb"
    params.output.file_name = name
  if (params.crystal_symmetry.unit_cell is None):
    params.crystal_symmetry.unit_cell = \
      command_line.symmetry.unit_cell()
  if (params.crystal_symmetry.space_group is None):
    params.crystal_symmetry.space_group = \
      command_line.symmetry.space_group_info()
  phil_params = master_params.format(python_object=params)
  phil_params.show()
  print "#phil __OFF__"
  #
  # Final checks.
  #
  if (params.reference.file_name is None):
    raise Sorry("Required file name is missing: reference.file_name")
  if (params.other.file_name is None):
    raise Sorry("Required file name is missing: other.file_name")
  if (params.output.file_name is None):
    raise Sorry("Required file name is missing: output.file_name")
  #
  # Processing of input PDB files.
  #
  pdb_objs = []
  sites_carts = []
  centers_of_mass = []
  for param_group in [params.reference, params.other]:
    pdb_obj = pdb.hierarchy.input(file_name=param_group.file_name)
    pdb_obj.atoms = pdb_obj.hierarchy.atoms()
    pdb_objs.append(pdb_obj)
    sites_carts.append(pdb_obj.atoms.extract_xyz())
    sites_sel = sites_carts[-1]
    if (param_group.atom_selection is not None):
      sel = pdb_obj.hierarchy.atom_selection_cache().selection(
        param_group.atom_selection)
      sites_sel = sites_sel.select(sel)
    print "Number of selected sites:", sites_sel.size()
    centers_of_mass.append(sites_sel.mean())
  #
  # Consolidation of crystal symmetries.
  #
  crystal_symmetry = command_line.symmetry
  for pdb_obj in pdb_objs:
    crystal_symmetry_from_pdb = pdb_obj.input.crystal_symmetry()
    if (crystal_symmetry_from_pdb is not None):
      crystal_symmetry = crystal_symmetry.join_symmetry(
        other_symmetry=crystal_symmetry_from_pdb,
        force=False)
  if (crystal_symmetry.unit_cell() is None):
    raise Sorry("Unknown unit cell parameters."
      "\n  Use --unit_cell or --symmetry to supply unit cell parameters.")
  if (crystal_symmetry.space_group_info() is None):
    raise Sorry("Unknown space group symmetry."
      "\n  Use --space_group or --symmetry to supply symmetry information.")
  crystal_symmetry.show_summary()
  #
  # Obtain transformation to reference setting.
  #   To ensure all allowed origin shifts are parallel to the basis vectors.
  #
  cb_op_to_ref = crystal_symmetry.change_of_basis_op_to_reference_setting()
  sym_ref = crystal_symmetry.change_basis(cb_op=cb_op_to_ref)
  #
  # Obtain allowed origin shifts.
  #   This is the most convenient interface. Essentially we just need
  #   sgtbx.structure_seminvariants.
  #
  match_symmetry = euclidean_model_matching.euclidean_match_symmetry(
    space_group_info=sym_ref.space_group_info(),
    use_k2l=False,
    use_l2n=False)
  #
  # Compute the symmetry operation which maps the center of mass of
  # "other" closest to the center of mass of "reference."
  #
  centers_frac = [
    sym_ref.unit_cell().fractionalize(cb_op_to_ref.c() * center_cart)
      for center_cart in centers_of_mass]
  dist_info = sgtbx.min_sym_equiv_distance_info(
    sym_ref.special_position_settings().sym_equiv_sites(centers_frac[0]),
    centers_frac[1],
    match_symmetry.continuous_shift_flags)
  sym_op = cb_op_to_ref.inverse().apply(dist_info.sym_op())
  print "Rotation in fractional space:", sym_op.r().as_xyz()
  sym_op = sym_op.as_rational().as_float() \
         + matrix.col(dist_info.continuous_shifts())
  print "Translation in fractional space: (%s)" % (
    ", ".join(["%.6g" % t for t in sym_op.t]))
  #
  centers_frac = [sym_ref.unit_cell().fractionalize(center_cart)
    for center_cart in centers_of_mass]
  sym_center_frac = sym_op * centers_frac[1]
  sym_center_cart = crystal_symmetry.unit_cell().orthogonalize(sym_center_frac)
  print "Centers of mass:"
  print "               Reference: (%s)" % ", ".join(["%8.2f" % v
    for v in centers_of_mass[0]])
  print "          Original other: (%s)" % ", ".join(["%8.2f" % v
    for v in centers_of_mass[1]])
  print "  Symmetry related other: (%s)" % ", ".join(["%8.2f" % v
    for v in sym_center_cart])
  print "Cartesian distance between centers of mass: %.4f" % dist_info.dist()
  #
  # Internal consistency check (in input setting).
  #
  assert approx_equal(crystal_symmetry.unit_cell().distance(
    centers_frac[0], sym_center_frac), dist_info.dist())
  #
  # Transform atomic coordinates of "other."
  #
  sites_frac_other = crystal_symmetry.unit_cell().fractionalize(
    sites_cart=sites_carts[1])
  sites_frac_other_superposed = sym_op * sites_frac_other
  sites_cart_other_superposed = crystal_symmetry.unit_cell().orthogonalize(
    sites_frac=sites_frac_other_superposed)
  #
  # Replace original coordinates with transformed coordinates.
  #
  pdb_objs[1].atoms.set_xyz(new_xyz=sites_cart_other_superposed)
  #
  # Write (selected) transformed coordinates.
  #
  pdb_hierarchy = pdb_objs[1].hierarchy
  if (params.output.atom_selection is not None):
    sel = pdb_hierarchy.atom_selection_cache().selection(
      params.output.atom_selection)
    pdb_hierarchy = pdb_hierarchy.select(atom_selection=sel)
  pdb_hierarchy.write_pdb_file(
    file_name=params.output.file_name,
    crystal_symmetry=crystal_symmetry,
    append_end=True,
    atoms_reset_serial_first_value=1)

if (__name__ == "__main__"):
  run(sys.argv[1:])
