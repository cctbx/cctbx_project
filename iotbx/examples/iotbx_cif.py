#
# Command to run this example:
#   iotbx.python iotbx_cif.py
#
# See also:
#   http://cctbx.sourceforge.net/iotbx_cif
#

def run():
  quartz_as_cif = """\
data_quartz
_space_group_name_H-M_alt         'P 62 2 2'
_cell_length_a                    5.01
_cell_length_b                    5.01
_cell_length_c                    5.47
_cell_angle_alpha                 90
_cell_angle_beta                  90
_cell_angle_gamma                 120
loop_
  _atom_site_label
  _atom_site_type_symbol
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  _atom_site_U_iso_or_equiv
   Si Si 0.500 0.500 0.333 0.200
   O O 0.197 -0.197 0.833 0.200
  """

  import iotbx.cif

  quartz_structure = iotbx.cif.reader(
    input_string=quartz_as_cif).build_crystal_structures()["quartz"]
  quartz_structure.show_summary().show_scatterers()
  print

  # Examine the site symmetry of each scatterer
  for scatterer in quartz_structure.scatterers():
    print "%s:" % scatterer.label, "%8.4f %8.4f %8.4f" % scatterer.site
    site_symmetry = quartz_structure.site_symmetry(scatterer.site)
    print "  point group type:", site_symmetry.point_group_type()
    print "  special position operator:", site_symmetry.special_op_simplified()

  # Let's use scattering factors from the International Tables
  quartz_structure.scattering_type_registry(table="it1992")

  # Now calculate some structure factors
  f_calc = quartz_structure.structure_factors(d_min=2).f_calc()
  f_calc_sq = f_calc.as_intensity_array()
  f_calc_sq.show_summary().show_array()

  # Output the intensities to a CIF file
  f_calc_sq.as_cif_simple(
    array_type="calc", data_name="quartz", out=open("quartz.hkl", "wb"))
  from iotbx.cif import model

  # Create an instance of model.cif, equivalent to a full CIF file
  cif = model.cif()

  # Create an instance of model.block, equivalent to a CIF data block
  cif_block = model.block()

  # Add the unit cell parameters to the cif_block
  unit_cell = quartz_structure.unit_cell()
  params = unit_cell.parameters()
  cif_block["_cell_length_a"] = params[0]
  cif_block["_cell_length_b"] = params[1]
  cif_block["_cell_length_c"] = params[2]
  cif_block["_cell_angle_alpha"] = params[3]
  cif_block["_cell_angle_beta"] = params[4]
  cif_block["_cell_angle_gamma"] = params[5]
  cif_block["_cell_volume"] = unit_cell.volume()

  # now we will create a CIF loop containing the space group symmetry operations
  space_group = quartz_structure.space_group()
  symop_loop = model.loop(header=("_space_group_symop_id",
                            "_space_group_symop_operation_xyz"))
  for symop_id, symop in enumerate(space_group):
    symop_loop.add_row((symop_id + 1, symop.as_xyz()))

  # add the symop_loop and space group items to the cif_block
  space_group_type = quartz_structure.space_group_info().type()
  cif_block["_space_group_crystal_system"] = space_group.crystal_system().lower()
  cif_block["_space_group_IT_number"] = space_group_type.number()
  cif_block["_space_group_name_H-M_alt"] = space_group_type.lookup_symbol()
  cif_block["_space_group_name_Hall"] = space_group_type.hall_symbol()
  cif_block.add_loop(symop_loop)

  # add cif_block to the cif object with the data block name "quartz"
  cif["quartz"] = cif_block

  # print the cif object to standard output
  print cif

  from iotbx.cif import validation

  cif_model = iotbx.cif.reader(input_string=quartz_as_cif).model()
  cif_model["quartz"]["_diffrn_radiation_probe"] = "xray"
  cif_model["quartz"]["_space_group_crystal_system"] = "Monoclinic"
  cif_model["quartz"]["_space_group_IT_number"] = "one hundred and eighty"
  symop_loop.add_column("_space_group_symop_sg_id", [1]*12)
  cif_model["quartz"].add_loop(symop_loop)

  cif_core_dic = validation.smart_load_dictionary(name="cif_core.dic")
  cif_model.validate(cif_core_dic, show_warnings=True)

if __name__ == "__main__":
  run()
