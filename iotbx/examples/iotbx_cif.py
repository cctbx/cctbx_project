from __future__ import division

import iotbx.pdb
import iotbx.cif
from iotbx.cif import model
from iotbx.cif import validation

import mmtbx.model


#
# Command to run this example:
#   iotbx.python iotbx_cif.py
#
# See also:
#   http://cctbx.sourceforge.net/iotbx_cif
#
#




# hard-coded mmCIF words (probably outdated)
# does not pass its own validation



quartz_as_cif = """\
data_quartz
_cell.angle_beta                  90.000
_cell.angle_gamma                 120.000
_cell.length_b                    5.010
_cell.length_c                    5.470
_cell.angle_alpha                 90.000
_cell.volume                      118.903
_cell.length_a                    5.010
_space_group.crystal_system       hexagonal
_space_group.name_H-M_alt         'P 62 2 2'
_space_group.IT_number            180
_space_group.name_Hall            ' P 62 2 (x,y,z+1/3)'
_symmetry.space_group_name_H-M    'P 62 2 2'
_symmetry.Int_Tables_number       180
_symmetry.space_group_name_Hall   ' P 62 2 (x,y,z+1/3)'

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.pdbx_PDB_model_num
   ATOM 1 SI . SI . 1 ? 1.25200 2.16900 1.82300 1.000 15.79000 SI ? A ? 1 1
   ATOM 2 O . O . 2 ? 1.48000 -0.85500 4.55800 1.000 15.79000 O ? A ? 2 1
  """

def run():
  inp = iotbx.pdb.input(source_info=None, lines=quartz_as_cif.split('\n'))
  modelm = mmtbx.model.manager(
      model_input=inp)
  quartz_structure = modelm.get_xray_structure()

  # Examine the site symmetry of each scatterer
  for scatterer in quartz_structure.scatterers():
    print "%s:" % scatterer.label, "%8.4f %8.4f %8.4f" % scatterer.site
    site_symmetry = quartz_structure.site_symmetry(scatterer.site)
    print "  point group type:", site_symmetry.point_group_type()
    print "  special position operator:", site_symmetry.special_op_simplified()

  # Let's use scattering factors from the International Tables
  modelm.setup_scattering_dictionaries(scattering_table="it1992")
  quartz_structure = modelm.get_xray_structure()

  # Now calculate some structure factors
  f_calc = quartz_structure.structure_factors(d_min=2).f_calc()
  f_calc_sq = f_calc.as_intensity_array()
  f_calc_sq.show_summary().show_array()

  # Output the intensities to a CIF file
  # This probably should be deprecated as well
  f_calc_sq.as_cif_simple(
    array_type="calc", data_name="quartz", out=open("quartz-intensities.cif", "wb"))

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

  cif_model = iotbx.cif.reader(input_string=quartz_as_cif).model()
  cif_model["quartz"]["_diffrn_radiation_probe"] = "xray"
  cif_model["quartz"]["_space_group_crystal_system"] = "Monoclinic"
  cif_model["quartz"]["_space_group_IT_number"] = "one hundred and eighty"
  symop_loop.add_column("_space_group_symop_sg_id", [1]*12)
  cif_model["quartz"].add_loop(symop_loop)

  pdbx_v50 = validation.smart_load_dictionary(name="mmcif_pdbx_v50.dic")
  print "validation"
  cif_model.validate(pdbx_v50, show_warnings=True)

if __name__ == "__main__":
  run()
