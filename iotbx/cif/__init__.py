#import boost.python
#ext = boost.python.import_ext("iotbx_cif_ext")

import libtbx.load_env
has_antlr3 = libtbx.env.has_module('antlr3')

from cctbx import adptbx, crystal
from cctbx.xray import structure
from iotbx.cif import model, builders

def python_reader(file_path=None, file_object=None, input_string=None,
                  builder=None):
  assert [file_path, file_object, input_string].count(None) == 2
  assert has_antlr3
  if builder is None:
    builder = builders.cif_model_builder()
  from iotbx.cif import cifLexer, cifParser
  import antlr3
  if file_object is not None:
    char_stream = antlr3.ANTLRInputStream(file_object)
  elif file_path is not None:
    char_stream = antlr3.ANTLRFileStream(file_path)
  else:
    char_stream = antlr3.ANTLRStringStream(input_string)
  lexer = cifLexer.cifLexer(char_stream)
  tokens = antlr3.CommonTokenStream(lexer)
  parser = cifParser.cifParser(tokens)
  parser.parse(builder=builder)
  return builder


class crystal_symmetry_as_cif_block:

  def __init__(self, crystal_symmetry):
    self.cif_block = model.block()
    sym_loop = model.loop(data={
      '_space_group_symop_operation_xyz':
      [s.as_xyz() for s in crystal_symmetry.space_group()],
      '_space_group_symop_id':
      range(1, len(crystal_symmetry.space_group())+1)})
    self.cif_block.add_loop(sym_loop)
    #
    sg_type = crystal_symmetry.space_group_info().type()
    sg = sg_type.group()
    self.cif_block['_space_group_crystal_system'] = sg.crystal_system()
    self.cif_block['_space_group_IT_number'] = sg_type.number()
    self.cif_block['_space_group_name_H-M_alt'] = sg_type.universal_hermann_mauguin_symbol()
    self.cif_block['_space_group_name_Hall'] = sg_type.hall_symbol()
    #
    uc = crystal_symmetry.unit_cell()
    a,b,c,alpha,beta,gamma = uc.parameters()
    self.cif_block['_cell_length_a'] = a
    self.cif_block['_cell_length_b'] = b
    self.cif_block['_cell_length_c'] = c
    self.cif_block['_cell_angle_alpha'] = alpha
    self.cif_block['_cell_angle_beta'] = beta
    self.cif_block['_cell_angle_gamma'] = gamma

class xray_structure_as_cif_block(crystal_symmetry_as_cif_block):

  def __init__(self, xray_structure):
    crystal_symmetry_as_cif_block.__init__(
      self, xray_structure.crystal_symmetry())
    scatterers = xray_structure.scatterers()
    atom_site_loop = model.loop(header=(
      '_atom_site_label', '_atom_site_type_symbol',
      '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z',
      '_atom_site_U_iso_or_equiv', '_atom_site_occupancy'))
    uc = xray_structure.unit_cell()
    for sc in scatterers:
      atom_site_loop.add_row((
        sc.label, sc.scattering_type, sc.site[0], sc.site[1], sc.site[2],
        sc.u_iso_or_equiv(uc), sc.occupancy))
    aniso_scatterers = scatterers.select(scatterers.extract_use_u_aniso())
    aniso_loop = model.loop(header=('_atom_site_aniso_label',
                                    '_atom_site_aniso_U_11',
                                    '_atom_site_aniso_U_22',
                                    '_atom_site_aniso_U_33',
                                    '_atom_site_aniso_U_12',
                                    '_atom_site_aniso_U_13',
                                    '_atom_site_aniso_U_23'))
    for sc in aniso_scatterers:
      u_cif = adptbx.u_star_as_u_cif(uc, sc.u_star)
      aniso_loop.add_row((
        sc.label, u_cif[0], u_cif[1], u_cif[2], u_cif[3], u_cif[4], u_cif[5]))

    self.cif_block.add_loop(atom_site_loop)
    self.cif_block.add_loop(aniso_loop)

class miller_array_as_cif_block(crystal_symmetry_as_cif_block):

  def __init__(self, array):
    crystal_symmetry_as_cif_block.__init__(self, array.crystal_symmetry())


def cctbx_data_structure_from_cif(
  file_object=None, file_path=None, data_structure_builder=None,
  reader=None, block_heading=None, **kwds):
  assert data_structure_builder is not None
  if reader is None:
    reader = python_reader
  cif_model = reader(file_path=file_path, file_object=file_object).model()
  if block_heading is not None:
    return data_structure_builder(cif_model[block_heading], **kwds)
  else:
    xs = None
    for block in cif_model.values():
      try:
        return data_structure_builder(block)
      except:
        continue
