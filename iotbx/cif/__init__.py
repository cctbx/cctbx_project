import libtbx.load_env
has_antlr3 = libtbx.env.has_module('antlr3')

if has_antlr3:
  import boost.python
  ext = boost.python.import_ext("iotbx_cif_ext")

from cctbx.array_family import flex
from cctbx import adptbx, covariance, crystal, sgtbx
from iotbx.cif import model, builders
from libtbx.containers import OrderedDict
from libtbx.utils import format_float_with_standard_uncertainty \
     as format_float_with_su
from libtbx.utils import Sorry
from scitbx import matrix

import math, sys

class reader:

  def __init__(self, file_path=None, file_object=None, input_string=None,
               builder=None, max_errors=50):
    assert [file_path, file_object, input_string].count(None) == 2
    assert has_antlr3
    self.file_path = file_path
    if builder is None:
      builder = builders.cif_model_builder()
    self.builder = builder
    if file_path is not None:
      if isinstance(file_path, unicode):
        file_object = open(file_path, 'rb')
      else:
        self.parser = ext.fast_reader(file_path, builder)
    if file_object is not None:
      input_string = file_object.read()
    if input_string is not None:
      self.parser = ext.fast_reader(input_string, builder)
    self.show_errors(max_errors)

  def model(self):
    return self.builder.model()

  def error_count(self):
    return self.parser.lexer_errors().size()\
           + self.parser.parser_errors().size()

  def show_errors(self, max_errors=50, out=None):
    if out is None: out = sys.stdout
    for msg in self.parser.lexer_errors()[:max_errors]:
      print >> out, msg
    for msg in self.parser.parser_errors()[:max_errors]:
      print >> out, msg

  def get_block_else_raise(self, data_block_name):
    block = self.model().get(key=data_block_name)
    if (block is None):
      if (self.file_path is None):
        msg = 'Unknown CIF data block name: "%s"' % data_block_name
      else:
        msg = 'Unknown CIF data block name "%s" in file: "%s"' % (
          data_block_name, self.file_path)
      raise RuntimeError(msg)
    return block

  def build_crystal_structure(self, data_block_name="global"):
    block = self.get_block_else_raise(data_block_name=data_block_name)
    return builders.crystal_structure_builder(cif_block=block).structure

  def build_miller_arrays(self, data_block_name="global"):
    block = self.get_block_else_raise(data_block_name=data_block_name)
    return builders.miller_array_builder(cif_block=block).arrays().values()

fast_reader = reader # XXX backward compatibility 2010-08-25

class crystal_symmetry_as_cif_block:

  def __init__(self, crystal_symmetry, cell_covariance_matrix=None):
    self.cif_block = model.block()
    sym_loop = model.loop(data=OrderedDict((
      ('_space_group_symop_id',
       range(1, len(crystal_symmetry.space_group())+1)),
      ('_space_group_symop_operation_xyz',
       [s.as_xyz() for s in crystal_symmetry.space_group()]))))
    self.cif_block.add_loop(sym_loop)
    #
    sg_type = crystal_symmetry.space_group_info().type()
    sg = sg_type.group()
    self.cif_block['_space_group_crystal_system'] = sg.crystal_system().lower()
    self.cif_block['_space_group_IT_number'] = sg_type.number()
    self.cif_block['_space_group_name_H-M_alt'] = sg_type.universal_hermann_mauguin_symbol()
    self.cif_block['_space_group_name_Hall'] = sg_type.hall_symbol()
    #
    uc = crystal_symmetry.unit_cell()
    params = list(uc.parameters())
    volume = uc.volume()
    if cell_covariance_matrix is not None:
      diag = cell_covariance_matrix.matrix_packed_u_diagonal()
      for i in range(6):
        if diag[i] > 0:
          params[i] = format_float_with_su(params[i], math.sqrt(diag[i]))
      d_v_d_params = matrix.row(uc.d_volume_d_params())
      vcv = matrix.sqr(
        cell_covariance_matrix.matrix_packed_u_as_symmetric())
      var_v = (d_v_d_params * vcv).dot(d_v_d_params)
      volume = format_float_with_su(volume, math.sqrt(var_v))
    a,b,c,alpha,beta,gamma = params
    self.cif_block['_cell_length_a'] = a
    self.cif_block['_cell_length_b'] = b
    self.cif_block['_cell_length_c'] = c
    self.cif_block['_cell_angle_alpha'] = alpha
    self.cif_block['_cell_angle_beta'] = beta
    self.cif_block['_cell_angle_gamma'] = gamma
    self.cif_block['_cell_volume'] = volume


class xray_structure_as_cif_block(crystal_symmetry_as_cif_block):
  sources = {
    "it1992": "International Tables Volume C Table 6.1.1.4 (pp. 500-502)",
    "wk1995": "Waasmaier & Kirfel (1995), Acta Cryst. A51, 416-431",
  }

  def __init__(self, xray_structure, covariance_matrix=None,
               cell_covariance_matrix=None):
    crystal_symmetry_as_cif_block.__init__(
      self, xray_structure.crystal_symmetry(),
      cell_covariance_matrix=cell_covariance_matrix)
    scatterers = xray_structure.scatterers()
    uc = xray_structure.unit_cell()
    if covariance_matrix is not None:
      param_map = xray_structure.parameter_map()
      covariance_diagonal = covariance_matrix.matrix_packed_u_diagonal()
      u_star_to_u_cif_linear_map_pow2 = flex.pow2(flex.double(
        uc.u_star_to_u_cif_linear_map()))
      u_star_to_u_iso_linear_form = matrix.row(
        uc.u_star_to_u_iso_linear_form())
    fmt = "%.6f"

    # _atom_site_* loop
    atom_site_loop = model.loop(header=(
      '_atom_site_label', '_atom_site_type_symbol',
      '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z',
      '_atom_site_U_iso_or_equiv', '_atom_site_adp_type',
      '_atom_site_occupancy'))
    fp_fdp_table = {}
    for i_seq, sc in enumerate(scatterers):
      # site
      if covariance_matrix is not None and sc.flags.grad_site():
        site = []
        for i in range(3):
          idx = param_map[i_seq].site
          if idx > -1:
            var = covariance_diagonal[idx+i]
          else: var = 0
          if var > 0:
            site.append(format_float_with_su(sc.site[i], math.sqrt(var)))
          else: site.append(fmt % sc.site[i])
      else:
        site = [fmt % sc.site[i] for i in range(3)]
      # u_eq
      if (covariance_matrix is not None and
          (sc.flags.grad_u_iso() or sc.flags.grad_u_aniso())):
        if sc.flags.grad_u_iso():
          u_iso_or_equiv = format_float_with_su(
            sc.u_iso, math.sqrt(covariance.variance_for_u_iso(
              i_seq, covariance_matrix, param_map)))
        else:
          cov = covariance.extract_covariance_matrix_for_u_aniso(
            i_seq, covariance_matrix, param_map).matrix_packed_u_as_symmetric()
          var = (u_star_to_u_iso_linear_form * matrix.sqr(cov)
                 ).dot(u_star_to_u_iso_linear_form)
          u_iso_or_equiv = format_float_with_su(
            sc.u_iso_or_equiv(uc), math.sqrt(var))
      else:
        u_iso_or_equiv = fmt % sc.u_iso_or_equiv(uc)
      if sc.flags.use_u_aniso():
        adp_type = 'Uani'
      else:
        adp_type = 'Uiso'
      atom_site_loop.add_row((
        sc.label, sc.scattering_type, site[0], site[1], site[2], u_iso_or_equiv,
        adp_type, fmt%sc.occupancy))
      fp_fdp_table.setdefault(sc.scattering_type, (sc.fp, sc.fdp))
    self.cif_block.add_loop(atom_site_loop)

    # _atom_site_aniso_* loop
    aniso_scatterers = scatterers.select(scatterers.extract_use_u_aniso())
    if aniso_scatterers.size():
      labels = list(scatterers.extract_labels())
      aniso_loop = model.loop(header=('_atom_site_aniso_label',
                                      '_atom_site_aniso_U_11',
                                      '_atom_site_aniso_U_22',
                                      '_atom_site_aniso_U_33',
                                      '_atom_site_aniso_U_12',
                                      '_atom_site_aniso_U_13',
                                      '_atom_site_aniso_U_23'))
      for sc in aniso_scatterers:
        u_cif = adptbx.u_star_as_u_cif(uc, sc.u_star)
        if covariance_matrix is not None:
          row = [sc.label]
          idx = param_map[labels.index(sc.label)].u_aniso
          if idx > -1:
            var = covariance_diagonal[idx:idx+6] * u_star_to_u_cif_linear_map_pow2
            for i in range(6):
              if var[i] > 0:
                row.append(
                  format_float_with_su(u_cif[i], math.sqrt(var[i])))
              else:
                row.append(fmt%u_cif[i])
          else:
            row = [sc.label] + [fmt%u_cif[i] for i in range(6)]
        else:
          row = [sc.label] + [fmt%u_cif[i] for i in range(6)]
        aniso_loop.add_row(row)
      self.cif_block.add_loop(aniso_loop)

    # _atom_type_* loop
    atom_type_loop = model.loop(header=('_atom_type_symbol',
                                        '_atom_type_scat_dispersion_real',
                                        '_atom_type_scat_dispersion_imag',
                                        '_atom_type_scat_source'))
    scattering_type_registry = xray_structure.scattering_type_registry()
    params = xray_structure.scattering_type_registry_params
    scat_source = self.sources.get(params.table)
    scattering_type_registry = xray_structure.scattering_type_registry()
    for atom_type, gaussian in scattering_type_registry.as_type_gaussian_dict().iteritems():
      if params.custom_dict and atom_type in params.custom_dict:
        scat_source = "Custom %i-Gaussian" %gaussian.n_terms()
      elif scat_source is None:
        scat_source = """\
%i-Gaussian fit: Grosse-Kunstleve RW, Sauter NK, Adams PD:
Newsletter of the IUCr Commission on Crystallographic Computing 2004, 3, 22-31."""
        scat_source = scat_source %gaussian.n_terms()
      fp, fdp = fp_fdp_table[atom_type]
      atom_type_loop.add_row((atom_type, "%.5f" %fp, "%.5f" %fdp, scat_source))
    self.cif_block.add_loop(atom_type_loop)


class miller_indices_as_cif_loop:

  def __init__(self, indices, prefix='_refln_'):
    self.refln_loop = model.loop(header=(
      '%sindex_h' %prefix, '%sindex_k' %prefix, '%sindex_l' %prefix))
    for hkl in indices:
      self.refln_loop.add_row(hkl)


class miller_arrays_as_cif_block(crystal_symmetry_as_cif_block,
                                 miller_indices_as_cif_loop):

  def __init__(self, array, array_type=None,
               column_name=None, column_names=None,
               miller_index_prefix='_refln_'):
    crystal_symmetry_as_cif_block.__init__(self, array.crystal_symmetry())
    miller_indices_as_cif_loop.__init__(
      self, array.indices(), prefix=miller_index_prefix)
    self.indices = array.indices()
    self.add_miller_array(array, array_type, column_name, column_names)
    self.cif_block.add_loop(self.refln_loop)

  def add_miller_array(self, array, array_type=None,
                       column_name=None, column_names=None):
    """
    Accepts a miller array, and one of array_type, column_name or column_names.
    """

    assert [array_type, column_name, column_names].count(None) == 2
    if array_type is not None:
      assert array_type in ('calc', 'meas')
    elif column_name is not None:
      column_names = [column_name]
    assert array.size() == self.indices.size()
    if array.is_complex_array():
      if column_names is None:
        column_names = ['_refln_F_%s' %array_type, '_refln_phase_%s' %array_type]
      else: assert len(column_names) == 2
      if '_A_' in column_names[0] and '_B_' in column_names[1]:
        data = [flex.real(array.data()).as_string(),
                 flex.imag(array.data()).as_string()]
      else:
        data = [flex.abs(array.data()).as_string(),
                 array.phases().data().as_string()]
    else:
      if array_type is not None:
        if array.is_xray_intensity_array():
          obs_ext = 'squared_'
        else: obs_ext = ''
        column_names = ['_refln_F_%s%s' %(obs_ext, array_type)]
        if array.sigmas() is not None:
          column_names.append('_refln_F_%ssigma' %obs_ext)
      if isinstance(array.data(), flex.std_string):
        data = [array.data()]
      else:
        data = [array.data().as_string()]
      if array.sigmas() is not None and len(column_names) == 2:
        data.append(array.sigmas().as_string())
    columns = OrderedDict(zip(column_names, data))
    for key in columns:
      assert key not in self.refln_loop
    self.refln_loop.add_columns(columns)

class distances_as_cif_loop(object):

  def __init__(self,
               pair_asu_table,
               site_labels,
               sites_frac=None,
               sites_cart=None,
               covariance_matrix=None,
               cell_covariance_matrix=None,
               parameter_map=None,
               eps=2e-16):
    assert [sites_frac, sites_cart].count(None) == 1
    fmt = "%.4f"
    asu_mappings = pair_asu_table.asu_mappings()
    space_group_info = sgtbx.space_group_info(group=asu_mappings.space_group())
    unit_cell = asu_mappings.unit_cell()
    self.loop = model.loop(header=(
      "_geom_bond_atom_site_label_1",
      "_geom_bond_atom_site_label_2",
      "_geom_bond_distance",
      "_geom_bond_site_symmetry_2"
    ))
    distances = crystal.calculate_distances(
      pair_asu_table, sites_frac,
      covariance_matrix=covariance_matrix,
      cell_covariance_matrix=cell_covariance_matrix,
      parameter_map=parameter_map)
    for d in distances:
      if site_labels[d.i_seq].startswith('H') or site_labels[d.j_seq].startswith('H'):
        continue
      if d.variance is not None and d.variance > eps:
        distance = format_float_with_su(d.distance, math.sqrt(d.variance))
      else:
        distance = fmt % d.distance
      sym_code = space_group_info.cif_symmetry_code(d.rt_mx_ji)
      if sym_code == "1": sym_code = "."
      self.loop.add_row((site_labels[d.i_seq],
                         site_labels[d.j_seq],
                         distance,
                         sym_code))
    self.distances = distances.distances
    self.pair_counts = distances.pair_counts

class angles_as_cif_loop(object):

  def __init__(self,
               pair_asu_table,
               site_labels,
               sites_frac=None,
               sites_cart=None,
               covariance_matrix=None,
               cell_covariance_matrix=None,
               parameter_map=None,
               eps=2e-16):
    assert [sites_frac, sites_cart].count(None) == 1
    fmt = "%.1f"
    asu_mappings = pair_asu_table.asu_mappings()
    space_group_info = sgtbx.space_group_info(group=asu_mappings.space_group())
    unit_cell = asu_mappings.unit_cell()
    self.loop = model.loop(header=(
      "_geom_angle_atom_site_label_1",
      "_geom_angle_atom_site_label_2",
      "_geom_angle_atom_site_label_3",
      "_geom_angle",
      "_geom_angle_site_symmetry_1",
      "_geom_angle_site_symmetry_3"
    ))
    angles = crystal.calculate_angles(
      pair_asu_table, sites_frac,
      covariance_matrix=covariance_matrix,
      cell_covariance_matrix=cell_covariance_matrix,
      parameter_map=parameter_map)
    for a in angles:
      i_seq, j_seq, k_seq = a.i_seqs
      if site_labels[i_seq].startswith('H') or site_labels[k_seq].startswith('H'):
        continue
      sym_code_ji = space_group_info.cif_symmetry_code(a.rt_mx_ji)
      sym_code_ki = space_group_info.cif_symmetry_code(a.rt_mx_ki)
      if sym_code_ji == "1": sym_code_ji = "."
      if sym_code_ki == "1": sym_code_ki = "."
      if a.variance is not None and a.variance > eps:
        angle = format_float_with_su(a.angle, math.sqrt(a.variance))
      else:
        angle = fmt % a.angle
      self.loop.add_row((site_labels[i_seq],
                         site_labels[j_seq],
                         site_labels[k_seq],
                         angle,
                         sym_code_ji,
                         sym_code_ki,
                         ))
    self.angles = angles.angles


def cctbx_data_structure_from_cif(
  file_object=None, file_path=None, data_structure_builder=None,
  block_heading=None, **kwds):
  assert data_structure_builder is not None
  cif_model = reader(file_path=file_path, file_object=file_object).model()
  if not len(cif_model):
    raise Sorry("No data block found in CIF")
  if block_heading is not None:
    return data_structure_builder(cif_model[block_heading], **kwds)
  else:
    errors = []
    for block in cif_model.values():
      try:
        return data_structure_builder(block)
      except Exception, e:
        errors.append(e)
        continue
    if errors:
      raise errors[0]
