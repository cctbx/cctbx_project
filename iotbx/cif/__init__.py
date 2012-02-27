"""
R. J. Gildea, L. J. Bourhis, O. V. Dolomanov, R. W. Grosse-Kunstleve,
H. Puschmann, P. D. Adams and J. A. K. Howard:
iotbx.cif: a comprehensive CIF toolbox.
J. Appl. Cryst. (2011). 44, 1259-1263.

http://dx.doi.org/10.1107/S0021889811041161

http://cctbx.sourceforge.net/iotbx_cif

"""

import boost.python
ext = boost.python.import_ext("iotbx_cif_ext")

from cctbx.array_family import flex
from cctbx import adptbx
from cctbx import covariance
from iotbx.cif import model, builders, geometry
from libtbx.containers import OrderedDict
from libtbx.utils import format_float_with_standard_uncertainty \
     as format_float_with_su
from libtbx.utils import Sorry
from libtbx.utils import flat_list
from libtbx import smart_open
from scitbx import matrix

import math, sys

distances_as_cif_loop = geometry.distances_as_cif_loop
angles_as_cif_loop = geometry.angles_as_cif_loop

class CifParserError(Sorry):
  __orig_module__ = __module__
  __module__ = Exception.__module__

class reader(object):

  def __init__(self,
               file_path=None,
               file_object=None,
               input_string=None,
               cif_object=None,
               builder=None,
               raise_if_errors=True,
               strict=True):
    assert [file_path, file_object, input_string].count(None) == 2
    self.file_path = file_path
    if builder is None:
      builder = builders.cif_model_builder(cif_object)
    else: assert cif_object is None
    self.builder = builder
    if file_path is not None:
      file_object = smart_open.for_reading(file_path)
    else:
      file_path = "memory"
    if file_object is not None:
      input_string = file_object.read()
    self.parser = ext.fast_reader(builder, input_string, file_path, strict)
    if raise_if_errors and len(self.parser.lexer_errors()):
      raise CifParserError(self.parser.lexer_errors()[0])
    if raise_if_errors and len(self.parser.parser_errors()):
      raise CifParserError(self.parser.parser_errors()[0])

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

  def build_crystal_structures(self, data_block_name=None):
    xray_structures = cctbx_data_structures_from_cif(
      cif_model=self.model(),
      file_path=self.file_path,
      data_block_name=data_block_name,
      data_structure_builder=builders.crystal_structure_builder).xray_structures
    if data_block_name is not None:
      return xray_structures[data_block_name]
    else:
      return xray_structures

  def build_miller_arrays(self,
                          data_block_name=None,
                          base_array_info=None):
    arrays = cctbx_data_structures_from_cif(
      cif_model=self.model(),
      file_path=self.file_path,
      data_block_name=data_block_name,
      data_structure_builder=builders.miller_array_builder,
      base_array_info=base_array_info).miller_arrays
    if data_block_name is not None:
      return arrays[data_block_name]
    else:
      return arrays

  def as_miller_arrays(self, data_block_name=None,
                       crystal_symmetry=None,
                       force_symmetry=False,
                       merge_equivalents=True,
                       base_array_info=None):
    if base_array_info is None:
      from cctbx import miller
      base_array_info = miller.array_info(
        source=self.file_path, source_type="cif")
    if data_block_name is not None:
      arrays = self.build_miller_arrays(
        data_block_name=data_block_name,
        base_array_info=base_array_info).values()
    else:
      arrays = flat_list([
        arrays.values() for arrays in
        self.build_miller_arrays(base_array_info=base_array_info).values()])
    other_symmetry=crystal_symmetry
    for i, array in enumerate(arrays):
      if crystal_symmetry is not None:
        crystal_symmetry = other_symmetry.join_symmetry(
          other_symmetry=array.crystal_symmetry(),
          force=force_symmetry)
        arrays[i] = array.customized_copy(crystal_symmetry=crystal_symmetry)
        arrays[i].set_info(array.info())
    return arrays

fast_reader = reader # XXX backward compatibility 2010-08-25

class crystal_symmetry_as_cif_block(object):

  def __init__(self, crystal_symmetry,
               cell_covariance_matrix=None,
               format="coreCIF"):
    self.format = format.lower()
    assert self.format in ("corecif", "mmcif")
    if self.format == "mmcif": self.separator = '.'
    else: self.separator = '_'
    self.cif_block = model.block()
    sg_prefix = '_space_group%s' %self.separator
    cell_prefix = '_cell%s' %self.separator
    sym_loop = model.loop(data=OrderedDict((
      (sg_prefix+'symop_id',
       range(1, len(crystal_symmetry.space_group())+1)),
      (sg_prefix+'symop_operation_xyz',
       [s.as_xyz() for s in crystal_symmetry.space_group()]))))
    self.cif_block.add_loop(sym_loop)
    #
    sg_type = crystal_symmetry.space_group_info().type()
    sg = sg_type.group()
    self.cif_block[sg_prefix+'crystal_system'] = sg.crystal_system().lower()
    self.cif_block[sg_prefix+'IT_number'] = sg_type.number()
    self.cif_block[sg_prefix+'name_H-M_alt'] = sg_type.lookup_symbol()
    self.cif_block[sg_prefix+'name_Hall'] = sg_type.hall_symbol()
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
    self.cif_block[cell_prefix+'length_a'] = a
    self.cif_block[cell_prefix+'length_b'] = b
    self.cif_block[cell_prefix+'length_c'] = c
    self.cif_block[cell_prefix+'angle_alpha'] = alpha
    self.cif_block[cell_prefix+'angle_beta'] = beta
    self.cif_block[cell_prefix+'angle_gamma'] = gamma
    self.cif_block[cell_prefix+'volume'] = volume


class xray_structure_as_cif_block(crystal_symmetry_as_cif_block):
  sources = {
    "it1992": "International Tables Volume C Table 6.1.1.4 (pp. 500-502)",
    "wk1995": "Waasmaier & Kirfel (1995), Acta Cryst. A51, 416-431",
  }
  inelatic_references = {
    "henke" : "Henke, Gullikson and Davis, At. Data and Nucl. Data Tables, 1993, 54, 2",
    "sasaki" : "Sasaki, KEK Report, 1989, 88-14, 1",
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
                                        '_atom_type_scat_source',
                                        '_atom_type_scat_dispersion_source'))
    scattering_type_registry = xray_structure.scattering_type_registry()
    params = xray_structure.scattering_type_registry_params
    scat_source = self.sources.get(params.table)
    disp_source = self.inelatic_references.get(
      xray_structure.inelastic_form_factors_source)
    # custom?
    if disp_source is None:
      disp_source = xray_structure.inelastic_form_factors_source
    if disp_source is None:
      disp_source = "undefined"
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
      atom_type_loop.add_row(
        (atom_type, "%.5f" %fp, "%.5f" %fdp, scat_source, disp_source))
    self.cif_block.add_loop(atom_type_loop)


class miller_indices_as_cif_loop(object):

  def __init__(self, indices, prefix='_refln', separator='_'):
    prefix += separator
    self.refln_loop = model.loop(header=(
      '%sindex_h' %prefix, '%sindex_k' %prefix, '%sindex_l' %prefix))
    for hkl in indices:
      self.refln_loop.add_row(hkl)


class miller_arrays_as_cif_block(crystal_symmetry_as_cif_block,
                                 miller_indices_as_cif_loop):

  def __init__(self, array, array_type=None,
               column_name=None, column_names=None,
               miller_index_prefix='_refln',
               format="coreCIF"):
    crystal_symmetry_as_cif_block.__init__(
      self, array.crystal_symmetry(), format=format)
    miller_indices_as_cif_loop.__init__(
      self, array.indices(), prefix=miller_index_prefix,
      separator=self.separator)
    self.prefix = miller_index_prefix + self.separator
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
        column_names = [self.prefix+'F_'+array_type,
                        self.prefix+'phase_'+array_type]
      else: assert len(column_names) == 2
      if (('_A_' in column_names[0] and '_B_' in column_names[1]) or
          ('.A_' in column_names[0] and '.B_' in column_names[1])):
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
        column_names = [self.prefix+'F_'+obs_ext+array_type]
        if array.sigmas() is not None:
          column_names.append(self.prefix+'F_'+obs_ext+'sigma')
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


class cctbx_data_structures_from_cif(object):
  def __init__(self,
               file_object=None,
               file_path=None,
               cif_model=None,
               data_structure_builder=None,
               data_block_name=None,
               base_array_info=None,
               **kwds):
    assert file_object is None or cif_model is None
    if data_structure_builder is None:
      data_structure_builders = (
        builders.miller_array_builder, builders.crystal_structure_builder)
    else:
      assert data_structure_builder in (
        builders.miller_array_builder, builders.crystal_structure_builder)
      data_structure_builders = (data_structure_builder,)

    self.xray_structures = OrderedDict()
    self.miller_arrays = OrderedDict()
    if cif_model is None:
      cif_model = reader(file_path=file_path, file_object=file_object).model()
    if not len(cif_model):
      raise Sorry("No data block found in CIF")
    if data_block_name is not None and not data_block_name in cif_model:
      if (file_path is None):
        msg = 'Unknown CIF data block name: "%s"' % data_block_name
      else:
        msg = 'Unknown CIF data block name "%s" in file: "%s"' % (
          data_block_name, file_path)
      raise RuntimeError(msg)
    errors = []
    for key, block in cif_model.items():
      if data_block_name is not None and key != data_block_name: continue
      for builder in data_structure_builders:
        if builder == builders.crystal_structure_builder:
          if '_atom_site_fract_x' in block or '_atom_site_Cartn_x' in block:
            self.xray_structures.setdefault(key, builder(block).structure)
        elif builder == builders.miller_array_builder:
          if base_array_info is not None:
            base_array_info = base_array_info.customized_copy(labels=[key])
          if '_refln_index_h' in block or '_refln.index_h' in block:
            self.miller_arrays.setdefault(
              key, builder(block, base_array_info=base_array_info).arrays())
