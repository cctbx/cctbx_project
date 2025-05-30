"""
Tools for reading and writing mmCIF files.

R. J. Gildea, L. J. Bourhis, O. V. Dolomanov, R. W. Grosse-Kunstleve,
H. Puschmann, P. D. Adams and J. A. K. Howard:
iotbx.cif: a comprehensive CIF toolbox.
J. Appl. Cryst. (2011). 44, 1259-1263.

https://doi.org/10.1107/S0021889811041161

"""
from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
from six.moves import range
from six.moves import zip
import six
ext = bp.import_ext("iotbx_cif_ext")

from cctbx.array_family import flex
from cctbx import miller
from iotbx.cif import model, builders, geometry
from libtbx.containers import OrderedDict
from libtbx.utils import Sorry
from libtbx.utils import flat_list
from libtbx.utils import detect_binary_file
from libtbx import smart_open

import sys

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
    self.original_arrays = None
    if file_path is not None:
      file_object = smart_open.for_reading(file_path)
    else:
      file_path = "memory"
    if file_object is not None:
      input_string = file_object.read()
      file_object.close()
    # check input_string for binary, and abort if necessary
    binary_detector = detect_binary_file()
    binary_detector.monitor_initial = min(
      len(input_string), binary_detector.monitor_initial)
    if binary_detector.is_binary_file(block=input_string):
      raise CifParserError("Binary file detected, aborting parsing.")
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
      print(msg, file=out)
    for msg in self.parser.parser_errors()[:max_errors]:
      print(msg, file=out)

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
    cctbxdat = cctbx_data_structures_from_cif(
      cif_model=self.model(),
      file_path=self.file_path,
      data_block_name=data_block_name,
      data_structure_builder=builders.miller_array_builder,
      base_array_info=base_array_info)

    self.original_arrays = cctbxdat.original_arrays
    if data_block_name is not None:
      return cctbxdat.miller_arrays[data_block_name]
    else:
      return cctbxdat.miller_arrays

  def as_miller_arrays(self, data_block_name=None,
                       crystal_symmetry=None,
                       force_symmetry=False,
                       merge_equivalents=True,
                       base_array_info=None,
                       anomalous=None):
    if base_array_info is None:
      base_array_info = miller.array_info(
        source=self.file_path, source_type="cif")
    if data_block_name is not None:
      arrays = list(self.build_miller_arrays(
        data_block_name=data_block_name,
        base_array_info=base_array_info).values())
    else:
      arrays = flat_list([
        list(arrays.values()) for arrays in
        self.build_miller_arrays(base_array_info=base_array_info).values()])
    other_symmetry=crystal_symmetry
    for i in range(len(arrays)):
      if crystal_symmetry is not None:
        crystal_symmetry_from_file = arrays[i].crystal_symmetry()
        crystal_symmetry = crystal_symmetry_from_file.join_symmetry(
          other_symmetry=other_symmetry,
          force=force_symmetry)
        arrays[i] = arrays[i].customized_copy(
          crystal_symmetry=crystal_symmetry, info=arrays[i].info())
      if anomalous is not None:
        arrays[i] = arrays[i].customized_copy(
          anomalous_flag=anomalous, info=arrays[i].info())
    return arrays

  def as_original_arrays(self):
    return self.original_arrays

fast_reader = reader # XXX backward compatibility 2010-08-25

def atom_type_cif_loop(xray_structure, format="mmcif"):
  format = format.lower()
  assert format in ("corecif", "mmcif")
  if format == "mmcif": separator = '.'
  else: separator = '_'

  sources = {
    "it1992": "International Tables Volume C Table 6.1.1.4 (pp. 500-502)",
    "wk1995": "Waasmaier & Kirfel (1995), Acta Cryst. A51, 416-431",
  }
  inelastic_references = {
    "henke" : "Henke, Gullikson and Davis, At. Data and Nucl. Data Tables, 1993, 54, 2",
    "sasaki" : "Sasaki, KEK Report, 1989, 88-14, 1",
  }

  scattering_type_registry = xray_structure.scattering_type_registry()
  unique_gaussians = scattering_type_registry.unique_gaussians_as_list()
  max_n_gaussians = max([gaussian.n_terms() for gaussian in unique_gaussians])
  max_n_gaussians = max(max_n_gaussians, 4) # Need for compliance with mmcif_pdbx_v50
  # _atom_type_* loop
  header = ['_atom_type%ssymbol' %separator,
            '_atom_type%sscat_dispersion_real' %separator,
            '_atom_type%sscat_dispersion_imag' %separator]
  header.extend(['_atom_type%sscat_Cromer_Mann_a%i' %(separator, i+1)
                 for i in range(max_n_gaussians)])
  header.extend(['_atom_type%sscat_Cromer_Mann_b%i' %(separator, i+1)
                 for i in range(max_n_gaussians)])
  header.extend(['_atom_type%sscat_Cromer_Mann_c' %separator,
                 '_atom_type%sscat_source' %separator,
                 '_atom_type%sscat_dispersion_source' %separator])
  atom_type_loop = model.loop(header=header)
  gaussian_dict = scattering_type_registry.as_type_gaussian_dict()
  scattering_type_registry = xray_structure.scattering_type_registry()
  params = xray_structure.scattering_type_registry_params
  fp_fdp_table = {}
  for sc in xray_structure.scatterers():
    fp_fdp_table.setdefault(sc.scattering_type, (sc.fp, sc.fdp))
  disp_source = inelastic_references.get(
    xray_structure.inelastic_form_factors_source)
  # custom?
  if disp_source is None:
    disp_source = xray_structure.inelastic_form_factors_source
  if disp_source is None:
    disp_source = "."
  for atom_type, gaussian in six.iteritems(scattering_type_registry.as_type_gaussian_dict()):
    scat_source = sources.get(params.table)
    if params.custom_dict and atom_type in params.custom_dict:
      scat_source = "Custom %i-Gaussian" %gaussian.n_terms()
    elif scat_source is None:
      scat_source = """\
%i-Gaussian fit: Grosse-Kunstleve RW, Sauter NK, Adams PD:
Newsletter of the IUCr Commission on Crystallographic Computing 2004, 3, 22-31."""
      scat_source = scat_source %gaussian.n_terms()
    if disp_source == ".":
      fp, fdp = ".", "."
    else:
      fp, fdp = fp_fdp_table[atom_type]
      fp = "%.5f" %fp
      fdp = "%.5f" %fdp
    row = [atom_type, fp, fdp]
    #gaussian = gaussian_dict[sc.scattering_type]
    gaussian_a = ["%.5f" %a for a in gaussian.array_of_a()]
    gaussian_b = ["%.5f" %a for a in gaussian.array_of_b()]
    gaussian_a.extend(["."]*(max_n_gaussians-gaussian.n_terms()))
    gaussian_b.extend(["."]*(max_n_gaussians-gaussian.n_terms()))
    row.extend(gaussian_a + gaussian_b)
    row.extend([gaussian.c(), scat_source, disp_source])
    atom_type_loop.add_row(row)

  return atom_type_loop


def miller_indices_as_cif_loop(indices, prefix='_refln_'):
    refln_loop = model.loop(header=(
      '%sindex_h' %prefix, '%sindex_k' %prefix, '%sindex_l' %prefix))
    for hkl in indices:
      refln_loop.add_row(hkl)
    return refln_loop


class miller_arrays_as_cif_block():

  def __init__(self, array, array_type=None,
               column_name=None, column_names=None,
               miller_index_prefix='_refln',
               format="mmcif"):
    wformat = format.lower()
    assert wformat in ("corecif", "mmcif")
    if wformat == "mmcif":
      separator = '.'
    else:
      separator = '_'
    self.cif_block = array.crystal_symmetry().as_cif_block(format=format)
    self.prefix = miller_index_prefix + separator
    self.indices = array.indices().deep_copy()
    self.refln_loop = None
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
                 array.phases(deg=True).data().as_string()]
    elif array.is_hendrickson_lattman_array():
      if column_names is None:
        column_names = [self.prefix+'HL_%s_iso' %abcd for abcd in 'ABCD']
      else: assert len(column_names) == 4
      data = [d.as_string() for d in array.data().as_abcd()]
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
      if array.anomalous_flag():
        if ((array.sigmas() is not None and len(column_names) == 4) or
            (array.sigmas() is None and len(column_names) == 2)):
          data = []
          asu, matches = array.match_bijvoet_mates()
          for anomalous_sign in ("+", "-"):
            sel = matches.pairs_hemisphere_selection(anomalous_sign)
            sel.extend(matches.singles_hemisphere_selection(anomalous_sign))
            if (anomalous_sign == "+"):
              indices = asu.indices().select(sel)
              hemisphere_column_names = column_names[:len(column_names)//2]
            else:
              indices = -asu.indices().select(sel)
              hemisphere_column_names = column_names[len(column_names)//2:]
            hemisphere_data = asu.data().select(sel)
            hemisphere_array = miller.array(miller.set(
              array.crystal_symmetry(), indices), hemisphere_data)
            if array.sigmas() is not None:
              hemisphere_array.set_sigmas(asu.sigmas().select(sel))
            if self.refln_loop is None:
              # then this is the first array to be added to the loop,
              # hack so we don't have both hemispheres of indices
              self.indices = indices
            self.add_miller_array(
              hemisphere_array, column_names=hemisphere_column_names)
          return
      if array.sigmas() is not None and len(column_names) == 2:
        data.append(array.sigmas().as_string())
    if not (self.indices.size() == array.indices().size() and
            self.indices.all_eq(array.indices())):
      from cctbx.miller import match_indices
      other_indices = array.indices().deep_copy()
      match = match_indices(self.indices, other_indices)
      if match.singles(0).size():
        # array is missing some reflections indices that already appear in the loop
        # therefore pad the data with '?' values
        other_indices.extend(self.indices.select(match.single_selection(0)))
        for d in data:
          d.extend(flex.std_string(['?']*(other_indices.size() - d.size())))
        for d in data:
          assert d.size() == other_indices.size()
        match = match_indices(self.indices, other_indices)
      if match.singles(1).size():
        # this array contains some reflections that are not already present in the
        # cif loop, therefore need to add rows of '?' values
        single_indices = other_indices.select(match.single_selection(1))
        self.indices.extend(single_indices)
        n_data_columns = len(self.refln_loop) - 3
        for hkl in single_indices:
          row = list(hkl) + ['?'] * n_data_columns
          self.refln_loop.add_row(row)
        match = match_indices(self.indices, other_indices)

      match = match_indices(self.indices, other_indices)
      perm = match.permutation()
      data = [d.select(perm) for d in data]

    if self.refln_loop is None:
      self.refln_loop = miller_indices_as_cif_loop(self.indices, prefix=self.prefix)
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
    self.original_arrays = OrderedDict()
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
    wavelengths = {}
    for key, block in cif_model.items():
      if data_block_name is not None and key != data_block_name: continue
      for builder in data_structure_builders:
        if builder == builders.crystal_structure_builder:
          if '_atom_site_fract_x' in block or '_atom_site_Cartn_x' in block:
            self.xray_structures.setdefault(key, builder(block).structure)
        elif builder == builders.miller_array_builder:
          block_wavelengths = builders.get_wavelengths(block)
          if (block_wavelengths is not None):
            wavelengths = block_wavelengths
          if base_array_info is not None:
            base_array_info = base_array_info.customized_copy(labels=[key])
          if ( '_refln_index_h' in block or '_refln.index_h' in block or
               '_diffrn_refln' in block
               ):
            b = builder(block, base_array_info=base_array_info,
                wavelengths=wavelengths)
            self.miller_arrays.setdefault( key, b.arrays())
            self.original_arrays.setdefault( key, b.origarrays())


# This defines the order that categories will appear in the CIF file
category_order = [
  '_cell',
  '_space_group',
  '_space_group_symop',
  '_symmetry',
  '_computing',
  '_software',
  '_em_software',
  '_citation',
  '_citation_author',
  '_reflns',
  '_reflns_shell',
  '_refine',
  '_refine_ls_restr',
  '_refine_ls_shell',
  '_pdbx_refine_tls',
  '_pdbx_refine_tls_group',
  '_struct_asym',
  '_struct_conf_type',
  '_struct_conf',
  '_struct_conn',
  '_struct_sheet',
  '_struct_sheet_order',
  '_struct_sheet_range',
  '_pdbx_struct_sheet_hbond',
  '_struct_ncs_ens',
  '_struct_ncs_dom',
  '_struct_ncs_ens_gen',
  '_struct_ncs_dom_lim',
  '_struct_ncs_oper',
  '_refine_ls_restr_ncs',
  '_entity',
  '_entity_poly',
  '_entity_poly_seq',
  '_atom_type',
  '_chem_comp',
  '_chem_comp_atom',
  '_atom_site',
]

def category_sort_function(key):
  key_category = key.split('.')[0]
  try:
    return category_order.index(key_category)
  except ValueError as e:
    # any categories we don't know about will end up at the end of the file
    return len(category_order)
