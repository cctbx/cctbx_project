from __future__ import absolute_import, division, print_function
from cctbx import adptbx, crystal, miller, sgtbx, uctbx, xray
from cctbx.array_family import flex
import iotbx.cif
from iotbx.cif import model
from libtbx.utils import Sorry
from libtbx.containers import OrderedDict, OrderedSet
import warnings
from six import string_types
from six.moves import range
import six, re
from six.moves import zip

# Refer to https://www.iucr.org/__data/iucr/cifdic_html/2/cif_mm.dic/index.html for definitons
# of elements and columns in a CIF file

class CifBuilderError(Sorry):
  __module__ = Exception.__module__

def CifBuilderWarning(message):
  warnings.showwarning(message, UserWarning, 'CifBuilderWarning', '')

class cif_model_builder(object):

  def __init__(self, cif_object=None):
    self._model = cif_object
    if self._model is None:
      self._model = model.cif()
    self._current_block = None
    self._current_save = None

  def add_data_block(self, data_block_heading):
    self._current_block = model.block()
    if data_block_heading.lower() == 'global_':
      block_name = data_block_heading
    else:
      block_name = data_block_heading[data_block_heading.find('_')+1:]
    self._model[block_name] = self._current_block

  def add_loop(self, header, columns):
    if self._current_save is not None:
      block = self._current_save
    else:
      block = self._current_block
    loop = model.loop()
    assert len(header) == len(columns)
    n_columns = len(columns)
    for i in range(n_columns):
      loop[header[i]] = columns[i]
    if loop.name() not in block.loops.keys():
      block.add_loop(loop)
    else:
      raise Sorry("Loop containing tags %s appears repeated" %', '.join(loop.keys()))

  def add_data_item(self, key, value):
    if self._current_save is not None:
      if key not in self._current_save:
        self._current_save[key] = value
      else:
        raise Sorry("Data item %s received multiple values" % key)
    elif self._current_block is not None:
      if key not in self._current_block:
        self._current_block[key] = value
      else:
        raise Sorry("Data item %s received multiple values" % key)
    else: # support for global_ blocks in non-strict mode
      pass

  def start_save_frame(self, save_frame_heading):
    assert self._current_save is None
    self._current_save = model.save()
    save_name = save_frame_heading[save_frame_heading.find('_')+1:]
    self._current_block[save_name] = self._current_save

  def end_save_frame(self):
    self._current_save = None

  def model(self):
    return self._model


class builder_base(object):

  __equivalents__ = {
    '_space_group_symop_operation_xyz': ('_symmetry_equiv_pos_as_xyz',
                                         '_space_group_symop.operation_xyz',
                                         '_symmetry_equiv.pos_as_xyz'),
    '_space_group_symop_id': ('_symmetry_equiv_pos_site_id',
                              '_space_group_symop.id',
                              '_symmetry_equiv.id'),
    '_space_group_name_Hall': ('_symmetry_space_group_name_Hall',
                               '_space_group.name_Hall',
                               '_symmetry.space_group_name_Hall'),
    '_space_group_name_H-M_alt': ('_symmetry_space_group_name_H-M',
                                  '_space_group.name_H-M_alt',
                                  '_symmetry.space_group_name_H-M'),
    '_space_group_IT_number': ('_symmetry_Int_Tables_number',
                                 '_symmetry.Int_Tables_number'
                                 '_space_group.IT_number'),
    '_cell_length_a': ('_cell.length_a',),
    '_cell_length_b': ('_cell.length_b',),
    '_cell_length_c': ('_cell.length_c',),
    '_cell_angle_alpha': ('_cell.angle_alpha',),
    '_cell_angle_beta': ('_cell.angle_beta',),
    '_cell_angle_gamma': ('_cell.angle_gamma',),
    '_cell_volume': ('_cell.volume',),
    '_refln_index_h': ('_refln.index_h',),
    '_refln_index_k': ('_refln.index_k',),
    '_refln_index_l': ('_refln.index_l',),
  }

  def get_cif_item(self, key, default=None):
    value = self.cif_block.get(key)
    if value is not None: return value
    for equiv in self.__equivalents__.get(key, []):
      value = self.cif_block.get(equiv)
      if value is not None: return value
    return default


class crystal_symmetry_builder(builder_base):

  def __init__(self, cif_block, strict=False):
    # The order of priority for determining space group is:
    #   sym_ops, hall symbol, H-M symbol, space group number
    self.cif_block = cif_block
    sym_ops = self.get_cif_item('_space_group_symop_operation_xyz')
    sym_op_ids = self.get_cif_item('_space_group_symop_id')
    space_group = None
    space_group_from_ops = None
    space_group_from_other = None
    # Symmetry from operators
    if sym_ops is not None:
      if isinstance(sym_ops, string_types):
        sym_ops = flex.std_string([sym_ops])
      if sym_op_ids is not None:
        if isinstance(sym_op_ids, string_types):
          sym_op_ids = flex.std_string([sym_op_ids])
        assert len(sym_op_ids) == len(sym_ops)
      self.sym_ops = {}
      space_group_from_ops = sgtbx.space_group()
      if isinstance(sym_ops, string_types): sym_ops = [sym_ops]
      for i, op in enumerate(sym_ops):
        try:
          s = sgtbx.rt_mx(op)
        except RuntimeError as e:
          str_e = str(e)
          if "Parse error: " in str_e:
            raise CifBuilderError("Error interpreting symmetry operator: %s" %(
              str_e.split("Parse error: ")[-1]))
          else:
            raise
        if sym_op_ids is None:
          sym_op_id = i+1
        else:
          try:
            sym_op_id = int(sym_op_ids[i])
          except ValueError as e:
            raise CifBuilderError("Error interpreting symmetry operator id: %s" %(
              str(e)))
        self.sym_ops[sym_op_id] = s
        space_group_from_ops.expand_smx(s)
    # Symmetry from other
    hall_symbol = self.get_cif_item('_space_group_name_Hall')
    hm_symbol = self.get_cif_item('_space_group_name_H-M_alt')
    sg_number = self.get_cif_item('_space_group_IT_number')
    if sg_number not in (None, '?'):
      try: space_group_from_other = sgtbx.space_group_info(number=sg_number).group()
      except Exception: pass
    if hall_symbol not in (None, '?'):
      try: space_group_from_other = sgtbx.space_group(hall_symbol)
      except Exception: pass
    if hm_symbol not in (None, '?'):
      try: space_group_from_other = sgtbx.space_group_info(symbol=hm_symbol).group()
      except Exception: pass
    # Check consistency.
    # Use space group equivalence operation in cctbx.sgtbx.space_group

    if (space_group_from_other is not None) and (
       space_group_from_ops is not None) and (not
         space_group_from_other == space_group_from_ops):
      ops1 = [o.as_xyz() for o in space_group_from_other.all_ops()]
      ops2 = [o.as_xyz() for o in space_group_from_ops.all_ops()]
      ops1.sort()
      ops2.sort()
      msg1 = "\n"+"\n".join(ops1)+"\n"
      msg2 = "\n"+"\n".join(ops2)
      raise CifBuilderError(
          "Inconsistent symmetry information found:%s ---vs---%s"%(msg1, msg2))
    for sg in [space_group_from_other, space_group_from_ops]:
      if(sg is not None):
        space_group = sg
        break
    if (space_group is None and strict):
      raise CifBuilderError(
        "No symmetry instructions could be extracted from the cif block")
    #
    items = [self.get_cif_item("_cell_length_"+s) for s in "abc"]
    for i, item in enumerate(items):
      if isinstance(item, flex.std_string):
        raise CifBuilderError(
          "Data item _cell_length_%s cannot be declared in a looped list"
          %("abc"[i]))
    for s in ["alpha", "beta", "gamma"]:
      item = self.get_cif_item("_cell_angle_"+s)
      if isinstance(item, flex.std_string):
        raise CifBuilderError(
          "Data item _cell_angle_%s cannot be declared in a looped list" %s)
      if (item == "?"):
        item = "90" # enumeration default for angles is 90 degrees
      items.append(item)
    ic = items.count(None)
    ic_question = items.count('?')
    if ic == 6 or ic_question > 0:
      if (strict):
        raise CifBuilderError(
          "Unit cell parameters not found in the cif file")
      unit_cell = None
    elif (ic == 0):
      try:
        vals = [float_from_string(s) for s in items]
      except ValueError:
        raise CifBuilderError("Invalid unit cell parameters are given")
      try:
        unit_cell = uctbx.unit_cell(vals)
      except RuntimeError as e:
        if "cctbx Error: Unit cell" in str(e):
          raise CifBuilderError(e)
        else:
          raise
    elif (space_group is not None):
      unit_cell = uctbx.infer_unit_cell_from_symmetry(
        [float_from_string(s) for s in items if s is not None], space_group)
    else:
      raise CifBuilderError(
        "Not all unit cell parameters are given in the cif file")
    if unit_cell is not None and space_group is not None:
      if not space_group.is_compatible_unit_cell(unit_cell):
        # try primitive setting
        space_group_input = space_group
        space_group = space_group.info().primitive_setting().group()
        if not space_group.is_compatible_unit_cell(unit_cell):
          raise CifBuilderError(
            "Space group is incompatible with unit cell parameters:\n" + \
            "  Space group: %s\n" %space_group_input.info() + \
            "  Unit cell: %s" %unit_cell)
    self.crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,
                                             space_group=space_group)

class crystal_structure_builder(crystal_symmetry_builder):

  def __init__(self, cif_block):
    # XXX To do: interpret _atom_site_refinement_flags
    crystal_symmetry_builder.__init__(self, cif_block, strict=True)
    atom_sites_frac = [
      as_double_or_none_if_all_question_marks(
        _, column_name='_atom_site_fract_%s' %axis)
      for _, axis in [(cif_block.get('_atom_site_fract_%s' %axis), axis)
                      for axis in ('x','y','z')]]
    if atom_sites_frac.count(None) == 3:
      atom_sites_cart = [as_double_or_none_if_all_question_marks(
        _, column_name='_atom_site_Cartn_%s' %axis)
                         for _ in [cif_block.get('_atom_site_Cartn_%s' %axis)
                                   for axis in ('x','y','z')]]
      if atom_sites_cart.count(None) != 0:
        raise CifBuilderError("No atomic coordinates could be found")
      atom_sites_cart = flex.vec3_double(*atom_sites_cart)
      # XXX do we need to take account of _atom_sites_Cartn_tran_matrix_ ?
      atom_sites_frac = self.crystal_symmetry.unit_cell().fractionalize(
        atom_sites_cart)
    else:
      if atom_sites_frac.count(None) != 0:
        raise CifBuilderError("No atomic coordinates could be found")
      atom_sites_frac = flex.vec3_double(*atom_sites_frac)
    labels = cif_block.get('_atom_site_label')
    type_symbol = cif_block.get('_atom_site_type_symbol')
    if type_symbol:
      type_symbol = flex.std_string(
        s.replace('0+', '').replace('0-', '') for s in type_symbol)
    U_iso_or_equiv = flex_double_else_none(
      cif_block.get('_atom_site_U_iso_or_equiv',
      cif_block.get('_atom_site_U_equiv_geom_mean')))
    if U_iso_or_equiv is None:
      B_iso_or_equiv = flex_double_else_none(
        cif_block.get('_atom_site_B_iso_or_equiv',
        cif_block.get('_atom_site_B_equiv_geom_mean')))
    adp_type = cif_block.get('_atom_site_adp_type')
    occupancy = flex_double_else_none(cif_block.get('_atom_site_occupancy'))
    scatterers = flex.xray_scatterer()
    atom_site_aniso_label = flex_std_string_else_none(
      cif_block.get('_atom_site_aniso_label'))
    if atom_site_aniso_label is not None:
      atom_site_aniso_label = atom_site_aniso_label
      adps = [cif_block.get('_atom_site_aniso_U_%i' %i)
              for i in (11,22,33,12,13,23)]
      have_Bs = False
      if adps.count(None) > 0:
        adps = [cif_block.get('_atom_site_aniso_B_%i' %i)
                for i in (11,22,33,12,13,23)]
        have_Bs = True
      if adps.count(None) == 6:
        adps = None
      elif adps.count(None) > 0:
        CifBuilderError("Some ADP items are missing")
      else:
        sel = None
        for adp in adps:
          f = (adp == "?")
          if (sel is None): sel = f
          else:             sel &= f
        sel = ~sel
        atom_site_aniso_label = atom_site_aniso_label.select(sel)
        try:
          adps = [flex.double(adp.select(sel)) for adp in adps]
        except ValueError as e:
          raise CifBuilderError("Error interpreting ADPs: " + str(e))
        adps = flex.sym_mat3_double(*adps)
    for i in range(len(atom_sites_frac)):
      kwds = {}
      if labels is not None:
        kwds.setdefault('label', str(labels[i]))
      if type_symbol is not None:
        kwds.setdefault('scattering_type', str(type_symbol[i]))
      if (atom_site_aniso_label is not None
          and adps is not None
          and labels is not None
          and labels[i] in atom_site_aniso_label):
        adp = adps[flex.first_index(atom_site_aniso_label, labels[i])]
        if have_Bs: adp = adptbx.b_as_u(adp)
        kwds.setdefault('u', adptbx.u_cif_as_u_star(
          self.crystal_symmetry.unit_cell(), adp))
      elif U_iso_or_equiv is not None:
        kwds.setdefault('u', float_from_string(U_iso_or_equiv[i]))
      elif B_iso_or_equiv is not None:
        kwds.setdefault('b', float_from_string(B_iso_or_equiv[i]))
      if occupancy is not None:
        kwds.setdefault('occupancy', float_from_string(occupancy[i]))
      scatterers.append(xray.scatterer(**kwds))
    scatterers.set_sites(atom_sites_frac)

    wvl_str = self.get_cif_item('_diffrn_radiation_wavelength')
    if not isinstance(wvl_str, str) and wvl_str is not None:
      wvl_str = wvl_str[0]
    wavelength = float_from_string(wvl_str) if (wvl_str and wvl_str!='?') else None

    self.structure = xray.structure(crystal_symmetry=self.crystal_symmetry,
                                    scatterers=scatterers,
                                    wavelength=wavelength)


class miller_array_builder(crystal_symmetry_builder):
  # Changes to this class should pass regression tests:
  # cctbx_project\mmtbx\regression\tst_cif_as_mtz_wavelengths.py
  # cctbx_project\iotbx\cif\tests\tst_lex_parse_build.py
  # phenix_regression\cif_as_mtz\tst_cif_as_mtz.py

  observation_types = {
    # known types of column data to be tagged as either amplitudes or intensities as per
    # https://www.iucr.org/__data/iucr/cifdic_html/2/cif_mm.dic/index.html
    '_refln.F_squared': xray.intensity(),
    '_refln_F_squared': xray.intensity(),
    '_refln.intensity': xray.intensity(),
    '_refln.I(+)': xray.intensity(),
    '_refln.I(-)': xray.intensity(),
    '_refln.F_calc': xray.amplitude(),
    '_refln.F_meas': xray.amplitude(),
    '_refln.FP': xray.amplitude(),
    '_refln.F-obs': xray.amplitude(),
    '_refln.Fobs': xray.amplitude(),
    '_refln.F-calc': xray.amplitude(),
    '_refln.Fcalc': xray.amplitude(),
    '_refln.pdbx_F_': xray.amplitude(),
    '_refln.pdbx_I_': xray.intensity(),
    '_refln.pdbx_anom_difference': xray.amplitude(),
  }

  def guess_observationtype(self, labl):
    for okey in self.observation_types.keys():
      if labl.startswith(okey):
        return self.observation_types[okey]
    return None

  def __init__(self, cif_block, base_array_info=None, wavelengths=None):
    crystal_symmetry_builder.__init__(self, cif_block)
    self._arrays = OrderedDict()
    self._origarrays = OrderedDict() # used for presenting raw data tables in HKLviewer
    basearraylabels = []
    if base_array_info is not None:
      self.crystal_symmetry = self.crystal_symmetry.join_symmetry(
        other_symmetry=base_array_info.crystal_symmetry_from_file,
      force=True)
      if base_array_info.labels:
        basearraylabels = base_array_info.labels
    if (wavelengths is None):
      wavelengths = {}
    if base_array_info is None:
      base_array_info = miller.array_info(source_type="cif")
    refln_containing_loops = self.get_miller_indices_containing_loops()
    for self.indices, refln_loop in refln_containing_loops:
      self.wavelength_id_array = None
      self.crystal_id_array = None
      self.scale_group_array = None
      wavelength_ids = [None]
      crystal_ids = [None]
      scale_groups = [None]
      for key, value in six.iteritems(refln_loop):
        # Get wavelength_ids, crystal_id, scale_group_code columns for selecting data of other
        # columns in self.get_selection() used by self.flex_std_string_as_miller_array()
        if (key.endswith('wavelength_id') or
            key.endswith('crystal_id') or
            key.endswith('scale_group_code')):
          data = as_int_or_none_if_all_question_marks(value, column_name=key)
          if data is None:
            continue
          counts = data.counts()
          if key.endswith('wavelength_id'):
            wavelength_ids = list(counts.keys())
          if len(counts) == 1: continue
          array = miller.array(
            miller.set(self.crystal_symmetry, self.indices).auto_anomalous(), data)
          if key.endswith('wavelength_id'):
            self.wavelength_id_array = array
            wavelength_ids = list(counts.keys())
          elif key.endswith('crystal_id'):
            self.crystal_id_array = array
            crystal_ids = list(counts.keys())
          elif key.endswith('scale_group_code'):
            self.scale_group_array = array
            scale_groups = list(counts.keys())
      labelsuffix = []
      wavelbl = []
      cryslbl = []
      scalegrplbl = []
      self._origarrays["HKLs"] = self.indices
      alllabels = list(sorted(refln_loop.keys()))
      remaininglabls = alllabels[:] # deep copy the list
      # Parse labels matching cif column conventions
      # https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/refln.html
      # and extract groups of labels or just single columns.
      # Groups corresponds to the map coefficients, phase and amplitudes,
      # amplitudes or intensities with sigmas and hendrickson-lattman columns.
      phaseamplabls, remaininglabls = self.get_phase_amplitude_labels(remaininglabls)
      mapcoefflabls, remaininglabls = self.get_mapcoefficient_labels(remaininglabls)
      HLcoefflabls, remaininglabls = self.get_HL_labels(remaininglabls)
      data_sig_obstype_labls, remaininglabls = self.get_FSigF_ISigI_labels(remaininglabls)
      for w_id in wavelength_ids:
        for crys_id in crystal_ids:
          for scale_group in scale_groups:
            # If reflection data files contain more than one crystal, wavelength or scalegroup
            # then add their id(s) as a suffix to data labels computed below. Needed for avoiding
            # ambuguity but avoid when not needed to make labels more human readable!
            if (len(wavelength_ids) > 1 or len(wavelengths) > 1) and w_id is not None:
              wavelbl = ["wavelength_id=%i" %w_id]
            if len(crystal_ids) > 1 and crys_id is not None:
              cryslbl = ["crystal_id=%i" %crys_id]
            if len(scale_groups) > 1 and scale_group is not None:
              scalegrplbl = ["scale_group_code=%i" %scale_group]
            labelsuffix = scalegrplbl + cryslbl + wavelbl
            jlablsufx = ""
            if len(labelsuffix):
              jlablsufx = "," + ",".join(labelsuffix)
            for mapcoefflabl in mapcoefflabls:
              A_array = refln_loop[ mapcoefflabl[0] ]
              B_array = refln_loop[ mapcoefflabl[1] ]
              # deselect any ? marks in the two arrays, assuming both A and B have the same ? marks
              selection = self.get_selection( A_array, wavelength_id=w_id,
               crystal_id=crys_id, scale_group_code=scale_group)
              A_array = A_array.select(selection)
              B_array = B_array.select(selection)
              # form the miller array with map coefficients
              data = flex.complex_double( flex.double(A_array), flex.double(B_array) )
              millarr = miller.array( miller.set(self.crystal_symmetry,
                                             self.indices.select(selection)
                                            ).auto_anomalous(), data)
              # millarr will be None for column data not matching w_id,crys_id,scale_group values
              if millarr is None: continue
              labl = basearraylabels + mapcoefflabl + labelsuffix
              millarr.set_info(base_array_info.customized_copy(labels= labl ,
                                                wavelength=wavelengths.get(w_id, None)))
              self._arrays[mapcoefflabl[0] + jlablsufx ] = millarr
            for phaseamplabl in phaseamplabls:
              amplitudestrarray = refln_loop[ phaseamplabl[0] ]
              phasestrarray = refln_loop[ phaseamplabl[1] ]
              millarr = self.flex_std_string_as_miller_array(
                amplitudestrarray, wavelength_id=w_id, crystal_id=crys_id,
                scale_group_code=scale_group)
              phasesmillarr = self.flex_std_string_as_miller_array(
                phasestrarray, wavelength_id=w_id, crystal_id=crys_id,
                scale_group_code=scale_group)
              # millarr will be None for column data not matching w_id,crys_id,scale_group values
              if millarr is None or phasesmillarr is None: continue
              if(not check_array_sizes(millarr, phasesmillarr, phaseamplabl[0], phaseamplabl[1])):
                continue
              phases = as_flex_double(phasesmillarr, phaseamplabl[1])
              millarr = millarr.phase_transfer(phases, deg=True)
              labl = basearraylabels + phaseamplabl + labelsuffix
              millarr.set_info(base_array_info.customized_copy(labels= labl ,
                                                wavelength=wavelengths.get(w_id, None)))
              self._arrays[phaseamplabl[0] +jlablsufx ] = millarr
            for datlabl,siglabl,otype in data_sig_obstype_labls:
              datastrarray = refln_loop[datlabl]
              millarr = self.flex_std_string_as_miller_array(
                datastrarray, wavelength_id=w_id, crystal_id=crys_id,
                scale_group_code=scale_group)
              # millarr will be None for column data not matching w_id,crys_id,scale_group values
              if millarr is None: continue
              millarr = as_flex_double(millarr, datlabl)
              datsiglabl = [datlabl]
              if siglabl:
                sigmasstrarray = refln_loop[siglabl]
                sigmas = self.flex_std_string_as_miller_array(
                  sigmasstrarray, wavelength_id=w_id, crystal_id=crys_id,
                  scale_group_code=scale_group)
                sigmas = as_flex_double(sigmas, siglabl)
                if check_array_sizes(millarr, sigmas, datlabl, siglabl):
                  millarr.set_sigmas(sigmas.data())
                  datsiglabl = [datlabl, siglabl]
                else:
                  sigmas.set_info(base_array_info.customized_copy(labels= [siglabl],
                                                    wavelength=wavelengths.get(w_id, None)))
                  self._arrays[ siglabl +jlablsufx ] = sigmas

              datsiglabl = basearraylabels + datsiglabl + labelsuffix
              millarr.set_info(base_array_info.customized_copy(labels= datsiglabl,
                                                wavelength=wavelengths.get(w_id, None)))
              if otype is not None:
                millarr.set_observation_type(otype)
              self._arrays[ datlabl +jlablsufx ] = millarr
            for hl_labels in HLcoefflabls:
              hl_values = [ cif_block.get(hl_key) for hl_key in hl_labels ]
              if hl_values.count(None) == 0:
                selection = self.get_selection(
                  hl_values[0], wavelength_id=w_id,
                  crystal_id=crys_id, scale_group_code=scale_group)
                hl_values = [as_double_or_none_if_all_question_marks(
                  hl.select(selection), column_name=lab)
                              for hl, lab in zip(hl_values, hl_labels)]
                # hl_values will be None for column data not matching w_id,crys_id,scale_group values
                if hl_values == [None,None,None,None]: continue
                millarr = miller.array(miller.set(
                  self.crystal_symmetry, self.indices.select(selection)
                  ).auto_anomalous(), flex.hendrickson_lattman(*hl_values))
                hlabels = basearraylabels + hl_labels + labelsuffix
                millarr.set_info(base_array_info.customized_copy(labels= hlabels,
                                                  wavelength=wavelengths.get(w_id, None)))
                self._arrays[ hl_labels[0] + jlablsufx ] = millarr
            # pick up remaining columns if any that weren't identified above
            for label in alllabels:
              if "index_" in label:
                continue
              datastrarray = refln_loop[label]
              if label in remaininglabls:
                labels = basearraylabels + [label]  + labelsuffix
                lablsufx = jlablsufx
                millarr = self.flex_std_string_as_miller_array(
                  datastrarray, wavelength_id=w_id, crystal_id=crys_id,
                  scale_group_code=scale_group)
                # millarr will be None for column data not matching w_id,crys_id,scale_group values
                if (label.endswith('wavelength_id') or
                 label.endswith('crystal_id') or # get full array if any of these labels, not just subsets
                 label.endswith('scale_group_code')):
                  millarr = self.flex_std_string_as_miller_array(
                    datastrarray, wavelength_id=None, crystal_id=None,
                    scale_group_code=None)
                  lablsufx = ""
                  labels = basearraylabels + [label]
                if millarr is None: continue
                otype = self.guess_observationtype(label)
                if otype is not None:
                  millarr.set_observation_type(otype)
                millarr.set_info(base_array_info.customized_copy(labels= labels,
                                                  wavelength=wavelengths.get(w_id, None)))
                self._arrays[ label + lablsufx ] = millarr
              origarr = self.flex_std_string_as_miller_array(
                  datastrarray, wavelength_id=w_id, crystal_id=crys_id,
                  scale_group_code=scale_group,
                  allowNaNs=True)
              newlabel = label.replace("_refln.", "")
              newlabel2 = newlabel.replace("_refln_", "")
              if origarr: # want only genuine miller arrays
                self._origarrays[newlabel2 + jlablsufx ] = origarr.data()
    # Convert any groups of I+,I-,SigI+,SigI- (or amplitudes) arrays into anomalous arrays
    # i.e. both friedel mates in the same array
    for key, array in six.iteritems(self._arrays.copy()):
      plus_key = ""
      if '_minus' in key:
        minus_key = key
        plus_key = key.replace('_minus', '_plus')
      elif '-' in key:
        minus_key = key
        plus_key = key.replace('-', '+')
      elif '_plus' in key:
        plus_key = key
        minus_key = key.replace('_plus', '_minus')
      elif '+' in key:
        plus_key = key
        minus_key = key.replace('+', '-')
      if plus_key in self._arrays and minus_key in self._arrays:
        plus_array = self._arrays.pop(plus_key)
        minus_array = self._arrays.pop(minus_key)
        minus_array = minus_array.customized_copy(
          indices=-minus_array.indices()).set_info(minus_array.info())
        array = plus_array.concatenate(
          minus_array, assert_is_similar_symmetry=False)
        array = array.customized_copy(anomalous_flag=True)
        array.set_info(minus_array.info().customized_copy(
          labels=list(
            OrderedSet(plus_array.info().labels+minus_array.info().labels))))
        array.set_observation_type(plus_array.observation_type())
        self._arrays.setdefault(key, array)
    if len(self._arrays) == 0:
      raise CifBuilderError("No reflection data present in cif block")
    # Sort the ordered dictionary to resemble the order of columns in the cif file
    # This is to avoid any F_meas arrays accidentally being put adjacent to
    # pdbx_anom_difference arrays in the self._arrays OrderedDict. Otherwise these
    # arrays may unintentionally be combined into a reconstructed anomalous amplitude
    # array when saving as an mtz file due to a problem in the iotbx/mtz module.
    # See http://phenix-online.org/pipermail/cctbxbb/2021-March/002289.html
    arrlstord = []
    arrlst = list(self._arrays)
    for arr in arrlst:
      for i,k in enumerate(refln_loop.keys()):
        if arr.split(",")[0] == k:
          arrlstord.append((arr, i))
    # arrlstord must have the same keys as in the self._arrays dictionary
    assert sorted(arrlst) == sorted([ e[0] for e in arrlstord] )
    sortarrlst = sorted(arrlstord, key=lambda arrord: arrord[1])
    self._ordarrays = OrderedDict()
    for sortkey,i in sortarrlst:
      self._ordarrays.setdefault(sortkey, self._arrays[sortkey])
    self._arrays = self._ordarrays


  def get_HL_labels(self, keys):
    lstkeys = list(keys) # cast into list if not a list
    HLquads = []
    alllabels = " ".join(lstkeys)
    """ Hendrickson-Lattmann labels could look like: 'HLAM', 'HLBM', 'HLCM', 'HLDM'
    or like 'HLanomA', 'HLanomB', 'HLanomC', 'HLanomD'
    Use a regular expression to group them accordingly
    """
    allmatches = re.findall(r"(\S*(HL(\S*)[abcdABCD](\S*)))", alllabels )
    HLtagslst = list(set([ (e[2], e[3]) for e in allmatches ]))
    usedkeys = []
    for m in HLtagslst:
      hllist = []
      for hm in allmatches:
        if m==(hm[2], hm[3]):
          hllist.append((hm[0], hm[1]))
      if len(hllist) == 4:
        HLquads.append([ e[0] for e in hllist])
        for e in hllist:
          usedkeys.append(e[0])
    remainingkeys = []
    for e in lstkeys:
      if e not in usedkeys:
        remainingkeys.append(e)
    return HLquads, remainingkeys


  def get_mapcoefficient_labels(self, keys):
    # extract map coeffficients labels from list of cif column labels
    # e.g. ( _refln.A_calc_au _refln.B_calc_au ) , ( _refln.A_calc _refln.B_calc )
    lstkeys = list(keys) # cast into list if not a list
    remainingkeys = lstkeys[:] # deep copy the list
    alllabels = " ".join(lstkeys)
    mapcoefflabels = []
    A_matches = re.findall(r"( (\s*_refln[\._]A_)(\S*) )", alllabels, re.VERBOSE ) # [('_refln.PHWT', '_refln.PH', 'WT'), ('_refln.PHDELWT', '_refln.PH', 'DELWT')]
    for label in lstkeys:
      for m in A_matches:
        Blabel = m[1].replace("A_","B_") + m[2]
        if Blabel == label:
          mapcoefflabels.append([ m[0], label])
          remainingkeys.remove(m[0])
          remainingkeys.remove(label)
    return mapcoefflabels, remainingkeys


  def get_phase_amplitude_labels(self, keys):
    # extract phase and amplitudes labels from list of cif column labels
    # e.g. ( _refln.F_calc _refln.phase_calc ) , ( _refln.FC_ALL _refln.PHIC_ALL ), ( _refln.FWT _refln.PHWT )
    lstkeys = list(keys) # cast into list if not a list
    remainingkeys = lstkeys[:] # deep copy the list
    alllabels = " ".join(lstkeys)
    phase_amplitudelabels = []
    PHmatches = re.findall(r"((\S*PH)([^I]\S*))", alllabels ) # [('_refln.PHWT', '_refln.PH', 'WT'), ('_refln.PHDELWT', '_refln.PH', 'DELWT')]
    for label in lstkeys:
      for m in PHmatches:
        PFlabel = m[1].replace("PH","F") + m[2]
        Flabel = m[1].replace("PH","") + m[2]
        if Flabel == label or PFlabel == label:
          phase_amplitudelabels.append([ label, m[0]])
          remainingkeys.remove(label)
          remainingkeys.remove(m[0])
    alllabels = " ".join(remainingkeys)
    PHImatches = re.findall(r"((\S*PHI)(\S*))", alllabels ) # [('_refln.PHIC', '_refln.PHI', 'C'), ('_refln.PHIC_ALL', '_refln.PHI', 'C_ALL')]
    for label in lstkeys:
      for m in PHImatches:
        PFlabel = m[1].replace("PHI","F") + m[2]
        Flabel = m[1].replace("PHI","") + m[2]
        if Flabel == label or PFlabel == label:
          phase_amplitudelabels.append([ label, m[0]])
          remainingkeys.remove(label)
          remainingkeys.remove(m[0])
    alllabels = " ".join(remainingkeys)
    PHDELmatches = re.findall(r"(((\S*)PH)([^I]\S*(WT)))", alllabels ) # [('_refln.PHDELWT', '_refln.PH', '_refln.', 'DELWT', 'WT')]
    for label in lstkeys:
      for m in PHDELmatches:
        Flabel = m[2] + m[3].replace("WT","FWT")
        if Flabel == label:
          phase_amplitudelabels.append([ label, m[0]])
          remainingkeys.remove(label)
          remainingkeys.remove(m[0])
    alllabels = " ".join(remainingkeys)
    phase_matches = re.findall(r"((\S*[\._])phase(\S*))", alllabels ) # [('_refln.phase_calc', '_refln.', '')]
    for label in lstkeys:
      for m in phase_matches:
        phaselabel = m[0]
        Flabl = m[1] + m[2]
        Flabel = m[1] + "F" + m[2]
        Faulabel = m[1] + "F" + m[2] + "_au"
        if Flabl in label or Flabel in label or Faulabel in label: # in case of _refln.F_calc_au and _refln.phase_calc
          if label in remainingkeys and m[0] in remainingkeys: # in case
            if (Flabel + "_sigma_au") in remainingkeys or (Flabel + "_sigma") in remainingkeys:
              continue # give priority to F_meas, F_meas_sigma or  F_meas_au, F_meas_sigma_au
            phase_amplitudelabels.append([ label, m[0]])
            remainingkeys.remove(label)
            remainingkeys.remove(m[0])
    return phase_amplitudelabels, remainingkeys


  def get_FSigF_ISigI_labels(self, keys):
    # extract amplitudea, sigmas or intensitiy, sigmas labels from list of cif column labels
    # e.g. ( _refln.F_meas_sigma_au _refln.F_meas), ( _refln.intensity_sigma _refln.intensity ) ,
    # ( _refln.pdbx_I_plus_sigma _refln.pdbx_I_plus )
    lstkeys = list(keys) # cast into list if not a list
    remainingkeys = lstkeys[:] # deep copy the list
    alllabels = " ".join(lstkeys)
    labelpairs = []
    sigma_matches = re.findall(r"((\S*[\._])SIG(\S*))", alllabels ) # catch label pairs like F(+),SIGF(+)
    for label in lstkeys:
      for m in sigma_matches:
        FIlabel = m[1] + m[2]
        if FIlabel == label:
          labelpairs.append([ label, m[0], self.guess_observationtype(label)])
          remainingkeys.remove(label)
          remainingkeys.remove(m[0])
    alllabels = " ".join(remainingkeys)
    sigma_matches = re.findall(r"((\S*)_sigma(_*\S*))", alllabels ) # [('_refln.F_meas_sigma_au', '_refln.F_meas', '_au'), ('_refln.intensity_sigma', '_refln.intensity', ''), ('_refln.pdbx_I_plus_sigma', '_refln.pdbx_I_plus', '')]
    for label in lstkeys:
      for m in sigma_matches:
        FIlabel = m[1] + m[2]
        if FIlabel == label:
          labelpairs.append([ label, m[0], self.guess_observationtype(label)])
          remainingkeys.remove(label)
          remainingkeys.remove(m[0])
    alllabels = " ".join(remainingkeys)
    # catch generic meas and sigma labels
    anymeas_matches = re.findall(r"((\S*)_meas(\S*))", alllabels ) + re.findall(r"((\S*)_calc(\S*))", alllabels )
    anysigma_matches = re.findall(r"((\S*)_sigma(\S*))", alllabels )
    for mmatch in anymeas_matches:
      for smatch in anysigma_matches:
        if mmatch[1]==smatch[1] and mmatch[2]==smatch[2]:
          remainingkeys.remove(mmatch[0])
          if smatch[0] in remainingkeys: # in case of say F_squared_calc, F_squared_meas, F_squared_sigma all being present
            remainingkeys.remove(smatch[0])
            labelpairs.append([ mmatch[0], smatch[0], self.guess_observationtype(mmatch[0])])
          else:
            labelpairs.append([ mmatch[0], None, self.guess_observationtype(mmatch[0])])
    return labelpairs, remainingkeys


  def get_miller_indices_containing_loops(self):
    loops = []
    for loop in self.cif_block.loops.values():
      for key in loop.keys():
        if 'index_h' not in key: continue
        hkl_str = [loop.get(key.replace('index_h', 'index_%s' %i)) for i in 'hkl']
        if hkl_str.count(None) > 0:
          raise CifBuilderError(
            "Miller indices missing from current CIF block (%s)"
            %key.replace('index_h', 'index_%s' %'hkl'[hkl_str.index(None)]))
        hkl_int = []
        for i,h_str in enumerate(hkl_str):
          try:
            h_int = flex.int(h_str)
          except ValueError as e:
            raise CifBuilderError(
              "Invalid item for Miller index %s: %s" % ("HKL"[i], str(e)))
          hkl_int.append(h_int)
        indices = flex.miller_index(*hkl_int)
        loops.append((indices, loop))
        break
    return loops


  def get_selection(self, value,
                    wavelength_id=None,
                    crystal_id=None,
                    scale_group_code=None,
                    allowNaNs = False):
    if allowNaNs:
      selection = flex.bool(value.size(), True)
    else:
      selection = ~((value == '.') | (value == '?'))
    if self.wavelength_id_array is not None and wavelength_id is not None:
      selection &= (self.wavelength_id_array.data() == wavelength_id)
    if self.crystal_id_array is not None and crystal_id is not None:
      selection &= (self.crystal_id_array.data() == crystal_id)
    if self.scale_group_array is not None and scale_group_code is not None:
      selection &= (self.scale_group_array.data() == scale_group_code)
    return selection


  def flex_std_string_as_miller_array(self, value,
                                      wavelength_id=None,
                                      crystal_id=None,
                                      scale_group_code=None,
                                      allowNaNs = False):
    # Create a miller_array object of only the data and indices matching the
    # wavelength_id, crystal_id and scale_group_code submitted or full array if these are None
    selection = self.get_selection(
      value, wavelength_id=wavelength_id,
      crystal_id=crystal_id, scale_group_code=scale_group_code,
      allowNaNs=allowNaNs)
    data = value.select(selection)
    try:
      data = flex.int(data)
      indices = self.indices.select(selection)
    except ValueError:
      try:
        data = flex.double(data)
        indices = self.indices.select(selection)
      except ValueError:
        # if flex.std_string return all values including '.' and '?'
        data = value
        indices = self.indices
    if data.size() == 0: return None
    return miller.array(
      miller.set(self.crystal_symmetry, indices).auto_anomalous(), data)


  def arrays(self):
    return self._arrays


  def origarrays(self):
    """
    return dictionary of raw data found in cif file cast into flex.double arrays
    or just string arrays as a fall back.
    """
    return self._origarrays


def as_flex_double(array, key):
  if isinstance(array.data(), flex.double):
    return array
  elif isinstance(array.data(), flex.int):
    return array.customized_copy(
      data=array.data().as_double()).set_info(array.info())
  else:
    try:
      flex.double(array.data())
    except ValueError as e:
      e_str = str(e)
      if e_str.startswith("Invalid floating-point value: "):
        i = e_str.find(":") + 2
        raise CifBuilderError("Invalid floating-point value for %s: %s"
                              %(key, e_str[i:].strip()))
      else:
        raise CifBuilderError(e_str)

def check_array_sizes(array1, array2, key1, key2):
  if array1.size() != array2.size():
    msg = "Miller arrays '%s' and '%s' are of different sizes" %(key1, key2)
    CifBuilderWarning(message=msg)
    return False
  return True

def none_if_all_question_marks_or_period(cif_block_item):
  if (cif_block_item is None): return None
  result = cif_block_item
  if (result.all_eq("?")): return None
  elif (result.all_eq(".")): return None
  return result

def as_int_or_none_if_all_question_marks(cif_block_item, column_name=None):
  strings = none_if_all_question_marks_or_period(cif_block_item)
  if (strings is None): return None
  try:
    return flex.int(strings)
  except ValueError as e:
    # better error message if column_name is given
    e_str = str(e)
    if column_name is not None and e_str.startswith(
      "Invalid integer value: "):
      i = e_str.find(":") + 2
      raise CifBuilderError("Invalid integer value for %s: %s"
                            %(column_name, e_str[i:].strip()))
    else:
      raise CifBuilderError(e_str)

def as_double_or_none_if_all_question_marks(cif_block_item, column_name=None):
  strings = none_if_all_question_marks_or_period(cif_block_item)
  if (strings is None): return None
  try:
    return flex.double(strings)
  except ValueError as e:
    # better error message if column_name is given
    e_str = str(e)
    if column_name is not None and e_str.startswith(
      "Invalid floating-point value: "):
      i = e_str.find(":") + 2
      raise CifBuilderError("Invalid floating-point value for %s: %s"
                            %(column_name, e_str[i:].strip()))
    else:
      raise CifBuilderError(e_str)

def flex_double(flex_std_string):
  try:
    return flex.double(flex_std_string)
  except ValueError as e:
    raise CifBuilderError(str(e))

def flex_double_else_none(cif_block_item):
  strings = none_if_all_question_marks_or_period(cif_block_item)
  if (strings is None): return None
  try:
    return flex.double(strings)
  except ValueError:
    pass
  return None

def flex_std_string_else_none(cif_block_item):
  if isinstance(cif_block_item, flex.std_string):
    return cif_block_item
  else:
    return None

def float_from_string(string):
  """a cif string may be quoted,
and have an optional esd in brackets associated with it"""
  if isinstance(string, float):
    return string
  return float(string.strip('\'').strip('"').split('(')[0])

def get_wavelengths(cif_block):
  for loop in cif_block.loops.values():
    for key in loop.keys():
      if ("_diffrn_radiation_wavelength." in key):
        wavelength_ids = loop.get("_diffrn_radiation_wavelength.id")
        wavelength_strs = loop.get("_diffrn_radiation_wavelength.wavelength")
        if (not None in [wavelength_ids, wavelength_strs]):
          wl_ = {}
          for wavelength_id,wavelength in zip(wavelength_ids,wavelength_strs):
            try :
              wl_id = int(wavelength_id)
              wl_[int(wavelength_id)] = float(wavelength)
            except ValueError :
              pass
          return wl_
        else :
          return None
  wavelength_id = cif_block.get("_diffrn_radiation_wavelength.id")
  wavelength_str = cif_block.get("_diffrn_radiation_wavelength.wavelength")
  if (not None in [wavelength_id, wavelength_str]):
    try :
      wl_id = int(wavelength_id)
      return { int(wavelength_id) : float(wavelength_str) }
    except ValueError :
      pass
  return None
