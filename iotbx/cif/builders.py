from cctbx import adptbx, crystal, miller, sgtbx, uctbx, xray
from cctbx.array_family import flex
import iotbx.cif
from iotbx.cif import model

try:
  import PyCifRW
  import PyCifRW.CifFile
except ImportError:
  PyCifRW = None

class cif_model_builder:

  def __init__(self):
    self._model = model.cif()
    self._current_block = None
    self._current_save = None

  def add_data_block(self, data_block_heading):
    self._current_block = model.block()
    block_name = data_block_heading[data_block_heading.find('_')+1:]
    self._model[block_name] = self._current_block

  def add_loop(self, header, data):
    if self._current_save is not None:
      block = self._current_save
    else:
      block = self._current_block
    loop = model.loop()
    n_columns = len(header)
    assert len(data) % n_columns == 0, "Wrong number of data items for loop"
    if n_columns == 1:
      columns = [data]
    else:
      columns = iotbx.cif.ext.looped_data_as_columns(data, n_columns)
    for i in range(n_columns):
      loop[header[i]] = columns[i]
    block.add_loop(loop)

  def add_data_item(self, key, value):
    if self._current_save is not None:
      self._current_save[key] = value
    else:
      self._current_block[key] = value

  def start_save_frame(self, save_frame_heading):
    assert self._current_save is None
    self._current_save = model.save()
    save_name = save_frame_heading[save_frame_heading.find('_')+1:]
    self._current_block[save_name] = self._current_save

  def end_save_frame(self):
    self._current_save = None

  def model(self):
    return self._model

if PyCifRW is not None:
  class PyCifRW_model_builder:

    def __init__(self):
      self._model = PyCifRW.CifFile.CifFile()
      self._current_block = None

    def add_data_block(self, data_block_heading):
      if self._current_block is not None:
        # add previous block to the model
        # PyCifRW seems to copy the block object rather than pass by reference
        self._model[self._current_block_name] = self._current_block
      self._current_block = PyCifRW.CifFile.CifBlock()
      self._current_block_name = data_block_heading[data_block_heading.find('_')+1:]

    def add_loop(self, header, data):
      loop = PyCifRW.CifFile.CifLoopBlock(dimension=1)
      for key in header:
        loop[key] = []
      self._makeloop(loop, [data])
      self._current_block.insert_loop(loop)

    def _makeloop(self, loopstructure, itemlists):
      # helper function borrowed from PyCifRW YappsStarParser_1_1.py
      if itemlists[-1] == []: itemlists.pop(-1)
      if loopstructure.dimension == 1 and loopstructure.loops == []:
        storage_iter = loopstructure.fast_load_iter()
      else:
        storage_iter = loopstructure.load_iter()
      nowloop = loopstructure
      for datalist in itemlists:
        for datavalue in datalist:
          try:
            nowloop,target = storage_iter.next()
          except StopIteration:
            print "StopIter at %s/%s" % (datavalue,datalist)
            raise StopIteration
          target.append(datavalue)
        nowloop.popout = True
        nowloop,blank = storage_iter.next()
      return loopstructure

    def add_data_item(self, key, value):
      self._current_block[str(key)] = str(value)

    def model(self):
      if self._current_block is not None:
        # add last block to the model
        self._model[self._current_block_name] = self._current_block
      return self._model

  def start_save_frame(self, save_frame_heading):
    pass

  def end_save_frame(self):
    pass



# Build cctbx data structures from either iotbx cif model, or PyCifRW model

class crystal_symmetry_builder:

  def __init__(self, cif_block):
    # The order of priority for determining space group is:
    #   sym_ops, hall symbol, H-M symbol, space group number
    sym_ops = cif_block.get('_space_group_symop_operation_xyz',
              cif_block.get('_symmetry_equiv_pos_as_xyz'))
    sym_op_ids = cif_block.get('_space_group_symop_id',
                 cif_block.get('_symmetry_equiv_pos_site_id'))
    if sym_ops is not None:
      if sym_op_ids is not None:
        assert len(sym_op_ids) == len(sym_ops)
      self.sym_ops = {}
      space_group = sgtbx.space_group()
      if isinstance(sym_ops, basestring): sym_ops = [sym_ops]
      for i, op in enumerate(sym_ops):
        s = sgtbx.rt_mx(op)
        if sym_op_ids is None:
          sym_op_id = i+1
        else:
          sym_op_id = int(sym_op_ids[i])
        self.sym_ops[sym_op_id] = s
        space_group.expand_smx(s)
    else:
      hall_symbol = cif_block.get('_space_group_name_Hall',
                    cif_block.get('_symmetry_space_group_name_Hall'))
      hm_symbol = cif_block.get('_space_group_name_H-M_alt',
                  cif_block.get('_symmetry_space_group_name_H-M'))
      sg_number = cif_block.get('_space_group_symop_sg_id',
                  cif_block.get('_symmetry_Int_Tables_number'))
      if hall_symbol is not None:
        space_group = sgtbx.space_group(hall_symbol)
      elif hm_symbol is not None:
        space_group = sgtbx.space_group_info(symbol=hm_symbol).group()
      elif sg_number is not None:
        space_group = sgtbx.space_group_info(number=sg_number).group()
      else:
        raise RuntimeError("No symmetry instructions are present in the cif block")
    try:
      cell_params = [float_from_string(
        cif_block['_cell_length_%s' %dim]) for dim in ('a','b','c')]
      cell_params.extend([float_from_string(
        cif_block['_cell_angle_%s' %angle]) for angle in ('alpha','beta','gamma')])
      unit_cell = uctbx.unit_cell(cell_params)
    except KeyError:
      raise RuntimeError("Not all unit cell parameters are given in the cif file")
    self.crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,
                                             space_group=space_group)

class crystal_structure_builder(crystal_symmetry_builder):

  def __init__(self, cif_block):
    # XXX To do: interpret _atom_site_refinement_flags
    crystal_symmetry_builder.__init__(self, cif_block)
    atom_sites_frac = [cif_block.get('_atom_site_fract_%s' %axis)
                       for axis in ('x','y','z')]
    if atom_sites_frac.count(None) == 3:
      atom_sites_cart = [cif_block.get('_atom_site_Cartn_%s' %axis)
                         for axis in ('x','y','z')]
      assert atom_sites_cart.count(None) == 0
      atom_sites_cart = flex.vec3_double(
        flex.double(atom_sites_cart[0]),
        flex.double(atom_sites_cart[1]),
        flex.double(atom_sites_cart[2]))
      # XXX do we need to take account of _atom_sites_Cartn_tran_matrix_ ?
      atom_sites_frac = self.crystal_symmetry.unit_cell().fractionalize(
        atom_sites_cart)
    else:
      assert atom_sites_frac.count(None) == 0
      atom_sites_frac = flex.vec3_double(
        flex.double(flex.std_string(atom_sites_frac[0])),
        flex.double(flex.std_string(atom_sites_frac[1])),
        flex.double(flex.std_string(atom_sites_frac[2])))
    labels = cif_block.get('_atom_site_label')
    type_symbol = cif_block.get('_atom_site_type_symbol')
    U_iso_or_equiv = flex_double_else_none(
      cif_block.get('_atom_site_U_iso_or_equiv'))
    B_iso_or_equiv = flex_double_else_none(
      cif_block.get('_atom_site_B_iso_or_equiv'))
    adp_type = cif_block.get('_atom_site_adp_type')
    occupancy = flex_double_else_none(cif_block.get('_atom_site_occupancy'))
    scatterers = flex.xray_scatterer()
    atom_site_aniso_label = cif_block.get('_atom_site_aniso_label')
    if atom_site_aniso_label is not None:
      atom_site_aniso_label = flex.std_string(atom_site_aniso_label)
      adps = [cif_block.get('_atom_site_aniso_U_%i' %i)
              for i in (11,22,33,12,13,23)]
      have_Bs = False
      if adps.count(None) > 0:
        adps = [cif_block.get('_atom_site_aniso_B_%i' %i)
                for i in (11,22,33,12,13,23)]
        have_Bs = True
      assert adps.count(None) == 0
      adps = [flex.double(flex.std_string(adp)) for adp in adps]
      adps = flex.sym_mat3_double(*adps)
    for i in range(len(atom_sites_frac)):
      kwds = {}
      if labels is not None:
        kwds.setdefault('label', str(labels[i]))
      if type_symbol is not None:
        kwds.setdefault('scattering_type', str(type_symbol[i]))
      if (atom_site_aniso_label is not None
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

    self.structure = xray.structure(crystal_symmetry=self.crystal_symmetry,
                                    scatterers=scatterers)

class miller_array_builder(crystal_symmetry_builder):

  observation_types = {
    '_refln_F_squared': xray.intensity(),
    '_refln_intensity': xray.intensity(),
    '_refln_F': xray.amplitude(),
    '_refln_A': None,
  }

  def __init__(self, cif_block, base_array_info=None):
    crystal_symmetry_builder.__init__(self, cif_block)
    self._arrays = {}
    if base_array_info is None:
      base_array_info = miller.array_info(source_type="cif")
    hkl = [flex.int(flex.std_string(cif_block.get('_refln_index_%s' %i)))
           for i in ('h','k','l')]
    indices = flex.miller_index(*hkl)

    phase_calc = cif_block.get('_refln_phase_calc')
    phase_meas = cif_block.get('_refln_phase_meas')
    for prefix in self.observation_types.keys():
      sigmas = cif_block.get('%s_sigma' %prefix)
      obs_type = self.observation_types[prefix]
      if sigmas is not None:
        sigmas = flex.double(flex.std_string(sigmas))
      for array_type in ('meas', 'calc'):
        label = '_'.join((prefix, array_type))
        if prefix == '_refln_A':
          data = [cif_block.get(label),
                  cif_block.get('_'.join((prefix.replace('A', 'B'), array_type)))]
          if data.count(None) == 0:
            data = flex.complex_double(
              flex.double(flex.std_string(data[0])),
              flex.double(flex.std_string(data[1])))
          else: continue
        else:
          data = cif_block.get(label)
          if data is not None:
            data = flex.double(flex.std_string(data))
          else: continue
        array = miller.array(
          miller.set(self.crystal_symmetry, indices).auto_anomalous(),
          data, sigmas)
        if obs_type is not xray.intensity():
          if array_type == 'calc' and phase_calc is not None:
            array = array.phase_transfer(
              flex.double(flex.std_string(phase_calc)), deg=True)
          elif array_type == 'meas' and phase_meas is not None:
            array = array.phase_transfer(
              flex.double(flex.std_string(phase_meas)), deg=True)
        array.set_observation_type(obs_type)
        labels = [label]
        if sigmas is not None:
          labels.append('%s_sigma' %prefix)
        array.set_info(base_array_info.customized_copy(labels=labels))
        self._arrays.setdefault(label, array)
    refln_loop = cif_block.get('_refln')
    if refln_loop is not None:
      for key, value in refln_loop.iteritems():
        if (key not in self._arrays and '_index_' not in key
            and not key.endswith('sigma')
            and not key.endswith('calc')
            and not key.endswith('meas')):
          try:
            data = flex.int(flex.std_string(value))
          except RuntimeError:
            try:
              data = flex.double(flex.std_string(value))
            except RuntimeError:
              data = value
          array = miller.array(
            miller.set(self.crystal_symmetry, indices).auto_anomalous(), data)
          self._arrays.setdefault(key, array)

    if len(self._arrays) == 0:
      raise RuntimeError("No reflection data present in cif block")

  def arrays(self):
    return self._arrays

def flex_double_else_none(strings):
  if strings is None: return None
  try:
    return flex.double(flex.std_string(strings))
  except RuntimeError, e:
    if 'bad lexical cast' in str(e): return None
    else: raise e

def float_from_string(string):
  """a cif string may be quoted,
and have an optional esd in brackets associated with it"""
  if isinstance(string, float):
    return string
  return float(string.strip('\'').strip('"').split('(')[0])
