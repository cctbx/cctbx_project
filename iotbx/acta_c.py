from cctbx import xray
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.array_family import flex
from PyCifRW.CifFile import CifFile

def cif_float(s):
  return float(s.split("(")[0])

def cif_as_miller_array(file_name):
  cif_file = CifFile(datasource=file_name, strict=0)
  assert len(cif_file.keys()) == 1
  cif_section = cif_file[cif_file.keys()[0]]
  unit_cell = uctbx.unit_cell([cif_float(cif_section[param])
    for param in [
      "_cell_length_a","_cell_length_b","_cell_length_c",
      "_cell_angle_alpha","_cell_angle_beta","_cell_angle_gamma"]])
  space_group = sgtbx.space_group()
  for xyz in cif_section["_symmetry_equiv_pos_as_xyz"]:
    space_group.expand_smx(xyz)
  crystal_symmetry = crystal.symmetry(
    unit_cell=unit_cell,
    space_group=space_group)
  for squared in ["_squared", ""]:
    try: refln_items = [cif_section[param] for param in [
      "_refln_index_h",
      "_refln_index_k",
      "_refln_index_l",
      "_refln_F%s_meas" % squared,
      "_refln_F%s_sigma" % squared]]
    except KeyError: refln_items = None
    else: break
  assert refln_items is not None
  refln_items.append(cif_section.get("_refln_observed_status", None))
  if (refln_items[-1] is None):
    refln_items[-1] = ["o"]*len(refln_items[0])
  for refln_item in refln_items[1:]:
    assert len(refln_item) == len(refln_items[0])
  indices = flex.miller_index()
  data = flex.double()
  sigmas = flex.double()
  for h,k,l,meas,sigma,status in zip(*refln_items):
    if (status != "o"): continue
    indices.append([int(i) for i in (h,k,l)])
    data.append(float(meas))
    sigmas.append(float(sigma))
  miller_array = miller.array(
    miller_set=miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=indices).auto_anomalous(),
    data=data,
    sigmas=sigmas)
  if (squared == ""):
    miller_array.set_observation_type_xray_amplitude()
  else:
    miller_array.set_observation_type_xray_intensity()
  miller_array.set_info(miller.array_info(
    source=file_name,
    labels=[s%squared for s in ["F%s_meas", "F%s_sigma"]]))
  return miller_array

def get_cif_column_string(
      cif_block, item_name, size=None, mandatory=False):
  result = flex.std_string()
  if (not cif_block.has_key(item_name)):
    assert not mandatory
    if (size is not None):
      result.resize(size)
  else:
    for value in cif_block[item_name]:
      result.append(value)
  if (size is not None):
    assert result.size() == size
  return result

def get_cif_column_double(
      cif_block, item_name, size=None, default_value=0, mandatory=False):
  result = flex.double()
  if (not cif_block.has_key(item_name)):
    assert not mandatory
    if (size is not None):
      result.resize(size, default_value)
  else:
    for value in cif_block[item_name]:
      result.append(cif_float(value))
  if (size is not None):
    assert result.size() == size
  return result

def get_cif_column_vec3_double(
      cif_block, item_names, size=None, mandatory=False):
  result = flex.vec3_double()
  if (not cif_block.has_key(item_names[0])):
    assert not mandatory
    if (size is not None):
      result.resize(size)
  else:
    for values in zip(*[cif_block[item_name] for item_name in item_names]):
      result.append([cif_float(value) for value in values])
  if (size is not None):
    assert result.size() == size
  return result

def cif_as_xray_structure(file_name, data_block_name):
  cif_file = CifFile(datasource=file_name, strict=0)
  cif_block = cif_file[data_block_name]
  unit_cell = uctbx.unit_cell([cif_float(cif_block[param])
    for param in [
      "_cell_length_a","_cell_length_b","_cell_length_c",
      "_cell_angle_alpha","_cell_angle_beta","_cell_angle_gamma"]])
  space_group_info = sgtbx.space_group_info(
    symbol=cif_block["_symmetry_space_group_name_H-M"])
  crystal_symmetry = crystal.symmetry(
    unit_cell=unit_cell,
    space_group_info=space_group_info)
  sites = get_cif_column_vec3_double(
    cif_block,
    ["_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"],
    mandatory=True)
  labels = get_cif_column_string(
    cif_block,
    "_atom_site_label",
    size=sites.size())
  type_symbols = get_cif_column_string(
    cif_block,
    "_atom_site_type_symbol",
    size=sites.size())
  u_iso_or_equiv = get_cif_column_double(
    cif_block,
    "_atom_site_U_iso_or_equiv",
    size=sites.size())
  occupancies = get_cif_column_double(
    cif_block,
    "_atom_site_occupancy",
    default_value=1,
    size=sites.size())
  result = xray.structure(crystal_symmetry=crystal_symmetry)
  for label,scattering_type,site,u_iso,occupancy in zip(
        labels, type_symbols, sites, u_iso_or_equiv, occupancies):
    if (scattering_type == ""):
      scattering_type = None
    scatterer = xray.scatterer(
      label=label,
      site=site,
      scattering_type=scattering_type,
      u=u_iso,
      occupancy=occupancy)
    result.add_scatterer(scatterer)
  return result
