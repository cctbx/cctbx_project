import iotbx.pdb
from iotbx.pdb import parser
from iotbx.pdb import cryst1_interpretation
from iotbx.pdb import scattering_type_interpr
from cctbx import xray
from cctbx import crystal
from cctbx.eltbx import tiny_pse
from cctbx import adptbx
from cctbx.array_family import flex
from cStringIO import StringIO

def xray_structure_as_pdb_file(self, remark=None, remarks=[],
                                     fractional_coordinates=False,
                                     res_name=None,
                                     connect=None):
  if (remark is not None):
    remarks.insert(0, remark)
  s = StringIO()
  for remark in remarks:
    print >> s, "REMARK", remark
  print >> s, "REMARK Number of scatterers:", self.scatterers().size()
  print >> s, "REMARK At special positions:", \
    self.special_position_indices().size()
  if (fractional_coordinates):
    print >> s, "REMARK Fractional coordinates"
  else:
    print >> s, "REMARK Cartesian coordinates"
  print >> s, iotbx.pdb.format_cryst1_record(crystal_symmetry=self)
  print >> s, iotbx.pdb.format_scale_records(unit_cell=self.unit_cell())
  if (res_name is not None):
    res_name_i = res_name.upper()
  serial = 0
  for scatterer in self.scatterers():
    serial += 1
    if (scatterer.anisotropic_flag):
      u_cart = adptbx.u_star_as_u_cart(self.unit_cell(), scatterer.u_star)
      u = adptbx.u_cart_as_u_iso(u_cart)
    else:
      u_cart = None
      u = scatterer.u_iso
    xyz = scatterer.site
    if (not fractional_coordinates):
      xyz = self.unit_cell().orthogonalize(xyz)
    label = scatterer.label.upper()
    element_symbol = scatterer.element_symbol()
    if (element_symbol is None): element_symbol = "Q"
    assert len(element_symbol) in (1,2)
    if (res_name is None):
      res_name_i = label[:3]
    print >> s, iotbx.pdb.format_atom_record(
      serial=serial,
      name=label[:4],
      resName=res_name_i,
      resSeq=serial,
      site=xyz,
      occupancy=scatterer.occupancy,
      tempFactor=adptbx.u_as_b(u),
      element=element_symbol.upper())
    if (u_cart is not None):
      print >> s, iotbx.pdb.format_anisou_record(
        serial=serial,
        name=label[:4],
        resName=res_name_i,
        resSeq=serial,
        u_cart=u_cart,
        element=element_symbol.upper())
  if (connect is not None):
    assert len(connect) == self.scatterers().size()
    i = 0
    for bonds in connect:
      i += 1
      l = "CONNECT%5d" % i
      for bond in bonds:
        l += "%5d" % (bond+1)
      print >> s, l
  print >> s, "END"
  return s.getvalue()

xray.structure.as_pdb_file = xray_structure_as_pdb_file

def from_pdb(file_name=None, file_iterator=None, pdb_records=None,
             crystal_symmetry=None, force_symmetry=False,
             ignore_atom_element_q=True,
             scan_atom_element_columns=True,
             fractional_coordinates=False,
             min_distance_sym_equiv=0.5,
             keep_scatterer_pdb_records=False):
  assert [file_name, file_iterator, pdb_records].count(None) == 2
  if (pdb_records is None):
    if (file_iterator is None):
      file_iterator = open(file_name)
    pdb_records = parser.collect_records(
      raw_records=file_iterator,
      ignore_master=True)
  cryst1_symmetry = None
  for record in pdb_records:
    if (record.record_name.startswith("CRYST1")):
      try:
        cryst1_symmetry = cryst1_interpretation.crystal_symmetry(
          cryst1_record=record)
      except: pass
      break
  if (cryst1_symmetry is not None):
    crystal_symmetry = cryst1_symmetry.join_symmetry(
      other_symmetry=crystal_symmetry,
      force=force_symmetry)
  assert crystal_symmetry is not None, "Unknown crystal symmetry."
  assert crystal_symmetry.unit_cell() is not None, "Unknown unit cell."
  assert crystal_symmetry.space_group_info() is not None,"Unknown space group."
  if (scan_atom_element_columns):
    have_useful_atom_element_columns = (
      scattering_type_interpr.scan_atom_element_columns(pdb_records)
        .n_uninterpretable == 0)
  else:
    have_useful_atom_element_columns = None
  structure = xray.structure(
    special_position_settings=crystal.special_position_settings(
      crystal_symmetry=crystal_symmetry,
      min_distance_sym_equiv=min_distance_sym_equiv))
  scatterer = None
  scatterer_pdb_record = None
  if (keep_scatterer_pdb_records):
    scatterer_pdb_records = []
  prev_record = None
  for record in pdb_records:
    if (record.record_name in ("ATOM", "HETATM")):
      if (scatterer is not None):
        try: structure.add_scatterer(scatterer)
        except RuntimeError, e:
          raise RuntimeError("%s%s" % (prev_record.error_prefix(), str(e)))
        if (keep_scatterer_pdb_records):
          scatterer_pdb_records.append(scatterer_pdb_record)
      if (ignore_atom_element_q
          and (record.element == " Q" or record.name[1] == "Q")):
        scatterer = None
        scatterer_pdb_record = None
      else:
        if (fractional_coordinates):
          site=record.coordinates
        else:
          site=structure.unit_cell().fractionalize(record.coordinates)
        scatterer = xray.scatterer(
          label="%s%03d" % (record.name, record.serial),
          site=site,
          b=record.tempFactor,
          occupancy=record.occupancy,
          scattering_type=scattering_type_interpr.from_pdb_atom_record(
            record, have_useful_atom_element_columns))
        scatterer_pdb_record = record
    elif (record.record_name == "SIGATM"):
      if (not prev_record.record_name in ("ATOM", "HETATM")):
        record.raise_FormatError(
          "SIGATM record without preceeding ATOM or HETATM record.")
      if (prev_record.raw[7:27] != record.raw[7:27]):
        record.raise_FormatError(
          "SIGATM record does not match preceeding %s record."
            % prev_record.record_name)
    elif (record.record_name == "ANISOU"):
      if (not prev_record.record_name in ("ATOM", "HETATM", "SIGATM")):
        record.raise_FormatError(
          "ANISOU record without preceeding ATOM or HETATM record.")
      if (prev_record.raw[7:27] != record.raw[7:27]):
        record.raise_FormatError(
          "ANISOU record does not match preceeding %s record."
            % prev_record.record_name)
      if (scatterer is not None):
        scatterer = scatterer.customized_copy(
          u=adptbx.u_cart_as_u_star(structure.unit_cell(), record.Ucart))
    elif (scatterer is not None):
      try: structure.add_scatterer(scatterer)
      except RuntimeError, e:
        raise RuntimeError("%s%s" % (prev_record.error_prefix(), str(e)))
      if (keep_scatterer_pdb_records):
        scatterer_pdb_records.append(scatterer_pdb_record)
      scatterer = None
      scatterer_pdb_record = None
    prev_record = record
  if (scatterer is not None):
    try: structure.add_scatterer(scatterer)
    except RuntimeError, e:
      raise RuntimeError("%s%s" % (prev_record.error_prefix(), str(e)))
    if (keep_scatterer_pdb_records):
      scatterer_pdb_records.append(scatterer_pdb_record)
  if (keep_scatterer_pdb_records):
    structure.scatterer_pdb_records = scatterer_pdb_records
  return structure
