from iotbx.pdb import parser
from iotbx.pdb import residue_info
from cctbx import xray
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import adptbx
from cctbx.array_family import flex
from scitbx.python_utils.math_utils import iround
from cStringIO import StringIO

def xray_structure_as_pdb_file(self, remark=None, remarks=[],
                                     fractional_coordinates=00000,
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
  # CRYST1
  #  7 - 15       Real(9.3)      a             a (Angstroms).
  # 16 - 24       Real(9.3)      b             b (Angstroms).
  # 25 - 33       Real(9.3)      c             c (Angstroms).
  # 34 - 40       Real(7.2)      alpha         alpha (degrees).
  # 41 - 47       Real(7.2)      beta          beta (degrees).
  # 48 - 54       Real(7.2)      gamma         gamma (degrees).
  # 56 - 66       LString        sGroup        Space group.
  # 67 - 70       Integer        z             Z value.
  print >> s, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %s" % (
    self.unit_cell().parameters()
    + (str(self.space_group_info()).replace(" ", ""),))
  # ATOM
  #  7 - 11  Integer       serial        Atom serial number.
  # 13 - 16  Atom          name          Atom name.
  # 17       Character     altLoc        Alternate location indicator.
  # 18 - 20  Residue name  resName       Residue name.
  # 22       Character     chainID       Chain identifier.
  # 23 - 26  Integer       resSeq        Residue sequence number.
  # 27       AChar         iCode         Code for insertion of residues.
  # 31 - 38  Real(8.3)     x             Orthogonal coordinates for X in
  #                                      Angstroms.
  # 39 - 46  Real(8.3)     y             Orthogonal coordinates for Y in
  #                                      Angstroms.
  # 47 - 54  Real(8.3)     z             Orthogonal coordinates for Z in
  #                                      Angstroms.
  # 55 - 60  Real(6.2)     occupancy     Occupancy.
  # 61 - 66  Real(6.2)     tempFactor    Temperature factor.
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
    assert serial < 100000
    atom_07_27 = ("%5d %-4s %-3s  %4d" % (
      serial, label[:4], label[:3], serial%10000),)
    print >> s, "ATOM  %s    %8.3f%8.3f%8.3f%6.2f%6.2f" % (
      atom_07_27 + xyz + (scatterer.occupancy, adptbx.u_as_b(u)))
    if (u_cart is not None):
      print >> s, "ANISOU%s%7d%7d%7d%7d%7d%7d" % (
        atom_07_27 + tuple([iround(x*10000) for x in u_cart]))
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

def from_pdb(file_name=None, file_iterator=None,
             crystal_symmetry=None, force_symmetry=00000,
             fractional_coordinates=00000):
  assert [file_name, file_iterator].count(None) == 1
  if (file_iterator is None):
    file_iterator = open(file_name)
  pdb_records = parser.collect_records(raw_records=file_iterator)
  unit_cell = None
  space_group_info = None
  for record in pdb_records:
    if (record.record_name.startswith("CRYST1")):
      try: unit_cell = uctbx.unit_cell(record.ucparams)
      except: pass
      try: space_group_info = sgtbx.space_group_info(symbol=record.sGroup)
      except: pass
      break
  crystal_symmetry=crystal.symmetry(
    unit_cell=unit_cell,
    space_group_info=space_group_info).join_symmetry(
    other_symmetry=crystal_symmetry,
    force=force_symmetry)
  assert crystal_symmetry.unit_cell() is not None
  assert crystal_symmetry.space_group_info() is not None
  structure = xray.structure(
    special_position_settings=crystal.special_position_settings(
      crystal_symmetry=crystal_symmetry))
  scatterers = flex.xray_scatterer()
  prev_record = None
  for record in pdb_records:
    if (record.record_name in ("ATOM", "HETATM")):
      if (fractional_coordinates):
        site=record.coordinates
      else:
        site=structure.unit_cell().fractionalize(record.coordinates)
      try:
        caasf = residue_info.get(
          residue_name=record.resName,
          atom_name=record.name).scattering_label
      except KeyError:
        caasf = ""
      scatterer = xray.scatterer(
        label="%s%03d" % (record.name, record.serial),
        site=site,
        b=record.tempFactor,
        occupancy=record.occupancy,
        caasf=caasf)
      scatterers.append(scatterer)
    elif (record.record_name == "ANISOU"):
      if (not prev_record.record_name in ("ATOM", "HETATM")):
        record.raise_FormatError(
          "ANISOU record without preceeding ATOM or HETATM record.")
      if (prev_record.raw[7:27] != record.raw[7:27]):
        record.raise_FormatError(
          "ANISOU record does not match preceeding %s record."
            % prev_record.record_name)
      scatterers[-1] = scatterer.copy(
        u=adptbx.u_cart_as_u_star(structure.unit_cell(), record.Ucart))
    prev_record = record
  structure.add_scatterers(scatterers)
  return structure
