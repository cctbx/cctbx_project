from cctbx import xray
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import adptbx
from phenix.foundation.Structure import PDBClass # XXX move to iotbx
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
  i = 0
  for scatterer in self.scatterers():
    i += 1
    assert not scatterer.anisotropic_flag, "Not implemented." # XXX
    xyz = scatterer.site
    if (not fractional_coordinates):
      xyz = self.unit_cell().orthogonalize(xyz)
    label = scatterer.label.upper()
    print >> s, "ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f" % ((
      i,
      label[:4],
      label[:3],
      i,
      ) + xyz + (
      scatterer.occupancy,
      adptbx.u_as_b(scatterer.u_iso)))
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

def from_pdb(file_name, crystal_symmetry=None, force_symmetry=00000):
  pdb_file = PDBClass(file_name)
  try:
    cryst1 = pdb_file.GetUnitCellRecord()
  except:
    cryst1 = None
  try:
    unit_cell = uctbx.unit_cell(cryst1.ucparams)
  except:
    unit_cell = None
  try:
    space_group_info = sgtbx.space_group_info(symbol=cryst1.sGroup)
  except:
    space_group_info = None
  else:
    if (space_group_info.group() is None):
      space_group_info = None
  crystal_symmetry=crystal.symmetry(
    unit_cell=unit_cell,
    space_group_info=space_group_info).join_symmetry(
    other_symmetry=crystal_symmetry,
    force=force_symmetry)
  structure = xray.structure(
    special_position_settings=crystal.special_position_settings(
      crystal_symmetry=crystal_symmetry))
  for atom in pdb_file.GetAtomList():
    scatterer = xray.scatterer(
      label="%s%03d" % (atom.name, atom.serial),
      site=structure.unit_cell().fractionalize(atom.coordinates),
      b=atom.tempFactor,
      occupancy=atom.occupancy)
    assert structure.unit_cell() is not None
    assert structure.space_group_info() is not None
    structure.add_scatterer(scatterer)
  return structure
