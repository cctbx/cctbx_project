from cctbx import xray
from cctbx import adptbx
from cStringIO import StringIO

def xray_structure_as_pdb_file(self, remark=None, remarks=[],
                                     fractional_coordinates=00000):
  if (remark != None):
    remarks.insert(0, remark)
  i = 0
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
    self.unit_cell().parameters() + (str(self.space_group_info()),))
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
  for scatterer in self.scatterers():
    assert not scatterer.anisotropic_flag, "Not implemented." # XXX
    xyz = scatterer.site
    if (not fractional_coordinates):
      xyz = self.unit_cell().orthogonalize(xyz)
    i += 1
    label = scatterer.label.upper()
    print >> s, "ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f" % ((
      i,
      label,
      label,
      i,
      ) + xyz + (
      scatterer.occupancy,
      adptbx.u_as_b(scatterer.u_iso)))
  print >> s, "END"
  return s.getvalue()

xray.structure.as_pdb_file = xray_structure_as_pdb_file
