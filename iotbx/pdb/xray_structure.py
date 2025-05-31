"""Convert an xray_structure to PDB file format"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb.hierarchy
from cctbx import adptbx
from six.moves import cStringIO as StringIO

def as_pdb_file(self,
      remark,
      remarks,
      fractional_coordinates,
      resname,
      connect):
  if (remark is not None):
    remarks.insert(0, remark)
  s = StringIO()
  for remark in remarks:
    print("REMARK", remark, file=s)
  print("REMARK Number of scatterers:", self.scatterers().size(), file=s)
  print("REMARK At special positions:", \
    self.special_position_indices().size(), file=s)
  if (fractional_coordinates):
    print("REMARK Fractional coordinates", file=s)
  else:
    print("REMARK Cartesian coordinates", file=s)
  print(iotbx.pdb.format_cryst1_record(crystal_symmetry=self), file=s)
  print(iotbx.pdb.format_scale_records(unit_cell=self.unit_cell()), file=s)
  atom = iotbx.pdb.hierarchy.atom_with_labels()
  if (resname is not None):
    atom.resname = resname.upper()
  serial = 0
  for scatterer in self.scatterers():
    serial += 1
    atom.serial = iotbx.pdb.hy36encode(width=5, value=serial)
    if (scatterer.flags.use_u_aniso_only()):
      atom.uij = adptbx.u_star_as_u_cart(self.unit_cell(), scatterer.u_star)
      atom.b = adptbx.u_as_b(adptbx.u_cart_as_u_iso(atom.uij))
    else:
      atom.uij_erase()
      atom.b = adptbx.u_as_b(scatterer.u_iso)
    if (fractional_coordinates):
      atom.xyz = scatterer.site
    else:
      atom.xyz = self.unit_cell().orthogonalize(scatterer.site)
    atom.occ = scatterer.occupancy
    label = scatterer.label.upper()
    atom.name = label[:4]
    if (resname is None):
      atom.resname = label[:3]
    element_symbol = scatterer.element_symbol()
    if (element_symbol is None): element_symbol = "Q"
    assert len(element_symbol) in (1,2)
    atom.element = element_symbol.upper()
    atom.resseq = iotbx.pdb.hy36encode(width=4, value=serial)
    print(atom.format_atom_record_group(), file=s)
  if (connect is not None):
    assert len(connect) == self.scatterers().size()
    i = 0
    for bonds in connect:
      i += 1
      l = "CONNECT%5d" % i
      for bond in bonds:
        l += "%5d" % (bond+1)
      print(l, file=s)
  print("END", file=s)
  return s.getvalue()
