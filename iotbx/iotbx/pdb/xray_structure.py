import iotbx.pdb
from cctbx import xray
from cctbx import adptbx
from cStringIO import StringIO

def as_pdb_file(self,
      remark,
      remarks,
      fractional_coordinates,
      res_name,
      connect):
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
