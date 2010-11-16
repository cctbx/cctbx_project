"""Parser for Kim Henrick's standard_geometry.cif file:
  http://www.ebi.ac.uk/~henrick/newpdb/standard_geometry.cif
"""

def get_cif_number(type, s):
  s = s.strip()
  if (s == "."): return None
  assert len(s) != 0
  assert s.find(" ") < 0
  return type(s)

def get_cif_float(s): return get_cif_number(float, s)
def get_cif_int(s): return get_cif_number(int, s)

class chem_comp_atom(object):

  __slots__ = [
    "atom_id",
    "alt_atom_id",
    "type_symbol",
    "charge",
    "pdbx_align",
    "pdbx_aromatic_flag",
    "pdbx_leaving_atom_flag",
    "pdbx_stereo_config",
    "ccp4_type_energy",
    "ccp4_partial_charge",
    "pdbx_model_cartn_x_ideal",
    "pdbx_model_cartn_y_ideal",
    "pdbx_model_cartn_z_ideal",
    "pdbx_ordinal",
    "ref_radii_asa_probe",
    "descriptor"]

  def site_ideal(self):
    return (self.pdbx_model_cartn_x_ideal,
            self.pdbx_model_cartn_y_ideal,
            self.pdbx_model_cartn_z_ideal)

class process_chem_comp_atom_buffer(object):

  __slots__ = ["comp_id", "atoms"]

  def __init__(self, buffer):
    comd_id = None
    self.atoms = []
    for line in buffer:
      if (comd_id is None):
        comp_id = self.comp_id = line[:3]
      else:
        assert line[:3] == comd_id
      assert line[3:5] == "  "
      a = chem_comp_atom()
      a.atom_id = line[5:9]
      assert line[9:12] == "   "
      a.alt_atom_id = line[12:16]
      assert line[16] == " "
      a.type_symbol = line[17:19]
      assert line[19] == " "
      a.charge = line[20]
      assert line[21] == " "
      a.pdbx_align = line[22]
      assert line[23] == " "
      a.pdbx_aromatic_flag = line[24]
      assert line[25] == " "
      a.pdbx_leaving_atom_flag = line[26]
      assert line[27] == " "
      a.pdbx_stereo_config = line[28]
      assert line[29:33] == "    "
      a.ccp4_type_energy = line[33:37]
      assert line[37:41] == "    "
      a.ccp4_partial_charge = get_cif_float(line[41:48])
      assert line[48:52] == "    "
      a.pdbx_model_cartn_x_ideal = get_cif_float(line[52:59])
      assert line[59] == " "
      a.pdbx_model_cartn_y_ideal = get_cif_float(line[60:67])
      assert line[67] == " "
      a.pdbx_model_cartn_z_ideal = get_cif_float(line[68:75])
      assert line[75] == " "
      a.pdbx_ordinal = get_cif_int(line[76:79])
      assert line[79] == " "
      a.ref_radii_asa_probe = get_cif_float(line[80:85])
      if (a.ref_radii_asa_probe is not None):
        assert line[85] == " "
        a.descriptor = line[86:].strip()
      else:
        assert line[80:85] == "  .  "
        a.descriptor = line[85:].strip()
      self.atoms.append(a)

  def count_non_hydrogen_atoms(self):
    result = 0
    for atom in self.atoms:
      if (atom.type_symbol not in [" H", "H "]):
        result += 1
    return result

def process_chem_comps(file_name):
  lines = open(file_name).read().splitlines()
  line_iter = iter(lines)
  for line in line_iter:
    if (line.rstrip() == "_chem_comp_atom.descriptor"):
      break
  else:
    raise RuntimeError("Unexpected end of file.")
  chem_comps = {}
  buffer = []
  for line in line_iter:
    if (line.rstrip() == "loop_"):
      break
    if (line.strip() == "#"):
      if (len(buffer) != 0):
        chem_comp = process_chem_comp_atom_buffer(buffer=buffer)
        chem_comps[chem_comp.comp_id] = chem_comp
        buffer = []
    else:
      buffer.append(line)
  else:
    raise RuntimeError("Unexpected end of file.")
  assert len(buffer) == 0
  return chem_comps


if __name__=="__main__":
  import sys
  process_chem_comps(sys.argv[1])
