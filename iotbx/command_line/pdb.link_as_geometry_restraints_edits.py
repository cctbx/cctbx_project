from iotbx import pdb
from libtbx.str_utils import show_string
from libtbx import adopt_init_args
import sys

def field_as_number(
      line,
      icols,
      target_type,
      strict,
      error_message,
      substitute_value=0):
  try: return target_type(line[icols[0]:icols[1]])
  except ValueError:
    if (strict): raise RuntimeError(
      error_message+"\n"
      + "  PDB input line: " + line.rstrip())
  return substitute_value

class atom_labels(object):

  def __init__(self,
        chain=None,
        resname=None,
        resseq=None,
        icode=None,
        name=None,
        altloc=None,
        segid=None):
    adopt_init_args(self, locals())

  def selection_string(self):
    result = []
    a = result.append
    a('chain %s' % show_string(self.chain))
    a('resname %s' % show_string(self.resname))
    a('resseq %s' % self.resseq)
    a('icode %s' % show_string(self.icode))
    a('name %s' % show_string(self.name))
    a('altloc %s' % show_string(self.altloc))
    return show_string(" and ".join(result))

class link_record(object):

  def __init__(self, line):
    # 13 - 16      Atom            name1       Atom name.
    # 17           Character       altLoc1     Alternate location indicator.
    # 18 - 20      Residue name    resName1    Residue name.
    # 22           Character       chainID1    Chain identifier.
    # 23 - 26      Integer         resSeq1     Residue sequence number.
    # 27           AChar           iCode1      Insertion code.
    # 31 - 40      distance (REFMAC extension: F10.5)
    # 43 - 46      Atom            name2       Atom name.
    # 47           Character       altLoc2     Alternate location indicator.
    # 48 - 50      Residue name    resName2    Residue name.
    # 52           Character       chainID2    Chain identifier.
    # 53 - 56      Integer         resSeq2     Residue sequence number.
    # 57           AChar           iCode2      Insertion code.
    # 60 - 65      SymOP           sym1        Symmetry operator for 1st atom.
    # 67 - 72      SymOP           sym2        Symmetry operator for 2nd atom.
    # 73 - 80      margin (REFMAC extension: _chem_link.id)
    line = line + " "*max(0,80-len(line))
    assert line[:6] == "LINK  "
    self.labels_pair = [atom_labels(
      name=line[12:16],
      altloc=line[16],
      resname=line[17:20],
      chain=line[21],
      resseq=field_as_number(
        line=line, icols=(22,26), target_type=int, strict=True,
        error_message="Serial number must be an integer:"),
      icode=line[26])]
    try: self.distance = float(line[30:40])
    except ValueError: self.distance = None
    self.labels_pair.append(atom_labels(
      name=line[42:46],
      altloc=line[46],
      resname=line[47:50],
      chain=line[51],
      resseq=field_as_number(
        line=line, icols=(52,56), target_type=int, strict=True,
        error_message="Serial number must be an integer:"),
      icode=line[56]))
    self.symops = [line[59:65], line[66:72]]
    self.margin = line[72:80]

  def as_geometry_restraints_edits(self):
    if (self.distance is None):
      distance_ideal = "None"
    else:
      distance_ideal = "%.6g" % self.distance
    return """\
  bond {
    action = *add delete change
    atom_selection_1 = %s
    atom_selection_2 = %s
    symmetry_operation = None
    distance_ideal = %s
    sigma = None
  }
""" % tuple([labels.selection_string() for labels in self.labels_pair]
            + [distance_ideal])

def run(args, command_name="iotbx.pdb.link_as_geometry_restraints_edits"):
  for file_name in args:
    section = pdb.input(file_name=file_name).connectivity_annotation_section()
    print "refinement.geometry_restraints.edits {"
    for line in section:
      if (not line.startswith("LINK  ")): continue
      link = link_record(line=line)
      assert link.symops[0].strip() == "" # not implemented
      assert link.symops[1].strip() == "" # not implemented
      sys.stdout.write(link.as_geometry_restraints_edits())
    print "}"

if (__name__ == "__main__"):
  run(sys.argv[1:])
