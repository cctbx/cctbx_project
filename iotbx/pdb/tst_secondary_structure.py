from libtbx.utils import Sorry, Usage
from iotbx import file_reader
from libtbx import test_utils
import libtbx.phil
from libtbx import adopt_init_args, group_args
import string, sys, os
from iotbx.pdb import secondary_structure as ss
import iotbx


def exercise_single():
  from scitbx.array_family import flex
  from libtbx import test_utils
  two_char_chain_records = """\
HELIX    1   1 THRA1   11  ASPA1   39  1                                  29
HELIX    2   2 GLUA1   46  ARGA1   73  1                                  28
HELIX    3   3 THRA1   93  ALAA1  120  1                                  28
HELIX    4   4 PROA1  124  HISA1  133  1                                  10
HELIX    5   5 LEUA1  135  ARGA1  154  1                                  20
HELIX    6   6 THRa1   11  TYRa1   37  1                                  27
HELIX    7   7 GLUa1   46  GLNa1   72  1                                  27
HELIX    8   8 THRa1   93  ALAa1  120  1                                  28
HELIX    9   9 PROa1  124  HISa1  133  1                                  10"""
  lines = flex.std_string()
  lines.extend(flex.split_lines(two_char_chain_records))
  secstr = ss.annotation.from_records(records=lines)
  ss_out = secstr.as_pdb_str()
  assert not test_utils.show_diff(ss_out, two_char_chain_records)
  print "OK"

def tst_pdb_file():
  if not libtbx.env.has_module("phenix_regression"):
    print "Skipping tst_pdb_file(): phenix_regression not available"
    return
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  pdb_in = file_reader.any_file(pdb_file, force_type="pdb")
  old_ss = pdb_in.file_object.input.secondary_structure_section()
  structure = pdb_in.file_object.input.extract_secondary_structure()
  new_ss = structure.as_pdb_str()
  old_ss = "\n".join(old_ss)
  assert not test_utils.show_diff(new_ss, old_ss)
  print "OK"

def tst_parsing_phil():
  phil_str = """\
secondary_structure {
  helix {
    selection = chain A and resseq 1:18
    helix_type = alpha pi *3_10 unknown
    }
  helix {
    selection = chain A and resseq 37:48
    helix_type = *alpha pi 3_10 unknown
    }
  helix {
    selection = chain A and resseq 57:65
    helix_type = *alpha pi 3_10 unknown
    }
  helix {
    selection = chain A and resseq 119:133
    helix_type = *alpha pi 3_10 unknown
    }
  helix {
    selection = chain A and resseq 134:136
    helix_type = alpha pi *3_10 unknown
    }
  helix {
    selection = chain A and resseq 138:152
    helix_type = *alpha pi 3_10 unknown
    }
  helix {
    selection = chain A and resseq 165:178
    helix_type = *alpha pi 3_10 unknown
    }
  helix {
    selection = chain A and resseq 181:191
    helix_type = *alpha pi 3_10 unknown
    }
  helix {
    selection = chain A and resseq 192:194
    helix_type = alpha pi *3_10 unknown
    }
  helix {
    selection = chain A and resseq 195:209
    helix_type = *alpha pi 3_10 unknown
    }
  helix {
    selection = chain A and resseq 216:225
    helix_type = *alpha pi 3_10 unknown
    }
  helix {
    selection = chain A and resseq 228:233
    helix_type = *alpha pi 3_10 unknown
    }
  helix {
    selection = chain A and resseq 235:251
    helix_type = *alpha pi 3_10 unknown
    }
  helix {
    selection = chain A and resseq 252:260
    helix_type = *alpha pi 3_10 unknown
    }
  helix {
    selection = chain A and resseq 263:275
    helix_type = *alpha pi 3_10 unknown
    }

  sheet {
    first_strand = chain A and resseq 13:14
    sheet_id = A
    strand {
      selection = chain A and resseq 27:30
      sense = antiparallel
      bond_start_current = chain A and resseq 29 and name O
      bond_start_previous = chain A and resseq 13 and name N
    }
    strand {
      selection = chain A and resseq 156:159
      sense = parallel
      bond_start_current = chain A and resseq 156 and name O
      bond_start_previous = chain A and resseq 28 and name N
    }
    strand {
      selection = chain A and resseq 51:54
      sense = parallel
      bond_start_current = chain A and resseq 51 and name O
      bond_start_previous = chain A and resseq 157 and name N
    }
    strand {
      selection = chain A and resseq 74:77
      sense = parallel
      bond_start_current = chain A and resseq 74 and name O
      bond_start_previous = chain A and resseq 52 and name N
    }
  }
}"""
  # copy-paste from pdb file
  equiv_pdb_str = """\
HELIX    1   1 ALA A   16  THR A   18  5                                   3
HELIX    2   2 ASP A   37  GLY A   48  1                                  12
HELIX    3   3 SER A   57  GLY A   65  1                                   9
HELIX    4   4 ASN A  119  PHE A  133  1                                  15
HELIX    5   5 PRO A  134  ARG A  136  5                                   3
HELIX    6   6 GLY A  138  ALA A  152  1                                  15
HELIX    7   7 ASP A  165  VAL A  178  1                                  14
HELIX    8   8 ASP A  181  ARG A  191  1                                  11
HELIX    9   9 SER A  192  ASP A  194  5                                   3
HELIX   10  10 SER A  195  GLN A  209  1                                  15
HELIX   11  11 ALA A  216  ALA A  225  1                                  10
HELIX   12  12 SER A  228  GLY A  233  1                                   6
HELIX   13  13 ARG A  235  GLY A  251  1                                  17
HELIX   14  14 SER A  252  ALA A  260  1                                   9
HELIX   15  15 SER A  263  LEU A  275  1                                  13
SHEET    1   A 5 ARG A  13  ASP A  14  0
SHEET    2   A 5 LEU A  27  SER A  30 -1  O  ARG A  29   N  ARG A  13
SHEET    3   A 5 VAL A 156  HIS A 159  1  O  VAL A 156   N  PHE A  28
SHEET    4   A 5 ASP A  51  ASP A  54  1  N  ALA A  51   O  LEU A 157
SHEET    5   A 5 ASP A  74  LEU A  77  1  O  HIS A  74   N  VAL A  52"""
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  pdb_inp = iotbx.pdb.input(source_info=None, file_name=pdb_file)
  pdb_h = pdb_inp.construct_hierarchy()
  old_ss = pdb_inp.secondary_structure_section()
  import mmtbx.command_line.geometry_minimization
  master_phil = mmtbx.command_line.geometry_minimization.master_params()
  helix_phil = iotbx.phil.parse(phil_str)
  working_params = master_phil.fetch(source=helix_phil).extract()
  phil_helices = working_params.secondary_structure.helix
  phil_sheets = working_params.secondary_structure.sheet
  annot = ss.annotation.from_phil(phil_helices, phil_sheets, pdb_h)
  new_ss = annot.as_pdb_str()
  assert old_ss == new_ss
  print "OK"

def exercise(args):
  exercise_single()
  tst_pdb_file()
  tst_parsing_phil()

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
