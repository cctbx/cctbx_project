from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import sys
from iotbx.pdb.tst_secondary_structure import pdb_1ywf_sample_strings, \
    get_annotation
import iotbx.pdb
from iotbx.pdb.utils import all_chain_ids
from iotbx.pdb.secondary_structure import annotation
from six.moves import cStringIO as StringIO


def tst_parsing_phil_single_helix():
  phil_hbond_portion = """\
  hbond {
    donor = chain A and resid 37 and name O
    acceptor = chain A and resid 41 and name N
    }
  hbond {
    donor = chain A and resid 38 and name O
    acceptor = chain A and resid 42 and name N
    }
  hbond {
    donor = chain A and resid 39 and name O
    acceptor = chain A and resid 43 and name N
    }
  hbond {
    donor = chain A and resid 40 and name O
    acceptor = chain A and resid 44 and name N
    }
  hbond {
    donor = chain A and resid 41 and name O
    acceptor = chain A and resid 45 and name N
    }
  hbond {
    donor = chain A and resid 42 and name O
    acceptor = chain A and resid 46 and name N
    }
  hbond {
    donor = chain A and resid 43 and name O
    acceptor = chain A and resid 47 and name N
    }
  hbond {
    donor = chain A and resid 44 and name O
    acceptor = chain A and resid 48 and name N
    }
  """
  phil_str1 = """\
secondary_structure.protein.helix {
  selection = chain A and resseq 37:48
  helix_type = *alpha pi 3_10 unknown
  %s
  }
  """ % phil_hbond_portion

  phil_str2 = """\
secondary_structure.protein.helix {
  selection = chain A and resseq 37:48
  helix_type = alpha *pi 3_10 unknown
  }
  """

  phil_str3 = """\
secondary_structure.protein.helix {
  selection = chain A and resseq 37:48
  helix_type = alpha pi *3_10 unknown
  }
  """

  phil_str4 = """\
secondary_structure.protein.helix {
  selection = chain A and resseq 37:48
  helix_type = alpha pi 3_10 *unknown
  }
  """

  phil_str5 = """\
secondary_structure.protein.helix {
  serial_number = 1
  selection = chain A and resseq 37:48
  helix_type = *alpha pi 3_10 unknown
  }
  """

  phil_str6 = """\
secondary_structure.protein.helix {
  serial_number = 2
  helix_identifier = A
  selection = chain A and resseq 37:48
  helix_type = *alpha pi 3_10 unknown
  }
  """

  phil_str7 = """\
secondary_structure.protein.helix {
  serial_number = 3
  helix_identifier = BB
  selection = chain A and resseq 37:48
  helix_type = *alpha pi 3_10 unknown
  }
  """

  phil_str8 = """\
secondary_structure.protein.helix {
  serial_number = 1
  helix_identifier = C12
  selection = chain A and resseq 37:48
  helix_type = *alpha pi 3_10 unknown
  }
  """

  phil_str9 = """\
secondary_structure.protein.helix {
  serial_number = 1
  helix_identifier = D123
  selection = chain A and resseq 37:48
  helix_type = *alpha pi 3_10 unknown
  }
  """


  result1 = """\
protein.helix {
  helix_identifier = 0
  selection = chain 'A' and resid   37  through   48
  helix_type = alpha
  hbond {
    donor = chain A and resid 37 and name O
    acceptor = chain A and resid 41 and name N
  }
  hbond {
    donor = chain A and resid 38 and name O
    acceptor = chain A and resid 42 and name N
  }
  hbond {
    donor = chain A and resid 39 and name O
    acceptor = chain A and resid 43 and name N
  }
  hbond {
    donor = chain A and resid 40 and name O
    acceptor = chain A and resid 44 and name N
  }
  hbond {
    donor = chain A and resid 41 and name O
    acceptor = chain A and resid 45 and name N
  }
  hbond {
    donor = chain A and resid 42 and name O
    acceptor = chain A and resid 46 and name N
  }
  hbond {
    donor = chain A and resid 43 and name O
    acceptor = chain A and resid 47 and name N
  }
  hbond {
    donor = chain A and resid 44 and name O
    acceptor = chain A and resid 48 and name N
  }
}"""
  result1_1 = """\
protein.helix {
  helix_identifier = 0
  selection = chain 'A' and resid   37  through   48
  helix_type = alpha
}"""
  result2_9 = """\
==================================================
protein.helix {
  helix_identifier = 0
  selection = chain 'A' and resid   37  through   48
  helix_type = pi
}
HELIX    0   0 ASP A   37  GLY A   48  3                                  12
==================================================
protein.helix {
  helix_identifier = 0
  selection = chain 'A' and resid   37  through   48
  helix_type = 3_10
}
HELIX    0   0 ASP A   37  GLY A   48  5                                  12
==================================================
protein.helix {
  helix_identifier = 0
  selection = chain 'A' and resid   37  through   48
  helix_type = unknown
}
HELIX    0   0 ASP A   37  GLY A   48  1                                  12
==================================================
protein.helix {
  serial_number = 1
  helix_identifier = 1
  selection = chain 'A' and resid   37  through   48
  helix_type = alpha
}
HELIX    1   1 ASP A   37  GLY A   48  1                                  12
==================================================
protein.helix {
  serial_number = 2
  helix_identifier = A
  selection = chain 'A' and resid   37  through   48
  helix_type = alpha
}
HELIX    2   A ASP A   37  GLY A   48  1                                  12
==================================================
protein.helix {
  serial_number = 3
  helix_identifier = BB
  selection = chain 'A' and resid   37  through   48
  helix_type = alpha
}
HELIX    3  BB ASP A   37  GLY A   48  1                                  12
==================================================
protein.helix {
  serial_number = 1
  helix_identifier = C12
  selection = chain 'A' and resid   37  through   48
  helix_type = alpha
}
HELIX    1 C12 ASP A   37  GLY A   48  1                                  12
==================================================
protein.helix {
  serial_number = 1
  helix_identifier = D123
  selection = chain 'A' and resid   37  through   48
  helix_type = alpha
}
HELIX    1 D12 ASP A   37  GLY A   48  1                                  12
"""

  annot, ss_from_file = get_annotation(
      phil_lines=phil_str1,
      pdb_lines=pdb_1ywf_sample_strings)
  assert annot.get_n_helices() == 1
  assert annot.get_n_sheets() == 0
  h = annot.helices[0]
  assert h.get_n_defined_hbonds() == 8
  assert annot.get_n_defined_hbonds() == 8
  res = h.as_restraint_group(show_hbonds=True)
  assert not test_utils.show_diff(res, result1,
      strip_trailing_whitespace=True)
  res = h.as_restraint_group(show_hbonds=False)
  assert not test_utils.show_diff(res, result1_1,
      strip_trailing_whitespace=True)
  out = StringIO()
  for ph_str in [phil_str2, phil_str3, phil_str4, phil_str5, phil_str6,
      phil_str7, phil_str8, phil_str9]:
    print("="*50, file=out)
    annot, ss_from_file = get_annotation(
        phil_lines=ph_str,
        pdb_lines=pdb_1ywf_sample_strings)
    assert annot.get_n_defined_hbonds() == 0
    assert annot.get_n_helices() == 1
    assert annot.get_n_sheets() == 0
    h = annot.helices[0]
    print(h.as_restraint_group(), file=out)
    print(h.as_pdb_str(), file=out)
  assert not test_utils.show_diff(out.getvalue(), result2_9,
      strip_trailing_whitespace=True)
  print("OK")

def tst_parsing_phil_single_sheet():
  phil_hbond_portion = """\
  hbond {
    donor = chain A and resid 29 and name O
    acceptor = chain A and resid 13 and name N
  }
  hbond {
    donor = chain A and resid 13 and name O
    acceptor = chain A and resid 29 and name N
  }
  hbond {
    donor = chain A and resid 53 and name O
    acceptor = chain A and resid 159 and name N
  }
  hbond {
    donor = chain A and resid 157 and name O
    acceptor = chain A and resid 53 and name N
  }
  hbond {
    donor = chain A and resid 51 and name O
    acceptor = chain A and resid 157 and name N
  }
  hbond {
    donor = chain A and resid 76 and name O
    acceptor = chain A and resid 54 and name N
  }
  hbond {
    donor = chain A and resid 52 and name O
    acceptor = chain A and resid 76 and name N
  }
  hbond {
    donor = chain A and resid 74 and name O
    acceptor = chain A and resid 52 and name N
  }
  """
  phil_str_1 = """\
  secondary_structure.protein.sheet {
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
    %s
  }
  """ % phil_hbond_portion
  phil_str_1_1 = """\
  secondary_structure.protein.sheet {
    first_strand = chain A and resseq 13:14
    sheet_id = AA
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
  """
  result1 = """
protein.sheet {
  sheet_id = "  A"
  first_strand = chain 'A' and resid   13  through   14
  strand {
    selection = chain 'A' and resid   27  through   30
    sense = antiparallel
    bond_start_current = chain 'A' and resid   29  and name O
    bond_start_previous = chain 'A' and resid   13  and name N
  }
  strand {
    selection = chain 'A' and resid  156  through  159
    sense = parallel
    bond_start_current = chain 'A' and resid  156  and name O
    bond_start_previous = chain 'A' and resid   28  and name N
  }
  strand {
    selection = chain 'A' and resid   51  through   54
    sense = parallel
    bond_start_current = chain 'A' and resid   51  and name O
    bond_start_previous = chain 'A' and resid  157  and name N
  }
  strand {
    selection = chain 'A' and resid   74  through   77
    sense = parallel
    bond_start_current = chain 'A' and resid   74  and name O
    bond_start_previous = chain 'A' and resid   52  and name N
  }
  hbond {
    donor = chain A and resid 29 and name O
    acceptor = chain A and resid 13 and name N
  }
  hbond {
    donor = chain A and resid 13 and name O
    acceptor = chain A and resid 29 and name N
  }
  hbond {
    donor = chain A and resid 53 and name O
    acceptor = chain A and resid 159 and name N
  }
  hbond {
    donor = chain A and resid 157 and name O
    acceptor = chain A and resid 53 and name N
  }
  hbond {
    donor = chain A and resid 51 and name O
    acceptor = chain A and resid 157 and name N
  }
  hbond {
    donor = chain A and resid 76 and name O
    acceptor = chain A and resid 54 and name N
  }
  hbond {
    donor = chain A and resid 52 and name O
    acceptor = chain A and resid 76 and name N
  }
  hbond {
    donor = chain A and resid 74 and name O
    acceptor = chain A and resid 52 and name N
  }
}"""
  result1_1 = """
protein.sheet {
  sheet_id = "  A"
  first_strand = chain 'A' and resid   13  through   14
  strand {
    selection = chain 'A' and resid   27  through   30
    sense = antiparallel
    bond_start_current = chain 'A' and resid   29  and name O
    bond_start_previous = chain 'A' and resid   13  and name N
  }
  strand {
    selection = chain 'A' and resid  156  through  159
    sense = parallel
    bond_start_current = chain 'A' and resid  156  and name O
    bond_start_previous = chain 'A' and resid   28  and name N
  }
  strand {
    selection = chain 'A' and resid   51  through   54
    sense = parallel
    bond_start_current = chain 'A' and resid   51  and name O
    bond_start_previous = chain 'A' and resid  157  and name N
  }
  strand {
    selection = chain 'A' and resid   74  through   77
    sense = parallel
    bond_start_current = chain 'A' and resid   74  and name O
    bond_start_previous = chain 'A' and resid   52  and name N
  }
}"""

  annot, ss_from_file = get_annotation(
      phil_lines=phil_str_1,
      pdb_lines=pdb_1ywf_sample_strings)
  assert annot.get_n_helices() == 0
  assert annot.get_n_sheets() == 1
  sh = annot.sheets[0]
  assert sh.get_n_defined_hbonds() == 8
  assert annot.get_n_defined_hbonds() == 8
  res = sh.as_restraint_group(show_hbonds=True)
  assert not test_utils.show_diff(res, result1,
      strip_trailing_whitespace=True)
  res = sh.as_restraint_group(show_hbonds=False)
  assert not test_utils.show_diff(res, result1_1,
      strip_trailing_whitespace=True)
  print("OK")


def tst_hybrid_ids():
  inp = iotbx.pdb.input(source_info=None, lines=pdb_1ywf_sample_strings)
  ss_ann = inp.extract_secondary_structure()
  ch_ids = all_chain_ids()[1:1000]
  chain_ids_dict = {'A': ch_ids}
  ss_ann.multiply_to_asu_2(chain_ids_dict)
  # print (ss_ann.as_pdb_str())
  ss_big = annotation.from_records(ss_ann.as_pdb_str().split('\n'))
  # print (ss_big.as_restraint_groups())
  ss_big.as_restraint_groups()

def exercise(args):
  tst_parsing_phil_single_helix()
  tst_parsing_phil_single_sheet()
  tst_hybrid_ids()


if (__name__ == "__main__"):
  exercise(sys.argv[1:])
