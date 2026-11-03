from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from mmtbx.secondary_structure import ss_validation
from libtbx.utils import null_out
from six.moves import cStringIO as StringIO

# An explicit mmCIF model: an 8-residue ideal poly-ALA alpha helix on chain
# "AAA" - a chain id longer than 2 characters, which only mmCIF can hold - with
# the helix annotation supplied directly as a _struct_conf record. Reading this
# reproduces exactly "a model supplied in mmcif format with chain ids longer
# than 2 characters"; no PDB->mmCIF reformatting is involved.
helix_in_cif = """\
data_default
_cell.length_a                    14.213
_cell.length_b                    13.872
_cell.length_c                    18.795
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.volume                      3705.674
_space_group.crystal_system       triclinic
_space_group.IT_number            1
_space_group.name_H-M_alt         'P 1'
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.space_group_name_Hall   ' P 1'
_symmetry.Int_Tables_number       1
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _struct_asym.id
   A

loop_
  _chem_comp.id
   ALA

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.auth_atom_id
  _atom_site.pdbx_PDB_model_num
   ATOM 1 N . ALA AAA 1 ? -1.20400 -0.51400 0.64300 1.000 0.00000 N ? A ? 1 N 1
   ATOM 2 CA . ALA AAA 1 ? 0.00000 0.00000 0.00000 1.000 0.00000 C ? A ? 1 CA 1
   ATOM 3 C . ALA AAA 1 ? 0.86600 -1.13400 -0.53700 1.000 0.00000 C ? A ? 1 C 1
   ATOM 4 O . ALA AAA 1 ? 1.43200 -1.03500 -1.62500 1.000 0.00000 O ? A ? 1 O 1
   ATOM 5 CB . ALA AAA 1 ? 0.80800 0.86000 0.97400 1.000 0.00000 C ? A ? 1 CB 1
   ATOM 6 N . ALA AAA 2 ? 0.96500 -2.21300 0.23400 1.000 0.00000 N ? A ? 2 N 1
   ATOM 7 CA . ALA AAA 2 ? 1.76100 -3.36800 -0.16200 1.000 0.00000 C ? A ? 2 CA 1
   ATOM 8 C . ALA AAA 2 ? 1.25800 -3.96300 -1.47200 1.000 0.00000 C ? A ? 2 C 1
   ATOM 9 O . ALA AAA 2 ? 2.04700 -4.36600 -2.32700 1.000 0.00000 O ? A ? 2 O 1
   ATOM 10 CB . ALA AAA 2 ? 1.75100 -4.43400 0.93600 1.000 0.00000 C ? A ? 2 CB 1
   ATOM 11 N . ALA AAA 3 ? -0.06200 -4.01600 -1.62500 1.000 0.00000 N ? A ? 3 N 1
   ATOM 12 CA . ALA AAA 3 ? -0.67300 -4.56200 -2.83000 1.000 0.00000 C ? A ? 3 CA 1
   ATOM 13 C . ALA AAA 3 ? -0.24500 -3.78200 -4.06800 1.000 0.00000 C ? A ? 3 C 1
   ATOM 14 O . ALA AAA 3 ? 0.01100 -4.36400 -5.12300 1.000 0.00000 O ? A ? 3 O 1
   ATOM 15 CB . ALA AAA 3 ? -2.19900 -4.56000 -2.71100 1.000 0.00000 C ? A ? 3 CB 1
   ATOM 16 N . ALA AAA 4 ? -0.16800 -2.46200 -3.93400 1.000 0.00000 N ? A ? 4 N 1
   ATOM 17 CA . ALA AAA 4 ? 0.22900 -1.60000 -5.04000 1.000 0.00000 C ? A ? 4 CA 1
   ATOM 18 C . ALA AAA 4 ? 1.62800 -1.94600 -5.53700 1.000 0.00000 C ? A ? 4 C 1
   ATOM 19 O . ALA AAA 4 ? 1.88800 -1.94900 -6.74000 1.000 0.00000 O ? A ? 4 O 1
   ATOM 20 CB . ALA AAA 4 ? 0.16900 -0.12800 -4.62700 1.000 0.00000 C ? A ? 4 CB 1
   ATOM 21 N . ALA AAA 5 ? 2.52900 -2.23900 -4.60200 1.000 0.00000 N ? A ? 5 N 1
   ATOM 22 CA . ALA AAA 5 ? 3.90200 -2.58700 -4.94300 1.000 0.00000 C ? A ? 5 CA 1
   ATOM 23 C . ALA AAA 5 ? 3.95400 -3.82700 -5.82700 1.000 0.00000 C ? A ? 5 C 1
   ATOM 24 O . ALA AAA 5 ? 4.74400 -3.89900 -6.76700 1.000 0.00000 O ? A ? 5 O 1
   ATOM 25 CB . ALA AAA 5 ? 4.73200 -2.81100 -3.67700 1.000 0.00000 C ? A ? 5 CB 1
   ATOM 26 N . ALA AAA 6 ? 3.10400 -4.80400 -5.52000 1.000 0.00000 N ? A ? 6 N 1
   ATOM 27 CA . ALA AAA 6 ? 3.05200 -6.04300 -6.28600 1.000 0.00000 C ? A ? 6 CA 1
   ATOM 28 C . ALA AAA 6 ? 2.70700 -5.77400 -7.74700 1.000 0.00000 C ? A ? 6 C 1
   ATOM 29 O . ALA AAA 6 ? 3.26600 -6.39400 -8.65300 1.000 0.00000 O ? A ? 6 O 1
   ATOM 30 CB . ALA AAA 6 ? 2.03600 -7.01200 -5.67800 1.000 0.00000 C ? A ? 6 CB 1
   ATOM 31 N . ALA AAA 7 ? 1.78300 -4.84500 -7.97200 1.000 0.00000 N ? A ? 7 N 1
   ATOM 32 CA . ALA AAA 7 ? 1.36100 -4.49300 -9.32200 1.000 0.00000 C ? A ? 7 CA 1
   ATOM 33 C . ALA AAA 7 ? 2.53300 -3.97600 -10.14900 1.000 0.00000 C ? A ? 7 C 1
   ATOM 34 O . ALA AAA 7 ? 2.65800 -4.29400 -11.33200 1.000 0.00000 O ? A ? 7 O 1
   ATOM 35 CB . ALA AAA 7 ? 0.24500 -3.44700 -9.28400 1.000 0.00000 C ? A ? 7 CB 1
   ATOM 36 N . ALA AAA 8 ? 3.39000 -3.17800 -9.52100 1.000 0.00000 N ? A ? 8 N 1
   ATOM 37 CA . ALA AAA 8 ? 4.55200 -2.61500 -10.19700 1.000 0.00000 C ? A ? 8 CA 1
   ATOM 38 C . ALA AAA 8 ? 5.47200 -3.71400 -10.72000 1.000 0.00000 C ? A ? 8 C 1
   ATOM 39 O . ALA AAA 8 ? 6.01400 -3.61300 -11.82100 1.000 0.00000 O ? A ? 8 O 1
   ATOM 40 CB . ALA AAA 8 ? 5.32600 -1.68600 -9.26000 1.000 0.00000 C ? A ? 8 CB 1
loop_
  _struct_conf.conf_type_id
  _struct_conf.id
  _struct_conf.pdbx_PDB_helix_id
  _struct_conf.beg_label_comp_id
  _struct_conf.beg_label_asym_id
  _struct_conf.beg_label_seq_id
  _struct_conf.pdbx_beg_PDB_ins_code
  _struct_conf.end_label_comp_id
  _struct_conf.end_label_asym_id
  _struct_conf.end_label_seq_id
  _struct_conf.pdbx_end_PDB_ins_code
  _struct_conf.pdbx_PDB_helix_class
  _struct_conf.details
  _struct_conf.pdbx_PDB_helix_length
  HELX_P  1  1  ALA  AAA  1  ?  ALA  AAA  8  ?  1  ?  8

loop_
  _struct_conf_type.id
  _struct_conf_type.criteria
  _struct_conf_type.reference
  HELX_P  ?  ?
"""

# Same idea, but the _struct_conf helix and the 2-strand _struct_sheet point at
# residues that are absent from the (3-residue) model, so both are reported as
# records "without corresponding atoms" - the other two as_pdb_str() call sites.
empty_annotations_in_cif = """\
data_default
_cell.length_a                    10.246
_cell.length_b                    11.422
_cell.length_c                    12.097
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.volume                      1415.710
_space_group.crystal_system       triclinic
_space_group.IT_number            1
_space_group.name_H-M_alt         'P 1'
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.space_group_name_Hall   ' P 1'
_symmetry.Int_Tables_number       1
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z

loop_
  _struct_asym.id
   A

loop_
  _chem_comp.id
   ALA

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.auth_atom_id
  _atom_site.pdbx_PDB_model_num
   ATOM 1 N . ALA AAA 1 ? -1.20400 -0.51400 0.64300 1.000 0.00000 N ? A ? 1 N 1
   ATOM 2 CA . ALA AAA 1 ? 0.00000 0.00000 0.00000 1.000 0.00000 C ? A ? 1 CA 1
   ATOM 3 C . ALA AAA 1 ? 0.86600 -1.13400 -0.53700 1.000 0.00000 C ? A ? 1 C 1
   ATOM 4 O . ALA AAA 1 ? 1.43200 -1.03500 -1.62500 1.000 0.00000 O ? A ? 1 O 1
   ATOM 5 CB . ALA AAA 1 ? 0.80800 0.86000 0.97400 1.000 0.00000 C ? A ? 1 CB 1
   ATOM 6 N . ALA AAA 2 ? 0.96500 -2.21300 0.23400 1.000 0.00000 N ? A ? 2 N 1
   ATOM 7 CA . ALA AAA 2 ? 1.76100 -3.36800 -0.16200 1.000 0.00000 C ? A ? 2 CA 1
   ATOM 8 C . ALA AAA 2 ? 1.25800 -3.96300 -1.47200 1.000 0.00000 C ? A ? 2 C 1
   ATOM 9 O . ALA AAA 2 ? 2.04700 -4.36600 -2.32700 1.000 0.00000 O ? A ? 2 O 1
   ATOM 10 CB . ALA AAA 2 ? 1.75100 -4.43400 0.93600 1.000 0.00000 C ? A ? 2 CB 1
   ATOM 11 N . ALA AAA 3 ? -0.06200 -4.01600 -1.62500 1.000 0.00000 N ? A ? 3 N 1
   ATOM 12 CA . ALA AAA 3 ? -0.67300 -4.56200 -2.83000 1.000 0.00000 C ? A ? 3 CA 1
   ATOM 13 C . ALA AAA 3 ? -0.24500 -3.78200 -4.06800 1.000 0.00000 C ? A ? 3 C 1
   ATOM 14 O . ALA AAA 3 ? 0.01100 -4.36400 -5.12300 1.000 0.00000 O ? A ? 3 O 1
   ATOM 15 CB . ALA AAA 3 ? -2.19900 -4.56000 -2.71100 1.000 0.00000 C ? A ? 3 CB 1
loop_
  _struct_conf.conf_type_id
  _struct_conf.id
  _struct_conf.pdbx_PDB_helix_id
  _struct_conf.beg_label_comp_id
  _struct_conf.beg_label_asym_id
  _struct_conf.beg_label_seq_id
  _struct_conf.pdbx_beg_PDB_ins_code
  _struct_conf.end_label_comp_id
  _struct_conf.end_label_asym_id
  _struct_conf.end_label_seq_id
  _struct_conf.pdbx_end_PDB_ins_code
  _struct_conf.pdbx_PDB_helix_class
  _struct_conf.details
  _struct_conf.pdbx_PDB_helix_length
  HELX_P  1  1  ALA  AAA  100  ?  ALA  AAA  110  ?  1  ?  11

loop_
  _struct_conf_type.id
  _struct_conf_type.criteria
  _struct_conf_type.reference
  HELX_P  ?  ?

loop_
  _struct_sheet.id
  _struct_sheet.type
  _struct_sheet.number_strands
  _struct_sheet.details
  A  ?  2  ?

loop_
  _struct_sheet_order.sheet_id
  _struct_sheet_order.range_id_1
  _struct_sheet_order.range_id_2
  _struct_sheet_order.offset
  _struct_sheet_order.sense
  A  1  2  ?  anti-parallel

loop_
  _struct_sheet_range.sheet_id
  _struct_sheet_range.id
  _struct_sheet_range.beg_label_comp_id
  _struct_sheet_range.beg_label_asym_id
  _struct_sheet_range.beg_label_seq_id
  _struct_sheet_range.pdbx_beg_PDB_ins_code
  _struct_sheet_range.end_label_comp_id
  _struct_sheet_range.end_label_asym_id
  _struct_sheet_range.end_label_seq_id
  _struct_sheet_range.pdbx_end_PDB_ins_code
  A  1  ALA  AAA  200  ?  ALA  AAA  205  ?
  A  2  ALA  AAA  210  ?  ALA  AAA  215  ?

loop_
  _pdbx_struct_sheet_hbond.sheet_id
  _pdbx_struct_sheet_hbond.range_id_1
  _pdbx_struct_sheet_hbond.range_id_2
  _pdbx_struct_sheet_hbond.range_1_label_atom_id
  _pdbx_struct_sheet_hbond.range_1_label_comp_id
  _pdbx_struct_sheet_hbond.range_1_label_asym_id
  _pdbx_struct_sheet_hbond.range_1_label_seq_id
  _pdbx_struct_sheet_hbond.range_1_PDB_ins_code
  _pdbx_struct_sheet_hbond.range_2_label_atom_id
  _pdbx_struct_sheet_hbond.range_2_label_comp_id
  _pdbx_struct_sheet_hbond.range_2_label_asym_id
  _pdbx_struct_sheet_hbond.range_2_label_seq_id
  _pdbx_struct_sheet_hbond.range_2_PDB_ins_code
  A  1  2  O  ALA  AAA  205  ?  N  ALA  AAA  210  ?
"""

def model_from_cif(cif_str):
  return mmtbx.model.manager(
      model_input=iotbx.pdb.input(source_info=None, lines=cif_str),
      log=null_out())

def get_validation_params(bad_hbond_cutoff):
  params = ss_validation.validate.get_default_params().ss_validation
  params.nproc = 1
  # Drive every detected H-bond into the "bad" bucket so that the helix is
  # flagged as a bad annotation and its record gets printed (the failing path).
  params.bad_hbond_cutoff = bad_hbond_cutoff
  return params

def tst_helix_with_long_chain_id():
  """A helix on a >2 char chain id must not crash ss_validation when its
  record is printed; it has to fall back to mmCIF output."""
  model = model_from_cif(helix_in_cif)
  log = StringIO()
  ss_validation.validate(
      model=model,
      params=get_validation_params(bad_hbond_cutoff=0.5),
      log=log)
  out = log.getvalue()
  print(out)
  print('='*80)
  assert "Bad annotation found:" in out, out
  # The helix does not fit in PDB format, so the printed record must be mmCIF.
  assert "_struct_conf" in out, out

def tst_empty_annotations_with_long_chain_id():
  """Empty helix/sheet records (no matching atoms) on a >2 char chain id must
  not crash ss_validation when they are reported; they have to fall back to
  mmCIF output. Exercises the other two as_pdb_str() call sites."""
  model = model_from_cif(empty_annotations_in_cif)
  log = StringIO()
  ss_validation.validate(
      model=model,
      params=get_validation_params(bad_hbond_cutoff=3.5),
      log=log)
  out = log.getvalue()
  print(out)
  assert "Helices without corresponding atoms" in out, out
  assert "Sheets without corresponding atoms" in out, out
  # Neither record fits PDB format, so both must be rendered as mmCIF.
  assert "_struct_conf" in out, out
  assert "_struct_sheet" in out, out

if (__name__ == "__main__"):
  tst_helix_with_long_chain_id()
  tst_empty_annotations_with_long_chain_id()
  print("OK")
