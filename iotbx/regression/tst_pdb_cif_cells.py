from __future__ import absolute_import, division, print_function
import iotbx.pdb

#
# Example from PDB ID 7c9c
#

t_pdb_str = """\
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000
ATOM      1  P    DT D   1     166.659 153.609 135.182  1.00 61.01           P
ATOM      2  OP1  DT D   1     166.715 152.154 135.422  1.00 61.01           O
ATOM      3  OP2  DT D   1     166.269 154.124 133.851  1.00 61.01           O
ATOM      4  O5'  DT D   1     165.696 154.272 136.272  1.00 61.01           O
ATOM      5  C5'  DT D   1     166.016 154.187 137.657  1.00 61.01           C
ATOM      6  C4'  DT D   1     166.635 155.480 138.154  1.00 61.01           C
ATOM      7  O4'  DT D   1     165.721 156.583 137.917  1.00 61.01           O
ATOM      8  C3'  DT D   1     166.950 155.515 139.650  1.00 61.01           C
ATOM      9  O3'  DT D   1     168.189 156.162 139.861  1.00 61.01           O
ATOM     10  C2'  DT D   1     165.801 156.336 140.211  1.00 61.01           C
ATOM     11  C1'  DT D   1     165.600 157.334 139.096  1.00 61.01           C
ATOM     12  N1   DT D   1     164.274 157.975 139.127  1.00 61.01           N
ATOM     13  C2   DT D   1     164.180 159.327 138.934  1.00 61.01           C
ATOM     14  O2   DT D   1     165.141 160.036 138.720  1.00 61.01           O
ATOM     15  N3   DT D   1     162.912 159.825 138.994  1.00 61.01           N
ATOM     16  C4   DT D   1     161.753 159.119 139.229  1.00 61.01           C
ATOM     17  O4   DT D   1     160.654 159.657 139.264  1.00 61.01           O
ATOM     18  C5   DT D   1     161.924 157.705 139.429  1.00 61.01           C
ATOM     19  C7   DT D   1     160.735 156.834 139.691  1.00 61.01           C
ATOM     20  C6   DT D   1     163.163 157.207 139.374  1.00 61.01           C
ATOM     21  P    DT D   2     168.926 156.085 141.287  1.00 56.21           P
ATOM     22  OP1  DT D   2     170.334 155.745 141.021  1.00 56.21           O
ATOM     23  OP2  DT D   2     168.132 155.242 142.199  1.00 56.21           O
ATOM     24  O5'  DT D   2     168.871 157.586 141.817  1.00 56.21           O
ATOM     25  C5'  DT D   2     169.222 158.639 140.942  1.00 56.21           C
ATOM     26  C4'  DT D   2     168.956 159.998 141.561  1.00 56.21           C
ATOM     27  O4'  DT D   2     167.552 160.326 141.452  1.00 56.21           O
ATOM     28  C3'  DT D   2     169.313 160.123 143.038  1.00 56.21           C
ATOM     29  O3'  DT D   2     170.046 161.312 143.239  1.00 56.21           O
ATOM     30  C2'  DT D   2     167.950 160.183 143.722  1.00 56.21           C
ATOM     31  C1'  DT D   2     167.118 160.874 142.665  1.00 56.21           C
ATOM     32  N1   DT D   2     165.663 160.625 142.790  1.00 56.21           N
ATOM     33  C2   DT D   2     164.787 161.668 142.631  1.00 56.21           C
ATOM     34  O2   DT D   2     165.142 162.807 142.394  1.00 56.21           O
ATOM     35  N3   DT D   2     163.471 161.332 142.759  1.00 56.21           N
ATOM     36  C4   DT D   2     162.955 160.085 143.024  1.00 56.21           C
ATOM     37  O4   DT D   2     161.753 159.880 143.123  1.00 56.21           O
ATOM     38  C5   DT D   2     163.923 159.037 143.179  1.00 56.21           C
ATOM     39  C7   DT D   2     163.475 157.639 143.468  1.00 56.21           C
ATOM     40  C6   DT D   2     165.218 159.353 143.056  1.00 56.21           C
"""
t_cif_str = """\
data_7C9C
_cell.angle_alpha                  90.00
_cell.angle_alpha_esd              ?
_cell.angle_beta                   90.00
_cell.angle_beta_esd               ?
_cell.angle_gamma                  90.00
_cell.angle_gamma_esd              ?
_cell.entry_id                     7C9C
_cell.details                      ?
_cell.formula_units_Z              ?
_cell.length_a                     1.00
_cell.length_a_esd                 ?
_cell.length_b                     1.00
_cell.length_b_esd                 ?
_cell.length_c                     1.00
_cell.length_c_esd                 ?
_cell.volume                       ?
_cell.volume_esd                   ?
_cell.Z_PDB                        ?
#
_symmetry.entry_id                         7C9C
_symmetry.cell_setting                     ?
_symmetry.Int_Tables_number                1
_symmetry.space_group_name_Hall            ?
_symmetry.space_group_name_H-M             'P 1'
_symmetry.pdbx_full_space_group_name_H-M   ?
#

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1    P  P     . DT  A 1 1   ? 166.659 153.609 135.182 1.00 61.01  ? 1   DT  D P     1
ATOM   2    O  OP1   . DT  A 1 1   ? 166.715 152.154 135.422 1.00 61.01  ? 1   DT  D OP1   1
ATOM   3    O  OP2   . DT  A 1 1   ? 166.269 154.124 133.851 1.00 61.01  ? 1   DT  D OP2   1
ATOM   4    O  "O5'" . DT  A 1 1   ? 165.696 154.272 136.272 1.00 61.01  ? 1   DT  D "O5'" 1
ATOM   5    C  "C5'" . DT  A 1 1   ? 166.016 154.187 137.657 1.00 61.01  ? 1   DT  D "C5'" 1
ATOM   6    C  "C4'" . DT  A 1 1   ? 166.635 155.480 138.154 1.00 61.01  ? 1   DT  D "C4'" 1
ATOM   7    O  "O4'" . DT  A 1 1   ? 165.721 156.583 137.917 1.00 61.01  ? 1   DT  D "O4'" 1
ATOM   8    C  "C3'" . DT  A 1 1   ? 166.950 155.515 139.650 1.00 61.01  ? 1   DT  D "C3'" 1
ATOM   9    O  "O3'" . DT  A 1 1   ? 168.189 156.162 139.861 1.00 61.01  ? 1   DT  D "O3'" 1
ATOM   10   C  "C2'" . DT  A 1 1   ? 165.801 156.336 140.211 1.00 61.01  ? 1   DT  D "C2'" 1
ATOM   11   C  "C1'" . DT  A 1 1   ? 165.600 157.334 139.096 1.00 61.01  ? 1   DT  D "C1'" 1
ATOM   12   N  N1    . DT  A 1 1   ? 164.274 157.975 139.127 1.00 61.01  ? 1   DT  D N1    1
ATOM   13   C  C2    . DT  A 1 1   ? 164.180 159.327 138.934 1.00 61.01  ? 1   DT  D C2    1
ATOM   14   O  O2    . DT  A 1 1   ? 165.141 160.036 138.720 1.00 61.01  ? 1   DT  D O2    1
ATOM   15   N  N3    . DT  A 1 1   ? 162.912 159.825 138.994 1.00 61.01  ? 1   DT  D N3    1
ATOM   16   C  C4    . DT  A 1 1   ? 161.753 159.119 139.229 1.00 61.01  ? 1   DT  D C4    1
ATOM   17   O  O4    . DT  A 1 1   ? 160.654 159.657 139.264 1.00 61.01  ? 1   DT  D O4    1
ATOM   18   C  C5    . DT  A 1 1   ? 161.924 157.705 139.429 1.00 61.01  ? 1   DT  D C5    1
ATOM   19   C  C7    . DT  A 1 1   ? 160.735 156.834 139.691 1.00 61.01  ? 1   DT  D C7    1
ATOM   20   C  C6    . DT  A 1 1   ? 163.163 157.207 139.374 1.00 61.01  ? 1   DT  D C6    1
ATOM   21   P  P     . DT  A 1 2   ? 168.926 156.085 141.287 1.00 56.21  ? 2   DT  D P     1
ATOM   22   O  OP1   . DT  A 1 2   ? 170.334 155.745 141.021 1.00 56.21  ? 2   DT  D OP1   1
ATOM   23   O  OP2   . DT  A 1 2   ? 168.132 155.242 142.199 1.00 56.21  ? 2   DT  D OP2   1
ATOM   24   O  "O5'" . DT  A 1 2   ? 168.871 157.586 141.817 1.00 56.21  ? 2   DT  D "O5'" 1
ATOM   25   C  "C5'" . DT  A 1 2   ? 169.222 158.639 140.942 1.00 56.21  ? 2   DT  D "C5'" 1
ATOM   26   C  "C4'" . DT  A 1 2   ? 168.956 159.998 141.561 1.00 56.21  ? 2   DT  D "C4'" 1
ATOM   27   O  "O4'" . DT  A 1 2   ? 167.552 160.326 141.452 1.00 56.21  ? 2   DT  D "O4'" 1
ATOM   28   C  "C3'" . DT  A 1 2   ? 169.313 160.123 143.038 1.00 56.21  ? 2   DT  D "C3'" 1
ATOM   29   O  "O3'" . DT  A 1 2   ? 170.046 161.312 143.239 1.00 56.21  ? 2   DT  D "O3'" 1
ATOM   30   C  "C2'" . DT  A 1 2   ? 167.950 160.183 143.722 1.00 56.21  ? 2   DT  D "C2'" 1
ATOM   31   C  "C1'" . DT  A 1 2   ? 167.118 160.874 142.665 1.00 56.21  ? 2   DT  D "C1'" 1
ATOM   32   N  N1    . DT  A 1 2   ? 165.663 160.625 142.790 1.00 56.21  ? 2   DT  D N1    1
ATOM   33   C  C2    . DT  A 1 2   ? 164.787 161.668 142.631 1.00 56.21  ? 2   DT  D C2    1
ATOM   34   O  O2    . DT  A 1 2   ? 165.142 162.807 142.394 1.00 56.21  ? 2   DT  D O2    1
ATOM   35   N  N3    . DT  A 1 2   ? 163.471 161.332 142.759 1.00 56.21  ? 2   DT  D N3    1
ATOM   36   C  C4    . DT  A 1 2   ? 162.955 160.085 143.024 1.00 56.21  ? 2   DT  D C4    1
ATOM   37   O  O4    . DT  A 1 2   ? 161.753 159.880 143.123 1.00 56.21  ? 2   DT  D O4    1
ATOM   38   C  C5    . DT  A 1 2   ? 163.923 159.037 143.179 1.00 56.21  ? 2   DT  D C5    1
ATOM   39   C  C7    . DT  A 1 2   ? 163.475 157.639 143.468 1.00 56.21  ? 2   DT  D C7    1
ATOM   40   C  C6    . DT  A 1 2   ? 165.218 159.353 143.056 1.00 56.21  ? 2   DT  D C6    1
"""

def exercise():
  pdb_inp = iotbx.pdb.input(lines=t_pdb_str, source_info=None)
  cif_inp = iotbx.pdb.input(lines=t_cif_str, source_info=None)
  pdb_cs = pdb_inp.crystal_symmetry()
  cif_cs = cif_inp.crystal_symmetry()
  # print(pdb_cs.unit_cell())
  # print(cif_cs.unit_cell())
  assert pdb_cs.unit_cell() == cif_cs.unit_cell()

if (__name__ == "__main__"):
  exercise()
