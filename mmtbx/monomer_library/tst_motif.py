from mmtbx.monomer_library import server
from libtbx.test_utils import show_diff
from cStringIO import StringIO
import sys

tmp_cif = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
tst tst 'test compound' g 1 2 l

data_comp_tst
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
  tst a b c 1.2
  tst d e f 2.3
  tst . . . .
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
  tst a b c 3.4 4.5
  tst . . . . .
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
  tst a b c 5.6 6.7
  tst . . . . .
loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
  tst a b c d e 1.2 2.3 7
  tst . . . . . . . .
loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
  tst f g h i j k
  tst . . . . . .
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
  tst p1 a 1.2
  tst p1 b 2.3
  tst p1 c 3.4
  tst p2 d 1.2
  tst p2 e 2.3
  tst p2 f 3.4
  tst . . .

###############################################################################

data_link_list
loop_
_chem_link.id
_chem_link.comp_id_1
_chem_link.mod_id_1
_chem_link.group_comp_1
_chem_link.comp_id_2
_chem_link.mod_id_2
_chem_link.group_comp_2
_chem_link.name
  tst_lnk c1 m1 g1 c2 m2 g2 'test link'

data_link_tst_lnk
loop_
_chem_link_bond.link_id
_chem_link_bond.atom_1_comp_id
_chem_link_bond.atom_id_1
_chem_link_bond.atom_2_comp_id
_chem_link_bond.atom_id_2
_chem_link_bond.type
_chem_link_bond.value_dist
_chem_link_bond.value_dist_esd
  tst_lnk 1 A 2 B single 2.3 1.2
  tst_lnk . . . . . . .
loop_
_chem_link_angle.link_id
_chem_link_angle.atom_1_comp_id
_chem_link_angle.atom_id_1
_chem_link_angle.atom_2_comp_id
_chem_link_angle.atom_id_2
_chem_link_angle.atom_3_comp_id
_chem_link_angle.atom_id_3
_chem_link_angle.value_angle
_chem_link_angle.value_angle_esd
  tst_lnk 1 C 2 D 3 E 13.8 1.8
  tst_lnk . . . . . . . .
loop_
_chem_link_tor.link_id
_chem_link_tor.id
_chem_link_tor.atom_1_comp_id
_chem_link_tor.atom_id_1
_chem_link_tor.atom_2_comp_id
_chem_link_tor.atom_id_2
_chem_link_tor.atom_3_comp_id
_chem_link_tor.atom_id_3
_chem_link_tor.atom_4_comp_id
_chem_link_tor.atom_id_4
_chem_link_tor.value_angle
_chem_link_tor.value_angle_esd
_chem_link_tor.period
  tst_lnk u 1 F 2 G 3 H 4 I 1.2 2.5 3
  tst_lnk . . . . . . . . . . . .
loop_
_chem_link_chir.link_id
_chem_link_chir.id
_chem_link_chir.atom_centre_comp_id
_chem_link_chir.atom_id_centre
_chem_link_chir.atom_1_comp_id
_chem_link_chir.atom_id_1
_chem_link_chir.atom_2_comp_id
_chem_link_chir.atom_id_2
_chem_link_chir.atom_3_comp_id
_chem_link_chir.atom_id_3
_chem_link_chir.volume_sign
  tst_lnk v 1 J 2 K 3 L 4 M positiv
  tst_lnk . . . . . . . . . .
loop_
_chem_link_plane.link_id
_chem_link_plane.plane_id
_chem_link_plane.atom_comp_id
_chem_link_plane.atom_id
_chem_link_plane.dist_esd
  tst_lnk p 1 N 0.1
  tst_lnk p 2 O 0.2
  tst_lnk p 3 P 0.3
  tst_lnk p 4 Q 0.4
  tst_lnk q 5 R 0.5
  tst_lnk q 6 S 0.6
  tst_lnk q 7 T 0.7
  tst_lnk q 8 U 0.8
  tst_lnk . . . .

###############################################################################

data_mod_list
loop_
_chem_mod.id
_chem_mod.name
_chem_mod.comp_id
_chem_mod.group_id
  tst_mod 'test modification' c g

data_mod_tst_mod
loop_
_chem_mod_atom.mod_id
_chem_mod_atom.function
_chem_mod_atom.atom_id
_chem_mod_atom.new_atom_id
_chem_mod_atom.new_type_symbol
_chem_mod_atom.new_type_energy
_chem_mod_atom.new_partial_charge
  tst_mod add    . A B C 1.2
  tst_mod add    . . . . .
  tst_mod change D E F G 2.3
  tst_mod change . . . . .
  tst_mod delete H . . . .
  tst_mod delete . . . . .
loop_
_chem_mod_bond.mod_id
_chem_mod_bond.function
_chem_mod_bond.atom_id_1
_chem_mod_bond.atom_id_2
_chem_mod_bond.new_type
_chem_mod_bond.new_value_dist
_chem_mod_bond.new_value_dist_esd
  tst_mod add    A B x 1.2 2.3
  tst_mod add    . . . . .
  tst_mod change C D y 3.4 4.5
  tst_mod change . . . . .
  tst_mod delete E F z 5.6 6.7
  tst_mod delete . . . . .
loop_
_chem_mod_angle.mod_id
_chem_mod_angle.function
_chem_mod_angle.atom_id_1
_chem_mod_angle.atom_id_2
_chem_mod_angle.atom_id_3
_chem_mod_angle.new_value_angle
_chem_mod_angle.new_value_angle_esd
  tst_mod add A B C 1.2 2.3
  tst_mod add . . . . .
  tst_mod change D E F 1.2 2.3
  tst_mod change . . . . .
  tst_mod delete G H I . .
  tst_mod delete . . . . .
loop_
_chem_mod_tor.mod_id
_chem_mod_tor.function
_chem_mod_tor.id
_chem_mod_tor.atom_id_1
_chem_mod_tor.atom_id_2
_chem_mod_tor.atom_id_3
_chem_mod_tor.atom_id_4
_chem_mod_tor.new_value_angle
_chem_mod_tor.new_value_angle_esd
_chem_mod_tor.new_period
  tst_mod add    x A B C D 1.2 2.3 5
  tst_mod add    . . . . . . . .
  tst_mod change y E F G H 3.4 5.5 7
  tst_mod change . . . . . . . .
  tst_mod delete . I J K L . . .
  tst_mod delete . . . . . . . .
loop_
_chem_mod_chir.mod_id
_chem_mod_chir.function
_chem_mod_chir.id
_chem_mod_chir.atom_id_centre
_chem_mod_chir.atom_id_1
_chem_mod_chir.atom_id_2
_chem_mod_chir.atom_id_3
_chem_mod_chir.new_volume_sign
  tst_mod add    x A B C D negativ
  tst_mod add    . . . . . .
  tst_mod change y E F G H positiv
  tst_mod change . . . . . .
  tst_mod delete . I J K L .
  tst_mod delete . . . . . .
loop_
_chem_mod_plane_atom.mod_id
_chem_mod_plane_atom.function
_chem_mod_plane_atom.plane_id
_chem_mod_plane_atom.atom_id
_chem_mod_plane_atom.new_dist_esd
  tst_mod add    p A 0.1
  tst_mod add    p . .
  tst_mod change p B 0.2
  tst_mod change p . .
  tst_mod delete p C 0.3
  tst_mod delete p . .
  tst_mod add    q D 0.4
  tst_mod add    q . .
  tst_mod change q E 0.5
  tst_mod change q . .
  tst_mod delete q F 0.6
  tst_mod delete q . .
"""

expected_out_tst_comp = """\
geometry_restraints.motif {
  id = "tst"
  description = "test compound"
  info = "file: tmp.cif"
  atom = [name scattering_type nonbonded_type partial_charge]
  atom = "a" "b" "c" 1.2
  atom = "d" "e" "f" 2.3
  atom = "" "" "" 0
  bond = [atom_name*2 type distance_ideal weight id]
  bond = "a" "b" "c" 3.4 0.0493827 ""
  bond = "" "" "" 0 0 ""
  angle = [atom_name*3 angle_ideal weight id]
  angle = "a" "b" "c" 5.6 0.0222767 ""
  angle = "" "" "" 0 0 ""
  dihedral = [atom_name*4 angle_ideal weight periodicity id]
  dihedral = "b" "c" "d" "e" 1.2 0.189036 7 "a"
  dihedral = "" "" "" "" 0 0 0 ""
  chirality = [atom_name*4 volume_sign both_signs volume_ideal weight id]
  chirality = "g" "h" "i" "j" "k" False 0 0 "f"
  chirality = "" "" "" "" "" False 0 0 ""
  planarity {
    id = "p1"
    atom = [name weight]
    atom = "a" 0.694444
    atom = "b" 0.189036
    atom = "c" 0.0865052
  }
  planarity {
    id = "p2"
    atom = [name weight]
    atom = "d" 0.694444
    atom = "e" 0.189036
    atom = "f" 0.0865052
  }
  planarity {
    id = ""
    atom = [name weight]
    atom = "" 0
  }
}
"""

expected_out_tst_lnk = """\
geometry_restraints.motif_manipulation {
  id = "tst_lnk"
  description = "test link"
  info = "file: tmp.cif"
  bond = add [(motif_id atom_name)*2 type distance_ideal weight id]
  bond = add "1" "A" "2" "B" 2.3 0.694444 ""
  bond = add "" "" "" "" 0 0 ""
  angle = add [(motif_id atom_name)*3 type angle_ideal weight id]
  angle = add "1" "C" "2" "D" "3" "E" 13.8 0.308642 ""
  angle = add "" "" "" "" "" "" 0 0 ""
  dihedral = add [(motif_id atom_name)*4 angle_ideal weight periodicity id]
  dihedral = add "1" "F" "2" "G" "3" "H" "4" "I" 1.2 0.16 3 "u"
  dihedral = add "" "" "" "" "" "" "" "" 0 0 0 ""
  chirality = add [(motif_id atom_name)*4 \\
                   volume_sign volume_ideal weight id]
  chirality = add "1" "J" "2" "K" "3" "L" "4" "M" \\
                  "positiv" 0 0 "v"
  chirality = add "" "" "" "" "" "" "" "" \\
                  "" 0 0 ""
  planarity {
    action = add
    motif_id = ""
    id = "p"
    atom = [motif_id name weight]
    atom = "1" "N" 100
    atom = "2" "O" 25
    atom = "3" "P" 11.1111
    atom = "4" "Q" 6.25
  }
  planarity {
    action = add
    motif_id = ""
    id = "q"
    atom = [motif_id name weight]
    atom = "5" "R" 4
    atom = "6" "S" 2.77778
    atom = "7" "T" 2.04082
    atom = "8" "U" 1.5625
  }
  planarity {
    action = add
    motif_id = ""
    id = ""
    atom = [motif_id name weight]
    atom = "" "" 0
  }
}
"""

expected_out_tst_mod = """\
geometry_restraints.motif_manipulation {
  id = "tst_mod"
  description = "test modification"
  info = "file: tmp.cif"
  atom = add [motif_id name scattering_type nonbonded_type partial_charge]
  atom = add "" "A" "B" "C" 1.2
  atom = add "" "" "" "" 0.0
  atom = change [motif_id motif_atom_name \\
                 name scattering_type nonbonded_type partial_charge]
  atom = change "" "D" \\
                "E" "F" "G" 2.3
  atom = change "" "" \\
                "" "" "" None
  atom = delete [motif_id motif_atom_name]
  atom = delete "" "H"
  atom = delete "" ""
  bond = add [(motif_id atom_name)*2 type distance_ideal weight id]
  bond = add "" "A" "" "B" 1.2 0.189036 ""
  bond = add "" "" "" "" 0 0 ""
  bond = change [(motif_id atom_name)*2 type distance_ideal weight id]
  bond = change "" "C" "" "D" 3.4 0.0493827 ""
  bond = change "" "" "" "" None None ""
  bond = delete [(motif_id atom_name)*2]
  bond = delete "" "E" "" "F"
  bond = delete "" "" "" ""
  angle = delete [(motif_id atom_name)*3]
  angle = delete "" "A" "" "B" "" "C"
  angle = delete "" "" "" "" "" ""
  angle = delete "" "D" "" "E" "" "F"
  angle = delete "" "" "" "" "" ""
  angle = delete "" "G" "" "H" "" "I"
  angle = delete "" "" "" "" "" ""
  dihedral = add [(motif_id atom_name)*4 angle_ideal weight periodicity id]
  dihedral = add "" "A" "" "B" "" "C" "" "D" 1.2 0.189036 5 "x"
  dihedral = add "" "" "" "" "" "" "" "" 0 0 0 ""
  dihedral = change [(motif_id atom_name)*4 angle_ideal weight periodicity id]
  dihedral = change "" "E" "" "F" "" "G" "" "H" 3.4 0.0330579 7 "y"
  dihedral = change "" "" "" "" "" "" "" "" None None None ""
  dihedral = delete [(motif_id atom_name)*4]
  dihedral = delete "" "I" "" "J" "" "K" "" "L"
  dihedral = delete "" "" "" "" "" "" "" ""
  chirality = add [(motif_id atom_name)*4 \\
                   volume_sign volume_ideal weight id]
  chirality = add "" "A" "" "B" "" "C" "" "D" \\
                  "negativ" 0 0 "x"
  chirality = add "" "" "" "" "" "" "" "" \\
                  "" 0 0 ""
  chirality = change [(motif_id atom_name)*4 \\
                      volume_sign volume_ideal weight id]
  chirality = change "" "E" "" "F" "" "G" "" "H" \\
                     "positiv" None None "y"
  chirality = change "" "" "" "" "" "" "" "" \\
                     "" None None ""
  chirality = delete [(motif_id atom_name)*4]
  chirality = delete "" "I" "" "J" "" "K" "" "L"
  chirality = delete "" "" "" "" "" "" "" ""
  planarity {
    action = change
    motif_id = ""
    id = "p"
    atom = add [motif_id name weight]
    atom = add "" "A" 100
    atom = add "" "" 0
    atom = change [motif_id name weight]
    atom = change "" "B" 25
    atom = change "" "" 0
    atom = delete [motif_id name]
    atom = delete "" "C"
    atom = delete "" ""
  }
  planarity {
    action = change
    motif_id = ""
    id = "q"
    atom = add [motif_id name weight]
    atom = add "" "D" 6.25
    atom = add "" "" 0
    atom = change [motif_id name weight]
    atom = change "" "E" 4
    atom = change "" "" 0
    atom = delete [motif_id name]
    atom = delete "" "F"
    atom = delete "" ""
  }
}
"""

def exercise():
  verbose = "--verbose" in sys.argv[1:]
  list_cif = server.mon_lib_list_cif()
  srv = server.server(list_cif=list_cif)
  open("tmp.cif", "w").write(tmp_cif)
  srv.process_cif(file_name="tmp.cif")
  comp_comp_id = srv.get_comp_comp_id_direct(comp_id="tst")
  motif = comp_comp_id.as_geometry_restraints_motif()
  out = StringIO()
  motif.show(out=out)
  if (verbose): sys.stdout.write(out.getvalue())
  assert not show_diff(out.getvalue(), expected_out_tst_comp)
  for link_link_id in srv.link_link_id_list:
    out = StringIO()
    link_link_id.as_geometry_restraints_motif_manipulation().show(out=out)
    if (verbose): sys.stdout.write(out.getvalue())
    if (link_link_id.chem_link.id == "tst_lnk"):
      assert not show_diff(out.getvalue(), expected_out_tst_lnk)
  for mod_mod_id in srv.mod_mod_id_list:
    out = StringIO()
    mod_mod_id.as_geometry_restraints_motif_manipulation().show(out=out)
    if (verbose): sys.stdout.write(out.getvalue())
    if (mod_mod_id.chem_mod.id == "tst_mod"):
      assert not show_diff(out.getvalue(), expected_out_tst_mod)
  print "OK"

if (__name__ == "__main__"):
  exercise()
