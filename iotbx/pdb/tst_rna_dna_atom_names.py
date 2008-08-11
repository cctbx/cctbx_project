import iotbx.pdb

def exercise_rna_dna_atom_names():
  aliases = [
    ("1D2", " H21", "G DG"),
    ("1D2'", " H2'", "ANY"),
    ("1D2*", " H2'", "ANY"),
    ("1D4", " H41", "C DC"),
    ("1D5'", " H5'", "ANY"),
    ("1D5*", " H5'", "ANY"),
    ("1D5M", " H71", "DT"),
    ("1D6", " H61", "A DA"),
    ("1H2", " H21", "G DG"),
    ("1H2'", " H2'", "ANY"),
    ("1H2*", " H2'", "ANY"),
    ("1H4", " H41", "C DC"),
    ("1H5'", " H5'", "ANY"),
    ("1H5*", " H5'", "ANY"),
    ("1H5M", " H71", "DT"),
    ("1H6", " H61", "A DA"),
    ("2D2", " H22", "G DG"),
    ("2D2'", "H2''", "DA DC DG DT"),
    ("2D2*", "H2''", "DA DC DG DT"),
    ("2D4", " H42", "C DC"),
    ("2D5'", "H5''", "ANY"),
    ("2D5*", "H5''", "ANY"),
    ("2D5M", " H72", "DT"),
    ("2D6", " H62", "A DA"),
    ("2DO'", "HO2'", "A C G U"),
    ("2DO*", "HO2'", "A C G U"),
    ("2DOP", "HOP2", "ANY"),
    ("2H2", " H22", "G DG"),
    ("2H2'", "H2''", "DA DC DG DT"),
    ("2H2*", "H2''", "DA DC DG DT"),
    ("2H4", " H42", "C DC"),
    ("2H5'", "H5''", "ANY"),
    ("2H5*", "H5''", "ANY"),
    ("2H5M", " H72", "DT"),
    ("2H6", " H62", "A DA"),
    ("2HO'", "HO2'", "A C G U"),
    ("2HO*", "HO2'", "A C G U"),
    ("2HOP", "HOP2", "ANY"),
    ("3D5M", " H73", "DT"),
    ("3DOP", "HOP3", "ANY"),
    ("3H5M", " H73", "DT"),
    ("3HOP", "HOP3", "ANY"),
    ("C1'", " C1'", "ANY"),
    ("C1*", " C1'", "ANY"),
    ("C2", " C2 ", "ANY"),
    ("C2'", " C2'", "ANY"),
    ("C2*", " C2'", "ANY"),
    ("C3'", " C3'", "ANY"),
    ("C3*", " C3'", "ANY"),
    ("C4", " C4 ", "ANY"),
    ("C4'", " C4'", "ANY"),
    ("C4*", " C4'", "ANY"),
    ("C5", " C5 ", "ANY"),
    ("C5'", " C5'", "ANY"),
    ("C5*", " C5'", "ANY"),
    ("C5M", " C7 ", "DT"),
    ("C6", " C6 ", "ANY"),
    ("C7", " C7 ", "DT"),
    ("C8", " C8 ", "A G DA DG"),
    ("D1", " H1 ", "G DG"),
    ("D1'", " H1'", "ANY"),
    ("D1*", " H1'", "ANY"),
    ("D2", " H2 ", "A DA"),
    ("D2'", " H2'", "ANY"),
    ("D2*", " H2'", "ANY"),
    ("D2''", "H2''", "DA DC DG DT"),
    ("D2'1", " H2'", "ANY"),
    ("D2*1", " H2'", "ANY"),
    ("D2'2", "H2''", "DA DC DG DT"),
    ("D2*2", "H2''", "DA DC DG DT"),
    ("D21", " H21", "G DG"),
    ("D22", " H22", "G DG"),
    ("D3", " H3 ", "U DT"),
    ("D3'", " H3'", "ANY"),
    ("D3*", " H3'", "ANY"),
    ("D3T", "HO3'", "ANY"),
    ("D4'", " H4'", "ANY"),
    ("D4*", " H4'", "ANY"),
    ("D41", " H41", "C DC"),
    ("D42", " H42", "C DC"),
    ("D5", " H5 ", "C U DC"),
    ("D5'", " H5'", "ANY"),
    ("D5*", "HO5'", "ANY"),
    ("D5''", "H5''", "ANY"),
    ("D5'1", " H5'", "ANY"),
    ("D5*1", " H5'", "ANY"),
    ("D5'2", "H5''", "ANY"),
    ("D5*2", "H5''", "ANY"),
    ("D5M1", " H71", "DT"),
    ("D5M2", " H72", "DT"),
    ("D5M3", " H73", "DT"),
    ("D5T", "HO5'", "ANY"),
    ("D6", " H6 ", "C U DC DT"),
    ("D61", " H61", "A DA"),
    ("D62", " H62", "A DA"),
    ("D71", " H71", "DT"),
    ("D72", " H72", "DT"),
    ("D73", " H73", "DT"),
    ("D8", " H8 ", "A G DA DG"),
    ("DO2'", "HO2'", "A C G U"),
    ("DO2*", "HO2'", "A C G U"),
    ("H1", " H1 ", "G DG"),
    ("H1'", " H1'", "ANY"),
    ("H1*", " H1'", "ANY"),
    ("H2", " H2 ", "A DA"),
    ("H2'", " H2'", "ANY"),
    ("H2*", " H2'", "ANY"),
    ("H2''", "H2''", "DA DC DG DT"),
    ("H2'1", " H2'", "ANY"),
    ("H2*1", " H2'", "ANY"),
    ("H2'2", "H2''", "DA DC DG DT"),
    ("H2*2", "H2''", "DA DC DG DT"),
    ("H21", " H21", "G DG"),
    ("H22", " H22", "G DG"),
    ("H3", " H3 ", "U DT"),
    ("H3'", " H3'", "ANY"),
    ("H3*", " H3'", "ANY"),
    ("H3T", "HO3'", "ANY"),
    ("H4'", " H4'", "ANY"),
    ("H4*", " H4'", "ANY"),
    ("H41", " H41", "C DC"),
    ("H42", " H42", "C DC"),
    ("H5", " H5 ", "C U DC"),
    ("H5'", " H5'", "ANY"),
    ("H5*", "HO5'", "ANY"),
    ("H5''", "H5''", "ANY"),
    ("H5'1", " H5'", "ANY"),
    ("H5*1", " H5'", "ANY"),
    ("H5'2", "H5''", "ANY"),
    ("H5*2", "H5''", "ANY"),
    ("H5M1", " H71", "DT"),
    ("H5M2", " H72", "DT"),
    ("H5M3", " H73", "DT"),
    ("H5T", "HO5'", "ANY"),
    ("H6", " H6 ", "C U DC DT"),
    ("H61", " H61", "A DA"),
    ("H62", " H62", "A DA"),
    ("H71", " H71", "DT"),
    ("H72", " H72", "DT"),
    ("H73", " H73", "DT"),
    ("H8", " H8 ", "A G DA DG"),
    ("HO2'", "HO2'", "A C G U"),
    ("HO2*", "HO2'", "A C G U"),
    ("HO3'", "HO3'", "ANY"),
    ("HO3*", "HO3'", "ANY"),
    ("HO5'", "HO5'", "ANY"),
    ("HO5*", "HO5'", "ANY"),
    ("HOP2", "HOP2", "ANY"),
    ("HOP3", "HOP3", "ANY"),
    ("N1", " N1 ", "ANY"),
    ("N2", " N2 ", "G DG"),
    ("N3", " N3 ", "ANY"),
    ("N4", " N4 ", "C DC"),
    ("N6", " N6 ", "A DA"),
    ("N7", " N7 ", "A G DA DG"),
    ("N9", " N9 ", "A G DA DG"),
    ("O1P", " OP1", "ANY"),
    ("O2", " O2 ", "C U DC DT"),
    ("O2'", " O2'", "A C G U"),
    ("O2*", " O2'", "A C G U"),
    ("O2P", " OP2", "ANY"),
    ("O3'", " O3'", "ANY"),
    ("O3*", " O3'", "ANY"),
    ("O3P", " OP3", "ANY"),
    ("O3T", " OP3", "ANY"),
    ("O4", " O4 ", "U DT"),
    ("O4'", " O4'", "ANY"),
    ("O4*", " O4'", "ANY"),
    ("O5'", " O5'", "ANY"),
    ("O5*", " O5'", "ANY"),
    ("O5T", " OP3", "ANY"),
    ("O6", " O6 ", "G DG"),
    ("OP1", " OP1", "ANY"),
    ("OP2", " OP2", "ANY"),
    ("OP3", " OP3", "ANY"),
    ("P", " P  ", "ANY")]
  for atom_name in [None, ""]:
    info = iotbx.pdb.rna_dna_atom_names_info(atom_name=atom_name)
    assert info.reference_name is None
    assert info.compatible_residue_names() == "None"
    assert not info.is_compatible_with(residue_name="")
    assert not info.is_hydrogen()
    assert not info.is_deuterium()
    assert not info.is_o2prime()
    assert not info.is_ho2prime()
    assert not info.is_h2primeprime()
    assert not info.is_in_phosphate_group()
    assert not info.is_op3_or_hop3()
    assert not info.is_ho5prime()
    assert not info.is_ho3prime()
    assert not info.change_h2primeprime_to_ho2prime()
    assert not info.change_ho5prime_to_hop3()
    info.change_to_unknown()
  for a,r,f in aliases:
    for lower in [False, True]:
      if (lower): work_name = "  " + a.lower() + " "
      else:       work_name = a
      info = iotbx.pdb.rna_dna_atom_names_info(atom_name=work_name)
      assert info.reference_name == r
      assert info.compatible_residue_names() == f
      assert not info.is_compatible_with(residue_name="")
      assert not info.is_compatible_with(residue_name="D")
      for n in f.replace("ANY", "A C G U DA DC DG DT").split():
        assert info.is_compatible_with(residue_name=n)
        assert not info.is_compatible_with(residue_name=n+"X")
      assert info.is_hydrogen() == (a.find("H") >= 0 or a.find("D") >= 0)
      assert info.is_deuterium() == (a.find("D") >= 0)
      assert info.is_o2prime() == (info.reference_name == " O2'")
      assert info.is_ho2prime() == (info.reference_name == "HO2'")
      assert info.is_h2primeprime() == (info.reference_name == "H2''")
      assert info.is_in_phosphate_group() == (info.reference_name in [
        " P  ", " OP1", " OP2", "HOP2", " OP3", "HOP3"])
      assert info.is_op3_or_hop3() == (info.reference_name in [
        " OP3", "HOP3"])
      if (not info.is_in_phosphate_group()):
        assert not info.is_op3_or_hop3()
      elif (info.is_op3_or_hop3()):
        assert info.is_in_phosphate_group()
      assert info.is_ho5prime() == (info.reference_name == "HO5'")
      assert info.is_ho3prime() == (info.reference_name == "HO3'")
      if (not info.is_h2primeprime()):
        assert not info.change_h2primeprime_to_ho2prime()
      if (not info.is_ho5prime()):
        assert not info.change_ho5prime_to_hop3()
      if (info.is_h2primeprime()):
        assert info.change_h2primeprime_to_ho2prime()
        assert info.reference_name == "HO2'"
        assert info.is_hydrogen()
        assert info.is_deuterium() == (a.find("D") >= 0)
        assert info.compatible_residue_names() == "A C G U"
        assert not info.is_o2prime()
        assert info.is_ho2prime()
        assert not info.is_h2primeprime()
        assert not info.is_in_phosphate_group()
        assert not info.is_op3_or_hop3()
        assert not info.is_ho5prime()
        assert not info.is_ho3prime()
      elif (info.is_ho5prime()):
        assert info.change_ho5prime_to_hop3()
        assert info.reference_name == "HOP3"
        assert info.is_hydrogen()
        assert info.is_deuterium() == (a.find("D") >= 0)
        assert info.compatible_residue_names() == "ANY"
        assert not info.is_o2prime()
        assert not info.is_ho2prime()
        assert not info.is_h2primeprime()
        assert info.is_in_phosphate_group()
        assert info.is_op3_or_hop3()
        assert not info.is_ho5prime()
        assert not info.is_ho3prime()
      info.change_to_unknown()
      assert info.reference_name is None
      assert info.compatible_residue_names() == "None"
  for a,r,f in aliases:
    info = iotbx.pdb.rna_dna_atom_names_info(atom_name=a+"X")
    assert info.reference_name is None
    assert info.compatible_residue_names() == "None"
  for a,r,f in aliases:
    for i in xrange(len(a)):
      info = iotbx.pdb.rna_dna_atom_names_info(atom_name=a[:i]+"X"+a[i+1:])
      assert info.reference_name is None
      assert info.compatible_residue_names() == "None"

mon_lib_rna_dna_names = {
  "AR": """
    P O1P O2P O5* C5* H5*1 H5*2 C4* H4* O4* C1* H1* N9 C8 H8 N7 C5 C4 N3
    C2 H2 N1 C6 N6 H61 H62 C2* H2* O2* HO2* C3* H3* O3*""".split(),
  "CR": """
    P O1P O2P O5* C5* H5*1 H5*2 C4* H4* O4* C1* H1* N1 C2 O2 N3 C4 N4 H41
    H42 C5 H5 C6 H6 C2* H2* O2* HO2* C3* H3* O3*""".split(),
  "GR": """
    P O1P O2P O5* C5* H5*1 H5*2 C4* H4* O4* C1* H1* N9 C8 H8 N7 C5 C4 N3 C2
    N2 H21 H22 N1 H1 C6 O6 C2* H2* O2* HO2* C3* H3* O3*""".split(),
  "UR": """
    P O1P O2P O5* C5* H5*1 H5*2 C4* H4* O4* C1* H1* N1 C2 O2 N3 H3 C4 O4 C5
    H5 C6 H6 C2* H2* O2* HO2* C3* H3* O3*""".split(),
  "AD": """
    P O1P O2P O5* C5* H5*1 H5*2 C4* H4* O4* C1* H1* N9 C8 H8 N7 C5 C4 N3 C2
    H2 N1 C6 N6 H61 H62 C2* H2*1 H2*2 C3* H3* O3*""".split(),
  "CD": """
    P O1P O2P O5* C5* H5*1 H5*2 C4* H4* O4* C1* H1* N1 C2 O2 N3 C4 N4 H41
    H42 C5 H5 C6 H6 C2* H2*1 H2*2 C3* H3* O3*""".split(),
  "GD": """
    P O1P O2P O5* C5* H5*1 H5*2 C4* H4* O4* C1* H1* N9 C8 H8 N7 C5 C4 N3 C2
    N2 H21 H22 N1 H1 C6 O6 C2* H2*1 H2*2 C3* H3* O3*""".split(),
  "TD": """
    P O1P O2P O5* C5* H5*1 H5*2 C4* H4* O4* C1* H1* N1 C2 O2 N3 H3 C4 O4 C5
    C5M H5M1 H5M2 H5M3 C6 H6 C2* H2*1 H2*2 C3* H3* O3*""".split()
}

def exercise_mon_lib_names():
  n_problems = 0
  for mon_code,mon_names in mon_lib_rna_dna_names.items():
    mon_names = mon_names + ["O3T", "HO5*", "HO3*"]
    q = ""
    if (mon_code[0] in ["A", "C", "G"]): q = "?"
    def interpretation(mon_names):
      result = iotbx.pdb.rna_dna_atom_names_interpretation(
        residue_name=q+mon_code[0],
        atom_names=mon_names)
      if (mon_code[1] == "R"):
        assert result.residue_name == mon_code[0]
      else:
        assert result.residue_name == "D"+mon_code[0]
      return result
    #
    interpreted = interpretation(mon_names=mon_names)
    interp_mon_names = interpreted.mon_lib_names()
    for inp,out in zip(mon_names, interp_mon_names):
      if (inp == "HO5*"):
        assert out == "HOP3" # added to monomer library
      else:
        assert out == inp
    #
    mon_names_without_phosphate = []
    for mon_name,info in zip(mon_names, interpreted.infos):
      if (info.is_in_phosphate_group()): continue
      mon_names_without_phosphate.append(mon_name)
    interpreted = interpretation(mon_names=mon_names_without_phosphate)
    interp_mon_names = interpreted.mon_lib_names()
    for inp,out in zip(mon_names_without_phosphate, interp_mon_names):
      assert out == inp

cns_rna_dna_names = {
  "GUA": """
    P O1P O2P O5' C5' H5' H5'' C4' H4' O4' C1' H1' N9 C4 N3 C2 N2 H21 H22
    N1 H1 C6 O6 C5 N7 C8 H8 C2' H2' O2' HO2' C3' H3' O3'""".split(),
  "ADE": """
    P O1P O2P O5' C5' H5' H5'' C4' H4' O4' C1' H1' N9 C4 N3 C2 H2 N1 C6
    N6 H61 H62 C5 N7 C8 H8 C2' H2' O2' HO2' C3' H3' O3'""".split(),
  "CYT": """
    P O1P O2P O5' C5' H5' H5'' C4' H4' O4' C1' H1' N1 C6 H6 C2 O2 N3 C4
    N4 H41 H42 C5 H5 C2' H2' O2' HO2' C3' H3' O3'""".split(),
  "THY": """
    P O1P O2P O5' C5' H5' H5'' C4' H4' O4' C1' H1' N1 C6 H6 C2 O2 N3 H3 C4
    O4 C5 C7 H71 H72 H73 C2' H2' O2' HO2' C3' H3' O3'""".split(),
  "URI": """
    P O1P O2P O5' C5' H5' H5'' C4' H4' O4' C1' H1' N1 C6 H6 C2 O2 N3 H3 C4
    O4 C5 H5 C2' H2' O2' HO2' C3' H3' O3'""".split()
}

def apply_cns_deox(atom_names):
  result = list(atom_names)
  result.remove("O2'")
  result.remove("HO2'")
  result.append("H2''")
  return result

def apply_cns_5pho(atom_names): # 5-terminus (with phosphate)
  result = list(atom_names)
  result.append("H5T") # bond O5T
  result.append("O5T") # bond P
  return result

def apply_cns_5ter(atom_names): # 5-terminus (without phosphate)
  result = list(atom_names)
  result.remove("P")
  result.remove("O1P")
  result.remove("O2P")
  result.append("H5T") # bond O5'
  return result

def apply_cns_3ter(atom_names): # 3-terminus (without phosphate)
  result = list(atom_names)
  result.append("H3T") # bond O3'
  return result

def exercise_cns_names():
  for residue_name,atom_names_base in cns_rna_dna_names.items():
    for deox in [False, True]:
      q = ""
      if (residue_name == "URI"):
        if (deox): continue
      elif (residue_name == "THY"):
        if (not deox): continue
      else:
        q = "?"
      if (deox):
        atom_names_base = apply_cns_deox(atom_names=atom_names_base)
      def do_nothing(atom_names): return atom_names
      for apply3 in [do_nothing, apply_cns_3ter]:
        atom_names = apply3(atom_names=atom_names_base)
        for apply5 in [do_nothing, apply_cns_5pho, apply_cns_5ter]:
          interpreted = iotbx.pdb.rna_dna_atom_names_interpretation(
            residue_name=q+residue_name[0],
            atom_names=apply5(atom_names=atom_names))
          if (not deox):
            assert interpreted.residue_name == residue_name[0]
          else:
            assert interpreted.residue_name == "D"+residue_name[0]
          assert interpreted.have_phosphate == (apply5 != apply_cns_5ter)
          assert interpreted.have_op3_or_hop3 == (apply5 == apply_cns_5pho)
          assert interpreted.have_ho3prime == (apply3 == apply_cns_3ter)
          assert interpreted.n_unexpected == 0
          assert interpreted.n_expected == len(interpreted.atom_names)
          for atom_name,info in zip(interpreted.atom_names, interpreted.infos):
            if (atom_name == "O1P"):
              assert info.reference_name == " OP1"
            elif (atom_name == "O2P"):
              assert info.reference_name == " OP2"
            elif (atom_name == "O5T"):
              assert info.reference_name == " OP3"
            elif (atom_name == "H5T"):
              if (apply5 == apply_cns_5ter):
                assert info.reference_name == "HO5'"
              else:
                assert info.reference_name == "HOP3"
            else:
              info.reference_name.strip() == atom_name

def exercise():
  exercise_rna_dna_atom_names()
  exercise_mon_lib_names()
  exercise_cns_names()
  print "OK"

if (__name__ == "__main__"):
  exercise()
