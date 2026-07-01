"""Tests for mmtbx.kinemage.validation"""

from __future__ import absolute_import, division, print_function
import os

# Small PDB with a short helix (3 residues), enough for backbone connections
# and basic validation
pdb_str = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       2.458   1.000   1.000  1.00 10.00           C
ATOM      3  C   ALA A   1       3.009   2.425   1.000  1.00 10.00           C
ATOM      4  O   ALA A   1       2.249   3.390   1.000  1.00 10.00           O
ATOM      5  CB  ALA A   1       2.982   0.231   2.207  1.00 10.00           C
ATOM      6  N   ALA A   2       4.316   2.498   1.000  1.00 10.00           N
ATOM      7  CA  ALA A   2       4.988   3.798   1.000  1.00 10.00           C
ATOM      8  C   ALA A   2       6.504   3.664   1.000  1.00 10.00           C
ATOM      9  O   ALA A   2       7.054   2.562   1.000  1.00 10.00           O
ATOM     10  CB  ALA A   2       4.538   4.637   2.196  1.00 10.00           C
ATOM     11  N   ALA A   3       7.142   4.835   1.000  1.00 10.00           N
ATOM     12  CA  ALA A   3       8.598   4.893   1.000  1.00 10.00           C
ATOM     13  C   ALA A   3       9.108   6.330   1.000  1.00 10.00           C
ATOM     14  O   ALA A   3       8.330   7.286   1.000  1.00 10.00           O
ATOM     15  CB  ALA A   3       9.137   4.148   2.215  1.00 10.00           C
END
"""

# PDB with altlocs to test the altloc bug fix (line 650 in original)
pdb_altloc_str = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       2.458   1.000   1.000  1.00 10.00           C
ATOM      3  C   ALA A   1       3.009   2.425   1.000  1.00 10.00           C
ATOM      4  O   ALA A   1       2.249   3.390   1.000  1.00 10.00           O
ATOM      5  CB AALA A   1       2.982   0.231   2.207  0.60 10.00           C
ATOM      6  CB BALA A   1       2.800   0.400   2.100  0.40 10.00           C
END
"""


# PDB with hets (SO4 ligand), ions (ZN), waters (HOH + WAT),
# and a single-atom SO4 (like in 1nxb) to test het/ion/water handling
pdb_het_str = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       2.458   1.000   1.000  1.00 10.00           C
ATOM      3  C   ALA A   1       3.009   2.425   1.000  1.00 10.00           C
ATOM      4  O   ALA A   1       2.249   3.390   1.000  1.00 10.00           O
ATOM      5  CB  ALA A   1       2.982   0.231   2.207  1.00 10.00           C
HETATM    6 ZN    ZN A 100      10.000  10.000  10.000  1.00 15.00          ZN
HETATM    7  S   SO4 A 200      15.000  15.000  15.000  1.00 20.00           S
HETATM    8  O1  SO4 A 200      16.200  15.500  15.200  1.00 20.00           O
HETATM    9  O2  SO4 A 200      14.200  16.000  15.500  1.00 20.00           O
HETATM   10  O3  SO4 A 200      14.800  13.800  15.800  1.00 20.00           O
HETATM   11  O4  SO4 A 200      15.300  15.200  13.600  1.00 20.00           O
HETATM   12  O   HOH A 301      20.000  20.000  20.000  1.00 25.00           O
HETATM   13  O   WAT A 302      22.000  22.000  22.000  1.00 25.00           O
HETATM   14  S   SO4 A 400      25.000  25.000  25.000  1.00 20.00           S
END
"""


def exercise_build_kinemage():
  """Test that _build_kinemage produces valid kinemage output with expected
  sections."""
  from mmtbx.kinemage.validation import (
    build_name_hash, _build_bond_hash, _build_kinemage)
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import monomer_library
  from cctbx import geometry_restraints
  from mmtbx.validation.rotalyze import rotalyze
  from mmtbx.validation.ramalyze import ramalyze
  from mmtbx.validation.cbetadev import cbetadev
  from mmtbx.validation.restraints import combined as restraints_combined
  from iotbx import pdb

  pdb_io = pdb.input(source_info=None, lines=pdb_str)
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    pdb_inp=pdb_io,
    substitute_non_crystallographic_unit_cell_if_necessary=True)

  hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  sites_cart = processed_pdb_file.all_chain_proxies.sites_cart
  geometry = processed_pdb_file.geometry_restraints_manager()
  xray_structure = processed_pdb_file.xray_structure()
  i_seq_name_hash = build_name_hash(pdb_hierarchy=hierarchy)
  flags = geometry_restraints.flags.flags(default=True)
  pair_proxies = geometry.pair_proxies(flags=flags, sites_cart=sites_cart)
  bond_proxies = pair_proxies.bond_proxies
  quick_bond_hash = _build_bond_hash(bond_proxies, i_seq_name_hash)

  rot_outliers = rotalyze(pdb_hierarchy=hierarchy, outliers_only=True)
  cb = cbetadev(pdb_hierarchy=hierarchy, outliers_only=True)
  rama = ramalyze(pdb_hierarchy=hierarchy, outliers_only=True)
  restraints_result = restraints_combined(
    pdb_hierarchy=hierarchy,
    xray_structure=xray_structure,
    geometry_restraints_manager=geometry,
    ignore_hd=True,
    outliers_only=True)

  kin_out = _build_kinemage(
    hierarchy=hierarchy,
    bond_hash=quick_bond_hash,
    i_seq_name_hash=i_seq_name_hash,
    pdbID="test",
    rot_outliers=rot_outliers,
    rama_result=rama,
    cb_result=cb,
    restraints_result=restraints_result,
    keep_hydrogens=False)

  # Check that critical kinemage sections are present
  assert "@kinemage 1" in kin_out, "Missing @kinemage header"
  assert "@group {test}" in kin_out, "Missing @group"
  assert "@subgroup" in kin_out, "Missing @subgroup"
  assert "@vectorlist {mc}" in kin_out, "Missing mainchain vectorlist"
  assert "master= {mainchain}" in kin_out, "Missing mainchain master"
  assert "@vectorlist {Calphas}" in kin_out, "Missing Calphas vectorlist"
  # Check that output is a string (not accidentally iterating chars)
  assert isinstance(kin_out, str)
  print("  exercise_build_kinemage: OK")


def exercise_altloc_handling():
  """Test that the altloc assignment bug (== vs =) is fixed by verifying
  that altloc PDB files produce valid kinemage output without errors."""
  from mmtbx.kinemage.validation import (
    build_name_hash, _build_bond_hash, get_kin_lots)
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import monomer_library
  from cctbx import geometry_restraints
  from iotbx import pdb

  pdb_io = pdb.input(source_info=None, lines=pdb_altloc_str)
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    pdb_inp=pdb_io,
    substitute_non_crystallographic_unit_cell_if_necessary=True)

  hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  sites_cart = processed_pdb_file.all_chain_proxies.sites_cart
  geometry = processed_pdb_file.geometry_restraints_manager()
  i_seq_name_hash = build_name_hash(pdb_hierarchy=hierarchy)
  flags = geometry_restraints.flags.flags(default=True)
  pair_proxies = geometry.pair_proxies(flags=flags, sites_cart=sites_cart)
  bond_proxies = pair_proxies.bond_proxies
  quick_bond_hash = _build_bond_hash(bond_proxies, i_seq_name_hash)

  for model in hierarchy.models():
    for chain in model.chains():
      kin_out = get_kin_lots(
        chain=chain,
        bond_hash=quick_bond_hash,
        i_seq_name_hash=i_seq_name_hash,
        pdbID="altloc_test",
        index=0)
      # Should produce output without errors
      assert isinstance(kin_out, str)
      # With altlocs, we should see altloc controls if present
      assert len(kin_out) > 0, "Empty kinemage output for altloc PDB"
  print("  exercise_altloc_handling: OK")


def exercise_make_multikin():
  """Test the make_multikin entry point end-to-end."""
  from mmtbx.kinemage.validation import make_multikin
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import monomer_library
  from iotbx import pdb
  import tempfile

  pdb_io = pdb.input(source_info=None, lines=pdb_str)
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    pdb_inp=pdb_io,
    substitute_non_crystallographic_unit_cell_if_necessary=True)

  outfile = tempfile.mktemp(suffix='.kin')
  try:
    result = make_multikin(
      f=outfile,
      processed_pdb_file=processed_pdb_file,
      pdbID="test_kin",
      keep_hydrogens=True)
    assert os.path.exists(result), "Output file not created"
    with open(result, 'r') as f:
      content = f.read()
    assert "@kinemage 1" in content, "Missing kinemage header in output file"
    assert "@group {test_kin}" in content, "Missing group in output file"
    assert len(content) > 100, "Output file too small"
  finally:
    if os.path.exists(outfile):
      os.remove(outfile)
  print("  exercise_make_multikin: OK")


def exercise_helper_functions():
  """Test the extracted helper functions individually."""
  from mmtbx.kinemage.validation import (
    midpoint, get_chain_color, get_ions, get_default_header, get_footer,
    _get_prev_connection, _draw_rna_virtual_backbone)

  # midpoint
  mid = midpoint([0,0,0], [2,2,2])
  assert mid == [1.0, 1.0, 1.0], "midpoint failed"

  # get_chain_color
  assert get_chain_color(0) == 'white'
  assert get_chain_color(7) == 'white'  # wraps around

  # get_ions
  ions = get_ions("{test} 1.0 2.0 3.0\n")
  assert "@spherelist" in ions

  # header and footer
  assert "@kinemage 1" in get_default_header()
  assert "@master" in get_footer()

  # _get_prev_connection
  prev_key, prev_xyz = _get_prev_connection(
    {'A': 'keyA', ' ': 'key_blank'},
    {'A': (1,2,3), ' ': (4,5,6)},
    'A')
  assert prev_key == 'keyA'
  prev_key, prev_xyz = _get_prev_connection(
    {' ': 'key_blank'},
    {' ': (4,5,6)},
    'A')
  assert prev_key == 'key_blank'
  prev_key, prev_xyz = _get_prev_connection({}, {}, 'A')
  assert prev_key is None

  # _draw_rna_virtual_backbone with missing keys should return empty
  vbb = _draw_rna_virtual_backbone(
    type('MockRG', (), {'resseq_as_int': lambda self: 5})(),
    {}, {}, {}, {}, {}, {})
  assert vbb == ""

  print("  exercise_helper_functions: OK")


def exercise_deleted_functions():
  """Verify that deleted functions are truly gone from the module."""
  import mmtbx.kinemage.validation as v
  assert not hasattr(v, 'add_fan'), "add_fan should have been deleted"
  assert not hasattr(v, 'add_spring'), "add_spring should have been deleted"
  assert not hasattr(v, 'get_residue_bonds'), "get_residue_bonds should have been deleted"
  assert not hasattr(v, 'get_angle_outliers'), "get_angle_outliers should have been deleted"
  assert not hasattr(v, 'get_bond_outliers'), "get_bond_outliers should have been deleted"
  assert not hasattr(v, 'rotamer_outliers'), "rotamer_outliers should have been deleted"
  assert not hasattr(v, 'rama_outliers'), "rama_outliers should have been deleted"
  assert not hasattr(v, 'pperp_outliers'), "pperp_outliers should have been deleted"
  print("  exercise_deleted_functions: OK")


def exercise_het_ion_water():
  """Test handling of heteroatoms (ligands), ions, and waters."""
  from mmtbx.kinemage.validation import (
    build_name_hash, _build_bond_hash, get_kin_lots)
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import monomer_library
  from cctbx import geometry_restraints
  from iotbx import pdb

  pdb_io = pdb.input(source_info=None, lines=pdb_het_str)
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    pdb_inp=pdb_io,
    substitute_non_crystallographic_unit_cell_if_necessary=True)

  hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  sites_cart = processed_pdb_file.all_chain_proxies.sites_cart
  geometry = processed_pdb_file.geometry_restraints_manager()
  i_seq_name_hash = build_name_hash(pdb_hierarchy=hierarchy)
  flags = geometry_restraints.flags.flags(default=True)
  pair_proxies = geometry.pair_proxies(flags=flags, sites_cart=sites_cart)
  bond_proxies = pair_proxies.bond_proxies
  quick_bond_hash = _build_bond_hash(bond_proxies, i_seq_name_hash)

  for model in hierarchy.models():
    for chain in model.chains():
      kin_out = get_kin_lots(
        chain=chain,
        bond_hash=quick_bond_hash,
        i_seq_name_hash=i_seq_name_hash,
        pdbID="het_test",
        index=0)

      # Ion (ZN) should appear as a spherelist
      assert "@spherelist {het M}" in kin_out, "Missing ion spherelist"
      assert "color= gray" in kin_out, "Ion spherelist should be gray"
      # The ZN key should appear in the ion spherelist
      assert "zn" in kin_out.lower(), "ZN ion not found in output"

      # SO4 ligand should have het bonds drawn (pink vectorlist)
      assert "@vectorlist {het}" in kin_out, "Missing het vectorlist"
      assert "color= pink" in kin_out, "Het bonds should be pink"
      # SO4 has S-O bonds that should be drawn
      assert "so4" in kin_out.lower(), "SO4 ligand not found in output"

      # HOH water should appear as a balllist
      assert "@balllist {water O}" in kin_out, "Missing water balllist"
      assert "peachtint" in kin_out, "Water balllist should be peachtint"

      # WAT water should ALSO appear as a water ball (not as a het)
      # This tests the fix: using res_class == "common_water" instead of
      # resname.lower() == 'hoh'
      water_lines = [l for l in kin_out.splitlines()
                     if "wat" in l.lower() and "water" not in l.lower()
                     and "hoh" not in l.lower()]
      # WAT should NOT appear in the het vectorlist
      het_section = ""
      in_het = False
      for line in kin_out.splitlines():
        if "@vectorlist {het}" in line:
          in_het = True
          continue
        if in_het and line.startswith("@"):
          in_het = False
        if in_het:
          het_section += line + "\n"
      assert "wat" not in het_section.lower(), \
        "WAT should be drawn as water ball, not as het bonds"

      # Single-atom SO4 (res 400) should appear as a sphere in the ion list,
      # not be invisible. This tests the fix for structures like 1nxb where
      # SO4 has only an S atom.
      ion_section = ""
      in_ion = False
      for line in kin_out.splitlines():
        if "@spherelist {het M}" in line:
          in_ion = True
        if in_ion:
          ion_section += line + "\n"
      # The single-atom SO4 should be in the ion spherelist alongside ZN
      assert " s  " in ion_section.lower() or "so4" in ion_section.lower(), \
        "Single-atom SO4 should appear in ion spherelist"

  print("  exercise_het_ion_water: OK")


def run():
  print("Testing mmtbx.kinemage.validation:")
  exercise_helper_functions()
  exercise_deleted_functions()
  exercise_altloc_handling()
  exercise_build_kinemage()
  exercise_make_multikin()
  exercise_het_ion_water()
  print("All tests passed.")


if __name__ == "__main__":
  run()
