"""Tests for mmtbx.kinemage.validation"""

from __future__ import absolute_import, division, print_function
import os
import sys
from libtbx import easy_run

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

# Longer PDB with 8 residues and HELIX/SHEET records for ribbon testing.
# Two-stranded beta sheet (res 1-3, 6-8) and a short loop (res 4-5).
# Coordinates are idealized to give proper CA distances (~3.8A).
pdb_ribbon_str = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
HELIX    1 H1  ALA A    4  ALA A    5  1                                   2
SHEET    1 S1  2 ALA A   1  ALA A   3  0
SHEET    2 S1  2 ALA A   6  ALA A   8 -1
ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       2.430   1.000   1.000  1.00 10.00           C
ATOM      3  C   ALA A   1       3.013   2.404   1.000  1.00 10.00           C
ATOM      4  O   ALA A   1       2.253   3.370   1.000  1.00 10.00           O
ATOM      5  CB  ALA A   1       2.930   0.250   2.230  1.00 10.00           C
ATOM      6  N   ALA A   2       4.315   2.500   1.000  1.00 10.00           N
ATOM      7  CA  ALA A   2       4.990   3.800   1.000  1.00 10.00           C
ATOM      8  C   ALA A   2       6.500   3.660   1.000  1.00 10.00           C
ATOM      9  O   ALA A   2       7.050   2.560   1.000  1.00 10.00           O
ATOM     10  CB  ALA A   2       4.540   4.640   2.200  1.00 10.00           C
ATOM     11  N   ALA A   3       7.140   4.840   1.000  1.00 10.00           N
ATOM     12  CA  ALA A   3       8.600   4.900   1.000  1.00 10.00           C
ATOM     13  C   ALA A   3       9.110   6.330   1.000  1.00 10.00           C
ATOM     14  O   ALA A   3       8.330   7.290   1.000  1.00 10.00           O
ATOM     15  CB  ALA A   3       9.140   4.150   2.220  1.00 10.00           C
ATOM     16  N   ALA A   4      10.410   6.400   1.000  1.00 10.00           N
ATOM     17  CA  ALA A   4      11.050   7.700   1.000  1.00 10.00           C
ATOM     18  C   ALA A   4      12.560   7.600   1.000  1.00 10.00           C
ATOM     19  O   ALA A   4      13.110   6.500   1.000  1.00 10.00           O
ATOM     20  CB  ALA A   4      10.550   8.540   2.200  1.00 10.00           C
ATOM     21  N   ALA A   5      13.200   8.770   1.000  1.00 10.00           N
ATOM     22  CA  ALA A   5      14.650   8.830   1.000  1.00 10.00           C
ATOM     23  C   ALA A   5      15.160  10.260   1.000  1.00 10.00           C
ATOM     24  O   ALA A   5      14.380  11.220   1.000  1.00 10.00           O
ATOM     25  CB  ALA A   5      15.190   8.080   2.220  1.00 10.00           C
ATOM     26  N   ALA A   6      16.460  10.340   1.000  1.00 10.00           N
ATOM     27  CA  ALA A   6      17.100  11.640   1.000  1.00 10.00           C
ATOM     28  C   ALA A   6      18.610  11.540   1.000  1.00 10.00           C
ATOM     29  O   ALA A   6      19.160  10.440   1.000  1.00 10.00           O
ATOM     30  CB  ALA A   6      16.600  12.480   2.200  1.00 10.00           C
ATOM     31  N   ALA A   7      19.250  12.710   1.000  1.00 10.00           N
ATOM     32  CA  ALA A   7      20.700  12.770   1.000  1.00 10.00           C
ATOM     33  C   ALA A   7      21.210  14.200   1.000  1.00 10.00           C
ATOM     34  O   ALA A   7      20.430  15.160   1.000  1.00 10.00           O
ATOM     35  CB  ALA A   7      21.240  12.020   2.220  1.00 10.00           C
ATOM     36  N   ALA A   8      22.510  14.270   1.000  1.00 10.00           N
ATOM     37  CA  ALA A   8      23.150  15.570   1.000  1.00 10.00           C
ATOM     38  C   ALA A   8      24.660  15.470   1.000  1.00 10.00           C
ATOM     39  O   ALA A   8      25.210  14.370   1.000  1.00 10.00           O
ATOM     40  CB  ALA A   8      22.650  16.410   2.200  1.00 10.00           C
END
"""

# PDB with altlocs to test the altloc bug fix (line 650 in original)
# PDB with two CYS residues forming a disulfide bond (SG-SG ~2.03A)
pdb_disulfide_str = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
SSBOND   1 CYS A    2    CYS A    4                          2.03
ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       2.458   1.000   1.000  1.00 10.00           C
ATOM      3  C   ALA A   1       3.009   2.425   1.000  1.00 10.00           C
ATOM      4  O   ALA A   1       2.249   3.390   1.000  1.00 10.00           O
ATOM      5  CB  ALA A   1       2.982   0.231   2.207  1.00 10.00           C
ATOM      6  N   CYS A   2       4.316   2.498   1.000  1.00 10.00           N
ATOM      7  CA  CYS A   2       4.988   3.798   1.000  1.00 10.00           C
ATOM      8  C   CYS A   2       6.504   3.664   1.000  1.00 10.00           C
ATOM      9  O   CYS A   2       7.054   2.562   1.000  1.00 10.00           O
ATOM     10  CB  CYS A   2       4.538   4.637   2.196  1.00 10.00           C
ATOM     11  SG  CYS A   2       5.100   6.350   2.100  1.00 10.00           S
ATOM     12  N   ALA A   3       7.142   4.835   1.000  1.00 10.00           N
ATOM     13  CA  ALA A   3       8.598   4.893   1.000  1.00 10.00           C
ATOM     14  C   ALA A   3       9.108   6.330   1.000  1.00 10.00           C
ATOM     15  O   ALA A   3       8.330   7.286   1.000  1.00 10.00           O
ATOM     16  CB  ALA A   3       9.137   4.148   2.215  1.00 10.00           C
ATOM     17  N   CYS A   4      10.410   6.400   1.000  1.00 10.00           N
ATOM     18  CA  CYS A   4      11.050   7.700   1.000  1.00 10.00           C
ATOM     19  C   CYS A   4      12.560   7.600   1.000  1.00 10.00           C
ATOM     20  O   CYS A   4      13.110   6.500   1.000  1.00 10.00           O
ATOM     21  CB  CYS A   4      10.550   8.540   2.200  1.00 10.00           C
ATOM     22  SG  CYS A   4       5.120   8.380   2.100  1.00 10.00           S
END
"""

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
  from mmtbx.kinemage import kin_vec
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


def exercise_ribbon_rendering():
  """Test the ribbon_rendering module classes and functions."""
  from mmtbx.kinemage.ribbon_rendering import (
      Crayon, Range, RibbonElement,
      build_secondary_structure_map, consolidate_sheets,
      spline_interpolate, generate_chain_ribbons)
  from scitbx.matrix import col
  from iotbx import pdb
  from libtbx.utils import null_out

  # --- Test Range class ---
  r = Range()
  assert r.type == 'COIL'
  assert r.strand == 1
  assert r.previous is None

  r = Range(ss_type='HELIX', chain_id='A', init_seq=10, end_seq=20)
  assert r.type == 'HELIX'
  assert r.initSeqNum == 10
  assert r.endSeqNum == 20

  # --- Test RibbonElement class ---
  re1 = RibbonElement()
  assert re1.type == 'COIL'
  re2 = RibbonElement(set_range=Range(ss_type='HELIX'))
  assert re2.type == 'HELIX'
  re3 = RibbonElement(set_range=Range(ss_type='TURN'))
  assert re3.type == 'COIL'  # TURN maps to COIL
  # same_sse
  assert not re1.same_sse(re2)
  r_shared = Range(ss_type='SHEET')
  re4 = RibbonElement(set_range=r_shared)
  re5 = RibbonElement(set_range=r_shared)
  assert re4.same_sse(re5)
  # sorting by strand
  re4.range.strand = 2
  re5.range = Range(ss_type='SHEET', strand=1)
  assert not re4 < re5  # strand 2 > strand 1

  # --- Test Crayon class ---
  c = Crayon(do_rainbow=False, color="red")
  assert c.get_kin_string() == "red"
  assert c.should_print()
  c2 = Crayon(do_rainbow=False)
  assert c2.get_kin_string() == ""

  # --- Test spline_interpolate ---
  # Straight line: 6 points along x-axis (need >= 4 for one spline segment)
  pts = [col((float(i), 0.0, 0.0)) for i in range(6)]
  result = spline_interpolate(pts, 4)
  assert len(result) > 0
  # All y and z coords should be ~0 for a straight line
  for pt in result:
    assert abs(pt[1]) < 1e-6
    assert abs(pt[2]) < 1e-6

  # --- Test build_secondary_structure_map with annotation ---
  hierarchy = pdb.input(source_info=None, lines=pdb_ribbon_str).construct_hierarchy()
  annotation = pdb.input(source_info=None,
      lines=pdb_ribbon_str).extract_secondary_structure()
  ss_map = build_secondary_structure_map(
      hierarchy, annotation=annotation, log=null_out())
  # Check that the map has entries for all 8 residues
  assert len(ss_map) == 8
  # Residues 1-3 should be SHEET, 4-5 should be HELIX, 6-8 should be SHEET
  assert ss_map[1].type == 'SHEET', "res 1 should be SHEET, got %s" % ss_map[1].type
  assert ss_map[3].type == 'SHEET', "res 3 should be SHEET, got %s" % ss_map[3].type
  assert ss_map[4].type == 'HELIX', "res 4 should be HELIX, got %s" % ss_map[4].type
  assert ss_map[5].type == 'HELIX', "res 5 should be HELIX, got %s" % ss_map[5].type
  assert ss_map[6].type == 'SHEET', "res 6 should be SHEET, got %s" % ss_map[6].type
  assert ss_map[8].type == 'SHEET', "res 8 should be SHEET, got %s" % ss_map[8].type

  # strand IDs should be ints after conversion
  for rng in ss_map.values():
    if rng.type == 'SHEET':
      assert isinstance(rng.strand, int), "strand should be int, got %s" % type(rng.strand)

  # --- Test consolidate_sheets ---
  consolidate_sheets(ss_map)
  # Sheet strands with strand==2 should have previous pointing to strand==1
  for rng in ss_map.values():
    if rng.type == 'SHEET' and rng.strand == 2 and rng.duplicateOf is None:
      assert rng.previous is not None, "strand 2 should have a previous strand"

  # --- Test build_secondary_structure_map without annotation (from_ca fallback) ---
  ss_map_auto = build_secondary_structure_map(
      hierarchy, annotation=None, log=null_out())
  assert len(ss_map_auto) == 8
  # All entries should be valid Range objects
  for rng in ss_map_auto.values():
    assert rng.type in ('COIL', 'HELIX', 'SHEET')

  # --- Test generate_chain_ribbons ---
  for chain in hierarchy.models()[0].chains():
    ribbon_out = generate_chain_ribbons(
        chain=chain, secondary_structure=ss_map,
        chain_id=chain.id, chain_color="white")
    assert len(ribbon_out) > 0, "generate_chain_ribbons produced no output"
    assert "@subgroup" in ribbon_out
    assert "@colorset" in ribbon_out
    # Should contain ribbon list content
    assert "master= {ribbon}" in ribbon_out

  print("  exercise_ribbon_rendering: OK")


def exercise_ribbon_in_kinemage():
  """Test that ribbon content appears within the main group in
  _build_kinemage output, not as a separate top-level group."""
  from mmtbx.kinemage.validation import (
      build_name_hash, _build_bond_hash, _build_kinemage,
      get_default_header, get_footer)
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import monomer_library
  from cctbx import geometry_restraints
  from mmtbx.validation.rotalyze import rotalyze
  from mmtbx.validation.ramalyze import ramalyze
  from mmtbx.validation.cbetadev import cbetadev
  from mmtbx.validation.restraints import combined as restraints_combined
  from iotbx import pdb
  from libtbx.utils import null_out

  # Check that header and footer have ribbon master controls
  header = get_default_header()
  footer = get_footer()
  assert "@master {ribbon} indent" in header, "Missing ribbon master in header"
  assert "@master {ribbon} off" in footer, "Missing ribbon master in footer"

  # Build full kinemage from the 8-residue PDB with SS records
  pdb_io = pdb.input(source_info=None, lines=pdb_ribbon_str)
  annotation = pdb_io.extract_secondary_structure()
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
      pdbID="ribbon_test",
      rot_outliers=rot_outliers,
      rama_result=rama,
      cb_result=cb,
      restraints_result=restraints_result,
      keep_hydrogens=False,
      ss_annotation=annotation)

  # Ribbon content should be present
  assert "master= {ribbon}" in kin_out, "No ribbon content in kinemage output"
  assert "@subgroup {ribbon" in kin_out, "Missing ribbon subgroup"
  assert "@colorset" in kin_out, "Missing ribbon colorset"

  # Ribbon should NOT be a separate top-level @group
  assert "@group {ribbon}" not in kin_out, \
      "Ribbon should be a subgroup within the main group, not a separate @group"

  # There should be exactly one top-level @group (the main model group)
  group_lines = [l for l in kin_out.splitlines()
                 if l.startswith("@group")]
  assert len(group_lines) == 1, \
      "Expected 1 @group, got %d: %s" % (len(group_lines), group_lines)
  assert "{ribbon_test}" in group_lines[0]

  print("  exercise_ribbon_in_kinemage: OK")


def exercise_deleted_functions():
  """Verify that deleted functions are truly gone from the module."""
  import mmtbx.kinemage.validation as v
  assert not hasattr(v, 'add_fan'), "add_fan should have been deleted"
  assert not hasattr(v, 'add_spring'), "add_spring should have been deleted"
  assert not hasattr(v, 'get_residue_bonds'), "get_residue_bonds should have been deleted"
  assert not hasattr(v, 'get_angle_outliers'), "get_angle_outliers should have been deleted"
  assert not hasattr(v, 'get_bond_outliers'), "get_bond_outliers should have been deleted"
  print("  exercise_deleted_functions: OK")


def exercise_build_name_hash():
  """Test that build_name_hash returns structured dicts (not raw strings)."""
  from mmtbx.kinemage.validation import build_name_hash
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import monomer_library
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
  name_hash = build_name_hash(pdb_hierarchy=hierarchy)

  # Should have entries for all atoms
  assert len(name_hash) == len(list(hierarchy.atoms()))

  # Each entry should be a dict with the expected keys
  for i_seq, info in name_hash.items():
    assert isinstance(info, dict), "Expected dict, got %s" % type(info)
    for key in ('name', 'altloc', 'resname', 'chain_id', 'resseq', 'icode'):
      assert key in info, "Missing key '%s' in name_hash entry" % key

  # Check specific values for the first atom (N of ALA A 1)
  first_atom = list(hierarchy.atoms())[0]
  info = name_hash[first_atom.i_seq]
  assert info['name'].strip() == 'N', "Expected 'N', got '%s'" % info['name'].strip()
  assert info['resname'].strip() == 'ALA', "Expected 'ALA', got '%s'" % info['resname'].strip()
  assert info['chain_id'].strip() == 'A', "Expected 'A', got '%s'" % info['chain_id'].strip()

  print("  exercise_build_name_hash: OK")


def exercise_same_residue():
  """Test the _same_residue helper function."""
  from mmtbx.kinemage.validation import _same_residue

  info_a = {'chain_id': 'A', 'resseq': '   1', 'icode': ' ',
            'name': ' N  ', 'altloc': ' ', 'resname': 'ALA'}
  info_b = {'chain_id': 'A', 'resseq': '   1', 'icode': ' ',
            'name': ' CA ', 'altloc': ' ', 'resname': 'ALA'}
  info_c = {'chain_id': 'A', 'resseq': '   2', 'icode': ' ',
            'name': ' N  ', 'altloc': ' ', 'resname': 'ALA'}
  info_d = {'chain_id': 'B', 'resseq': '   1', 'icode': ' ',
            'name': ' N  ', 'altloc': ' ', 'resname': 'ALA'}
  info_e = {'chain_id': 'A', 'resseq': '   1', 'icode': 'A',
            'name': ' N  ', 'altloc': ' ', 'resname': 'ALA'}

  # Same residue (different atom names)
  assert _same_residue(info_a, info_b) == True
  # Different resseq
  assert _same_residue(info_a, info_c) == False
  # Different chain
  assert _same_residue(info_a, info_d) == False
  # Different icode
  assert _same_residue(info_a, info_e) == False

  print("  exercise_same_residue: OK")


def exercise_disulfide_bonds():
  """Test _build_ss_bond_list and SS bond drawing in get_kin_lots."""
  from mmtbx.kinemage.validation import (
    build_name_hash, _build_bond_hash, _build_ss_bond_list, get_kin_lots)
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import monomer_library
  from cctbx import geometry_restraints
  from iotbx import pdb

  pdb_io = pdb.input(source_info=None, lines=pdb_disulfide_str)
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

  # Test _build_ss_bond_list
  ss_bonds = _build_ss_bond_list(bond_proxies, i_seq_name_hash)
  assert len(ss_bonds) >= 1, "Expected at least 1 SS bond, got %d" % len(ss_bonds)
  # Each entry should be a tuple of two i_seqs
  for i_seq_0, i_seq_1 in ss_bonds:
    info_0 = i_seq_name_hash[i_seq_0]
    info_1 = i_seq_name_hash[i_seq_1]
    assert info_0['name'].strip() == 'SG', "Expected SG, got %s" % info_0['name'].strip()
    assert info_1['name'].strip() == 'SG', "Expected SG, got %s" % info_1['name'].strip()
    assert info_0['resname'].strip() == 'CYS'
    assert info_1['resname'].strip() == 'CYS'

  # Test that get_kin_lots includes SS bond vectorlist when ss_bonds provided
  quick_bond_hash = _build_bond_hash(bond_proxies, i_seq_name_hash)
  for model in hierarchy.models():
    for chain in model.chains():
      kin_out = get_kin_lots(
        chain=chain,
        bond_hash=quick_bond_hash,
        i_seq_name_hash=i_seq_name_hash,
        pdbID="ss_test",
        index=0,
        ss_bonds=ss_bonds,
        sites_cart=sites_cart)
      assert "@vectorlist {SS}" in kin_out, "Missing SS vectorlist"
      assert "color= yellow" in kin_out, "SS bonds should be yellow"

  # Test that get_kin_lots without ss_bonds does NOT include SS vectorlist
  for model in hierarchy.models():
    for chain in model.chains():
      kin_out_no_ss = get_kin_lots(
        chain=chain,
        bond_hash=quick_bond_hash,
        i_seq_name_hash=i_seq_name_hash,
        pdbID="ss_test",
        index=0)
      assert "@vectorlist {SS}" not in kin_out_no_ss, \
        "SS vectorlist should not appear without ss_bonds"

  print("  exercise_disulfide_bonds: OK")


def exercise_draw_residue_bonds():
  """Test _draw_residue_bonds sorting into mc/sc/het categories."""
  from mmtbx.kinemage.validation import (
    build_name_hash, _build_bond_hash, _draw_residue_bonds)
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import monomer_library
  from cctbx import geometry_restraints
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
  i_seq_name_hash = build_name_hash(pdb_hierarchy=hierarchy)
  flags = geometry_restraints.flags.flags(default=True)
  pair_proxies = geometry.pair_proxies(flags=flags, sites_cart=sites_cart)
  bond_proxies = pair_proxies.bond_proxies
  quick_bond_hash = _build_bond_hash(bond_proxies, i_seq_name_hash)

  mc_atoms = ["N", "CA", "C", "O", "OXT",
              "P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'", "O4'", "C1'",
              "C3'", "O3'", "C2'", "O2'"]

  # Test with the first residue
  for model in hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        for conformer in residue_group.conformers():
          for residue in conformer.residues():
            key_hash = {}
            xyz_hash = {}
            het_hash = {}
            iseq_altloc = {}
            for atom in residue.atoms():
              key = "%s %s %s" % (atom.name.lower(), residue.resname.lower(),
                                  chain.id)
              key_hash[atom.name.strip()] = key
              xyz_hash[atom.name.strip()] = atom.xyz
              iseq_altloc[atom.i_seq] = ' '

            drawn_bonds = []
            result = _draw_residue_bonds(
              residue, quick_bond_hash, i_seq_name_hash, key_hash,
              xyz_hash, het_hash, iseq_altloc, ' ',
              mc_atoms, False, drawn_bonds)

            # Result should be a dict with expected keys
            assert isinstance(result, dict)
            for key in ('mc', 'sc', 'mc_h', 'sc_h', 'het', 'het_h'):
              assert key in result, "Missing key '%s' in result" % key

            # For ALA, we should have mainchain bonds (N-CA, CA-C, C-O)
            # and sidechain bonds (CA-CB)
            assert len(result['mc']) > 0, "Expected mainchain bonds for ALA"
            assert len(result['sc']) > 0, "Expected sidechain bonds for ALA (CA-CB)"
            # No hydrogens in input, so no H bonds
            assert result['mc_h'] == ''
            assert result['sc_h'] == ''
            # ALA is not a het
            assert result['het'] == ''
            assert result['het_h'] == ''

            # Only test first residue
            print("  exercise_draw_residue_bonds: OK")
            return

  print("  exercise_draw_residue_bonds: OK")


def exercise_track_amino_acid_atom():
  """Test _track_amino_acid_atom backbone tracking for inter-residue bonds."""
  from mmtbx.kinemage.validation import _track_amino_acid_atom
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import monomer_library
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
  chain = list(hierarchy.models())[0].chains()[0]
  residue_groups = list(chain.residue_groups())
  assert len(residue_groups) == 3

  # Simulate tracking through the 3 residues
  prev_C_key = {}
  prev_C_xyz = {}
  prev_CA_key = {}
  prev_CA_xyz = {}
  prev_resid = None

  for rg in residue_groups:
    cur_C_xyz = {}
    cur_C_key = {}
    cur_CA_xyz = {}
    cur_CA_key = {}
    mc_parts = []
    ca_parts = []
    for ag in rg.atom_groups():
      for atom in ag.atoms():
        key = "key_%s_%s" % (atom.name.strip(), rg.resseq_as_int())
        _track_amino_acid_atom(
          atom, key, ' ', rg, prev_resid,
          cur_C_xyz, cur_C_key, cur_CA_xyz, cur_CA_key,
          prev_C_key, prev_C_xyz, prev_CA_key, prev_CA_xyz,
          mc_parts, ca_parts)

    if rg.resseq_as_int() > 1:
      # After processing residue 2 or 3, there should be inter-residue
      # connections if the previous residue had C/CA atoms
      assert len(mc_parts) > 0 or len(ca_parts) > 0, \
        "Expected inter-residue connections for residue %d" % rg.resseq_as_int()

    prev_C_key = cur_C_key
    prev_C_xyz = cur_C_xyz
    prev_CA_key = cur_CA_key
    prev_CA_xyz = cur_CA_xyz
    prev_resid = rg.resid()

  print("  exercise_track_amino_acid_atom: OK")


def exercise_track_rna_dna_atom():
  """Test _track_rna_dna_atom with mock RNA-like atoms."""
  from mmtbx.kinemage.validation import _track_rna_dna_atom

  # Create mock atom and residue_group objects
  class MockAtom:
    def __init__(self, name, xyz, i_seq=0):
      self.name = name
      self.xyz = xyz
      self.i_seq = i_seq

  class MockResidueGroup:
    def __init__(self, resseq):
      self._resseq = resseq
    def resseq_as_int(self):
      return self._resseq
    def resid(self):
      return "%4d " % self._resseq

  # Track a P atom for residue 2
  cur_O3_xyz = {}
  cur_O3_key = {}
  prev_O3_key = {' ': 'prev_O3_key'}
  prev_O3_xyz = {' ': (1.0, 2.0, 3.0)}
  p_hash_key = {}
  p_hash_xyz = {}
  c1_hash_key = {}
  c1_hash_xyz = {}
  c4_hash_key = {}
  c4_hash_xyz = {}
  mc_parts = []

  rg = MockResidueGroup(2)
  p_atom = MockAtom(' P  ', (4.0, 5.0, 6.0))

  _track_rna_dna_atom(
    p_atom, 'p_key_2', ' ', rg, '   1 ',
    cur_O3_xyz, cur_O3_key,
    prev_O3_key, prev_O3_xyz,
    p_hash_key, p_hash_xyz,
    c1_hash_key, c1_hash_xyz,
    c4_hash_key, c4_hash_xyz,
    mc_parts)

  # P atom should be stored in p_hash
  assert 2 in p_hash_key, "P atom not stored in p_hash_key"
  assert p_hash_key[2] == 'p_key_2'
  assert p_hash_xyz[2] == (4.0, 5.0, 6.0)
  # Should have generated an O3'-P connection
  assert len(mc_parts) == 1, "Expected 1 O3'-P connection, got %d" % len(mc_parts)

  # Track C1' and C4' atoms
  c1_atom = MockAtom(" C1'", (7.0, 8.0, 9.0))
  c4_atom = MockAtom(" C4'", (10.0, 11.0, 12.0))

  _track_rna_dna_atom(
    c1_atom, 'c1_key_2', ' ', rg, '   1 ',
    cur_O3_xyz, cur_O3_key,
    prev_O3_key, prev_O3_xyz,
    p_hash_key, p_hash_xyz,
    c1_hash_key, c1_hash_xyz,
    c4_hash_key, c4_hash_xyz,
    mc_parts)
  _track_rna_dna_atom(
    c4_atom, 'c4_key_2', ' ', rg, '   1 ',
    cur_O3_xyz, cur_O3_key,
    prev_O3_key, prev_O3_xyz,
    p_hash_key, p_hash_xyz,
    c1_hash_key, c1_hash_xyz,
    c4_hash_key, c4_hash_xyz,
    mc_parts)

  assert 2 in c1_hash_key
  assert 2 in c4_hash_key

  # Track O3' atom
  o3_atom = MockAtom(" O3'", (13.0, 14.0, 15.0))
  _track_rna_dna_atom(
    o3_atom, 'o3_key_2', ' ', rg, '   1 ',
    cur_O3_xyz, cur_O3_key,
    prev_O3_key, prev_O3_xyz,
    p_hash_key, p_hash_xyz,
    c1_hash_key, c1_hash_xyz,
    c4_hash_key, c4_hash_xyz,
    mc_parts)

  assert ' ' in cur_O3_key
  assert cur_O3_key[' '] == 'o3_key_2'

  print("  exercise_track_rna_dna_atom: OK")


def exercise_make_multikin_with_ribbons():
  """Test that make_multikin output includes ribbon content."""
  from mmtbx.kinemage.validation import make_multikin
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import monomer_library
  from iotbx import pdb
  import tempfile

  # Use the 8-residue PDB with SS records for reliable ribbon generation
  pdb_io = pdb.input(source_info=None, lines=pdb_ribbon_str)
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
      pdbID="ribbon_mk",
      keep_hydrogens=True)
    assert os.path.exists(result), "Output file not created"
    with open(result, 'r') as f:
      content = f.read()

    # Ribbon content should be present in the output
    assert "master= {ribbon}" in content, "No ribbon content in make_multikin output"
    assert "@subgroup {ribbon" in content, "Missing ribbon subgroup"

    # Should also have standard kinemage elements
    assert "@kinemage 1" in content
    assert "@group {ribbon_mk}" in content

    # Ribbon master controls
    assert "@master {ribbon} indent" in content, "Missing ribbon master in header"
    assert "@master {ribbon} off" in content, "Missing ribbon master off in footer"
  finally:
    if os.path.exists(outfile):
      os.remove(outfile)
  print("  exercise_make_multikin_with_ribbons: OK")


def exercise_make_multikin_with_disulfide():
  """Test that make_multikin includes disulfide bonds."""
  from mmtbx.kinemage.validation import make_multikin
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import monomer_library
  from iotbx import pdb
  import tempfile

  pdb_io = pdb.input(source_info=None, lines=pdb_disulfide_str)
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
      pdbID="ss_mk",
      keep_hydrogens=True)
    assert os.path.exists(result), "Output file not created"
    with open(result, 'r') as f:
      content = f.read()

    # Disulfide bonds should be drawn
    assert "@vectorlist {SS}" in content, "Missing SS vectorlist in output"
    assert "color= yellow" in content, "SS bonds should be yellow"
  finally:
    if os.path.exists(outfile):
      os.remove(outfile)
  print("  exercise_make_multikin_with_disulfide: OK")


def run():
  print("Testing mmtbx.kinemage.validation:")
  exercise_helper_functions()
  exercise_deleted_functions()
  exercise_build_name_hash()
  exercise_same_residue()
  exercise_altloc_handling()
  exercise_build_kinemage()
  exercise_draw_residue_bonds()
  exercise_track_amino_acid_atom()
  exercise_track_rna_dna_atom()
  exercise_disulfide_bonds()
  exercise_make_multikin()
  exercise_make_multikin_with_ribbons()
  exercise_make_multikin_with_disulfide()
  exercise_ribbon_rendering()
  exercise_ribbon_in_kinemage()
  print("All tests passed.")


if __name__ == "__main__":
  run()
