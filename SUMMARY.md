# Change Summary

## 2026-02-22 — 4449349c
Replace `pdb_label_columns()` column slicing with structured dict in `build_name_hash()` for CIF compatibility. Fix `pperp_outliers()` key format bug. Deduplicate bond hash logic in `cablam.py`.

## 2026-02-22 — 62c754a4
Add disulfide bond (SS) drawing to kinemage output. Detects CYS SG-SG pairs from geometry restraints and draws yellow vectors in the sidechain master group, matching old MolProbity/prekin format.

## 2026-02-22 — c7c67d52
Fix Rama outlier kinemage markup to draw from CA midpoints instead of actual CA coordinates, matching old MolProbity/prekin style.

## 2026-02-22 — 744c6fce
Refactor `mmtbx/kinemage/validation.py`: fix altloc bug, string comparison bugs, duplicate geometry deviations bug, duplicate footer line; delete dead code (add_fan, add_spring, get_residue_bonds, get_angle_outliers, get_bond_outliers); consolidate make_multikin/export into shared _build_kinemage; update make_probe_dots to probe2 API; extract helpers from get_kin_lots; add tst_validation.py.

## 2026-02-23 — dc7bbdfb
Add ribbon rendering to MolProbity kinemage output. Extract ribbon code from `programs/ribbons.py` into new `kinemage/ribbon_rendering.py` library; integrate into `_build_kinemage()` as off-by-default subgroup; add 15 tests covering ribbons, disulfide bonds, build_name_hash dict format, _same_residue, _draw_residue_bonds, backbone tracking, and make_multikin integration.

## 2026-02-23 — 177d54c9
Fix water and ligand classification in kinemage output. Use `res_class == "common_water"` instead of `resname == 'hoh'` so WAT/DOD waters render correctly. Make het bond drawing the fallback in `_draw_residue_bonds` so `common_saccharide` (NAG, etc.) gets bonds drawn.

## 2026-02-23 — e1d0e07b
Draw single-atom het residues (e.g., lone S from SO4 in 1nxb) as gray spheres instead of silently dropping them. Moved single-atom check after the water check to avoid misclassifying HOH.

## 2026-02-23 — 65a6cb7c
Fix errant ribbon triangle when a chain starts with a SHEET residue (visible in 1ubq). Rename ribbon master controls from `{protein}`/`{nucleic acid}` to `{protein ribbon}`/`{NA ribbon}` to avoid confusion with molecular representation masters.

## 2026-02-26 — 05e7174f
Fix reduce2 to emit PDB v3 atom names for RNA/DNA hydrogens. The monomer library returns old-convention names (H4*, H5*1, HO2*) but PDB v3 uses primes (H4', H5', HO2'). Added translation in `reduce_hydrogen.py` using the existing mapping from `iotbx.pdb`, with a fallback `*`→`'` for modified nucleotides. Updated tests and added `tst_reduce2_atom_names.py`.

## 2026-02-27 — 08b6f9a7
Redesign `calculate_overall_residue_quality_score` as a continuous triage priority metric. Replace additive integer penalties with severity-based combination rule (worst + 0.25 * rest). Use continuous values for clashes, C-beta, and bonds/angles. Scores now have 1 decimal place.

## 2026-02-27 — ef2b6618
Differentiate chiral outlier types in residue quality score. Split single chirality penalty into three severity levels: handedness swap (10.0), tetrahedral geometry outlier (5.0), pseudochiral naming error (1.5). Uses outlier_type() classification from restraints.py.

## 2026-03-01
Fix probe2/reduce2 crash on structures with atoms near crystallographic symmetry mates. `getBondedNeighborLists` now processes ASU bond proxies (j_sym==0), which contain covalent bonds the restraints engine expresses as ASU proxies instead of simple proxies when an atom is near a symmetry element. Fixes "Found Hydrogen with no neighbors" error (e.g. PDB 3q9v GLN A 179 HE21).

## 2026-03-04
Speed up clashscore2 ~2x by making the second probe2 run (full VDW output for Coot) optional via `save_probe_output` flag (default False). Also use set for clash dedup and iterate by residue group for water classification.
