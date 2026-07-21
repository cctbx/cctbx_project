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

## 2026-03-04
Speed up probe2 ~23% by replacing `rvec3` matrix operations with direct tuple arithmetic in hot loops (distance checks, dot/spike location), and caching DotInfo sort keys to avoid repeated `_makeName` string formatting during sorting.

## 2026-03-05
Fix missing probe dots in kinemage output. Two bugs in `make_probe_dots`: (1) passed top-level probe2 params to Optimizer instead of `params.probe` sub-params, (2) used nonexistent `DataManager.set_default_output_dir()` instead of `CCTBXParser`. Both errors were silently swallowed.

## 2026-03-10
Add `make_probe_dots_from_model()` to reuse an already-hydrogenated model for probe dots (avoids re-running reduce2). Add `probe_dots_kin` param to `_build_kinemage()` to accept pre-computed dots. Save hydrogenated model on `clashscore2` for downstream kinemage use.

## 2026-03-21
Tune residue score: use non-linear power curve (exponent 1.3) for clash severity, cap clash-count bonus at 4.0 via log2, raise twisted-peptide severity from 10 to 15. Add `tst_utils.py` with tests for severity functions and ranking invariants.

## 2026-03-28
Pass `useNeutronDistances=nuclear` to the reduce2 Optimizer call in `clashscore2.py`, so neutron-distance models use correct hydrogen placement during clash analysis.

## 2026-04-16
Add `build_kinemage_from_model()` high-level helper in `mmtbx/kinemage/validation.py` that takes an `mmtbx.model.manager` and handles the bond_hash / ss_bonds / validator / ss-annotation / probe-dots plumbing internally. Callers may inject pre-run validator results or probe dots. Add two regression tests covering the H and non-H paths.

## 2026-04-16
Add `include_cablam_wheels` and `plain_coils` toggle parameters to kinemage output. Cablam wheels (purple/magenta peptide-plane score wedges) are now off by default in `_build_kinemage`; gate the wheel rendering in `cablamalyze.as_kinemage(include_wheels=...)`. `plain_coils` suppresses the rear deadblack coil halo. Both flags threaded through all entry points (`build_kinemage_from_model`, `make_multikin`, `export_molprobity_result_as_kinemage`). Add toggle tests.

## 2026-04-24
Pass `restraint_objects` from parent model manager to sub-models in `make_probe_dots_from_model` so probe2 gets correct bond topology for CCD-restrained residues (e.g. AS).

## 2026-04-24
Add suitealyze RNA suite outlier markup to kinemage output. Thread `suite_result` through `_build_kinemage`, `build_kinemage_from_model`, `make_multikin`, and `export_molprobity_result_as_kinemage`.

## 2026-05-11
Fix ribbon rendering: nucleic acid ribbons now render as SHEET (flat with arrowheads, matching Prekin's BETA behavior) instead of thin coil lines. Raise minimum contiguous segment to 3 residues (matching Prekin). Filter degenerate short HELIX/SHEET ribbon elements to COIL.

## 2026-07-21
Grade RNA sugar-pucker and suite problems by severity instead of scoring every outlier the same. Add `_rna_suite_severity` and `_rna_pucker_severity` in `mmtbx/validation/utils.py`: a delta outlier (wrong sugar pucker) weighs 5.0, an epsilon-only problem is a milder backbone issue, and suite outliers scale with how far the conformer sits from any named cluster.

## 2026-07-21
Extend sugar-pucker validation to DNA. Lift the `is_rna` gate in `rna_validate.rna_puckers` (that check was only "has an O2'", so it excluded all DNA), add DNA-specific delta and O3'-perp thresholds to `rna_sugar_pucker_analysis` selected by a new `is_dna` argument, and admit DNA-only structures at the top-level gate in `validation/molprobity/__init__.py`. RNA behaviour and the pdb_interpretation restraints path are unchanged. Retuned thresholds alone were not enough: DNA needs a contradiction rule plus a P-perp confidence margin, because roughly a third of DNA sugars are conformationally intermediate where a binary call carries no information. Flags 0.44% of DNA against the old rule's 24.5%, with 2.5x enrichment on PDB-REDO movement where the old rule scores 0.92x, i.e. chance.

## 2026-07-21
Add test coverage for the RNA severity functions and the DNA pucker path. New `tst_rna_sugar_pucker_analysis_dna.py` covers the `is_dna` argument at the evaluate() level: that `is_dna=False` is the default and leaves the RNA constants untouched (the pdb_interpretation restraints path depends on this), that a delta contradicting P-perp is flagged, that a delta merely *between* the two windows is not (the case the DNA rule exists for, which the RNA thresholds get wrong), that the confidence margin suppresses a contradiction when P-perp is near its own cutoff, and that clearing a `dna_*` threshold falls back to the RNA value rather than disabling the test. New `tst_rna_validate_dna.py` covers the lifted gate end to end, `include_dna=False` restoring the old behaviour exactly, and the guard that skips a DNA residue with an incomplete sugar instead of raising KeyError. Three exercises added to `tst_utils.py` for `_rna_suite_severity` and `_rna_pucker_severity`, including an ordering check so the tiers cannot be retuned out of order. Both new files registered in `mmtbx/run_tests.py`.

