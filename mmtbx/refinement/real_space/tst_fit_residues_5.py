from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space.fit_residues
import mmtbx.refinement.real_space
from scitbx.array_family import flex
from mmtbx.rotamer.rotamer_eval import RotamerEval
rotamer_eval = RotamerEval()

pdb_answer = """\
CRYST1   16.402   13.234   20.833  90.00  90.00  90.00 P 1
ATOM      1  N   LYS A 407       9.045   6.777  10.558  1.00 10.36           N
ATOM      2  CA  LYS A 407       7.962   7.610  10.037  1.00 10.78           C
ATOM      3  C   LYS A 407       7.986   8.907  10.843  1.00  8.79           C
ATOM      4  O   LYS A 407       7.084   9.187  11.629  1.00 10.53           O
ATOM      5  CB  LYS A 407       6.606   6.911  10.148  1.00 30.00           C
ATOM      6  CG  LYS A 407       6.444   5.717   9.220  1.00 30.00           C
ATOM      7  CD  LYS A 407       7.791   5.206   8.737  1.00 30.00           C
ATOM      8  CE  LYS A 407       8.932   6.044   9.288  1.00 30.00           C
ATOM      9  NZ  LYS A 407       8.484   6.932  10.396  1.00 30.00           N
END
"""

pdb_poor = pdb_answer

pdb_for_map = pdb_answer

def exercise(use_map, i_pdb, d_min = 1.5, resolution_factor = 0.1):
  """
  Best fitting residue is a rotamer outlier and a clash. Don't use
  vdW clashes.
  Use map: residue does not move (no better fit to map).
  Do not use map: residue switches to nearest non-outlier.
  """
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor,
    i_pdb             = i_pdb,
    d_min             = d_min,
    residues          = ["LYS"],
    resolution_factor = resolution_factor,
    pdb_for_map       = pdb_for_map)
  #
  found = False
  for model in t.ph_answer.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        conformers = residue_group.conformers()
        for conformer in residue_group.conformers():
          residue = conformer.only_residue()
          fl = rotamer_eval.evaluate_residue(residue = residue)
          if(residue.resname=="LYS"):
            assert fl.strip().upper()=="OUTLIER"
            found=True
  assert found
  #
  if(use_map): map_data = t.target_map
  else:        map_data = None
  result = mmtbx.refinement.real_space.fit_residues.run(
    pdb_hierarchy     = t.ph_poor,
    vdw_radii         = None,
    crystal_symmetry  = t.crystal_symmetry,
    map_data          = map_data,
    backbone_sample   = False,
    rotatable_hd      = t.rotatable_hd,
    rotamer_manager   = t.rotamer_manager,
    sin_cos_table     = t.sin_cos_table,
    mon_lib_srv       = t.mon_lib_srv)
  result.show()
  result.pdb_hierarchy.write_pdb_file(file_name = "refined_%s.pdb"%str(i_pdb),
    crystal_symmetry = t.crystal_symmetry)
  #
  if(use_map):
    mmtbx.refinement.real_space.check_sites_match(
      ph_answer  = t.ph_answer,
      ph_refined = result.pdb_hierarchy,
      tol        = 1.e-6)
  else:
    s1 = t.ph_answer.atoms().extract_xyz()
    s2 = result.pdb_hierarchy.atoms().extract_xyz()
    dist = flex.max(flex.sqrt((s1 - s2).dot()))
    assert dist>3.57, dist
  #
  found = False
  for model in result.pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        conformers = residue_group.conformers()
        for conformer in residue_group.conformers():
          residue = conformer.only_residue()
          fl = rotamer_eval.evaluate_residue(residue = residue)
          if(residue.resname=="LYS"):
            if(use_map): assert fl=="OUTLIER"
            else:        assert fl=="mmpt"
            found=True
  assert found

if(__name__ == "__main__"):
  t0 = time.time()
  for i_pdb, use_map in enumerate([True, False]):
    print("use_map", use_map, "-"*20)
    exercise(use_map=use_map, i_pdb=i_pdb)
  print("Time: %6.4f"%(time.time()-t0))
