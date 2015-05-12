
from __future__ import division
from libtbx.utils import null_out
import libtbx.load_env
import warnings
import os
import sys

def exercise_assembly (verbose=False) :
  # This test subclasses the main class in sliding_window to substitute
  # fake annealing results with predefined coordinate offsets.  The results
  # will then be processed and assembled as they were actual results.
  from mmtbx.command_line import build_alternate_conformations
  from mmtbx.building.alternate_conformations import sliding_window
  from mmtbx.building import alternate_conformations as alt_confs
  from mmtbx import building
  from scitbx.array_family import flex
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/CypA_refine_3.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/3k0n.mtz",
    test=os.path.isfile)
  assert (not None in [pdb_file, mtz_file])
  class fragment_driver (sliding_window.fragment_refinement_driver) :
    def refine_window (O, window) :
      processed_pdb_file = O.get_processed_pdb_file()
      assert (processed_pdb_file is not None)
      hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
      pdb_atoms = hierarchy.atoms()
      sites_cart = pdb_atoms.extract_xyz()
      dxyz = [ None ] * 3
      dxyz[0] = flex.vec3_double(pdb_atoms.size(), (1.0,0.2,0.3))
      dxyz[1] = flex.vec3_double(pdb_atoms.size(), (-1.0,-0.1,0.0))
      dxyz[2] = dxyz[1]
      trials = []
      if (window.residue_id_str in [" A PHE 113 ", " A SER 110 "]) :
        for i in range(3) :
          trials.append(alt_confs.trial_result(
            sites_cart=(sites_cart + dxyz[i]).select(window.selection),
            min_fofc=4.0,
            mean_fofc=4.5,
            rmsd=1.2,
            max_dev=1.3,
            cc=0.99))
      elif (window.residue_id_str in [" A ASN 108 "]) :
        trials.append(alt_confs.trial_result(
          sites_cart=(sites_cart + dxyz[0]).select(window.selection),
          min_fofc=4.0,
          mean_fofc=4.5,
          rmsd=1.2,
          max_dev=1.3,
          cc=0.99))
      e = sliding_window.ensemble(window=window, sites_trials=trials)
      n_keep = e.filter_trials(
        sites_cart=O.sites_cart,
        min_rmsd=O.params.min_rmsd,
        min_dev=O.min_required_deviation)
      if (n_keep > 0) :
        return e
      return None
  class program_driver (build_alternate_conformations.build_and_refine) :
    def refine (O, title=None) : pass
    def rejoin (O) : return False
    def build_conformers (O, stop_if_none=None) :
      driver = fragment_driver(
        fmodel=O.fmodel,
        pdb_hierarchy=O.pdb_hierarchy,
        processed_pdb_file=O.processed_pdb_file,
        params=O.params.sliding_window,
        mp_params=O.params,
        out=O.out)
      O.pdb_hierarchy = driver.assemble()
      O.processed_pdb_file = None
  args = [ pdb_file, mtz_file, "nproc=2", "high_resolution=2.0",
    "create_dir=False", "output.prefix=tst_sliding_window" ]
  out = null_out()
  if (verbose) :
    out = sys.stdout
  hierarchy = build_alternate_conformations.run(args=args,
    driver_class=program_driver,
    out=out)
  for residue_group in building.iter_residue_groups(hierarchy) :
    n_atom_groups = len(residue_group.atom_groups())
    if (n_atom_groups == 3) :
      assert (106 <= residue_group.resseq_as_int() <= 115)
    else :
      assert (n_atom_groups == 1)

if (__name__ == "__main__") :
  if (not libtbx.env.has_module("phenix_regression")) :
    warnings.warn("phenix_regression missing, skipping test")
  else :
    exercise_assembly(verbose=("--verbose" in sys.argv))
    print "OK"
