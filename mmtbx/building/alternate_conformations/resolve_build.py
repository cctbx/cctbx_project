
# XXX requires solve_resolve (inline import)

from __future__ import absolute_import, division, print_function
from mmtbx import building
from libtbx.str_utils import make_sub_header
from libtbx.utils import null_out
import time
import os

resolve_build_params_str = """
build_new_loop = False
  .type = bool
n_random_loop = 200
  .type = int
anneal = False
  .type = bool
"""

class resolve_builder(building.box_build_refine_base):
  def __init__(self, params, *args, **kwds):
    self.params = params
    building.box_build_refine_base.__init__(self, *args, **kwds)
    self.box_map_coeffs = self.box.box_map_coefficients(
      d_min=self.d_min,
      resolution_factor=self.resolution_factor)
    self.box.write_pdb_file("box_start.pdb")
    self.box.write_ccp4_map()
    self.run_resolve()

  def run_resolve(self):
    from solve_resolve.resolve_python import resolve_in_memory
    from iotbx import pdb
    from scitbx.array_family import flex
    make_sub_header("RESOLVE build", out=self.out)
    mean_density_start = self.mean_density_at_sites()
    cc_start = self.cc_model_map()
    sites_start = self.get_selected_sites(hydrogens=False)
    t1 = time.time()
    pdb_inp = self.box_selected_hierarchy.as_pdb_input()
    inp_hierarchy = pdb_inp.construct_hierarchy()
    chain = inp_hierarchy.only_model().only_chain()
    first_resseq = chain.residue_groups()[0].resseq_as_int()
    seq = "".join(chain.only_conformer().as_sequence(substitute_unknown='A'))
    resolve_args = [
      "start_chain 1 %d" % first_resseq,
      "extend_only",
      "skip_hetatm",
      "no_merge_ncs_copies",
      "no_optimize_ncs",
      "i_ran_seed %d" % int(time.time() % os.getpid()),
    ]
    if (self.params.build_new_loop) : # XXX not really working...
      n_res = len(chain.residue_groups())
      assert (n_res >= 3)
      k = 0
      for residue_group in chain.residue_groups()[1:-1] :
        print("  removing residue group %s %s" % \
          (chain.id, residue_group.resid()), file=self.out)
        chain.remove_residue_group(residue_group)
      resolve_args.extend([
        "loop_only",
        "build_outside_model",
        "no_sub_segments",
        "n_random_loop %d" % self.params.n_random_loop,
        "loop_length %d" % (n_res - 2),
        "rms_random_loop 0.3",
        "rho_min_main_low 0.5",
        "rho_min_main_base 0.5",
        "n_internal_start 0",
      ])
    else :
      resolve_args.extend([
        "rebuild_in_place",
        "replace_existing",
        "richardson_rotamers",
        "min_z_value_rho -3.0",
        "delta_phi   20.00",
        "dist_cut_base 3.0",
        "n_random_frag 0",
        "group_ca_length 4",
        "group_length 2",
      ])
    out = null_out()
    if (self.debug):
      out = self.out
    cmn = resolve_in_memory.run(
      map_coeffs=self.box_map_coeffs,
      pdb_inp=inp_hierarchy.as_pdb_input(),
      build=True,
      input_text="\n".join(resolve_args),
      chain_type="PROTEIN",
      seq_file_as_string=seq,
      out=out)
    new_pdb_input = pdb.input(
      source_info='string',
      lines=flex.split_lines(cmn.atom_db.pdb_out_as_string))
    new_hierarchy = new_pdb_input.construct_hierarchy()
    print("  %d atoms rebuilt" % len(new_hierarchy.atoms()), file=self.out)
    new_hierarchy.write_pdb_file("resolve.pdb")
    selection_moved = flex.size_t()
    sites_new = flex.vec3_double()
    for atom in new_hierarchy.atoms():
      id_str = atom.id_str()
      if (not id_str in self.atom_id_mapping):
        raise KeyError("Atom ID %s not recognized in RESOLVE model." % id_str)
      i_seq = self.atom_id_mapping[id_str]
      selection_moved.append(i_seq)
      sites_new.append(atom.xyz)
    sites_cart_selected = self.box_selected_hierarchy.atoms().extract_xyz()
    sites_cart_selected.set_selected(selection_moved, sites_new)
    self.box_selected_hierarchy.atoms().set_xyz(sites_cart_selected)
    sites_cart_box = self.box.xray_structure_box.sites_cart()
    sites_cart_box.set_selected(self.selection_in_box, sites_cart_selected)
    self.box.xray_structure_box.set_sites_cart(sites_cart_box)
    self.box.pdb_hierarchy_box.atoms().set_xyz(sites_cart_box)
    t2 = time.time()
    print("  RESOLVE time: %.1fs" % (t2-t1), file=self.out)
    selection_rebuilt = self.selection_in_box.select(selection_moved)
    minimize_sel = flex.bool(self.n_sites_box, False).set_selected(
      self.selection_in_box, True).set_selected(selection_rebuilt, False)
    # atoms present in the selection but not in the RESOLVE model (usually
    # hydrogen atoms) need to be minimized to follow the rebuilt sites
    if (minimize_sel.count(True) > 0):
      print("  Performing geometry minimzation on unbuilt sites", file=self.out)
      self.geometry_minimization(
        selection=minimize_sel,
        nonbonded=False)
    self.box.write_pdb_file("box_resolve.pdb")
    # two alternatives here: restrain other atoms tightly, and minimize the
    # entire box, or restrain selected atoms loosely, and refine only those
    self.restrain_atoms(
      selection=self.others_in_box,
      reference_sigma=0.02)
    if (self.params.anneal):
      self.anneal(start_temperature=2500)
    else :
      self.real_space_refine(selection=self.selection_all_box)
    self.box.write_pdb_file("box_refined.pdb")
    self.box.write_ccp4_map()
    mean_density_end = self.mean_density_at_sites()
    cc_end = self.cc_model_map()
    print("  mean density level: start=%.2fsigma  end=%.2fsigma" \
      % (mean_density_start, mean_density_end), file=self.out)
    print("  model-map CC: start=%.3f  end=%.3f" % (cc_start,
      cc_end), file=self.out)
    sites_final = self.get_selected_sites(hydrogens=False)
    print("  rmsd to starting model: %.3f Angstrom" % \
      sites_final.rms_difference(sites_start), file=self.out)
    t3 = time.time()
    print("  Total build and refine time: %.1fs" % (t3-t1), file=self.out)
