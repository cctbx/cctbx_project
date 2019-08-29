# LIBTBX_SET_DISPATCHER_NAME mmtbx.screen_conformations

from __future__ import absolute_import, division, print_function
import iotbx.phil
from libtbx.utils import Sorry
from libtbx import easy_run
import time
import sys

master_phil = iotbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
selection = None
  .type = atom_selection
coot = False
  .type = bool
nproc = Auto
  .type = int
include scope mmtbx.building.alternate_conformations.density_sampling.master_params
""", process_includes=True)

def run(args, out=None ):
  if (out is None) : out = sys.stdout
  from mmtbx.building.alternate_conformations import density_sampling
  import mmtbx.maps.utils
  import mmtbx.utils
  from mmtbx.monomer_library import server
  import iotbx.pdb.hierarchy
  get_class = iotbx.pdb.common_residue_names_get_class
  mon_lib_srv = server.server()
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    process_pdb_file=False,
    scattering_table="n_gaussian")
  params = cmdline.params
  working_phil = master_phil.format(python_object=params)
  master_phil.fetch_diff(source=working_phil).show(out=out)
  fmodel = cmdline.fmodel
  hierarchy = cmdline.pdb_hierarchy
  sele_cache = hierarchy.atom_selection_cache()
  assert (params.selection is not None)
  selection = sele_cache.selection(params.selection)
  assert (selection.count(True) > 0)
  have_results = False
  pdb_atoms = hierarchy.atoms()
  sites_cart = pdb_atoms.extract_xyz()
  ensembles = []
  t1 = time.time()
  for chain in hierarchy.only_model().chains():
    prev_residue = next_residue = None
    residue_groups = chain.residue_groups()
    n_groups = len(residue_groups)
    for i_res, residue_group in enumerate(residue_groups):
      if (i_res < n_groups - 1):
        next_residue = residue_groups[i_res+1].atom_groups()[0]
      i_seqs = residue_group.atoms().extract_i_seq()
      if (selection.select(i_seqs).all_eq(True)):
        atom_group = residue_group.only_atom_group()
        if (get_class(atom_group.resname) != "common_amino_acid"):
          continue
        confs = density_sampling.screen_residue(
          residue=atom_group,
          prev_residue=prev_residue,
          next_residue=next_residue,
          sites_cart=sites_cart,
          fmodel=fmodel,
          mon_lib_srv=mon_lib_srv,
          params=params,
          map_file_name="map_coeffs.mtz",
          out=out)
        ensemble = []
        for conf in confs :
          conf_atoms = conf.set_atom_sites(pdb_atoms)
          pdb_lines = []
          for atom in conf_atoms :
            pdb_lines.append(atom.format_atom_record())
          ensemble.append("\n".join(pdb_lines))
        ensembles.append(ensemble)
      prev_residue = residue_group.atom_groups()[0]
  if (len(ensembles) == 0):
    raise Sorry("No alternate conformations found.")
  t2 = time.time()
  print("search time: %.1fs" % (t2-t1), file=out)
  for i_ens, ensemble in enumerate(ensembles):
    ensemble_hierarchy = iotbx.pdb.hierarchy.root()
    for k, model_str in enumerate(ensemble):
      input = iotbx.pdb.hierarchy.input(pdb_string=model_str)
      model = input.hierarchy.only_model().detached_copy()
      model.id = str(k+1)
      ensemble_hierarchy.append_model(model)
    f = open("ensemble_%d.pdb" % (i_ens+1), "w")
    f.write(ensemble_hierarchy.as_pdb_string())
    f.close()
    print("wrote ensemble_%d.pdb" % (i_ens+1))
  if (params.coot):
    easy_run.call("coot --pdb %s --auto map_coeffs.mtz --pdb ensemble.pdb" %
      params.input.pdb.file_name[0])

if (__name__ == "__main__"):
  run(sys.argv[1:])
