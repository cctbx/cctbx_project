
from libtbx.utils import Sorry
from cStringIO import StringIO
import os
import sys

morph_params_str = """
morph {
  pdb_file = None
    .type = path
    .multiple = True
  cif_file = None
    .type = path
    .multiple = True
    .help = Restraint (.cif) files for ligands etc.
  output_directory = None
    .type = path
    .help = Output directory (defaults to current working directory)
  output_prefix = morph
    .type = str
  serial_format = %03d
    .type = str
  frames = 30
    .type = ints
    .help = Number of frames for each transition (list of integers, or just \
            a single integer that applies to each transition)
  pause = 0
    .type = int
    .help = Number of frames to pause at each starting structure
  delete_waters = True
    .type = bool
  verbose = False
    .type = bool
  fitting {
    align_atoms = name CA
      .type = str
      .help = Selection of atoms to use in alignment ("None" to disable)
    reference_structure = 0
      .type = int
      .help = Index (C array notation) of model to use as reference in alignment
    sieve_fit = True
      .type = bool
      .help = Perform sieve-fitting of structures
  }
  minimization {
    minimize = True
      .type = bool
    n_min_steps = 60
      .type = int
    exclude_dihedrals = False
      .type = bool
    interpolate_dihedrals = False
      .type = bool
    dihedral_weight = 0.1
      .type = float
  }
}
pdb_interpretation {
  include scope mmtbx.monomer_library.pdb_interpretation.master_params
}
"""

def morph_models (params, out=None) :
  assert len(params.morph.pdb_file) > 1
  assert (len(params.morph.frames) == (len(params.morph.pdb_file) - 1))
  if (out is None) : out = sys.stdout
  from mmtbx.monomer_library import pdb_interpretation, server
  from iotbx import file_reader
  pdb_hierarchies = []
  for pdb_file in params.morph.pdb_file :
    pdb_in = file_reader.any_file(pdb_file, force_type="pdb")
    pdb_in.check_file_type("pdb")
    hierarchy = pdb_in.file_object.construct_hierarchy()
    pdb_hierarchies.append(hierarchy)
  new_pdb = homogenize_structures(
    pdb_hierarchies=pdb_hierarchies,
    delete_waters=params.morph.delete_waters)
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  for cif_file in params.morph.cif_file :
    print "Loading CIF file %s" % cif_file
    cif_object = server.read_cif(file_name=cif_file)
    mon_lib_serv.process_cif_object(cif_object=cif_object, file_name=cif_file)
    ener_lib.process_cif_object(cif_object=cif_object, file_name=cif_file)
  if (params.morph.minimization.interpolate_dihedrals) :
    params.pdb_interpretation.peptide_link.discard_psi_phi = False
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params.pdb_interpretation,
    pdb_inp=new_pdb[0],
    substitute_non_crystallographic_unit_cell_if_necessary=True)
  all_chain_proxies = processed_pdb_file.all_chain_proxies
  static_coords = [ all_chain_proxies.pdb_hierarchy.atoms().extract_xyz() ]
  for pdb_inp in new_pdb[1:] :
    sites = pdb_inp.atoms().extract_xyz()
    static_coords.append(sites)
  if (params.morph.fitting.align_atoms is not None) :
    print >> out, "Superposing on initial structure..."
    selection_cache = all_chain_proxies.pdb_hierarchy.atom_selection_cache()
    selection = selection_cache.selection(params.morph.fitting.align_atoms)
    if (selection.count(True) == 0) :
      raise Sorry("No atoms in alignment selection!")
    i_ref = params.morph.fitting.reference_structure
    sites_fixed = static_coords[i_ref]
    j = 0
    while (j < len(static_coords)) :
      if (j != i_ref) :
        sites_moving = static_coords[j]
        if (params.morph.fitting.sieve_fit) :
          from scitbx.math import superpose
          lsq_fit = superpose.sieve_fit(
            sites_fixed=sites_fixed,
            sites_moving=sites_moving,
            selection=selection)
          sites_moving_new = lsq_fit.other_sites_best_fit()
        else :
          sites_moving_new = fit_sites(
            sites_fixed=sites_fixed,
            sites_moving=sites_moving,
            selection=selection)
        static_coords[j] = sites_moving_new
      j += 1
  print "Ready to morph"
  morphs = []
  restraints_manager = processed_pdb_file.geometry_restraints_manager()
  for i in range(len(params.morph.pdb_file) - 1) :
    morph = adiabatic_mapping(
      pdb_hierarchy=all_chain_proxies.pdb_hierarchy,
      restraints_manager=restraints_manager,
      start_coords = static_coords[i],
      end_coords = static_coords[i+1],
      params = params.morph.minimization,
      nsteps = params.morph.frames[i])
    morphs.append(morph)
  serial = 1
  if (params.morph.output_directory is not None) :
    output_base = os.path.join(params.morph.output_directory,
      params.morph.output_prefix)
  else :
    output_base = params.morph.output_prefix
  for i, morph in enumerate(morphs) :
    serial = morph.write_pdb_files(
      output_base=output_base,
      serial=serial,
      serial_format=params.morph.serial_format,
      pause=params.morph.pause,
      pause_at_end=(i == (len(morphs) - 1)))
  f = open("%s.pml" % output_base, "w")
  for i in range(1, serial) :
    print >> f, "load %s_%04d.pdb, morph" % (output_base, i)
  f.close()
  print >> out, "PyMOL script is %s.pml" % output_base

def homogenize_structures (pdb_hierarchies,
                           delete_heteroatoms=False,
                           delete_waters=True) :
  """
  Eliminate atoms not present in all models.  This ignores residue names, so
  corresponding mainchain atoms should be preserved in most cases.
  """
  import iotbx.pdb
  chain_ids = []
  residues_by_chain = {}
  for hierarchy in pdb_hierarchies :
    if (len(hierarchy.models()) > 1) :
      raise Sorry("Multiple model PDB files not supported.")
  atom_sets = []
  for hierarchy in pdb_hierarchies :
    atom_set = set([])
    for atom in hierarchy.atoms() :
      labels = atom.fetch_labels()
      if (delete_waters) and (labels.resname in ["HOH", "WAT"]) :
        continue
      atom_info = (labels.chain_id, labels.resid(), atom.name, labels.altloc)
      atom_set.add(atom_info)
    atom_sets.append(atom_set)
  common_set = atom_sets[0].intersection(*(atom_sets[1:]))
  n_atoms = len(common_set)
  if (n_atoms == 0) :
    raise RuntimeError("No atoms left in structure.")
  print "%d atoms in common." % len(common_set)
  atom_lists = []
  for hierarchy in pdb_hierarchies :
    atom_set = set([])
    atoms = hierarchy.atoms()
    i = 0
    while i < len(atoms) :
      atom = atoms[i]
      labels = atom.fetch_labels()
      atom_info = (labels.chain_id, labels.resid(), atom.name, labels.altloc)
      if (not atom_info in common_set) :
        del atoms[i]
      else :
        i += 1
    atom_lists.append(atoms)
  index_lookup = {}
  atoms_out = [[]]
  for i, atom in enumerate(atom_lists[0]) :
    labels = atom.fetch_labels()
    atom_info = (labels.chain_id, labels.resid(), atom.name, labels.altloc)
    index_lookup[atom_info] = i
    atoms_out[0].append(labels.format_atom_record())
  for atoms in atom_lists[1:] :
    assert (len(atoms) == len(atom_lists[0]))
    sorted_atoms = [None] * len(atoms)
    for atom in atoms :
      labels = atom.fetch_labels()
      atom_info = (labels.chain_id, labels.resid(), atom.name, labels.altloc)
      i = index_lookup[atom_info]
      sorted_atoms[i] = labels.format_atom_record()
    assert (len(sorted_atoms) == len(atoms))
    assert (not None in sorted_atoms)
    atoms_out.append(sorted_atoms)
  pdb_out = []
  for atoms_strings in atoms_out :
    pdb_in = iotbx.pdb.input(source_info=None, lines=atoms_strings)
    pdb_out.append(pdb_in)
  return pdb_out

def fit_sites (sites_fixed, sites_moving, selection) : # TODO
  """
  Simple least-squares superposition of sites on reference structure
  """
  from scitbx.math import superpose
  sites_fixed_aln = sites_fixed.select(selection)
  sites_moving_aln = sites_moving.select(selection)
  lsq_fit_obj = superpose.least_squares_fit(
    reference_sites=sites_fixed_aln,
    other_sites=sites_moving_aln)
  sites_moved = lsq_fit_obj.r.elems * sites_moving + lsq_fit_obj.t.elems
  return sites_moved

class morph (object) :
  """
  Container for interpolated sites.
  """
  def __init__ (self, pdb_hierarchy) :
    self.pdb_hierarchy = pdb_hierarchy.deep_copy()
    self._frames = []

  def add_frame (self, sites_cart) :
    self._frames.append(sites_cart)

  def _write_pdb (self, file_name) :
    f = open(file_name, "w")
    f.write(self.pdb_hierarchy.as_pdb_string())
    f.close()
    print "  wrote %s" % os.path.basename(file_name)

  def write_pdb_files (self, output_base, serial, serial_format="%04d",
      pause=0, pause_at_end=False) :
    file_format = "%s_%s.pdb" % (output_base, serial_format)
    k = serial
    if (pause != 0) :
      for j in range(pause) :
        self.pdb_hierarchy.atoms().set_xyz(self._frames[0])
        file_name = file_format % k
        self._write_pdb(file_name)
        k += 1
    for sites in self._frames :
      self.pdb_hierarchy.atoms().set_xyz(sites)
      file_name = file_format % k
      self._write_pdb(file_name)
      k += 1
    if (pause_at_end) and (pause != 0) :
      for j in range(pause) :
        self.pdb_hierarchy.atoms().set_xyz(self._frames[-1])
        file_name = file_format % k
        self._write_pdb(file_name)
        k += 1
    return k

def adiabatic_mapping (pdb_hierarchy,
                       restraints_manager,
                       start_coords,
                       end_coords,
                       params,
                       nsteps=10,
                       verbose=False,
                       out=None) :
  """
  Linear interpolation with energy minimization.  The number of minimizer
  cycles should be kept low to prevent each "frame" from snapping back to the
  starting coordinates.  Interpolating the dihedral restraint target angles
  may help remediate this problem.
  """
  from mmtbx.command_line import geometry_minimization
  from cctbx import geometry_restraints
  import scitbx.lbfgs
  assert (start_coords.size() == end_coords.size())
  grm = restraints_manager
  grm_flags = geometry_restraints.flags.flags(
    default=True,
    dihedral=(not params.exclude_dihedrals))
  term_params = scitbx.lbfgs.termination_parameters(
    max_iterations=params.n_min_steps)
  dihedral_steps = []
  n = nsteps - 1
  if (params.interpolate_dihedrals) :
    for proxy in grm.dihedral_proxies :
      sites1 = []
      sites2 = []
      for k in range(4) :
        sites1.append(start_coords[proxy.i_seqs[k]])
        sites2.append(end_coords[proxy.i_seqs[k]])
      dihe_start = geometry_restraints.dihedral(
        sites=sites1,
        angle_ideal=0,
        weight=1).angle_model
      dihe_end = geometry_restraints.dihedral(
        sites=sites2,
        angle_ideal=0,
        weight=1).angle_model
      delta = (dihe_end - dihe_start) / n
      dihedral_steps.append(delta)
      proxy.weight = 0.1
  m = morph(pdb_hierarchy)
  m.add_frame(start_coords)
  current_xyz = start_coords
  final_xyz = end_coords
  while n > 0 :
    print >> out, "Interpolation step %d" % (nsteps - n)
    if (verbose) :
      print >> out, ""
    new_xyz = current_xyz.deep_copy()
    dxyz = (final_xyz - new_xyz) * (1. / n)
    new_xyz += dxyz
    if (params.minimize) :
      if (params.interpolate_dihedrals) :
        for k, proxy in enumerate(grm.dihedral_proxies) :
          proxy.angle_ideal += dihedral_steps[k]
      minimized = geometry_minimization.lbfgs(
        sites_cart=new_xyz,
        geometry_restraints_manager=grm,
        geometry_restraints_flags=grm_flags,
        lbfgs_termination_params=term_params)
      if (verbose) :
        n_iter = minimized.minimizer.iter()
        print >> out, ""
        print >> out, "Number of minimization iterations:", n_iter
        print >> out, ""
        print >> out, "Energies at end of minimization:"
        minimized.final_target_result.show(f=out)
    print >> out, "RMS coordinate change: %.3f" % current_xyz.rms_difference(new_xyz)
    m.add_frame(new_xyz)
    n -= 1
    current_xyz = new_xyz
  m.add_frame(final_xyz)
  return m

def run (args, out=None) :
  import iotbx.phil
  if (out is None) : out = sys.stdout
  usage_str = "mmtbx.interpolate model1.pdb model2.pdb [...modelX.pdb] [options]"
  if (len(args) == 0) or ("--help" in args) :
    params_out = StringIO()
    master_phil = iotbx.phil.parse(morph_params_str)
    master_phil.show(out=params_out)
    print >> out, """
mmtbx.interpolate - simple morphing with energy minimization

Limitations:
  1. Any atoms not present in all structures will be deleted.
  2. Equivalent atoms must have identical chain IDs and residue numbers;
     however, residue names do not have to be the same.
  3. This method usually doesn't work very well on conformational changes that
     involve major rearrangements and/or refolding.

Usage:
  %s

Full parameter list:
  %s
""" % (usage_str, params_out.getvalue())
    return None
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=morph_params_str,
    pdb_file_def="morph.pdb_file",
    cif_file_def="morph.cif_file")
  working_phil = cmdline.work
  params = working_phil.extract()
  eff_out = open("%s.eff" % params.morph.output_prefix, "w")
  print >> out, "Writing effective parameters to %s.eff" % \
    params.morph.output_prefix
  working_phil.show(out=eff_out)
  #working_phil.show(out=sys.stdout)
  eff_out.close()
  morph_models(params, out=out)

if __name__ == "__main__" :
  run(sys.argv[1:])

#---end
