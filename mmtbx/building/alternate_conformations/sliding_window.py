
# TODO handle existing alternate conformations more sensibly

from __future__ import absolute_import, division, print_function
from mmtbx.building import alternate_conformations as alt_confs
import mmtbx.building
from libtbx import adopt_init_args, Auto, slots_getstate_setstate
from libtbx.utils import null_out
import string
import sys

master_params_str = """
window_radius = 2
  .type = int(value_min=1, value_max=4)
min_rmsd = 0.5
  .type = float
min_deviation = Auto
  .type = float
exclude_resnames = None
  .type = strings
include scope mmtbx.building.alternate_conformations.real_space_annealing.master_params_str
prefilter {
  rotameric_only = False
    .type = bool
  use_difference_map = True
    .type = bool
  sampling_radius = 2.5
    .type = float
}
assembly {
  primary_occupancy = 0.5
    .type = float
}
"""

class fragment_refinement_driver(object):
  def __init__(self,
      fmodel,
      params,
      mp_params,
      pdb_hierarchy,
      processed_pdb_file,
      selection=None,
      cif_objects=(),
      verbose=True,
      debug=False,
      out=None):
    if (out is None) : out = sys.stdout
    adopt_init_args(self, locals())
    self.asynchronous_output = False
    from mmtbx.rotamer import rotamer_eval
    from scitbx.array_family import flex
    assert (processed_pdb_file is not None) or (len(pdb_file_names) > 0)
    assert (0 < self.params.window_radius <= 4)
    self.pdb_hierarchy = pdb_hierarchy
    self.processed_pdb_file = processed_pdb_file
    self.get_processed_pdb_file(log=out)
    self.sites_cart = self.pdb_hierarchy.atoms().extract_xyz().deep_copy()
    if (self.selection is None):
      self.selection = flex.bool(self.sites_cart.size(), True)
    self.min_required_deviation = self.params.min_deviation
    if (self.min_required_deviation is Auto):
      self.min_required_deviation = fmodel.f_obs().d_min() / 2
    self._ensembles = []
    self.nproc_1 = self.nproc_2 = 1
    two_fofc_map, fofc_map = mmtbx.building.get_difference_maps(fmodel=fmodel)
    windows = []
    r = rotamer_eval.RotamerEval(data_version="8000")
    exclude_resnames = []
    if (params.exclude_resnames is not None):
      exclude_resnames = [ n.upper() for n in params.exclude_resnames ]
    for chain in self.pdb_hierarchy.only_model().chains():
      if not chain.is_protein():
        continue
      residues = chain.residue_groups()
      fragments = alt_confs.fragment_single_conformer_chain(residues)
      for fragment_residues in fragments :
        start = params.window_radius
        end = - params.window_radius
        for i_res, residue in enumerate(fragment_residues[start:end]):
          j_res = i_res + start
          atom_groups = residue.atom_groups()
          main_conf = atom_groups[0]
          if (main_conf.resname.upper() in exclude_resnames):
            continue
          residue_id = main_conf.id_str()
          ag_i_seqs = main_conf.atoms().extract_i_seq()
          if (not self.selection.select(ag_i_seqs).all_eq(True)):
            continue
          if (len(atom_groups) != 1):
            if (self.verbose):
              print("  residue %s already has multiple conformations"%\
                residue_id, file=out)
            continue
          ag_i_seqs_no_hd = flex.size_t()
          for atom in main_conf.atoms():
            if (atom.element.strip() not in ["H","D"]):
              ag_i_seqs_no_hd.append(atom.i_seq)
          # XXX this is probably not optimal; what should I do about the
          # adjacent residues?  it would be good to check Ramachandran plot too
          if (self.params.prefilter.rotameric_only):
            n_outliers = alt_confs.score_rotamers(hierarchy=hierarchy,
              selection=ag_i_seqs)
            if (n_outliers > 0):
              if (self.verbose):
                print("  residue %s is a rotamer outlier" % residue_id, file=out)
              continue
          if (self.params.prefilter.use_difference_map):
            map_stats = building.local_density_quality(
              fofc_map=fofc_map,
              two_fofc_map=two_fofc_map,
              atom_selection=ag_i_seqs_no_hd,
              xray_structure=fmodel.xray_structure,
              radius=self.params.prefilter.sampling_radius)
            if ((map_stats.number_of_atoms_in_difference_holes() == 0) and
                (map_stats.fraction_of_nearby_grid_points_above_cutoff()==0)):
              if (self.verbose):
                print("  no difference density for %s" % residue_id, file=out)
              continue
          window_selection = flex.size_t()
          offset = - self.params.window_radius
          while (offset <= self.params.window_radius):
            adjacent_group = fragment_residues[j_res+offset].atom_groups()[0]
            window_selection.extend(adjacent_group.atoms().extract_i_seq())
            offset += 1
          windows.append(residue_window(
            residue_id_str=residue_id,
            selection=window_selection,
            residue_selection=ag_i_seqs_no_hd,
            sites_reference=self.sites_cart.select(selection),
            window_radius=self.params.window_radius))
    if (len(windows) == 0):
      raise Sorry("No peptide segments meeting the filtering criteria could "+
        "be extracted from the selected atoms.")
    else :
      print("%d fragments will be refined." % len(windows), file=out)
    if (self.mp_params.nproc == 1):
      pass
    elif (self.mp_params.technology == "multiprocessing"):
      if (self.params.n_trials == 1) and (len(self.params.partial_occupancy) == 1):
        # only one refinement per window, so parallelize residue iteration
        self.nproc_1 = self.mp_params.nproc
      else :
        # multiple refinements per window, so parallelize at that level
        # FIXME actually, this needs to be smarter - if the number of
        # available processors is greater than the number of refinements per
        # window, it will be more efficient to parallelize the window loop
        self.nproc_2 = self.mp_params.nproc
    else :
      # queuing system, so we can only parallelize residue iteration
      self.nproc_1 = self.mp_params.nproc
      self.out = null_out()
      self.processed_pdb_file = None
    print("", file=out)
    alt_confs.print_trial_header(out)
    ensembles = []
    if (self.nproc_1 == 1):
      self.asynchronous_output = True
      for window in windows :
        ens = self.refine_window(window)
        ensembles.append(ens)
    else :
      ensembles = easy_mp.parallel_map(
        func=self.refine_window,
        iterable=windows,
        processes=self.nproc_1,
        qsub_command=mp_params.qsub_command,
        method=mp_params.technology)
    self._ensembles = [ e for e in ensembles if (e is not None) ]
    # XXX reassert order
    print("", file=out)
    if (len(self._ensembles) == 0):
      print("WARNING: no ensembles passed filtering step", file=out)
      print("", file=out)
    self._ensembles.sort(lambda a,b: a.selection[0] < b.selection[0])
    self.processed_pdb_file = processed_pdb_file
    if (debug):
      for k, ens in enumerate(filtered):
        pdb_out = ens.dump_pdb_file(
          pdb_hierarchy=pdb_hierarchy,
          crystal_symmetry=fmodel.f_obs())
        print("wrote %s" % pdb_out, file=out)

  def n_ensembles(self):
    return len(self._ensembles)

  def get_processed_pdb_file(self, log=None):
    if (self.processed_pdb_file is None):
      if (log is None):
        log = null_out()
      import mmtbx.building
      self.processed_pdb_file = mmtbx.building.reprocess_pdb(
        pdb_hierarchy=self.pdb_hierarchy,
        crystal_symmetry=self.crystal_symmetry,
        cif_objects=self.cif_objects,
        log=log)
    return self.processed_pdb_file

  def refine_window(self, window):
    from mmtbx.building.alternate_conformations import real_space_annealing
    processed_pdb_file = self.get_processed_pdb_file()
    assert (processed_pdb_file is not None)
    hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
    #print window.residue_id_str
    log = null_out()
    if (self.debug):
      log = sys.stdout
    refinements = real_space_annealing.refine_into_difference_density(
      fmodel=self.fmodel,
      pdb_hierarchy=hierarchy.deep_copy(),
      processed_pdb_file=processed_pdb_file,
      selection=window.selection,
      selection_score=window.residue_selection,
      params=self.params,
      nproc=self.nproc_2,
      out=null_out()).get_filtered_trials(log=log)
    result = ensemble(
      window=window,
      sites_trials=refinements)
    n_keep = result.filter_trials(
      sites_cart=self.sites_cart,
      min_rmsd=self.params.min_rmsd,
      min_dev=self.min_required_deviation)
    if (n_keep > 0):
      if (self.asynchronous_output):
        print(result, file=self.out)
      return result
    return None

  def group_ensembles(self):
    groups = []
    current_group = last_ensemble = None
    for ens in self._ensembles :
      # this will check for overlap
      if ((last_ensemble is None) or
          (alt_confs.get_selection_gap(last_ensemble.selection,
                                       ens.selection) > 0)):
        current_group = [ ens ]
        groups.append(current_group)
      else :
        current_group.append(ens)
      last_ensemble = ens
    return [ ensemble_group(ensembles) for ensembles in groups ]

  def assemble(self, out=None):
    if (out is None) : out = self.out
    return assemble_fragments(
      pdb_hierarchy=self.get_processed_pdb_file().all_chain_proxies.pdb_hierarchy,
      groups=self.group_ensembles(),
      out=out)

class residue_window(slots_getstate_setstate):
  """
  Container for information about a sliding window along a protein chain.
  Uniquely identified by ID of central residue.
  """
  __slots__ = ["residue_id_str", "selection", "residue_selection",
    "sites_reference", "window_radius"]
  def __init__(self, **kwds):
    kwds = dict(kwds)
    for name in self.__slots__ :
      setattr(self, name, kwds[name])

  def __str__(self):
    return "window: %s [%d atoms in %d residues]" % (self.residue_id_str,
      len(self.selection), (2*self.window_radius)+1)

class ensemble(residue_window):
  """
  Container for results of refinement of a single sliding window.  The
  hierarchy will be a multi-model ensemble if more than one trial was run.
  """
  __slots__ = ["residue_id_str", "selection", "residue_selection",
    "sites_reference", "window_radius", "sites_trials"]
  def __init__(self,
      window,
      sites_trials):
    assert (isinstance(sites_trials, list))
    residue_window.__init__(self,
      residue_id_str=window.residue_id_str,
      selection=window.selection,
      residue_selection=window.residue_selection,
      sites_reference=window.sites_reference,
      window_radius=window.window_radius,
      sites_trials=sites_trials)

  def show_summary(self, out=None, comprehensive=False):
    if (out is None) : out = sys.stdout
    print(str(self), file=out)
    rmsds = []
    ccs = []
    for trial in self.sites_trials :
      rmsds.append(trial.rmsd)
      ccs.append(trial.cc)
    print("  %d conformation(s) generated" % len(self.sites_trials), file=out)
    if (len(rmsds) > 0):
      print("    max. rmsd:  %7.3f" % max(rmsds), file=out)
      print("    min. rmsd:  %7.3f" % min(rmsds), file=out)
      print("    best CC:    %7.3f" % max(ccs), file=out)
      if (comprehensive):
        for k, trial in enumerate(self.sites_trials):
          print("    trial %d:" % (k+1), file=out)
          trial.show_summary(out=out, prefix="      ")

  def __str__(self):
    lines = []
    prefix = "%s  " % (self.residue_id_str)
    for i, trial in enumerate(self.sites_trials):
      lines.append(prefix + ("%4d  " % (i+1)) + str(trial))
    return "\n".join(lines)

  def as_pdb_hierarchy(self, pdb_hierarchy):
    if (len(self.sites_trials) == 0):
      return None
    import iotbx.pdb.hierarchy
    root = iotbx.pdb.hierarchy.root()
    for k, trial in enumerate(self.sites_trials):
      new_hierarchy = pdb_hierarchy.deep_copy()
      pdb_atoms = new_hierarchy.atoms()
      sites_cart = pdb_atoms.extract_xyz()
      sites_cart.set_selected(self.selection, trial.sites_cart)
      pdb_atoms.set_xyz(sites_cart)
      new_model = new_hierarchy.only_model().detached_copy()
      new_model.id = str(k+1)
      root.append_model(new_model)
    return root

  def dump_pdb_file(self, pdb_hierarchy, file_prefix="ensemble",
      crystal_symmetry=None,):
    window_id = "_".join(self.residue_id_str.split())
    file_name = "%s_%s.pdb" % (file_prefix, window_id)
    new_hierarchy = self.as_pdb_hierarchy(pdb_hierarchy)
    if (new_hierarchy is not None):
      f = open(file_name, "w")
      f.write(new_hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry))
      f.close()
      return file_name
    return None

  def max_rmsd(self):
    rmsds = [ t.rmsd for t in self.sites_trials ]
    return max(rmsds)

  def max_deviation(self):
    dev = [ t.max_dev for t in self.sites_trials ]
    return max(dev)

  def filter_trials(self, sites_cart, min_rmsd, min_dev, cluster=True):
    """
    Reduce the trials to only those which are significantly different from
    the input conformation and each other.
    """
    filtered_trials = []
    for trial in self.sites_trials :
      if (trial.rmsd >= min_rmsd) and (trial.max_dev >= min_dev):
        filtered_trials.append(trial)
    # FIXME this is pretty dumb logic...
    clusters = []
    cluster_sites = []
    for trial in filtered_trials :
      sites_trial = sites_cart.deep_copy()
      sites_trial.set_selected(self.selection, trial.sites_cart)
      sites_compare = sites_trial.select(self.residue_selection)
      if (len(clusters) == 0):
        clusters.append(trial)
        cluster_sites.append(sites_compare)
        continue
      k = 0
      while k < len(clusters):
        other = clusters[k]
        other_sites = cluster_sites[k]
        for other, other_sites in zip(clusters, cluster_sites):
          if ((other_sites.max_distance(sites_compare) < min_dev) or
              (other_sites.rms_difference(sites_compare) < min_rmsd)):
            if (other.cc < trial.cc) : # or should I pick the best mFo-DFc?
              clusters[k] = trial
              cluster_sites[k] = sites_compare
            else :
              pass
            break
        else :
          clusters.append(trial)
          cluster_sites.append(sites_compare)
        k += 1
    self.sites_trials = clusters
    return self.n_confs()

  def n_confs(self):
    return len(self.sites_trials)

class ensemble_group(object):
  """
  Collection of ensembles which overlap in some way.  The assumption is that
  each ensemble is unique to a specific central residue, which is enforced in
  the code.
  """
  def __init__(self, ensembles):
    self.ensembles = ensembles

  def n_confs_max(self):
    return max([ e.n_confs() for e in self.ensembles ])

  def n_confs_min(self):
    return min([ e.n_confs() for e in self.ensembles ])

  def show_summary(self, out=sys.stdout):
    print("Fragment (%d ensembles):" % len(self.ensembles), file=out)
    print("   %s --> %s" % (self.ensembles[0].residue_id_str,
      self.ensembles[-1].residue_id_str), file=out) # FIXME assumes sorting!
    print("   %d conformations" % self.n_confs_max(), file=out)

  def get_central_residues_selection(self):
    from scitbx.array_family import flex
    selection = flex.size_t()
    for ens in self.ensembles :
      assert (len(selection.intersection(ens.residue_selection)) == 0)
      if (len(selection) > 0):
        if (ens.residue_selection[0] <= selection[-1]):
          raise RuntimeError("Residue %s is out of order" % ens.residue_id_str)
      selection.extend(ens.residue_selection)
    return selection

  def get_entire_selection(self):
    from scitbx.array_family import flex
    selection = []
    for ens in self.ensembles :
      selection.extend(ens.selection)
    return flex.size_t(sorted(list(set(selection))))

  def combine_ensembles(self, pdb_hierarchy):
    """
    Iterate over all residues covered by the group, and pick the alternate
    conformers out for each residue to combine into a continous fragment.
    If additional conformers are needed, the previously existing atom group
    will be duplicated as necessary.
    """
    from mmtbx import building
    from scitbx.array_family import flex
    hierarchy = pdb_hierarchy.deep_copy()
    pdb_atoms = hierarchy.atoms()
    pdb_atoms.reset_i_seq()
    n_atoms = len(pdb_atoms)
    covered_selection = flex.bool(pdb_atoms.size(), False)
    sites_cart = pdb_atoms.extract_xyz().deep_copy()
    residue_groups = []
    n_max = self.n_confs_max()
    selection = self.get_entire_selection()
    central_selection = self.get_central_residues_selection()
    sub_hierarchy = hierarchy.deep_copy().select(selection)
    sub_hierarchy.atoms().reset_i_seq()
    sub_iselection = pdb_atoms.extract_i_seq().select(selection)
    sub_central_selection = flex.bool(pdb_atoms.size(),
      False).set_selected(central_selection, True).select(sub_iselection)
    fragment_atom_groups = []
    fragment_i_seqs = []
    for residue_group in mmtbx.building.iter_residue_groups(sub_hierarchy):
      new_atom_groups = []
      n_confs = 0
      atom_group = residue_group.only_atom_group()
      ag_new_isel = atom_group.atoms().extract_i_seq()
      ag_old_isel = sub_iselection.select(ag_new_isel)
      fragment_i_seqs.append(ag_old_isel)
      original_sel = flex.bool(n_atoms, False).set_selected(ag_old_isel, True)
      def new_atom_group(trial, selection):
        sites_new = sites_cart.deep_copy()
        sites_new.set_selected(selection, trial.sites_cart)
        new_ag = atom_group.detached_copy()
        new_ag.atoms().set_xyz(sites_new.select(original_sel))
        #for k, atom in enumerate(new_ag.atoms()):
        #  atom.i_seq = ag_old_isel[k]
        return new_ag
      if (sub_central_selection.select(ag_new_isel).all_eq(True)):
        # residue is the center of the window for a single ensemble, so use
        # the coordinates from that enemble
        found_residue = False
        for ens in self.ensembles :
          if (original_sel.select(ens.residue_selection).all_eq(True)):
            assert (not found_residue)
            for trial in self.sites_trials :
              new_atom_groups.append(new_atom_group(trial,
                ens.selection))
              n_confs += 1
            found_residue = True
        assert (found_residue)
      else :
        # adjacent residue only, so we take the best n_max out of all
        # conformations (padding with the first conf if necessary)
        all_trials = []
        ens_selection = None
        for ens in self.ensembles :
          ens_selection_ = flex.bool(n_atoms, False).set_selected(
            ens.selection, True)
          if (ens_selection_.select(ag_old_isel).all_eq(True)):
            all_trials.extend(ens.sites_trials)
            ens_selection = ens.selection
            break
        if (len(all_trials) == 0):
          raise RuntimeError("No matching atoms found for %s" %
            atom_group.id_str())
        # XXX sort by maximum deviation, or rmsd?
        all_trials.sort(lambda a,b: cmp(b.max_dev, a.max_dev))
        k = 0
        while (n_confs < n_max) and (k < len(all_trials)):
          trial = all_trials[k]
          new_atom_groups.append(new_atom_group(trial, ens_selection))
          n_confs += 1
          k += 1
      # fill in the remaining conformations with copies of the first
      while (n_confs < n_max):
        new_ag = atom_group.detached_copy()
        new_atom_groups.append(new_ag)
        n_confs += 1
      fragment_atom_groups.append(new_atom_groups)
    return fragment_atom_groups, fragment_i_seqs

def assemble_fragments(
    pdb_hierarchy,
    groups,
    primary_occupancy=0.5,
    out=None):
  if (out is None) : out = sys.stdout
  pdb_hierarchy = pdb_hierarchy.deep_copy()
  pdb_hierarchy.atoms().reset_i_seq()
  all_atom_groups_and_i_seqs = []
  for group in groups :
    group.show_summary(out=out)
    atom_groups, isels = group.combine_ensembles(pdb_hierarchy)
    all_atom_groups_and_i_seqs.extend(list(zip(atom_groups, isels)))
  n_placed = 0
  for residue_group in mmtbx.building.iter_residue_groups(pdb_hierarchy):
    if (len(residue_group.atom_groups()) > 1) : continue
    rg_isel = residue_group.atoms().extract_i_seq()
    for atom_groups, ag_isel in all_atom_groups_and_i_seqs :
      if (ag_isel.all_eq(rg_isel)):
        main_conf = residue_group.only_atom_group()
        main_conf.altloc = "A"
        for atom in main_conf.atoms():
          atom.occ = primary_occupancy
        for k, atom_group in enumerate(atom_groups):
          for atom in atom_group.atoms():
            atom.occ = (1.0 - primary_occupancy) / len(atom_groups)
          atom_group.altloc = string.uppercase[k+1]
          residue_group.append_atom_group(atom_group)
        n_placed += 1
        break
  assert (n_placed == len(all_atom_groups_and_i_seqs))
  return pdb_hierarchy
