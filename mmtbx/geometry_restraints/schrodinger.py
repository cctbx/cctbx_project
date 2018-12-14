from __future__ import division, print_function

import os
from copy import copy
from math import sqrt
from subprocess import Popen
from time import sleep

import numpy as np

import cctbx.geometry_restraints.manager
import cctbx.geometry_restraints.energies
import cctbx.geometry_restraints.flags
from cctbx.array_family import flex
from libtbx.utils import Sorry


master_phil_str = """
  use_schrodinger = False
    .help = Use Schrodinger force field during refinement.
    .type = bool
  maestro_file = None
    .help = Schrodinger Maestro file containing structure providing context for energies and gradients.
    .type = path
  forcefield = *default OPLS3e 2005 2005X
    .help = Force field to use during refinement.
    .type = choice
  scaling_factor = 1
    .help = Scaling factor of schrodinger energies and gradients.
    .type = float
  escale = 1
    .help = Scaling factor of the schrodinger energies.
    .type = float
  gscale = 1
    .help = Scaling factor of the schrodinger gradients.
    .type = float
  use_symmetry = False
    .help = Use crystal symmetry during energy and gradient calculation.
    .type = bool
  selection = None
    .help = Use Schrodinger force field only for selected region.
    .type = atom_selection
  maestro_asl = None
    .help = ASL for given Maestro file. The selected structure should be identical or a superset against selection.
    .type = str
  flags {
    bond = True
      .type = bool
    angle = True
      .type = bool
    dihedral = True
      .type = bool
    nonbonded = True
      .type = bool
    sgb = True
      .type = bool
    pipack = True
      .type = bool
    selfcont = True
      .type = bool
    hbond = True
      .type = bool
    lipo = True
      .type = bool
  }
  mix {
    bond = False
      .type = bool
    angle = False
      .type = bool
    dihedral = False
      .type = bool
    chirality = False
      .type = bool
    planarity = False
      .type = bool
    nonbonded = False
      .type = bool
    parallelity = False
      .type = bool
    bond_similarity = False
      .type = bool
    ramachandran = False
      .type = bool
    reference_coordinate = True
      .type = bool
    reference_dihedral = True
      .type = bool
    ncs_dihedral = True
      .type = bool
    den_restraints = True
      .type = bool
  }
  log_file = None
    .help = File name for log.
    .type = path
  debug = True
    .help = Debug manager by writing out structures that are evaluated.
    .type = bool
"""

# TODO
"""
  add_water = True
    .help = Waters not present in original Maestro file will be added.
    .type = bool
"""


class schrodinger_manager(cctbx.geometry_restraints.manager.manager):
  """Geometry restrainst manager using Schrodinger forcefield."""

  SLEEP_TIMESTEP = 0.01
  FORCEFIELD_SERVER = 'primex_phenix_server.py'
  RESIDUALS_TO_SUBTRACT = [
      "bond", "angle", "dihedral", "chirality", "planarity", "nonbonded",
      "parallelity", 'ramachandran', 'bond_similarity']
  RESIDUALS_TO_KEEP = [
      "reference_coordinate", "reference_dihedral", "ncs_dihedral",
      "den_restraints"]
  RESIDUALS_TO_MIX = [
      "bond", "angle", "dihedral", "nonbonded"]
  UNIQUE_IDS_USED = []


  def __init__(self, pdb_hierarchy, params, directory='.', use_symmetry=True,
               cleanup=False, grm=None, **kwargs):
    """
    pdb_hierarchy : iotbx.hierarchy.root
    params : phil params
    directory : Path
    use_symmetry : bool
    cleanup : bool
    grm : cctbx.geometry_restraints.manager.manager
    """

    self.grm = grm
    if grm is not None:
      # Hack to initialize manager from
      # cctbx.geometry_restraints.manager.manager instance
      for k, v in grm.__dict__.iteritems():
        setattr(self, k, v)
    else:
      super(schrodinger_manager, self).__init__(**kwargs)

    self.directory = directory
    self.params = params
    self.pdb_hierarchy = pdb_hierarchy
    # Only use symmetry when the user requested it via the command line, and
    # PHENIX has not internally rejected it as in real space refinement.
    self.use_symmetry = use_symmetry and self.params.schrodinger.use_symmetry
    self._sel_cache = self.pdb_hierarchy.atom_selection_cache()
    self._selection = self.params.schrodinger.selection

    # Define some file names to read and write later on
    # Files to write
    self._coor_fname = os.path.join(self.directory, 'COORDINATES.npy')
    self._stop_fname = os.path.join(self.directory, 'STOP')
    self._conformer_fname = os.path.join(self.directory, "CONFORMER.pdb")
    self._selection_fname = os.path.join(self.directory, "SELECTION.pdb")
    self._shift_fname = os.path.join(self.directory, "SHIFT.npy")
    self._server_id_fname = os.path.join(self.directory, 'SERVER_ID')
    self._flags_fname = os.path.join(self.directory, 'FLAGS')
    # Files to read
    self._energies_fname = os.path.join(self.directory, 'ENERGIES.npy')
    self._gradients_fname = os.path.join(self.directory, 'GRADIENTS.npy')
    self._running_fname = os.path.join(self.directory, 'RUNNING')
    self._error_fname = os.path.join(self.directory, 'ERROR')

    # Get selections belonging to each conformer
    self._conformer_selections = {}
    self._get_conformer_selections()
    # TODO Check conformer consistency. Each conformer should have the same
    # constitution.

    if cleanup:
      # Check if there are left over files from a previous run, if so remove them.
      self.stop_server(wait=1, feel_sorry=False)
      self.cleanup()

    # Start force field server
    # Give a unique ID to the server so we know when its asking for energies.
    self.server_id = len(self.UNIQUE_IDS_USED)
    self.UNIQUE_IDS_USED.append(self.server_id)

  def _get_conformer_selections(self):
    altloc_indices = self.pdb_hierarchy.altloc_indices()
    for altloc in altloc_indices:
      if altloc.strip() == '' and len(altloc_indices) > 1:
        continue
      # Selection that will be affected by force field.
      sel_str = 'all'
      if altloc:
        sel_str = "(altloc ' ' or altloc {altloc})".format(altloc=altloc)
      conf_isel = self._sel_cache.selection(sel_str).iselection()
      ff_isel = self._get_selection(altloc=altloc).iselection()
      self._conformer_selections[altloc] = (conf_isel, ff_isel)

  def update_schrodinger_flags(self, bond=True, angle=True, dihedral=True,
      nonbonded=True, sgb=True, pipack=True, hbond=True, selfcont=True,
      lipo=True, symmetry=0):
    d = {'include_bond': bond,
         'include_angl': angle,
         'include_tor': dihedral,
         'include_nb': nonbonded,
         'include_sgb': sgb,
         'include_pipack': pipack,
         'include_hb': hbond,
         'include_selfcont': selfcont,
         'include_lipo': lipo,
         'new_sym_mode': symmetry}
    with open(self._flags_fname, 'w') as f:
      for key, value in d.iteritems():
        line = '{key}\t{value}\n'.format(key=key, value=value)
        f.write(line)
    total_sleep = 0
    while os.path.exists(self._flags_fname):
      sleep(self.SLEEP_TIMESTEP)
      total_sleep += self.SLEEP_TIMESTEP
      if total_sleep > 20:
        raise Sorry("Flag file not read.")

  def energies_sites(self, 
      sites_cart, 
      flags=None,
      custom_nonbonded_function=None,
      compute_gradients=False,
      gradients=None,
      disable_asu_cache=False,
      normalization=False,
      external_energy_function=None,
      extension_objects=[],
      site_labels=None,
      ):

    # Calculate initial energies instance, so the end result has access to bond
    # deviation methods etc.
    energies = super(schrodinger_manager, self).energies_sites(
        sites_cart,
        flags=flags,
        custom_nonbonded_function=custom_nonbonded_function,
        compute_gradients=compute_gradients,
        gradients=gradients,
        disable_asu_cache=disable_asu_cache,
        normalization=False,
        external_energy_function=external_energy_function,
        extension_objects=extension_objects,
        site_labels=site_labels,
    )
    # Calculate the average gradient length, so we can scale the OPLS3e
    # ones later on.
    #if compute_gradients:
    #  g_length = sqrt(flex.mean(energies.gradients.dot()))

    # Check if there is at least one atom that should be accounted for by the
    # forcefield.
    iselection = self._get_selection().iselection()
    nselected = iselection.size()
    if nselected == 0:
      return energies

    opts = self.params.schrodinger
    # Subtract energy and gradient part that is covered by OPLS3e force field.
    # First determine which energy terms are being replaced. Terms that will be
    # exchanged are set to True in the flags instance.
    flags_subtract = copy(flags)
    if flags_subtract is None:
      flags_subtract = cctbx.geometry_restraints.flags.flags(default=True)
    for residual_name in self.RESIDUALS_TO_KEEP:
      setattr(flags_subtract, residual_name, False)
    for residual_name in self.RESIDUALS_TO_SUBTRACT:
      try:
        mix = getattr(opts.mix, residual_name)
      except AttributeError:
        continue
      if mix:
        setattr(flags_subtract, residual_name, False)

    # In case the whole structure is represented by the force field, we can
    # recalculate the terms using the original manager. Else we need to
    # subtract residuals we want to keep globally and some other terms locally.
    global_forcefield = iselection.size() == self.pdb_hierarchy.atoms_size()
    if global_forcefield:
      manager_subtract = super(schrodinger_manager, self)
      sites_cart_subtract = sites_cart
    else:
      manager_subtract = super(schrodinger_manager, self).select(iselection=iselection)
      sites_cart_subtract = sites_cart.select(iselection)

    energies_subtract = manager_subtract.energies_sites(
      sites_cart_subtract,
      flags=flags_subtract,
      custom_nonbonded_function=custom_nonbonded_function,
      compute_gradients=compute_gradients,
      gradients=None,
      disable_asu_cache=disable_asu_cache,
      normalization=False,
      external_energy_function=external_energy_function,
      extension_objects=extension_objects,
      site_labels=site_labels,
    )

    # Subtract residuals
    for residual_name in self.RESIDUALS_TO_SUBTRACT:
      residual_name_full = residual_name + "_residual_sum"
      residual = getattr(energies_subtract, residual_name_full)
      current = getattr(energies, residual_name_full)
      energy = current - residual
      setattr(energies, residual_name_full, energy)
      energies.residual_sum -= residual
    # Subtract gradients
    if compute_gradients:
      for n, i in enumerate(iselection):
        energies.gradients[i] -= np.asarray(energies_subtract.gradients[n])

    # Check if we need to start up the server. Reasons are 
    # 1) server is not running.
    # 2) symmetry requirement has changed.
    start_server = False
    if not os.path.isfile(self._running_fname):
      start_server = True
    else:
      with open(self._running_fname) as f:
        for line in f:
          if line.startswith('SYMMETRY'):
            current_symmetry = bool(int(line.split()[1]))
            if self.use_symmetry != current_symmetry:
              self.stop_server()
              start_server = True

    if start_server:
      self.start_server()
      f = opts.flags
      mix = opts.mix
      bond = f.bond and not mix.bond
      angle = f.angle and not mix.angle
      dihedral = f.dihedral and not mix.dihedral
      nonbonded = f.nonbonded  and not mix.nonbonded
      self.update_schrodinger_flags(
        bond=bond, angle=angle, dihedral=dihedral, nonbonded=nonbonded,
        sgb=f.sgb, pipack=f.pipack, hbond=f.hbond, selfcont=f.selfcont, lipo=f.lipo,
        symmetry=int(self.use_symmetry))

    # Check if this manager is properly connected to the running server, i.e.
    # that the selection / structure of this manager is represented by the
    # server.
    is_active_manager = self.is_active_manager()
    # Make it the active manager by writing id to file
    if not is_active_manager:
      with open(self._server_id_fname, 'w') as f:
        f.write("{:d}".format(self.server_id))

    sites_cart = np.asarray(sites_cart)
    nconformers = len(self._conformer_selections)
    energy_scale = opts.scaling_factor * opts.escale / nconformers
    gradients_scale = opts.scaling_factor * opts.gscale / nconformers

    for altloc, (conf_isel, ff_isel) in self._conformer_selections.iteritems():
      # Save conformer and selection pdb file for server to read and match
      if not is_active_manager:
        conformer = self.pdb_hierarchy.select(conf_isel)
        ff_part = self.pdb_hierarchy.select(ff_isel)
        crystal_symmetry = None
        ff_part.write_pdb_file(self._selection_fname)
        conformer.write_pdb_file(self._conformer_fname, crystal_symmetry=crystal_symmetry)
        total_sleep = 0
        while os.path.exists(self._conformer_fname):
          sleep(self.SLEEP_TIMESTEP)
          total_sleep += self.SLEEP_TIMESTEP
          if total_sleep > 20:
            raise Sorry("Conformer file is not being read. Shutting down.")

      # Save coordinates to file
      conformer_sites_cart = sites_cart[conf_isel]
      np.save(self._coor_fname, conformer_sites_cart)

      # Wait for the energy server to write out energies and gradients
      total_sleep = 0
      while True:
        if (os.path.isfile(self._energies_fname) and 
            os.path.isfile(self._gradients_fname)):
          # Get energy of conformer
          try:
            prime_energies = np.load(self._energies_fname) * energy_scale
          except IOError:
            sleep(self.SLEEP_TIMESTEP)
            total_sleep += self.SLEEP_TIMESTEP
            if total_sleep > 10:
              break
            continue
          if opts.debug:
            print('PRIME ENERGY ({}): {:.7f}'.format(altloc, prime_energies[0]))
          # In case no gradients were requested, gradients is None
          if compute_gradients:
            try:
              prime_gradients = np.load(self._gradients_fname)
              prime_gradients *= gradients_scale
            except IOError:
              total_sleep += self.SLEEP_TIMESTEP
              sleep(self.SLEEP_TIMESTEP)
              if total_sleep > 10:
                break
              continue

          # TODO add term specific energies from PrimeX. For now just give each
          # term a part of the energy
          energies.residual_sum += prime_energies[0]
          for residual_name in ["bond", "angle", "dihedral", "nonbonded"]:
            residual_name_full = residual_name + "_residual_sum"
            current = getattr(energies, residual_name_full)
            setattr(energies, residual_name_full, current + prime_energies[0] / 4.0)
          if compute_gradients:
            for n, i in enumerate(ff_isel):
              energies.gradients[i] += prime_gradients[n]

          # Remove files we don't need anymore
          for fname in self._energies_fname, self._gradients_fname:
            os.remove(fname)
          break
        elif os.path.isfile(self._error_fname):
            with open(self._error_fname) as f:
                error = f.read()
            self.cleanup()
            raise Sorry(error)
        else:
          if not os.path.isfile(self._running_fname):
            raise Sorry("SERVER NOT RUNNING. ABORTING.")
        sleep(self.SLEEP_TIMESTEP)

    energies.normalization = normalization
    energies.finalize_target_and_gradients()
    if opts.debug:
      print('ENERGY AFTER FINALIZING:', energies.target)
    return energies
  
  def is_active_manager(self):
    """Check if current manager is coupled to the external server"""
    is_active_manager = False
    if os.path.exists(self._server_id_fname):
      with open(self._server_id_fname) as f:
        server_id = int(f.read())
      if server_id == self.server_id:
        is_active_manager = True
    return is_active_manager

  def shift_cart(self, shift_cart):
      """Shift Cartesian coordinates of the whole structure present in the maestro_file."""
      shift_cart = np.asarray(shift_cart, dtype=np.float64)
      # Check if a previous shift file is still there, if so, wait
      totsleep = 0
      np.save(self._shift_fname, shift_cart)
      while os.path.isfile(self._shift_fname):
          sleep(self.SLEEP_TIMESTEP)
          totsleep += self.SLEEP_TIMESTEP
          if totsleep > 10:
              raise Sorry("SHIFT FILE IS NOT BEING PROCESSED.")
      sleep(self.SLEEP_TIMESTEP)

  def select(self, selection=None, iselection=None):
    manager_new = super(
        schrodinger_manager, self).select(selection=selection, iselection=iselection)

    if selection is None and iselection is None:
      pdb_hierarchy = self.pdb_hierarchy
    elif selection is None:
      pdb_hierarchy = self.pdb_hierarchy.select(iselection)
    else:
      pdb_hierarchy = self.pdb_hierarchy.select(selection)
    schrodinger_manager_new = schrodinger_manager(
        pdb_hierarchy, self.params, 
        directory=self.directory, grm=manager_new)
    return schrodinger_manager_new

  def discard_symmetry(self, new_unit_cell):
    manager_new = super(schrodinger_manager, self).discard_symmetry(new_unit_cell)
    schrodinger_manager_new = schrodinger_manager(
        self.pdb_hierarchy, self.params, directory=self.directory,
        use_symmetry=False, grm=manager_new)
    return schrodinger_manager_new

  def start_server(self):
    """Startup external PrimeServer to obtain energies and gradients."""

    opts = self.params.schrodinger
    schrodinger_dir = os.environ["SCHRODINGER"]
    run_call = os.path.join(schrodinger_dir, 'run')
    cmd = [run_call, '-FROM', 'psp', self.FORCEFIELD_SERVER, 
           opts.maestro_file]
    if opts.forcefield != 'default':
        cmd += ['-ff', opts.forcefield]
    if self.use_symmetry:
      cmd.append('-use_symmetry')
    if opts.maestro_asl is not None:
      cmd += ['-asl', self.params.schrodinger.maestro_asl]
    if opts.log_file is not None:
      cmd += ['-log', self.params.schrodinger.log_file]
    if self.params.schrodinger.debug:
      cmd += ['-debug']

    # Remove CCTBX / PHENIX from the PYTHONPATH when calling Schrodinger python
    env = os.environ.copy()
    env["PYTHONPATH"] = ""

    # Start external forcefield server
    self._process = Popen(cmd, env=env)
    while not os.path.isfile(self._running_fname):
      sleep(self.SLEEP_TIMESTEP)

  def stop_server(self, wait=10, feel_sorry=True):
    """Stop running server."""
    with open(self._stop_fname, 'w') as f:
      pass
    tot_sleep_time = 0
    while os.path.isfile(self._running_fname):
      sleep(self.SLEEP_TIMESTEP)
      tot_sleep_time += self.SLEEP_TIMESTEP
      if tot_sleep_time > wait:
        if feel_sorry:
          raise Sorry("Can't stop server")
        else:
          break
    self.cleanup()

  def is_running(self):
    if self._process.poll() is None:
      return True
    else:
      self._process.communicate()
      return False

  def cleanup(self):
    all_names = ('coor stop conformer selection shift '
                 'server_id energies gradients running flags error').split()
    for name in all_names:
      fn = getattr(self, '_' + name + '_fname')
      try:
        os.remove(fn)
      except OSError:
        pass

  def _get_selection(self, altloc=''):
    """Return selection of structure that is using the Schrodinger force field."""

    sel_str = ""
    lsel_str = self.params.schrodinger.selection
    if lsel_str is not None:
      sel_str += lsel_str

    if altloc:
      if sel_str:
        sel_str += " and "
      sel_str += "(altloc ' ' or altloc {})".format(altloc)

    if not sel_str:
      sel_str = "all"
    selection = self._sel_cache.selection(sel_str)
    return selection

  # The following methods output a new manager instance, and thus should be
  # wrapped.
  def reduce_for_tardy(self,
        tardy_tree,
        omit_bonds_with_slack_greater_than=0,
        include_den_restraints=False):
    raise NotImplementedError

  def new_including_isolated_sites(self,
        n_additional_sites,
        model_indices=None,
        conformer_indices=None,
        sym_excl_indices=None,
        donor_acceptor_excl_groups=None,
        site_symmetry_table=None,
        nonbonded_types=None,
        nonbonded_charges=None):
    raise NotImplementedError

