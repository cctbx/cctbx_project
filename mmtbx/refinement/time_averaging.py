from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys, copy
from libtbx.test_utils import approx_equal
from cctbx import xray
from mmtbx.dynamics import cartesian_dynamics
from iotbx import pdb
from cctbx import adptbx
from mmtbx import utils


class run(object):

  def __init__(self, fmodels, model, target_weights, log, time_step = 0.0005,
                     n_cycles = 5, temperature = 300., n_steps = 200,
                     eq_cycles = None):
    tx = 1.
    f_calc_average = None
    if 0: fmodels.fmodel_xray().xray_structure.set_b_iso(value=1.0)
    self.xray_structures = []
    self.pdb_hierarchy = model.pdb_hierarchy
    twr = target_weights.xyz_weights_result
    for macro_cycle in xrange(n_cycles):
      print >> log, "Cycle number: ", macro_cycle
      print >> log, "start R=%6.4f" % fmodels.fmodel_xray().r_work()
      xrs_start = fmodels.fmodel_xray().xray_structure.deep_copy_scatterers()
      cd_manager = cartesian_dynamics.cartesian_dynamics(
        structure                   = fmodels.fmodel_xray().xray_structure,
        restraints_manager          = model.restraints_manager,
        temperature                 = temperature,
        n_steps                     = n_steps,
        time_step                   = time_step,
        fmodel                      = fmodels.fmodel_xray(),
        shift_update                = None,
        xray_target_weight          = twr.wx * twr.wx_scale,
        chem_target_weight          = 1,
        xray_structure_last_updated = model.xray_structure,
        log = log)
      xrs_final = cd_manager.xray_structure_last_updated
      if(eq_cycles is not None and macro_cycle > eq_cycles):
        self.xray_structures.append(xrs_final)
      result = xrs_start.distances(other = xrs_final)
      fmodels.update_xray_structure(
        xray_structure = cd_manager.xray_structure_last_updated,
        update_f_calc  = True)
      print >> log, "final R=%6.4f" % fmodels.fmodel_xray().r_work(), \
        xrs_start.distances(other = xrs_final).min_max_mean().as_tuple()
      f_calc = fmodels.fmodel_xray().f_calc()
      if(f_calc_average is None):
        f_calc_average = f_calc
      else:
        a_prime = math.exp(-cd_manager.n_steps*time_step/tx)
        f_calc_data = f_calc.data()
        f_calc_average_data = f_calc_average.data()
        f_calc_average_data = a_prime * f_calc_average_data + (1.-a_prime) * f_calc_data
        f_calc_average = f_calc_average.array(data = f_calc_average_data)
      fmodels.fmodel_xray().update(f_calc = f_calc_average)
      print >> log, "final R=%6.4f" % fmodels.fmodel_xray().r_work()
    # Modify occupancies in-place
    for i_model in xrange(len(self.xray_structures)):
      self.xray_structures[i_model].set_occupancies(
        value = (100./len(self.xray_structures))/100)
    #
    self.write_pdb(out = open("ta_multi_model.pdb", "w"))
    #
    fmodel = utils.fmodel_simple(
      xray_structures = self.xray_structures,
      f_obs = fmodels.fmodel_xray().f_obs,
      r_free_flags = fmodels.fmodel_xray().r_free_flags)
    print >> log, "FINAL Rwork=%6.4f Rfree=%6.4f" % (
      fmodel.r_work(), fmodel.r_free())

  def write_pdb(self, out, max_models = 100, num_models = 10):
    assert num_models <= 100
    crystal_symmetry = self.xray_structures[0].crystal_symmetry()
    print >> out, pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry)
    print >> out, pdb.format_scale_records(
      unit_cell = crystal_symmetry.unit_cell())
    atoms_reset_serial = True
    np = min(len(self.xray_structures)/100, 100)
    if(num_models < len(self.xray_structures)):
      num_models = len(self.xray_structures)
    for i_model, xrs in enumerate(self.xray_structures):
      if(i_model%np == 0):
        print >> out, "MODEL %8d"%(i_model+1)
        # XXX Copy from utils.py
        scatterers = xrs.scatterers()
        sites_cart = xrs.sites_cart()
        u_isos = xrs.extract_u_iso_or_u_equiv()
        occupancies = scatterers.extract_occupancies()
        u_carts = scatterers.extract_u_cart_plus_u_iso(xrs.unit_cell())
        scat_types = scatterers.extract_scattering_types()
        pdb_atoms = self.pdb_hierarchy.atoms()
        for j_seq,atom in enumerate(pdb_atoms):
          atom.xyz = sites_cart[j_seq]
          atom.occ = occupancies[j_seq]
          atom.b = adptbx.u_as_b(u_isos[j_seq])
          e = scat_types[j_seq]
          if (len(e) > 1 and "+-0123456789".find(e[1]) >= 0):
            atom.element = "%2s" % e[:1]
            atom.charge = "%-2s" % e[1:]
          elif (len(e) > 2):
            atom.element = "%2s" % e[:2]
            atom.charge = "%-2s" % e[2:]
          else:
            atom.element = "%2s" % e
            atom.charge = "  "
          if (scatterers[j_seq].flags.use_u_aniso()):
            atom.uij = u_carts[j_seq]
          else:
            atom.uij = (-1,-1,-1,-1,-1,-1)
        if (atoms_reset_serial):
          atoms_reset_serial_first_value = 1
        else:
          atoms_reset_serial_first_value = None
        out.write(self.pdb_hierarchy.as_pdb_string(
          append_end=False,
          atoms_reset_serial_first_value=atoms_reset_serial_first_value))
        #
      print >> out, "ENDMDL"
    print >> out, "END"
