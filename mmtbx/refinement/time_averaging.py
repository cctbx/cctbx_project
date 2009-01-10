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

  def __init__(self, fmodels,
                     model,
                     target_weights,
                     log,
                     time_step = 0.0005,
                     number_of_macro_cycles = 11,
                     temperature = 300.,
                     n_steps = 10,
                     eq_cycles = 0,
                     b_target = 1.0,
                     b_target_cycles = 10,
                     tx = 2.0):
    fmodel = fmodels.fmodel_xray()
    assert fmodel.xray_structure is model.xray_structure
    xray_structure_last_updated = model.xray_structure.deep_copy_scatterers()
    xrs_start = fmodels.fmodel_xray().xray_structure.deep_copy_scatterers()
    xray_gradient = None
    reset_velocities = True
    f_calc_average = None
    wx = target_weights.xyz_weights_result.wx * \
      target_weights.xyz_weights_result.wx_scale
    #
    if 1:
      b_isos = fmodel.xray_structure.scatterers().extract_u_iso()/adptbx.b_as_u(1)
      assert b_target_cycles != 0
      b_decs = (b_isos - b_target)/b_target_cycles
    #
    self.xray_structures = []
    self.pdb_hierarchy = model.pdb_hierarchy
    for macro_cycle in xrange(number_of_macro_cycles):
      print >> log
      if(macro_cycle == 0):
        print >> log, "mc: %d r_work=%6.4f r_free=%6.4f"%(
          macro_cycle, fmodel.r_work(), fmodel.r_free())
      assert fmodel.xray_structure is model.xray_structure
      cd_manager = cartesian_dynamics.cartesian_dynamics(
        structure                   = model.xray_structure,
        restraints_manager          = model.restraints_manager,
        temperature                 = temperature,
        n_steps                     = n_steps,
        time_step                   = time_step,
        initial_velocities_zero_fraction = 0,
        fmodel                      = fmodel,
        xray_target_weight          = wx,
        chem_target_weight          = target_weights.xyz_weights_result.w,
        xray_structure_last_updated = xray_structure_last_updated,
        shift_update                = 0.0,
        xray_gradient               = xray_gradient,
        reset_velocities            = reset_velocities,
        update_f_calc               = False,
        log                         = log)
      reset_velocities = False
      xray_structure_last_updated = \
        cd_manager.xray_structure_last_updated.deep_copy_scatterers()
      result = xrs_start.distances(other =
        model.xray_structure).min_max_mean().as_tuple()
      xray_gradient = cd_manager.xray_gradient
      assert model.xray_structure is cd_manager.structure
      fmodel.update_xray_structure(
        xray_structure = model.xray_structure,
        update_f_calc  = True,
        update_f_mask  = True)
      if(macro_cycle > 0):
        print >> log, "mc: %d r_work=%6.4f r_free=%6.4f deviation (max,mean)=%5.3f %5.3f"%(
          macro_cycle, fmodel.r_work(), fmodel.r_free(),result[1], result[2])
      #
      if(eq_cycles is not None):
        if(macro_cycle > eq_cycles):
          self.xray_structures.append(model.xray_structure.deep_copy_scatterers())
      else:
        self.xray_structures.append(model.xray_structure.deep_copy_scatterers())
      #
      f_calc = fmodels.fmodel_xray().f_calc()
      if(f_calc_average is None):
        f_calc_average = f_calc
      else:
        a_prime = math.exp(-cd_manager.n_steps*time_step/tx)
        f_calc_data = f_calc.data()
        f_calc_average_data = f_calc_average.data()
        f_calc_average_data = a_prime * f_calc_average_data + (1.-a_prime) * f_calc_data
        f_calc_average = f_calc_average.array(data = f_calc_average_data)
      fmodel.update(f_calc = f_calc_average)
      #
      if 1:
        b_isos = fmodel.xray_structure.scatterers().extract_u_iso()/adptbx.b_as_u(1)
        b_isos = b_isos - b_decs
        sel = b_isos <= b_target
        b_isos = b_isos.set_selected(sel, b_target)
        fmodels.fmodel_xray().xray_structure.set_b_iso(values = b_isos)
      #
      if(macro_cycle > 0):
        print >> log, "mc: %d r_work=%6.4f r_free=%6.4f deviation (max,mean)=%5.3f %5.3f"%(
          macro_cycle, fmodel.r_work(), fmodel.r_free(),result[1], result[2])
    # Modify occupancies in-place or reset to zero
    for i_model in xrange(len(self.xray_structures)):
      self.xray_structures[i_model].set_occupancies(
        value = (100./len(self.xray_structures))/100)
    #
    self.write_pdb(out = open("ta_multi_model.pdb", "w"), log = log)
    #
    if 1:
      fmodel_ = utils.fmodel_simple(
        xray_structures = self.xray_structures,
        f_obs = fmodel.f_obs,
        r_free_flags = fmodel.r_free_flags)
      print >> log, "FINAL Rwork=%6.4f Rfree=%6.4f" % (
        fmodel_.r_work(), fmodel_.r_free())


  def write_pdb(self, log, out, max_models = 100, num_models = 10):
    assert num_models <= 100
    crystal_symmetry = self.xray_structures[0].crystal_symmetry()
    print >> out, pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry)
    print >> out, pdb.format_scale_records(
      unit_cell = crystal_symmetry.unit_cell())
    atoms_reset_serial = True
    #
    if(num_models > len(self.xray_structures)):
      num_models = len(self.xray_structures)
    np = min(len(self.xray_structures)/num_models, max_models)
    cntr = 0
    for i_model, xrs in enumerate(self.xray_structures):
      if(i_model%np == 0):
        cntr += 1
        print >> out, "MODEL %8d"%cntr
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
