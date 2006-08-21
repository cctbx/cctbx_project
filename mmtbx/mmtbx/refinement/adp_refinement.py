from cctbx.array_family import flex
import iotbx.phil
from mmtbx.refinement import minimization_individual_adp
import mmtbx.refinement.group_b
from mmtbx.tls import tools
from mmtbx.refinement import print_statistics
import scitbx.lbfgs
from cctbx import xray
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry

from scitbx.python_utils.misc import user_plus_sys_time

time_adp_refinement_py = 0.0
time_set_flags         = 0.0

def show_times(out = None):
  if(out is None): out = sys.stdout
  total = time_adp_refinement_py + \
          time_set_flags
  if(total > 0.01):
     print >> out, "ADP refinement:"
     print >> out, "  time spent in adp_refinement.py          = %-7.2f" % time_adp_refinement_py
     print >> out, "  set_flags                                = %-7.2f" % time_set_flags
  return total

group_adp_master_params = iotbx.phil.parse("""\
  number_of_macro_cycles   = 3
    .type = int
  max_number_of_iterations = 25
    .type = int
  one_residue_one_group    = True
    .type = bool
  convergence_test         = False
    .type = bool
  selection                = None
    .type=str
    .multiple=True
""")

tls_master_params = iotbx.phil.parse("""\
  selection                   = None
    .type=str
    .multiple=True
  one_residue_one_group       = None
    .type = bool
  refine_T                    = True
    .type = bool
  refine_L                    = True
    .type = bool
  refine_S                    = True
    .type = bool
  number_of_macro_cycles      = 2
    .type = int
  max_number_of_iterations    = 25
    .type = int
  start_tls_value             = None
    .type = float
  run_finite_differences_test = False
    .type = bool
  eps                         = 1.e-6
    .type = float
  write_out_as {
    tls_and_adp_local = True
      .type = bool
    tls_zero_and_adp_total = False
      .type = bool
  }
""")

individual_adp_master_params = iotbx.phil.parse("""\
  convert_to_isotropic = True
    .type = bool
  iso {
    max_number_of_iterations = 25
      .type = int
    automatic_randomization_if_all_equal = True
      .type = bool
    shake_biso               = None
      .type = float
    set_biso                 = None
      .type = float
    set_biso_random          = False
      .type = bool
    set_biso_to_wilson_b     = False
      .type = bool
    scaling {
      scale_max       = 3.0
        .type = float
      scale_min       = 10.0
        .type = float
    }
  }
  aniso {
    selection = None
      .type=str
  }
""")

adp_restraints_master_params = iotbx.phil.parse("""\
  iso {
    use_u_local_only = False
      .type = bool
    sphere_radius = 5.0
      .type = float
    distance_power = 1.69
      .type = float
    average_power = 1.03
      .type = float
    wilson_b_weight_auto = False
      .type = bool
    wilson_b_weight = None
      .type = float
    plain_pairs_radius = 5.0
      .type = float
    refine_ap_and_dp = False
      .type = bool
    b_iso_max = None
      .type = float
  }
""")

class manager(object):
  def __init__(
            self,
            fmodel,
            model,
            anisotropic_flags,
            group_adp_selections      = None,
            group_adp_params          = group_adp_master_params.extract(),
            tls_selections            = None,
            tls_params                = tls_master_params.extract(),
            individual_adp_params     = individual_adp_master_params.extract(),
            adp_restraints_params     = adp_restraints_master_params.extract(),
            wilson_b                  = None,
            refine_adp_individual     = None,
            refine_adp_group          = None,
            refine_tls                = None,
            tan_b_iso_max             = None,
            restraints_manager        = None,
            wu_individual             = None,
            wx                        = None,
            macro_cycle               = None,
            log                       = None,
            fmodel_neutron            = None,
            wx_individual_neutron     = None,
            neutron_scattering_dict   = None,
            xray_scattering_dict      = None,
            wxnu_scale                = None):
    global time_adp_refinement_py
    timer = user_plus_sys_time()
    if(log is None): log = sys.stdout
    tan_u_iso = False
    param = 0
    if(tan_b_iso_max > 0.0):
       tan_u_iso = True
       param = int(tan_b_iso_max)
    if(macro_cycle == 1):
       offset = True
    else:
       offset = False

    if(refine_tls):
       print_statistics.make_sub_header(text = "TLS refinement",
                                        out  = log)
       if(refine_adp_individual or refine_adp_group):
          if(macro_cycle == 1):
             gbr_selections = []
             for s in tls_selections:
                 gbr_selections.append(s.iselection())
          else:
             gbr_selections = group_adp_selections
          group_b_manager = mmtbx.refinement.group_b.manager(
             fmodel                   = fmodel,
             selections               = gbr_selections,
             convergence_test         = group_adp_params.convergence_test,
             max_number_of_iterations = 50,
             number_of_macro_cycles   = 1,
             log                      = log)

       set_flags(xray_structure        = fmodel.xray_structure,
                 anisotropic_flags     = anisotropic_flags,
                 refine_adp_individual = refine_adp_individual,
                 refine_adp_group      = refine_adp_group,
                 tls_selections        = tls_selections,
                 refine_tls            = refine_tls,
                 tan_u_iso             = tan_u_iso,
                 param                 = param,
                 offset                = offset)
       model.show_groups(tls = True, out = log)
       current_target_name = fmodel.target_name
       fmodel.update(target_name = "ls_wunit_k1")
       self.tls_refinement_manager = tools.tls_refinement(
          fmodel                      = fmodel,
          model                       = model,
          selections                  = tls_selections,
          refine_T                    = tls_params.refine_T,
          refine_L                    = tls_params.refine_L,
          refine_S                    = tls_params.refine_S,
          number_of_macro_cycles      = tls_params.number_of_macro_cycles,
          max_number_of_iterations    = tls_params.max_number_of_iterations,
          start_tls_value             = tls_params.start_tls_value,
          run_finite_differences_test = tls_params.run_finite_differences_test,
          eps                         = tls_params.eps,
          out                         = log,
          macro_cycle = macro_cycle)
       fmodel.update(target_name = current_target_name)
       fmodel.update_xray_structure(
                   xray_structure = self.tls_refinement_manager.fmodel.xray_structure,
                   update_f_calc  = True,
                   out            = log)
       model.xray_structure = fmodel.xray_structure

    if(refine_adp_individual):
       print_statistics.make_sub_header(text= "Individual ADP refinement",
                                        out = log)
       set_flags(xray_structure     = fmodel.xray_structure,
              anisotropic_flags     = anisotropic_flags,
              refine_adp_individual = refine_adp_individual,
              refine_adp_group      = refine_adp_group,
              tls_selections        = tls_selections,
              refine_tls            = refine_tls,
              tan_u_iso             = tan_u_iso,
              param                 = param,
              offset                = offset)
       lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
           max_iterations = individual_adp_params.iso.max_number_of_iterations)
       fmodel.xray_structure.approx_equal(other = model.xray_structure)
       self.minimized = minimization_individual_adp.lbfgs(
               restraints_manager       = restraints_manager,
               fmodel                   = fmodel,
               model                    = model,
               lbfgs_termination_params = lbfgs_termination_params,
               wx                       = wx,
               wu                       = wu_individual,
               wilson_b                 = wilson_b,
               tan_b_iso_max            = tan_b_iso_max,
               iso_restraints           = adp_restraints_params.iso,
               verbose                  = 0,
               fmodel_neutron           = fmodel_neutron,
               wx_individual_neutron    = wx_individual_neutron,
               neutron_scattering_dict  = neutron_scattering_dict,
               xray_scattering_dict     = xray_scattering_dict,
               wxnu_scale               = wxnu_scale)
       self.minimized.show(text = "LBFGS minimization", out  = log)
       fmodel.update_xray_structure(
                                xray_structure = self.minimized.xray_structure,
                                update_f_calc  = True,
                                out            = log)
       model.xray_structure = fmodel.xray_structure
       fmodel.xray_structure.approx_equal(other = model.xray_structure)
    if(refine_adp_group):
       print_statistics.make_sub_header(text= "group isotropic ADP refinement",
                                        out = log)
       set_flags(xray_structure        = fmodel.xray_structure,
              anisotropic_flags     = anisotropic_flags,
              refine_adp_individual = refine_adp_individual,
              refine_adp_group      = refine_adp_group,
              tls_selections        = tls_selections,
              refine_tls            = refine_tls,
              tan_u_iso             = tan_u_iso,
              param                 = param,
              offset                = offset)
       group_b_manager = mmtbx.refinement.group_b.manager(
          fmodel                   = fmodel,
          selections               = group_adp_selections,
          convergence_test         = group_adp_params.convergence_test,
          max_number_of_iterations = group_adp_params.max_number_of_iterations,
          number_of_macro_cycles   = group_adp_params.number_of_macro_cycles,
          log                      = log)
    time_adp_refinement_py += timer.elapsed()


def set_flags(xray_structure,
              anisotropic_flags,
              refine_adp_individual,
              refine_adp_group,
              tls_selections,
              refine_tls,
              tan_u_iso,
              param,
              offset):
  global time_set_flags
  timer = user_plus_sys_time()
  #XXX n_tls_groups Python loops over scatterers: highly inefficient.
  xrs = xray_structure
  ueq1 = xray_structure.extract_u_iso_or_u_equiv()
  xray.set_scatterer_grad_flags(scatterers = xrs.scatterers()) #set all false
  if(refine_adp_individual and refine_tls):
     if(anisotropic_flags.count(True) > 0):
        raise Sorry(
           "Simultaneous TLS and individual anisotropic ADP refinement "\
           "is not implemented.")
     u_iso = xrs.scatterers().extract_u_iso()
     if((u_iso < 0.0).count(True) > 0):
        raise Sorry("All u_iso must be > 0 at this step.")
     tools.split_u(xray_structure, tls_selections, offset)
     for sc in xrs.scatterers():
         sc.flags.set_use_u_iso(True)
         sc.flags.set_grad_u_iso(True)
         sc.flags.set_tan_u_iso(tan_u_iso)
         sc.flags.param = param
         assert sc.u_iso != -1.0
         if(sc.u_star != (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)):
            sc.flags.set_use_u_aniso(True)
  elif(refine_tls):
     for tls_selection in tls_selections:
         xray_structure.convert_to_anisotropic(selection = tls_selection)
  elif(refine_adp_individual):
     for anisotropic_flag, sc in zip(anisotropic_flags, xrs.scatterers()):
         if(anisotropic_flag):
            sc.flags.set_use_u_aniso(True)
            sc.flags.set_grad_u_aniso(True)
            sc.flags.set_use_u_iso(False)
            sc.flags.set_grad_u_iso(False)
            assert sc.u_star != (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)
         else:
            sc.flags.set_use_u_aniso(False)
            sc.flags.set_grad_u_aniso(False)
            sc.flags.set_use_u_iso(True)
            sc.flags.set_grad_u_iso(True)
            sc.flags.set_tan_u_iso(tan_u_iso)
            sc.flags.param = param
            assert sc.u_iso != -1.0
  elif(refine_adp_group):
     #XXX Do this convertion only for gbr selected atoms.
     xray_structure.convert_to_isotropic()
  else:
     raise Sorry("None of refine_adp_individual, refine_adp_group, "\
                 "refine_tls is defined.")
  u_iso = xray_structure.extract_u_iso_or_u_equiv()
  assert (u_iso < 0.0).count(True) == 0
  ueq2 = xray_structure.extract_u_iso_or_u_equiv()
  assert approx_equal(ueq1, ueq2)
  time_set_flags += timer.elapsed()
