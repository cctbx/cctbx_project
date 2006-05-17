from cctbx.array_family import flex
import iotbx.phil
from mmtbx.refinement import minimization_individual_adp
import mmtbx.refinement.group_b
from mmtbx.tls import tools
from mmtbx.refinement import print_statistics
import scitbx.lbfgs
from cctbx import xray
from libtbx.utils import Sorry

group_adp_master_params = iotbx.phil.parse("""\
  number_of_macro_cycles   = 5
    .type = int
  max_number_of_iterations = 50
    .type = int
  one_residue_one_group    = True
    .type = bool
  convergence_test         = True
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
  number_of_macro_cycles      = 3
    .type = int
  max_number_of_iterations    = 25
    .type = int
  start_tls_value             = None
    .type = float
  run_finite_differences_test = False
    .type = bool
  eps                         = 1.e-6
    .type = float
""")

individual_adp_master_params = iotbx.phil.parse("""\
  iso {
    max_number_of_iterations = 25
      .type = int
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
            log                       = None):
    if(log is None): log = sys.stdout
    if(refine_adp_individual):
       print_statistics.make_sub_header(text= "individual ADP refinement",
                                        out = log)

       tan_u_iso = False
       param = 0
       if(tan_b_iso_max > 0.0):
          tan_u_iso = True
          param = int(tan_b_iso_max)
       #XXX Inefficient code, happens at least N-macro-cycles times
       xray.set_scatterer_grad_flags(scatterers = fmodel.xray_structure.scatterers())
       if(not refine_tls):
          for anisotropic_flag, sc in zip(anisotropic_flags,
                                          fmodel.xray_structure.scatterers()):
              if(anisotropic_flag):
                 sc.flags.set_use_u_aniso(True)
                 sc.flags.set_grad_u_aniso(True)
                 sc.flags.set_use_u_iso(False)
                 sc.flags.set_grad_u_iso(False)
              else:
                 sc.flags.set_use_u_aniso(False)
                 sc.flags.set_grad_u_aniso(False)
                 sc.flags.set_use_u_iso(True)
                 sc.flags.set_grad_u_iso(True)
                 sc.flags.set_tan_u_iso(tan_u_iso)
                 sc.flags.param = param
       else:
          #print "OK1"
          for tls_selection in tls_selections:
              for anisotropic_flag, sc, tls_si in zip(anisotropic_flags,
                                              fmodel.xray_structure.scatterers(), tls_selection):
                  if(tls_si and sc.u_star != (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0)):
                     sc.flags.set_use_u_aniso(True)
                     sc.flags.set_use_u_iso(True)
                     if(not anisotropic_flag):
                       sc.flags.set_grad_u_aniso(False)
                       sc.flags.set_grad_u_iso(True)
                     else:
                       sc.flags.set_grad_u_aniso(True)
                       sc.flags.set_grad_u_iso(False)
                     #sc.flags.set_tan_u_iso(tan_u_iso)
                     #sc.flags.param = param
                  else:
                     sc.flags.set_use_u_aniso(False)
                     sc.flags.set_use_u_iso(True)
                     if(not anisotropic_flag):
                       sc.flags.set_grad_u_aniso(False)
                       sc.flags.set_grad_u_iso(True)
                     else:
                       sc.flags.set_grad_u_aniso(True)
                       sc.flags.set_grad_u_iso(False)
                     sc.flags.set_tan_u_iso(False)
                     sc.flags.param = 0
                  if(sc.u_iso <= 0): sc.flags.set_use_u_iso(False)
       ##############
       c1=0;c2=0;c3=0;c4=0;c5=0
       for sc in fmodel.xray_structure.scatterers():
         c1 += sc.flags.use_u_iso()
         c2 += sc.flags.use_u_aniso()
         c3 += sc.flags.grad_u_iso()
         c4 += sc.flags.grad_u_aniso()
         c5 += sc.flags.tan_u_iso()
       #print "use_u_iso()    = ", c1
       #print "use_u_aniso()  = ", c2
       #print "grad_u_iso()   = ", c3
       #print "grad_u_aniso() = ", c4
       #print "tan_u_iso()    = ", c5
       ##############
       lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
           max_iterations = individual_adp_params.iso.max_number_of_iterations)
       self.minimized = minimization_individual_adp.lbfgs(
            restraints_manager       = restraints_manager,
            fmodel                   = fmodel,
            lbfgs_termination_params = lbfgs_termination_params,
            wx                       = wx,
            wu                       = wu_individual,
            wilson_b                 = wilson_b,
            tan_b_iso_max            = tan_b_iso_max,
            iso_restraints           = adp_restraints_params.iso,
            verbose                  = 0)
       self.minimized.show(text = "LBFGS minimization", out  = log)
       fmodel.update_xray_structure(
                               xray_structure           = self.minimized.xray_structure,
                               update_f_calc            = True,
                               out                      = log)
    if(refine_adp_group):
       print_statistics.make_sub_header(text= "group isotropic ADP refinement",
                                        out = log)
       group_b_manager = mmtbx.refinement.group_b.manager(
          fmodel                   = fmodel,
          selections               = group_adp_selections,
          convergence_test         = group_adp_params.convergence_test,
          max_number_of_iterations = group_adp_params.max_number_of_iterations,
          number_of_macro_cycles   = group_adp_params.number_of_macro_cycles,
          log                      = log,
          tan_b_iso_max            = tan_b_iso_max)
    if(refine_tls):
       print_statistics.make_sub_header(text = "TLS refinement",
                                        out  = log)
       tls_refinement_manager = tools.tls_refinement(
          fmodel                      = fmodel,
          selections                  = tls_selections,
          refine_T                    = tls_params.refine_T,
          refine_L                    = tls_params.refine_L,
          refine_S                    = tls_params.refine_S,
          number_of_macro_cycles      = tls_params.number_of_macro_cycles,
          max_number_of_iterations    = tls_params.max_number_of_iterations,
          start_tls_value             = tls_params.start_tls_value,
          run_finite_differences_test = tls_params.run_finite_differences_test,
          eps                         = tls_params.eps)
