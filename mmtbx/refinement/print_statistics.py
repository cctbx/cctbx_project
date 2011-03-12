from cctbx.array_family import flex
import scitbx.math.euler_angles
from scitbx import matrix
from libtbx.utils import format_cpu_times, getenv_bool
from libtbx import adopt_init_args
import sys, time
from libtbx import str_utils
from libtbx.str_utils import prefix_each_line_suffix
from libtbx import introspection
from stdlib import math
from cctbx import xray
import cctbx.xray.structure_factors.global_counters
from mmtbx import max_lik
from itertools import count


enable_show_process_info = getenv_bool(
  "MMTBX_PRINT_STATISTICS_ENABLE_SHOW_PROCESS_INFO")

time_collect_and_process = 0.0

def show_times(out = None):
  if(out is None): out = sys.stdout
  total = time_collect_and_process
  if(total > 0.01):
     print >> out, "Collect and process                      = %-7.2f" % time_collect_and_process
  return total

def show_process_info(out):
  print >> out, "\\/"*39
  introspection.virtual_memory_info().show_if_available(out=out, show_max=True)
  xray.structure_factors.global_counters.show(out=out)
  print >> out, format_cpu_times()
  print >> out, "/\\"*39
  out.flush()

def make_header(line, out=None):
  if (out is None): out = sys.stdout
  if (enable_show_process_info):
    show_process_info(out=out)
  str_utils.make_header(line, out=out, header_len=80)

def make_sub_header(text, out=None):
  if (out is None): out = sys.stdout
  str_utils.make_sub_header(text, out=out, header_len=80)

def macro_cycle_header(macro_cycle, number_of_macro_cycles, out=None):
  if (out is None): out = sys.stdout
  #show_process_info(out=out)
  header_len = 80
  macro_cycle = str(macro_cycle)
  number_of_macro_cycles = str(number_of_macro_cycles)
  macro_cycle_str = len(macro_cycle)
  number_of_macro_cycles_str = len(number_of_macro_cycles)
  line_len = len(" REFINEMENT MACRO_CYCLE "+macro_cycle+" OF "+\
             number_of_macro_cycles)+1
  fill_len = header_len - line_len
  fill_rl = fill_len/2
  fill_r = fill_rl
  fill_l = fill_rl
  if (fill_rl*2 != fill_len): fill_r +=1
  str1 = "\n"+"*"*(fill_l-1)+" REFINEMENT MACRO_CYCLE "+macro_cycle+" OF "
  str2 = number_of_macro_cycles+" "+"*"*(fill_r)+"\n"
  out_string = str1+str2
  print >> out, out_string
  out.flush()

def show_histogram(data,
                   n_slots,
                   out=None,
                   prefix=""):
    if (out is None): out = sys.stdout
    print >> out, prefix
    histogram = flex.histogram(data    = data,
                               n_slots = n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print >> out, "%7.3f - %7.3f: %d" % (low_cutoff, high_cutoff, n)
      low_cutoff = high_cutoff
    out.flush()
    return histogram

def show_4_histograms(h_1, h_2, h_3, h_4, n_slots, out):
  format = "|%5.3f-%5.3f: %4d|%5.3f-%5.3f: %4d|%6.2f-%6.2f: %4d|%6.2f-%6.2f: %4d  |"
  lc_1 = h_1.data_min()
  lc_2 = h_2.data_min()
  lc_3 = h_3.data_min()
  lc_4 = h_4.data_min()
  s_1 = enumerate(h_1.slots())
  s_2 = enumerate(h_2.slots())
  s_3 = enumerate(h_3.slots())
  s_4 = enumerate(h_4.slots())
  for (i_1,n_1),(i_2,n_2),(i_3,n_3),(i_4,n_4) in zip(s_1,s_2,s_3,s_4):
    hc_1 = h_1.data_min() + h_1.slot_width() * (i_1+1)
    hc_2 = h_2.data_min() + h_2.slot_width() * (i_2+1)
    hc_3 = h_3.data_min() + h_3.slot_width() * (i_3+1)
    hc_4 = h_4.data_min() + h_4.slot_width() * (i_4+1)
    outs = (lc_1,hc_1,n_1,lc_2,hc_2,n_2,lc_3,hc_3,n_3,lc_4,hc_4,n_4)
    print >> out, format % outs
    lc_1 = hc_1
    lc_2 = hc_2
    lc_3 = hc_3
    lc_4 = hc_4
  out.flush()


def write_pdb_header(out,
                     xray_structure,
                     monitor):
   remark = "REMARK "
   print >> out, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f" % \
       xray_structure.unit_cell().parameters(), str(xray_structure.space_group_info())
   monitor.print_items(f = out, remark = remark)
   out.flush()

def show_rigid_body_rotations_and_translations(
      out,
      prefix,
      frame,
      euler_angle_convention,
      rotations,
      translations):
  assert euler_angle_convention in ["xyz", "zyz"]
  euler_angles_as_matrix = getattr(
    scitbx.math.euler_angles, euler_angle_convention+"_matrix")
  print >> out, prefix_each_line_suffix(
    prefix=prefix+frame, lines_as_one_string=
"                            rotation (deg)                 translation (A)   "
"\n"
"                         %s              total           xyz          total "
      % euler_angle_convention, suffix=frame)
  for i,r,t in zip(count(1), rotations, translations):
    r = list(r)
    r.reverse()
    r_total = abs(scitbx.math.r3_rotation_axis_and_angle_from_matrix(
      r=euler_angles_as_matrix(*r)).angle(deg=True))
    t_total = abs(matrix.col(t))
    print >> out, (prefix + frame +
      " group %4d: %8.3f %8.3f %8.3f %7.2f  %6.2f %6.2f %6.2f %6.2f "
        % tuple([i] + r + [r_total] + list(t) + [t_total])
      + frame).rstrip()
  out.flush()

class refinement_monitor(object):
  def __init__(self, params,
                     model_ref = None,
                     out=None,
                     short=False,
                     neutron_refinement = None,
                     call_back_handler=None,
                     is_neutron_monitor=False):
    adopt_init_args(self, locals())
    if (self.out is None): self.out = sys.stdout
    self.model_ini = None
    self.steps = []
    self.wilson_b = None
    self.bond_start = None
    self.angle_start= None
    self.bond_final = None
    self.angle_final= None
    self.r_works         = []
    self.r_frees         = []
    self.targets_w       = []
    self.targets_t       = []
    self.k_sols          = []
    self.b_sols          = []
    self.b_anisos        = []
    self.wxcs            = []
    self.wxus            = []
    self.wxc_scales      = []
    self.wxu_scales      = []
    self.angles_xc       = []
    self.angles_xu       = []
    self.alphas          = []
    self.betas           = []
    self.foms            = []
    self.phers           = []
    self.as_ave          = []
    self.as_max          = []
    self.bs_ave          = []
    self.bs_max          = []
    self.cs_ave          = []
    self.cs_max          = []
    self.ds_ave          = []
    self.ds_max          = []
    self.ps_ave          = []
    self.ps_max          = []
    self.rs_ave          = []
    self.rs_min          = []
    self.targets_c       = []
    self.gs_c_norm       = []
    self.wcs             = []
    self.bs_iso_max_a    = []
    self.bs_iso_min_a    = []
    self.bs_iso_ave_a    = []
    self.bs_iso_max_p    = []
    self.bs_iso_min_p    = []
    self.bs_iso_ave_p    = []
    self.bs_iso_max_s    = []
    self.bs_iso_min_s    = []
    self.bs_iso_ave_s    = []
    self.tus             = []
    self.wus             = []
    self.ds_ic_min       = []
    self.ds_ic_max       = []
    self.ds_ic_ave       = []
    self.ds_rc_min       = []
    self.ds_rc_max       = []
    self.ds_rc_ave       = []
    self.dopts_ic        = []
    self.dopts_rc        = []
    self.n_zeros         = []
    self.n_zeros_p       = []
    self.k1s_w           = []
    self.k1s_t           = []
    self.k3s_w           = []
    self.k3s_t           = []
    self.scale_ml        = []
    self.n_solv          = []
    self.rigid_body_shift_accumulator             = None

  def collect(self, model,
                    fmodel,
                    step,
                    tan_b_iso_max,
                    wilson_b = None,
                    target_weights = None,
                    rigid_body_shift_accumulator = None):
    global time_collect_and_process
    t1 = time.time()
    if(self.model_ini is None): self.model_ini = model.deep_copy()
    if(wilson_b is not None): self.wilson_b = wilson_b
    self.steps.append(step)
    ###
    self.r_works         .append(fmodel.r_work()                  )
    self.r_frees         .append(fmodel.r_free()                  )
    t_r = fmodel.target_functor()(compute_gradients=False)
    self.targets_w       .append(t_r.target_work())
    self.targets_t       .append(t_r.target_test())
    self.k_sols          .append(fmodel.k_sol_b_sol()[0]          )
    self.b_sols          .append(fmodel.k_sol_b_sol()[1]          )
    self.b_anisos        .append(fmodel.b_cart()                  )
    if(target_weights is not None):
       self.wxcs            .append(target_weights.wx_xyz()          )
       self.wxus            .append(target_weights.wx_adp()          )
       self.wxc_scales      .append(target_weights.holder.wxc_scale  )
       self.wxu_scales      .append(target_weights.holder.wxu_scale  )
       angle_xc = target_weights.holder.angle_xc
       angle_xu = target_weights.holder.angle_xu
       if(angle_xc is not None): angle_xc = angle_xc * 180 / math.pi
       else:                     angle_xc = 0.0
       if(angle_xu is not None): angle_xu = angle_xu * 180 / math.pi
       else:                     angle_xu = 0.0
       self.angles_xc       .append(angle_xc)
       self.angles_xu       .append(angle_xu)
    alpha, beta = fmodel.alpha_beta()
    self.alphas          .append(flex.mean_default(alpha.data(),0)          )
    self.betas           .append(flex.mean_default(beta.data(),0)           )
    self.foms            .append(flex.mean_default(fmodel.figures_of_merit(),0))
    self.phers           .append(flex.mean_default(fmodel.phase_errors(),0) )
    geom = model.geometry_statistics(ignore_hd = True,
                                     ignore_side_chain = True,
                                     )
    self.main_chain_geometry_statistics = geom
    if False:
      print "Main chain stats"
      geom.show()
    geom = model.geometry_statistics(ignore_hd = not self.neutron_refinement)
    if(geom is not None):
      self.as_ave          .append(geom.a_mean                      )
      self.as_max          .append(geom.a_max                       )
      self.bs_ave          .append(geom.b_mean                      )
      self.bs_max          .append(geom.b_max                       )
      self.cs_ave          .append(geom.c_mean                      )
      self.cs_max          .append(geom.c_max                       )
      self.ds_ave          .append(geom.d_mean                      )
      self.ds_max          .append(geom.d_max                       )
      self.ps_ave          .append(geom.p_mean                      )
      self.ps_max          .append(geom.p_max                       )
      self.rs_ave          .append(geom.n_mean                      )
      self.rs_min          .append(geom.n_min                       )
      self.targets_c       .append(geom.target                      )
      self.gs_c_norm       .append(geom.norm_of_gradients           )
    if(target_weights is not None):
       self.wcs             .append(target_weights.wc()              )
    #
    hd_sel = None
    if(not self.neutron_refinement and not self.is_neutron_monitor):
      hd_sel = model.xray_structure.hd_selection()
    b_isos = model.xray_structure.extract_u_iso_or_u_equiv() * math.pi**2*8
    if(hd_sel is not None): b_isos = b_isos.select(~hd_sel)
    self.bs_iso_max_a    .append(flex.max_default( b_isos, 0)     )
    self.bs_iso_min_a    .append(flex.min_default( b_isos, 0)     )
    self.bs_iso_ave_a    .append(flex.mean_default(b_isos, 0)     )
    s_sel = model.solvent_selection()
    if(hd_sel is not None): s_sel = s_sel.select(~hd_sel)
    if(s_sel.count(True) > 0):
       b_isos_s = b_isos.select(s_sel)
       b_isos_p = b_isos.select(~s_sel)
       self.bs_iso_max_p .append(flex.max_default( b_isos_p, 0)   )
       self.bs_iso_min_p .append(flex.min_default( b_isos_p, 0)   )
       self.bs_iso_ave_p .append(flex.mean_default(b_isos_p, 0)   )
       self.bs_iso_max_s .append(flex.max_default( b_isos_s, 0)   )
       self.bs_iso_min_s .append(flex.min_default( b_isos_s, 0)   )
       self.bs_iso_ave_s .append(flex.mean_default(b_isos_s, 0)   )
    adp = model.adp_statistics()
    #self.tus             .append(adp.target_adp_iso               )
    if(target_weights is not None):
       self.wus             .append(target_weights.wu()              )
    rem = relative_errors(
              model                  = model,
              model_ini              = self.model_ini,
              model_ref              = self.model_ref,
              compute_optimal_errors = False)
    self.ds_ic_min       .append(rem.d_ic_min                     )
    self.ds_ic_max       .append(rem.d_ic_max                     )
    self.ds_ic_ave       .append(rem.d_ic_ave                     )
    if(self.model_ref is not None):
       self.ds_rc_min    .append(rem.d_rc_min                     )
       self.ds_rc_max    .append(rem.d_rc_max                     )
       self.ds_rc_ave    .append(rem.d_rc_ave                     )
    if(fmodel.target_attributes().pseudo_ml):
       n_zero = fmodel.f_star_w_star_obj().number_of_f_star_zero()
       percent = n_zero*100./fmodel.f_obs.data().size()
       self.n_zeros      .append(n_zero                           )
       self.n_zeros_p    .append(percent                          )
    self.k1s_w           .append(fmodel.scale_k1_w()              )
    self.k1s_t           .append(fmodel.scale_k1_t()              )
    self.k3s_w           .append(fmodel.scale_k3_w()              )
    self.k3s_t           .append(fmodel.scale_k3_t()              )


    scale_ml = None
    if (fmodel.target_name=="twin_lsq_f"):
      scale_ml = 1
    else:
      if(fmodel.alpha_beta_params.method == "calc"):
         if(fmodel.alpha_beta_params.fix_scale_for_calc_option == None):
            scale_ml = fmodel.scale_ml()
         else:
            scale_ml =fmodel.alpha_beta_params.fix_scale_for_calc_option
      if(fmodel.alpha_beta_params.method == "est"):
         scale_ml = 1.0

    self.scale_ml        .append(scale_ml                         )
    self.n_solv          .append(model.number_of_ordered_solvent_molecules())
    if([self.bond_start,self.angle_start].count(None) == 2):
       if(len(self.bs_ave)>0):
         self.bond_start  = self.bs_ave[0]
         self.angle_start = self.as_ave[0]
    try:
      self.bond_final  = self.bs_ave[len(self.bs_ave)-1]
      self.angle_final = self.as_ave[len(self.as_ave)-1]
    except:
      if(len(self.bs_ave)>0):
        self.bond_final  = self.bs_ave[0]
        self.angle_final = self.as_ave[0]
    ###
    xrs = model.xray_structure
    ###
    if(rigid_body_shift_accumulator is not None):
       self.rigid_body_shift_accumulator = rigid_body_shift_accumulator
    ###
    t2 = time.time()
    time_collect_and_process += (t2 - t1)
    self.call_back(model, fmodel)

  def call_back (self, model, fmodel, method="monitor_collect") :
    if self.call_back_handler is not None and callable(self.call_back_handler) :
      self.call_back_handler(self, model, fmodel, method)

  def show(self, out=None, remark=""):
    global time_collect_and_process
    t1 = time.time()
    if(out is None): out = self.out
    separator = "-"*72
    if(not self.short):
       #
       if(self.rigid_body_shift_accumulator is not None):
          print >> out, remark + "Information about total rigid body shift of selected groups:"
          show_rigid_body_rotations_and_translations(
            out=out,
            prefix=remark,
            frame=" ",
            euler_angle_convention
              =self.rigid_body_shift_accumulator.euler_angle_convention,
            rotations=self.rigid_body_shift_accumulator.rotations,
            translations=self.rigid_body_shift_accumulator.translations)
    #
    if(self.short):
       print >> out, remark
       print >> out, remark + "*********** REFINEMENT STATISTICS STEP BY STEP: NEUTRON DATA ***********"
    else:
       print >> out, remark + "****************** REFINEMENT STATISTICS STEP BY STEP ******************"
       print >> out, remark + "leading digit, like 1_, means number of macro-cycle                     "
       print >> out, remark + "0    : statistics at the very beginning when nothing is done yet        "
       if(self.params.main.bulk_solvent_and_scale):
          print >> out, remark + "1_bss: bulk solvent correction and/or (anisotropic) scaling             "
       if("individual_sites" in self.params.refine.strategy):
          print >> out, remark + "1_xyz: refinement of coordinates                                        "
       if("individual_adp" in self.params.refine.strategy):
          print >> out, remark + "1_adp: refinement of ADPs (Atomic Displacement Parameters)              "
       if(self.params.main.simulated_annealing):
          print >> out, remark + "1_sar: simulated annealing refinement of x,y,z                          "
       if(self.params.main.ordered_solvent):
          print >> out, remark + "1_wat: ordered solvent update (add / remove)                            "
       if("rigid_body" in self.params.refine.strategy):
          print >> out, remark + "1_rbr: rigid body refinement                                            "
       if("group_adp" in self.params.refine.strategy):
          print >> out, remark + "1_gbr: group B-factor refinement                                        "
       if("occupancies" in self.params.refine.strategy):
          print >> out, remark + "1_occ: refinement of occupancies                                        "
       print >> out, remark + separator
    #
    #
    a,b,c,d,e,f,g,h,i,j = [None,]*10
    print >> out, remark + \
      " R-factors, x-ray target values and norm of gradient of x-ray target"
    print >> out, remark + \
      " stage     r-work r-free  xray_target_w  xray_target_t"
    format = remark + "%9s  %6.4f %6.4f %14.6e %14s"
    for a,b,c,d,e in zip(self.steps,
                           self.r_works,
                           self.r_frees,
                           self.targets_w,
                           self.targets_t):
        if (e is None): e = "None"
        else:           e = "%14.6e" % e
        print >> out, format % (a,b,c,d,e)
    print >> out, remark + separator
    #
    #
    a,b,c,d,e,f,g,h,i,j = [None,]*10
    if(len(self.wxcs) > 0):
       print >> out, remark + " Weights for target T = Exray * wxc * wxc_scale + Echem * wc and"
       print >> out, remark + " angles between gradient vectors, eg. (d_Exray/d_sites, d_Echem/d_sites)"
       print >> out, remark + \
          " stage              wxc          wxu  wxc_sc  wxu_sc /_gxc,gc /_gxu,gu"
       format = remark + "%9s  %12.4e %12.4e%8.3f%8.3f%9.3f%9.3f"
       for a,b,c,d,e,f,g in zip(self.steps,
                                self.wxcs,
                                self.wxus,
                                self.wxc_scales,
                                self.wxu_scales,
                                self.angles_xc,
                                self.angles_xu):
           print >> out, format % (a,b,c,d,e,f,g)
       print >> out, remark + separator
    else:
       assert len(self.wxcs)+len(self.wxus)+len(self.wxc_scales)+\
              len(self.wxu_scales)+len(self.angles_xc)+len(self.angles_xu) == 0
    #
    #
    a,b,c,d,e,f,g,h,i,j = [None,]*10
    print >> out, remark + \
     " stage     k_sol   b_sol     b11     b22     b33     b12     b13     b23"
    format = remark + "%9s  %5.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f"
    for a,b,c,d in zip(self.steps, self.k_sols, self.b_sols, self.b_anisos):
        print >> out, format % (a,b,c,d[0],d[1],d[2],d[3],d[4],d[5])
    print >> out, remark + separator
    #
    #
    a,b,c,d,e,f,g,h,i,j = [None,]*10
    print >> out, remark + " stage     <pher>    fom   alpha         beta"
    format = remark + "%9s %7.3f %6.4f %7.4f %12.3f"
    for a,b,c,d,e in zip(self.steps,
                         self.phers,
                         self.foms,
                         self.alphas,
                         self.betas):
        print >> out, format % (a,b,c,d,e)
    print >> out, remark + separator
    #
    #
    if(not self.short and len(self.as_ave)>0):
       a,b,c,d,e,f,g,h,i,j = [None,]*10
       if(len(self.wcs) > 0):
          print >> out, remark + \
          " stage       angl   bond   chir   dihe   plan   repu  geom_target     wc"
          format = remark + "%9s %7.3f%7.3f%7.3f%7.3f%7.3f%7.3f %12.4e%7.2f"
          for a,b,c,d,e,f,g,h,i in zip(self.steps,
                                       self.as_ave,
                                       self.bs_ave,
                                       self.cs_ave,
                                       self.ds_ave,
                                       self.ps_ave,
                                       self.rs_ave,
                                       self.targets_c,
                                       self.wcs):
              print >> out, format % (a,b,c,d,e,f,g,h,i)
       else:
          print >> out, remark + \
          " stage       angl   bond   chir   dihe   plan   repu  geom_target"
          format = remark + "%9s %7.3f%7.3f%7.3f%7.3f%7.3f%7.3f %12.4e"
          for a,b,c,d,e,f,g,h in zip(self.steps,
                                       self.as_ave,
                                       self.bs_ave,
                                       self.cs_ave,
                                       self.ds_ave,
                                       self.ps_ave,
                                       self.rs_ave,
                                       self.targets_c):
              print >> out, format % (a,b,c,d,e,f,g,h)

       print >> out, remark + separator
    #
    #
    if(not self.short and len(self.as_ave)>0):
       a,b,c,d,e,f,g,h,i,j = [None,]*10
       print >> out, remark + "                      Maximal deviations:"
       print >> out, remark + \
               " stage       angl   bond   chir   dihe   plan   repu       |grad|"
       format = remark + "%9s %7.3f%7.3f%7.3f%7.3f%7.3f%7.3f %12.4e"
       for a,b,c,d,e,f,g,h in zip(self.steps,
                                  self.as_max,
                                  self.bs_max,
                                  self.cs_max,
                                  self.ds_max,
                                  self.ps_max,
                                  self.rs_min,
                                  self.gs_c_norm):
           print >> out, format % (a,b,c,d,e,f,g,h)
       print >> out, remark + separator
       #
       #
       a,b,c,d,e,f,g,h,i,j = [None,]*10
       if(len(self.bs_iso_max_s) == 0):
          print >> out, remark + " stage      b_max   b_min   b_ave"
       else:
          print >> out, remark + "           |-----overall-----|---macromolecule----|------solvent-------|"
          print >> out, remark + "  stage    b_max  b_min  b_ave  b_max  b_min  b_ave  b_max  b_min  b_ave"
       if(len(self.bs_iso_max_s) > 0):
          format = remark + "%9s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f"
          for a,b,c,d,e,f,g,h,i,j in zip(self.steps,
                                         self.bs_iso_max_a,
                                         self.bs_iso_min_a,
                                         self.bs_iso_ave_a,
                                         self.bs_iso_max_p,
                                         self.bs_iso_min_p,
                                         self.bs_iso_ave_p,
                                         self.bs_iso_max_s,
                                         self.bs_iso_min_s,
                                         self.bs_iso_ave_s):
              print >> out, format % (a,b,c,d,e,f,g,h,i,j)
       if(len(self.bs_iso_max_s) == 0):
          format = remark + "%9s %7.2f %7.2f %7.2f"
          for a,b,c,d in zip(self.steps,
                                         self.bs_iso_max_a,
                                         self.bs_iso_min_a,
                                         self.bs_iso_ave_a):
              print >> out, format % (a,b,c,d)
       print >> out, remark + separator
    #
    #
    if(not self.short):
       a,b,c,d,e,f,g,h,i,j = [None,]*10
       if(len(self.ds_rc_min) > 0 and len(self.dopts_ic) == 0):
          if(self.model_ref is not None):
             print >> out, remark + "Stage        Deviation of refined model from"
             print >> out, remark + "                     start model        reference model"
             print >> out, remark + "                  max    min   mean    max    min   mean"
             format = remark + "%11s   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f"
       if(len(self.ds_rc_min) == 0 and len(self.dopts_ic) == 0):
          print >> out, remark + " stage        Deviation of refined"
          print >> out, remark + "              model from start model"
          print >> out, remark + "               max      min     mean"
          format = remark + "%9s %8.3f %8.3f %8.3f"
       if(len(self.ds_rc_min) > 0 and len(self.dopts_ic) > 0 and len(self.dopts_rc) > 0):
          if(self.model_ref is not None):
             print >> out, remark + "Stage                      Deviation of refined model from"
             print >> out, remark + "                 Start model                         Reference model"
             print >> out, remark + "                   max      min     mean      opt      max      min     mean      opt"
             format = remark + "%11s   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f"
       if(len(self.ds_rc_min) == 0 and len(self.dopts_ic) > 0):
             print >> out, remark + "Stage                Deviation of refined"
             print >> out, remark + "                     model from start model"
             print >> out, remark + "                  max      min     mean    opt"
             format = remark + "%11s %8.3f %8.3f %8.3f %8.3f"
       if(len(self.ds_rc_min) == 0 and len(self.dopts_ic) == 0):
          for a,b,c,d in zip(self.steps,
                             self.ds_ic_max,
                             self.ds_ic_min,
                             self.ds_ic_ave):
              print >> out, format % (a,b,c,d)
       if(len(self.ds_rc_min) > 0 and len(self.dopts_ic) == 0):
          for a,b,c,d,e,f,g in zip(self.steps,
                                   self.ds_ic_max,
                                   self.ds_ic_min,
                                   self.ds_ic_ave,
                                   self.ds_rc_max,
                                   self.ds_rc_min,
                                   self.ds_rc_ave):
              print >> out, format % (a,b,c,d,e,f,g)
       if(len(self.ds_rc_min) == 0 and len(self.dopts_ic) > 0):
          for a,b,c,d,e in zip(self.steps,
                                   self.ds_ic_max,
                                   self.ds_ic_min,
                                   self.ds_ic_ave,
                                   self.dopts_ic):
              print >> out, format % (a,b,c,d,e)
       if(len(self.ds_rc_min) > 0 and len(self.dopts_ic) > 0 and len(self.dopts_rc) > 0):
          for a,b,c,d,e,f,g,h,i in zip(self.steps,
                                       self.ds_ic_max,
                                       self.ds_ic_min,
                                       self.ds_ic_ave,
                                       self.dopts_ic,
                                       self.ds_rc_max,
                                       self.ds_rc_min,
                                       self.ds_rc_ave,
                                       self.dopts_rc):
              print >> out, format % (a,b,c,d,e,f,g,h,i)
       print >> out, remark + separator
    #
    #
    a,b,c,d,e,f,g,h,i,j = [None,]*10
    if(len(self.n_zeros) > 0):
       print >> out, remark + "Fobs* statistics"
       print >> out, remark + "Stage        number_of_f_star_zero %_from_total"
       format = remark + "%11s                 %6d %5.2f"
       for a,b,c in zip(self.steps, self.n_zeros, self.n_zeros_p):
         print >> out, format % (a, b ,c)
    #
    #
    if 0:
      a,b,c,d,e,f,g,h,i,j = [None,]*10
      print >> out, remark + \
                        " stage         k1_w     k1_t     k3_w     k3_t scale_ml"
      format = remark + "%9s  %8.4f %8.4f %8.4f %8.4f %8.4f"
      for a,b,c,d,e,f in zip(self.steps,
                             self.k1s_w,
                             self.k1s_t,
                             self.k3s_w,
                             self.k3s_t,
                             self.scale_ml):
          print >> out, format % (a,b,c,d,e,f)
      print >> out, remark + separator
    #
    #
    if(not self.short):
       a,b,c,d,e,f,g,h,i,j = [None,]*10
       proceed = False
       n_solv_0 = self.n_solv[0]
       for n in self.n_solv:
           if(n != n_solv_0):
              proceed = True
              break
       if(proceed):
          print >> out, remark + " stage  number of ordered solvent"
          format = remark + "%11s  %15d"
          for a,b in zip(self.steps, self.n_solv):
              print >> out, format % (a,b)
          print >> out, remark + separator
    #
    #
#    if(not self.short):
#       a,b,c,d,e,f,g,h,i,j = [None,]*10
#       if(len(self.wus) > 0):
#          print >> out, remark + " stage      adp_target      wu"
#          format = remark + "%9s %12.4e %7.3f"
#          for a,b,c in zip(self.steps,self.tus,self.wus):
#              print >> out, format % (a,b,c)
#       else:
#          print >> out, remark + " stage      adp_target"
#          format = remark + "%9s %12.4e"
#          for a,b in zip(self.steps,self.tus):
#              print >> out, format % (a,b)
#       print >> out, remark + separator
    out.flush()
    #
    #
    t2 = time.time()
    time_collect_and_process += (t2 - t1)

class relative_errors(object):
  def __init__(self, model, model_ini, model_ref, compute_optimal_errors):
    #XXX fix this later
    if(model.xray_structure.scatterers().size() !=
                                 model_ini.xray_structure.scatterers().size()):
       xray_structure   = model.xray_structure_macromolecule()
       xray_structure_0 = model_ini.xray_structure_macromolecule()
    else:
       xray_structure   = model.xray_structure
       xray_structure_0 = model_ini.xray_structure
    array_of_distances_between_each_atom_ini_curr = \
                              flex.sqrt(xray_structure.difference_vectors_cart(
                                                       xray_structure_0).dot())
    self.d_ic_min = \
             flex.min_default( array_of_distances_between_each_atom_ini_curr,0)
    self.d_ic_max = \
             flex.max_default( array_of_distances_between_each_atom_ini_curr,0)
    self.d_ic_ave = \
             flex.mean_default(array_of_distances_between_each_atom_ini_curr,0)
    if(model_ref is not None):
       #XXX fix this later
       if(model.xray_structure.scatterers().size() !=
                                 model_ref.xray_structure.scatterers().size()):
          xray_structure   = model.xray_structure_macromolecule()
          xray_structure_ref = model_ref.xray_structure_macromolecule()
       else:
          xray_structure   = model.xray_structure
          xray_structure_ref = model_ref.xray_structure
       array_of_distances_between_each_atom_ini_ref = \
                             flex.sqrt(xray_structure.difference_vectors_cart(
                                                    xray_structure_ref).dot())
       self.d_rc_min = \
             flex.min_default( array_of_distances_between_each_atom_ini_ref,0)
       self.d_rc_max = \
             flex.max_default( array_of_distances_between_each_atom_ini_ref,0)
       self.d_rc_ave = \
             flex.mean_default(array_of_distances_between_each_atom_ini_ref,0)
    if(compute_optimal_errors):
       assert model.xray_structure.scatterers().size() == \
                                   model_ini.xray_structure.scatterers().size()
       xrs = model.xray_structure
       xrs_ini = model_ini.xray_structure
       self.dopt_ic = max_lik.sasha_error_calculator(
                                  r1f  = xrs_ini.sites_frac(),
                                  r1c  = xrs_ini.sites_cart(),
                                  r2f  = xrs.sites_frac(),
                                  r2c  = xrs.sites_cart(),
                                  lab1 = xrs_ini.scatterers().extract_labels(),
                                  lab2 = xrs.scatterers().extract_labels(),
                                  uc   = xrs_ini.unit_cell(),
                                  sg   = xrs_ini.space_group(),
                                  rad  = 3.0).doptimal()
       if(model_ref is not None):
          xrs = model.xray_structure
          xrs_ref = model_ref.xray_structure
          assert xrs.scatterers().size() == xrs_ref.scatterers().size()
          self.dopt_rc = max_lik.sasha_error_calculator(
                                  r1f  = xrs_ref.sites_frac(),
                                  r1c  = xrs_ref.sites_cart(),
                                  r2f  = xrs.sites_frac(),
                                  r2c  = xrs.sites_cart(),
                                  lab1 = xrs_ref.scatterers().extract_labels(),
                                  lab2 = xrs.scatterers().extract_labels(),
                                  uc   = xrs_ref.unit_cell(),
                                  sg   = xrs_ref.space_group(),
                                  rad  = 3.0).doptimal()
