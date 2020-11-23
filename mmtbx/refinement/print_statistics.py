from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import scitbx.math.euler_angles
from scitbx import matrix
from libtbx.utils import format_cpu_times, getenv_bool
from libtbx import adopt_init_args, slots_getstate_setstate
import sys, time
from libtbx import str_utils
from libtbx.str_utils import prefix_each_line_suffix, format_value
from libtbx import introspection
import math
from cctbx import xray
import cctbx.xray.structure_factors.global_counters
from libtbx import easy_pickle
from itertools import count
from libtbx import group_args
from six.moves import zip

enable_show_process_info = getenv_bool(
  "MMTBX_PRINT_STATISTICS_ENABLE_SHOW_PROCESS_INFO")

time_collect_and_process = 0.0

def show_times(out = None):
  if(out is None): out = sys.stdout
  total = time_collect_and_process
  if(total > 0.01):
     print("Collect and process                      = %-7.2f" % time_collect_and_process, file=out)
  return total

def show_process_info(out):
  print("\\/"*39, file=out)
  introspection.virtual_memory_info().show_if_available(out=out, show_max=True)
  xray.structure_factors.global_counters.show(out=out)
  print(format_cpu_times(), file=out)
  print("/\\"*39, file=out)
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
  fill_rl = fill_len//2
  fill_r = fill_rl
  fill_l = fill_rl
  if (fill_rl*2 != fill_len): fill_r +=1
  str1 = "\n"+"*"*(fill_l-1)+" REFINEMENT MACRO_CYCLE "+macro_cycle+" OF "
  str2 = number_of_macro_cycles+" "+"*"*(fill_r)+"\n"
  out_string = str1+str2
  print(out_string, file=out)
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
  print(prefix_each_line_suffix(
    prefix=prefix+frame, lines_as_one_string=
"                            rotation (deg)                 translation (A)   "
"\n"
"                         %s              total           xyz          total "
      % euler_angle_convention, suffix=frame), file=out)
  for i,r,t in zip(count(1), rotations, translations):
    r = list(r)
    r.reverse()
    r_total = abs(scitbx.math.r3_rotation_axis_and_angle_from_matrix(
      r=euler_angles_as_matrix(*r)).angle(deg=True))
    t_total = abs(matrix.col(t))
    print((prefix + frame +
      " group %4d: %8.3f %8.3f %8.3f %7.2f  %6.2f %6.2f %6.2f %6.2f "
        % tuple([i] + r + [r_total] + list(t) + [t_total])
      + frame).rstrip(), file=out)
  out.flush()

# these are the steps we actually want to display in the GUI
show_actions = {
   "bss" : "bss",
   "sol" : "sol",
   "ion" : "ion",
   "rbr" : "rbr",
   "realsrl" : "rsrl",
   "realsrg" : "rsrg",
   "den" : "den",
   "tardy" : "SA",
   "sacart" : "SA",
   "xyzrec" : "xyz",
   "adp" : "adp",
   "occ" : "occ",
   "fp_fdp" : "anom",
}

class refinement_monitor(object):
  __arrays__ = [
    "steps",
    "r_works",
    "r_frees",
    "geom",
    "bs_iso_max_a",
    "bs_iso_min_a",
    "bs_iso_ave_a",
    "n_solv",
    "shifts",
  ]

  def __init__(
        self,
        params,
        out=None,
        neutron_refinement = None,
        call_back_handler=None,
        is_neutron_monitor=False):
    adopt_init_args(self, locals())
    if (self.out is None): self.out = sys.stdout
    self.wilson_b = None
    self.bond_start = None
    self.angle_start= None
    self.bond_final = None
    self.angle_final= None
    self.rigid_body_shift_accumulator = None
    self.sites_cart_start = None
    for name in self.__arrays__ :
      setattr(self, name, [])
    self.is_amber_monitor = False
    self.geom = group_args(bonds=[], angles=[])

  def dump_statistics(self, file_name):
    stats = {}
    for name in self.__arrays__ :
      stats[name] = getattr(self, name)
    easy_pickle.dump(file_name, stats)

  def collect(self, model,
                    fmodel,
                    step,
                    wilson_b = None,
                    rigid_body_shift_accumulator = None):
    global time_collect_and_process
    t1 = time.time()
    if(self.sites_cart_start is None):
      self.sites_cart_start = model.get_sites_cart()
    sites_cart_curr = model.get_sites_cart()
    if(sites_cart_curr.size()==self.sites_cart_start.size()):
      self.shifts.append(
        flex.mean(flex.sqrt((self.sites_cart_start-sites_cart_curr).dot())))
    else: self.shifts.append("n/a")
    if(wilson_b is not None): self.wilson_b = wilson_b
    self.steps.append(step)
    self.r_works.append(fmodel.r_work())
    self.r_frees.append(fmodel.r_free())
    use_amber = False
    if hasattr(self.params, "amber"): # loaded amber scope
      use_amber = self.params.amber.use_amber
      self.is_amber_monitor=use_amber
    use_afitt = False
    if hasattr(self.params, "afitt"): # loaded amber scope
      use_afitt = self.params.afitt.use_afitt
    general_selection = None
    if use_afitt:
      from mmtbx.geometry_restraints import afitt
      general_selection = afitt.get_non_afitt_selection(
        model.restraints_manager,
        model.get_sites_cart(),
        model.get_hd_selection(),
        None)
    geom = model.geometry_statistics()
    if(geom is not None):
      self.geom.bonds.append(geom.bond().mean)
      self.geom.angles.append(geom.angle().mean)
    hd_sel = None
    if(not self.neutron_refinement and not self.is_neutron_monitor):
      hd_sel = model.get_hd_selection()
    b_isos = model.get_xray_structure().extract_u_iso_or_u_equiv() * math.pi**2*8
    if(hd_sel is not None): b_isos = b_isos.select(~hd_sel)
    self.bs_iso_max_a.append(flex.max_default( b_isos, 0))
    self.bs_iso_min_a.append(flex.min_default( b_isos, 0))
    self.bs_iso_ave_a.append(flex.mean_default(b_isos, 0))
    self.n_solv.append(model.number_of_ordered_solvent_molecules())
    if(len(self.geom.bonds)>0):
      if([self.bond_start,self.angle_start].count(None) == 2):
        if(len(self.geom.bonds)>0):
          self.bond_start  = self.geom.bonds[0]
          self.angle_start = self.geom.angles[0]
      if(len(self.geom.bonds)>0):
        self.bond_final  = self.geom.bonds[len(self.geom.bonds)-1]
        self.angle_final = self.geom.angles[len(self.geom.angles)-1]
      elif(len(self.geom)==1):
        self.bond_final  = self.geom.bonds[0]
        self.angle_final = self.geom.angles[0]
    if(rigid_body_shift_accumulator is not None):
      self.rigid_body_shift_accumulator = rigid_body_shift_accumulator
    t2 = time.time()
    time_collect_and_process += (t2 - t1)
    self.call_back(model, fmodel, method=step)

  def call_back(self, model, fmodel, method="monitor_collect"):
    if self.call_back_handler is not None and callable(self.call_back_handler):
      self.call_back_handler(self, model, fmodel, method)

  def show(self, out=None, remark=""):
    global time_collect_and_process
    t1 = time.time()
    max_step_len = max([len(s) for s in self.steps])
    if(out is None): out = self.out
    separator = "-"*72
    if(self.rigid_body_shift_accumulator is not None):
      print(remark + "Information about total rigid body shift of selected groups:", file=out)
      show_rigid_body_rotations_and_translations(
        out=out,
        prefix=remark,
        frame=" ",
        euler_angle_convention
          =self.rigid_body_shift_accumulator.euler_angle_convention,
        rotations=self.rigid_body_shift_accumulator.rotations,
        translations=self.rigid_body_shift_accumulator.translations)
    #
    print(remark + "****************** REFINEMENT STATISTICS STEP BY STEP ******************", file=out)
    print(remark + "leading digit, like 1_, means number of macro-cycle                     ", file=out)
    print(remark + "0    : statistics at the very beginning when nothing is done yet        ", file=out)
    if(self.params.main.bulk_solvent_and_scale):
       print(remark + "1_bss: bulk solvent correction and/or (anisotropic) scaling             ", file=out)
    if("individual_sites" in self.params.refine.strategy):
       print(remark + "1_xyz: refinement of coordinates                                        ", file=out)
    if("individual_adp" in self.params.refine.strategy):
       print(remark + "1_adp: refinement of ADPs (Atomic Displacement Parameters)              ", file=out)
    if(self.params.main.simulated_annealing):
       print(remark + "1_sar: simulated annealing refinement of x,y,z                          ", file=out)
    if(self.params.main.ordered_solvent):
       print(remark + "1_wat: ordered solvent update (add / remove)                            ", file=out)
    if("rigid_body" in self.params.refine.strategy):
       print(remark + "1_rbr: rigid body refinement                                            ", file=out)
    if("group_adp" in self.params.refine.strategy):
       print(remark + "1_gbr: group B-factor refinement                                        ", file=out)
    if("occupancies" in self.params.refine.strategy):
       print(remark + "1_occ: refinement of occupancies                                        ", file=out)
    print(remark + separator, file=out)
    #
    has_bonds_angles=True
    if len(self.geom.bonds):
      print(remark + \
        " stage r-work r-free bonds angles b_min b_max b_ave n_water shift", file=out)
      format = remark + "%s%ds"%("%",max_step_len)+\
        " %6.4f %6.4f %5.3f %6.3f %5.1f %5.1f %5.1f %3d %s"
    else:
      print(remark + \
        " stage       r-work r-free b_min b_max b_ave n_water shift", file=out)
      format = remark + "%s%ds"%("%",max_step_len)+\
        " %6.4f %6.4f %5.1f %5.1f %5.1f %3d %s"
    for a,b,c, d1,d2, e,f,g,h,i in zip(self.steps,
                                   self.r_works,
                                   self.r_frees,
                                   self.geom.bonds,
                                   self.geom.angles,
                                   self.bs_iso_min_a,
                                   self.bs_iso_max_a,
                                   self.bs_iso_ave_a,
                                   self.n_solv,
                                   self.shifts):
        if(type(1.)==type(i)): i = "     "+str("%5.3f"%i)
        else: i = "%9s"%i
        if has_bonds_angles:
          print(format % (a,b,c,d1,d2,e,f,g,h,i), file=out)
        else:
          print(format % (a,b,c,e,f,g,h,i), file=out)
    print(remark + separator, file=out)
    out.flush()
    #
    t2 = time.time()
    time_collect_and_process += (t2 - t1)

  def format_stats_for_phenix_gui(self):
    steps = []
    r_works = []
    r_frees = []
    as_ave = []
    bs_ave = []
    for i_step, label in enumerate(self.steps):
      label = label.replace(":", "")
      fields = label.split("_")
      if (len(fields) < 2):
        steps.append(label)
      else :
        cycle = fields[0]
        action = "_".join(fields[1:])
        action_label = show_actions.get(action, None)
        if (action_label is None) : continue
        steps.append(cycle + "_" + action_label)
      r_works.append(self.r_works[i_step])
      r_frees.append(self.r_frees[i_step])
      if (self.geom is not None) and (len(self.geom.bonds) != 0):
        as_ave.append(self.geom.angles[i_step])
        bs_ave.append(self.geom.bonds[i_step])
      else :
        as_ave.append(None)
        bs_ave.append(None)
    return stats_table(
      steps=steps,
      r_works=r_works,
      r_frees=r_frees,
      as_ave=as_ave,
      bs_ave=bs_ave,
      neutron_flag=self.is_neutron_monitor)

  def show_current_r_factors_summary(self, out, prefix=""):
    if (len(self.steps) == 0):
      return
    last_step = self.steps[-1].replace(":", "")
    fields = last_step.split("_")
    action_label = None
    if (len(fields) >= 2):
      cycle = fields[0]
      action = "_".join(fields[1:])
      action_label = show_actions.get(action, None)
      if (action_label is None):
        return
      action_label = cycle + "_" + action_label
    if (action_label is None):
      if (len(self.steps) == 1):
        action_label = "start"
      else :
        action_label = "end"
    print("%s%-6s  r_work=%s  r_free=%s" % (prefix,
      action_label,
      format_value("%.4f", self.r_works[-1]),
      format_value("%.4f", self.r_frees[-1])), file=out)

class stats_table(slots_getstate_setstate):
  __slots__ = [
    "steps",
    "r_works",
    "r_frees",
    "as_ave",
    "bs_ave",
    "neutron_flag",
  ]
  def __init__(self, **kwds):
    for attr in self.__slots__ :
      setattr(self, attr, kwds.get(attr))

# we need something simpler for the Phenix GUI...
class coordinate_shifts(object):
  def __init__(self, hierarchy_start, hierarchy_end):
    from scitbx.array_family import flex
    self.hierarchy_shifted = hierarchy_end.deep_copy()
    atoms_shifted = self.hierarchy_shifted.atoms()
    coords_start = {}
    for atom in hierarchy_start.atoms():
      id_str = atom.fetch_labels().id_str()
      coords_start[id_str] = atom.xyz
    def get_distance(xyz1, xyz2):
      x1,y1,z1 = xyz1
      x2,y2,z2 = xyz2
      return math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
    for i_seq, atom in enumerate(atoms_shifted):
      id_str = atom.fetch_labels().id_str()
      if (id_str in coords_start):
        atom.b = get_distance(coords_start[id_str], atom.xyz)
      else :
        atom.b = -1.0

  def get_shifts(self):
    return self.hierarchy_shifted.atoms().extract_b()

  def min_max_mean(self):
    shifts = self.hierarchy_shifted.atoms().extract_b()
    shifts = shifts.select(shifts >= 0)
    return shifts.min_max_mean()

  def save_pdb_file(self, file_name):
    f = open(file_name, "w")
    f.write(self.hierarchy_shifted.as_pdb_string())
    f.close()

class trajectory_output(object):
  """
  Callback object for saving the intermediate results of refinement as a stack
  of PDB and MTZ files.  Equivalent to the interactivity with Coot in the
  Phenix GUI, but intended for command-line use and demonstrative purposes.
  """
  def __init__(self, file_base="refine", filled_maps=True, log=sys.stdout,
      verbose=True):
    adopt_init_args(self, locals())
    self._i_trajectory = 0

  def __call__(self, monitor, model, fmodel, method="monitor_collect"):
    import iotbx.map_tools
    self._i_trajectory += 1
    file_base = "%s_traj_%d" % (self.file_base, self._i_trajectory)
    pdb_hierarchy = model.get_hierarchy()
    two_fofc_map_coeffs = fmodel.map_coefficients(map_type="2mFo-DFc",
      fill_missing=self.filled_maps)
    fofc_map_coeffs = fmodel.map_coefficients(map_type="mFo-DFc")
    iotbx.map_tools.write_map_coeffs(
      fwt_coeffs=two_fofc_map_coeffs,
      delfwt_coeffs=fofc_map_coeffs,
      file_name=file_base+".mtz")
    f = open(file_base + ".pdb", "w")
    f.write(pdb_hierarchy.as_pdb_string(
      crystal_symmetry=model.get_xray_structure()))
    f.close()
    print("wrote model to %s.pdb" % file_base, file=self.log)
    print("wrote map coefficients to %s.mtz" % file_base, file=self.log)

class annealing_callback(object):
  def __init__(self, model, monitor):
    self.model = model
    self.monitor = monitor

  def __call__(self, fmodel):
    self.monitor.call_back(model=self.model,
      fmodel=fmodel,
      method="anneal")
