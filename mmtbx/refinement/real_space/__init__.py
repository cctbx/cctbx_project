from __future__ import division
from __future__ import print_function
from libtbx import adopt_init_args, group_args
from scitbx.array_family import flex
from scitbx.matrix import rotate_point_around_axis
import time, sys
from cctbx import maptbx
import mmtbx.utils
from mmtbx.rotamer.rotamer_eval import RotamerEval
import iotbx.pdb
from cctbx import miller
from libtbx.str_utils import format_value
import mmtbx.model.statistics
import libtbx.load_env
from mmtbx.utils import rotatable_bonds
from cctbx.eltbx import tiny_pse
from cctbx import eltbx
from libtbx.test_utils import approx_equal
from mmtbx.maps.correlation import five_cc
import mmtbx.model
from cctbx.eltbx import tiny_pse
from cctbx import eltbx

# Test utils (common to all tests in this folder).
# Perhaps to move to mmtbx/regression.

import scitbx.math
from libtbx.utils import null_out

def setup_test(pdb_answer, pdb_poor, i_pdb, d_min, resolution_factor,
               pdb_for_map = None):
  rotamer_manager = mmtbx.idealized_aa_residues.rotamer_manager.load()
  sin_cos_table = scitbx.math.sin_cos_table(n=10000)
  #
  pip = mmtbx.model.manager.get_default_pdb_interpretation_params()
  pip.pdb_interpretation.link_distance_cutoff=999
  # answer
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_answer)
  model_answer = mmtbx.model.manager(model_input=pdb_inp, process_input=True,
    log=null_out(), pdb_interpretation_params=pip)
  with open("answer_%s.pdb"%str(i_pdb), "w") as the_file:
    the_file.write(model_answer.model_as_pdb())
  #
  model_ = model_answer.deep_copy()
  if(pdb_for_map is not None):
    pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_for_map)
    model_ = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  #
  xrs_answer = model_.get_xray_structure()
  f_calc = xrs_answer.structure_factors(d_min = d_min).f_calc()
  fft_map = f_calc.fft_map(resolution_factor=resolution_factor)
  fft_map.apply_sigma_scaling()
  target_map = fft_map.real_map_unpadded()
  mtz_dataset = f_calc.as_mtz_dataset(column_root_label = "FCmap")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "answer_%s.mtz"%str(i_pdb))
  # poor
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_poor)
  model_poor = mmtbx.model.manager(model_input=pdb_inp, log=null_out(),
    pdb_interpretation_params=pip)
  with open("poor_%s.pdb"%str(i_pdb), "w") as the_file:
    the_file.write(model_poor.model_as_pdb())
  #
  return group_args(
    rotamer_manager  = rotamer_manager,
    sin_cos_table    = sin_cos_table,
    target_map       = target_map,
    xrs_poor         = model_poor.get_xray_structure(),
    ph_answer        = model_answer.get_hierarchy(),
    vdw              = model_answer.get_vdw_radii(),
    mon_lib_srv      = model_answer.get_mon_lib_srv(),
    ph_poor          = model_poor.get_hierarchy(),
    crystal_symmetry = f_calc.crystal_symmetry(),
    model_poor       = model_poor)

def check_sites_match(ph_answer, ph_refined, tol):
  s1 = flex.vec3_double()
  s2 = flex.vec3_double()
  for a1,a2 in zip(ph_answer.atoms(), ph_refined.atoms()):
    if((not a1.element.strip().upper() in ["H","D"]) and
       (not a2.element.strip().upper() in ["H","D"])):
      s1.append(a1.xyz)
      s2.append(a2.xyz)
  dist = flex.max(flex.sqrt((s1 - s2).dot()))
  print(dist, tol)
  assert dist < tol,  [dist, tol]

# End test utils

class rsr_model(object):
  def __init__(self,
               model,
               target_map_object=None):
    adopt_init_args(self, locals())
    self.unit_cell = self.model.crystal_symmetry().unit_cell()
    self.sites_cart_start = self.model.get_sites_cart()
    self.s1 = self.sites_cart_start.deep_copy()
    self.states_collector = mmtbx.utils.states(
      pdb_hierarchy  = self.model.get_hierarchy().deep_copy(),
      xray_structure = self.model.get_xray_structure().deep_copy_scatterers(),
      counter        = 1)
    self.states_collector.add(sites_cart = self.model.get_sites_cart())
    #
    self.five_cc = None
    self.rmsd_b = None
    self.rmsd_a = None
    self.dist_from_start = 0
    self.dist_from_previous = 0
    self.stats_evaluations = []
    self.cc_mask=None
    self.cc_box=None
    #
    self.initialize()

  def initialize(self):
    five_cc_o = five_cc(
      map               = self.target_map_object.map_data,
      xray_structure    = self.model.get_xray_structure(),
      d_min             = self.target_map_object.d_min,
      compute_cc_box    = True,
      compute_cc_image  = False,
      compute_cc_mask   = True,
      compute_cc_volume = False,
      compute_cc_peaks  = False).result
    self.cc_mask = five_cc_o.cc_mask
    self.cc_box  = five_cc_o.cc_box
    # XXX use model.statistics or better method of model!
    if(self.model.restraints_manager_available()):
      es = self.model.get_restraints_manager().geometry.energies_sites(
        sites_cart = self.model.get_sites_cart())
      self.rmsd_a = es.angle_deviations()[2]
      self.rmsd_b = es.bond_deviations()[2]

  def show(self, prefix="", log=None):
    if(log is None): log = sys.stdout
    print("%smodel-to-map fit, CC_mask: %-6.4f"%(prefix, self.cc_mask), file=log)
    print("%smoved from start:          %-6.4f"%(prefix, self.dist_from_start), file=log)
    gs = self.model.geometry_statistics()
    result = None
    if(gs is not None):
      gs.show(prefix=prefix, log=log, uppercase=False)
      result = gs.result()
    self.stats_evaluations.append(group_args(
      cc       = self.cc_mask,
      geometry = result))

  def update(self, xray_structure=None, sites_cart=None):
    assert [xray_structure, sites_cart].count(None)!=0
    s0 = self.sites_cart_start
    if(xray_structure is not None):
      s2 = xray_structure.sites_cart()
      self.model.set_xray_structure(xray_structure = xray_structure)
    elif(sites_cart is not None):
      s2 = sites_cart
      self.model.set_sites_cart(sites_cart=s2)
    else:
      s2 = self.model.get_sites_cart()
    self.dist_from_start    = flex.mean(flex.sqrt((s0 - s2).dot()))
    self.dist_from_previous = flex.mean(flex.sqrt((self.s1 - s2).dot()))
    self.initialize()
    self.states_collector.add(sites_cart = s2)
    self.s1 = self.model.get_sites_cart() # must be last

def flatten(l):
  if l is None: return None
  return sum(([x] if not (isinstance(x, list) or isinstance(x, flex.size_t)) else flatten(x) for x in l), [])

def need_sidechain_fit(
      residue,
      rotamer_evaluator,
      mon_lib_srv,
      unit_cell,
      f_map,
      outliers_only=False,
      fdiff_map=None,
      small_f_map=0.9):
  """
  Important: maps assumed to be sigma-scaled!
  """
  get_class = iotbx.pdb.common_residue_names_get_class
  assert get_class(residue.resname) == "common_amino_acid"
  if(residue.resname.strip().upper() in ["ALA", "GLY"]): return False
  ### If it is rotamer OUTLIER
  if(rotamer_evaluator.evaluate_residue(residue)=="OUTLIER"):
    return True
  if(outliers_only): return False
  ###
  cl = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
    residue         = residue,
    mon_lib_srv     = mon_lib_srv,
    backbone_sample = False).clusters
  if(len(cl)==0): return False
  # service functions
  def anal(x):
    for i,e in enumerate(x):
      if(e<0): return True
      r=None
      if(i+1<len(x)):
        e1=abs(x[i])
        e2=abs(x[i+1])
        if(e1>e2):
          if(e2!=0):
            r = e1/e2
        else:
          if(e1!=0):
            r = e2/e1
      if(r is not None and r>3): return True
    return False
  def anal2(x):
    for i,e in enumerate(x):
      if(e<-3.0): return True
    return False
  def anal3(x): return (flex.double(x)>=small_f_map).count(True)==len(x)
  #
  last = cl[0].vector[len(cl[0].vector)-1]
  vector = flatten(cl[0].vector)
  bs = residue.atoms().extract_b()
  #
  weights = []
  for el in residue.atoms().extract_element():
    std_lbl = eltbx.xray_scattering.get_standard_label(
      label=el, exact=True, optional=True)
    weights.append(tiny_pse.table(std_lbl).weight())
  #
  side_chain_sel = flex.size_t()
  main_chain_sel = flex.size_t()
  for i_seq, a in enumerate(list(residue.atoms())):
    if(a.name.strip().upper() not in ["N","CA","C","O","CB"]):
      side_chain_sel.append(i_seq)
    elif(a.name.strip().upper() in ["N","CA","C"]):
      main_chain_sel.append(i_seq)
  #
  sites_frac = unit_cell.fractionalize(residue.atoms().extract_xyz())
  #
  mv = []
  mv_orig = []
  if(fdiff_map is not None): diff_mv = []
  mv2 = flex.double()
  for v_ in vector:
    sf = sites_frac[v_]
    f_map_epi = f_map.eight_point_interpolation(sf)
    mv.append(     f_map_epi/weights[v_]*bs[v_])
    mv2.append(    f_map_epi/weights[v_])
    mv_orig.append(f_map_epi)
    if(fdiff_map is not None):
      diff_mv.append(fdiff_map.value_at_closest_grid_point(sf))
  f  = anal(mv)
  if(fdiff_map is not None): f2 = anal2(diff_mv)
  f3 = anal3(mv_orig)
  # main vs side chain
  mvbb = flex.double()
  mvbb_orig = flex.double()
  for mcs in main_chain_sel:
    sf = sites_frac[mcs]
    f_map_epi = f_map.eight_point_interpolation(sf)
    mvbb.append(     f_map_epi/weights[mcs])
    mvbb_orig.append(f_map_epi)
  f4 = flex.min(mvbb_orig)<small_f_map or flex.mean(mvbb)<flex.mean(mv2)
  c_id = "none"
  if(residue.parent() is not None):
    c_id = residue.parent().parent().id.strip()
  id_str = "%s_%s_%s"%(c_id, residue.resname.strip(), residue.resid().strip())
  #
  result = False
  if(fdiff_map is not None):
    if((f or f2 and not f3) and not f4): result = True
    else: result = False
  else:
    if((f       and not f3) and not f4): result = True
    else: result = False
  return result

class cluster(object):
  def __init__(self,
               axis,
               atoms_to_rotate,
               atom_names=None,
               vector=None,
               selection=None):
    adopt_init_args(self, locals())
    self.vector_flat = None

  def get_vector_flat(self):
    if(self.vector is not None):
      if(self.vector_flat is None):
        self.vector_flat = flex.size_t(flatten(self.vector))
    return self.vector_flat

  def show(self):
    if(self.atom_names is None): return
    cl = self
    an = self.atom_names
    print(cl.axis, ",".join([an[i].strip() for i in cl.axis]), \
        cl.atoms_to_rotate, \
        ",".join([an[i].strip() for i in cl.atoms_to_rotate]), "<>",\
        ",".join([an[i].strip() for i in cl.selection]), "<>",\
        ",".join([an[i].strip() for i in cl.get_vector_flat()]))

class aa_residue_axes_and_clusters(object):
  def __init__(self,
               residue,
               mon_lib_srv,
               backbone_sample):
    self.clusters               = []
    atoms                       = residue.atoms()
    atoms_as_list               = list(atoms)
    atom_names                  = atoms.extract_name()
    self.weights                = flex.double()
    self.clash_eval_selection   = flex.size_t()
    self.clash_eval_h_selection = flex.bool(len(atoms_as_list), False)
    self.rsr_eval_selection     = flex.size_t()
    # Backbone sample
    backrub_axis  = []
    backrub_atoms_to_rotate = []
    backrub_atoms_to_evaluate = []
    counter = 0 # XXX DOES THIS RELY ON ORDER?
    for atom in atoms:
      an = atom.name.strip().upper()
      ae = atom.element.strip().upper()
      if(ae in ["H","D"]):
        self.clash_eval_h_selection[counter]=True
      if(an in ["N", "C"]):
        backrub_axis.append(counter)
      else:
        backrub_atoms_to_rotate.append(counter)
      if(an in ["CA", "O", "CB"]):
        backrub_atoms_to_evaluate.append(counter)
      if(not an in ["CA", "O", "CB", "C", "N", "HA", "H"]):
        self.clash_eval_selection.append(counter)
        if(not ae in ["H","D"]):
          self.rsr_eval_selection.append(counter)
      std_lbl = eltbx.xray_scattering.get_standard_label(
        label=ae, exact=True, optional=True)
      self.weights.append(tiny_pse.table(std_lbl).weight())
      #
      counter += 1
    #
    if(backbone_sample):
      if(len(backrub_axis)==2 and len(backrub_atoms_to_evaluate)>0):
        self.clusters.append(cluster(
          axis            = flex.size_t(backrub_axis),
          atom_names      = atom_names,
          atoms_to_rotate = flex.size_t(backrub_atoms_to_rotate),
          selection       = flex.size_t(backrub_atoms_to_evaluate)))
    self.axes_and_atoms_aa_specific = \
      rotatable_bonds.axes_and_atoms_aa_specific(
        residue = residue, mon_lib_srv = mon_lib_srv)
    if(self.axes_and_atoms_aa_specific is not None):
      for i_aa, aa in enumerate(self.axes_and_atoms_aa_specific):
        if(i_aa == len(self.axes_and_atoms_aa_specific)-1):
          selection = flex.size_t(aa[1])
        else:
          selection = flex.size_t([aa[1][0]])
        # Exclude pure H or D rotatable groups
        elements_to_rotate = flex.std_string()
        for etr in aa[1]:
          elements_to_rotate.append(atoms_as_list[etr].element.strip())
        c_H = elements_to_rotate.count("H")
        c_D = elements_to_rotate.count("D")
        etr_sz = elements_to_rotate.size()
        if(c_H==etr_sz or c_D==etr_sz or c_H+c_D==etr_sz):
          continue
        #
        self.clusters.append(cluster(
          axis            = flex.size_t(aa[0]),
          atom_names      = atom_names,
          atoms_to_rotate = flex.size_t(aa[1]),
          selection       = flex.size_t(selection)))
      vector_selections = []
      if(len(self.clusters)>0):
        for i_aa, aa in enumerate(self.axes_and_atoms_aa_specific):
          for aa_ in aa[0]:
            if(not aa_ in vector_selections):
              vector_selections.append(aa_)
        vector_selections.append(
          self.clusters[len(self.clusters)-1].atoms_to_rotate)
        for cl in self.clusters:
          cl.vector = vector_selections

class residue_monitor(object):
  def __init__(self,
               residue,
               id_str,
               selection_all,
               selection_sidechain=None,
               selection_backbone=None,
               selection_c=None,
               selection_n=None,
               map_cc_sidechain=None,
               map_cc_backbone=None,
               map_cc_all=None,
               rotamer_status=None):
    adopt_init_args(self, locals())

  def format_info_string(self):
    return "%7s %6s    %6s     %6s %9s"%(
      self.id_str,
      format_value("%6.3f",self.map_cc_all),
      format_value("%6.3f",self.map_cc_backbone),
      format_value("%6.3f",self.map_cc_sidechain),
      self.rotamer_status)

class structure_monitor(object):
  def __init__(self,
               pdb_hierarchy,
               xray_structure,
               target_map_object=None,
               geometry_restraints_manager=None):
    adopt_init_args(self, locals())
    self.unit_cell = self.xray_structure.unit_cell()
    self.xray_structure = xray_structure.deep_copy_scatterers()
    self.xray_structure_start = xray_structure.deep_copy_scatterers()
    self.states_collector = mmtbx.utils.states(
      pdb_hierarchy  = self.pdb_hierarchy,
      xray_structure = self.xray_structure,
      counter        = 1)
    self.states_collector.add(sites_cart = self.xray_structure.sites_cart())
    self.rotamer_manager = RotamerEval()
    self.assert_pdb_hierarchy_xray_structure_sync()
    #
    self.five_cc = None
    self.map_cc_whole_unit_cell = None
    self.map_cc_around_atoms = None
    self.map_cc_per_atom = None
    self.rmsd_b = None
    self.rmsd_a = None
    self.dist_from_start = 0
    self.dist_from_previous = 0
    self.number_of_rotamer_outliers = 0
    self.residue_monitors = None
    self.stats_evaluations = []
    #
    self.initialize()

  def assert_pdb_hierarchy_xray_structure_sync(self):
    return #XXX
    sc1 = self.xray_structure.sites_cart()
    sc2 = self.pdb_hierarchy.atoms().extract_xyz()
    assert approx_equal(sc1, sc2, 1.e-3)

  def initialize(self):
    self.assert_pdb_hierarchy_xray_structure_sync()
    # residue monitors
    self.residue_monitors = []
    backbone_atoms = ["N","CA","C","O","CB"]
    get_class = iotbx.pdb.common_residue_names_get_class
    sites_cart = self.xray_structure.sites_cart()
    current_map = self.compute_map(xray_structure = self.xray_structure)
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          conformers = residue_group.conformers()
          if(len(conformers)>1): continue
          for conformer in residue_group.conformers():
            residue = conformer.only_residue()
            id_str="%s%s%s"%(chain.id,residue.resname,residue.resseq.strip())
            if(get_class(residue.resname) == "common_amino_acid"):
              residue_i_seqs_backbone  = flex.size_t()
              residue_i_seqs_sidechain = flex.size_t()
              residue_i_seqs_all       = flex.size_t()
              residue_i_seqs_c         = flex.size_t()
              residue_i_seqs_n         = flex.size_t()
              for atom in residue.atoms():
                an = atom.name.strip()
                bb = an in backbone_atoms
                residue_i_seqs_all.append(atom.i_seq)
                if(bb): residue_i_seqs_backbone.append(atom.i_seq)
                else:   residue_i_seqs_sidechain.append(atom.i_seq)
                if(an == "C"): residue_i_seqs_c.append(atom.i_seq)
                if(an == "N"): residue_i_seqs_n.append(atom.i_seq)
              sca = sites_cart.select(residue_i_seqs_all)
              scs = sites_cart.select(residue_i_seqs_sidechain)
              scb = sites_cart.select(residue_i_seqs_backbone)
              if(scs.size()==0): ccs = None
              else: ccs = self.map_cc(sites_cart=scs, other_map = current_map)
              if(sca.size()==0): cca = None
              else: cca = self.map_cc(sites_cart=sca, other_map = current_map)
              if(scb.size()==0): ccb = None
              else: ccb = self.map_cc(sites_cart=scb, other_map = current_map)
              self.residue_monitors.append(residue_monitor(
                residue             = residue,
                id_str              = id_str,
                selection_sidechain = residue_i_seqs_sidechain,
                selection_backbone  = residue_i_seqs_backbone,
                selection_all       = residue_i_seqs_all,
                selection_c         = residue_i_seqs_c,
                selection_n         = residue_i_seqs_n,
                map_cc_sidechain    = ccs,
                map_cc_backbone     = ccb,
                map_cc_all          = cca,
                rotamer_status= self.rotamer_manager.evaluate_residue(residue)))
            else:
              residue_i_seqs_all = residue.atoms().extract_i_seq()
              sca = sites_cart.select(residue_i_seqs_all)
              cca = self.map_cc(sites_cart=sca, other_map = current_map)
              self.residue_monitors.append(residue_monitor(
                residue       = residue,
                id_str        = id_str,
                selection_all = residue_i_seqs_all,
                map_cc_all    = cca))
    # globals
    self.five_cc = five_cc(
        map            = self.target_map_object.map_data,
        xray_structure = self.xray_structure,
        d_min          = self.target_map_object.d_min)
    self.map_cc_whole_unit_cell = self.map_cc(other_map = current_map)
    self.map_cc_around_atoms = self.map_cc(other_map = current_map,
      sites_cart = sites_cart)
    self.map_cc_per_atom = self.map_cc(other_map = current_map,
      sites_cart = sites_cart, per_atom = True)
    if(self.geometry_restraints_manager is not None):
      es = self.geometry_restraints_manager.energies_sites(sites_cart=sites_cart)
      self.rmsd_a = es.angle_deviations()[2]
      self.rmsd_b = es.bond_deviations()[2]
    self.dist_from_start = flex.mean(self.xray_structure_start.distances(
      other = self.xray_structure))
    self.number_of_rotamer_outliers = 0
    for r in self.residue_monitors:
      if(r.rotamer_status == "OUTLIER"):
        self.number_of_rotamer_outliers += 1
    self.assert_pdb_hierarchy_xray_structure_sync()

  def compute_map(self, xray_structure):
    self.assert_pdb_hierarchy_xray_structure_sync()
    mc = self.target_map_object.miller_array.structure_factors_from_scatterers(
      xray_structure = xray_structure).f_calc()
    fft_map = miller.fft_map(
      crystal_gridding     = self.target_map_object.crystal_gridding,
      fourier_coefficients = mc)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded()

  def map_cc_histogram_per_atom(self, radius=2, n_slots=10):
    self.assert_pdb_hierarchy_xray_structure_sync()
    from mmtbx.maps import correlation
    current_map = self.compute_map(xray_structure = self.xray_structure)
    return correlation.histogram_per_atom(
      map_1      = current_map,
      map_2      = self.target_map_object.map_data,
      sites_cart = self.xray_structure.sites_cart(),
      unit_cell  = self.xray_structure.unit_cell(),
      radius     = radius,
      n_slots    = n_slots)

  def map_cc(self, other_map, sites_cart=None, atom_radius=2, per_atom=False):
    self.assert_pdb_hierarchy_xray_structure_sync()
    from mmtbx.maps import correlation
    if(sites_cart is not None):
      if(per_atom):
        result = correlation.from_map_map_atoms_per_atom(
          map_1      = other_map,
          map_2      = self.target_map_object.map_data,
          sites_cart = sites_cart,
          unit_cell  = self.xray_structure.unit_cell(),
          radius     = atom_radius)
      else:
        result = correlation.from_map_map_atoms(
          map_1      = other_map,
          map_2      = self.target_map_object.map_data,
          sites_cart = sites_cart,
          unit_cell  = self.xray_structure.unit_cell(),
          radius     = atom_radius)
    else:
      result = correlation.from_map_map(
        map_1 = other_map,
        map_2 = self.target_map_object.map_data)
    return result

  def show(self, prefix="", log=None):
    self.assert_pdb_hierarchy_xray_structure_sync()
    if(log is None): log = sys.stdout
    fmt = """%s CC_mask:                   %-6.3f
%s CC_volume:                 %-6.3f
%s CC_peaks:                  %-6.3f
%s rmsd (bonds):              %-s
%s rmsd (angles):             %-s
%s Dist. moved from start:    %-6.3f
%s Dist. moved from previous: %-6.3f
%s All-atom clashscore        %-s
%s Ramachandran plot:
%s   outliers:                %-s %%
%s   allowed:                 %-s %%
%s   favored:                 %-s %%
%s Omega angle:
%s   cis-proline:             %-s %%
%s   twisted proline:         %-s %%
%s   cis-general:             %-s %%
%s   twisted-general:         %-s %%
%s CaBLAM analysis:
%s   outliers:                %-s %%
%s   disfavored:              %-s %%
%s   ca outliers:             %-s %%
%s Rotamer outliers:          %-s %%
%s C-beta deviations:         %-s %%
"""
    mso = None
    try:
      if self.geometry_restraints_manager is not None and False:
        # XXX False at the end is intentional, because currently I want to
        # disable this 'if' branch. Reason is - nothing from extended
        # model_statistics (with GRM) is being used, so no reason to spend
        # time calculating statistics over various restraints.
        mso = mmtbx.model.statistics.geometry(
          pdb_hierarchy      = self.pdb_hierarchy,
          molprobity_scores  = libtbx.env.has_module("probe"),
          restraints_manager = self.geometry_restraints_manager)
      else:
        mso = mmtbx.model.statistics.geometry_no_grm(
          pdb_hierarchy      = self.pdb_hierarchy,
          molprobity_scores  = libtbx.env.has_module("probe"))
    except Exception:
      # some part of validation failed
      pass
    self.stats_evaluations.append(
        group_args(
          cc = group_args(
              cc_mask   = self.five_cc.result.cc_mask,
              cc_volume = self.five_cc.result.cc_volume,
              cc_peaks  = self.five_cc.result.cc_peaks),
          geometry = mso,
          rmsd_a = self.rmsd_a,
          rmsd_b = self.rmsd_b))
    if mso is not None and self.five_cc is not None:
      print(fmt%(
        # prefix, self.map_cc_whole_unit_cell,
        # prefix, self.map_cc_around_atoms,
        prefix, self.five_cc.cc_mask,
        prefix, self.five_cc.cc_volume,
        prefix, self.five_cc.cc_peaks,
        prefix, format_value("%-6.2f", self.rmsd_b).strip(),
        prefix, format_value("%-6.2f", self.rmsd_a).strip(),
        prefix, self.dist_from_start,
        prefix, self.dist_from_previous,
        prefix, format_value("%-6.2f", mso.clashscore),
        prefix,
        prefix, format_value("%-5.2f", mso.ramachandran_outliers),
        prefix, format_value("%-5.2f", mso.ramachandran_allowed),
        prefix, format_value("%-5.2f", mso.ramachandran_favored),
        prefix,
        prefix, format_value("%-5.2f", mso.cis_proline),
        prefix, format_value("%-5.2f", mso.twisted_proline),
        prefix, format_value("%-5.2f", mso.cis_general),
        prefix, format_value("%-5.2f", mso.twisted_general),
        prefix,
        prefix, format_value("%-5.2f", mso.cablam_outliers),
        prefix, format_value("%-5.2f", mso.cablam_disfavored),
        prefix, format_value("%-5.2f", mso.cablam_ca_outliers),
        prefix, format_value("%6.2f", mso.rotamer_outliers).strip(),
        prefix, format_value("%-5.2f", mso.c_beta_dev_percent)), file=log)

  def show_residues(self, map_cc_all=0.8, map_cc_sidechain=0.8, log=None):
    self.assert_pdb_hierarchy_xray_structure_sync()
    if(log is None): log = sys.stdout
    header_printed = True
    for r in self.residue_monitors:
      i1=r.map_cc_all < map_cc_all
      i2=r.rotamer_status == "OUTLIER"
      i4=r.map_cc_sidechain is not None and r.map_cc_sidechain<map_cc_sidechain
      if([i1,i2,i4].count(True)>0):
        if(header_printed):
          print("Residue     CC        CC         CC   Rotamer", file=log)
          print("     id    all  backbone  sidechain        id", file=log)
          header_printed = False
        print(r.format_info_string(), file=log)

  def update(self, xray_structure, accept_as_is=True):
    if(not accept_as_is):
      current_map = self.compute_map(xray_structure = xray_structure)
      sites_cart  = xray_structure.sites_cart()
      sites_cart_ = self.xray_structure.sites_cart()
      for r in self.residue_monitors:
        sca = sites_cart.select(r.selection_all)
        scs = sites_cart.select(r.selection_sidechain)
        scb = sites_cart.select(r.selection_backbone)
        map_cc_all       = self.map_cc(sites_cart = sca, other_map = current_map)
        map_cc_sidechain = self.map_cc(sites_cart = scs, other_map = current_map)
        map_cc_backbone  = self.map_cc(sites_cart = scb, other_map = current_map)
        flag = map_cc_all      >= r.map_cc_all and \
               map_cc_backbone >= r.map_cc_backbone and \
               map_cc_sidechain>= r.map_cc_sidechain
        if(flag):
          residue_sites_cart_new = sites_cart.select(r.selection_all)
          sites_cart_ = sites_cart_.set_selected(r.selection_all,
            residue_sites_cart_new)
      xray_structure = xray_structure.replace_sites_cart(sites_cart_)
    # re-initialize monitor
    self.dist_from_previous = flex.mean(self.xray_structure.distances(
      other = xray_structure))
    self.xray_structure = xray_structure
    self.pdb_hierarchy.adopt_xray_structure(xray_structure)
    self.initialize()
    self.states_collector.add(sites_cart = xray_structure.sites_cart())
    self.assert_pdb_hierarchy_xray_structure_sync()

def selection_around_to_negate(
      xray_structure,
      selection_within_radius,
      iselection,
      selection_good=None,
      iselection_backbone=None,
      iselection_n_external=None,
      iselection_c_external=None):
  # takes ~0.002 seconds
  if([selection_good,iselection_backbone].count(None)==0):
    selection_backbone = flex.bool(selection_good.size(), iselection_backbone)
    selection_good = selection_good.set_selected(selection_backbone, True)
  sel_around = xray_structure.selection_within(
    radius    = selection_within_radius,
    selection = flex.bool(xray_structure.scatterers().size(), iselection))
  if(selection_good is not None):
    ssb = flex.bool(selection_good.size(), iselection)
    sel_around_minus_self = sel_around.set_selected(ssb, False)
  else:
    sel_around_minus_self = flex.size_t(tuple(
      set(sel_around.iselection()).difference(set(iselection))))
  if(selection_good is not None):
    negate_selection = sel_around_minus_self & selection_good
  else:
    negate_selection = sel_around_minus_self
  if(iselection_n_external is not None and iselection_n_external.size()>0):
    negate_selection[iselection_n_external[0]]=False
  if(iselection_c_external is not None and iselection_c_external.size()>0):
    negate_selection[iselection_c_external[0]]=False
  return negate_selection

def negate_map_around_selected_atoms_except_selected_atoms(
      xray_structure,
      map_data,
      negate_selection,
      atom_radius):
  # XXX time and memory inefficient
  sites_cart_p1 = xray_structure.select(negate_selection).expand_to_p1(
    sites_mod_positive=True).sites_cart()
  around_atoms_selections = maptbx.grid_indices_around_sites(
    unit_cell  = xray_structure.unit_cell(),
    fft_n_real = map_data.focus(),
    fft_m_real = map_data.all(),
    sites_cart = sites_cart_p1,
    site_radii = flex.double(sites_cart_p1.size(), atom_radius))
  return maptbx.negate_selected_in_place(map_data=map_data,
    selection=around_atoms_selections)

class score2(object):
  def __init__(self,
               unit_cell,
               target_map,
               residue,
               vector = None,
               selection=None):
    adopt_init_args(self, locals())
    self.target = None
    self.sites_cart = None
    self.i_seqs = []
    self.weights = flex.double()
    for el in self.residue.atoms().extract_element():
      std_lbl = eltbx.xray_scattering.get_standard_label(
        label=el, exact=True, optional=True)
      self.weights.append(tiny_pse.table(std_lbl).weight())
    self.occ = self.residue.atoms().extract_occ()
    self.vector_flat = None
    if(vector is not None):
      self.vector_flat = flex.size_t(flatten(self.vector))
    self.sites_cart = self.residue.atoms().extract_xyz()
    if(selection is None): selection = self.vector_flat
    self.target = maptbx.real_space_target_simple(
      unit_cell   = self.unit_cell,
      density_map = self.target_map,
      sites_cart  = self.sites_cart,
      selection   = selection)

  def update(self, sites_cart, selection=None):
    if(selection is None): selection = self.vector_flat
    target = maptbx.real_space_target_simple(
      unit_cell   = self.unit_cell,
      density_map = self.target_map,
      sites_cart  = sites_cart,
      selection   = selection)
    if(target > self.target):
      self.sites_cart = sites_cart
      self.target = target

class score(object):
  def __init__(self,
               unit_cell,
               target_map,
               residue,
               rotamer_eval = None,
               vector = None):
    adopt_init_args(self, locals())
    self.target = None
    self.sites_cart = None
    self.i_seqs = []
    self.weights = flex.double()
    for el in self.residue.atoms().extract_element():
      std_lbl = eltbx.xray_scattering.get_standard_label(
        label=el, exact=True, optional=True)
      self.weights.append(tiny_pse.table(std_lbl).weight())
    self.occ = self.residue.atoms().extract_occ()
    self.vector_flat = flatten(self.vector)

  def compute_target(self, sites_cart, selection=None):
    sites_frac = self.unit_cell.fractionalize(sites_cart)
    result = 0
    vals = []
    if(selection is None): i_seqs = self.vector_flat
    else:                  i_seqs = selection
    for i_seq in i_seqs:
      vals.append(self.target_map.eight_point_interpolation(sites_frac[i_seq])/
        self.weights[i_seq]/self.occ[i_seq])
    #
    sz = len(vals)
    if(sz>3):
      deltas = []
      for i in xrange(sz):
        if(i+1<sz and i>1):
          e1=abs(vals[i])
          e2=abs(vals[i+1])
          r=e1/e2
          deltas.append(r)
      if(max(deltas)>5 or min(deltas)<1./5):
        return 0
    return sum(vals)

  def update(self, sites_cart, selection=None, tmp=None):
    target = self.compute_target(sites_cart = sites_cart, selection=selection)
    assert self.target is not None
    if(target > self.target):
      self.residue.atoms().set_xyz(sites_cart)
      fl = self.rotamer_eval is None or \
        self.rotamer_eval.evaluate_residue(residue = self.residue) != "OUTLIER"
      if(fl):
        self.target = target
        self.sites_cart = sites_cart

  def reset(self, sites_cart, selection=None):
    self.target = self.compute_target(sites_cart = sites_cart,
      selection = selection)
    self.sites_cart = sites_cart

def torsion_search(clusters, scorer, sites_cart, start, stop, step):
  def generate_range(start, stop, step):
    assert abs(start) <= abs(stop)
    inc = start
    result = []
    while abs(inc) <= abs(stop):
      result.append(inc)
      inc += step
    return result
  for i_cl, cl in enumerate(clusters):
    if(i_cl == 0): sites_cart_start = sites_cart.deep_copy()
    else:          sites_cart_start = scorer.sites_cart.deep_copy()
    scorer.reset(sites_cart=sites_cart_start, selection=cl.selection)
    sites_cart_ = scorer.sites_cart.deep_copy()
    for angle_deg in generate_range(start=start, stop=stop, step=step):
      xyz_moved = sites_cart_.deep_copy()
      for atom in cl.atoms_to_rotate:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = sites_cart_[cl.axis[0]],
          axis_point_2 = sites_cart_[cl.axis[1]],
          point        = sites_cart_[atom],
          angle        = angle_deg, deg=True)
        xyz_moved[atom] = new_xyz
      scorer.update(sites_cart = xyz_moved, selection = cl.selection)
  return scorer

def torsion_search_nested(
      clusters,
      scorer,
      sites_cart):
  n_angles = len(clusters)
  print(n_angles)
  if(n_angles == 3):
    r1 = [-3,-7,-9]
    r2 = [3,7,9]
  elif(n_angles == 4):
    r1 = [-5,-5,-10,-10]
    r2 = [5,5,10,10]
  else: return
  nested_loop = flex.nested_loop(begin=r1, end=r2, open_range=False)
  selection = clusters[0].atoms_to_rotate
  scorer.reset(sites_cart = sites_cart, selection = selection)
  for angles in nested_loop:
    xyz_moved = sites_cart.deep_copy()
    for i, angle in enumerate(angles):
      cl = clusters[i]
      for atom in cl.atoms_to_rotate:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = xyz_moved[cl.axis[0]],
          axis_point_2 = xyz_moved[cl.axis[1]],
          point        = xyz_moved[atom],
          angle        = angle, deg=True)
        xyz_moved[atom] = new_xyz
    scorer.update(sites_cart = xyz_moved, selection = selection)
  return scorer

class score3(object):
  def __init__(self,
               unit_cell,
               target_map,
               residue,
               rotamer_eval):
    adopt_init_args(self, locals())
    self.target = None
    self.sites_cart = self.residue.atoms().extract_xyz()
    self.status = self.rotamer_eval.evaluate_residue(residue = self.residue)

  def compute_target(self, sites_cart, selection=None):
    if(selection is not None):
      return maptbx.real_space_target_simple(
        unit_cell   = self.unit_cell,
        density_map = self.target_map,
        sites_cart  = sites_cart,
        selection   = selection)
    else:
      return maptbx.real_space_target_simple(
        unit_cell   = self.unit_cell,
        density_map = self.target_map,
        sites_cart  = sites_cart)

  def update(self, sites_cart, selection=None):
    target = self.compute_target(sites_cart=sites_cart, selection=selection)
    assert self.target is not None
    den = (abs(self.target)+abs(target))*100
    num = abs(abs(self.target)-abs(target))*2
    if(den==0): second_cond = False
    else:       second_cond = num/den<5.
    if(target > self.target):
      self.residue.atoms().set_xyz(sites_cart)
      fl = self.rotamer_eval is None or \
        self.rotamer_eval.evaluate_residue(residue = self.residue) != "OUTLIER"
      if(fl):
        self.target = target
        self.sites_cart = sites_cart
    elif(second_cond):
      fl = self.rotamer_eval is None or \
        self.rotamer_eval.evaluate_residue(residue = self.residue) == "OUTLIER"
      if(fl):
        self.residue.atoms().set_xyz(sites_cart)
        fl = self.rotamer_eval is None or \
          self.rotamer_eval.evaluate_residue(residue = self.residue) != "OUTLIER"
        if(fl):
          self.target = target
          self.sites_cart = sites_cart

  def reset(self, sites_cart, selection=None):
    self.target = self.compute_target(sites_cart = sites_cart,
      selection = selection)
    self.sites_cart = sites_cart
