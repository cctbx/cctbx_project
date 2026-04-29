from __future__ import absolute_import, division, print_function
from mmtbx.utils import rotatable_bonds
from scitbx.matrix import rotate_point_around_axis
from cctbx.array_family import flex
from cctbx import maptbx
import mmtbx.model
import iotbx.pdb
from libtbx.utils import Sorry

import boost_adaptbx.boost.python as bp
from six.moves import range

hydrogens_master_params_str = """
refine = individual riding *Auto
  .type = choice
  .help = Choice for refinement: riding model or full (H is refined as \
          other atoms, useful at very high resolutions only)
  .short_caption = Hydrogen refinement model
  .expert_level=1
force_riding_adp = None
  .type = bool
optimize_scattering_contribution = True
  .type = bool
contribute_to_f_calc = True
  .type = bool
  .help = Add H contribution to Xray (Fcalc) calculations
  .short_caption=Include hydrogens in Fcalc
  .expert_level=1
high_resolution_limit_to_include_scattering_from_h = 1.6
  .type = float
  .short_caption = High-resolution limit to include scattering from H
  .expert_level=2
real_space_optimize_x_h_orientation = True
  .type = bool
  .short_caption = Optimize X-H orientation in real-space
  .expert_level = 1
xh_bond_distance_deviation_limit = 0.0
  .type = float
  .help = Idealize XH bond distances if deviation from ideal is greater \
          than xh_bond_distance_deviation_limit
  .short_caption = X-H bond distance deviation limit
  .expert_level=2
"""

def shortcut(residue, first, log):
  """
  Avoid TARDY, huge speedup for huge structures!
  """
  d = {}
  for a in residue.atoms():
    name = a.name.strip().upper()
    if(name[0]=="D"):
      name = list(name)
      name[0]="H"
      name = "".join(name)
    d[name] = a.i_seq
  rn = residue.resname.strip().upper()
  result = []
  try:
    if(first):
      try: # H1,H2,H3 are missing way to often, also due to reduce bug
        axis  = [d["CA"], d["N"]]
        atoms = [d["H1"], d["H2"], d["H3"]]
        result.append([axis, atoms])
      except KeyError: pass
    if  (rn == "MET"):
      axis  = [d['SD'],  d['CE']           ]
      atoms = [d['HE1'], d['HE2'], d['HE3']]
      result.append([axis, atoms])
    elif(rn == "LYS"):
      axis  = [d['CE'],  d['NZ']           ]
      atoms = [d['HZ1'], d['HZ2'], d['HZ3']]
      result.append([axis, atoms])
    elif(rn == "SER"):
      axis  = [d['CB'], d['OG']]
      atoms = [d['HG']]
      result.append([axis, atoms])
    elif(rn == "THR"):
      axis  = [d['CB'],   d['CG2']            ]
      atoms = [d['HG21'], d['HG22'], d['HG23']]
      result.append([axis, atoms])
      axis  = [d['CB'], d['OG1']]
      atoms = [d['HG1']         ]
      result.append([axis, atoms])
    elif(rn == "LEU"):
      axis  = [d['CG'],   d['CD1']            ]
      atoms = [d['HD11'], d['HD12'], d['HD13']]
      result.append([axis, atoms])
      axis  = [d['CG'],   d['CD2']            ]
      atoms = [d['HD21'], d['HD22'], d['HD23']]
      result.append([axis, atoms])
    elif(rn == "CYS"):
      if("HG" in d.keys()): # not a disulfide bridge
        axis  = [d['CB'], d['SG']]
        atoms = [d['HG']         ]
        result.append([axis, atoms])
    elif(rn == "VAL"):
      axis  = [d['CB'],   d['CG1']            ]
      atoms = [d['HG11'], d['HG12'], d['HG13']]
      result.append([axis, atoms])
      axis  = [d['CB'],   d['CG2']            ]
      atoms = [d['HG21'], d['HG22'], d['HG23']]
      result.append([axis, atoms])
    elif(rn == "TYR"):
      axis  = [d['CZ'], d['OH']]
      atoms = [d['HH']         ]
      result.append([axis, atoms])
    elif(rn == "ALA"):
      axis  = [d['CA'],  d['CB']           ]
      atoms = [d['HB1'], d['HB2'], d['HB3']]
      result.append([axis, atoms])
    elif(rn == "ILE"):
      axis  = [d['CB'],   d['CG2']            ]
      atoms = [d['HG21'], d['HG22'], d['HG23']]
      result.append([axis, atoms])
      axis  = [d['CG1'],  d['CD1']            ]
      atoms = [d['HD11'], d['HD12'], d['HD13']]
      result.append([axis, atoms])
    elif(rn == "MSE"):
      axis  =  [d['SE'],  d['CE']           ]
      atoms =  [d['HE1'], d['HE2'], d['HE3']]
      result.append([axis, atoms])
    else:
      if(len(result)>0): return result
      else:              return None
  except KeyError:
    m="Residue %s %s is missing expected H atoms. Skipping."%(
      residue.resname, str(residue.resseq))
    print(m, file=log)
    return None
  return result

def rotatable(pdb_hierarchy, mon_lib_srv, restraints_manager, log,
              use_shortcut=True):
  """
  General tool to identify rotatable H, such as C-O-H, C-H3, in any molecule.
  """
  result = []
  def analyze_group_aa_specific(g, atoms, psel):
    result = []
    for gi in g:
      assert len(gi[0])==2 # because this is axis
      assert len(gi[1])>0  # because these are atoms rotating about this axis
      # condition 1: axis does not contain H or D
      a1, a2 = atoms[gi[0][0]], atoms[gi[0][1]]
      e1 = a1.element.strip().upper()
      e2 = a2.element.strip().upper()
      #
      condition_00_ = psel[a1.i_seq] and psel[a2.i_seq]
      condition_00 = flex.bool([condition_00_])
      for gi1_ in gi[1]:
        condition_00.append(condition_00_ and psel[atoms[gi1_].i_seq])
      condition_00 = condition_00.all_eq(True)
      if condition_00: continue
      #
      condition_1 = [e1,e2].count("H")==0 and [e1,e2].count("D")==0
      # condition 2: all atoms to rotate are H or D
      condition_2 = True
      rot_atoms = []
      for gi1i in gi[1]:
        if(not atoms[gi1i].element.strip().upper() in ["H","D"]):
          condition_2 = False
          break
      rot_atoms = []
      axis = None
      if(condition_1 and condition_2 and not condition_00):
        axis = [a1.i_seq, a2.i_seq]
        for gi1i in gi[1]:
          rot_atoms.append(atoms[gi1i].i_seq)
        result.append([axis, rot_atoms])
    if(len(result)>0): return result
    else: return None
  def analyze_group_general(g, atoms, bps, psel):
    result = []
    for gi in g:
      condition_1, condition_2, condition_3 = None,None,None
      assert len(gi[0])==2 # because this is axis
      assert len(gi[1])>0  # because these are atoms rotating about this axis
      # condition 1: axis does not contain H or D
      a1, a2 = atoms[gi[0][0]], atoms[gi[0][1]]
      e1 = a1.element.strip().upper()
      e2 = a2.element.strip().upper()
      #
      condition_00_ = psel[a1.i_seq] and psel[a2.i_seq]
      condition_00 = flex.bool([condition_00_])
      for gi1_ in gi[1]:
        condition_00.append(condition_00_ and psel[atoms[gi1_].i_seq])
      condition_00 = condition_00.all_eq(True)
      if condition_00: continue
      #
      condition_1 = [e1,e2].count("H")==0 and [e1,e2].count("D")==0
      s1 = set(gi[1])
      if(condition_1):
        # condition 2: all atoms to rotate are H or D
        condition_2 = True
        for gi1i in gi[1]:
          if(not atoms[gi1i].element.strip().upper() in ["H","D"]):
            condition_2 = False
            break
        if(condition_2):
          # condition 3: one of axis atoms is terminal (bonded to another axis
          #              atom and hydrogens
          condition_3 = False
          for gia in gi[0]:
            bonds_involved_into = []
            for bp in bps:
              if(gia in bp.i_seqs):
                for i_seq in bp.i_seqs:
                  if(atoms[i_seq].element.strip().upper() in ["H","D"]):
                    bonds_involved_into.append(i_seq)
            s2 = set(bonds_involved_into)
            s = list(s1 & s2)
            if(len(s)>0): condition_3 = True
          #
          if(condition_1 and condition_2 and condition_3):
          #if(condition_1 and condition_2 and condition_3 and not condition_00):
            axis = [a1.i_seq, a2.i_seq]
            rot_atoms = []
            in_plane = False
            for i in bonds_involved_into:
              if(psel[atoms[i].i_seq]): in_plane = True
              rot_atoms.append(atoms[i].i_seq)
            if(not in_plane): result.append([axis, rot_atoms])
    if(len(result)>0): return result
    else: return None
  def helper_1(residue, mon_lib_srv, log, result, psel):
    fr = rotatable_bonds.axes_and_atoms_aa_specific(
      residue=residue, mon_lib_srv=mon_lib_srv,
      remove_clusters_with_all_h=False, log=log)
    if(fr is not None):
      r = analyze_group_aa_specific(g=fr, atoms=residue.atoms(), psel=psel)
      if(r is not None):
        for r_ in r: result.append(r_)
  def helper_2(atoms, restraints_manager):
    elements = atoms.extract_element()
    names = atoms.extract_name()
    # create tardy_model
    iselection = atoms.extract_i_seq()
    sites_cart = atoms.extract_xyz()
    masses     = [1]*sites_cart.size()
    labels     = list(range(sites_cart.size()))
    grm_i = restraints_manager.select(iselection)
    bps, asu = grm_i.geometry.get_all_bond_proxies(
      sites_cart = sites_cart)
    edge_list = []
    for bp in bps: edge_list.append(bp.i_seqs)
    fixed_vertex_lists = []
    tmp_r = []
    # try all possible edges (bonds) as potential fixed vertices and
    # accept only non-redundant
    for bp in bps:
      tardy_tree = scitbx.graph.tardy_tree.construct(
        sites              = sites_cart,
        edge_list          = edge_list,
        fixed_vertex_lists = [bp.i_seqs])
      tardy_model = scitbx.rigid_body.tardy_model(
        labels        = labels,
        sites         = sites_cart,
        masses        = masses,
        tardy_tree    = tardy_tree,
        potential_obj = None)
      fr = rotatable_bonds.axes_and_atoms_aa_specific(
        residue=residue, mon_lib_srv=mon_lib_srv,
        remove_clusters_with_all_h=False, log=None,
        tardy_model = tardy_model)
      if(fr is not None):
        r = analyze_group_general(g=fr, atoms=atoms,bps=bps,psel=psel)
        if(r is not None and len(r)>0):
          for r_ in r:
            if(not r_ in tmp_r):
              if(not r_ in tmp_r): tmp_r.append(r_)
    for r in tmp_r:
      if(not r in result):
        result.append(r)
  #
  if(restraints_manager is not None):
    psel = flex.bool(pdb_hierarchy.atoms().size(), False)
    for p in restraints_manager.geometry.planarity_proxies:
      for i in p.i_seqs:
        psel[i] = True
  # very handy for debugging: do not remove
  #NAMES = pdb_hierarchy.atoms().extract_name()
  #
  get_class = iotbx.pdb.common_residue_names_get_class
  import scitbx.graph.tardy_tree
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      residue_groups = chain.residue_groups()
      n_residues = len(residue_groups)
      for i_rg, residue_group in enumerate(residue_groups):
        first = i_rg == 0
        last  = i_rg+1 == n_residues
        first_or_last = first or last
        conformers = residue_group.conformers()
        for conformer in residue_group.conformers():
          for residue in conformer.residues():
            if(residue.resname.strip().upper() == "PRO"): continue
            atoms = residue.atoms()
            if(get_class(name=residue.resname)=="common_water" and
               len(atoms)==1):
                 continue
            hd = False
            for a in residue.atoms():
              if(a.element_is_hydrogen()):
                hd=True
                break
            if(not hd): continue
            if(get_class(name=residue.resname)=="common_amino_acid" and not last):
              if(use_shortcut):
                r_ = shortcut(residue=residue, first=first, log=log)
                if(r_ is not None): result.extend(r_)
              else:
                helper_1(residue, mon_lib_srv, log, result, psel)
                if(first):
                  helper_2(atoms, restraints_manager)
            elif(get_class(name=residue.resname)=="common_amino_acid"):
              helper_1(residue, mon_lib_srv, log, result, psel)
              helper_2(atoms, restraints_manager)
            elif((restraints_manager is not None)):
              helper_2(atoms, restraints_manager)
  # very handy for debugging: do not remove
  #for r_ in result:
  #  print "  analyze_group:", r_, \
  #    [NAMES[i] for i in r_[0]], [NAMES[i] for i in r_[1]], residue.resname
  return result

def count_rotatable(selections):
  result = 0
  for s in selections:
    result += len(s[1])
  return result

class map_manager(object):

  def __init__(self, fmodel, map_type):
    self.fmodel = fmodel
    self.map_type = map_type
    cs = fmodel.f_obs().crystal_symmetry()
    self.unit_cell = cs.unit_cell()
    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell        = self.unit_cell,
      space_group_info = cs.space_group_info(),
      symmetry_flags   = maptbx.use_space_group_symmetry,
      step             = 0.25)
    self.omit_map = None
    self.size = self.fmodel.xray_structure.scatterers().size()

  def update_omit_map(self, omit_selection):
    fmodel = self.fmodel.deep_copy()
    xrs = fmodel.xray_structure
    xrs.set_occupancies(value=0, selection = omit_selection)
    fmodel.update_xray_structure(
      xray_structure = xrs,
      update_f_calc  = True,
      update_f_mask  = False)
    mc = fmodel.electron_density_map().map_coefficients(
      map_type   = "mFobs-DFmodel",
      isotropize = False,
      exclude_free_r_reflections = True)
    fft_map = mc.fft_map(crystal_gridding = self.crystal_gridding)
    fft_map.apply_sigma_scaling()
    #fft_map.apply_volume_scaling()
    self.omit_map = fft_map.real_map_unpadded()

  def score(self, sites_cart):
    return maptbx.real_space_target_simple(
      unit_cell   = self.unit_cell,
      density_map = self.omit_map,
      sites_cart  = sites_cart)

def map_statistics(model, fmodel):
  result = flex.double()
  mm = map_manager(fmodel = fmodel, map_type = "mFobs-DFmodel")
  for m in model.get_hierarchy().models():
    for c in m.chains():
      for rg in c.residue_groups():
        h_sel      = flex.size_t()
        sites_cart = flex.vec3_double()
        for atom in rg.atoms():
          if not atom.element_is_hydrogen(): continue
          h_sel.append(atom.i_seq)
          sites_cart.append(atom.xyz)
        if h_sel.size()==0: continue
        mm.update_omit_map(omit_selection = h_sel)
        for site_cart in sites_cart:
          score = mm.score(sites_cart = flex.vec3_double([site_cart]))
          score = min(3.5, max(0, score))
          result.append(score)
  return result

def fit_rotatable2(model, fmodel):
  """
  Slow, OMIT map based fitting of rotatable H.
  """
  # Reset mask params if needed
  mask_params = fmodel.mask_params
  mpih = mask_params.ignore_hydrogens
  if(mpih):
    mask_params.ignore_hydrogens=False
    fmodel.update(mask_params=mask_params)
  # Set X-H lengths based on data type
  scattering_table = model.get_scattering_table()
  if scattering_table is None:
    raise Sorry("scattering_table must be set.")
  use_neutron_distances = False
  if scattering_table in ["neutron", "electron"]:
    use_neutron_distances = True
  model.set_hydrogen_bond_length(use_neutron_distances = use_neutron_distances)
  model.idealize_h_riding()
  #
  fmodel.update_xray_structure(
    xray_structure = model.get_xray_structure(),
    update_f_calc  = True,
    update_f_mask  = True)
  fmodel.update_all_scales()
  #
  # FIT
  mm = map_manager(fmodel = fmodel, map_type = "mFobs-DFmodel")
  get_class = iotbx.pdb.common_residue_names_get_class
  sites_cart = model.get_sites_cart()

  #for m in model.get_hierarchy().models():
  #  for c in m.chains():
  #    first=True
  #    for r in c.residues():

  for m in model.get_hierarchy().models():
    for c in m.chains():
      first=True
      for residue_group in c.residue_groups():
        conformers = residue_group.conformers()
        for conformer in conformers:
          r = conformer.only_residue()

          if not get_class(r.resname)=="common_amino_acid": continue
          s = mmtbx.hydrogens.shortcut(residue=r, first=first, log=model.log)
          first=False
          if s is None: continue
          # print(r.resname, s)
          omit_selection = flex.size_t()
          for cl in s:
            omit_selection.extend(flex.size_t(cl[1]))
          mm.update_omit_map(omit_selection = omit_selection)
          for cl in s:
            axis = cl[0]
            a1, a2               = sites_cart[axis[0]], sites_cart[axis[1]]
            sel_to_rotate        = flex.size_t(cl[1])
            sites_cart_to_rotate = sites_cart.select(sel_to_rotate)
            score_start          = mm.score(sites_cart = sites_cart_to_rotate)
            score_best           = score_start
            sites_cart_moved     = flex.vec3_double(sites_cart_to_rotate.size())
            sites_cart_best      = None
            assert len(sel_to_rotate) in [1,3]
            if len(sel_to_rotate)==1: stop=360
            else:                     stop=60
            for angle in range(0, stop, 1):
              sites_cart_moved = flex.vec3_double(sites_cart_to_rotate.size())
              for isite, site_cart in enumerate(sites_cart_to_rotate):
                site_cart_rotated = rotate_point_around_axis(
                  axis_point_1 = a1,
                  axis_point_2 = a2,
                  point        = site_cart,
                  angle        = angle,
                  deg          = True)
                sites_cart_moved[isite] = site_cart_rotated
              score = mm.score(sites_cart = sites_cart_moved)
              if score > score_best:
                score_best = score
                sites_cart_best = sites_cart_moved.deep_copy()
            if sites_cart_best is not None:
              sites_cart = sites_cart.set_selected(sel_to_rotate, sites_cart_best)
  model.set_sites_cart(sites_cart)
  fmodel.xray_structure.set_sites_cart(model.get_sites_cart())
  if(mpih):
    mask_params.ignore_hydrogens=False
    fmodel.update(mask_params=mask_params)
  fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)
  fmodel.update_all_scales()

def fit_rotatable(
      pdb_hierarchy,
      xray_structure,
      map_data,
      rotatable_h_selection):
  unit_cell = xray_structure.unit_cell()
  sites_cart = xray_structure.sites_cart()
  scatterers = xray_structure.scatterers()
  for sel_ in rotatable_h_selection:
    ed_val = -1
    angle = 0.
    angular_step = 1
    axis = sel_[0]
    points_i_seqs = sel_[1]
    sites_frac_best = flex.vec3_double(len(points_i_seqs))
    while angle <= 360:
      sites_frac_tmp  = flex.vec3_double(len(points_i_seqs))
      ed_val_ = 0
      for i_seq, point_i_seq in enumerate(points_i_seqs):
        site_cart_new = rotate_point_around_axis(
          axis_point_1 = sites_cart[axis[0]],
          axis_point_2 = sites_cart[axis[1]],
          point        = sites_cart[point_i_seq],
          angle        = angle,
          deg          = True)
        site_frac_new = unit_cell.fractionalize(site_cart_new)
        ed_val_ += abs(maptbx.eight_point_interpolation(map_data,site_frac_new))
        sites_frac_tmp[i_seq] = site_frac_new
      if(ed_val_ > ed_val):
        ed_val = ed_val_
        sites_frac_best = sites_frac_tmp.deep_copy()
      angle += angular_step
    for i_seq, point_i_seq in enumerate(points_i_seqs):
      scatterers[point_i_seq].site = sites_frac_best[i_seq]
  pdb_hierarchy.adopt_xray_structure(xray_structure)

def run_fit_rotatable(
      fmodel,
      ref_model,
      angular_step,
      log = None,
      use_h_omit_map = False,
      map_type="2mFo-DFc"):
  pdb_hierarchy = ref_model.get_hierarchy()
  xrs = fmodel.xray_structure
  rotatable_h_selection = rotatable(
    pdb_hierarchy      = pdb_hierarchy,
    mon_lib_srv        = ref_model.get_mon_lib_srv(),
    restraints_manager = ref_model.get_restraints_manager(),
    log                = log)
  rotatable_h_selection_1d = []
  for s in rotatable_h_selection:
    rotatable_h_selection_1d.extend(s[1])
  rotatable_h_selection_1d = flex.size_t(rotatable_h_selection_1d)
  rotatable_h_selection_1d_bool = flex.bool(xrs.scatterers().size(),
    rotatable_h_selection_1d)
  if(log is not None):
    print("Real-space grid search fit H (or D) atoms:", file=log)
    print("  start:  r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(),
      fmodel.r_free()), file=log)
  if(use_h_omit_map):
    xrs_omit = fmodel.xray_structure.select(~rotatable_h_selection_1d_bool)
    fmodel.update_xray_structure(xray_structure = xrs_omit, update_f_calc= True)
    if(log is not None):
      print("  H omit: r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(),
        fmodel.r_free()), file=log)
  fft_map = fmodel.electron_density_map().fft_map(
    resolution_factor = 1./4.,
    map_type          = map_type,
    symmetry_flags    = maptbx.use_space_group_symmetry)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  fit_rotatable(
    pdb_hierarchy  = pdb_hierarchy,
    xray_structure = xrs,
    rotatable_h_selection=rotatable_h_selection,
    map_data       = map_data)
  fmodel.update_xray_structure(xray_structure = xrs, update_f_calc=True)
  ref_model.xray_structure = xrs
  if(log is not None):
    print("  final:  r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(),
      fmodel.r_free()), file=log)
