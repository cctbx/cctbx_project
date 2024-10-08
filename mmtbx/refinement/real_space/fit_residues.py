from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx import adopt_init_args
import iotbx.pdb
import mmtbx.refinement.real_space
import mmtbx.refinement.real_space.fit_residue
import sys, time
from cctbx import crystal, maptbx
from mmtbx.maps import correlation
import boost_adaptbx.boost.python as bp
cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")
fit_ext = bp.import_ext("mmtbx_rotamer_fit_ext")

import scitbx.math
from mmtbx.rotamer.rotamer_eval import RotamerEval
rotamer_manager = RotamerEval()

# XXX Keep for debugging
def write_map_file(crystal_symmetry, map_data, file_name):
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
    file_name   = file_name,
    unit_cell   = crystal_symmetry.unit_cell(),
    space_group = crystal_symmetry.space_group(),
    map_data    = map_data,
    labels      = flex.std_string([""]))

negate_map_table = {
  #"ala": False,
  "asn": 5,
  "asp": 5,
  "cys": 5,
  "gln": 6,
  "glu": 6,
  #"gly": False,
  "his": 6,
  "ile": 5,
  "leu": 5,
  "met": 6,
  "mse": 6,
  "phe": 6.5,
  #"pro": 5,#False,
  "ser": 5,
  "thr": 5,
  "trp": 7.4,
  "tyr": 7.7,
  "val": 5,
  "arg": 10,
  "lys": 8,
  "pro": 3
}

def get_map(fmodel, crystal_gridding, map_type):
  coeffs = fmodel.electron_density_map().map_coefficients(
    map_type     = map_type,
    fill_missing = False,
    isotropize   = False)
  fft_map = coeffs.fft_map(crystal_gridding = crystal_gridding)
  fft_map.apply_sigma_scaling()
  result = fft_map.real_map_unpadded()
  return result

def fit_altlocs_with_masking(model, fmodel, sin_cos_table=None,
                             rotamer_manager=None):
  log = model.log
  if rotamer_manager is None:
    rotamer_manager = mmtbx.idealized_aa_residues.rotamer_manager.load(
      residues = None, rotamers = "favored_allowed")
  if sin_cos_table is None:
    sin_cos_table = scitbx.math.sin_cos_table(n=10000)
  #
  crystal_symmetry = fmodel.xray_structure.crystal_symmetry()
  cg = maptbx.crystal_gridding(
    unit_cell        = crystal_symmetry.unit_cell(),
    space_group_info = crystal_symmetry.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = 0.3)
  #
  map_data = get_map(fmodel=fmodel, crystal_gridding=cg, map_type='2mFo-DFc')
  #
  protein_selection = model.selection("protein")
  sites_cart = model.get_sites_cart()
  sites_cart_start = sites_cart.deep_copy()
  cis = model.get_hierarchy().get_conformer_indices()
  indices = cis.conformer_indices
  vals = list(set(cis.index_altloc_mapping.values()))
  for k,v in zip(cis.index_altloc_mapping.keys(), vals):
    if v == 0: continue
    sel  = indices == 0
    sel |= indices == v
    sel = sel & protein_selection
    model_conf = model.select(sel)
    # Compute mask for other conformers
    model_other = model.select(~sel)
    ss = "not (name CA or name O or name N or name C)"
    sssel = model_other.selection(string=ss)
    model_other = model_other.select(sssel)
    model_other = model_other.remove_hydrogens()
    mask_p1 = mmtbx.masks.mask_from_xray_structure(
      xray_structure           = model_other.get_xray_structure(),
      p1                       = True,
      for_structure_factors    = True,
      solvent_radius           = 0,
      shrink_truncation_radius = 0,
      atom_radius              = 0.5,
      n_real                   = map_data.accessor().all(),
      in_asu                   = False).mask_data
    maptbx.unpad_in_place(map=mask_p1)
    #
    map_data_masked = map_data * mask_p1
    #
    # FOR DEBUGGING
    #
    #with open("%s%d.pdb"%(k,v),"w") as fo:
    #  fo.write(model_conf.model_as_pdb())
    #write_map_file(crystal_symmetry = model.crystal_symmetry(),
    #               map_data         = map_data_masked,
    #               file_name        = "%s%d.map"%(k,v))
    #
    result = mmtbx.refinement.real_space.fit_residues.run(
      pdb_hierarchy     = model_conf.get_hierarchy(),
      vdw_radii         = model_conf.get_vdw_radii(),
      crystal_symmetry  = model_conf.crystal_symmetry(),
      map_data          = map_data_masked,
      backbone_sample   = False,
      rotatable_hd      = model_conf.rotatable_hd_selection(iselection=False),
      rotamer_manager   = rotamer_manager,
      sin_cos_table     = sin_cos_table,
      mon_lib_srv       = model_conf.get_mon_lib_srv(),
      fit_altlocs_method= "masking",
      log               = log)
    #
    sites_cart = sites_cart.set_selected(
      sel, result.pdb_hierarchy.atoms().extract_xyz())
    model_conf.set_sites_cart(result.pdb_hierarchy.atoms().extract_xyz())
    #
    # FOR DEBUGGING
    #
    #with open("%s%d_fitted.pdb"%(k,v),"w") as fo:
    #  fo.write(model_conf.model_as_pdb())
    #
    model.set_sites_cart(sites_cart)
    fmodel.xray_structure.set_sites_cart(model.get_sites_cart())
    fmodel.update_xray_structure(update_f_calc=True)
    #
    map_data = get_map(fmodel=fmodel, crystal_gridding=cg, map_type='2mFo-DFc')
  #
  # FOR DEBUGGING
  #
  #with open("fitted1.pdb","w") as fo:
  #  fo.write(model.model_as_pdb())
  #print("pre-FINAL:", fmodel.r_work())
  #
  sites_cart = model.get_sites_cart()
  uc = model.crystal_symmetry().unit_cell()
  for m in model.get_hierarchy().models():
    for chain in m.chains():
      for residue_group in chain.residue_groups():
        conformers = residue_group.conformers()
        rs = []
        for conformer in residue_group.conformers():
          residue = conformer.only_residue()
          sel = residue.atoms().extract_i_seq()

          if len(residue_group.conformers()) > 1:
            t1 = maptbx.real_space_target_simple(
              unit_cell   = uc,
              density_map = map_data,
              sites_cart  = residue.atoms().extract_xyz())
            mmtbx.refinement.real_space.fit_residue.tune_up(
               target_map           = map_data,
               residue              = residue,
               mon_lib_srv          = model.get_mon_lib_srv(),
               rotamer_manager      = rotamer_manager.rotamer_evaluator,
               unit_cell            = uc,
               monitor              = None,
               torsion_search_start = -20,
               torsion_search_stop  = 20,
               torsion_search_step  = 1)
            t2 = maptbx.real_space_target_simple(
              unit_cell   = uc,
              density_map = map_data,
              sites_cart  = residue.atoms().extract_xyz())
            if t2>t1:
              sites_cart = sites_cart.set_selected(sel, residue.atoms().extract_xyz())
            #print("             tuneup:", t1, t2)

          ts = maptbx.real_space_target_simple(
            unit_cell   = uc,
            density_map = map_data,
            sites_cart  = sites_cart_start.select(sel))
          te = maptbx.real_space_target_simple(
            unit_cell   = uc,
            density_map = map_data,
            sites_cart  = sites_cart.select(sel))
          f = ""
          if te > ts: f = "<<< better"
          r = (ts-te)/(ts+te)*100*2
          label = "%2s %5s %s%4s"%(
            chain.id, residue.resseq, conformer.altloc, residue.resname)
          #print(label, "%6.2f %6.2f %6.2f"%(ts, te, r), f)
          rs.append(r)
        #
        if len(rs)>0:
          got_better = 0
          got_worse  = 0
          for r in rs:
            if r<-1 or r<-1 : got_better+=1
            if r> 1 or r> 1 : got_worse +=1
          if got_better>0 and got_worse==0:
            #print(" "*30, "<<< fit accepted")
            for conformer in residue_group.conformers():
              sel = conformer.atoms().extract_i_seq()
              sel = flex.bool(model.size(), sel)
              sites_cart_start = sites_cart_start.set_selected(sel, sites_cart)
  #
  model.set_sites_cart(sites_cart_start)
  fmodel.xray_structure.set_sites_cart(model.get_sites_cart())
  fmodel.update_xray_structure(update_f_calc=True)
  #print("FINAL:", fmodel.r_work())
  #
  # FOR DEBUGGING
  #
  #with open("fitted2.pdb","w") as fo:
  #  fo.write(model.model_as_pdb())

class run(object):
  def __init__(self,
               pdb_hierarchy,
               crystal_symmetry,
               rotamer_manager,
               sin_cos_table,
               mon_lib_srv,
               filter_by_cc = False,
               d_min = None,
               map_data_scale = 2.5,
               diff_map_data_threshold = -2.5,
               exclude_selection=None,
               rotatable_hd=None,
               bselection=None,
               map_data=None,
               vdw_radii=None,
               backbone_sample=False, # XXX
               diff_map_data=None,
               fit_altlocs_method=None,
               log = None):
    adopt_init_args(self, locals())
    assert fit_altlocs_method in [None, 'sampling', 'masking']
    # Final filtering by CC_mask requires the map and resolution
    if self.filter_by_cc:
      assert self.d_min is not None
      assert self.map_data is not None
    #
    self.pdb_hierarchy_start = self.pdb_hierarchy.deep_copy()
    self.processed = 0
    self.total_time_residue_loop = 0
    t0 = time.time()
    self.atoms = self.pdb_hierarchy.atoms()
    self.cmv = None
    if(self.map_data is not None):
      self.cmv = mmtbx.refinement.real_space.common_map_values(
        pdb_hierarchy = self.pdb_hierarchy,
        unit_cell     = self.crystal_symmetry.unit_cell(),
        map_data      = self.map_data)
    self.mes = []
    if(self.bselection is None):
      o = mmtbx.refinement.real_space.side_chain_fit_evaluator(
        pdb_hierarchy           = self.pdb_hierarchy,
        crystal_symmetry        = self.crystal_symmetry,
        rotamer_evaluator       = self.rotamer_manager.rotamer_evaluator,
        map_data                = self.map_data,
        diff_map_data           = self.diff_map_data,
        map_data_scale          = self.map_data_scale,
        exclude_selection       = self.exclude_selection,
        diff_map_data_threshold = self.diff_map_data_threshold)
      self.mes.extend(o.mes)
      self.bselection = o.sel_all() # or all_possible() ?
    if(self.log is None): self.log = sys.stdout
    self.special_position_settings = crystal.special_position_settings(
      crystal_symmetry = self.crystal_symmetry)
    # Even better would be to pass it here. Ideally just use model
    asc = self.pdb_hierarchy.atom_selection_cache()
    self.selection_water_as_set = \
      set(asc.selection(string = "water").iselection())
    self.target_map = map_data
    self.mes.append("outliers start: %d"%self.count_outliers())
    #
    self.loop(function = self.one_residue_iteration)
    #
    self.mes.append("outliers final: %d"%self.count_outliers())
    #
    self.mes.append("residues processed: %d"%self.processed)
    if(self.processed > 0):
      self.mes.append("average time/residue: %6.4f"%(
        self.total_time_residue_loop/self.processed))
    self.mes.append("time to fit residues: %6.4f"%(time.time()-t0))
    # Final chack by CCmask per residue
    if self.filter_by_cc:
      self._filter_by_cc()

  def _filter_by_cc(self):
    def __get_map_from_hierarchy(ph):
      xrs = ph.extract_xray_structure(
        crystal_symmetry = self.crystal_symmetry)
      hd_sel = xrs.hd_selection()
      xrs = xrs.select(~hd_sel)
      fc = xrs.structure_factors(d_min=self.d_min).f_calc()
      crystal_gridding = maptbx.crystal_gridding(
        unit_cell             = xrs.unit_cell(),
        space_group_info      = xrs.space_group_info(),
        pre_determined_n_real = self.map_data.accessor().all(),
        symmetry_flags        = maptbx.use_space_group_symmetry)
      fft_map = fc.fft_map(crystal_gridding = crystal_gridding)
      return fft_map.real_map_unpadded()
    # XXX This may be a memory bottleneck (3 maps in memory!)
    map_calc_start = __get_map_from_hierarchy(ph=self.pdb_hierarchy_start)
    map_calc_final = __get_map_from_hierarchy(ph=self.pdb_hierarchy)
    #
    for m1,m2 in zip(self.pdb_hierarchy_start.models(),
                     self.pdb_hierarchy.models()):
      for c1,c2 in zip(m1.chains(), m2.chains()):
        for rg1, rg2 in zip(c1.residue_groups(), c2.residue_groups()):
          confs1 = rg1.conformers()
          confs2 = rg2.conformers()
          if(len(confs1)>1): continue
          for con1, con2 in zip(confs1, confs2):
            res1 = con1.only_residue()
            res2 = con2.only_residue()
            if(not self.bselection[res1.atoms()[0].i_seq]): continue
            atoms1, atoms2 = res1.atoms(), res2.atoms()
            dist = flex.sqrt((atoms1.extract_xyz() -
                              atoms2.extract_xyz()).dot())
            moved = False
            if(flex.max(dist) > 0.5): moved = True
            if(moved):
              sites_cart_start = flex.vec3_double()
              sites_cart_final = flex.vec3_double()
              for a1, a2 in zip(atoms1, atoms2):
                if a1.element_is_hydrogen(): continue
                sites_cart_start.append(a1.xyz)
                sites_cart_final.append(a2.xyz)
              cc_start = correlation.from_map_map_atoms(
                map_1      = map_calc_start,
                map_2      = self.map_data,
                sites_cart = sites_cart_start,
                unit_cell  = self.crystal_symmetry.unit_cell(),
                radius     = 2)
              cc_final = correlation.from_map_map_atoms(
                map_1      = map_calc_final,
                map_2      = self.map_data,
                sites_cart = sites_cart_start,
                unit_cell  = self.crystal_symmetry.unit_cell(),
                radius     = 2)
              if(cc_final < cc_start and abs(cc_final-cc_start)>0.02):
                re = self.rotamer_manager.rotamer_evaluator
                re_start = re.evaluate_residue(res1)
                re_final = re.evaluate_residue(res2)
                k = "REVERT: %s %s %s cc_start: %6.4f (%s) cc_final: %6.4f (%s)"%(
                  c1.id, str(res1.resseq), res1.resname, cc_start, re_start,
                  cc_final, re_final)
                print(k, file=self.log)
                rg2.atoms().set_xyz(rg1.atoms().extract_xyz())

  def show(self, prefix=""):
    for m in self.mes:
      print("%s%s"%(prefix,m), file=self.log)
    self.log.flush()

  def _selection_around_minus_self(self, residue, radius):
    if(self.special_position_settings is None): return None
    if(self.vdw_radii is None): return None
    residue_i_selection = flex.size_t()
    for a in residue.atoms():
      if(not a.name.strip() in ["N", "C", "O"]):
        residue_i_selection.append(a.i_seq)
    sites_cart = self.atoms.extract_xyz()
    #
    residue_b_selection = flex.bool(sites_cart.size(), residue_i_selection)
    selection_around_residue = self.special_position_settings.pair_generator(
      sites_cart      = sites_cart,
      distance_cutoff = radius
        ).neighbors_of(primary_selection = residue_b_selection).iselection()

    residue_i_selection = residue.atoms().extract_i_seq()

    selection_around_residue_minus_residue = flex.size_t(
      list(set(selection_around_residue).difference(
        set(residue_i_selection)).difference(self.selection_water_as_set)))
    # exclude rotatable H
    selection_around_residue_minus_residue_minus_rotatableH = flex.size_t()
    for s in selection_around_residue_minus_residue:
      if(not self.rotatable_hd[s]):
        selection_around_residue_minus_residue_minus_rotatableH.append(s)
    return selection_around_residue_minus_residue_minus_rotatableH

  def get_nonbonded_bumpers(self, residue, radius):
    #
    # Symmetry-related atoms are treated differently and less comprihensively
    # See fit_residue.py in "def loop(..)"
    #
    if(self.special_position_settings is None): return None
    if(self.vdw_radii is None): return None
    selection_around_residue_minus_residue_minus_rotatableH = \
      self._selection_around_minus_self(residue=residue, radius=radius)
    #
    radii = flex.double()
    sites_cart = flex.vec3_double()
    for i in selection_around_residue_minus_residue_minus_rotatableH:
      atom = self.atoms[i]
      an = atom.name.strip()
      rad = mmtbx.refinement.real_space.get_radius(
        atom = atom, vdw_radii = self.vdw_radii)
      good = True
      if(self.diff_map_data is not None):
        sf = self.crystal_symmetry.unit_cell().fractionalize(atom.xyz)
        mv = self.diff_map_data.eight_point_interpolation(sf)
        if(mv<-2.): good = False
      if(residue.resname=="CYS" and atom.parent().resname=="CYS" and an=="SG"):
        continue
      if(self.map_data is not None):
        key = "%s_%s_%s"%(
          atom.parent().parent().parent().id, atom.parent().resname, an)
        sf = self.crystal_symmetry.unit_cell().fractionalize(atom.xyz)
        mv = self.map_data.eight_point_interpolation(sf)
        if(mv < self.cmv[key]/3): good = False
      if(good):
        radii.append(rad)
        sites_cart.append(atom.xyz)
    return fit_ext.fixed(sites_cart = sites_cart, radii = radii)

  def one_residue_iteration(self, residue):
    t0 = time.time()
    negate_rad = negate_map_table[residue.resname.strip().lower()]
    if(not negate_rad): return
    xyzrad_bumpers = self.get_nonbonded_bumpers(
      residue=residue, radius=negate_rad)
    self.processed +=1
    if self.fit_altlocs_method=="masking":
      trust_map_values_real=False
    else:
      trust_map_values_real=True
    mmtbx.refinement.real_space.fit_residue.run(
      residue           = residue,
      vdw_radii         = self.vdw_radii,
      xyzrad_bumpers    = xyzrad_bumpers,
      backbone_sample   = self.backbone_sample,
      unit_cell         = self.crystal_symmetry.unit_cell(),
      target_map        = self.target_map,
      target_map_for_cb = self.target_map,
      mon_lib_srv       = self.mon_lib_srv,
      rotamer_manager   = self.rotamer_manager,
      rotatable_hd      = self.rotatable_hd,
      sin_cos_table     = self.sin_cos_table,
      cmv               = self.cmv,
      trust_map_values_real = trust_map_values_real,
      log               = self.log)
    self.total_time_residue_loop += (time.time()-t0)

  def loop(self, function):
    get_class = iotbx.pdb.common_residue_names_get_class
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          conformers = residue_group.conformers()
          # DEAL WITH ALTLOCS
          if(len(conformers)>1 and self.fit_altlocs_method=='sampling'):
            n = len(conformers)
            overall = {}
            for conformer in conformers:
              residue = conformer.only_residue()

              #if int(residue.resseq) != 30: continue

              #
              rn = residue.resname.strip().lower()
              if rn in ["ala", "pro", "gly"]: continue
              negate_rad = negate_map_table[residue.resname.strip().lower()]
              xrs = self.pdb_hierarchy.extract_xray_structure(
                crystal_symmetry = self.crystal_symmetry)
              negate_selection = mmtbx.refinement.real_space.selection_around_to_negate(
                xray_structure          = xrs,
                selection_within_radius = negate_rad,
                iselection              = residue_group.atoms().extract_i_seq())

              target_map_work = mmtbx.refinement.real_space.\
                negate_map_around_selected_atoms_except_selected_atoms(
                  xray_structure   = xrs,
                  map_data         = self.target_map.deep_copy(),
                  negate_selection = negate_selection,
                  atom_radius      = 1.5)

              if(get_class(residue.resname) != "common_amino_acid"): continue
              sel = residue.atoms().extract_i_seq()
              residue_dc = residue.standalone_copy()
              residue_dc.atoms().reset_i_seq()
              o = mmtbx.refinement.real_space.fit_residue.find_all_conformers(
                residue     = residue_dc,
                map_data    = target_map_work, #self.target_map,
                mon_lib_srv = self.mon_lib_srv,
                unit_cell   = self.crystal_symmetry.unit_cell(),
                threshold   = 1.0,
                rotamer_evaluator = rotamer_manager)
              for k,v in o.unique.items():
                overall.setdefault(k, []).append(v)
              #sites_cart_conformers = o.sorted_by_map_value(n=n)
              #if sites_cart_conformers is not None:
              ##  for sites_cart, conformer in zip(sites_cart_conformers, conformers):
              ##    residue = conformer.only_residue()
              ##    residue.atoms().set_xyz( sites_cart )
              #  residue.atoms().set_xyz( sites_cart_conformers[0] )
            unique = {}
            for k,v in overall.items():
              vbest  = -1.e9
              ssbest = None
              for it in v:
                vv, ss = it
                if vv > vbest:
                  vbest = vv
                  ssbest = ss
              unique[k] = [vbest, ssbest]
            if unique != {}:
              print("OVERALL unique")
              code=flex.std_string()
              vals=flex.double()
              for k, v in unique.items():
                print("   ", k, v)
                code.append(k)
                vals.append(v[0])
              sel = flex.sort_permutation(vals, reverse=True)
              print(list(vals.select(sel)))



          # SINGLE CONFORMATION
          else:
            for conformer in residue_group.conformers():
              residue = conformer.only_residue()
              if(not self.bselection[residue.atoms()[0].i_seq]): continue
              xyz_start = residue.atoms().extract_xyz()
              function(residue = residue)
              # Check for symmetry clash
              sels = self._selection_around_minus_self(residue=residue, radius=1.5)
              if(sels is not None and sels.size()>0):
                print("   revert: symmetry clash", file=self.log)
                residue.atoms().set_xyz(xyz_start)
                # XXX
                # XXX For debugging
                # XXX
                #atoms=self.pdb_hierarchy.atoms()
                #for s in sels:
                #  atom = atoms[s]
                #  key = "%s_%s_%s"%(
                #    atom.parent().parent().parent().id, atom.parent().resname, atom.name)
                #  print(key, residue.resname, "LOOK-"*10)

  def count_outliers(self):
    o = mmtbx.refinement.real_space.side_chain_fit_evaluator(
      pdb_hierarchy     = self.pdb_hierarchy,
      crystal_symmetry  = self.crystal_symmetry,
      rotamer_evaluator = self.rotamer_manager.rotamer_evaluator)
    return o.cntr_outliers

# XXX
# XXX Looks obsolete. Not used anywhere.
# XXX
# XXX Where is this used? Looks like severe duplication!
#class fix_outliers(object):
#  def __init__(self,
#               pdb_hierarchy,
#               rotamer_manager,
#               sin_cos_table,
#               mon_lib_srv,
#               f_map=None,
#               fdiff_map=None,
#               unit_cell=None,
#               accept_only_if_max_shift_is_smaller_than=None,
#               log = None):
#    adopt_init_args(self, locals())
#    assert [f_map, fdiff_map, unit_cell].count(None) in [0,3]
#    ac = accept_only_if_max_shift_is_smaller_than
#    if(self.log is None): self.log = sys.stdout
#    get_class = iotbx.pdb.common_residue_names_get_class
#    for model in self.pdb_hierarchy.models():
#      for chain in model.chains():
#        for residue_group in chain.residue_groups():
#          conformers = residue_group.conformers()
#          if(len(conformers)>1): continue
#          for conformer in residue_group.conformers():
#            residue = conformer.only_residue()
#            id_str="%s%s%s"%(chain.id,residue.resname,residue.resseq.strip())
#            if(get_class(residue.resname) == "common_amino_acid"):
#              # Idealize rotamer: move to nearest rotameric state
#              re = self.rotamer_manager.rotamer_evaluator
#              if(re.evaluate_residue(residue)=="OUTLIER"):
#                go=True
#                if([f_map, fdiff_map, unit_cell].count(None)==0):
#                  go = mmtbx.refinement.real_space.need_sidechain_fit(
#                    residue           = residue,
#                    rotamer_evaluator = self.rotamer_manager.rotamer_evaluator,
#                    mon_lib_srv       = self.mon_lib_srv,
#                    unit_cell         = self.unit_cell,
#                    f_map             = f_map,
#                    fdiff_map         = fdiff_map)
#                if(go):
#                  mmtbx.refinement.real_space.fit_residue.run(
#                    residue         = residue,
#                    backbone_sample = False,
#                    mon_lib_srv     = self.mon_lib_srv,
#                    rotamer_manager = self.rotamer_manager,
#                    sin_cos_table   = self.sin_cos_table,
#                    accept_only_if_max_shift_is_smaller_than = ac,
#                    log = self.log)
