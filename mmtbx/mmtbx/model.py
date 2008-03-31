from cctbx.array_family import flex
import math, time
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal, not_approx_equal
import sys, random
from stdlib import math
from cctbx import xray
from cctbx import adptbx
import mmtbx.restraints
from iotbx import pdb
from cctbx import geometry_restraints
from cctbx.geometry_restraints.lbfgs import lbfgs as cctbx_geometry_restraints_lbfgs
import scitbx.lbfgs
from libtbx.utils import Sorry, user_plus_sys_time
from mmtbx.tls import tools
from cctbx import adp_restraints
from mmtbx import ias
from mmtbx import utils
from mmtbx import model_statistics
from mmtbx.solvent import ordered_solvent
import iotbx.pdb


time_model_show = 0.0

class xh_connectivity_table(object):
  # XXX need angle information as well
  def __init__(self, geometry, xray_structure):
    bond_proxies_simple = geometry.geometry.pair_proxies(sites_cart =
      xray_structure.sites_cart()).bond_proxies.simple
    self.table = []
    scatterers = xray_structure.scatterers()
    for proxy in bond_proxies_simple:
      i_seq, j_seq = proxy.i_seqs
      i_x, i_h = None, None
      if(scatterers[i_seq].element_symbol() in ["H", "D"]):
        i_h = i_seq
        i_x = j_seq
        site_x = scatterers[i_x].site
        site_h = scatterers[i_h].site
        const_vect = flex.double(site_h)-flex.double(site_x)
        distance_model = xray_structure.unit_cell().distance(site_x, site_h)
        self.table.append([i_x, i_h, const_vect, proxy.distance_ideal,
                           distance_model])
      if(scatterers[j_seq].element_symbol() in ["H", "D"]):
        i_h = j_seq
        i_x = i_seq
        site_x = scatterers[i_x].site
        site_h = scatterers[i_h].site
        const_vect = flex.double(site_h)-flex.double(site_x)
        distance_model = xray_structure.unit_cell().distance(site_x, site_h)
        self.table.append([i_x, i_h, const_vect, proxy.distance_ideal,
                           distance_model])

class manager(object):
  def __init__(self, xray_structure,
                     atom_attributes_list,
                     processed_pdb_files_srv = None,
                     restraints_manager = None,
                     ias_xray_structure = None,
                     refinement_flags = None,
                     ias_manager = None,
                     wilson_b = None,
                     tls_groups = None,
                     anomalous_scatterer_groups = None,
                     log = None):
    self.log = log
    self.processed_pdb_files_srv = processed_pdb_files_srv
    self.restraints_manager = restraints_manager
    self.xray_structure = xray_structure
    self.xray_structure_initial = self.xray_structure.deep_copy_scatterers()
    self.atom_attributes_list = atom_attributes_list
    self.refinement_flags = refinement_flags
    self.wilson_b = wilson_b
    self.tls_groups = tls_groups
    if(anomalous_scatterer_groups is not None and
      len(anomalous_scatterer_groups) == 0):
      anomalous_scatterer_groups = None
    self.anomalous_scatterer_groups = anomalous_scatterer_groups
    # IAS related, need a real cleaning!
    self.ias_manager = ias_manager
    self.ias_xray_structure = ias_xray_structure
    self.use_ias = False
    self.ias_selection = None

  def xh_connectivity_table(self):
    result = None
    if(self.restraints_manager is not None):
      if(self.xray_structure.hd_selection().count(True) > 0):
        xray_structure = self.xray_structure
        if(self.ias_selection and self.ias_selection.count(True) > 0):
          xray_structure = self.xray_structure.select(~self.ias_selection)
        result = xh_connectivity_table(
          geometry       = self.restraints_manager,
          xray_structure = xray_structure).table
    return result

  def idealize_h(self, xh_bond_distance_deviation_limit=0, show=True): # XXX _limit is not used
    if(self.xray_structure.hd_selection().count(True) > 0):
      sol_hd = self.solvent_selection().set_selected(
        ~self.xray_structure.hd_selection(), False)
      mac_hd = self.xray_structure.hd_selection().set_selected(
        self.solvent_selection(), False)
      selx = ~self.xray_structure.hd_selection()
      if(self.ias_selection is not None):
        selx.set_selected(self.ias_selection, False)
      sites_cart_mac_before = self.xray_structure.sites_cart().select(selx)
      xhd = flex.double()
      for t in self.xh_connectivity_table():
        xhd.append(abs(t[-1]-t[-2]))
        if(abs(t[3]-t[4]) < xh_bond_distance_deviation_limit):
          sol_hd[t[1]] = False
          mac_hd[t[1]] = False
      if(show):
        print >> self.log, \
        "X-H deviation from ideal before regularization (bond): mean=%6.3f max=%6.3f"%\
        (flex.mean(xhd), flex.max(xhd))
      for sel_pair in [(mac_hd, False), (sol_hd, True)]*2:
        if(sel_pair[0].count(True) > 0):
          sel = sel_pair[0]
          if(self.ias_selection is not None and self.ias_selection.count(True) > 0):
            sel = sel.select(~self.ias_selection)
          self.geometry_minimization(
            selection = sel,
            bond      = True,
            nonbonded = sel_pair[1],
            angle     = True,
            dihedral  = True,
            chirality = True,
            planarity = True)
      sites_cart_mac_after = self.xray_structure.sites_cart().select(selx)
      assert approx_equal(flex.max(sites_cart_mac_before.as_double() -
        sites_cart_mac_after.as_double()), 0)
      xhd = flex.double()
      for t in self.xh_connectivity_table():
        xhd.append(abs(t[-1]-t[-2]))
      if(show):
        print >> self.log,\
        "X-H deviation from ideal after  regularization (bond): mean=%6.3f max=%6.3f"%\
        (flex.mean(xhd), flex.max(xhd))

  def add_to_aal(self, next_to_i_seq, attr):
    self.atom_attributes_list = self.atom_attributes_list[:next_to_i_seq+1] + \
      [attr] + self.atom_attributes_list[next_to_i_seq+1:]

  def extact_water(self):
    top = {}
    solvent_sel = self.solvent_selection()
    get_class = iotbx.pdb.common_residue_names_get_class
    for i_seq, aal in enumerate(self.atom_attributes_list):
      if(get_class(name = aal.resName) == "common_water"):
        assert solvent_sel[i_seq]
        top.setdefault(aal.chainID,      {}).\
            setdefault(aal.residue_id(), []).append([i_seq,aal])
      else:
        assert not solvent_sel[i_seq]
    water_residues = []
    for v0 in top.values():
      for v1 in v0.values():
        water_residues.append(v1)
        assert len(v1) <= 3
    return water_residues

  def renumber_water(self):
    water_residues = self.extact_water()
    for i_seq, water_residue in enumerate(water_residues):
      assert len(water_residue) <= 3
      for wra in water_residue:
        wra[1].resSeq = i_seq

  def add_hydrogens(self, element = "H"):
    water_residues = self.extact_water()
    result = []
    n_built = 0
    for water_residue in water_residues:
      if(len(water_residue) == 2):
        n_built += 1
        o_attr = None
        h_attr = None
        for wr in water_residue:
          if(wr[1].element.strip() == "O"): o_attr = wr
          else: h_attr = wr
        assert [o_attr, h_attr].count(None) == 0
        h_name = h_attr[1].name.strip()
        if(len(h_name) == 1):
          atom_name = h_name+"1"
        elif(len(h_name) == 2):
          if(h_name[0].isdigit()):
            if(int(h_name[0]) == 1): atom_name = h_name[1]+"2"
            elif(int(h_name[0]) == 2): atom_name = h_name[1]+"1"
            else: raise RuntimeError
          elif(h_name[1].isdigit()):
            if(int(h_name[1]) == 1): atom_name = h_name[0]+"2"
            elif(int(h_name[1]) == 2): atom_name = h_name[0]+"1"
            else: raise RuntimeError
          else: raise RuntimeError
        else: raise RuntimeError
        h_attr = create_atom_attr(other_attr     = o_attr[1],
                                  i_seq          = o_attr[0],
                                  xray_structure = self.xray_structure,
                                  element        = h_attr[1].element,
                                  atom_name      = atom_name)
        result.append([o_attr[0], h_attr])
      elif(len(water_residue) == 1):
        n_built += 1
        water_residue = water_residue[0]
        assert water_residue[1].element.strip() == "O"
        for atom_name in (element+"1", element+"2"):
          h_attr = create_atom_attr(other_attr     = water_residue[1],
                                    i_seq          = water_residue[0],
                                    xray_structure = self.xray_structure,
                                    element        = element,
                                    atom_name      = atom_name)
          result.append([water_residue[0], h_attr])
    #
    result.sort()
    result.reverse()
    next_to_i_seqs = []
    print >> self.log, "Number of H added:", len(result)
    if(len(result) > 0):
      for res in result:
        i_seq, attr = res[0], res[1]
        next_to_i_seqs.append(i_seq)
        scatterer = xray.scatterer(
          label           = attr.name,
          scattering_type = attr.element,
          site            = self.xray_structure.unit_cell().fractionalize(attr.coordinates),
          u               = adptbx.b_as_u(attr.tempFactor),
          occupancy       = attr.occupancy)
        self.xray_structure.add_scatterer(
          scatterer = scatterer,
          insert_at_index = i_seq+1)
        self.add_to_aal(next_to_i_seq = i_seq, attr = attr)
      if(self.refinement_flags is not None):
        self.refinement_flags.add(next_to_i_seqs   = next_to_i_seqs,
                                  sites_individual = True)
      # XXX very inefficient: re-process PDB from scratch and create restraints
      raw_records = self.write_pdb_file()
      processed_pdb_file, pdb_inp = self.processed_pdb_files_srv.\
        process_pdb_files(raw_records = raw_records)
      new_xray_structure = processed_pdb_file.xray_structure(
        show_summary = False).deep_copy_scatterers()
      new_atom_attributes_list = \
        processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list[:]
      assert approx_equal(new_xray_structure.sites_cart(),
        self.xray_structure.sites_cart(), 0.001)
      assert len(new_atom_attributes_list) == len(self.atom_attributes_list)
      for attr1, attr2 in zip(self.atom_attributes_list, new_atom_attributes_list):
        assert attr1.name.strip() == attr2.name.strip()
        assert attr1.element.strip() == attr2.element.strip()
      self.atom_attributes_list = new_atom_attributes_list[:]
      # XXX now we gonna loose old grm (with all NCS, edits, etc...)
      if(self.restraints_manager.ncs_groups is not None):
        raise Sorry("Hydrogen building is not compatible with NCS refinement.")

      sctr_keys = \
        self.xray_structure.scattering_type_registry().type_count_dict().keys()
      has_hd = "H" in sctr_keys or "D" in sctr_keys
      geometry = processed_pdb_file.geometry_restraints_manager(
        show_energies      = False,
        plain_pairs_radius = self.restraints_manager.geometry.plain_pairs_radius,
        edits              = None, #self.params.geometry_restraints.edits, XXX
        assume_hydrogens_all_missing = not has_hd)
      # self.remove_selected_geometry_restraints(manager = geometry) XXX this is lost too
      new_restraints_manager = mmtbx.restraints.manager(
        geometry      = geometry,
        normalization = self.restraints_manager.normalization)
      self.restraints_manager = new_restraints_manager
      self.idealize_h()

  def geometry_minimization(self,
                            max_number_of_iterations = 500,
                            number_of_macro_cycles   = 3,
                            selection                = None,
                            bond                     = False,
                            nonbonded                = False,
                            angle                    = False,
                            dihedral                 = False,
                            chirality                = False,
                            planarity                = False):
    assert max_number_of_iterations+number_of_macro_cycles > 0
    assert [bond,nonbonded,angle,dihedral,chirality,planarity].count(False) < 6
    from mmtbx.command_line import geometry_minimization
    import scitbx.lbfgs
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = max_number_of_iterations)
    geometry_restraints_flags = geometry_restraints.flags.flags(
      bond      = bond,
      nonbonded = nonbonded,
      angle     = angle,
      dihedral  = dihedral,
      chirality = chirality,
      planarity = planarity)
    for i in xrange(number_of_macro_cycles):
      sites_cart = self.xray_structure.sites_cart()
      sites_cart_orig = sites_cart.deep_copy()
      if(self.ias_selection is not None and self.ias_selection.count(True) > 0):
        sites_cart = sites_cart.select(~self.ias_selection)
      minimized = geometry_minimization.lbfgs(
        sites_cart                  = sites_cart,
        geometry_restraints_manager = self.restraints_manager.geometry,
        geometry_restraints_flags   = geometry_restraints_flags,
        lbfgs_termination_params    = lbfgs_termination_params,
        sites_cart_selection        = selection)
      if(self.ias_selection is not None):
        for i_seq, ias_s in enumerate(self.ias_selection): # assumes that IAS appended to the back
          if(not ias_s):
            sites_cart_orig[i_seq] = sites_cart[i_seq]
      else:
        sites_cart_orig = sites_cart
      self.xray_structure.set_sites_cart(sites_cart = sites_cart_orig)

  def extract_ncs_groups(self):
    result = None
    if(self.restraints_manager.ncs_groups is not None):
      result = self.restraints_manager.ncs_groups.extract_ncs_groups(
        sites_cart = self.xray_structure.sites_cart())
    return result

  def deep_copy(self):
    return self.select(selection = flex.bool(
      self.xray_structure.scatterers().size(), True))

  def add_ias(self, fmodel=None, ias_params=None, file_name=None,
                                                             build_only=False):
    if(self.ias_manager is not None):
       self.remove_ias()
       fmodel.update_xray_structure(xray_structure = self.xray_structure,
                                    update_f_calc = True)
    print >> self.log, ">>> Adding IAS.........."
    self.old_refinement_flags = None
    if not build_only: self.use_ias = True
    self.ias_manager = ias.manager(
                    geometry             = self.restraints_manager.geometry,
                    atom_attributes_list = self.atom_attributes_list,
                    xray_structure       = self.xray_structure,
                    fmodel               = fmodel,
                    params               = ias_params,
                    file_name            = file_name,
                    log                  = self.log)
    if(not build_only):
      self.ias_xray_structure = self.ias_manager.ias_xray_structure
      ias_size = self.ias_xray_structure.scatterers().size()
      tail = flex.bool(ias_size, True)
      tail_false = flex.bool(ias_size, False)
      self.ias_selection = flex.bool(
                      self.xray_structure.scatterers().size(),False)
      self.ias_selection.extend(tail)
      self.xray_structure.concatenate_inplace(other = self.ias_xray_structure)
      print >> self.log, "Scattering dictionary for combined xray_structure:"
      self.xray_structure.scattering_type_registry().show()
      self.xray_structure_initial.concatenate_inplace(
                                           other = self.ias_xray_structure)
      if(self.refinement_flags is not None):
         self.old_refinement_flags = self.refinement_flags.deep_copy()
         # define flags
         ssites = flex.bool(self.ias_xray_structure.scatterers().size(), False)
         sadp = flex.bool(self.ias_xray_structure.scatterers().size(), False)
         # XXX set occ refinement ONLY for involved atoms
         # XXX now it refines only occupancies of IAS !!!
         occupancy_flags = []
         ms = self.ias_selection.count(False)
         for i in range(1, self.ias_selection.count(True)+1):
           occupancy_flags.append([flex.size_t([ms+i-1])])
         # set flags
         self.refinement_flags.inflate(
           sites_individual     = ssites,
           s_occupancies        = occupancy_flags,
           adp_individual_iso   = sadp,
           adp_individual_aniso = sadp)
         # adjust flags
         self.refinement_flags.sites_individual.set_selected(self.ias_selection, False)
         self.refinement_flags.sites_individual.set_selected(~self.ias_selection, True)
         self.refinement_flags.adp_individual_aniso.set_selected(self.ias_selection, False)
         self.refinement_flags.adp_individual_iso.set_selected(self.ias_selection, True)

         #occs = flex.double(self.xray_structure.scatterers().size(), 0.9)
         #self.xray_structure.scatterers().set_occupancies(occs, ~self.ias_selection)
         # D9
         sel = self.xray_structure.scatterers().extract_scattering_types() == "D9"
         self.xray_structure.convert_to_anisotropic(selection = sel)
         self.refinement_flags.adp_individual_aniso.set_selected(sel, True)
         self.refinement_flags.adp_individual_iso.set_selected(sel, False)
    # add to aal:
    i_seq = 0
    for sc in self.ias_xray_structure.scatterers():
      i_seq += 1
      new_atom_name = sc.label.strip()
      if(len(new_atom_name) < 4): new_atom_name = " " + new_atom_name
      while(len(new_atom_name) < 4): new_atom_name = new_atom_name+" "
      new_attr = pdb.atom.attributes(name        = new_atom_name,
                                     resName     = "IAS",
                                     element     = sc.element_symbol(),
                                     is_hetatm   = True,
                                     resSeq      = i_seq)
      self.atom_attributes_list.append(new_attr)


  def remove_ias(self):
    print >> self.log, ">>> Removing IAS..............."
    self.use_ias = False
    if(self.ias_manager is not None):
       self.ias_manager = None
    if(self.old_refinement_flags is not None):
       self.refinement_flags = self.old_refinement_flags.deep_copy()
       self.old_refinement_flags = None
    if(self.ias_selection is not None):
       self.xray_structure.select_inplace(selection = ~self.ias_selection)
       n_non_ias = self.ias_selection.count(False)
       self.ias_selection = None
       self.xray_structure.scattering_type_registry().show()
       self.atom_attributes_list = self.atom_attributes_list[:n_non_ias]

  def show_rigid_bond_test(self, out=None):
    if(out is None): out = sys.stdout
    bond_proxies_simple = \
            self.restraints_manager.geometry.pair_proxies().bond_proxies.simple
    scatterers = self.xray_structure.scatterers()
    unit_cell = self.xray_structure.unit_cell()
    rbt_array = flex.double()
    for proxy in bond_proxies_simple:
        i_seqs = proxy.i_seqs
        i,j = proxy.i_seqs
        atom_i = self.atom_attributes_list[i]
        atom_j = self.atom_attributes_list[j]
        if(atom_i.element.strip() not in ["H","D"] and
                                      atom_j.element.strip() not in ["H","D"]):
           sc_i = scatterers[i]
           sc_j = scatterers[j]
           if(sc_i.flags.use_u_aniso() and sc_j.flags.use_u_aniso()):
              p = adp_restraints.rigid_bond_pair(sc_i.site,
                                                 sc_j.site,
                                                 sc_i.u_star,
                                                 sc_j.u_star,
                                                 unit_cell)
              rbt_value = p.delta_z()*10000.
              rbt_array.append(rbt_value)
              print >> out, "%s %s %10.3f"%(atom_i.name, atom_j.name, rbt_value)
    print >> out, "RBT values (*10000):"
    print >> out, "  mean = %.3f"%flex.mean(rbt_array)
    print >> out, "  max  = %.3f"%flex.max(rbt_array)
    print >> out, "  min  = %.3f"%flex.min(rbt_array)


  def restraints_manager_energies_sites(self,
                                        geometry_flags    = None,
                                        compute_gradients = False,
                                        gradients         = None,
                                        disable_asu_cache = False):
    sites_cart = self.xray_structure.sites_cart()
    if(self.use_ias and self.ias_selection is not None and
       self.ias_selection.count(True) > 0):
      sites_cart = sites_cart.select(~self.ias_selection)
    return self.restraints_manager.energies_sites(
      sites_cart        = sites_cart,
      geometry_flags    = geometry_flags,
      compute_gradients = compute_gradients,
      gradients         = gradients,
      disable_asu_cache = disable_asu_cache)

  def solvent_selection(self):
    water = ordered_solvent.water_ids()
    result = flex.bool()
    get_class = iotbx.pdb.common_residue_names_get_class
    for a in self.atom_attributes_list:
      element = (a.element).strip().upper()
      resName = (a.resName).strip()
      name    = (a.name).strip()
      if((element in water.element_types) and
         (name in water.atom_names) and \
         (get_class(name = resName) == "common_water")):
        result.append(True)
      else: result.append(False)
    return result

  def xray_structure_macromolecule(self):
    sel = self.solvent_selection()
    if(self.use_ias): sel = sel | self.ias_selection
    result = self.xray_structure.select(~sel)
    return result

  def select(self, selection):
    # XXX ignores IAS
    new_atom_attributes_list = []
    for attr, sel in zip(self.atom_attributes_list, selection):
      if(sel): new_atom_attributes_list.append(attr)
    new_refinement_flags = None
    if(self.refinement_flags is not None):
      new_refinement_flags = self.refinement_flags.select(selection)
    new_restraints_manager = None
    if(self.restraints_manager is not None):
      new_restraints_manager = self.restraints_manager.select(
        selection = selection)
      new_restraints_manager.geometry.pair_proxies(sites_cart =
        self.xray_structure.sites_cart().select(selection)) # XXX is it necessary ?
    new = manager(
      processed_pdb_files_srv    = self.processed_pdb_files_srv,
      restraints_manager         = new_restraints_manager,
      xray_structure             = self.xray_structure.select(selection),
      atom_attributes_list       = new_atom_attributes_list,
      refinement_flags           = new_refinement_flags,
      tls_groups                 = self.tls_groups, # XXX not selected, potential bug
      anomalous_scatterer_groups = self.anomalous_scatterer_groups,
      log                        = self.log)
    new.xray_structure_initial = \
      self.xray_structure_initial.deep_copy_scatterers()
    new.xray_structure.scattering_type_registry()
    return new

  def number_of_ordered_solvent_molecules(self):
    return self.solvent_selection().count(True)

  def show_groups(self, rigid_body = None, tls = None,
                        out = None, text="Information about rigid groups"):
    global time_model_show
    timer = user_plus_sys_time()
    selections = None
    if(rigid_body is not None):
       selections = self.refinement_flags.sites_rigid_body
    if(tls is not None): selections = self.refinement_flags.adp_tls
    if(self.refinement_flags.sites_rigid_body is None and
                                 self.refinement_flags.adp_tls is None): return
    assert selections is not None
    if (out is None): out = sys.stdout
    print >> out
    line_len = len("| "+text+"|")
    fill_len = 80 - line_len-1
    upper_line = "|-"+text+"-"*(fill_len)+"|"
    print >> out, upper_line
    next = "| Total number of atoms = %-6d  Number of rigid groups = %-3d                |"
    natoms_total = self.xray_structure.scatterers().size()
    print >> out, next % (natoms_total, len(selections))
    print >> out, "| group: start point:                        end point:                       |"
    print >> out, "|               x      B  atom   residue <>        x      B  atom   residue   |"
    next = "| %5d: %8.3f %6.2f %5s %4s %4s <> %8.3f %6.2f %5s %4s %4s   |"
    sites = self.xray_structure.sites_cart()
    b_isos = self.xray_structure.extract_u_iso_or_u_equiv() * math.pi**2*8
    n_atoms = 0
    for i_seq, selection in enumerate(selections):
        try:
          i_selection = selection.iselection()
          n_atoms += i_selection.size()
        except:
          i_selection = selection
          n_atoms += i_selection.size()
        start = i_selection[0]
        final = i_selection[i_selection.size()-1]
        first = self.atom_attributes_list[start]
        last  = self.atom_attributes_list[final]
        print >> out, next % (i_seq+1, sites[start][0], b_isos[start],
          first.name, first.resName, first.resSeq, sites[final][0],
          b_isos[final], last.name, last.resName, last.resSeq)
    print >> out, "|"+"-"*77+"|"
    print >> out
    out.flush()
    time_model_show += timer.elapsed()

  def remove_solvent(self):
    result = self.select(selection = ~self.solvent_selection())
    return result

  def show_occupancy_statistics(self, out=None, text=""):
    global time_model_show
    timer = user_plus_sys_time()
    # XXX make this more complete and smart
    if(out is None): out = sys.stdout
    print >> out, "|-"+text+"-"*(80 - len("| "+text+"|") - 1)+"|"
    occ = self.xray_structure.scatterers().extract_occupancies()
    occ_min = flex.min(occ)
    occ_max = flex.max(occ)
    n_zeros = (occ < 0.1).count(True)
    percent_small = n_zeros * 100 / occ.size()
    n_large = (occ > 2.0).count(True)
    if(occ_min < 0.0):
       raise Sorry("There are atoms with negative occupancies. Check input "\
                   "PDB file.")
    if(percent_small > 30.0):
       print >> out, "| *** WARNING: there more than 30 % of atoms with small occupancy (< 0.1) *** |"
    if(n_large > 0):
       print >> out, "| *** WARNING: there are some atoms with large occupancy (> 2.0) ***          |"
    if(abs(occ_max-occ_min) >= 0.01):
       print >> out, "| occupancies: max = %-6.2f min = %-6.2f number of "\
                     "occupancies < 0.1 = %-6d |"%(occ_max,occ_min,n_zeros)
    else:
       print >> out, "| occupancies: max = %-6.2f min = %-6.2f number of "\
                     "occupancies < 0.1 = %-6d |"%(occ_max,occ_min,n_zeros)
    print >> out, "|"+"-"*77+"|"
    out.flush()
    time_model_show += timer.elapsed()

  def write_pdb_file(self, out = None, selection = None, xray_structure = None):
    return utils.write_pdb_file(
      xray_structure       = self.xray_structure,
      atom_attributes_list = self.atom_attributes_list,
      selection            = selection,
      out                  = out)

  def add_solvent(self, solvent_xray_structure,
                        atom_name    = "O",
                        residue_name = "HOH",
                        chain_id     = None,
                        refine_occupancies = False,
                        refine_adp = None):
    assert refine_adp is not None
    if(refine_adp == "isotropic"):
      solvent_xray_structure.convert_to_isotropic()
    elif(refine_adp == "anisotropic"):
      occ = solvent_xray_structure.scatterers().extract_occupancies()
      occ_sel_iso = (occ < 0.5).iselection()
      solvent_xray_structure.convert_to_anisotropic(selection = ~occ_sel_iso)
      solvent_xray_structure.convert_to_isotropic(selection = occ_sel_iso)
    else: raise RuntimeError
    ms = self.xray_structure.scatterers().size() #
    self.xray_structure = \
      self.xray_structure.concatenate(solvent_xray_structure)
    occupancy_flags = None
    if(refine_occupancies):
      occupancy_flags = []
      for i in range(1, solvent_xray_structure.scatterers().size()+1):
        occupancy_flags.append([flex.size_t([ms+i-1])])
    if(self.refinement_flags.individual_sites):
      ssites = flex.bool(solvent_xray_structure.scatterers().size(), True)
    else: ssites = None
    if(self.refinement_flags.adp_individual_iso):
      sadp_iso = solvent_xray_structure.use_u_iso()
    else: sadp_iso = None
    if(self.refinement_flags.adp_individual_aniso):
      sadp_aniso = solvent_xray_structure.use_u_aniso()
    else: sadp_aniso = None
    self.refinement_flags.inflate(
      sites_individual       = ssites,
      adp_individual_iso     = sadp_iso,
      adp_individual_aniso   = sadp_aniso,
      s_occupancies          = occupancy_flags)
    new_atom_name = atom_name.strip()
    if(len(new_atom_name) < 4): new_atom_name = " " + new_atom_name
    while(len(new_atom_name) < 4): new_atom_name = new_atom_name+" "
    #
    water_resseqs = flex.size_t()
    get_class = iotbx.pdb.common_residue_names_get_class
    for aattr in self.atom_attributes_list:
      if(get_class(name = aattr.resName) == "common_water"):
        water_resseqs.append(int(aattr.resSeq))
    i_seq = flex.max_default(water_resseqs, 0)
    #
    for sc in solvent_xray_structure.scatterers():
        i_seq += 1
        new_attr = pdb.atom.attributes(name        = new_atom_name,
                                       resName     = residue_name,
                                       chainID     = chain_id,
                                       element     = sc.element_symbol(),
                                       is_hetatm   = True,
                                       resSeq      = i_seq)
        self.atom_attributes_list.append(new_attr)
    geometry = self.restraints_manager.geometry
    number_of_new_solvent = solvent_xray_structure.scatterers().size()
    if(geometry.model_indices is None):
       model_indices = None
    else:
       model_indices = flex.size_t(number_of_new_solvent, 0)
    if(geometry.conformer_indices is None):
       conformer_indices = None
    else:
       conformer_indices = flex.size_t(number_of_new_solvent, 0)
    geometry = geometry.new_including_isolated_sites(
           n_additional_sites  = number_of_new_solvent,
           model_indices       = model_indices,
           conformer_indices   = conformer_indices,
           site_symmetry_table = solvent_xray_structure.site_symmetry_table(),
           nonbonded_types     = flex.std_string(number_of_new_solvent, "OH2"))
    self.restraints_manager = mmtbx.restraints.manager(
                         geometry      = geometry,
                         ncs_groups    = self.restraints_manager.ncs_groups,
                         normalization = self.restraints_manager.normalization)
    if (self.restraints_manager.ncs_groups is not None):
      self.restraints_manager.ncs_groups.register_additional_isolated_sites(
        number=number_of_new_solvent)
    self.restraints_manager.geometry.update_plain_pair_sym_table(
                                 sites_frac = self.xray_structure.sites_frac())
    assert len(self.atom_attributes_list) == \
                                        self.xray_structure.scatterers().size()


  def scale_adp(self, scale_max, scale_min):
    b_isos = self.xray_structure.extract_u_iso_or_u_equiv() * math.pi**2*8
    b_isos_mean = flex.mean(b_isos)
    max_b_iso = b_isos_mean * scale_max
    min_b_iso = b_isos_mean / scale_min
    sel_outliers_max = b_isos > max_b_iso
    sel_outliers_min = b_isos < min_b_iso
    b_isos.set_selected(sel_outliers_max, max_b_iso)
    b_isos.set_selected(sel_outliers_min, min_b_iso)
    self.xray_structure.set_b_iso(values = b_isos)

  def geometry_statistics(self):
    sites_cart = self.xray_structure.sites_cart()
    if(self.use_ias): sites_cart = sites_cart.select(~self.ias_selection)
    return model_statistics.geometry(
      sites_cart         = sites_cart,
      restraints_manager = self.restraints_manager)

  def show_geometry_statistics(self, message = "", out = None):
    global time_model_show
    if(out is None): out = self.log
    timer = user_plus_sys_time()
    result = self.geometry_statistics()
    result.show(message = message, out = out)
    time_model_show += timer.elapsed()
    return result

  def adp_statistics(self):
    return model_statistics.adp(model = self)

  def show_adp_statistics(self,
                          prefix         = "",
                          padded         = False,
                          pdb_deposition = False,
                          out            = None):
    global time_model_show
    if(out is None): out = self.log
    timer = user_plus_sys_time()
    result = self.adp_statistics()
    result.show(out = out, prefix = prefix, padded = padded,
      pdb_deposition = pdb_deposition)
    time_model_show += timer.elapsed()
    return result

  def energies_adp(self, iso_restraints, compute_gradients, use_hd):
    assert self.refinement_flags is not None
    n_aniso = 0
    if(self.refinement_flags.adp_individual_aniso is not None):
      n_aniso = self.refinement_flags.adp_individual_aniso.count(True)
    if(n_aniso == 0):
      energies_adp_iso = self.restraints_manager.energies_adp_iso(
        xray_structure    = self.xray_structure,
        parameters        = iso_restraints,
        use_u_local_only  = iso_restraints.use_u_local_only,
        use_hd            = use_hd,
        compute_gradients = compute_gradients)
      target = energies_adp_iso.target
    else:
      energies_adp_aniso = self.restraints_manager.energies_adp_aniso(
        xray_structure    = self.xray_structure,
        compute_gradients = compute_gradients)
      target = energies_adp_aniso.target
    u_iso_gradients = None
    u_aniso_gradients = None
    if(compute_gradients):
      if(n_aniso == 0):
        u_iso_gradients = energies_adp_iso.gradients
      else:
        u_aniso_gradients = energies_adp_aniso.gradients_aniso_star
        u_iso_gradients = energies_adp_aniso.gradients_iso
    class result(object):
      def __init__(self):
        self.target = target
        self.u_iso_gradients = u_iso_gradients
        self.u_aniso_gradients = u_aniso_gradients
    return result()

  def set_refine_individual_sites(self, selection = None):
    self.xray_structure.scatterers().flags_set_grads(state=False)
    if(selection is None):
      selection = self.refinement_flags.sites_individual
    self.xray_structure.scatterers().flags_set_grad_site(
      iselection = selection.iselection())

  def set_refine_individual_adp(self, selection_iso = None,
                                      selection_aniso = None):
    self.xray_structure.scatterers().flags_set_grads(state=False)
    if(selection_iso is None):
      selection_iso = self.refinement_flags.adp_individual_iso
    if(selection_iso is not None):
      self.xray_structure.scatterers().flags_set_grad_u_iso(
        iselection = selection_iso.iselection())
    if(selection_aniso is None):
      selection_aniso = self.refinement_flags.adp_individual_aniso
    if(selection_aniso is not None):
      self.xray_structure.scatterers().flags_set_grad_u_aniso(
        iselection = selection_aniso.iselection())

def create_atom_attr(other_attr, i_seq, xray_structure, element, atom_name):
  x,y,z = xray_structure.sites_cart()[i_seq]
  r = random.randrange(10000, 11000)/10000. * random.choice([-1,1])
  if(other_attr.altLoc is None): altLoc = " "
  else: altLoc = other_attr.altLoc
  if(other_attr.chainID is None): chainID = " "
  else: chainID = other_attr.chainID
  if(other_attr.resSeq is None): resSeq = 1
  else: resSeq = other_attr.resSeq
  if(other_attr.iCode is None): iCode = " "
  else: iCode = other_attr.iCode
  if(other_attr.segID is None): segID = "    "
  else: segID = other_attr.segID
  if(other_attr.charge is None): charge = "  "
  else: charge = other_attr.charge
  pdb_record = pdb.format_atom_record(
    record_name = other_attr.record_name(),
    serial      = 0,
    name        = atom_name,
    altLoc      = altLoc,
    resName     = other_attr.resName,
    chainID     = chainID,
    resSeq      = resSeq,
    iCode       = iCode,
    site        = (x+r,y,z),
    occupancy   = 1.0,
    tempFactor  = adptbx.u_as_b(xray_structure.extract_u_iso_or_u_equiv()[i_seq]),
    segID       = segID,
    element     = element,
    charge      = charge)
  attr = pdb.atom.attributes()
  attr.set_from_ATOM_record(pdb.parser.pdb_record(pdb_record))
  return attr
