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
from mmtbx import dbe
from mmtbx import utils
from mmtbx import model_statistics


time_model_show = 0.0

class xh_connectivity_table(object):
  def __init__(self, geometry, xray_structure):
    bond_proxies_simple = geometry.geometry.pair_proxies(sites_cart =
      xray_structure.sites_cart()).bond_proxies.simple
    self.table = []
    scatterers = xray_structure.scatterers()
    for proxy in bond_proxies_simple:
      i_seq, j_seq = proxy.i_seqs
      i_x, i_h = None, None
      if(scatterers[i_seq].element_symbol() == "H"):
        i_h = i_seq
        i_x = j_seq
        const_vect = flex.double(scatterers[i_h].site)- \
          flex.double(scatterers[i_x].site)
        self.table.append([i_x, i_h, const_vect, proxy.distance_ideal])
      if(scatterers[j_seq].element_symbol() == "H"):
        i_h = j_seq
        i_x = i_seq
        const_vect = flex.double(scatterers[i_h].site)- \
          flex.double(scatterers[i_x].site)
        self.table.append([i_x, i_h, const_vect, proxy.distance_ideal])

class manager(object):
  def __init__(self, xray_structure,
                     atom_attributes_list,
                     restraints_manager = None,
                     ias_xray_structure = None,
                     refinement_flags = None,
                     dbe_manager = None,
                     wilson_b = None,
                     tls_groups = None,
                     anomalous_scatterer_groups = None,
                     log = None):
    self.log = log
    self.restraints_manager = restraints_manager
    self.xray_structure = xray_structure.deep_copy_scatterers()
    self.xray_structure_initial = self.xray_structure.deep_copy_scatterers()
    self.atom_attributes_list = atom_attributes_list[:]
    self.refinement_flags = refinement_flags
    self.solvent_selection = self._solvent_selection()
    self.solvent_selection_ini = self._solvent_selection()
    self.wilson_b = wilson_b
    self.tls_groups = tls_groups
    if(anomalous_scatterer_groups is not None and
      len(anomalous_scatterer_groups) == 0):
      anomalous_scatterer_groups = None
    self.anomalous_scatterer_groups = anomalous_scatterer_groups
    # DBE related, need a real cleaning!
    self.dbe_manager = dbe_manager
    self.ias_xray_structure = ias_xray_structure
    self.use_dbe = False
    self.dbe_selection = None
    self.use_dbe_true_ = None
    self.use_dbe_false_ = None
    self.inflated = False
    self.dbe_added = False
    ###
    self.xh_connectivity = None
    if(restraints_manager is not None):
      if(xray_structure.hd_selection().count(True) > 0):
        self.xh_connectivity = xh_connectivity_table(
          geometry       = restraints_manager,
          xray_structure = xray_structure)
    ###
    if(self.refinement_flags is not None and [self.refinement_flags,
                                self.refinement_flags.adp_tls].count(None)==0):
       tlsos = tools.generate_tlsos(
                                selections     = self.refinement_flags.adp_tls,
                                xray_structure = self.xray_structure,
                                value          = 0.0)
       self.tls_groups.tlsos = tlsos

  def idealize_h(self, xh_bond_distance_deviation_limit):
    assert self.xh_connectivity is not None
    scatterers = self.xray_structure.scatterers()
    unit_cell = self.xray_structure.unit_cell()
    for bond in self.xh_connectivity.table:
      i_x = bond[0]
      i_h = bond[1]
      dist = unit_cell.distance(scatterers[i_x].site, scatterers[i_h].site)
      const_vect = bond[2]
      if(abs(bond[3] - dist) > xh_bond_distance_deviation_limit):
        scatterers[i_h].site = list(flex.double(scatterers[i_x].site) +
          flex.double(const_vect))

  def extract_ncs_groups(self):
    result = None
    if(self.restraints_manager.ncs_groups is not None):
      result = self.restraints_manager.ncs_groups.extract_ncs_groups(
        sites_cart = self.xray_structure.sites_cart())
    return result

  def deep_copy(self):
    new = manager(restraints_manager    = self.restraints_manager,
                  xray_structure        = self.xray_structure, # XXX not a deep copy
                  atom_attributes_list  = self.atom_attributes_list,
                  refinement_flags      = self.refinement_flags,
                  tls_groups            = self.tls_groups,
                  anomalous_scatterer_groups = self.anomalous_scatterer_groups,
                  log                   = self.log)
    selection = flex.bool(self.xray_structure.scatterers().size(), True)
    # XXX not a deep copy
    if(self.restraints_manager is not None):
       new.restraints_manager = mmtbx.restraints.manager(
            geometry      = self.restraints_manager.geometry.select(selection),
            ncs_groups    = self.restraints_manager.ncs_groups,
            normalization = self.restraints_manager.normalization)
       new.restraints_manager.geometry.pair_proxies(sites_cart =
                                              self.xray_structure.sites_cart())
    new.xray_structure       = self.xray_structure.deep_copy_scatterers()
    new.xray_structure_initial   = self.xray_structure_initial.deep_copy_scatterers()
    new.atom_attributes_list = self.atom_attributes_list[:]
    new.solvent_selection    = self.solvent_selection.deep_copy()
    if(self.refinement_flags is not None):
       new.refinement_flags = self.refinement_flags.select(selection)
    return new

  def add_dbe(self, fmodel=None, dbe_params=None, file_name=None,
                                                             build_only=False):
    if(self.dbe_manager is not None):
       self.remove_dbe()
       fmodel.update_xray_structure(xray_structure = self.xray_structure,
                                    update_f_calc = True)
    print >> self.log, ">>> Adding BDE.........."
    self.old_refinement_flags = None
    if not build_only: self.use_dbe = True
    self.dbe_manager = dbe.manager(
                    geometry             = self.restraints_manager.geometry,
                    atom_attributes_list = self.atom_attributes_list,
                    xray_structure       = self.xray_structure,
                    fmodel               = fmodel,
                    params               = dbe_params,
                    file_name            = file_name,
                    log                  = self.log)
    if(not build_only):
      self.ias_xray_structure = self.dbe_manager.ias_xray_structure
      dbe_size = self.ias_xray_structure.scatterers().size()
      tail = flex.bool(dbe_size, True)
      tail_false = flex.bool(dbe_size, False)
      self.dbe_selection = flex.bool(
                      self.xray_structure.scatterers().size(),False)
      self.dbe_selection.extend(tail)
      self.solvent_selection.extend(tail_false)
      self.xray_structure.concatenate_inplace(other = self.ias_xray_structure)
      print >> self.log, "Scattering dictionary for combined xray_structure:"
      self.xray_structure.scattering_type_registry().show()
      self.xray_structure_initial.concatenate_inplace(
                                           other = self.ias_xray_structure)
      if(self.refinement_flags is not None):
         self.old_refinement_flags = self.refinement_flags.deep_copy()
         self.refinement_flags.inflate(
              size = self.ias_xray_structure.scatterers().size(), flag = True)

         # refining simultaneously make convergence slower
         #self.refinement_flags.sites_individual[0].set_selected(
         #                                               self.dbe_selection, True)
         #self.refinement_flags.sites_individual[0].set_selected(
         #                                             ~self.dbe_selection, False)

         if 1: # XXX turned on
           self.refinement_flags.sites_individual[0].set_selected(
                                                          self.dbe_selection, False)
           self.refinement_flags.sites_individual[0].set_selected(
                                                        ~self.dbe_selection, True)

         if 1:
           # XXX 1) DO THIS WHEN INFLATING. THIS WILL NOT PROPERLY INFLATE IF THERE ARE ALT.CONF.
           # XXX 2) set to refine q fro only participating in IAS atoms, and not all
           self.refinement_flags.individual_occupancies = True
           self.refinement_flags.occupancies_individual = \
             [[i] for i in flex.bool(
               self.xray_structure.scatterers().size(), True).iselection()]

         if 0:
           self.refinement_flags.individual_occupancies = True
           self.refinement_flags.occupancies_individual = [flex.bool(
                                    self.xray_structure.scatterers().size(), True)]
           self.refinement_flags.occupancies_individual[0].set_selected(
                                                          ~self.dbe_selection, False)
           self.refinement_flags.occupancies_individual[0].set_selected(
                                        self.xray_structure.hd_selection(), True)


         #self.xray_structure.convert_to_anisotropic(selection = self.dbe_selection)
         self.refinement_flags.adp_individual_aniso[0].set_selected(
                                                     self.dbe_selection, False) # False
         self.refinement_flags.adp_individual_iso[0].set_selected(
                                                      self.dbe_selection, True) # True
         #occs = flex.double(self.xray_structure.scatterers().size(), 0.9)
         #self.xray_structure.scatterers().set_occupancies(occs, ~self.dbe_selection)
         # D9
         sel = self.xray_structure.scatterers().extract_scattering_types() == "D9"
         self.xray_structure.convert_to_anisotropic(selection = sel)
         self.refinement_flags.adp_individual_aniso[0].set_selected(sel, True)
         self.refinement_flags.adp_individual_iso[0].set_selected(sel, False)

  def remove_dbe(self):
    print >> self.log, ">>> Removing IAS..............."
    self.use_dbe = False
    if(self.dbe_manager is not None):
       self.dbe_manager = None
    if(self.old_refinement_flags is not None):
       self.refinement_flags = self.old_refinement_flags.deep_copy()
       self.old_refinement_flags = None
    if(self.dbe_selection is not None):
       self.xray_structure.select_inplace(selection = ~self.dbe_selection)
       n_non_ias = self.dbe_selection.count(False)
       self.solvent_selection = self.solvent_selection[:n_non_ias]
       self.dbe_selection = None
       self.xray_structure.scattering_type_registry().show()


  def write_dbe_pdb_file(self, out = None):
    if(out is None):
       out = sys.stdout
    print >> out, pdb.format_cryst1_record(
      crystal_symmetry = self.xray_structure.crystal_symmetry())
    print >> out, pdb.format_scale_records(
      unit_cell = self.xray_structure.crystal_symmetry().unit_cell())
    sites_cart = self.xray_structure.select(self.dbe_selection).sites_cart()
    scatterers = self.xray_structure.select(self.dbe_selection).scatterers()
    u_carts = scatterers.extract_u_cart_or_u_cart_plus_u_iso(self.xray_structure.unit_cell())
    u_isos      = self.xray_structure.extract_u_iso_or_u_equiv()
    for i_seq, sc in enumerate(self.xray_structure.select(self.dbe_selection).scatterers()):
        print >> out, pdb.format_atom_record(
                                    record_name = "HETATM",
                                    serial      = i_seq+1,
                                    name        = sc.label,
                                    resName     = "DBE",
                                    resSeq      = i_seq+1,
                                    site        = sites_cart[i_seq],
                                    occupancy   = sc.occupancy,
                                    tempFactor  = adptbx.u_as_b(u_isos[i_seq]),
                                    element     = sc.label)
        if(scatterers[i_seq].flags.use_u_aniso()):
           print >> out, pdb.format_anisou_record(
                                    serial      = i_seq+1,
                                    name        = sc.label,
                                    resName     = "DBE",
                                    resSeq      = i_seq+1,
                                    u_cart      = u_carts[i_seq],
                                    element     = sc.label)

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
    if(self.use_dbe and self.dbe_selection is not None and
       self.dbe_selection.count(True) > 0):
      sites_cart = sites_cart.select(~self.dbe_selection)
    return self.restraints_manager.energies_sites(
      sites_cart        = sites_cart,
      geometry_flags    = geometry_flags,
      compute_gradients = compute_gradients,
      gradients         = gradients,
      disable_asu_cache = disable_asu_cache)

  def _solvent_selection(self):
    labels = self.xray_structure.scatterers().extract_labels()
    res_name_tags = ["HOH","SOL","SOLV","WAT","DOD","TIP3"]
    atom_name_tags = ["O","OH2","H","H1","H2","D"]
    element_tags = ["O","H","D",""]
    result = flex.bool()
    for a in self.atom_attributes_list:
        element = (a.element).strip()
        resName = (a.resName).strip()
        name    = (a.name).strip()
        if((element in element_tags) and (name in atom_name_tags) and \
           (resName in res_name_tags)):
           result.append(True)
        else:
           result.append(False)
    return result

  def xray_structure_macromolecule(self):
    sel = self.solvent_selection
    if(self.use_dbe): sel = sel | self.dbe_selection
    result = self.xray_structure.select(~sel)
    return result

  def update(self, selection):
    self.xray_structure.select_inplace(selection)
    new_atom_attributes_list = []
    new_solvent_selection = flex.bool()
    for attr, solsel, sel in zip(self.atom_attributes_list,
                                 self.solvent_selection,
                                 selection):
        if(sel):
           new_atom_attributes_list.append(attr)
           new_solvent_selection.append(solsel)
    assert len(new_atom_attributes_list) == \
                                        self.xray_structure.scatterers().size()
    self.atom_attributes_list = new_atom_attributes_list
    self.solvent_selection = new_solvent_selection
    self.xray_structure.scattering_type_registry()
    if(self.restraints_manager is not None):
      self.restraints_manager = self.restraints_manager.select(
        selection=selection)
    if(self.refinement_flags is not None):
       self.refinement_flags = self.refinement_flags.select(selection)

  def number_of_ordered_solvent_molecules(self):
    return self.solvent_selection.count(True)

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
    self.update(selection = ~self.solvent_selection)

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

  def write_pdb_file(self, out, selection = None, xray_structure = None):
    utils.write_pdb_file(xray_structure       = self.xray_structure,
                         atom_attributes_list = self.atom_attributes_list,
                         selection            = selection,
                         out                  = out)

  def remove_atoms(self, atom_type         = None,
                         leave_only_labels = None):
    assert [atom_type, leave_only_labels].count(None) == 1
    if(atom_type is not None):
       remove_atoms_selection = (
         self.xray_structure.scatterers().extract_scattering_types()
           != atom_type)
       if (remove_atoms_selection.all_eq(True)): return
    if(leave_only_labels is not None):
       remove_atoms_selection = flex.bool(len(self.atom_attributes_list), False)
       for i_seq, atom in enumerate(self.atom_attributes_list):
           for item in leave_only_labels:
               if(item == atom.name.strip()):
                  remove_atoms_selection[i_seq] = True
    self.update(selection = remove_atoms_selection)

  def atoms_selection(self, scattering_type = None):
    scattering_types = \
                   self.xray_structure.scatterers().extract_scattering_types()
    return (scattering_type == scattering_types)


  def remove_atom_with_i_seqs(self, i_seq = None, i_seqs = None):
    assert [i_seq, i_seqs].count(None) == 1
    remove_atom_selection = flex.bool(len(self.atom_attributes_list), True)
    if(i_seq is not None):
       remove_atom_selection[i_seq] = False
    if(i_seqs is not None):
       for i_seq_i in i_seqs:
           remove_atom_selection[i_seq_i] = False
    self.update(selection = remove_atom_selection)

  def add_solvent(self, solvent_xray_structure,
                        solvent_selection,
                        atom_name    = "O",
                        residue_name = "HOH",
                        chain_id     = None):
    self.xray_structure = \
                        self.xray_structure.concatenate(solvent_xray_structure)
    self.refinement_flags.inflate(
      size        = solvent_xray_structure.scatterers().size(),
      use_u_iso   = solvent_xray_structure.use_u_iso(),
      use_u_aniso = solvent_xray_structure.use_u_aniso())
    new_atom_name = atom_name.strip()
    if(len(new_atom_name) < 4): new_atom_name = " " + new_atom_name
    while(len(new_atom_name) < 4): new_atom_name = new_atom_name+" "
    i_seq = 0
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
    self.solvent_selection = solvent_selection


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

  def geometry_minimization(self,
                            max_number_of_iterations = 100,
                            number_of_macro_cycles   = 100):
    raise RuntimeError("Not implemented.")
    if(max_number_of_iterations == 0 or number_of_macro_cycles == 0): return
    sso_start = stereochemistry_statistics(
                          xray_structure         = self.xray_structure,
                          restraints_manager     = self.restraints_manager,
                          use_dbe                = self.use_dbe,
                          dbe_selection          = self.dbe_selection,
                          text                   = "start")
    sites_cart = self.xray_structure.sites_cart()
    first_target_value = None
    for macro_cycles in xrange(1,number_of_macro_cycles+1):
        minimized = cctbx_geometry_restraints_lbfgs(
          sites_cart                  = sites_cart,
          geometry_restraints_manager = self.restraints_manager.geometry,
          lbfgs_termination_params    = scitbx.lbfgs.termination_parameters(
                                    max_iterations = max_number_of_iterations))
        if(first_target_value is None):
           first_target_value = minimized.first_target_value
    self.xray_structure = \
                 self.xray_structure.replace_sites_cart(new_sites = sites_cart)
    sso_end = stereochemistry_statistics(
                          xray_structure         = self.xray_structure,
                          restraints_manager     = self.restraints_manager,
                          use_dbe                = self.use_dbe,
                          dbe_selection          = self.dbe_selection,
                          text                   = "final")
    assert approx_equal(first_target_value, sso_start.target)
    assert approx_equal(minimized.final_target_value, sso_end.target)
    sso_start.show(out = self.log)
    sso_end.show(out = self.log)

  def geometry_statistics(self):
    sites_cart = self.xray_structure.sites_cart()
    if(self.use_dbe): sites_cart = sites_cart.select(~self.dbe_selection)
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

  def energies_adp(self, iso_restraints, compute_gradients):
    assert self.refinement_flags is not None
    n_aniso = 0
    for sel in self.refinement_flags.adp_individual_aniso:
      n_aniso += sel.count(True)
    if(n_aniso == 0):
      energies_adp_iso = self.restraints_manager.energies_adp_iso(
        xray_structure    = self.xray_structure,
        parameters        = iso_restraints,
        use_u_local_only  = iso_restraints.use_u_local_only,
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
