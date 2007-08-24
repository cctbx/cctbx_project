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


time_model_show = 0.0


class manager(object):
  def __init__(self, restraints_manager,
                     restraints_manager_ini,
                     xray_structure,
                     atom_attributes_list,
                     ias_xray_structure = None,
                     refinement_flags = None,
                     dbe_manager = None,
                     wilson_b = None,
                     tls_groups = None,
                     anomalous_scatterer_groups = None,
                     log = None):
    self.log = log
    self.restraints_manager = restraints_manager
    self.restraints_manager_ini = restraints_manager_ini
    self.xray_structure = xray_structure.deep_copy_scatterers()
    self.xray_structure_ini = self.xray_structure.deep_copy_scatterers()
    self.crystal_symmetry = self.xray_structure.crystal_symmetry()
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
    # DBE related
    self.dbe_manager = dbe_manager
    self.ias_xray_structure = ias_xray_structure
    self.use_dbe = False
    self.dbe_selection = None
    self.use_dbe_true_ = None
    self.use_dbe_false_ = None
    self.inflated = False
    self.dbe_added = False

    if(self.refinement_flags is not None and [self.refinement_flags,
                                self.refinement_flags.adp_tls].count(None)==0):
       tlsos = tools.generate_tlsos(
                                selections     = self.refinement_flags.adp_tls,
                                xray_structure = self.xray_structure,
                                value          = 0.0)
       self.tls_groups.tlsos = tlsos

  def extract_ncs_groups(self):
    result = None
    if(self.restraints_manager.ncs_groups is not None):
      result = self.restraints_manager.ncs_groups.extract_ncs_groups(
        sites_cart = self.xray_structure.sites_cart())
    return result

  def deep_copy(self):
    new = manager(restraints_manager    = self.restraints_manager,
                  restraints_manager_ini= self.restraints_manager_ini,
                  xray_structure        = self.xray_structure,
                  atom_attributes_list  = self.atom_attributes_list,
                  refinement_flags      = self.refinement_flags,
                  tls_groups            = self.tls_groups,
                  anomalous_scatterer_groups = self.anomalous_scatterer_groups,
                  log                   = self.log)
    selection = flex.bool(self.xray_structure.scatterers().size(), True)
    # XXX not a deep copy
    if(self.restraints_manager is not None):
       new.restraints_manager_ini = self.restraints_manager_ini
       new.restraints_manager = mmtbx.restraints.manager(
            geometry      = self.restraints_manager.geometry.select(selection),
            ncs_groups    = self.restraints_manager.ncs_groups,
            normalization = self.restraints_manager.normalization)
       new.restraints_manager.geometry.pair_proxies(sites_cart =
                                              self.xray_structure.sites_cart())
    new.xray_structure       = self.xray_structure.deep_copy_scatterers()
    new.xray_structure_ini   = self.xray_structure_ini.deep_copy_scatterers()
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
      self.xray_structure_ini.concatenate_inplace(
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
           self.refinement_flags.individual_occupancies = True
           self.refinement_flags.occupancies_individual = [flex.bool(
                                self.xray_structure.scatterers().size(), True)]

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
       self.xray_structure_ini.select_inplace(selection = ~self.dbe_selection)
       n_non_ias = self.dbe_selection.count(False)
       self.solvent_selection = self.solvent_selection[:n_non_ias]
       self.dbe_selection = None
       self.xray_structure.scattering_type_registry().show()

  def write_dbe_pdb_file(self, out = None):
    if(out is None):
       out = sys.stdout
    print >> out, pdb.format_cryst1_record(
                                      crystal_symmetry = self.crystal_symmetry)
    print >> out, pdb.format_scale_records(
                                 unit_cell = self.crystal_symmetry.unit_cell())
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
                                        sites_cart           = None,
                                        geometry_flags       = None,
                                        compute_gradients    = False,
                                        gradients            = None,
                                        disable_asu_cache    = False,
                                        lock_for_line_search = False):
    if(sites_cart is None): sites_cart = self.xray_structure.sites_cart()
    if(self.use_dbe and self.dbe_selection is not None and
                                           self.dbe_selection.count(True) > 0):
       sites_cart = sites_cart.select(~self.dbe_selection)
    return self.restraints_manager.energies_sites(
                                   sites_cart           = sites_cart,
                                   geometry_flags       = geometry_flags,
                                   compute_gradients    = compute_gradients,
                                   gradients            = gradients,
                                   disable_asu_cache    = disable_asu_cache,
                                   lock_for_line_search = lock_for_line_search)

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
    next = "| %5d: %8.3f %6.2f %5s %4s %4d <> %8.3f %6.2f %5s %4s %4d   |"
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
                        chain_id     = "S"):
    resSeqs = flex.int()
    for aa in self.atom_attributes_list:
        if(aa.chainID.strip() == chain_id.strip()):
           resSeqs.append(int(aa.resSeq))
    if(resSeqs.size() > 0):
       i_seq = flex.max(resSeqs)
    else:
       i_seq = 0
    self.xray_structure = \
                        self.xray_structure.concatenate(solvent_xray_structure)
    self.refinement_flags.inflate(
                             size = solvent_xray_structure.scatterers().size())
    new_atom_name = atom_name.strip()
    if(len(new_atom_name) < 4): new_atom_name = " " + new_atom_name
    while(len(new_atom_name) < 4): new_atom_name = new_atom_name+" "
    #i_seq = solvent_xray_structure.scatterers().size()
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
    if(max_number_of_iterations == 0 or number_of_macro_cycles == 0): return
    sso_start = stereochemistry_statistics(
                          xray_structure         = self.xray_structure,
                          xray_structure_ref     = self.xray_structure_ini,
                          restraints_manager     = self.restraints_manager,
                          restraints_manager_ref = self.restraints_manager_ini,
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
                          xray_structure_ref     = self.xray_structure_ini,
                          restraints_manager     = self.restraints_manager,
                          restraints_manager_ref = self.restraints_manager_ini,
                          use_dbe                = self.use_dbe,
                          dbe_selection          = self.dbe_selection,
                          text                   = "final")
    assert approx_equal(first_target_value, sso_start.target)
    assert approx_equal(minimized.final_target_value, sso_end.target)
    sso_start.show(out = self.log)
    sso_end.show(out = self.log)


  def geometry_statistics(self, other = None,
                                show = False,
                                text = "",
                                short = True):
    global time_model_show
    timer = user_plus_sys_time()
    if(other is not None):
       stereochemistry_statistics_obj = stereochemistry_statistics(
                             xray_structure         = self.xray_structure,
                             xray_structure_ref     = other.xray_structure,
                             restraints_manager     = self.restraints_manager,
                             restraints_manager_ref = other.restraints_manager,
                             use_dbe                = self.use_dbe,
                             dbe_selection          = self.dbe_selection,
                             text                   = text)
    else:
       stereochemistry_statistics_obj = stereochemistry_statistics(
                          xray_structure         = self.xray_structure,
                          xray_structure_ref     = self.xray_structure_ini,
                          restraints_manager     = self.restraints_manager,
                          restraints_manager_ref = self.restraints_manager_ini,
                          use_dbe                = self.use_dbe,
                          dbe_selection          = self.dbe_selection,
                          text                   = text)
    if(show): stereochemistry_statistics_obj.show(out = self.log,short = short)
    time_model_show += timer.elapsed()
    return stereochemistry_statistics_obj

  def geometry_statistics_simple(self, show = False, prefix = "", out = None):
    global time_model_show
    timer = user_plus_sys_time()
    if(out == None): out = self.log
    if(self.dbe_selection is not None and self.dbe_selection.count(True) > 0):
       xray_structure = self.xray_structure.select(~self.dbe_selection)
       solvent_selection = self.solvent_selection.select(~self.dbe_selection)
    else:
       xray_structure = self.xray_structure
       solvent_selection = self.solvent_selection
    geometry_statistics_obj = geometry_statistics(
                    xray_structure     = xray_structure.deep_copy_scatterers(),
                    restraints_manager = self.restraints_manager,
                    prefix             = prefix,
                    out                = out)
    if(show): geometry_statistics_obj.show()
    time_model_show += timer.elapsed()
    return geometry_statistics_obj

  def adp_statistics(self, show     = False,
                           prefix   = "",
                           padded   = False,
                           out      = None):
    global time_model_show
    timer = user_plus_sys_time()
    if(out == None): out = self.log
    if(self.dbe_selection is not None and self.dbe_selection.count(True) > 0):
       xray_structure = self.xray_structure.select(~self.dbe_selection)
       solvent_selection = self.solvent_selection.select(~self.dbe_selection)
    else:
       xray_structure = self.xray_structure
       solvent_selection = self.solvent_selection
    adp_statistics_obj = adp_statistics(
                     xray_structure    = xray_structure.deep_copy_scatterers(),
                     solvent_selection = solvent_selection.deep_copy(),
                     wilson_b          = self.wilson_b,
                     prefix            = prefix,
                     padded            = padded,
                     out               = out)
    if(show): adp_statistics_obj.show()
    time_model_show += timer.elapsed()
    return adp_statistics_obj

class geometry_statistics(object):
  def __init__(self,
               xray_structure,
               restraints_manager,
               prefix        = "",
               out           = None):
    if(out is None): out = sys.stdout
    self.out = out
    self.prefix = prefix
    self.xray_structure = xray_structure
    self.restraints_manager = restraints_manager
    self.a_target, self.a_mean, self.a_max, self.a_min = 0.,0.,0.,0.
    self.b_target, self.b_mean, self.b_max, self.b_min = 0.,0.,0.,0.
    self.c_target, self.c_mean, self.c_max, self.c_min = 0.,0.,0.,0.
    self.d_target, self.d_mean, self.d_max, self.d_min = 0.,0.,0.,0.
    self.p_target, self.p_mean, self.p_max, self.p_min = 0.,0.,0.,0.
    self.n_target, self.n_mean, self.n_max, self.n_min = 0.,0.,0.,0.
    self.target = 0.
    energies_sites = restraints_manager.energies_sites(
                               sites_cart        = xray_structure.sites_cart(),
                               compute_gradients = False)
    self.esg      = energies_sites.geometry
    self.b_target = self.esg.bond_residual_sum
    self.a_target = self.esg.angle_residual_sum
    self.d_target = self.esg.dihedral_residual_sum
    self.c_target = self.esg.chirality_residual_sum
    self.n_target = self.esg.nonbonded_residual_sum
    self.p_target = self.esg.planarity_residual_sum
    self.b_min, self.b_max, self.b_mean = self.esg.bond_deviations()
    self.n_min, self.n_max, self.n_mean = self.esg.nonbonded_deviations()
    self.a_min, self.a_max, self.a_mean = self.esg.angle_deviations()
    self.d_min, self.d_max, self.d_mean = self.esg.dihedral_deviations()
    self.c_min, self.c_max, self.c_mean = self.esg.chirality_deviations()
    self.p_min, self.p_max, self.p_mean = self.esg.planarity_deviations()
    self.target               = self.esg.target # XXX not normalized, I presume ?!
    self.number_of_restraints = self.esg.number_of_restraints # XXX not normalized !

  def show_bond_angle_nonbonded_histogram(self, n_slots = 10, out = None):
    if(out is None): out = self.out
    prefix = self.prefix
    rmg = self.restraints_manager.geometry
    bond_deltas = flex.abs(geometry_restraints.bond_deltas(
                         sites_cart         = self.xray_structure.sites_cart(),
                         sorted_asu_proxies = rmg.pair_proxies().bond_proxies))
    angle_deltas = geometry_restraints.angle_deltas(
                                 sites_cart = self.xray_structure.sites_cart(),
                                 proxies    = rmg.angle_proxies)
    nonbonded_distances = self.esg.nonbonded_distances()
    h_1 = flex.histogram(data = flex.abs(bond_deltas),  n_slots = n_slots)
    h_2 = flex.histogram(data = flex.abs(angle_deltas), n_slots = n_slots)
    h_3 = flex.histogram(data = flex.abs(nonbonded_distances), n_slots=n_slots)
    print >> out, prefix+"|------------------------------------------------------------|"
    print >> out, prefix+"|        Histogram of deviations from ideal values for       |"
    print >> out, prefix+"|Bonds             |Angles                |Nonbonded contacts|"
    format = "|%5.3f-%5.3f:%6d|%7.3f-%7.3f:%6d|%5.3f-%5.3f:%6d|"
    lc_1 = h_1.data_min()
    lc_2 = h_2.data_min()
    lc_3 = h_3.data_min()
    s_1 = enumerate(h_1.slots())
    s_2 = enumerate(h_2.slots())
    s_3 = enumerate(h_3.slots())
    for (i_1,n_1),(i_2,n_2),(i_3,n_3) in zip(s_1, s_2, s_3):
        hc_1 = h_1.data_min() + h_1.slot_width() * (i_1+1)
        hc_2 = h_2.data_min() + h_2.slot_width() * (i_2+1)
        hc_3 = h_3.data_min() + h_3.slot_width() * (i_3+1)
        print >> out, prefix+format % (lc_1,hc_1,n_1, lc_2,hc_2,n_2, lc_3,hc_3,n_3)
        lc_1 = hc_1
        lc_2 = hc_2
        lc_3 = hc_3
    out.flush()
    print >> out, prefix+"|------------------------------------------------------------|"
    out.flush()

  def show(self, out = None):
    if(out is None): out = self.out
    prefix = self.prefix
    fmt1 = "|%7.3f%8.3f%7.3f|%12.3f|              |"
    fmt2 = "|%7.3f%8.3f%7.3f|%12.3f|%12.3f  |"
    print >> out, prefix+"|-Geometry statistics----------------------------------------|"
    print >> out, prefix+"|Type     | Deviation from ideal |   Targets  |Target (sum)  |"
    print >> out, prefix+"|         |  rmsd     max    min |            |              |"
    print >> out, prefix+"|bond     "+fmt1%(self.b_mean,self.b_max,self.b_min,self.b_target)
    print >> out, prefix+"|angle    "+fmt1%(self.a_mean,self.a_max,self.a_min,self.a_target)
    print >> out, prefix+"|chirality"+fmt2%(self.c_mean,self.c_max,self.c_min,self.c_target,self.target)
    print >> out, prefix+"|planarity"+fmt1%(self.p_mean,self.p_max,self.p_min,self.p_target)
    print >> out, prefix+"|dihedral "+fmt1%(self.d_mean,self.d_max,self.d_min,self.d_target)
    print >> out, prefix+"|nonbonded"+fmt1%(self.n_mean,self.n_max,self.n_min,self.n_target)
    print >> out, prefix+"|------------------------------------------------------------|"
    self.show_bond_angle_nonbonded_histogram(out = out)
    out.flush()

class adp_statistics(object):
  def __init__(self,
               xray_structure,
               solvent_selection,
               wilson_b      = None,
               prefix        = "",
               padded        = False,
               out           = None):
    if(out is None): out = sys.stdout
    self.out = out
    self.wilson_b = wilson_b
    self.prefix = prefix
    self.padded = padded
    eps = math.pi**2*8
    scattering_types = xray_structure.scatterers().extract_scattering_types()
    hd_selection = (scattering_types == "D") | (scattering_types == "H")
    # XXX too generous but safe and quick
    xs_all = xray_structure.deep_copy_scatterers()
    xs_hyd = xs_all.select(hd_selection)
    xs_all = xs_all.select(~hd_selection)
    solvent_selection = solvent_selection.select(~hd_selection)
    xs_sol = xs_all.select(solvent_selection)
    xs_mac = xs_all.select(~solvent_selection)
    #
    self.b_isos_all = xs_all.extract_u_iso_or_u_equiv() * eps
    b_isos_sol = xs_sol.extract_u_iso_or_u_equiv() * eps
    b_isos_mac = xs_mac.extract_u_iso_or_u_equiv() * eps
    b_isos_hyd = xs_hyd.extract_u_iso_or_u_equiv() * eps
    #
    self.b_min_all, self.b_max_all, self.b_mean_all = \
                              self._min_max_mean_default(self.b_isos_all, None)
    self.b_min_sol, self.b_max_sol, self.b_mean_sol = \
                                   self._min_max_mean_default(b_isos_sol, None)
    self.b_min_mac, self.b_max_mac, self.b_mean_mac = \
                                   self._min_max_mean_default(b_isos_mac, None)
    self.b_min_hyd, self.b_max_hyd, self.b_mean_hyd = \
                                   self._min_max_mean_default(b_isos_hyd, None)
    #
    sel_aniso_all = xs_all.use_u_aniso()
    sel_aniso_sol = xs_sol.use_u_aniso()
    sel_aniso_mac = xs_mac.use_u_aniso()
    sel_aniso_hyd = xs_hyd.use_u_aniso()
    sel_iso_all   = xs_all.use_u_iso()
    sel_iso_sol   = xs_sol.use_u_iso()
    sel_iso_mac   = xs_mac.use_u_iso()
    sel_iso_hyd   = xs_hyd.use_u_iso()
    self.n_aniso_all = sel_aniso_all.count(True)
    self.n_aniso_sol = sel_aniso_sol.count(True)
    self.n_aniso_mac = sel_aniso_mac.count(True)
    self.n_aniso_hyd = sel_aniso_hyd.count(True)
    self.n_iso_all = sel_iso_all.count(True)
    self.n_iso_sol = sel_iso_sol.count(True)
    self.n_iso_mac = sel_iso_mac.count(True)
    self.n_iso_hyd = sel_iso_hyd.count(True)
    #
    uc = xs_all.unit_cell()
    self.a_all = xs_all.scatterers().anisotropy(unit_cell = uc).select(sel_aniso_all)
    a_sol      = xs_sol.scatterers().anisotropy(unit_cell = uc).select(sel_aniso_sol)
    a_mac      = xs_mac.scatterers().anisotropy(unit_cell = uc).select(sel_aniso_mac)
    a_hyd      = xs_hyd.scatterers().anisotropy(unit_cell = uc).select(sel_aniso_hyd)
    self.a_min_all, self.a_max_all, self.a_mean_all = \
                                   self._min_max_mean_default(self.a_all, None)
    self.a_min_sol, self.a_max_sol, self.a_mean_sol = \
                                        self._min_max_mean_default(a_sol, None)
    self.a_min_mac, self.a_max_mac, self.a_mean_mac = \
                                        self._min_max_mean_default(a_mac, None)
    self.a_min_hyd, self.a_max_hyd, self.a_mean_hyd = \
                                        self._min_max_mean_default(a_hyd, None)

  def _min_max_mean_default(self, values, default):
    return flex.min_default(values, default), \
           flex.max_default(values, default), \
           flex.mean_default(values, default)

  def _fmtl(self,a,b,c,d,e,f,g,h, pad):
    n = []
    for item in [a,b,c,d,e,f,g,h]:
      if(item is None): item = str(item)
      elif(str(item).count(".")): item = str("%8.2f"%item).strip()
      else: item = str("%6d"%item).strip()
      n.append(item)
    fmt = "%-6s %-6s %-8s%-8s%-8s %-6s%-7s%-7s"+pad+"|"
    return fmt%(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7])

  def _histogram(self, values, out, pad_l,pad_r):
    result = []
    histogram = flex.histogram(data = values, n_slots = 10)
    # n_slots must be even
    low_cutoff_1 = histogram.data_min()
    for (i_1,n_1) in enumerate(histogram.slots()):
      high_cutoff_1 = histogram.data_min() + histogram.slot_width()*(i_1+1)
      result.append([i_1,low_cutoff_1,high_cutoff_1,n_1])
      low_cutoff_1 = high_cutoff_1
    result1, result2 = result[:len(result)/2], result[len(result)/2:]
    fmt = "|"+pad_l+"   %d:%10.3f -%8.3f:%5d   |   %d:%10.3f -%8.3f:%5d    "+pad_r+"|"
    for r1, r2 in zip(result1, result2):
      print >> out, self.prefix+fmt%(r1[0],r1[1],r1[2],r1[3], r2[0],r2[1],r2[2],r2[3])
    print >> out, self.prefix+"|"+pad_l+"                            =>continue=>                              "+pad_r+"|"

  def show(self, out = None, padded = None):
    prefix = self.prefix
    if(padded is None): padded = self.padded
    if(out is None): out = self.out
    pad_l, pad_r = "", ""
    if(padded):
       pad_l = " "*3
       pad_r = " "*4
    #
    if(self.wilson_b is not None):
       bw = str("%9.3f"%self.wilson_b).strip()
       s1 = "|-ADP statistics (Wilson B = "+bw+")"
    else:
       s1 = "|-ADP statistics"
    if(padded): s1 = s1 + "-"*( 79-len(s1+"|") ) + "|"
    else: s1 = s1 + "-"*( 72-len(s1+"|") ) + "|"
    print >> out, prefix+s1
    #
    print >> out, prefix+"|"+pad_l+" Atom    | Number of   | Isotropic or equivalent| Anisotropy lmin/max "+pad_r+"|"
    print >> out, prefix+"|"+pad_l+" type    |iso    aniso | min     max     mean   | min   max    mean   "+pad_r+"|"
    print >> out, prefix+"|"+pad_l+" - - - - |- - - - - - -| - - - - - - - - - - - -| - - - - - - - - - - "+pad_r+"|"
    #
    print >> out, prefix+"|"+pad_l+" Solv+Mac: "+self._fmtl(self.n_iso_all,self.n_aniso_all,self.b_min_all,self.b_max_all,self.b_mean_all,self.a_min_all,self.a_max_all,self.a_mean_all,pad_r)
    print >> out, prefix+"|"+pad_l+" Sol.    : "+self._fmtl(self.n_iso_sol,self.n_aniso_sol,self.b_min_sol,self.b_max_sol,self.b_mean_sol,self.a_min_sol,self.a_max_sol,self.a_mean_sol,pad_r)
    print >> out, prefix+"|"+pad_l+" Mac.    : "+self._fmtl(self.n_iso_mac,self.n_aniso_mac,self.b_min_mac,self.b_max_mac,self.b_mean_mac,self.a_min_mac,self.a_max_mac,self.a_mean_mac,pad_r)
    print >> out, prefix+"|"+pad_l+" Hyd.    : "+self._fmtl(self.n_iso_hyd,self.n_aniso_hyd,self.b_min_hyd,self.b_max_hyd,self.b_mean_hyd,self.a_min_hyd,self.a_max_hyd,self.a_mean_hyd,pad_r)
    print >> out, prefix+"|"+pad_l+" "+"- "*34+pad_r+" |"
    #
    name = "|"+pad_l+"    Distribution of isotropic (or equivalent) ADP for non-H atoms:    "+pad_r+"|"
    print >> out, self.prefix+name
    print >> out, self.prefix+"|"+pad_l+" Bin#      value range     #atoms | Bin#      value range     #atoms  "+pad_r+"|"
    self._histogram(values = self.b_isos_all, out = out, pad_l=pad_l,pad_r=pad_r)
    #
    if([self.a_min_all,self.a_max_all,self.a_mean_all].count(None) == 0 and
       abs(self.a_min_all+self.a_max_all+self.a_mean_all-3) > 0.001):
       print >> out, prefix+"|"+pad_l+" "+"- "*34+pad_r+" |"
       name = "|"+pad_l+"                     Distribution of anisotropy:                      "+pad_r+"|"
       print >> out, self.prefix+name
       print >> out, self.prefix+"|"+pad_l+" Bin#      value range     #atoms | Bin#      value range     #atoms  "+pad_r+"|"
       self._histogram(values = self.a_all, out = out, pad_l=pad_l,pad_r=pad_r)
    #
    if(padded): print >> out, prefix+"|"+"-"*77+"|"
    else: print >> out, prefix+"|"+"-"*70+"|"
    out.flush()

# XXX This below will retire soon: too complex, hard to maintain
class stereochemistry_statistics(object):
  def __init__(self,
               xray_structure,
               xray_structure_ref,
               restraints_manager,
               restraints_manager_ref,
               use_dbe,
               dbe_selection,
               text=""):
    adopt_init_args(self, locals())
    self.a_target, self.a_mean, self.a_max, self.a_min = 0.,0.,0.,0.
    self.b_target, self.b_mean, self.b_max, self.b_min = 0.,0.,0.,0.
    self.c_target, self.c_mean, self.c_max, self.c_min = 0.,0.,0.,0.
    self.d_target, self.d_mean, self.d_max, self.d_min = 0.,0.,0.,0.
    self.p_target, self.p_mean, self.p_max, self.p_min = 0.,0.,0.,0.
    self.r_target, self.r_mean, self.r_max, self.r_min = 0.,0.,0.,0.
    self.target = 0.
    self.delta_model_mean = 0.
    self.delta_model_max = 0.
    self.delta_model_min = 0.
    if(self.use_dbe):
       self.energies_sites_ref = self.restraints_manager_ref.energies_sites(
         sites_cart        = self.xray_structure_ref.sites_cart().select(~self.dbe_selection),
         compute_gradients = True)
       self.energies_sites = self.restraints_manager.energies_sites(
         sites_cart        = self.xray_structure.sites_cart().select(~self.dbe_selection),
         compute_gradients = True)
    else:
       self.energies_sites_ref = self.restraints_manager_ref.energies_sites(
                         sites_cart        = self.xray_structure_ref.sites_cart(),
                         compute_gradients = True)
       self.energies_sites = self.restraints_manager.energies_sites(
                             sites_cart        = self.xray_structure.sites_cart(),
                             compute_gradients = True)
    self.esg      = self.energies_sites.geometry
    self.esg_ref  = self.energies_sites_ref.geometry
    self.b_target = self.esg.bond_residual_sum
    self.a_target = self.esg.angle_residual_sum
    self.d_target = self.esg.dihedral_residual_sum
    self.c_target = self.esg.chirality_residual_sum
    self.r_target = self.esg.nonbonded_residual_sum
    self.p_target = self.esg.planarity_residual_sum
    self.b_min, self.b_max, self.b_mean = self.esg.bond_deviations()
    self.r_min, self.r_max, self.r_mean = self.esg.nonbonded_deviations()
    self.a_min, self.a_max, self.a_mean = self.esg.angle_deviations()
    self.d_min, self.d_max, self.d_mean = self.esg.dihedral_deviations()
    self.c_min, self.c_max, self.c_mean = self.esg.chirality_deviations()
    self.p_min, self.p_max, self.p_mean = self.esg.planarity_deviations()
    self.target               = self.esg.target
    self.gradients            = self.esg.gradients
    self.number_of_restraints = self.esg.number_of_restraints
    if(self.number_of_restraints > 0):
       self.target_normalized    = self.target / self.number_of_restraints
       self.gradients_normalized = \
                                self.gradients * (1./self.number_of_restraints)
    self.get_model_diff()

  def get_model_diff(self):
    if(self.xray_structure.scatterers().size() != \
                                  self.xray_structure_ref.scatterers().size()):
       self.delta_model_mean = None
       self.delta_model_max  = None
       self.delta_model_min  = None
    else:
       array_of_distances_between_each_atom = flex.sqrt(
                 self.xray_structure.difference_vectors_cart(
                                                self.xray_structure_ref).dot())
       self.delta_model_mean = flex.mean_default(
                                       array_of_distances_between_each_atom, 0)
       self.delta_model_max = flex.max_default(
                                       array_of_distances_between_each_atom, 0)
       self.delta_model_min = flex.min_default(
                                       array_of_distances_between_each_atom, 0)

  def show_bond_angle_histogram(self, n_slots = 10, out=None):
    if (out is None): out = sys.stdout
    rmg_1 = self.restraints_manager_ref.geometry
    rmg_2 = self.restraints_manager.geometry
    bond_deltas_1 = flex.abs(geometry_restraints.bond_deltas(
                     sites_cart         = self.xray_structure_ref.sites_cart(),
                     sorted_asu_proxies = rmg_1.pair_proxies().bond_proxies))
    bond_deltas_2 = flex.abs(geometry_restraints.bond_deltas(
                        sites_cart         = self.xray_structure.sites_cart(),
                        sorted_asu_proxies = rmg_2.pair_proxies().bond_proxies))
    angle_deltas_1 = geometry_restraints.angle_deltas(
                             sites_cart = self.xray_structure_ref.sites_cart(),
                             proxies    = rmg_1.angle_proxies)
    angle_deltas_2 = geometry_restraints.angle_deltas(
                                 sites_cart = self.xray_structure.sites_cart(),
                                 proxies    = rmg_2.angle_proxies)
    nonbonded_distances_1 = self.esg_ref.nonbonded_distances()
    nonbonded_distances_2 = self.esg.nonbonded_distances()
    h_1 = flex.histogram(data    = flex.abs(bond_deltas_1),
                         n_slots = n_slots)
    h_2 = flex.histogram(other   = h_1,
                         data    = flex.abs(bond_deltas_2))
    h_3 = flex.histogram(data    = flex.abs(angle_deltas_1),
                         n_slots = n_slots)
    h_4 = flex.histogram(other   = h_3,
                         data    = flex.abs(angle_deltas_2))
    h_5 = flex.histogram(data    = flex.abs(nonbonded_distances_1),
                         n_slots = n_slots)
    h_6 = flex.histogram(other   = h_5,
                         data    = flex.abs(nonbonded_distances_2))
    print >> out, "|-----------------------------------------------------------------------------|"
    print >> out, "|                 Histograms for start / current models of                    |"
    print >> out, "|                       |                           |                         |"
    print >> out, "|        deviations from ideal values for           |                         |"
    print >> out, "|                       |                           |                         |"
    print >> out, "| bonds                 |angles                     |nonbonded contacts       |"
    print >> out, "|                       |                           |                         |"
    show_6_histograms(h_1, h_2, h_3, h_4, h_5, h_6, n_slots = n_slots, out=out)
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def show(self, out=None, short=True):
    if (out is None): out = sys.stdout
    line_len = len("|"+self.text+"|")
    fill_len = 80 - line_len-1
    print >> out
    print >> out, "|"+self.text+"-"*(fill_len)+"|"
    print >> out, "|Type| Deviation from ideal |   Targets  ||Target (sum)|| Deviation of start  |"
    print >> out, "|    |  rmsd     max    min |            ||            || model from current  |"
    print >> out, "|bond|%7.3f%8.3f%7.3f|%12.3f||            ||  mean   max    min  |"%\
        (self.b_mean,self.b_max,self.b_min,self.b_target)
    print >> out, "|angl|%7.3f%8.3f%7.3f|%12.3f||            ||                     |"%\
        (self.a_mean,self.a_max,self.a_min,self.a_target)
    if([self.delta_model_mean,self.delta_model_max,self.delta_model_min].count(None) == 0):
       print >> out, "|chir|%7.3f%8.3f%7.3f|%12.3f||%12.3f||%7.3f%7.3f%7.3f|"%\
              (self.c_mean,self.c_max,self.c_min,self.c_target,self.target,\
               self.delta_model_mean,self.delta_model_max,self.delta_model_min)
    else:
       print >> out, "|chir|%7.3f%8.3f%7.3f|%12.3f||%12.3f||      undefined      |"%\
                  (self.c_mean,self.c_max,self.c_min,self.c_target,self.target)
    print >> out, "|plan|%7.3f%8.3f%7.3f|%12.3f||            ||                     |"%\
        (self.p_mean,self.p_max,self.p_min,self.p_target)
    print >> out, "|dihe|%7.3f%8.3f%7.3f|%12.3f||            ||                     |"%\
        (self.d_mean,self.d_max,self.d_min,self.d_target)
    print >> out, "|repu|%7.3f%8.3f%7.3f|%12.3f||            ||                     |"%\
        (self.r_mean,self.r_max,self.r_min,self.r_target)
    print >> out, "|"+"-"*77+"|"
    ng = self.energies_sites.ncs_groups
    if (ng is not None and len(ng.rms_with_respect_to_averages) > 0):
      print >> out, "| %-75s |" % (
        "RMS differences with respect to the average:")
      for i_group,rms in enumerate(ng.rms_with_respect_to_averages):
        l = ["|   NCS group %2d:" % (i_group+1)]
        if (rms is None):
          l.append("None")
        else:
          l.append("size = %2d min = %9.6f max = %9.6f mean = %9.6f" % (
            rms.size(), flex.min(rms), flex.max(rms), flex.mean(rms)))
        l = " ".join(l)
        if   (len(l) <= 77): l += " "*(77-len(l)) + " |"
        elif (len(l) == 78): l += "|"
        print >> out, l
      print >> out, "|"+"-"*77+"|"
    if(not short):
       self.show_bond_angle_histogram(out=out)
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

def show_6_histograms(h_1, h_2, h_3, h_4, h_5, h_6, n_slots, out):
  format = "|%4.2f-%4.2f:%6d/%6d|%6.2f-%6.2f:%6d/%6d|%4.2f-%4.2f:%6d/%6d  |"
  lc_1 = h_1.data_min()
  lc_2 = h_2.data_min()
  lc_3 = h_3.data_min()
  lc_4 = h_4.data_min()
  lc_5 = h_5.data_min()
  lc_6 = h_6.data_min()
  s_1 = enumerate(h_1.slots())
  s_2 = enumerate(h_2.slots())
  s_3 = enumerate(h_3.slots())
  s_4 = enumerate(h_4.slots())
  s_5 = enumerate(h_5.slots())
  s_6 = enumerate(h_6.slots())
  for (i_1,n_1),(i_2,n_2),(i_3,n_3),(i_4,n_4),(i_5,n_5),(i_6,n_6) in zip(s_1,
                                                          s_2,s_3,s_4,s_5,s_6):
    hc_1 = h_1.data_min() + h_1.slot_width() * (i_1+1)
    hc_2 = h_2.data_min() + h_2.slot_width() * (i_2+1)
    hc_3 = h_3.data_min() + h_3.slot_width() * (i_3+1)
    hc_4 = h_4.data_min() + h_4.slot_width() * (i_4+1)
    hc_5 = h_5.data_min() + h_5.slot_width() * (i_5+1)
    hc_6 = h_6.data_min() + h_6.slot_width() * (i_6+1)
    assert lc_1 == lc_2 and hc_1 == hc_2
    assert lc_3 == lc_4 and hc_3 == hc_4
    assert lc_5 == lc_6 and hc_5 == hc_6
    outs = (lc_1,hc_1,n_1,n_2, lc_3,hc_3,n_3,n_4, lc_5,hc_5,n_5,n_6)
    print >> out, format % outs
    lc_1 = hc_1
    lc_2 = hc_2
    lc_3 = hc_3
    lc_4 = hc_4
    lc_5 = hc_5
    lc_6 = hc_6
  out.flush()
