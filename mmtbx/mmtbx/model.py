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
from libtbx.utils import Sorry



class manager(object):
  def __init__(self, processed_pdb_file,
                     rigid_body_selections = None,
                     group_b_selections    = None,
                     log = None):
    self.log = log
    self.restraints_manager = None
    self.processed_pdb_file = processed_pdb_file
    self.xray_structure = \
                     processed_pdb_file.xray_structure().deep_copy_scatterers()
    self.xray_structure_ini = self.xray_structure.deep_copy_scatterers()
    self.crystal_symmetry = \
                  processed_pdb_file.all_chain_proxies.stage_1.crystal_symmetry
    self.atom_attributes_list = \
           processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list[:]
    self.solvent_selection = self._solvent_selection()
    self.solvent_selection_ini = self._solvent_selection()
    self.locked = False
    self.rigid_body_selections = rigid_body_selections
    self.group_b_selections = group_b_selections
    if(self.rigid_body_selections is not None):
    #XXX BUG
       dim = self.xray_structure.scatterers().size()
       for sel in self.rigid_body_selections:
           assert sel.size() == dim


  def setup_restraints_manager(self,
                               plain_pairs_radius = 5.0,
                               normalization      = True,
                               show_energies      = False):
    self.restraints_manager = mmtbx.restraints.manager(
           geometry      = self.processed_pdb_file.geometry_restraints_manager(
                                      show_energies      = show_energies,
                                      plain_pairs_radius = plain_pairs_radius),
           normalization = normalization)
    if(not self.locked):
       selection = flex.bool(self.xray_structure_ini.scatterers().size(), True)
       self.restraints_manager_ini = mmtbx.restraints.manager(
            geometry      = self.restraints_manager.geometry.select(selection),
            ncs_groups    = self.restraints_manager.ncs_groups,
            normalization = self.restraints_manager.normalization)
       self.locked = True

  def deep_copy(self):
    new = manager(processed_pdb_file    = self.processed_pdb_file,
                  rigid_body_selections = None,
                  group_b_selections    = None,
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
    new.locked               = self.locked
    if(self.rigid_body_selections is not None):
       new.rigid_body_selections = []
       for item in self.rigid_body_selections:
           new.rigid_body_selections.append(item.deep_copy())
    if(self.group_b_selections is not None):
       new.group_b_selections = []
       for item in self.group_b_selections:
           new.group_b_selections.append(item.deep_copy())
    return new

  def _solvent_selection(self):
    labels = self.xray_structure.scatterers().extract_labels()
    res_name_tags = ["HOH","SOL","SOLV","WAT","DOD"]
    atom_name_tags = ["O","OH2","H","H1","H2"]
    element_tags = ["O","H"]
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
    return self.xray_structure.select(~self.solvent_selection)

  def update(self, selection):
    self.xray_structure.select_inplace(selection)
    new_atom_attributes_list = []
    rigid_body_selections = []
    group_b_selections    = []
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
    if(self.rigid_body_selections is not None):
       for s in self.rigid_body_selections:
           rigid_body_selections.append(s.select(selection))
    self.rigid_body_selections = rigid_body_selections
    if(self.group_b_selections is not None):
       for s in self.group_b_selections:
           group_b_selections.append(s.select(selection))
    self.group_b_selections = group_b_selections
    self.solvent_selection = new_solvent_selection
    self.xray_structure.scattering_type_registry()
    self.restraints_manager = mmtbx.restraints.manager(
            geometry      = self.restraints_manager.geometry.select(selection),
            ncs_groups    = self.restraints_manager.ncs_groups,
            normalization = self.restraints_manager.normalization)

  def number_of_ordered_solvent_molecules(self):
    return self.solvent_selection.count(True)

  def show_rigid_body_groups_info(self, out,
                                        text="Information about rigid groups"):
    if(self.rigid_body_selections is None): return
    if (out is None): out = sys.stdout
    print >> out
    line_len = len("| "+text+"|")
    fill_len = 80 - line_len-1
    upper_line = "|-"+text+"-"*(fill_len)+"|"
    print >> out, upper_line
    next = "| Total number of atoms = %-6d  Number of rigid groups = %-3d                |"
    natoms_total = self.xray_structure.scatterers().size()
    print >> out, next % (natoms_total, len(self.rigid_body_selections))
    print >> out, "| group: start point:                        end point:                       |"
    print >> out, "|               x      B  atom   residue <>        x      B  atom   residue   |"
    next = "| %5d: %8.3f %6.2f %5s %4s %4d <> %8.3f %6.2f %5s %4s %4d   |"
    sites = self.xray_structure.sites_cart()
    b_isos = self.xray_structure.extract_u_iso_or_u_equiv() * math.pi**2*8
    n_atoms = 0
    for i_seq, selection in enumerate(self.rigid_body_selections):
        n_atoms += selection.count(True)
        i_selection = selection.iselection()
        start = i_selection[0]
        final = i_selection[i_selection.size()-1]
        first = self.atom_attributes_list[start]
        last  = self.atom_attributes_list[final]
        print >> out, next % (i_seq+1, sites[start][0], b_isos[start],
          first.name, first.resName, first.resSeq, sites[final][0],
          b_isos[final], last.name, last.resName, last.resSeq)
    if(n_atoms != natoms_total):
       print >> out, "|                                                                             |"
       print >> out, "|                 *** Error in rigid groups definition ***                    |"
       print >> out, "|                                                                             |"
       print >> out, "| Total number of atoms in specified rigid groups does not equal to the total |"
       print >> out, "| number of atoms in the model:                                               |"
       print >> out, "| Atoms in model        = %-7d                                             |"%natoms_total
       print >> out, "| Atoms in rigid groups = %-7d                                             |"%n_atoms
       print >> out, "|"+"-"*77+"|"
       raise Sorry("Error in rigid groups definition")
    print >> out, "|"+"-"*77+"|"
    print >> out
    out.flush()


  def remove_solvent(self):
    self.update(selection = ~self.solvent_selection)

  def show_occupancy_statistics(self, out=None, text=""):
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

  def write_pdb_file(self, out, crystal_symmetry = None):
    if(crystal_symmetry is None):
       crystal_symmetry = self.crystal_symmetry
    if(crystal_symmetry is not None):
       print >> out, pdb.format_cryst1_record(
                                         crystal_symmetry = crystal_symmetry)
       print >> out, pdb.format_scale_records(
                                    unit_cell = crystal_symmetry.unit_cell())
    sites_cart  = self.xray_structure.sites_cart()
    occupancies = self.xray_structure.scatterers().extract_occupancies()
    u_isos      = self.xray_structure.extract_u_iso_or_u_equiv()
    for i_seq,atom in enumerate(self.atom_attributes_list):
        if(atom.name is None): name = "    "
        else: name = atom.name
        if(atom.altLoc is None): altLoc = " "
        else: altLoc = atom.altLoc
        if(atom.chainID is None): chainID = " "
        else: chainID = atom.chainID
        if(atom.resSeq is None): resSeq = 1
        else: resSeq = atom.resSeq
        if(atom.iCode is None): iCode = " "
        else: iCode = atom.iCode
        if(atom.segID is None): segID = "    "
        else: segID = atom.segID
        if(atom.element is None): element = "  "
        else: element = atom.element
        if(atom.charge is None): charge = "  "
        else: charge = atom.charge
        print >> out, pdb.format_atom_record(
                                    record_name = atom.record_name(),
                                    serial      = i_seq,
                                    name        = name,
                                    altLoc      = altLoc,
                                    resName     = atom.resName,
                                    chainID     = chainID,
                                    resSeq      = resSeq,
                                    iCode       = iCode,
                                    site        = sites_cart[i_seq],
                                    occupancy   = occupancies[i_seq],
                                    tempFactor  = adptbx.u_as_b(u_isos[i_seq]),
                                    segID       = segID,
                                    element     = element,
                                    charge      = charge)
    print >> out, "END"

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

  def remove_atom_with_i_seqs(self, i_seq = None, i_seqs = None):
    assert [i_seq, i_seqs].count(None) == 1
    remove_atom_selection = flex.bool(len(self.atom_attributes_list), True)
    if(i_seq is not None):
       remove_atom_selection[i_seq] = False
    if(i_seqs is not None):
       for i_seq_i in i_seqs:
           remove_atom_selection[i_seq_i] = False
    self.update(selection = remove_atom_selection)

  def replace_solvent(self, xray_structure,
                            solvent_selection,
                            atom_name    = None,
                            residue_name = None,
                            chain_id     = None):
    self.remove_solvent()
    mac = xray_structure.select(~solvent_selection)
    sol = xray_structure.select(solvent_selection)
    assert mac.scatterers().size() == self.xray_structure.scatterers().size()
    atom_atom_distances = \
              flex.sqrt(mac.difference_vectors_cart(self.xray_structure).dot())
    assert approx_equal(flex.mean_default(atom_atom_distances,0), 0)
    self.xray_structure = self.xray_structure.concatenate(sol)
    # XXX TODO NCS restraints
    # XXX RALF: throw exception if self.reduced_solvent_selection affects NCS
    geometry = self.restraints_manager.geometry
    number_of_new_solvent = sol.scatterers().size()
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
           site_symmetry_table = sol.site_symmetry_table(),
           nonbonded_types     = flex.std_string(number_of_new_solvent, "OH2"))
    self.restraints_manager = mmtbx.restraints.manager(
                         geometry      = geometry,
                         ncs_groups    = self.restraints_manager.ncs_groups,
                         normalization = self.restraints_manager.normalization)
    self.solvent_selection = solvent_selection
    if(atom_name is None):
       new_atom_name = " O  "
    else:
       new_atom_name = atom_name.strip()
       if(len(new_atom_name) < 4): new_atom_name = " " + new_atom_name
       while(len(new_atom_name) < 4):
          new_atom_name = new_atom_name+" "
    for i_seq, sc in enumerate(sol.scatterers()):
        if(residue_name is None): residue_name = "HOH"
        if(chain_id is None):     chain_id = "S"
        new_attr = pdb.atom.attributes(name        = new_atom_name,
                                       resName     = residue_name,
                                       chainID     = chain_id,
                                       element     = sc.element_symbol(),
                                       is_hetatm   = True,
                                       resSeq      = i_seq)
        self.atom_attributes_list.append(new_attr)

  def build_hydrogens(self):
    pass

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
                          text                   = "final")
    assert approx_equal(first_target_value, sso_start.target)
    assert approx_equal(minimized.final_target_value, sso_end.target)
    sso_start.show(out = self.log)
    sso_end.show(out = self.log)


  def geometry_statistics(self, other = None, show = False, text = ""):
    if(other is not None):
       stereochemistry_statistics_obj = stereochemistry_statistics(
                             xray_structure         = self.xray_structure,
                             xray_structure_ref     = other.xray_structure,
                             restraints_manager     = self.restraints_manager,
                             restraints_manager_ref = other.restraints_manager,
                             text                   = text)
    else:
       stereochemistry_statistics_obj = stereochemistry_statistics(
                          xray_structure         = self.xray_structure,
                          xray_structure_ref     = self.xray_structure_ini,
                          restraints_manager     = self.restraints_manager,
                          restraints_manager_ref = self.restraints_manager_ini,
                          text                   = text)
    if(show): stereochemistry_statistics_obj.show(out = self.log)
    return stereochemistry_statistics_obj

  def adp_statistics(self, iso_restraints,
                           other    = None,
                           wilson_b = None,
                           tan_b_iso_max= None,
                           show     = False,
                           text     = ""):
    if(other is not None):
       adp_statistics_obj = adp_statistics(
                             xray_structure         = self.xray_structure,
                             xray_structure_ref     = other.xray_structure,
                             restraints_manager     = self.restraints_manager,
                             restraints_manager_ref = other.restraints_manager,
                             iso_restraints         = iso_restraints,
                             wilson_b               = wilson_b,
                             tan_b_iso_max              = tan_b_iso_max,
                             text                   = text)
    else:
       adp_statistics_obj = adp_statistics(
                          xray_structure         = self.xray_structure,
                          xray_structure_ref     = self.xray_structure_ini,
                          restraints_manager     = self.restraints_manager,
                          restraints_manager_ref = self.restraints_manager_ini,
                          iso_restraints         = iso_restraints,
                          tan_b_iso_max              = tan_b_iso_max,
                          wilson_b               = wilson_b,
                          text                   = "")
    if(show): adp_statistics_obj.show(out = self.log)
    return adp_statistics_obj



class adp_statistics(object):
  def __init__(self,
               xray_structure,
               xray_structure_ref,
               restraints_manager,
               restraints_manager_ref,
               iso_restraints,
               wilson_b = None,
               tan_b_iso_max = None,
               text=""):
    adopt_init_args(self, locals())
    self.text = self.text + "ADP statistics"
    energies_adp_iso = restraints_manager.energies_adp_iso(
                                            xray_structure    = xray_structure,
                                            parameters        = iso_restraints,
                                            wilson_b          = wilson_b,
                                            tan_b_iso_max     = tan_b_iso_max,
                                            compute_gradients = True)
    energies_adp_iso_ref = restraints_manager_ref.energies_adp_iso(
                                        xray_structure    = xray_structure_ref,
                                        parameters        = iso_restraints,
                                        wilson_b          = wilson_b,
                                        tan_b_iso_max     = tan_b_iso_max,
                                        compute_gradients = True)
    eps = math.pi**2*8
    self.b_isos = xray_structure.extract_u_iso_or_u_equiv() * eps
    self.b_isos_ref = xray_structure_ref.extract_u_iso_or_u_equiv() * eps
    # XXX TODO NCS restraints
    # XXX RALF/PAVEL show NCS statistics
    self.b_iso_min  = flex.min_default(self.b_isos, 0)
    self.b_iso_max  = flex.max_default(self.b_isos, 0)
    self.b_iso_mean = flex.mean_default(self.b_isos, 0)
    self.b_iso_min_ref  = flex.min_default(self.b_isos_ref, 0)
    self.b_iso_max_ref  = flex.max_default(self.b_isos_ref, 0)
    self.b_iso_mean_ref = flex.mean_default(self.b_isos_ref, 0)
    self.target_adp_iso = energies_adp_iso.target
    self.grad_adp_iso = energies_adp_iso.gradients
    self.target_adp_iso_ref = energies_adp_iso_ref.target
    self.grad_adp_iso_ref = energies_adp_iso_ref.gradients
    self.norm_of_grad_adp_iso = self.grad_adp_iso.norm()
    self.norm_of_grad_adp_iso_ref = self.grad_adp_iso_ref.norm()
    self.n_zero = (self.b_isos < 0.5).count(True)
    self.n_zero_ref = (self.b_isos_ref < 0.5).count(True)

  def show(self, out=None):
    if (out is None): out = sys.stdout
    line_len = len("| "+self.text+"|")
    fill_len = 80 - line_len-1
    ends1 = " "*36 + "|"
    ends2 = " "*32 + "|"
    p = " "
    v = "|"
    print >> out
    print >> out, "|-"+self.text+"-"*(fill_len)+"|"
    if(self.wilson_b is not None):
       print >> out, "| Wilson B =%8.3f | target = %12.5f | norm of gradient "\
                  "= %11.5f |"%(self.wilson_b,self.target_adp_iso,
                                                     self.norm_of_grad_adp_iso)
    else:
       print >> out, "| Wilson B = %7s | target = %12.5f | norm of gradient "\
                  "= %11.5f |"%(str(None),self.target_adp_iso,
                                                     self.norm_of_grad_adp_iso)
    print >> out, "|"+"-"*77+"|"
    print >> out, "|"+12*" "+"Reference model"+" "*10+"|"+" "*12+\
                  "Current model"+" "*14+"|"
    print >> out, "| "+"- "*38+"|"
    print >> out, "|                           Isotropic B-factors:        "\
                  "                      |"
    print >> out, "|     min      max      mean          |     min      max "\
                  "     mean            |"
    print >> out, "| %7.3f  %7.3f   %7.3f   "%\
                  (self.b_iso_min_ref,self.b_iso_max_ref,self.b_iso_mean_ref),\
                  "      | %7.3f  %7.3f   %7.3f            |"%\
                  (self.b_iso_min,self.b_iso_max,self.b_iso_mean)
    print >> out, "| number of B < 0.5:%6d"%self.n_zero_ref,"           | ",\
                  "number of B < 0.5:%6d             |"%self.n_zero
    print >> out, "| "+"- "*38+"|"
    print >> out, "|                     Distribution of isotropic B-factors:"\
                  "                    |"
    histogram_1 = flex.histogram(data = self.b_isos_ref, n_slots = 10)
    low_cutoff_1 = histogram_1.data_min()
    histogram_2 = flex.histogram(data = self.b_isos, n_slots = 10)
    low_cutoff_2 = histogram_2.data_min()
    for (i_1,n_1),(i_2,n_2) in zip(enumerate(histogram_1.slots()),
                                   enumerate(histogram_2.slots())):
      high_cutoff_1 = histogram_1.data_min() + histogram_1.slot_width()*(i_1+1)
      high_cutoff_2 = histogram_2.data_min() + histogram_2.slot_width()*(i_2+1)
      print >> out, "|  %9.3f -%9.3f:%8d      |    %9.3f -%9.3f:%8d      |" % \
             (low_cutoff_1,high_cutoff_1,n_1,low_cutoff_2,high_cutoff_2,n_2)
      low_cutoff_1 = high_cutoff_1
      low_cutoff_2 = high_cutoff_2
    print >> out, "|"+"-"*77+"|"
    print >> out
    out.flush()

class stereochemistry_statistics(object):
  def __init__(self,
               xray_structure,
               xray_structure_ref,
               restraints_manager,
               restraints_manager_ref,
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
    bond_deltas_1 = geometry_restraints.bond_deltas(
                     sites_cart         = self.xray_structure_ref.sites_cart(),
                     sorted_asu_proxies = rmg_1.pair_proxies().bond_proxies)
    bond_deltas_2 = geometry_restraints.bond_deltas(
                        sites_cart         = self.xray_structure.sites_cart(),
                        sorted_asu_proxies = rmg_2.pair_proxies().bond_proxies)
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
    print >> out, "|                        |                          |                         |"
    print >> out, "|        deviations from ideal values for           |                         |"
    print >> out, "|                        |                          |                         |"
    print >> out, "| bonds                  | angles                   | nonbonded contacts      |"
    print >> out, "|                        |                          |                         |"
    show_6_histograms(h_1, h_2, h_3, h_4, h_5, h_6, n_slots = n_slots, out=out)
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def show(self, out=None):
    if (out is None): out = sys.stdout
    line_len = len("| "+self.text+"|")
    fill_len = 80 - line_len-1
    print >> out, "| "+self.text+"-"*(fill_len)+"|"
    print >> out, "|Type| Deviation from ideal |   Targets  ||Target (sum)|| Deviation of start  |"
    print >> out, "|    |  mean     max    min |            ||            || model from current  |"
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
        print >> out, \
          "|   NCS group %2d: min = %9.6f max = %9.6f mean = %9.6f %s|" % (
            i_group+1, flex.min(rms), flex.max(rms), flex.mean(rms), " "*11)
      print >> out, "|"+"-"*77+"|"
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
  format = "|%5.3f-%5.3f: %5d/%5d|%6.2f-%6.2f: %5d/%5d|%5.2f-%5.2f: %5d/%5d |"
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
