from __future__ import absolute_import, division, print_function
import six
import mmtbx.model
import iotbx.pdb
import boost_adaptbx.boost.python as bp
from libtbx.utils import null_out
from cctbx.array_family import flex

aa_codes = [ # XXX PVA: Temporary, ad hoc, remove later!
"ALA",
"ARG",
"ASN",
"ASP",
"CYS",
"GLN",
"GLU",
"GLY",
"HIS",
"ILE",
"LEU",
"LYS",
"MET",
"MSE",
"PHE",
"PRO",
"SER",
"THR",
"TRP",
"TYR",
"VAL"
]

ext = bp.import_ext("cctbx_geometry_restraints_ext")

def mon_lib_query(residue, mon_lib_srv):
    get_func = getattr(mon_lib_srv, "get_comp_comp_id", None)
    if (get_func is not None): return get_func(comp_id=residue)
    return mon_lib_srv.get_comp_comp_id_direct(comp_id=residue)

def exclude_h_on_SS(model):
  rm = model.get_restraints_manager()
  bond_proxies_simple, asu = rm.geometry.get_all_bond_proxies(
    sites_cart = model.get_sites_cart())
  els = model.get_hierarchy().atoms().extract_element()
  ss_i_seqs = []
  all_proxies = [p for p in bond_proxies_simple]
  for proxy in asu:
    all_proxies.append(proxy)
  for proxy in all_proxies:
    if(  isinstance(proxy, ext.bond_simple_proxy)): i,j=proxy.i_seqs
    elif(isinstance(proxy, ext.bond_asu_proxy)):    i,j=proxy.i_seq,proxy.j_seq
    else: assert 0 # never goes here
    if([els[i],els[j]].count("S")==2): # XXX may be coordinated if metal edits used
      ss_i_seqs.extend([i,j])
  sel_remove = flex.size_t()
  for proxy in all_proxies:
    if(  isinstance(proxy, ext.bond_simple_proxy)): i,j=proxy.i_seqs
    elif(isinstance(proxy, ext.bond_asu_proxy)):    i,j=proxy.i_seq,proxy.j_seq
    else: assert 0 # never goes here
    if(els[i] in ["H","D"] and j in ss_i_seqs): sel_remove.append(i)
    if(els[j] in ["H","D"] and i in ss_i_seqs): sel_remove.append(j)
  return model.select(~flex.bool(model.size(), sel_remove))


def exclude_h_on_coordinated_S(model): # XXX if edits used it should be like in exclude_h_on_SS
  rm = model.get_restraints_manager().geometry
  els = model.get_hierarchy().atoms().extract_element()
  # Find possibly coorinated S
  exclusion_list = ["H","D","T","S","O","P","N","C","SE"]
  sel_s = []
  for proxy in rm.pair_proxies().nonbonded_proxies.simple:
    i,j = proxy.i_seqs
    if(els[i] == "S" and not els[j] in exclusion_list): sel_s.append(i)
    if(els[j] == "S" and not els[i] in exclusion_list): sel_s.append(j)
  # Find H attached to possibly coordinated S
  bond_proxies_simple, asu = rm.get_all_bond_proxies(
    sites_cart = model.get_sites_cart())
  sel_remove = flex.size_t()
  for proxy in bond_proxies_simple:
    i,j = proxy.i_seqs
    if(els[i] in ["H","D"] and j in sel_s): sel_remove.append(i)
    if(els[j] in ["H","D"] and i in sel_s): sel_remove.append(j)
  return model.select(~flex.bool(model.size(), sel_remove))

def add(model,
        use_neutron_distances=False,
        adp_scale=1,
        exclude_water=True,
        protein_only=False,
        stop_for_unknowns=True,
        remove_first=True):
  if(remove_first):
    model = model.select(~model.get_hd_selection())
  pdb_hierarchy = model.get_hierarchy()
  mon_lib_srv = model.get_mon_lib_srv()
  get_class = iotbx.pdb.common_residue_names_get_class
  """
  for pmodel in pdb_hierarchy.models():
    for chain in pmodel.chains():
      for residue_group in chain.residue_groups():
        for conformer in residue_group.conformers():
          for residue in conformer.residues():
            print list(residue.atoms().extract_name())
  """
  #XXX This breaks for 1jxt, residue 2, TYR
  for chain in pdb_hierarchy.only_model().chains():
    for rg in chain.residue_groups():
      for ag in rg.atom_groups():
        #print list(ag.atoms().extract_name())
        if(get_class(name=ag.resname) == "common_water"): continue
        if(protein_only and
           not ag.resname.strip().upper() in aa_codes): continue
        actual = [a.name.strip().upper() for a in ag.atoms()]
        mlq = mon_lib_query(residue=ag.resname, mon_lib_srv=mon_lib_srv)
        expected_h   = []
        for k, v in six.iteritems(mlq.atom_dict()):
          if(v.type_symbol=="H"): expected_h.append(k)
        missing_h = list(set(expected_h).difference(set(actual)))
        if 0: print(ag.resname, missing_h)
        new_xyz = ag.atoms().extract_xyz().mean()
        hetero = ag.atoms()[0].hetero
        for mh in missing_h:
          # TODO: this should be probably in a central place
          if len(mh) < 4: mh = (' ' + mh).ljust(4)
          a = (iotbx.pdb.hierarchy.atom()
            .set_name(new_name=mh)
            .set_element(new_element="H")
            .set_xyz(new_xyz=new_xyz)
            .set_hetero(new_hetero=hetero))
          ag.append_atom(a)
  pdb_hierarchy.atoms().reset_serial()
  #pdb_hierarchy.sort_atoms_in_place()
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  p.pdb_interpretation.use_neutron_distances = use_neutron_distances
  p.pdb_interpretation.proceed_with_excessive_length_bonds=True
  #p.pdb_interpretation.restraints_library.cdl=False # XXX this triggers a bug !=360
  ro = model.get_restraint_objects()
  model = mmtbx.model.manager(
    model_input               = None,
    pdb_hierarchy             = pdb_hierarchy,
    build_grm                 = True,
    stop_for_unknowns         = stop_for_unknowns,
    crystal_symmetry          = model.crystal_symmetry(),
    restraint_objects         = ro,
    pdb_interpretation_params = p,
    log                       = null_out())
#  # Remove lone H
#  sel_h = model.get_hd_selection()
#  sel_isolated = model.isolated_atoms_selection()
#  sel_lone = sel_h & sel_isolated
#  model = model.select(~sel_lone)
  # Only keep H which have been parameterized in riding H procedure
  sel_h = model.get_hd_selection()
  model.setup_riding_h_manager()
  sel_h_in_para = flex.bool(
                [bool(x) for x in model.riding_h_manager.h_parameterization])
  sel_h_not_in_para = sel_h_in_para.exclusive_or(sel_h)
  model = model.select(~sel_h_not_in_para)

  model = exclude_h_on_SS(model = model)
  model = exclude_h_on_coordinated_S(model = model)
  # Reset occupancies, ADPs and idealize
  model.reset_adp_for_hydrogens(scale = adp_scale)
  model.reset_occupancy_for_hydrogens_simple()
  model.idealize_h_riding()
  return model
