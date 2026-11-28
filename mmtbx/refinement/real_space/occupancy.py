from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")
import mmtbx.refinement.real_space.adp as dependency
import mmtbx.refinement.occupancies
import mmtbx.refinement.refinement_flags
from libtbx.utils import Sorry
from scitbx.array_family import flex

max_rounds = 10

def raise_ncs_and_custom_occ_selections_not_supported(occ_params, model):
  custom_occ_selections = len(occ_params.constrained_group) + \
                          len(occ_params.individual) + \
                          len(occ_params.remove_selection)
  # Safeguard against using NCS constraints and custom occupancy selections
  if(custom_occ_selections>0 and model.ncs_constraints_present()):
    msg="NCS constraints and custom occupancy selections are not supported."
    raise Sorry(msg)

class ncs_aware_refinement(object):
  def __init__(self, map_model_manager, d_min, atom_radius, occupancy_params,
               nproc=1, log = None):
    self.mmm              = map_model_manager
    self.nproc            = nproc
    self.d_min            = d_min
    self.atom_radius      = atom_radius
    self.log              = log
    self.occupancy_params = occupancy_params
    # Determine if need to do occ refinement
    f1 = self.occupancy_params.constrained_group
    if(len(f1)==1):
     if(len(f1[0].selection)==1):
       if(f1[0].selection[0] is None):
         f1 = []
    f2 = self.occupancy_params.individual
    if(len(f2)==1):
      if(len(f2)==1):
        if(f2[0] is None):
          f2 = []
    f3 = self.occupancy_params.remove_selection
    if(len(f3)==1):
      if(len(f3)==1):
        if(f3[0] is None):
          f3 = []
    custom_occ_selections = len(f1) + len(f2) + len(f3)
    # Selections with empty ([None]) entries removed; used from here on
    self._constrained_group, self._individual, self._remove_selection = f1,f2,f3
    # Safeguard against using NCS constraints and custom occupancy selections
    raise_ncs_and_custom_occ_selections_not_supported(
      occ_params = self.occupancy_params, model = self.mmm.model())
    self.selections = mmtbx.refinement.occupancies.occupancy_selections(
      model                              = self.mmm.model(),
      add_water                          = False,
      other_constrained_groups           = f1,
      other_individual_selection_strings = f2,
      remove_selection                   = f3,
      as_flex_arrays                     = True,
      constrain_correlated_3d_groups     = False,
      log                                = self.log)
    proceed = True
    if(self.selections is None or len(self.selections)==0):
      proceed = False
    if(proceed):
      #
      ncs_groups = self.mmm.model().get_ncs_groups()
      if(ncs_groups is None or len(ncs_groups)==0):
        values = self.run_one()
        self.mmm.model().set_occupancies(values = values)
      else:
        values = self.mmm.model().get_occ()
        for i, g in enumerate(ncs_groups):
          values_g = self.run_one(selection = g.master_iselection)
          values = values.set_selected(g.master_iselection, values_g)
          for j, c in enumerate(g.copies):
            values = values.set_selected(c.iselection, values_g)
        self.mmm.model().set_occupancies(values = values)

  def run_one(self, selection=None):
    """
    Refine occupancies of the whole model (or of the selected part, an NCS
    master) against the map in one calculation. nproc is not used (see ADP
    refinement): there is one code path and the result does not depend on it.
    """
    model = self.mmm.model()
    if(selection is not None): model = model.select(selection)
    return self.run_one_one(args = [model,])

  def run_one_one(self, args):
    model = args[0].deep_copy()
    # selections for refinable occupancies
    selections = mmtbx.refinement.occupancies.occupancy_selections(
      model                              = model,
      add_water                          = False,
      other_constrained_groups           = self._constrained_group,
      other_individual_selection_strings = self._individual,
      remove_selection                   = self._remove_selection,
      as_flex_arrays                     = True,
      constrain_correlated_3d_groups     = False,
      log                                = self.log)
    if(selections is None or len(selections)==0):
      return model.get_occ()
    #
    fmodel = dependency.map_and_model_to_fmodel(
      map_data       = self.mmm.map_data().deep_copy(),
      xray_structure = model.get_xray_structure(),
      atom_radius    = self.atom_radius,
      reset_adp      = False,
      d_min          = self.d_min)
    model.set_xray_structure(xray_structure = fmodel.xray_structure)
    #
    refla = mmtbx.refinement.refinement_flags.manager(
      occupancies=True, s_occupancies = selections)
    model.set_refinement_flags(flags = refla)
    if(self.log is not None):
      print("r_start: %6.4f"%fmodel.r_work(), file=self.log)
    # Rounds of refinement; the data are renewed before each (map outside
    # the mask around atoms comes from the present model, see
    # adp.fill_outside_mask), until occupancies stop changing: at least
    # three rounds (as before), max_rounds at most.
    for cycle in range(max_rounds):
      occ_before = fmodel.xray_structure.scatterers().extract_occupancies()
      if(cycle > 2):
        if(flex.max(flex.abs(occ_before - occ_previous)) < 0.01): break
      occ_previous = occ_before
      dependency.update_f_obs(fmodel = fmodel)
      o = mmtbx.refinement.occupancies.manager(
        fmodels                  = mmtbx.fmodels(fmodel_xray = fmodel),
        model                    = model,
        max_number_of_iterations = 50,
        number_of_macro_cycles   = 1,
        occupancy_max            = 1,
        occupancy_min            = 0,
        exclude_hd               = False,
        log                      = self.log)
    if(self.log is not None):
      print("r_final: %6.4f"%fmodel.r_work(), file=self.log)
    #
    return fmodel.xray_structure.scatterers().extract_occupancies()
