from __future__ import absolute_import, division, print_function

import time
from six.moves import cStringIO as StringIO
import iotbx.pdb
#import mmtbx.model

import libtbx.load_env
import cctbx.geometry_restraints.process_nonbonded_proxies as pnp
from cctbx import adptbx
from iotbx import phil
from cctbx.array_family import flex
from libtbx import group_args
from libtbx.str_utils import make_sub_header
from libtbx.utils import null_out
#from libtbx import easy_run
#from cctbx import miller
#from mmtbx import monomer_library
from mmtbx import real_space_correlation
#from mmtbx import map_tools
from elbow.command_line.ready_set import model_interface as ready_set_model_interface
from six.moves import zip


master_params_str = """
validate_ligands {

place_hydrogens = True
  .type = bool
  .help = Add H atoms with ready_set.

nproc = 1
  .type = int

}
"""


def master_params():
  return phil.parse(master_params_str, process_includes = False)


# =============================================================================


class manager(dict):
  '''
  Manager class for ALL ligands
  '''
  def __init__(self,
               model,
               fmodel,
               params,
               log=None):
    self.model = model
    self.params = params
    self.log   = log
    self.fmodel = fmodel
    #
    self.nproc = params.nproc

    if fmodel is not None:
      self.two_fofc_map = fmodel.map_coefficients(map_type="2mFo-DFc",
        fill_missing=False,
        merge_anomalous=True).fft_map(
          resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
      self.fofc_map = fmodel.map_coefficients(map_type="mFo-DFc",
        fill_missing=False,
        merge_anomalous=True).fft_map(
          resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
      self.fmodel_map = fmodel.map_coefficients(map_type="Fmodel",
        fill_missing=False,
        merge_anomalous=True).fft_map(
          resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
    else:
      self.two_fofc_map = None
      self.fofc_map = None
      self.fmodel_map = None


  def parallel_populate(self, args):
    from libtbx import easy_mp
    def _run(lr, func):
      func = getattr(lr, func)
      rc = func()
      return rc
    funcs = []
    ligand_results = []
    inputs = []
    for id_tuple, rg, altloc, conformer_isel in args:
      id_list = list(id_tuple)
      id_list.append(altloc)
      lr = ligand_result(
        model = self.model,
        fmodel = self.fmodel,
        readyset_model = self.readyset_model,
        two_fofc_map = self.two_fofc_map,
        fofc_map = self.fofc_map,
        fmodel_map = self.fmodel_map,
        isel = conformer_isel,
        id_list = id_list)
      ligand_results.append(lr)
      for attr, func in lr._result_attrs.items():
        funcs.append([lr, func])
        inputs.append([id_tuple, altloc, func])

    results = []
    t0=time.time()
    for i, (args, res, err_str) in enumerate(easy_mp.multi_core_run(
      _run,
      funcs,
      self.nproc,
      )):
      results.append([args, res, err_str])
      if self.nproc>1:
        print('\n  Returning Selection : %s AltLoc : %s Func : %s' % tuple(inputs[i]))
        print('  Cumulative time: %6.2f (s)' % (time.time()-t0))
      if err_str:
        print('Error output from %s' % args)
        print(err_str)
        print('_'*80)

    i=0
    for lr in ligand_results:
      for attr, func in lr._result_attrs.items():
        setattr(lr, attr, results[i][1])
        i+=1
    return ligand_results


  def run(self):
    self.get_readyset_model_with_grm()
    args = []
    def _generate_ligand_isel():
      #done = []
      ph = self.model.get_hierarchy()
      ligand_isel_dict = self.get_ligands(ph = ph)
      for id_tuple, isel in ligand_isel_dict.items():
        for rg in ph.select(isel).residue_groups():
          ligand_dict = {}
          for conformer in rg.conformers():
            altloc = conformer.altloc
            #key = (id_tuple, altloc)
            #print (key)
            #if key in done: continue
            #done.append(key)
            conformer_isel = conformer.atoms().extract_i_seq()
            yield id_tuple, rg, altloc, conformer_isel
    for id_tuple, rg, altloc, conformer_isel in _generate_ligand_isel():
      args.append([id_tuple, rg, altloc, conformer_isel])
    results = self.parallel_populate(args)
    for lr, (id_tuple, rg, altloc, conformer_isel) in zip(results,
                                                          _generate_ligand_isel()):
      ligand_dict = self.setdefault(id_tuple, {})
      ligand_dict[altloc] = lr


  def get_readyset_model_with_grm(self):
    '''
    Run ready set by default (to get ligand cif and/or add H atoms)
    TODO: Once it it refactored, make running it optional
          Complicated case is when cif is input for one ligand, but missing
          for another ligand
    '''
    self.readyset_model = None
    #if not self.model.has_hd():
    if (self.params.place_hydrogens):
      params=["add_h_to_water=False",
               "optimise_final_geometry_of_hydrogens=False",
               "--silent"]
    else:
      params=["add_h_to_water=False",
              "hydrogens=False",
              "optimise_final_geometry_of_hydrogens=False",
              "--silent"]
    assert (libtbx.env.has_module(name="reduce"))
    assert (libtbx.env.has_module(name="elbow"))
    # runs reduce internally
    self.readyset_model = ready_set_model_interface(
        model  = self.model,
        params = params)
    self.readyset_model.set_log(null_out())




  def get_ligands(self, ph):
    '''
    Store ligands as list of iselections --> better way? Careful if H will be
    added at some point!
    '''
    ligand_isel_dict = {}
    get_class = iotbx.pdb.common_residue_names_get_class
    exclude = ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
               "modified_rna_dna", "ccp4_mon_lib_rna_dna", "common_water",
                "common_element"]
    for model in ph.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          for resname in rg.unique_resnames():
            if (not get_class(name=resname) in exclude):
              iselection = rg.atoms().extract_i_seq()
              id_tuple = (model.id, chain.id, rg.resseq)
              ligand_isel_dict[id_tuple] = iselection
    return ligand_isel_dict


  def show_ligand_counts(self):
    make_sub_header(' Ligands in input model ', out=self.log)
    for id_tuple, ligand_dict in self.items():
      for altloc, lr in ligand_dict.items():
        print(lr.id_str)


  def show_adps(self):
    '''
    Show results for ADPs of ligand and surrounding atoms
    '''
    make_sub_header(' ADPs ', out=self.log)
    pad1 = ' '*18
    print(pad1, "min   max    mean   n_iso   n_aniso", file=self.log)
    for id_tuple, ligand_dict in self.items():
      for altloc, lr in ligand_dict.items():
        adps = lr.get_adps()
        print(lr.id_str.ljust(14), '%7s%7s%7s%7s%7s' %
          (round(adps.b_min,1), round(adps.b_max,1), round(adps.b_mean,1),
           adps.n_iso, adps.n_aniso), file = self.log)
        if (adps.b_mean_within is not None):
          print('neighbors'.ljust(14), '%7s%7s%7s' %
            (round(adps.b_min_within,1), round(adps.b_max_within,1),
             round(adps.b_mean_within,1) ), file = self.log)


  def show_ligand_occupancies(self):
    '''
    Show results for ligand occupancies
    '''
    make_sub_header(' Occupancies ', out=self.log)
    pad1 = ' '*20
    print('If three values: min, max, mean, otherwise the same occupancy for entire ligand.', \
      file=self.log)
    for id_tuple, ligand_dict in self.items():
      for altloc, lr in ligand_dict.items():
        occs = lr.get_occupancies()
        if (occs.occ_min == occs.occ_max):
          print(lr.id_str.ljust(16), occs.occ_min, file = self.log)
        else:
          print(lr.id_str.ljust(16), '%s   %s   %s' %
            (occs.occ_min, occs.occ_max, occs.occ_mean), file = self.log)


  def show_ccs(self):
    '''
    Show results for correlation coefficients
    '''
    if self.fmodel is None: return
    make_sub_header(' Correlation coefficients ', out=self.log)
    for id_tuple, ligand_dict in self.items():
      for altloc, lr in ligand_dict.items():
        ccs = lr.get_ccs()
        cc_two_fofc = round(ccs.cc_two_fofc, 2)
        cc_fofc = round(ccs.cc_fofc, 2)
        fofc_min  = round(ccs.fofc_min, 2)
        fofc_max  = round(ccs.fofc_max, 2)
        fofc_mean = round(ccs.fofc_mean, 2)
        print(lr.id_str.ljust(16),
          cc_two_fofc, cc_fofc, fofc_min, fofc_max, fofc_mean, file = self.log)


  def show_nonbonded_overlaps(self):
    '''
    Print results for overlaps
    '''
    for id_tuple, ligand_dict in self.items():
      for altloc, lr in ligand_dict.items():
        clashes_result = lr.get_overlaps()
        print(clashes_result.clashes_str, file=self.log)


# =============================================================================


class ligand_result(object):
  '''
  Class storing info per ligand
  '''
  def __init__(self,
               model,
               fmodel,
               readyset_model,
               two_fofc_map,
               fofc_map,
               fmodel_map,
               isel,
               id_list):
    self.model = model
    self.fmodel = fmodel
    self.isel = isel
    self.readyset_model = readyset_model
    self.two_fofc_map = two_fofc_map
    self.fofc_map = fofc_map
    self.fmodel_map = fmodel_map

    # results
    self._result_attrs = {'_occupancies' : 'get_occupancies',
                          '_adps'        : 'get_adps',
                          '_overlaps' : 'get_overlaps',
                          '_ccs'         : 'get_ccs'
    }
    for attr, func in self._result_attrs.items():
      setattr(self, attr, None)
      assert hasattr(self, func)
    # to be used internally
    self._ph = self.model.get_hierarchy()
    self._atoms = self._ph.select(self.isel).atoms()
    self._xrs = self.model.get_xray_structure()
    self._xrs_ligand = self.model.select(isel).get_xray_structure()
    self._get_id_str(id_list=id_list)


  def __repr__(self):
    outl = 'ligand %s\n' % self.id_str
    for attr in self._result_attrs:
      outl += '  %s : %s\n' % (attr, getattr(self, attr))
    return outl


  def _get_id_str(self, id_list):
    rg_ligand = self._ph.select(self.isel).only_residue_group()
    resname = ",".join(rg_ligand.unique_resnames())
    self.id_str = " ".join([id_list[0], id_list[1], resname, id_list[2], id_list[3]])
    self.sel_str = " ".join(['chain', id_list[1], 'and resseq', id_list[2]])
    if (id_list[0] != ''):
      self.sel_str = " ".join(['model', id_list[0], 'and', self.sel_str])
    if (id_list[3] != ''):
      self.sel_str = " ".join([self.sel_str, 'and altloc', id_list[3]])
    self.id_str = self.id_str.strip()


  def get_ccs(self):
    # still a stub
    if self.fmodel is None: return
    manager = real_space_correlation.selection_map_statistics_manager(
      atom_selection    = self.isel,
      xray_structure    = self._xrs,
      fft_m_real        = self.two_fofc_map.all(),
      fft_n_real        = self.two_fofc_map.focus(),
      exclude_hydrogens = True)
    stats_two_fofc = manager.analyze_map(
      map       = self.two_fofc_map,
      model_map = self.fmodel_map,
      min       = 1.5)
    stats_fofc = manager.analyze_map(
      map       = self.fofc_map,
      model_map = self.fmodel_map,
      min       = -3.0)

#    params = real_space_correlation.master_params().extract()
#
#    results = real_space_correlation.simple(
#      fmodel        = self.fmodel,
#      pdb_hierarchy = self._ph.select(self.isel),
#      params        = None,
#      show_results  = True,
#      log           = None)
#    params.map_2.type = 'mFobs-DFc'
#    sel_bool = self._ph.atom_selection_cache().selection(
#      string = self.sel_str)
#    self.fmodel.update_xray_structure(
#      xray_structure      = self._xrs.select(~sel_bool),
#      update_f_calc       = True,
#      update_f_mask       = True,
#      force_update_f_mask = True)
#    fmodel = mmtbx.f_model.manager(
#     f_obs          = self.fmodel.f_obs(),
#     r_free_flags   = self.fmodel.r_free_flags(),
#     xray_structure = self._xrs.select(~sel_bool))
#
#    fmodel.update_all_scales()
#    results_fofc = real_space_correlation.simple(
#      fmodel        = fmodel,
#      pdb_hierarchy = self._ph.select(self.isel),
#      params        = params,
#      show_results  = True,
#      log           = None)

#    mc_diff = map_tools.electron_density_map(
#      fmodel = self.fmodel).map_coefficients(
#        map_type         = "mFo-DFc",
#        isotropize       = True,
#        fill_missing     = False)
#    crystal_gridding =
#    fft_map = miller.fft_map(
#      crystal_gridding     = crystal_gridding,
#      fourier_coefficients = mc_diff)
#    fft_map.apply_sigma_scaling()
#    map_data = fft_map.real_map_unpadded()
#    box = mmtbx.utils.extract_box_around_model_and_map(
#      xray_structure = self._xrs.select(self.isel),
#      map_data       = map_data,
#      box_cushion    = 2.1)

    self._ccs = group_args(
      cc_two_fofc = stats_two_fofc.cc,
      cc_fofc = stats_fofc.cc,
      two_fofc_min = stats_two_fofc.min,
      two_fofc_max = stats_two_fofc.max,
      two_fofc_mean = stats_two_fofc.mean,
      fofc_min = stats_fofc.min,
      fofc_max = stats_fofc.max,
      fofc_mean = stats_fofc.mean,
      n_below_two_fofc_cutoff = stats_two_fofc.n_below_min,
      n_below_fofc_cutoff = stats_fofc.n_below_min
      )
    return self._ccs


  def get_adps(self):
    if self._adps is None:
      b_isos  = self._xrs_ligand.extract_u_iso_or_u_equiv() * adptbx.u_as_b(1.)
      n_iso   = self._xrs_ligand.use_u_iso().count(True)
      n_aniso = self._xrs_ligand.use_u_aniso().count(True)
      n_zero  = (b_isos < 0.01).count(True)
      # TODO: what number as cutoff?
      n_above_100 = (b_isos > 100).count(True)
      isel_above_100 = (b_isos > 100).iselection()
      b_min, b_max, b_mean = b_isos.min_max_mean().as_tuple()

      #if this selection is used somewhere else, it might be better to do it outside
      within_radius = 3.0 #TODO should this be a parameter?
      sel_within_str = '(within (%s, %s)) and (protein or water)' % (within_radius, self.sel_str)
      isel_within = self.model.iselection(sel_within_str)
      xrs_within = self._xrs.select(isel_within)
      b_isos_within = xrs_within.extract_u_iso_or_u_equiv() * adptbx.u_as_b(1.)
      b_min_within, b_max_within, b_mean_within = b_isos_within.min_max_mean().as_tuple()

      self._adps = group_args(
        n_iso          = n_iso,
        n_aniso        = n_aniso,
        n_zero         = n_zero,
        n_above_100    = n_above_100,
        isel_above_100 = isel_above_100,
        b_min          = b_min,
        b_max          = b_max,
        b_mean         = b_mean,
        b_min_within   = b_min_within,
        b_max_within   = b_max_within,
        b_mean_within  = b_mean_within
        )
    return self._adps


  def get_occupancies(self):
    if self._occupancies is None:
      eps = 1.e-6
      occ = self._atoms.extract_occ()
      mmm = occ.min_max_mean()
      occ_min = mmm.min
      occ_max = mmm.max
      occ_mean = mmm.mean
      negative_count = (occ<0).count(True)
      negative_isel = (occ<0).iselection()
      zero_count = (flex.abs(occ)<eps).count(True)
      zero_isel = (flex.abs(occ)<eps).iselection()
      less_than_dot9_isel = (occ<0.9).iselection()

      self._occupancies = group_args(
      occ_min             = occ_min,
      occ_max             = occ_max,
      occ_mean            = occ_mean,
      negative_count      = negative_count,
      negative_isel       = negative_isel,
      zero_count          = zero_count,
      zero_isel           = zero_isel,
      less_than_dot9_isel = less_than_dot9_isel
      )

    return self._occupancies


  def get_overlaps(self):
    '''
    Obtain overlaps involving ligands
    '''
    # A model with H atoms is necessary to process overlaps
    if not self.readyset_model.has_hd():
      return None
    elif self._overlaps is None:
      # sel_within could be done in __init__?
      within_radius = 3.0
      # TODO clashes with other ligands?
      sel_within_str = '%s or (residues_within (%s, %s)) and (protein or water)' \
        % (self.sel_str, within_radius, self.sel_str)
      sel_within = self.readyset_model.selection(sel_within_str)
      isel_ligand_within = self.readyset_model.select(sel_within).iselection(self.sel_str)
      #sel = flex.bool([True]*len(sel_within))
      model_within = self.readyset_model.select(sel_within)

      processed_nbps = pnp.manager(model = model_within)
      clashes = processed_nbps.get_clashes()
      clashes_dict   = clashes._clashes_dict

      ligand_clashes_dict = dict()
      for iseq_tuple, record in clashes_dict.items():
        if (iseq_tuple[0] in isel_ligand_within or
            iseq_tuple[1] in isel_ligand_within):
          ligand_clashes_dict[iseq_tuple] = record

      ligand_clashes = pnp.clashes(
                      clashes_dict = ligand_clashes_dict,
                      model        = model_within)

      string_io = StringIO()
      ligand_clashes.show(log=string_io, show_clashscore=False)
      results = ligand_clashes.get_results()

      self._overlaps = group_args(
        n_clashes      = results.n_clashes,
        clashscore     = results.clashscore,
        n_clashes_sym  = results.n_clashes_sym,
        clashscore_sym = results.clashscore_sym,
        clashes_str    = string_io.getvalue(),
        clashes_dict   = clashes._clashes_dict)

    return self._overlaps
