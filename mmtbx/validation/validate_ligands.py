from __future__ import division, print_function

import time, os
import iotbx.pdb
import mmtbx.model
import cctbx.geometry_restraints.nonbonded_overlaps as nbo
from cctbx import adptbx
from iotbx import phil
from cctbx.array_family import flex
from libtbx import group_args
from libtbx.str_utils import make_sub_header
from libtbx.utils import null_out
from libtbx import easy_run
from mmtbx import monomer_library

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
# Manager class for ALL ligands

class manager(dict):

  def __init__(self,
               model,
               params,
               model_fn,
               log=None):
    self.model = model
    self.params = params
    self.log   = log
    #
    self.nproc = params.nproc
    # TODO:  tidy up filename requirement once ready set is refactored
    self.model_fn = model_fn

  # ---------------------------------------------------------------------------

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
        readyset_model = self.readyset_model,
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

  # ---------------------------------------------------------------------------

  def run(self):
    args = []
    self.get_readyset_model_with_grm()

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

  # ---------------------------------------------------------------------------

  def get_readyset_model_with_grm(self):
    # TODO: user should have possibility to use existing H atoms,
    # in this case, how to run readyset (because if might be still necessary to get
    # ligand cif file)
    self.readyset_model = None
    fn_pdb, fn_cif = self.run_ready_set(file_name = self.model_fn)
    if fn_cif is not None:
      cif_object = monomer_library.server.read_cif(file_name=fn_cif)
      cif_objects = [(fn_cif, cif_object)]
    else:
      cif_objects = []
    if fn_pdb is not None:
      pdb_inp = iotbx.pdb.input(file_name=fn_pdb)
      self.readyset_model = mmtbx.model.manager(
        model_input = pdb_inp,
        restraint_objects = cif_objects,
        build_grm   = True,
        log         = null_out())

  # ---------------------------------------------------------------------------

  def run_ready_set(self, file_name):
    print('\nRunning ready_set to get ligand cif_files...')
    import libtbx.load_env
    has_ready_set = libtbx.env.has_module(name="phenix")
    if not has_ready_set:
      raise_sorry('phenix.ready_set could not be detected on your system.')
    if (self.params.place_hydrogens):
      cmd = "phenix.ready_set {} --silent".format(file_name)
      print('Placing H atoms...')
    else:
      cmd = "phenix.ready_set {} hydrogens=False --silent".format(file_name)
      print('H atoms are NOT added to the model...')
    print(cmd)
    out = easy_run.fully_buffered(cmd)
    if (out.return_code != 0):
      msg_str = "ready_set crashed - dumping stderr:\n%s"
      raise RuntimeError(msg_str % ( "\n".join(out.stderr_lines)))
    fn_pdb = file_name.replace('.pdb','.updated.pdb')
    fn_cif = file_name.replace('.pdb','.ligands.cif')
    if (os.path.isfile(fn_cif)):
      print('Ligand(s) cif file created: ', fn_cif)
    else:
      fn_cif = None
    if (os.path.isfile(fn_pdb)):
      print('Updated model file: ', fn_pdb)
    else:
      fn_pdb = None
    return fn_pdb, fn_cif

  # ---------------------------------------------------------------------------

  def get_ligands(self, ph):
    # Store ligands as list of iselections --> better way? Careful if H will be
    # added at some point!
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

  # ---------------------------------------------------------------------------

  def show_ligand_counts(self):
    make_sub_header(' Ligands in input model ', out=self.log)
    for id_tuple, ligand_dict in self.items():
      for altloc, lr in ligand_dict.items():
        print(lr.id_str)

  # ---------------------------------------------------------------------------

  def show_adps(self):
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

  # ---------------------------------------------------------------------------

  def show_ligand_occupancies(self):
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

  # ---------------------------------------------------------------------------

  def show_nonbonded_overlaps(self):
    make_sub_header(' Nonbonded overlaps', out=self.log)
    overlaps_found = False
    out_str = '{:>16}-{:>16}     {:>6.3f}       {:>6.3f}   {:^10}'
    for id_tuple, ligand_dict in self.items():
      for altloc, lr in ligand_dict.items():
        nbo = lr.get_nonbonded_overlaps()
        nbo_proxies = nbo.result.nb_overlaps_proxies_all
        out_list = []
        argmented_counts = [0,0]
        def _adjust_count(d):
          d=(abs(d)-0.4)*10
          return d+1
        if (len(nbo_proxies) > 0):
          overlaps_found = True
          for data in nbo_proxies:
            d = list(data)
            rec_list = [x.replace('pdb=','') for x in d[0]]
            rec_list = [x.replace('"','') for x in rec_list]
            #rec_list.extend(d[1:3])
            rec_list.append(d[3])
            rec_list.append(d[4])
            rec_list.append('1'*bool(d[5]) + ' '*(not bool(d[5])))
            ptr = 0
            if rec_list[4].strip(): ptr=1
            #argmented_counts[ptr] += _adjust_count(rec_list[2])
            out_list.append(out_str.format(*rec_list))
          out_string = '\n'.join(out_list)
    if overlaps_found:
      lbl_str = '{:^33} {:^11} {:^11} {:<10}'
      print(lbl_str.format(*["Overlapping atoms","model distance", "vdW distance",
                     "sym overlap"]), file=self.log)
      print('-'*73, file=self.log)
      print(out_string, file=self.log)
    else:
      print('No clashes found', file=self.log)
        #print argmented_counts
#        nbo.show(
#          log=self.log,
#          nbo_type='all',
#          normalized_nbo=False)

# =============================================================================
# Class storing info per ligand

class ligand_result(object):

  def __init__(self, model, readyset_model, isel, id_list):
    self.model = model
    self.isel = isel
    self.readyset_model = readyset_model

    # results
    self._result_attrs = {'_occupancies' : 'get_occupancies',
                          '_adps'        : 'get_adps',
                          '_nb_overlaps' : 'get_nonbonded_overlaps'
    }
    for attr, func in self._result_attrs.items():
      setattr(self, attr, None)
      assert hasattr(self, func)
    # to be used internally
    self._ph = self.model.get_hierarchy()
    self._atoms = self._ph.select(self.isel).atoms()
    self._xrs = self.model.select(isel).get_xray_structure()
    self._get_id_str(id_list=id_list)

  # ---------------------------------------------------------------------------

  def __repr__(self):
    outl = 'ligand %s\n' % self.id_str
    for attr in self._result_attrs:
      outl += '  %s : %s\n' % (attr, getattr(self, attr))
    return outl

  # ---------------------------------------------------------------------------

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

  # ---------------------------------------------------------------------------

  def get_adps(self):
    if self._adps is None:
      b_isos = self._xrs.extract_u_iso_or_u_equiv() * adptbx.u_as_b(1.)
      n_iso   = self._xrs.use_u_iso().count(True)
      n_aniso = self._xrs.use_u_aniso().count(True)
      n_zero = (b_isos < 0.01).count(True)
      # TODO: what number as cutoff?
      n_above_100 = (b_isos > 100).count(True)
      isel_above_100 = (b_isos > 100).iselection()
      b_min, b_max, b_mean = b_isos.min_max_mean().as_tuple()

      #if this selection is used somewhere else, it might be better to do it outside
      within_radius = 3.0 #TODO should this be a parameter?
      sel_within_str = '(within (%s, %s)) and (protein or water)' % (within_radius, self.sel_str)
      isel_within = self.model.iselection(sel_within_str)
      xrs_within = self.model.select(isel_within).get_xray_structure()
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

  # ---------------------------------------------------------------------------

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

  # ---------------------------------------------------------------------------

  def get_nonbonded_overlaps(self):
    # this function should use a model with H atoms!!
    if self.readyset_model is None:
      return None
    if self._nb_overlaps is None:
      # sel_within could be done in __init__?
      within_radius = 3.0
      sel_within_str = '%s or (within (%s, %s)) and (protein or water)' \
        % (self.sel_str, within_radius, self.sel_str)
      #isel_within = self.readyset_model.iselection(sel_within_str)
      sel_within = self.readyset_model.selection(sel_within_str)
      sel_ligand_within = self.readyset_model.select(sel_within).selection(self.sel_str)
      xrs_within = self.readyset_model.select(sel_within).get_xray_structure()
      #xrs = self.model.get_xray_structure()
      geometry = self.readyset_model.get_restraints_manager().select(sel_within).geometry

      #sel = flex.bool([True]*len(sel_within))

      self._nb_overlaps = nbo.info(
        geometry_restraints_manager=geometry,
        macro_molecule_selection=~sel_ligand_within, #TODO should be only ligand, not sel_within?
        sites_cart=xrs_within.sites_cart(),
        site_labels=xrs_within.scatterers().extract_labels(),
        hd_sel=xrs_within.hd_selection())

    return self._nb_overlaps

