from __future__ import absolute_import, division, print_function
import time, sys
import re
from six.moves import cStringIO as StringIO
import iotbx.pdb
from mmtbx import map_tools
from cctbx import maptbx
import cctbx.geometry_restraints.process_nonbonded_proxies as pnp
from cctbx import adptbx
from iotbx import phil
from cctbx.array_family import flex
from libtbx import group_args
from cctbx import miller
from libtbx.str_utils import make_sub_header
import mmtbx.maps.polder
import mmtbx.maps.correlation

master_params_str = """
validate_ligands {

nproc = 1
  .type = int

}
"""

def master_params():
  return phil.parse(master_params_str, process_includes = False)

# =============================================================================

class manager(list):
  '''
  Manager class for all ligands
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

  # ----------------------------------------------------------------------------

  def parallel_populate(self, args):
    from libtbx import easy_mp
    def _run(lr, func):
      func = getattr(lr, func)
      rc = func()
      return rc
    funcs = []
    ligand_results = []
    inputs = []
    for ligand_isel, sel_str in args:
      lr = ligand_result(
        model       = self.model,
        fmodel      = self.fmodel,
        ligand_isel = ligand_isel,
        sel_str     = sel_str)
      ligand_results.append(lr)
      for attr, func in lr._result_attrs.items():
        funcs.append([lr, func])
        inputs.append(func)

    results = []
    t0=time.time()
    for i, (args, res, err_str) in enumerate(easy_mp.multi_core_run(
      _run,
      funcs,
      self.params.nproc,
      )):
      results.append([args, res, err_str])
      if self.params.nproc>1:
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

  # ----------------------------------------------------------------------------

  def run(self):
    args = []
    for ligand_isel, sel_str in self.generate_ligand_iselections():
      args.append([ligand_isel, sel_str])
    results = self.parallel_populate(args)
    for r in results:
      self.append(r)

  # ----------------------------------------------------------------------------

  def generate_ligand_iselections(self):
    '''
    Get iselections for each ligand
    '''
    ph = self.model.get_hierarchy()
    get_class = iotbx.pdb.common_residue_names_get_class
    exclude = ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
               "modified_rna_dna", "ccp4_mon_lib_rna_dna", "common_water",
                "common_element"]
    for model in ph.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          for conformer in rg.conformers():
            residue = conformer.only_residue()
            resname = residue.resname
            if (not get_class(name=resname) in exclude):
              iselection = residue.atoms().extract_i_seq()
              sel_str = 'chain %s and resseq %s and resname %s ' % (chain.id,
                rg.resseq_as_int(), resname)
              if conformer.altloc:
                sel_str = sel_str + ' and altloc %s' % conformer.altloc
              #print(sel_str)
              yield iselection, sel_str

  # ----------------------------------------------------------------------------

  def show_ligand_counts(self):
    make_sub_header(' Ligands in input model ', out=self.log)
    for ligand_result in self:
      names = []
      for _a in ligand_result._atoms_ligand:
        names.append(_a.name.strip())
      print(ligand_result.id_str, '    ' + ', '.join(names), file=self.log)


  # ----------------------------------------------------------------------------

  def show_table(self, out):
      '''
      Print summary table
      '''
      lab_row1 =  ['','','', '% bad', '','', '', '', '', '']
      lab_row2 =  ['','','CC', 'map values', '','', 'ADPs', 'occupancies',
        'bond', 'angle']
      lab_row3 =  ['ligand', 'suspicious', '2Fo-Fc', 'Fo-Fc','clashes','H-bonds',\
        'min   max   mean   owab', 'min   max   mean', ' rmsz  outliers ',
        ' rmsz  outliers ']
      lab1_str = '{:^14}|{:^12}|{:^9}|{:^12}|{:^9}|{:^9}|{:^28}|{:^21}|{:^17}|{:^17}|'
      print('\n' + lab1_str.format(*lab_row1), file=out)
      print(lab1_str.format(*lab_row2), file=out)
      print(lab1_str.format(*lab_row3), file=out)
      print('-'*157, file=out)
      for lr in self:
        ccs     = lr.get_ccs()
        clashes = lr.get_overlaps()
        adps    = lr.get_adps()
        owab    = lr.get_owab()
        occs    = lr.get_occupancies()
        map_vals = lr.get_map_values()
        rmsds   = lr.get_rmsds()
        is_suspicious = lr.check_if_suspicious()
        n_atoms = lr._atoms_ligand_noH.size()
        if is_suspicious: check='***'
        else:             check =''

        #line = [lr.id_str,check, ccs.cc_2fofc, clashes.n_clashes,
        #  round(adps.b_min,1), round(adps.b_max,1), round(adps.b_mean,1)]
        #print(table_str.format(*line), file=self.log)

        val_id      = f"{lr.id_str:^14}"
        val_check   = f"{check:^12}"
        val_cc      = f"{ccs.cc_2fofc:^9.2f}" if ccs is not None else f"{'':^9}"
        val_mapvals = f"{round((map_vals.fofc_map_values<=-3).count(True)/n_atoms, 2):^12}" if map_vals is not None else f"{'':^12}"
        val_clash   = f"{clashes.n_clashes:^9}" if clashes.n_clashes != 0 else f"{'-':^9}"
        val_hbonds  = f"{clashes.n_hbonds:^9}" if clashes.n_hbonds != 0 else f"{'-':^9}"
        val_b_min   = f"{round(adps.b_min,1):^7}"
        val_b_max   = f"{round(adps.b_max,1):^7}"
        val_b_mean  = f"{round(adps.b_mean,1):^7}"
        val_owab    = f"{round(owab,1):^7}"
        val_o_min   = f"{round(occs.occ_min,1):^7}" if occs.occ_min != occs.occ_mean else f"{'':^7}"
        val_o_max   = f"{round(occs.occ_max,1):^7}" if occs.occ_max != occs.occ_mean else f"{'':^7}"
        val_o_mean  = f"{round(occs.occ_mean,1):^7}"
        val_bond_str = str(round(rmsds.bond_rmsz,2)) + '  ' + str(rmsds.bond_n_outliers) + \
          '(' + str(rmsds.bond_n) + ')'
        val_bond  = f"{val_bond_str:^17}"
        val_angle_str = str(round(rmsds.angle_rmsz,2)) + '  ' + str(rmsds.angle_n_outliers) + \
          '(' + str(rmsds.angle_n) + ')'
        val_angle  = f"{val_angle_str:^17}"

        row = f"{val_id}|{val_check}|{val_cc}|{val_mapvals}|{val_clash}|\
{val_hbonds}|{val_b_min}{val_b_max}{val_b_mean}{val_owab}|\
{val_o_min}{val_o_max}{val_o_mean}|{val_bond}|{val_angle}|"
        print(row, file=out)
        #
        val_sites = f"{'sites':^14}"
        _b_min_within = f"{round(adps.b_min_within,1):^7}" if adps.b_min_within is not None else f"{'':^7}"
        _b_max_within = f"{round(adps.b_max_within,1):^7}" if adps.b_max_within is not None else f"{'':^7}"
        _b_mean_within = f"{round(adps.b_mean_within,1):^7}" if adps.b_mean_within is not None else f"{'':^7}"
        sites_row = f"{val_sites}|{'':^12}|{'':^9}|{'':^12}|{'':^9}|{'':^9}|\
{_b_min_within}{_b_max_within}{_b_mean_within}\
{'':^7}|{'':^7}{'':^7}{'':^7}|"
        print(sites_row, file=out)

  # def show_adps(self):
  #   '''
  #   Show results for ADPs of ligand and surrounding atoms
  #   '''
  #   make_sub_header(' ADPs ', out=self.log)
  #   pad1 = ' '*18
  #   print(pad1, "min   max    mean   n_iso   n_aniso", file=self.log)
  #   for id_tuple, ligand_dict in self.items():
  #     for altloc, lr in ligand_dict.items():
  #       adps = lr.get_adps()
  #       print(lr.id_str.ljust(14), '%7s%7s%7s%7s%7s' %
  #         (round(adps.b_min,1), round(adps.b_max,1), round(adps.b_mean,1),
  #          adps.n_iso, adps.n_aniso), file = self.log)
  #       if (adps.b_mean_within is not None):
  #         print('neighbors'.ljust(14), '%7s%7s%7s' %
  #           (round(adps.b_min_within,1), round(adps.b_max_within,1),
  #            round(adps.b_mean_within,1) ), file = self.log)


  # def show_ligand_occupancies(self):
  #   '''
  #   Show results for ligand occupancies
  #   '''
  #   make_sub_header(' Occupancies ', out=self.log)
  #   pad1 = ' '*20
  #   print('If three values: min, max, mean, otherwise the same occupancy for entire ligand.', \
  #     file=self.log)
  #   for id_tuple, ligand_dict in self.items():
  #     for altloc, lr in ligand_dict.items():
  #       occs = lr.get_occupancies()
  #       if (occs.occ_min == occs.occ_max):
  #         print(lr.id_str.ljust(16), occs.occ_min, file = self.log)
  #       else:
  #         print(lr.id_str.ljust(16), '%s   %s   %s' %
  #           (occs.occ_min, occs.occ_max, occs.occ_mean), file = self.log)


  # def show_ccs(self):
  #   '''
  #   Show results for correlation coefficients
  #   '''
  #   if self.fmodel is None: return
  #   make_sub_header(' Correlation coefficients ', out=self.log)
  #   for id_tuple, ligand_dict in self.items():
  #     for altloc, lr in ligand_dict.items():
  #       ccs = lr.get_ccs()
  #       cc_two_fofc = round(ccs.cc_two_fofc, 2)
  #       cc_fofc = round(ccs.cc_fofc, 2)
  #       fofc_min  = round(ccs.fofc_min, 2)
  #       fofc_max  = round(ccs.fofc_max, 2)
  #       fofc_mean = round(ccs.fofc_mean, 2)
  #       print(lr.id_str.ljust(16),
  #         cc_two_fofc, cc_fofc, fofc_min, fofc_max, fofc_mean, file = self.log)


  # def show_nonbonded_overlaps(self):
  #   '''
  #   Print results for overlaps
  #   '''
  #   for id_tuple, ligand_dict in self.items():
  #     for altloc, lr in ligand_dict.items():
  #       clashes_result = lr.get_overlaps()
  #       print(clashes_result.clashes_str, file=self.log)

  def show_sites_within(self):
    make_sub_header(' Sites within 3 A', out=self.log)
    for lr in self:
      adps    = lr.get_adps()
      isel_within_noH = adps.isel_within_noH
      ph_within = lr._ph.select(isel_within_noH)
      print(lr.id_str, file=self.log)
      for rg in ph_within.residue_groups():
        for c in rg.conformers():
          print('    ' + c.only_residue().id_str().split('"')[1], file=self.log)

# =============================================================================

class ligand_result(object):
  '''
  Class that stores validation info per ligand
  '''
  def __init__(self,
               model,
               fmodel,
               ligand_isel,
               sel_str):
    self.model       = model
    self.fmodel      = fmodel
    self.ligand_isel = ligand_isel
    self.sel_str     = sel_str

    # results
    self._result_attrs = {
      '_occupancies'   : 'get_occupancies',
      '_adps'          : 'get_adps',
      '_owab'          : 'get_owab',
      '_overlaps'      : 'get_overlaps',
      '_rmsds'         : 'get_rmsds',
      '_ccs'           : 'get_ccs',
      '_is_suspicious' : 'check_if_suspicious',
      '_map_values'    : 'get_map_values',
      #'_polder_ccs'  : 'get_polder_ccs',
    }

    self._set_internals()

    for attr, func in self._result_attrs.items():
      setattr(self, attr, None)
      assert hasattr(self, func)



  # ----------------------------------------------------------------------------

  def __repr__(self):
    outl = 'ligand %s\n' % self.id_str
    for attr in self._result_attrs:
      outl += '  %s : %s\n' % (attr, getattr(self, attr))
    return outl

  # ----------------------------------------------------------------------------

  def check_if_suspicious(self):
    '''
    If ligand metrics fulfil certain criteria, it is flagged as suspicious
    '''
    if self._is_suspicious is not None:
      return self._is_suspicious
    #
    self._is_suspicious = False
    # get info
    n_atoms = self._atoms_ligand_noH.size()
    occs = self.get_occupancies()
    adps = self.get_adps()
    clashes = self.get_overlaps()
    ccs = None
    map_vals = None

    if self.fmodel is not None:
      ccs = self.get_ccs()
      map_vals = self.get_map_values()
    # apply criteria for metrics
    if ccs is not None:
      if ccs.cc_2fofc < 0.5:
        self._is_suspicious = True
    if adps.b_mean_within is not None:
      if adps.b_mean > 4 * adps.b_mean_within:
        self._is_suspicious = True
    if occs.occ_mean == 0.:
      self._is_suspicious = True
    if clashes.n_clashes > 0.5 * n_atoms:
      self._is_suspicious = True
    if map_vals is not None:
      if (map_vals.fofc_map_values<-3).count(True) >= 0.5 * n_atoms:
        self._is_suspicious = True
    #
    return self._is_suspicious

  # ----------------------------------------------------------------------------

  def get_rmsds(self):

    if self._rmsds is not None:
      return self._rmsds
    # ligand without H atoms
    model_ligand = self.model.select(self.ligand_isel_noH)
    stats = model_ligand.geometry_statistics()

    bond = stats.bond(origin_id=0)
    #print('number bond outliers', len(bond.outliers))
    #print('bond rmsd',bond.mean)
    bond_z = stats.bond(origin_id=0, return_rmsZ=True)
    #print('bond rmsz',bond_z.mean)

    angle = stats.angle(origin_id=0)
    #print('number angle outliers', len(angle.outliers))
    #print('angle rmsd',angle.mean)
    angle_z = stats.angle(origin_id=0, return_rmsZ=True)
    #print('angle rmsz',angle_z.mean)

    dihedral = stats.dihedral()
    #print('number dihedral outliers', len(dihedral.outliers))
    #print('dihedral rmsd',dihedral.mean)
    dihedral_z = stats.dihedral(return_rmsZ=True)

    plane = stats.planarity()
    #print('plane rmsd', plane.mean)

    # Note that geo.bond.outliers will be empty list, need to supply origin_id
    # in stats.bond() to get list of outliers
    # So the stats.result() object does not have the list of outliers
    # geo = stats.result()

    self._rmsds = group_args(
      bond_rmsd  = bond.mean,
      bond_rmsz  = bond_z.mean,
      bond_n     = bond.n,
      bond_n_outliers = len(bond.outliers),
      angle_rmsd = angle.mean,
      angle_rmsz = angle_z.mean,
      angle_n     = angle.n,
      angle_n_outliers = len(angle.outliers),
      #chirality_rmsd = geo.chirality.mean,
      planarity_rmsd = plane.mean,
      dihedral_rmsd = dihedral.mean,
      dihedral_rmsz = dihedral_z.mean,
      dihedral_n     = dihedral.n,
      dihedral_n_outliers = len(dihedral.outliers),)
    return self._rmsds

  # ----------------------------------------------------------------------------

  def get_occupancies(self):
    '''
    Get occupancies of non-H atoms
    '''
    if self._occupancies is not None:
      return self._occupancies
    eps = 1.e-6
    occ = self._atoms_ligand_noH.extract_occ()
    occ_mmm = occ.min_max_mean()

    self._occupancies = group_args(
      occ_min             = occ_mmm.min,
      occ_max             = occ_mmm.max,
      occ_mean            = occ_mmm.mean,
      negative_count      = (occ<0).count(True),
      negative_isel       = (occ<0).iselection(),
      zero_count          = (flex.abs(occ)<eps).count(True),
      zero_isel           = (flex.abs(occ)<eps).iselection(),
      less_than_dot9_isel = (occ<0.9).iselection()
      )

    return self._occupancies

  # ----------------------------------------------------------------------------

  def get_adps(self):
    '''
    Get isotropic B-factors of non-H atoms
    '''
    if self._adps is not None:
      return self._adps
    b_isos  = self._xrs_ligand_noH.extract_u_iso_or_u_equiv() * adptbx.u_as_b(1.)
    n_iso   = self._xrs_ligand_noH.use_u_iso().count(True)
    n_aniso = self._xrs_ligand_noH.use_u_aniso().count(True)
    n_zero  = (b_isos < 0.01).count(True)
    #n_above_100 = (b_isos > 100).count(True)
    #isel_above_100 = (b_isos > 100).iselection()
    b_min, b_max, b_mean = b_isos.min_max_mean().as_tuple()

    within_radius = 3.0 #TODO should this be a parameter?
    # if ligand has alternative conformation, ignore it
    #if 'altloc' in self.sel_str:
    #  import re
    #  s = re.sub(r'altloc \w\s*(and\s*)?', '', self.sel_str)
    #  _sel_str = s.strip()
    #else:
    #  _sel_str = self.sel_str
    _sel_str = self.sel_str

    sel_within_str_noH = '(residues_within (%s, %s)) and protein and not water \
    and not (element H or element D) and not (%s)' % \
    (within_radius, self.sel_str, _sel_str)
    #print(sel_within_str_noH)
    isel_within_noH = self.model.iselection(sel_within_str_noH)
    xrs_within_noH = self._xrs.select(isel_within_noH)
    b_isos_within = xrs_within_noH.extract_u_iso_or_u_equiv() * adptbx.u_as_b(1.)
    b_min_within, b_max_within, b_mean_within = b_isos_within.min_max_mean().as_tuple()

    self._adps = group_args(
      n_iso           = n_iso,
      n_aniso         = n_aniso,
      n_zero          = n_zero,
      b_min           = b_min,
      b_max           = b_max,
      b_mean          = b_mean,
      b_min_within    = b_min_within,
      b_max_within    = b_max_within,
      b_mean_within   = b_mean_within,
      isel_within_noH = isel_within_noH
      )

    return self._adps

  # ----------------------------------------------------------------------------

  def get_owab(self):
    '''
    Compute occupancy weighted average B-factor (owab)
    '''
    if self._owab is not None:
      return self._owab
    #
    eps = 1.e-6
    occ = self._atoms_ligand_noH.extract_occ()
    b_isos  = self._xrs_ligand_noH.extract_u_iso_or_u_equiv() * adptbx.u_as_b(1.)
    sum_occ = flex.sum(occ)
    # yes, sum_occ=0 really happens: 2ace, 1lmc, 1lrl, 1v2u, 2ou9
    if sum_occ < eps:
      sum_occ = 1.e-6
    owab = flex.sum(occ*b_isos)/sum_occ
    self._owab = owab

    return self._owab

  # ----------------------------------------------------------------------------

  def _set_internals(self):
    self._ph = self.model.get_hierarchy()
    self._atoms_ligand = self._ph.select(self.ligand_isel).atoms()
    self._xrs = self.model.get_xray_structure()
    self._xrs_ligand = self._xrs.select(self.ligand_isel)

    for rg in self._ph.select(self.ligand_isel).residue_groups():
      for c in rg.conformers():
        _resname = c.only_residue().resname

    self.resname = _resname.strip()

    _id_str = self._atoms_ligand[0].id_str()
    if _id_str.startswith("model"):
      _id_str = _id_str.split('pdb="')[1].split('"')[0]
    else:
      _id_str = _id_str.split('"')[1]
#    altloc = _id_str[4]
#    resseq = _id_str[10:14]
#    chain  = _id_str[8:10]
#    self.sel_str = " ".join(['chain', chain, 'and resseq', resseq, 'and resname', _resname])
#    if (altloc != ' '):
#      self.sel_str = " ".join(['altloc', altloc, 'and', self.sel_str])
    _id_str = _id_str.strip().split(' ')
    self.id_str = " ".join(_id_str[1:]).strip()
    #
    _noH = ' and not (element H or element D)'
    self.ligand_isel_noH = self.model.iselection(self.sel_str + _noH)

    self._xrs_ligand_noH = self._xrs.select(self.ligand_isel_noH)
    #self._xrs_ligand_noH = \
    #  self.model.select(self.ligand_isel_noH).get_xray_structure()
    self._atoms_ligand_noH = self._ph.select(self.ligand_isel_noH).atoms()

  # ----------------------------------------------------------------------------

  def get_polder_ccs(self):
    if self._polder_ccs is not None:
      return self._polder_ccs

    from mmtbx.maps.polder import master_params_str as polder_params_str
    _params = iotbx.phil.parse(
      input_string=polder_params_str, process_includes=True).extract()

    polder_object = mmtbx.maps.polder.compute_polder_map(
      f_obs            = self.fmodel.f_obs(),
      r_free_flags     = self.fmodel.r_free_flags(),
      model            = self.model,
      params           = _params.polder,
      selection_string = self.sel_str )

    polder_object.validate()
    polder_object.run()
    r = polder_object.get_results()
    vr = r.validation_results
    print('Map 1: calculated Fobs with ligand')
    print('Map 2: calculated Fobs without ligand')
    print('Map 3: real Fobs data')
    print('CC(1,2): %6.4f' % vr.cc12)
    print('CC(1,3): %6.4f' % vr.cc13)
    print('CC(2,3): %6.4f' % vr.cc23)

    self._polder_ccs = group_args(
      cc12 = vr.cc12,
      cc13 = vr.cc13,
      cc23 = vr.cc23
      )

    return self._polder_ccs

  # ----------------------------------------------------------------------------

  def get_ccs(self):
    if self.fmodel is None:
      return
    if self._ccs is not None:
      return self._ccs

    cs = self.fmodel.f_obs().crystal_symmetry()
    crystal_gridding = maptbx.crystal_gridding(
      unit_cell        = cs.unit_cell(),
      space_group_info = cs.space_group_info(),
      symmetry_flags   = maptbx.use_space_group_symmetry,
      step             = 0.6)

    # Dfmodel map including ligand
    #fmodel.update_xray_structure(
    #  xray_structure = self.model.get_xray_structure(),
    #  update_f_calc=True
    #)
    m2 = self.compute_maps(
      fmodel           = self.fmodel,
      crystal_gridding = crystal_gridding,
      map_type         = "DFmodel")

    fmodel = self.fmodel.deep_copy()
    sel = self.model.selection(string=self.sel_str)
    # 2mFo-DFc map without ligand
    fmodel.update_xray_structure(
      xray_structure = fmodel.xray_structure.select(~sel),
      update_f_calc=True
    )
    #fmodel.show_short(show_k_mask=True, log=None, prefix="")
    m1 = self.compute_maps(
      fmodel           = fmodel,
      crystal_gridding = crystal_gridding,
      map_type         = "2mFo-DFc")

    sites_cart = self.model.get_sites_cart().select(sel)
    sel = maptbx.grid_indices_around_sites(
      unit_cell  = cs.unit_cell(),
      fft_n_real = m1.focus(),
      fft_m_real = m1.all(),
      sites_cart = sites_cart,
      site_radii = flex.double(sites_cart.size(), 1.0))
    m1 = m1.set_selected(m1<0, 0)
    m2 = m2.set_selected(m1<0, 0)
    cc = flex.linear_correlation(
      x=m1.select(sel).as_1d(),
      y=m2.select(sel).as_1d()).coefficient()

    self._ccs = group_args(
      cc_2fofc = cc,
       )
    return self._ccs

  # ----------------------------------------------------------------------------

  def compute_maps(self, fmodel, crystal_gridding, map_type):
    map_coefficients = map_tools.electron_density_map(
      fmodel = fmodel).map_coefficients(
        map_type         = map_type,
        isotropize       = True,
        fill_missing     = False)
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = map_coefficients)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded()

  # ----------------------------------------------------------------------------

  def get_map_values(self):
    if self.fmodel is None:
      return
    if self._map_values is not None:
      return self._map_values

    fofc_map_values = flex.double()

#   # get map coefficients
#    mc = map_tools.electron_density_map(fmodel = self.fmodel).map_coefficients(
#        map_type         = "mFo-DFc",
#        isotropize       = True,
#        fill_missing     = False)
#    # TODO d_min and cg should be accessible in class
#    d_min = self.fmodel.f_obs().d_min()
#    cg = self.fmodel.f_obs().crystal_gridding(
#      d_min             = d_min,
#      symmetry_flags    = maptbx.use_space_group_symmetry,
#      resolution_factor = 0.25)
#    # get map_data
#    map_fo = miller.fft_map(
#        crystal_gridding     = cg,
#        fourier_coefficients = mc)
#    map_fo.apply_sigma_scaling()
#    map_data_fofc = map_fo.real_map_unpadded()

    # could make crystal gridding class attribute?
    cs = self.fmodel.f_obs().crystal_symmetry()
    crystal_gridding = maptbx.crystal_gridding(
      unit_cell        = cs.unit_cell(),
      space_group_info = cs.space_group_info(),
      symmetry_flags   = maptbx.use_space_group_symmetry,
      step             = 0.6)
    # mFo-DFc map with ligand present
    m_fofc = self.compute_maps(
      fmodel           = self.fmodel,
      crystal_gridding = crystal_gridding,
      map_type         = "2mFo-DFc")

    # get fo-fc map values at atom centers: fo-fc
    unit_cell = self.model.crystal_symmetry().unit_cell()
    for site_cart, _a in zip(self._atoms_ligand_noH.extract_xyz(), self._atoms_ligand_noH):
      site_frac = unit_cell.fractionalize(site_cart)
      map_val = m_fofc.eight_point_interpolation(site_frac)
      #print(_a.id_str(), map_val)
      fofc_map_values.append(map_val)

    self._map_values = group_args(
      fofc_map_values = fofc_map_values,
       )

    return self._map_values

  # ----------------------------------------------------------------------------

  def get_overlaps(self):
    '''
    Obtain overlaps involving ligands
    '''
    # A model with H atoms is necessary to process overlaps
    if not self.model.has_hd():
      return None
    if self._overlaps is not None:
      return self._overlaps

    within_radius = 3.0

    # TODO clashes with other ligands?
    sel_within_str = '%s or (residues_within (%s, %s)) and (protein or water)' \
      % (self.sel_str, within_radius, self.sel_str)
    #print(sel_within_str)

    sel_within = self.model.selection(sel_within_str)
    #print(sel_within.count(True))
    model_within = self.model.select(sel_within)
    isel_ligand_within = model_within.iselection(self.sel_str)


    ##isel_ligand_within = sel_within.iselection()
    #isel_ligand_within = self.model.select(sel_within).iselection(self.sel_str)
    ##sel = flex.bool([True]*len(sel_within))
    #model_within = self.model.select(sel_within)
    # debug
    #_id_str = self.id_str.replace(" ", "_")
    #fn = "site_%s.pdb" % _id_str
    #clean_filename = re.sub(r"\s+", "_", fn)
    #f = open(clean_filename,"w")
    #f.write(model_within.model_as_pdb())
    #f.close()
    # debug end

    processed_nbps = pnp.manager(model = model_within)
    clashes = processed_nbps.get_clashes()
    hbonds = processed_nbps.get_hbonds()

    clashes_dict   = clashes._clashes_dict
    hbonds_dict = hbonds._hbonds_dict

    ligand_clashes_dict = {}
    for iseq_tuple, record in clashes_dict.items():
      if (iseq_tuple[0] in isel_ligand_within or
          iseq_tuple[1] in isel_ligand_within):
        ligand_clashes_dict[iseq_tuple] = record

    ligand_clashes = pnp.clashes(
                    clashes_dict = ligand_clashes_dict,
                    model        = model_within)

    ligand_hbonds_dict = {}
    for iseq_tuple, record in hbonds_dict.items():
      if (iseq_tuple[0] in isel_ligand_within or
          iseq_tuple[1] in isel_ligand_within):
        ligand_hbonds_dict[iseq_tuple] = record

    ligand_hbonds = pnp.hbonds(
                    hbonds_dict  = ligand_hbonds_dict,
                    model        = model_within)

    results_hbonds = ligand_hbonds.get_results()

    #clashes.show(log=sys.stdout)
    #hbonds.show(log=sys.stdout)

    #string_io = StringIO()
    #ligand_clashes.show(log=string_io, show_clashscore=False)
    #print(string_io.getvalue())

    results = ligand_clashes.get_results()

    self._overlaps = group_args(
      n_clashes      = results.n_clashes,
      clashscore     = results.clashscore,
      n_clashes_sym  = results.n_clashes_sym,
      #clashscore_sym = results.clashscore_sym,
      #clashes_str    = string_io.getvalue(),
      #clashes_dict   = clashes._clashes_dict,
      n_hbonds = results_hbonds.n_hbonds)

    return self._overlaps
