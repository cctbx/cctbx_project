from __future__ import absolute_import, division, print_function
import sys, math
from cctbx.array_family import flex
from libtbx import adopt_init_args
from mmtbx import utils
import scitbx.math
from cctbx import adptbx
from cctbx import geometry_restraints

import six
from six.moves import range

def selector(hierarchy, ignore_hd=True):
  ensemble_asc = hierarchy.atom_selection_cache()
  ensemble_str = 'not resname HOH'
  if ignore_hd:
    ensemble_str += ' and not element H'
  ensemble_sel = ensemble_asc.selection(ensemble_str)
  ensemble_hierarchy = hierarchy.select(ensemble_sel)
  return ensemble_hierarchy

def get_sites_carts(hierarchys, ignore_hd=True):
  sites_carts = []
  for n, hierarchy in enumerate(hierarchys):
    ensemble_hierarchy = selector(hierarchy, ignore_hd=ignore_hd)
    sites_carts.append(ensemble_hierarchy.atoms().extract_xyz())
  return sites_carts

def get_selected_hierarchys(hierarchys, ignore_hd=True):
  sh = []
  for n, hierarchy in enumerate(hierarchys):
    ensemble_hierarchy = selector(hierarchy, ignore_hd=ignore_hd)
    sh.append(ensemble_hierarchy)
  return sh

def ensemble_mean_hierarchy(hierarchys,
                            ignore_hd=True,
                            verbose=False,
                            out=None,
                            ):
  sites_carts = get_sites_carts(hierarchys, ignore_hd=ignore_hd)
  mean_sites_cart = sites_carts[0].deep_copy()
  for n, sites_cart in enumerate(sites_carts):
    if not n: continue
    mean_sites_cart += sites_cart
  mean_sites_cart = mean_sites_cart * float(1./len(sites_carts))
  if verbose:
    print('Mean Structure Stats', file=out)
    for i, sites_cart in enumerate(sites_carts):
      print('  RMS to mean %3d : %0.2f' % (i+1,
                                           mean_sites_cart.rms_difference(sites_cart)),
            file=out)
  mean_hierarchy = hierarchys[0].deep_copy()
  mean_hierarchy = selector(mean_hierarchy, ignore_hd=ignore_hd)
  mean_hierarchy.atoms().set_xyz(mean_sites_cart)
  return mean_hierarchy

def closest_to_mean(hierarchys, mean_hierarchy, ignore_hd=True, verbose=False):
  sites_carts = get_sites_carts(hierarchys, ignore_hd=ignore_hd)
  mean_sites_cart = mean_hierarchy.atoms().extract_xyz()
  min_sites_cart = None
  min_rms = 1e9
  min_index = None
  for i, sites_cart in enumerate(sites_carts):
    rms = mean_sites_cart.rms_difference(sites_cart)
    if rms<min_rms:
      min_rms=rms
      min_sites_cart = sites_cart
      min_index = i
  assert min_sites_cart
  close_hierarchy = hierarchys[0].deep_copy()
  close_hierarchy = selector(close_hierarchy, ignore_hd=ignore_hd)
  close_hierarchy.atoms().set_xyz(min_sites_cart)
  if verbose:
    print('Closest to Mean\n  %3d : %0.2f (rms)' % (min_index+1, min_rms))
  return close_hierarchy, min_index

def get_centroid_hierarchy(hierarchys, ignore_hd=True, verbose=False):
  sites_carts = get_sites_carts(hierarchys, ignore_hd=ignore_hd)
  rmss = {}
  for i, sites_cart1 in enumerate(sites_carts):
    for j, sites_cart2 in enumerate(sites_carts):
      if i==j: break
      rms = sites_cart2.rms_difference(sites_cart1)
      rmss.setdefault(i, {})
      rmss[i][j]=rms
      rmss.setdefault(j, {})
      rmss[j][i]=rms
  min_sites_cart = None
  min_rms = 1e9
  min_index = None
  max_sites_cart = None
  max_rms = -1e9
  max_index = None
  for key, item in rmss.items():
    s = sum(item.values())
    if s<min_rms:
      min_rms=s
      min_sites_cart = sites_carts[key]
      min_index = key
    if s>max_rms:
      max_rms=s
      max_sites_cart = sites_carts[key]
      max_index = key
  assert min_sites_cart
  centroid_hierarchy = hierarchys[0].deep_copy()
  centroid_hierarchy = selector(centroid_hierarchy, ignore_hd=ignore_hd)
  centroid_hierarchy.atoms().set_xyz(min_sites_cart)
  least_hierarchy = hierarchys[0].deep_copy()
  least_hierarchy = selector(least_hierarchy, ignore_hd=ignore_hd)
  least_hierarchy.atoms().set_xyz(min_sites_cart)
  if verbose:
    print('Centroid Structure\n  %3d : %0.2f (average rms)' % (min_index+1, min_rms/len(rmss)))
    print('Extreme Structure\n  %3d : %0.2f (average rms)' % (max_index+1, max_rms/len(rmss)))
  return centroid_hierarchy, least_hierarchy, min_index, max_index

def dist2(xyz1, xyz2):
  d2=0
  for i in range(3):
    d2+=(xyz2[i]-xyz1[i])**2
  return d2

def get_rmsf_B_factor_per_residue_per_atom(hierarchys,
                                           reference,
                                           ignore_hd=True,
                                           verbose=False):
  sites_carts = get_sites_carts(hierarchys, ignore_hd=ignore_hd)
  sel_hierarchys = get_selected_hierarchys(hierarchys, ignore_hd=ignore_hd)
  ref_sites_cart = reference.atoms().extract_xyz()
  diff = []
  tempFactor = {}
  RMSF = {}
  for j, atom in enumerate(reference.atoms()):
    diff.append(0.)
    for sites_cart in sites_carts:
      # differences between the states and reference_state for each atom.
      d2 = dist2(atom.xyz, sites_cart[j])
      diff[j] += d2
    diff[j] = math.sqrt(diff[j]/len(hierarchys))
    key_atom = sel_hierarchys[0].atoms()[j]
    key = ((key_atom.parent().parent().parent().id,
            key_atom.parent().parent().resseq,
            key_atom.parent().resname,
            key_atom.name
           ))
    tempFactor.setdefault(key, [])

    tempFactor[key].append(diff[j]*8*math.pi**2)

    key = ((key_atom.parent().parent().parent().id,
            key_atom.parent().parent().resseq,
            key_atom.parent().resname,
           ))
    RMSF.setdefault(key, [])
    RMSF[key].append(diff[j])

  for key, item in RMSF.items():
    RMSF[key] =  sum(item)/len(item)
  return tempFactor, RMSF, diff

class manager(object):
  def __init__(self,
               ensemble_obj):
    adopt_init_args(self, locals())

  def print_atom_statisitics(self,
                             selection_string = 'name CA',
                             selection_bool_array = None,
                             model = None):
    if model == None:
      model = self.ensemble_obj.model
    pdb_hierarchy = model.get_hierarchy()
    pdb_atoms = pdb_hierarchy.atoms()
    if selection_bool_array == None:
      selection_bool_array = pdb_hierarchy.atom_selection_cache().selection(selection_string)
    xrs = model.get_xray_structure().deep_copy_scatterers()
    xrs.convert_to_isotropic()
    b_iso_atoms = xrs.scatterers().extract_u_iso()/adptbx.b_as_u(1)
    q_atoms = xrs.scatterers().extract_occupancies()
    for i_seq, x in enumerate(selection_bool_array):
      if x:
        atom_info = pdb_atoms[i_seq].fetch_labels()
        if self.ensemble_obj.er_data.ke_protein_running == None:
          print(' {0:6} {1:6} {2:6} {3:6} {4:6d} {5:6.3f} {6:6.3f} '.format(
                   atom_info.name,
                   atom_info.resseq,
                   atom_info.resname,
                   atom_info.chain_id,
                   i_seq,
                   b_iso_atoms[i_seq],
                   q_atoms[i_seq]), file=self.ensemble_obj.log)
        else:
          print(' {0:6} {1:6} {2:6} {3:6} {4:6d} {5:6.3f} {6:6.3f} {7:6.3f}'.format(
                   atom_info.name,
                   atom_info.resseq,
                   atom_info.resname,
                   atom_info.chain_id,
                   i_seq,
                   b_iso_atoms[i_seq],
                   q_atoms[i_seq],
                   self.ensemble_obj.er_data.ke_protein_running[i_seq]), file=self.ensemble_obj.log)

  #Print KE stats for simulation/model validation
  def kinetic_energy_stats(self):
    if self.ensemble_obj is not None:
      if self.ensemble_obj.er_data.ke_protein_running is not None:
        utils.print_header(
          line ="Non-solvent KE Statistics | MC : "+str(self.ensemble_obj.macro_cycle),
          out  = self.ensemble_obj.log)
        ke_basic_stats = scitbx.math.basic_statistics(self.ensemble_obj.er_data.ke_protein_running)
        print('  {0:<11}  {1:>12} {2:>12} {3:>12} {4:>12} {5:>12} '.format(
                                    '','min','max','mean', 'sdev', 'skew'), file=self.ensemble_obj.log)
        print('  KE MC {0:<5}: {1:12.3f} {2:12.3f} {3:12.3f} {4:12.3f} {5:12.3f}'.format(
                                    self.ensemble_obj.macro_cycle,
                                    ke_basic_stats.min,
                                    ke_basic_stats.max,
                                    ke_basic_stats.mean,
                                    ke_basic_stats.biased_standard_deviation,
                                    ke_basic_stats.skew), file=self.ensemble_obj.log)
        ke_atom_number_tuple_list = []
        ke_list_histo = []
        for n, ke in enumerate(self.ensemble_obj.er_data.ke_protein_running):
          ke_atom_number_tuple_list.append( (n, ke) )
          ke_list_histo.append(ke)
        assert len(ke_atom_number_tuple_list) == len(self.ensemble_obj.er_data.ke_protein_running)
        sorted_by_ke_ke_atom_number_tuple_list = \
          sorted(ke_atom_number_tuple_list, key=lambda ke:ke[-1])
        ke_list_histo = sorted(ke_list_histo)
        #Lowest KE Atoms
        pdb_atoms = self.ensemble_obj.pdb_hierarchy().atoms()
        print("\nNon-solvent atoms lowest KE : ", file=self.ensemble_obj.log)
        print('  {0:3} : {1:>44} {2:>12} {3:>12}'.format(
                'rank',
                'KE',
                'dmean/sdev',
                '%cum freq'), file=self.ensemble_obj.log)
        low_five_percent = (len(ke_atom_number_tuple_list) * 0.05)
        cntr = 0
        lowest_range = min(25, int(0.5 * len(ke_atom_number_tuple_list) ) )
        while cntr < lowest_range:
          atom_info = pdb_atoms[sorted_by_ke_ke_atom_number_tuple_list[cntr][0]].fetch_labels()
          assert atom_info.i_seq == sorted_by_ke_ke_atom_number_tuple_list[cntr][0]
          print(' {0:5} : {1:6} {2:6} {3:6} {4:6} {5:6} {6:9.3f} {7:12.3f} {8:12.1f}'.format(
                cntr+1,
                sorted_by_ke_ke_atom_number_tuple_list[cntr][0],
                atom_info.name,
                atom_info.resname,
                atom_info.chain_id,
                atom_info.resseq,
                sorted_by_ke_ke_atom_number_tuple_list[cntr][1],
                (sorted_by_ke_ke_atom_number_tuple_list[cntr][1]-ke_basic_stats.mean)\
                  / ke_basic_stats.biased_standard_deviation,
                100 * (float(cntr)/float(len(ke_atom_number_tuple_list))) ), file=self.ensemble_obj.log)
          cntr+=1
        #Highest KE Atoms
        print("\nNon-solvent atoms highest KE : ", file=self.ensemble_obj.log)
        print('  {0:3} : {1:>44} {2:>12} {3:>12}'.format(
                'rank',
                'KE',
                'dmean/sdev',
                '%cum freq'), file=self.ensemble_obj.log)
        cntr = len(ke_atom_number_tuple_list) - min(25, int(0.5 * len(ke_atom_number_tuple_list) ) )
        while cntr < len(ke_atom_number_tuple_list):
          atom_info = pdb_atoms[sorted_by_ke_ke_atom_number_tuple_list[cntr][0]].fetch_labels()
          assert atom_info.i_seq == sorted_by_ke_ke_atom_number_tuple_list[cntr][0]
          print(' {0:5} : {1:6} {2:6} {3:6} {4:6} {5:6} {6:9.3f} {7:12.3f} {8:12.1f}'.format(
                cntr+1,
                sorted_by_ke_ke_atom_number_tuple_list[cntr][0],
                atom_info.name,
                atom_info.resname,
                atom_info.chain_id,
                atom_info.resseq,
                sorted_by_ke_ke_atom_number_tuple_list[cntr][1],
                (sorted_by_ke_ke_atom_number_tuple_list[cntr][1]-ke_basic_stats.mean)\
                  / ke_basic_stats.biased_standard_deviation,
                100 * (float(cntr)/float(len(ke_atom_number_tuple_list))) ), file=self.ensemble_obj.log)
          cntr+=1
        #XXX Add print stats by for <ke>/residue

        #Histogram
        bin_list, bin_range_list = self.bin_generator_equal_range(
            array          = ke_list_histo[:-int(0.1 * (len(ke_list_histo) ) )],
            number_of_bins = 50)
        bin_range_list[-1][1] = max(ke_list_histo)
        self.bivariate_histogram(
          bin_array      = ke_list_histo,
          value_array    = ke_list_histo,
          name           = 'KE Histogram',
          bin_list       = bin_list,
          bin_range_list = bin_range_list)
        print("|"+"-"*77+"|\n", file=self.ensemble_obj.log)

  def bin_generator_equal_range(self, array, number_of_bins = 10):
    array_max   = max(array)
    array_min   = min(array)
    array_range = array_max - array_min
    bin_range   = array_range / (number_of_bins)
    bin_list = []
    bin_range_list = []
    val = array_min
    for x in range(number_of_bins):
      bin_list.append(val)
      bin_range_list.append([val, val+bin_range])
      val += bin_range
    return bin_list, bin_range_list

  def bin_generator_equal_size(self, array, number_of_bins = 8, return_bin_ave = False):
    sorted_array = sorted(array)
    bin_size = int(len(sorted_array) / number_of_bins+1)
    bin_list = []
    bin_range_list = []
    bin_ave_calc = []
    bin_average_array = flex.double()
    for n, x in enumerate(sorted_array):
      bin_ave_calc.append(x)
      if n%bin_size==0:
        bin_list.append(x)
        if n > 1:
          bin_average_array.append( sum(bin_ave_calc)/len(bin_ave_calc) )
          bin_ave_calc = []
        if len(bin_list) == number_of_bins:
          bin_range_list.append([sorted_array[n], sorted_array[-1]])
        else:
          bin_range_list.append([sorted_array[n],sorted_array[n+bin_size-1]])
    bin_average_array.append( sum(bin_ave_calc)/len(bin_ave_calc) )
    if return_bin_ave:
      return bin_list, bin_range_list, bin_average_array
    else:
      return bin_list, bin_range_list

  def bivariate_histogram(self, bin_array,
                                value_array,
                                name,
                                bin_list,
                                bin_range_list,
                                verbose = True,
                                return_bin_ave = False):
    if verbose:
      print("\nBivariate histogram "+name+" MC: "+str(self.ensemble_obj.macro_cycle), file=self.ensemble_obj.log)
    assert len(value_array) == len(bin_array)
    binned_array = [[]*1 for i in range(len(bin_list))]
    assert len(bin_list) == len(binned_array)
    for x in range(len(bin_array)):
      for bin_int, bin_value in enumerate(bin_list):
        flag = False
        if bin_array[x] < bin_value:
          binned_array[bin_int-1].append( value_array[x] )
          flag = True
          break
      if not flag:
        binned_array[-1].append( value_array[x] )
    if return_bin_ave:
      return_bin_array = flex.double()
    if verbose:
      print("{0:>20} | {1:>6} {2:>10} {3:>10} {4:>10} {5:>10}"\
        .format('Bin Range',
                'Freq',
                'Mean',
                'Min',
                'Max',
                'StdDev'), file=self.ensemble_obj.log)
    for n,bin_array in enumerate(binned_array):
      if len(bin_array) > 0:
        bin_array = flex.double(bin_array)
        if return_bin_ave:
          return_bin_array.append( sum(bin_array)/len(bin_array) )
        if verbose:
          print("{0:3d} {1:7.2f} -{2:7.2f} | {3:6d} {4:10.3f} {5:10.3f} {6:10.3f} {7:10.3f}".format(n+1,
                     bin_range_list[n][0],
                     bin_range_list[n][1],
                     len(bin_array),
                     sum(bin_array)/len(bin_array),
                     min(bin_array),
                     max(bin_array),
                     scitbx.math.basic_statistics(bin_array).biased_standard_deviation), file=self.ensemble_obj.log)
      else:
        if verbose:
          print("{0:3d} {1:7.2f} -{2:7.2f} | {3:6d} {4:10.3f} {5:10.3f} {6:10.3f} {7:10.3f}".format(n+1,
                     bin_range_list[n][0],
                     bin_range_list[n][1],
                     len(bin_array),
                     0,
                     0,
                     0,
                     0), file=self.ensemble_obj.log)
    if return_bin_ave:
      return return_bin_array

  def ensemble_reduction(self,
                         rfree_tolerance = 0.0025):
    #Reduces number of models to minimum required to reproduce Rfree
    utils.print_header("Ensemble reducer", out = self.ensemble_obj.log)
    self.ensemble_obj.show_overall(message        = "Full simulation fmodel final",
                             fmodel_running = False)
    final_rfree = self.ensemble_obj.fmodel_total.r_free()
    final_rwork = self.ensemble_obj.fmodel_total.r_work()

    # XXX no b_iso - how to apply this???
#    print >> self.ensemble_obj.log, "\nApply B_iso to all model in ensemble"
#    shift_b_iso  = self.ensemble_obj.fmodel_total.b_iso()
#    print >> self.ensemble_obj.log, 'Shift B_iso : {0:8.3f}'.format(shift_b_iso)
#    for x in self.ensemble_obj.er_data.xray_structures:
#      x.shift_us(b_shift = shift_b_iso)

    total_number_xrs = len(self.ensemble_obj.er_data.xray_structures)
    print("\nReduce ensemble with equal distribution though trajectory :", file=self.ensemble_obj.log)
    print("Rfree tolerance (%) : ", rfree_tolerance * 100, file=self.ensemble_obj.log)
    print('\n {0:>12} {1:>8} {2:>8} {3:>8}'\
        .format('Num','Rwork','Rfree','k1'), file=self.ensemble_obj.log)
    target_rfree = final_rfree
    final_div    = None
    for div_int in [1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000]:
      if div_int <= total_number_xrs:
        self.fmodel_ens = self.ensemble_obj.fmodel_total.deep_copy()
        cntr = 0.0
        fcalc_total = None
        fmask_total = None

  #      self.fmodel_ens.update(k_sols  = self.ensemble_obj.fmodel_total.k_sols(),
  #                             b_sol   = self.ensemble_obj.fmodel_total.b_sol(),
  #                             b_cart  = self.ensemble_obj.fmodel_total.b_cart() )

        for x in range(total_number_xrs):
          if x%int(div_int) == 0:
            #Apply back trace of Biso here...
            self.fmodel_ens.update_xray_structure(
              xray_structure      = self.ensemble_obj.er_data.xray_structures[x],
              update_f_calc       = True,
              update_f_mask       = True,
              force_update_f_mask = True)
            if fcalc_total == None:
              fcalc_total = self.fmodel_ens.f_calc().data().deep_copy()
              fmask_total = self.fmodel_ens.f_masks()[0].data().deep_copy()
              cntr = 1
            else:
              fcalc_total += self.fmodel_ens.f_calc().data().deep_copy()
              fmask_total += self.fmodel_ens.f_masks()[0].data().deep_copy()
              cntr += 1
          if x == total_number_xrs-1:
            self.fmodel_ens.update(
              f_calc = self.ensemble_obj.copy_ma.array(data = (fcalc_total / cntr)),
              f_mask = self.ensemble_obj.copy_ma.array(data = (fmask_total / cntr)) )
            self.fmodel_ens.update_all_scales(
              log    = self.ensemble_obj.log,
              remove_outliers=False,
              params = self.ensemble_obj.bsp)
            if cntr < 4:
              break
            print("Ens: {0:8d} {1:8.3f} {2:8.3f} {3:8.3f}"\
              .format(cntr,
                      self.fmodel_ens.r_work(),
                      self.fmodel_ens.r_free(),
                      self.fmodel_ens.scale_k1()
                      ), file=self.ensemble_obj.log)
            if self.fmodel_ens.r_free() < (target_rfree + rfree_tolerance):
              final_div    = div_int
              final_f_calc = self.ensemble_obj.copy_ma.array(data = (fcalc_total / cntr))
              final_f_mask = self.ensemble_obj.copy_ma.array(data = (fmask_total / cntr))
              if self.fmodel_ens.r_free() < target_rfree:
                target_rfree = self.fmodel_ens.r_free()

    if final_div == None:
      print("Warning pdb ensemble does not contain sufficent models and missrepresents simulation.  Simulation Rfree: {0:2.3f} %".format(100*(final_rfree)), file=self.ensemble_obj.log)
    else:
      #Update fmodel_total
      self.ensemble_obj.fmodel_total.update(f_calc = final_f_calc,
                                      f_mask = final_f_mask)
      self.ensemble_obj.fmodel_total.update_all_scales(
        log    = self.ensemble_obj.log,
        remove_outliers=False,
        params = self.ensemble_obj.bsp)
      #Parse arrays for output PDB
      copy_ed_data_xray_structures = []
      copy_pdb_hierarchys = []
      copy_ed_data_ke_pdb = []
      for x in range(len(self.ensemble_obj.er_data.xray_structures)):
        if x%int(final_div) == 0:
          copy_ed_data_xray_structures.append(self.ensemble_obj.er_data.xray_structures[x])
          copy_pdb_hierarchys.append(self.ensemble_obj.er_data.pdb_hierarchys[x])
          copy_ed_data_ke_pdb.append(self.ensemble_obj.er_data.ke_pdb[x])
      self.ensemble_obj.er_data.xray_structures          = copy_ed_data_xray_structures
      self.ensemble_obj.er_data.pdb_hierarchys           = copy_pdb_hierarchys
      self.ensemble_obj.er_data.ke_pdb                   = copy_ed_data_ke_pdb
      print("Final pdb ensemble contains {0:3d} models".format(len(self.ensemble_obj.er_data.xray_structures)), file=self.ensemble_obj.log)
      assert len(self.ensemble_obj.er_data.xray_structures) == len(self.ensemble_obj.er_data.pdb_hierarchys)
      assert len(self.ensemble_obj.er_data.xray_structures) == len(self.ensemble_obj.er_data.ke_pdb)

      print("|"+"-"*77+"|\n", file=self.ensemble_obj.log)

  def ensemble_rmsf_stats( self,
                           ensemble_hierarchys,
                           transfer_b_factors=True,
                           ignore_hd = True,
                           max_print=10,
                           verbose = False,
                           out = None,
                           ):
    if (out is None): out = sys.stdout
    if verbose:
      utils.print_header("Ensemble mean and centroid geometry statistics", out = out)
    ensemble_size = len(ensemble_hierarchys)
    self.mean_hierarchy = ensemble_mean_hierarchy( ensemble_hierarchys,
                                                   ignore_hd=ignore_hd,
                                                   verbose=verbose,
                                                   )
    close_hierarchy, self.closest_to_mean_index = closest_to_mean(
      ensemble_hierarchys,
      self.mean_hierarchy,
      ignore_hd=ignore_hd,
      verbose=verbose,
      )
    self.centroid_hierarchy, least_hierarchy, self.centroid_index, self.least_index = \
      get_centroid_hierarchy( ensemble_hierarchys,
                              ignore_hd=ignore_hd,
                              verbose=verbose,
                              )
    self.tempFactor, self.per_residue, self.per_atom = \
      get_rmsf_B_factor_per_residue_per_atom(
        ensemble_hierarchys,
        self.centroid_hierarchy,
        # mean_sites_cart,
        ignore_hd=ignore_hd,
        verbose=verbose,
        )

    if verbose:
      print('Per residue rmsf', file=out)
      for i, (key, item) in enumerate(self.per_residue.items()):
        print('  %5d : %s %0.2f' % (i,key,item), file=out)
        if i>=max_print: break

      print('B-factor', file=out)
      for i, (key, item) in enumerate(self.tempFactor.items()):
        print('  %5d : %s %7.2f' % (i,key,item[0]), file=out)
        if i>=max_print: break

      print('Per atom rmsf', file=out)
      for i, atom in enumerate(ensemble_hierarchys[0].atoms()):
        print('  %5d : %s %0.2f' % (i, atom.quote(), self.per_atom[i]), file=out)
        if i>=max_print: break

    if transfer_b_factors:
      atoms = self.centroid_hierarchy.atoms()
      occupancies = atoms.extract_occ()
      occupancies *= len(ensemble_hierarchys)
      atoms.set_occ(occupancies)
      for i, (key, item) in enumerate(self.tempFactor.items()):
        atom = atoms[i]
        atom.b = item[0]

    return self

  def write_mean_hierarchy(self, filename, crystal_symmetry):
    if not hasattr(self, 'mean_hierarchy'):
      assert 0, 'need to run ensemble_rmsf_stats'
    self.mean_hierarchy.write_pdb_file(
      filename,
      crystal_symmetry = crystal_symmetry,
    )

  def write_centroid_hierarchy(self, filename, crystal_symmetry):
    if not hasattr(self, 'centroid_hierarchy'):
      assert 0, 'need to run ensemble_rmsf_stats'
    self.centroid_hierarchy.write_pdb_file(
      filename,
      crystal_symmetry = crystal_symmetry,
    )

  def ensemble_mean_geometry_stats(self,
                                   restraints_manager,
                                   xray_structure,
                                   ensemble_xray_structures,
                                   ignore_hd = True,
                                   verbose = False,
                                   out = None,
                                   return_pdb_string = False):
    if (out is None): out = sys.stdout
    if verbose:
      utils.print_header("Ensemble mean geometry statistics", out = out)
    ensemble_size = len(ensemble_xray_structures)
    print("Ensemble size : ", ensemble_size, file=out)

    # Dictionaries to store deltas
    ensemble_bond_deltas = {}
    ensemble_angle_deltas = {}
    ensemble_chirality_deltas = {}
    ensemble_planarity_deltas = {}
    ensemble_dihedral_deltas = {}

    # List to store rmsd of each model
    structures_bond_rmsd = flex.double()
    structures_angle_rmsd = flex.double()
    structures_chirality_rmsd = flex.double()
    structures_planarity_rmsd = flex.double()
    structures_dihedral_rmsd = flex.double()

    # Remove water and hd atoms from global restraints manager
    selection = flex.bool()
    for sc in xray_structure.scatterers():
      if sc.label.find('HOH') > -1:
        selection.append(True)
      else:
        selection.append(False)
    if ignore_hd:
      hd_selection = xray_structure.hd_selection()
      assert hd_selection.size() == selection.size()
      for n in range(hd_selection.size()):
        if hd_selection[n] or selection[n]:
          selection[n] = True
    restraints_manager = restraints_manager.select(selection = ~selection)

    # Get all deltas
    for n, structure in enumerate(ensemble_xray_structures):
      if verbose:
        print("\nModel : ", n+1, file=out)
      sites_cart = structure.sites_cart()
      # Remove water and hd atoms from individual structures sites cart
      selection = flex.bool()
      for sc in structure.scatterers():
        if sc.label.find('HOH') > -1:
          selection.append(True)
        else:
          selection.append(False)
      if ignore_hd:
        hd_selection = structure.hd_selection()
        assert hd_selection.size() == selection.size()
        for n in range(hd_selection.size()):
          if hd_selection[n] or selection[n]:
            selection[n] = True
      sites_cart = sites_cart.select(~selection)
      assert sites_cart is not None
      site_labels = None
      energies_sites = restraints_manager.energies_sites(
          sites_cart        = sites_cart,
          compute_gradients = False)

      # Rmsd of individual model
      bond_rmsd = energies_sites.geometry.bond_deviations()[2]
      angle_rmsd = energies_sites.geometry.angle_deviations()[2]
      chirality_rmsd = energies_sites.geometry.chirality_deviations()[2]
      planarity_rmsd = energies_sites.geometry.planarity_deviations()[2]
      dihedral_rmsd = energies_sites.geometry.dihedral_deviations()[2]

      structures_bond_rmsd.append(bond_rmsd)
      structures_angle_rmsd.append(angle_rmsd)
      structures_chirality_rmsd.append(chirality_rmsd)
      structures_planarity_rmsd.append(planarity_rmsd)
      structures_dihedral_rmsd.append(dihedral_rmsd)

      if verbose:
        print("  Model RMSD", file=out)
        print("    bond      : %.6g" % bond_rmsd, file=out)
        print("    angle     : %.6g" % angle_rmsd, file=out)
        print("    chirality : %.6g" % chirality_rmsd, file=out)
        print("    planarity : %.6g" % planarity_rmsd, file=out)
        print("    dihedral  : %.6g" % dihedral_rmsd, file=out)

      # Bond
      pair_proxies = restraints_manager.geometry.pair_proxies(flags=None, sites_cart=sites_cart)
      assert pair_proxies is not None
      if verbose:
        pair_proxies.bond_proxies.show_histogram_of_deltas(
          sites_cart  = sites_cart,
          n_slots     = 10,
          f           = out)
      for proxy in pair_proxies.bond_proxies.simple:
        bond_simple_proxy = geometry_restraints.bond(
            sites_cart = sites_cart,
            proxy      = proxy)
        if proxy.i_seqs in ensemble_bond_deltas:
          ensemble_bond_deltas[proxy.i_seqs][0]+=bond_simple_proxy.delta
          ensemble_bond_deltas[proxy.i_seqs][1]+=1
        else:
          ensemble_bond_deltas[proxy.i_seqs] = [bond_simple_proxy.delta, 1]
        if verbose:
          print("bond simple :", proxy.i_seqs, file=out)
          print("  distance_ideal : %.6g" % proxy.distance_ideal, file=out)
          print("  distance_model : %.6g" % bond_simple_proxy.distance_model, file=out)
          print("  detla          : %.6g" % bond_simple_proxy.delta, file=out)
      if (pair_proxies.bond_proxies.asu.size() > 0):
        asu_mappings = pair_proxies.bond_proxies.asu_mappings()
        for proxy in pair_proxies.bond_proxies.asu:
          rt_mx = asu_mappings.get_rt_mx_ji(pair=proxy)
          bond_asu_proxy = geometry_restraints.bond(
              sites_cart   = sites_cart,
              asu_mappings = asu_mappings,
              proxy        = proxy)
          proxy_i_seqs = (proxy.i_seq, proxy.j_seq)
          if proxy_i_seqs in ensemble_bond_deltas:
            ensemble_bond_deltas[proxy_i_seqs][0]+=bond_asu_proxy.delta
            ensemble_bond_deltas[proxy_i_seqs][1]+=1
          else:
            ensemble_bond_deltas[proxy_i_seqs] = [bond_asu_proxy.delta, 1]
          if verbose:
            print("bond asu :", (proxy.i_seq, proxy.j_seq), rt_mx, file=out)
            print("  distance_ideal : %.6g" % proxy.distance_ideal, file=out)
            print("  distance_model : %.6g" % bond_asu_proxy.distance_model, file=out)
            print("  delta          : %.6g" % bond_asu_proxy.delta, file=out)

      # Angle
      if verbose:
        restraints_manager.geometry.angle_proxies.show_histogram_of_deltas(
            sites_cart  = sites_cart,
            n_slots     = 10,
            f           = out)
      for proxy in restraints_manager.geometry.angle_proxies:
        angle_proxy = geometry_restraints.angle(
            sites_cart = sites_cart,
            proxy      = proxy)
        if proxy.i_seqs in ensemble_angle_deltas:
          ensemble_angle_deltas[proxy.i_seqs][0]+=angle_proxy.delta
          ensemble_angle_deltas[proxy.i_seqs][1]+=1
        else:
          ensemble_angle_deltas[proxy.i_seqs] = [angle_proxy.delta, 1]
        if verbose:
          print("angle : ", proxy.i_seqs, file=out)
          print("  angle_ideal   : %.6g" % proxy.angle_ideal, file=out)
          print("  angle_model   : %.6g" % angle_proxy.angle_model, file=out)
          print("  delta         : %.6g" % angle_proxy.delta, file=out)

      # Chirality
      if verbose:
        restraints_manager.geometry.chirality_proxies.show_histogram_of_deltas(
            sites_cart  = sites_cart,
            n_slots     = 10,
            f           = out)
      for proxy in restraints_manager.geometry.chirality_proxies:
        chirality_proxy = geometry_restraints.chirality(
            sites_cart = sites_cart,
            proxy      = proxy)
        if proxy.i_seqs in ensemble_chirality_deltas:
          ensemble_chirality_deltas[proxy.i_seqs][0]+=chirality_proxy.delta
          ensemble_chirality_deltas[proxy.i_seqs][1]+=1
        else:
          ensemble_chirality_deltas[proxy.i_seqs] = [chirality_proxy.delta, 1]
        if verbose:
          print("chirality : ", proxy.i_seqs, file=out)
          print("  chirality_ideal : %.6g" % proxy.volume_ideal, file=out)
          print("  chirality_model : %.6g" % chirality_proxy.volume_model, file=out)
          print("  chirality       : %.6g" % chirality_proxy.delta, file=out)

      # Planarity
      for proxy in restraints_manager.geometry.planarity_proxies:
        planarity_proxy = geometry_restraints.planarity(
            sites_cart = sites_cart,
            proxy      = proxy)
        proxy_i_seqs = []
        for i_seq in proxy.i_seqs:
          proxy_i_seqs.append(i_seq)
        proxy_i_seqs = tuple(proxy_i_seqs)
        if proxy_i_seqs in ensemble_planarity_deltas:
          ensemble_planarity_deltas[proxy_i_seqs][0]+=planarity_proxy.rms_deltas()
          ensemble_planarity_deltas[proxy_i_seqs][1]+=1
        else:
          ensemble_planarity_deltas[proxy_i_seqs] = [planarity_proxy.rms_deltas(), 1]
        if verbose:
          print("planarity : ", proxy_i_seqs, file=out)
          print("  planarity rms_deltas : %.6g" % planarity_proxy.rms_deltas(), file=out)

      # Dihedral
      if verbose:
        restraints_manager.geometry.dihedral_proxies.show_histogram_of_deltas(
            sites_cart  = sites_cart,
            n_slots     = 10,
            f           = out)
      for proxy in restraints_manager.geometry.dihedral_proxies:
        dihedral_proxy = geometry_restraints.dihedral(
            sites_cart = sites_cart,
            proxy      = proxy)
        if proxy.i_seqs in ensemble_dihedral_deltas:
          ensemble_dihedral_deltas[proxy.i_seqs][0]+=dihedral_proxy.delta
          ensemble_dihedral_deltas[proxy.i_seqs][1]+=1
        else:
          ensemble_dihedral_deltas[proxy.i_seqs] = [dihedral_proxy.delta, 1]
        if verbose:
          print("dihedral : ", proxy.i_seqs, file=out)
          print("  dihedral_ideal  : %.6g" % proxy.angle_ideal, file=out)
          print("  periodicity     : %.6g" % proxy.periodicity, file=out)
          print("  dihedral_model  : %.6g" % dihedral_proxy.angle_model, file=out)
          print("  delta           : %.6g" % dihedral_proxy.delta, file=out)

    # Calculate RMSDs for ensemble model
    # Bond
    mean_bond_delta = flex.double()
    for proxy, info in six.iteritems(ensemble_bond_deltas):
      # assert info[1] == ensemble_size
      if info[1]!=ensemble_size:
        print('skipping bond RMSD calns of ensemble %s' % info, file=out)
        continue
      mean_delta = info[0] / info[1]
      mean_bond_delta.append(mean_delta)
    bond_delta_sq = mean_bond_delta * mean_bond_delta
    ensemble_bond_rmsd = math.sqrt(flex.mean_default(bond_delta_sq, 0))

    # Angle
    mean_angle_delta = flex.double()
    for proxy, info in six.iteritems(ensemble_angle_deltas):
      assert info[1] == ensemble_size
      mean_delta = info[0] / info[1]
      mean_angle_delta.append(mean_delta)
    angle_delta_sq = mean_angle_delta * mean_angle_delta
    ensemble_angle_rmsd = math.sqrt(flex.mean_default(angle_delta_sq, 0))

    # Chirality
    mean_chirality_delta = flex.double()
    for proxy, info in six.iteritems(ensemble_chirality_deltas):
      assert info[1] == ensemble_size
      mean_delta = info[0] / info[1]
      mean_chirality_delta.append(mean_delta)
    chirality_delta_sq = mean_chirality_delta * mean_chirality_delta
    ensemble_chirality_rmsd = math.sqrt(flex.mean_default(chirality_delta_sq, 0))

    # Planarity
    mean_planarity_delta = flex.double()
    for proxy, info in six.iteritems(ensemble_planarity_deltas):
      assert info[1] == ensemble_size
      mean_delta = info[0] / info[1]
      mean_planarity_delta.append(mean_delta)
    planarity_delta_sq = mean_planarity_delta * mean_planarity_delta
    ensemble_planarity_rmsd = math.sqrt(flex.mean_default(planarity_delta_sq, 0))

    # Dihedral
    mean_dihedral_delta = flex.double()
    for proxy, info in six.iteritems(ensemble_dihedral_deltas):
      assert info[1] == ensemble_size
      mean_delta = info[0] / info[1]
      mean_dihedral_delta.append(mean_delta)
    dihedral_delta_sq = mean_dihedral_delta * mean_dihedral_delta
    ensemble_dihedral_rmsd = math.sqrt(flex.mean_default(dihedral_delta_sq, 0))

    # Calculate <structure rmsd>
    assert ensemble_size == structures_bond_rmsd
    assert ensemble_size == structures_angle_rmsd
    assert ensemble_size == structures_chirality_rmsd
    assert ensemble_size == structures_planarity_rmsd
    assert ensemble_size == structures_dihedral_rmsd
    structure_bond_rmsd_mean = structures_bond_rmsd.min_max_mean().mean
    structure_angle_rmsd_mean = structures_angle_rmsd.min_max_mean().mean
    structure_chirality_rmsd_mean = structures_chirality_rmsd.min_max_mean().mean
    structure_planarity_rmsd_mean = structures_planarity_rmsd.min_max_mean().mean
    structure_dihedral_rmsd_mean = structures_dihedral_rmsd.min_max_mean().mean

    # Show summary
    utils.print_header("Ensemble RMSD summary", out = out)
    print("  RMSD (mean delta per restraint)", file=out)
    print("    bond      : %.6g" % ensemble_bond_rmsd, file=out)
    print("    angle     : %.6g" % ensemble_angle_rmsd, file=out)
    print("    chirality : %.6g" % ensemble_chirality_rmsd, file=out)
    print("    planarity : %.6g" % ensemble_planarity_rmsd, file=out)
    print("    dihedral  : %.6g" % ensemble_dihedral_rmsd, file=out)
    print("  RMSD (mean RMSD per structure)", file=out)
    print("    bond      : %.6g" % structure_bond_rmsd_mean, file=out)
    print("    angle     : %.6g" % structure_angle_rmsd_mean, file=out)
    print("    chirality : %.6g" % structure_chirality_rmsd_mean, file=out)
    print("    planarity : %.6g" % structure_planarity_rmsd_mean, file=out)
    print("    dihedral  : %.6g" % structure_dihedral_rmsd_mean, file=out)
    if ignore_hd:
      print("\n  Calculated excluding H/D", file=out)
    else:
      print("\n  Calculated including H/D", file=out)

    if return_pdb_string:
      ens_geo_pdb_string  = "REMARK   3"
      ens_geo_pdb_string += "\nREMARK   3  NUMBER STRUCTURES IN ENSEMBLE : {0:5d}".format(ensemble_size)
      if ignore_hd:
        ens_geo_pdb_string += "\nREMARK   3  RMS DEVIATIONS FROM IDEAL VALUES (EXCLUDING H/D)"
      else:
        ens_geo_pdb_string += "\nREMARK   3  RMS DEVIATIONS FROM IDEAL VALUES (INCLUDING H/D)"
      ens_geo_pdb_string += "\nREMARK   3  RMSD (MEAN DELTA PER RESTRAINT)"
      ens_geo_pdb_string += "\nREMARK   3    BOND      : {0:5.3f}".format(ensemble_bond_rmsd)
      ens_geo_pdb_string += "\nREMARK   3    ANGLE     : {0:5.3f}".format(ensemble_angle_rmsd)
      ens_geo_pdb_string += "\nREMARK   3    CHIRALITY : {0:5.3f}".format(ensemble_chirality_rmsd)
      ens_geo_pdb_string += "\nREMARK   3    PLANARITY : {0:5.3f}".format(ensemble_planarity_rmsd)
      ens_geo_pdb_string += "\nREMARK   3    DIHEDRAL  : {0:5.2f}".format(ensemble_dihedral_rmsd)
      ens_geo_pdb_string += "\nREMARK   3  RMSD (MEAN RMSD PER STRUCTURE)"
      ens_geo_pdb_string += "\nREMARK   3    BOND      : {0:5.3f}".format(structure_bond_rmsd_mean)
      ens_geo_pdb_string += "\nREMARK   3    ANGLE     : {0:5.3f}".format(structure_angle_rmsd_mean)
      ens_geo_pdb_string += "\nREMARK   3    CHIRALITY : {0:5.3f}".format(structure_chirality_rmsd_mean)
      ens_geo_pdb_string += "\nREMARK   3    PLANARITY : {0:5.3f}".format(structure_planarity_rmsd_mean)
      ens_geo_pdb_string += "\nREMARK   3    DIHEDRAL  : {0:5.2f}".format(structure_dihedral_rmsd_mean)
      ens_geo_pdb_string += "\nREMARK   3"
      return ens_geo_pdb_string

if __name__ == '__main__':
  from iotbx import pdb

  def get_pdb_hierarchies(pdb_file_names):
    pdb_hierarchys = []
    if len(pdb_file_names)==1:
      ensemble_filename = pdb_file_names[0]
      pdb_inp = pdb.input(ensemble_filename)
      pdb_hierarchy = pdb_inp.construct_hierarchy()
      for i, model in enumerate(pdb_hierarchy.models()):
        pdb_hierarchys.append(pdb_hierarchy.deep_copy())
        for j in range(len(pdb_hierarchy.models())-1,-1,-1):
          if j!=i:
            pdb_hierarchys[-1].remove_model(j)
        # pdb_hierarchys[-1].write_pdb_file('test_%s.pdb' % i)
    else:
      for i, pdb_file_name in enumerate(pdb_file_names):
        print(i,pdb_file_name)
        pdb_inp = pdb.input(pdb_file_name)
        pdb_hierarchy = pdb_inp.construct_hierarchy()
        pdb_hierarchys.append(pdb_hierarchy)
    return pdb_inp, pdb_hierarchys

  erm = manager(None)
  print('Ensemble Refinement Manager')
  ensemble_filenames = sys.argv[1:]
  print('ensemble_filename', ensemble_filenames)
  pdb_input, pdb_hierarchys = get_pdb_hierarchies(ensemble_filenames)
  print('Number of models : %s' % (len(pdb_hierarchys)))
  erm.ensemble_rmsf_stats(pdb_hierarchys)
  erm.write_centroid_hierarchy('centroid.pdb', pdb_input.crystal_symmetry())
