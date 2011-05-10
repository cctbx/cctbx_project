from cctbx.array_family import flex
from libtbx import adopt_init_args
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
from mmtbx import utils
import scitbx.math

class manager(object):
  def __init__(self,
               ta_obj):
    adopt_init_args(self, locals())

  #Print KE stats for simulation/model validation
  def kinetic_energy_stats(self):
    if self.ta_obj is not None:
      if self.ta_obj.tad.ke_protein_running is not None:
        utils.print_header(
          line ="Non-solvent KE Statistics | MC : "+str(self.ta_obj.macro_cycle),
          out  = self.ta_obj.log)
        ke_basic_stats = scitbx.math.basic_statistics(self.ta_obj.tad.ke_protein_running)
        print >> self.ta_obj.log, '  {0:<11}  {1:>12} {2:>12} {3:>12} {4:>12} {5:>12} '.format(
                                    '','min','max','mean', 'sdev', 'skew')
        print >> self.ta_obj.log, '  KE MC {0:<5}: {1:12.3f} {2:12.3f} {3:12.3f} {4:12.3f} {5:12.3f}'.format(
                                    self.ta_obj.macro_cycle,
                                    ke_basic_stats.min,
                                    ke_basic_stats.max,
                                    ke_basic_stats.mean,
                                    ke_basic_stats.biased_standard_deviation,
                                    ke_basic_stats.skew)
        ke_atom_number_tuple_list = []
        ke_list_histo = []
        for n, ke in enumerate(self.ta_obj.tad.ke_protein_running):
          ke_atom_number_tuple_list.append( (n, ke) )
          ke_list_histo.append(ke)
        assert len(ke_atom_number_tuple_list) == len(self.ta_obj.tad.ke_protein_running)
        sorted_by_ke_ke_atom_number_tuple_list = \
          sorted(ke_atom_number_tuple_list, key=lambda ke:ke[-1])
        ke_list_histo = sorted(ke_list_histo)
        #Lowest KE Atoms
        pdb_atoms = self.ta_obj.pdb_hierarchy.atoms()
        print >> self.ta_obj.log, "\nNon-solvent atoms lowest KE : "
        print >> self.ta_obj.log, \
              '  {0:3} : {1:>44} {2:>12} {3:>12}'.format(
                'rank',
                'KE',
                'dmean/sdev',
                '%cum freq')
        low_five_percent = (len(ke_atom_number_tuple_list) * 0.05)
        cntr = 0
        lowest_range = min(25, int(0.5 * len(ke_atom_number_tuple_list) ) )
        while cntr < lowest_range:
          atom_info = pdb_atoms[sorted_by_ke_ke_atom_number_tuple_list[cntr][0]].fetch_labels()
          assert atom_info.i_seq == sorted_by_ke_ke_atom_number_tuple_list[cntr][0]
          print >> self.ta_obj.log, \
              ' {0:5} : {1:6} {2:6} {3:6} {4:6} {5:6} {6:9.3f} {7:12.3f} {8:12.1f}'.format(
                cntr+1,
                sorted_by_ke_ke_atom_number_tuple_list[cntr][0],
                atom_info.name,
                atom_info.resname,
                atom_info.chain_id,
                atom_info.resseq,
                sorted_by_ke_ke_atom_number_tuple_list[cntr][1],
                (sorted_by_ke_ke_atom_number_tuple_list[cntr][1]-ke_basic_stats.mean)\
                  / ke_basic_stats.biased_standard_deviation,
                100 * (float(cntr)/float(len(ke_atom_number_tuple_list))) )
          cntr+=1
        #Highest KE Atoms
        print >> self.ta_obj.log, "\nNon-solvent atoms highest KE : "
        print >> self.ta_obj.log, \
              '  {0:3} : {1:>44} {2:>12} {3:>12}'.format(
                'rank',
                'KE',
                'dmean/sdev',
                '%cum freq')
        cntr = len(ke_atom_number_tuple_list) - min(25, int(0.5 * len(ke_atom_number_tuple_list) ) )
        while cntr < len(ke_atom_number_tuple_list):
          atom_info = pdb_atoms[sorted_by_ke_ke_atom_number_tuple_list[cntr][0]].fetch_labels()
          assert atom_info.i_seq == sorted_by_ke_ke_atom_number_tuple_list[cntr][0]
          print >> self.ta_obj.log, \
              ' {0:5} : {1:6} {2:6} {3:6} {4:6} {5:6} {6:9.3f} {7:12.3f} {8:12.1f}'.format(
                cntr+1,
                sorted_by_ke_ke_atom_number_tuple_list[cntr][0],
                atom_info.name,
                atom_info.resname,
                atom_info.chain_id,
                atom_info.resseq,
                sorted_by_ke_ke_atom_number_tuple_list[cntr][1],
                (sorted_by_ke_ke_atom_number_tuple_list[cntr][1]-ke_basic_stats.mean)\
                  / ke_basic_stats.biased_standard_deviation,
                100 * (float(cntr)/float(len(ke_atom_number_tuple_list))) )
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
        print >> self.ta_obj.log, "|"+"-"*77+"|\n"

  def bin_generator_equal_range(self, array, number_of_bins = 10):
    array_max   = max(array)
    array_min   = min(array)
    array_range = array_max - array_min
    bin_range   = array_range / (number_of_bins)
    bin_list = []
    bin_range_list = []
    val = array_min
    for x in xrange(number_of_bins):
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
      print >> self.ta_obj.log, "\nBivariate histogram "+name+" MC: "+str(self.ta_obj.macro_cycle)
    assert len(value_array) == len(bin_array)
    binned_array = [[]*1 for i in xrange(len(bin_list))]
    assert len(bin_list) == len(binned_array)
    for x in xrange(len(bin_array)):
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
      print >> self.ta_obj.log, "{0:>20} | {1:>6} {2:>10} {3:>10} {4:>10} {5:>10}"\
        .format('Bin Range',
                'Freq',
                'Mean',
                'Min',
                'Max',
                'StdDev')
    for n,bin_array in enumerate(binned_array):
      if len(bin_array) > 0:
        bin_array = flex.double(bin_array)
        if return_bin_ave:
          return_bin_array.append( sum(bin_array)/len(bin_array) )
        if verbose:
          print >> self.ta_obj.log, "{0:3d} {1:7.2f} -{2:7.2f} | {3:6d} {4:10.3f} {5:10.3f} {6:10.3f} {7:10.3f}".format(n+1,
                     bin_range_list[n][0],
                     bin_range_list[n][1],
                     len(bin_array),
                     sum(bin_array)/len(bin_array),
                     min(bin_array),
                     max(bin_array),
                     scitbx.math.basic_statistics(bin_array).biased_standard_deviation)
      else:
        if verbose:
          print >> self.ta_obj.log, "{0:3d} {1:7.2f} -{2:7.2f} | {3:6d} {4:10.3f} {5:10.3f} {6:10.3f} {7:10.3f}".format(n+1,
                     bin_range_list[n][0],
                     bin_range_list[n][1],
                     len(bin_array),
                     0,
                     0,
                     0,
                     0)
    if return_bin_ave:
      return return_bin_array

  def ensemble_reduction(self):
    #Reduces number of models to minimum required to reproduce Rfree
    utils.print_header("Time-averaging ensemble reducer", out = self.ta_obj.log)
    self.ta_obj.show_overall(message        = "Full simulation fmodel final",
                             fmodel_running = False)
    final_rfree = self.ta_obj.fmodel_total.r_free()
    final_rwork = self.ta_obj.fmodel_total.r_work()
    print >> self.ta_obj.log, "\nApply B_iso to all model in ensemble"
    shift_b_iso  = self.ta_obj.fmodel_total.b_iso()
    print >> self.ta_obj.log, 'Shift B_iso : {0:8.3f}'.format(shift_b_iso)
    for x in self.ta_obj.tad.xray_structures:
      x.shift_us(b_shift = shift_b_iso)
    total_number_xrs = len(self.ta_obj.tad.xray_structures)
    div_int =  1
    print >> self.ta_obj.log, "\nReduce ensemble with equal distribution though trajectory :"
    print >> self.ta_obj.log, '\n {0:>12} {1:>8} {2:>8} {3:>8} {4:>8} {5:>8} {6:>8}'\
        .format('Num','Rwork','Rfree','k1','biso','ksol','bsol')
    target_rfree = final_rfree
    final_div    = None
    while div_int <= total_number_xrs:
      self.fmodel_ens = self.ta_obj.fmodel_total.deep_copy()
      cntr = 0.0
      fcalc_total = None
      fmask_total = None
      self.fmodel_ens.update(k_sol  = self.ta_obj.fmodel_total.k_sol(),
                             b_sol  = self.ta_obj.fmodel_total.b_sol(),
                             b_cart = self.ta_obj.fmodel_total.b_cart() )
      for x in xrange(total_number_xrs):
        if x%int(div_int) == 0:
          #Apply back trace of Biso here...
          self.fmodel_ens.update_xray_structure(
            xray_structure      = self.ta_obj.tad.xray_structures[x],
            update_f_calc       = True,
            update_f_mask       = True,
            force_update_f_mask = True)
          if fcalc_total == None:
            fcalc_total = self.fmodel_ens.f_calc().data().deep_copy()
            fmasks = self.fmodel_ens.shell_f_masks()
            assert len(fmasks) == 1
            fmask_total = fmasks[0].data().deep_copy()
            cntr = 1
          else:
            fcalc_total += self.fmodel_ens.f_calc().data().deep_copy()
            fmasks = self.fmodel_ens.shell_f_masks()
            assert len(fmasks) == 1
            fmask_total += fmasks[0].data().deep_copy()
            cntr += 1
        if x == total_number_xrs-1:
          self.fmodel_ens.update(
            f_calc = self.ta_obj.copy_ma.array(data = (fcalc_total / cntr)),
            f_mask = self.ta_obj.copy_ma.array(data = (fmask_total / cntr)) )
          self.fmodel_ens.update_solvent_and_scale(
                                   verbose       = 0,
                                   out           = self.ta_obj.log,
                                   params        = self.ta_obj.bsp,
                                   optimize_mask = False)
          print >> self.ta_obj.log, "Ens: {0:8d} {1:8.3f} {2:8.3f} {3:8.3f} {4:8.3f} {5:8.3f} {6:8.3f}"\
            .format(cntr,
                    100*self.fmodel_ens.r_work(),
                    100*self.fmodel_ens.r_free(),
                    self.fmodel_ens.scale_k1(),
                    self.fmodel_ens.b_iso(),
                    self.fmodel_ens.k_sol(),
                    self.fmodel_ens.b_sol() )

          #Tolerate 0.1% reduction in Rfree
          if self.fmodel_ens.r_free() < (target_rfree+0.001):
            final_div    = div_int
            final_f_calc = self.ta_obj.copy_ma.array(data = (fcalc_total / cntr))
            final_f_mask = self.ta_obj.copy_ma.array(data = (fmask_total / cntr))
            if self.fmodel_ens.r_free() < target_rfree:
              target_rfree = self.fmodel_ens.r_free()
      div_int = div_int * 2.0

    if final_div == None:
      print >> self.ta_obj.log, "Warning pdb ensemble does not contain sufficent models and missrepresents simulation.  Simulation Rfree: {0:2.3f} %".format(100*(final_rfree))
    else:
      #Update fmodel_total
      self.ta_obj.fmodel_total.update(f_calc = final_f_calc,
                                      f_mask = final_f_mask)
      self.ta_obj.fmodel_total.update_solvent_and_scale(
                                verbose       = 0,
                                out           = self.ta_obj.log,
                                params        = self.ta_obj.bsp,
                                optimize_mask = False)
      #Parse arrays for output PDB
      copy_tad_xray_structures = []
      copy_pdb_hierarchys      = []
      copy_tad_ke_pdb          = []
      copy_tad_geo_a_b_c_d_n_p = []
      copy_tad_geo_a_b_c_d_n_p_n_res = []
      for x in xrange(len(self.ta_obj.tad.xray_structures)):
        if x%int(final_div) == 0:
          copy_tad_xray_structures.append(self.ta_obj.tad.xray_structures[x])
          copy_pdb_hierarchys.append(self.ta_obj.tad.pdb_hierarchys[x])
          copy_tad_ke_pdb.append(self.ta_obj.tad.ke_pdb[x])
          copy_tad_geo_a_b_c_d_n_p.append(self.ta_obj.tad.geo_a_b_c_d_n_p[x])
          copy_tad_geo_a_b_c_d_n_p_n_res.append(self.ta_obj.tad.geo_a_b_c_d_n_p_n_res[x])
      self.ta_obj.tad.xray_structures        = copy_tad_xray_structures
      self.ta_obj.tad.pdb_hierarchys         = copy_pdb_hierarchys
      self.ta_obj.tad.ke_pdb                 = copy_tad_ke_pdb
      self.ta_obj.tad.geo_a_b_c_d_n_p        = copy_tad_geo_a_b_c_d_n_p
      self.ta_obj.tad.geo_a_b_c_d_n_p_n_res  = copy_tad_geo_a_b_c_d_n_p_n_res
      print >> self.ta_obj.log, "Final pdb ensemble contains {0:3d} models".format(len(self.ta_obj.tad.xray_structures))
      assert len(self.ta_obj.tad.xray_structures) == len(self.ta_obj.tad.pdb_hierarchys)
      assert len(self.ta_obj.tad.xray_structures) == len(self.ta_obj.tad.ke_pdb)
      assert len(self.ta_obj.tad.xray_structures) == len(self.ta_obj.tad.geo_a_b_c_d_n_p)
      assert len(self.ta_obj.tad.xray_structures) == len(self.ta_obj.tad.geo_a_b_c_d_n_p_n_res)
      print >> self.ta_obj.log, "|"+"-"*77+"|\n"

  def ensemble_to_fmodel(self):
    #TST function for different recombination methods
    utils.print_header("Time-averaging ensemble pruner", out = self.ta_obj.log)
    total_number_xrs = len(self.ta_obj.tad.xray_structures)
    print >> self.ta_obj.log, "Total number of stored structures : ", total_number_xrs
    self.ta_obj.fmodel_ens = self.ta_obj.fmodel_running.deep_copy()
    print >> self.ta_obj.log, "\nRandom distribution though trajectory (with enrichment): \n"
    print >> self.ta_obj.log, "20% structures with greatest individual Rfree discarded"

    indi_rfree = flex.double()
    # For use with best structure ens below
    xrs_rfree_tuple_list = []
    for x in xrange(total_number_xrs):
      self.ta_obj.fmodel_ens.update_xray_structure(
        xray_structure = self.ta_obj.tad.xray_structures[x],
        update_f_calc  = True,
        update_f_mask  = True,
        force_update_f_mask = True)
      self.ta_obj.fmodel_ens.update_solvent_and_scale(verbose = 0, out = self.ta_obj.log, optimize_mask = False)
      print >> self.ta_obj.log, x, self.ta_obj.fmodel_ens.r_work(), self.ta_obj.fmodel_ens.r_free()
      indi_rfree.append(self.ta_obj.fmodel_ens.r_free())
      xrs_rfree_tuple_list.append( (self.ta_obj.tad.xray_structures[x],self.ta_obj.fmodel_ens.r_free()) )

    sorted_indi_rfree = sorted(indi_rfree)
    print >> self.ta_obj.log, "80% cutoff : ", sorted_indi_rfree[int(total_number_xrs*0.8)]
    best_xrs_sel = flex.bool(indi_rfree < sorted_indi_rfree[int(total_number_xrs*0.8)])

    best_xrs = []
    for x in xrange(total_number_xrs):
      if best_xrs_sel[x]:
        best_xrs.append(self.ta_obj.tad.xray_structures[x])

    print >> self.ta_obj.log, "Total number xrs    : ", len(self.ta_obj.tad.xray_structures)
    print >> self.ta_obj.log, "Selected number xrs : ", len(best_xrs)

    for n in xrange(20):
      random_sel = random.sample(xrange(len(best_xrs)), total_number_xrs/4)
      print >> self.ta_obj.log, "Random selection (enriched) : ", random_sel
      cntr = 0
      fcalc_total = None
      fmask_total = None
      for x in xrange(len(best_xrs)):
        if x in random_sel:
          self.ta_obj.fmodel_ens.update_xray_structure(xray_structure = best_xrs[x],
                                                update_f_calc  = True,
                                                update_f_mask  = True,
                                                force_update_f_mask = True)
          if fcalc_total == None:
            fcalc_total = self.ta_obj.fmodel_ens.f_calc().data().deep_copy()
            fmasks = self.ta_obj.fmodel_ens.shell_f_masks()
            assert len(fmasks) == 1
            fmask_total = fmasks[0].data().deep_copy()
            cntr = 1
          else:
            fcalc_total += self.ta_obj.fmodel_ens.f_calc().data().deep_copy()
            fmasks = self.ta_obj.fmodel_ens.shell_f_masks()
            assert len(fmasks) == 1
            fmask_total += fmasks[0].data().deep_copy()
            cntr += 1
        if x == len(best_xrs)-1:
          self.ta_obj.fmodel_ens.update(
            f_calc = self.ta_obj.copy_ma.array(data = (fcalc_total / cntr)),
            f_mask = self.ta_obj.copy_ma.array(data = (fmask_total / cntr)) )
          self.ta_obj.fmodel_ens.update_solvent_and_scale(verbose = 0, out = self.ta_obj.log, optimize_mask = False)
          print >> self.ta_obj.log, "ERE: {0:8d} {1:8.3f} {2:8.3f} {3:8.3f} {4:8.3f} {5:8.3f} {6:8.3f}"\
            .format(cntr,
                    100*self.ta_obj.fmodel_ens.r_work(),
                    100*self.ta_obj.fmodel_ens.r_free(),
                    self.ta_obj.fmodel_ens.scale_k1(),
                    self.ta_obj.fmodel_ens.b_iso(),
                    self.ta_obj.fmodel_ens.k_sol(),
                    self.ta_obj.fmodel_ens.b_sol() )

    print  >> self.ta_obj.log, "\nBest ens : \n"
    x = total_number_xrs
    div_list = []
    while x > 1:
      div_list.append(int(x))
      x = x / 2

    sorted_xrs_rfree_tuple_list = sorted(xrs_rfree_tuple_list, key=lambda rfree:rfree[-1])

    for num in div_list:
      cntr = 0
      fcalc_total = None
      fmask_total = None
      while cntr <= num:
        if cntr < num:
          self.ta_obj.fmodel_ens.update_xray_structure(
            xray_structure = sorted_xrs_rfree_tuple_list[cntr][0],
            update_f_calc  = True,
            update_f_mask  = True,
            force_update_f_mask = True)
          print  >> self.ta_obj.log, self.ta_obj.fmodel_ens.r_free()
          if fcalc_total == None:
            fcalc_total = self.ta_obj.fmodel_ens.f_calc().data().deep_copy()
            fmasks = self.ta_obj.fmodel_ens.shell_f_masks()
            assert len( fmasks) ==1
            fmask_total = fmasks[0].data().deep_copy()
            cntr = 1
          else:
            fcalc_total += self.ta_obj.fmodel_ens.f_calc().data().deep_copy()
            fmasks = self.ta_obj.fmodel_ens.shell_f_masks()
            assert len(fmasks) == 1
            fmask_total += fmasks[0].data().deep_copy()
            cntr += 1
        else:
          self.ta_obj.fmodel_ens.update(
            f_calc = self.ta_obj.copy_ma.array(data = (fcalc_total / cntr)),
            f_mask = self.ta_obj.copy_ma.array(data = (fmask_total / cntr)) )
          self.ta_obj.fmodel_ens.update_solvent_and_scale(verbose = 0, out = self.ta_obj.log, optimize_mask = False)
          print >> self.ta_obj.log, "EBS: {0:8d} {1:8.3f} {2:8.3f} {3:8.3f} {4:8.3f} {5:8.3f} {6:8.3f}"\
            .format(cntr,
                    100*self.ta_obj.fmodel_ens.r_work(),
                    100*self.ta_obj.fmodel_ens.r_free(),
                    self.ta_obj.fmodel_ens.scale_k1(),
                    self.ta_obj.fmodel_ens.b_iso(),
                    self.ta_obj.fmodel_ens.k_sol(),
                    self.ta_obj.fmodel_ens.b_sol() )
          cntr += 1

    print >> self.ta_obj.log, "\nRandom distribution though trajectory: \n"
    for n in xrange(20):
      random_sel = random.sample(xrange(total_number_xrs), total_number_xrs/4)
      print >> self.ta_obj.log, "Random selection :", random_sel
      cntr = 0
      fcalc_total = None
      fmask_total = None
      for x in xrange(total_number_xrs):
        if x in random_sel:
          self.ta_obj.fmodel_ens.update_xray_structure(
            xray_structure = self.ta_obj.tad.xray_structures[x],
            update_f_calc  = True,
            update_f_mask  = True,
            force_update_f_mask = True)
          if fcalc_total == None:
            fcalc_total = self.ta_obj.fmodel_ens.f_calc().data().deep_copy()
            fmasks = self.ta_obj.fmodel_ens.shell_f_masks()
            assert len( fmasks) ==1
            fmask_total = fmasks[0].data().deep_copy()
            cntr = 1
          else:
            fcalc_total += self.ta_obj.fmodel_ens.f_calc().data().deep_copy()
            fmasks = self.ta_obj.fmodel_ens.shell_f_masks()
            assert len(fmasks) ==1
            fmask_total += fmasks[0].data().deep_copy()
            cntr += 1
        if x == total_number_xrs-1:
          self.ta_obj.fmodel_ens.update(
            f_calc = self.ta_obj.copy_ma.array(data = (fcalc_total / cntr)),
            f_mask = self.ta_obj.copy_ma.array(data = (fmask_total / cntr)) )
          self.ta_obj.fmodel_ens.update_solvent_and_scale(verbose = 0, out = self.ta_obj.log, optimize_mask = False)
          print >> self.ta_obj.log, "ERS: {0:8d} {1:8.3f} {2:8.3f} {3:8.3f} {4:8.3f} {5:8.3f} {6:8.3f}"\
            .format(cntr,
                    100*self.ta_obj.fmodel_ens.r_work(),
                    100*self.ta_obj.fmodel_ens.r_free(),
                    self.ta_obj.fmodel_ens.scale_k1(),
                    self.ta_obj.fmodel_ens.b_iso(),
                    self.ta_obj.fmodel_ens.k_sol(),
                    self.ta_obj.fmodel_ens.b_sol() )

    print >> self.ta_obj.log, "\nEqual distribution though trajectory: \n"
    div_int =  1
    print >> self.ta_obj.log, '\n {0:>12} {1:>8} {2:>8} {3:>8} {4:>8} {5:>8} {6:>8}'\
        .format('Num','Rwork','Rfree','k1','biso','ksol','bsol')
    while div_int < total_number_xrs:
      cntr = 0
      fcalc_total = None
      fmask_total = None
      for x in xrange(total_number_xrs):
        if x%int(div_int) == 0:
          self.ta_obj.fmodel_ens.update_xray_structure(
            xray_structure = self.ta_obj.tad.xray_structures[x],
            update_f_calc  = True,
            update_f_mask  = True,
            force_update_f_mask = True)
          if fcalc_total == None:
            fcalc_total = self.ta_obj.fmodel_ens.f_calc().data().deep_copy()
            fmasks = self.ta_obj.fmodel_ens.shell_f_masks()
            assert len(fmasks) ==1
            fmask_total = fmasks[0].data().deep_copy()
            cntr = 1
          else:
            fcalc_total += self.ta_obj.fmodel_ens.f_calc().data().deep_copy()
            fmasks = self.ta_obj.fmodel_ens.shell_f_masks()
            assert len(fmasks) == 1
            fmask_total += fmasks[0].data().deep_copy()
            cntr += 1
        if x == total_number_xrs-1:
          self.ta_obj.fmodel_ens.update(
            f_calc = self.ta_obj.copy_ma.array(data = (fcalc_total / cntr)),
            f_mask = self.ta_obj.copy_ma.array(data = (fmask_total / cntr)) )
          self.ta_obj.fmodel_ens.update_solvent_and_scale(verbose = 0, out = self.ta_obj.log, optimize_mask = False)
          print >> self.ta_obj.log, "EES: {0:8d} {1:8.3f} {2:8.3f} {3:8.3f} {4:8.3f} {5:8.3f} {6:8.3f}"\
            .format(cntr,
                    100*self.ta_obj.fmodel_ens.r_work(),
                    100*self.ta_obj.fmodel_ens.r_free(),
                    self.ta_obj.fmodel_ens.scale_k1(),
                    self.ta_obj.fmodel_ens.b_iso(),
                    self.ta_obj.fmodel_ens.k_sol(),
                    self.ta_obj.fmodel_ens.b_sol() )
      div_int = div_int * 2.0
    print >> self.ta_obj.log, "|"+"-"*77+"|\n"

  def chebyshev_smoothing(self,
                          x_data,
                          y_data,
                          x_miller,
                          n_chebyshev_terms = 10,
                          ):
    max_y = max(y_data)
    reparam_y = -flex.log(y_data)
    fit_lsq = chebyshev_lsq_fit.chebyshev_lsq_fit(
      n_terms = n_chebyshev_terms,
      x_obs   = x_data,
      y_obs   = reparam_y,
      w_obs   = None)
    min_h = flex.min(x_miller)
    max_h = flex.max(x_miller)
    cheb_pol = chebyshev_polynome(
      n_chebyshev_terms,
      min_h,
      max_h,
      fit_lsq.coefs)
    def reverse_reparam(values):
      return flex.exp(-values)
    y_fitted = reverse_reparam(cheb_pol.f(x_data))
    fitted_miller_array = reverse_reparam(cheb_pol.f(x_miller))
    fitted_miller_array = self.ta_obj.fmodel_running.f_obs.array(data=fitted_miller_array)
    return fitted_miller_array
