'''
Author      : Uervirojnangkoorn, M.
Created     : 12/1/2014
Description : Optimizer main module.
'''
from __future__ import absolute_import, division, print_function
import numpy as np
import math, random
from cctbx.array_family import flex
from datetime import datetime
from six.moves import range

class sisa_optimizer(object):
  '''
  Setup project environment, sort reflections and return selected indices,
  then perform microcycle optimization (macrocyle is done in the command file).
  '''
  def __init__(self):
    '''
    Constructor
    '''
    self.phi_for_hl=flex.double(range(0,360,1))

  def calc_pdf_cdf_from_hl(self, hl_given):
    cdf_from_pdf_set=[]
    phi_for_hl_radian = self.phi_for_hl * math.pi / 180
    for i_refl in range(len(hl_given)):
      pdf_from_hl=[0]*len(phi_for_hl_radian)
      cdf_from_pdf=[0]*len(phi_for_hl_radian)

      for i_phi in range(len(phi_for_hl_radian)):
        phi_given_radian=phi_for_hl_radian[i_phi]
        pdf_from_hl[i_phi]=math.exp((hl_given[i_refl][0]*math.cos(phi_given_radian))+
                                    (hl_given[i_refl][1]*math.sin(phi_given_radian))+
                                    (hl_given[i_refl][2]*math.cos(2*phi_given_radian))+
                                    (hl_given[i_refl][3]*math.sin(2*phi_given_radian)));
        cdf_from_pdf[i_phi]=sum(pdf_from_hl[0:i_phi+1]);

      #scale the cdf
      max_cdf_from_pdf=max(cdf_from_pdf);
      for j_phi in range(len(cdf_from_pdf)):
        cdf_from_pdf[j_phi]/=max_cdf_from_pdf;

      cdf_from_pdf_set.append(cdf_from_pdf);

    return cdf_from_pdf_set

  def setup_map_coeff(self, miller_arrays, indices_selected, phi_selected, fom_selected):

    flex_fp = miller_arrays[0].data()
    flex_sigmas = miller_arrays[0].sigmas()
    flex_phi_new = miller_arrays[1].data()
    flex_fom_new = miller_arrays[2].data()
    for i in range(len(indices_selected)):
      flex_fom_new[indices_selected[i]] = fom_selected[i]
      flex_phi_new[indices_selected[i]] = phi_selected[i]

    flex_fp_w_new = flex_fp * flex_fom_new
    miller_array_fp_w = miller_arrays[0].customized_copy(data=flex_fp_w_new)

    map_coeff = miller_array_fp_w.phase_transfer(phase_source=flex_phi_new, deg=True)

    return map_coeff

  def calcskew(self, map_coeff, iparams):
    real_map = map_coeff.fft_map().real_map()

    ed_limit = iparams.fit_params.ed_sigma_thres*real_map.sample_standard_deviation()

    ed_limit_up=ed_limit
    ed_limit_dn=ed_limit*-1

    #truncate the electron density beyond sigma limit
    real_map.set_selected(real_map>ed_limit_up,ed_limit_up)
    real_map.set_selected(real_map<ed_limit_dn,ed_limit_dn)

    skew = flex.mean(flex.pow(real_map,3))/pow(flex.mean_sq(real_map),3/2)

    return skew

  def pickbestidv(self, list_idv, list_fit, lo_sigma_fit, up_sigma_fit):
    list_phis_good = []
    mean_fit = sum(list_fit)/len(list_fit)
    std_fit = np.std(list_fit)
    fit_thres_lo = mean_fit + (lo_sigma_fit * std_fit)
    fit_thres_hi = mean_fit + (up_sigma_fit * std_fit)

    for i_idv in range(len(list_idv)):
      if list_fit[i_idv] > fit_thres_lo and list_fit[i_idv] <= fit_thres_hi:
        list_phis_good.append(list_idv[i_idv])

    #find centroid phase for phis_good
    from mmtbx.sisa.optimize.mod_util import util_handler
    uth = util_handler()
    flex_phi_bar, dummy = uth.calcphibar(list_phis_good)

    return flex_phi_bar, mean_fit, std_fit, len(list_phis_good)

  def calc_stats(self, miller_arrays, indices_selected, phi_selected, fom_selected, iparams):
    map_coeff = self.setup_map_coeff(miller_arrays, indices_selected,
                                      phi_selected, fom_selected)
    skew = self.calcskew(map_coeff, iparams)

    fp_selected = flex.double([miller_arrays[0].data()[inds] for inds in indices_selected])
    phic_selected = flex.double([miller_arrays[4].data()[inds] for inds in indices_selected])

    mpe_phi = 0
    mapcc_phi = 0
    if iparams.hklrefin is not None and \
      np.sum(phic_selected) > 0:
      from mmtbx.sisa.optimize.mod_util import util_handler
      uth = util_handler()
      mapcc_phi, mpe_phi = uth.calcphicc(fp_selected,
                                           [1]*len(fom_selected),
                                           fom_selected,
                                           phic_selected,
                                           phi_selected,
                                           True)
    return skew, mapcc_phi, mpe_phi

  def run_optimize(self, micro_cycle_no, stack_no, miller_arrays, indices_selected, cdf_set, iparams):

    #start ga
    from mmtbx.sisa.optimize.mod_ga import ga_handler
    gah = ga_handler()

    from mmtbx.sisa.optimize.mod_util import util_handler
    uth = util_handler()

    list_phis_intcycle = []
    list_fit_intcycle = []

    #initial fist population and calculate their fitness
    ga_idv_length = len(indices_selected)
    cur_pop=gah.initialse(iparams.ga_params.pop_size, ga_idv_length, cdf_set, self.phi_for_hl)
    cur_fit=[0]*len(cur_pop)

    #setup new fp, phic (if given), fom, and new fom
    fp_selected = flex.double([miller_arrays[0].data()[inds] for inds in indices_selected])
    phib_selected = flex.double([miller_arrays[1].data()[inds] for inds in indices_selected])
    fom_selected = flex.double([miller_arrays[2].data()[inds] for inds in indices_selected])

    fom_selected_new = fom_selected + 0.2

    #calculate initial mapcc and mpe
    skew, mapcc, mpe = self.calc_stats(miller_arrays, indices_selected,
                                        phib_selected, fom_selected, iparams)

    txt_prn_out = "Starting stack %2.0f: microcycle %2.0f (initial skew=%6.2f, mapcc=%6.2f, mpe=%6.2f)\n"%(\
      stack_no+1, micro_cycle_no+1, skew, mapcc, mpe*180/math.pi)
    print(txt_prn_out)
    txt_pop_hist_out=""
    for i_idv in range(len(cur_pop)):
      map_coeff = self.setup_map_coeff(miller_arrays, indices_selected,
            flex.double(cur_pop[i_idv]), fom_selected_new)
      cur_fit[i_idv] = self.calcskew(map_coeff, iparams)

      #collect the first population
      list_phis_intcycle.append(cur_pop[i_idv])
      list_fit_intcycle.append(cur_fit[i_idv])

    txt_prn_tmp = 'gen'.center(5)+'<skew>'.center(7)+'std_skew'.center(8)+'n_accidv'.center(10)+ \
      'skew'.center(6)+'mapcc'.center(7)+'mpe'.center(5)+'mapccp'.center(7)+'mpep'.center(6)+'time_spent (min)'.center(16)+'\n'
    print(txt_prn_tmp)
    txt_prn_out += txt_prn_tmp
    #Set up population map
    #[[7,1,5],
    #[2,3,0],
    #[4,6,9]]
    map_width=int(math.sqrt(iparams.ga_params.pop_size))
    map_1D=random.sample(range(iparams.ga_params.pop_size), iparams.ga_params.pop_size)

    map_2D=[]
    for i_width in range(map_width):
      map_2D.append(map_1D[int(map_width*i_width):int(map_width*(i_width+1))])

    map_visit_order=random.sample(range(iparams.ga_params.pop_size), iparams.ga_params.pop_size)

    #start a generation
    conv_gen = iparams.ga_params.max_gen
    for i_gen in range(iparams.ga_params.max_gen):
      time_gen_start=datetime.now()

      #calculate crossover rate for this generation
      ga_ratio_cross=iparams.ga_params.crossover_start_rate + \
        (math.pow(i_gen/iparams.ga_params.max_gen, iparams.ga_params.crossover_slope) * \
        (iparams.ga_params.crossover_end_rate-iparams.ga_params.crossover_start_rate))

      num_point_cross=int(round(ga_ratio_cross * ga_idv_length))
      for i_idv in range(len(map_visit_order)):
        mom_x=int(math.fmod(map_visit_order[i_idv], map_width))
        mom_y=int(math.floor(map_visit_order[i_idv]/map_width))
        mom_id=map_2D[mom_x][mom_y]
        mom=cur_pop[mom_id]
        mom_fit=cur_fit[mom_id]


        #draw an xmap around mom and grab the idv_id stored in map_2D
        mate_candidate_id=[]
        xmap_width = (iparams.ga_params.xmap_radius * 2) + 1
        for i_xmap_y in range(xmap_width):
          for i_xmap_x in range(xmap_width):
            i_xmap_tmp_x=mom_x+i_xmap_x-iparams.ga_params.xmap_radius
            i_xmap_tmp_y=mom_y+i_xmap_y-iparams.ga_params.xmap_radius
            if i_xmap_tmp_x >= xmap_width:
              i_xmap_tmp_x-= xmap_width
            if i_xmap_tmp_y >= xmap_width:
              i_xmap_tmp_y-= xmap_width
            if ~(map_2D[i_xmap_tmp_x][i_xmap_tmp_y]==mom_id):
              mate_candidate_id.append(map_2D[i_xmap_tmp_x][i_xmap_tmp_y])

        mate_candiate_id_random_order=random.sample(range(len(mate_candidate_id)), iparams.ga_params.num_sel_mate)

        tmp_mate_fit_set=[]
        for i_mate in range(iparams.ga_params.num_sel_mate):
          tmp_mate_fit_set.append(
                [mate_candidate_id[mate_candiate_id_random_order[i_mate]],
                cur_fit[mate_candidate_id[mate_candiate_id_random_order[i_mate]]]])

        tmp_mate_fit_sort=sorted(tmp_mate_fit_set, key=lambda fit: fit[1],reverse=True)

        dad_id=tmp_mate_fit_sort[0][0]
        dad=cur_pop[dad_id]

        #perform ga operator - perform under probability
        #otherwise keeps the mom
        if random.random() < iparams.ga_params.prob_of_cross:
          child1,child2,cross_template = gah.crossover(mom, dad, ga_ratio_cross)
          child1 = gah.mutation(
                  child1,
                  iparams.ga_params.prob_of_mut,
                  iparams.ga_params.num_point_mut,
                  cdf_set,
                  self.phi_for_hl)
          child2 = gah.mutation(
                  child2,
                  iparams.ga_params.prob_of_mut,
                  iparams.ga_params.num_point_mut,
                  cdf_set,
                  self.phi_for_hl)

          #recalculate fitness
          map_coeff = self.setup_map_coeff(miller_arrays, indices_selected,
            flex.double(child1), fom_selected_new)
          child1_fit = self.calcskew(map_coeff, iparams)

          map_coeff = self.setup_map_coeff(miller_arrays, indices_selected,
            flex.double(child2), fom_selected_new)
          child2_fit = self.calcskew(map_coeff, iparams)

          if child1_fit >= child2_fit:
            child_sel=child1[:]
            child_fit_sel=child1_fit
          else:
            child_sel=child2[:]
            child_fit_sel=child2_fit

          if child_fit_sel > cur_fit[mom_id]:
            cur_pop[mom_id]=child_sel[:]
            cur_fit[mom_id]=child_fit_sel


      ''' collect some stats '''
      #1. record some stats
      time_gen_end=datetime.now()
      time_gen_spent=time_gen_end-time_gen_start


      #2. mapcc and mpe among population
      n_idv_pick=int(round(0.05*iparams.ga_params.pop_size))
      mpe_idv_pick=[0]*n_idv_pick
      mapcc_idv_pick=[0]*n_idv_pick
      id_idv_pick=random.sample(range(iparams.ga_params.pop_size), n_idv_pick)
      for i in range(n_idv_pick):
        i_idv_pick = id_idv_pick[i]
        mpe_to_others = [0] * iparams.ga_params.pop_size
        mapcc_to_others = [0] * iparams.ga_params.pop_size
        for j in range(iparams.ga_params.pop_size):
          mapcc_to_others[j], mpe_to_others[j] = uth.calcphicc(
                fp_selected,
                fom_selected_new,
                fom_selected_new,
                cur_pop[i_idv_pick],
                cur_pop[j],
                True)

        mpe_idv_pick[i]=sum(mpe_to_others)/iparams.ga_params.pop_size
        mapcc_idv_pick[i]=sum(mapcc_to_others)/iparams.ga_params.pop_size

      mpe_avg_gen = sum(mpe_idv_pick)/len(mpe_idv_pick)
      mapcc_avg_gen = sum(mapcc_idv_pick)/len(mapcc_idv_pick)

      #3. collect the population in generation
      for i_idv in range(len(cur_pop)):
        list_phis_intcycle.append(cur_pop[i_idv])
        list_fit_intcycle.append(cur_fit[i_idv])


      #4. calculate averaged phis among available phis now
      flex_phis_intcycle, \
      mean_fit_intcycle, \
      std_fit_intcycle, \
      num_good_idv_intcycle = self.pickbestidv(
                                list_phis_intcycle,
                                list_fit_intcycle,
                                iparams.ga_params.skew_sigma_sel_lo,
                                iparams.ga_params.skew_sigma_sel_hi)

      flex_phis_fit_intcycle, mapcc_phis_intcycle, mpe_phis_intcycle = self.calc_stats(\
              miller_arrays, indices_selected, flex_phis_intcycle, fom_selected_new, iparams)

      txt_prn_tmp = '%2.0f %7.2f %8.2f %6.0f %8.2f %5.2f %6.2f %6.2f %6.2f %10.2f\n'%(i_gen+1, \
        mean_fit_intcycle, std_fit_intcycle, num_good_idv_intcycle, flex_phis_fit_intcycle, \
        mapcc_phis_intcycle, mpe_phis_intcycle*180/math.pi, \
        mapcc_avg_gen, mpe_avg_gen*180/math.pi, time_gen_spent.seconds/60)
      print(txt_prn_tmp)
      txt_prn_out += txt_prn_tmp
      #check termination
      if mapcc_avg_gen >= 0.9:
        conv_gen = i_gen
        break




    txt_prn_out += 'Population converged at generation '+str(conv_gen+1)+'\n'

    #complete one internal cycle
    if iparams.flag_log_verbose_on:
      #for phis in list_phis_intcycle:
      #  print phis
      """
      file_name_pop_hist_out = iparams.project_name+'/'+iparams.run_name+"/log_pop_hist_step_"+str(stack_no+1)+"_intcycle_"+str(micro_cycle_no+1)+".log"
      output = open(file_name_pop_hist_out, 'w')
      output.write(list(list_phis_intcycle))
      output.close()
      """


    return flex_phis_intcycle, fom_selected_new, flex_phis_fit_intcycle, txt_prn_out
