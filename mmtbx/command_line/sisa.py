# LIBTBX_SET_DISPATCHER_NAME phenix.sisa
'''
Author      : Uervirojnangkoorn, M.
Created     : 12/1/2014
Description : Commands linked to sisa libraries.
'''
from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx.easy_mp import pool_map
import math
import sys,os
from six.moves import range

def read_input(args):
  from mmtbx.sisa.optimize.mod_input import process_input
  iparams, txt_out_input = process_input(args)
  return iparams, txt_out_input

def sisa_optimize_mproc(micro_cycle_no, stack_no, miller_arrays, indices_selected, cdf_set, iparams):
  from mmtbx.sisa.optimize.mod_optimize import sisa_optimizer
  somer = sisa_optimizer()
  result = somer.run_optimize(micro_cycle_no, stack_no, miller_arrays, indices_selected, cdf_set, iparams)
  return result

def update_miller_arrays(miller_arrays, indices_selected, phis_selected, foms_selected):
  flex_phib = miller_arrays[1].data()
  flex_fomb = miller_arrays[2].data()
  for i in range(len(indices_selected)):
    flex_phib[indices_selected[i]] = phis_selected[i]
    flex_fomb[indices_selected[i]] = foms_selected[i]

  miller_arrays_out = []
  miller_arrays_out.append(miller_arrays[0])
  miller_arrays_out.append(miller_arrays[1].customized_copy(data=flex_phib))
  miller_arrays_out.append(miller_arrays[2].customized_copy(data=flex_fomb))
  miller_arrays_out.append(miller_arrays[3])
  miller_arrays_out.append(miller_arrays[4])

  return miller_arrays_out

if __name__=="__main__":
  txt_out = ''

  iparams, txt_out_input = read_input(sys.argv[:1])
  txt_out += txt_out_input

  from mmtbx.sisa.optimize.mod_mtz import mtz_handler
  mtzh = mtz_handler()
  miller_arrays, fp_sort_index_stacks, txt_out_format = mtzh.format_miller_arrays(iparams)
  print(txt_out_format)
  txt_out += txt_out_format

  for i in range(iparams.n_macro_cycles):
    txt_out += 'Macrocycle no. %4.0f\n'%(i+1)
    print('Macrocycle no. %4.0f\n'%(i+1))
    for j in range(len(fp_sort_index_stacks)):
      #select the index group
      i_sel = fp_sort_index_stacks[j]

      #generate cdf_set for selected reflections
      from mmtbx.sisa.optimize.mod_optimize import sisa_optimizer
      somer = sisa_optimizer()
      hl_selected = flex.hendrickson_lattman([miller_arrays[3].data()[ii_sel] for ii_sel in i_sel])
      cdf_set = somer.calc_pdf_cdf_from_hl(hl_selected)

      def sisa_optimize_mproc_wrapper(arg):
        return sisa_optimize_mproc(arg, j, miller_arrays, i_sel, cdf_set, iparams)

      sisa_optimize_results = pool_map(
              args=range(iparams.n_micro_cycles),
              func=sisa_optimize_mproc_wrapper,
              processes=iparams.n_processors)

      list_phis = []
      foms_sum = None
      list_skews = []
      for result in sisa_optimize_results:
        if result is not None:
          phis, foms, skews, txt_out_optim = result
          list_phis.append(phis)
          list_skews.append(skews)
          if foms_sum is None:
            foms_sum = foms[:]
          else:
            foms_sum += foms[:]
          print(txt_out_optim)
          txt_out += txt_out_optim

      #calculate centroid of list_phis
      phis_averaged, skews_averaged, skews_std, n_phis_selected = somer.pickbestidv(
                                list_phis,
                                list_skews,
                                -99, 99)
      foms_averaged = foms_sum/len(list_phis)

      skew_phis_averaged, mapcc_phis_averaged, mpe_phis_averaged = somer.calc_stats(\
              miller_arrays, i_sel, phis_averaged, foms_averaged, iparams)

      txt_out_tmp = 'Averaged phis skew=%6.2f mapcc=%6.2f mpe=%6.2f'%( \
        skew_phis_averaged, mapcc_phis_averaged, mpe_phis_averaged*180/math.pi)
      print(txt_out_tmp)
      txt_out += txt_out_tmp


      #update miller_arrays
      miller_arrays = update_miller_arrays(miller_arrays,
                                           i_sel,
                                           phis_averaged, foms_averaged)

      #output mtz for optimized stack n
      file_name_out = iparams.project_name + '/' + iparams.run_name + '/' + \
        'sisa_cycle_'+str(i+1)+'_stack_'+str(j+1)+'.mtz'
      mtzh.write_mtz(miller_arrays, file_name_out)

  f = open(iparams.project_name + '/' + iparams.run_name +'/log.txt', 'w')
  f.write(txt_out)
  f.close()

  print('Sisa done.')

  if iparams.autodm:
    print('Proceed with automatic density modification...(your density-modified map will be AutoBuild_run_n_/overall_best_denmod_map_coeffs.mtz.')
    cmd='phenix.autobuild data=' + file_name_out + ' seq_file=' + str(iparams.seq_file) + \
        ' maps_only=True n_cycle_build_max=1 n_cycle_rebuild_max=0' + \
        ' input_ha_file=' + str(iparams.ha_file) + ' model=' + str(iparams.model_file)
    print('Running: '+cmd)
    os.system(cmd)
