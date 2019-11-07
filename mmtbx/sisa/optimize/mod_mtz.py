'''
Author      : Uervirojnangkoorn, M.
Created     : 12/1/2014
Description : Handling mtz file.
'''
from __future__ import absolute_import, division, print_function
import sys
import numpy as np
from iotbx import reflection_file_reader
from cctbx.array_family import flex
from cctbx import miller
from cctbx import crystal
from six.moves import range

class mtz_handler(object):
  '''
  Handling mtz files.
  '''
  def __init__(self):
    '''
    Constructor
    '''

  def format_miller_arrays(self, iparams):
    '''
    Read in mtz file and format to miller_arrays_out object with
    index[0] --> FP, SIGFP
    index[1] --> PHIB
    index[2] --> FOM
    index[3] --> HLA, HLB, HLC, HLD
    index[4] --> optional PHIC
    '''
    #readin reflection file
    reflection_file = reflection_file_reader.any_reflection_file(iparams.data)

    file_content=reflection_file.file_content()
    column_labels=file_content.column_labels()
    col_name=iparams.column_names.split(',')

    miller_arrays=reflection_file.as_miller_arrays()
    flex_centric_flags = miller_arrays[0].centric_flags().data()
    crystal_symmetry = crystal.symmetry(
        unit_cell=miller_arrays[0].unit_cell(), space_group=miller_arrays[0].space_group())

    #grab all required columns
    flag_fp_found = 0
    flag_phib_found = 0
    flag_fom_found = 0
    flag_hl_found = 0
    ind_miller_array_fp = 0
    ind_miller_array_phib = 0
    ind_miller_array_fom = 0
    ind_miller_array_hl = 0
    for i in range(len(miller_arrays)):
      label_string = miller_arrays[i].info().label_string()
      labels=label_string.split(',')
      #only look at first index string
      if labels[0]==col_name[0]:
        #grab FP, SIGFP
        flex_fp_all=miller_arrays[i].data()
        flex_sigmas_all=miller_arrays[i].sigmas()
        flag_fp_found=1
        ind_miller_array_fp = i
      elif labels[0]==col_name[2]:
        #grab PHIB
        flex_phib_all=miller_arrays[i].data()
        flag_phib_found=1
        ind_miller_array_phib = i
      elif labels[0]==col_name[3]:
        #grab FOM
        flex_fom_all=miller_arrays[i].data()
        flag_fom_found=1
        ind_miller_array_fom = i
      elif labels[0]==col_name[4]:
        #grab HLA,HLB,HLC,HLD
        flex_hl_all=miller_arrays[i].data()
        flag_hl_found=1
        ind_miller_array_hl = i

    if flag_hl_found==1 and flag_phib_found == 0:
      #calculate PHIB and FOM from HL
      miller_array_phi_fom = miller_arrays[ind_miller_array_hl].phase_integrals()
      flex_phib_all = miller_array_phi_fom.phases(deg=True).data()
      flex_fom_all = miller_array_phi_fom.amplitudes().data()
      flag_phib_found = 1
      flag_fom_found = 1

    if flag_fp_found==0 or flag_phib_found==0 or flag_fom_found==0 or flag_hl_found==0:
      print("couldn't find all required columns")
      sys.exit()

    miller_indices_sel = miller_arrays[ind_miller_array_fp].indices()
    print('No. reflections for read-in miller arrays - indices:%6.0f fp:%6.0f phib:%6.0f fom:%6.0f HL:%6.0f)'%( \
          len(miller_indices_sel), len(flex_fp_all), len(flex_phib_all), len(flex_fom_all), len(flex_hl_all)))

    miller_indices = flex.miller_index()
    flex_fp = flex.double()
    flex_sigmas = flex.double()
    flex_phib = flex.double()
    flex_fom = flex.double()
    flex_hl = flex.hendrickson_lattman()
    #format all miller arrays to the same length
    for miller_index in miller_indices_sel:
      fp_cn, phib_cn, fom_cn, hl_cn = (0,0,0,0)

      matches = miller.match_multi_indices(
                    miller_indices_unique=flex.miller_index([miller_index]),
                    miller_indices=miller_arrays[ind_miller_array_fp].indices())
      if len(matches.pairs()) > 0:
        fp_cn = 1
        fp = flex_fp_all[matches.pairs()[0][1]]
        sigmas = flex_sigmas_all[matches.pairs()[0][1]]

      matches = miller.match_multi_indices(
                    miller_indices_unique=flex.miller_index([miller_index]),
                    miller_indices=miller_arrays[ind_miller_array_phib].indices())
      if len(matches.pairs()) > 0:
        phib_cn = 1
        phib = flex_phib_all[matches.pairs()[0][1]]

      matches = miller.match_multi_indices(
                    miller_indices_unique=flex.miller_index([miller_index]),
                    miller_indices=miller_arrays[ind_miller_array_fom].indices())
      if len(matches.pairs()) > 0:
        fom_cn = 1
        fom = flex_fom_all[matches.pairs()[0][1]]

      matches = miller.match_multi_indices(
                    miller_indices_unique=flex.miller_index([miller_index]),
                    miller_indices=miller_arrays[ind_miller_array_hl].indices())
      if len(matches.pairs()) > 0:
        hl_cn = 1
        hl = flex_hl_all[matches.pairs()[0][1]]

      if (fp_cn + phib_cn + fom_cn + hl_cn) == 4:
        miller_indices.append(miller_index)
        flex_fp.append(fp)
        flex_sigmas.append(sigmas)
        flex_phib.append(phib)
        flex_fom.append(fom)
        flex_hl.append(hl)



    print('No. reflections after format - indices:%6.0f fp:%6.0f phib:%6.0f fom:%6.0f HL:%6.0f)'%( \
          len(miller_indices), len(flex_fp), len(flex_phib), len(flex_fom), len(flex_hl)))

    flex_hla = flex.double()
    flex_hlb = flex.double()
    flex_hlc = flex.double()
    flex_hld = flex.double()
    for i in range(len(flex_hl)):
      data_hl_row=flex_hl[i]
      flex_hla.append(data_hl_row[0])
      flex_hlb.append(data_hl_row[1])
      flex_hlc.append(data_hl_row[2])
      flex_hld.append(data_hl_row[3])
    '''
    Read benchmark MTZ (PHICalc) for MPE calculation
    '''
    flex_phic = flex.double([0]*len(flex_fp))
    if iparams.hklrefin is not None:
      reflection_file = reflection_file_reader.any_reflection_file(iparams.hklrefin)
      miller_arrays_bench=reflection_file.as_miller_arrays()
      flex_phic_raw = None
      for i in range(len(miller_arrays_bench)):
        label_string = miller_arrays_bench[i].info().label_string()
        labels=label_string.split(',')
        #only look at first index string
        if labels[0] == iparams.column_phic:
          #grab PHIC
          if miller_arrays_bench[i].is_complex_array():
            flex_phic_raw = miller_arrays_bench[i].phases(deg=True).data()
          else:
            flex_phic_raw = miller_arrays_bench[i].data()
          miller_indices_phic = miller_arrays_bench[i].indices()

      if flex_phic is not None:
        matches = miller.match_multi_indices(
                  miller_indices_unique=miller_indices,
                  miller_indices=miller_indices_phic)

        flex_phic = flex.double([flex_phic_raw[pair[1]] for pair in matches.pairs()])

    #format miller_arrays_out
    miller_set=miller.set(
            crystal_symmetry=crystal_symmetry,
            indices=miller_indices,
            anomalous_flag=False)
    miller_array_out = miller_set.array(
            data=flex_fp,
            sigmas=flex_sigmas).set_observation_type_xray_amplitude()

    #check if Wilson B-factor is applied
    flex_fp_for_sort = flex_fp[:]
    if iparams.flag_apply_b_factor:
      try:
        #get wilson_plot
        from mmtbx.scaling import xtriage
        from libtbx.utils import null_out
        xtriage_args = [
          iparams.data,
          "",
          "",
          "log=tst_xtriage_1.log"
        ]
        result = xtriage.run(args=xtriage_args, out=null_out())
        ws = result.wilson_scaling

        print('Wilson K=%6.2f B=%6.2f'%(ws.iso_p_scale, ws.iso_b_wilson))
        sin_theta_over_lambda_sq = miller_array_out.two_theta(wavelength=iparams.wavelength) \
                                    .sin_theta_over_lambda_sq().data()
        wilson_expect = flex.exp(-2 * ws.iso_b_wilson * sin_theta_over_lambda_sq)
        flex_fp_for_sort = wilson_expect * flex_fp
      except Exception:
        print('Error calculating Wilson scale factors. Continue without applying B-factor.')


    flex_d_spacings = miller_array_out.d_spacings().data()

    mtz_dataset = miller_array_out.as_mtz_dataset(column_root_label="FP")

    for data,lbl,typ in [(flex_phib, "PHIB", "P"),
        (flex_fom, "FOMB", "W"),
        (flex_hla,"HLA","A"),
        (flex_hlb,"HLB","A"),
        (flex_hlc,"HLC","A"),
        (flex_hld,"HLD","A"),
        (flex_phic,"PHIC","P")]:
        mtz_dataset.add_miller_array(miller_array_out.array(data=data),
            column_root_label=lbl,
            column_types=typ)

    miller_arrays_out = mtz_dataset.mtz_object().as_miller_arrays()

    '''
    getting sorted indices for the selected reflections in input mtz file
    list_fp_sort_index: stores indices of sorted FP in descending order
    '''
    import operator
    fp_sort_index= [i for (i,j) in sorted(enumerate(flex_fp_for_sort), key=operator.itemgetter(1))]
    fp_sort_index.reverse()

    """
    for i in range(100):
      print miller_indices[fp_sort_index[i]], flex_d_spacings[fp_sort_index[i]], flex_fp[fp_sort_index[i]], flex_sigmas[fp_sort_index[i]], wilson_expect[fp_sort_index[i]]

    exit()
    """

    #calculate sum of fp^2 from percent_f_squared
    flex_fp_squared = flex_fp ** 2
    f_squared_per_stack = (iparams.percent_f_squared * np.sum(flex_fp_squared))/100
    fp_sort_index_stacks = []
    sum_fp_now, i_start = (0,0)
    for i in range(len(fp_sort_index)):
      i_sel = fp_sort_index[i_start:i+1]
      sum_fp_now = np.sum([flex_fp_squared[ii_sel] for ii_sel in i_sel])
      if sum_fp_now >= f_squared_per_stack:
        fp_sort_index_stacks.append(fp_sort_index[i_start:i+1])
        i_start = i+1
        if len(fp_sort_index_stacks) == iparams.n_stacks:
          break

    txt_out = 'stack_no sum(f_squared) %total  n_refl\n'
    for i in range(len(fp_sort_index_stacks)):
      sum_fp = np.sum([flex_fp_squared[ii_sel] for ii_sel in fp_sort_index_stacks[i]])
      txt_out += '%6.0f %14.2f %8.2f %6.0f\n'%(i+1, sum_fp, \
        (sum_fp/np.sum(flex_fp_squared))*100, len(fp_sort_index_stacks[i]))

    return miller_arrays_out, fp_sort_index_stacks, txt_out

  def write_mtz(self, miller_arrays, file_name_out):
    crystal_symmetry = crystal.symmetry(
        unit_cell=miller_arrays[0].unit_cell(), space_group=miller_arrays[0].space_group())

    #get data from miller_arrays
    flex_fp = miller_arrays[0].data()
    flex_sigmas = miller_arrays[0].sigmas()
    flex_phib = miller_arrays[1].data()
    flex_fom = miller_arrays[2].data()
    flex_hl = miller_arrays[3].data()

    #format hla,b,c,d
    flex_hla = flex.double()
    flex_hlb = flex.double()
    flex_hlc = flex.double()
    flex_hld = flex.double()
    for i in range(len(flex_hl)):
      data_hl_row=flex_hl[i]
      flex_hla.append(data_hl_row[0])
      flex_hlb.append(data_hl_row[1])
      flex_hlc.append(data_hl_row[2])
      flex_hld.append(data_hl_row[3])

    #format miller_arrays_out
    miller_set=miller.set(
            crystal_symmetry=crystal_symmetry,
            indices=miller_arrays[0].indices(),
            anomalous_flag=False)
    miller_array_out = miller_set.array(
            data=flex_fp,
            sigmas=flex_sigmas).set_observation_type_xray_amplitude()

    mtz_dataset = miller_array_out.as_mtz_dataset(column_root_label="FP")

    for data,lbl,typ in [(flex_phib, "PHIB", "P"),
        (flex_fom, "FOMB", "W"),
        (flex_hla,"HLA","A"),
        (flex_hlb,"HLB","A"),
        (flex_hlc,"HLC","A"),
        (flex_hld,"HLD","A")]:
        mtz_dataset.add_miller_array(miller_array_out.array(data=data),
            column_root_label=lbl,
            column_types=typ)

    mtz_dataset.mtz_object().write(file_name=file_name_out)
