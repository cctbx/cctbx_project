# LIBTBX_SET_DISPATCHER_NAME cxi.brehm_diederichs
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
from __future__ import absolute_import, division, print_function

import iotbx.phil
from cctbx.crystal import symmetry
from libtbx.utils import Usage # implicit import
from libtbx.utils import multi_out
import sys
from xfel.cxi.util import is_odd_numbered # implicit import
from xfel.command_line.cxi_merge import master_phil
#-----------------------------------------------------------------------
from xfel.command_line.cxi_xmerge import xscaling_manager

extra_phil = """
asymmetric = 3
  .type = int(value_min=0, value_max=3)
show_plot = True
  .type = bool
save_plot = False
  .type = bool
"""

#-----------------------------------------------------------------------
def run(args):
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil + "\n" + extra_phil).show()
  work_params = phil.work.extract()
  log = open("%s_%s_merging.log" % (work_params.output.prefix,work_params.scaling.algorithm), "w")
  out = multi_out()
  out.register("log", log, atexit_send_to=None)
  out.register("stdout", sys.stdout)

  print("Target unit cell and space group:", file=out)
  print("  ", work_params.target_unit_cell, file=out)
  print("  ", work_params.target_space_group, file=out)

  miller_set = symmetry(
      unit_cell=work_params.target_unit_cell,
      space_group_info=work_params.target_space_group
    ).build_miller_set(
      anomalous_flag=not work_params.merge_anomalous,
      d_min=work_params.d_min)
  from xfel.merging.general_fcalc import random_structure
  i_model = random_structure(work_params)

# ---- Augment this code with any special procedures for x scaling
  scaler = xscaling_manager(
    miller_set=miller_set,
    i_model=i_model,
    params=work_params,
    log=out)
  scaler.read_all_mysql()
  print("finished reading the database")
  sg = miller_set.space_group()

  hkl_asu = scaler.observations["hkl_id"]
  imageno = scaler.observations["frame_id"]
  intensi = scaler.observations["i"]
  lookup = scaler.millers["merged_asu_hkl"]
  origH = scaler.observations["H"]
  origK = scaler.observations["K"]
  origL = scaler.observations["L"]

  from cctbx.array_family import flex

  print("# observations from the database",len(scaler.observations["hkl_id"]))
  hkl = flex.miller_index(flex.select(lookup,hkl_asu))
  from cctbx import miller

  hkl_list = miller_set.customized_copy(indices = hkl)

  ARRAY = miller.array(miller_set = hkl_list, data = intensi)
  LATTICES = miller.array(miller_set = hkl_list, data = imageno)

  from cctbx.merging.brehm_diederichs import run_multiprocess, run

  L = (ARRAY, LATTICES) # tuple(data,lattice_id)
  from libtbx import easy_pickle
  presort_file = work_params.output.prefix+"_intensities_presort.pickle"
  print("pickling these intensities to", presort_file)
  easy_pickle.dump(presort_file,L)

    ######  INPUTS #######
    #       data = miller array: ASU miller index + intensity (sigmas not implemented yet)
    #       lattice_id = flex double: assignment of each miller index to a lattice number
    ######################
  if work_params.nproc < 5:
    print("Sorting the lattices with 1 processor")
    result = run(L, asymmetric=work_params.asymmetric, nproc=1, verbose=True,
                 show_plot=work_params.show_plot, save_plot=work_params.save_plot)
  else:
    print("Sorting the lattices with %d processors"%work_params.nproc)
    result = run_multiprocess(L, asymmetric=work_params.asymmetric, nproc=work_params.nproc, verbose=False,
                              show_plot=work_params.show_plot, save_plot=work_params.save_plot)
  for key in result.keys():
    print(key,len(result[key]))

  # 2) pickle the postsort (reindexed) ARRAY, LATTICES XXX not done yet; not clear if needed

  reverse_lookup = {}
  frame_id_list = list(scaler.frames_mysql["frame_id"])
  for key in result.keys():
    for frame in result[key]:
      frame_idx = frame_id_list.index(frame)
      reverse_lookup[scaler.frames_mysql["unique_file_name"][frame_idx]] = key

  lookup_file = work_params.output.prefix+"_lookup.pickle"
  reverse_lookup_file = work_params.output.prefix+"_reverse_lookup.pickle"
  easy_pickle.dump(lookup_file, result)
  easy_pickle.dump(reverse_lookup_file, reverse_lookup)


if (__name__ == "__main__"):
  import sys
  result = run(args = sys.argv[1:])
