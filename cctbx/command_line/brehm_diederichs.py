#!/usr/bin/env python
#
# brehm_diederichs.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import os
import iotbx.phil
from libtbx.phil import command_line
from cctbx import sgtbx
from cctbx.array_family import flex
from iotbx.reflection_file_reader import any_reflection_file


master_phil_scope = iotbx.phil.parse("""
asymmetric = 3
  .type = int(value_min=0, value_max=3)
nproc = 1
  .type = int(value_min=1)
show_plot = True
  .type = bool
save_plot = False
  .type = bool
suffix = "_reindexed"
  .type = str
""")


def run(args):
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil, files = cmd_line.process_and_fetch(
    args=args, custom_processor="collect_remaining")
  working_phil.show()
  params = working_phil.extract()

  miller_array_all = None
  lattice_ids = None
  space_group = None
  file_name_dict = {}
  lattice_id = -1

  for file_name in files:
    lattice_id += 1
    #print "lattice_id: %i" %(lattice_id)
    reader = any_reflection_file(file_name)

    as_miller_arrays = reader.as_miller_arrays(merge_equivalents=False)
    #for ma in as_miller_arrays: print ma.info().labels
    intensities = [ma for ma in as_miller_arrays
                   if ma.info().labels == ['I', 'SIGI']][0]
    intensities = intensities.customized_copy(anomalous_flag=True).set_info(
      intensities.info())
    intensities.set_observation_type_xray_intensity()
    #intensities.crystal_symmetry().show_summary()
    #print intensities.info().labels
    if space_group is None:
      space_group = intensities.space_group()
    else:
      assert intensities.space_group() == space_group
    assert reader.file_type() == 'ccp4_mtz'

    file_name_dict[lattice_id] = file_name

    ids = intensities.customized_copy(
      data=flex.double(intensities.size(), lattice_id), sigmas=None)
    assert ids.size() == intensities.size()
    if miller_array_all is None:
      miller_array_all = intensities
      lattice_ids = ids
    else:
      miller_array_all = miller_array_all.customized_copy(
        indices=miller_array_all.indices().concatenate(intensities.indices()),
        data=miller_array_all.data().concatenate(intensities.data()),
        sigmas=miller_array_all.sigmas().concatenate(intensities.sigmas()))
      lattice_ids = lattice_ids.customized_copy(
        indices=lattice_ids.indices().concatenate(ids.indices()),
        data=lattice_ids.data().concatenate(ids.data()))
    assert miller_array_all.size() == lattice_ids.size()

    intensities = intensities.map_to_asu()
    intensities = intensities.customized_copy(anomalous_flag=True)
    intensities_p1 = intensities.expand_to_p1().merge_equivalents().array()
    intensities = intensities_p1.customized_copy(
      crystal_symmetry=intensities.crystal_symmetry())

  L = (miller_array_all, lattice_ids)
  L[0].crystal_symmetry().show_summary()
  from cctbx.merging import brehm_diederichs
  if params.nproc == 1:
    result_sets = brehm_diederichs.run(
      L, asymmetric=params.asymmetric, nproc=1, show_plot=params.show_plot,
      save_plot=params.save_plot)
  else:
    result_sets = brehm_diederichs.run_multiprocess(
      L, asymmetric=params.asymmetric, nproc=params.nproc,
      show_plot=params.show_plot, save_plot=params.save_plot)

  out_file = open('reindex.txt', 'wb')

  for reindexing_op, wedges in result_sets.iteritems():
    cb_op = sgtbx.change_of_basis_op(reindexing_op)
    for wedge in wedges:
      file_name = file_name_dict[wedge]
      if out_file is not None:
        print(file_name, cb_op.as_hkl(), file=out_file)
      basename = os.path.basename(file_name)
      out_name = os.path.splitext(basename)[0] + params.suffix + ".mtz"
      reader = any_reflection_file(file_name)
      assert reader.file_type() == 'ccp4_mtz'
      mtz_object = reader.file_content()
      if not cb_op.is_identity_op():
        print("reindexing %s" %file_name)
        mtz_object.change_basis_in_place(cb_op)
      mtz_object.write(out_name)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
