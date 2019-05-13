from __future__ import division
from __future__ import print_function
from six.moves import range
#-*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.print_pickle
#

"""
Simple utility for printing the contents of a cctbx.xfel pickle file
"""

import sys, os
from iotbx.detectors.cspad_detector_formats import detector_format_version as detector_format_function
from iotbx.detectors.cspad_detector_formats import reverse_timestamp
from cctbx import sgtbx # import dependency
from cctbx.array_family import flex

def keywise_printout(data):
  for key in data:
    if key == 'ACTIVE_AREAS':
      print(int(len(data[key])/4), "active areas, first one: ", list(data[key][0:4]))
    elif key == 'observations':
      print(key, data[key], "Showing unit cell/spacegroup:")
      obs = data[key][0]
      uc = obs.unit_cell()
      uc.show_parameters()
      obs.space_group().info().show_summary()
      d = uc.d(obs.indices())
      print("Number of observations:", len(obs.indices()))
      print("Max resolution: %f"%flex.min(d))
      print("Mean I/sigma:", flex.mean(obs.data())/flex.mean(obs.sigmas()))
      print("I/sigma > 1 count:", (obs.data()/obs.sigmas() > 1).count(True))
      print("I <= 0:", len(obs.data().select(obs.data() <= 0)))

      from cctbx.crystal import symmetry
      sym = symmetry(unit_cell = uc, space_group = obs.space_group())
      mset = sym.miller_set(indices = obs.indices(), anomalous_flag=False)
      binner = mset.setup_binner(n_bins=20)
      acceptable_resolution_bins = []
      binned_avg_i_sigi = []
      for i in binner.range_used():
        d_max, d_min = binner.bin_d_range(i)
        sel = (d <= d_max) & (d > d_min)
        sel &= obs.data() > 0
        intensities = obs.data().select(sel)
        sigmas = obs.sigmas().select(sel)
        n_refls = len(intensities)
        avg_i = flex.mean(intensities) if n_refls > 0 else 0
        avg_i_sigi = flex.mean(intensities / sigmas) if n_refls > 0 else 0
        acceptable_resolution_bins.append(avg_i_sigi >= 1.0)

      acceptable_resolution_bins = [acceptable_resolution_bins[i] if False not in acceptable_resolution_bins[:i+1] else False
                                    for i in range(len(acceptable_resolution_bins))]
      best_res = None
      for i, ok in zip(binner.range_used(), acceptable_resolution_bins):
        d_max, d_min = binner.bin_d_range(i)
        if ok:
          best_res = d_min
        else:
          break
      if best_res is None:
        print("Highest resolution with I/sigI >= 1.0: None")
      else:
        print("Highest resolution with I/sigI >= 1.0: %f"%d_min)

    elif key == 'mapped_predictions':
      print(key, data[key][0][0], "(only first shown of %d)"%len(data[key][0]))
    elif key == 'correction_vectors' and data[key] is not None and data[key][0] is not None:
      if data[key][0] is None:
        print(key, "None")
      else:
        print(key, data[key][0][0], "(only first shown)")
    elif key == "DATA":
      print(key,"len=%d max=%f min=%f dimensions=%s"%(data[key].size(),flex.max(data[key]),flex.min(data[key]),str(data[key].focus())))
    elif key == "WAVELENGTH":
      print("WAVELENGTH", data[key], ", converted to eV:", 12398.4187/data[key])
    elif key == "fuller_kapton_absorption_correction":
      print(key, data[key])
      if doplots:
        c = data[key][0]
        hist = flex.histogram(c, n_slots=30)
        from matplotlib import pyplot as plt
        plt.scatter(hist.slot_centers(), hist.slots())
        plt.show()

        obs = data['observations'][0]
        preds = data['mapped_predictions'][0]
        p1 = preds.select(c == 1.0)
        p2 = preds.select((c != 1.0) & (c <= 1.5))
        plt.scatter(preds.parts()[1], preds.parts()[0], c='g')
        plt.scatter(p1.parts()[1], p1.parts()[0], c='b')
        plt.scatter(p2.parts()[1], p2.parts()[0], c='r')
        plt.show()

    else:
      print(key, data[key])

def generate_streams_from_path(tar_or_other):
    from tarfile import ReadError
    import tarfile
    try:
      T = tarfile.open(name=tar_or_other, mode='r')
      K = T.getmembers()
      NT = len(K)
      for nt in range(NT):
        k = os.path.basename(K[nt].path)
        fileIO = T.extractfile(member=K[nt])
        yield fileIO,K[nt].path
    except ReadError as e:
      fileIO = open(tar_or_other,'rb')
      yield fileIO,tar_or_other

def generate_data_from_streams(args, verbose=False):
  from six.moves import cPickle as pickle
  for path in args:
    if not os.path.isfile(path):
      if verbose: print("Not a file:", path)
      continue

    # interpret the object as a tar of pickles, a pickle, or none of the above
    for fileIO,path in generate_streams_from_path(path):
      try:
        data = pickle.load(fileIO)
        if not isinstance(data, dict):
          if verbose: print("\nNot a dictionary pickle",path)
          continue
        else:
          if verbose: print("\nPrinting contents of", path)
          data["path"] = path
          yield data

      except pickle.UnpicklingError as e:
        if verbose: print("\ndoesn't unpickle",path)
      except EOFError as e:
        if verbose: print("\nEOF error",path)

if __name__=="__main__":
  args = sys.argv[1:]
  if "--break" in args:
    args.remove("--break")
    dobreak = True
  else:
    dobreak = False

  if "--plots" in args:
    args.remove("--plots")
    doplots = True
  else:
    doplots = False

  for data in generate_data_from_streams(args, verbose=True):
    if 'TIMESTAMP' in data:
      # this is how FormatPYunspecified guesses the address
      if not "DETECTOR_ADDRESS" in data:
        # legacy format; try to guess the address
        LCLS_detector_address = 'CxiDs1-0|Cspad-0'
        if "DISTANCE" in data and data["DISTANCE"] > 1000:
          # downstream CS-PAD detector station of CXI instrument
          LCLS_detector_address = 'CxiDsd-0|Cspad-0'
      else:
        LCLS_detector_address = data["DETECTOR_ADDRESS"]

      detector_format_version = detector_format_function(
        LCLS_detector_address, reverse_timestamp(data['TIMESTAMP'])[0])
      print("Detector format version:", detector_format_version)
      image_pickle = True
    else:
      print("Not an image pickle")
      image_pickle = False

    keywise_printout(data)

    if image_pickle:
      import dxtbx
      image = dxtbx.load(path)
      tile_manager = image.detectorbase.get_tile_manager(image.detectorbase.horizons_phil_cache)
      tiling = tile_manager.effective_tiling_as_flex_int(reapply_peripheral_margin = True)
      print(int(len(tiling)/4), "translated active areas, first one: ", list(tiling[0:4]))

    if dobreak:
      print("Entering break. The pickle is loaded in the variable 'data'")
      try:
        from IPython import embed
      except ImportError:
        import pdb; pdb.set_trace()
      else:
        embed()
