from __future__ import absolute_import, division, print_function
from six.moves import range
import h5py
import copy
from scitbx.array_family import flex
from matplotlib import pyplot as plt

if __name__=="__main__":
  import sys
  from xfel.vonHamos import fit_3_gaussian

  fname = sys.argv[1]
  f = h5py.File(fname, 'r')

  histograms = f['mydata/histograms']
  cspad = f['mydata/cspad_sum']
  image = cspad[:,:,1]
  fig=True
  if fig:
    plt.figure()
    plt.imshow(image,interpolation="none")
    plt.title("Raw ADU count, no histogram fitting")
    plt.colorbar()
    plt.clim(0, 25000)

  inelastic = copy.deepcopy(image)
  single_process = False

  if single_process:
    for islow in range(185):
      for jfast in range(388):
        histo1 = histograms[islow*388+jfast,:]
        nphotons = fit_3_gaussian.test_fit(histo1,plot=False)
        print(islow*388+jfast, "of %d, # photons= %d"%(len(histograms),nphotons))
        inelastic[islow,jfast]=nphotons
  else:
    # same thing, multiprocess
    from libtbx import easy_mp
    import libtbx.introspection
    nproc = libtbx.introspection.number_of_processors(return_value_if_unknown=4)

    def run_fast(args):
      islow = args[0]
      result_val = flex.int(388)
      for jfast in range(388):
        histo1 = histograms[islow*388+jfast,:]
        nphotons = fit_3_gaussian.test_fit(histo1,plot=False)
        print(islow*388+jfast, "of %d, # photons= %d"%(len(histograms),nphotons))
        result_val[jfast]=nphotons
      return {islow:result_val}

    iterable = [(i,) for i in range(185)]
    results = easy_mp.parallel_map(
      func=run_fast,
      iterable=iterable,
      processes=nproc,
      method="multiprocessing",
      preserve_order=True
    )

    for result in results:
      for key in result.keys():
        values = result[key]
        for i in range(len(values)):
          inelastic[key,i] = values[i]

  # save fitting result:
  # some heuristics to output a proper filename
  if fname.find('_offline')>=0:
    filename = fname.replace('_offline','_fitted')
  else:
    filename = fname.replace('.h5','_fitted.h5')
  ofile = h5py.File(filename, 'w')
  group = ofile.create_group("mydata")
  group.create_dataset("cspad_photons", data=inelastic)
  ofile.close()
  print('*** fitted image saved as', filename)

  if fig:
    plt.figure()
    plt.imshow(inelastic,interpolation="none")
    plt.title("Dark + inelastic + elastic histogram fitting")
    plt.colorbar()
    plt.clim(0, 220)
    plt.show()
