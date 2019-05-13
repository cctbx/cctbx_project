from __future__ import division
from __future__ import print_function
from six.moves import range
import h5py, os, sys
import numpy as np
from libtbx.phil import parse
from xfel.euxfel.agipd_cxigeom2nexus import agipd_cxigeom2nexus

message = '''Program to calculate the maximum, mean and/or stdev image from an EU-XFEL run. Applicable only for the AGIPD detector. Needs a) cxi file b) what mode to use c) geometry file for generation of the master file d) detector distance (optional)'''

phil_scope = parse("""
  cxi_file = None
    .type = str
    .help = cheetah file used to read in image data (.cxi format)
  composite_mode = *max mean stdev
    .type = str
    .multiple = True
    .help = specify what statistical metric(s) should be used for output image
  geom_file = None
    .type = str
    .help = geometry file to be read in for AGIPD detector (.geom).
  detector_distance = None
    .type = float
    .help = AGIPD Detector distance
""")




class composite_image_writer(object):
  ''' class to write composite image'''
  def __init__(self, args):
    self.params_from_phil(args)
    n_frames = None

  def params_from_phil(self,args):
    user_phil = []
    for arg in args:
      if os.path.isfile(arg):
        user_phil.append(parse(file_name=arg))
      else:
        try:
          user_phil.append(parse(arg))
        except Exception, e:
          raise Sorry("Unrecognized argument: %s"%arg)
    self.params = phil_scope.fetch(sources=user_phil).extract()

  def copy_attributes(self,src, dest):
    for attr in src.attrs:
      dest.attrs[attr] = src.attrs[attr]

  def recursive_copy(self,src, dest, mode='max', n_frames=None):
    print(src, type(src))
    self.copy_attributes(src, dest)

    assert n_frames is not None, 'Need to provide n_frames'
    assert type(src) in [h5py._hl.group.Group, h5py._hl.files.File]
    for key in src:
      if type(src[key]) in [h5py._hl.group.Group, h5py._hl.files.File]:
        dest_child = dest.create_group(key)
        self.recursive_copy(src[key], dest_child,mode=mode,n_frames=n_frames)
      elif type(src[key]) is h5py._hl.dataset.Dataset:
        dataset = src[key]
        print(key, dataset.shape)
        if dataset.shape == (n_frames, 8192, 128):
          if dataset.name != "/entry_1/data_1/data":
            print("Skipping data block", dataset.name)
            continue
          print('=====================================')
          if mode == 'max':
            dmax = np.zeros((8192, 128))
            dmax[:] = -np.inf
            for i in range(n_frames):
              frame = dataset[i]
              dmax = np.maximum(dmax, frame)
            result = dmax.reshape(1 ,8192, 128)
          elif mode == 'mean':
            dsum = np.zeros((8192, 128))
            for i in range(n_frames):
              frame = dataset[i]
              dsum += frame
            result = (dsum/n_frames).reshape(1 ,8192, 128)
          elif mode == 'stdev':
            dsum = np.zeros((8192, 128))
            dsumsq = np.zeros((8192, 128))
            for i in range(n_frames):
              frame = dataset[i]
              dsum += frame
              dsumsq += frame*frame
            result = dsumsq - dsum*dsum
            result = np.sqrt(result)/n_frames

          new_dataset = dest.create_dataset(os.path.basename(dataset.name), data = result)
        else:
          new_dataset = dest.create_dataset(os.path.basename(dataset.name), data = dataset)
        self.copy_attributes(dataset, new_dataset)
      else:
        assert False, type(src)

if __name__ == '__main__':

  image_writer = composite_image_writer(sys.argv[1:])
  handle = h5py.File(image_writer.params.cxi_file, 'r')
  n_frames = len(handle['entry_1/instrument_1/detector_1/distance'])
  for mode in image_writer.params.composite_mode:
    outfile = os.path.splitext(image_writer.params.cxi_file)[0]+"_"+mode+".h5"
    output_handle = h5py.File(outfile, 'w')
    image_writer.recursive_copy(handle, output_handle,mode=mode,n_frames=n_frames)
    nexus_helper_str = ['cxi_file='+outfile] \
                     + ['geom_file='+image_writer.params.geom_file]
    if image_writer.params.detector_distance is not None:
      nexus_helper_str += ['detector_distance='+str(image_writer.params.detector_distance)]
    nexus_helper = agipd_cxigeom2nexus(nexus_helper_str)
    nexus_helper.create_nexus_master_file()
