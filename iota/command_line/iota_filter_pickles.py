from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME iota.filter_pickles

import h5py
import os, shutil, argparse, glob, ntpath
import numpy as np
from six.moves import range, zip, map


def main(hdf5Filename, data, destDir, powder):
  """
  ideally, if pulse is on, you should expect to have 1 - 1 -1 for all events
  otherwise 0 -1 - 1. Currently, we have more than one pickle read out for
  each rotatation angle. Solution is picking the pickle with maximum value for Ipm3.
  """
  f = h5py.File(hdf5Filename)

  evr92 = f['evr92'] #pulse
  evr95 = f['evr95'] #
  evr97 = f['evr97'] #rayonix
  rayonixTime = f['rayonix_time']
  rot_wrap = range(len(rayonixTime))
  #rotation angles willl not be recorded if the motor does not moved
  if 'sam_rot_wrap' in f:
    rot_wrap = f['sam_rot_wrap']
  ipm3 = f['ipm3']

  goodRayonixTime = {}
  goodRayonixIpm3 = {}
  for _rt, _evr92, _evr95, _evr97, rot, _ipm3 in zip(rayonixTime, evr92, evr95, evr97, rot_wrap, ipm3):
    if (_evr95 == 1 and _evr97 == 1) or powder:
      if (_evr92,_evr95, _evr97, rot) in goodRayonixTime:
        #found more than one
        goodRayonixTime[(_evr92,_evr95, _evr97, rot)].append(_rt)
        goodRayonixIpm3[(_evr92,_evr95, _evr97, rot)].append(_ipm3)
      else:
        goodRayonixTime[(_evr92,_evr95, _evr97, rot)] = [_rt]
        goodRayonixIpm3[(_evr92,_evr95, _evr97, rot)] = [_ipm3]

  goodKey = []
  goodRayonixTimeList = []
  goodRayonixIpm3List = []
  for key, val in goodRayonixTime.items():
    if len(val) > 1:
      goodIndex = np.argmax(goodRayonixIpm3[key])
    else:
      goodIndex = 0
    goodKey.append(key)
    goodRayonixTimeList.append(goodRayonixTime[key][goodIndex])
    goodRayonixIpm3List.append(goodRayonixIpm3[key][goodIndex])

  isort = np.argsort(goodRayonixTimeList)
  sortedRayonixTime = np.array(goodRayonixTimeList)[isort]
  sortedRayonixIpm3 = np.array(goodRayonixIpm3List)[isort]
  sortedKey = np.array(goodKey)[isort]

  #make destination directory
  destDirOn = destDir + '/on'
  destDirOff = destDir + '/off'
  try:
    os.makedirs(destDirOn)
    os.makedirs(destDirOff)
  except OSError:
    print ("WARNING: Problem creating folder. If the folder exists, it will get overwritten.")
    pass

  txtOn = ''
  txtOff = ''
  txtError = ''
  cnSuccess, cnErrors = (0, 0)
  for _key, _rt, _ipm3  in zip(sortedKey, sortedRayonixTime, sortedRayonixIpm3):
    filename_pattern = 'shot-'+str(_rt)+'*'
    #copy to a good directory
    try:
      _evr92, _evr95, _evr97, _rot = tuple(_key)
      txtOut = ','.join(map(str, list(_key)))
      txtOut += ',%d, %f\n'%(_rt,_ipm3)
      #determine output folder
      if _evr92 == 0:
        destDirSelected = destDirOff
        txtOff += txtOut
      else:
        destDirSelected = destDirOn
        txtOn += txtOut
      #copy any files that match with the pattern
      for filename in glob.glob(data + '/' + filename_pattern):
        shutil.copy(filename, destDirSelected+'/'+ntpath.basename(filename))
        print (filename, _evr92, destDirSelected, 'success')
      cnSuccess += 1
    except IOError as err:
      print (err)
      txtError += ','.join(map(str, list(_key)))
      txtError += ',%d, %f\n'%(_rt,_ipm3)
      cnErrors += 1
      pass
  #write out log file
  f = open(destDirOn+'/shotangles.txt','w')
  f.write(txtOn)
  f.close()
  f = open(destDirOff+'/shotangles.txt','w')
  f.write(txtOff)
  f.close()
  f = open(destDir+'/errors.txt','w')
  f.write(txtError)
  print ("Successfully moved %d files"%(cnSuccess))
  print ("%d files not found"%(cnErrors))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    description='Filers out pickle files read out from xtc with insignificant events'
  )
  parser.add_argument(
    'hdf5Filename',
    metavar='HDF5 Filename',
    help='HDF5 Filename'
  )
  parser.add_argument(
    'data',
    metavar='DATA',
    help='Path to all pickles'
  )
  parser.add_argument(
    'destDir',
    metavar='DESTINATION FOLDER',
    help='Path to destination folder'
  )
  parser.add_argument('--powder', dest='powder', action='store_const',
                    const=True, default=False,
                    help='Flag for powder patterns')
  args = parser.parse_args()
  if args.powder:
    print ("Powder diffraction. Ignore event codes")
  main(args.hdf5Filename, args.data, args.destDir, args.powder)
