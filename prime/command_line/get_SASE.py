from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.get_SASE
"""
Author      : Uervirojnangkoorn, M.
Created     : 7/29/2016
Description : prime.get_SASE experiment_name run_no
"""
import argparse
from psana import *
import cPickle as pickle

def main(exp, run_from, run_to):
  e_dict = {}
  for run in range(run_from, run_to):
    ds = DataSource('exp=%s:run=%d'%(exp, run))
    feeSpec = Detector('FEE-SPEC0')
    config = ds.env().configStore()
    events = ds.events()
    for evt in events:
      control_data = config.get(ControlData.ConfigV3)
      pvLabels = control_data.pvLabels()
      fname_value = ""
      for label in pvLabels:
        if label.name() == "filename":
          fname_value = label.value()
      feeSpec_evt = feeSpec.get(evt)
      if feeSpec_evt is not None:
        spectrum = feeSpec_evt.hproj()
        if fname_value.find('test')==-1 and fname_value.strip() != '' and fname_value.find('image')==-1:
          print fname_value, spectrum
          e_dict[fname_value] = spectrum
  #write out energy-fname pickle
  pickle.dump(e_dict, open(exp+'_'+str(run_from)+'_'+str(run_to)+'_sase.pickle', "wb"))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    description='DAQ Control IOC info XTC dumper'
  )
  parser.add_argument(
    'exp',
    metavar='EXP',
    help='the experiment number'
  )
  parser.add_argument(
    'run_from',
    metavar='RUN FROM',
    type=int,
    help='the starting run number'
  )
  parser.add_argument(
    'run_to',
    metavar='RUN TO',
    type=int,
    help='the ending run number'
  )
  args = parser.parse_args()
  main(args.exp, args.run_from, args.run_to)
