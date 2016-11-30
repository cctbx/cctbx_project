from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.check_wavelength
"""
Author      : Uervirojnangkoorn, M.
Created     : 7/29/2016
Description : prime.check_wavelength /path/to/image/pickles sase.pickle
"""
import argparse
import cPickle as pickle
import os
import numpy as np
import matplotlib.pyplot as plt

def calculate_centroid(sase_pickle_file, start_ev, end_ev):
  sase_pickle = pickle.load(open(sase_pickle_file,'rb'))
  n_pixels = 0
  sase_centroids = {}
  sase_lambda_centroids = {}
  bin_range = []
  for key in sase_pickle.keys():
    n_pixels = len(sase_pickle[key])
    break
  if n_pixels:
    ev_per_pixel = (end_ev - start_ev)/n_pixels
    bin_range = start_ev + (np.array(range(n_pixels))*ev_per_pixel)
    for key in sase_pickle.keys():
      sase_centroids[key] = np.sum(bin_range * np.array(sase_pickle[key]))/np.sum(np.array(sase_pickle[key]))
      sase_lambda_centroids[key] = 12398.4187/sase_centroids[key]
  else:
    raise('Error reading sase pickle.')
    exit(1)
  return sase_pickle, sase_centroids, sase_lambda_centroids, bin_range

def compare_wavelength(img_pickle_dir, sase_spectrum, sase_centroids, sase_lambda_centroids, bin_range):
  #read sase pickle
  img_set = {}
  for img_fname in os.listdir(img_pickle_dir):
    if img_fname.endswith('.pickle'):
        img_pickle = pickle.load(open(os.path.join(img_pickle_dir, img_fname), 'rb'))
        img_set[img_fname] = 12398.4187/img_pickle['WAVELENGTH']
        img_fname_only = img_fname.split('.')[0]
        if img_fname_only in sase_centroids:
          print img_fname, img_set[img_fname], img_pickle['WAVELENGTH'], sase_centroids[img_fname_only], sase_lambda_centroids[img_fname_only], img_set[img_fname]-sase_centroids[img_fname_only]
          """
          plt.plot(bin_range, sase_spectrum[img_fname_only], c='g', linewidth=2)
          plt.scatter(sase_centroids[img_fname_only], 0, s=20, marker='x', color='r')
          plt.scatter(img_set[img_fname], 0, s=20, marker='x', color='b')
          plt.title(img_fname+' header: %8.4f calculated: %8.4f' %(img_set[img_fname], sase_centroids[img_fname_only]))
          plt.show()
          """

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Compare image-header wavelength with the calculated centroid of energy spectrum')
  parser.add_argument('img_pickle_dir', metavar='dir', help='Path to image pickles')
  parser.add_argument('sase_pickle_file', metavar='sase', help='Pickle file storing sase spectrum (see prime.get_SASE)')
  parser.add_argument('start_ev', metavar='start (eV)', type=int, help='Lower limit of the energy range (eV)')
  parser.add_argument('end_ev', metavar='end (ev)', type=int, help='Upper limit of the energy range (eV)')
  parser.add_argument('mode', metavar='mode', help='c for comparison, w for write-out only ')

  args = parser.parse_args()
  sase_spectrum, sase_centroids, sase_lambda_centroids, bin_range = calculate_centroid(args.sase_pickle_file, args.start_ev, args.end_ev)
  if args.mode == 'c':
    compare_wavelength(args.img_pickle_dir, sase_spectrum, sase_centroids, sase_lambda_centroids, bin_range)
  sase_pickle_filename = args.sase_pickle_file.split('.')[0]
  pickle.dump(sase_lambda_centroids, open(sase_pickle_filename+'_lambda.pickle', 'wb'))
