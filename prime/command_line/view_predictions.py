from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.view_predictions
import argparse
import cPickle as pickle
import matplotlib.pyplot as plt

def main(img_fname, int_fname):
  """
  read image pickle and integration pickle then show predictions on the top
  of the diffraction pattern.
  """

  img_pickle = pickle.load(open(img_fname,"rb"))
  int_pickle = pickle.load(open(int_fname,"rb"))
  obs = int_pickle["observations"][0]
  obs_asu = obs.map_to_asu()
  mapped_predictions = int_pickle['mapped_predictions'][0]
  mm_predictions = mapped_predictions
  spot_pred_x_mm = [pred[0]-int_pickle['xbeam'] for pred in mm_predictions]
  spot_pred_y_mm = [pred[1]-int_pickle['ybeam'] for pred in mm_predictions]
  #zdata = img_pickle['DATA'].set_selected(img_pickle['DATA']<=800, 0)
  #zdata = zdata.set_selected(zdata > 800, 1)
  #data = zdata.as_numpy_array()
  #plt.hist(data.ravel(), bins=256, range=(0,3500))
  #imgplot = plt.imshow(data)
  #plt.colorbar()
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.scatter(spot_pred_x_mm, spot_pred_y_mm)
  shown_indices = [(-10,-4,-3), (10,4,3), (-10,-7,-2), (10,7,2), (-12,-8,-2), (12,8,2)]
  for ind, ind_asu, I, x, y in zip(obs.indices(), obs_asu.indices(), obs.data(), spot_pred_x_mm, spot_pred_y_mm):
    if ind_asu in shown_indices:
      ax.text(x, y, str(ind)+' %8.1f'%(I), fontsize=10)
  ax.axis([0,4096,0,4096])
  plt.show()


if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description='View integration predictions on the top of a diffraction pattern'
  )
  parser.add_argument(
    'img',
    metavar='IMG',
    help='Diffraction pattern as a pickle file'
  )
  parser.add_argument(
    'int',
    metavar='INT',
    help='Integration result pickle'
  )
  args = parser.parse_args()
  main(args.img, args.int)
