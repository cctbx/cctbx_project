from __future__ import division
from libtbx import easy_pickle as ep
import scipy as sp
import matplotlib.pyplot as plt

def make_png(image_pickle, integration_pickle, file_name=None, res=300):
  """ Write a png file visualizing integration results.
  :param image_pickle: path to image pickle file.
  :param integration_pickle: path to integration pickle file.
  :param res: resolution of output file in dpi (6x6 image size).
  :file_name: desired output file name. Defaults to the integration_pickle name.
  """

  if file_name is None:
    import os
    # Change extension of `image_pickle` to .png
    file_name = os.path.splitext(integration_pickle)[0] + ".png"

  # Load image pickle, and convert to image
  img_dict = ep.load(image_pickle)
  img_data = img_dict['DATA'].as_numpy_array()
  image = sp.misc.toimage(img_data, high=img_data.max())

  # Load integration pickle, and get coordinates of predictions
  int_d = ep.load(integration_pickle)
  predictions = int_d["mapped_predictions"][0]
  pred_coords = predictions.as_double().as_numpy_array() \
                           .reshape(2, len(predictions), order='F')

  # Create the figure
  fig = plt.figure(figsize=(6, 6))

  # 1st set of axes for the image
  ax = fig.add_subplot(1,1,1)
  ax.set_xlim(0, len(img_data[1]))
  ax.set_ylim(0, len(img_data[0]))
  ax.set_aspect('equal')
  ax.imshow(image, origin=None)

  # 2nd set of axes for the predictions
  ax2 = fig.add_axes(ax.get_position(), frameon=False)
  ax2.set_xlim(0, len(img_data[1]))
  ax2.set_ylim(0, len(img_data[0]))
  ax2.set_aspect('equal')

  ax2.plot(pred_coords[1], pred_coords[0], ms=2, linestyle='.', marker='o',
           alpha=0.3, markeredgecolor='r', markerfacecolor='none',
           markeredgewidth=0.3)

  plt.savefig(file_name, dpi=res, format='png')




