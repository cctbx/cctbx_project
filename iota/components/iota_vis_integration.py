from __future__ import absolute_import, division, print_function

'''
Author      : Zeldin O.B., Lyubimov A.Y.
Created     : 10/12/2014
Last Changed: 04/13/2015
Description : Creates a PNG file visualizing integration results.
'''

from libtbx import easy_pickle as ep
import numpy as np
import matplotlib as mpl

def make_png(image_pickle, integration_pickle, file_name=None, res=600,
             show_spots=True):
  """ Write a png file visualizing integration results.
  :param image_pickle: path to image pickle file.
  :param integration_pickle: path to integration pickle file.
  :param res: resolution of output file in dpi (6x6 image size).
  :param show_spots: Uses the `obsspot` field from the `correction_vector` key to plot spot positions. Fails silently if `correction_vector` does not exist. Use option indexing.verbose_cv=True in cxi.index to make sure that this is created.
  :file_name: desired output file name. Defaults to the integration_pickle name.
  """

  if file_name is None:
    import os
    # Change extension of `image_pickle` to .png
    file_name = os.path.splitext(integration_pickle)[0] + ".png"

  # Load image pickle, and convert to image
  img_dict = ep.load(image_pickle)
  img_data = img_dict['DATA'].as_numpy_array()

  # Load integration pickle, and get coordinates of predictions
  int_d = ep.load(integration_pickle)
  predictions = int_d["mapped_predictions"][0]
  pred_coords = predictions.as_double().as_numpy_array() \
                           .reshape(2, len(predictions), order='F')

  if show_spots and 'correction_vectors' in int_d:
    spot_coords = [x['obsspot'] for x in int_d['correction_vectors'][0]]
    spot_coords = np.array(spot_coords).transpose()
    plot_spots = True
  else:
    plot_spots = False

  # Get some other useful info from integration pickle
  point_group = int_d['pointgroup']
  unit_cell = int_d['current_orientation'][0].unit_cell().parameters()
  unit_cell = ', '.join("{:.1f}".format(u) for u in unit_cell)
  mosaicity = int_d['mosaicity']

  # Create the figure
  import matplotlib.pyplot as plt
  fig = plt.figure(figsize=(6, 6))

  # 1st set of axes for the image
  ax = fig.add_subplot(1,1,1)
  ax.set_xlim(0, len(img_data[1]))
  ax.set_ylim(0, len(img_data[0]))
  ax.set_aspect('equal')

  # To display image properly, need to set img limit to 0.01th percentile
  # a few maxed-out pixels will make the entire image appear blank white
  clim = (img_data.min(), np.percentile(img_data, 99.99))
  ax.imshow(img_data, origin=None, cmap='Greys', clim=clim)

  # 2nd set of axes for the predictions
  ax2 = fig.add_axes(ax.get_position(), frameon=False)  # superimposed axes
  ax2.set_xlim(0, len(img_data[1]))
  ax2.set_ylim(0, len(img_data[0]))
  ax2.set_aspect('equal')

  ax2.plot(pred_coords[1], pred_coords[0], ms=2, linestyle='.', marker='o',
           alpha=0.3, markeredgecolor='r', markerfacecolor='none',
           markeredgewidth=0.3)

  if plot_spots:
    ax2.plot(spot_coords[1], spot_coords[0], ms=2, linestyle='.', marker='d',
             alpha=0.3, markeredgecolor='b', markerfacecolor='none',
             markeredgewidth=0.3, linewidth=1)


  plt.title("Unit cell: {} ({}) \nNominal mosaicity: {}" \
            .format(point_group, unit_cell, mosaicity))
  plt.savefig(file_name, dpi=res, format='png')


def cv_png(image_pickle, integration_pickle, file_name=None, res=600,
             show_spots=True):
  """ Write a png file visualizing the correction vectors.
  :param image_pickle: path to image pickle file.
  :param integration_pickle: path to integration pickle file.
  :param res: resolution of output file in dpi (6x6 image size).
  :file_name: desired output file name. Defaults to the integration_pickle name.
  """

  cmap = mpl.cm.jet

  if file_name is None:
    import os
    # Change extension of `image_pickle` to .png
    file_name = os.path.splitext(integration_pickle)[0] + ".png"

  # Load image pickle, and convert to image
  img_dict = ep.load(image_pickle)
  img_data = img_dict['DATA'].as_numpy_array()

  # Load integration pickle, and get coordinates of predictions
  int_d = ep.load(integration_pickle)
  spot_coords = [x['obsspot'] for x in int_d['correction_vectors'][0]]
  spot_coords = np.array(spot_coords).transpose()
  pred_coords = [x['predspot'] for x in int_d['correction_vectors'][0]]
  pred_coords = np.array(pred_coords).transpose()
  corr_vecs = pred_coords - spot_coords
  corr_factors = np.sqrt(abs(corr_vecs[0]**2 - corr_vecs[1]**2))

  # Get some other useful info from integration pickle
  point_group = int_d['pointgroup']
  unit_cell = int_d['current_orientation'][0].unit_cell().parameters()
  unit_cell = ', '.join("{:.1f}".format(u) for u in unit_cell)
  mosaicity = int_d['mosaicity']

  # Create the figure
  import matplotlib.pyplot as plt
  fig = plt.figure(figsize=(6, 6))

  # 1st set of axes for the image
  ax = fig.add_subplot(1,1,1)
  plt.subplots_adjust(right=0.8)
  ax.set_xlim(0, len(img_data[1]))
  ax.set_ylim(0, len(img_data[0]))
  ax.set_aspect('equal')

  # To display image properly, need to set img limit to 0.01th percentile
  # a few maxed-out pixels will make the entire image appear blank white
  clim = (img_data.min(), np.percentile(img_data, 99.99))
  ax.imshow(img_data, origin=None, cmap='Greys', clim=clim)

  # 2nd set of axes for the predictions
  ax2 = fig.add_axes(ax.get_position(), frameon=False)  # superimposed axes
  ax2.set_xlim(0, len(img_data[1]))
  ax2.set_ylim(0, len(img_data[0]))
  ax2.set_aspect('equal')


  norm = mpl.colors.Normalize(vmin=0, vmax=round(max(corr_factors)))
  ax2.scatter(spot_coords[1], spot_coords[0], c=norm(corr_factors),
              cmap=cmap, alpha=0.5, s=5, linewidths=0)

  cax = plt.axes([0.85, 0.15, 0.05, 0.69])
  cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
  cbar.set_label('Spot-prediction distance [px]')

  ax.set_title("Prediction-offset values\nfor spotfinder results")
  plt.savefig(file_name, dpi=res, format='png')
