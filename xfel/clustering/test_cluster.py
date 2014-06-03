from __future__ import division
__author__ = 'zeldin'

import logging
from xfel.clustering.cluster import Cluster
import matplotlib.pyplot as plt

FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)

#0. Get data
path_to_some_integration_pickles = r"/Volumes/Data/XFEL/SFC_May_2014_XPP/" \
                                   r"processed_data/SFC_CRYO/Int_3.7_ShortCell"
all_longcell = Cluster.from_directories([path_to_some_integration_pickles],
                                        'SFC_longcell')
logging.info("data imported")

#1. Set up mega-plot
plt.figure(figsize=(20,15))
orr_axes = [plt.subplot2grid((3,3), (0,0)),
            plt.subplot2grid((3,3), (0,1)),
            plt.subplot2grid((3,3), (0,2))]
clust_ax = plt.subplot2grid((3,3), (1,0), colspan=3, rowspan=2)
all_longcell.visualise_orientational_distribution(orr_axes, cbar=False)
_, big_cluster = all_longcell.ab_cluster(5000, ax=clust_ax)

plt.show()


#
#best_data = big_cluster.total_intensity_filter(res=3.7,
#                                               completeness_threshold=0.95,
#                                               plot=True)
#
#print big_cluster.info
#big_cluster.dump_file_list()
