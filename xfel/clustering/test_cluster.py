
from __future__ import division

import logging

from xfel.clustering.cluster import Cluster

FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
""" A test cluster script for playing with the class """
path_to_some_integration_pickles = r"/Volumes/Data/XFEL/SFC_May_2014_XPP/" \
                                   r"processed_data/SFC_CRYO/Int_3.7_LongCell"

all_longcell = Cluster.from_directories([path_to_some_integration_pickles],
                                        'SFC_longcell')
logging.info("data imported")
#clust = test_cluster.point_group_filter('P222')
all_longcell.visualise_orientational_distribution()
#big_cluster = all_longcell.ab_cluster(5000, plot=False, max_only=True)
#
#best_data = big_cluster.total_intensity_filter(res=3.7,
#                                               completeness_threshold=0.95,
#                                               plot=True)
#
#print big_cluster.info
#big_cluster.dump_file_list()
