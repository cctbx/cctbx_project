
from __future__ import absolute_import, division, print_function
from .Cluster import Cluster

import logging
FORMAT = '%(levelname)s %(module)s.%(funcName)s: %(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMAT)
""" A test cluster script for playing with the class """
path_to_some_integration_pickles = r"/usr/local/Dropbox/Stanford_Postdoc/CODING/Fraser_test_data"

test_cluster = Cluster.from_directories([path_to_some_integration_pickles],
                                        'test_script')
logging.info("data imported")
clust = test_cluster.point_group_filer('P222')
clust = clust.ab_cluster(80, plot=False)
big_cluster = max(clust, key=lambda x: len(x.members))

best_data = big_cluster.total_intensity_filter(res=8, completeness_threshold=0.5, plot=False)

print(best_data.info)
best_data.dump_file_list()
