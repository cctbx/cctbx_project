import numpy as np


def is_outlier(points, thresh=3.5):
    """http://stackoverflow.com/a/22357811/2077270"""
    if len(points.shape) == 1:
        points = points[:, None]
    if points.size:
        median = np.median(points, axis=0)
    else:
        median = np.nan
    diff = np.sum((points - median) ** 2, axis=-1)
    diff = np.sqrt(diff)
    if diff.size:
        med_abs_deviation = np.median(diff)
    else:
        med_abs_deviation = np.nan
    if med_abs_deviation == 0:
        return np.zeros(points.shape[0], bool)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

