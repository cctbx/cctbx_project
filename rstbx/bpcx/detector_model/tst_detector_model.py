from __future__ import absolute_import, division, print_function
from scitbx import matrix
from rstbx.bpcx import sensor

def work_detector_model():
    # trivial test of instantiation

    i = matrix.col((1, 0, 0))
    j = matrix.col((0, 1, 0))
    k = matrix.col((0, 0, 1))

    d1 = i
    d2 = -j
    lim = (0,50)

    panel0 = sensor(matrix.col((-126, 144, -193)), d1, d2, lim, lim)

    panel0.set_origin(k)
    panel0.set_frame(k, i, j)

    assert(k.dot(matrix.col(panel0.origin)) > 0.999)

    print('OK')
    return

if __name__ == '__main__':
    work_detector_model()
