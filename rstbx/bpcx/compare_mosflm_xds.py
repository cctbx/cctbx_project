from __future__ import absolute_import, division, print_function
# Compare intensities and positions of measurements from XDS INTEGRATE (which
# are not necessarily correctly LP corrected) and from Mosflm via sortmtz and
# scala to sum partials but not merge or scale e.g.
#
# #!/bin/sh
#
# sortmtz hklin thau_2_001.mtz hklout sorted.mtz << eof
# H K L M/ISYM BATCH
# eof
#
# scala hklin sorted.mtz hklout summed.mtz << eof
# run 1 batch 1 to 1800
# scales constant
# cycles 0
# sdcorrection noadjust norefine both 1.0 0.0 0.0
# output unmerged
# eof
#
# N.B. all comparisons done with data in memory, so careful with the sizes of
# example data sets to use. Also assumes X, Y definitions between Mosflm and
# XDS are inverted, and computes deviations in terms of pixels and images.
# This requires reading of the batch headers. Try also comparing with pixels
# and degrees (perhaps small differences in postrefinement are affecting
# things?)

import sys
import math
from iotbx import mtz
from cctbx.array_family import flex
from annlib_ext import AnnAdaptor as ann_adaptor
from cctbx.sgtbx import space_group, space_group_symbols
from cctbx.miller import map_to_asu
from six.moves import range
from six.moves import zip

def meansd(values):
    assert(len(values) > 3)

    mean = sum(values) / len(values)
    var = sum([(v - mean) * (v - mean) for v in values]) / (len(values) - 1)

    return mean, math.sqrt(var)

def cc(a, b):

    assert(len(a) == len(b))

    ma, sa = meansd(a)
    mb, sb = meansd(b)

    r = (1 / (len(a) - 1)) * sum([((a[j] - ma) / sa) * ((b[j] - mb) / sb)
                                  for j in range(len(a))])

    return r

def get_phi_range(mtz_object):
    bs = mtz_object.batches()

    rs = [(b.phiend() - b.phistt()) for b in bs]

    return meansd(rs)[0]

def get_hkl_xyz_isigi(mtz_file):
    m = mtz.object(mtz_file)

    r = get_phi_range(m)

    sg = m.space_group()

    xdet_col = m.get_column('XDET')
    ydet_col = m.get_column('YDET')
    rot_col = m.get_column('ROT')
    i_col = m.get_column('I')
    sigi_col = m.get_column('SIGI')

    mi_s = m.extract_miller_indices()
    map_to_asu(sg.type(), False, mi_s)

    xdet_s = xdet_col.extract_values(not_a_number_substitute = 0.0)
    ydet_s = ydet_col.extract_values(not_a_number_substitute = 0.0)
    rot_s = rot_col.extract_values(not_a_number_substitute = 0.0)
    i_s = i_col.extract_values(not_a_number_substitute = 0.0)
    sigi_s = sigi_col.extract_values(not_a_number_substitute = 0.0)

    hkl_xyz_isigi = []

    # N.B. swapping X, Y here

    for j in range(mi_s.size()):
        hkl_xyz_isigi.append(((mi_s[j][0], mi_s[j][1], mi_s[j][2]),
                              (ydet_s[j], xdet_s[j], rot_s[j] / r),
                              (i_s[j], sigi_s[j])))

    return hkl_xyz_isigi

def read_xds_integrate(xds_integrate_file):

    # N.B. in here need to reduce indices to asymmetric unit

    sg_num = 0
    r = 0.0

    for record in open(xds_integrate_file):
        if not '!' in record[:1]:
            break
        if 'SPACE_GROUP_NUMBER' in record:
            sg_num = int(record.split()[-1])
        if 'OSCILLATION_RANGE' in record:
            r = float(record.split()[-1])

    if not sg_num:
        raise RuntimeError('spacegroup missing')

    if not r:
        raise RuntimeError('rotation missing')

    sg = space_group(space_group_symbols(sg_num).hall())

    observations = []

    hkls = flex.miller_index()
    xyzs = []
    isigmas = []

    for record in open(xds_integrate_file):
        if '!' in record[:1]:
            continue
        values = record.split()
        hkls.append( [int(h) for h in values[:3]] )
        xyzs.append((float(values[5]), float(values[6]), float(values[7])))
        isigmas.append([float(x) for x in values[3:5]])

    map_to_asu(sg.type(), False, hkls)

    for hkl, xyz, isigma in zip(hkls, xyzs, isigmas):
        observations.append((hkl, xyz, isigma))

    return observations

def main(mtz_file, xds_integrate_file):
    mos_hkl_xyz_isigi = get_hkl_xyz_isigi(mtz_file)

    print('Read %d observations from %s' % (len(mos_hkl_xyz_isigi), mtz_file))

    xds_hkl_xyz_isigi = read_xds_integrate(xds_integrate_file)

    print('Read %d observations from %s' % \
          (len(xds_hkl_xyz_isigi), xds_integrate_file))

    # treat XDS as reference, mosflm as query (arbitrary)

    reference = flex.double()
    query = flex.double()

    for hkl, xyz, isigi in xds_hkl_xyz_isigi:
        reference.append(xyz[0])
        reference.append(xyz[1])
        reference.append(xyz[2])

    for hkl, xyz, isigi in mos_hkl_xyz_isigi:
        query.append(xyz[0])
        query.append(xyz[1])
        query.append(xyz[2])

    ann = ann_adaptor(data = reference, dim = 3, k = 1)
    ann.query(query)

    i_s_mos = []
    i_s_xds = []

    for j in range(len(mos_hkl_xyz_isigi)):
        c = ann.nn[j]
        if xds_hkl_xyz_isigi[c][0] == mos_hkl_xyz_isigi[j][0]:
            i_s_mos.append(mos_hkl_xyz_isigi[j][2][0])
            i_s_xds.append(xds_hkl_xyz_isigi[c][2][0])

    print('Matched %d observations' % len(i_s_mos))

    print(cc(i_s_mos, i_s_xds))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
