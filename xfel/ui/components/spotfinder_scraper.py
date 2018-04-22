from __future__ import division

from xfel.util.reflection_length import ReflectionsRadialLengthsFromFiles
from dials.array_family import flex
import os

def get_spot_length_stats(outdir, ref_stats=None):
  def join(filename):
    return os.path.join(outdir, str(filename))
  files = map(join, os.listdir(outdir))
  datablocks = sorted([f for f in files if 'datablock' in f])
  strong = sorted([f for f in files if 'strong' in f])
  if ref_stats:
    timestamps, n_strong = ref_stats
    ts_strong = timestamps.select(n_strong > 0)
    assert len(strong) == len(datablocks), \
      "Mismatch between n_strong: %d, n_datablocks: %d" % (len(strong), len(datablocks))
    assert len(strong) == len(timestamps), \
      "Mismatch between n_strong: %d, n_timestamps: %d" % (len(strong), len(timestamps))
    strong_spot_timestamps = flex.double()
    strong_spot_lengths = flex.double()
    for i in xrange(len(datablocks)):
      ts, s, db = timestamps[i], strong[i], datablocks[i]
      lengths = ReflectionsRadialLengthsFromFiles([s, db]).get_spot_lengths_px()
      for length in lengths:
        strong_spot_lengths.append(length)
        strong_spot_timestamps.append(ts)
    return strong_spot_timestamps, strong_spot_lengths
  else:
    strong_spot_lengths = ReflectionsRadialLengthsFromFiles(
      datablocks + strong).get_spot_lengths_px()
    return (flex.double(), strong_spot_lengths)
