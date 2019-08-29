from __future__ import absolute_import, division, print_function
from six.moves import range

from xfel.util.reflection_length import ReflectionsRadialLengthsFromFiles
from dials.array_family import flex
import os
from six.moves import zip
from six.moves import map

def get_spot_length_stats(outdir, ref_stats=None):
  def join(filename):
    return os.path.join(outdir, str(filename))

  files = map(join, os.listdir(outdir))
  datablocks = sorted([f for f in files if 'datablock' in f])
  strong = sorted([f for f in files if 'strong' in f])
  strong_spot_timestamps = flex.double()
  strong_spot_lengths = flex.double()
  strong_spot_intensities = flex.double()

  def process_pair(refl, block, ts=None):
    measurer = ReflectionsRadialLengthsFromFiles([refl, block])
    lengths = measurer.get_spot_lengths_px()
    intensities = measurer.get_intensities()
    assert len(lengths) == len(intensities)
    for i, length in enumerate(lengths):
      strong_spot_lengths.append(length)
      strong_spot_intensities.append(intensities[i])
      if ts:
        strong_spot_timestamps.append(ts)

  if ref_stats:
    timestamps, n_strong = ref_stats
    # Can't use reference timestamps at all if there are fewer datablocks than timestamps.
    assert len(datablocks) == len(timestamps), \
      "Mismatch between n_strong: %d, n_timestamps: %d" % (len(strong), len(timestamps))
    # In case of unanticipated errors during spotfinding, strong.pickle may not be written.
    # Get a list strong of the same length as datablocks with placeholder None values.
    if len(strong) != len(datablocks):
      strong = [db.split("datablock.json")[0] + "strong.pickle"
        if db.split("datablock.json")[0] + "strong.pickle" in strong
        else None
        for db in datablocks]
    for i in range(len(timestamps)):
      ts, s, db = timestamps[i], strong[i], datablocks[i]
      if s is not None:
        process_pair(s, db, ts=ts)
  else:
    if len(strong) != len(datablocks):
      for pickle in strong:
        matching_json = pickle.split("_strong.pickle")[0] + "_datablock.json"
        if matching_json in datablocks:
          process_pair(pickle, matching_json)
    else: # Looks redundant but doesn't depend on file naming conventions
      for s, db in zip(strong, datablocks):
        process_pair(s, db)
  return strong_spot_timestamps, strong_spot_lengths, strong_spot_intensities
