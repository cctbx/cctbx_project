from __future__ import absolute_import, division, print_function
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.find_highres_shots
#
import os, glob
import libtbx.load_env
from dials.util import show_mail_on_error
from dials.util.options import ArgumentParser
from libtbx.phil import parse
from libtbx import easy_pickle
from scitbx.array_family import flex
import six

help_message = """
Script to find high resolution shots from XFEL data
Example usage:
cxi.find_highres_shots data=r000[6-9]/001/integration data=r001[1-4]/001/integration \
  min_good_reflections=100 min_IsigI=10

This looks at each integration pickle and selects the strong reflections (I/sigI >= 10).  Then,
we throw out all images with less than 100 strong reflections.  For the images that are left,
we find the 10 images whose strong reflections diffract to the highest resolution.  The
command prints out a cctbx.image_viewer command you can directly run to see these images.
"""

phil_str = '''
  data = None
    .type = str
    .multiple = True
    .help = Paths to integration dirs
  max_images = 10
    .type = int
    .help = Find at most this many images
  min_resolution = None
    .type = float
    .help = Mininmum resolution to accept
  min_IsigI = 1.0
    .type = float
    .help = Minimum accepted value for I/sigI for determining number of good reflections
  min_good_reflections = None
    .type = int
    .help = Only accept images with at least this many measurements with I/sigI >= min_IsigI
  max_to_examine = None
    .type = int
    .help = Only examine this many images
'''

phil_scope = parse(phil_str)

class Script(object):
  """ Script to find good images """
  def __init__(self):
    """ Set up the option parser. Arguments come from the command line or a phil file """
    self.usage = "%s data=path"%(libtbx.env.dispatcher_name)
    self.parser = ArgumentParser(
      usage = self.usage,
      phil = phil_scope,
      epilog = help_message)

  def run(self):
    """ Find the high res images """

    params, options = self.parser.parse_args(
      show_diff_phil=True)

    assert params.max_images is not None
    assert params.min_IsigI is not None

    final_data = {}
    worst_resolution = None
    worst_image = None
    n_tested = 0

    all_paths = []
    for item in params.data:
      all_paths.extend(glob.glob(item))

    # loop through the integration pickles
    for path in all_paths:
      for filename in os.listdir(path):
        if not os.path.splitext(filename)[1] == ".pickle":
          continue

        filepath = os.path.join(path, filename)
        try:
          data = easy_pickle.load(filepath)
        except Exception as e:
          print("Couldn't read", filepath)
          continue

        if not 'observations' in data:
          print("Not an integration pickle", filepath)
          continue

        n_tested += 1

        # Find the critical parameters for this image
        obs = data['observations'][0]
        uc = obs.unit_cell()
        obs_strong = obs.select((obs.data()/obs.sigmas() >= params.min_IsigI))
        strong_measurements = len(obs_strong.indices())
        d = uc.d(obs_strong.indices())
        resolution = flex.min(d)

        if params.min_resolution is not None and resolution > min_resolution:
          print("Rejecting", filepath, "because resolution", resolution, "is too low")
          continue

        if params.min_good_reflections is not None and strong_measurements < params.min_good_reflections:
          print("Rejecting", filepath, "because number of measurements where I/sigmaI >=", params.min_IsigI, strong_measurements, "is too low")
          continue

        # Save the image if it's good enough
        if len(final_data) < params.max_images:
          print("Accepting", filepath, "not yet reached maximum number good images")
          final_data[filepath] = data, resolution
          if worst_resolution is None or resolution > worst_resolution:
            worst_resolution = resolution
            worst_image = filepath
        elif resolution < worst_resolution:
          print("Accepting", filepath, "better than the worst image")
          final_data.pop(worst_image)
          final_data[filepath] = data, resolution

          worst_resolution = 0
          for key, entry in six.iteritems(final_data):
            data, resolution = entry
            if resolution > worst_resolution:
              worst_resolution = resolution
              worst_image = key
        else:
          print("Rejecting", filepath)

        if params.max_to_examine is not None and n_tested > params.max_to_examine: break
      if params.max_to_examine is not None and n_tested > params.max_to_examine: break

    # Results
    print("Final best images and resolutions:")
    for key in final_data:
      print(key, final_data[key][1])

    # Show a viewing command
    cmd = "cctbx.image_viewer "

    for key in final_data:
      dirname = os.path.dirname(key).replace('integration', 'out')
      basename = os.path.basename(key)

      # tsl: time stamp long, tss: time stamp short
      if basename.endswith("_00000.pickle"):
        tsl = basename.lstrip('int-').split('_00000.pickle')[0]
      elif basename.endswith(".pickle"):
        tsl = os.path.splitext(basename.lstrip('int-'))[0]

      tss = ''.join([c for c in tsl if c in '1234567890'])

      cmd += (os.path.join(dirname, "idx-%s.pickle "%tss))

    print("Use this command to display these images:")
    print(cmd.rstrip())

if __name__ == "__main__":
  with show_mail_on_error():
    script = Script()
    script.run()
