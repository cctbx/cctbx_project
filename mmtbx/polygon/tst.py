from __future__ import absolute_import, division, print_function
import iotbx.phil
from mmtbx import polygon
import libtbx, os
import libtbx.load_env

params1 = iotbx.phil.parse("""\
polygon {
  keys_to_show = *r_work *r_free *bond_rmsd *angle_rmsd *adp_mean_all
}
""")

params3 = iotbx.phil.parse("""\
polygon {
  keys_to_show = %s
}
""" % " ".join([ "*%s" % k for k in polygon.keys_to_show ]))

def example_1():
  # show selected characteristics, apply pre-defined filters
  pr, unused_definitions = polygon.master_params.fetch(sources = [params1],
    track_unused_definitions = True)
  polygon.polygon(params = pr.extract())

def example_2():
  #
  #pr, unused_definitions = polygon.master_params.fetch(sources = [params2],
  #  track_unused_definitions = True, d_min=1.0)
  polygon.polygon(d_min=1.0)

def example_3():
  # show selected characteristics, apply default selection
  # d_min comes from the structure that we are comparing against database
  pr, unused_definitions = polygon.master_params.fetch(sources = [params1],
    track_unused_definitions = True)
  polygon.polygon(params = pr.extract(),
                  d_min  = 2.0)

if (__name__ == "__main__"):
  file_name = libtbx.env.find_in_repositories(
      relative_path = "chem_data/polygon_data/all_mvd.pickle",
      test = os.path.isfile)
  if(file_name is None):
    print("Skip POLYGON test: database file is not available.")
  else:
    print("\nEXAMPLE 1:")
    example_1()
    print("\nEXAMPLE 2:")
    example_2()
    print("\nEXAMPLE 3:")
    example_3()
