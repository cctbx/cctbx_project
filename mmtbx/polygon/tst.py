from scitbx.array_family import flex
import iotbx.phil
from mmtbx import polygon
import libtbx, os

params1 = iotbx.phil.parse("""\
polygon {
  keys_to_show = *rfact *rfree *bndav *angav *adpav
  filter
  {
    key = *dhigh
    value_min = 2
    value_max = 2.5
  }
  filter
  {
    key = *natoms
    value_min = 0
    value_max = 50000
  }
  filter
  {
    key = *data
    target_value = fobs_x
  }
}
""")

params2 = iotbx.phil.parse("""\
polygon {
  keys_to_show = *rfact *rfree *bndav *angav *adpav
}
""")

def example_1():
  # show selected characteristics, apply pre-defined filters
  pr, unused_definitions = polygon.master_params.fetch(sources = [params1],
    track_unused_definitions = True)
  polygon.polygon(params = pr.extract())

def example_2():
  # show selected characteristics for the whole database
  pr, unused_definitions = polygon.master_params.fetch(sources = [params2],
    track_unused_definitions = True)
  polygon.polygon(params = pr.extract())

def example_3():
  # show selected characteristics, apply default selection
  # d_min comes from the structure that we are comparing against database
  pr, unused_definitions = polygon.master_params.fetch(sources = [params2],
    track_unused_definitions = True)
  polygon.polygon(params = pr.extract(),
                  d_min  = 2.0)


if (__name__ == "__main__"):
  file_name = libtbx.env.find_in_repositories(
      relative_path = "chem_data/polygon_data/pdb_2009-01-14_ord.txt",
      test = os.path.isfile)
  if(file_name is None):
    print "Skip POLYGON test: database file is not available."
    return
  print "\nEXAMPLE 1:"
  example_1()
  print "\nEXAMPLE 2:"
  example_2()
  print "\nEXAMPLE 3:"
  example_3()
