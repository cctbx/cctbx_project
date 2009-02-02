from scitbx.array_family import flex
from libtbx import group_args, smart_open
import iotbx.phil
from mmtbx import polygon
import libtbx, os


def run():
  c1 = group_args(key = "dhigh",  value_min = 2,    value_max = 2.5,  target_value  = None)
  c2 = group_args(key = "natoms", value_min = 0,    value_max = 50000, target_value = None)
  c3 = group_args(key = "data",   value_min = None, value_max = None, target_value  = "fobs_x")
  filter_keys_and_targets = [c1,c2, c3]

  database_file = libtbx.env.find_in_repositories(
    relative_path = "chem_data/polygon_data/pdb_2009-01-14_ord.txt",
    test          = os.path.isfile)

  polygon.polygon(file_name               = database_file,
                  selected_keys           = ["rfact", "rfree", "angav", "bndav"],
                  filter_keys_and_targets = filter_keys_and_targets,
                  max_reject_fraction     = 0.1,
                  n_histogram_slots       = 10)


if (__name__ == "__main__"):
  run()
