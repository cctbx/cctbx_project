"""
the parameters should have this scope
  omit {
    perform_omit = True
    fraction = 0.15
    max_number = 1e5
    number_of_sets = 100
    root_name = 'omit_'
  }
"""
from __future__ import absolute_import, division, print_function
from six.moves import range

class random_omit_data(object):
  def __init__(self,
               miller_array,
               parameters):
    self.miller_array = miller_array
    self.parameters = parameters

  def write_datasets(self):

    comfile = open(self.parameters.root_name+"exec.com", "w")

    for nth in range(self.parameters.number_of_sets):
      file_name = self.parameters.root_name + str(nth)+".mtz"
      print("__REPLACE_1__%s__REPLACE_2__ > %s.log"%(file_name,nth), file=comfile)
      tmp_select = self.miller_array.generate_r_free_flags(
        fraction = self.parameters.fraction ,
        max_free =  self.parameters.max_number,
        use_lattice_symmetry = False)
      tmp_miller = self.miller_array.select( ~tmp_select.data() )
      tmp_mtz_dataset = tmp_miller.as_mtz_dataset(
        column_root_label="Fdelta")
      tmp_mtz_dataset.mtz_object().write(
        file_name=file_name)
