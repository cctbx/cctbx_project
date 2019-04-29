from __future__ import division, print_function

from datetime import datetime

class perf_handler():

  def __init__(self, module_name):
    self.module_name = module_name
    self.t_st = datetime.now()

  def get_elapsed_times(self):
    self.t_en = datetime.now()
    self.t_spent =  self.t_en - self.t_st
    print(self.module_name, self.t_st.strftime("%d/%m/%Y %H:%M:%S"), self.t_en.strftime("%d/%m/%Y %H:%M:%S"), self.t_spent.microseconds)
