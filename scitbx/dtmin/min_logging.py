from __future__ import print_function
from __future__ import division

DEF_TAB = 3

SILENCE = 0 #
LOGFILE = 1 #
VERBOSE = 2 #
FURTHER = 3 #
TESTING = 4 #

class Logging:
  def initial_statistics(self): pass

  def current_statistics(self): pass

  def final_statistics(self): pass

  def __init__(self):
    self.output_level = 0 #default all output

  ### printing methods

  def log_output(self, depth, output_level, string, add_return):
    if output_level <= self.output_level:
      print(depth*DEF_TAB*' ' + string, end='\n' if add_return else '')

  def log_tab_printf(self, t, output_level, text, fformat): #format is already taken
    if output_level <= self.output_level:
      print((t*DEF_TAB*' ' + text) % fformat, end = "")

  def log_blank(self, output_level):
    self.log_output(0,output_level,"",True)

  def log_underline(self, output_level, string):
    self.log_blank(output_level)
    self.log_tab(1,output_level, string)
    self.log_line(1,output_level, len(string),'-')

  def log_line(self, t, output_level, len, c):
    self.log_tab(t, output_level,len*c)

  def log_protocol(self, output_level, macrocycle_protocol):
    self.log_underline(output_level, "Protocol:")
    if (macrocycle_protocol[0] == "off"):
      self.log_tab(1,output_level,"No parameters to refine this macrocycle")
      self.log_blank(output_level)
    else:
      for i in range(len(macrocycle_protocol)):
        self.log_tab(2,output_level,macrocycle_protocol[i].upper() + " ON")
      self.log_blank(output_level)

  def log_tab(self, t, output_level, text, add_return=True):
    self.log_output(t, output_level, text, add_return)

  def log_ellipsis_start(self, output_level, text):
    self.log_blank(output_level)
    self.log_output(1,output_level,text+"...",True)

  def log_ellipsis_end(self, output_level):
    self.log_output(2,output_level,"...Done",True)
    self.log_blank(output_level)

  def log_parameters(self, output_level, macrocycle_parameter_names):
    self.log_tab(1,output_level, "Parameters:")
    for i in range(len(macrocycle_parameter_names)):
      self.log_tab(1,output_level, "Refined Parameter #:" + str(i+1) + " " + macrocycle_parameter_names[i])

  def log_hessian(self, output_level, hessian, macrocycle_parameter_names):
    max_dim = 100 #very generous but not infinite
    (n_rows, n_cols) = hessian.all()
    assert(n_rows == n_cols)
    self.log_tab(1,output_level,"Matrix")
    if n_rows > max_dim:
      self.log_tab(1,output_level,"Matrix is too large to write to log: size = "+str(n_rows))
    for i in range(min(n_rows, max_dim-1)):
      hess = ""
      line = ""
      for j in range(min(n_cols, max_dim-1)):
        hess = ("%+3.1e " % hessian[i,j])
        if hessian[i,j] == 0:
          line += " ---0--- "
        else:
          line += hess
      line = line[:-1] # cut the last space off for neatness
      if n_rows >= max_dim:
        self.log_tab_printf(1,output_level,"%3d [%s etc...]\n",(i+1,line))
      else:
        self.log_tab_printf(1,output_level,"%3d [%s]\n",(i+1,line))
    if n_rows >= max_dim:
      self.log_tab_printf(0,output_level," etc...\n",())
    self.log_blank(output_level)
    self.log_tab(1,output_level,"Matrix Diagonals")
    for i in range(n_rows):
      self.log_tab_printf(1,output_level,"%3d [% 9.4g] %s\n",
          (i+1,hessian[i,i],
          macrocycle_parameter_names[i] if (i < len(macrocycle_parameter_names)) else ""))
    self.log_blank(output_level)

  def log_vector(self, output_level, what, vec, macrocycle_parameter_names):
    self.log_tab(1,output_level,what)
    for i in range(len(vec)):
      self.log_tab_printf(1,output_level,"%3d [% 9.4g] %s\n",
          (i+1,vec[i],
          macrocycle_parameter_names[i] if (i < len(macrocycle_parameter_names)) else ""))
    self.log_blank(output_level)
