from  cctbx.array_family import flex
import sys
import os
import string

class plot_data:
  def __init__(self,
               plot_title,
               x_label,
               y_label,
               x_data,
               y_data,
               y_legend,
               comments):
    self.plot_title = plot_title

    self.x_label = x_label
    self.y_label = y_label
    self.x_data = x_data
    self.comments = comments
    self.domain_flag = 'N'

    ## The x_data is an (flex) array of numbers
    ## The Y_data should be an list of (flex) arrays
    ## The legends should be an array of strings
    self.y_legend = []
    self.y_data = []
    if y_data is not None:
      assert ( len(y_data)==len(self.x_data) )
      self.y_data.append(y_data)
      self.y_legend.append(y_legend)
      if flex.min(y_data) < 0:
        self.domain_flag = 'A'

  def add_data(self,
               y_data,
               y_legend):
    assert ( len(y_data)==len(self.x_data) )
    self.y_data.append(y_data)
    self.y_legend.append(y_legend)
    if self.domain_flag =='N':
      if flex.min(y_data) < 0:
        self.domain_flag = 'A'


def plot_data_loggraph(plot_data,output):
  ## First we need to print the header information
  output.write('$TABLE: %s:\n'%(plot_data.plot_title) )
  output.write('$GRAPHS\n')
  output.write(':%s \n' %(plot_data.comments))
  index_string = ''
  for ii in range(len(plot_data.y_data)+1):
    index_string += '%d,'%(ii+1)
  output.write(':%s:%s:\n'%(plot_data.domain_flag,index_string[:-1]))
  output.write('$$\n')
  ## replace spaces for loggraph with underscores
  tmp_legend = plot_data.x_label
  spaces = 0
  spaces = string.find(tmp_legend,' ')
  if spaces>0:
    tmp_legend = tmp_legend.replace(' ','_')
  label_string = '%s'%(tmp_legend)
  for ii in range(len(plot_data.y_data)):
    tmp_legend = plot_data.y_legend[ii]
    ## loggraph does not like spaces in the legend names
    ## lets replace them with underscores
    spaces = 0
    spaces = string.find(tmp_legend,' ')
    if spaces>0:
      tmp_legend = tmp_legend.replace(' ','_')
    label_string += '   %s'%( tmp_legend )
  output.write('%s   $$ \n'%(label_string))
  output.write('$$\n')
  for ii in range(len(plot_data.x_data)):
    data_string = '%f'%(plot_data.x_data[ii])
    for jj in range(len(plot_data.y_data)):
      data_string +='   %f'%(plot_data.y_data[jj][ii])
    output.write('%s\n'%(data_string))
  output.write('$$\n')
