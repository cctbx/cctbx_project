from __future__ import division
from cctbx import sgtbx

class polar_manager(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''


    def get_polar_and_rotxy(self, selected_pickle_file, polar_info):
      file_op = open(polar_info,'r')
      result_op = file_op.read()
      pickle_info_data=result_op.split("\n")

      polar_hkl = 'None'
      rotx = 0
      roty = 0
      for line in pickle_info_data:
        if line.find(selected_pickle_file) > 0:
          data = line.split(' ')
          polar_hkl = data[3]
          k = float(data[4])
          G = float(data[5])
          rotx = float(data[6])
          roty = float(data[7])

          break

      return polar_hkl, k, G, rotx, roty

    def get_pickle_filename(self, selected_frame_id, polar_info):
      file_op = open(polar_info,'r')
      result_op = file_op.read()
      pickle_info_data=result_op.split("\n")

      pickle_filename = ''
      for line in pickle_info_data:
        line_data = line.split(' ')
        if line_data[1] == str(selected_frame_id):
          pickle_filename = line_data[0]
          break

      return pickle_filename

    def get_pickle_filename_and_polar(self, selected_frame_id, polar_info):
      file_op = open(polar_info,'r')
      result_op = file_op.read()
      pickle_info_data=result_op.split("\n")

      pickle_filename = ''
      polar_hkl = ''
      for line in pickle_info_data:
        line_data = line.split(' ')
        if line_data[1] == str(selected_frame_id):
          pickle_filename = line_data[0]
          polar_hkl = line_data[2]
          break

      return pickle_filename, polar_hkl

    def convert_miller_polar(self, miller_array, polar_hkl):
      cb_op = sgtbx.change_of_basis_op(polar_hkl)
      new_miller_array = None
      if miller_array is not None:
        new_miller_array = miller_array.change_basis(cb_op)

      miller_index = miller_array.indices()
      new_miller_index = new_miller_array.indices()
      flex_iobs=miller_array.data()
      flex_new_iobs=new_miller_array.data()

      """
      print "Change basis for this pickle. Top 10 miller indices:"
      for i in range(10):
        print miller_index[i], flex_iobs[i], new_miller_index[i], flex_new_iobs[i]
      """

      return new_miller_array
