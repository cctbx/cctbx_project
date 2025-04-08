"""Tools for manipulating ranges
"""
from __future__ import absolute_import, division, print_function
from six.moves import range
from six.moves import zip
def range_parser( txt ):
  splitter = " "
  if "," in txt:
    splitter = ","
  txt_ranges = txt.split(splitter)
  ranges = []
  for range_ in txt_ranges:
    tmp = range_.split("-")
    assert len(tmp)<=2
    if len(tmp)==2:
      ranges.append( [int(tmp[0]),int(tmp[1])]  )
    if len( tmp) == 1:
      ranges.append( [ int(tmp[0]) ] )
  return ranges

def range_to_list(range_):
  result = []
  for item in range_:
    if len(item)==1:
      result.append( item[0] )
    if len(item)==2:
      for ii in range(item[0],item[1]+1):
        result.append( ii )
  return result



class serial_file_name_handler(object):
  def __init__(self, base_name, wildcard="#"):
    self.user_base_name = base_name
    self.wildcard = wildcard
    self.base, self.extension, self.n = self.process_user_input()

  def process_user_input(self):
    at_serial = False
    at_extension = False
    p1=""
    p3=""
    n=0
    for ii in self.user_base_name:
      if not at_serial:
        if (ii == self.wildcard):
          at_serial=True
        if not at_extension:
          if not at_serial:
            p1+=ii
      if at_serial:
        n += 1
        if (ii!=self.wildcard):
          at_serial=False
          at_extension=True
      if at_extension:
        p3 += ii
    if len(p3)==0:
      n = n
    else:
      n = n-1
    return p1,p3,n


  def name_from_number(self, serial_id=1):
    id = str(serial_id)
    return self.base+(self.n-len(id))*"0"+id+self.extension

  def names_from_range_list(self, range_list):
    result = []
    for id in range_list:
      result.append( self.name_from_number( id ) )
    return result

  def names_from_range(self, range_txt):
    result = []
    range_list = range_to_list( range_parser( range_txt ) )
    result = self.names_from_range_list( range_list )
    return result



def tst_file_names():
  base="lysozyme_1_####.img"
  obj = serial_file_name_handler(base)
  assert obj.base=="lysozyme_1_"
  assert obj.extension == ".img"
  assert obj.n == 4
  assert obj.name_from_number( 23 ) == "lysozyme_1_0023.img"

  base="lysozyme_1.####"
  obj = serial_file_name_handler(base)
  assert obj.base=="lysozyme_1."
  assert obj.extension == ""
  assert obj.n == 4
  assert obj.name_from_number( 23 ) == "lysozyme_1.0023"

  base = "######.lysozyme_1.img"
  obj = serial_file_name_handler(base)
  assert obj.base==""
  assert obj.extension == ".lysozyme_1.img"
  assert obj.n == 6
  assert obj.name_from_number( 23 ) == "000023.lysozyme_1.img"
  lst = obj.names_from_range_list( [1,2,3] )
  assert lst[0] == "000001.lysozyme_1.img"
  assert lst[1] == "000002.lysozyme_1.img"
  assert lst[2] == "000003.lysozyme_1.img"
  assert len(lst)==3


  base = "######.lysozyme_1.img"
  obj = serial_file_name_handler(base)
  range_txt = "1-3"
  name_list =  obj.names_from_range(range_txt)
  assert name_list[0] == '000001.lysozyme_1.img'
  assert name_list[1] == '000002.lysozyme_1.img'
  assert name_list[2] == '000003.lysozyme_1.img'

def tst_ranges():
  range_txt="1-9,10,11,13-16"
  result = range_to_list( range_parser( range_txt ) )
  tst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16]
  for ii, jj in zip(result,tst):
    assert(ii==jj)
  range_txt="1-9 10 11 13-16"
  result = range_to_list( range_parser( range_txt ) )
  tst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16]
  for ii, jj in zip(result,tst):
    assert(ii==jj)

if __name__ == "__main__":
  tst_ranges()
  tst_file_names()
  print("OK")
