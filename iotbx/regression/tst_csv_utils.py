from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx.test_utils import open_tmp_file, Exception_expected
from six.moves import zip

def exercise():
  exercise_writer()
  exercise_reader()

def exercise_writer():
  from iotbx import csv_utils

  x = (1,2,3,4,5)
  y = (6,7,8,9,10)
  f = open_tmp_file()
  field_names = ('x','y')
  csv_utils.writer(f, (x,y), field_names=field_names)
  f.close()
  f = open(f.name, 'r')
  content = [l.strip() for l in f.readlines()]

  text = ['x,y']
  text += ['%s,%s' %(row[0],row[1]) for row in zip(x,y)]
  assert content == text
  f.close()

  x = (1,2,3,4,5)
  y = (6,7,8,9,10)
  f = open_tmp_file()
  csv_utils.writer(f, (x,y), delimiter=';')
  f.close()
  f = open(f.name, 'r')
  content = [l.strip() for l in f.readlines()]
  text = ['%s;%s' %(row[0],row[1]) for row in zip(x,y)]
  assert content == text
  f.close()

  x = flex.int(x)
  y = flex.int(y)
  f = open_tmp_file()
  csv_utils.writer(f, (x,y), field_names=field_names)
  f.close()
  f = open(f.name, 'r')
  content = [l.strip() for l in f.readlines()]
  text = ['x,y']
  text += ['%s,%s' %(row[0],row[1]) for row in zip(x,y)]
  assert content == text
  f.close()

  y.append(11)
  f = open_tmp_file()
  try:
    csv_utils.writer(f, (x,y), field_names=field_names)
  except AssertionError:
    pass
  else:
    raise Exception_expected
  f.close()


def exercise_reader():
  from iotbx import csv_utils

  x = (1,2,3,4,5)
  y = (6,7,8,9,10)
  f = open_tmp_file()
  field_names = ('x','y')
  csv_utils.writer(f, (x,y), field_names=field_names,delimiter=';')
  f.close()
  f = open(f.name, 'r')
  a = csv_utils.reader(f, data_type=int, field_names=True,delimiter=';')
  f.close()
  assert tuple(a.data[0]) == x
  assert tuple(a.data[1]) == y

  x = (1,2,3,4,5)
  y = (1.1,2.2,3.3,4.4,5.5)
  f = open_tmp_file()
  csv_utils.writer(f, (x,y))
  f.close()
  f = open(f.name, 'r')
  data_type_list = (int, float)
  a = csv_utils.reader(f, data_type_list=data_type_list)
  f.close()
  assert tuple(a.data[0]) == x
  assert tuple(a.data[1]) == y

  f = open(f.name, 'r')
  data_type_list = (int, float)
  try:
    a = csv_utils.reader(f, data_type=int,
                         data_type_list=data_type_list)
    # Can't pass data_type AND data_type_list
  except AssertionError:
    pass
  else:
    raise Exception_expected
  f.close()

  f = open(f.name, 'r')
  a = csv_utils.reader(f)
  f.close()
  assert list(a.data[0]) == [str(i) for i in x]
  assert list(a.data[1]) == [str(i) for i in y]

def run():
  exercise()
  print("OK")

if __name__ == '__main__':
  run()
