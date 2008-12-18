from scitbx.array_family import flex
from libtbx.test_utils import Exception_expected
import tempfile
import sys

def exercise():
  v = sys.version_info
  if (v[0] == 2 and v[1] <= 2):
    print "Skipping iotbx.csv_utils tests: csv extension not available"
    return
  exercise_writer()
  exercise_reader()

def exercise_writer():
  from iotbx import csv_utils

  x = (1,2,3,4,5)
  y = (6,7,8,9,10)
  filename = tempfile.mktemp()
  f = open(filename, 'w')
  field_names = ('x','y')
  csv_utils.writer(f, (x,y), field_names=field_names)
  f.close()
  f = open(filename, 'r')
  content = f.readlines()
  text = ['x,y\r\n']
  text += ['%s,%s\r\n' %(row[0],row[1]) for row in zip(x,y)]
  assert content == text
  f.close()

  x = (1,2,3,4,5)
  y = (6,7,8,9,10)
  filename = tempfile.mktemp()
  f = open(filename, 'w')
  csv_utils.writer(f, (x,y), delimiter=';')
  f.close()
  f = open(filename, 'r')
  content = f.readlines()
  text = ['%s;%s\r\n' %(row[0],row[1]) for row in zip(x,y)]
  assert content == text

  x = flex.int(x)
  y = flex.int(y)
  f = open(filename, 'w')
  csv_utils.writer(f, (x,y), field_names=field_names)
  f.close()
  f = open(filename, 'r')
  content = f.readlines()
  text = ['x,y\r\n']
  text += ['%s,%s\r\n' %(row[0],row[1]) for row in zip(x,y)]
  assert content == text
  f.close()

  y.append(11)
  f = open(filename, 'w')
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
  filename = tempfile.mktemp()
  f = open(filename, 'w')
  field_names = ('x','y')
  csv_utils.writer(f, (x,y), field_names=field_names,delimiter=';')
  f.close()
  f = open(filename, 'r')
  a = csv_utils.reader(f, data_type=int, field_names=True,delimiter=';')
  f.close()
  assert tuple(a.data[0]) == x
  assert tuple(a.data[1]) == y

  x = (1,2,3,4,5)
  y = (1.1,2.2,3.3,4.4,5.5)
  filename = tempfile.mktemp()
  f = open(filename, 'w')
  csv_utils.writer(f, (x,y))
  f.close()
  f = open(filename, 'r')
  data_type_list = (int, float)
  a = csv_utils.reader(f, data_type_list=data_type_list)
  f.close()
  assert tuple(a.data[0]) == x
  assert tuple(a.data[1]) == y

  f = open(filename, 'r')
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

  f = open(filename, 'r')
  a = csv_utils.reader(f)
  f.close()
  assert list(a.data[0]) == [str(i) for i in x]
  assert list(a.data[1]) == [str(i) for i in y]

def run():
  exercise()
  print "OK"

if __name__ == '__main__':
  run()
