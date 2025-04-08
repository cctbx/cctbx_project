"""
Tools to read and write csv formatted files
"""

from __future__ import absolute_import, division, print_function
import csv
from six.moves import range

if 'unix' not in csv.list_dialects():
  csv.register_dialect('unix', lineterminator='\n', quoting=csv.QUOTE_NONE)

class writer(object):
  def __init__(self,
               file_object,
               fields,
               field_names=None,
               delimiter=','):
    for field in fields:
      assert len(field) == len(fields[0])
    iter_object = self._iter_rows(fields)
    writer = csv.writer(file_object, delimiter=delimiter,
                        dialect='unix', quoting=csv.QUOTE_NONE)
    if field_names:
      writer.writerow(field_names)
    writer.writerows(iter_object)

  def _iter_rows(self, args):
    for i in range(len(args[0])):
      yield tuple([arg[i] for arg in args])

class reader(object):
  def __init__(self,
               file_object,
               data_type=None,
               data_type_list=None,
               field_names=False,
               delimiter=','):
    """Supply a single data_type which is to be applied to all fields,
    or supply a data_type_list containing the type of each field.
    String is the default type if none is supplied.
    """
    assert data_type is None or data_type_list is None
    if data_type is None and data_type_list is None:
      data_type = str

    reader = csv.reader(file_object, delimiter=delimiter)
    data = []

    for i_row,row in enumerate(reader):
      if i_row == 0:
        n_data = len(row)
        if data_type_list is not None:
          assert len(data_type_list) == n_data
        else:
          data_type_list = [data_type] * n_data
        for i in range(n_data):
          if field_names:
            data.append([])
          else:
            data.append([data_type_list[i](row[i])])
        n_data = len(data)
      else:
        for i in range(n_data):
          data[i].append(data_type_list[i](row[i]))
    self.data = data
