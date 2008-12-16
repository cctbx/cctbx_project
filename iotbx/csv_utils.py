import csv

class writer(object):
  def __init__(self,
               file_object,
               fields,
               field_names=None,
               delimiter=','):
    for field in fields:
      assert len(field) == len(fields[0])
    iter_object = self._iter_rows(fields)
    writer = csv.writer(file_object, iter_object, delimiter=delimiter)
    if field_names:
      writer.writerow(field_names)
    writer.writerows(iter_object)

  def _iter_rows(self, args):
    for i in range(len(args[0])):
      yield tuple([arg[i] for arg in args])

class reader(object):
  def __init__(self,
               file_object,
               data_type=str,
               field_names=False,
               delimiter=','):
    """data_type can be a single type which is then applied to all fields,
    or a list of types which apply to each field.
    """
    reader = csv.reader(file_object, delimiter=delimiter)
    data = []

    for row in reader:
      if reader.line_num == 1:
        n_data = len(row)
        if type(data_type) in (list,tuple):
          assert len(data_type) == n_data
          data_type_list = data_type
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
