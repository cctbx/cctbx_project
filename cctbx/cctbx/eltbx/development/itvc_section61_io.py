import cctbx.eltbx.gaussian_fit
from cctbx.eltbx.gaussian_fit import international_tables_stols
from cctbx.array_family import flex
from scitbx.python_utils.misc import line_feeder
from scitbx.python_utils.misc import adopt_init_args

class table6111_entry:

  def __init__(self, element, atomic_number, method, table_y):
    adopt_init_args(self, locals())
    self.atomic_symbol = element
    if (element == "Cval"):
      self.atomic_symbol = "C"
    elif (element == "Sival"):
      self.atomic_symbol = "Si"
    else:
      for sign in ["+", "-"]:
        i = element.find(sign)
        if (i > 0):
          self.element = element.replace(sign,"") + sign
          self.atomic_symbol = element[:i]
          break

class table6111:

  def __init__(self, file_name):
    self.file_name = file_name
    self.entries = {}
    self.elements = []

  def enter_entry(self, entry):
    self.entries[entry.element] = entry
    self.elements.append(entry.element)

  def enter_block(self, elements, atomic_numbers, methods, value_rows):
    i_column = -1
    for element,atomic_number,method in zip(elements,atomic_numbers,methods):
      i_column += 1
      table_y = flex.double()
      for value_row in value_rows:
        table_y.append(value_row[i_column])
      self.enter_entry(table6111_entry(element,atomic_number,method,table_y))

def read_table6111(file_name):
  tab = table6111(file_name)
  lf = line_feeder(open(file_name))
  while 1:
    line = lf.next()
    if (lf.eof): break
    if (line.startswith("Element")):
      elements = line.split()[1:]
      line = lf.next()
      assert line.lstrip().startswith("Z"), line
      atomic_numbers = [int(z) for z in line.split()[1:]]
      assert len(atomic_numbers) == len(elements), line
      line = lf.next()
      assert line.startswith("Method"), line
      methods = line.split()[1:]
      assert len(methods) == len(elements), line
      line = lf.next()
      assert line.find("sin") > 0, line
      stols = flex.double()
      value_rows = []
      while 1:
        line = lf.next()
        assert not lf.eof
        if (len(line.strip()) == 0): continue
        raw_value_row = line.rstrip().split("\t")
        assert len(raw_value_row) == len(elements) + 1, line
        stols.append(float(raw_value_row[0]))
        assert stols[-1] == international_tables_stols[stols.size()-1], line
        value_row = []
        for value in raw_value_row[1:]:
          if (len(value.strip()) == 0):
            value_row.append(0)
          else:
            try: value_row.append(float(value))
            except ValueError, e: raise ValueError(line)
        value_rows.append(value_row)
        if (stols.size() == 62):
          assert stols[-1] == 6, line
          break
      tab.enter_block(elements, atomic_numbers, methods, value_rows)
  return tab
