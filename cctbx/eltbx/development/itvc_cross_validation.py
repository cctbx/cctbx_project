from __future__ import absolute_import, division, print_function
from cctbx.eltbx.development import itvc_section61_io
from cctbx.eltbx.gaussian_fit import international_tables_stols
from libtbx.option_parser import OptionParser
from six.moves import zip

def show_differences_if_any(label, stols, y0, y1):
  if (not y0.all_eq(y1)):
    print("Different:", label)
    for stol,y0i,y1i in zip(stols,y0,y1):
      if (y0i != y1i):
        print(stol, y0i, y1i)

def run(file_names):
  assert len(file_names) == 2
  tabs = []
  for file_name in file_names:
    tabs.append(itvc_section61_io.read_table6111(file_name))
  all_labels = {}
  for tab in tabs:
    for label in tab.elements: all_labels[label] = 1
  for label in all_labels:
    e0 = tabs[0].entries.get(label, None)
    e1 = tabs[1].entries.get(label, None)
    if ([e0,e1].count(None) == 0):
      assert e0.atomic_number == e1.atomic_number, \
        (label, e0.atomic_number, e1.atomic_number)
      assert e0.method == e1.method, \
        (label, e0.method, e1.method)
      min_size = min(e0.table_y.size(), e1.table_y.size())
      show_differences_if_any(
        label=label,
        y0=e0.table_y[:min_size],
        y1=e1.table_y[:min_size],
        stols=international_tables_stols[:min_size])
    for tab in tabs:
      ei = tab.entries.get(label, None)
      if (ei is not None and ei.element != ei.atomic_symbol):
        ee = tab.entries.get(ei.atomic_symbol, None)
        if ([ee,ei].count(None) == 0):
          assert ee.table_y.size() == 62
          assert ei.table_y.size() == 62
          show_differences_if_any(
            label=label,
            y0=ee.table_y[-6:],
            y1=ei.table_y[-6:],
            stols=international_tables_stols[-6:])

def main():
  parser = OptionParser(
    usage="usage: python %prog [options] file_name_1 file_name_2")
  (options, args) = parser.parse_args()
  if (len(args) != 2):
    parser.print_help()
    return
  run(file_names=args)

if (__name__ == "__main__"):
  main()
