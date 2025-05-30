"""Read a CBF file"""

from __future__ import absolute_import, division, print_function
import iotbx.cif
import iotbx.cif.validation
import libtbx.load_env
import sys, os
import six
op = os.path

def run(args):
  dic_path = libtbx.env.find_in_repositories(
    relative_path="cbflib/doc/cif_img_1.5.4_28Jul07.dic",
    test=op.isfile)
  if (dic_path is not None):
    cif_dic = iotbx.cif.validation.smart_load_dictionary(file_path=dic_path)
  else:
    cif_dic = None
  for file_name in args:
    raw = open(file_name, "rb").read()
    pattern = "--CIF-BINARY-FORMAT-SECTION--"
    i = raw.find(pattern)
    if (i >= 0):
      j = raw.rfind(pattern+"--")
    sliced = raw[:i] + "\n" + raw[j:]
    cif = iotbx.cif.fast_reader(input_string=sliced)
    model = cif.model()
    for k, v in six.iteritems(model):
      print(v["_array_data.data"])
    if (cif_dic is not None):
      model.validate(cif_dic)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
