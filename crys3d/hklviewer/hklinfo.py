
from __future__ import absolute_import, division, print_function

import sys
from iotbx.gui_tools.reflections import ArrayInfo
from iotbx.reflection_file_reader import any_reflection_file

def run(data_file, philparams, log = sys.stdout):
  """
  Print a table of properties of miller_array objects from a reflection file to stdout formatted like

7 Miller arrays in this dataset:
 Labels        |       Type      |  #HKLs  |     min,max data       |Anomalous|Sym.uniq.|Data compl.|
  R-free-flags |     R-free flag |   30451 |        0.0,         1.0|   False |    True |    0.9999 |
  FOBS,SIGFOBS |       Amplitude |   30443 |     3.9975,      1557.4|   False |    True |   0.99964 |
  IOBS,SIGIOBS |       Intensity |   30451 |    -484.25,  2.1485e+05|   False |    True |    0.9999 |
 I(+),SIGI(+),
  I(-),SIGI(-) |       Intensity |   57347 |    -727.01,  2.1485e+05|    True |    True |    0.9984 |
 F(+),SIGF(+),
  F(-),SIGF(-) |       Amplitude |   57338 |     3.9975,      1557.4|    True |    True |   0.99824 |
          DANO |  Floating-point |   30451 |    -102.58,      194.38|   False |    True |    0.9999 |
         SIGDP |       Amplitude |   30450 |        0.0,      85.081|   False |    True |   0.99987 |

  """

  hkl_file = any_reflection_file(data_file)
  arrays = hkl_file.as_miller_arrays(merge_equivalents=philparams.merge_equivalents)

  print("%d Miller arrays in this dataset:" %len(arrays))
  delimiter = philparams.delimiter
  array_info_format_tpl=[]
  for i,array in enumerate(arrays):
    arrayinfo = ArrayInfo(array,philparams.wrap_labels)
    info_fmt, headerstr, infostr = arrayinfo.get_selected_info_columns_from_phil(philparams)
    if i==0:
      print(headerstr)
    print(infostr)

