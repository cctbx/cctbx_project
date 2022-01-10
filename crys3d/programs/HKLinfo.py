from __future__ import absolute_import, division, print_function

try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import os

from iotbx.gui_tools.reflections import ArrayInfo
from iotbx.reflection_file_reader import any_reflection_file


class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s reflectionfile
or
  %(prog)s reflectionfile philinput.txt

where reflectionfile can be any of the conventional reflection file formats
such as .mtz, .cif, .sca or .hkl file. This will print a table to the screen
listing properties of the reflection data arrays present in the file. The
name of the properties to be listed can be shown by typing:
%(prog)s --show-defaults
and noting which PHIL parameters in the scope \"selected_info\" are set to True.
These can be changed either by specifying these on the command line or by
entering the assignments into a text file, say \"philinput.txt\" that is
submitted on the command line together with the name of the reflection file.
""" % locals()

  datatypes = ['miller_array', 'phil' ]
  master_phil_str ="""
merge_equivalents = False
  .type = bool
  .caption = "merging symmetry equivalent reflections into unique wedge in reciprocal space"
  .short_caption = "merging reflections into a symmetry unique wedge"
reconstruct_amplitudes = False
  .type = bool
  .caption = "Convert mean amplitudes, F,SIGF, and anomalous differences, DANO,SIGDANO, into anomalous arrays, F(+),SIGF(+),F(-),SIGF(-)"
  .short_caption = "Turn mean amplitudes and anomalous differences into anomalous arrays"
""" + ArrayInfo.arrayinfo_phil_str

  def validate(self):
    pass

  def run(self):
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
    data_file = self.data_manager.get_miller_array_names()[0]
    hkl_file = any_reflection_file(data_file)
    arrays = hkl_file.as_miller_arrays(merge_equivalents=self.params.merge_equivalents,
                                       reconstruct_amplitudes=self.params.reconstruct_amplitudes)

    print("%d Miller arrays in this dataset:" %len(arrays))
    delimiter = self.params.delimiter
    array_info_format_tpl=[]
    for i,array in enumerate(arrays):
      arrayinfo = ArrayInfo(array,self.params.wrap_labels)
      info_fmt, headerstr, infostr = arrayinfo.get_selected_info_columns_from_phil(self.params)
      if i==0:
        print(headerstr)
      print(infostr)
