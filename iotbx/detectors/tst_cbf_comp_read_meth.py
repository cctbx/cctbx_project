from __future__ import absolute_import, division, print_function
# test comparative read methods for cbf library
# default method provided by CBFlib-0.8.1; second method is optimized read
import os
from libtbx.development.timers import Profiler
cases = ["ana/myo2_3_00001.cbf",
"insulin2009/insu_before_1_00001.cbf",
"insulin2009/insu_before_1_00262.cbf",
"insulin2009/insu_before_1_00093.cbf",
"insulin2009/insu_before_1_00112.cbf",
"insulin2009/insu_before_1_00200.cbf",
"pilatus.web.psi.ch_insulin_0.2/DATA/DATASETS/insulin_0.2/run2_1_00052.cbf",
"pilatus.web.psi.ch_insulin_0.2/DATA/DATASETS/insulin_0.2/run2_1_00010.cbf",
"pilatus.web.psi.ch_insulin_0.2/DATA/DATASETS/insulin_0.2/run2_1_00149.cbf",
"pilatus.web.psi.ch_insulin_0.2/DATA/DATASETS/insulin_0.2/run2_1_00167.cbf",
"ribosome/images/colD55A_13_1_00172.cbf",
"ribosome/images/colD55A_13_2_00276.cbf",
"ribosome/images/colD55A_13_1_00227.cbf",
"ribosome/images/colD55A_13_2_00071.cbf",
"ribosome/images/colD55A_13_1_00336.cbf",
"ribosome/images/colD55A_13_2_00225.cbf",
"ribosome/images/colD55A_13_1_00022.cbf",
"ribosome/images/colD55A_13_2_00085.cbf",
"ribosome/G817_1_00002.cbf",
"ribosome/G817_2_00001.cbf",
"pilatus_2M_graeme/insulin_1_fixed_0043.img",
"pilatus_2M_graeme/insulin_1_fixed_0289.img",
]

def generate_paths(dirpath):
  for item in cases:
    file = os.path.join(dirpath,item)
    yield file

def test_all(dirpath,timer=False):
  for file in generate_paths(dirpath):
    from iotbx.detectors.pilatus_minicbf import PilatusImage
    if timer: print (os.path.basename(file))
    P = PilatusImage(file)
    if timer: G = Profiler("cbflib no-opt    read")
    P.read(algorithm="cbflib")
    read1 = P.linearintdata
    if timer: G = Profiler("cbflib optimized read")
    P.linearintdata = None #won't read again without resetting first
    P.read(algorithm="cbflib_optimized")
    read2 = P.linearintdata
    if timer: G = Profiler("buffer-based read")
    P.linearintdata = None #won't read again without resetting first
    P.read(algorithm="buffer_based")
    read3 = P.linearintdata
    if timer: del G
    expected_image_size = {"Pilatus-6M":(2527,2463),
                           "Pilatus-2M":(1679,1475),
                           "Pilatus-300K":(619,487)}[P.vendortype]
    assert read1.accessor().focus() == read2.accessor().focus() == expected_image_size
    from cbflib_adaptbx import assert_equal
    #print "Equality of arrays from two decompress methods", assert_equal(read1,read2), "\n"
    assert assert_equal(read1,read2)
    assert assert_equal(read1,read3)

if __name__=="__main__":
  test_all(dirpath=".",timer=False)
  print ("OK")
