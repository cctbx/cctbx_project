import os, shutil
import generate_boost_array_bpl
import generate_vector_algebra
generate_boost_array_bpl.run()
generate_vector_algebra.run()
print "Copying to", "../cctbx/basic/boost_array_bpl.h"
shutil.copy("boost_array_bpl.h", "../cctbx/basic/boost_array_bpl.h")
os.unlink("boost_array_bpl.h")
print "Copying to", "../cctbx/vector/algebra.h"
shutil.copy("algebra.h", "../cctbx/vector/algebra.h")
os.unlink("algebra.h")
