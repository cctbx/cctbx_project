import os, shutil
import generate_carray_bpl
import generate_vector_algebra
generate_carray_bpl.run()
generate_vector_algebra.run()
print "Copying to", "../cctbx/carray_bpl.h"
shutil.copy("carray_bpl.h", "../cctbx/carray_bpl.h")
os.unlink("carray_bpl.h")
print "Copying to", "../cctbx/vector/algebra.h"
shutil.copy("algebra.h", "../cctbx/vector/algebra.h")
os.unlink("algebra.h")
