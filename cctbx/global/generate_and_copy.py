import os, shutil
import generate_array_bpl
import generate_vector_algebra
generate_array_bpl.run()
generate_vector_algebra.run()
print "Copying to", "../cctbx/array_bpl.h"
shutil.copy("array_bpl.h", "../cctbx/array_bpl.h")
os.unlink("array_bpl.h")
print "Copying to", "../cctbx/vector/algebra.h"
shutil.copy("algebra.h", "../cctbx/vector/algebra.h")
os.unlink("algebra.h")
