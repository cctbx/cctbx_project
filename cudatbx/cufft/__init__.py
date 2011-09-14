
import atexit

def real_to_complex_3d_in_place (data) :
  import boost.python
  ext = boost.python.import_ext("cudatbx_cufft_ext")
  real_type = type(data).__name__
  if (real_type == "double") :
    return ext.real_to_complex_3d_in_place_dp(data)
  else :
    assert (real_type == "float")
    return ext.real_to_complex_3d_in_place_sp(data)

def complex_to_complex_3d_in_place (data, direction) :
  import boost.python
  ext = boost.python.import_ext("cudatbx_cufft_ext")
  complex_type = type(data).__name__
  if (complex_type == "complex_double") :
    return ext.complex_to_complex_3d_in_place_dp(data, direction)
  else :
    assert (complex_type == "complex_float")
    return ext.complex_to_complex_3d_in_place_sp(data, direction)

def complex_to_real_3d_in_place (data, n) :
  import boost.python
  ext = boost.python.import_ext("cudatbx_cufft_ext")
  complex_type = type(data).__name__
  if (complex_type == "complex_double") :
    return ext.complex_to_real_3d_in_place_dp(data, n)
  else :
    assert (complex_type == "complex_float")
    return ext.complex_to_real_3d_in_place_sp(data, n)

def clean_up () :
  import boost.python
  ext = boost.python.import_ext("cudatbx_cufft_ext")
  ext.clean_up()

atexit.register(clean_up)
