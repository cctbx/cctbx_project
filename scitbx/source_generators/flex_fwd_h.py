from scitbx.source_generators.utils import join_open
from scitbx.source_generators.utils import write_this_is_auto_generated
import libtbx.load_env

this = "scitbx.source_generators.flex_fwd_h"

motivation = """\
/* The declarations in this file facilitate cross-module functionality
   on platforms that do not support comparison of type expressions
   across dynamically loaded library boundaries. On such platforms
   Boost.Python uses type_id::name() for comparing type expressions.
   For a given type, with some compilers (e.g. some EDG based
   compilers) the result of type_id::name() depends on the first type
   expression encountered in a translation unit. To ensure that
   type_id::name() produces the same result in all translation
   units, this file should be included at the top of all Boost.Python
   extension modules that involve the types in the function
   signatures below.
 */
"""

common_code = """\
#ifndef SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_FWD_H
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_FWD_H

#include <boost/python/detail/prefix.hpp>

#if defined(BOOST_MSVC) && BOOST_MSVC == 1600 && defined(_WIN64)
#pragma optimize("g", off)
#endif

%s#endif // SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_FWD_H
"""

full_code = """\
#define SCITBX_ARRAY_FAMILY_BOOST_PYTHON_FLEX_FWD_FULL_CODE

%s
#include <complex>
#include <vector>
#include <set>

namespace scitbx { namespace boost_python {

  struct misc_fwd
  {
    friend void f(std::vector<unsigned> const&);
    friend void f(std::set<unsigned> const&);
    friend void f(std::vector<std::set<unsigned> > const&);
  };

}} // namespace scitbx::boost_python

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/vec3.h>

#include <boost_adaptbx/std_pair_fwd.h>

#if defined(__sgi) && !defined(__GNUC__)

namespace scitbx { namespace af { namespace boost_python {

  struct misc_fwd
  {
    friend void f(af::tiny<unsigned, 2> const&);
    friend void f(af::tiny<unsigned, 3> const&);
    friend void f(af::tiny<unsigned, 4> const&);
    friend void f(af::tiny<double, 3> const&);
  };

  template <typename T>
  struct flex_fwd
  {
    friend void f(shared_plain<T> const&);
    friend void f(shared<T> const&);
    friend void f(versa<T> const&);
    friend void f(versa<T, flex_grid<> > const&);
    friend void f(versa<T, c_grid<2> > const&);
    friend void f(versa<T, c_grid<3> > const&);
    friend void f(versa<T, c_grid_padded<2> > const&);
    friend void f(versa<T, c_grid_padded<3> > const&);
    friend void f(ref<T> const&);
    friend void f(ref<T, flex_grid<> > const&);
    friend void f(ref<T, c_grid<2> > const&);
    friend void f(ref<T, c_grid<3> > const&);
    friend void f(ref<T, c_grid_padded<2> > const&);
    friend void f(ref<T, c_grid_padded<3> > const&);
    friend void f(const_ref<T> const&);
    friend void f(const_ref<T, flex_grid<> > const&);
    friend void f(const_ref<T, c_grid<2> > const&);
    friend void f(const_ref<T, c_grid<3> > const&);
    friend void f(const_ref<T, c_grid_padded<2> > const&);
    friend void f(const_ref<T, c_grid_padded<3> > const&);
  };

  inline void
  flex_fwd_types()
  {
    flex_fwd<bool>();
    flex_fwd<int>();
    flex_fwd<long>();
    flex_fwd<std::size_t>();
    flex_fwd<float>();
    flex_fwd<double>();
    flex_fwd<std::complex<double> >();
    flex_fwd<std::string>();
    flex_fwd<vec3<double> >();
    flex_fwd<tiny<std::size_t, 2> >();
    flex_fwd<sym_mat3<double> >();

    // for shared_ext.cpp
    flex_fwd<std::vector<std::size_t> >();
    flex_fwd<std::set<std::size_t> >();
  }

}}} // namespace scitbx::af::boost_python

#endif // defined(__sgi) && !defined(__GNUC__)

"""

def write(this, target_dir, common_code, full_code):
  f = join_open(target_dir, "flex_fwd.h", "w")
  write_this_is_auto_generated(f, this)
  if (libtbx.env.build_options.write_full_flex_fwd_h):
    code = full_code % motivation
  else:
    code = ""
  f.write(common_code % code)

def run(target_dir):
  write(this, target_dir, common_code, full_code)

if (__name__ == "__main__"):
  run(".")
