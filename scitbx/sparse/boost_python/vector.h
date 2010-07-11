#ifndef SCITBX_SPARSE_BOOST_PYTHON_VECTOR_H
#define SCITBX_SPARSE_BOOST_PYTHON_VECTOR_H

#include <boost/python/dict.hpp>
#include <scitbx/sparse/vector.h>

namespace scitbx { namespace sparse { namespace boost_python {

  template <typename T, template<class> class C>
  struct vector_from_dict
  {
    typedef vector<T, C> wt;
    typedef typename wt::index_type index_type;

    static wt make_on_stack(index_type n, boost::python::dict d) {
      wt result(n);
      fill(result, d);
      return result;
    }

    static wt *make_on_heap(index_type n, boost::python::dict d) {
      wt *result = new vector<T, C>(n);
      fill(*result, d);
      return result;
    }

    static void fill(wt &result, boost::python::dict d) {
      using namespace boost::python;
      list key = d.keys();
      index_type nz = len(key);
      for (index_type l=0; l<nz; ++l) {
        object k = key[l];
        index_type i = extract<index_type>(k);
        T x = extract<T>(d[k]);
        result[i] = x;
      }
      result.compact();
    }
  };

}}}

#endif // GUARD
