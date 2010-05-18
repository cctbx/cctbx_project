#ifndef SCITBX_SPARSE_BOOST_PYTHON_VECTOR_H
#define SCITBX_SPARSE_BOOST_PYTHON_VECTOR_H

#include <boost/python/dict.hpp>
#include <scitbx/sparse/vector.h>

namespace scitbx { namespace sparse { namespace boost_python {

  template <typename T>
  struct vector_from_dict : vector<T>
  {
    typedef vector<T> wt;
    typedef typename wt::index_type index_type;

    vector_from_dict (index_type n, boost::python::dict d)
        : wt(n)
    {
      using namespace boost::python;
      list key = d.keys();
      index_type nz = len(key);
      for (index_type l=0; l<nz; ++l) {
        object k = key[l];
        index_type i = extract<index_type>(k);
        T x = extract<T>(d[k]);
        (*this)[i] = x;
      }
      this->compact();
    }
  };

}}}

#endif // GUARD
