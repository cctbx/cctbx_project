#ifndef SCITBX_SPARSE_IO_H
#define SCITBX_SPARSE_IO_H

#include <ostream>
#include <iomanip>
#include <scitbx/error.h>
#include <scitbx/sparse/vector.h>
#include <scitbx/sparse/matrix.h>
#include <scitbx/sparse/lu_factorization.h>
#include <scitbx/array_family/simple_io.h>

namespace scitbx { namespace sparse {

template<class T, template<class> class C>
struct vector_dense_display_t
{
  vector_dense_display_t(const sparse::vector<T, C>& v) : content(v) {}
  const sparse::vector<T, C>& content;
};

template<class T, template<class> class C>
vector_dense_display_t<T, C> dense_display(const sparse::vector<T, C>& v) {
  return vector_dense_display_t<T, C>(v);
}

template<class T, template<class> class C>
std::ostream& operator<<(std::ostream& o,
                         const sparse::vector_dense_display_t<T, C>& disp)
{
  typedef typename sparse::vector<T, C>::index_type index_type;
  typedef typename sparse::vector<T, C>::index_difference_type
          index_difference_type;
  typedef typename sparse::vector<T, C>::const_iterator const_iterator;
  std::streamsize width = o.width();
  sparse::vector<T, C> const& v = disp.content;
  v.compact();
  o << std::setw(0) << "{ ";
  index_difference_type last_non_zero = -1;
  bool first = true;
  for (const_iterator p  = v.begin(); p != v.end(); p++) {
    if (!first) o << ", ";
    else first = false;
    for (index_type i=1; i < p.index() - last_non_zero; i++) {
      o << std::setw(width) << "0" << ", ";
    }
    last_non_zero = p.index();
    o << std::setw(width) << *p;
  }
  index_difference_type trailing_zeroes_ = v.size() - (last_non_zero+1);
  if (trailing_zeroes_ > 0) {
    index_type trailing_zeroes = trailing_zeroes_;
    if (trailing_zeroes < v.size()) o << ", ";
    for (index_type i=1; i < trailing_zeroes; i++) {
      o << std::setw(width) << "0" << ", ";
    }
    o << std::setw(width) << "0";
  }
  o << " }";
  return o;
}

template<class T, template<class> class C>
struct vector_compressed_display_t {
  vector_compressed_display_t(const sparse::vector<T, C>& v) : content(v) {}
  const sparse::vector<T, C>& content;
};

template<class T, template<class> class C>
vector_compressed_display_t<T, C>
compressed_display(const sparse::vector<T, C>& v) {
  return vector_compressed_display_t<T, C>(v);
}

template<class T, template<class> class C>
std::ostream& operator<<(std::ostream& o,
                         const sparse::vector_compressed_display_t<T, C>& disp)
{
  sparse::vector<T, C> const& v = disp.content;
  typedef typename sparse::vector<T, C>::index_type size_type;
  typedef typename sparse::vector<T, C>::const_iterator const_iterator;
  o << "{ ";
  bool first = true;
  for (const_iterator p  = v.begin(); p != v.end(); p++) {
    if (!first) o << ", ";
    else first = false;
    o << p.index() << ": " << *p;
   }
  o << " }";
  return o;
}

template<class T>
struct matrix_dense_display_t {
  matrix_dense_display_t(const sparse::matrix<T>& m) : content(m) {}
  const sparse::matrix<T>& content;
};

template<class T>
matrix_dense_display_t<T> dense_display(const sparse::matrix<T>& m) {
  return matrix_dense_display_t<T>(m);
}

template<class T>
std::ostream& operator<<(std::ostream& o,
                         const sparse::matrix_dense_display_t<T>& disp)
{
  typedef typename sparse::matrix<T>::index_type index_type;
  typedef typename sparse::matrix<T>::const_row_iterator const_row_iterator;
  const sparse::matrix<T>& m = disp.content;
  sparse::matrix<T> mt = m.transpose();
  std::streamsize width = o.width();
  o << std::setw(0) << "{\n";
  for (index_type j=0; j < mt.n_cols(); j++) {
    o << std::setw(width) << dense_display(mt.col(j));
    if (j != mt.n_cols()-1) o << ",";
    o << "\n";
  }
  o << "}\n";
  return o;
}

template<class M>
std::ostream& operator << (std::ostream& o,
                           sparse::gilbert_peierls_lu_factorization<M>& lu)
{
  std::streamsize width = o.width();
  o << "A = \n" << std::setw(width) << dense_display(lu.factored());
  o << "\nL = \n" << std::setw(width) << dense_display(lu.l());
  o << "\nU = \n" << std::setw(width) << dense_display(lu.u());
  o << "\nrows permutation = " << lu.rows_permutation().const_ref();
  return o;
}


}}




#endif
