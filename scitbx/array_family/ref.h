#ifndef SCITBX_ARRAY_FAMILY_REF_H
#define SCITBX_ARRAY_FAMILY_REF_H

#include <scitbx/error.h>
#include <scitbx/math/traits.h>
#include <scitbx/array_family/error.h>
#include <scitbx/array_family/accessors/trivial.h>
#include <scitbx/array_family/detail/ref_helpers.h>
#include <algorithm>
#include <boost/scoped_array.hpp>

namespace scitbx { namespace af {

  template <typename ElementType,
            typename AccessorType = trivial_accessor>
  class const_ref
  {
    public:
      SCITBX_ARRAY_FAMILY_TYPEDEFS

      typedef AccessorType accessor_type;
      typedef typename accessor_type::index_type index_type;
      typedef typename accessor_type::index_value_type index_value_type;

      const_ref() {}

      const_ref(const ElementType* begin, accessor_type const& accessor)
      : begin_(begin), accessor_(accessor)
      {
        init();
      }

      // convenience constructors
      const_ref(const ElementType* begin, index_value_type const& n0)
      : begin_(begin), accessor_(n0)
      {
        init();
      }

      const_ref(const ElementType* begin, index_value_type const& n0,
                                          index_value_type const& n1)
      : begin_(begin), accessor_(n0, n1)
      {
        init();
      }

      const_ref(const ElementType* begin, index_value_type const& n0,
                                          index_value_type const& n1,
                                          index_value_type const& n2)
      : begin_(begin), accessor_(n0, n1, n2)
      {
        init();
      }

      accessor_type const& accessor() const { return accessor_; }
      size_type size() const { return size_; }

      /// Matrix ref interface (only make sense for relevant accessors)
      //@{
      bool is_square() const { return accessor_.is_square(); }
      std::size_t n_rows() const { return accessor_.n_rows(); }
      std::size_t n_columns() const { return accessor_.n_columns(); }

      bool is_diagonal(bool require_square=true) const;
      //@}


      const ElementType* begin() const { return begin_; }
      const ElementType* end() const { return end_; }
      const_reverse_iterator rbegin() const {
        return const_reverse_iterator(end_);
      }
      const_reverse_iterator rend() const {
        return const_reverse_iterator(begin_);
      }
      ElementType const& front() const { return begin_[0]; }
      ElementType const& back() const { return end_[-1]; }

      ElementType const&
      operator[](size_type i) const { return begin_[i]; }

      ElementType const&
      at(size_type i) const
      {
        if (i >= size_) throw_range_error();
        return begin_[i];
      }

      const_ref<ElementType>
      as_1d() const
      {
        return const_ref<ElementType>(begin_, size_);
      }

      value_type const&
      operator()(index_type const& i) const
      {
        return begin_[accessor_(i)];
      }

      // Convenience operator()
      value_type const& operator()(index_value_type const& i0) const
      {
        return begin_[accessor_(i0)];
      }

      value_type const& operator()(index_value_type const& i0,
                                   index_value_type const& i1) const
      {
        return begin_[accessor_(i0, i1)];
      }

      value_type const& operator()(index_value_type const& i0,
                                   index_value_type const& i1,
                                   index_value_type const& i2) const
      {
        return begin_[accessor_(i0, i1, i2)];
      }

      bool all_eq(const_ref const& other) const;

      bool all_eq(ElementType const& other) const;

      bool all_ne(const_ref const& other) const;

      bool all_ne(ElementType const& other) const;

      bool all_lt(const_ref const& other) const;

      bool all_lt(ElementType const& other) const;

      bool all_gt(const_ref const& other) const;

      bool all_gt(ElementType const& other) const;

      bool all_le(const_ref const& other) const;

      bool all_le(ElementType const& other) const;

      bool all_ge(const_ref const& other) const;

      bool all_ge(ElementType const& other) const;

      bool
      all_approx_equal(
        const_ref const& other,
        typename scitbx::math::abs_traits<
          ElementType>::result_type const& tolerance) const;

      bool
      all_approx_equal(
        ElementType const& other,
        typename scitbx::math::abs_traits<
          ElementType>::result_type const& tolerance) const;

      bool
      all_approx_equal_relatively(
        const_ref const& other,
        typename scitbx::math::abs_traits<
          ElementType>::result_type const& relative_error) const;

      bool
      all_approx_equal_relatively(
        ElementType const& other,
        typename scitbx::math::abs_traits<
          ElementType>::result_type const& relative_error) const;

    protected:
      void
      init()
      {
        size_ = accessor_.size_1d();
        end_ = begin_ + size_;
      }

      const ElementType* begin_;
      accessor_type accessor_;
      size_type size_;
      const ElementType* end_;
  };

  template<typename ElementType, class AccessorType>
  bool
  const_ref<ElementType, AccessorType>::is_diagonal(bool require_square) const {
    if (require_square && !is_square()) return false;
    for (index_value_type ir=0;ir<n_rows();ir++)
      for (index_value_type ic=0;ic<n_columns();ic++)
        if (ir != ic && (*this)(ir,ic)) return false;
    return true;
  }


  template <class E> class expression;


  template <typename ElementType,
            typename AccessorType = trivial_accessor>
  class ref : public const_ref<ElementType, AccessorType>
  {
    public:
      SCITBX_ARRAY_FAMILY_TYPEDEFS

      typedef const_ref<ElementType, AccessorType> base_class;
      typedef AccessorType accessor_type;
      typedef typename accessor_type::index_type index_type;
      typedef typename accessor_type::index_value_type index_value_type;

      ref() {}

      ref(ElementType* begin, accessor_type accessor)
      : base_class(begin, accessor)
      {}

      // convenience constructors
      ref(ElementType* begin, index_value_type const& n0)
      : base_class(begin, n0)
      {}

      ref(ElementType* begin, index_value_type const& n0,
                              index_value_type const& n1)
      : base_class(begin, n0, n1)
      {}

      ref(ElementType* begin, index_value_type const& n0,
                              index_value_type const& n1,
                              index_value_type const& n2)
      : base_class(begin, n0, n1, n2)
      {}

      template <class E>
      ref &operator=(expression<E> const &e);

      ElementType*
      begin() const { return const_cast<ElementType*>(this->begin_); }

      ElementType*
      end() const { return const_cast<ElementType*>(this->end_); }

      reverse_iterator rbegin() const {
        return reverse_iterator(this->end());
      }

      reverse_iterator rend() const {
        return reverse_iterator(this->begin());
      }

      ElementType&
      front() const { return begin()[0]; }

      ElementType&
      back() const { return end()[-1]; }

      ElementType&
      operator[](size_type i) const { return begin()[i]; }

      ElementType&
      at(size_type i) const
      {
        if (i >= this->size_) throw_range_error();
        return begin()[i];
      }

      ref const&
      fill(ElementType const& x) const
      {
        std::fill(begin(), end(), x);
        return *this;
      }

      ref<ElementType>
      as_1d() const
      {
        return ref<ElementType>(this->begin(), this->size_);
      }

      value_type&
      operator()(index_type const& i) const
      {
        return begin()[this->accessor_(i)];
      }

      // Convenience operator()
      value_type&
      operator()(index_value_type const& i0) const
      {
        return begin()[this->accessor_(i0)];
      }

      value_type&
      operator()(index_value_type const& i0,
                 index_value_type const& i1) const
      {
        return begin()[this->accessor_(i0, i1)];
      }

      value_type&
      operator()(index_value_type const& i0,
                 index_value_type const& i1,
                 index_value_type const& i2) const
      {
        return begin()[this->accessor_(i0, i1, i2)];
      }

      /// Matrix ref interface (only make sense for relevant accessors)
      //@{

      //! Swaps two rows in place.
      void
      swap_rows(index_value_type const& i1, index_value_type const& i2) const
      {
        std::swap_ranges(&(*this)(i1,0), &(*this)(i1+1,0), &(*this)(i2,0));
      }

      //! Swaps two columns in place.
      void
      swap_columns(index_value_type const& i1,
                   index_value_type const& i2) const
      {
        for(index_value_type ir=0;ir<this->n_rows();ir++) {
          std::swap((*this)(ir,i1), (*this)(ir,i2));
        }
      }

      //! Sets diagonal matrix.
      /*! Off-diagonal elements are set to zero.
       */
      void set_diagonal(ElementType const& d, bool require_square=true) const {
        SCITBX_ASSERT(!require_square || this->is_square());
        this->fill(0);
        index_value_type m = this->n_rows(), n = this->n_columns();
        for(index_value_type i=0; i < std::min(m,n); i++) (*this)(i,i) = d;
      }

      //! Sets diagonal matrix.
      /*! Off-diagonal elements are set to zero.
       */
      void set_diagonal(af::const_ref<ElementType> const& d,
                        bool require_square=true) const
      {
        SCITBX_ASSERT(!require_square || this->is_square());
        SCITBX_ASSERT(this->n_rows() >= d.size());
        SCITBX_ASSERT(this->n_columns() >= d.size());
        this->fill(0);
        index_value_type m = this->n_rows(), n = this->n_columns();
        for(index_value_type i=0; i < d.size(); i++) (*this)(i,i) = d[i];
      }

      //! Sets identity matrix.
      /*! Off-diagonal elements are set to zero.
       */
      void set_identity(bool require_square=true) const {
        set_diagonal(1, require_square);
      }

      /// Efficiently transpose a square matrix in-place
      void transpose_square_in_place() const {
        SCITBX_ASSERT(this->is_square());
        for (index_value_type ir=0;ir<this->n_rows();ir++)
          for (index_value_type ic=ir+1;ic<this->n_columns();ic++)
            std::swap((*this)(ir, ic), (*this)(ic, ir));
      }

      /// Transpose a matrix in-place.
      /** It involves a copy if it is not squared. */
      void transpose_in_place();

      //@}
  };

  template <typename ElementType, class AccessorType>
  void ref<ElementType, AccessorType>::transpose_in_place() {
    if (this->is_square()) {
      this->transpose_square_in_place();
    }
    else {
      boost::scoped_array<ElementType> mt_buffer(new ElementType[this->size()]);
      ref mt(mt_buffer.get(), this->n_columns(), this->n_rows());
      for (index_value_type ir=0;ir<this->n_rows();ir++)
        for (index_value_type ic=0;ic<this->n_columns();ic++)
          mt(ic, ir) = (*this)(ir, ic);
      std::copy(mt.begin(), mt.end(), this->begin());
      this->accessor_ = mt.accessor();
      this->init();
    }
  }

  template <typename ArrayType>
  const_ref<typename ArrayType::value_type>
  make_const_ref(ArrayType const& a)
  {
    typedef typename ArrayType::value_type value_type;
    typedef const_ref<value_type> return_type;
    typedef typename return_type::accessor_type accessor_type;
    return return_type(
      (a.size() == 0 ? 0 : &(*(a.begin()))),
      accessor_type(a.size()));
  }

  template <typename ArrayType>
  ref<typename ArrayType::value_type>
  make_ref(ArrayType& a)
  {
    typedef typename ArrayType::value_type value_type;
    typedef ref<value_type> return_type;
    typedef typename return_type::accessor_type accessor_type;
    return return_type(
      (a.size() == 0 ? 0 : &(*(a.begin()))),
      accessor_type(a.size()));
  }

  /// Wrapper for expression template
  /** This class is just a wrapper, which shall be inherited using the CRTP:

      class foo : public af::expression<foo>;

      Each member function of expression<E> forward to E's member function
      with the same name. Thus E shall implement those member functions
      which makes sense to its concept.

      See scitbx::sparse::matrix_times_dense_vector for a real-life example
      and scitbx/sparse/tests/tst_sparse for the natural syntax made
      possible by this mechanism
   */
  template <class E>
  class expression
  {
  public:
    E const &heir() const { return static_cast<E const &>(*this); }

    /// Total number of elements
    std::size_t size() const { return heir().size(); }

    /// Assign the elements of the expression to the memory referred to by x
    template <class ElementType, class AccessorType>
    void assign_to(af::ref<ElementType, AccessorType> const &x) const
    {
      heir().assign_to(x);
    }
  };

  template <class ElementType, class AccessorType>
  template <class E>
  ref<ElementType, AccessorType>
  &ref<ElementType, AccessorType>::operator=(expression<E> const &e) {
    e.assign_to(*this);
    return *this;
  }

  }} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_REF_H
