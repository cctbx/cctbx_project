#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_PACKED_MATRIX_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_PACKED_MATRIX_H

namespace scitbx { namespace af {

  /// Conversion from packed size to matrix dimension
  inline
  unsigned
  dimension_from_packed_size(std::size_t packed_size)
  {
    unsigned n = static_cast<unsigned>(
      (std::sqrt(1.0+8.0*static_cast<double>(packed_size))-1.0)/2.0 + 0.5);
    SCITBX_ASSERT(n*(n+1)/2 == packed_size);
    return n;
  }

  /// Accessor for the upper diagonal of a square matrix packed by row
  class packed_u_accessor
  {
    public:
      typedef unsigned index_value_type;
      typedef af::tiny_plain<unsigned, 2> index_type;

      /// Construct an accessor to a n x n matrix
      explicit
      packed_u_accessor(unsigned n_=0) : n(n_) {}

      unsigned n_columns() const { return n; }
      unsigned n_rows() const { return n; }
      bool is_square() const { return true; }

      /// The size of the storage for the whole upper diagonal
      std::size_t
      size_1d() const { return n*(n+1)/2; }

      /// The index in the storage array for element (i,j) of the matrix
      /** Precondition: i <= j
          Not enforced for efficiency
      */
      unsigned
      operator()(unsigned i, unsigned j) const
      {
        return i*(n-1) - i*(i-1)/2 + j;
      }

      unsigned n;
  };


  /// Accessor for the lower diagonal of a square matrix packed by row
  class packed_l_accessor
  {
  public:
    typedef unsigned index_value_type;
    typedef af::tiny_plain<unsigned, 2> index_type;

    /// Construct an accessor to a n x n matrix
    explicit
    packed_l_accessor(unsigned n_=0) : n(n_) {}

    unsigned n_columns() const { return n; }
    unsigned n_rows() const { return n; }
    bool is_square() const { return true; }

    /// The size of the storage for the whole lower diagonal
    std::size_t
    size_1d() const { return n*(n+1)/2; }

    /// The index in the storage array for element (i,j) of the matrix
    /** Precondition: i >= j
     Not enforced for efficiency
     */
    unsigned
    operator()(unsigned i, unsigned j) const
    {
      return i*(i+1)/2 + j;
    }

    unsigned n;
  };

}}

#endif // GUARD
