// Included from tiny.h, small.h, shared.h, versa.h
// DO NOT INCLUDE THIS FILE DIRECTLY.

      bool all_eq(base_class const& other) const
      {
        return this->const_ref().all_eq(other.const_ref());
      }

      bool all_eq(ElementType const& other) const
      {
        return this->const_ref().all_eq(other);
      }

      bool all_ne(base_class const& other) const
      {
        return this->const_ref().all_ne(other.const_ref());
      }

      bool all_ne(ElementType const& other) const
      {
        return this->const_ref().all_ne(other);
      }

      bool all_lt(base_class const& other) const
      {
        return this->const_ref().all_lt(other.const_ref());
      }

      bool all_lt(ElementType const& other) const
      {
        return this->const_ref().all_lt(other);
      }

      bool all_gt(base_class const& other) const
      {
        return this->const_ref().all_gt(other.const_ref());
      }

      bool all_gt(ElementType const& other) const
      {
        return this->const_ref().all_gt(other);
      }

      bool all_le(base_class const& other) const
      {
        return this->const_ref().all_le(other.const_ref());
      }

      bool all_le(ElementType const& other) const
      {
        return this->const_ref().all_le(other);
      }

      bool all_ge(base_class const& other) const
      {
        return this->const_ref().all_ge(other.const_ref());
      }

      bool all_ge(ElementType const& other) const
      {
        return this->const_ref().all_ge(other);
      }

      bool all_approx_equal(
        base_class const& other,
        ElementType const& tolerance) const
      {
        return this->const_ref().all_approx_equal(other.const_ref(),tolerance);
      }

      bool all_approx_equal(
        ElementType const& other,
        ElementType const& tolerance) const
      {
        return this->const_ref().all_approx_equal(other, tolerance);
      }
