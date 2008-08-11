// Included from small_plain.h and shared_plain.h
// DO NOT INCLUDE THIS FILE DIRECTLY.

      void assign(size_type const& sz, ElementType const& x) {
        if (sz > capacity()) {
          clear();
          reserve(sz);
          std::uninitialized_fill_n(begin(), sz, x);
          m_set_size(sz);
        }
        else if (sz > size()) {
          std::fill(begin(), end(), x);
          std::uninitialized_fill(end(), begin() + sz, x);
          m_set_size(sz);
        }
        else {
          std::fill_n(begin(), sz, x);
          erase(begin() + sz, end());
        }
      }

      void assign(const ElementType* first, const ElementType* last)
      {
        size_type sz = last - first;
        if (sz > capacity()) {
          clear();
          reserve(sz);
          std::uninitialized_copy(first, last, begin());
          m_set_size(sz);
        }
        else if (sz > size()) {
          const ElementType* mid = first + size() ;
          std::copy(first, mid, begin());
          std::uninitialized_copy(mid, last, end());
          m_set_size(sz);
        }
        else {
          std::copy(first, last, begin());
          erase(begin() + sz, end());
        }
      }

      void push_back(ElementType const& x) {
        if (size() < capacity()) {
          new (end()) ElementType(x);
          m_incr_size(1);
        }
        else {
          m_insert_overflow(end(), size_type(1), x, true);
        }
      }

      void pop_back() {
        m_decr_size(1);
        typedef typename has_trivial_destructor<ElementType>::value htd;
        detail::destroy_array_element(end(), htd());
      }

      ElementType* insert(ElementType* pos, ElementType const& x) {
        size_type n = pos - begin();
        if (size() == capacity()) {
          m_insert_overflow(pos, size_type(1), x, false);
        }
        else {
          if (pos == end()) {
            new (end()) ElementType(x);
            m_incr_size(1);
          }
          else {
            new (end()) ElementType(*(end() - 1));
            m_incr_size(1);
            ElementType x_copy = x;
            std::copy_backward(pos, end() - 2, end() - 1);
            *pos = x_copy;
          }
        }
        return begin() + n;
      }

      // restricted to ElementType
      void insert(ElementType* pos,
                  const ElementType* first, const ElementType* last)
      {
        size_type n = last - first;
        if (n == 0) return;
        if (size() + n > capacity()) {
          m_insert_overflow(pos, first, last);
        }
        else {
          size_type n_move_up = end() - pos;
          ElementType* old_end = end();
          if (n_move_up > n) {
            std::uninitialized_copy(end() - n, end(), end());
            m_incr_size(n);
            std::copy_backward(pos, old_end - n, old_end);
            std::copy(first, last, pos);
          }
          else {
            const ElementType* mid = first + n_move_up;
            std::uninitialized_copy(mid, last, end());
            m_incr_size(n - n_move_up);
            std::uninitialized_copy(pos, old_end, end());
            m_incr_size(n_move_up);
            std::copy(first, mid, pos);
          }
        }
      }

      // non-std
      void extend(const ElementType* first, const ElementType* last)
      {
        // pos == end();
        size_type n = last - first;
        if (n == 0) return;
        if (size() + n > capacity()) {
          m_insert_overflow(end(), first, last);
        }
        else {
          std::uninitialized_copy(first, last, end());
          m_incr_size(n);
        }
      }

      void insert(ElementType* pos, size_type const& n, ElementType const& x) {
        if (n == 0) return;
        if (size() + n > capacity()) {
          m_insert_overflow(pos, n, x, false);
        }
        else {
          ElementType x_copy = x;
          size_type n_move_up = end() - pos;
          ElementType* old_end = end();
          if (n_move_up > n) {
            std::uninitialized_copy(end() - n, end(), end());
            m_incr_size(n);
            std::copy_backward(pos, old_end - n, old_end);
            std::fill_n(pos, n, x_copy);
          }
          else {
            std::uninitialized_fill_n(end(), n - n_move_up, x_copy);
            m_incr_size(n - n_move_up);
            std::uninitialized_copy(pos, old_end, end());
            m_incr_size(n_move_up);
            std::fill(pos, old_end, x_copy);
          }
        }
      }

      ElementType* erase(ElementType* pos) {
        if (pos + 1 != end()) {
          std::copy(pos + 1, end(), pos);
        }
        pop_back();
        return pos;
      }

      ElementType* erase(ElementType* first, ElementType* last) {
        ElementType* i = std::copy(last, end(), first);
        typedef typename has_trivial_destructor<ElementType>::value htd;
        detail::destroy_array_elements(i, end(), htd());
        m_decr_size(last - first);
        return first;
      }

      void resize(size_type const& new_size, ElementType const& x) {
        if (new_size < size())  {
          erase(begin() + new_size, end());
        }
        else {
          insert(end(), new_size - size(), x);
        }
      }

      void resize(size_type const& new_size) {
        resize(new_size, ElementType());
      }

      void clear() {
        erase(begin(), end());
      }

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
      // non-std
      template <typename OtherElementType>
      void
      assign(const OtherElementType* first, const OtherElementType* last) {
        clear();
        size_type sz = last - first;
        reserve(sz);
        uninitialized_copy_typeconv(first, last, begin());
        m_set_size(sz);
      }
#endif

      // non-std
      template <typename OtherArrayType>
      void
      assign(OtherArrayType const& other) {
        if (other.size() == 0) {
          clear();
        }
        else {
          assign(&*(other.begin()), (&*(other.begin()))+other.size());
        }
      }
