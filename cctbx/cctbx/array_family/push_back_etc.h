// $Id$
// Included from small_plain.h and shared_plain.h
// DO NOT INCLUDE THIS FILE DIRECTLY.

      void clear() { this->resize(0); }

      void push_back(const ElementType& x) {
        this->auto_resize(this->size()+1, x);
      }
      void pop_back() { this->resize(this->size()-1); }

      ElementType*
      insert(ElementType* pos, size_type n, const ElementType& x) {
        ElementType* old_end = this->end();
        this->resize(this->size()+n);
        std::copy_backward(pos, old_end, this->end());
        std::fill(pos, pos+n, x);
        return pos;
      }

      ElementType*
      insert(ElementType* pos, const ElementType& x) {
        return insert(pos, 1, x);
      }

      template <typename OtherElementType>
      ElementType*
      insert(
        ElementType* pos,
        const OtherElementType* first,
        const OtherElementType* last) {
        size_type n = last - first;
        if (n > 0) {
          ElementType* old_end = this->end();
          this->resize(this->size()+n);
          std::copy_backward(pos, old_end, this->end());
          copy_typeconv(first, last, pos);
        }
        return pos;
      }

      ElementType* erase(ElementType* first, ElementType* last) {
        size_type n = last - first;
        std::copy(last, this->end(), first);
        this->resize(this->size()-n);
        return first;
      }

      ElementType* erase(ElementType* pos) {
        return erase(pos, pos+1);
      }

      void assign(size_type n, const ElementType& x = ElementType()) {
        this->resize(n);
        std::fill(this->begin(), this->end(), x);
      }

      template <typename OtherElementType>
      void assign(
        const OtherElementType* first,
        const OtherElementType* last) {
        this->resize(last - first);
        copy_typeconv(first, last, this->begin());
      }
