// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_UTILS_H
#define CCTBX_SGTBX_UTILS_H

#include <string>
#include <algorithm>
#include <cctbx/sgtbx/basic.h>

namespace cctbx { namespace sgtbx {

  //! class for communicating string parsing errors.
  /*! This class is used by functions such as
      sgtbx::SpaceGroup::ParseHallSymbol()
      or a constructor of class sgtbx::RTMx to communitcate errors
      that are detected during the interpretation of an input string.<br>
      Intended use:<pre>
      using namespace sgtbx;
      parse_string HSym("P x");
      SpaceGroup SgOps;
      try {
        SgOps.ParseHallSymbol(HSym);
      }
      catch (const error& e) {
        std::cout << e.what() << std::endl;
        std::cout << HSym.string() << std::endl;
        for (int i = 0; i < HSym.where(); i++) std::cout << "_";
        std::cout << "^" << std::endl;
      }
      </pre>
      This will produce the output:<pre>
      cctbx Error: Improper symbol for rotational order.
      P x
      __^
      </pre>
   */
  class parse_string {
    public:
      //! Initialize the parse_string with the input string.
      explicit
      parse_string(const std::string& str = "")
        : s(str), pos(0), marked_pos(0) { }
      //! Get back the input string.
      const char* string() const { return s.c_str(); }
      //! Index of last character accessed by the parsing algorithm.
      int where() const { return pos; }
      //! For internal use only.
      char operator()() const { return s[pos]; }
      //! For internal use only.
      const char* peek() { return s.c_str() + pos; }
      //! For internal use only.
      char get() {
        if (pos >= s.size()) return '\0';
        return s[pos++];
      }
      //! For internal use only.
      void skip(int n = 1) { pos += n; }
      //! For internal use only.
      void set_mark() { marked_pos = pos; }
      //! For internal use only.
      void go_to_mark() { pos = marked_pos; }
    private:
      std::string s;
      int pos;
      int marked_pos;
  };

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  template <class AnyType>
  inline void
  swap(AnyType* a, AnyType* b, int n)
  {
    for(int i=0;i<n;i++) {
      AnyType tmp = a[i];
      a[i] = b[i];
      b[i] = tmp;
    }
  }

  template <class AnyType> // XXX use vecrefnd
  class Array2DAdaptor {
    public:
      Array2DAdaptor(AnyType *M, int mc) : m_M(M), m_mc(mc) { }
      AnyType& operator()(int ir, int ic) { return m_M[ir * m_mc + ic]; }
    private:
      AnyType*  m_M;
      int m_mc;
  };

  int ChangeBaseFactor(const int *Old, int OldBF, int *New, int NewBF, int n);

  int SignHemisphere(const Vec3& v);

  class CmpiVect {
    private:
      int m_n;
    public:
      CmpiVect(int n) : m_n(n) { }
      bool operator()(const int *a, const int *b) const;
  };

  template <std::size_t MaxN>
  class loop_n_from_m
  {
    public:
      loop_n_from_m() : m_m(0), m_n(0), m_over(1) {}
      loop_n_from_m(std::size_t m, std::size_t n) : m_m(m), m_n(n), m_over(0) {
        cctbx_assert(m_m >= m_n);
        cctbx_assert(MaxN >= m_n);
        for(std::size_t i=0;i<m_n;i++) m_current[i] = i;
      }
      bool incr()
      {
        if (m_n > 0) {
          std::size_t p, l;
          p = l = m_n - 1;
          for (;;) {
                m_current[p]++;
            if (m_current[p] + l == m_m + p) {
              if (p == 0) break;
              p--;
            }
            else if (p < l) {
              m_current[p + 1] = m_current[p];
                        p++;
            }
            else {
              return true;
            }
          }
        }
        m_over++;
        return false;
      }
      std::size_t m() { return m_m; }
      std::size_t n() { return m_n; }
      std::size_t operator[](std::size_t i) { return m_current[i]; }
      std::size_t over() const { return m_over; }
    private:
      std::size_t m_m;
      std::size_t m_n;
      std::size_t m_current[MaxN];
      std::size_t m_over;
  };

  template <typename ArrayType>
  class NestedLoop
  {
    public:
      NestedLoop() : m_over(1) {}
      NestedLoop(const ArrayType& end)
        : m_begin(end), m_end(end), m_current(end), m_over(0) {
        std::fill(m_begin.begin(), m_begin.end(), 0);
        m_current = m_begin;
      }
      NestedLoop(const ArrayType& begin, const ArrayType& end)
        : m_begin(begin), m_end(end), m_current(begin), m_over(0) {
        cctbx_assert(m_begin.size() == m_end.size());
      }
      bool incr()
      {
        for (std::size_t i = m_current.size(); i != 0;) {
          i--;
          m_current[i]++;
          if (m_current[i] < m_end[i]) return true;
          m_current[i] = m_begin[i];
        }
        m_over++;
        return false;
      }
      const ArrayType& begin() const { return m_begin; }
      const ArrayType& end() const { return m_end; }
      const ArrayType& operator()() { return m_current; }
      std::size_t over() const { return m_over; }
      const ArrayType& next() { incr(); return m_current; }
    private:
      ArrayType m_begin;
      ArrayType m_end;
      ArrayType m_current;
      std::size_t m_over;
  };

#endif // DOXYGEN_SHOULD_SKIP_THIS

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_UTILS_H
