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
#include <cctbx/sgtbx/basic.h>

namespace sgtbx {

  //! class for communicating string parsing errors.
  /*! This class is used by functions such as sgtbx::SgOps::ParseHallSymbol()
      or a constructor of class sgtbx::RTMx to communitcate errors
      that are detected during the interpretation of an input string.<br>
      Intended use:<pre>
      using namespace sgtbx;
      parse_string HSym("P x");
      SgOps sgops;
      try {
        sgops.ParseHallSymbol(HSym);
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
    private:
      std::string s;
      int pos;
      int marked_pos;
    public:
      //! Initialize the parse_string with the input string.
      inline parse_string(const std::string& str)
        : s(str), pos(0), marked_pos(0) { }
      //! Get back the input string.
      inline const char* string() const { return s.c_str(); }
      //! Index of last character accessed by the parsing algorithm.
      inline int where() const { return pos; }
      //! For internal use only.
      inline char operator()() const { return s[pos]; }
      //! For internal use only.
      inline const char* peek() { return s.c_str() + pos; }
      //! For internal use only.
      inline char get() {
        if (pos >= s.size()) return '\0';
        return s[pos++];
      }
      //! For internal use only.
      inline void skip(int n = 1) { pos += n; }
      //! For internal use only.
      inline void set_mark() { marked_pos = pos; }
      //! For internal use only.
      inline void go_to_mark() { pos = marked_pos; }
  };

  template <class T>
  inline void swap(T* a, T* b, int n)
  {
    for(int i=0;i<n;i++) {
      T tmp = a[i];
      a[i] = b[i];
      b[i] = tmp;
    }
  }

  template <class T>
  class Array2DAdaptor {
    private:
      T*  m_M;
      int m_mc;
    public:
      inline Array2DAdaptor(T *M, int mc) : m_M(M), m_mc(mc) { }
      inline T& operator()(int ir, int ic) { return m_M[ir * m_mc + ic]; }
  };

  int ChangeBaseFactor(const int *Old, int OldBF, int *New, int NewBF, int n);

  int SignHemisphere(const Vec3& v);

  class CmpiVect {
    private:
      int m_n;
    public:
      inline CmpiVect(int n) : m_n(n) { }
      bool operator()(const int *a, const int *b) const;
  };

  int NextOf_n_from_m(int m, int n, int *ix);

  template <typename ArrayType>
  class NestedLoop
  {
    public:
      NestedLoop() : m_over(1) {}
      NestedLoop(const ArrayType& min, const ArrayType& max)
        : m_min(min), m_max(max), m_cur(m_min), m_over(0) {}
      bool incr()
      {
        std::size_t i = m_cur.size() - 1;
        for (;;) {
          m_cur[i]++;
          if (m_cur[i] <= m_max[i]) return true;
          m_cur[i] = m_min[i];
          if (i == 0) {
            m_over++;
            break;
          }
          i--;
        }
        return false;
      }
      inline const ArrayType& min() const { return m_min; }
      inline const ArrayType& max() const { return m_max; }
      inline const ArrayType& operator()() { return m_cur; }
      inline std::size_t over() const { return m_over; }
      inline const ArrayType& next() { incr(); return m_cur; }
    private:
      ArrayType m_min;
      ArrayType m_max;
      ArrayType m_cur;
      std::size_t m_over;
  };

} // namespace sgtbx

#endif // CCTBX_SGTBX_UTILS_H
