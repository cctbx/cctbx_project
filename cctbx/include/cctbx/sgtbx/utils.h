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

#include <cstddef>
#include <string>
#include <cctbx/sgtbx/basic.h>

namespace cctbx { namespace sgtbx {

  //! class for communicating string parsing errors.
  /*! This class is used by functions such as
      cctbx::sgtbx::SpaceGroup::ParseHallSymbol()
      or a constructor of class cctbx::sgtbx::RTMx to communitcate errors
      that are detected during the interpretation of an input string.<br>
      Intended use:<pre>
      using namespace cctbx::sgtbx;
      parse_string HSym("P x");
      SpaceGroup SgOps;
      try {
        SgOps.ParseHallSymbol(HSym);
      }
      catch (const cctbx::error& e) {
        std::cout << e.what() << std::endl;
        std::cout << HSym.string() << std::endl;
        for (std::size_t i = 0; i < HSym.where(); i++) std::cout << "_";
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
      std::size_t where() const { return pos; }
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
      void skip(std::size_t n = 1) { pos += n; }
      //! For internal use only.
      void set_mark() { marked_pos = pos; }
      //! For internal use only.
      void go_to_mark() { pos = marked_pos; }
    private:
      std::string s;
      std::size_t pos;
      std::size_t marked_pos;
  };

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  template <class AnyType>
  inline void
  swap(AnyType* a, AnyType* b, std::size_t n)
  {
    for(std::size_t i=0;i<n;i++) {
      AnyType tmp = a[i];
      a[i] = b[i];
      b[i] = tmp;
    }
  }

  int ChangeBaseFactor(const int *Old, int OldBF, int *New, int NewBF, int n);

  int SignHemisphere(const af::int3& v);

  class CmpiVect {
    private:
      int m_n;
    public:
      CmpiVect(int n) : m_n(n) { }
      bool operator()(const int *a, const int *b) const;
  };

#endif // DOXYGEN_SHOULD_SKIP_THIS

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_UTILS_H
