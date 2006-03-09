#ifndef IOTBX_PDB_SMALL_STR_H
#define IOTBX_PDB_SMALL_STR_H

#include <cstring>
#include <ctype.h> // cannot use <cctype> since MIPSpro 7.3.1.2 defines macros

namespace iotbx { namespace pdb {

  template <unsigned N>
  struct small_str
  {
    char elems[N+1];

    small_str() { elems[0] = '\0'; }

    small_str(char c)
    {
      elems[0] = c;
      elems[1] = '\0';
    }

    small_str(const char* s)
    {
      replace_with(s);
    }

    small_str(
      const char* s_data,
      unsigned s_size,
      unsigned i_begin=0,
      char pad_with='\0')
    {
      unsigned j = 0;
      while (i_begin < s_size) {
        if (j == N) {
          elems[j] = '\0';
          return;
        }
        elems[j++] = s_data[i_begin++];
      }
      if (pad_with != '\0') {
        while (j < N) {
          elems[j++] = ' ';
        }
      }
      elems[j] = '\0';
    }

    void
    replace_with(const char* s)
    {
      unsigned i = 0;
      while(i<N) {
        if (*s == '\0') break;
        elems[i++] = *s++;
      }
      elems[i] = '\0';
    }

    static
    unsigned
    capacity() { return N; }

    unsigned
    size() const { return std::strlen(elems); }

    unsigned
    stripped_size() const
    {
      const char* e = elems;
      do {
        if (*e == '\0') return 0;
      }
      while (isspace(*e++));
      unsigned i = 0;
      for(unsigned j=0;e[j]!='\0';j++) {
        if (!isspace(e[j])) i = j;
      }
      return i+1;
    }

    bool
    operator==(small_str const& other) const
    {
      return (std::strcmp(elems, other.elems) == 0);
    }

    bool
    operator!=(small_str const& other) const
    {
      return (std::strcmp(elems, other.elems) != 0);
    }

    bool
    operator<(small_str const& other) const
    {
      return (std::strcmp(elems, other.elems) < 0);
    }

    bool
    operator>(small_str const& other) const
    {
      return (std::strcmp(elems, other.elems) > 0);
    }

    bool
    operator<=(small_str const& other) const
    {
      return (std::strcmp(elems, other.elems) <= 0);
    }

    bool
    operator>=(small_str const& other) const
    {
      return (std::strcmp(elems, other.elems) >= 0);
    }
  };

  typedef small_str<1> str1;
  typedef small_str<2> str2;
  typedef small_str<4> str4;
  typedef small_str<6> str6;

}} // namespace iotbx::pdb

#endif // IOTBX_PDB_SMALL_STR_H
