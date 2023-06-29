#ifndef IOTBX_PDB_SMALL_STR_H
#define IOTBX_PDB_SMALL_STR_H

#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <ctype.h> // cannot use <cctype> since MIPSpro 7.3.1.2 defines macros

namespace iotbx { namespace pdb {

  inline
  unsigned
  stripped_size(const char* s)
  {
    for(;;) {
      if (*s == '\0') return 0;
      if (!isspace(*s)) break;
      s++;
    }
    unsigned i = 0;
    for(unsigned j=1;s[j]!='\0';j++) {
      if (!isspace(s[j])) i = j;
    }
    return i+1U;
  }

  inline
  unsigned
  rstripped_size(const char* s)
  {
    unsigned i = static_cast<unsigned>(std::strlen(s));
    if (i == 0) return 0;
    for(;;) {
      i--;
      if (!isspace(s[i])) {
        return i+1U;
      }
      if (i == 0) {
        return 0;
      }
    }
  }

  inline
  void
  copy_left_justified(
    char* dest,
    unsigned dest_size,
    const char *src,
    unsigned src_size,
    char pad_with)
  {
    unsigned i = 0;
    if (src != 0) {
      unsigned n = (dest_size < src_size ? dest_size : src_size);
      while(i<n) {
        char c = src[i];
        if (c == '\0') break;
        dest[i++] = c;
      }
    }
    while(i<dest_size) dest[i++] = pad_with;
  }

  inline
  void
  copy_right_justified(
    char* dest,
    unsigned dest_size,
    const char *src,
    unsigned src_size,
    char pad_with)
  {
    unsigned i = 0;
    if (src == 0) {
      while(i<dest_size) dest[i++] = pad_with;
    }
    else {
      unsigned n = (dest_size < src_size ? dest_size : src_size);
      while(i<n) {
        if (src[i] == '\0') break;
        i++;
      }
      while(i<dest_size) {
        *dest++ = pad_with;
        i++;
      }
      i = 0;
      while(i<n) {
        char c = src[i];
        if (c == '\0') break;
        dest[i++] = c;
      }
    }
  }

  // inline
  // void
  // copy_right_justified_std_string(
  //   std::string dest,
  //   const char *src,
  //   unsigned src_size,
  //   char pad_with)
  // {

  // }

  enum small_str_no_init_t { small_str_no_init };

  template <unsigned N>
  struct small_str
  {
    char elems[N+1U];

    small_str(small_str_no_init_t) {}

    small_str() { elems[0] = '\0'; }

    small_str(char c)
    {
      elems[0] = c;
      elems[1] = '\0';
    }

    small_str(const char* s, bool truncate_to_fit=false)
    {
      replace_with(s, truncate_to_fit);
    }

    small_str(
      const char* s_data,
      unsigned s_size,
      unsigned i_begin=0,
      char pad_with='\0');

    void
    replace_with(const char* s, bool truncate_to_fit=false);

    static
    unsigned
    capacity() { return N; }

    unsigned
    size() const { return static_cast<unsigned>(std::strlen(elems)); }

    unsigned
    stripped_size() const { return pdb::stripped_size(elems); }

    small_str<N>
    strip() const;

    unsigned
    rstripped_size() const { return pdb::rstripped_size(elems); }

    small_str<N> &
    upper_in_place()
    {
      for(char* e=elems; *e != '\0'; e++) {
        *e = toupper(*e);
      }
      return *this;
    }

    small_str<N> &
    replace_in_place(char old, char new_)
    {
      for(char* e=elems; *e != '\0'; e++) {
        if (*e == old) *e = new_;
      }
      return *this;
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

    void
    copy_left_justified(char* dest, unsigned dest_size, char pad_with) const
    {
      pdb::copy_left_justified(dest, dest_size, elems, N, pad_with);
    }

    void
    copy_right_justified(char* dest, unsigned dest_size, char pad_with) const
    {
      pdb::copy_right_justified(dest, dest_size, elems, N, pad_with);
    }

    template <unsigned OtherN>
    small_str<N+OtherN>
    concatenate(
      small_str<OtherN> const& other) const;
  };

  template <unsigned N>
  small_str<N>::small_str(
      const char* s_data,
      unsigned s_size,
      unsigned i_begin,
      char pad_with)
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

  template <unsigned N>
  void
  small_str<N>::replace_with(const char* s, bool truncate_to_fit)
  {
    if (s == 0) s = "";
    unsigned i = 0;
    while (i < N) {
      elems[i++] = *s;
      if (*s++ == '\0') return;
    }
    elems[i] = '\0';
    if (!truncate_to_fit && *s != '\0') {
      for(i=1U; s[i] != '\0'; i++);
      char buf[128];
      std::snprintf(buf, sizeof(buf),
        "string is too long for target variable"
        " (maximum length is %u character%s, %u given).",
          capacity(), (capacity() == 1U ? "" : "s"), N+i);
      throw std::invalid_argument(buf);
    }
  }

  template <unsigned N>
  small_str<N>
  small_str<N>::strip() const
  {
    const char* e = elems;
    for(;;) {
      if (*e == '\0') return small_str<N>();
      if (!isspace(*e)) break;
      e++;
    }
    unsigned i = 0;
    for(unsigned j=1;e[j]!='\0';j++) {
      if (!isspace(e[j])) i = j;
    }
    small_str<N> result(small_str_no_init);
    std::memcpy(result.elems, e, ++i);
    result.elems[i] = '\0';
    return result;
  }

  template <unsigned N>
  template <unsigned OtherN>
  small_str<N+OtherN>
  small_str<N>::concatenate(
    small_str<OtherN> const& other) const
  {
    small_str<N+OtherN> result((small_str_no_init));
    char *r = result.elems;
    const char *c = elems;
    while(*c != '\0') *r++ = *c++;
    c = other.elems;
    while(*c != '\0') *r++ = *c++;
    *r = '\0';
    return result;
  }

  typedef small_str<1> str1;
  typedef small_str<2> str2;
  typedef small_str<3> str3;
  typedef small_str<4> str4;
  typedef small_str<5> str5;
  typedef small_str<6> str6;
  typedef small_str<7> str7;
  typedef small_str<8> str8;

}} // namespace iotbx::pdb

#endif // IOTBX_PDB_SMALL_STR_H
