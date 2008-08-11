#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

namespace {

  /* Slightly modified code from tcsh-6.13.00/sh.glob.c
     http://www.tcsh.org/
   */

  inline
  int
  globcharcoll(int c1, int c2)
  {
      return (c1 - c2);
  }

  /* t_pmatch():
   *      Return 2 on exact match,
   *      Return 1 on substring match.
   *      Return 0 on no match.
   *      *estr will point to the end of the longest exact or substring match.
   */
  int
  t_pmatch(const char *string, const char *pattern, const char **estr,
           char escape_char='\0')
  {
      char    stringc, patternc, patternc_esc;
      int     match, negate_range;
      char    rangec;
      const char *oestr, *pestr;

      for (;; ++string) {
          stringc = *string;
          /*
           * apollo compiler bug: switch (patternc = *pattern++) dies
           */
          patternc = patternc_esc = *pattern++;
          if (   escape_char != '\0'
              && patternc == escape_char) {
            patternc = *pattern++;
          }
          switch (patternc_esc) {
          case '\0':
              *estr = string;
              return (stringc == 0 ? 2 : 1);
          case '?':
              if (stringc == 0)
                  return (0);
              *estr = string;
              break;
          case '*':
              if (!*pattern) {
                  while (*string) string++;
                  *estr = string;
                  return (2);
              }
              oestr = *estr;
              pestr = 0;

              do {
                  switch(t_pmatch(string, pattern, estr, escape_char)) {
                  case 0:
                      break;
                  case 1:
                      pestr = *estr;
                      break;
                  case 2:
                      return 2;
                  default:
                      throw std::runtime_error("Internal error: " __FILE__);
                      /* Cannot happen */
                  }
                  *estr = string;
              }
              while (*string++);

              if (pestr) {
                  *estr = pestr;
                  return 1;
              }
              else {
                  *estr = oestr;
                  return 0;
              }

          case '[':
              match = 0;
              if ((negate_range = (*pattern == '^')) != 0)
                  pattern++;
              while ((rangec = *pattern++) != '\0') {
                  if (rangec == ']')
                      break;
                  if (match)
                      continue;
                  if (rangec == '-'
                      && *(pattern-2) != '[' && *pattern  != ']') {
                      match = (globcharcoll(stringc, *pattern) <= 0 &&
                      globcharcoll(*(pattern-2), stringc) <= 0);
                      pattern++;
                  }
                  else
                      match = (stringc == rangec);
              }
              if (rangec == 0)
                  return (-1);
              if (match == negate_range)
                  return (0);
              *estr = string;
              break;
          default:
              if (patternc != stringc)
                  return (0);
              *estr = string;
              break;
          }
      }
  }

  bool
  is_match(
    std::string const& string,
    std::string const& pattern,
    char escape_char='\0')
  {
    const char* estr;
    return (t_pmatch(string.c_str(), pattern.c_str(), &estr, escape_char)
            == 2);
  }

  BOOST_PYTHON_FUNCTION_OVERLOADS(is_match_overloads, is_match, 2, 3)

  void init_module()
  {
    using namespace boost::python;
    def("is_match", is_match, is_match_overloads(
      (arg_("string"), arg_("pattern"), arg_("escape_char")='\0')));
  }

} // namespace <anonymous>

BOOST_PYTHON_MODULE(iotbx_wildcard_ext)
{
  init_module();
}
