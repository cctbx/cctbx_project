#ifndef SCITBX_MATRIX_AS_XYZ_H
#define SCITBX_MATRIX_AS_XYZ_H

#include <scitbx/rational.h>
#include <scitbx/error.h>

namespace scitbx { namespace matrix {

  template <typename IntType>
  std::string
  rational_as_xyz(
    int n_rows,
    int n_columns,
    const IntType* r_num,
    IntType r_den,
    const IntType* t_num,
    IntType t_den,
    bool decimal,
    bool t_first,
    const char* letters_xyz,
    const char* separator)
  {
    SCITBX_ASSERT(letters_xyz != 0 && strlen(letters_xyz) == n_columns);
    SCITBX_ASSERT(separator != 0);
    std::string result;
    for (int i = 0; i < n_rows; i++) {
      std::string R_term;
      if (r_num != 0) {
        for (int j = 0; j < n_columns; j++) {
          boost::rational<IntType> R_frac(r_num[i*n_columns+j], r_den);
          if (R_frac != 0) {
            if (R_frac > 0) {
              if (!R_term.empty()) {
                R_term += "+";
              }
            }
            else {
              R_term += "-";
              R_frac *= -1;
            }
            if (R_frac != 1) {
              R_term += format(R_frac, decimal) + "*";
            }
            R_term += letters_xyz[j];
          }
        }
      }
      if (i != 0) result += separator;
      if (t_num == 0) {
        if (R_term.empty()) result += "0";
        else                result += R_term;
      }
      else {
        boost::rational<IntType> T_frac(t_num[i], t_den);
        if (T_frac == 0) {
          if (R_term.empty()) result += "0";
          else                result += R_term;
        }
        else if (R_term.empty()) {
          result += format(T_frac, decimal);
        }
        else if (t_first) {
          result += format(T_frac, decimal);
          if (R_term[0] != '-') result += "+";
          result += R_term;
        }
        else {
          result += R_term;
          if (T_frac > 0) result += "+";
          result += format(T_frac, decimal);
        }
      }
    }
    return result;
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_AS_XYZ_H
