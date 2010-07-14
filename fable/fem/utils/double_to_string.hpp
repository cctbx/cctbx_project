#ifndef FEM_UTILS_DOUBLE_TO_STRING_HPP
#define FEM_UTILS_DOUBLE_TO_STRING_HPP

#include <fem/utils/char.hpp>

namespace fem { namespace utils {

  struct double_to_string_scientific_notation
  {
    static const int buffer_capacity = 256;
    char buffer[buffer_capacity];

    double_to_string_scientific_notation() {}

    double_to_string_scientific_notation(
      double value,
      int w, // "w" part of Fortran "w.d" format
      int d, // "d" part of Fortran "w.d" format
      int exp_scale, // Fortran "P" scaling
      char e_or_d)
    {
      ASSERTBX(w > 0);
      ASSERTBX(w < buffer_capacity);
      ASSERTBX(d >= 0);
      int nd = d;
      int na = d;
      if (exp_scale == 0) {
        // pass
      }
      else if (exp_scale > 0) {
        nd++;
        na -= exp_scale - 1;
      }
      else if (exp_scale < 0) {
        nd += exp_scale;
      }
      if (value == 0) {
        if (w < 1 + d + 4 || d == 0 || nd <= 0) {
          std::memset(buffer, '*', w);
          return;
        }
        int i_buf = w;
        buffer[--i_buf] = '0';
        buffer[--i_buf] = '0';
        buffer[--i_buf] = '+';
        buffer[--i_buf] = e_or_d;
        for(int i=0;i<d;i++) {
          buffer[--i_buf] = '0';
        }
        buffer[--i_buf] = '.';
        if (i_buf != 0) {
          buffer[--i_buf] = '0';
          while (i_buf != 0) {
            buffer[--i_buf] = ' ';
          }
        }
        return;
      }
      //
      struct pow_tab_entry {
        int i;
        double g;
        double l;
      };
      static const pow_tab_entry pow_tab[] = {
        {256, 1e256, 1e-256},
        {128, 1e128, 1e-128},
        { 64, 1e064, 1e-064},
        { 32, 1e032, 1e-032},
        { 16, 1e016, 1e-016},
        {  8, 1e008, 1e-008},
        {  4, 1e004, 1e-004},
        {  2, 1e002, 1e-002},
        {  1, 1e001, 1e-001}};
      bool v_negative = (value < 0);
      if (v_negative) value = -value;
      int required_w = (v_negative ? 2 : 1) + nd + 4;
      if (w < required_w || d == 0 || nd <= 0 || na < 0) {
        std::memset(buffer, '*', w);
        return;
      }
      int iexp = 0;
      if (value > 1) {
        for(int i_t=0;i_t<9;i_t++) {
          pow_tab_entry const& t = pow_tab[i_t];
          if (value > t.g) {
            iexp += t.i;
            value *= t.l;
          }
        }
        iexp++;
        value *= 0.1;
      }
      else if (value < 1) {
        for(int i_t=0;i_t<9;i_t++) {
          pow_tab_entry const& t = pow_tab[i_t];
          if (value < t.l) {
            iexp -= t.i;
            value *= t.g;
          }
        }
      }
      double multiplier = std::pow(10., nd); // XXX
      long ival = static_cast<long>(value * multiplier + 0.5);
      long jval = 1;
      for(int j=1;j<nd;j++) jval *= 10;
      if (ival == jval*10) {
        ival = jval;
        iexp++;
      }
      iexp -= exp_scale;
      //
      int i_buf = w;
      bool e_negative = (iexp < 0);
      if (e_negative) iexp = -iexp;
      while (iexp != 0) {
        buffer[--i_buf] = int_as_digit(iexp % 10);
        iexp /= 10;
      }
      while (i_buf > w-2) {
        buffer[--i_buf] = '0';
      }
      buffer[--i_buf] = (e_negative ? '-' : '+');
      if (i_buf != w-4) buffer[--i_buf] = e_or_d;
      //
      for(int i=0;i<na;i++) {
        buffer[--i_buf] = int_as_digit(ival % 10);
        ival /= 10;
      }
      buffer[--i_buf] = '.';
      if (ival == 0) {
        if (i_buf != 0) {
          buffer[--i_buf] = '0';
        }
      }
      else {
        do {
          buffer[--i_buf] = int_as_digit(ival % 10);
          ival /= 10;
        }
        while (ival != 0);
      }
      if (i_buf == 0) {
        if (v_negative) buffer[0] = '-';
      }
      else {
        buffer[--i_buf] = (v_negative ? '-' : ' ');
        while (i_buf != 0) {
          buffer[--i_buf] = ' ';
        }
      }
    }
  };

}} // namespace fem::utils

#endif // GUARD
