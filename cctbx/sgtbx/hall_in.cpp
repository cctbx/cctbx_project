/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/rotation_matrices.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>

namespace cctbx { namespace sgtbx {

  namespace {

    namespace hall {

      bool is_space(char c)
      {
        if (c == '\0') return false;
        if (c == '_') return true;
        return isspace(c);
      }

      bool is_char(char c)
      {
        if (c == '\0') return false;
        return !is_space(c);
      }

    } // namespace hall

    int get_abs_order(char c)
    {
      if (c == '1') return 1;
      if (c == '2') return 2;
      if (c == '3') return 3;
      if (c == '4') return 4;
      if (c == '6') return 6;
      return 0;
    }

    int get_screw(char c)
    {
      if (c == '1') return 1;
      if (c == '2') return 2;
      if (c == '3') return 3;
      if (c == '4') return 4;
      if (c == '5') return 5;
      return 0;
    }

    char get_reference_axis(char c)
    {
      c = tolower(c);
      if (   c == 'x'
          || c == 'y'
          || c == 'z') return c;
      return '\0';
    }

    char get_direction_code(char c)
    {
      if (   c == '\''
          || c ==  '"'
          || c ==  '*') return c;
      if (   c == ','
          || c == '.') return '\'';
      if (   c == ';'
          || c == ':') return  '"';
      return '\0';
    }

    tr_vec const& get_translation(char symbol)
    {
      struct tr_map
      {
        tr_map(char s, tr_vec const& v) : symbol(s), vec(v) {}
        char symbol;
        tr_vec vec;
      };
      static const tr_map hall_translations[] =
      {
        tr_map('a', tr_vec_12(6, 0, 0)),
        tr_map('b', tr_vec_12(0, 6, 0)),
        tr_map('c', tr_vec_12(0, 0, 6)),
        tr_map('n', tr_vec_12(6, 6, 6)),
        tr_map('u', tr_vec_12(3, 0, 0)),
        tr_map('v', tr_vec_12(0, 3, 0)),
        tr_map('w', tr_vec_12(0, 0, 3)),
        tr_map('d', tr_vec_12(3, 3, 3)),
      };
      const std::size_t n_hall_translations =
        sizeof hall_translations / sizeof (*hall_translations);
      symbol = tolower(symbol);
      for(std::size_t i=0;i<n_hall_translations;i++) {
        if (hall_translations[i].symbol == symbol) {
          return hall_translations[i].vec;
        }
      }
      static tr_vec null(0);
      return null;
    }

    const rot_mx get_rot_mx(bool improper, int abs_order,
                            char reference_axis, char direction_code)
    {
      struct tab_rot_mx_entry
      {
        int order;
        char direction_code;
        const rot_mx* r;
      };
      using namespace tables::rotation_matrices;
      const tab_rot_mx_entry tab_rot_mx[] = {
        { 1, '\0', &r_1_000 },
        { 2, '\0', &r_2_001 },
        { 2, '\'', &r_2_1b0 },
        { 2,  '"', &r_2_110 },
        { 3, '\0', &r_3_001 },
        { 3,  '*', &r_3_111 },
        { 4, '\0', &r_4_001 },
        { 6, '\0', &r_6_001 },
      };
      const std::size_t n = (sizeof tab_rot_mx / sizeof (*tab_rot_mx));
      for(std::size_t i=0;i<n;i++) {
        if (   tab_rot_mx[i].order == abs_order
            && tab_rot_mx[i].direction_code == direction_code) {
          rot_mx r;
          if (!improper)  r =  (*tab_rot_mx[i].r);
          else            r = -(*tab_rot_mx[i].r);
          if (reference_axis == 'x') return r_3_111 * r * r_3i111;
          if (reference_axis == 'y') return r_3i111 * r * r_3_111;
          return r;
        }
      }
      throw CCTBX_INTERNAL_ERROR();
    }

    tr_vec
    parse_short_cb_op(parse_string& hall_symbol,
                      const char* stop_chars,
                      int t_den)
    {
      CCTBX_ASSERT(t_den % 12 == 0);
      tr_vec result(t_den);
      for (std::size_t row = 0; row < 3; row++) {
        while (hall::is_space(hall_symbol())) hall_symbol.skip();
        if (row && hall_symbol() == ',') {
          hall_symbol.skip();
          while (hall::is_space(hall_symbol())) hall_symbol.skip();
        }
        if (strchr(stop_chars, hall_symbol())) return tr_vec(0);
        int i = 1;
        int n = sscanf(hall_symbol.peek(), "%d%n", &result[row], &i);
        hall_symbol.skip(i - 1);
        if (n != 1) return tr_vec(0);
        hall_symbol.skip();
        result[row] *= (t_den / 12);
      }
      return result;
    }

  } // namespace <anonymous>

  std::size_t
  space_group::parse_hall_symbol_cb_op(
    parse_string& hall_symbol,
    change_of_basis_op& cb_op,
    bool pedantic,
    bool no_centring_type_symbol)
  {
    std::size_t n_added_mx = 0;

    // Interpret the lattice type code.
    if (!no_centring_type_symbol) {
      while (hall::is_space(hall_symbol())) hall_symbol.skip();
      if (hall_symbol() == '-') {
        expand_inv(tr_vec(t_den()));
        hall_symbol.skip();
        n_added_mx++;
      }
      if (hall_symbol() == '\0') throw error("Lattice type not specified.");
      n_added_mx += expand_conventional_centring_type(hall_symbol());
      hall_symbol.skip();
    }

    const char char_after_lattice_type_symbol = hall_symbol();
    while (hall::is_space(hall_symbol())) hall_symbol.skip();
    if (hall_symbol() == '\0' || hall_symbol() == '(') {
      if (pedantic) throw error("Matrix symbol expected.");
      if (hall_symbol() == '\0') return n_added_mx;
    }
    if (   !no_centring_type_symbol
        && pedantic
        && !hall::is_space(char_after_lattice_type_symbol))
      throw error("Space expected after lattice type symbol.");

    // Loop over the matrix symbols.
    int  i_mx_symbol = 0;
    int  first_abs_order = 0;
    char first_reference_axis  = '\0';
    while (hall_symbol() != '\0' && hall_symbol() != '(')
    {
      bool improper = false;
      int abs_order =  0;
      int screw = 0;
      char reference_axis = '\0';
      char direction_code = '\0';
      rot_mx smx_r;
      tr_vec smx_t(t_den());

      if (hall_symbol() == '-') {
        improper = true;
        hall_symbol.skip();
        if (!hall::is_char(hall_symbol())) {
          throw error("Incomplete matrix symbol.");
        }
      }
           abs_order = get_abs_order(hall_symbol());
      if (!abs_order) {
        throw error("Improper symbol for rotational order.");
      }
      hall_symbol.skip();

          screw = get_screw(hall_symbol());
      if (screw) {
        if (screw >= abs_order) {
          throw error("Improper screw translation.");
        }
        hall_symbol.skip();
      }

      while (hall::is_char(hall_symbol()))
      {
        if (  reference_axis == '\0') {
              reference_axis = get_reference_axis(hall_symbol());
          if (reference_axis != '\0') {
            if (    abs_order == 1
                || (abs_order == 3 && direction_code == '*')) {
              throw error("Inconsistent matrix symbol.");
            }
            hall_symbol.skip();
            continue;
          }
        }
        else if (get_reference_axis(hall_symbol()) != '\0') {
          throw error("Multiple axis symbols.");
        }

        if (  direction_code == '\0') {
              direction_code = get_direction_code(hall_symbol());
          if (direction_code != '\0') {
            if (   !(abs_order == 2 && (  direction_code ==  '"'
                                       || direction_code == '\''))
                && !(abs_order == 3 && direction_code == '*')) {
              throw error("Inconsistent matrix symbol.");
            }
            if (screw) {
              throw error("Screw translation for non-principal direction.");
            }
            hall_symbol.skip();
            continue;
          }
        }
        else if (get_direction_code(hall_symbol()) != '\0') {
          throw error("Multiple axis symbols.");
        }

        tr_vec const& hall_translation = get_translation(hall_symbol());
        if (hall_translation.is_valid()) {
          smx_t = (smx_t
                   + hall_translation.new_denominator(t_den())).mod_positive();
          hall_symbol.skip();
          continue;
        }

        if (hall_symbol() == '(') {
          if (pedantic) {
            throw error("Space expected before change-of-basis operator.");
          }
          break;
        }

        throw error("Malformed matrix symbol.");
      }

      if (reference_axis == '\0') {
        if      (i_mx_symbol == 0) {
          if (     abs_order != 1
              && !(abs_order == 3 && direction_code == '*'))
            reference_axis = 'z';
        }
        else if (i_mx_symbol == 1) {
          if      (abs_order == 2) {
            if      (first_abs_order == 2 || first_abs_order == 4) {
              if (direction_code == '\0') {
                reference_axis = 'x';
              }
            }
            else if (first_abs_order == 3 || first_abs_order == 6) {
              if (direction_code == '\0') {
                direction_code = '\'';
              }
              reference_axis = first_reference_axis;
            }
          }
          else if (   abs_order == 3
                   && (first_abs_order == 2 || first_abs_order == 4)
                   && direction_code == '\0') {
            direction_code = '*';
          }
        }
        else if (i_mx_symbol == 2) {
          if (abs_order == 3 && direction_code == '\0') {
            direction_code = '*';
          }
        }
      }

      if (reference_axis == '\0' && (   direction_code ==  '"'
                              || direction_code == '\'')) {
        reference_axis = 'z';
      }

      if (reference_axis == '\0' && abs_order != 1 && direction_code != '*') {
        throw error("Need explicit axis symbol.");
      }

      smx_r = get_rot_mx(improper, abs_order, reference_axis, direction_code);

      if (screw) {
        std::size_t i;
        switch (reference_axis) {
          case 'x': i = 0; break;
          case 'y': i = 1; break;
          default:  i = 2; break;
        }
        CCTBX_ASSERT(smx_t.den() * screw % abs_order == 0);
        smx_t[i] += smx_t.den() * screw / abs_order;
      }

      expand_smx(rt_mx(smx_r, smx_t));

      if (i_mx_symbol == 0) {
        first_abs_order = abs_order;
        first_reference_axis  = reference_axis;
      }
      i_mx_symbol++;

      if (improper || abs_order != 1) {
        n_added_mx++;
      }

      while (hall::is_space(hall_symbol())) hall_symbol.skip();
    }

    // Interpret the change-of-basis operator.
    if (hall_symbol() == '(') {
      hall_symbol.skip();
      hall_symbol.set_mark();
      rt_mx v = rt_mx(parse_short_cb_op(hall_symbol, ")", cb_t_den), cb_r_den);
      if (!v.is_valid()) {
        hall_symbol.go_to_mark();
        v = rt_mx(hall_symbol, ")", cb_r_den, cb_t_den);
        if (!v.is_valid()) {
          throw error("Malformed change-of-basis operator.");
        }
      }
      while (hall::is_space(hall_symbol())) hall_symbol.skip();
      if (hall_symbol() != ')') {
        throw error(
          "Closing parenthesis expected after change-of-basis operator");
      }
      try {
        cb_op = change_of_basis_op(v);
      }
      catch (error const&) {
        throw error("Change-of-basis operator is not invertible.");
      }
      hall_symbol.skip();
    }

    while (hall::is_space(hall_symbol())) hall_symbol.skip();
    if (hall_symbol() != '\0') {
      throw error("Unexpected extra character.");
    }

    return n_added_mx;
  }

  std::size_t
  space_group::parse_hall_symbol(
    parse_string& hall_symbol,
    bool pedantic,
    bool no_centring_type_symbol)
  {
    change_of_basis_op cb_op(0, 0);
    std::size_t n_added_mx = parse_hall_symbol_cb_op(
      hall_symbol, cb_op, pedantic, no_centring_type_symbol);
    if (cb_op.is_valid()) {
      space_group new_sg = change_basis(cb_op);
      *this = new_sg;
    }
    return n_added_mx;
  }

  space_group::space_group(
    parse_string& hall_symbol,
    bool pedantic,
    bool no_centring_type_symbol,
    bool no_expand,
    int t_den)
  : no_expand_(no_expand)
  {
    reset(t_den);
    parse_hall_symbol(hall_symbol, pedantic, no_centring_type_symbol);
  }

  space_group::space_group(
    std::string const& hall_symbol,
    bool pedantic,
    bool no_centring_type_symbol,
    bool no_expand,
    int t_den)
  : no_expand_(no_expand)
  {
    reset(t_den);
    parse_string hall_symbol_ps(hall_symbol);
    parse_hall_symbol(hall_symbol_ps, pedantic, no_centring_type_symbol);
  }

  space_group::space_group(
    const char* hall_symbol,
    bool pedantic,
    bool no_centring_type_symbol,
    bool no_expand,
    int t_den)
  : no_expand_(no_expand)
  {
    reset(t_den);
    parse_string hall_symbol_ps(hall_symbol);
    parse_hall_symbol(hall_symbol_ps, pedantic, no_centring_type_symbol);
  }

  space_group::space_group(
    space_group_symbols const& symbols,
    int t_den)
  : no_expand_(false)
  {
    reset(t_den);
    parse_string hall_symbol_ps(symbols.hall());
    parse_hall_symbol(hall_symbol_ps, true);
  }

}} // namespace cctbx::sgtbx
