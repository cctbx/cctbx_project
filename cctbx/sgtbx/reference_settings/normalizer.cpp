/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored parts of cctbx/sgtbx/reference.h
               and normalizers.cpp (rwgk)
 */

#include <cctbx/sgtbx/reference_settings.h>
#include <cctbx/sgtbx/space_group.h>

namespace cctbx { namespace sgtbx { namespace reference_settings {

namespace normalizer {

  addl_generators const&
  addl_generators_table(std::size_t i)
  {
    static const addl_generators table[] = {
                {    0,   0 },
// BEGIN_COMPILED_IN_REFERENCE_DATA
      /*   1 */ { "-1",   0 }, // Note: the affine normalizers for
      /*   2 */ {    0,   0 }, //   space groups 1 to 15 CANNOT be
      /*   3 */ { "-1",   0 }, //   represented by rotation symbols.
      /*   4 */ { "-1",   0 },
      /*   5 */ { "-1",   0 },
      /*   6 */ { "-1",   0 },
      /*   7 */ { "-1",   0 },
      /*   8 */ { "-1",   0 },
      /*   9 */ { "-1",   0 },
      /*  10 */ {    0,   0 },
      /*  11 */ {    0,   0 },
      /*  12 */ {    0,   0 },
      /*  13 */ {    0,   0 },
      /*  14 */ {    0,   0 },
      /*  15 */ {    0,   0 },
      /*  16 */ { "-1",   "3* -2'" },  // Space groups 16 to 74:
      /*  17 */ { "-1",   "-2'w" },    //   the L2N column is for the
      /*  18 */ { "-1",   "-2'" },     //   affine normalizer only
      /*  19 */ { "-1",   "3* -2'd" },
      /*  20 */ { "-1",   "-2'w" },
      /*  21 */ { "-1",   "-2'" },
      /*  22 */ { "-1",   "3* -2'" },
      /*  23 */ { "-1",   "3* -2'" },
      /*  24 */ { "-1",   "3* -2'd" },
      /*  25 */ { "-1",   "-2'" },
      /*  26 */ { "-1",   0 },
      /*  27 */ { "-1",   "-2'" },
      /*  28 */ { "-1",   0 },
      /*  29 */ { "-1",   0 },
      /*  30 */ { "-1",   0 },
      /*  31 */ { "-1",   0 },
      /*  32 */ { "-1",   "-2'" },
      /*  33 */ { "-1",   0 },
      /*  34 */ { "-1",   "-2'" },
      /*  35 */ { "-1",   "-2'" },
      /*  36 */ { "-1",   0 },
      /*  37 */ { "-1",   "-2'" },
      /*  38 */ { "-1",   0 },
      /*  39 */ { "-1",   0 },
      /*  40 */ { "-1",   0 },
      /*  41 */ { "-1",   0 },
      /*  42 */ { "-1",   "-2'" },
      /*  43 */ { "-1uv", "-2'" },
      /*  44 */ { "-1",   "-2'" },
      /*  45 */ { "-1",   "-2'" },
      /*  46 */ { "-1",   0 },
      /*  47 */ {    0,   "3* -2'" },
      /*  48 */ {    0,   "3* -2'" },
      /*  49 */ {    0,   "-2'" },
      /*  50 */ {    0,   "-2'" },
      /*  51 */ {    0,   0 },
      /*  52 */ {    0,   0 },
      /*  53 */ {    0,   0 },
      /*  54 */ {    0,   0 },
      /*  55 */ {    0,   "-2'" },
      /*  56 */ {    0,   "-2'" },
      /*  57 */ {    0,   0 },
      /*  58 */ {    0,   "-2'" },
      /*  59 */ {    0,   "-2'" },
      /*  60 */ {    0,   0 },
      /*  61 */ {    0,   "3*" },
      /*  62 */ {    0,   0 },
      /*  63 */ {    0,   0 },
      /*  64 */ {    0,   0 },
      /*  65 */ {    0,   "-2'" },
      /*  66 */ {    0,   "-2'" },
      /*  67 */ {    0,   "-2'buv" },
      /*  68 */ {    0,   "-2'buv" },
      /*  69 */ {    0,   "3* -2'" },
      /*  70 */ {    0,   "3* -2'" },
      /*  71 */ {    0,   "3* -2'" },
      /*  72 */ {    0,   "-2'" },
      /*  73 */ {    0,   "3* -2'd" },
      /*  74 */ {    0,   "-2'bd" },
      /*  75 */ { "-1",   "-2'" },
      /*  76 */ {    0,   "2\"" },
      /*  77 */ { "-1",   "-2'" },
      /*  78 */ {    0,   "2\"" },
      /*  79 */ { "-1",   "-2'" },
      /*  80 */ { "-1a",  "2\"" },
      /*  81 */ { "-1",   "-2'" },
      /*  82 */ { "-1",   "-2'" },
      /*  83 */ {    0,   "-2'" },
      /*  84 */ {    0,   "-2'" },
      /*  85 */ {    0,   "-2'" },
      /*  86 */ {    0,   "-2'" },
      /*  87 */ {    0,   "-2'" },
      /*  88 */ {    0,   "-2'd" }, // Space group 88 is the only space
      /*  89 */ { "-1",   0 },      // group where the L2N column is
      /*  90 */ { "-1",   0 },      // different for the two origin choices.
      /*  91 */ {    0,   0 },
      /*  92 */ {    0,   0 },
      /*  93 */ { "-1",   0 },
      /*  94 */ { "-1",   0 },
      /*  95 */ {    0,   0 },
      /*  96 */ {    0,   0 },
      /*  97 */ { "-1",   0 },
      /*  98 */ { "-1aw", 0 },
      /*  99 */ { "-1",   0 },
      /* 100 */ { "-1",   0 },
      /* 101 */ { "-1",   0 },
      /* 102 */ { "-1",   0 },
      /* 103 */ { "-1",   0 },
      /* 104 */ { "-1",   0 },
      /* 105 */ { "-1",   0 },
      /* 106 */ { "-1",   0 },
      /* 107 */ { "-1",   0 },
      /* 108 */ { "-1",   0 },
      /* 109 */ { "-1a",  0 },
      /* 110 */ { "-1a",  0 },
      /* 111 */ { "-1",   0 },
      /* 112 */ { "-1",   0 },
      /* 113 */ { "-1",   0 },
      /* 114 */ { "-1",   0 },
      /* 115 */ { "-1",   0 },
      /* 116 */ { "-1",   0 },
      /* 117 */ { "-1",   0 },
      /* 118 */ { "-1",   0 },
      /* 119 */ { "-1",   0 },
      /* 120 */ { "-1",   0 },
      /* 121 */ { "-1",   0 },
      /* 122 */ { "-1aw", 0 },
      /* 123 */ {    0,   0 },
      /* 124 */ {    0,   0 },
      /* 125 */ {    0,   0 },
      /* 126 */ {    0,   0 },
      /* 127 */ {    0,   0 },
      /* 128 */ {    0,   0 },
      /* 129 */ {    0,   0 },
      /* 130 */ {    0,   0 },
      /* 131 */ {    0,   0 },
      /* 132 */ {    0,   0 },
      /* 133 */ {    0,   0 },
      /* 134 */ {    0,   0 },
      /* 135 */ {    0,   0 },
      /* 136 */ {    0,   0 },
      /* 137 */ {    0,   0 },
      /* 138 */ {    0,   0 },
      /* 139 */ {    0,   0 },
      /* 140 */ {    0,   0 },
      /* 141 */ {    0,   0 },
      /* 142 */ {    0,   0 },
      /* 143 */ { "-1",   "2 -2'" },
      /* 144 */ {    0,   "2 2\"" },
      /* 145 */ {    0,   "2 2\"" },
      /* 146 */ { "-1",   "-2\"" },
      /* 147 */ {    0,   "2 -2'" },
      /* 148 */ {    0,   "-2\"" },
      /* 149 */ { "-1",   "2" },
      /* 150 */ { "-1",   "2" },
      /* 151 */ {    0,   "2" },
      /* 152 */ {    0,   "2" },
      /* 153 */ {    0,   "2" },
      /* 154 */ {    0,   "2" },
      /* 155 */ { "-1",   0 },
      /* 156 */ { "-1",   "2" },
      /* 157 */ { "-1",   "2" },
      /* 158 */ { "-1",   "2" },
      /* 159 */ { "-1",   "2" },
      /* 160 */ { "-1",   0 },
      /* 161 */ { "-1",   0 },
      /* 162 */ {    0,   "2" },
      /* 163 */ {    0,   "2" },
      /* 164 */ {    0,   "2" },
      /* 165 */ {    0,   "2" },
      /* 166 */ {    0,   0 },
      /* 167 */ {    0,   0 },
      /* 168 */ { "-1",   "-2'" },
      /* 169 */ {    0,   "2\"" },
      /* 170 */ {    0,   "2\"" },
      /* 171 */ {    0,   "2\"" },
      /* 172 */ {    0,   "2\"" },
      /* 173 */ { "-1",   "-2'" },
      /* 174 */ { "-1",   "-2'" },
      /* 175 */ {    0,   "-2'" },
      /* 176 */ {    0,   "-2'" },
      /* 177 */ { "-1",   0 },
      /* 178 */ {    0,   0 },
      /* 179 */ {    0,   0 },
      /* 180 */ {    0,   0 },
      /* 181 */ {    0,   0 },
      /* 182 */ { "-1",   0 },
      /* 183 */ { "-1",   0 },
      /* 184 */ { "-1",   0 },
      /* 185 */ { "-1",   0 },
      /* 186 */ { "-1",   0 },
      /* 187 */ { "-1",   0 },
      /* 188 */ { "-1",   0 },
      /* 189 */ { "-1",   0 },
      /* 190 */ { "-1",   0 },
      /* 191 */ {    0,   0 },
      /* 192 */ {    0,   0 },
      /* 193 */ {    0,   0 },
      /* 194 */ {    0,   0 },
      /* 195 */ { "-1",   "-2'" },
      /* 196 */ { "-1",   "-2'" },
      /* 197 */ { "-1",   "-2'" },
      /* 198 */ { "-1",   "-2'd" },
      /* 199 */ { "-1",   "-2'd" },
      /* 200 */ {    0,   "-2'" },
      /* 201 */ {    0,   "-2'" },
      /* 202 */ {    0,   "-2'" },
      /* 203 */ {    0,   "-2'" },
      /* 204 */ {    0,   "-2'" },
      /* 205 */ {    0,   0 },
      /* 206 */ {    0,   "-2'd" },
      /* 207 */ { "-1",   0 },
      /* 208 */ { "-1",   0 },
      /* 209 */ { "-1",   0 },
      /* 210 */ { "-1d",  0 },
      /* 211 */ { "-1",   0 },
      /* 212 */ {    0,   0 },
      /* 213 */ {    0,   0 },
      /* 214 */ { "-1",   0 },
      /* 215 */ { "-1",   0 },
      /* 216 */ { "-1",   0 },
      /* 217 */ { "-1",   0 },
      /* 218 */ { "-1",   0 },
      /* 219 */ { "-1",   0 },
      /* 220 */ { "-1",   0 },
      /* 221 */ {    0,   0 },
      /* 222 */ {    0,   0 },
      /* 223 */ {    0,   0 },
      /* 224 */ {    0,   0 },
      /* 225 */ {    0,   0 },
      /* 226 */ {    0,   0 },
      /* 227 */ {    0,   0 },
      /* 228 */ {    0,   0 },
      /* 229 */ {    0,   0 },
      /* 230 */ {    0,   0 }
// END_COMPILED_IN_REFERENCE_DATA
    };
    return table[i];
  }

  af::shared<rt_mx>
  get_addl_generators(
    int sg_number,
    bool flag_affine,
    bool flag_k2l,
    bool flag_l2n)
  {
    CCTBX_ASSERT(0 < sg_number && sg_number <= 230);
    af::shared<rt_mx> result;
    for(int i_type=0;i_type<2;i_type++) {
      const char* hall_mx_symbol = 0;
      if      (i_type == 0 && flag_k2l)
        hall_mx_symbol = addl_generators_table(sg_number).k2l;
      else if (i_type == 1 && flag_l2n && (sg_number >= 75 || flag_affine))
        hall_mx_symbol = addl_generators_table(sg_number).l2n;
      if (hall_mx_symbol) {
        space_group sg_addl_g(true); // no_expand: no group multiplication
        parse_string ps(hall_mx_symbol);
        std::size_t n_added_mx = sg_addl_g.parse_hall_symbol(ps, true, true);
        CCTBX_ASSERT(n_added_mx > 0);
        CCTBX_ASSERT(sg_addl_g.n_ltr() == 1);
        if (sg_addl_g.is_centric()) {
          result.push_back(sg_addl_g(0, 1, 0));
        }
        for (std::size_t i=1;i<sg_addl_g.n_smx();i++) {
          result.push_back(sg_addl_g.smx(i));
        }
      }
    }
    return result;
  }

  void
  get_monoclinic_affine_trial_ranges(
    rot_mx const& cb_mx_r,
    int& r00,
    int& r22)
  {
    /* International Tables Volume A, chapter 15, tables 15.3.3 & 15.3.4.

     M.C = n00*c00 + n02*c20,  n00*c01 + n02*c21,  n00*c02 + n02*c22,
           c10,                c11,                c12,
           n20*c00 + n22*c20,  n20*c01 + n22*c21,  n20*c02 + n22*c22

     Determine trial range for n00 and n20:
       max(lcm(c00, c20) / c00,
           lcm(c01, c21) / c01,
           lcm(c02, c22) / c02)
     Determine trial range for n02 and n22:
       max(lcm(c00, c20) / c20,
           lcm(c01, c21) / c21,
           lcm(c02, c22) / c22)
     */
    using scitbx::fn::absolute;
    r00 = 1;
    r22 = 1;
    for(std::size_t i=0;i<3;i++) {
      int l = boost::lcm(cb_mx_r[i], cb_mx_r[6 + i]);
      if (cb_mx_r[i]) {
        int n = absolute(l / cb_mx_r[i]);
        if (r00 < n) r00 = n;
      }
      if (cb_mx_r[i + 6]) {
        int n = absolute(l / cb_mx_r[6 + i]);
        if (r22 < n) r22 = n;
      }
    }
    r00++;
    r22++;
  }

  bool
  check_monoclinic_affine_restrictions(
    int sg_number,
    rot_mx const& r)
  {
    // International Tables Volume A, chapter 15, tables 15.3.3 & 15.3.4.
    switch (sg_number) {
      case  3:
      case  4:
      case  6:
      case 10:
      case 11: /* M2 */
        break;

      case  5:
      case  8:
      case 12: /* M4 */
      case  9:
      case 15: /* M6 or M12 */
        if (r[0] % (2 * r.den()) == 0) return false;
        if (r[6] % (2 * r.den()) != 0) return false;
        if (r[8] % (2 * r.den()) == 0) return false;
        break;

      case  7:
      case 13:
      case 14: /* M5 */
        if (r[0] % (2 * r.den()) == 0) return false;
        if (r[2] % (2 * r.den()) != 0) return false;
        if (r[8] % (2 * r.den()) == 0) return false;
        break;

      default:
        throw CCTBX_INTERNAL_ERROR();
    }
    return true;
  }

} // namespace normalizer

}}} // namespace cctbx::sgtbx::reference_settings
