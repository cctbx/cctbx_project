#ifndef CCTBX_SGTBX_GROUP_CODES_H
#define CCTBX_SGTBX_GROUP_CODES_H

namespace cctbx { namespace sgtbx {

  namespace crystal_system {

    enum code
    {
      undefined,
      triclinic,
      monoclinic,
      orthorhombic,
      tetragonal,
      trigonal,
      hexagonal,
      cubic
    };

    const char* label(code const& c);

  } // namespace crystal_system

  namespace matrix_group {

    class code
    {
      public:
        code(crystal_system::code const& x, int l, int p, int m)
        : x_(x), l_(l), p_(p), m_(m)
        {}

        bool operator==(code const& rhs) const
        {
          return    x_ == rhs.x_
                 && l_ == rhs.l_
                 && p_ == rhs.p_
                 && m_ == rhs.m_;
        }

        bool operator!=(code const& rhs) const
        {
          return !(*this == rhs);
        }

        int index() const { return m_; }

        code point_group_type() const
        {
          return code(x_, l_ - p_, 0, m_ + p_);
        }

        code laue_group_type() const
        {
          return code(x_, 0, 0, m_ + l_);
        }

        crystal_system::code crystal_system() const
        {
          return x_;
        }

        const char* label() const;

      private:
        crystal_system::code x_;
        int l_;
        int p_;
        int m_;
    };

    static const code undefined   (crystal_system::undefined,     0,  0,  0);
    static const code unknown     (crystal_system::undefined,     0,  0,  1);
    static const code code_1      (crystal_system::triclinic,     1,  0,  2);
    static const code code_1b     (crystal_system::triclinic,     0,  0,  3);
    static const code code_2      (crystal_system::monoclinic,    2,  0,  4);
    static const code code_m      (crystal_system::monoclinic,    1,  0,  5);
    static const code code_2_m    (crystal_system::monoclinic,    0,  0,  6);
    static const code code_222    (crystal_system::orthorhombic,  2,  0,  7);
    static const code code_mm2    (crystal_system::orthorhombic,  1,  0,  8);
    static const code code_mmm    (crystal_system::orthorhombic,  0,  0,  9);
    static const code code_4      (crystal_system::tetragonal,    2,  0, 10);
    static const code code_4b     (crystal_system::tetragonal,    1,  0, 11);
    static const code code_4_m    (crystal_system::tetragonal,    0,  0, 12);
    static const code code_422    (crystal_system::tetragonal,    4,  0, 13);
    static const code code_4mm    (crystal_system::tetragonal,    3,  0, 14);
    static const code code_4b2m   (crystal_system::tetragonal,    2,  1, 15);
    static const code code_4bm2   (crystal_system::tetragonal,    1,  0, 16);
    static const code code_4_mmm  (crystal_system::tetragonal,    0,  0, 17);
    static const code code_3      (crystal_system::trigonal,      1,  0, 18);
    static const code code_3b     (crystal_system::trigonal,      0,  0, 19);
    static const code code_321    (crystal_system::trigonal,      8,  2, 20);
    static const code code_312    (crystal_system::trigonal,      7,  1, 21);
    static const code code_32     (crystal_system::trigonal,      6,  0, 22);
    static const code code_3m1    (crystal_system::trigonal,      5,  2, 23);
    static const code code_31m    (crystal_system::trigonal,      4,  1, 24);
    static const code code_3m     (crystal_system::trigonal,      3,  0, 25);
    static const code code_3bm1   (crystal_system::trigonal,      2,  2, 26);
    static const code code_3b1m   (crystal_system::trigonal,      1,  1, 27);
    static const code code_3bm    (crystal_system::trigonal,      0,  0, 28);
    static const code code_6      (crystal_system::hexagonal,     2,  0, 29);
    static const code code_6b     (crystal_system::hexagonal,     1,  0, 30);
    static const code code_6_m    (crystal_system::hexagonal,     0,  0, 31);
    static const code code_622    (crystal_system::hexagonal,     4,  0, 32);
    static const code code_6mm    (crystal_system::hexagonal,     3,  0, 33);
    static const code code_6b2m   (crystal_system::hexagonal,     2,  1, 34);
    static const code code_6bm2   (crystal_system::hexagonal,     1,  0, 35);
    static const code code_6_mmm  (crystal_system::hexagonal,     0,  0, 36);
    static const code code_23     (crystal_system::cubic,         1,  0, 37);
    static const code code_m3b    (crystal_system::cubic,         0,  0, 38);
    static const code code_432    (crystal_system::cubic,         2,  0, 39);
    static const code code_4b3m   (crystal_system::cubic,         1,  0, 40);
    static const code code_m3bm   (crystal_system::cubic,         0,  0, 41);

  } // namespace matrix_group

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_GROUP_CODES_H
