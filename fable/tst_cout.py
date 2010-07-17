from fable import cout
from libtbx.test_utils import \
  Exception_expected, show_diff, anchored_block_show_diff as absd
import libtbx.load_env
from cStringIO import StringIO
import os
op = os.path

def head_off(i): return i + 5
def tail_off(i): return -(i + 12) - 1

def exercise_simple(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/valid", test=op.isdir)
  def get(
        file_name,
        top_unit_name=None,
        data_specializations=True,
        arr_nd_size_max=None):
    if (verbose):
      print "exercise_simple:", file_name
    file_names = [op.join(t_dir, file_name)]
    common_report_stringio = StringIO()
    return cout.process(
      file_names=file_names,
      top_unit_name=top_unit_name,
      data_specializations=data_specializations,
      fem_do_safe=False,
      arr_nd_size_max=arr_nd_size_max,
      common_report_stringio=common_report_stringio)
  #
  assert not show_diff(get("add_reals.f"), """\
#include <fem.hpp> // Fortran EMulation library of fable module

namespace placeholder_please_replace {

using namespace fem::major_types;

using fem::common;

void
program_prog(
  int argc,
  char const* argv[])
{
  if (argc != 1) {
    throw std::runtime_error("Unexpected command-line arguments.");
  }
  float a = fem::float0;
  float b = fem::float0;
  float c = a + b;
}

} // namespace placeholder_please_replace

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    placeholder_please_replace::program_prog);
}
""")
  #
  assert not absd(get("add_real_integer.f"), tail_off(1), """\
  float a = fem::float0;
  int i = fem::int0;
  float c = a + i;
""")
  #
  assert not absd(get("logical_a_or_b.f"), tail_off(1), """\
  bool a = fem::bool0;
  bool b = fem::bool0;
  bool c = a || b;
""")
  #
  assert not absd(get("add_dp_integer.f"), tail_off(1), """\
  double a = fem::double0;
  int i = fem::int0;
  double c = a + i;
""")
  #
  assert not absd(get("add_strings.f"), tail_off(8), """\
  fem::str<3> a = "x\\"z";
  fem::str<4> b = "i\\\\'l";
  fem::str<7> c = a + b;
""")
  #
  lines = get("real_array_sum.f")
  assert not absd(lines, tail_off(1), """\
  arr<float> a(dimension(2), fem::fill0);
  float sum_a = a(1) + a(2);
""")
  lines = get("real_array_sum.f", arr_nd_size_max=2)
  assert not absd(lines, tail_off(2), """\
  arr_1d<2, float> a(fem::fill0);
""")
  #
  assert not absd(get("write_star.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  fem::str<1> c = "x";
  write(6, star), c;
  write(6, star), "i is zero.";
  bool l = fem::bool0;
  write(6, star), l;
  int i = fem::int0;
  write(6, star), i;
  fem::integer_star_8 j = fem::zero<fem::integer_star_8>();
  write(6, star), j;
  float r = fem::float0;
  write(6, star), r;
  double d = fem::double0;
  write(6, star), d;
  write(6, star), 1.e111;
  write(6, star), -1.e111;
  write(6, star);
  write(6, star), c, c, c, c, c, c, c, c, c, c, c, c;
  write(6, star), "i is ", "zero", ".";
  write(6, star), l, l, l, l, l, l, l, l, l, l, l, l;
  write(6, star), i, i, i, i, i, i, i, i, i, i, i, i;
  write(6, star), j, j, j, j, j, j, j, j, j, j, j, j;
  write(6, star), r, r, r, r, r, r, r, r, r, r, r, r;
  write(6, star), d, d, d, d, d, d, d, d, d, d, d, d;
  fem::str<1> s1 = "x";
  fem::str<2> s2 = "yz";
  write(6, star), s1, s1;
  write(6, star), s1, s2;
  write(6, star), s2, s1;
  write(6, star), s2, s2;
  write(6, star), s1, 12;
  write(6, star), s2, 34;
  write(6, star), 56, s1;
  write(6, star), 78, s2;
  write(6, star), "aBcD ", 12, " eFgHi ", 345;
""")
  #
  assert not absd(get("tab_syndrome.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int i = 1;
  write(6, star), i;
  i = 3;
  write(6, star), i;
  i = 5;
  write(6, star), i;
  i = 7;
  write(6, star), i;
  i = 123456;
  write(6, star), i;
""")
  #
  assert not absd(get("ops.f"), tail_off(1), """\
  bool la = fem::bool0;
  bool lb = !la;
  bool lc = la && lb;
  lc = la || lb;
  float a = fem::float0;
  float b = fem::float0;
  float c = a + b;
  c = a - b;
  c = a * b;
  c = a / b;
  fem::str<2> sa = "x";
  fem::str<3> sb = "abc";
  fem::str<5> sc = sa + sb;
""")
  #
  assert not absd(get("mod_integers.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  write(6, star), fem::mod(13, 5);
""")
  #
  assert not absd(get("do_enddo.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int i = fem::int0;
  FEM_DO(i, 1, 2) {
    write(6, star), i;
  }
  {
    int fem_do_last = 2 * 3;
    FEM_DO(i, 1, fem_do_last) {
      write(6, star), i;
    }
  }
  int j = fem::int0;
  FEM_DOSTEP(j, 3, 5, 2) {
    write(6, star), j;
  }
""")
  #
  assert not absd(get("if_endif.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int i = fem::int0;
  if (i == 0) {
    write(6, star), "i is zero.";
  }
""")
  #
  assert not absd(get("if_else_endif.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int i = fem::int0;
  FEM_DO(i, 0, 1) {
    if (i == 0) {
      write(6, star), "i is zero.";
    }
    else {
      write(6, star), "i is not zero.";
    }
  }
""")
  #
  assert not absd(get("if_elseif_else_endif.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int i = fem::int0;
  FEM_DO(i, 0, 2) {
    if (i == 0) {
      write(6, star), "i is zero.";
    }
    else if (i == 1) {
      write(6, star), "i is one.";
    }
    else {
      write(6, star), "i is not zero or one.";
    }
  }
""")
  #
  assert not absd(get("if_write.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int i = fem::int0;
  FEM_DO(i, 3, 4) {
    if (i < 4) {
      write(6, star), "i is less than four.";
    }
    if (i >= 4) {
      write(6, star), "i is greater than or equal to four.";
    }
  }
""")
  #
  assert not absd(get("assign_to_array_elements.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  arr<float> abc(dimension(2), fem::fill0);
  abc(1) = 10;
  abc(2) = 20;
  write(6, star), abc(1), abc(2);
""")
  #
  assert not absd(get("parameter_n.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int i = fem::int0;
  const int n = 2;
  FEM_DO(i, 1, n) {
    write(6, star), i;
  }
""")
  #
  assert not absd(get("do_ix_num.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int ix = fem::int0;
  const int num = 2;
  FEM_DO(ix, 1, num) {
    write(6, star), ix;
  }
""")
  #
  assert not absd(get("scopes_1.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int ix = fem::int0;
  int ix_sum = fem::int0;
  FEM_DO(ix, 1, 2) {
    ix_sum += ix;
  }
  write(6, star), ix_sum;
""")
  #
  assert not absd(get("scopes_2.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int ix = fem::int0;
  int ix_sum = fem::int0;
  FEM_DO(ix, 1, 3) {
    if (ix == 1) {
      write(6, star), "ix is one.";
    }
    else {
      ix_sum += ix;
    }
    write(6, star), ix_sum;
  }
  int ix_sum_sq = fem::int0;
  FEM_DO(ix, 2, 3) {
    ix_sum_sq += ix * ix;
  }
  write(6, star), ix_sum_sq;
""")
  #
  assert not absd(get("scopes_3.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  const int n_data = 10;
  float d_max = fem::float0;
  int i = fem::int0;
  arr<float> data(dimension(n_data), fem::fill0);
  float d = fem::float0;
  if (n_data <= 100) {
    write(6, star), "branch 1.";
  }
  else {
    d_max = 0;
    FEM_DO(i, 1, n_data) {
      d = data(i);
      if (d_max < d) {
        d_max = d;
      }
    }
  }
""")
  #
  assert not absd(get("scopes_4.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  const int n_data = 10;
  float d_max = fem::float0;
  int i = fem::int0;
  arr<float> data(dimension(n_data), fem::fill0);
  arr<float> d(dimension(1), fem::fill0);
  if (n_data <= 100) {
    write(6, star), "branch 1.";
  }
  else {
    d_max = 0;
    FEM_DO(i, 1, n_data) {
      d(1) = data(i);
      if (d_max < d(1)) {
        d_max = d(1);
      }
    }
  }
""")
  #
  assert not absd(get("arr_float_2.f"), tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int i = fem::int0;
  int j = fem::int0;
  arr<float, 2> data(dimension(2, 3), fem::fill0);
  FEM_DO(i, 1, 2) {
    FEM_DO(j, 1, 3) {
      data(i, j) = i + 10 * j;
    }
  }
  FEM_DO(j, 1, 3) {
    FEM_DO(i, 1, 2) {
      write(6, star), data(i, j);
    }
  }
""")
  lines = get("arr_float_2.f", arr_nd_size_max=6)
  assert not absd(lines, tail_off(11), """\
  arr_2d<2, 3, float> data(fem::fill0);
""")
  lines = get("arr_float_2.f", arr_nd_size_max=5)
  assert not absd(lines, tail_off(11), """\
  arr<float, 2> data(dimension(2, 3), fem::fill0);
""")
  #
  lines = get("subroutine_1.f")
  assert not absd(lines, head_off(3), """\
void
sub1(
  common& cmn)
{
  common_write write(cmn);
  write(6, star), "output from sub1.";
}

void
sub2(
  common& cmn)
{
  common_write write(cmn);
  write(6, star), "output from sub2.";
}
""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  common_write write(cmn);
  write(6, star), "first line in prog.";
  sub1(cmn);
  sub1(cmn);
  sub2(cmn);
  sub2(cmn);
  write(6, star), "last line in prog.";
""")
  #
  lines = get("subroutine_2.f")
  assert not absd(lines, head_off(3), expected="""\
void
sub1(
  common& cmn,
  int& i)
{
  common_write write(cmn);
  i = 3;
  write(6, star), "sub1", i;
  i = 7;
}

void
sub2(
  common& cmn,
  int& i,
  int& j)
{
  common_write write(cmn);
  j = 4;
  write(6, star), "sub2", i, j;
  i = 8;
  j = 5;
}
""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  common_write write(cmn);
  write(6, star), "first line in prog.";
  int i = fem::int0;
  sub1(cmn, i);
  write(6, star), "prog", i;
  int j = fem::int0;
  sub2(cmn, i, j);
  write(6, star), "prog", i, j;
  write(6, star), "last line in prog.";
""")
  #
  lines = get("subroutine_3.f")
  assert not absd(lines, head_off(3), """\
void
sub1(
  int& num)
{
  num = 9;
}

void
sub2(
  arr_ref<int> nums,
  int const& nums_size)
{
  nums(dimension(star));
  int i = fem::int0;
  FEM_DO(i, 1, nums_size) {
    nums(i) = i * 10;
  }
}

void
sub3(
  arr_ref<int> nums,
  int const& nums_size)
{
  nums(dimension(star));
  int i = fem::int0;
  FEM_DO(i, 1, nums_size) {
    nums(i) = i * 20;
  }
}
""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int num = fem::int0;
  sub1(num);
  write(6, star), "num after sub1:", num;
  int n = 1;
  sub2(num, n);
  write(6, star), "num after sub2", num;
  arr<int> nums(dimension(2), fem::fill0);
  sub1(nums);
  write(6, star), "nums after sub1:", nums(1), nums(2);
  n = 2;
  sub2(nums, n);
  write(6, star), "nums after sub2:", nums(1), nums(2);
  sub3(nums, n);
  write(6, star), "nums after sub3:", nums(1), nums(2);
""")
  #
  assert not absd(get("combine_decl_1.f"), tail_off(1), """\
  arr<int> vals(dimension(2), fem::fill0);
  write(6, star), vals(1);
  arr<float> abc(dimension(4), fem::fill0);
  write(6, star), abc(3);
""")
  #
  assert not absd(get("implied_program.f"), tail_off(0), """\
//C1
void
program_unnamed(
  int argc,
  char const* argv[])
{
  if (argc != 1) {
    throw std::runtime_error("Unexpected command-line arguments.");
  }
  common cmn;
  common_write write(cmn);
  int num = fem::int0;
  write(6, star), num;
}
//C2
""")
  #
  lines = get("implied_trailing_program.f")
  assert not absd(lines, head_off(3), expected="""\
void
sub(
  common& cmn)
{
  common_write write(cmn);
  write(6, star), "write sub";
}
""")
  assert not absd(lines, tail_off(1), """\
  sub(cmn);
""")
  #
  lines = get("common_0.f")
  assert not absd(lines, head_off(0), expected="""\

struct common_com
{
  int num;

  common_com() :
    num(fem::int0)
  {}
};

""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  common_write write(cmn);
  write(6, star), cmn.num;
""")
  #
  lines = get("common_1.f")
  assert not absd(lines, head_off(1), expected="""\
struct common_com
{
  int num;

  common_com() :
    num(fem::int0)
  {}
};

struct common :
  fem::common,
  common_com
{};

void
sub(
  common& cmn)
{
  common_write write(cmn);
  write(6, star), cmn.num;
}
""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  sub(cmn);
  cmn.num = 7;
  sub(cmn);
""")
  #
  lines = get("common_2.f")
  assert not absd(lines, head_off(1), expected="""\
struct common_com
{
  int num;
  float val;

  common_com() :
    num(fem::int0),
    val(fem::float0)
  {}
};

struct common :
  fem::common,
  common_com
{};

void
sub(
  common& cmn)
{
  common_write write(cmn);
  write(6, star), cmn.num;
  write(6, star), cmn.val;
}
""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  cmn.num = 3;
  cmn.val = 9;
  sub(cmn);
""")
  #
  lines = get("common_3.f")
  assert not absd(lines, head_off(1), expected="""\
struct common_com
{
  arr<int> vals;

  common_com() :
    vals(dimension(2), fem::fill0)
  {}
};

struct common :
  fem::common,
  common_com
{};

void
sub(
  common& cmn)
{
  // COMMON com
  arr_ref<int> vals(cmn.vals, dimension(2));
  //
  vals(1) = vals(2) + 1;
}
""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  common_write write(cmn);
  // COMMON com
  arr_ref<int> vals(cmn.vals, dimension(2));
  //
  vals(2) = 4;
  sub(cmn);
  write(6, star), vals(1);
""")
  #
  lines = get("common_4.f")
  assert not absd(lines, head_off(1), expected="""\
struct common_com
{
  arr<int> n;

  common_com() :
    n(dimension(2), fem::fill0)
  {}
};

struct common :
  fem::common,
  common_com
{};

void
sub(
  common& cmn,
  int const& num)
{
  // COMMON com
  arr_ref<int> n(cmn.n, dimension(2));
  //
  n(1) = num + 1;
  n(2) = num + 3;
}
""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  common_write write(cmn);
  // COMMON com
  arr_cref<int> n(cmn.n, dimension(2));
  //
  write(6, star), n(1), n(2);
  sub(cmn, 5);
  write(6, star), n(2), n(1);
""")
  #
  for file_name in ["save_0.f", "save_1.f"]:
    lines = get(file_name)
    assert not absd(lines, head_off(0), expected="""\

struct common :
  fem::common
{
  fem::cmn_sve program_prog_sve;
};

struct program_prog_save
{
  int num;

  program_prog_save() :
    num(fem::int0)
  {}
};

""")
    assert not absd(lines, tail_off(1), """\
  common cmn;
  FEM_CMN_SVE(program_prog);
  common_write write(cmn);
  write(6, star), sve.num;
""")
  #
  lines = get("save_2.f")
  assert not absd(lines, head_off(0), expected="""\

struct common :
  fem::common
{
  fem::cmn_sve sub_sve;
};

struct sub_save
{
  int num;

  sub_save() :
    num(fem::int0)
  {}
};

void
sub(
  common& cmn)
{
  FEM_CMN_SVE(sub);
  common_write write(cmn);
  // SAVE
  int& num = sve.num;
  //
  write(6, star), num;
  num++;
}

""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  sub(cmn);
  sub(cmn);
""")
  #
  lines = get("conv_recipe.f")
  assert not absd(lines, head_off(0), expected="""\

struct common_abc
{
  float a;
  float b;
  float c;

  common_abc() :
    a(fem::float0),
    b(fem::float0),
    c(fem::float0)
  {}
};

struct common :
  fem::common,
  common_abc
{
  fem::cmn_sve show_resolution_sve;
};

struct show_resolution_save
{
  float ass;
  float bss;
  float css;
  bool first;

  show_resolution_save() :
    ass(fem::float0),
    bss(fem::float0),
    css(fem::float0),
    first(fem::bool0)
  {}
};

//C cctbx_project/compcomm/newsletter09/conv_recipe.py, svn rev. 9983
//C
void
show_resolution(
  common& cmn,
  int const& h,
  int const& k,
  int const& l)
{
  FEM_CMN_SVE(show_resolution);
  common_write write(cmn);
  // COMMON abc
  float& a = cmn.a;
  float& b = cmn.b;
  float& c = cmn.c;
  //
  // SAVE
  float& ass = sve.ass;
  float& bss = sve.bss;
  float& css = sve.css;
  bool& first = sve.first;
  //
  if (is_called_first_time) {
    first = true;
  }
  if (first) {
    first = false;
    if (a <= 0 || b <= 0 || c <= 0) {
      write(0, "(1x,a)"), "invalid unit cell constants.";
      FEM_STOP(0);
    }
    ass = 1 / (a * a);
    bss = 1 / (b * b);
    css = 1 / (c * c);
  }
  float dss = h * h * ass + k * k * bss + l * l * css;
  if (dss == 0) {
    write(6, "(3(1x,i3),1x,a)"), h, k, l, "    infinity";
  }
  else {
    write(6, "(3(1x,i3),1x,f12.6)"), h, k, l, fem::sqrt(1 / dss);
  }
}

""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  cmn.a = 11.0f;
  cmn.b = 12.0f;
  cmn.c = 13.0f;
  show_resolution(cmn, 0, 0, 0);
  show_resolution(cmn, 1, 2, 3);
""")
  #
  assert not absd(get("read_star_integer.f"), tail_off(1), """\
  common_read read(cmn);
  common_write write(cmn);
  int num = fem::int0;
  read(5, star), num;
  write(6, star), num + 1;
""")
  #
  assert not absd(get("write_implied_do_1.f"), tail_off(1), """\
  int i = fem::int0;
  {
    write_loop wloop(cmn, 6, "(2i4)");
    FEM_DO(i, 1, 2) {
      wloop, nums(i);
    }
  }
  int j = fem::int0;
  {
    write_loop wloop(cmn, 6, "(3i4)");
    FEM_DO(i, 1, 2) {
      wloop, nums(i);
      FEM_DO(j, 1, 2) {
        wloop, nums(j);
      }
    }
  }
  {
    write_loop wloop(cmn, 6, "(5i4)");
    FEM_DO(i, 1, 2) {
      wloop, nums(i);
      FEM_DO(j, 1, 2) {
        wloop, nums(j);
      }
      FEM_DO(j, 1, 2) {
        wloop, nums(j);
      }
    }
  }
  int k = fem::int0;
  {
    write_loop wloop(cmn, 6, "(7i4)");
    FEM_DO(i, 1, 2) {
      wloop, nums(i);
      FEM_DO(j, 1, 2) {
        wloop, nums(j);
        FEM_DO(k, 1, 2) {
          wloop, nums(k);
        }
      }
    }
  }
""")
  #
  lines = get("write_extra_parentheses.f")
  assert not absd(lines, tail_off(35), """\
int
i1(
  int const& i)
{
  int return_value = fem::int0;
  return_value = 7 * i;
  return return_value;
}

int
i2(
  int const& i,
  int const& j)
{
  int return_value = fem::int0;
  return_value = i * 8 + j;
  return return_value;
}

int
i3(
  int const& i,
  int const& j,
  int const& k)
{
  int return_value = fem::int0;
  return_value = i * 9 + j * 29 + k;
  return return_value;
}
""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  common_write write(cmn);
  int i = 3;
  int j = 4;
  int k = 5;
  int l = 6;
  int m = 7;
  write(6, "(1i5)"), -i * 3;
  write(6, "(1i5)"), -i * 3;
  write(6, "(2i5)"), -i * 3, -j * 4;
  write(6, "(3i5)"), -i * 3, -j * 4, -k * 7;
  write(6, "(3i5)"), -i * 3, -j * 4, -k * 7;
  write(6, "(3i5)"), -i * 3, -j * 4, -k * 7;
  write(6, "(4i5)"), -i * 3, -j * 4, -k * 7, -l * 8;
  write(6, "(4i5)"), -i * 3, -j * 4, -k * 7, -l * 8;
  write(6, "(4i5)"), -i * 3, -j * 4, -k * 7, -l * 8;
  write(6, "(5i5)"), -i * 3, -j * 4, -k * 7, -l * 8, -m * 9;
  write(6, "(4i5)"), i1(-i * 3), -j * 4, -k * 7, -l * 8;
  write(6, "(4i5)"), i1(-i * 3), -j * 4, -k * 7, -l * 8;
  write(6, "(4i5)"), i1(-i * 3), -j * 4, -k * 7, -l * 8;
  write(6, "(4i5)"), i1(-i * 3), -j * 4, -k * 7, -l * 8;
  write(6, "(4i5)"), i1(-i * 3), -j * 4, -k * 7, -l * 8;
  write(6, "(2i5)"), i1(-i * 3), i3(-j * 4, -k * 7, -l * 8);
  write(6, "(2i5)"), i1(-i * 3), i3(-j * 4, -k * 7, -l * 8);
  write(6, "(2i5)"), i1(-i * 3), i3(-j * 4, -k * 7, -l * 8);
""")
  #
  for file_name in [
        "function_two_returns_1.f",
        "function_two_returns_2.f"]:
    lines = get(file_name)
  assert not absd(lines, tail_off(14), """\
int
fun(
  int const& i)
{
  int return_value = fem::int0;
  if (i < 3) {
    return_value = i * 4;
    return return_value;
  }
  return_value = -i;
  return return_value;
}
""")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  common_write write(cmn);
  write(6, star), fun(2);
  write(6, star), fun(3);
""")
  #
  assert not absd(get("do_label.f"), tail_off(1), """\
  int num = fem::int0;
  FEM_DO(num, 1, 2) {
    write(6, star), num;
  }
  int num2 = fem::int0;
  FEM_DO(num, 1, 2) {
    FEM_DO(num2, 3, 4) {
      write(6, star), num, num2;
    }
  }
  FEM_DO(num, 1, 2) {
    FEM_DO(num2, 3, 4) {
      write(6, star), num, num2;
      write(6, star), num * 10, num2 * 10;
    }
  }
""")
  #
  assert not absd(get("label_format.f"), tail_off(1), """\
  int num = 102;
  write(6, "(2(i4,1x,2i5,/,'_^'))"), num, num + 1, num + 2, num + 3,
    num + 4, num + 5;
""")
  #
  assert not absd(get("write_internal_file.f"), tail_off(1), """\
  int num = -2;
  fem::str<5> buf = fem::char0;
  write(buf, "(i3)"), num;
  write(6, "('num = (',a,')')"), buf;
""")
  #
  assert not absd(get("open_chain.f"), tail_off(1), """\
  cmn.io.open(10, "fable_tmp_661de075");
  cmn.io.close(10);
  cmn.io.open(10, "fable_tmp_661de075")
    .form("formatted")
    .status("unknown");
  cmn.io.close(10);
  try {
    cmn.io.open(10, "fable_tmp_661de075")
      .access("sequential")
      .form("formatted")
      .status("new");
  }
  catch (fem::io_err const&) {
    goto statement_10;
  }
  goto statement_20;
  statement_10:
  write(6, "(a)"), "open err statement";
  statement_20:
  try {
    cmn.io.close(10)
      .status("keep");
  }
  catch (fem::io_err const&) {
    goto statement_30;
  }
  goto statement_40;
  statement_30:
  write(6, "(a)"), "close err statement";
  statement_40:
  write(6, "(a)"), "Done.";
""")
  #
  assert not absd(get("goto_spaghetti.f"), tail_off(1), """\
  int i = fem::int0;
  int j = fem::int0;
  write(6, star), "start";
  goto statement_20;
  statement_10:
  i = 3;
  write(6, star), "stmt 10";
  goto statement_30;
  statement_20:
  write(6, star), "stmt 20";
  i = 2;
  statement_30:
  write(6, star), "stmt 30", i;
  if (i == 2) {
    goto statement_10;
  }
  FEM_DO(j, 1, 2) {
    if (j == 2) {
      goto statement_40;
    }
    write(6, star), "loop j is", j;
    statement_40:;
  }
  write(6, star), "end";
""")
  #
  lines = get("sub_nums_size_capacity.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  common& cmn,
  arr_cref<int> nums,
  int const& nums_size,
  int const& nums_capacity)
{
  nums(dimension(nums_capacity));
  common_write write(cmn);
  int i = fem::int0;
  FEM_DO(i, 1, nums_size) {
    write(6, "(i3)"), nums(i);
  }
}
""")
  assert not absd(lines, tail_off(1), """\
  arr<int> nums(dimension(3), fem::fill0);
  nums(1) = 12;
  nums(2) = 34;
  sub(cmn, nums, 2, 3);
""")
  #
  lines = get("passing_arrays.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  common& cmn,
  int const& num,
  arr_cref<int> nums1,
  arr_cref<int, 2> nums2)
{
  nums1(dimension(6));
  nums2(dimension(2, 3));
  common_write write(cmn);
  write(6, "(i1)"), num;
  int i = fem::int0;
  FEM_DO(i, 1, 6) {
    write(6, "(i2)"), nums1(i);
  }
  int j = fem::int0;
  FEM_DO(i, 1, 2) {
    FEM_DO(j, 1, 3) {
      write(6, "(i2)"), nums2(i, j);
    }
  }
}
""")
  assert not absd(lines, tail_off(1), """\
  arr<int> nums(dimension(6), fem::fill0);
  nums(1) = 3;
  int i = fem::int0;
  FEM_DO(i, 2, 6) {
    nums(i) = nums(i - 1) + i;
  }
  sub(cmn, nums, nums, nums);
""")
  #
  lines = get("unused_args.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  int const& /* num */,
  arr_cref<int> /* nums1 */,
  arr_cref<int, 2> /* nums2 */)
{
}
""")
  assert not absd(lines, tail_off(1), """\
  sub(1, 2, 3);
""")
  #
  lines = get("passing_arrays_2.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  arr_ref<int> nums1,
  arr_ref<int, 2> nums2)
{
  nums1(dim1(0, 1));
  nums2(dim1(2, 4).dim2(-1, 2));
  nums1(0) = 23;
  nums1(1) = 45;
  nums2(4, -1) = 67;
  int i = fem::int0;
  int j = fem::int0;
  FEM_DO(i, 2, 4) {
    FEM_DO(j, 0, 2) {
      nums2(i, j) = i * 10 + j;
    }
  }
}
""")
  assert not absd(lines, tail_off(1), """\
  arr<int> nums(dimension(12), fem::fill0);
  sub(nums, nums);
  int i = fem::int0;
  {
    write_loop wloop(cmn, 6, "(6i3)");
    FEM_DO(i, 1, 12) {
      wloop, nums(i);
    }
  }
""")
  #
  lines = get("array_origin.f")
  assert not absd(lines, head_off(14), """\
  arr<int> nums1(dim1(0, 1), fem::fill0);
  int k = fem::int0;
  arr<int, 2> nums2(dim1(0, 1).dim2(-1, 2), fem::fill0);
  int j = fem::int0;
  arr<int, 3> nums3(dim1(0, 1).dim2(3).dim3(-1, 2), fem::fill0);
""")
  lines = get("array_origin.f", arr_nd_size_max=24)
  assert not absd(lines, head_off(14), """\
  arr_1d<2, int> nums1(dim1(0, 1), fem::fill0);
  int k = fem::int0;
  arr_2d<2, 4, int> nums2(dim1(0, 1).dim2(-1, 2), fem::fill0);
  int j = fem::int0;
  arr_3d<2, 3, 4, int> nums3(dim1(0, 1).dim2(3).dim3(-1, 2), fem::fill0);
""")
  #
  lines = get("power.f")
  assert not absd(lines, tail_off(1), """\
  const float val = fem::pow(1.2f, 3.4f);
  write(6, "(f5.3)"), val;
  float x = 1.2f + fem::pow(3.4f, 5.6f) / 7.8f;
  write(6, "(f5.1)"), x;
  x = fem::pow((1.2f + 3.4f), 5.6f) / 7.8f;
  write(6, "(f5.1)"), x;
  x = fem::pow((1.2f + 3.4f), (5.6f / 7.8f));
  write(6, "(f5.3)"), x;
  x = -fem::pow2(1.3f);
  write(6, "(f5.2)"), x;
  x = fem::pow2((-1.3f));
  write(6, "(f4.2)"), x;
  x = ((-1.4f));
  write(6, "(f5.2)"), x;
  x = fem::pow3((-1.5f));
  write(6, "(f5.2)"), x;
  x = fem::pow4((-1.6f));
  write(6, "(f4.2)"), x;
""")
  #
  lines = get("stop_bare.f")
  assert not absd(lines, tail_off(2), """\
    if (i == 2) {
      FEM_STOP(0);
    }
    write(6, "(a,i2)"), "iteration", i;
""")
  #
  lines = get("stop_integer.f")
  assert not absd(lines, tail_off(2), """\
    if (i == 2) {
      FEM_STOP(2345);
    }
    write(6, "(a,i2)"), "iteration", i;
""")
  #
  lines = get("stop_string.f")
  assert not absd(lines, tail_off(2), """\
    if (i == 2) {
      FEM_STOP("Break");
    }
    write(6, "(a,i2)"), "iteration", i;
""")
  #
  lines = get("passing_strings.f")
  assert not absd(lines, head_off(3), """\
void
sub1(
  common& cmn,
  str_cref str)
{
  common_write write(cmn);
  write(6, "(i1)"), fem::len(str);
  write(6, "(a)"), str;
}
""")
  assert not absd(lines, tail_off(1), """\
  fem::str<2> str2 = "Pq";
  fem::str<3> str3 = "rSt";
  sub1(cmn, str2);
  sub1(cmn, str3);
  sub2(cmn, "a");
""")
  #
  lines = get("write_internal_file_2.f")
  assert not absd(lines, tail_off(1), """\
  arr<int> nums(dimension(2), fem::fill0);
  nums(1) = -2;
  nums(2) = 3;
  fem::str<8> buf = fem::char0;
  int i = fem::int0;
  {
    write_loop wloop(buf, "(2i3)");
    FEM_DO(i, 1, 2) {
      wloop, nums(i);
    }
  }
  write(6, "('nums = (',a,')')"), buf;
""")
  #
  lines = get("open_write_read.f")
  assert not absd(lines, head_off(19), """\
  cmn.io.open(10, buf)
    .form("formatted");
  int num = fem::int0;
  read(10, "(i6)"), num;
  cmn.io.close(10);
""")
  #
  lines = get("write_internal_file_3.f")
  assert not absd(lines, tail_off(1), """\
  int num = -2;
  fem::str<6> buf = "AbCdEf";
  write(buf(2, 4), "(i3)"), num;
  write(6, "('num = (',a,')')"), buf;
""")
  #
  lines = get("open_write_read_2.f")
  assert not absd(lines, tail_off(1), """\
  cmn.io.open(10, "fable_tmp_7895777d")
    .form("unformatted")
    .status("unknown");
  write(10, fem::unformatted), -123;
  cmn.io.close(10);
  cmn.io.open(10, "fable_tmp_7895777d")
    .form("unformatted")
    .status("old");
  int num = fem::int0;
  read(10, fem::unformatted), num;
  cmn.io.close(10);
  if (num !=  - 123) {
    write(6, "(a)"), "FAILURE int", num;
  }
  else {
    write(6, "(a)"), "OK";
  }
""")
  #
  lines = get("read_err.f")
  assert not absd(lines, head_off(15), """\
  try {
    read(10, "(i1)"), num;
  }
  catch (fem::io_err const&) {
    goto statement_10;
  }
  write(6, "(a)"), "FAILURE exercise_file_fmt";
  goto statement_20;
  statement_10:
  write(6, "(a)"), "success exercise_file_fmt";
  statement_20:;
""")
  #
  lines = get("read_end.f")
  assert not absd(lines, tail_off(6), """\
  try {
    read(10, "(i1)"), num1, num2;
  }
  catch (fem::read_end const&) {
    goto statement_10;
  }
""")
  #
  lines = get("read_end_err.f")
  assert not absd(lines, tail_off(6), """\
  try {
    read(10, "(i1)"), num1, num2;
  }
  catch (fem::read_end const&) {
    goto statement_10;
  }
  catch (fem::io_err const&) {
    goto statement_20;
  }
""")
  #
  lines = get("write_err.f")
  assert not absd(lines, tail_off(7), """\
  try {
    write(10, "(i1)"), num;
  }
  catch (fem::io_err const&) {
    goto statement_10;
  }
""")
  #
  lines = get("write_read_end_err_implied_do.f")
  assert not absd(lines, tail_off(25), """\
  try {
    write_loop wloop(cmn, 10, "(2i3)");
    FEM_DO(i, 8, 9) {
      wloop, i + 23;
    }
  }
  catch (fem::io_err const&) {
    goto statement_10;
  }
  statement_10:
""")
  assert not absd(lines, tail_off(8), """\
  try {
    read_loop rloop(cmn, 10, "(2i3)");
    FEM_DO(i, 1, 2) {
      rloop, nums(i);
    }
  }
  catch (fem::read_end const&) {
    goto statement_20;
  }
  catch (fem::io_err const&) {
    goto statement_30;
  }
  statement_20:
  statement_30:
""")
  #
  lines = get("goto_last_do.f")
  assert not absd(lines, head_off(10), """\
  if (num == 1) {
    num = 2;
    goto statement_10;
  }
  num = 3;
  statement_10:
  FEM_DO(i, 1, num) {
    write(6, "(i2)"), i;
  }
""")
  #
  lines = get("inquire.f")
  assert not absd(lines, head_off(9), """\
  FEM_DOSTEP(i, fem::len(s), 1, -1) {
""")
  assert not absd(lines, head_off(33), """\
  cmn.io.inquire_unit(10)
    .name(cvar);
""")
  assert not absd(lines, tail_off(3), """\
  try {
    cmn.io.inquire_file("fable_tmp_5d70aa2a")
      .exist(lvar);
  }
  catch (fem::io_err const&) {
    goto statement_10;
  }
  write(6, "(a,l1,a)"), "(", lvar, ")";
  goto statement_20;
  statement_10:
""")
  #
  lines = get("data_type_star.f")
  assert not absd(lines, tail_off(1), """\
  fem::logical_star_1 l1 = false;
  fem::integer_star_2 i2 = 4;
  fem::integer_star_4 i4 = 8;
  fem::integer_star_8 i8 = fem::zero<fem::integer_star_8>();
  if (i2 * 2 == i4) {
    i8 = 16;
    if (i4 * 2 == i8) {
      write(6, "(a)"), "OK integers";
      l1 = true;
    }
  }
  if (!l1) {
    write(6, "(a)"), "FAILURE integers";
  }
  fem::real_star_4 r4 = 3.14f;
  fem::real_star_8 r8 = 6.28f;
  if (fem::abs(r4 * 2 - r8) < 1.e-5f) {
    write(6, "(a)"), "OK reals";
  }
  else {
    write(6, "(a)"), "FAILURE reals";
  }
""")
  #
  lines = get("integer_star_2_array.f")
  assert not absd(lines, tail_off(1), """\
  arr<fem::integer_star_2> nums(dimension(2), fem::fill0);
  nums(1) = 9;
  nums(2) = -6;
  int num_sum = nums(1) + nums(2);
  write(6, "(i1)"), num_sum;
""")
  #
  lines = get("power_2.f")
  assert not absd(lines, tail_off(8), """\
  vals(1) = 1.2f;
  vals(2) = fem::pow2(vals(1));
  vals(3) = fem::pow(2, vals(2));
  vals(4) = fem::pow(vals(2), vals(3));
""")
  #
  lines = get("subroutine_5.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  common& cmn,
  str_cref str)
{
  common_write write(cmn);
  if (str(1, 1) == " ") {
    write(6, "(a)"), "str starts with a blank";
  }
  else {
    write(6, "(a)"), "str does not start with a blank";
  }
}
""")
  #
  lines = get("string_compare.f")
  assert not absd(lines, tail_off(1), """\
  str2 = " y";
  write(6, "(a,2l1)"), "p", str2 == " y", str2 != " y";
  write(6, "(a,2l1)"), "q", str2 == " y ", str2 != " y ";
  write(6, "(a,2l1)"), "r", str2 == " yz", str2 != " yz";
  write(6, "(a,2l1)"), "s", " y" == str2, " y" != str2;
  //C
""")
  #
  lines = get("decl_before_if.f")
  assert not absd(lines, head_off(9), """\
  int num = fem::int0;
  if (num_max > 41) {
    num = 41;
  }
  else {
    num = num_max;
  }
""")
  #
  lines = get("const_analysis_1.f", top_unit_name="prog")
  assert not absd(lines, head_off(3), """\
void
sub1(
  int& num)
{
  num = 12;
}

void
sub2(
  int& num)
{
  sub1(num);
}
""")
  #
  lines = get("const_analysis_2.f", top_unit_name="prog")
  assert not absd(lines, head_off(3), """\
void
sub1(
  common& cmn,
  int& num)
{
  common_read read(cmn);
  read(5, "(i2)"), num;
}

void
sub2(
  common& cmn,
  int& num)
{
  sub1(cmn, num);
}
""")
  #
  lines = get("double_literal.f")
  assert not absd(lines, tail_off(1), """\
  write(6, star), 1.2e2, 3.4e2f, 5.6e2;
""")
  #
  lines = get("read_implied_do_1.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  common& cmn,
  arr_ref<int> nums)
{
  nums(dimension(star));
  common_read read(cmn);
  int i = fem::int0;
  {
    read_loop rloop(cmn, 5, "(2i3)");
    FEM_DO(i, 1, 2) {
      rloop, nums(i);
    }
  }
}
""")
  #
  lines = get("dim_with_parameter.f")
  assert not absd(lines, tail_off(4), """\
  const int num = 2;
  arr<int> nums1(dimension(num), fem::fill0);
  nums1(1) = 12;
  nums1(2) = 34;
  arr<int> nums2(dimension(num), fem::fill0);
""")
  #
  lines = get("data_10.f", data_specializations=False)
  assert not absd(lines, head_off(7), """\
struct program_prog_save
{
  arr<int> num;

  program_prog_save() :
    num(dimension(2), fem::fill0)
  {}
};
""")
  assert not absd(lines, tail_off(1), """\
  const int one = 1;
  if (is_called_first_time) {
    fem::data_values data((values, 12, 34));
    {
      int fem_do_last = one + one;
      FEM_DO(ind, 1, fem_do_last) {
        data, num(ind);
      }
    }
  }
  write(6, star), num(1), num(2);
""")
  #
  lines = get("data_16.f")
  assert not absd(lines, head_off(7), """\
struct program_unnamed_save
{
  static const int num = 2;

  arr<float> vals;

  program_unnamed_save() :
    vals(dimension(num), fem::fill0)
  {}
};

const int program_unnamed_save::num;
""")
  assert not absd(lines, tail_off(8), """\
  if (is_called_first_time) {
    fem::data((values, num*datum(1.2f))), vals;
  }
""")
  #
  lines = get("data_22.f")
  assert not absd(lines, head_off(7), """\
struct program_prog_save
{
  arr<int> nums;
  arr<fem::str<2> > s2s;

  program_prog_save() :
    nums(dimension(2), fem::fill0),
    s2s(dimension(2), fem::fill0)
  {}
};
""")
  assert not absd(lines, tail_off(5), """\
    fem::data_values data((values, 12, "Xy", 34, "Za"));
    FEM_DO(i, 1, 2) {
      data, nums(i), s2s(i);
    }
""")
  #
  lines = get("data_23.f")
  assert not absd(lines, tail_off(7), """\
    fem::data_values data((values, 12, "Xy", 34, "Za", 56, "cD", 78, "eF"));
    FEM_DO(j, 1, 2) {
      FEM_DO(i, 1, 2) {
        data, nums(i, j), s2s(i, j);
      }
    }
""")
  #
  lines = get("data_24.f")
  assert not absd(lines, tail_off(7), """\
    fem::data_values data((values, 12, "Xy", 34, "Za", 56, "cD", 78, "eF"));
    FEM_DO(j, 1, 2) {
      FEM_DOSTEP(i, 1, 2, 2) {
        data, nums(i, j), s2s(i, j);
      }
    }
    FEM_DO(j, 1, 2) {
      FEM_DOSTEP(i, 2, 2, 2) {
        data, nums(i, j), s2s(i, j);
      }
    }
""")
  #
  lines = get("data_25.f", data_specializations=False)
  assert not absd(lines, tail_off(2), """\
  const fem::str<2> s12 = "xy";
  const fem::str<2> s34 = "ab";
  if (is_called_first_time) {
    fem::data((values, s12)), s4(1, 2);
    fem::data((values, s34)), s4(3, 4);
  }
""")
  #
  lines = get("data_26.f", data_specializations=False)
  assert not absd(lines, head_off(55), """\
    fem::data((values, 1, 2, 3)), num1, num2, num3;
""")
  #
  lines = get("data_27.f", data_specializations=False)
  assert not absd(lines, tail_off(3), """\
    fem::data((values, "A")), s2s(1)(1, 1);
    fem::data((values, "b")), s2s(2)(2, 2);
    fem::data((values, "C")), s2s(1)(2, 2);
    fem::data((values, "d")), s2s(2)(1, 1);
""")
  #
  lines = get("data_28.f", data_specializations=False)
  assert not absd(lines, tail_off(11), """\
    fem::data_values data((values, 1, 2, 3, 4));
    FEM_DO(i, 1, 2) {
      data, nums1(i);
    }
    FEM_DO(i, 1, 2) {
      data, nums2(i);
    }
""")
  #
  lines = get("data_29.f", data_specializations=False)
  assert not absd(lines, tail_off(3), """\
    {
      fem::data_values data((values, 12, 34));
      FEM_DOSTEP(i, 1, 4, 2) {
        data, nums(i);
      }
    }
    {
      fem::data_values data((values, 56, 78));
      FEM_DOSTEP(i, 2, 4, 2) {
        data, nums(i);
      }
    }
""")
  #
  lines = get("parameter_save.f")
  assert not absd(lines, tail_off(1), """\
  FEM_CMN_SVE(program_prog);
  common_write write(cmn);
  const int size = 2;
  arr<int> nums_local(dimension(size), fem::fill0);
  write(6, star), nums_local;
  arr_cref<int> nums_save(sve.nums_save, dimension(size));
  write(6, star), nums_save;
""")
  #
  lines = get("function_write.f")
  assert not absd(lines, tail_off(1), """\
  FEM_DO(i, 1, 3) {
    num = fun(cmn, num);
  }
""")
  #
  lines = get("parameters_recursive.f")
  assert not absd(lines, tail_off(5), """\
  const int num1 = 2;
  const int num2 = num1 + 3;
  arr<int> nums(dimension(num2), fem::fill0);
""")
  #
  lines = get("strings_size_dim_data.f")
  assert not absd(lines, head_off(7), """\
struct program_prog_save
{
  static const int base = 4;
  static const int size = base - 1;
  static const int dim = base - 2;

  arr<fem::str<size> > strings;

  program_prog_save() :
    strings(dimension(dim), fem::fill0)
  {}
};

const int program_prog_save::base;
const int program_prog_save::size;
const int program_prog_save::dim;
""")
  #
  lines = get("format_escape.f")
  assert not absd(lines, tail_off(1), """\
  write(6, "('Text with \\"quote\\" \\\\ and \\\\ backslashes.')");
""")
  #
  lines = get("main_cmn_indirect.f")
  assert not absd(lines, tail_off(1), """\
  common cmn;
  sub_main(cmn);
""")
  lines = get("main_cmn_indirect.f", top_unit_name="sub_main")
  assert not absd(lines, -4, """\
  write(6, "(a)"), "sub_main";
""")
  #
  lines = get("subroutine_7.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  int const& num1,
  fem::star_type const& /* UNHANDLED: star argument */,
  int const& num2,
  fem::star_type const& /* UNHANDLED: star argument */)
{
  if (num1 + num2 == 2) {
    return;
  }
  if (num1 == 3) {
    return;
  }
  if (num1 == 4) {
    return;
  }
  if (num1 == 5) {
    return;
  }
}
""")
  assert not absd(lines, tail_off(1), """\
  FEM_DO(i, 0, 5) {
    sub(i, star /* 10 UNHANDLED */, 1, star /* 20 UNHANDLED */);
    write(6, "(a)"), "regular";
    goto statement_30;
    write(6, "(a)"), "goto 10";
    goto statement_30;
    write(6, "(a)"), "goto 20";
    statement_30:;
  }
""")
  #
  lines = get("parameters_recursive_2.f")
  assert not absd(lines, tail_off(1), """\
  const int num1 = 1;
  const int num2 = num1 + 1;
  arr<int> nums(dimension(num2), fem::fill0);
  write(6, star), nums, num1;
""")
  #
  lines = get("data_30.f", data_specializations=False)
  assert not absd(lines, tail_off(11), """\
  if (is_called_first_time) {
    fem::data((values, 12, 34)), nums;
  }
  arr<float> data(dimension(2), fem::fill0);
""")
  #
  lines = get("string_array.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  common& cmn,
  str_arr_cref<> strings)
{
  strings(dimension(star));
  common_write write(cmn);
  int i = fem::int0;
  FEM_DO(i, 1, 2) {
    write(6, "(a)"), strings(i);
  }
}
""")
  #
  lines = get("string_sub_array_passing.f")
  assert not absd(lines, head_off(3), """\
void
sub1(
  common& cmn,
  str_arr_cref<> strs1)
""")
  assert not absd(lines, head_off(16), """\
void
sub2(
  str_arr_ref<> strs1)
""")
  #
  lines = get("write_pow.f")
  assert not absd(lines, tail_off(1), """\
  write(6, star), fem::pow3(2);
""")
  #
  lines = get("subroutine_8.f")
  assert not absd(lines, head_off(7), """\
struct sub_save
{
  int i;

  sub_save() :
    i(fem::int0)
  {}
};
""")
  assert not absd(lines, head_off(16), """\
void
sub(
  common& cmn,
  int const& sz,
  arr_cref<int> nums)
{
  FEM_CMN_SVE(sub);
  nums(dimension(sz));
  common_write write(cmn);
  // SAVE
  int& i = sve.i;
  //
  FEM_DO(i, 1, sz) {
    write(6, star), nums(i);
  }
}
""")
  #
  lines = get("subroutine_9.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  common& cmn,
  int const& n2,
  arr_cref<int, 2> nums)
{
  const int n1 = 2;
  nums(dimension(n1, n2));
  common_write write(cmn);
  write(6, star), nums;
}
""")
  #
  lines = get("parameter_save_common.f")
  assert not absd(lines, head_off(1), """\
struct common_cmn
{
  static const int ld = 2;

  arr<int> nums_cmn;

  common_cmn() :
    nums_cmn(dimension(ld), fem::fill0)
  {}
};

const int common_cmn::ld;

struct common :
  fem::common,
  common_cmn
{
  fem::cmn_sve program_prog_sve;
};

struct program_prog_save
{
  static const int ld = 2;

  arr<int> nums_sve;

  program_prog_save() :
    nums_sve(dimension(ld), fem::fill0)
  {}
};

const int program_prog_save::ld;
""")
  #
  lines = get("equivalence_01.f")
  assert not absd(lines, tail_off(5), """\
  local_equivalences loc_equivalences;
  {
    using fem::mbr; // member
    mbr<int> num1;
    mbr<int> num2;
    loc_equivalences.allocate(),
      equivalence(num1, num2)
        .align<1>()
         .with<2>()
    ;
  }
  int& num1 = loc_equivalences.bind<int>();
  int& num2 = loc_equivalences.bind<int>();
""")
  #
  lines = get("subroutine_write_iunit.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  common& cmn,
  int const& iunit)
{
  common_write write(cmn);
  write(iunit, "(a)"), "ABc";
}
""")
  #
  lines = get("common_variants.f")
  assert not absd(lines, head_off(1), """\
struct common :
  fem::common
{
  fem::variant_core common_scr;
  fem::cmn_sve sub1a_sve;
  fem::cmn_sve sub1b_sve;
  fem::cmn_sve sub2a_sve;
  fem::cmn_sve sub2b_sve;
  fem::cmn_sve sub3_sve;
  fem::cmn_sve sub4_sve;
};

struct sub1a_save
{
  fem::variant_bindings scr_bindings;
};

void
sub1a(
  common& cmn)
{
  FEM_CMN_SVE(sub1a);
  common_variant scr(cmn.common_scr, sve.scr_bindings);
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<int> i;
      mbr<int> j(dimension(2));
      scr.allocate(), i, j;
    }
  }
  int& i = scr.bind<int>();
  arr_ref<int> j(scr.bind<int>(), dimension(2));
  i = 12;
  j(1) = 34;
  j(2) = 65;
}
""")
  assert not absd(lines, tail_off(45), """\
  /* int const& i */ scr.bind<int>();
  int& j = scr.bind<int>();
""")
  assert not absd(lines, tail_off(21), """\
  /* arr_cref<int> i( */ scr.bind<int>() /* , dimension(2)) */ ;
  int& j = scr.bind<int>();
""")
  #
  lines = get("common_equivalence_1.f")
  assert not absd(lines, tail_off(7), """\
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<int> nums(dimension(2));
      mbr<int> numse(dimension(4));
      mbr<int> numx;
      scr.allocate(),
        equivalence(nums, numse)
          .align<1>()
           .with<2>(),
        numx
      ;
    }
  }
  arr_ref<int> nums(scr.bind<int>(), dimension(2));
  arr_ref<int> numse(scr.bind<int>(), dimension(4));
  int const& numx = scr.bind<int>();
""")
  #
  lines = get("common_equivalence_2.f")
  assert not absd(lines, tail_off(14), """\
      scr.allocate(),
        equivalence(nc, nl)
          .align<1>(arr_index(2))
           .with<2>(arr_index(1))
      ;
""")
  #
  lines = get("common_equivalence_3.f")
  assert not absd(lines, tail_off(19), """\
      mbr<int> nums3(dimension(4));
      mbr<int> nums1(dimension(5));
      mbr<int> nums2(dimension(3));
      scr.allocate(),
        equivalence(nums3, nums1, nums2)
          .align<2>(arr_index(2))
           .with<3>(arr_index(3))
           .with<1>(arr_index(4))
      ;
""")
  assert not absd(lines, tail_off(14), """\
  arr_ref<int> nums3(scr.bind<int>(), dimension(4));
  arr_ref<int> nums1(scr.bind<int>(), dimension(5));
  arr_ref<int> nums2(scr.bind<int>(), dimension(3));
""")
  #
  lines = get("common_equivalence_4.f")
  assert not absd(lines, tail_off(19), """\
      mbr<int> nums1(dimension(6));
      mbr<int> nums2(dimension(3));
      mbr<int> nums3(dimension(5));
      scr.allocate(),
        equivalence(nums1, nums2, nums3)
          .align<1>(arr_index(3))
           .with<2>(arr_index(2))
          .align<1>(arr_index(6))
           .with<3>(arr_index(4))
      ;
""")
  #
  lines = get("common_equivalence_5.f")
  assert not absd(lines, head_off(24), """\
      mbr<int> inside(dimension(2));
      mbr<int> data(dimension(4));
      scr.allocate(),
        equivalence(inside, data)
          .align<2>()
           .with<1>()
      ;
""")
  assert not absd(lines, head_off(34), """\
  arr_ref<int> data(scr.bind<int>(), dimension(4));
""")
  assert not absd(lines, tail_off(19), """\
      mbr<int> inside(dimension(3));
      scr.allocate(), inside;
""")
  #
  lines = get("equivalence_09.f")
  assert not absd(lines, head_off(8), """\
struct program_prog_save
{
  fem::variant_bindings scr_bindings;
  fem::variant_core_and_bindings save_equivalences;
};
""")
  assert not absd(lines, tail_off(6), """\
  common_variant scr(cmn.common_scr, sve.scr_bindings);
  save_equivalences sve_equivalences(sve.save_equivalences);
  if (is_called_first_time) {
    using fem::mbr; // member of variant common or equivalence
    {
      mbr<int> nc;
      mbr<int> nce;
      scr.allocate(),
        equivalence(nc, nce)
          .align<1>()
           .with<2>()
      ;
    }
    {
      mbr<int> ns;
      mbr<int> nse;
      sve_equivalences.allocate(),
        equivalence(ns, nse)
          .align<1>()
           .with<2>()
      ;
    }
  }
  local_equivalences loc_equivalences;
  {
    using fem::mbr; // member
    mbr<int> nl;
    mbr<int> nle;
    loc_equivalences.allocate(),
      equivalence(nl, nle)
        .align<1>()
         .with<2>()
    ;
  }
  int& nc = scr.bind<int>();
  int& nce = scr.bind<int>();
  int& ns = sve_equivalences.bind<int>();
  int& nse = sve_equivalences.bind<int>();
  int& nl = loc_equivalences.bind<int>();
  int& nle = loc_equivalences.bind<int>();
""")
  #
  lines = get("equivalence_10.f")
  assert not absd(lines, tail_off(9), """\
  local_equivalences loc_equivalences;
  {
    using fem::mbr; // member
    mbr<int> nl(dimension(2));
    mbr<int> nle(dimension(2));
    loc_equivalences.allocate(),
      equivalence(nl, nle)
        .align<1>()
         .with<2>()
    ;
  }
  arr_ref<int> nc(scr.bind<int>(), dimension(2));
  arr_ref<int> nce(scr.bind<int>(), dimension(2));
  arr_ref<int> ns(sve_equivalences.bind<int>(), dimension(2));
  arr_ref<int> nse(sve_equivalences.bind<int>(), dimension(2));
  arr_ref<int> nl(loc_equivalences.bind<int>(), dimension(2));
  arr_ref<int> nle(loc_equivalences.bind<int>(), dimension(2));
""")
  #
  lines = get("equivalence_05.f")
  assert not absd(lines, tail_off(8), """\
    loc_equivalences.allocate(),
      equivalence(s1, s2)
        .align<1>(str_index(1, 1))
         .with<2>(str_index(2, 2))
    ;
""")
  #
  lines = get("equivalence_06.f")
  assert not absd(lines, head_off(20), """\
    loc_equivalences.allocate(),
      equivalence(s1, s2, s3, s4)
        .align<1>(arr_index(1)(1, 1))
         .with<2>(arr_index(2)(3, 3))
""")
  assert not absd(lines, head_off(32), """\
  /* str_arr_ref<> s3( */ loc_equivalences.bind_str() /* , dimension(2)) */ ;
  /* str_ref s4 */ loc_equivalences.bind_str();
""")
  #
  lines = get("equivalence_repeated.f")
  assert not absd(lines, tail_off(12), """\
  const int itwo = 2;
  const int ione = 1;
  local_equivalences loc_equivalences;
  {
    using fem::mbr; // member
    mbr<int> nums1(dimension(2));
    mbr<int> nums2(dimension(itwo));
    loc_equivalences.allocate(),
      equivalence(nums1, nums2)
        .align<1>(arr_index(1))
         .with<2>(arr_index(ione))
        .align<1>(arr_index(ione))
         .with<2>(arr_index(1))
    ;
  }
""")
  #
  lines = get("equivalence_data.f", data_specializations=False)
  assert not absd(lines, head_off(8), """\
struct program_prog_save
{
  fem::variant_bindings scr_bindings;
};
""")
  assert not absd(lines, tail_off(3), """\
  if (is_called_first_time) {
    fem::data((values, 12, 34, 56)), numse;
  }
""")
  #
  lines = get("common_name_clash.f")
  assert not absd(lines, head_off(29), """\
void
sub1init(
  common& cmn)
{
  // COMMON cmn1
  int& num2 = static_cast<common_cmn1&>(cmn).num2;
  //
  cmn.num1 = 12;
  num2 = 34;
}

void
sub2init(
  common& cmn)
{
  // COMMON cmn2
  arr_ref<int> num2(static_cast<common_cmn2&>(cmn).num2, dimension(2));
  //
  num2(1) = 56;
  num2(2) = 78;
  cmn.num3 = 90;
}
""")
  #
  lines = get("external_arg_simple.f")
  assert not absd(lines, head_off(3), """\
typedef void (*show1_function_pointer)(common&, int const&);

void
show1(
  common& cmn,
  int const& i)
{
  common_write write(cmn);
  write(6, star), 10 + i;
}
""")
  assert not absd(lines, head_off(25), """\
void
show(
  common& cmn,
  show1_function_pointer which,
  int const& i)
{
  which(cmn, i);
}
""")
  #
  lines = get("external_arg_function.f")
  assert not absd(lines, head_off(30), """\
void
sub1(
  common& cmn,
  fun_function_pointer func)
{
  FEM_CMN_SVE(sub1);
  common_write write(cmn);
  // SAVE
  int& i = sve.i;
  //
  i = func(cmn, i);
  write(6, star), i;
}
""")
  #
  lines = get("function_calls_with_its_name_as_arg.f")
  assert not absd(lines, head_off(13), """\
int
ifun(
  common& cmn,
  int const& iarg)
{
  int return_value = fem::int0;
  return_value = 12;
  sub(cmn, return_value, iarg);
  return return_value;
}
""")
  assert not absd(lines, tail_off(1), """\
  write(6, star), ifun(cmn, 34);
""")
  #
  lines = get("function_no_arg.f")
  assert not absd(lines, tail_off(1), """\
  write(6, star), ifun();
  write(6, star), jfun(cmn);
""")
  #
  lines = get("if_arithmetic.f")
  assert not absd(lines, head_off(9), """\
  switch (fem::if_arithmetic(iarg - 2)) {
    case -1: goto statement_10;
    case  0: goto statement_20;
    default: goto statement_30;
  }
""")
  #
  lines = get("if_spaghetti.f")
  assert not absd(lines, head_off(22), """\
      statement_10:
      if (i == j) {
        goto statement_14;
      }
""")
  #
  lines = get("common_name_clash_2.f")
  assert not absd(lines, head_off(68), """\
  // COMMON cmn2
  arr_cref<int> num2(static_cast<common_cmn2&>(cmn).num2, dimension(2));
  int& num3 = cmn.num3;
  //
  int i = fem::int0;
  FEM_DO(i, 1, 2) {
    write(6, star), i, num2, num3;
  }
  FEM_DO(i, 3, 4) {
    write(6, star), i, num2, num3;
  }
""")
  #
  lines = get("dependency_cycle.f")
  assert not absd(lines, head_off(1), """\
/* Dependency cycles: 1
     sub1 sub2
 */
""")
  assert not absd(lines, head_off(7), """\
// forward declaration (dependency cycle)
void sub1(common&, int const&);
""")
  #
  lines = get("common_name_clash_3.f")
  assert not absd(lines, head_off(68), """\
  arr_cref<int> num2(static_cast<common_cmn2&>(cmn).num2, dimension(2));
  int& num3 = cmn.num3;
  //
  int j = fem::int0;
  int i = fem::int0;
  j = 0;
""")
  #
  lines = get("external_arg_non_const.f")
  assert not absd(lines, head_off(3), """\
typedef void (*exch_imp_function_pointer)(common&, arr_cref<int>, arr_ref<int>);

void
exch_imp(
  common& cmn,
  arr_cref<int> nc,
  arr_ref<int> nm)
""")
  #
  lines = get("parameter_for_arg_and_cmn_dim.f")
  assert not absd(lines, head_off(1), """\
struct common_scr
{
  static const int isz = 2;
""")
  assert not absd(lines, head_off(24), """\
  const int isz = 2;
  nums_arg(dimension(isz));
""")
  #
  lines = get("character_1_array_passing.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  common& cmn,
  str_arr_cref<> strs1)
{
  strs1(dimension(2));
  common_write write(cmn);
  write(6, star), strs1;
}
""")
  assert not absd(lines, tail_off(2), """\
  arr<fem::str<1> > strs1(dimension(2), fem::fill0);
  strs1(1) = "X";
  strs1(2) = "y";
""")
  #
  lines = get("do_variable_passed.f")
  assert not absd(lines, head_off(3), """\
void
sub(
  common& cmn,
  int& iarg)
{
  common_write write(cmn);
  FEM_DO(iarg, 1, 2) {
    write(6, star), iarg + 13;
  }
}
""")
  #
  lines = get("intrinsics_extra.f")
  assert not absd(lines, tail_off(2), """\
  fem::str<9> d = fem::char0;
  fem::date(d);
  write(6, "(a)"), d;
  fem::str<8> t = fem::char0;
  fem::time(t);
  write(6, "(a)"), t;
  fem::str<70> e = fem::char0;
  fem::getenv(" PATH ", e);
  write(6, "(a)"), e;
  float tm = fem::float0;
  fem::cpu_time(tm);
""")
  #
  lines = get("blockdata_unnamed.f", data_specializations=False)
  assert not absd(lines, head_off(21), """\
void
blockdata_unnamed(
  common& cmn)
{
  FEM_CMN_SVE(blockdata_unnamed);
  if (is_called_first_time) {
    fem::data((values, 3)), cmn.i;
  }
}
""")
  assert not absd(lines, tail_off(3), """\
  common cmn;
  blockdata_unnamed(cmn);
""")
  #
  lines = get("read_end_empty.f")
  assert not absd(lines, tail_off(8), """\
  statement_10:
  try {
    read(5, "()");
  }
""")
  #
  lines = get("rewind.f")
  assert not absd(lines, tail_off(18), """\
  cmn.io.rewind(1);
  read(1, "(i3)"), num;
  write(6, star), num;
  try {
    cmn.io.rewind(1);
  }
  catch (fem::io_err const&) {
    goto statement_10;
  }
  goto statement_20;
  statement_10:
""")
  #
  lines = get("read_rec_iostat.f")
  assert not absd(lines, tail_off(1), """\
  read(11, fem::unformatted).rec(21).iostat(ios), num;
  try {
    read(12, fem::unformatted).iostat(ios), num;
  }
  catch (fem::read_end const&) {
    goto statement_20;
  }
  catch (fem::io_err const&) {
    goto statement_10;
  }
  {
    read_loop rloop(cmn, 13, fem::unformatted);
    rloop.rec(23).iostat(ios);
    FEM_DO(i, 1, 2) {
      rloop, nums(i);
    }
  }
  try {
    read_loop rloop(cmn, 14, fem::unformatted);
    rloop.iostat(ios);
    FEM_DO(i, 1, 2) {
      rloop, nums(i);
    }
  }
  catch (fem::read_end const&) {
    goto statement_40;
  }
  catch (fem::io_err const&) {
    goto statement_30;
  }
  statement_10:
  statement_20:
  statement_30:
  statement_40:;
""")
  #
  lines = get("goto_computed.f")
  assert not absd(lines, tail_off(5), """\
      else {
        switch (i) {
          case 1: goto statement_10;
          case 2: goto statement_20;
          default: break;
        }
      }
      statement_10:
      write(6, star), "statement 10", j;
      goto statement_30;
      statement_20:
""")
  #
  lines = get("unformatted_experiments.f")
  assert not absd(lines, head_off(9), """\
  cmn.io.open(1, fem::file_not_specified)
    .form("unformatted")
    .status("unknown");
""")
  #
  lines = get("data_31.f", data_specializations=False)
  assert not absd(lines, tail_off(13), """\
    {
      fem::data_values data;
      data.values, 1, 2, 3, 4, 5, 6, 7, 8;
      data.values, 9, 10, 11, 12, 13, 14, 15, 16;
      data.values, 17, 18, 19, 20, 21, 22, 23, 24;
      data.values, 25, 26, 27, 28, 29, 30, 31, 32;
      data, nums1;
    }
    {
      fem::data_values data;
      data.values, 2, 3, 4, 5, 6, 7, 8, 9;
      data.values, 10, 11, 12, 13, 14, 15, 16, 17;
      data.values, 18, 19, 20, 21, 22, 23, 24, 25;
      data.values, 26, 27, 28, 29, 30, 31, 32, 33;
      data.values, 34;
      FEM_DO(i, 1, 33) {
        data, nums2(i);
      }
    }
""")
  #
  lines = get("data_32.f")
  assert not absd(lines, tail_off(12), """\
    num = -34;
    str = "YuIo";
    {
      static const int values[] = {
        +12, -34
      };
      fem::data_of_type<int>(FEM_VALUES_AND_SIZE),
        nums;
    }
    {
      static const int values[] = {
        -56, +78
      };
      fem::data_of_type<int> data(FEM_VALUES_AND_SIZE);
      FEM_DO(i, 1, 2) {
        data, numsi(i);
      }
    }
    {
      static const char* values[] = {
        "Cde", "FgH"
      };
      fem::data_of_type_str(FEM_VALUES_AND_SIZE),
        strs;
    }
    {
      static const char* values[] = {
        "IjkL", "MNOp"
      };
      fem::data_of_type_str data(FEM_VALUES_AND_SIZE);
      FEM_DO(i, 1, 2) {
        data, strsi(i);
      }
    }
    numj = 91;
    numsj(1) = 23;
    strsj(1) = "Hjklo";
    numsj(2) = 45;
    strsj(2) = "ASdfg";
""")
  assert not absd(lines, head_off(27), """\
    static const int values[] = {
      -24, +35
    };
    fem::data_of_type<int>(FEM_VALUES_AND_SIZE),
      nums;
""")
  #
  lines = get("const_expressions.f", arr_nd_size_max=6)
  assert not absd(lines, tail_off(10), """\
  arr_2d<n2 - 5, n3 - 48, int> nums1(fem::fill0);
""")
  assert not absd(lines, tail_off(3), """\
  const int n6 = fem::pow2(n1);
  arr<int> nums3(dimension(n6), fem::fill0);
""")
  lines = get("const_expressions.f", arr_nd_size_max=-6)
  assert not absd(lines, tail_off(10), """\
  arr_2d<n2 - 5, n3 - 48, int> nums1(fem::no_fill0);
""")
  #
  lines = get("common_save_members.f")
  assert not absd(lines, tail_off(9), """\
  // COMMON globals
  int& ci = cmn.ci;
  fem::str<8>& cc = cmn.cc;
  arr_ref<int> cai(cmn.cai, dimension(2));
  str_arr_ref<1> cas(cmn.cas, dimension(2));
  //
  // SAVE
  int& i = sve.i;
  arr_ref<int> sai(sve.sai, dimension(2));
  str_arr_ref<1> sas(sve.sas, dimension(2));
  fem::str<5>& sc = sve.sc;
  int& si = sve.si;
  //
  si = 12;
  ci = 34;
  sc = "WeRtY";
  cc = "uIoPqWeR";
  FEM_DO(i, 1, 2) {
    sai(i) = i + 37;
    cai(i) = i + 41;
  }
  sas(1) = "xYz";
  sas(2) = "EfG";
  cas(1) = "uvWx";
  cas(2) = "PqrS";
""")
  #
  lines = get("subroutine_4.f")
  assert not absd(lines, head_off(3), """\
//C1
//C c2
void
sub1(
  str_cref letter,
  int& num)
{
  //C3
  if (letter(1, 1) == "x") {
    num += 10;
  }
  //C4
}
//C c5

//C
//C6
void
sub2(
  str_cref letter,
  int& num)
{
  //C7
  sub1(letter, num);
  if (letter(1, 1) == "x") {
    num++;
  }
  else {
    num += 2;
  }
  //C8
}

//C
//C9
""")
  #
  lines = get("comments.f")
  assert not absd(lines, head_off(3), """\
//C
//C1
//Cc2
void
program_prog(
  int argc,
  char const* argv[])
{
  if (argc != 1) {
    throw std::runtime_error("Unexpected command-line arguments.");
  }
  common cmn;
  common_write write(cmn);
  int i = fem::int0;
  arr<int> nums(dimension(2), fem::fill0);
  //C3
  //C c4
  //C5
  //Cc6
  //C7
  //C c8
  FEM_DO(i, 1, 2) {
    //C9
    //Cc10
    //C
    //C12
    //C
    //C c13
    //Cc14
    //C
    //C c15
    nums(i) = i + 47;
    //C16
    //Cc17
  }
  //C
  //C c18
  try {
    write(6, star), nums;
  }
  catch (fem::io_err const&) {
    goto statement_10;
  }
  //C19
  //Cc20
  goto statement_20;
  //C21
  //C c22
  statement_10:
  FEM_STOP("write error");
  //C23
  //Cc24
  statement_20:;
  //C25
}
//C  c26
//C27
//C
""")
  #
  lines = get("long_lines.f")
  assert not absd(lines, tail_off(15), """\
  write(6, star), numbers(1), numbers(2), numbers(3), numbers(4),
    numbers(5), numbers(6), numbers(7), numbers(8);
  write(6, star), numbers(1), numbers(2), numbers(3), numbers(4),
    numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
    numbers(10);
  write(6, star), numbers(1), numbers(2), numbers(3), numbers(4),
    numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
    numbers(10), numbers(11), numbers(12);
  write(6, star), numbers(1), numbers(2), numbers(3), numbers(4),
    numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
    numbers(10), numbers(11), numbers(12), numbers(13), numbers(14);
  write(6, star), numbers(1), numbers(2), numbers(3), numbers(4),
    numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
    numbers(10), numbers(11), numbers(12), numbers(13), numbers(14),
    numbers(15), numbers(16);
  write(6, star), numbers(1), numbers(2), numbers(3), numbers(4),
    numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
    numbers(10), numbers(11), numbers(12), numbers(13), numbers(14),
    numbers(15), numbers(16), numbers(17), numbers(18);
  write(6, star), numbers(1), numbers(2), numbers(3), numbers(4),
    numbers(5), numbers(6), numbers(7), numbers(8), numbers(9),
    numbers(10), numbers(11), numbers(12), numbers(13), numbers(14),
    numbers(15), numbers(16), numbers(17), numbers(18), numbers(19),
    numbers(20);
  write(6, "(a)"),
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz)!@#$%^&*(`"
    "~-_+=[{]}\\\\|;:'\\",<.>/?";
  write(6, "(a)"),
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz)!@#$%^&*\\\\"
    "`~-_+=[{]}(|;:'\\",<.>/?";
  write(6, "(a)"),
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz)!@#$%^&*("
    "\\\\~-_+=[{]}`|;:'\\",<.>/?";
  fem::str<127> s =
    "qwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnmqwertyuioasdfghjkzxcv"
    "bnmqwerjkdfghjkertyjkxcghidfbndtyuiklmbvftyuiknbvdtyuh";
  write(6, "(a)"), s;
  write(6,
    "(/,' Sorry, your unit cell, range of hkl, size of map,',"
    "' and resolution will require ',/,' redimensioning of the program.',/,/,"
    "' This is quite easy:  You need to edit the source file ',"
    "' for the program',/,' and increase the value of \\"base_size\\" from ',i2,"
    "' to ',' a larger value',/,/,"
    "'  Then recompile the program and try again.',/,/,"
    "' If you do not have the source code, then you can obtain',/,"
    "' a version with a larger dimension from ',/,' our web site.',/)"),
    12;
  write(6,
    "(' first = ',f8.4,/,/,' second:              ',f5.2,/,"
    "' third:               ',f5.2,/,' fourth:              ',f5.2,/,"
    "' fifth:               ',f5.2,/,' sixth:               ',f5.1)"),
    1.2f, 3.4f, 5.6f, 7.8f, 9.1f, 2.3f;
""")
  assert not absd(lines, tail_off(1), """\
  if (nnnnn1 < 0 || nnnnn2 < 0 || nnnnn3 < 0 || nnnnn4 < 0 ||
      nnnnn5 < 0 || nnnnn6 <= 0) {
    write(6, "(a)"), "or ok";
  }
  if (nnnnn1 == 0 && nnnnn2 == 0 && nnnnn3 == 0 && nnnnn4 == 0 &&
      nnnnn5 == 0 && nnnnn6 <= 0) {
    write(6, "(a)"), "and ok";
  }
""")
  #
  lines = get("format_used_twice.f")
  assert not absd(lines, tail_off(1), """\
  static const char* format_10 = "(i2)";
  write(6, format_10), 12;
  write(6, "(i3)"), 345;
  write(6, format_10), 67;
""")
  #
  lines = get("blockdata_named.f")
  assert "\n".join(lines).find("Missing function implementation") < 0
  assert not absd(lines, head_off(1), """\
struct common_com
""")
  #
  lines = get("do_while.f")
  assert not absd(lines, tail_off(1), """\
  int i = 123;
  while (i < 169) {
    write(6, star), i;
    i += 45;
  }
  while (i < 281) {
    write(6, star), i;
    i += 67;
  }
""")

def exercise_syntax_error(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/syntax_error", test=op.isdir)
  from fable.read import Error
  def fail(file_name):
    if (verbose):
      print "exercise_syntax_error:", file_name
    cout.process(file_names=[op.join(t_dir, file_name)])
  try:
    fail("bad_open_err_label.f")
  except Error, e:
    assert str(e).startswith("Invalid statement label:")
    assert str(e).endswith("""\
  |      open(1, file=name, err=1.3)|
--------------------------------^""")
  else: raise Exception_expected
  try:
    fail("power_no_base.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      x = **3.4|
-------------^""")
  else: raise Exception_expected
  try:
    fail("power_no_exponent.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      x = 1.2**|
----------------^""")
  else: raise Exception_expected

def exercise_semantic_error(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/semantic_error", test=op.isdir)
  from fable import SemanticError
  def fail(file_name):
    if (verbose):
      print "exercise_semantic_error:", file_name
    cout.process(file_names=[op.join(t_dir, file_name)])
  try:
    fail("assignment_to_parameter.f")
  except SemanticError, e:
    assert str(e).startswith("Assignment to PARAMETER n:")
    assert str(e).endswith("""\
  |      n = 1|
---------^""")
  else: raise Exception_expected
  try:
    fail("inquire_no_unit_no_file.f")
  except SemanticError, e:
    assert str(e).startswith("Missing UNIT or FILE in INQUIRE statement:")
    assert str(e).endswith("""\
  |      inquire(exist=lexist)|
---------^""")
  else: raise Exception_expected
  try:
    fail("inquire_both_unit_and_file.f")
  except SemanticError, e:
    assert str(e).startswith(
      "Conflicting UNIT vs. FILE in INQUIRE statement"
      " (exactly one is needed):")
    assert str(e).endswith("""\
  |      inquire(10, file='fable_tmp')|
---------^""")
  else: raise Exception_expected
  try:
    fail("recursion_in_declaration.f")
  except SemanticError, e:
    assert str(e).startswith("Recursion in declaration:")
    assert str(e).endswith("""\
  |      dimension nums(nums)|
------------------------^""")
  else: raise Exception_expected

def exercise_unsupported(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/unsupported", test=op.isdir)
  def get(file_name):
    if (verbose):
      print "exercise_unsupported:", file_name
    return cout.process(file_names=[op.join(t_dir, file_name)])
  #
  assert not absd(get("goto_into_loop.f"), tail_off(1), """\
  int i = fem::int0;
  FEM_DO_SAFE(i, 1, 2) {
    statement_10:
    write(6, star), i;
  }
  goto statement_10;
""")

def exercise_dynamic_parameters(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/valid", test=op.isdir)
  def get(file_name, dynamic_parameters):
    if (verbose):
      print "exercise_dynamic_parameter:", file_name
    file_names = [op.join(t_dir, file_name)]
    return cout.process(
      file_names=file_names,
      top_unit_name="prog",
      dynamic_parameters=dynamic_parameters)
  #
  lines = get("dynamic_parameters_1.f", [
    cout.dynamic_parameter_props(
      name="root_size", ctype="int", default="3")])
  assert not absd(lines, head_off(0), """\

struct dynamic_parameters_t
{
  int root_size;

  dynamic_parameters_t(
    int argc,
    char const* argv[])
  :
    root_size(3)
  {
    fem::dynamic_parameters_from_argv(argc, argv, 1)
      .reset_if_given(root_size)
    ;
  }
};

struct common :
  fem::common
{
  dynamic_parameters_t dynamic_parameters;

  common(
    dynamic_parameters_t const& dynamic_parameters_)
  :
    dynamic_parameters(dynamic_parameters_)
  {}
};

void
sub(
  common& cmn,
  arr_ref<int> nums)
{
  const int root_size = cmn.dynamic_parameters.root_size;
""")
  assert not absd(lines, tail_off(7), """\
  common cmn(dynamic_parameters_t(argc, argv));
""")
  #
  lines = get("dynamic_parameters_2.f", [
    cout.dynamic_parameter_props(
      name="nums_size", ctype="int", default="2")])
  assert not absd(lines, head_off(17), """\
struct common_com
{
  const int nums_size;
  arr<int> nums;

  common_com(
    dynamic_parameters_t const& dynamic_parameters)
  :
    nums_size(dynamic_parameters.nums_size),
    nums(dimension(nums_size), fem::fill0)
  {}
};

struct common :
""")
  assert not absd(lines, head_off(39), """\
    common_com(dynamic_parameters_),
""")
  #
  lines = get("dynamic_parameters_3.f", [
    cout.dynamic_parameter_props(
      name="base_size", ctype="int", default="3")])
  assert not absd(lines, head_off(17), """\
struct common_com
{
  const int base_size;
  const int nums_size;
  arr<int> nums;

  common_com(
    dynamic_parameters_t const& dynamic_parameters)
  :
    base_size(dynamic_parameters.base_size),
    nums_size(base_size * 2),
    nums(dimension(nums_size), fem::fill0)
  {}
};

""")
  #
  lines = get("dynamic_parameters_4.f", [
    cout.dynamic_parameter_props(
      name="base_size", ctype="int", default="3")])
  assert not absd(lines, head_off(30), """\
struct sub_save
{
  const int base_size;
  arr<int> nums;

  sub_save(
    dynamic_parameters_t const& dynamic_parameters)
  :
    base_size(dynamic_parameters.base_size),
    nums(dimension(base_size * 2), fem::fill0)
  {}
};

""")
  assert not absd(lines, head_off(48), """\
  FEM_CMN_SVE_DYNAMIC_PARAMETERS(sub);
""")
  #
  lines = get("dynamic_parameters_5.f", [
    cout.dynamic_parameter_props(
      name="base_size", ctype="int", default="3")])
  assert not absd(lines, head_off(34), """\
  const int base_size = cmn.dynamic_parameters.base_size;
  nums(dimension(base_size * 2));
""")

def exercise_common_equivalence_simple(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/valid", test=op.isdir)
  def get(file_name, common_names, expected_common_report=None):
    if (verbose):
      print "exercise_common_equivalence_simple:", file_name
    file_names = [op.join(t_dir, file_name)]
    common_report_stringio = StringIO()
    lines = cout.process(
      file_names=file_names,
      top_unit_name="prog",
      common_equivalence_simple=set(common_names.split(",")),
      common_report_stringio=common_report_stringio)
    if (expected_common_report is None):
      assert common_report_stringio.getvalue() == ""
    else:
      assert not show_diff(
        common_report_stringio.getvalue(),
        expected_common_report)
    return lines
  #
  for i in [1,2]:
    lines = get("common_equivalence_simple_%d.f" % i, "info")
    assert not absd(lines, tail_off(2), """\
  common cmn;
  common_write write(cmn);
  // COMMON info
  arr_ref<int> nums(cmn.nums, dimension(2));
  //
  int& n1 = nums(1); // SIMPLE EQUIVALENCE
  n1 = 12;
  int& n2 = nums(2); // SIMPLE EQUIVALENCE
  n2 = 34;
""")
  #
  lines = get("common_equivalence_simple_3.f", "tab")
  assert not absd(lines, head_off(1), """\
struct common_tab
{
  int na;
  int nb_memory[2];
  int nc_memory[1-0+1];
  int nd_memory[(2-(-1)+1) * 3];

  arr_ref<int> nb;
  arr_ref<int> nc;
  arr_ref<int, 2> nd;

  common_tab() :
    na(fem::int0),
    nb(*nb_memory, dimension(2), fem::fill0),
    nc(*nc_memory, dim1(0, 1), fem::fill0),
    nd(*nd_memory, dim1(-1, 2).dim2(3), fem::fill0)
  {}
};
""")
  assert not absd(lines, tail_off(5), """\
  arr_ref<int> nums(cmn.na, dimension(17)); // SIMPLE EQUIVALENCE
""")
  #
  lines = get("common_equivalence_simple_4.f", "first",
    expected_common_report="""\
Name clash: n2 in COMMONs: first, second

""")
  assert not absd(lines, tail_off(6), """\
  arr_ref<int> nums(cmn.n1, dimension(3)); // SIMPLE EQUIVALENCE
""")
  assert not absd(lines, tail_off(2), """\
  int& m2 = n2; // SIMPLE EQUIVALENCE
""")
  #
  lines = get("common_equivalence_simple_5.f", "all")
  assert not absd(lines, tail_off(2), """\
  arr_ref<int> m1a(n1(1), dimension(2)); // SIMPLE EQUIVALENCE
  write(6, star), m1a;
  arr_ref<int> m1b(n1, dimension(2)); // SIMPLE EQUIVALENCE
  write(6, star), m1b;
  arr_ref<int> m1c(n1(1), dimension(2)); // SIMPLE EQUIVALENCE
  write(6, star), m1c;
  arr_cref<int> m2(n2, dimension(6)); // SIMPLE EQUIVALENCE
  write(6, star), m2;
  arr_cref<int> m2a(n2(1, 1), dimension(6)); // SIMPLE EQUIVALENCE
  write(6, star), m2a;
  arr_cref<int> m2b(n2, dimension(6)); // SIMPLE EQUIVALENCE
  write(6, star), m2b;
  arr_cref<int> m2c(n2(1, 1), dimension(6)); // SIMPLE EQUIVALENCE
""")
  #
  lines = get("common_equivalence_simple_6.f", "com")
  assert not absd(lines, head_off(3), """\
  fem::str<3> s3_memory[2];
  fem::str<8> s8;

  str_arr_ref<1> s3;
""")
  assert not absd(lines, tail_off(24), """\
  str_ref s6(s3, 6); // SIMPLE EQUIVALENCE
""")
  assert not absd(lines, tail_off(20), """\
  str_arr_ref<1> s2(s3, 2, dimension(3)); // SIMPLE EQUIVALENCE
""")
  assert not absd(lines, tail_off(10), """\
  str_ref s8e(cmn.s8, 8); // SIMPLE EQUIVALENCE
""")
  assert not absd(lines, tail_off(8), """\
  str_arr_ref<1> s4(cmn.s8, 4, dimension(2)); // SIMPLE EQUIVALENCE
""")
  assert not absd(lines, tail_off(4), """\
  str_arr_ref<1> s1(s3(2), 1, dimension(5)); // SIMPLE EQUIVALENCE
""")

def run(args):
  assert args in [[], ["--verbose"]]
  verbose = (len(args) != 0)
  exercise_simple(verbose=verbose)
  exercise_syntax_error(verbose=verbose)
  exercise_semantic_error(verbose=verbose)
  exercise_unsupported(verbose=verbose)
  exercise_dynamic_parameters(verbose=verbose)
  exercise_common_equivalence_simple(verbose=verbose)
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
