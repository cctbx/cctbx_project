from libtbx.test_utils import approx_equal
from libtbx.utils import Usage
from libtbx import easy_run
from libtbx import easy_pickle
from libtbx.path import full_command_path
import platform
import sys, os
op = os.path

__this_script__ = "cctbx_project/compcomm/newsletter09/sf_times.py"

fortran_template = r"""C %(this_script)s

      subroutine cos_wrapper(result, arg)
      REAL result
      REAL arg
      result = COS(arg)
      return
      end

      subroutine exp_wrapper(result, arg)
      REAL result
      REAL arg
      result = EXP(arg)
      return
      end

      subroutine sf(abcss, n_scatt, xyz, b_iso, n_refl, hkl, f_calc)
      implicit none
      REAL abcss(3)
      integer n_scatt
      REAL xyz(3, *)
      REAL b_iso(*)
      integer n_refl
      integer hkl(3, *)
      REAL f_calc(2, *)
      integer i_refl, i_scatt, j, h
      REAL phi, cphi, sphi, dss, ldw, dw, a, b
      DO i_refl=1,n_refl
        a = 0
        b = 0
        DO i_scatt=1,n_scatt
          phi = 0
          DO j=1,3
            phi = phi + hkl(j,i_refl) * xyz(j,i_scatt)
          enddo
          phi = phi * 2 * 3.1415926535897931
          call cos_wrapper(cphi, phi)
          call cos_wrapper(sphi, phi - 3.1415926535897931*0.5)
          dss = 0
          DO j=1,3
            h = hkl(j,i_refl)
            dss = dss + h*h * abcss(j)
          enddo
          ldw = -0.25 * dss * b_iso(i_scatt)
          call exp_wrapper(dw, ldw)
          a = a + dw * cphi
          b = b + dw * sphi
        enddo
        f_calc(1, i_refl) = a
        f_calc(2, i_refl) = b
      enddo
      return
      end

      program run
      implicit none
      REAL abcss(3)
      integer n_scatt
      parameter(n_scatt=%(n_scatt)s)
      REAL xyz(3, n_scatt)
      REAL b_iso(n_scatt)
      integer n_refl
      parameter(n_refl=%(n_refl)s)
      integer hkl(3, n_refl)
      REAL f_calc(2, n_refl)
      integer i, j, jr
      REAL a, b, max_a, max_b
      abcss(1) = 1/(11.0*11.0)
      abcss(2) = 1/(12.0*12.0)
      abcss(3) = 1/(13.0*13.0)
      jr = 0
      DO i=1,n_scatt
        DO j=1,3
          jr = mod(jr*1366+150889, 714025)
          xyz(j,i) = (mod(jr, 20000) - 10000) / 10000.0
        enddo
      enddo
      DO i=1,n_scatt
        jr = mod(jr*1366+150889, 714025)
        b_iso(i) = mod(jr, 10000) / 100.0
      enddo
      if (n_scatt .le. 10) then
        DO i=1,n_scatt
          write(6, '(4(1x,f9.6))')
     &      xyz(1,i), xyz(2,i), xyz(3, i), b_iso(i)
        enddo
      endif
      DO i=1,n_refl
        DO j=1,3
          jr = mod(jr*1366+150889, 714025)
          hkl(j,i) = mod(jr, 10) - 5
        enddo
      enddo
      call sf(abcss, n_scatt, xyz, b_iso, n_refl, hkl, f_calc)
      if (n_refl .le. 100) then
        DO i=1,n_refl
          write(6, '(3(1x,i3),1x,f12.6,1x,f12.6)')
     &      hkl(1,i), hkl(2,i), hkl(3,i),
     &      f_calc(1,i), f_calc(2,i)
        enddo
      else
        max_a = 0
        max_b = 0
        DO i=1,n_refl
          a = f_calc(1,i)
          b = f_calc(2,i)
          if (max_a .lt. a) max_a = a
          if (max_b .lt. b) max_b = b
        enddo
        write(6, '(2(1x,f12.6))') max_a, max_b
      endif
      end
"""

cpp_template = r"""// %(this_script)s

#include <vector>
#include <cstdio>
#include <cmath>
#include <cstddef>

#define DO1(i,n) for(i=1;i<=n;i++)

template <typename T>
struct dim1
{
  std::vector<T> data;
  dim1(int n) : data(n) {}
  T& operator()(int i) { return data[i-1]; }
};

template <typename T>
struct dim2
{
  int n1;
  std::vector<T> data;
  dim2(int n1_, int n2) : n1(n1_), data(n1*n2) {}
  T& operator()(int i, int j) { return data[i-1+(j-1)*n1]; }
};

typedef dim2<int> int2d;
typedef dim1<float> real1d;
typedef dim2<float> real2d;

    void
    cos_wrapper(float& result, float const& arg)
    {
      result = std::cos(arg);
    }

    void
    exp_wrapper(float& result, float const& arg)
    {
      result = std::exp(arg);
    }

    void
    sf(real1d& abcss,
       int n_scatt, real2d& xyz, real1d& b_iso,
       int n_refl, int2d& hkl, real2d& f_calc)
    {
      int i_refl, i_scatt, j, h;
      float phi, cphi, sphi, dss, ldw, dw, a, b;
      DO1(i_refl, n_refl) {
        a = 0;
        b = 0;
        DO1(i_scatt, n_scatt) {
          phi = 0;
          DO1(j, 3) {
            phi = phi + hkl(j,i_refl) * xyz(j,i_scatt);
          }
          phi = phi * 2 * 3.1415926535897931f;
          cos_wrapper(cphi, phi);
          cos_wrapper(sphi, phi - 3.1415926535897931f*0.5f);
          dss = 0;
          DO1(j, 3) {
            h = hkl(j,i_refl);
            dss = dss + h*h * abcss(j);
          }
          ldw = -0.25f * dss * b_iso(i_scatt);
          exp_wrapper(dw, ldw);
          a = a + dw * cphi;
          b = b + dw * sphi;
        }
        f_calc(1, i_refl) = a;
        f_calc(2, i_refl) = b;
      }
    }

    int
    main()
    {
      real1d abcss(3);
      int n_scatt;
      n_scatt = %(n_scatt)s;
      real2d xyz(3, n_scatt);
      real1d b_iso(n_scatt);
      int n_refl;
      n_refl = %(n_refl)s;
      int2d hkl(3, n_refl);
      real2d f_calc(2, n_refl);
      int i, j, jr;
      float a, b, max_a, max_b;
      abcss(1) = 1/(11.0f*11.0f);
      abcss(2) = 1/(12.0f*12.0f);
      abcss(3) = 1/(13.0f*13.0f);
      jr = 0;
      DO1(i, n_scatt) {
        DO1(j, 3) {
          jr = (jr*1366+150889) %% 714025;
          xyz(j,i) = (jr %% 20000 - 10000) / 10000.0f;
        }
      }
      DO1(i, n_scatt) {
        jr = (jr*1366+150889) %% 714025;
        b_iso(i) = (jr %% 10000) / 100.0f;
      }
      if (n_scatt <= 10) {
        DO1(i, n_scatt) {
          std::printf(" %%9.6f %%9.6f %%9.6f %%9.6f\n",
            xyz(1,i), xyz(2,i), xyz(3, i), b_iso(i));
        }
      }
      DO1(i, n_refl) {
        DO1(j, 3) {
          jr = (jr*1366+150889) %% 714025;
          hkl(j,i) = jr %% 10 - 5;
        }
      }
      sf(abcss, n_scatt, xyz, b_iso, n_refl, hkl, f_calc);
      if (n_refl <= 100) {
        DO1(i, n_refl) {
          std::printf(" %%3d %%3d %%3d %%12.6f %%12.6f\n",
            hkl(1,i), hkl(2,i), hkl(3,i),
            f_calc(1,i), f_calc(2,i));
        }
      }
      else {
        max_a = 0.0f;
        max_b = 0.0f;
        DO1(i, n_refl) {
          a = f_calc(1,i);
          b = f_calc(2,i);
          if (max_a < a) max_a = a;
          if (max_b < b) max_b = b;
        }
        std::printf(" %%12.6f %%12.6f\n", max_a, max_b);
      }
      return 0;
    }
"""

def compare_with_cctbx_structure_factors(n_scatt, n_refl, output_lines):
  from cctbx import xray
  from cctbx import miller
  from cctbx import crystal
  from cctbx.array_family import flex
  crystal_symmetry = crystal.symmetry(
    unit_cell=(11,12,13,90,90,90),
    space_group_symbol="P1")
  scatterers = flex.xray_scatterer()
  miller_indices = flex.miller_index()
  f_calc = flex.complex_double()
  for line in output_lines:
    flds = line.split()
    assert len(flds) in [4,5]
    if (len(flds) == 4):
      x,y,z,b_iso = [float(s) for s in flds]
      scatterers.append(
        xray.scatterer(site=(x,y,z), b=b_iso, scattering_type="const"))
    else:
      miller_indices.append([int(s) for s in flds[:3]])
      f_calc.append(complex(float(flds[3]), float(flds[4])))
  assert scatterers.size() == n_scatt
  assert miller_indices.size() == n_refl
  xs = xray.structure(
    crystal_symmetry=crystal_symmetry,
    scatterers=scatterers)
  fc = miller_array = miller.set(
    crystal_symmetry=crystal_symmetry,
    indices=miller_indices,
    anomalous_flag=False).array(data=f_calc)
  fc2 = fc.structure_factors_from_scatterers(
    xray_structure=xs,
    algorithm="direct",
    cos_sin_table=False).f_calc()
  for f1,f2 in zip(fc.data(), fc2.data()):
    assert approx_equal(f1, f2, eps=1e-5)

def build_run(a_out, n_scatt, n_refl, build_cmd, check_max_a_b):
  if (op.isfile("a.out")):
    os.remove("a.out")
  assert not op.isfile("a.out")
  if (a_out is None):
    easy_run.fully_buffered(command=build_cmd).raise_if_errors_or_output()
  else:
    open("a.out", "w").write(a_out)
    os.chmod("a.out", 0755)
  assert op.isfile("a.out")
  run_cmd = "/usr/bin/time  -p ./a.out"
  buffers = easy_run.fully_buffered(command=run_cmd)
  assert len(buffers.stderr_lines) == 3
  if (n_scatt <= 10 and n_refl <= 100):
    assert len(buffers.stdout_lines) == n_scatt + n_refl
  else:
    assert len(buffers.stdout_lines) == 1
    max_a, max_b = [float(s) for s in buffers.stdout_lines[0].split()]
  if (check_max_a_b):
    if (n_scatt == 2000 and n_refl == 20000):
      assert approx_equal(max_a, 35.047157, eps=1e-4)
      assert approx_equal(max_b, 25.212738, eps=1e-4)
    elif (n_scatt == 100 and n_refl == 1000):
      assert approx_equal(max_a,  4.493645, eps=1e-4)
      assert approx_equal(max_b, 10.515532, eps=1e-4)
    elif (n_scatt <= 10 and n_refl <= 100):
      compare_with_cctbx_structure_factors(
        n_scatt=n_scatt,
        n_refl=n_refl,
        output_lines=buffers.stdout_lines)
    else:
      raise RuntimeError, (max_a, max_b)
  utime = float(buffers.stderr_lines[1].split()[1])
  return utime

def fortran_write_build_run(
      a_out, n_scatt, n_refl, real, build_cmd, replace_cos, replace_exp):
  if (a_out is None):
    this_script = __this_script__
    srctxt = fortran_template % vars()
    if (replace_cos):
      srctxt = srctxt.replace(
        "COS(arg)",
        "arg / (abs(arg)+1.0)")
    if (replace_exp):
      srctxt = srctxt.replace(
        "EXP(arg)",
        "max(0.0, 1.0 - arg*arg)")
    srctxt = srctxt.replace("REAL", real)
    open("tmp.f", "w").write(srctxt)
    if (not op.isfile("sf_f.f")):
      open("sf_f.f", "w").write(srctxt)
  return build_run(
    a_out=a_out,
    n_scatt=n_scatt,
    n_refl=n_refl,
    build_cmd=build_cmd+" tmp.f",
    check_max_a_b=(not (replace_cos or replace_exp)))

def cpp_write_build_run(
      a_out, n_scatt, n_refl, real, build_cmd, replace_cos, replace_exp):
  if (a_out is None):
    this_script = __this_script__
    srctxt = cpp_template % vars()
    if (replace_cos):
      srctxt = srctxt.replace(
        "std::cos(arg)",
        "arg / (std::abs(arg)+1.0f)")
    if (replace_exp):
      srctxt = srctxt.replace(
        "std::exp(arg)",
        "std::max(float(0), 1.0f - arg*arg)")
    srctxt = srctxt.replace("float", real)
    open("tmp.cpp", "w").write(srctxt)
    if (not op.isfile("sf_cpp.cpp")):
      open("sf_cpp.cpp", "w").write(srctxt)
  return build_run(
    a_out=a_out,
    n_scatt=n_scatt,
    n_refl=n_refl,
    build_cmd=build_cmd+" tmp.cpp",
    check_max_a_b=(not (replace_cos or replace_exp)))

def run_combinations(
      compiler_versions,
      all_utimes,
      a_out_archive,
      write_a_out_archive,
      n_scatt,
      n_refl,
      compiler_build_opts_list,
      real_list,
      write_build_run):
  for compiler,build_opts in compiler_build_opts_list:
    if (write_a_out_archive):
      have_compiler = (full_command_path(command=compiler) is not None)
      if (not have_compiler):
        compiler_version = "n/a"
      else:
        compiler_version = easy_run.fully_buffered(
          command=compiler+" --version",
          join_stdout_stderr=True).stdout_lines[0]
      compiler_versions.append(compiler_version)
      a_out = None
    build_cmd = " ".join([compiler, build_opts])
    print build_cmd
    utimes = []
    for real in real_list:
      print "  %s" % real
      for replace_cos in [False, True]:
        print "    replace_cos", replace_cos
        for replace_exp in [False, True]:
          print "      replace_exp", replace_exp
          sys.stdout.flush()
          a_out_key = (build_cmd, replace_cos, replace_exp)
          if (not write_a_out_archive):
            compiler_version, a_out = a_out_archive.get(a_out_key, (
              "n/a", None))
          if (compiler_version != "n/a"):
            utime = write_build_run(
              a_out=a_out,
              n_scatt=n_scatt,
              n_refl=n_refl,
              real=real,
              build_cmd=build_cmd,
              replace_cos=replace_cos,
              replace_exp=replace_exp)
            print "        %4.2f" % utime
            if (write_a_out_archive):
              a_out_archive[a_out_key] = (
                compiler_version,
                open("a.out", "rb").read())
          else:
            utime = -1.0
            print "        n/a"
          utimes.append(utime)
          sys.stdout.flush()
    all_utimes.append(utimes)

def gcc_static_is_available():
  tmp_cpp = r"""\
#include <iostream>
#include <cmath>
int
main()
{
  std::cout << static_cast<int>(std::cos(0.0)) << std::endl;
  return 0;
}
"""
  open("tmp.cpp", "w").write(tmp_cpp)
  buffers = easy_run.fully_buffered(command="g++ -O tmp.cpp")
  buffers.raise_if_errors_or_output()
  buffers = easy_run.fully_buffered(command="g++ -static -O tmp.cpp")
  if (len(buffers.stderr_lines) != 0):
    return False
  return True

def usage():
  raise Usage(
    "cctbx.python sf_times.py unit_test|quick|production|a_out_archive.pickle")

def run(args):
  if (len(args) != 1): usage()
  build_platform = platform.platform()
  build_node = platform.node()
  compiler_versions = []
  gcc_static = "-static "
  a_out_archive = {}
  write_a_out_archive = True
  if (args[0] == "unit_test"):
    n_scatt, n_refl = 10, 100
  elif (args[0] == "quick"):
    n_scatt, n_refl = 100, 1000
  elif (args[0] == "production"):
    n_scatt, n_refl = 2000, 20000
  elif (op.isfile(args[0])):
    n_scatt, n_refl, \
    build_platform, build_node, \
    compiler_versions, gcc_static, a_out_archive = easy_pickle.load(
      file_name=args[0])
    write_a_out_archive = False
  else:
    usage()
  if (write_a_out_archive and not gcc_static_is_available()):
    gcc_static = ""
  all_utimes = []
  if (1):
    run_combinations(
      compiler_versions,
      all_utimes,
      a_out_archive=a_out_archive,
      write_a_out_archive=write_a_out_archive,
      n_scatt=n_scatt,
      n_refl=n_refl,
      compiler_build_opts_list=[
        ("ifort", "-O"),
        ("gfortran", gcc_static+"-O -ffast-math"),
        ("g77", gcc_static+"-O -ffast-math")],
      real_list=["real*4", "real*8"],
      write_build_run=fortran_write_build_run)
  if (1):
    run_combinations(
      compiler_versions,
      all_utimes,
      a_out_archive=a_out_archive,
      write_a_out_archive=write_a_out_archive,
      n_scatt=n_scatt,
      n_refl=n_refl,
      compiler_build_opts_list=[
        ("icpc", "-static -O"),
        ("g++", gcc_static+"-O -ffast-math")],
      real_list=["float", "double"],
      write_build_run=cpp_write_build_run)
  if (write_a_out_archive and len(a_out_archive) != 0):
    print "Writing file: a_out_archive.pickle"
    easy_pickle.dump(
      file_name="a_out_archive.pickle",
      obj=(
        n_scatt,
        n_refl,
        build_platform,
        build_node,
        compiler_versions,
        gcc_static,
        a_out_archive))
  print
  print "current_platform:", platform.platform()
  print "current_node:", platform.node()
  print "build_platform:", build_platform
  print "build_node:", build_node
  print 'gcc_static: "%s"' % gcc_static
  for compiler_version in compiler_versions:
    print "compiler:", compiler_version
  print "n_scatt * n_refl: %d * %d" % (n_scatt, n_refl)
  for utimes in all_utimes:
    print " ".join(["%6.2f" % u for u in utimes])

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
