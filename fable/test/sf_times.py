from libtbx.test_utils import approx_equal
from libtbx.utils import Usage
from libtbx import easy_run
import libtbx.load_env
import platform
import time
import sys, os
op = os.path

__this_script__ = "cctbx_project/fable/test/sf_times.py"
       # based on  cctbx_project/compcomm/newsletter09/sf_times.py

setup_dir = "/net/cci/setup/Linux"

ifort_versions = ["intel121.sh", "intel111.sh", "ifort91.sh"]

icc_versions = [
  "intel121.sh",
  "intel111.sh",
  "icc101.sh",
  "icc91.sh"]

gcc_versions = [
  "gcc-4.6.1_fc8.sh",
  "gcc-4.5.3_fc8.sh",
  "gcc-4.4.6_fc8.sh",
  "gcc-4.3.6_fc8.sh",
  "gcc-4.2.4_fc8.sh"]

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

def build_run(
      setup_cmd, ld_preload_flag, n_scatt, n_refl, build_cmd, check_max_a_b):
  if (op.isfile("a.out")):
    os.remove("a.out")
  assert not op.isfile("a.out")
  print build_cmd
  buffers = easy_run.fully_buffered(command=build_cmd)
  msg = buffers.format_errors_if_any()
  if (msg is not None):
    if (0):
      print build_cmd
      print
      print msg
      print
      STOP()
    return None
  assert op.isfile("a.out")
  run_cmd = setup_cmd
  if (ld_preload_flag):
    run_cmd += 'env LD_PRELOAD='\
    '"/net/marbles/raid1/rwgk/dist/opt_resources/linux64/libimf.so:"'\
    '"/net/marbles/raid1/rwgk/dist/opt_resources/linux64/libirc.so" '
  run_cmd += '/usr/bin/time -p ./a.out'
  buffers = easy_run.fully_buffered(command=run_cmd)
  if (len(buffers.stderr_lines) != 3):
    print "v"*79
    print "\n".join(buffers.stderr_lines)
    print "^"*79
    raise RuntimeError(
      "Unexpected number of output lines"
      " (3 expected; acutal output see above).")
  if (n_scatt == 0):
    pass
  elif (n_scatt <= 10 and n_refl <= 100):
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
      if (libtbx.env.has_module(name="cctbx")):
        compare_with_cctbx_structure_factors(
          n_scatt=n_scatt,
          n_refl=n_refl,
          output_lines=buffers.stdout_lines)
    else:
      raise RuntimeError, (max_a, max_b)
  utime = float(buffers.stderr_lines[1].split()[1])
  return utime

def finalize_cpp_build_cmd(source_cpp):
  from fable import simple_compilation
  comp_env = simple_compilation.environment()
  return comp_env.assemble_include_search_paths(no_quotes=False) \
    + " " + source_cpp

def write_build_run(
      setup_cmd, ld_preload_flag, n_scatt, n_refl, real, lang, build_cmd,
      replace_cos, replace_exp):
  this_script = __this_script__
  for_txt = fortran_template % vars()
  if (replace_cos):
    for_txt = for_txt.replace(
      "COS(arg)",
      "arg / (abs(arg)+1.0)")
  if (replace_exp):
    for_txt = for_txt.replace(
      "EXP(arg)",
      "max(0.0, 1.0 - arg*arg)")
  for_txt = for_txt.replace("REAL", real)
  open("tmp.f", "w").write(for_txt)
  from fable import cout
  cpp_txt = cout.process(
    file_names=["tmp.f"],
    namespace="sf_test",
    fem_do_safe=False,
    inline_all=True)
  open("tmp.cpp", "w").write("\n".join(cpp_txt)+"\n")
  if (lang.lower() == "f"):
    build_cmd += " tmp.f"
  elif (lang.lower() == "c"):
    build_cmd += finalize_cpp_build_cmd("tmp.cpp")
  else:
    raise RuntimeError('Unknown lang: "%s"' % lang)
  return build_run(
    setup_cmd=setup_cmd,
    ld_preload_flag=ld_preload_flag,
    n_scatt=n_scatt,
    n_refl=n_refl,
    build_cmd=build_cmd,
    check_max_a_b=(not (replace_cos or replace_exp)))

def run_combinations(
      compiler_versions,
      all_utimes,
      n_scatt,
      n_refl,
      compiler_build_opts_list,
      real_list):
  for lang,setup_sh_list,compiler,build_opts in compiler_build_opts_list:
    for setup_sh in setup_sh_list:
      if (setup_sh is None):
        setup_cmd = ""
      else:
        setup_cmd = ". %s/%s; " % (setup_dir, setup_sh)
      compiler_version = easy_run.fully_buffered(
        command=setup_cmd+compiler+" --version",
        join_stdout_stderr=True).stdout_lines[0]
      if (lang in ["f", "c"]):
        ld_preload_flags = [False, True]
      else:
        ld_preload_flags = [False]
      for ld_preload_flag in ld_preload_flags:
        iml = ["", " Intel Math Lib"][int(ld_preload_flag)]
        compiler_versions.append(compiler_version + iml)
        build_cmd = " ".join([setup_cmd+compiler, build_opts])
        print build_cmd
        utimes = []
        if (n_scatt != 0):
          for real in real_list:
            print "  %s" % real
            for replace_cos in [False, True]:
              print "    replace_cos", replace_cos
              for replace_exp in [False, True]:
                print "      replace_exp", replace_exp
                sys.stdout.flush()
                if (compiler_version != "n/a"):
                  utime = write_build_run(
                    setup_cmd=setup_cmd,
                    ld_preload_flag=ld_preload_flag,
                    n_scatt=n_scatt,
                    n_refl=n_refl,
                    real=real,
                    lang=lang,
                    build_cmd=build_cmd,
                    replace_cos=replace_cos,
                    replace_exp=replace_exp)
                  if (utime is not None):
                    print "        %4.2f" % utime
                  else:
                    utime = -1.0
                    print "        err"
                else:
                  utime = -1.0
                  print "        n/a"
                utimes.append(utime)
                sys.stdout.flush()
        else:
          if (lang.lower() == "f"):
            f_source = libtbx.env.find_in_repositories(
              relative_path="lapack_fem/dsyev_test.f",
              test=op.isfile,
              optional=False)
            build_cmd_compl = build_cmd + " " + f_source
          else:
            cpp_source = libtbx.env.find_in_repositories(
              relative_path="lapack_fem/dsyev_test.cpp",
              test=op.isfile,
              optional=False)
            build_cmd_compl = build_cmd + finalize_cpp_build_cmd(cpp_source)
          utime = build_run(
            setup_cmd=setup_cmd,
            ld_preload_flag=ld_preload_flag,
            n_scatt=n_scatt,
            n_refl=n_refl,
            build_cmd=build_cmd_compl,
            check_max_a_b=False)
          if (utime is None):
            print "err"
            utime = -1.0
          else:
            print "utime: %.2f" % utime
          utimes.append(utime)
        all_utimes.append((utimes, build_cmd + iml))

def usage():
  raise Usage("fable.python sf_times.py unit_test|quick|production")

def run(args):
  if (len(args) != 1): usage()
  t_start = time.time()
  build_platform = platform.platform()
  build_node = platform.node()
  compiler_versions = []
  if (args[0] == "unit_test"):
    n_scatt, n_refl = 10, 100
  elif (args[0] == "quick"):
    n_scatt, n_refl = 100, 1000
  elif (args[0] == "production"):
    n_scatt, n_refl = 2000, 20000
  elif (args[0] == "dsyev"):
    n_scatt, n_refl = 0, 0
    del icc_versions[0]
    del icc_versions[-1]
  else:
    usage()
  gcc_sh = gcc_versions + [None]
  icc_sh = icc_versions
  if (args[0] == "quick"):
    gcc_sh = gcc_sh[:2]
    icc_sh = icc_sh[:1]
  all_utimes = []
  run_combinations(
    compiler_versions,
    all_utimes,
    n_scatt=n_scatt,
    n_refl=n_refl,
    compiler_build_opts_list=[
      ("F", ifort_versions, "ifort", "-O"),
      ("f", gcc_sh, "gfortran", "-O -ffast-math"),
      ("f", gcc_sh, "gfortran", "-O -ffast-math -march=native"),
      ("C", icc_sh, "icpc", "-O"),
      ("c", gcc_sh, "g++", "-O -ffast-math"),
      ("c", gcc_sh, "g++", "-O -ffast-math -march=native"),
      ("c", [None], "clang++",
        "-O -U__GXX_WEAK__ -Wno-logical-op-parentheses -ffast-math"),
      ("c", [None], "clang++",
        "-O -U__GXX_WEAK__ -Wno-logical-op-parentheses -ffast-math"
        " -march=native")],
    real_list=["real*4", "real*8"])
  print
  print "current_platform:", platform.platform()
  print "current_node:", platform.node()
  print "build_platform:", build_platform
  print "build_node:", build_node
  for compiler_version in compiler_versions:
    print "compiler:", compiler_version
  if (n_scatt != 0):
    print "n_scatt * n_refl: %d * %d" % (n_scatt, n_refl)
    print '''\
"s" or "d": single-precision or double-precision floating-point variables
"E" or "e": using the library exp(arg) function or "max(0.0, 1.0 - arg*arg)"
"C" or "c": using the library cos(arg) function or "arg / (abs(arg)+1.0)"'''
    print "  sEC    seC    sEc    sec    dEC    deC    dEc    dec"
  else:
    print "dsyev times:"
  useful_utimes = []
  for utimes,build_cmd in all_utimes:
    if (max(utimes) != -1.0):
      print " ".join(["%6.2f" % u for u in utimes]), build_cmd
      useful_utimes.append((utimes,build_cmd))
  if (len(useful_utimes) > 1):
    print "Relative to first:"
    for utimes,build_cmd in useful_utimes:
      print " ".join(["%6.2f" % (u/max(u0,0.01))
        for u,u0 in zip(utimes,useful_utimes[0][0])]), build_cmd
  print "Wall clock time: %.2f s" % (time.time()-t_start)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
