from scitbx import matrix
from mmtbx.tls import tls_as_xyz

def print_step(m):
  print "-"*80
  tls_as_xyz.print_step(m)


def exercise_01():
  # getTLS_test001.mes
  T = matrix.sym(sym_mat3=[0.01, 0.16, 0.64, 0,0,0])
  L = matrix.sym(sym_mat3=[0.01, 0.04, 0.09, 0,0,0])
  S = matrix.sqr([0,0,0,0,0,0,0,0,0])
  print_step("Input (getTLS_test001.mes):")
  print "  T:\n", T
  print "  L:\n", L
  print "  S:\n", S
  tls_as_xyz.decompose_tls(T=T, L=L, S=S)

def exercise_02():
  # getTLS_test002.mes
  T = matrix.sym(sym_mat3=[0.325, 0.325, 0.16, -0.315, 0, 0])
  L = matrix.sym(sym_mat3=[0.01, 0.04, 0.09, 0,0,0])
  S = matrix.sqr([0,0,0,0,0,0,0,0,0])
  print_step("Input (getTLS_test002.mes):")
  print "  T:\n", T
  print "  L:\n", L
  print "  S:\n", S
  tls_as_xyz.decompose_tls(T=T, L=L, S=S)

def exercise_03():
  # getTLS_test003.mes
  T = matrix.sym(sym_mat3=[0.11000001, 0.34999996, 0.34999996,
                          -0.05, -0.05, -0.28999996])
  L = matrix.sym(sym_mat3=[0.01, 0.04, 0.09, 0,0,0])
  S = matrix.sqr([0,0,0,0,0,0,0,0,0])
  print_step("Input (getTLS_test003.mes):")
  print "  T:\n", T
  print "  L:\n", L
  print "  S:\n", S
  tls_as_xyz.decompose_tls(T=T, L=L, S=S)

def exercise_04():
  # getTLS_test004.mes
  T = matrix.sym(sym_mat3=[0.18685097, 0.45522982, 0.16791928,
                           -0.13412423, 0.03046584, -0.25211182])
  L = matrix.sym(sym_mat3=[0.01, 0.04, 0.09, 0,0,0])
  S = matrix.sqr([0,0,0,0,0,0,0,0,0])
  print_step("Input (getTLS_test004.mes):")
  print "  T:\n", T
  print "  L:\n", L
  print "  S:\n", S
  tls_as_xyz.decompose_tls(T=T, L=L, S=S)

def exercise_11():
  # getTLS_test011.mes
  T = matrix.sym(sym_mat3=[0.32500002, 0.32500002, 0.16000001,
                           -0.315, 0,0])
  L = matrix.sym(sym_mat3=[0.05, 0.05, 0.04, -0.04, 0,0])
  S = matrix.sqr([0,0,0,0,0,0,0,0,0])
  print_step("Input (getTLS_test011.mes):")
  print "  T:\n", T
  print "  L:\n", L
  print "  S:\n", S
  tls_as_xyz.decompose_tls(T=T, L=L, S=S)


if (__name__ == "__main__"):
  exercise_01()
  exercise_02()
  exercise_03()
  exercise_04()
  exercise_11()
