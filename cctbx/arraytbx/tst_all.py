import os

tests = (
"tst_af_1",
"tst_af_2",
"tst_af_3",
"tst_af_4",
"tst_mat3",
"tst_vec3",
"tst_sym_mat3",
"python tst_flex.py",
"python tst_flex_utils.py",
)

for tst in tests:
  os.system(tst)
