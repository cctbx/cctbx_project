from __future__ import absolute_import, division, print_function

from mmtbx.validation.ramalyze import find_region_max_value
from mmtbx.validation.ramalyze import RAMA_GENERAL, RAMA_GLYCINE, RAMA_CISPRO, \
    RAMA_TRANSPRO, RAMA_PREPRO, RAMA_ILE_VAL

def exercise_1():
  # ((-63.0, -43.0), 1.0)
  print(find_region_max_value(RAMA_GENERAL, -100,0))
  assert find_region_max_value(RAMA_GENERAL, -100,0) == ((-63.0, -43.0), 1.0)
  assert find_region_max_value(RAMA_GENERAL, -120, 40) == ((-63.0, -43.0), 1.0)
  assert find_region_max_value(RAMA_GENERAL, -45, -55) == ((-63.0, -43.0), 1.0)
  assert find_region_max_value(RAMA_GENERAL, 0, 0) == None

  # ((-115.0, 131.0), 0.57068)
  assert find_region_max_value(RAMA_GENERAL, -80, 100) == ((-115.0, 131.0), 0.57068)
  assert find_region_max_value(RAMA_GENERAL, -160, 179) == ((-115.0, 131.0), 0.57068)
  assert find_region_max_value(RAMA_GENERAL, -80, 60) == ((-115.0, 131.0), 0.57068)
  assert find_region_max_value(RAMA_GENERAL, -120, -179) == ((-115.0, 131.0), 0.57068)

  # ((53.0, 43.0), 0.323004)
  assert find_region_max_value(RAMA_GENERAL, 60, 40) == ((53.0, 43.0), 0.323004)
  assert find_region_max_value(RAMA_GENERAL, 75, 0) == ((53.0, 43.0), 0.323004)
  assert find_region_max_value(RAMA_GENERAL, 60, 60) == ((53.0, 43.0), 0.323004)

  # ((53.0, -127.0), 0.0246619)
  assert find_region_max_value(RAMA_GENERAL, 53, -127) == ((53.0, -127.0), 0.0246619)
  assert find_region_max_value(RAMA_GENERAL, 54, -128) == ((53.0, -127.0), 0.0246619)
  assert find_region_max_value(RAMA_GENERAL, 53, -130) == ((53.0, -127.0), 0.0246619)

  print("==================================================")
  # ((-63.0, -41.0), 1.0)
  assert find_region_max_value(RAMA_GLYCINE, -80, -20) == ((-63.0, -41.0), 1.0)
  assert find_region_max_value(RAMA_GLYCINE, -80, 60) == ((-63.0, -41.0), 1.0)
  assert find_region_max_value(RAMA_GLYCINE, -80, 65) == ((-63.0, -41.0), 1.0)
  # ((63.0, 41.0), 1.0)
  assert find_region_max_value(RAMA_GLYCINE, 80, 0) == ((63.0, 41.0), 1.0)
  assert find_region_max_value(RAMA_GLYCINE, 80, -60) == ((63.0, 41.0), 1.0)
  # ((79.0, -173.0), 0.553852)
  assert find_region_max_value(RAMA_GLYCINE, 100, -160) == ((79.0, -173.0), 0.553852)
  assert find_region_max_value(RAMA_GLYCINE, 120, 160) == ((79.0, -173.0), 0.553852)
  assert find_region_max_value(RAMA_GLYCINE, -100, 140) == ((79.0, -173.0), 0.553852)
  assert find_region_max_value(RAMA_GLYCINE, -100, -160) == ((79.0, -173.0), 0.553852)

  print("==================================================")
  # ((-89.0, 5.0), 0.701149)
  assert find_region_max_value(RAMA_CISPRO, -80, 0) == ((-89.0, 5.0), 0.701149)
  assert find_region_max_value(RAMA_CISPRO, -60, -20) == ((-89.0, 5.0), 0.701149)
  assert find_region_max_value(RAMA_CISPRO, -100, 40) == ((-89.0, 5.0), 0.701149)
  # ((-75.0, 155.0), 1.0)
  assert find_region_max_value(RAMA_CISPRO, -80, 140) == ((-75.0, 155.0), 1.0)
  assert find_region_max_value(RAMA_CISPRO, -80, -178) == ((-75.0, 155.0), 1.0)

  print("==================================================")
  # ((-57.0, -37.0), 0.99566)
  assert find_region_max_value(RAMA_TRANSPRO, -60, -20) == ((-57.0, -37.0), 0.99566)
  assert find_region_max_value(RAMA_TRANSPRO, -80, 0) == ((-57.0, -37.0), 0.99566)
  assert find_region_max_value(RAMA_TRANSPRO, -40, -40) == ((-57.0, -37.0), 0.99566)
  # ((-81.0, 65.0), 0.0896269)
  assert find_region_max_value(RAMA_TRANSPRO, -80, 60) == ((-81.0, 65.0), 0.0896269)
  # ((-59.0, 143.0), 1.0)
  assert find_region_max_value(RAMA_TRANSPRO, -60, 140) == ((-59.0, 143.0), 1.0)
  assert find_region_max_value(RAMA_TRANSPRO, -80, -179) == ((-59.0, 143.0), 1.0)

  print("==================================================")
  # ((-70.1, 149.0), 0.9619998)
  assert find_region_max_value(RAMA_PREPRO, -120, 140) == ((-67.0, 147.0), 0.992025)
  assert find_region_max_value(RAMA_PREPRO, -120, 60) == ((-67.0, 147.0), 0.992025)
  assert find_region_max_value(RAMA_PREPRO, -160, 80) == ((-67.0, 147.0), 0.992025)
  assert find_region_max_value(RAMA_PREPRO, -160, 160) == ((-67.0, 147.0), 0.992025)
  # ((-57.0, -45.0), 1.0)
  assert find_region_max_value(RAMA_PREPRO, -60, -40) == ((-57.0, -45.0), 1.0)
  assert find_region_max_value(RAMA_PREPRO, -45, -55) == ((-57.0, -45.0), 1.0)
  # ((49.0, 57.0), 0.185259)
  assert find_region_max_value(RAMA_PREPRO, 49, 57) == ((49.0, 57.0), 0.185259)
  assert find_region_max_value(RAMA_PREPRO, 60, 60) == ((49.0, 57.0), 0.185259)
  assert find_region_max_value(RAMA_PREPRO, 55, 55) == ((49.0, 57.0), 0.185259)

  print("==================================================")
  # ((-63.0, -45.0), 1.0)
  assert find_region_max_value(RAMA_ILE_VAL, -60, -40) == ((-63.0, -45.0), 1.0)
  assert find_region_max_value(RAMA_ILE_VAL, -120, -60) == ((-63.0, -45.0), 1.0)
  assert find_region_max_value(RAMA_ILE_VAL, -120, 20) == ((-63.0, -45.0), 1.0)
  assert find_region_max_value(RAMA_ILE_VAL, -80, 0) == ((-63.0, -45.0), 1.0)
  # ((-121.0, 129.0), 0.76163)
  assert find_region_max_value(RAMA_ILE_VAL, -100, 140) == ((-121.0, 129.0), 0.76163)
  assert find_region_max_value(RAMA_ILE_VAL, -160, 140) == ((-121.0, 129.0), 0.76163)
  assert find_region_max_value(RAMA_ILE_VAL, -60, 140) == ((-121.0, 129.0), 0.76163)
  assert find_region_max_value(RAMA_ILE_VAL, -130, -179) == ((-121.0, 129.0), 0.76163)



if __name__ == '__main__':
  exercise_1()
