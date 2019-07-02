///
/// Generated code. DO NOT EDIT
///
/// Generator: cctbx/sgtbx/direct_space_asu/proto/generate_cpp_asu_table.py
///
/// Dependencies: cctbx/sgtbx/direct_space_asu/reference_table.py
///

////////////////
// sgtbx/direct_space_asu/reference_table.py 5e11890e01d068e93a1273a8a8bea2a4
// sgtbx/direct_space_asu/proto/generate_cpp_asu_table.py c67cc76278e002b08b66746710aedcb9
// sgtbx/direct_space_asu/short_cuts.py bf76fe6d609a622908e94403d24b40d7
////////////////

#include "reference_table.h"

namespace cctbx { namespace sgtbx { namespace asu {

namespace {

typedef sg_vec3 vvv;

facet_collection::pointer  asu_001 ()   //  Hall:  P 1
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_002 ()   //  Hall:  -P 1
{
  return facet_collection_asu(
      x0(y0(z2) & y2(z2))
    & x2(y0(z2) & y2(z2))
    & y0
    & +y1
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_003 ()   //  Hall:  P 2y
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(x2)
    & z2(x2)
  );
}

facet_collection::pointer  asu_004 ()   //  Hall:  P 2yb
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(x0(+y2) & x2(+y2))
    & z2(x0(+y2) & x2(+y2))
  );
}

facet_collection::pointer  asu_005 ()   //  Hall:  C 2y
{
  return facet_collection_asu(
      x0(z2)
    & x2(z2)
    & y0
    & +y2
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_006 ()   //  Hall:  P -2y
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & y2
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_007 ()   //  Hall:  P -2yc
{
  return facet_collection_asu(
      x0
    & +x1
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_008 ()   //  Hall:  C -2y
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & y4(+x2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_009 ()   //  Hall:  C -2yc
{
  return facet_collection_asu(
      x0
    & +x1
    & y0(+z2)
    & y4(+z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_010 ()   //  Hall:  -P 2y
{
  return facet_collection_asu(
      x0(z2)
    & x2(z2)
    & y0
    & y2
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_011 ()   //  Hall:  -P 2yb
{
  return facet_collection_asu(
      x0
    & +x1
    & y0(z0(x2) & z2(x2))
    & y4
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_012 ()   //  Hall:  -C 2y
{
  return facet_collection_asu(
      x0(z2)
    & x2(z2)
    & y0
    & y4(x4(z2))
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_013 ()   //  Hall:  -P 2yc
{
  return facet_collection_asu(
      x0(z0(y2) & z4)
    & x2(z0(y2) & z4)
    & y0
    & +y1
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_014 ()   //  Hall:  -P 2ybc
{
  return facet_collection_asu(
      x0(y0(z2))
    & +x1
    & y0(x2(z2))
    & y4(+z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_015 ()   //  Hall:  -C 2yc
{
  return facet_collection_asu(
      x0(z4)
    & x2(z4)
    & y0
    & +y2
    & z0(y4(x4))
    & z2(-y4(x4))
  );
}

facet_collection::pointer  asu_016 ()   //  Hall:  P 2 2
{
  return facet_collection_asu(
      x0(z2)
    & x2(z2)
    & y0(z2)
    & y2(z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_017 ()   //  Hall:  P 2c 2
{
  return facet_collection_asu(
      x0(-z4 & z1*3/4)
    & x2(-z4 & z1*3/4)
    & y0(z2)
    & y2(z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_018 ()   //  Hall:  P 2 2ab
{
  return facet_collection_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_019 ()   //  Hall:  P 2ac 2ab
{
  return facet_collection_asu(
      x0
    & +x2
    & y0(-z2)
    & y2(z2)
    & z0(+y2)
    & +z1
  );
}

facet_collection::pointer  asu_020 ()   //  Hall:  C 2c 2
{
  return facet_collection_asu(
      x0(z4)
    & x2(-z4)
    & y0
    & y2(-z0)
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_021 ()   //  Hall:  C 2 2
{
  return facet_collection_asu(
      x0(z2)
    & x4(y4)
    & y0(z2)
    & y2(z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_022 ()   //  Hall:  F 2 2
{
  return facet_collection_asu(
      x0(z2)
    & x4(-z4 & z1*3/4)
    & y0(z2)
    & y4(-z4 & z1*3/4)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_023 ()   //  Hall:  I 2 2
{
  return facet_collection_asu(
      x0
    & x2(-y0)
    & y0
    & y2(-z0)
    & z0
    & z2(-x0)
  );
}

facet_collection::pointer  asu_024 ()   //  Hall:  I 2b 2c
{
  return facet_collection_asu(
      x0(y4)
    & x2(y4)
    & y0(z4)
    & y2(z4)
    & z0(x4)
    & z2(x4)
  );
}

facet_collection::pointer  asu_025 ()   //  Hall:  P 2 -2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_026 ()   //  Hall:  P 2c -2
{
  return facet_collection_asu(
      x0
    & x2
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_027 ()   //  Hall:  P 2 -2c
{
  return facet_collection_asu(
      x0(+z2)
    & x2(+z2)
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_028 ()   //  Hall:  P 2 -2a
{
  return facet_collection_asu(
      x0(y2)
    & x4
    & y0
    & +y1
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_029 ()   //  Hall:  P 2c -2ac
{
  return facet_collection_asu(
      x0(+z2)
    & x4(+z2)
    & y0
    & +y1
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_030 ()   //  Hall:  P 2 -2bc
{
  return facet_collection_asu(
      x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_031 ()   //  Hall:  P 2ac -2
{
  return facet_collection_asu(
      x0
    & x2
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_032 ()   //  Hall:  P 2 -2ab
{
  return facet_collection_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_033 ()   //  Hall:  P 2c -2n
{
  return facet_collection_asu(
      x0
    & +x2
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_034 ()   //  Hall:  P 2 -2n
{
  return facet_collection_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_035 ()   //  Hall:  C 2 -2
{
  return facet_collection_asu(
      x0
    & x4(y4)
    & y0
    & y2
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_036 ()   //  Hall:  C 2c -2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & +y2
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_037 ()   //  Hall:  C 2 -2c
{
  return facet_collection_asu(
      x0(+z2)
    & x4(y4)
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_038 ()   //  Hall:  A 2 -2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_039 ()   //  Hall:  A 2 -2b
{
  return facet_collection_asu(
      x0(+z2)
    & x2(+z2)
    & y0(+z2)
    & y4
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_040 ()   //  Hall:  A 2 -2a
{
  return facet_collection_asu(
      x0(+z2)
    & x4
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_041 ()   //  Hall:  A 2 -2ab
{
  return facet_collection_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_042 ()   //  Hall:  F 2 -2
{
  return facet_collection_asu(
      x0
    & x4(+z2)
    & y0
    & y4(+z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_043 ()   //  Hall:  F 2 -2d
{
  return facet_collection_asu(
      x0
    & x4(-y0(+z2))
    & y0
    & +y4
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_044 ()   //  Hall:  I 2 -2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_045 ()   //  Hall:  I 2 -2c
{
  return facet_collection_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_046 ()   //  Hall:  I 2 -2a
{
  return facet_collection_asu(
      x0(y2)
    & x4
    & y0
    & +y1
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_047 ()   //  Hall:  -P 2 2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & z2
  );
}

facet_collection::pointer  asu_048 ()   //  Hall:  -P 2ab 2bc
{
  return facet_collection_asu(
      x0(-y0(z2))
    & x4(-z4 & z1*3/4)
    & ~y4(-z4 & z1*3/4)
    & y4(-z4 & z1*3/4)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_049 ()   //  Hall:  -P 2 2c
{
  return facet_collection_asu(
      x0(z4)
    & x2(z4)
    & y0(z4)
    & y2(z4)
    & z0
    & z2
  );
}

facet_collection::pointer  asu_050 ()   //  Hall:  -P 2ab 2b
{
  return facet_collection_asu(
      x0(-y2)
    & x4(-y4 & y1*3/4)
    & y0
    & +y1
    & z0(-y4 & y1*3/4)
    & z2(-y4 & y1*3/4)
  );
}

facet_collection::pointer  asu_051 ()   //  Hall:  -P 2a 2a
{
  return facet_collection_asu(
      x0(z2)
    & x4
    & y0
    & y2
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_052 ()   //  Hall:  -P 2a 2bc
{
  return facet_collection_asu(
      x0
    & +x1
    & y0(-x4 & x1*3/4)
    & y4(z4)
    & z0(-x2)
    & z2(-x2)
  );
}

facet_collection::pointer  asu_053 ()   //  Hall:  -P 2ac 2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & +y1
    & z0(y2)
    & z4(x4)
  );
}

facet_collection::pointer  asu_054 ()   //  Hall:  -P 2a 2ac
{
  return facet_collection_asu(
      x0(-z4)
    & x2(z4)
    & y0(-x4)
    & y2(-x4)
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_055 ()   //  Hall:  -P 2 2ab
{
  return facet_collection_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & z2
  );
}

facet_collection::pointer  asu_056 ()   //  Hall:  -P 2ab 2ac
{
  return facet_collection_asu(
      x0(y2(-z0))
    & x4(-y4 & y1*3/4)
    & y0
    & +y1
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_057 ()   //  Hall:  -P 2c 2b
{
  return facet_collection_asu(
      x0(-y2)
    & x2(-y2)
    & y0
    & +y1
    & z0(-y4 & y1*3/4)
    & z4
  );
}

facet_collection::pointer  asu_058 ()   //  Hall:  -P 2 2n
{
  return facet_collection_asu(
      x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & z2
  );
}

facet_collection::pointer  asu_059 ()   //  Hall:  -P 2ab 2a
{
  return facet_collection_asu(
      x0(-y0(z2))
    & x4
    & ~y4
    & y4
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_060 ()   //  Hall:  -P 2n 2ab
{
  return facet_collection_asu(
      x0(z4)
    & x2(-z4)
    & y0
    & y2(-x0(-z0))
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_061 ()   //  Hall:  -P 2ac 2ab
{
  return facet_collection_asu(
      x0
    & x2(-y0(-z0))
    & y0
    & +y2
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_062 ()   //  Hall:  -P 2ac 2n
{
  return facet_collection_asu(
      x0
    & x2(-y0(-z0))
    & y0(+z2)
    & y4
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_063 ()   //  Hall:  -C 2c 2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & +y2
    & z0(y4(x4))
    & z4
  );
}

facet_collection::pointer  asu_064 ()   //  Hall:  -C 2ac 2
{
  return facet_collection_asu(
      x0
    & x4(z4)
    & y0
    & +y2
    & z0(y4)
    & z2(+y4)
  );
}

facet_collection::pointer  asu_065 ()   //  Hall:  -C 2 2
{
  return facet_collection_asu(
      x0
    & x4(y4)
    & y0
    & y2
    & z0
    & z2
  );
}

facet_collection::pointer  asu_066 ()   //  Hall:  -C 2 2c
{
  return facet_collection_asu(
      x0(z4)
    & x4(y4)
    & y0(z4)
    & y2(z4)
    & z0
    & z2
  );
}

facet_collection::pointer  asu_067 ()   //  Hall:  -C 2a 2
{
  return facet_collection_asu(
      x0
    & x2
    & y0(x4)
    & y4
    & z0(x4)
    & z2(x4)
  );
}

facet_collection::pointer  asu_068 ()   //  Hall:  -C 2a 2ac
{
  return facet_collection_asu(
      x0(z4)
    & x2(z4)
    & y0(x4)
    & y4(z4)
    & z0(+x2 & y4(x4))
    & +z2
  );
}

facet_collection::pointer  asu_069 ()   //  Hall:  -F 2 2
{
  return facet_collection_asu(
      x0
    & x4(z4)
    & y0
    & y4(z4)
    & z0
    & z2
  );
}

facet_collection::pointer  asu_070 ()   //  Hall:  -F 2uv 2vw
{
  return facet_collection_asu(
      x0(-y0(z2))
    & x8(-z8 & z1*5/8)
    & ~y8(-z1*3/8 & z1*7/8)
    & y8(-z8 & z1*5/8)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_071 ()   //  Hall:  -I 2 2
{
  return facet_collection_asu(
      x0
    & x4(y4(z4))
    & y0
    & y2
    & z0
    & z2
  );
}

facet_collection::pointer  asu_072 ()   //  Hall:  -I 2 2c
{
  return facet_collection_asu(
      x0(z4)
    & x4(y4(z4))
    & y0(z4)
    & y2(z4)
    & z0
    & z2
  );
}

facet_collection::pointer  asu_073 ()   //  Hall:  -I 2b 2c
{
  return facet_collection_asu(
      x0(y4)
    & x4(z4(y4))
    & y0(z4)
    & y2(-z4)
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_074 ()   //  Hall:  -I 2b 2
{
  return facet_collection_asu(
      x0
    & x4(-z4 & z1*3/4)
    & y0(z2)
    & y4
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_075 ()   //  Hall:  P 4
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_076 ()   //  Hall:  P 4w
{
  return facet_collection_asu(
      x0(+z4)
    & x2(+z4)
    & y0(+z1*3/4)
    & y2(+z1*3/4)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_077 ()   //  Hall:  P 4c
{
  return facet_collection_asu(
      x0(+z2)
    & x2(+z2)
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_078 ()   //  Hall:   P 4cw
{
  return facet_collection_asu(
      x0(+-z1*3/4)
    & x2(+-z1*3/4)
    & y0(+-z4)
    & y2(+-z4)
    & z1
    & +z0
  );
}

facet_collection::pointer  asu_079 ()   //  Hall:  I 4
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_080 ()   //  Hall:  I 4bw
{
  return facet_collection_asu(
      x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0
    & +z4
  );
}

facet_collection::pointer  asu_081 ()   //  Hall:  P -4
{
  return facet_collection_asu(
      x0(-y0(z2))
    & x2
    & y0
    & y2(-x2(z2))
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_082 ()   //  Hall:  I -4
{
  return facet_collection_asu(
      x0(z0(-y0))
    & x2(-y0(z4))
    & y0
    & y2(-x0(z4))
    & z0
    & z2(-y0)
  );
}

facet_collection::pointer  asu_083 ()   //  Hall:  -P 4
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z2
  );
}

facet_collection::pointer  asu_084 ()   //  Hall:  -P 4c
{
  return facet_collection_asu(
      x0(-y0(z4))
    & x2
    & y0
    & y2(-x2(z4))
    & z0
    & z2
  );
}

facet_collection::pointer  asu_085 ()   //  Hall:  -P 4a
{
  return facet_collection_asu(
      ~x4(-~y4)
    & x4(z0(-~y4) & z2(-~y4))
    & ~y4
    & y4(-x4)
    & z0(-y0(-x0))
    & z2(-y0(-x0))
  );
}

facet_collection::pointer  asu_086 ()   //  Hall:  -P 4bc
{
  return facet_collection_asu(
      ~x4(-~y4(z4))
    & x4(z0(-~y4) & z2(+-~y4))
    & ~y4
    & y4(-x4(z4))
    & z0(-y0(-x0))
    & z2(-y0(-x0))
  );
}

facet_collection::pointer  asu_087 ()   //  Hall:  -I 4
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z4(y4(x4) & x2(-y0))
  );
}

facet_collection::pointer  asu_088 ()   //  Hall:  -I 4ad
{
  return facet_collection_asu(
      x0
    & x4
    & y0(-x0(z2) | -x4(+z4))
    & y4(-x0(-z8 & z1*5/8))
    & z0
    & +z1
  );
}

facet_collection::pointer  asu_089 ()   //  Hall:  P 4 2
{
  return facet_collection_asu(
      x0(p0)
    & x2
    & y0
    & y2(-x2)
    & z0(p0)
    & z2(p0)
  );
}

facet_collection::pointer  asu_090 ()   //  Hall:  P 4ab 2ab
{
  return facet_collection_asu(
      x0
    & x2(-y0)
    & y0
    & y2(-x0)
    & z0(p0)
    & z2(p0)
  );
}

facet_collection::pointer  asu_091 ()   //  Hall:  P 4w 2c
{
  return facet_collection_asu(
      x0(z8(-y0))
    & +x1
    & y0
    & +y1
    & z0(x2)
    & z8(m1)
  );
}

facet_collection::pointer  asu_092 ()   //  Hall:  P 4abw 2nw
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z8(-y2)
  );
}

facet_collection::pointer  asu_093 ()   //  Hall:  P 4c 2
{
  return facet_collection_asu(
      x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0(y2)
    & z4(-p0 & m1)
  );
}

facet_collection::pointer  asu_094 ()   //  Hall:  P 4n 2n
{
  return facet_collection_asu(
      x0(-y0)
    & x2(z2(-y2))
    & y0(z2(-x0))
    & +y2
    & z0(p0)
    & z2(p0)
  );
}

facet_collection::pointer  asu_095 ()   //  Hall:   P 4cw 2c
{
  return facet_collection_asu(
      x1(z8(-y0))
    & +x0
    & y0
    & +y1
    & z0(-x2)
    & z8(p0)
  );
}

facet_collection::pointer  asu_096 ()   //  Hall:  P 4nw 2abw
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z8(-x2)
  );
}

facet_collection::pointer  asu_097 ()   //  Hall:  I 4 2
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0(p0)
    & z4(m2)
  );
}

facet_collection::pointer  asu_098 ()   //  Hall:  I 4bw 2bw
{
  return facet_collection_asu(
      x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0(m1 & -p0)
    & z8(-y4 & y1*3/4)
  );
}

facet_collection::pointer  asu_099 ()   //  Hall:  P 4 -2
{
  return facet_collection_asu(
      x0
    & y2
    & z0
    & +z1
    & -p0
  );
}

facet_collection::pointer  asu_100 ()   //  Hall:  P 4 -2ab
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0
    & +z1
    & m2
  );
}

facet_collection::pointer  asu_101 ()   //  Hall:  P 4c -2c
{
  return facet_collection_asu(
      x0(+z2)
    & y2(+z2)
    & z0
    & +z1
    & -p0
  );
}

facet_collection::pointer  asu_102 ()   //  Hall:  P 4n -2n
{
  return facet_collection_asu(
      x0(+z2)
    & y2(+z2)
    & z0
    & +z1
    & -p0
  );
}

facet_collection::pointer  asu_103 ()   //  Hall:  P 4 -2c
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_104 ()   //  Hall:  P 4 -2n
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_105 ()   //  Hall:  P 4c -2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_106 ()   //  Hall:  P 4c -2ab
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & +y2
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_107 ()   //  Hall:  I 4 -2
{
  return facet_collection_asu(
      x0
    & y2
    & z0
    & +z2
    & -p0
  );
}

facet_collection::pointer  asu_108 ()   //  Hall:  I 4 -2c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0
    & +z2
    & m2
  );
}

facet_collection::pointer  asu_109 ()   //  Hall:  I 4bw -2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & +z4
  );
}

facet_collection::pointer  asu_110 ()   //  Hall:  I 4bw -2c
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & +y2
    & z0
    & +z4
  );
}

facet_collection::pointer  asu_111 ()   //  Hall:  P -4 2
{
  return facet_collection_asu(
      x0(z2)
    & y2(z2)
    & z0
    & +z1
    & -p0
  );
}

facet_collection::pointer  asu_112 ()   //  Hall:  P -4 2c
{
  return facet_collection_asu(
      x0(z4 & z0(-y0))
    & x2(z4)
    & y0(z4)
    & y2(z4 & z0(-x2))
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_113 ()   //  Hall:  P -4 2ab
{
  return facet_collection_asu(
      x0(-y0(z2))
    & y0
    & z0
    & +z1
    & m2
  );
}

facet_collection::pointer  asu_114 ()   //  Hall:  P -4 2n
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2(-z0))
    & z0
    & +z2
  );
}

facet_collection::pointer  asu_115 ()   //  Hall:  P -4 -2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & y2
    & z0(p0)
    & z2(p0)
  );
}

facet_collection::pointer  asu_116 ()   //  Hall:  P -4 -2c
{
  return facet_collection_asu(
      x0(y2)
    & x2(y2 & z0(-y2))
    & y0
    & +y1
    & z0(y2 & y0(-x0))
    & z4(m1 & -p0)
  );
}

facet_collection::pointer  asu_117 ()   //  Hall:  P -4 -2ab
{
  return facet_collection_asu(
      x0(z0(-y0) & z2(-y0))
    & x2(-y0)
    & y0
    & +y2
    & z0(m2)
    & z2(m2)
  );
}

facet_collection::pointer  asu_118 ()   //  Hall:  P -4 -2n
{
  return facet_collection_asu(
      x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0(y2(-x2) & x0(-y0))
    & z4(~p2 & -m2)
  );
}

facet_collection::pointer  asu_119 ()   //  Hall:  I -4 -2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & y2
    & z0(p0)
    & z4(m2)
  );
}

facet_collection::pointer  asu_120 ()   //  Hall:  I -4 -2c
{
  return facet_collection_asu(
      x0(z0(-y0))
    & x2(-y0)
    & y0
    & +y2
    & z0(m2)
    & z4(p0)
  );
}

facet_collection::pointer  asu_121 ()   //  Hall:  I -4 2
{
  return facet_collection_asu(
      x0
    & y2(-x0(z4))
    & z0
    & z2(-x0)
    & -p0
  );
}

facet_collection::pointer  asu_122 ()   //  Hall:  I -4 2bw
{
  return facet_collection_asu(
      x0(y2 & z0(-y0))
    & x2(y2)
    & y0
    & +y1
    & z0(y2(-x2))
    & z8(-y4 & y1*3/4)
  );
}

facet_collection::pointer  asu_123 ()   //  Hall:  -P 4 2
{
  return facet_collection_asu(
      x0
    & y2
    & z0
    & z2
    & -p0
  );
}

facet_collection::pointer  asu_124 ()   //  Hall:  -P 4 2c
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z4(p0)
  );
}

facet_collection::pointer  asu_125 ()   //  Hall:  -P 4a 2b
{
  return facet_collection_asu(
      ~x4(-~y4)
    & ~y4
    & z0(p0)
    & z2(p0)
    & -m0
  );
}

facet_collection::pointer  asu_126 ()   //  Hall:  -P 4a 2bc
{
  return facet_collection_asu(
      ~x4(-~y4)
    & x4
    & ~y4
    & y4(-x4)
    & z0(-y0(-x0) & x4(-~y4))
    & z4(p0)
  );
}

facet_collection::pointer  asu_127 ()   //  Hall:  -P 4 2ab
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0
    & z2
    & m2
  );
}

facet_collection::pointer  asu_128 ()   //  Hall:  -P 4 2n
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z4(m2)
  );
}

facet_collection::pointer  asu_129 ()   //  Hall:  -P 4a 2a
{
  return facet_collection_asu(
      ~x4
    & y4
    & z0(-m0)
    & z2(-m0)
    & -p0
  );
}

facet_collection::pointer  asu_130 ()   //  Hall:  -P 4a 2ac
{
  return facet_collection_asu(
      ~x4(-~y4)
    & x4(z0(-~y4))
    & ~y4
    & y4(-x4)
    & z0(-y0(-x0))
    & z4(-m0)
  );
}

facet_collection::pointer  asu_131 ()   //  Hall:  -P 4c 2
{
  return facet_collection_asu(
      x0
    & x2
    & y0
    & y2
    & z0
    & z4(p0)
  );
}

facet_collection::pointer  asu_132 ()   //  Hall:  -P 4c 2c
{
  return facet_collection_asu(
      x0(z4)
    & y2(z4)
    & z0
    & z2
    & -p0
  );
}

facet_collection::pointer  asu_133 ()   //  Hall:  -P 4ac 2b
{
  return facet_collection_asu(
      ~x4
    & x4(-z0 | -~y4)
    & ~y4
    & +y4
    & z0(-y0(-x0))
    & z4(p0)
  );
}

facet_collection::pointer  asu_134 ()   //  Hall:  -P 4ac 2bc
{
  return facet_collection_asu(
      ~x4(z4)
    & ~y4(z4)
    & z0(p0)
    & z2(p0)
    & -m0
  );
}

facet_collection::pointer  asu_135 ()   //  Hall:  -P 4c 2ab
{
  return facet_collection_asu(
      x0(z4(-y0))
    & x2(-y0)
    & y0
    & +y2
    & z0
    & z4(m2)
  );
}

facet_collection::pointer  asu_136 ()   //  Hall:  -P 4n 2n
{
  return facet_collection_asu(
      x0(z4)
    & y2(z4(-x0))
    & z0
    & z2
    & -p0
  );
}

facet_collection::pointer  asu_137 ()   //  Hall:  -P 4ac 2a
{
  return facet_collection_asu(
      ~x4
    & x4
    & ~y4
    & y4
    & z0(-y0(-x0))
    & z4(-m0)
  );
}

facet_collection::pointer  asu_138 ()   //  Hall:  -P 4ac 2ac
{
  return facet_collection_asu(
      ~x4
    & y4(-~x4(z4))
    & z0(-m0)
    & z2(-m0 & ~x4(-y4))
    & -p0
  );
}

facet_collection::pointer  asu_139 ()   //  Hall:  -I 4 2
{
  return facet_collection_asu(
      x0
    & y2
    & z0
    & z4(m2)
    & -p0
  );
}

facet_collection::pointer  asu_140 ()   //  Hall:  -I 4 2c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0
    & z4(p0)
    & m2
  );
}

facet_collection::pointer  asu_141 ()   //  Hall:  -I 4bd 2
{
  return facet_collection_asu(
      x0
    & x2
    & ~y4
    & y4
    & z0(-y0)
    & z8(-p4)
  );
}

facet_collection::pointer  asu_142 ()   //  Hall:  -I 4bd 2c
{
  return facet_collection_asu(
      x0(z8(-~y4) & z0(-y0))
    & x2(-~y4)
    & ~y4
    & +y4
    & z0(x4)
    & z8(m1/4)
  );
}

facet_collection::pointer  asu_143 ()   //  Hall:  P 3
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0
    & +z1
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_144 ()   //  Hall:  P 31
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0
    & +z3
  );
}

facet_collection::pointer  asu_145 ()   //  Hall:   P 32
{
  return facet_collection_asu(
      y0
    & +y1
    & x0
    & +x1
    & z0
    & +z3
  );
}

facet_collection::pointer  asu_146 ()   //  Hall:  R 3
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0
    & +z3
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_147 ()   //  Hall:  -P 3
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z2(p0(-y0))
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_148 ()   //  Hall:  -R 3
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z6(-h0(x3) | -k0(-y0 | -m1))
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_149 ()   //  Hall:  P 3 2
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(-h0 | -k0)
    & z2(-h0 | -k0)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_150 ()   //  Hall:  P 3 2"
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(p0)
    & z2(p0)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_151 ()   //  Hall:  P 31 2 (x,y,z+1/3)
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(-h0 | -h1)
    & z6(-k0 | -k1)
  );
}

facet_collection::pointer  asu_152 ()   //  Hall:  P 31 2"
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(-p0)
    & z6(-p0)
  );
}

facet_collection::pointer  asu_153 ()   //  Hall:  P 32 2 (x,y,z+1/6)
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(-h0 | -h1)
    & z6(x0(-y0) & m1)
  );
}

facet_collection::pointer  asu_154 ()   //  Hall:   P 32 2"
{
  return facet_collection_asu(
      y0
    & +y1
    & x0
    & +x1
    & z0(p0)
    & z6(p0)
  );
}

facet_collection::pointer  asu_155 ()   //  Hall:  R 3 2"
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(p0)
    & z6(x3 & ~p3)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_156 ()   //  Hall:  P 3 -2"
{
  return facet_collection_asu(
      z0
    & +z1
    & h0
    & m1
    & k0
  );
}

facet_collection::pointer  asu_157 ()   //  Hall:  P 3 -2
{
  return facet_collection_asu(
      y0
    & z0
    & +z1
    & k1
    & m1(y3)
    & p0
  );
}

facet_collection::pointer  asu_158 ()   //  Hall:  P 3 -2"c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0
    & +z2
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_159 ()   //  Hall:  P 3 -2c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0
    & +z2
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_160 ()   //  Hall:  R 3 -2"
{
  return facet_collection_asu(
      z0
    & +z3
    & h0
    & m1
    & k0
  );
}

facet_collection::pointer  asu_161 ()   //  Hall:  R 3 -2"c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0
    & +z6
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_162 ()   //  Hall:  -P 3 2
{
  return facet_collection_asu(
      y0
    & z0(-h0)
    & z2(-h0)
    & k1
    & m1(y3)
    & p0
  );
}

facet_collection::pointer  asu_163 ()   //  Hall:  -P 3 2c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z4(-h0 | -k0)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_164 ()   //  Hall:  -P 3 2"
{
  return facet_collection_asu(
      y0(z2)
    & z0
    & +z1
    & k1
    & -h0
  );
}

facet_collection::pointer  asu_165 ()   //  Hall:  -P 3 2"c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z4(p0)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_166 ()   //  Hall:  -R 3 2"
{
  return facet_collection_asu(
      z0(p0)
    & z6(x3)
    & h0
    & m1
    & k0
  );
}

facet_collection::pointer  asu_167 ()   //  Hall:  -R 3 2"c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z12(y3 & p3)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_168 ()   //  Hall:  P 6
{
  return facet_collection_asu(
      y0
    & z0
    & +z1
    & k1
    & m1(y3)
    & p0(-y0)
  );
}

facet_collection::pointer  asu_169 ()   //  Hall:  P 61
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0
    & +z6
  );
}

facet_collection::pointer  asu_170 ()   //  Hall:   P 65
{
  return facet_collection_asu(
      y0
    & +y1
    & x0
    & +x1
    & z0
    & +z6
  );
}

facet_collection::pointer  asu_171 ()   //  Hall:  P 62
{
  return facet_collection_asu(
      x1(y2)
    & y0(x2)
    & z0
    & +z3
    & p0(y2)
  );
}

facet_collection::pointer  asu_172 ()   //  Hall:   P 64
{
  return facet_collection_asu(
      y0(-x2)
    & x1(-y2)
    & z0
    & +z3
    & p0(-x2)
  );
}

facet_collection::pointer  asu_173 ()   //  Hall:  P 6c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0
    & +z2
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_174 ()   //  Hall:  P -6
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0
    & z2
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_175 ()   //  Hall:  -P 6
{
  return facet_collection_asu(
      y0
    & z0
    & z2
    & k1
    & m1(y3)
    & p0(-y0)
  );
}

facet_collection::pointer  asu_176 ()   //  Hall:  -P 6c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(p0(-y0))
    & z4
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_177 ()   //  Hall:  P 6 2
{
  return facet_collection_asu(
      y0
    & z0(-h0)
    & z2(-h0)
    & k1
    & m1(y3)
    & p0(-y0)
  );
}

facet_collection::pointer  asu_178 ()   //  Hall:  P 61 2 (x,y,z+5/12)
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z12(-h0 | -h1)
  );
}

facet_collection::pointer  asu_179 ()   //  Hall:  P 65 2 (x,y,z+1/12)
{
  return facet_collection_asu(
      x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z12(m1 & x0(-y0))
  );
}

facet_collection::pointer  asu_180 ()   //  Hall:  P 62 2 (x,y,z+1/3)
{
  return facet_collection_asu(
      x1(y2)
    & y0(x2)
    & z0(k1)
    & z6(-h0)
    & p0(y2)
  );
}

facet_collection::pointer  asu_181 ()   //  Hall:   P 64 2 (x,y,z+1/6)
{
  return facet_collection_asu(
      y0(p2)
    & p0(-y2)
    & z6(-m1)
    & z0(k1)
    & x1(p2)
  );
}

facet_collection::pointer  asu_182 ()   //  Hall:  P 6c 2c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(p0)
    & z4(-h0 | -k0)
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_183 ()   //  Hall:  P 6 -2
{
  return facet_collection_asu(
      y0
    & z0
    & +z1
    & k1
    & -h0
  );
}

facet_collection::pointer  asu_184 ()   //  Hall:  P 6 -2c
{
  return facet_collection_asu(
      y0
    & z0
    & +z2
    & k1
    & m1(y3)
    & p0(-y0)
  );
}

facet_collection::pointer  asu_185 ()   //  Hall:  P 6c -2
{
  return facet_collection_asu(
      y0
    & z0
    & +z2
    & k1
    & m1(y3)
    & p0
  );
}

facet_collection::pointer  asu_186 ()   //  Hall:  P 6c -2c
{
  return facet_collection_asu(
      y0(+z2)
    & z0
    & +z1
    & k1
    & -h0
  );
}

facet_collection::pointer  asu_187 ()   //  Hall:  P -6 2
{
  return facet_collection_asu(
      z0
    & z2
    & h0
    & m1
    & k0
  );
}

facet_collection::pointer  asu_188 ()   //  Hall:  P -6c 2
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(-h0 | -k0)
    & z4
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_189 ()   //  Hall:  P -6 -2
{
  return facet_collection_asu(
      y0
    & z0
    & z2
    & k1
    & m1(y3)
    & p0
  );
}

facet_collection::pointer  asu_190 ()   //  Hall:  P -6c -2c
{
  return facet_collection_asu(
      x0(-y0)
    & y0
    & z0(p0)
    & z4
    & k1
    & m1(-h1 | -k1)
    & h1
  );
}

facet_collection::pointer  asu_191 ()   //  Hall:  -P 6 2
{
  return facet_collection_asu(
      y0
    & z0
    & z2
    & k1
    & -h0
  );
}

facet_collection::pointer  asu_192 ()   //  Hall:  -P 6 2c
{
  return facet_collection_asu(
      y0
    & z0
    & z4(-h0)
    & k1
    & m1(y3)
    & p0(-y0)
  );
}

facet_collection::pointer  asu_193 ()   //  Hall:  -P 6c 2
{
  return facet_collection_asu(
      y0
    & z0(-h0)
    & z4
    & k1
    & m1(y3)
    & p0
  );
}

facet_collection::pointer  asu_194 ()   //  Hall:  -P 6c 2c
{
  return facet_collection_asu(
      z0(p0)
    & z4
    & h0
    & m1
    & k0
  );
}

facet_collection::pointer  asu_195 ()   //  Hall:  P 2 2 3
{
  return facet_collection_asu(
      z0(y2 & x2)
    & m1(-y2)
    & cut(vvv(0,1,-1),0)(cut(vvv(-1,0,1),0))
    & cut(vvv(1,0,-1),0)
  );
}

facet_collection::pointer  asu_196 ()   //  Hall:  F 2 2 3
{
  return facet_collection_asu(
      p0(m2)
    & cut(vvv(-1,0,-1),r1/2)(cut(vvv(0,-1,1),0))
    & cut(vvv(-1,0,1),r1/2)(cut(vvv(0,-1,-1),0))
    & cut(vvv(0,1,1),0)
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_197 ()   //  Hall:  I 2 2 3
{
  return facet_collection_asu(
      z0(x2)
    & p0(cut(vvv(0,-1,1),0))
    & +m1
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_198 ()   //  Hall:  P 2ac 2ab 3
{
  return facet_collection_asu(
      x0(-y0)
    & x2
    & y2(+z0 & x2(+z2))
    & cut(vvv(-1,0,1),r1/2)(m2)
    & cut(vvv(1,0,-1),0)(p0)
    & cut(vvv(0,1,1),0)
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_199 ()   //  Hall:  I 2b 2c 3
{
  return facet_collection_asu(
      x2(-y4)
    & y2(-z4)
    & z0(x4)
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0)(+x2))
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_200 ()   //  Hall:  -P 2 2 3
{
  return facet_collection_asu(
      x2
    & y2
    & z0
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_201 ()   //  Hall:  -P 2ab 2bc 3
{
  return facet_collection_asu(
      ~z4(x4)
    & p0(cut(vvv(0,-1,1),0)(-x0))
    & m2(cut(vvv(0,-1,1),0)(x2))
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_202 ()   //  Hall:  -F 2 2 3
{
  return facet_collection_asu(
      z0
    & p0(x4)
    & cut(vvv(-1,0,-1),r1/2)(cut(vvv(0,-1,1),0))
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_203 ()   //  Hall:  -F 2uv 2vw 3
{
  return facet_collection_asu(
      p0(cut(vvv(0,-1,1),0)(-x0))
    & (m1/4)(cut(vvv(0,-1,1),0) | cut(vvv(0,-1,-1),-r1/4)(~z4))
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,1,1),r1/4)
  );
}

facet_collection::pointer  asu_204 ()   //  Hall:  -I 2 2 3
{
  return facet_collection_asu(
      x2
    & z0
    & p0(cut(vvv(0,-1,1),0)(x4))
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_205 ()   //  Hall:  -P 2ac 2ab 3
{
  return facet_collection_asu(
      x2(-z0(cut(vvv(0,-1,1),0)))
    & y2(cut(vvv(0,-1,1),0))
    & z0
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_206 ()   //  Hall:  -I 2b 2c 3
{
  return facet_collection_asu(
      z0(x4)
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(-1,0,-1),r1/2)
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,-1,-1),r1/2)(cut(vvv(-1,0,1),0))
  );
}

facet_collection::pointer  asu_207 ()   //  Hall:  P 4 2 3
{
  return facet_collection_asu(
      z0(x2)
    & p0
    & m1(-p0)
    & cut(vvv(0,1,-1),0)(x2)
  );
}

facet_collection::pointer  asu_208 ()   //  Hall:  P 4n 2 3
{
  return facet_collection_asu(
      cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(1,0,1),0)(cut(vvv(0,-1,-1),0))
    & cut(vvv(-1,0,1),r1/2)(y4)
    & cut(vvv(-1,0,-1),r1/2)(y4)
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,1,1),0)
    & cut(vvv(0,-1,1),r1/2)(-x4)
    & cut(vvv(0,-1,-1),r1/2)(-x4)
  );
}

facet_collection::pointer  asu_209 ()   //  Hall:  F 4 2 3
{
  return facet_collection_asu(
      p0(z0)
    & m2(z0)
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,1,1),0)(cut(vvv(0,-1,1),0))
  );
}

facet_collection::pointer  asu_210 ()   //  Hall:  F 4d 2 3
{
  return facet_collection_asu(
      y8(cut(vvv(1,0,1),-r1/4))
    & z8(m1/4)
    & p0(cut(vvv(-1,0,1),0))
    & m2(cut(vvv(1,0,1),-r1/2))
    & cut(vvv(0,1,1),0)(z0)
    & cut(vvv(1,0,-1),0)
    & cut(vvv(-1,0,-1),r1/2)
  );
}

facet_collection::pointer  asu_211 ()   //  Hall:  I 4 2 3
{
  return facet_collection_asu(
      z0(p0)
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(-1,0,-1),r1/2)(y4)
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,-1,-1),r1/2)(-x4)
  );
}

facet_collection::pointer  asu_212 ()   //  Hall:  P 4acd 2ab 3
{
  return facet_collection_asu(
      cut(vvv(-1,0,1),r1/2)
    & cut(vvv(0,1,1),0)(cut(vvv(1,0,-1),-r1/2))
    & cut(vvv(0,-1,-1),r1/2)(cut(vvv(-2,1,1),0))
    & cut(vvv(2,-1,-1),0)(x8)
    & cut(vvv(-1,2,-1),0)(y8)
    & cut(vvv(-2,1,-1),r1/2)(-x1*3/8)
  );
}

facet_collection::pointer  asu_213 ()   //  Hall:   P 4bd 2ab 3
{
  return facet_collection_asu(
      cut(vvv(0,1,-1),0)
    & -p0(cut(vvv(0,-1,1),0))
    & ~p2(cut(vvv(-1,1,-2),0))
    & cut(vvv(1,-1,2),0)(z8)
    & cut(vvv(-2,-1,-1),r1*3/2)(-x1*3/8)
    & cut(vvv(-1,-1,-2),r1*3/2)(-z1*3/8)
  );
}

facet_collection::pointer  asu_214 ()   //  Hall:  I 4bd 2c 3
{
  return facet_collection_asu(
      x8(cut(vvv(0,-1,-1),r1/4))
    & y8(cut(vvv(-1,0,-1),r1/4))
    & ~y8(cut(vvv(-1,0,1),-r1/4))
    & cut(vvv(-1,0,1),0)(cut(vvv(0,1,-1),0))
    & cut(vvv(0,-1,1),0)
    & cut(vvv(0,1,-1),r1/4)(-y0)
    & cut(vvv(1,-1,1),r1/8)(~p4)
  );
}

facet_collection::pointer  asu_215 ()   //  Hall:  P -4 2 3
{
  return facet_collection_asu(
      z0(x2)
    & p0
    & m1
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_216 ()   //  Hall:  F -4 2 3
{
  return facet_collection_asu(
      p0
    & m2
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,1,1),0)
  );
}

facet_collection::pointer  asu_217 ()   //  Hall:  I -4 2 3
{
  return facet_collection_asu(
      x2(-z0(y4))
    & z0
    & p0
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_218 ()   //  Hall:  P -4n 2 3
{
  return facet_collection_asu(
      x2(-z0(y4))
    & y2(-z0(x4))
    & z0
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_219 ()   //  Hall:  F -4a 2 3
{
  return facet_collection_asu(
      p0
    & m2(-p0(z0))
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,1,1),0)(-p0 | cut(vvv(0,-1,1),0)(x4))
  );
}

facet_collection::pointer  asu_220 ()   //  Hall:  I -4bd 2c 3
{
  return facet_collection_asu(
      -x4(-z0(-y1*3/8))
    & x2
    & -y4(-x2(-z8))
    & y2(-z4)
    & z0
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_221 ()   //  Hall:  -P 4 2 3
{
  return facet_collection_asu(
      x2
    & z0
    & p0
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_222 ()   //  Hall:  -P 4a 2bc 3
{
  return facet_collection_asu(
      (x1*3/4)(z4(y2) | cut(vvv(0,-1,1),0))
    & -z4
    & p0(cut(vvv(0,-1,1),0)(z2))
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_223 ()   //  Hall:  -P 4n 2 3
{
  return facet_collection_asu(
      z0
    & cut(vvv(1,0,-1),0)(cut(vvv(0,-1,1),0))
    & cut(vvv(-1,0,-1),r1/2)(y4)
    & cut(vvv(0,1,-1),0)
    & cut(vvv(0,-1,-1),r1/2)(x4)
  );
}

facet_collection::pointer  asu_224 ()   //  Hall:  -P 4bc 2bc 3
{
  return facet_collection_asu(
      p0
    & cut(vvv(-1,0,-1),r1)(y2)
    & cut(vvv(-1,0,1),r1/2)(y2)
    & cut(vvv(0,1,1),-r1/2)
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_225 ()   //  Hall:  -F 4 2 3
{
  return facet_collection_asu(
      z0
    & p0
    & m2
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_226 ()   //  Hall:  -F 4a 2 3
{
  return facet_collection_asu(
      z0
    & p0
    & m2(-p0)
    & cut(vvv(0,1,-1),0)(x4)
  );
}

facet_collection::pointer  asu_227 ()   //  Hall:  -F 4vw 2vw 3
{
  return facet_collection_asu(
      -y0(cut(vvv(1,0,1),0))
    & p0
    & m1/4
    & cut(vvv(0,1,1),r1/4)
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_228 ()   //  Hall:  -F 4ud 2vw 3
{
  return facet_collection_asu(
      -y0(cut(vvv(-1,0,1),r1/4))
    & p0(cut(vvv(0,-1,1),0))
    & m1/4
    & cut(vvv(0,1,1),r1/4)(cut(vvv(0,-1,1),0)(x8))
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_229 ()   //  Hall:  -I 4 2 3
{
  return facet_collection_asu(
      z0
    & p0
    & cut(vvv(-1,0,-1),r1/2)(y4)
    & cut(vvv(0,1,-1),0)
  );
}

facet_collection::pointer  asu_230 ()   //  Hall:  -I 4bd 2c 3
{
  return facet_collection_asu(
      x8(cut(vvv(0,1,-1),r1/4) & cut(vvv(0,-1,-1),r1/4))
    & ~x8(y0(-z4))
    & y8(cut(vvv(1,0,1),-r1/4))
    & ~y8(cut(vvv(1,0,-1),r1/4))
    & z4(y0)
    & cut(vvv(-1,0,1),0)
    & cut(vvv(1,0,1),0)(-z0)
    & cut(vvv(0,-1,1),0)(cut(vvv(1,0,-1),0))
    & cut(vvv(0,1,1),0)
  );
}

} // end of unnamed namespace

asu_func asu_table[230] = {
  asu_001, asu_002, asu_003, asu_004, asu_005, asu_006, asu_007, asu_008,
  asu_009, asu_010, asu_011, asu_012, asu_013, asu_014, asu_015, asu_016,
  asu_017, asu_018, asu_019, asu_020, asu_021, asu_022, asu_023, asu_024,
  asu_025, asu_026, asu_027, asu_028, asu_029, asu_030, asu_031, asu_032,
  asu_033, asu_034, asu_035, asu_036, asu_037, asu_038, asu_039, asu_040,
  asu_041, asu_042, asu_043, asu_044, asu_045, asu_046, asu_047, asu_048,
  asu_049, asu_050, asu_051, asu_052, asu_053, asu_054, asu_055, asu_056,
  asu_057, asu_058, asu_059, asu_060, asu_061, asu_062, asu_063, asu_064,
  asu_065, asu_066, asu_067, asu_068, asu_069, asu_070, asu_071, asu_072,
  asu_073, asu_074, asu_075, asu_076, asu_077, asu_078, asu_079, asu_080,
  asu_081, asu_082, asu_083, asu_084, asu_085, asu_086, asu_087, asu_088,
  asu_089, asu_090, asu_091, asu_092, asu_093, asu_094, asu_095, asu_096,
  asu_097, asu_098, asu_099, asu_100, asu_101, asu_102, asu_103, asu_104,
  asu_105, asu_106, asu_107, asu_108, asu_109, asu_110, asu_111, asu_112,
  asu_113, asu_114, asu_115, asu_116, asu_117, asu_118, asu_119, asu_120,
  asu_121, asu_122, asu_123, asu_124, asu_125, asu_126, asu_127, asu_128,
  asu_129, asu_130, asu_131, asu_132, asu_133, asu_134, asu_135, asu_136,
  asu_137, asu_138, asu_139, asu_140, asu_141, asu_142, asu_143, asu_144,
  asu_145, asu_146, asu_147, asu_148, asu_149, asu_150, asu_151, asu_152,
  asu_153, asu_154, asu_155, asu_156, asu_157, asu_158, asu_159, asu_160,
  asu_161, asu_162, asu_163, asu_164, asu_165, asu_166, asu_167, asu_168,
  asu_169, asu_170, asu_171, asu_172, asu_173, asu_174, asu_175, asu_176,
  asu_177, asu_178, asu_179, asu_180, asu_181, asu_182, asu_183, asu_184,
  asu_185, asu_186, asu_187, asu_188, asu_189, asu_190, asu_191, asu_192,
  asu_193, asu_194, asu_195, asu_196, asu_197, asu_198, asu_199, asu_200,
  asu_201, asu_202, asu_203, asu_204, asu_205, asu_206, asu_207, asu_208,
  asu_209, asu_210, asu_211, asu_212, asu_213, asu_214, asu_215, asu_216,
  asu_217, asu_218, asu_219, asu_220, asu_221, asu_222, asu_223, asu_224,
  asu_225, asu_226, asu_227, asu_228, asu_229, asu_230
};

}}}
