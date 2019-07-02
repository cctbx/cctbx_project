from __future__ import absolute_import, division, print_function
from cctbx.sgtbx.direct_space_asu import direct_space_asu
from cctbx.sgtbx.direct_space_asu.short_cuts import *
from six.moves import range

change_of_basis_tab = {
   78: ( 76, "a,b,-c+1"),
   95: ( 91, "-a+1,b,c"),
  145: (144, "b,a,c"),
  154: (152, "b,a,c"),
  170: (169, "b,a,c"),
  172: (171, "-b+1,-a+1,c"),
  181: (180, "b+1,-a-b,-c+1/6"),
  213: (212, "-b,c+1/2,a-1/2")
}

def apply_change_of_basis(target_sg_no):
  source_sg_no, cb_expr = change_of_basis_tab[target_sg_no]
  return get_asu(source_sg_no).change_basis(cb_expr)

def asu_001(): # P 1
  return (direct_space_asu('P 1')
    & x0
    & +x1
    & y0
    & +y1
    & z0
    & +z1
  )

def asu_002(): # P -1
  return (direct_space_asu('-P 1')
    & x0(y0(z2) & y2(z2))
    & x2(y0(z2) & y2(z2))
    & y0
    & +y1
    & z0
    & +z1
  )

def asu_003(): # P 1 2 1
  return (direct_space_asu('P 2y')
    & x0
    & +x1
    & y0
    & +y1
    & z0(x2)
    & z2(x2)
  )

def asu_004(): # P 1 21 1
  return (direct_space_asu('P 2yb')
    & x0
    & +x1
    & y0
    & +y1
    & z0(x0(+y2) & x2(+y2))
    & z2(x0(+y2) & x2(+y2))
  )

def asu_005(): # C 1 2 1
  return (direct_space_asu('C 2y')
    & x0(z2)
    & x2(z2)
    & y0
    & +y2
    & z0
    & +z1
  )

def asu_006(): # P 1 m 1
  return (direct_space_asu('P -2y')
    & x0
    & +x1
    & y0
    & y2
    & z0
    & +z1
  )

def asu_007(): # P 1 c 1
  return (direct_space_asu('P -2yc')
    & x0
    & +x1
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  )

def asu_008(): # C 1 m 1
  return (direct_space_asu('C -2y')
    & x0
    & +x1
    & y0
    & y4(+x2)
    & z0
    & +z1
  )

def asu_009(): # C 1 c 1
  return (direct_space_asu('C -2yc')
    & x0
    & +x1
    & y0(+z2)
    & y4(+z2)
    & z0
    & +z1
  )

def asu_010(): # P 1 2/m 1
  return (direct_space_asu('-P 2y')
    & x0(z2)
    & x2(z2)
    & y0
    & y2
    & z0
    & +z1
  )

def asu_011(): # P 1 21/m 1
  return (direct_space_asu('-P 2yb')
    & x0
    & +x1
    & y0(z0(x2) & z2(x2))
    & y4
    & z0
    & +z1
  )

def asu_012(): # C 1 2/m 1
  return (direct_space_asu('-C 2y')
    & x0(z2)
    & x2(z2)
    & y0
    & y4(x4(z2))
    & z0
    & +z1
  )

def asu_013(): # P 1 2/c 1
  return (direct_space_asu('-P 2yc')
    & x0(z0(y2) & z4)
    & x2(z0(y2) & z4)
    & y0
    & +y1
    & z0
    & +z2
  )

def asu_014(): # P 1 21/c 1
  return (direct_space_asu('-P 2ybc')
    & x0(y0(z2))
    & +x1
    & y0(x2(z2))
    & y4(+z2)
    & z0
    & +z1
  )

def asu_015(): # C 1 2/c 1
  return (direct_space_asu('-C 2yc')
    & x0(z4)
    & x2(z4)
    & y0
    & +y2
    & z0(y4(x4))
    & z2(-y4(x4))
  )

def asu_016(): # P 2 2 2
  return (direct_space_asu('P 2 2')
    & x0(z2)
    & x2(z2)
    & y0(z2)
    & y2(z2)
    & z0
    & +z1
  )

def asu_017(): # P 2 2 21
  return (direct_space_asu('P 2c 2')
    & x0(-z4 & z1*3/4)
    & x2(-z4 & z1*3/4)
    & y0(z2)
    & y2(z2)
    & z0
    & +z1
  )

def asu_018(): # P 21 21 2
  return (direct_space_asu('P 2 2ab')
    & x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z1
  )

def asu_019(): # P 21 21 21
  return (direct_space_asu('P 2ac 2ab')
    & x0
    & +x2
    & y0(-z2)
    & y2(z2)
    & z0(+y2)
    & +z1
  )

def asu_020(): # C 2 2 21
  return (direct_space_asu('C 2c 2')
    & x0(z4)
    & x2(-z4)
    & y0
    & y2(-z0)
    & z0
    & +z2
  )

def asu_021(): # C 2 2 2
  return (direct_space_asu('C 2 2')
    & x0(z2)
    & x4(y4)
    & y0(z2)
    & y2(z2)
    & z0
    & +z1
  )

def asu_022(): # F 2 2 2
  return (direct_space_asu('F 2 2')
    & x0(z2)
    & x4(-z4 & z1*3/4)
    & y0(z2)
    & y4(-z4 & z1*3/4)
    & z0
    & +z1
  )

def asu_023(): # I 2 2 2
  return (direct_space_asu('I 2 2')
    & x0
    & x2(-y0)
    & y0
    & y2(-z0)
    & z0
    & z2(-x0)
  )

def asu_024(): # I 21 21 21
  return (direct_space_asu('I 2b 2c')
    & x0(y4)
    & x2(y4)
    & y0(z4)
    & y2(z4)
    & z0(x4)
    & z2(x4)
  )

def asu_025(): # P m m 2
  return (direct_space_asu('P 2 -2')
    & x0
    & x2
    & y0
    & y2
    & z0
    & +z1
  )

def asu_026(): # P m c 21
  return (direct_space_asu('P 2c -2')
    & x0
    & x2
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  )

def asu_027(): # P c c 2
  return (direct_space_asu('P 2 -2c')
    & x0(+z2)
    & x2(+z2)
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  )

def asu_028(): # P m a 2
  return (direct_space_asu('P 2 -2a')
    & x0(y2)
    & x4
    & y0
    & +y1
    & z0
    & +z1
  )

def asu_029(): # P c a 21
  return (direct_space_asu('P 2c -2ac')
    & x0(+z2)
    & x4(+z2)
    & y0
    & +y1
    & z0
    & +z1
  )

def asu_030(): # P n c 2
  return (direct_space_asu('P 2 -2bc')
    & x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0
    & +z2
  )

def asu_031(): # P m n 21
  return (direct_space_asu('P 2ac -2')
    & x0
    & x2
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  )

def asu_032(): # P b a 2
  return (direct_space_asu('P 2 -2ab')
    & x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z1
  )

def asu_033(): # P n a 21
  return (direct_space_asu('P 2c -2n')
    & x0
    & +x2
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  )

def asu_034(): # P n n 2
  return (direct_space_asu('P 2 -2n')
    & x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z1
  )

def asu_035(): # C m m 2
  return (direct_space_asu('C 2 -2')
    & x0
    & x4(y4)
    & y0
    & y2
    & z0
    & +z1
  )

def asu_036(): # C m c 21
  return (direct_space_asu('C 2c -2')
    & x0
    & x2
    & y0
    & +y2
    & z0
    & +z2
  )

def asu_037(): # C c c 2
  return (direct_space_asu('C 2 -2c')
    & x0(+z2)
    & x4(y4)
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  )

def asu_038(): # A m m 2
  return (direct_space_asu('A 2 -2')
    & x0
    & x2
    & y0
    & y2
    & z0
    & +z2
  )

def asu_039(): # A b m 2
  return (direct_space_asu('A 2 -2b')
    & x0(+z2)
    & x2(+z2)
    & y0(+z2)
    & y4
    & z0
    & +z1
  )

def asu_040(): # A m a 2
  return (direct_space_asu('A 2 -2a')
    & x0(+z2)
    & x4
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  )

def asu_041(): # A b a 2
  return (direct_space_asu('A 2 -2ab')
    & x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z2
  )

def asu_042(): # F m m 2
  return (direct_space_asu('F 2 -2')
    & x0
    & x4(+z2)
    & y0
    & y4(+z2)
    & z0
    & +z1
  )

def asu_043(): # F d d 2
  return (direct_space_asu('F 2 -2d')
    & x0
    & x4(-y0(+z2))
    & y0
    & +y4
    & z0
    & +z1
  )

def asu_044(): # I m m 2
  return (direct_space_asu('I 2 -2')
    & x0
    & x2
    & y0
    & y2
    & z0
    & +z2
  )

def asu_045(): # I b a 2
  return (direct_space_asu('I 2 -2c')
    & x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z2
  )

def asu_046(): # I m a 2
  return (direct_space_asu('I 2 -2a')
    & x0(y2)
    & x4
    & y0
    & +y1
    & z0
    & +z2
  )

def asu_047(): # P m m m
  return (direct_space_asu('-P 2 2')
    & x0
    & x2
    & y0
    & y2
    & z0
    & z2
  )

def asu_048(): # P n n n :2
  return (direct_space_asu('-P 2ab 2bc')
    & x0(-y0(z2))
    & x4(-z4 & z1*3/4)
    & ~y4(-z4 & z1*3/4)
    & y4(-z4 & z1*3/4)
    & z0
    & +z1
  )

def asu_049(): # P c c m
  return (direct_space_asu('-P 2 2c')
    & x0(z4)
    & x2(z4)
    & y0(z4)
    & y2(z4)
    & z0
    & z2
  )

def asu_050(): # P b a n :2
  return (direct_space_asu('-P 2ab 2b')
    & x0(-y2)
    & x4(-y4 & y1*3/4)
    & y0
    & +y1
    & z0(-y4 & y1*3/4)
    & z2(-y4 & y1*3/4)
  )

def asu_051(): # P m m a
  return (direct_space_asu('-P 2a 2a')
    & x0(z2)
    & x4
    & y0
    & y2
    & z0
    & +z1
  )

def asu_052(): # P n n a
  return (direct_space_asu('-P 2a 2bc')
    & x0
    & +x1
    & y0(-x4 & x34)
    & y4(z4)
    & z0(-x2)
    & z2(-x2)
  )

def asu_053(): # P m n a
  return (direct_space_asu('-P 2ac 2')
    & x0
    & x2
    & y0
    & +y1
    & z0(y2)
    & z4(x4)
  )

def asu_054(): # P c c a
  return (direct_space_asu('-P 2a 2ac')
    & x0(-z4)
    & x2(z4)
    & y0(-x4)
    & y2(-x4)
    & z0
    & +z2
  )

def asu_055(): # P b a m
  return (direct_space_asu('-P 2 2ab')
    & x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & z2
  )

def asu_056(): # P c c n
  return (direct_space_asu('-P 2ab 2ac')
    & x0(y2(-z0))
    & x4(-y4 & y1*3/4)
    & y0
    & +y1
    & z0
    & +z2
  )

def asu_057(): # P b c m
  return (direct_space_asu('-P 2c 2b')
    & x0(-y2)
    & x2(-y2)
    & y0
    & +y1
    & z0(-y4 & y1*3/4)
    & z4
  )

def asu_058(): # P n n m
  return (direct_space_asu('-P 2 2n')
    & x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & z2
  )

def asu_059(): # P m m n :2
  return (direct_space_asu('-P 2ab 2a')
    & x0(-y0(z2))
    & x4
    & ~y4
    & y4
    & z0
    & +z1
  )

def asu_060(): # P b c n
  return (direct_space_asu('-P 2n 2ab')
    & x0(z4)
    & x2(-z4)
    & y0
    & y2(-x0(-z0))
    & z0
    & +z2
  )

def asu_061(): # P b c a
  return (direct_space_asu('-P 2ac 2ab')
    & x0
    & x2(-y0(-z0))
    & y0
    & +y2
    & z0
    & +z2
  )

def asu_062(): # P n m a
  return (direct_space_asu('-P 2ac 2n')
    & x0
    & x2(-y0(-z0))
    & y0(+z2)
    & y4
    & z0
    & +z1
  )

def asu_063(): # C m c m
  return (direct_space_asu('-C 2c 2')
    & x0
    & x2
    & y0
    & +y2
    & z0(y4(x4))
    & z4
  )

def asu_064(): # C m c a
  return (direct_space_asu('-C 2ac 2')
    & x0
    & x4(z4)
    & y0
    & +y2
    & z0(y4)
    & z2(+y4)
  )

def asu_065(): # C m m m
  return (direct_space_asu('-C 2 2')
    & x0
    & x4(y4)
    & y0
    & y2
    & z0
    & z2
  )

def asu_066(): # C c c m
  return (direct_space_asu('-C 2 2c')
    & x0(z4)
    & x4(y4)
    & y0(z4)
    & y2(z4)
    & z0
    & z2
  )

def asu_067(): # C m m a
  return (direct_space_asu('-C 2a 2')
    & x0
    & x2
    & y0(x4)
    & y4
    & z0(x4)
    & z2(x4)
  )

def asu_068(): # C c c a :2
  return (direct_space_asu('-C 2a 2ac')
    & x0(z4)
    & x2(z4)
    & y0(x4)
    & y4(z4)
    & z0(+x2 & y4(x4))
    & +z2
  )

def asu_069(): # F m m m
  return (direct_space_asu('-F 2 2')
    & x0
    & x4(z4)
    & y0
    & y4(z4)
    & z0
    & z2
  )

def asu_070(): # F d d d :2
  return (direct_space_asu('-F 2uv 2vw')
    & x0(-y0(z2))
    & x8(-z8 & z1*5/8)
    & ~y8(-z1*3/8 & z1*7/8)
    & y8(-z8 & z1*5/8)
    & z0
    & +z1
  )

def asu_071(): # I m m m
  return (direct_space_asu('-I 2 2')
    & x0
    & x4(y4(z4))
    & y0
    & y2
    & z0
    & z2
  )

def asu_072(): # I b a m
  return (direct_space_asu('-I 2 2c')
    & x0(z4)
    & x4(y4(z4))
    & y0(z4)
    & y2(z4)
    & z0
    & z2
  )

def asu_073(): # I b c a
  return (direct_space_asu('-I 2b 2c')
    & x0(y4)
    & x4(z4(y4))
    & y0(z4)
    & y2(-z4)
    & z0
    & +z2
  )

def asu_074(): # I m m a
  return (direct_space_asu('-I 2b 2')
    & x0
    & x4(-z4 & z1*3/4)
    & y0(z2)
    & y4
    & z0
    & +z1
  )

def asu_075(): # P 4
  return (direct_space_asu('P 4')
    & x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z1
  )

def asu_076(): # P 41 (enantiomorph of 76)
  return (direct_space_asu('P 4w')
    & x0(+z4)
    & x2(+z4)
    & y0(+z1*3/4)
    & y2(+z1*3/4)
    & z0
    & +z1
  )

def asu_077(): # P 42
  return (direct_space_asu('P 4c')
    & x0(+z2)
    & x2(+z2)
    & y0(+z2)
    & y2(+z2)
    & z0
    & +z1
  )

def asu_078(): # P 43 (enantiomorph of 78)
  return apply_change_of_basis(78)

def asu_079(): # I 4
  return (direct_space_asu('I 4')
    & x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z2
  )

def asu_080(): # I 41
  return (direct_space_asu('I 4bw')
    & x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0
    & +z4
  )

def asu_081(): # P -4
  return (direct_space_asu('P -4')
    & x0(-y0(z2))
    & x2
    & y0
    & y2(-x2(z2))
    & z0
    & +z1
  )

def asu_082(): # I -4
  return (direct_space_asu('I -4')
    & x0(z0(-y0))
    & x2(-y0(z4))
    & y0
    & y2(-x0(z4))
    & z0
    & z2(-y0)
  )

def asu_083(): # P 4/m
  return (direct_space_asu('-P 4')
    & x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z2
  )

def asu_084(): # P 42/m
  return (direct_space_asu('-P 4c')
    & x0(-y0(z4))
    & x2
    & y0
    & y2(-x2(z4))
    & z0
    & z2
  )

def asu_085(): # P 4/n :2
  return (direct_space_asu('-P 4a')
    & ~x4(-~y4)
    & x4(z0(-~y4) & z2(-~y4))
    & ~y4
    & y4(-x4)
    & z0(-y0(-x0))
    & z2(-y0(-x0))
  )

def asu_086(): # P 42/n :2
  return (direct_space_asu('-P 4bc')
    & ~x4(-~y4(z4))
    & x4(z0(-~y4) & z2(+-~y4))
    & ~y4
    & y4(-x4(z4))
    & z0(-y0(-x0))
    & z2(-y0(-x0))
  )

def asu_087(): # I 4/m
  return (direct_space_asu('-I 4')
    & x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z4(y4(x4) & x2(-y0))
  )

def asu_088(): # I 41/a :2
  return (direct_space_asu('-I 4ad')
    & x0
    & x4
    & y0(-x0(z2) | -x4(+z4))
    & y4(-x0(-z8 & z1*5/8))
    & z0
    & +z1
  )

def asu_089(): # P 4 2 2
  return (direct_space_asu('P 4 2')
    & x0(p0)
    & x2
    & y0
    & y2(-x2)
    & z0(p0)
    & z2(p0)
  )

def asu_090(): # P 4 21 2
  return (direct_space_asu('P 4ab 2ab')
    & x0
    & x2(-y0)
    & y0
    & y2(-x0)
    & z0(p0)
    & z2(p0)
  )

def asu_091(): # P 41 2 2 (enantiomorph of 95)
  return (direct_space_asu('P 4w 2c')
    & x0(z8(-y0))
    & +x1
    & y0
    & +y1
    & z0(x2)
    & z8(m1)
  )

def asu_092(): # P 41 21 2 (enantiomorph of 96)
  return (direct_space_asu('P 4abw 2nw')
    & x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z8(-y2)
  )

def asu_093(): # P 42 2 2
  return (direct_space_asu('P 4c 2')
    & x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0(y2)
    & z4(-p0 & m1)
  )

def asu_094(): # P 42 21 2
  return (direct_space_asu('P 4n 2n')
    & x0(-y0)
    & x2(z2(-y2))
    & y0(z2(-x0))
    & +y2
    & z0(p0)
    & z2(p0)
  )

def asu_095(): # P 43 2 2 (enantiomorph of 91)
  return apply_change_of_basis(95)

def asu_096(): # P 43 21 2 (enantiomorph of 92)
  # cannot be superimposed with enantiomorphic asu 92
  return (direct_space_asu('P 4nw 2abw')
    & x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z8(-x2)
  )

def asu_097(): # I 4 2 2
  return (direct_space_asu('I 4 2')
    & x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0(p0)
    & z4(m2)
  )

def asu_098(): # I 41 2 2
  return (direct_space_asu('I 4bw 2bw')
    & x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0(m1 & -p0)
    & z8(-y4 & y1*3/4)
  )

def asu_099(): # P 4 m m
  return (direct_space_asu('P 4 -2')
    & x0
    & y2
    & z0
    & +z1
    & -p0
  )

def asu_100(): # P 4 b m
  return (direct_space_asu('P 4 -2ab')
    & x0(-y0)
    & y0
    & z0
    & +z1
    & m2
  )

def asu_101(): # P 42 c m
  return (direct_space_asu('P 4c -2c')
    & x0(+z2)
    & y2(+z2)
    & z0
    & +z1
    & -p0
  )

def asu_102(): # P 42 n m
  return (direct_space_asu('P 4n -2n')
    & x0(+z2)
    & y2(+z2)
    & z0
    & +z1
    & -p0
  )

def asu_103(): # P 4 c c
  return (direct_space_asu('P 4 -2c')
    & x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z2
  )

def asu_104(): # P 4 n c
  return (direct_space_asu('P 4 -2n')
    & x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z2
  )

def asu_105(): # P 42 m c
  return (direct_space_asu('P 4c -2')
    & x0
    & x2
    & y0
    & y2
    & z0
    & +z2
  )

def asu_106(): # P 42 b c
  return (direct_space_asu('P 4c -2ab')
    & x0(-y0)
    & x2
    & y0
    & +y2
    & z0
    & +z2
  )

def asu_107(): # I 4 m m
  return (direct_space_asu('I 4 -2')
    & x0
    & y2
    & z0
    & +z2
    & -p0
  )

def asu_108(): # I 4 c m
  return (direct_space_asu('I 4 -2c')
    & x0(-y0)
    & y0
    & z0
    & +z2
    & m2
  )

def asu_109(): # I 41 m d
  return (direct_space_asu('I 4bw -2')
    & x0
    & x2
    & y0
    & y2
    & z0
    & +z4
  )

def asu_110(): # I 41 c d
  return (direct_space_asu('I 4bw -2c')
    & x0(-y0)
    & x2
    & y0
    & +y2
    & z0
    & +z4
  )

def asu_111(): # P -4 2 m
  return (direct_space_asu('P -4 2')
    & x0(z2)
    & y2(z2)
    & z0
    & +z1
    & -p0
  )

def asu_112(): # P -4 2 c
  return (direct_space_asu('P -4 2c')
    & x0(z4 & z0(-y0))
    & x2(z4)
    & y0(z4)
    & y2(z4 & z0(-x2))
    & z0
    & +z2
  )

def asu_113(): # P -4 21 m
  return (direct_space_asu('P -4 2ab')
    & x0(-y0(z2))
    & y0
    & z0
    & +z1
    & m2
  )

def asu_114(): # P -4 21 c
  return (direct_space_asu('P -4 2n')
    & x0(-y0)
    & x2
    & y0
    & y2(-x2(-z0))
    & z0
    & +z2
  )

def asu_115(): # P -4 m 2
  return (direct_space_asu('P -4 -2')
    & x0
    & x2
    & y0
    & y2
    & z0(p0)
    & z2(p0)
  )

def asu_116(): # P -4 c 2
  return (direct_space_asu('P -4 -2c')
    & x0(y2)
    & x2(y2 & z0(-y2))
    & y0
    & +y1
    & z0(y2 & y0(-x0))
    & z4(m1 & -p0)
  )

def asu_117(): # P -4 b 2
  return (direct_space_asu('P -4 -2ab')
    & x0(z0(-y0) & z2(-y0))
    & x2(-y0)
    & y0
    & +y2
    & z0(m2)
    & z2(m2)
  )

def asu_118(): # P -4 n 2
  return (direct_space_asu('P -4 -2n')
    & x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0(y2(-x2) & x0(-y0))
    & z4(~p2 & -m2)
  )

def asu_119(): # I -4 m 2
  return (direct_space_asu('I -4 -2')
    & x0
    & x2
    & y0
    & y2
    & z0(p0)
    & z4(m2)
  )

def asu_120(): # I -4 c 2
  return (direct_space_asu('I -4 -2c')
    & x0(z0(-y0))
    & x2(-y0)
    & y0
    & +y2
    & z0(m2)
    & z4(p0)
  )

def asu_121(): # I -4 2 m
  return (direct_space_asu('I -4 2')
    & x0
    & y2(-x0(z4))
    & z0
    & z2(-x0)
    & -p0
  )

def asu_122(): # I -4 2 d
  return (direct_space_asu('I -4 2bw')
    & x0(y2 & z0(-y0))
    & x2(y2)
    & y0
    & +y1
    & z0(y2(-x2))
    & z8(-y4 & y1*3/4)
  )

def asu_123(): # P 4/m m m
  return (direct_space_asu('-P 4 2')
    & x0
    & y2
    & z0
    & z2
    & -p0
  )

def asu_124(): # P 4/m c c
  return (direct_space_asu('-P 4 2c')
    & x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z4(p0)
  )

def asu_125(): # P 4/n b m :2
  return (direct_space_asu('-P 4a 2b')
    & ~x4(-~y4)
    & ~y4
    & z0(p0)
    & z2(p0)
    & -m0
  )

def asu_126(): # P 4/n n c :2
  return (direct_space_asu('-P 4a 2bc')
    & ~x4(-~y4)
    & x4
    & ~y4
    & y4(-x4)
    & z0(-y0(-x0) & x4(-~y4))
    & z4(p0)
  )

def asu_127(): # P 4/m b m
  return (direct_space_asu('-P 4 2ab')
    & x0(-y0)
    & y0
    & z0
    & z2
    & m2
  )

def asu_128(): # P 4/m n c
  return (direct_space_asu('-P 4 2n')
    & x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & z4(m2)
  )

def asu_129(): # P 4/n m m :2
  return (direct_space_asu('-P 4a 2a')
    & ~x4
    & y4
    & z0(-m0)
    & z2(-m0)
    & -p0
  )

def asu_130(): # P 4/n c c :2
  return (direct_space_asu('-P 4a 2ac')
    & ~x4(-~y4)
    & x4(z0(-~y4))
    & ~y4
    & y4(-x4)
    & z0(-y0(-x0))
    & z4(-m0)
  )

def asu_131(): # P 42/m m c
  return (direct_space_asu('-P 4c 2')
    & x0
    & x2
    & y0
    & y2
    & z0
    & z4(p0)
  )

def asu_132(): # P 42/m c m
  return (direct_space_asu('-P 4c 2c')
    & x0(z4)
    & y2(z4)
    & z0
    & z2
    & -p0
  )

def asu_133(): # P 42/n b c :2
  return (direct_space_asu('-P 4ac 2b')
    & ~x4
    & x4(-z0 | -~y4)
    & ~y4
    & +y4
    & z0(-y0(-x0))
    & z4(p0)
  )

def asu_134(): # P 42/n n m :2
  return (direct_space_asu('-P 4ac 2bc')
    & ~x4(z4)
    & ~y4(z4)
    & z0(p0)
    & z2(p0)
    & -m0
  )

def asu_135(): # P 42/m b c
  return (direct_space_asu('-P 4c 2ab')
    & x0(z4(-y0))
    & x2(-y0)
    & y0
    & +y2
    & z0
    & z4(m2)
  )

def asu_136(): # P 42/m n m
  return (direct_space_asu('-P 4n 2n')
    & x0(z4)
    & y2(z4(-x0))
    & z0
    & z2
    & -p0
  )

def asu_137(): # P 42/n m c :2
  return (direct_space_asu('-P 4ac 2a')
    & ~x4
    & x4
    & ~y4
    & y4
    & z0(-y0(-x0))
    & z4(-m0)
  )

def asu_138(): # P 42/n c m :2
  return (direct_space_asu('-P 4ac 2ac')
    & ~x4
    & y4(-~x4(z4))
    & z0(-m0)
    & z2(-m0 & ~x4(-y4))
    & -p0
  )

def asu_139(): # I 4/m m m
  return (direct_space_asu('-I 4 2')
    & x0
    & y2
    & z0
    & z4(m2)
    & -p0
  )

def asu_140(): # I 4/m c m
  return (direct_space_asu('-I 4 2c')
    & x0(-y0)
    & y0
    & z0
    & z4(p0)
    & m2
  )

def asu_141(): # I 41/a m d :2
  return (direct_space_asu('-I 4bd 2')
    & x0
    & x2
    & ~y4
    & y4
    & z0(-y0)
    & z8(-p4)
  )

def asu_142(): # I 41/a c d :2
  return (direct_space_asu('-I 4bd 2c')
    & x0(z8(-~y4) & z0(-y0))
    & x2(-~y4)
    & ~y4
    & +y4
    & z0(x4)
    & z8(m4)
  )

def asu_143(): # P 3
  return (direct_space_asu('P 3')
    & x0(-y0)
    & y0
    & z0
    & +z1
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_144(): # P 31 (enantiomorph of 145)
  return (direct_space_asu('P 31')
    & x0
    & +x1
    & y0
    & +y1
    & z0
    & +z3
  )

def asu_145(): # P 32 (enantiomorph of 144)
  return apply_change_of_basis(145)

def asu_146(): # R 3 :H
  return (direct_space_asu('R 3')
    & x0(-y0)
    & y0
    & z0
    & +z3
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_147(): # P -3
  return (direct_space_asu('-P 3')
    & x0(-y0)
    & y0
    & z0(p0(-y0))
    & z2(p0(-y0))
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_148(): # R -3 :H
  return (direct_space_asu('-R 3')
    & x0(-y0)
    & y0
    & z0(p0(-y0))
    & z6(-h0(x3) | -k0(-y0 | -m1))
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_149(): # P 3 1 2
  return (direct_space_asu('P 3 2')
    & x0(-y0)
    & y0
    & z0(-h0 | -k0)
    & z2(-h0 | -k0)
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_150(): # P 3 2 1
  return (direct_space_asu('P 3 2"')
    & x0(-y0)
    & y0
    & z0(p0)
    & z2(p0)
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_151(): # P 31 1 2 (enantiomorph of 153)
  return (direct_space_asu('P 31 2 (0 0 4)')
    & x0
    & +x1
    & y0
    & +y1
    & z0(-h0 | -h1)
    & z6(-k0 | -k1)
  )

def asu_152(): # P 31 2 1 (enantiomorph of 154)
  return (direct_space_asu('P 31 2"')
    & x0
    & +x1
    & y0
    & +y1
    & z0(-p0)
    & z6(-p0)
  )

def asu_153(): # P 32 1 2 (enantiomorph of 151)
  # cannot be superimposed with enantiomorphic asu 151
  return (direct_space_asu('P 32 2 (0 0 2)')
    & x0
    & +x1
    & y0
    & +y1
    & z0(-h0 | -h1)
    & z6(x0(-y0) & m1)
 )

def asu_154(): # P 32 2 1 (enantiomorph of 152)
  return apply_change_of_basis(154)

def asu_155(): # R 3 2 :H
  return (direct_space_asu('R 3 2"')
    & x0(-y0)
    & y0
    & z0(p0)
    & z6(x3 & ~p3)
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_156(): # P 3 m 1
  return (direct_space_asu('P 3 -2"')
    & z0
    & +z1
    & h0
    & m1
    & k0
  )

def asu_157(): # P 3 1 m
  return (direct_space_asu('P 3 -2')
    & y0
    & z0
    & +z1
    & k1
    & m1(y3)
    & p0
  )

def asu_158(): # P 3 c 1
  return (direct_space_asu('P 3 -2"c')
    & x0(-y0)
    & y0
    & z0
    & +z2
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_159(): # P 3 1 c
  return (direct_space_asu('P 3 -2c')
    & x0(-y0)
    & y0
    & z0
    & +z2
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_160(): # R 3 m :H
  return (direct_space_asu('R 3 -2"')
    & z0
    & +z3
    & h0
    & m1
    & k0
  )

def asu_161(): # R 3 c :H
  return (direct_space_asu('R 3 -2"c')
    & x0(-y0)
    & y0
    & z0
    & +z6
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_162(): # P -3 1 m
  return (direct_space_asu('-P 3 2')
    & y0
    & z0(-h0)
    & z2(-h0)
    & k1
    & m1(y3)
    & p0
  )

def asu_163(): # P -3 1 c
  return (direct_space_asu('-P 3 2c')
    & x0(-y0)
    & y0
    & z0(p0(-y0))
    & z4(-h0 | -k0)
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_164(): # P -3 m 1
  return (direct_space_asu('-P 3 2"')
    & y0(z2)
    & z0
    & +z1
    & k1
    & -h0
  )

def asu_165(): # P -3 c 1
  return (direct_space_asu('-P 3 2"c')
    & x0(-y0)
    & y0
    & z0(p0(-y0))
    & z4(p0)
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_166(): # R -3 m :H
  return (direct_space_asu('-R 3 2"')
    & z0(p0)
    & z6(x3)
    & h0
    & m1
    & k0
  )

def asu_167(): # R -3 c :H
  return (direct_space_asu('-R 3 2"c')
    & x0(-y0)
    & y0
    & z0(p0(-y0))
    & z12(y3 & p3)
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_168(): # P 6
  return (direct_space_asu('P 6')
    & y0
    & z0
    & +z1
    & k1
    & m1(y3)
    & p0(-y0)
  )

def asu_169(): # P 61 (enantiomorph of 170)
  return (direct_space_asu('P 61')
    & x0
    & +x1
    & y0
    & +y1
    & z0
    & +z6
  )

def asu_170(): # P 65 (enantiomorph of 169)
  return apply_change_of_basis(170)

def asu_171(): # P 62 (enantiomorph of 172)
  return (direct_space_asu('P 62')
    & x1(y2)
    & y0(x2)
    & z0
    & +z3
    & p0(y2)
  )

def asu_172(): # P 64 (enantiomorph of 171)
  return apply_change_of_basis(172)

def asu_173(): # P 63
  return (direct_space_asu('P 6c')
    & x0(-y0)
    & y0
    & z0
    & +z2
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_174(): # P -6
  return (direct_space_asu('P -6')
    & x0(-y0)
    & y0
    & z0
    & z2
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_175(): # P 6/m
  return (direct_space_asu('-P 6')
    & y0
    & z0
    & z2
    & k1
    & m1(y3)
    & p0(-y0)
  )

def asu_176(): # P 63/m
  return (direct_space_asu('-P 6c')
    & x0(-y0)
    & y0
    & z0(p0(-y0))
    & z4
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_177(): # P 6 2 2
  return (direct_space_asu('P 6 2')
    & y0
    & z0(-h0)
    & z2(-h0)
    & k1
    & m1(y3)
    & p0(-y0)
  )

def asu_178(): # P 61 2 2 (enantiomorph of 179)
  return (direct_space_asu('P 61 2 (0 0 5)')
    & x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z12(-h0 | -h1)
  )

def asu_179(): # P 65 2 2 (enantiomorph of 178)
  # cannot be superimposed with enantiomorphic asu 178
  return (direct_space_asu('P 65 2 (0 0 1)')
    & x0
    & +x1
    & y0
    & +y1
    & z0(p0)
    & z12(m1 & x0(-y0))
  )

def asu_180(): # P 62 2 2 (enantiomorph of 181)
  return (direct_space_asu('P 62 2 (0 0 4)')
    & x1(y2)
    & y0(x2)
    & z0(k1)
    & z6(-h0)
    & p0(y2)
  )

def asu_181(): # P 64 2 2 (enantiomorph of 180)
  result = apply_change_of_basis(181)
  assert result.hall_symbol == " P 64 2 (x,y,z+1/6)"
  result.hall_symbol = "P 64 2 (0 0 2)" # Int. Tab. Vol. B compatibility
  return result

def asu_182(): # P 63 2 2
  return (direct_space_asu('P 6c 2c')
    & x0(-y0)
    & y0
    & z0(p0)
    & z4(-h0 | -k0)
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_183(): # P 6 m m
  return (direct_space_asu('P 6 -2')
    & y0
    & z0
    & +z1
    & k1
    & -h0
  )

def asu_184(): # P 6 c c
  return (direct_space_asu('P 6 -2c')
    & y0
    & z0
    & +z2
    & k1
    & m1(y3)
    & p0(-y0)
  )

def asu_185(): # P 63 c m
  return (direct_space_asu('P 6c -2')
    & y0
    & z0
    & +z2
    & k1
    & m1(y3)
    & p0
  )

def asu_186(): # P 63 m c
  return (direct_space_asu('P 6c -2c')
    & y0(+z2)
    & z0
    & +z1
    & k1
    & -h0
  )

def asu_187(): # P -6 m 2
  return (direct_space_asu('P -6 2')
    & z0
    & z2
    & h0
    & m1
    & k0
  )

def asu_188(): # P -6 c 2
  return (direct_space_asu('P -6c 2')
    & x0(-y0)
    & y0
    & z0(-h0 | -k0)
    & z4
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_189(): # P -6 2 m
  return (direct_space_asu('P -6 -2')
    & y0
    & z0
    & z2
    & k1
    & m1(y3)
    & p0
  )

def asu_190(): # P -6 2 c
  return (direct_space_asu('P -6c -2c')
    & x0(-y0)
    & y0
    & z0(p0)
    & z4
    & k1
    & m1(-h1 | -k1)
    & h1
  )

def asu_191(): # P 6/m m m
  return (direct_space_asu('-P 6 2')
    & y0
    & z0
    & z2
    & k1
    & -h0
  )

def asu_192(): # P 6/m c c
  return (direct_space_asu('-P 6 2c')
    & y0
    & z0
    & z4(-h0)
    & k1
    & m1(y3)
    & p0(-y0)
  )

def asu_193(): # P 63/m c m
  return (direct_space_asu('-P 6c 2')
    & y0
    & z0(-h0)
    & z4
    & k1
    & m1(y3)
    & p0
  )

def asu_194(): # P 63/m m c
  return (direct_space_asu('-P 6c 2c')
    & z0(p0)
    & z4
    & h0
    & m1
    & k0
  )

def asu_195(): # P 2 3
  return (direct_space_asu('P 2 2 3')
    & z0(y2 & x2)
    & m1(-y2)
    & zy0(-zx0)
    & zx0
  )

def asu_196(): # F 2 3
  return (direct_space_asu('F 2 2 3')
    & p0(m2)
    & ~xz2(-zy0)
    & zx2(yz0)
    & -yz0
    & zy0
  )

def asu_197(): # I 2 3
  return (direct_space_asu('I 2 2 3')
    & z0(x2)
    & p0(-zy0)
    & +m1
    & zy0
  )

def asu_198(): # P 21 3
  return (direct_space_asu('P 2ac 2ab 3')
    & x0(-y0)
    & x2
    & y2(+z0 & x2(+z2))
    & zx2(m2)
    & zx0(p0)
    & -yz0
    & zy0
  )

def asu_199(): # I 21 3
  return (direct_space_asu('I 2b 2c 3')
    & x2(-y4)
    & y2(-z4)
    & z0(x4)
    & zx0(-zy0(+x2))
    & zy0
  )

def asu_200(): # P m -3
  return (direct_space_asu('-P 2 2 3')
    & x2
    & y2
    & z0
    & zx0(-zy0)
    & zy0
  )

def asu_201(): # P n -3 :2
  return (direct_space_asu('-P 2ab 2bc 3')
    & ~z4(x4)
    & p0(-zy0(-x0))
    & m2(-zy0(x2))
    & zy0
  )

def asu_202(): # F m -3
  return (direct_space_asu('-F 2 2 3')
    & z0
    & p0(x4)
    & ~xz2(-zy0)
    & zy0
  )

def asu_203(): # F d -3 :2
  return (direct_space_asu('-F 2uv 2vw 3')
    & p0(-zy0(-x0))
    & m4(-zy0 | -yz4(~z4))
    & zy0
    & yz4
  )

def asu_204(): # I m -3
  return (direct_space_asu('-I 2 2 3')
    & x2
    & z0
    & p0(-zy0(x4))
    & zy0
  )

def asu_205(): # P a -3
  return (direct_space_asu('-P 2ac 2ab 3')
    & x2(-z0(-zy0))
    & y2(-zy0)
    & z0
    & zx0(-zy0)
    & zy0
  )

def asu_206(): # I a -3
  return (direct_space_asu('-I 2b 2c 3')
    & z0(x4)
    & zx0(-zy0)
    & ~xz2
    & zy0
    & ~yz2(-zx0)
  )

def asu_207(): # P 4 3 2
  return (direct_space_asu('P 4 2 3')
    & z0(x2)
    & p0
    & m1(-p0)
    & zy0(x2)
  )

def asu_208(): # P 42 3 2
  return (direct_space_asu('P 4n 2 3')
    & zx0(-zy0)
    & -xz0(yz0)
    & zx2(y4)
    & ~xz2(y4)
    & zy0
    & -yz0
    & zy2(-x4)
    & ~yz2(-x4)
  )

def asu_209(): # F 4 3 2
  return (direct_space_asu('F 4 2 3')
    & p0(z0)
    & m2(z0)
    & zy0
    & -yz0(-zy0)
  )

def asu_210(): # F 41 3 2
  return (direct_space_asu('F 4d 2 3')
    & y8(-~xz4)
    & z8(m4)
    & p0(-zx0)
    & m2(-~xz2)
    & -yz0(z0)
    & zx0
    & ~xz2
  )

def asu_211(): # I 4 3 2
  return (direct_space_asu('I 4 2 3')
    & z0(p0)
    & zx0(-zy0)
    & ~xz2(y4)
    & zy0
    & ~yz2(-x4)
  )

def asu_212(): # P 43 3 2 (enantiomorph of 213)
  return (direct_space_asu('P 4acd 2ab 3')
    & zx2
    & -yz0(-zx2)
    & ~yz2(tx0)
    & -tx0(x8)
    & ty0(y8)
    & tz2(-x1*3/8)
  )

def asu_213(): # P 41 3 2 (enantiomorph of 212)
  return apply_change_of_basis(213)

def asu_214(): # I 41 3 2
  return (direct_space_asu('I 4bd 2c 3')
    & x8(~yz4)
    & y8(~xz4)
    & ~y8(-~zx1/4)
    & -zx0(zy0)
    & -zy0
    & ~zy4(-y0)
    & dy8(~p4)
  )

def asu_215(): # P -4 3 m
  return (direct_space_asu('P -4 2 3')
    & z0(x2)
    & p0
    & m1
    & zy0
  )

def asu_216(): # F -4 3 m
  return (direct_space_asu('F -4 2 3')
    & p0
    & m2
    & zy0
    & -yz0
  )

def asu_217(): # I -4 3 m
  return (direct_space_asu('I -4 2 3')
    & x2(-z0(y4))
    & z0
    & p0
    & zy0
  )

def asu_218(): # P -4 3 n
  return (direct_space_asu('P -4n 2 3')
    & x2(-z0(y4))
    & y2(-z0(x4))
    & z0
    & zx0(-zy0)
    & zy0
  )

def asu_219(): # F -4 3 c
  return (direct_space_asu('F -4a 2 3')
    & p0
    & m2(-p0(z0))
    & zy0
    & -yz0(-p0 | -zy0(x4))
  )

def asu_220(): # I -4 3 d
  return (direct_space_asu('I -4bd 2c 3')
    & -x4(-z0(-y1*3/8))
    & x2
    & -y4(-x2(-z8))
    & y2(-z4)
    & z0
    & zx0(-zy0)
    & zy0
  )

def asu_221(): # P m -3 m
  return (direct_space_asu('-P 4 2 3')
    & x2
    & z0
    & p0
    & zy0
  )

def asu_222(): # P n -3 n :2
  return (direct_space_asu('-P 4a 2bc 3')
    & x34(z4(y2) | -zy0)
    & -z4
    & p0(-zy0(z2))
    & zy0
  )

def asu_223(): # P m -3 n
  return (direct_space_asu('-P 4n 2 3')
    & z0
    & zx0(-zy0)
    & ~xz2(y4)
    & zy0
    & ~yz2(x4)
  )

def asu_224(): # P n -3 m :2
  return (direct_space_asu('-P 4bc 2bc 3')
    & p0
    & ~xz1(y2)
    & zx2(y2)
    & -~yz2
    & zy0
  )

def asu_225(): # F m -3 m
  return (direct_space_asu('-F 4 2 3')
    & z0
    & p0
    & m2
    & zy0
  )

def asu_226(): # F m -3 c
  return (direct_space_asu('-F 4a 2 3')
    & z0
    & p0
    & m2(-p0)
    & zy0(x4)
  )

def asu_227(): # F d -3 m :2
  return (direct_space_asu('-F 4vw 2vw 3')
    & -y0(-xz0)
    & p0
    & m4
    & yz4
    & zy0
  )

def asu_228(): # F d -3 c :2
  return (direct_space_asu('-F 4ud 2vw 3')
    & -y0(zx1/4)
    & p0(-zy0)
    & m4
    & yz4(-zy0(x8))
    & zy0
  )

def asu_229(): # I m -3 m
  return (direct_space_asu('-I 4 2 3')
    & z0
    & p0
    & ~xz2(y4)
    & zy0
  )

def asu_230(): # I a -3 d
  return (direct_space_asu('-I 4bd 2c 3')
    & x8(~zy4 & ~yz4)
    & ~x8(y0(-z4))
    & y8(-~xz4)
    & ~y8(~zx1/4)
    & z4(y0)
    & -zx0
    & -xz0(-z0)
    & -zy0(zx0)
    & -yz0
  )

def get_asu(space_group_number):
  return eval("asu_%03d" % space_group_number)()

if (__name__ == "__main__"):
  for i in range(1,231):
    get_asu(i).show_summary()
