from __future__ import absolute_import, division, print_function
expected_output_test1  = """Left cosets of :
  subgroup  H: P 3
  and group G: P 6 2 2

  Coset number :     0   (all operators from H)

               x,y,z                 h,k,l   Rotation:    1 ; direction:  (0, 0, 0) ; screw/glide:    (0,0,0)
            -y,x-y,z              k,-h-k,l   Rotation:    3 ; direction:  (0, 0, 1) ; screw/glide:    (0,0,0)
           -x+y,-x,z              -h-k,h,l   Rotation:    3 ; direction:  (0, 0, 1) ; screw/glide:    (0,0,0)

  Coset number :     1   (H+coset[1] = P 6)

             -x,-y,z               -h,-k,l   Rotation:    2 ; direction:  (0, 0, 1) ; screw/glide:    (0,0,0)
             x-y,x,z              h+k,-h,l   Rotation:    6 ; direction:  (0, 0, 1) ; screw/glide:    (0,0,0)
            y,-x+y,z              -k,h+k,l   Rotation:    6 ; direction:  (0, 0, 1) ; screw/glide:    (0,0,0)

  Coset number :     2   (H+coset[2] = P 3 2 1)

           x-y,-y,-z             h,-h-k,-l   Rotation:    2 ; direction:  (1, 0, 0) ; screw/glide:    (0,0,0)
          -x,-x+y,-z             -h-k,k,-l   Rotation:    2 ; direction:  (0, 1, 0) ; screw/glide:    (0,0,0)
              y,x,-z                k,h,-l   Rotation:    2 ; direction:  (1, 1, 0) ; screw/glide:    (0,0,0)

  Coset number :     3   (H+coset[3] = P 3 1 2)

            -y,-x,-z              -k,-h,-l   Rotation:    2 ; direction: (-1, 1, 0) ; screw/glide:    (0,0,0)
           -x+y,y,-z             -h,h+k,-l   Rotation:    2 ; direction:  (1, 2, 0) ; screw/glide:    (0,0,0)
            x,x-y,-z             h+k,-k,-l   Rotation:    2 ; direction:  (2, 1, 0) ; screw/glide:    (0,0,0)
"""

expected_output_test2 ="""Coset decomposition not successfull.
Group P 1 2 1 might not be a subgroup of P 3
Sorry.....
"""

class tmp_out:
  def __init__(self):
    self.store = """"""
  def write(self, txt):
    self.store+=txt



from cctbx.sgtbx import show_cosets

def test():
  test1 = tmp_out()
  show_cosets.run( "P3", "P622", out=test1 )
  assert test1.store == expected_output_test1

  test2 = tmp_out()
  show_cosets.run( "P2" , "P3", out=test2)
  assert test2.store == expected_output_test2
  print("OK")

if __name__ == "__main__":
  test()
