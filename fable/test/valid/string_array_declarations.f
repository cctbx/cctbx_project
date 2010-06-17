      subroutine sub(a,b,c,d,e,f,g,h,i,j,k,l,m)
      character a
      character b*(*)
      character*(*) c
      character d*(*)(*)
      character e*(*)(*)
      character f*(*)(*)
      character g*(*)(*)
      character h*(*)
      character i*(*)
      character j*(*)(*)
      character k*(*)(*)
      character l*(*)(*)
      character m*(*)(2,*)
      return
      end

      program prog
      integer z
      parameter(z=3)
      character a
      character*3 b
      character c*3
      character d(2)
      character*3 e(2)
      character f*3(2)
      character g(2)*3
      character*(z) h
      character i*(z)
      character*(z) j(2)
      character k*(z)(2)
      character l(2)*(z)
      character m(2,4)*(z)
      call sub(a,b,c,d,e,f,g,h,i,j,k,l,m)
      end
