      program prog
      implicit none
      byte vbyte
      logical vlogical
c!fem logical*1 vlogical1
c!fem logical*2 vlogical2
c!fem logical*4 vlogical4
c!fem logical*8 vlogical8
      integer vinteger
      integer*1 vinteger1
      integer*2 vinteger2
      integer*4 vinteger4
      integer*8 vinteger8
      real vreal
      real*4 vreal4
      real*8 vreal8
      double precision vdp
c!fem real*16 vreal16
      complex vcomplex
      complex*8 vcomplex8
      complex*16 vcomplex16
      double complex vdc
c!fem complex*32 vcomplex32
      integer i
c
      vbyte = x'FF'
      write(6, *) vbyte .eq. -1
      vbyte = x'80'
      write(6, *) vbyte .eq. -128
      vbyte = x'7F'
      write(6, *) vbyte .eq. 127
c
      vlogical = .false.
      write(6, *) vlogical
      vlogical = .true.
      write(6, *) vlogical
c!fem vlogical1 = .false.
c!fem write(6, *) vlogical1
c!fem vlogical1 = .true.
c!fem write(6, *) vlogical1
c!fem vlogical2 = .false.
c!fem write(6, *) vlogical2
c!fem vlogical2 = .true.
c!fem write(6, *) vlogical2
c!fem vlogical4 = .false.
c!fem write(6, *) vlogical4
c!fem vlogical4 = .true.
c!fem write(6, *) vlogical4
c!fem vlogical8 = .false.
c!fem write(6, *) vlogical8
c!fem vlogical8 = .true.
c!fem write(6, *) vlogical8
c
      vinteger = x'FFFFFFFF'
      write(6, *) vinteger
      vinteger = x'80000000'
      write(6, *) vinteger
      vinteger = x'7FFFFFFF'
      write(6, *) vinteger
      vinteger1 = x'FF'
      write(6, *) vinteger1
      vinteger1 = x'80'
      write(6, *) vinteger1
      vinteger1 = x'7F'
      write(6, *) vinteger1
      vinteger2 = x'FFFF'
      write(6, *) vinteger2
      vinteger2 = x'8000'
      write(6, *) vinteger2
      vinteger2 = x'7FFF'
      write(6, *) vinteger2
      vinteger4 = x'FFFFFFFF'
      write(6, *) vinteger4
      vinteger4 = x'80000000'
      write(6, *) vinteger4
      vinteger4 = x'7FFFFFFF'
      write(6, *) vinteger4
      vinteger8 = x'FFFFFFFF'
      write(6, *) vinteger8
      vinteger8 = x'80000000'
      write(6, *) vinteger8
      vinteger8 = x'7FFFFFFF'
      write(6, *) vinteger8
c
      vreal = 1.17550e-38
      write(6, *) vreal
      vreal = 3.40281e+38
      write(6, *) vreal
      vreal4 = 1.17550e-38
      write(6, *) vreal4
      vreal4 = 3.40281e+38
      write(6, *) vreal4
      vreal8 = 2.22508d-308
      write(6, *) vreal8
      vreal8 = 1.79768d+308
      write(6, *) vreal8
      vdp = 2.22508d-308
      write(6, *) vdp
      vdp = 1.79768d+308
      write(6, *) vdp
c!fem vreal16 = 1
c!fem do i=1,16
c!fem   vreal16 = vreal16 * 2.22508d-308
c!fem enddo
c!fem write(6, *) vreal16
c!fem vreal16 = 1
c!fem do i=1,16
c!fem   vreal16 = vreal16 * 1.79768d+308
c!fem enddo
c!fem write(6, *) vreal16
c
      vcomplex = cmplx(1., 2.e10)
      write(6, *) vcomplex
      vcomplex8 = cmplx(-3.e10, -4.)
      write(6, *) vcomplex8
      vcomplex16 = dcmplx(5.d0, 6.d10)
      write(6, *) vcomplex16
      vdc = dcmplx(-7.d10, -8.d0)
      write(6, *) vdc
c!fem vcomplex32 = dcmplx(9.d0, 1.d10)
c!fem write(6, *) vcomplex32
c
      end
