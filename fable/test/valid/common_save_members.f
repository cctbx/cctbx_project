      program prog
      save
      integer si
      character*5 sc
      integer sai(2)
      character*3 sas(2)
      integer ci
      character*8 cc
      integer cai(2)
      character*4 cas(2)
      common /globals/ ci, cc, cai, cas
      si = 12
      ci = 34
      sc = 'WeRtY'
      cc = 'uIoPqWeR'
      do i=1,2
        sai(i) = i + 37
        cai(i) = i + 41
      enddo
      sas(1) = 'xYz'
      sas(2) = 'EfG'
      cas(1) = 'uvWx'
      cas(2) = 'PqrS'
      write(6, *) si
      write(6, *) sc
      write(6, *) sai
      write(6, *) sas
      write(6, *) ci
      write(6, *) cc
      write(6, *) cai
      write(6, *) cas
      end
