      program prog
      integer e
      double precision val
      write(6, '(a)') '1.0'
      do e=-3,1
        write(6, '(d10.3)') 1.0d0*10**e
      enddo
      write(6, '(a)') '0.999'
      do e=-3,1
        write(6, '(d10.3)') 0.999d0*10**e
      enddo
      write(6, '(a)') '0.9999'
      do e=-3,1
        write(6, '(d10.3)') 0.9999d0*10**e
      enddo
      write(6, '(a)') '1.23'
      do e=-3,1
        write(6, '(d10.3)') 1.23d0*10**e
      enddo
      write(6, '(a)') '1.23**100'
      do e=-3,1
        write(6, '(d15.3)') 1.23d0*1.d100*10.d0**e
      enddo
      do i=1,2
        if (i .eq. 1) then
          val = 1.234d0
          write(6, '(a)') '1.234 with -3p through 5p'
        else
          val = 0.d0
          write(6, '(a)') '0 with -3p through 5p'
        endif
        write(6, '(i2,1x,-3p,d10.3)') -3, val
        write(6, '(i2,1x,-2p,d10.3)') -2, val
        write(6, '(i2,1x,-1p,d10.3)') -1, val
        write(6, '(i2,1x,0pd10.3)') 0, val
        write(6, '(i2,1x,1p,d10.3)') 1, val
        write(6, '(i2,1x,2p,d10.3)') 2, val
        write(6, '(i2,1x,3p,d10.3)') 3, val
        write(6, '(i2,1x,4p,d10.3)') 4, val
        write(6, '(i2,1x,5p,d10.3)') 5, val
      enddo
      write(6, '(''a'',1x,d1.0)') 0.d0
      write(6, '(''b'',1x,d10.0)') 0.d0
      write(6, '(''c'',1x,d6.1)') 0.d0
      write(6, '(''d'',1x,d6.1)') 1.d0
      write(6, '(''e'',1x,d6.1)') -1.d0
      write(6, '(''f'',1x,d7.1)') 0.d0
      write(6, '(''g'',1x,d7.1)') 1.d0
      write(6, '(''h'',1x,d7.1)') -1.d0
      write(6, '(''i'',1x,1p,d7.1)') -1.d0
      write(6, '(''j'',1x,d8.1)') 0.d0
      write(6, '(''k'',1x,d8.1)') 1.d0
      write(6, '(''l'',1x,d8.1)') -1.d0
      write(6, '(''m'',1x,d8.2)') -1.d0
      write(6, '(''n'',1x,d8.3)') -1.d0
      write(6, '(''o'',1x,1p,d8.1)') -1.d0
      write(6, '(''p'',1x,1p,d8.2)') -1.d0
      write(6, '(''q'',1x,2p,d8.1)') -1.d0
      write(6, '(''r'',1x,3p,d8.1)') -1.d0
c
      val = -1.d0/2097152.d0
      write(6, '(0p,d30.24)') val
      write(6, '(0p,d30.23)') val
      write(6, '(1p,d30.23)') val
      write(6, '(10p,d30.23)') val
      write(6, '(20p,d30.23)') val
      write(6, '(24p,d30.23)') val
      write(6, '(25p,d30.23)') val
c
      write(6, '(0p,e9.2)') 1.234
      write(6, '(1p,e9.2)') 1.234
      write(6, '(2p,e9.2)') 1.234
      write(6, '(3p,e9.2)') 1.234
      write(6, '(4p,e9.2)') 1.234
      end
