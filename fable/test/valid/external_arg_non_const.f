      program prog
      external exch_imp
      call sub(exch_imp)
      end

      subroutine sub(exch)
      external exch
      dimension nc(2)
      dimension nm(2)
      do i=1,2
        nc(i) = 13*i
      enddo
      call sub2(nc, nm, exch)
      write(6, *) nm
      end

      subroutine sub2(nc, nm, exch)
      dimension nc(2)
      dimension nm(2)
      call exch(nc, nm)
      end

      subroutine exch_imp(nc, nm)
      dimension nc(2)
      dimension nm(2)
      write(6, *) nc
      do i=1,2
        nm(i) = nc(i) + 19
      enddo
      end
