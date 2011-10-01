      program data
      common/misc/ skip(10,2)
      double precision skip
      integer*4 i
      do i=1,10
        write(*,*) skip(i,1), skip(i,2)
      enddo
      end

      block data
      common/misc/ skip(10,2)
      double precision skip
      data ((skip(i,j),i=1,10),j=1,2)/0.,3.,8*0.,1.d10,1.d10,2.,7*0./
      end
