      program prog
      logical l1, l2, l3
      do i=0,1
        l1 = i .ne. 0
        do j=0,1
          l2 = j .ne. 0
          do k=0,1
            l3 = k .ne. 0
            write(6, *)
     &        .not. l1,
     &        l1 .and. l2,
     &        l1 .or. l2,
     &        l1 .eqv. l2,
     &        l1 .neqv. l2,
     &        .not. l1 .and. l2,
     &        .not. l1 .or. l2,
     &        (l1 .and. l2) .or. l3,
     &        l1 .or. (l2 .and. l3)
          enddo
        enddo
      enddo
      end
