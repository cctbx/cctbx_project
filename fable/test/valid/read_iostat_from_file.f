      program prog
      character*(1) buf1
      character*(1) buf2
      character*(2) buf3
      character*(3) buf33
      character*(2) buf4(2)
      integer*4 k, l, m(2)
      double precision x
      integer ios, i

      open(10, file='fable_tmp_quuqu3ee', status='unknown')
      write(10, '(a2)') "12"
      close(10)

      open(10, file='fable_tmp_quuqu3ee')
      read(10, fmt='(a1)', iostat=ios) buf1
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, '"', buf1, '"'
      close(10)

      open(10, file='fable_tmp_quuqu3ee')
      read(10, fmt='(a2)', iostat=ios) buf2
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, '"', buf2, '"'
      close(10)

      open(10, file='fable_tmp_quuqu3ee')
      read(10, fmt='(a3)', iostat=ios) buf3
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, '"', buf3, '"'
      close(10)

      open(10, file='fable_tmp_quuqu3ee')
      read(10, fmt='(a3)', iostat=ios) buf33
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, '"', buf33, '"'
      close(10)

      open(10, file='fable_tmp_quuqu3ee')
      read(10, fmt='(a3)', iostat=ios) (buf4(i), i=1,2)
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0
      close(10)

      open(10, file='fable_tmp_quuqu3ee')
      read(10, fmt='(i3)', iostat=ios) k
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, k
      close(10)

      open(10, file='fable_tmp_quuqu3ee')
      read(10, fmt='(i3,i3)', iostat=ios) k, l
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, k, l
      close(10)

      open(10, file='fable_tmp_quuqu3ee')
      read(10, fmt='(i3)', iostat=ios) k, l
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0
      close(10)

      open(10, file='fable_tmp_quuqu3ee')
      read(10, fmt='(i3)', iostat=ios) (m(i), i=1,2)
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0
      close(10)

      open(10, file='fable_tmp_quuqu3ee')
      read(10, fmt='(d6.3)', iostat=ios) x
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0, x
      close(10)

      open(10, file='fable_tmp_quuqu3ee', status='unknown')
      write(10, '(a2)') "0x"
      close(10)

      open(10, file='fable_tmp_quuqu3ee')
      read(10, fmt='(i3)', iostat=ios) k
      write(*,*) ios .lt. 0, ios .eq. 0, ios .gt. 0
      close(10)

      end
