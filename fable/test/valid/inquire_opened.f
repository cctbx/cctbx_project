      program prog
C mkdir files
C touch files/direct
C ln files/direct files/hard_link
C ln -s files/direct files/soft_link
C ln -s files files_link
C ifort: T F F F F
C gfortran: T T T T T
      logical lvar
      open(1, file='files/direct')
      inquire(file='files/direct', opened=lvar)
      write(6, *) lvar
      inquire(file='./files/direct', opened=lvar)
      write(6, *) lvar
      inquire(file='files/hard_link', opened=lvar)
      write(6, *) lvar
      inquire(file='files/soft_link', opened=lvar)
      write(6, *) lvar
      inquire(file='files_link/direct', opened=lvar)
      write(6, *) lvar
      end
