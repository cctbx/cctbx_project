      program prog
      character*3 hols(4)
      data hols /1hX,2hYz,3HPqR,4hSTuv/
   10 format(1Ha,2hcD,3heFg,4hHijK,17hLMnO
     &PqRstUvWxyZ@#,3h   ,1h$)
      write(6, 10)
      write(6, '(4(''['',a,'']''))') hols
      call show(1hx, 1)
      call show(2hUs, 2)
      call show(3hPdW, 3)
      call show(19hrTiTGBrDYtATTSwDkSw, 19)
      end

      subroutine show(hol, l)
      character hol*(l)
      write(6, '(3a)') '>', hol(1:l), '<'
      end
