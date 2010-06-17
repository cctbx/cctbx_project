C Moving sub2 before sub1 changes the result:
C   Result with sub1 code first:
C          1           2         203
C         11          22         233
C        263
C         21          42        1263
C   Result with sub2 code first:
C          1           2           3
C         11          22          33
C         63
C         21          42        1063

      subroutine sub1
      save num2
      common /cmn/ num3
      data num1, num2, num3 /1, 2, 3/
      write(6, *) num1, num2, num3
      num1 = num1 + 10
      num2 = num2 + 20
      num3 = num3 + 30
      end

      subroutine sub2
      common /cmn/ num3
      data num3 /203/
      write(6, *) num3
      num3 = num3 + 1000
      end

      program prog
      call sub1
      call sub1
      call sub2
      call sub1
      end
