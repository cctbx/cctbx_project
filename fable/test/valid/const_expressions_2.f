      program prog
      parameter(n1=2)
      parameter(n2=n1*4)
      parameter(n3=(n2+n1)*5)
      parameter(n4=(n3-n2)**3)
      parameter(n5f=n2*(n4+1)/1.99)
      write(6, *) n5f
      end
