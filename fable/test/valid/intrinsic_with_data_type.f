      program prog
      intrinsic nint
      character*4 nint ! apparently ignored
      write(6, *) nint(4.5)
      end
