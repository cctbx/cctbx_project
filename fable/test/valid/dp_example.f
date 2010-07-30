      program dp_example
      integer array_size
      parameter(array_size=10)
      dimension array(array_size)
      write(6, *) 'I can store up to', array_size, 'numbers'
      do i=1,array_size
        array(i) = i
      enddo
      end
