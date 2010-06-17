      program prog
      dimension ns2(2), ns3(3), ns4(4), ns5(5)
      equivalence (ns2(2), ns3(1))
      equivalence (ns5(1), ns4(4))
      save
      do i=1,2
        ns2(i) = 20+i
      enddo
      do i=1,3
        ns3(i) = 30+i
      enddo
      do i=1,4
        ns4(i) = 40+i
      enddo
      do i=1,5
        ns5(i) = 50+i
      enddo
      write(6, *) ns2
      write(6, *) ns3
      write(6, *) ns4
      write(6, *) ns5
      end
