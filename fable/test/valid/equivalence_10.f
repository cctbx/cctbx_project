      program prog
      common /scr/ nc(2)
      dimension nce(2)
      dimension ns(2), nse(2)
      save ns
      dimension nl(2), nle(2)
      equivalence (nc, nce), (ns, nse), (nl, nle)
      nc(1) = 56
      nc(2) = -56
      ns(1) = 34
      ns(2) = -34
      nl(1) = 12
      nl(2) = -12
      write(6, *) nc, ns, nl
      write(6, *) nse, nle, nce
      end
