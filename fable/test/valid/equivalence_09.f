      program prog
      common /scr/ nc
      save ns
      equivalence (nc, nce), (ns, nse), (nl, nle)
      nc = 56
      ns = 34
      nl = 12
      write(6, *) nc, ns, nl
      write(6, *) nse, nle, nce
      end
