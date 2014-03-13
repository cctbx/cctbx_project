from __future__ import division
rdl_database = {
  "TRP" : {
    "default" : { # values from restraints library
        ("CA", "CB", "CG" ) : (113.6, 1.9),
        ("CB", "CG", "CD1") : (126.9, 1.5),
        },
    "m95" : {
        ("CA", "CB", "CG" ) : (115.6, 1.9),
        ("CB", "CG", "CD1") : (124.9, 1.5),
        },
    "t90" : {
        ("CA", "CB", "CG" ) : (115.6, 1.9),
        ("CB", "CG", "CD1") : (124.9, 1.5),
        },
      },
  "MET" : {
    "default" : { # values from restaints library
      ("CG", "SD", "CE") : (100.9, 2.2),
        },
    "default" : {
      ("CG", "SD", "CE") : (102.9, 2.2),
        },
      },
  "GLN" : {
    "default" : { # values from restaints library
      ("CA", "CB", "CG") : (114.1, 2.0),
        },
    "pp0?" : {
      ("CA", "CB", "CG") : (116.040, 2.0),
        },
      },
  }
