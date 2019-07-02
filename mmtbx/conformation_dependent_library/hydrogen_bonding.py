from __future__ import absolute_import, division, print_function

hb_database = { # donor acceptor
  ("F", "F") : 38.6, #(161.5 kJ/mol or 38.6 kcal/mol)
  ("O", "N") :  6.9,
  ("O", "O") :  5.0,
  ("N", "N") :  3.1,
  ("N", "O") :  1.9,
  #("HOH", "OH+") : 4.3,
  }

def run():
  print(hb_database)

if __name__=="__main__":
  run()#sys.argv[1])
