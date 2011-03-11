from stdlib import math as smath

expected_bond_lengths_by_element_pair = {
  # based on elbow/chemistry/BondLengths.py rev. 42
  # max of averages, rounded to one decimal
('H', 'H'): 0.0,
('AL', 'F'): 1.8,
('AS', 'C'): 2.0,
('AS', 'O'): 1.7,
('B', 'C'): 1.6,
('B', 'O'): 1.5,
('BR', 'C'): 2.0,
('C', 'C'): 1.5,
('C', 'CL'): 1.8,
('C', 'F'): 1.3,
('C', 'H'): 1.1,
('C', 'HG'): 2.3,
('C', 'N'): 1.4,
('C', 'O'): 1.4,
('C', 'P'): 1.7,
('C', 'S'): 1.7,
('C', 'SE'): 1.9,
('CO', 'N'): 2.0,
('CU', 'N'): 2.1,
('CU', 'O'): 1.8,
('F', 'O'): 1.8,
('FE', 'FE'): 2.2,
('FE', 'N'): 2.0,
('FE', 'O'): 2.0,
('FE', 'S'): 2.2,
('H', 'N'): 1.0,
('H', 'O'): 1.0,
('H', 'S'): 1.0,
('HG', 'O'): 2.3,
('MG', 'N'): 2.0,
('MG', 'O'): 2.2,
('N', 'N'): 1.3,
('N', 'NI'): 2.1,
('N', 'O'): 1.4,
('N', 'P'): 1.6,
('N', 'RU'): 1.8,
('N', 'S'): 1.6,
('O', 'O'): 1.4,
('O', 'P'): 1.6,
('O', 'S'): 1.5,
('O', 'U'): 1.8,
('O', 'V'): 2.0,
('O', 'W'): 2.0,
('P', 'S'): 1.7,
('S', 'S'): 2.0}


def build_edge_list(sites_cart, elements,slop=0.2):
  result = []
  for ii in range(len(sites_cart)):
    x1,y1,z1 = sites_cart[ii]
    for jj in range(ii+1,len(sites_cart)):
      x2,y2,z2 = sites_cart[jj]
      x2 = x2-x1
      y2 = y2-y1
      z2 = z2-z1
      dd = smath.sqrt( x2*x2+y2*y2+z2*z2 )
      expected_dist =  expected_bond_lengths_by_element_pair.get( (elements[ii], elements[jj]), False )
      if not expected_dist:
        expected_dist = expected_bond_lengths_by_element_pair.get( (elements[jj], elements[ii]), False )
      if not expected_dist:
        expected_dist = 1.7

      if dd <= expected_dist+slop:
        result.append( (ii,jj) )
      print dd, expected_dist

  return result


def tst_build_edge_list():
  sites_cart = [ (0,0,0), (0,0,1.53), (0,0,3.06) ]
  elements  =  [ 'C', 'O', 'N' ]
  tmp = build_edge_list(sites_cart, elements)
  assert (0,1) in tmp
  assert (1,2) in tmp
  assert (0,2) not in tmp

if __name__ == "__main__":
  tst_build_edge_list()
  print "OK"
