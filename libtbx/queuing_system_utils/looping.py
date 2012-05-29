def isodd(n):
  return n & 1 and True or False

def both_same_odd_or_even(n,m):
  if isodd(n)==isodd(m):
    return True
  return False

def generate_balanaced_list(n, only=None):
  for i in range(n):
    if only is not None and only!=i: continue
    for j in range(n):
      if j>i:
        if both_same_odd_or_even(i+1,j+1):
          pass
        else:
          continue
      elif j<i:
        if not both_same_odd_or_even(i+1, j+1):
          pass
        else:
          continue
      else:
        continue
      yield i,j
    if only==i: break

def generate_balanaced_list_from_list(l, only=None):
  if only is not None:
    assert only in l
    for i, item in enumerate(l):
      if only==item: break
    only = i
  for i,j in generate_balanaced_list(len(l), only=only):
    yield l[i], l[j]

if __name__=="__main__":
  print 'test generate_balanaced_list'
  pairs = {}
  for i,j in generate_balanaced_list(7):
    print i,j
    pairs.setdefault(i, [])
    assert j not in pairs[i]
    pairs[i].append(j)
    pairs.setdefault(j, [])
    assert i not in pairs[j]
    pairs[j].append(i)
  print pairs
  for key in pairs:
    assert len(pairs[key])==6
    assert key not in pairs[key]
  print "OK"
  for k, (i,j) in enumerate(generate_balanaced_list(7, only=0)):
    print k,i,j
  assert k==2
  print 'test generate_balanaced_list_from_list'
  l=list("abcdef")
  for i,j in generate_balanaced_list_from_list(l):
    print i,j
  print '-'*10
  for i,j in generate_balanaced_list_from_list(l, only=l[1]):
    print i,j
  print "OK"
