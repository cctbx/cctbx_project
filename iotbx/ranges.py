def range_parser( txt ):
  splitter = " "
  if "," in txt:
    splitter = ","
  txt_ranges = txt.split(splitter)
  ranges = []
  for range in txt_ranges:
    tmp = range.split("-")
    assert len(tmp)<=2
    if len(tmp)==2:
      ranges.append( [int(tmp[0]),int(tmp[1])]  )
    if len( tmp) == 1:
      ranges.append( [ int(tmp[0]) ] )
  return ranges

def range_to_list(range):
  result = []
  for item in range:

    if len(item)==1:
      result.append( item[0] )
    if len(item)==2:
      for ii in xrange(item[0],item[1]+1):
        result.append( ii )
  return result


def tst_ranges():
  range_txt="1-9,10,11,13-16"
  result = range_to_list( range_parser( range_txt ) )
  tst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16]
  for ii, jj in zip(result,tst):
    assert(ii==jj)

  range_txt="1-9 10 11 13-16"
  result = range_to_list( range_parser( range_txt ) )
  tst = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16]
  for ii, jj in zip(result,tst):
    assert(ii==jj)

if __name__ == "__main__":
  tst_ranges()
  print "OK"
