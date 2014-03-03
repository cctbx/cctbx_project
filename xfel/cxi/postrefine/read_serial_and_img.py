from __future__ import division
#read and extract serial no. (of the pickle file) and the image name.

import sys

file_a = open(sys.argv[1],'r')
data_a=file_a.read().split("\n")
search_word = 'DISTL'
line_no_now = 0
n_zeros = 5
for line in data_a:
  if line.find(search_word) > 0:
    data = data_a[line_no_now-1].split(' ')
    serial = data[0]
    img_file = data[1]
    print 'int-data_'+serial.zfill(n_zeros)+'.pickle', img_file
  line_no_now += 1
