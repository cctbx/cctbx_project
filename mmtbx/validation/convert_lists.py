import os, sys
import csv

def run(filename, list_type="csv"):
  print 'filename',filename
  print 'list_type',list_type
  if list_type=="csv":
    with open(filename, 'rb') as csvfile:
      data = {}
      data.setdefault("data", [])
      spamreader = csv.reader(csvfile) #, delimiter=' ', quotechar='|')
      for i, row in enumerate(spamreader):
        if i==0:
          headers = row
          data["headers"] = headers
        else:
          data["data"].append(row)
      def _cmp_pdb_id(l1, l2):
        if l1[1]<l2[1]: return -1
        return 1
      data["data"].sort(_cmp_pdb_id)
      for row in data["data"]:
        for i, item in enumerate(row):
          if i in [3,4,5,6,7,8]:
            row[i]=float(row[i])
          elif i in [9]:
            row[i]=int(row[i])
      outl =  "data = {"
      outl += '\n  "headers" : ['
      for header in data["headers"]:
        outl += '\n    "%s",' % header
      outl += '\n    ],'
      outl += '\n  "data" : ['
      for row in data["data"]:
        outl += '\n    %s,' % row
      outl += '\n    ],'
      outl += "\n  }"
      print outl
      f=file("top8000.py", "wb")
      f.write(outl)
      f.close()


if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
