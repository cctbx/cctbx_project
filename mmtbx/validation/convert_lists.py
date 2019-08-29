from __future__ import absolute_import, division, print_function
import os, sys
import csv

def run(filename,
        output,
        list_type="csv",
        ):
  print('filename',filename)
  print('list_type',list_type)
  if list_type=="csv":
    with open(filename, 'rb') as csvfile:
      data = {}
      data.setdefault("data", [])
      spamreader = csv.reader(csvfile) #, delimiter=' ', quotechar='|')
      for i, row in enumerate(spamreader):
        if i==0:
          for i, item in enumerate(row):
            row[i] = row[i].strip()
            row[i] = row[i].replace(" ", "_")
            row[i] = row[i].replace("%", "percent")
            row[i] = row[i].replace("#", "number")
            row[i] = row[i].lower()
          data["headers"] = row
        else:
          data["data"].append(row)
      def _cmp_pdb_id(l1, l2):
        if l1[1]<l2[1]: return -1
        return 1
      data["data"].sort(_cmp_pdb_id)
      def _guess_type(item):
        try:
          item=int(item)
          return int
        except ValueError : pass
        try:
          item=float(item)
          return float
        except ValueError: pass
        return str
      funcs = {}
      for row in data["data"]:
        for i, item in enumerate(row):
          funcs.setdefault(i, {})
          t = _guess_type(item)
          funcs[i].setdefault(t, 0)
          funcs[i][t]+=1
      for i in funcs:
        for k in funcs[i].keys():
          funcs[i][funcs[i][k]] = k
          del funcs[i][k]
      for row in data["data"]:
        for i, item in enumerate(row):
          if type(item)==type(""):
            if item.lower() in ["yes"]: row[i]=True
            elif item.lower() in ["no"]: row[i]=False
            elif item.lower() in ["none"]: row[i]=None
      for i in funcs:
        func = max(funcs[i].keys())
        func = funcs[i][func]
        for row in data["data"]:
          for j, item in enumerate(row):
            if i!=j: continue
            if type(item)!=type(""): continue
            try: row[j]=func(item)
            except ValueError as e:
              row[j]=float(item)
            if type(row[j])==type(""): row[j]=row[j].strip()
      #
      # filter
      #
      remove = []
      for i, header in enumerate(data["headers"]):
        if header in ["cluster_id"]:
          remove.append(i)
      remove.reverse()
      #
      # output
      #
      outl =  "data = {"
      outl += '\n  "headers" : ['
      for i, header in enumerate(data["headers"]):
        if i in remove: continue
        outl += '\n    "%s",' % header
      outl += '\n    ],'
      outl += '\n  "data" : ['
      for row in data["data"]:
        for i in remove:
          del row[i]
        outl += '\n    %s,' % row
      outl += '\n    ],'
      outl += "\n  }"
      outl += '''

if __name__=="__main__":
  print data["headers"]
  print
  for item in data["data"]:
    print item
  print
'''
      print(outl)
      f=open("%s.py" % output, "w")
      f.write(outl)
      f.close()

      os.system("python %s.py" % output)


if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
