from __future__ import absolute_import, division, print_function
import pandas as pd

def run(reference_file, test_file):
  refdata = pd.read_pickle(reference_file)
  cosdata = pd.read_pickle(test_file)

  #inner join
  merged_inner = pd.merge(left=refdata, right=cosdata, left_on='experiment', right_on='experiment')

  same = []
  isame = 0
  for idx in range(len(merged_inner["experiment"])):
    is_same = merged_inner['coset_x'][idx] == merged_inner['coset_y'][idx]
    if is_same: isame +=1
    same.append(is_same)
  merged_inner['same'] = same

  pd.set_option('display.max_rows', 300)
  # What's the size of the output data?
  print(merged_inner.shape)
  print(merged_inner)

  print("File %s, %d aligned experiments"%(reference_file, len(refdata["experiment"])))
  print("File %s, %d aligned experiments"%(test_file, len(cosdata["experiment"])))
  print("The common experiment set consists of %d"%len(merged_inner["experiment"]))
  A = len(merged_inner["same"])
  B = len(merged_inner[merged_inner["same"]==True])
  print("%d of %d common alignments agree on coset assignment, %0.2f%%"%(B,A,100.*B/A))

  all_expt = len(set(list(refdata["experiment"])).union(set(list(cosdata["experiment"]))))
  print("%d matches out of a total of %d experiments, or %.2f%%"%(B, all_expt, 100.*B/all_expt))
  return A, B/A

if __name__=="__main__":
  import sys
  reference_file = sys.argv[1]
  test_file = sys.argv[2]
  run(reference_file, test_file)
