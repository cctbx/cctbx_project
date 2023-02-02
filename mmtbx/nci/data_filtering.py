from __future__ import division
from scitbx.array_family import flex
from libtbx import group_args
import matplotlib.pyplot as plt
from libtbx import easy_pickle

def filter(fn, mode):
  fo = open(fn,"r")
  sh = flex.double()
  kh = flex.double()
  for l in fo.readlines():
    key, resolution, rama_fav, rama_out, Z, c_beta, rot_out, clashes,\
    bonds, angles, t1o_skew, t1o_kurtosis, dhao_skew, dhao_kurtosis, Zh,Zs,Zl, rw,rf = l.split()
    #
    resolution    = float(resolution   )
    rama_fav      = float(rama_fav     )
    rama_out      = float(rama_out     )
    c_beta        = float(c_beta       )
    rot_out       = float(rot_out      )
    clashes       = float(clashes      )
    bonds         = float(bonds        )
    angles        = float(angles       )
    t1o_skew      = float(t1o_skew     )
    t1o_kurtosis  = float(t1o_kurtosis )
    dhao_skew     = float(dhao_skew    )
    dhao_kurtosis = float(dhao_kurtosis)
    #
    if fn.count("all")>0   and Z  == "None": continue
    if fn.count("alpha")>0 and Zh == "None": continue
    if fn.count("beta")>0  and Zs == "None": continue
    if resolution is None: continue
    #
    if(mode=="t1o"):
      S = t1o_skew
      K = t1o_kurtosis
    elif(mode=="dhao"):
      S = dhao_skew
      K = dhao_kurtosis
    else: assert 0
    #
    ### HIGH RES BEST MODELS
    f0 = (fn.count("all")>0   and float(Z)  > -2.)
    f1 = (fn.count("alpha")>0 and float(Zh) > -2.)
    f2 = (fn.count("beta")>0  and float(Zs) > -2.)
    good = rama_out<1 and rama_fav>95 and \
      c_beta<0.1 and rot_out<2 and \
      clashes<10 and bonds<0.03 and angles<3 and (f0 or f1 or f2) and\
      rw != "None" and float(rw)<0.25 and resolution<1.5
    if(good):
      sh.append(S)
      kh.append(K)
  fo.close()
  return group_args(sh=sh, kh=kh)

def run():
  """
  This script generates filtered_raw_data.pkl from raw text files named all, alpha, beta.
  Every line of these files should contain comma-separated values of:
  pdb_id, resolution, rama_fav, rama_out, Z, c_beta, rot_out, clashes,
    bonds, angles, t1o_skew, t1o_kurtosis, dhao_skew, dhao_kurtosis, Zh,Zs,Zl, rw,rf

  Files themselves are not in the repository because of their size (~80Mb)
  """
  fig = plt.figure(figsize=(10,15))
  i=1
  filtered_results = {"all" : {}, "alpha": {}, "beta":{} }
  for ss in ["all", "alpha", "beta"]:
    for iseq, mode in enumerate(["t1o", "dhao"]):
      ra = filter(ss, mode=mode)
      filtered_results[ss][mode] = (ra.sh, ra.kh)
      # Plotting figure to make sure the filtered sets are OK
  #     ax = plt.subplot(int("42%d"%i))
  #     ax.scatter(ra.sh,   ra.kh,   s=3,   c='black',   edgecolors='none')
  #     i+=1
  #     plt.subplots_adjust(wspace=0.1, hspace=0.125)
  # fig.savefig("for_contours.pdf")

  # Writing results into file.
  # easy_pickle.dump(obj=filtered_results, file_name="filtered_raw_data.pkl")

if(__name__ == "__main__"):
  run()
