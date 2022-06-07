from __future__ import absolute_import, division, print_function
#
# bayesian_estimator.py
# tct 053008, updated 2014-10-29
#
# Creates a Bayesian estimator from a list of records containing y and {x_i}
# Then you can use it to estimate values of y from a new  {x_i}
#
# For typical use just copy and edit exercise_2 below
#
# For use with data that may have missing entries for some predictors, try
#  exercise_group below
#
# Parameters to set:  number of bins, minimum observations in each bin,
#  number of bins to smooth over.   Normally you want 20 or more observations
#  in each bin when you are training the estimator. You also want enough bins
#  so that adjacent bins have similar values compared to the SD of the
#  predictor variables (i.e., the range of values of x_1 in bin 3 should overlap
#  the range of values of x_1 in bin 4)  (no big jumps).
#  If you have too few observations to satisify these criteria, then you can
#  smooth the bins (it rarely hurts in any event).
#
# purpose: take a list of N measurements x_ij of M types x_i, and
#  N perfect values p_j, #   and generate an estimator that will estimate
#  y from one set of measurements of the M types.
#
# Method: assume that for any given value of y, the values of the
# M measurements will have Gaussian distributions around their means u_i(y)
# NOTE: Current version assumes flat priors.
#
# Here is how we can calculate p({x_i},y) any time we want for
# any {x_i} and any y:

#  From our assumption:
#       p({x_i},y)=p({x_i}|y)p(y)   # definition
#  If our distributions are Gaussian then we can represent the joint
#    probability of all the measurements (given y) of p({x_i}|y) as:
#    p({x_i}|y) = exp{-0.5(x-u)S-1(x-u)T } / (2pi det(S)
#  where x is {x_i}, u is {u_i} (the vector of means) and S is the covariance
#  matrix. We can calculate u and S from the members of our list of
#  measurements that led to the value y.
#
# Now given a set of measurements {x_i} we use Bayes' rule to estimate y:
#
# p(y)= p({x_i},y)/ integral_over_y[ p({x_i},y ]
#
# Implementation details:

# Need a set of bins n_bins_i for each variable x_i  and for y as we have
# no model for  how x_i and y are related (we want this to be general).

# If a bin has too few data, just pull data from nearby bins until there are
#   at least min_in_bin data.

# May be necessary to smooth the data in the bins if it is sparse.
# Average the means vector and the covariance matrices over smooth_bins if so


import os,sys,math,random
from scitbx import matrix
from cctbx.array_family import flex
from copy import deepcopy
from libtbx.utils import Sorry,null_out
from six.moves import zip
from six.moves import range

pi=2.*math.atan2(1.,0.)

class bayesian_estimator:
  def __init__(self,out=None,verbose=False,skip_covariance=True,
     minimum_records=20,trim_predictor_values=True):
    if out is None:
      self.out=sys.stdout
    else:
      self.out=out
    self.verbose=verbose
    self.trim_predictor_values=trim_predictor_values # limit to range of training set
    self.skip_covariance=skip_covariance
    self.minimum_records=minimum_records
    self.ok=True

  def create_estimator(self,record_list=None,
    n_bin=None,range_low=None,range_high=None,
    min_in_bin=20,use_flat_y_dist=False,
    smooth_bins=10):

    if not n_bin:
       raise Sorry(
     "Please set values for n_bin,eg., n_bin=100, min_in_bin=5, smooth_bins=5")
    self.record_list=record_list# vectors with measurements. First one is y
    self.n_bin=n_bin # bins of y
    self.smooth_bins=smooth_bins # smoothing the bins
    self.use_flat_y_dist=use_flat_y_dist # allow user to ignore dist of y
    self.min_in_bin=min_in_bin # minimum records in a bin
    self.range_low=range_low  # range_low of y
    self.range_high=range_high  # range_high of y
    self.record_list_range_low=None
    self.record_list_range_high=None
    # get set up:
    ok=self.check_estimator_inputs() # check and print out what we are doing..
    if not ok:
      self.ok=False
      return  # give up

    self.record_list.sort()  # sort on y (first variable in each record)
                             # so it is easy to bin
    # Split by bins of y
    self.get_range()  # decide on range of y if user didn't set it.
    self.split_into_bins() # split y and measurements x_i into self.n_bin bins
    self.merge_bins_if_nec() # if a bin has fewer than min_in_bin...take its
                             # neighbors too..
    # ready to create estimator
    self.create_estimator_lists()
    self.setup_p_y()

    # Now we can calculate p({x_i}|y) based on our means and covariance:
    #  prob=self.get_p_xx(values,bin)

    # and to get  y_bar and SD using p(y|{x_i}) we just sum up
    #  y_bar,sd=self.get_y_bar_sd_given_x({x_i})

  def apply_estimator(self,record_list_or_values,single_value=False):
    if not self.ok: return None,None

    if single_value:
      return self.get_y_bar_sd_given_x(record_list_or_values)

    predicted_data=flex.double()
    predicted_sd_data=flex.double()
    for values in record_list_or_values:
      y_bar,sd=self.get_y_bar_sd_given_x(values)
      predicted_data.append(y_bar)
      predicted_sd_data.append(sd)
    return predicted_data,predicted_sd_data

  def get_y_bar_sd_given_x(self,x_list):
    # return weighted value of y and SD of this
    p_y=flex.double()
    for bin in self.bins:
      p_y.append(self.get_p_xx(x_list,bin)*self.get_p_y(bin))
    sum=flex.sum(p_y)
    if sum<=0.0: sum=1.
    p_y=p_y/sum  # normalized probability
    weighted_bin=flex.sum(p_y*self.bin_mid_double)
    weighted_sq_bin=flex.sum(p_y*self.bin_mid_double*self.bin_mid_double)
    sd=weighted_sq_bin-weighted_bin**2
    if sd<0.: sd=0.
    sd=math.sqrt(sd)
    return weighted_bin,sd

  def get_p_xx(self,x_list,bin):
    if self.trim_predictor_values:
      x_vector=flex.double()
      for lower,higher,value in zip(self.record_list_range_low,
         self.record_list_range_high,x_list):
        x_vector.append(max(lower,min(higher,value)))
    else:
      x_vector=flex.double(x_list)

    means_vector,cov_matrix_inv,determinant=self.get_mean_cov_inv_determinant(bin)
    if means_vector is None or cov_matrix_inv is None or determinant is None:
      return 0.0 #
    delta_vect=x_vector-means_vector
    delta_row_matrix=matrix.rec(delta_vect,[len(x_list),1])
    delta_col_matrix=delta_row_matrix.transpose()
    inner=delta_col_matrix*cov_matrix_inv*delta_row_matrix
    v=inner.as_list_of_lists()[0][0]
    if v < -64.: v=-64.
    elif v> 64.: v=64.
    prob=math.exp(-0.5*v)/(2.*pi*determinant)
    return prob

  def get_p_y(self,bin):
    return self.bin_p_y[bin]

  def setup_p_y(self):
    sum=0.
    for bin in self.bins:
      if self.use_flat_y_dist:
        sum=sum+1.
      else:
        sum=sum+self.bin_orig_number_of_records[bin]
    self.bin_p_y=[]
    if sum==0: sum=1.
    for bin in self.bins:
      if self.use_flat_y_dist:
        self.bin_p_y.append(1./sum)
      else:
        self.bin_p_y.append(self.bin_orig_number_of_records[bin]/sum)

  def create_estimator_lists(self):
    self.bin_means_list=[]
    self.bin_cov_list=[]

    self.bin_cov_inv_list=[]
    self.bin_determinant_list=[]
    for bin in self.bins:
       records=self.get_records_in_bin(bin,
           self.bin_start_list,self.bin_after_end_list)
       if self.verbose:
         print("Working on bin ",bin," with ",len(records)," records", file=self.out)
       # now get vector of means and covariance matrix for these records.
       means_vector,cov_matrix=self.get_mean_cov_from_records(records)
       self.bin_means_list.append(means_vector)
       self.bin_cov_list.append(cov_matrix)

    # now smooth values of cov_matrix and means_list if desired
    self.smooth_cov_and_means()
    for bin in self.bins:
       cov_matrix=self.bin_cov_list[bin]
       try:
         self.bin_cov_inv_list.append(cov_matrix.inverse())
         self.bin_determinant_list.append(cov_matrix.determinant())
       except Exception as e:  # sets result to zero later if we do this...
         self.bin_cov_inv_list.append(None)
         self.bin_determinant_list.append(None)

  def smooth_cov_and_means(self): # smooth with window of smooth_bins
                                  # only smooth symmetrically, so ends not done
    if self.smooth_bins<2: return
    new_bin_cov_list=self.bin_cov_list
    new_bin_means_list=self.bin_means_list
    for bin in self.bins:
      local_window=self.smooth_bins//2
      if bin < local_window: local_window=bin
      elif self.n_bin-1-bin < local_window: local_window=self.n_bin-1-bin
      new_bin_cov_list[bin]=self.smooth_list(self.bin_cov_list,bin,local_window,
        is_matrix=True)
      new_bin_means_list[bin]=\
        self.smooth_list(self.bin_means_list,bin,local_window)
      if self.verbose:
        print("Bin ",bin," smoothed with window of ",local_window, file=self.out)
    self.bin_cov_list=new_bin_cov_list
    self.bin_means_list=new_bin_means_list

  def smooth_list(self,bin_value_list,bin,local_window,is_matrix=False):
    # smooth each entry in value_list[bin] with local window.
    new_entry_list=[]
    for i in range(len(bin_value_list[bin])): # one entry at a time
      sum=0.
      sumn=0.
      for local_bin in range(max(0,bin-local_window),
        min(bin+local_window+1,len(bin_value_list))): # XXX necessary?
        sumn+=1.
        if is_matrix:
          sum+=bin_value_list[local_bin].elems[i]
        else:
          sum+=bin_value_list[local_bin][i]
      if sumn==0.: sumn=1.
      new_entry_list.append(sum/sumn)
    if is_matrix:
      new_entry_list=matrix.rec(flex.double(new_entry_list),[self.n,self.n])
    else:
      new_entry_list=flex.double(new_entry_list)

    if self.verbose:
      print('smoothing bin ',bin,' value=',bin_value_list[bin],\
       'new_value=',new_entry_list, file=self.out)
    return new_entry_list



  def get_mean_cov_from_records(self,records):
    # get the vector of means and covariance matr
    # 070308: If self.skip_covariance=True then set covariance terms to zero.
    # rewrite our data into a list of flex.double() vectors
    if records is None or len(records)<2: return None,None

    list_of_vectors=[]
    for i in range(len(records[0])-1):
     list_of_vectors.append(flex.double())
    for record in records:
      values=record[1:]  # skip the y-value at the start
      for value,vector in zip(values,list_of_vectors):
        vector.append(value)

    list_of_means=[]
    list_of_indices=[]
    i=0
    for x in list_of_vectors:
      i+=1
      list_of_means.append(flex.mean(x))
      list_of_indices.append(i)
    if self.verbose:
      print("Means: ",list_of_means, file=self.out)
    m=[]
    for i,x,mean_x in zip(list_of_indices,list_of_vectors,list_of_means):
      for j,y,mean_y in zip(list_of_indices,list_of_vectors,list_of_means):
        covar=flex.mean((x-mean_x)*(y-mean_y))
        if self.verbose:
          print("covar: ",i,j,covar, file=self.out)
        if self.skip_covariance and i!=j:
          covar=0.0
          if self.verbose:
            print("NOTE: Covariance set to zero as skip_covariance=True", file=self.out)
        m.append(covar)
    mx=matrix.rec(m,2*[len(list_of_vectors)])
    return flex.double(list_of_means),mx

  def get_mean_cov_inv_determinant(self,bin):
    if bin<0 or bin>self.n_bin-1: return None,None
    return self.bin_means_list[bin],self.bin_cov_inv_list[bin], \
      self.bin_determinant_list[bin]

  def check_estimator_inputs(self):
    assert self.record_list is not None
    assert len(self.record_list)>0
    self.n=len(self.record_list[0])-1   # number of indep variables
    print("Creating Bayesian estimator from ",self.n," measurement vectors", file=self.out)
    print("and one vector of perfect values", file=self.out)
    assert self.n>0 # need at least one meas and one obs
    print("Total of ",len(self.record_list)," values of each", file=self.out)
    if len(self.record_list)<self.minimum_records:
      return False
    for v in self.record_list:
       assert len(v)==len(self.record_list[0])
       assert type(v)==type([1,2,3])

    if self.n_bin is None:
      self.n_bin=len(self.record_list)//50
      if self.n_bin < 2: self.n_bin=2
    print("Number of bins: ",self.n_bin, file=self.out)
    assert self.n_bin>0
    return True

  def get_range(self):
    if self.range_low is None or self.range_high is None:
      self.range_low=self.record_list[0][0]
      self.range_high=self.record_list[-1][0]
    print("Range of y: ",self.range_low,self.range_high, file=self.out)
    self.full_range=self.range_high-self.range_low
    assert self.full_range > 0.
    self.delta_range=self.full_range/self.n_bin
    print("Delta range: ",self.delta_range, file=self.out)
    self.bins=range(self.n_bin)

    # 2014-11-13 also get range for x-values
    n=len(self.record_list[0])-1
    self.record_list_range_low=n*[None]
    self.record_list_range_high=n*[None]
    for record in self.record_list:
      for i in range(n):
        ii=i+1
        if self.record_list_range_low[i] is None or \
             record[ii]<self.record_list_range_low[i]:
             self.record_list_range_low[i]=record[ii]
        if self.record_list_range_high[i] is None or \
             record[ii]>self.record_list_range_high[i]:
             self.record_list_range_high[i]=record[ii]
    print("Range for predictor variables:", file=self.out)
    for i in range(n):
      print("%d:  %7.2f - %7.2f " %(
          i,self.record_list_range_low[i],self.record_list_range_high[i]), file=self.out)

  def split_into_bins(self):
    self.bin_start_list=self.n_bin*[None]  # first record in this bin
    self.bin_after_end_list=self.n_bin*[None] # (last record in this bin)+1
    i_record=-1
    for record in self.record_list:
      i_record+=1
      i_bin=self.get_bin(record[0])  # get bin based on y
      if self.bin_start_list[i_bin] is None:
         self.bin_start_list[i_bin]=i_record
      if self.bin_after_end_list[i_bin] is None or \
             self.bin_after_end_list[i_bin] < i_record+1:
         self.bin_after_end_list[i_bin]=i_record+1

    # save original number of records in each bin
    self.bin_orig_number_of_records=[]
    self.bin_mid_list=[]
    for bin in self.bins:
       records=self.get_records_in_bin(bin,
           self.bin_start_list,self.bin_after_end_list)
       self.bin_orig_number_of_records.append(len(records))
       self.bin_mid_list.append(self.get_mid(bin))
    self.bin_mid_double=flex.double(self.bin_mid_list)

    if self.verbose:
      for bin,n in zip(self.bins,self.bin_orig_number_of_records):
        print("Records in bin ",bin,": ",n, file=self.out)

  def merge_bins_if_nec(self):
    # if a bin has fewer than min_in_bin...take neighbors
    self.new_bin_start_list=deepcopy(self.bin_start_list)
    self.new_bin_after_end_list=deepcopy(self.bin_after_end_list)
    for bin in self.bins:
       n=len(self.get_records_in_bin(bin,
           self.bin_start_list,self.bin_after_end_list))
       if n>=self.min_in_bin: continue

       offset=0
       while n < self.min_in_bin and offset < self.n_bin:
         offset+=1
         low=bin-offset
         high=bin+offset
         if low<0: low=0
         if high>self.n_bin-1: high=self.n_bin-1
         new_start=self.bin_start_list[low]
         new_after_end=self.bin_after_end_list[high]
         if new_start is not None:
            self.new_bin_start_list[bin]=new_start
         if new_after_end is not None:
            self.new_bin_after_end_list[bin]=new_after_end

         n=len(self.get_records_in_bin(bin,
           self.new_bin_start_list,self.new_bin_after_end_list))
       if self.verbose:
         print("Reset bin ",bin,\
            " to contain records ",self.new_bin_start_list[bin],\
            ' to ',self.new_bin_after_end_list[bin], file=self.out)
       if self.new_bin_start_list[bin] is None or  \
          self.new_bin_after_end_list[bin] is None:
          raise AssertionError("Sorry, need more data or smaller value of "+\
             "min_in_bin (currently "+str(self.min_in_bin)+")")
       if self.verbose:
         print("New number of records: ", \
         -self.new_bin_start_list[bin] +self.new_bin_after_end_list[bin], file=self.out)
    self.bin_start_list=self.new_bin_start_list
    self.bin_after_end_list=self.new_bin_after_end_list


  def get_records_in_bin(self,i_bin,bin_start_list,bin_after_end_list):
    # return the records in bin i_bin
    if i_bin<0 or i_bin>self.n_bin : return []
    start=bin_start_list[i_bin]
    after_end=bin_after_end_list[i_bin]
    if start is None or after_end is None: return []
    return self.record_list[start:after_end]


  def get_mid(self,i_bin):  # midpoint in y of bin i_bin
    return self.range_low+self.delta_range*float(i_bin)

  def get_bin(self,value):
    a_bin= float(self.n_bin)*(value - self.range_low)/self.full_range
    i_bin=int(a_bin)
    # bin 0 goes from range_low to range_low+(range_high-range_low)/self.n_bin
    if i_bin>self.n_bin-1: i_bin=self.n_bin-1
    if i_bin<0: i_bin=0
    return i_bin

  def generate_data(self,sig_list,n=1000,missing=True,non_linear=False):
    perfect_data=flex.double()
    data_record_list=[]
    record_list=[]
    for i in range(n):
      u=float(i)*1000./float(n)
      record=[]
      pp=u/100.
      if non_linear:
        r=pp/10
        r=r**2+.5*r**3+2.*r**4
      else:
        r=pp
      if missing and u>199. and u<275.: continue
      if missing and u>0 and u<50.:continue
      record.append(pp)
      perfect_data.append(pp)
      s1=random.gauss(0.,sig_list[0])
      for s in sig_list[1:]:
        ss=random.gauss(0.,s)
        record.append(r+s1+ss)
      record_list.append(record)
      data_record_list.append(record[1:])
    print("Set up ",len(record_list),\
         " y values and measurements of a,b ", file=self.out)
    return perfect_data,data_record_list,record_list

  def get_average_data(self,records):
    avg=flex.double()
    for values in records:
      v=flex.double(values)
      avg.append(flex.mean(v))
    return avg

  def exercise(self,iseed=39771):

    print("\nTESTING runs of Bayesian estimator with and without missing ", file=self.out)
    print("data and with and without use of covariance.\n", file=self.out)
    # set up our arrays
    # our model is a and b are related to p, with correlated error s1 and
    #   independent errors s2 and s3
    # a=p+s1+s2
    # b=p+s1+s3
    import random
    random.seed(iseed)
    result_list=[]
    for missing,non_linear in zip([True,False],[False,True]):
     for skip_covar in [False,True]:
      print("Missing data: ",missing,"   Non-linear: ",non_linear," Skip-covariance: ",skip_covar, file=self.out)
      record_list=[]
      sig_list=[.50,0.2,0.8]
      perfect_data,records,record_list=\
         self.generate_data(sig_list,n=1000,missing=missing,
         non_linear=non_linear)
      self.skip_covariance=skip_covar
      self.create_estimator(record_list=record_list,n_bin=100,
         use_flat_y_dist=True,
         min_in_bin=5,
         smooth_bins=0)

      more_perfect_data,more_records,more_record_list=\
         self.generate_data(sig_list,n=200,missing=missing,
         non_linear=non_linear)
      predicted_data,predicted_sd_data=\
         self.apply_estimator(more_records)

      c=flex.linear_correlation(more_perfect_data,predicted_data)
      cc=c.coefficient()

      average_data=self.get_average_data(more_records)
      c=flex.linear_correlation(more_perfect_data,average_data)
      cc_avg=c.coefficient()
      print(" CC: ",cc," just averaging: ",cc_avg, file=self.out)
      print(file=self.out)
      result_list.append(cc_avg)
      if 0:
        f=open('test'+str(non_linear)+'.out','w')
        for perf,pred,sd,values in zip(
          more_perfect_data,predicted_data,predicted_sd_data,more_records):
          print(perf,pred,sd, end=' ', file=f)
        for value in values: print(value, end=' ', file=f)
        print(file=f)
        f.close()
    return result_list


  def exercise_2(self,iseed=712771,out=sys.stdout):

   # You can just cut and paste this into a new python script and run it
   # Then edit it for your purpose.

   print("\nTESTING Bayesian estimator on sample randomized data\n", file=self.out)
   import math,random
   from mmtbx.scaling.bayesian_estimator import bayesian_estimator
   from cctbx.array_family import flex

   # create some training data. First value=target, rest are predictor variables
   record_list=[]
   random.seed(iseed)
   for i in range(1000):
    y=float(i)/100.
    x1=y+random.gauss(1.,.5)
    x2=y+random.gauss(-2.,.8)
    record_list.append([y,x1,x2])

   run=bayesian_estimator(out=out)
   run.create_estimator(
         record_list=record_list,
         n_bin=100,
         min_in_bin=5,
         smooth_bins=5)

   # Now create some observed data.
   # (just predictor variables, save target y for comparison)
   obs_record_list=[]
   target_y_list=[]  # just for comparison
   for i in range(200):
    y=float(i)/20.
    x1=y+random.gauss(1.,.5)
    x2=y+random.gauss(-2.,.8)
    obs_record_list.append([x1,x2]) # note no y here!
    target_y_list.append(y)  # just for comparison

   # get the predicted values:
   predicted_data,predicted_sd_data=\
       run.apply_estimator(obs_record_list)


   if 0:
     f=open("out.dat","w")
     for observed_record,predicted_value,y in zip (
         obs_record_list,predicted_data,target_y_list):
        print(observed_record,predicted_value,y, file=f)
     print("Data, predicted y, actual y written to out.dat", file=out)
     f.close()

   cc=flex.linear_correlation(flex.double(target_y_list),flex.double(predicted_data)).coefficient()
   print("CC: ",cc, file=self.out)
   return cc
#

def get_table_as_list(lines=None,text="",file_name=None,record_list=None,
   info_list=None,info_items=None,
   data_items=None,target_variable=None,
   start_column_header=None,
   select_only_complete=False, minus_one_as_none=True,
   d_min_as_data=False,out=sys.stdout):

  if not record_list:
     record_list=[]  # data records: [target,predictor1,predictor2...]
     info_list=[]    # info records [key, resolution]
  if not data_items:
     data_items=[]
     info_items=[]
  skip_columns=[]
  if not record_list:
    n=None
    first=True
    n_info=0
    if file_name:
      with open(file_name) as f:
        text=f.read()
    if not lines:
      lines=text.splitlines()
    for line in lines:
      if line.startswith("#"): continue # skip comments
      if not line.replace(" ",""): continue # skip empty lines
      if line.lstrip().startswith("Total"): continue # clean up file format
      spl=line.split()
      if first and not spl[0].lstrip().startswith('0'):
        print("Input data file columns:%s" %(line), file=out)
        # figure out if some of the columns are "key" and "d_min"
        all_info_items=[]
        all_items=line.split()
        if start_column_header:
           new_all_items=[]
           started=False
           i=0
           for item in all_items:
             if not started and item==start_column_header:
               started=True
             if started:
                new_all_items.append(item)
             else:
                skip_columns.append(i)
             i+=1
           all_items=new_all_items
           print("Skipped columns: %s" %(str(skip_columns)), file=out)
           print("Columns to read: %s" %(str(all_items)), file=out)
           if not started:
             raise Sorry("Sorry the header %s "% (start_column_header) +
               "was not found in the line %s" %(line))


        if all_items[0].lower()=='key':
          all_info_items.append(all_items[0])
          all_items=all_items[1:]
        if not d_min_as_data and \
          all_items[0].lower() in  ['d_min','dmin','res','resolution']:
          all_info_items.append('d_min') # standardize
          all_items=all_items[1:]
        target=all_items[0]
        all_items=all_items[1:]
        n_info=len(all_info_items)
        if not data_items:
          data_items=all_items
          info_items=all_info_items
        if not target_variable:
          target_variable=target
        first=False
        continue
      if skip_columns:
        spl=spl[len(skip_columns):]
      xx=[]
      info=spl[:n_info]
      for i in range(1,n_info):
        info[i]=float(info[i])
      for x in spl[n_info:]:
        if x.lower() in ["-1","-1.000","None","none"]:
          xx.append(None)
        else:
          xx.append(float(x))
      if n is None: n=len(xx)
      assert len(xx)==n
      if xx[0] is not None:
        record_list.append(xx)
        info_list.append(info)

  if select_only_complete:
    new_list=[]
    new_info_list=[]
    for x,y in zip(record_list,info_list):
      if not None in x:
        new_list.append(x)
        new_info_list.append(y)
    record_list=new_list
    info_list=new_info_list
  return record_list,info_list,target_variable,data_items,info_items

class estimator_group:
  # holder for a group of estimators based on the same data but using different
  # subsets for prediction.  Useful in case the request for estimation is
  # missing some of the predictor variables.

  # Also useful for using only data to a certain resolution cutoff in the
  # prediction. This is turned on if resolution_cutoffs is supplied and a
  #   resolution is supplied for apply_estimator()
  # if a resolution is requested that is outside the range of resolution_cutoffs
  #  then a predictor based on the nearest resolution is used.

  def __init__(self,
        n_bin=100,
        min_in_bin=5,
        smooth_bins=5,
        resolution_cutoffs=None,
        skip_covariance=True,
        minimum_records=20,
        verbose=False,
        out=sys.stdout):

    self.estimator_dict={}
    self.keys=[]
    self.combinations={} # keyed on resolution
    self.variable_names=[]
    self.info_names=[]
    self.n_bin=n_bin
    self.min_in_bin=min_in_bin
    self.smooth_bins=smooth_bins
    self.out=out
    self.verbose=verbose
    if self.verbose:
      self.verbose_out=out
    else:
      self.verbose_out=null_out()
    self.training_records=[]
    self.training_file_name='None'
    self.skip_covariance=skip_covariance
    self.minimum_records=minimum_records

    if resolution_cutoffs:
       resolution_cutoffs.sort()
       resolution_cutoffs.reverse() # highest first
    else:
       resolution_cutoffs=[0]
    self.resolution_cutoffs=resolution_cutoffs

  def show_summary(self):
    # show summary of this estimator
    print("\nSummary of Bayesian estimator:", file=self.out)
    print("\nVariable names:", end=' ', file=self.out)
    for x in self.variable_names: print(x, end=' ', file=self.out)
    print(file=self.out)
    print("Resolution cutoffs:", end=' ', file=self.out)
    for resolution_cutoff in self.resolution_cutoffs:
       print("%5.2f" %(resolution_cutoff), end=' ', file=self.out)
    print(file=self.out)

    if self.combinations:
      print("\nCombination groups:", file=self.verbose_out)
      for resolution_cutoff in self.resolution_cutoffs:
        print("Resolution cutoff: %5.2f A" %(resolution_cutoff), end=' ', file=self.verbose_out)
        for x in self.combinations[resolution_cutoff]:
          print("(", end=' ', file=self.verbose_out)
          for id in x:
            print("%s" %(self.variable_names[id]), end=' ', file=self.verbose_out)
          print(")  ", end=' ', file=self.verbose_out)
        print(file=self.verbose_out)

    print("Bins: %d   Smoothing bins: %d  Minimum in bins: %d " %(
        self.n_bin,self.smooth_bins,self.min_in_bin), file=self.verbose_out)

    print("Data records used in training estimator: %d\n" %(
      len(self.training_records)), file=self.out)
    print("Data file source of training data: %s\n" %(
      self.training_file_name), file=self.out)

  def add_estimator(self,estimator,resolution_cutoff=None,combination=None):
    key=self.get_key(combination=combination,
        resolution_cutoff=resolution_cutoff)
    self.estimator_dict[key]=estimator
    self.keys.append(key)
    if not resolution_cutoff in self.combinations:
      self.combinations[resolution_cutoff]=[]
    self.combinations[resolution_cutoff].append(combination)

  def set_up_estimators(self,
      record_list=None,
      data_items=None,
      file_name=None,
      text=None,
      select_only_complete=False,
      minimum_complete=True):

    if not record_list:
      record_list,info_list,target_variable,data_items,info_items=\
       self.get_record_list(
        file_name=file_name,
        text=text,
        select_only_complete=select_only_complete)

    n=len(record_list[0])
    assert data_items is not None


    self.training_records=record_list

    assert data_items
    self.variable_names=data_items
    self.info_names=info_items
    print("Total of %d records in data file" %(len(record_list)), file=self.out)
    if len(record_list)<1: return
    n=len(record_list[0])


    # now split in all possible ways and get estimators for each combination
    # if self.resolution_cutoffs is set, do it with each resolution cutoff
    # 2014-12-18 if minimum_complete=True, require that there are
    #  self.minimum_records complete ones
    for resolution_cutoff in self.resolution_cutoffs:
      print("\nResolution cutoff of %5.2f A" %(resolution_cutoff), file=self.verbose_out)
      local_record_list,local_info_list=self.select_records(
        record_list=record_list,
        info_list=info_list,resolution_cutoff=resolution_cutoff,
        minimum_complete=minimum_complete)
      if len(local_record_list)<1:
        print("No records selected for this resolution cutoff", file=self.verbose_out)
        self.delete_resolution_cutoff(resolution_cutoff)
        continue

      print("%d records selected for this resolution cutoff" %(
       len(local_record_list)), file=self.verbose_out)
      n_predictors=n-1
      print("\nPredictor variables: %d" %(n_predictors), file=self.verbose_out)
      for pv in data_items:
        print("%s" %(pv), end=' ', file=self.verbose_out)
      print(file=self.verbose_out)
      if len(data_items)!=n_predictors:
        raise Sorry("Number of predictor variables (%d) must be the same as " %(
            len(data_items)) +
          "the number \nof predictor variables in data records(%d)" %(
            n_predictors))

      combinations=self.get_combinations(list(range(n_predictors)))
      for combination in combinations:
       estimator=self.get_estimator(
          record_list=local_record_list,columns=combination)
       if estimator is not None:
         key=self.get_key(combination=combination,
           resolution_cutoff=resolution_cutoff)
         self.add_estimator(estimator,resolution_cutoff=resolution_cutoff,
           combination=combination)
    for resolution_cutoff in self.resolution_cutoffs: # make sure there is an estimator
      if not resolution_cutoff in self.combinations:
        self.delete_resolution_cutoff(resolution_cutoff)

  def get_record_list(self,
     file_name=None,
     text=None,
     select_only_complete=None):
      record_list=[]
      data_items=[]
      if text:
        print("Setting up Bayesian estimator using supplied data", file=self.out)
        lines=text.splitlines()
      else:
        if not file_name:
          import libtbx.load_env
          file_name=libtbx.env.find_in_repositories(
            relative_path=os.path.join("mmtbx","scaling","cc_ano_data.dat"),
            test=os.path.isfile)
        if not file_name:
          raise Sorry("Need file with training data")
        print("\nSetting up Bayesian estimator using data in %s\n" %(file_name), file=self.out)
        with open(file_name) as f:
          lines=f.readlines()
        self.training_file_name=file_name
      return get_table_as_list(
         lines=lines,
         select_only_complete=select_only_complete,out=self.out)

  def delete_resolution_cutoff(self,resolution_cutoff):
    new_cutoffs=[]
    for rc in self.resolution_cutoffs:
      if not rc==resolution_cutoff:
        new_cutoffs.append(rc)
    self.resolution_cutoffs=new_cutoffs

  def select_records(self,
      record_list=None,
      info_list=None,
      next_cutoff=None,
      minimum_complete=None,
      resolution_cutoff=None,
      first=True):
    # select those records that are >= resolution_cutoff and not >= the
    # next-highest cutoff in the list (if any)
    if not self.info_names or not resolution_cutoff or \
         not self.resolution_cutoffs or len(self.resolution_cutoffs)==1:
      return record_list,info_list

    if self.info_names[-1]!='d_min' and len(self.resolution_cutoffs)>1:
      raise Sorry("Cannot select records on resolution unless the database"+
        "(%s) has resolution information" %(str(self.training_file_name)))
    if first and next_cutoff is None:
      for rc in self.resolution_cutoffs:
        if (next_cutoff is None or rc <next_cutoff) and rc > resolution_cutoff:
          next_cutoff=rc

    new_records=[]
    new_info=[]
    count=0
    for record,info in zip(record_list,info_list):
      if info[-1]>=resolution_cutoff and \
        (next_cutoff is None or info[-1]<next_cutoff):
         new_records.append(record)
         new_info.append(info)
         if minimum_complete:
            if not (None in record):
              count+=1
    if next_cutoff is not None and next_cutoff<info_list[-1][-1] and (
       len(new_records)<self.minimum_records or
       minimum_complete and count < self.minimum_records):
      # try again with bigger range
      next_next_cutoff=None
      for rc in self.resolution_cutoffs:
        if (next_next_cutoff is None or rc <next_next_cutoff) \
          and rc > next_cutoff:
          next_next_cutoff=rc
      new_records,new_info=self.select_records(
          record_list=record_list,info_list=info_list,
          next_cutoff=next_next_cutoff,resolution_cutoff=resolution_cutoff,
          minimum_complete=minimum_complete,first=False)
    return new_records,new_info

  def get_key(self,combination=None,resolution_cutoff=None):
    if resolution_cutoff is None: resolution_cutoff=0
    return str(combination)+"_"+str(resolution_cutoff)

  def variable_names(self):
    return self.variable_names

  def get_estimator(self,record_list=[],
      columns=[]):

    from mmtbx.scaling.bayesian_estimator import bayesian_estimator
    from cctbx.array_family import flex

    from copy import deepcopy
    record_list=deepcopy(record_list)

    column_chooser=[]
    n=max(columns)
    for x in range(n+1):
      if x in columns:
        column_chooser.append(True)
      else:
        column_chooser.append(False)

    selected_list=[]
    for xx in record_list:
      new_x=[]
      have_all=True
      for x,keep in zip(xx,[True]+column_chooser):
         if keep:
           new_x.append(x)
           if x is None:
             have_all=False
      if have_all:
        selected_list.append(new_x)

    if not selected_list:
      return None

    estimator=bayesian_estimator(out=null_out(),
         skip_covariance=self.skip_covariance)
    estimator.create_estimator(
          record_list=selected_list,
         n_bin=self.n_bin,
         min_in_bin=self.min_in_bin,
         smooth_bins=self.smooth_bins)
    if not estimator.ok: return None

    return estimator

  def get_combinations(self,all_together):
    if len(all_together)<1: return []
    elif len(all_together)==1: return [all_together]
    else:  # return all sub_combinations
      combinations=[all_together]
      for n in range(len(all_together)):
        for x in self.get_combinations(all_together[:n]+all_together[n+1:]):
          if not x in combinations:
            combinations.append(x)
      return combinations

  def apply_estimators(self,value_list=None,data_items=None,
    resolution=None):
    assert len(value_list)==len(self.variable_names)
    if data_items != self.variable_names:
      print("WARNING: data items do not match working variables:", file=self.out)
      print("Item: variable name", file=self.out)
      for data_item,variable in zip(data_items,self.variable_names):
         print("%s: %s " %(data_item,variable), file=self.out)
      raise Sorry("Data items do not match working variables. Use verbose to see more info.")
    if len(self.resolution_cutoffs)>1 and resolution is None:
      raise Sorry("Must supply resolution if resolution_cutoffs is set")

    # choose which resolution cutoff to use
    rc=self.resolution_cutoffs[-1] # take highest-res  no matching resolution
    for resolution_cutoff in self.resolution_cutoffs:
      if resolution and resolution >= resolution_cutoff:
        rc=resolution_cutoff
        break

    # figure out which data are present and select the appropriate estimator
    columns=[]
    values=[]
    for x,n in zip(value_list,range(len(value_list))):
      if x is not None:
        columns.append(n)
        values.append(x)
    key=self.get_key(combination=columns,resolution_cutoff=rc)
    if not key in self.keys:
      return None,None
    estimator=self.estimator_dict[key]

    return estimator.apply_estimator(values,single_value=True)

prediction_values="""cc_perfect cc      skew    e
0.361   None    0.0124  1.04
0.476   0.382   0.0455  0.536
0.672   0.752   0.1546  0.14
0.47    0.321   0.0743  0.636
0.317   0.145   0.0299  0.795
0.046   0.026   -0.0088 0.978
0.055   None    -0.0121 1.007
0.391   None    0.0252  0.749
0.407   None    0.025   0.711
0.635   0.381   0.1601  0.495
0.684   0.481   0.1846  0.453
0.684   0.481   0.1846  0.453
0.057   0.008   -0.0002 0.938
0.349   0.189   -0.0004 0.781
0.482   0.375   0.0405  0.607
0.031   -0.046  0.0002  1.027
0.195   0.035   0.0073  0.914
0.268   0.12    0.0052  0.856
0.417   None    0.0427  0.657
0.524   None    0.0279  0.576
0.527   None    0.0256  0.751
0.284   0.152   0.0144  0.923
0.425   0.376   0.0803  0.655
0.13    -0.051  -0.0042 0.984
0.152   -0.029  0.0056  0.991
0.217   0.045   0.0117  0.969
0.177   0.029   -0.0084 0.979
0.266   0.088   0.0075  0.905
0.492   0.191   0.0105  0.664
0.309   0.021   -0.0021 0.874
0.62    0.354   0.1432  0.659
0.784   0.608   0.192   0.515
0.67    0.499   0.147   0.588
0.082   -0.059  0.0053  0.952
0.11    0.023   0.0163  0.96
0.174   0.061   0.0028  0.914
0.108   0.023   0.0072  0.989
0.572   None    0.0384  0.351
0.404   0.16    0.0306  0.818
0.625   0.464   0.1383  0.616
0.306   0.099   0.0056  0.791
0.21    0.147   0.0348  0.745
0.213   0.123   0.0195  0.76
0.391   0.303   0.0248  0.812
0.214   0.062   -0.0023 0.961
0.053   0.014   0.0154  1.016
0.066   None    -0.0295 0.919
0.34    None    0.0443  0.943
0.404   None    0.0448  0.856
0.371   0.275   0.0263  0.688
0.06    -0.08   0.005   0.997
0.426   0.22    0.047   0.589
0.178   0.034   0.0031  0.973
0.22    0.071   0.012   0.949
0.045   -0.02   0.0011  0.844
0.318   0.088   0.0155  0.965
0.026   -0.033  -0.0089 0.922
0.25    0.044   -0.0053 0.822
0.161   0.014   0.0194  0.998
0.352   0.168   0.0434  0.801
0.296   0.115   0.0156  0.817
0.231   0.1     0.0079  0.803
0.192   0.03    -0.0028 0.94
0.041   -0.022  -0.0091 1.005
0.236   0.057   0.0161  0.987
0.521   0.315   0.017   0.539
0.633   0.549   0.0266  0.314
0.175   0.089   0.0219  0.843
0.251   0.139   0.0024  0.972
0.175   0.068   0.0132  1.009
0.141   0.038   0.0224  0.93
0.342   0.141   0.0378  0.799
0.288   0.074   -0.0065 0.94
0.049   -0.074  -0.0049 1.013
0.326   0.166   0.0252  0.862
0.28    None    0.0173  0.952
0.045   None    -0.019  0.896
0.26    0.088   0.0022  0.835
0.103   -0.007  0.031   0.926
0.139   0.042   0.0038  0.87
0.302   0.226   0.0464  0.803
0.236   0.067   0.0255  0.987
0.255   0.11    0.0031  0.925
0.071   -0.026  -0.0097 0.996
0.508   0.352   0.0556  0.676
0.472   0.262   0.0668  0.695
0.427   0.172   0.0241  0.744
0.206   0.039   0.0038  0.913
0.194   0.039   -0.0038 0.92
0.006   -0.009  -0.0089 0.981
0.049   0.001   -0.005  0.993
0.257   0.139   0.0064  0.853
0.368   0.336   0.0599  0.727
0.29    0.186   0.0179  0.834
0.175   0.033   0.0212  0.867
0.326   0.107   0.0176  0.844
0.41    0.192   0.0192  0.762
0.341   0.111   0.0106  0.822
0.486   0.229   0.0675  0.783
0.497   0.262   0.08    0.802
0.385   0.134   0.0202  0.768
0.605   0.498   0.064   0.468
0.41    0.202   0.0293  0.757
0.506   0.273   0.0651  0.63
0.276   0.078   0.0391  0.97
0.286   0.098   0.0099  0.888
0.372   0.127   0.0349  0.795
0.311   0.087   0.0406  0.859
0.353   0.189   0.0057  0.709
0.35    0.184   0.0121  0.714
0.278   0.119   0.0122  0.829
0.459   0.151   0.0412  0.693
0.508   0.199   0.0454  0.658
0.593   0.294   0.0599  0.568
0.473   0.185   0.0471  0.789
0.604   0.352   0.0732  0.588
0.658   0.419   0.093   0.507
0.274   0.057   0.0103  0.849
0.282   0.071   0.0125  0.812
0.449   0.319   0.0284  0.647
0.528   0.455   0.0389  0.51
0.434   0.28    0.0264  0.648
0.351   0.19    0.0151  0.76
0.375   0.22    0.0135  0.726
0.298   0.161   0.0146  0.824
0.394   0.278   0.0231  0.702
0.29    0.147   0.0083  0.853
0.324   0.104   0.0141  0.952
0.371   0.15    0.032   1
0.463   0.201   0.028   0.791
0.461   0.292   0.0129  0.644
0.202   0.03    -0.0009 0.947
0.259   0.066   -0.0132 0.903
0.225   0.037   0.0121  0.924
0.323   0.197   0.0218  0.72
0.344   0.287   0.0503  0.619
0.302   0.204   0.0149  0.751
0.266   0.059   0.0206  0.892
0.252   0.059   0.0153  0.907
0.319   0.12    0.0252  0.833
0.273   0.064   0.0197  0.874
0.13    -0.008  0.0236  0.972
0.378   0.129   0.0451  0.804
0.478   0.301   0.0579  0.635
0.285   0.093   0.0091  0.857
0.317   0.133   0.032   0.8
0.263   0.177   0.0086  0.818
0.114   0.037   0.0031  0.931
0.31    0.081   0.0121  0.843
0.432   0.19    0.0482  0.716
0.303   0.075   0.0212  0.861
0.337   0.089   0.0264  0.847
0.322   0.086   0.0074  0.839
0.381   0.136   0.0425  0.773
0.474   0.211   0.0415  0.675
0.317   0.109   0.017   0.807
0.286   0.083   0.0024  0.833
0.39    0.175   0.0534  0.741
0.473   0.251   0.0671  0.587
0.185   0.058   -0.0086 0.951
0.272   0.138   -0.0133 0.843
0.282   0.145   0.0079  0.882
0.269   0.111   0.0149  0.862
0.406   0.171   0.0179  0.709
0.53    0.345   0.0473  0.531
0.217   0.147   0.0109  0.92
0.315   0.272   0.0307  0.819
0.295   0.189   0.0065  0.852
0.331   0.249   0.0061  0.796
0.196   0.017   0.0052  0.968
0.198   0.018   0.0122  0.982
0.502   None    0.0185  0.49
0.286   None    0.0138  0.993
0.247   None    -0.0123 1.033
0.341   0.158   0.0175  0.88
0.253   0.134   0.0064  0.903
0.368   0.209   0.0115  0.78
0.479   0.282   0.0399  0.658
0.618   0.474   0.1005  0.437
0.642   0.506   0.0807  0.427
0.34    0.301   0.0188  0.82
0.331   0.263   0.0254  0.846
0.374   0.354   0.0414  0.767
0.08    None    0.0179  0.801
0.071   None    0.0774  0.579
0.301   None    0.0168  0.721
0.317   None    0.02    0.845
0.045   None    -0.0017 0.972
0.272   None    0.0002  0.587
0.346   None    0.0127  0.307
0.24    None    0.0224  0.619
0.223   None    0.0286  0.611
0.183   None    -0.005  0.568
0.423   None    0.0686  0.441
0.487   None    0.0429  0.452
0.375   None    0.052   0.817
0.424   None    0.0113  0.5
0.368   None    0.0121  0.637
0.319   None    0.0026  0.594
0.317   None    0.0276  0.486
0.39    None    0.0287  0.568
0.526   None    0.0595  0.401
0.274   None    0.0148  0.648
0.481   None    0.014   0.337
0.261   None    0.0028  0.835
0.324   None    0.0513  0.659
0.347   None    0.0107  0.709
0.304   None    0.0254  0.603
0.385   None    0.0326  0.295
0.037   None    0.0397  0.991
0.305   None    0.0259  0.755
0.448   None    0.0457  0.699
0.247   None    -0.0111 0.739
0.105   None    0.0012  0.661
0.491   None    0.0492  0.486
0.338   None    -0.0031 0.648
0.385   None    0.0491  0.809
0.218   None    -0.002  0.679
0.544   None    0.0488  0.477
0.345   None    0.0295  0.643"""


pred_data="""cc_perfect    cc skew      e
0.023   0.0163  0.96
None    0.0028  0.914
0.023   0.0072  None
None    None 0.351
0.16    None  None
0.464   0.1383  0.616
0.099   0.0056  0.791
0.147   None 0.745
0.123   0.0195  0.76
None  None  None
0.303   0.0248  0.812"""

target_values=[0.476,0.672,0.47,0.317,0.046,0.635,0.684,0.684,0.057]

def exercise_group(out=sys.stdout):

  print("\nTESTING exercise_group(): group of predictors with same ", file=out)
  print("data but some missing entries.\n", file=out)

  import libtbx.load_env
  file_name=libtbx.env.find_in_repositories(
      relative_path=os.path.join("mmtbx","scaling","cc_ano_data.dat"),
      test=os.path.isfile)

  estimators=estimator_group(resolution_cutoffs=[0.,1.,2.,3.,4.,5.],
     out=out)
  estimators.set_up_estimators(file_name=file_name)
  estimators.show_summary()

  print("Running estimator now", file=out)
  prediction_values_as_list,info_values_as_list,\
     dummy_target_variable,dummy_data_items,dummy_info_items=\
    get_table_as_list(file_name=file_name, select_only_complete=False,out=out)

  # run through all prediction_values data
  from cctbx.array_family import flex
  target_list=flex.double()
  result_list=flex.double()
  for target_and_value,info in zip(
      prediction_values_as_list,info_values_as_list):
    resolution=info[-1]
    y,sd=estimators.apply_estimators(value_list=target_and_value[1:],
      data_items=estimators.variable_names,resolution=resolution)
    if y is None:
      raise Sorry("Estimator failed")
    else:
      target_list.append(target_and_value[0])
      result_list.append(y)
  cc1=flex.linear_correlation(target_list,result_list).coefficient()
  print("Prediction CC: %7.3f " %(cc1), file=out)

  # and using another dataset included here
  print(file=out)
  prediction_values_as_list,info_values_as_list,target_variable,\
     data_items,info_items=get_table_as_list(
     text=prediction_values,select_only_complete=False,out=out)

  estimators=estimator_group(out=out)
  estimators.set_up_estimators(
    text=prediction_values,select_only_complete=True,
    data_items=data_items)
  estimators.show_summary()

  all_value_list,all_info_list,target_variable,data_items,info_items=\
      get_table_as_list(text=pred_data,out=out)

  target_list=flex.double()
  result_list=flex.double()

  for value_list,target in zip(all_value_list,target_values):
    y,sd=estimators.apply_estimators(
     value_list=value_list,data_items=data_items)
    if y is not None:
      print("Y SD Target diff : %7.3f  %7.3f   %7.3f  %7.3f" %(
        y,sd,target,y-target), file=out)
      target_list.append(target)
      result_list.append(y)
    else:
      print("No data:  %7.3f " %(target), file=out)
  cc2=flex.linear_correlation(target_list,result_list).coefficient()
  print("Correlation: %6.3f" %(cc2), file=out)
  return cc1,cc2


# cross-validation of estimates of cc* from skew, cc_half, esqr
# 2014-10-28 tt
#

def get_suitable_target_and_values(target_and_values_list,test_values):
  suitable_test_values=[]
  for test_value in test_values:
    if test_value is not None:
      suitable_test_values.append(test_value)

  suitable_target_and_values_list=[]
  for target_and_values in target_and_values_list:
    new_target_and_values=[target_and_values[0]]
    for test_value,value in zip(test_values,target_and_values[1:]):
      if test_value is not None:
         new_target_and_values.append(value)
    if not None in new_target_and_values:
      suitable_target_and_values_list.append(new_target_and_values)
  return suitable_test_values,suitable_target_and_values_list

def jacknife(target_and_values_list,i,
    skip_covariance=True,
    skip_jacknife=False,out=sys.stdout):
  # take out entry i and get predictor and then predict value for entry i

  # if test entry does not have values for all the predictor variables then
  #  use only the ones that it has
  values=target_and_values_list[i][1:] # everything except target
  if skip_jacknife: # include everything
    tv_list=deepcopy(target_and_values_list)
  else:  # usual
    tv_list=deepcopy(target_and_values_list[:i])+\
            deepcopy(target_and_values_list[i+1:])
  # select only work entries that have all the necessary values
  suitable_test_values,suitable_target_and_values_list=\
     get_suitable_target_and_values(tv_list,values)
  if not suitable_test_values: return None,None

  target=target_and_values_list[i][0]
  estimator=bayesian_estimator(
    skip_covariance=skip_covariance,
    minimum_records=1,
    out=null_out())
  estimator.create_estimator(
        record_list=suitable_target_and_values_list,
        n_bin=100,
        min_in_bin=5,
        smooth_bins=5,
  )
  est,sd=estimator.apply_estimator(suitable_test_values,single_value=True)
  return target,est

def run_jacknife(args=None,no_jump=True,
     file_name=None,record_list=None,info_list=None,out=sys.stdout):

  if args is None: args=[]
  verbose=('verbose' in args)
  if 'testing' in args:
    print("\nTESTING jacknife run of bayesian estimator.\n", file=out)
  skip_covariance=(not 'include_covariance' in args)

  if file_name is not None or record_list is not None:
    # already have it
    pass
  elif args and os.path.isfile(args[0]):
    file_name=args[0]
    print("Setting up jacknife Bayesian predictor with data from %s" %(
      file_name), file=out)
  else:
    import libtbx.load_env
    file_name=libtbx.env.find_in_repositories(
      relative_path=os.path.join("mmtbx","scaling","cc_ano_data.dat"),
      test=os.path.isfile)
    print("Setting up jacknife Bayesian predictor with data from %s" %(
      file_name), file=out)
  if not record_list and (not file_name or not os.path.isfile(file_name)):
    raise Sorry("Unable to find the file '%s' " %(str(file_name)))

  if record_list is None:
    record_list,info_list,target_variable,data_items,info_items=\
       get_table_as_list(
         file_name=file_name,
         select_only_complete=False,
         out=out)
    print("Size of data list: %d" %(len(record_list)), file=out)
  if info_list is None:
    info_list=len(record_list)*["None"]

  # Set up estimator using all but one entry and predict it from the others

  if no_jump:
    n_jump=1
  elif no_jump is None and 'no_jump' in args:
    n_jump=1
  else:
    n_jump=20
    print("Testing every %d'th entry in jacknife" %(n_jump), file=out)
  from cctbx.array_family import flex
  target_list=flex.double()
  result_list=flex.double()
  for i in range(len(record_list)):
    if n_jump*(i//n_jump)!=i: continue
    target,value=jacknife(record_list,i,skip_covariance=skip_covariance,
        out=out)
    if target is None or value is None:
       raise Sorry("estimator failed to return a result")
    target_list.append(target)
    result_list.append(value)
    info=info_list[i]
    if verbose:
      for x in info:
        print(x, end=' ', file=out)
      print("%7.3f  %7.3f " %(target,value), file=out)
  cc=flex.linear_correlation(
    target_list,result_list).coefficient()
  print("CC: %7.3f " %(cc), file=out)
  return cc,target_list,result_list

if __name__=="__main__":
  args=sys.argv[1:]
  if not 'testing' in args: args.append('testing')
  if 'run_jacknife' in args:
    run_jacknife(args=args,no_jump=None)
  elif 'exercise_group' in args:
    exercise_group()
  else:
    run_jacknife(args=args,no_jump=None)
    exercise_group()
    run=bayesian_estimator()
    run.exercise()
    run.exercise_2()
