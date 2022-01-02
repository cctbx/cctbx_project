from __future__ import absolute_import, division, print_function
import math,sys
from six.moves import range

# breakpoint.py
# tool to analyze a plot where Y decreases with X and is best fit by 2 lines.
# Purpose of the tool is to identify the breakpoint in X where the slope changes
# Optionally convert Y to log10(Y)
# Optionally identify value of X where the first line has value lower than the second
#   line by an amount given by offset (or log10(offset))

input_target_offset=1
input_min_ratio=2.
input_use_log=True
input_value_dict={
0:6.49,
16.79:6.11,
33.57:5.92,
50.36:5.71,
67.15:5.15,
83.93:3.87,
100.72:3.35,
117.51:3,
134.29:2.78,
151.08:2.64,
167.87:2.52,
184.65:2.43,
}

def get_low_high(value_dict=None):
  x_low=None
  x_high=None
  for x in value_dict.keys():
    if x_low is None or value_dict[x] < value_dict[x_low]: x_low=x
    if x_high is None or value_dict[x] > value_dict[x_high]: x_high=x
  return x_low,x_high

def get_log(value_dict=None):
  for x in value_dict.keys():
    value_dict[x]=math.log10(value_dict[x])
  return value_dict

def get_closest_x(value=None,value_dict=None):
  keys = list(value_dict.keys())
  best_dist=None
  best_x=None
  for x in keys:
    dist=abs(x-value)
    if best_dist is None or dist< best_dist:
      best_x=x
      best_dist = dist
  return best_x

def get_y_calc(x=None,x1=None,y1=None,x2=None,y2=None):
  if x2==x1:
    return y1
  else:
    return y1+(x-x1)*(y2-y1)/(x2-x1)

def get_y_calc_from_multiple_lines(x=None,
     value_dict=None,x_midway=None,y_midway=None,min_x=None,max_x=None):
    if x < x_midway:
      y_calc=get_y_calc(x=x,
        x1=min_x,y1=value_dict[min_x],
        x2=x_midway,y2=y_midway)
    else:
      y_calc=get_y_calc(x=x,
        x1=x_midway,y1=y_midway,
        x2=max_x,y2=value_dict[max_x])

    return y_calc

def get_delta_from_multiple_lines(x=None,y=None,
     value_dict=None,x_midway=None,y_midway=None,min_x=None,max_x=None):

  return y-get_y_calc_from_multiple_lines(x=x,
     value_dict=value_dict,x_midway=x_midway,y_midway=y_midway,min_x=min_x,max_x=max_x)

def score_breakpoint(value_dict=None,x_midway=None,y_midway=None,min_x=None,max_x=None):
  import math
  score=0.
  score_n=0.
  for x in value_dict.keys():
    delta=get_delta_from_multiple_lines(x=x,y=value_dict[x],
     value_dict=value_dict,x_midway=x_midway,y_midway=y_midway,min_x=min_x,max_x=max_x)
    score+=delta**2
    score_n+=1
  if score_n>0:
    score=math.sqrt(score/score_n)
  return score

def find_breakpoint(value_dict=None,decreasing=True,min_ratio=1,
   use_log=True,target_offset=1.,sufficient_ratio=200,
    min_value=1,out=sys.stdout):
  assert decreasing #  only set up for decreasing values
  x_low,x_high=get_low_high(value_dict=value_dict)
  if decreasing:
    min_x=x_high
    max_x=x_low
  print("Low,high:",x_low,x_high,":",value_dict[x_low],value_dict[x_high], file=out)

  if decreasing and x_low < x_high:
    print("Not decreasing", file=out)
    return None # not decreasing

  ratio=value_dict[x_high]/max(1,value_dict[x_low])
  if not use_log:
    if ratio >-50 and ratio < 50:
      ratio=10**ratio
    else:
      from libtbx.utils import Sorry
      raise Sorry("Failed to use exponential (value too large): %s" %(ratio))

  if min_ratio and ratio <min_ratio:
    print("Ratio too low (min %7.2f  found %7.2f)" %(min_ratio,ratio), file=out)
    return None  # it is not enough to tell

  from copy import deepcopy
  value_dict=deepcopy(value_dict)

  # remove all entries outside of x_high to x_low and make sure all >=1
  for x in list(value_dict.keys()):
    if value_dict[x]<min_value: value_dict[x]=min_value
    if x < min_x: del value_dict[x]
    if x > max_x: del value_dict[x]

  # find highest x with value > sufficient_ratio*value_dict[max_x]
  lowest_x_keep=None
  for x in value_dict.keys():
    rr=value_dict[x]/value_dict[max_x]
    if not use_log: rr=10**rr
    if rr > sufficient_ratio and (
       lowest_x_keep is None or x > lowest_x_keep):
      lowest_x_keep=x
  if lowest_x_keep is not None:
    print("Removing all values of x < %7.2f" %(lowest_x_keep), file=out)
    for x in list(value_dict.keys()):
      if x < lowest_x_keep:
        del value_dict[x]
    min_x=lowest_x_keep
    x_high=min_x


  if use_log:
    value_dict=get_log(value_dict=value_dict)
    target_offset=math.log10(target_offset)

  ratio=value_dict[x_high]/max(1,value_dict[x_low])
  y_midway=value_dict[x_low]*ratio**0.5
  print("Low,high:",x_low,x_high,":",value_dict[x_low],value_dict[x_high], file=out)
  x_midway=get_closest_x(value=y_midway,value_dict=value_dict)

  y_midway=value_dict[x_midway]
  print("Midway:",y_midway," x-midway:",x_midway,\
    score_breakpoint(value_dict=value_dict,x_midway=x_midway,y_midway=y_midway,min_x=min_x,max_x=max_x), file=out)

  # XXX removing during 2to3 migration, no clear function?
  # keys=list(value_dict.keys())
  # keys.sort()

  # make a little grid search for refine x_midway, y_midway to minimize distance to either of the two lines
  delta_x=0.1*(max_x-min_x)/len(value_dict)
  delta_y=0.1*score_breakpoint(
    value_dict=value_dict,x_midway=x_midway,y_midway=y_midway,min_x=min_x,max_x=max_x)
  nx=41
  ny=41
  best_score=None
  best_xm=None
  best_ym=None
  for i in range(nx):
    xm=x_midway+(i-nx/2)*delta_x
    for j in range(ny):
      ym=y_midway+(j-ny/2)*delta_y
      score=score_breakpoint(value_dict=value_dict,x_midway=xm,y_midway=ym,min_x=min_x,max_x=max_x)
      if score is None: continue
      if best_score is None or score<best_score:
         best_score=score
         best_xm=xm
         best_ym=ym
  if best_score is not None:
    x_midway=best_xm
    y_midway=best_ym
    print("\n Best fit: x_midway: %7.2f  y_midway: %7.2f " %(
        x_midway,y_midway), file=out)

  keys=list(value_dict.keys())
  keys.sort()
  print("\n   X       Y-obs     Y-fit    Delta", file=out)

  for x in keys:
     delta=get_delta_from_multiple_lines(x=x,y=value_dict[x],
        value_dict=value_dict,x_midway=x_midway,y_midway=y_midway,min_x=min_x,max_x=max_x)
     pred_y=value_dict[x]-delta
     print("%7.2f  %7.2f  %7.2f  %7.2f " %(x,value_dict[x],pred_y,delta), file=out)

  # And identify value of x where value extrapolated from x<x_midway becomes less than
  #   y_midway minus target_offset (typically 1=1 log unit)
  if y_midway == value_dict[min_x]:
    dxdy=0
  else:
    dxdy=(x_midway-min_x)/(y_midway-value_dict[min_x])
  delta_x=-1.*dxdy*target_offset
  target_x=x_midway+delta_x
  print("Target X-value at breakpoint: %7.2f " %(target_x), file=out)
  return target_x

def get_value_dict(file_name=None):
  dd={}
  for line in open(file_name).readlines():
    if line.strip().startswith("X"): continue
    spl=line.split()
    if len(spl)<2: continue
    dd[float(spl[0])]=float(spl[1])
  return dd

def run(value_dict=None,use_log=None,target_offset=None,min_ratio=None):
  return find_breakpoint(value_dict=value_dict,use_log=use_log,
   target_offset=target_offset,
      min_ratio=min_ratio,
   )

if __name__=="__main__":
  import sys,os
  args=sys.argv[1:]
  if args and os.path.isfile(args[0]):
    value_dict=get_value_dict(file_name=args[0])
    bb=run(value_dict=value_dict,use_log=input_use_log,
      target_offset=input_target_offset,
      min_ratio=input_min_ratio,
      )
  else:
    bb=run(value_dict=input_value_dict,use_log=input_use_log,
      target_offset=input_target_offset,
      min_ratio=input_min_ratio,
       )
