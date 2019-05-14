from __future__ import absolute_import, division, print_function
from six.moves import range
import time

class phil_validation:
  def __init__(self,param):

    self.param = param
    self.application_level_validation()

  def application_level_validation(self):
    pass

class file_table:
  def __init__(self,param,query,enforce80=False,enforce81=False):
    import urllib2
    auth_handler = urllib2.HTTPBasicAuthHandler()
    auth_handler.add_password(realm="Webservices Auth",
                          uri="https://pswww.slac.stanford.edu",
                          user=param.web.user,
                          passwd=param.web.password)
    opener = urllib2.build_opener(auth_handler)
    # ...and install it globally so it can be used with urlopen.
    urllib2.install_opener(opener)
    R = urllib2.urlopen(query)
    if R.getcode() != 200:
      print("Status",R.getcode())
    import xml.etree.ElementTree
    X = xml.etree.ElementTree.XML(R.read())
    #from IPython import embed; embed()#help(X)
    self.runs = []
    self.items = []
    self.times = []
    for child in X:
      if child.tag=="issued":
        for change in child:
          if change.tag == "change":
            for token in change:
              if token.tag == "run":
                #from IPython import embed; embed()
                self.runs.append(int(token.text))
              if token.tag == "item":
                self.items.append(token.text)
              if token.tag == "time":
                self.times.append(token.text)
    self.unixtimes = []
    for item in self.times:
      self.unixtimes.append(  time.mktime(time.strptime(item[:19],"%Y-%m-%dT%H:%M:%S"))  )
    self.rundict = {}
    for i in range(len(self.runs)):
      if enforce80: #assume the required FEE spectrometer data is in stream 80
        if self.items[i].find("-s80-") < 0: continue
      if enforce81: #assume the required FEE spectrometer data is in stream 81
        if self.items[i].find("-s81-") < 0: continue
      if self.rundict.get(self.runs[i],0)==0:
        self.rundict[self.runs[i]]=dict(items=[],unixtimes=[])
      self.rundict[self.runs[i]]["items"].append(self.items[i])
      self.rundict[self.runs[i]]["unixtimes"].append(self.unixtimes[i])

  def get_runs(self,filter=None):
    values = []
    for key in self.rundict:
      if filter is None or (filter[0]<=key and key<=filter[1]):
        values.append(dict(run=key,time=max(self.rundict[key]["unixtimes"])))
    return values
  pass

class application:
  def __init__(self,param):
    self.param = param
    query = self.get_query1()
    FT = file_table(param,query)
    runs = FT.get_runs(filter = self.param.runlimits)
    # finally sort them by time
    runs.sort ( key = lambda A : A["time"])

    # Now prepare for the simulation
    data_timespan = runs[-1]["time"] - runs[0]["time"]
    simulation_timespan = data_timespan / self.param.speedup.factor
    print("Simulation duration %5.2f sec"%simulation_timespan)
    import time
    begin = time.time()
    runptr = 0
    while runptr<len(runs):
      time.sleep(1)
      dataclock = runs[runptr]["time"] - runs[0]["time"]
      if dataclock < self.param.speedup.factor * ( time.time() - begin):
        print("Run %d, time %s"%(runs[runptr]["run"],time.asctime(time.localtime(runs[runptr]["time"]))))
        runptr+=1

  def get_query1(self):
    return "https://pswww.slac.stanford.edu/ws-auth/dataexport/placed?exp_name=%s"%(
      self.param.data[0])
