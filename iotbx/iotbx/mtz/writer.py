from cctbx.array_family import flex
import math

def _DoubleOrComplex(writer,label,datatype,item_miller,item_data):
  if isinstance(item_data,flex.double):
    writer.addColumn(label,datatype,item_miller,item_data)
  elif isinstance(item_data,flex.complex_double):
    writer.addColumn(label,datatype,item_miller,flex.abs(item_data))
    writer.addColumn("phi_"+label,"P",item_miller,flex.arg(item_data)*180./math.pi)

def _columnCombinations(writer,label,datatype,carry_miller,carry_data):
  if len(carry_miller)==1:  # anomalous flag is False
    _DoubleOrComplex(writer,label,datatype,carry_miller[0],carry_data[0])
  elif len(carry_miller)==2:  # anomalous data
    _DoubleOrComplex(writer,label+"+",datatype,carry_miller[0],carry_data[0])
    _DoubleOrComplex(writer,label+"-",datatype,carry_miller[1],carry_data[1])

def add_miller_array(self, miller_array, label_data, label_sigmas=None):

  assert (miller_array.sigmas() == None) == (label_sigmas == None)
  anomalous_flag = miller_array.anomalous_flag()
  miller_indices = miller_array.indices()
  data = miller_array.data()
  sigmas = miller_array.sigmas()

  if (anomalous_flag):

    carry_miller = [flex.miller_index(),flex.miller_index()]
    if isinstance(data,flex.double):
      carry_data = [flex.double(),flex.double()]
    elif isinstance(data,flex.complex_double):
      carry_data = [flex.complex_double(),flex.complex_double()]
    if sigmas:
      carry_sigma= [flex.double(),flex.double()]

    asu, matches = miller_array.match_bijvoet_mates()
    for i,j in matches.pairs():
      for h,k in zip(asu.indices()[i], asu.indices()[j]): assert h == -k
      carry_miller[0].push_back(asu.indices()[i])
      carry_miller[1].push_back(asu.indices()[i])
      carry_data[0].push_back(asu.data()[i])
      carry_data[1].push_back(asu.data()[j])
      if sigmas:
        carry_sigma[0].push_back(asu.sigmas()[i])
        carry_sigma[1].push_back(asu.sigmas()[j])

    reciprocal_space_asu = asu.space_group_info().reciprocal_space_asu()
    for i in matches.singles():
      h = asu.indices()[i]
      f = asu.data()[i]
      if sigmas:
        g = asu.sigmas()[i]
      s = reciprocal_space_asu.which(h)
      assert s != 0
      if (reciprocal_space_asu.which(h) > 0):
        carry_miller[0].push_back(h)
        carry_data[0].push_back(f)
        if sigmas:
          carry_sigma[0].push_back(g)
      else:
        carry_miller[1].push_back((-h[0],-h[1],-h[2]))
        carry_data[1].push_back(f)
        if sigmas:
          carry_sigma[1].push_back(g)
    _columnCombinations(self,label_data,"G",carry_miller,carry_data)
    if sigmas:
      _columnCombinations(self,label_sigmas,"L",carry_miller,carry_sigma)

  else:
    D = miller_array.data()
    _columnCombinations(self,label_data,"F",[miller_array.indices()],[miller_array.data()])
    if sigmas:
      _columnCombinations(self,label_sigmas,"Q",[miller_array.indices()],[miller_array.sigmas()])
