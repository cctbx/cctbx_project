from cctbx.array_family import flex

def _DoubleOrComplex(writer,label,datatype,item_miller,item_data):
  if isinstance(item_data,flex.double):
    writer.addColumn(label,datatype,item_miller,item_data)
  elif isinstance(item_data,flex.complex_double):
    writer.addColumn(label,datatype,item_miller,flex.abs(item_data))
    writer.addColumn("phi_"+label,"P",item_miller,flex.arg(item_data, True))

def _columnCombinations(writer,label,datatype,carry_miller,carry_data):
  if len(carry_miller)==1:  # anomalous flag is False
    _DoubleOrComplex(writer,label,datatype,carry_miller[0],carry_data[0])
  elif len(carry_miller)==2:  # anomalous data
    _DoubleOrComplex(writer,label+"+",datatype,carry_miller[0],carry_data[0])
    _DoubleOrComplex(writer,label+"-",datatype,carry_miller[1],carry_data[1])

def add_miller_array(self, miller_array, label_data, label_sigmas=None):
  assert (miller_array.sigmas() == None) == (label_sigmas == None)
  if (miller_array.anomalous_flag()):
    carry_miller = []
    carry_data = []
    carry_sigma = []
    asu, matches = miller_array.match_bijvoet_mates()
    for i,s in ((0,"+"),(1,"-")):
      sel = matches.pairs_hemisphere_selection(s) \
          | matches.singles_hemisphere_selection(s)
      if (i == 0):
        carry_miller.append(asu.indices().select(sel))
      else:
        carry_miller.append(-asu.indices().select(sel))
      carry_data.append(asu.data().select(sel))
      if (asu.sigmas() != None):
        carry_sigma.append(asu.sigmas().select(sel))
    _columnCombinations(self,label_data,"G",carry_miller,carry_data)
    if (asu.sigmas()):
      _columnCombinations(self,label_sigmas,"L",carry_miller,carry_sigma)
  else:
    _columnCombinations(self,label_data,"F",
      [miller_array.indices()],[miller_array.data()])
    if (miller_array.sigmas()):
      _columnCombinations(self,label_sigmas,"Q",
        [miller_array.indices()],[miller_array.sigmas()])

def Hall_to_CCP4_Converter(spacegroup):
  symbol = spacegroup.type().lookup_symbol().replace(" ","")
  return symbol

def setSpaceGroup(self, spacegroup, symbol=None):
  if symbol==None:
    symbol = Hall_to_CCP4_Converter(spacegroup)
  assert not " " in symbol
  assert len(symbol)<=10
  self.ll_setSpaceGroup(spacegroup,symbol)
  
