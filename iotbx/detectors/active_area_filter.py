class active_area_filter:
  NEAR = 2
  def __init__(self,IT):
    from scitbx import matrix
    self.IT = IT
    self.centers = []
    from annlib_ext import AnnAdaptor

    reference = flex.double()

    for i in xrange(len(IT)//4):
      UL = matrix.col((float(IT[4*i]),float(IT[4*i+1])))
      LR = matrix.col((float(IT[4*i+2]),float(IT[4*i+3])))
      center = (UL+LR)/2.
      reference.append(center[0])
      reference.append(center[1])
    self.adapt = AnnAdaptor(data=reference,dim=2,k=self.NEAR)
  def __call__(self,predictions,hkllist,pxlsz):
    query = flex.double()
    for pred in predictions:
      query.append(pred[0]/pxlsz)
      query.append(pred[1]/pxlsz)
    self.adapt.query(query)
    selection = flex.bool()
    self.tile_id = flex.int()
    assert len(self.adapt.nn)==len(predictions)*self.NEAR
    for p in xrange(len(predictions)):
      is_in_active_area = False
      for n in xrange(self.NEAR):
        itile = self.adapt.nn[p*self.NEAR+n]
        if self.IT[4*itile]<predictions[p][0]/pxlsz<self.IT[4*itile+2] and\
           self.IT[4*itile+1]<predictions[p][1]/pxlsz<self.IT[4*itile+3]:
          is_in_active_area = True;break
      if is_in_active_area:
        self.tile_id.append(itile)
      selection.append(is_in_active_area)
    assert selection.count(True) == len(self.tile_id)
    return predictions.select(selection),hkllist.select(selection)

