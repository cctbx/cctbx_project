class reflnlist:

  def __init__(self, file):
    line_1 = file.read_line().split()
    assert len(line_1) in (3,4)
    line_1 = [int(field) for field in line_1]
    if (len(line_1) == 3):
      line_1.append(0)
    self.NumInts = line_1[0]
    self.NumFloats = line_1[1]
    self.NumStrings = line_1[2]
    self.NumInfoLines = line_1[3]
    self.InfoLines = []
    for i in xrange(self.NumInfoLines):
      self.InfoLines.append(file.read_line().strip())
    self.column_info = []
    for i in xrange(self.NumInts + self.NumFloat + self.NumStrings):
      column_name = file.read_line().strip()
      assert len(column_name) > 1
      assert column_name[0] in "nfs"
