try:
  from itertools import count
except:
  class count:

    def __init__(self, firstval=0):
      self.val = firstval

    def next(self):
      result = self.val
      self.val += 1
      return result

    def __iter__(self):
      return self
