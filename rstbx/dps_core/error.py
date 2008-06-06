from labelit.exception import AutoIndexError
class NoAutoIndex(AutoIndexError):
  def __init__(self):
    AutoIndexError.__init__(self,"(couldn't find 3 good basis vectors)")
