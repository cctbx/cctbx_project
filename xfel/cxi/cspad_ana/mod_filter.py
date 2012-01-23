from xfel.cxi.cspad_ana import cspad_tbx

class mod_filter(object):

  def __init__(self, timestamps_path):
    """ Pass a file containing a list of timestamps, or list of file names that
    contain timestamps. Only images whose timestamps appear in this list are
    filtered through to downwind module.
    """
    f = open(timestamps_path, 'rb')
    self.timestamps = f.read()
    f.close()
    self.nshots = 0

  def beginjob(self, evt, env):
    pass

  def event(self, evt, env):
    self.nshots += 1
    timestamp = cspad_tbx.evt_timestamp(evt)

    if timestamp not in self.timestamps:
      (evt.put(True, "skip_event"))
      #print "skipping shot %i" %self.nshots
      return
    #print "*** using shot %i ***" %self.nshots

  def endjob(self, env):
    pass
