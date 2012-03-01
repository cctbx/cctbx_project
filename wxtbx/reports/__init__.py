
import wx
import sys

class PDBLinkMixin (object) :
  """
  Subclass this along with wx.Frame, etc. to link to the PDB from a list of
  results.  get_pdb_id_for_viewing() must be re-implemented in subclasses.
  """
  def create_pdb_buttons (self, panel, sizer) :
    szr2 = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(szr2)
    btn1 = wx.Button(panel, -1, "Download selected structure")
    szr2.Add(btn1, 0, wx.ALL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnDownload, btn1)
    btn2 = wx.Button(panel, -1, "View PDB web page")
    szr2.Add(btn2, 0, wx.ALL, 5)
    self.Bind(wx.EVT_BUTTON, self.OnViewPDB, btn2)

  def get_pdb_id_for_viewing (self) :
    raise NotImplementedError()

  def OnViewPDB (self, event) :
    pdb_id = self.get_pdb_id_for_viewing()
    if (pdb_id is not None) :
      url = "http://www.rcsb.org/pdb/explore/explore.do?structureId=%s" %pdb_id
      if (sys.platform == "darwin") :
        from wxtbx import browser
        if (getattr(self, "web_frame", None) is None) :
          self.web_frame = browser.browser_frame(self, -1, "RCSB PDB")
          self.web_frame.SetHomepage("http://www.rcsb.org")
          self.Bind(wx.EVT_WINDOW_DESTROY, self.OnCloseWWW, self.web_frame)
          self.web_frame.Show()
        self.web_frame.Raise()
        self.web_frame.Open(url)
        self.web_frame.Refresh()
      else :
        from libtbx import easy_run
        easy_run.call("env LD_LIBRARY_PATH='' firefox %s" % url)

  def OnDownload (self, event) :
    pdb_id = self.get_pdb_id_for_viewing()
    if (pdb_id is not None) :
      def download (args) :
        import mmtbx.command_line.fetch_pdb
        return mmtbx.command_line.fetch_pdb.run2(args=args)
      from wxtbx.process_control import download_file_basic
      download_file_basic(
        window=self,
        dl_func=download,
        args=[pdb_id])

  def OnCloseWWW (self, event) :
    self.web_frame = None
