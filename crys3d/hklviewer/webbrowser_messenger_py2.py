from __future__ import absolute_import, division, print_function
import traceback
from libtbx.utils import Sorry, to_str
import threading
import time

from websocket_server import WebsocketServer

class WBmessenger(object):
  def __init__(self, viewerparent ):
    try:
      self.parent = viewerparent
      self.ProcessBrowserMessage = self.parent.ProcessBrowserMessage
      self.websockport = self.parent.websockport
      self.sleeptime = self.parent.sleeptime
      self.mprint = self.parent.mprint
      self.parent.lastviewmtrx
      self.browserisopen = False
      self.msgqueue = []
      self.msgdelim = ":\n"
      self.ishandling = False
      self.websockclient = None
      self.isterminating = False
      self.was_disconnected = None
      self.mywebsock = None
      self.websockeventloop = None
    except Exception as e:
      print( to_str(e) + "\n" + traceback.format_exc(limit=10))


  def Sleep(self, t):
    time.sleep(t)


  def OnWebsocketClientMessage(self, client, server, message):
    self.ProcessBrowserMessage(message)


  def StartWebsocket(self):
    self.server = WebsocketServer(self.websockport, host='127.0.0.1')
    if not self.server:
      raise Sorry("Could not connect to web browser")
    self.server.set_fn_new_client(self.OnConnectWebsocketClient)
    self.server.set_fn_client_left(self.OnDisconnectWebsocketClient)
    self.server.set_fn_message_received(self.OnWebsocketClientMessage)
    self.wst = threading.Thread(target=self.server.run_forever)
    self.wst.daemon = True
    self.wst.start()
    self.msgqueuethrd = threading.Thread(target = self.WebBrowserMsgQueue )
    self.msgqueuethrd.daemon = True
    self.msgqueuethrd.start()


  def StopWebsocket(self):
    try:
      if self.websockclient: # might not have been created if program is closed before a data set is shown
        self.websockclient['handler'].send_text(u"", opcode=0x8)
    except Exception as e:
      self.mprint( to_str(e) + "\n" + traceback.format_exc(limit=10), verbose=0)
    self.mprint("Shutting down Websocket listening thread", verbose=1)
    self.server.shutdown()
    self.parent.javascriptcleaned = True
    self.msgqueuethrd.join()
    self.mprint("Shutting down WebsocketServer", verbose=1)
    self.wst.join()
    self.isterminating = True


  def AddToBrowserMsgQueue(self, msgtype, msg=""):
    self.msgqueue.append( (msgtype, msg) )


  def WebBrowserMsgQueue(self):
    try:
      while True:
        nwait = 0.0
        time.sleep(self.sleeptime)
        if self.parent.javascriptcleaned:
          self.mprint("Shutting down WebBrowser message queue", verbose=1)
          return
        if len(self.msgqueue):
          pendingmessagetype, pendingmessage = self.msgqueue[0]
          gotsent = self.send_msg_to_browser(pendingmessagetype, pendingmessage)
          while not self.browserisopen:  #self.websockclient:
            time.sleep(self.sleeptime)
            nwait += self.sleeptime
            if nwait > self.parent.handshakewait or self.parent.javascriptcleaned or not self.viewerparams.scene_id is not None:
              return
          if gotsent:
            self.msgqueue.remove( self.msgqueue[0] )
          #if self.was_disconnected:
          #  nwait2 = 0.0
          #  while nwait2 < self.parent.handshakewait:
          #    nwait2 += self.sleeptime
          #  self.ReloadNGL()
# if the html content is huge the browser will be unresponsive until it has finished
# reading the html content. This may crash this thread. So try restarting this thread until
# browser is ready
    except Exception as e:
      self.mprint( str(e) + ", Restarting WebBrowserMsgQueue\n" \
                          + traceback.format_exc(limit=10), verbose=2)
      self.websockclient = None
      self.WebBrowserMsgQueue()


  def OnConnectWebsocketClient(self, client, server):
    self.websockclient = client
    self.mprint( "Browser connected:" + str( self.websockclient ), verbose=1 )
    if self.was_disconnected:
      self.was_disconnected = False
    if self.parent.lastviewmtrx and self.parent.viewerparams.scene_id is not None:
      self.parent.set_volatile_params()
      self.mprint( "Reorienting client after refresh:" + str( self.websockclient ), verbose=2 )
      self.AddToBrowserMsgQueue("ReOrient", self.parent.lastviewmtrx)
    else:
      self.parent.SetAutoView()


  def OnDisconnectWebsocketClient(self, client, server):
    self.mprint( "Browser disconnected:" + str( client ), verbose=1 )
    self.was_disconnected = True


  def send_msg_to_browser(self, msgtype, msg=""):
    message = u"" + msgtype + self.msgdelim + str(msg)
    if self.websockclient:
      nwait = 0.0
      while not ("Ready" in self.parent.lastmsg or "tooltip_id" in self.parent.lastmsg \
        or "CurrentViewOrientation" in self.parent.lastmsg or "AutoViewSet" in self.parent.lastmsg \
        or "ReOrient" in self.parent.lastmsg or "JavaScriptCleanUp" in self.parent.lastmsg or self.websockclient is None):
        time.sleep(self.sleeptime)
        nwait += self.sleeptime
        if nwait > 2.0 and self.browserisopen:
          self.mprint("ERROR: No handshake from browser!", verbose=0 )
          self.mprint("failed sending " + msgtype, verbose=1)
          self.mprint("Reopening webpage again", verbose=0)
          break
    if self.browserisopen and self.websockclient is not None:
      try:
        self.server.send_message(self.websockclient, message )
        return True
      except Exception as e:
        self.mprint( str(e) + "\n" + traceback.format_exc(limit=10), verbose=1)
        self.websockclient = None
        return False
    else:
      return self.parent.OpenBrowser()
