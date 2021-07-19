from __future__ import absolute_import, division, print_function

from six.moves import range
class MultiQueue(object):

  def __init__(self):

    self.queue_for = {}


  def create(self, name):

    from six.moves import queue
    self.queue_for[ name ] = queue.Queue()


  def remove(self, name):

    del self.queue_for[ name ]


  def put(self, name, item, block = True, timeout = None):

    self.queue_for[ name ].put( item, block = block, timeout = timeout )


  def get(self, name, block = True, timeout = None):

    return self.queue_for[ name ].get( block = block, timeout = timeout )


  def put_nowait(self, name, item):

    return self.queue_for[ name ].put_nowait( item )


  def get_nowait(self, name):

    return self.queue_for[ name ].get_nowait()


class Queue(object):

  def __init__(self, server):

    import socket
    import os
    self.identifier = "%s-%s-%s" % ( socket.getfqdn(), os.getpid(), id( self ) )

    self.server = server
    self.server.multiqueue.create( self.identifier )


  def shutdown(self):

    self.server.multiqueue.remove( self.identifier )
    self.server = None


  def put(self, value, block = True, timeout = None):

    self.server.multiqueue.put(
      self.identifier,
      value,
      block = block,
      timeout = timeout,
      )


  def get(self, block = True, timeout = None):

    return self.server.multiqueue.get(
      self.identifier,
      block = block,
      timeout = timeout,
      )


  def put_nowait(self, value):

    return self.server.multiqueue.put_nowait( self.identifier, value )


  def get_nowait(self):

    return self.server.multiqueue.get_nowait( self.identifier )


class Manager(object):

  def __init__(self, manager):

    self.manager = manager
    self.multiqueue = manager.get() # caching


  def Queue(self):

    return Queue( server = self )


  def __getstate__(self):

    result = self.__dict__.copy()
    result[ "address" ] = self.manager.address
    result[ "authkey" ] = str( self.manager._authkey )
    del result[ "manager" ]
    del result[ "multiqueue" ]
    return result


  def __setstate__(self, result):

    manager = self.get_client_manager(
      address = result[ "address" ],
      authkey = result[ "authkey" ],
      )
    result[ "manager" ] = manager
    result[ "multiqueue" ] = manager.get()
    del result[ "address" ]
    del result[ "authkey" ]
    self.__dict__ = result


  @classmethod
  def Server(cls, port = 0, keylength = 16):

    from multiprocessing.managers import BaseManager

    class QManager(BaseManager):

      pass

    multiqueue = MultiQueue()
    QManager.register( "get", lambda: multiqueue )

    import socket, string
    import random

    manager = QManager(
      address = ( socket.getfqdn(), port ),
      authkey = "".join(
        random.choice( string.ascii_letters ) for i in list(range( keylength))
        ),
      )
    manager.start()

    return cls( manager = manager )


  @classmethod
  def Client(cls, server, port, authkey):

    return cls(
      manager = cls.get_client_manager(
        address = (server, port ),
        autkey = authkey,
        )
      )


  @staticmethod
  def get_client_manager(address, authkey):

    from multiprocessing.managers import BaseManager

    class QManager(BaseManager):

      pass

    QManager.register( "get" )

    manager = QManager( address = address, authkey = authkey )
    manager.connect()

    return manager


class QFactory(object):
  """
  Creator pattern for socket queue, also include destruction
  """

  def __init__(self, port = 0, keylength = 16):

    self.manager = Manager.Server( port = port, keylength = keylength )


  def create(self):

    return self.manager.Queue()


  @staticmethod
  def destroy(queue):

    queue.shutdown()
