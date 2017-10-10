from __future__ import division
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
import http.client

def get_spotfinder_url(host,port):
  testurl = "%s:%d"%(host,port)
  Connection = http.client.HTTPConnection(testurl,strict=True)
  Connection.putrequest('POST',"/spotfinder")
  Connection.putheader('Expect',"200")
  Connection.endheaders()
  Connection.send("")
  Response = Connection.getresponse()
  print(Response.status,Response.reason)
  Connection.close()

get_spotfinder_url("localhost",8125)
