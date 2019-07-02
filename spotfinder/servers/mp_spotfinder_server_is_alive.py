from __future__ import absolute_import, division, print_function
import httplib

def get_spotfinder_url(host,port):
  testurl = "%s:%d"%(host,port)
  Connection = httplib.HTTPConnection(testurl,strict=True)
  Connection.putrequest('POST',"/spotfinder")
  Connection.putheader('Expect',"200")
  Connection.endheaders()
  Connection.send("")
  Response = Connection.getresponse()
  print(Response.status,Response.reason)
  Connection.close()

get_spotfinder_url("localhost",8125)
