from __future__ import absolute_import, division, print_function
from crys3d.regression import tests_HKLviewer
import asyncio, os.path, websockets, socket


global socket_connected
socket_connected = False

def find_free_port():
  import socket
  s = socket.socket()
  s.bind(('', 0))      # Bind to a free port provided by the host.
  port = s.getsockname()[1]
  s.close()
  return port


async def handler(websocket, path):
# WS server example
  while True:
    name = await websocket.recv()
    print(f"{name}")
    greeting = f"Hello {name}!"
    await websocket.send(greeting)
    print(greeting)
    if name=="Goodbye":
      global socket_connected
      await websocket.close()
      socket_connected = True
      return
    await asyncio.sleep(0.2)



websock_htmlstr = """
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<html><head><meta charset="utf-8" /></head>
<body>
<div id='mytext'></div>
<script>
document.getElementById('mytext').innerHTML = "Hoping to connect to localhost via websocket..."
var portnumber = %s;
//const mysocket = new WebSocket('ws://localhost:424242'); // testing connection failure
const mysocket = new WebSocket('ws://localhost:' + String(portnumber));
mysocket.addEventListener('open', function (event) {
  mysocket.send('Connection Established');
  mysocket.send('Goodbye');
  document.getElementById('mytext').innerHTML = "Websocket connection established to localhost:" + String(portnumber)
});
mysocket.onmessage = function(e) { console.log(e)};
mysocket.onopen = function(e) { console.log(e)  };

</script></body></html>

"""

def write_websocktest_html(port):
  with open("websocket_test.html","w") as f:
    f.write(websock_htmlstr %port)
  myurl = "file:///" + os.path.abspath( "websocket_test.html" )
  myurl = myurl.replace("\\", "/")
  _, webctrl = tests_HKLviewer.FindFireFox()
  #import webbrowser
  assert webctrl.open(myurl)


async def closing_time():
  dt = 0.2; t=0; maxtime = 60
  while t < maxtime:
    await asyncio.sleep(dt)
    t += dt
    global socket_connected
    if socket_connected:
      print("OK")
      asyncio.get_event_loop().call_soon(asyncio.get_event_loop().stop)
      return

  print('Timed out trying to connect to webbrowser')
  asyncio.get_event_loop().call_soon(asyncio.get_event_loop().stop)


if __name__ == '__main__':
  port = find_free_port()
  write_websocktest_html(port)
  print("Websockets server on localhost port %s waiting for browser connection." %port)

  tasks = asyncio.gather(
    websockets.serve(handler, "localhost", port),
    closing_time()
  )

  evl = asyncio.get_event_loop()
  evl.run_until_complete(tasks)
  evl.run_forever()
  assert socket_connected
