from __future__ import absolute_import, division, print_function
from crys3d.regression import tests_HKLviewer
import asyncio, os.path, websockets, socket

#async def run():
#  await tests_HKLviewer.exercise_websocket()

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
      asyncio.get_event_loop().stop()
      return
    await asyncio.sleep(0.2)


async def waithandler(websocket, path):
  # Wait for at most 5 second
  try:
    await asyncio.wait_for(handler(websocket, path), timeout=5.0)
  except asyncio.TimeoutError:
    print('Timed out trying to connect to webbrowser')
    asyncio.get_event_loop().stop()




websock_htmlstr = """
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<html><head><meta charset="utf-8" /></head>
<body><script>

var portnumber = %s;
const mysocket = new WebSocket('ws://localhost:' + String(portnumber));
mysocket.addEventListener('open', function (event) {
  mysocket.send('Connection Established');
  mysocket.send('Goodbye');
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


def main():
  port = find_free_port()
  write_websocktest_html(port)
  print("Websockets server on localhost port %s waiting for browser connection." %port)
  start_server = websockets.serve(waithandler, "localhost", port)
  asyncio.get_event_loop().run_until_complete(start_server)
  asyncio.get_event_loop().run_forever()




if __name__ == '__main__':
  main()
  print("OK")
