"""qttbx.viewers.chimerax: centering ChimeraX on the residues of a Phenix atom
selection. The Phenix GUI sends the selection of a clicked validation row;
ChimeraX cannot center on the x,y,z that Coot uses, only on atoms it can
name, so the selection becomes a residue spec, each residue kept to its copy
near the row's x,y,z. Residues, not atoms: Probe clash atoms are often
hydrogens that are not in the model ChimeraX shows. Probe clash dots are
drawn by the PhenixUI bundle's "probe" format.

No ChimeraX needed: a local HTTP server stands in for the REST server. It
decodes requests as ChimeraX's does (parse_qs on the query) and answers in
plain text with no charset, as ChimeraX does."""

import os
import subprocess
import sys
import tempfile
import threading
import time
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import parse_qs, urlparse

from libtbx.utils import Sorry, format_cpu_times
from qttbx.viewers import residue_selection_terms
from qttbx.viewers.chimerax import ChimeraXViewer, chimerax_residue_spec


def exercise_residue_spec():
  """Selections from validation results, and the forms older Phenix
  versions stored in results that restored jobs still carry."""
  for selection, expected in [
      # mmtbx.validation results
      ("chain A and resid 42", "/A:42"),
      ("chain A and resid 42A", "/A:42A"),
      ("chain A and resid 42 and name CB", "/A:42"),
      ("(chain A and resid 2 and name CA) or (chain A and resid 2A and name N)",
       "/A:2 | /A:2A"),
      ("(chain A and resid 7 and name CA) or (chain A and resid 7 and name C)",
       "/A:7"),
      ("chain B and resid A000", "/B:10000"),
      ("chain D and resid -3", "/D:-3"),
      ("chain C and resid 1 and name O5'", "/C:1"),
      ("chain C and resid 2 and name O5\\*", "/C:2"),
      ("CHAIN A AND RESID 5", "/A:5"),
      # a blank chain ID: ":5" alone would match residue 5 of every chain
      ("chain ' ' and resid 5", '(:5 & ::chain_id=" ")'),
      # missing atoms, sequence validation (stored by older versions);
      # resseq ignores the insertion code, and so does 42-42
      ("chain 'A' and resseq 42", "/A:42-42"),
      ("chain 'A' and resid   42 ", "/A:42"),
      ("chain 'A' and resid   42A", "/A:42A"),
      # barbed wire, sequence viewer ranges; ,20-20 keeps 20A, as resseq does
      ("(chain A and resseq 10:20)", "/A:10-20,20-20"),
      ("chain 'A' and (resseq 10:20 or resseq 25)",
       "/A:10-20,20-20 | /A:25-25"),
      # resid 10:20 stops at 20 itself, as 10-20 does; resid 42:42 is 42
      ("chain A and resid 10:20", "/A:10-20"),
      ("chain A and resid 42:42", "/A:42"),
      ("chain A and resid -12:-3", "/A:-12--3"),
      # a negative number, a lowercase insertion code (ChimeraX reads it as
      # a capital from PDB files, as it is from mmCIF), hybrid-36 with an
      # insertion code
      ("chain A and resid -12A", "/A:-12A"),
      ("chain A and resid 44a", "/A:44a,44A"),
      ("chain A and resid A000B", "/A:10000B"),
      # and binds tighter than or
      ("chain A and resid 1 or chain B and resid 2", "/A:1 | /B:2"),
      ("chain A and resseq 42 and icode B", "/A:42B"),
      # chain IDs that are not letters and digits: /A-2:5 matches nothing
      ("chain A-2 and resid 5", '(:5 & ::chain_id="A-2")'),
      ]:
    spec = chimerax_residue_spec(selection)
    assert spec == expected, (selection, spec, expected)


def exercise_residue_spec_unsupported():
  """Anything that does not name residues gives None, so the caller leaves
  ChimeraX alone instead of framing the whole structure or a wrong part."""
  for selection in [
      None, "", "all", "not water", "chain A", "segid ABCD",
      "resid 5 and resid 6", "chain A and within(5, resid 3)",
      " A  42  ALA",                 # an id_str, not a selection
      "chain A or resid 5", "chain A and resid 5 and not name CA",
      "chain A and resid", "(chain A and resid 5",
      "chain A and resid 5 and segid ABCD",
      "chain A and resid 5 and water",
      # left over after a whole selection
      "chain A and resid 5)", "resid 5 chain B",
      "chain A and resid 42AB", "chain A and resid 42*",
      # an insertion code cannot bound a range
      "chain A and resseq 10:20 and icode B", "chain A and resid 10A:20",
      ]:
    spec = chimerax_residue_spec(selection)
    assert spec is None, (selection, spec)
    assert residue_selection_terms(selection) is None, selection


_DROP = object()


class _Recorder(BaseHTTPRequestHandler):
  """reply(command) gives the text ChimeraX would log, (seconds, text) to
  answer late (or when released, see _release), or _DROP to close the
  connection unanswered; probe files are read while their open command
  runs, as ChimeraX does."""
  commands = []
  opened = []
  reply = staticmethod(lambda command: "ok")
  status_delay = 0  # seconds before answering cmdline.html
  release = threading.Event()
  answering = threading.Event()
  def do_GET(self):
    _Recorder.answering.set()
    try:
      self.answer()
    finally:
      _Recorder.answering.clear()
  def answer(self):
    if urlparse(self.path).path.endswith("cmdline.html"):
      _Recorder.release.wait(_Recorder.status_delay)
    query = parse_qs(urlparse(self.path).query)
    text = ""
    for command in query.get("command", []):
      _Recorder.commands.append(command)
      for part in command.split(";"):
        if part.strip().startswith("open") and "format probe" in part:
          with open(part.split('"')[1], encoding="utf-8") as f:
            _Recorder.opened.append(f.read())
      reply = _Recorder.reply(command)
      if reply is _DROP:
        self.close_connection = True
        return
      if isinstance(reply, tuple):
        _Recorder.release.wait(reply[0])
        reply = reply[1]
      text += reply
    try:
      self.send_response(200)
      self.send_header("Content-Type", "text/plain")
      self.end_headers()
      self.wfile.write(text.encode("utf-8"))
    except (BrokenPipeError, ConnectionResetError):
      pass  # the client stopped waiting
  def log_message(self, *args):
    pass


class _ChimeraX(object):
  """A viewer talking to _Recorder, which answers one request at a time,
  as ChimeraX's REST server does."""
  def __enter__(self):
    self.server = HTTPServer(("127.0.0.1", 0), _Recorder)
    # a short poll keeps shutdown quick
    threading.Thread(target=self.server.serve_forever, args=(0.05,),
                     daemon=True).start()
    self.viewer = ChimeraXViewer()
    self.viewer.url = "http://127.0.0.1:%d/" % self.server.server_address[1]
    del _Recorder.commands[:]
    del _Recorder.opened[:]
    _Recorder.reply = staticmethod(lambda command: "ok")
    _Recorder.status_delay = 0
    _Recorder.release = threading.Event()
    return self.viewer
  def __exit__(self, *args):
    _Recorder.release.set()
    self.server.shutdown()
    self.server.server_close()


def _release():
  """Let a late answer go now that the client stopped waiting for it, and
  wait until the server is free (it answers one request at a time)."""
  _Recorder.release.set()
  for _ in range(500):
    if not _Recorder.answering.is_set():
      break
    time.sleep(0.01)
  assert not _Recorder.answering.is_set(), "the server is stuck"
  _Recorder.release = threading.Event()


def _clip_list(near=None, far=None):
  """What "clip list" logs for near and far planes at these depths along the
  line of sight (each command wraps across lines and repeats as a link)."""
  lines = []
  for name, depth, sign in (("near", near, -1), ("far", far, 1)):
    if depth is not None:
      command = ("clip %s 0 axis 0.000,0.000,%.3f position\n1.153,2.259,%.4f"
                 " coordinateSystem scene" % (name, sign, -depth))
      lines.append("[%s](cxcmd:%s)" % (command, command))
  if not lines:
    return "Clipping is off\n"
  return "Using %d clip planes:\n%s\n" % (len(lines), "\n".join(lines))


_QUERY = "info models #9999; clip list"
_ZONE = "#9999 :<6"


_MARKER_INFO = 'model id #9999 type MarkerSet name "Phenix zoom point"\n'


def _session(state="", near=(), bracket=11.0):
  """ChimeraX's replies to a click: the state query, then the marker
  request (a residue spec in near has a copy within 6 A), then the view,
  where view puts the planes bracket A apart. echo writes to the Log panel
  of the ChimeraX window, not to the reply."""
  def reply(command):
    if command == _QUERY:
      return state
    text = ""
    for part in command.split(";"):
      part = part.strip()
      if part == "info models #9999":
        text += _MARKER_INFO
      elif part.startswith("info residues ("):
        spec = part[len("info residues ("):part.index(") & (")]
        if spec in near:
          text += "residue id #1/A:1 name ALA index 3\n"
      elif part == "clip list":
        text += _clip_list(near=-bracket / 2, far=bracket / 2)
    return text
  return reply


def _requests():
  return [[c.strip() for c in command.split(";")]
          for command in _Recorder.commands]


def _view(target):
  return ["show %s atoms" % target, "view %s clip false" % target,
          "zoom pixelSize 0.04"]


def _probe(xyz, *specs):
  """The marker request: the marker, then each residue's copies near it,
  each answer after the marker's info line"""
  commands = ["marker #9999 position %.3f,%.3f,%.3f radius 0.1" % xyz,
              'rename #9999 "Phenix zoom point"']
  for spec in specs:
    commands += ["info models #9999",
                 "info residues (%s) & (%s)" % (spec, _ZONE)]
  return commands


def exercise_zoom_to_selection():
  """Show the residues' atoms (view frames displayed atoms only) and view
  them, each kept to its copy within 6 A of the row's x,y,z: another model
  with the same chain and residue IDs, or a segid model whose rows name a
  chain that is somewhere else, is left out. Each residue is a part of its
  own in parentheses (a zone takes in everything after it). The point is a
  temporary marker, named, and closed by its own request even when view
  fails. Clipping is off and stays off; the user's ChimeraX selection is
  not touched."""
  with _ChimeraX() as viewer:
    _Recorder.reply = staticmethod(_session(near=("/A:2", "/A:2A")))
    viewer.zoom_to_selection(
      "(chain A and resid 2 and name CA) or (chain A and resid 2A and name N)",
      xyz=(1.0, -2.5, 3.25))
    requests = _requests()
    assert len(requests) == 4, requests
    assert _Recorder.commands[0] == _QUERY, requests
    assert requests[1] == _probe((1.0, -2.5, 3.25), "/A:2", "/A:2A"), requests
    assert requests[2] == _view(
      "((/A:2) & (%s)) | ((/A:2A) & (%s))" % (_ZONE, _ZONE)), requests
    assert requests[3] == ["close #9999"], requests
    assert not any(c.split()[0] in ("select", "echo")
                   for r in requests for c in r)


def exercise_zoom_residues_not_near():
  """A residue without a copy near the x,y,z is viewed wherever it is, so
  a pick of residues far apart keeps them all, and a model moved in
  ChimeraX is still centered. Ranges are not looked for."""
  for selection, near, probed, target in [
      ("chain A and (resseq 42 or resseq 47)", ("/A:47-47",),
       ("/A:42-42", "/A:47-47"),
       "(/A:42-42) | ((/A:47-47) & (%s))" % _ZONE),
      ("chain A and resid 42", (), ("/A:42",), "(/A:42)"),
      ("chain A and (resseq 10:20 or resseq 25)",
       ("/A:10-20,20-20", "/A:25-25"), ("/A:25-25",),
       "(/A:10-20,20-20) | ((/A:25-25) & (%s))" % _ZONE),
      ]:
    with _ChimeraX() as viewer:
      _Recorder.reply = staticmethod(_session(near=near))
      viewer.zoom_to_selection(selection, xyz=(1.0, 2.0, 3.0))
      requests = _requests()
      assert requests[1] == _probe((1.0, 2.0, 3.0), *probed), requests
      assert requests[2] == _view(target), (selection, requests)
      assert requests[3] == ["close #9999"], requests


def exercise_zoom_without_marker():
  """No marker: ranges (sequence viewer, barbed wire) are not narrowed to
  the residues around their centroid, rows without coordinates are viewed
  anywhere, and a model of the user's numbered #9999 is not touched. A
  marker left behind by a click that timed out is closed by the next."""
  ours = 'model id #9999 type MarkerSet name "Phenix zoom point"\n'
  for selection, xyz, state, view in [
      ("(chain A and resseq 10:20)", (1.0, 2.0, 3.0), "",
       _view("(/A:10-20,20-20)")),
      ("chain A and resid 42", None, "", _view("(/A:42)")),
      ("chain A and resid 42", (1.0, 2.0, 3.0),
       "model id #9999 type AtomicStructure name mine.pdb\n",
       _view("(/A:42)")),
      ("chain A and resid 42", None, ours, ["close #9999"] + _view("(/A:42)")),
      ]:
    with _ChimeraX() as viewer:
      _Recorder.reply = staticmethod(_session(state=state))
      viewer.zoom_to_selection(selection, xyz=xyz)
      requests = _requests()
      assert len(requests) == 2, (selection, requests)
      assert requests[1] == view, (selection, requests)
  with _ChimeraX() as viewer:
    _Recorder.reply = staticmethod(_session(state=ours, near=("/A:42",)))
    viewer.zoom_to_selection("chain A and resid 42", xyz=(1.0, 2.0, 3.0))
    requests = _requests()
    assert requests[1] == ["close #9999"] + _probe((1.0, 2.0, 3.0), "/A:42"), \
      requests
    assert requests[2] == _view("((/A:42) & (%s))" % _ZONE), requests


def exercise_zoom_keeps_clip_planes():
  """view resets or removes the near and far planes; they are put back
  without cutting into the residues. view (clip true) places both planes
  just clear of the residues. Two planes stay that way unless the user's
  were farther apart; then they are put back that far apart, centered on
  the residues (clip's default position is where the planes are, not the
  residues). One plane alone: the other is turned off."""
  target = "((/A:42) & (%s))" % _ZONE
  view = ["show %s atoms" % target, "view %s" % target]
  zoom = "zoom pixelSize 0.04"
  for state, expected in [
      # 10 A apart: thinner than view's 11 A
      (_clip_list(near=-5, far=5), [view + [zoom, "clip list"]]),
      # 14 A apart
      (_clip_list(near=-2, far=12), [
        view + [zoom, "clip list"],
        ["clip near -7.00 far 7.00 position %s" % target]]),
      (_clip_list(near=0), [view + ["clip far off", zoom]]),
      (_clip_list(far=0), [view + ["clip near off", zoom]]),
      ]:
    with _ChimeraX() as viewer:
      _Recorder.reply = staticmethod(_session(state=state, near=("/A:42",),
                                              bracket=11))
      viewer.zoom_to_selection("chain A and resid 42", xyz=(1.0, 2.0, 3.0))
      requests = _requests()
      assert requests[2:-1] == expected, (state, requests)
      assert requests[-1] == ["close #9999"], requests


def exercise_zoom_timeouts():
  """ChimeraX stops answering (busy with something else): the click gives
  up after one wait instead of sending the rest. A marker left behind is
  closed by the next click."""
  def late(step):
    def reply(command):
      text = _session(near=("/A:42",))(command)
      if (step == "query" and command == _QUERY or
          step == "marker" and command.startswith("marker") or
          step == "view" and command.startswith("show")):
        return (60, text)
      return text
    return reply
  for step, sent in [("query", 1), ("marker", 2), ("view", 3)]:
    with _ChimeraX() as viewer:
      _Recorder.reply = staticmethod(late(step))
      start = time.time()
      response = viewer.zoom_to_selection("chain A and resid 42",
                                          xyz=(1.0, 2.0, 3.0), timeout=0.3)
      elapsed = time.time() - start
      _release()
      assert response is None, step
      assert elapsed < 10, (step, elapsed)
      assert len(_Recorder.commands) == sent, (step, _requests())


def exercise_zoom_to_unsupported_selection():
  """A selection without residues sends nothing."""
  with _ChimeraX() as viewer:
    assert viewer.zoom_to_selection("all", xyz=(0, 0, 0)) is None
    assert viewer.zoom_to_selection(None) is None
    assert _Recorder.commands == [], _Recorder.commands


def exercise_zoom_without_chimerax():
  """ChimeraX closed by the user: the click does nothing and raises
  nothing."""
  with _ChimeraX() as viewer:
    pass
  assert viewer.zoom_to_selection("chain A and resid 42", (1, 2, 3)) is None


_PROBE = (
  ":1->2:wc: A   1 ASP  CA  : A  43 VAL  O   :0.262:0.407:46.294:12.067:8.673:0.000:0.0044:C:O:46.294:12.067:8.673:21.97:11.10\n"
  ":1->2:so: A   5 LEU HD11 : A  20 ILE HG12 :-0.31:-0.05:1.0:2.0:3.0:0.000:0.0044:C:C:1.0:2.0:3.0:21.97:11.10\n"
  ":1->2:hb: A   7 SER  OG  : A  30 ASP OD1  :0.100:0.200:1.0:2.0:3.0:0.000:0.0044:O:O:1.0:2.0:3.0:21.97:11.10\n"
  ":1->2:bo: A   9 ARG HH11 : A  12 GLU HG2  :-0.52:-0.10:1.0:2.0:3.0:0.000:0.0044:N:C:1.0:2.0:3.0:21.97:11.10\n")


def _files(tmp, probe_text=_PROBE):
  """A model in a directory with a non-ASCII name, another copy in a
  directory with a space, and a probe file."""
  model = os.path.join(tmp, "café", "model.pdb")
  other = os.path.join(tmp, "dir with space", "model.pdb")
  for f in (model, other):
    os.makedirs(os.path.dirname(f))
    open(f, "w").close()
  probe = os.path.join(tmp, "probe.txt")
  with open(probe, "w") as f:
    f.write(probe_text)
  return model, other, probe


def _chimerax_log(model, other, dots_on=None):
  """ChimeraX's replies: info models lists the open structures (and the
  dots, once loaded); a successful probe open logs its message."""
  def reply(command):
    text = ""
    if "format probe" in command:
      text += "Added probe dots to model.pdb #%s\n" % dots_on
    if "info models" in command:
      text += ("model id #1 type AtomicStructure filename %s\n"
               "model id #2 type AtomicStructure filename \"%s\"\n"
               "model id #3 type AtomicStructure filename %s\n"
               "model id #4 type AtomicStructure filename %sx\n"
               % (model, other, model, model))
      if dots_on is not None and "format probe" in command:
        text += 'model id #%s.1 type Model name "Probe markups"\n' % dots_on
    return text
  return reply


def exercise_load_probe_dots():
  """The bad overlaps (the clashes the table lists) go onto the newest copy
  of the model file open in ChimeraX, found under a non-ASCII path; its old
  dots are closed first. The temporary probe file is removed afterwards."""
  with tempfile.TemporaryDirectory() as tmp, _ChimeraX() as viewer:
    model, other, probe = _files(tmp)
    _Recorder.reply = staticmethod(_chimerax_log(model, other, dots_on="3"))
    assert viewer.load_probe_dots(probe, model) is not None
    assert _Recorder.commands[0].startswith("info models"), \
      _Recorder.commands
    commands = [c.strip() for c in _Recorder.commands[1].split(";")]
    assert commands[0] == 'close #3 & ##name="Probe markups"', commands
    words = commands[1].split('"')
    assert words[0] == "open " and \
      words[2] == " format probe structureModel #3", commands
    assert commands[2].startswith("info models"), commands
    assert _Recorder.opened == [
      "".join(l + "\n" for l in _PROBE.splitlines() if ":bo:" in l)
      ], _Recorder.opened
    assert not os.path.exists(words[1]), "temporary probe file left behind"
    # the copy in a directory with a space, reported quoted
    del _Recorder.commands[:]
    _Recorder.reply = staticmethod(_chimerax_log(model, other, dots_on="2"))
    viewer.load_probe_dots(probe, other)
    assert "structureModel #2" in _Recorder.commands[1], _Recorder.commands
    # the files of one model (e.g. the PDB and mmCIF files a refinement
    # wrote): whichever is open
    del _Recorder.commands[:]
    viewer.load_probe_dots(probe, [os.path.join(tmp, "model.cif"), other])
    assert "structureModel #2" in _Recorder.commands[1], _Recorder.commands


def exercise_load_probe_dots_fails():
  """ChimeraX without a probe reader (no PhenixUI 1.7.0 or newer) answers
  with an error and HTTP 200: the user is told to install it, not left
  without dots."""
  with tempfile.TemporaryDirectory() as tmp, _ChimeraX() as viewer:
    model, other, probe = _files(tmp)
    log = _chimerax_log(model, other)
    def reply(command):
      if "format probe" in command:
        return 'Invalid "format" argument: Should be one of ...\n'
      return log(command)
    _Recorder.reply = staticmethod(reply)
    try:
      viewer.load_probe_dots(probe, model)
    except Sorry as e:
      assert "PhenixUI" in str(e) and "Invalid \"format\"" in str(e), str(e)
    else:
      raise AssertionError("no Sorry")
    # any other failure is passed on without blaming PhenixUI
    def reply(command):
      if "format probe" in command:
        return "Permission denied: /tmp/x.probe\n"
      return log(command)
    _Recorder.reply = staticmethod(reply)
    try:
      viewer.load_probe_dots(probe, model)
    except Sorry as e:
      assert "PhenixUI" not in str(e) and "Permission denied" in str(e), \
        str(e)
    else:
      raise AssertionError("no Sorry")


def exercise_load_probe_dots_nothing_to_do():
  """The model is not open in ChimeraX, or the probe file is gone: nothing
  is opened. No bad overlaps (e.g. after fixing them and re-running Probe):
  the old dots are closed and nothing new is opened."""
  with tempfile.TemporaryDirectory() as tmp, _ChimeraX() as viewer:
    model, other, probe = _files(tmp)
    _Recorder.reply = staticmethod(lambda command: "")
    assert viewer.load_probe_dots(probe, model) is None
    assert len(_Recorder.commands) == 1, _Recorder.commands
    _Recorder.reply = staticmethod(_chimerax_log(model, other))
    del _Recorder.commands[:]
    assert viewer.load_probe_dots(probe + ".gone", model) is None
    assert _Recorder.commands == [], _Recorder.commands
  with tempfile.TemporaryDirectory() as tmp, _ChimeraX() as viewer:
    model, other, probe = _files(tmp, probe_text="".join(
      l + "\n" for l in _PROBE.splitlines() if ":bo:" not in l))
    _Recorder.reply = staticmethod(_chimerax_log(model, other))
    viewer.load_probe_dots(probe, model)
    assert _Recorder.commands[1] == 'close #3 & ##name="Probe markups"', \
      _Recorder.commands
    assert _Recorder.opened == [], _Recorder.opened


def exercise_load_probe_dots_no_answer():
  """ChimeraX drops the request: the temporary file is removed. ChimeraX
  is still reading it when the wait runs out: the file is left for it."""
  with tempfile.TemporaryDirectory() as tmp, _ChimeraX() as viewer:
    model, other, probe = _files(tmp)
    log = _chimerax_log(model, other)
    _Recorder.reply = staticmethod(
      lambda command: _DROP if "format probe" in command else log(command))
    assert viewer.load_probe_dots(probe, model) is None
    dots_file = _Recorder.commands[1].split('"')[3]
    assert not os.path.exists(dots_file), "temporary probe file left behind"
    del _Recorder.commands[:]
    _Recorder.reply = staticmethod(
      lambda command: (60, "") if "format probe" in command else log(command))
    assert viewer.load_probe_dots(probe, model, timeout=0.5) is None
    _release()
    dots_file = _Recorder.commands[1].split('"')[3]
    assert os.path.exists(dots_file), "removed while ChimeraX reads it"
    os.remove(dots_file)
    assert viewer.load_probe_dots(probe, None) is None


def exercise_busy_while_loading():
  """ChimeraX answers one request at a time, so while dots load a click or
  another load would wait for them: both are skipped, and is_alive answers
  without asking ChimeraX."""
  with tempfile.TemporaryDirectory() as tmp, _ChimeraX() as viewer:
    model, other, probe = _files(tmp)
    log = _chimerax_log(model, other, dots_on="3")
    _Recorder.reply = staticmethod(
      lambda command: (60, log(command)) if "format probe" in command
      else log(command))
    loading = threading.Thread(target=viewer.load_probe_dots,
                               args=(probe, model))
    loading.start()
    while len(_Recorder.commands) < 2:
      time.sleep(0.01)
    start = time.time()
    assert viewer.zoom_to_selection("chain A and resid 42", (1, 2, 3)) is None
    assert viewer.load_probe_dots(probe, model) is None
    assert viewer.is_alive()
    assert time.time() - start < 10, "waited for the load"
    _release()
    loading.join()
    assert len(_Recorder.commands) == 2, _Recorder.commands


def exercise_open_file_names():
  """Model and map files are opened by name, quoted: a space, '+' or '#'
  in the path would otherwise split or cut it. ChimeraX reads file names
  without escapes, so a Windows path is sent as it is."""
  with _ChimeraX() as viewer:
    paths = ["/a dir/model.pdb", "/a+b/map.ccp4", "/a#b/model.cif",
             "C:\\Users\\x y\\new\\model.pdb"]
    for path in paths:
      viewer.load_model(path)
    viewer.load_map(paths[1])
    expected = ['open "%s"' % path for path in paths + paths[1:2]]
    assert _Recorder.commands == expected, _Recorder.commands


def exercise_rest_timeouts():
  """ChimeraX answers one request at a time, so a long command keeps it
  busy. Opening files, other commands and the status check wait a limited
  time and return; a status check that times out means busy, not closed."""
  with _ChimeraX() as viewer:
    viewer._check_status()
    assert viewer._connected
    _Recorder.status_delay = 60
    start = time.time()
    viewer._check_status(timeout=0.3)
    elapsed = time.time() - start
    _release()
    assert viewer._connected, "a busy ChimeraX counted as closed"
    assert elapsed < 10, elapsed
    _Recorder.status_delay = 0
    _Recorder.reply = staticmethod(lambda command: (60, "late"))
    for send in (lambda: viewer.send_command(["version"], timeout=0.3),
                 lambda: viewer._load_file("/a/model.pdb", timeout=0.3)):
      start = time.time()
      response = send()
      elapsed = time.time() - start
      _release()
      assert response is None
      assert elapsed < 10, elapsed
  with _ChimeraX() as viewer:
    pass
  viewer._check_status()
  assert not viewer._connected, "a closed ChimeraX counted as running"


def exercise_close_busy_viewer():
  """Closing a busy ChimeraX (on exit from Phenix) waits a limited time;
  ChimeraX quits when it gets to the request."""
  with _ChimeraX() as viewer:
    _Recorder.reply = staticmethod(lambda command: (60, ""))
    viewer._connected = True
    viewer.process = subprocess.Popen(
      [sys.executable, "-c", "import time; time.sleep(60)"])
    try:
      start = time.time()
      viewer.close_viewer(timeout=0.3)
      elapsed = time.time() - start
      _release()
      assert elapsed < 10, elapsed
      assert _Recorder.commands == ["quit"], _Recorder.commands
      try:
        viewer.process.wait(timeout=1)
      except subprocess.TimeoutExpired:
        pass
      else:
        raise AssertionError("killed ChimeraX")
    finally:
      viewer.process.kill()
      viewer.process.wait()


def exercise_is_running():
  """is_running looks at the ChimeraX process, not its REST server, which
  may be busy for seconds."""
  viewer = ChimeraXViewer()
  assert not viewer.is_running()
  viewer.process = subprocess.Popen(
    [sys.executable, "-c", "import time; time.sleep(60)"])
  try:
    assert viewer.is_running()
  finally:
    viewer.process.kill()
    viewer.process.wait()
  assert not viewer.is_running()


if __name__ == "__main__":
  exercise_residue_spec()
  exercise_residue_spec_unsupported()
  exercise_zoom_to_selection()
  exercise_zoom_residues_not_near()
  exercise_zoom_without_marker()
  exercise_zoom_keeps_clip_planes()
  exercise_zoom_timeouts()
  exercise_zoom_to_unsupported_selection()
  exercise_zoom_without_chimerax()
  exercise_load_probe_dots()
  exercise_load_probe_dots_fails()
  exercise_load_probe_dots_nothing_to_do()
  exercise_load_probe_dots_no_answer()
  exercise_busy_while_loading()
  exercise_open_file_names()
  exercise_rest_timeouts()
  exercise_close_busy_viewer()
  exercise_is_running()
  print(format_cpu_times())
  print("OK")
