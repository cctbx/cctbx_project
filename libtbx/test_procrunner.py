from __future__ import absolute_import, division, print_function

import copy
import libtbx.procrunner
import mock
import os
import pytest

@pytest.mark.skipif(libtbx.procrunner.dummy, reason='procrunner class set to dummy mode')
@mock.patch('libtbx.procrunner._NonBlockingStreamReader')
@mock.patch('libtbx.procrunner.time')
@mock.patch('libtbx.procrunner.subprocess')
@mock.patch('libtbx.procrunner.Pipe')
def test_run_command_aborts_after_timeout(mock_pipe, mock_subprocess, mock_time, mock_streamreader):
  mock_pipe.return_value = mock.Mock(), mock.Mock()
  mock_process = mock.Mock()
  mock_process.returncode = None
  mock_subprocess.Popen.return_value = mock_process
  task = ['___']

  with pytest.raises(RuntimeError):
    libtbx.procrunner.run_process(task, -1, False)

  assert mock_subprocess.Popen.called
  assert mock_process.terminate.called
  assert mock_process.kill.called


@pytest.mark.skipif(libtbx.procrunner.dummy, reason='procrunner class set to dummy mode')
@mock.patch('libtbx.procrunner._NonBlockingStreamReader')
@mock.patch('libtbx.procrunner.subprocess')
def test_run_command_runs_command_and_directs_pipelines(mock_subprocess, mock_streamreader):
  (mock_stdout, mock_stderr) = (mock.Mock(), mock.Mock())
  mock_stdout.get_output.return_value = mock.sentinel.proc_stdout
  mock_stderr.get_output.return_value = mock.sentinel.proc_stderr
  (stream_stdout, stream_stderr) = (mock.sentinel.stdout, mock.sentinel.stderr)
  mock_process = mock.Mock()
  mock_process.stdout = stream_stdout
  mock_process.stderr = stream_stderr
  mock_process.returncode = 99
  command = ['___']
  def streamreader_processing(*args, **kwargs):
    return {(stream_stdout,): mock_stdout, (stream_stderr,): mock_stderr}[args]
  mock_streamreader.side_effect = streamreader_processing
  mock_subprocess.Popen.return_value = mock_process

  expected = {
    'stderr': mock.sentinel.proc_stderr,
    'stdout': mock.sentinel.proc_stdout,
    'exitcode': mock_process.returncode,
    'command': command,
    'runtime': mock.ANY,
    'timeout': False,
    'time_start': mock.ANY,
    'time_end': mock.ANY
  }

  actual = libtbx.procrunner.run_process(command, 0.5, False,
               callback_stdout=mock.sentinel.callback_stdout, callback_stderr=mock.sentinel.callback_stderr)

  assert mock_subprocess.Popen.called
  assert mock_subprocess.Popen.call_args[1]['env'] == os.environ
  mock_streamreader.assert_has_calls([mock.call(stream_stdout, output=mock.ANY, debug=mock.ANY, notify=mock.ANY, callback=mock.sentinel.callback_stdout),
                                      mock.call(stream_stderr, output=mock.ANY, debug=mock.ANY, notify=mock.ANY, callback=mock.sentinel.callback_stderr)],
                                     any_order=True)
  assert not mock_process.terminate.called
  assert not mock_process.kill.called
  assert actual == expected


@pytest.mark.skipif(libtbx.procrunner.dummy, reason='procrunner class set to dummy mode')
@mock.patch('libtbx.procrunner.subprocess')
def test_default_process_environment_is_parent_environment(mock_subprocess):
  mock_subprocess.Popen.side_effect = NotImplementedError() # cut calls short
  with pytest.raises(NotImplementedError):
    libtbx.procrunner.run_process(mock.Mock(), -1, False)
  assert mock_subprocess.Popen.call_args[1]['env'] == os.environ


@pytest.mark.skipif(libtbx.procrunner.dummy, reason='procrunner class set to dummy mode')
@mock.patch('libtbx.procrunner.subprocess')
def test_pass_custom_environment_to_process(mock_subprocess):
  mock_subprocess.Popen.side_effect = NotImplementedError() # cut calls short
  mock_env = { 'key': mock.sentinel.key }
  # Pass an environment dictionary
  with pytest.raises(NotImplementedError):
    libtbx.procrunner.run_process(mock.Mock(), -1, False, environment=copy.copy(mock_env))
  assert mock_subprocess.Popen.call_args[1]['env'] == mock_env


@pytest.mark.skipif(libtbx.procrunner.dummy, reason='procrunner class set to dummy mode')
@mock.patch('libtbx.procrunner.subprocess')
def test_pass_custom_environment_to_process_and_add_another_value(mock_subprocess):
  mock_subprocess.Popen.side_effect = NotImplementedError() # cut calls short
  mock_env1 = { 'keyA': mock.sentinel.keyA }
  mock_env2 = { 'keyB': str(mock.sentinel.keyB) }
  # Pass an environment dictionary
  with pytest.raises(NotImplementedError):
    libtbx.procrunner.run_process(mock.Mock(), -1, False, environment=copy.copy(mock_env1), environment_override=copy.copy(mock_env2))
  mock_env_sum = copy.copy(mock_env1)
  mock_env_sum.update(mock_env2)
  assert mock_subprocess.Popen.call_args[1]['env'] == mock_env_sum


@pytest.mark.skipif(libtbx.procrunner.dummy, reason='procrunner class set to dummy mode')
@mock.patch('libtbx.procrunner.subprocess')
def test_use_default_process_environment_and_add_another_value(mock_subprocess):
  mock_subprocess.Popen.side_effect = NotImplementedError() # cut calls short
  mock_env2 = { 'keyB': str(mock.sentinel.keyB) }
  with pytest.raises(NotImplementedError):
    libtbx.procrunner.run_process(mock.Mock(), -1, False, environment_override=copy.copy(mock_env2))
  random_environment_variable = list(os.environ)[0]
  random_environment_value = os.getenv(random_environment_variable)
  assert random_environment_variable and random_environment_variable != mock_env2.keys()[0]
  assert mock_subprocess.Popen.call_args[1]['env'][mock_env2.keys()[0]] == mock_env2.values()[0]
  assert mock_subprocess.Popen.call_args[1]['env'][random_environment_variable] == os.getenv(random_environment_variable)


@pytest.mark.skipif(libtbx.procrunner.dummy, reason='procrunner class set to dummy mode')
@mock.patch('libtbx.procrunner.subprocess')
def test_use_default_process_environment_and_override_a_value(mock_subprocess):
  mock_subprocess.Popen.side_effect = NotImplementedError() # cut calls short
  random_environment_variable = list(os.environ)[0]
  random_environment_value = os.getenv(random_environment_variable)
  with pytest.raises(NotImplementedError):
    libtbx.procrunner.run_process(mock.Mock(), -1, False, environment_override={ random_environment_variable: 'X' + random_environment_value })
  assert mock_subprocess.Popen.call_args[1]['env'][random_environment_variable] == 'X' + random_environment_value


@mock.patch('libtbx.procrunner.select')
def test_nonblockingstreamreader_can_read(mock_select):
  import time
  class _stream:
    def __init__(self):
      self.data = ""
      self.closed = False
    def write(self, string):
      self.data = self.data + string
    def read(self, n):
      if self.closed:
        return ""
      if self.data == "":
        time.sleep(0.01)
        return ""
      if (len(self.data) < n):
        data = self.data
        self.data = ""
      else:
        data = self.data[:n]
        self.data = self.data[n:]
      return data
    def close(self):
      self.closed=True
  teststream = _stream()

  def select_replacement(rlist, wlist, xlist, timeout):
    assert teststream in rlist
    if teststream.closed:
      return ([teststream], [], [])
    if teststream.data == "":
      return ([], [], [])
    return ([teststream], [], [])
  mock_select.select = select_replacement

  streamreader = libtbx.procrunner._NonBlockingStreamReader(teststream, output=False)
  assert not streamreader.has_finished()
  time.sleep(0.1)
  testdata = "abc\n" * 1024
  teststream.write(testdata)
  time.sleep(0.2)
  teststream.close()
  time.sleep(0.1)

  assert streamreader.has_finished()
  output = streamreader.get_output()
  assert len(output) == len(testdata)
  assert output == testdata


def test_lineaggregator_aggregates_data():
  callback = mock.Mock()
  aggregator = libtbx.procrunner._LineAggregator(callback=callback)

  aggregator.add('some')
  aggregator.add('string')
  callback.assert_not_called()
  aggregator.add("\n")
  callback.assert_called_once_with('somestring')
  callback.reset_mock()
  aggregator.add('more')
  aggregator.add('stuff')
  callback.assert_not_called()
  aggregator.flush()
  callback.assert_called_once_with('morestuff')
