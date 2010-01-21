
import libtbx.phil
from libtbx import easy_pickle, thread_utils
import sys, os

process_master_phil = libtbx.phil.parse("""
run_file = None
  .type = path
prefix = None
  .type = str
tmp_dir = None
  .type = path
debug = False
  .type = bool
""")

class detached_process_server (object) :
  def __init__ (self, run_file, status_file_prefix="", tmp_dir=None,
      d) :
    target = easy_pickle.load(run_file)
    assert isinstance(target, detached_process_driver)
    if tmp_dir is None :
      tmp_dir = target.output_dir
    else :
      assert os.path.isdir(tmp_dir)
    prefix = os.path.join(tmp_dir, status_file_prefix)
    self.libtbx_process = thread_utils.process_with_gui_callbacks(
      target=target,
      callback_stdout=self.callback_stdout,
      callback_final=self.callback_final,
      callback_err=self.callback_err,
      callback_abort=self.callback_abort,
      callback_other=self.callback_other,
      buffer_stdout=True)
    start_file = os.path.join(self.tmp_dir, self.prefix + ".libtbx_start")
    f = open(start_file, "w")
    f.write("1")
    f.close()
    prefix
    self.stdout_file = prefix + ".libtbx_stdout"
    self._stdout = open(self.stdout_file, "w")
    self.error_file =  prefix + ".libtbx_error"
    self.abort_file =  prefix + ".libtbx_abort"
    self.result_file = prefix + ".libtbx_result"
    self.status_file = prefix + ".libtbx_status"
    self.other_file =  prefix + ".libtbx_other"
    self._accumulated_callbacks = []
    self.libtbx_process.start()

  def callback_stdout (self, data) :
    self._stdout.write(data)
    self._stdout.flush()

  def callback_error (self, error, traceback_info) :
    easy_pickle.dump(self.error_file, (error, traceback_info))

  def callback_abort (self) :
    easy_pickle.dump(self.abort_file, True)

  def callback_final (self, result) :
    easy_pickle.dump(self.result_file, result)
    self._stdout.flush()
    self._stdout.close()

  def callback_other (self, status) :
    if not status.cached :
      return
    if status.accumulate :
      self._accumulated_callbacks.append(status)
      easy_pickle.dump(self.status_file, self._accumulated_callbacks)
    else :
      easy_pickle.dump(self.other_file, status)

class detached_process_driver (object) :
  def __init__ (self, output_dir, target) :
    adopt_init_args(self, locals())

  def __call__ (self, args, kwds, child_conn) :
    os.chdir(self.output_dir)
    import libtbx.callbacks
    callback = libtbx.callbacks.piped_callback(child_conn)
    libtbx.call_back.register_handler(callback)
    result = self.target(args, kwds, child_conn)
    return result

def run (args) :
  user_phil = []
  for arg in args :
    if os.path.isfile(arg) :
      file_name = os.path.abspath(arg)
      user_phil.append(libtbx.phil.parse("run_file = \"%s\"" % file_name))
    else :
      try :
        arg_phil = libtbx.phil.parse(arg)
      except RuntimeError :
        pass
      else :
        user_phil.append(arg_phil)
  params = master_phil.fetch(sources=user_phil).extract()
  detached_process_server(run_file=params.run_file,
    status_file_prefix=params.prefix,
    tmp_dir=params.tmp_dir)

#---end
