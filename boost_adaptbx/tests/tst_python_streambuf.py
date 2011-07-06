import boost.python
from boost.python import streambuf, ostream
ext = boost.python.import_ext("boost_adaptbx_python_streambuf_test_ext")
import StringIO
import cStringIO
from libtbx.test_utils import open_tmp_file, Exception_expected
from libtbx.option_parser import option_parser
import libtbx.object_oriented_patterns as oop
import os


class file_object_debug_proxy(oop.proxy):

  def __init__(self, *args, **kwds):
    super(file_object_debug_proxy, self).__init__(*args, **kwds)
    self.seek_call_log = []
    self.write_call_log = []

  def seek(self, off, mode):
    self.seek_call_log.append((off, mode))
    self.subject.seek(off, mode)

  def write(self, s):
    self.write_call_log.append(s)
    self.subject.write(s)


class io_test_case(object):
  phrase = "Coding should be fun"
  #         01234567890123456789

  def run(self):
    m = streambuf.default_buffer_size
    for n in xrange(50, 0, -1):
      streambuf.default_buffer_size = n
      self.exercise_read_failure()
      self.exercise_write_failure()
      self.exercise_read()
      self.exercise_write()
      self.exercise_seek_and_read()
      self.exercise_partial_read()
      self.exercise_write_and_seek()
    streambuf.default_buffer_size = m

  def exercise_read(self):
    self.create_file_object(mode='r')
    words = ext.test_read(streambuf(self.file_object), "read")
    assert words == "Coding, should, be, fun, [ fail, eof ]"
    self.file_object.close()

  def exercise_partial_read(self):
    self.create_file_object(mode='r')
    words = ext.test_read(streambuf(self.file_object), "partial read")
    assert words == "Coding, should, "
    trailing = self.file_object.read()
    assert  trailing == " be fun"
    self.file_object.close()

  def exercise_read_failure(self):
    self.create_file_object(mode='r')
    self.file_object.close()
    try:
      ext.test_read(streambuf(self.file_object), "read")
    except ValueError:
      pass
    else:
      raise Exception_expected
    self.file_object.close()

  def exercise_write(self):
    self.create_file_object(mode='w')
    report = ext.test_write(ostream(self.file_object), "write")
    assert report == ''
    assert self.file_content() == "2 times 1.6 equals 3.2"
    self.file_object.close()

  def exercise_seek_and_read(self):
    self.create_instrumented_file_object(mode='r')
    words = ext.test_read(streambuf(self.file_object), "read and seek")
    assert words == "should, should, uld, ding, fun, [ eof ]"
    n = streambuf.default_buffer_size
    soughts = self.file_object.seek_call_log
    # stringent tests carefully crafted to make sure the seek-in-buffer
    # optimisation works as expected
    # C.f. the comment in the C++ function actual_write_test
    assert soughts[-1] == (-4,2)
    soughts = soughts[:-1]
    if n >= 14:
      assert soughts == []
    else:
      assert soughts[0] == (6,0)
      soughts = soughts[1:]
      if 8 <= n <= 13:
        assert len(soughts) == 1 and self.only_seek_cur(soughts)
      elif n == 7:
        assert len(soughts) == 2 and self.only_seek_cur(soughts)
      else:
        assert soughts[0] == (6,0)
        soughts = soughts[1:]
        assert self.only_seek_cur(soughts)
        if n == 4: assert len(soughts) == 1
        else: assert len(soughts) == 2
    self.file_object.close()

  def exercise_write_and_seek(self):
    self.create_instrumented_file_object(mode='w')
    report = ext.test_write(ostream(self.file_object), "write and seek (cur)")
    assert report == ''
    expected = '1000 times 1000 equals 1000000'
    assert self.file_content() == expected
    assert self.file_object.tell() == 9
    if streambuf.default_buffer_size >= 30:
      assert self.file_object.write_call_log == [ expected ]
    self.file_object.close()

  def only_seek_cur(cls, seek_calls):
    return [ whence for off, whence in seek_calls if whence != 1 ] == []

  def create_instrumented_file_object(self, mode):
    self.create_file_object(mode)
    self.file_object = file_object_debug_proxy(self.file_object)


class stringio_test_case(io_test_case):

  stringio_type = StringIO.StringIO

  def exercise_write_failure(self):
    pass

  def create_file_object(self, mode):
    if mode == 'r':
      self.file_object = self.stringio_type(self.phrase)
    elif mode == 'w':
      self.file_object = self.stringio_type()
    else:
      raise NotImplementedError("Internal error in the test code")

  def file_content(self):
    return self.file_object.getvalue()


class cstringio_test_case(stringio_test_case):

  stringio_type = cStringIO.StringIO

  def exercise_write_failure(self):
    self.create_file_object(mode='r')
    try:
      ext.test_write(ostream(self.file_object), "write")
    except ValueError, err:
      # That Python file object has no 'write' attribute
      assert str(err).find("write") > -1
    except RuntimeError, err:
      # Redhat 8.0: basic_ios::clear(iostate) caused exception
      assert str(err).find("clear") > -1
    else:
      raise Exception_expected


class mere_file_test_case(io_test_case):

  def exercise_write_failure(self):
    import platform
    if (platform.platform().find("redhat-8.0") >= 0):
      return # avoid Abort
    self.create_file_object(mode='r')
    try:
      ext.test_write(streambuf(self.file_object), "write")
    except IOError, err:
      pass
    else:
      raise Exception_expected
    self.file_object.close()

  def create_file_object(self, mode):
    f = open_tmp_file()
    if mode.find('r') > -1:
      f.write(self.phrase)
    f.close()
    self.file_object = open(f.name, mode)

  def file_content(self):
    i = self.file_object.tell()
    self.file_object.flush()
    result = open(self.file_object.name).read()
    self.file_object.seek(i, 0)
    return result


def time_it(path, buffer_size):
  if (buffer_size is None):
    buffer_size = streambuf.default_buffer_size
  print "Buffer is %i bytes" % buffer_size
  path = os.path.expanduser(path)
  input = open(path, 'r')
  inp_buf = streambuf(python_file_obj=input, buffer_size=buffer_size)
  ext.time_read(input.name, inp_buf)
  output = open_tmp_file()
  out_buf = streambuf(python_file_obj=output, buffer_size=buffer_size)
  ext.time_write(output.name, out_buf)

def run(args):
  options = (option_parser()
              .option(None, '--time_on_file',
                      metavar="PATH",
                      help="time reading and writing."
                           "The file to read shall be a hkl file, i.e "
                           "each line has a format like "
                           "'int int int double double'. "
                           "The end shall be marked by 0 0 0")
              .option(None, '--buffer_size', type='int',
                      metavar="INT")
              ).process(args).options
  for i_trial in xrange(3):
    stringio_test_case().run()
  for i_trial in xrange(3):
    cstringio_test_case().run()
  for i_trial in xrange(3):
    mere_file_test_case().run()
  if options.time_on_file:
    time_it(options.time_on_file, options.buffer_size)

  print 'OK'

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
