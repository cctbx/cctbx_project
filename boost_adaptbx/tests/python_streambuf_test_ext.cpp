#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <boost_adaptbx/python_streambuf.h>
#include <boost/timer/timer.hpp>
#include <fstream>

namespace boost_adaptbx { namespace python { namespace {

  template <class StreamType>
  void append_status(StreamType const &s, std::string &result) {
    if (!s.good()) result += "[ ";
    if (s.bad()) result += "bad, ";
    if (s.fail()) result += "fail, ";
    if (s.eof()) result += "eof";
    if (!s.good()) result += " ]";
  }

  std::string actual_read_test(std::istream &is, std::string const &what) {
    // Coding should be fun
    // 012345678901234567890
    std::string word, result;
    if (what == "read") {
      while (is >> word) {
        result += word + ", ";
      }
    }
    else if (what == "read and seek") {
      is.seekg(6);
      is >> word; result += word + ", "; // should
      is.seekg(6, std::ios_base::beg);
      is >> word; result += word + ", "; // should
      is.seekg(-3, std::ios_base::cur);
      is >> word; result += word + ", "; // uld
      is.seekg(-11, std::ios_base::cur);
      is >> word; result += word + ", "; // ding
      is.seekg(-4, std::ios_base::end);
      is >> word; result += word + ", "; // fun
    }
    else if (what == "partial read") {
      is >> word; result += word + ", ";
      is >> word; result += word + ", ";
    }
    append_status(is, result);
    return result;
  }

  std::string actual_write_test(std::ostream &os, std::string const &what) {
    std::string result;
    if (what == "write") {
      os << 2 << " times " << 1.6 << " equals " << 3.2;
    }
    else if (what.find("write and seek") == 0) {
      os << 1000 << " timEs " << 5555 << " equalS " << 1000000;
      // 1000 timEs 5555 equalS 1700000
      // 0123456789012345678901234567890
      if (what.find("(cur)") != std::string::npos) {
        os.seekp(-19, std::ios_base::cur);
        os << 1000;
        os.seekp(6, std::ios_base::cur);
        os << "s";
        os.seekp(-14, std::ios_base::cur);
        os << "e";
      }
    }
    append_status(os, result);
    return result;
  }

  std::string test_read(streambuf& input,
                        std::string const &what)
  {
    streambuf::istream is(input);
    return actual_read_test(is, what);
  }

  std::string test_write(streambuf& output,
                         std::string const &what)
  {
    streambuf::ostream os(output);
    return actual_write_test(os, what);
  }

  double work_for_time_read(std::istream &is) {
    int h,k,l,foo;
    double f,s;
    double result = 0;
    while (is >> h >> k >> l) {
      if (h == 0 && k == 0 && l == 0) break;
      is >> f >> s >> foo;
      result += (h + k + l)*s*f;
    }
    return result;
  }

  void work_for_time_write(std::ostream &os) {
    int i=1, j=1;
    for (int n=0; n < 1000000; ++n) {
      os << j << " ";
      int j_ = (i + j) % 65536;
      i = j;
      j = j_;
    }
  }

  void time_read(char const *path, streambuf& input) {
    boost::timer::auto_cpu_timer t;
    streambuf::istream is(input);
    work_for_time_read(is);
    double py_t = t.elapsed().wall;
    std::ifstream std_is(path);
    t.start();
    work_for_time_read(std_is);
    double py_cpp = t.elapsed().wall;
    std::cout << "- Reading -\nPython adaptor: " << py_t;
    std::cout << "\nPure C++: " << py_cpp;
    if (py_t > py_cpp) {
      std::cout << "\noverhead: " << (py_t - py_cpp)/py_t*100 << " %";
    }
    std::cout << "\n\n";
  }

  void time_write(char const *path, streambuf& output) {
    boost::timer::auto_cpu_timer t;
    streambuf::ostream os(output);
    work_for_time_write(os);
    double py_t = t.elapsed().wall;
    std::ofstream std_os(path);
    t.start();
    work_for_time_write(std_os);
    double py_cpp = t.elapsed().wall;
    std::cout << "- Writing -\nPython adaptor: " << py_t;
    std::cout << "\nPure C++: " << py_cpp;
    if (py_t > py_cpp) {
      std::cout << "\noverhead: " << (py_t - py_cpp)/py_t*100 << " %";
    }
    std::cout << "\n\n";
  }

  void
  call_with_stderr_stdout_do_nothing(
    std::ostream& err_stream,
    std::ostream& out_stream)
  {}

  void
  wrap_all()
  {
    using namespace boost::python;
    def("test_read", test_read);
    def("test_write", actual_write_test);
    def("test_write", test_write);
    def("time_read", time_read);
    def("time_write", time_write);
    def("call_with_stderr_stdout_do_nothing",
         call_with_stderr_stdout_do_nothing);
  }

}}} // namespace boost_adaptbx::python::<anonymous>

BOOST_PYTHON_MODULE(boost_adaptbx_python_streambuf_test_ext)
{
  boost_adaptbx::python::wrap_all();
}
