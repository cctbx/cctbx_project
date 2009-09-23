#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost_adaptbx/python_file_stream.h>
#include <boost/timer.hpp>
#include <fstream>

namespace boost_adaptbx { namespace file_conversion {

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

  std::string test_read(python_file_buffer const &input,
                        std::string const &what)
  {
    istream is(&input);
    std::string result = actual_read_test(is, what);
    is.sync();
    return result;
  }

  std::string test_write(python_file_buffer const &output,
                         std::string const &what)
  {
    ostream os(&output);
    std::string result = actual_write_test(os, what);
    os.flush();
    return result;
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

  void time_read(char const *path, python_file_buffer const &input) {
    istream is(&input);
    boost::timer t;
    work_for_time_read(is);
    double py_t = t.elapsed();
    std::ifstream std_is(path);
    t.restart();
    work_for_time_read(std_is);
    double py_cpp = t.elapsed();
    std::cout << "- Reading -\nPython bridge:" << py_t;
    std::cout << "\nPure C++:" << py_cpp;
    std::cout << "\noverhead: " << (py_t - py_cpp)/py_t*100 << " %";
    std::cout << "\n\n";
  }

  void time_write(char const *path, python_file_buffer const &output) {
    ostream os(&output);
    boost::timer t;
    work_for_time_write(os);
    double py_t = t.elapsed();
    std::ofstream std_os(path);
    t.restart();
    work_for_time_write(std_os);
    double py_cpp = t.elapsed();
    std::cout << "- Writing -\nPython bridge:" << py_t;
    std::cout << "\nPure C++:" << py_cpp;
    std::cout << "\noverhead: " << (py_t - py_cpp)/py_t*100 << " %";
    std::cout << "\n\n";
  }

}} // boost_adaptbx::file_conversion

BOOST_PYTHON_MODULE(python_file_test_ext)
{
  using namespace boost::python;
  using namespace boost_adaptbx::file_conversion;
  def("test_read", test_read);
  def("test_write", test_write);
  def("time_read", time_read);
  def("time_write", time_write);
}
