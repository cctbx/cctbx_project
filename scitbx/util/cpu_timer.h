#ifndef SCITBX_CPU_TIMER_H
#define SCITBX_CPU_TIMER_H

#include <string>
#include <sstream>

#include <boost/date_time/posix_time/posix_time.hpp>

namespace scitbx { namespace util {

class cpu_timer
{
public:

  cpu_timer() : start_(boost::posix_time::microsec_clock::local_time()) {}

  void format(std::ostream &stream)
  {
    stream << this->elapsed() << "s wall";
  }

  std::string format()
  {
    std::ostringstream stream;
    this->format(stream);
    return stream.str();
  }

  void restart()
  {
    start_ = boost::posix_time::microsec_clock::local_time();
  }

  double elapsed() const
  {
    boost::posix_time::ptime t=boost::posix_time::microsec_clock::local_time();
    boost::posix_time::time_duration tdif = t-start_;
    return tdif.total_microseconds()*1.e-6;
  }

private:
  boost::posix_time::ptime start_;
};

}}
#endif
