#include <iostream>
#include "time_logger.h"

TimeLogger::TimeLogger(const char* fn_name):m_fn_name(fn_name)
{
    m_start = clock();
}

TimeLogger::~TimeLogger()
{
    clock_t end = clock();
    float elapsed_s = (float(end - m_start))/CLOCKS_PER_SEC;

    if( m_cumulative_time.count(m_fn_name) == 0 )
      m_cumulative_time[m_fn_name] = 0.0;

    m_cumulative_time[m_fn_name] += elapsed_s;

    std::cout << "TimeLogger " << m_fn_name << " elapsed: " << elapsed_s << "s; Cumulative: " << m_cumulative_time[m_fn_name] << "\n";
}

std::map<std::string,float> TimeLogger::m_cumulative_time;
