#pragma once

#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <map>
#include <string>

class TimeLogger
{
public:
    TimeLogger(const char* fn_name);
    ~TimeLogger();
private:
    clock_t m_start;
    std::string m_fn_name;
    static std::map<std::string,float> m_cumulative_time;
};
