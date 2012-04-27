#ifndef ERSATZ_NEXUZ_H
#define ERSATZ_NEXUZ_H

/*
 * ersatz_nexus
 *
 * A c++ class to wrap up HDF5 calls to create a simple NeXus-compatible HDF5
 * file.
 *
 */

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <hdf5.h>

using namespace std;

class ersatz_nexus
{
private:
  hid_t fid;
public:
  ersatz_nexus(const string & name);
  ~ersatz_nexus();
  hid_t create_group(const string & name);
  void create_attribute(const hid_t gid,
                        const string & name);
  void add_int_data(const hid_t gid,
                    const string & name,
                    const int rank,
                    const hsize_t * dimensions,
                    const int * data);
  void add_double_data(const hid_t gid,
                       const string & name,
                       const int rank,
                       const hsize_t * dimensions,
                       const double * data,
                       const bool signal);
};

string ersatz_nexus_file_name(const int index);
string group_name(const char * prefix,
                  const int & h,
                  const int & k,
                  const int & l);

#endif
