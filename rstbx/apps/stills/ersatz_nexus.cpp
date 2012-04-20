#include "ersatz_nexus.h"

ersatz_nexus::ersatz_nexus(const string & name)
{
  hid_t fapl;

  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
  fid = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  H5Pclose(fapl);
}

ersatz_nexus::~ersatz_nexus()
{
  H5Fclose(fid);
  fid = 0;
}

hid_t ersatz_nexus::create_group(const string & name)
{
  hid_t gid;

  gid = H5Gcreate(fid, (char *) name.c_str(), H5P_DEFAULT, H5P_DEFAULT, 
		  H5P_DEFAULT);

  return gid;
}

void ersatz_nexus::create_attribute(const hid_t gid,
                                    const string & name)
{
  hid_t atts, atttype, attid;

  atts = H5Screate(H5S_SCALAR);
  atttype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atttype, strlen(name.c_str()));
  attid = H5Acreate(gid, "NX_class", atttype, atts, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attid, atttype, (char *) name.c_str());
  H5Sclose(atts);
  H5Tclose(atttype);
  H5Aclose(attid);
}

void ersatz_nexus::add_int_data(const hid_t gid,
                                const string & name,
                                const int rank,
                                const hsize_t * dimensions,
                                const int * data)
{
  hid_t datatype, dataspace, dataprop, dataid;
  hid_t atts, atttype, attid;
  int signal_flag;

  dataspace = H5Screate_simple(rank, dimensions, dimensions);
  datatype = H5Tcopy(H5T_NATIVE_INT);
  dataprop = H5Pcreate(H5P_DATASET_CREATE);
  dataid = H5Dcreate(gid, (char *) name.c_str(), datatype, dataspace,
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataid, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Pclose(dataprop);

  /* add the signal flag as value 1, no idea what this does */

  atts = H5Screate(H5S_SCALAR);
  atttype = H5Tcopy(H5T_NATIVE_INT);
  H5Tset_size(atttype, 1);
  attid = H5Acreate(dataid, "signal", atttype, atts, H5P_DEFAULT, H5P_DEFAULT);
  signal_flag = 1;
  H5Awrite(attid, atttype, & signal_flag);
  H5Sclose(atts);
  H5Tclose(atttype);
  H5Aclose(attid);
  H5Dclose(dataid);
}

void ersatz_nexus::add_double_data(const hid_t gid,
                                   const string & name,
                                   const int rank,
                                   const hsize_t * dimensions,
                                   const double * data)
{
  hid_t datatype, dataspace, dataprop, dataid;
  hid_t atts, atttype, attid;
  int signal_flag;

  dataspace = H5Screate_simple(rank, dimensions, dimensions);
  datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
  dataprop = H5Pcreate(H5P_DATASET_CREATE);
  dataid = H5Dcreate(gid, (char *) name.c_str(), datatype, dataspace,
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataid, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Pclose(dataprop);

  /* add the signal flag as value 1, no idea what this does */

  atts = H5Screate(H5S_SCALAR);
  atttype = H5Tcopy(H5T_NATIVE_INT);
  H5Tset_size(atttype, 1);
  attid = H5Acreate(dataid, "signal", atttype, atts, H5P_DEFAULT, H5P_DEFAULT);
  signal_flag = 1;
  H5Awrite(attid, atttype, & signal_flag);
  H5Sclose(atts);
  H5Tclose(atttype);
  H5Aclose(attid);
  H5Dclose(dataid);
}

string ersatz_nexus_file_name(const int index)
{
    ostringstream filename;
    filename << "ersatz_" << setfill('0') << setw(5) << index << ".nxs";
    return filename.rdbuf()->str();
}

string group_name(const char * prefix,
                  const int & h,
                  const int & k,
                  const int & l)
{
    ostringstream group_name;
    group_name << prefix <<
      setw(4) << h <<
      setw(4) << k <<
      setw(4) << l;
    return group_name.rdbuf()->str();
}

int main_test(int argc,
              char ** argv)
{
  int * data;
  int i, j, k, rank, nx, ny;
  hid_t gid;
  hsize_t dim[2];

  char data_name[100];

  cout << ersatz_nexus_file_name(192) << endl;
  cout << group_name("group", 1, -2, 3) << endl;

  /* really crappy error trapping */
  if (argc < 2) return 1;

  ersatz_nexus en(argv[1]);

  nx = 100;
  ny = 100;

  data = (int *) malloc (sizeof(int) * nx * ny);

  dim[0] = ny; dim[1] = nx;
  rank = 2;

  gid = en.create_group("profile");
  en.create_attribute(gid, "NXentry");

  for (k = 0; k < 100; k ++) {
    for (i = 0; i < ny; i++) {
      for (j = 0; j < nx; j++) {
        data[i * nx + j] = (i + j) % (k + 1);
      }
    }

    sprintf(data_name, "/profile/peak%03d", k);

    gid = en.create_group(data_name);
    en.create_attribute(gid, "NXdata");

    sprintf(data_name, "data%03d", k);

    en.add_int_data(gid, data_name, rank, dim, data);
  }

  free(data);

  return 0;
}
