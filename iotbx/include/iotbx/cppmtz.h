#ifndef CPPMTZ_H
#define CPPMTZ_H

#include <iostream>
#include <exception>
#include <string>

#include <scitbx/array_family/shared.h>
#include <cctbx/miller.h>
#include <cctbx/hendrickson_lattman.h>
#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/space_group.h>
#include <sys/stat.h> // a cure for errors in irix compile
#include <cmtzlib.h>
#include <ccp4_array.h>

/* Correction 1:  in mtzdata.h, changed CCP4File to CCP4:CCP4File*/

namespace af = scitbx::af;

namespace iotbx{ namespace mtz {

class Error : public std::exception {
private:
  std::string s;
public:
  inline Error(std::string s):s(s){}
  virtual const char* what() const throw() {return s.c_str();}
  virtual ~Error() throw() {}
};

struct Column {
  CMtz::MTZCOL* p_col;
  Column(CMtz::MTZCOL* c);
  bool ccp4_isnan(int i)
  {
    return CCP4::ccp4_utils_isnan((union float_uint_uchar *)&p_col->ref[i]);
  }
  inline float lookup(const int& i){return p_col->ref[i];}
  std::string label();
  std::string type();
};
inline Column::Column(CMtz::MTZCOL* c):p_col(c){
  if (!p_col) {throw Error("Request for a non-existent column");}}

struct Dataset {
  CMtz::MTZSET* p_set;
  Dataset(CMtz::MTZSET* s);
  std::string dataset_name();
  double wavelength();
  int ncolumns();
  Column getColumn(const int&);
};
inline Dataset::Dataset(CMtz::MTZSET* s):p_set(s){
  if (!p_set) {throw Error("Request for a non-existent dataset");}}

struct Crystal {
  CMtz::MTZXTAL* p_xtal;
  Crystal(CMtz::MTZXTAL* c);
  std::string crystal_name();
  std::string project_name();
  cctbx::uctbx::unit_cell UnitCell();
  int ndatasets();
  Dataset getDataset(const int&);
};
inline Crystal::Crystal(CMtz::MTZXTAL* c):p_xtal(c){
  if (!p_xtal) {throw Error("Request for a non-existent crystal");}}

struct observation_arrays
{
  observation_arrays() {}

  observation_arrays(std::size_t size)
  {
    indices.reserve(size);
    data.reserve(size);
    sigmas.reserve(size);
  }

  af::shared<cctbx::miller::index<> > indices;
  af::shared<double> data;
  af::shared<double> sigmas;
};

class Mtz {
private:
  CMtz::MTZ* mtz;
public:
  Mtz(std::string s);
  ~Mtz();

  // Information identified with whole mtz
  std::string title();
  std::string SpaceGroup();
  cctbx::sgtbx::space_group getSgtbxSpaceGroup();
  int nsym() const;
  int& size();
  int& ncrystals() const;
  af::shared<std::string> columns() const;
  af::shared<std::string> history();
  Column getColumn(std::string s);
  Crystal getCrystal(const int&) const;
  Crystal columnToCrystal(std::string) const;

  void printHeader(int);
  void printHeaderAdv(int);

  // Information indexed by crystal number
  cctbx::uctbx::unit_cell UnitCell(const int& xtal);
  int ndatasets(const int& xtal);

  // Information indexed by crystal and dataset
  int& ncolumns(const int& xtal, const int& set);

  af::double2
  max_min_resolution();

  af::shared<cctbx::miller::index<> > MIx();

  af::shared<cctbx::miller::index<> >
  valid_indices(
    std::string const& column_label);

  af::shared<double>
  valid_values(
    std::string const& column_label);

  af::shared<cctbx::miller::index<> >
  valid_indices_anomalous(
    std::string const& column_label_plus,
    std::string const& column_label_minus);

  af::shared<double>
  valid_values_anomalous(
    std::string const& column_label_plus,
    std::string const& column_label_minus);

  observation_arrays
  valid_delta_anomalous(
    std::string const& column_label_f_data,
    std::string const& column_label_f_sigmas,
    std::string const& column_label_d_data,
    std::string const& column_label_d_sigmas);

  af::shared<std::complex<double> >
  valid_complex(
    std::string const& column_label_ampl,
    std::string const& column_label_phi);

  af::shared<std::complex<double> >
  valid_complex_anomalous(
    std::string const& column_label_ampl_plus,
    std::string const& column_label_phi_plus,
    std::string const& column_label_ampl_minus,
    std::string const& column_label_phi_minus);

  af::shared<cctbx::hendrickson_lattman<> >
  valid_hl(
    std::string const& column_label_a,
    std::string const& column_label_b,
    std::string const& column_label_c,
    std::string const& column_label_d);

  af::shared<int>
  valid_integers(
    std::string const& column_label);
};
}} //namespaces
#endif /* CPPMTZ_H*/
