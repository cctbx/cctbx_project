#ifndef CPPMTZ_H
#define CPPMTZ_H

#include <iostream>
#include <exception>
#include <string>

#include <sys/stat.h> // a cure for errors in irix compile
#include "cmtzlib.h"
#include <scitbx/array_family/shared.h>
#include <cctbx/miller.h>

/* Correction 1:  in mtzdata.h, changed CCP4File to CCP4:CCP4File*/

namespace af = scitbx::af;

namespace iotbx{ namespace mtz {

struct Foo {
  af::shared<std::string> value();
};

class Error : public std::exception {
private:
  std::string s;
public:
  inline Error(std::string s):s(s){}
  virtual const char* what() const throw();
  virtual ~Error() throw();
};
inline const char* Error::what() const throw() {return s.c_str();}
inline Error::~Error() throw() {}

struct Column {
  CMtz::MTZCOL* p_col;
  Column(CMtz::MTZCOL* c);
  inline double lookup(const int& i){return p_col->ref[i];}
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
  af::shared<double> UnitCell();
  int ndatasets();
  Dataset getDataset(const int&);
};
inline Crystal::Crystal(CMtz::MTZXTAL* c):p_xtal(c){
  if (!p_xtal) {throw Error("Request for a non-existent crystal");}}

class Mtz {
private:
  CMtz::MTZ* mtz;
public:
  Mtz(std::string s);
  ~Mtz();

  // Information identified with whole mtz
  std::string title();
  std::string SpaceGroup();
  int& size();
  int& ncrystals();
  af::shared<std::string> columns();
  af::shared<std::string> history();
  Column getColumn(std::string s);
  af::shared<double> getShared(std::string s);
  Crystal getCrystal(const int&);

  void printHeader(int);
  void printHeaderAdv(int);

  // Information indexed by crystal number
  af::shared<double> UnitCell(const int& xtal);
  int ndatasets(const int& xtal);

  // Information indexed by crystal and dataset
  int& ncolumns(const int& xtal, const int& set);

  af::shared< cctbx::miller::index<> > MIx();
};
}} //namespaces
#endif /* CPPMTZ_H*/
