#include "iotbx/cppmtz.h"
#include <scitbx/constants.h>
#include <scitbx/mat3.h>

namespace iotbx { namespace mtz {

Mtz::Mtz(std::string s):mtz(CMtz::MtzGet(s.c_str(),1)){
    if (!mtz) {throw Error("Mtz read failed"); }
}
Mtz::~Mtz(){
  CMtz::MtzFree(mtz);}

std::string Mtz::title(){return std::string(mtz->title);}
std::string Mtz::SpaceGroup(){return mtz->mtzsymm.spcgrpname;}

cctbx::sgtbx::space_group
Mtz::getSgtbxSpaceGroup(){
  cctbx::sgtbx::space_group sg;
  for  (int i = 0; i < mtz->mtzsymm.nsym; ++i) {
    scitbx::mat3<double> r_double;
    scitbx::vec3<double> t_double;
    for (int p=0;p<3;p++) {
      for (int q=0;q<3;q++) {
        r_double(p,q) = mtz->mtzsymm.sym[i][p][q];
      }
      t_double[p] = mtz->mtzsymm.sym[i][p][3];
    }
    sg.expand_smx(cctbx::sgtbx::rt_mx(r_double,t_double));
  }
  return sg;
}

int Mtz::nsym() const {return mtz->mtzsymm.nsym;}
int& Mtz::size(){return mtz->nref;}
int& Mtz::ncrystals() const {return mtz->nxtal;}
int Mtz::ndatasets(const int& xtal) {return mtz->xtal[xtal]->nset;}
int& Mtz::ncolumns(const int& xtal, const int& set) {return mtz->xtal[xtal]->set[set]->ncol;}

af::shared<std::string> Mtz::history(){
   af::shared<std::string> answer;
   for (int i = 0; i < mtz->histlines; ++i) {
     char* buffer=new char[MTZRECORDLENGTH+1];
     char* src= mtz->hist+ MTZRECORDLENGTH*i;
     char* dst= buffer;
     for (;dst-buffer<MTZRECORDLENGTH;) {*dst++=*src++;}
     buffer[MTZRECORDLENGTH]='\0';
     answer.push_back(std::string(buffer));
     delete [] buffer;
   }
   return answer;
}

af::shared<std::string> Mtz::columns() const {
   af::shared<std::string> answer;
   /* Loop over crystals */
    for (int i = 0; i < mtz->nxtal; ++i) {
   /* Loop over datasets for each crystal */
     for (int j = 0; j < mtz->xtal[i]->nset; ++j) {
   /* Loop over columns for each dataset */
      for (int k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
       if (mtz->xtal[i]->set[j]->col[k]->active) {
         answer.push_back(mtz->xtal[i]->set[j]->col[k]->label);
       }
      }
     }
    }
   return answer;
}

Column Mtz::getColumn(std::string s){
  return Column(CMtz::MtzColLookup(mtz,s.c_str()));
}

Crystal Mtz::getCrystal(const int& xtalid) const{
  return Crystal(CMtz::MtzIxtal(mtz, xtalid));
}

Crystal Mtz::columnToCrystal(std::string s) const{
  af::shared<std::string> cols = this->columns();
  bool foundlabel=false;
  for (int i = 0; i<cols.size(); ++i) {
    if (cols[i]==s) {
      foundlabel=true;
      break;
    }
  }
  if (!foundlabel) throw iotbx::mtz::Error("no such column label");

  for (int i = 0; i < this->ncrystals(); ++i) {
    iotbx::mtz::Crystal cryst = this->getCrystal(i);
    for (int j = 0; j < cryst.ndatasets(); ++j) {
      iotbx::mtz::Dataset data = cryst.getDataset(j);
      for (int k = 0; k < data.ncolumns(); ++k) {
        if (data.getColumn(k).label() == s) {
          return cryst;}
      }
    }
  }
  throw iotbx::mtz::Error("no such column label");//Never get here but make compiler happy
}

af::double2
Mtz::max_min_resolution()
{
  Column H(getColumn("H"));
  Column K(getColumn("K"));
  Column L(getColumn("L"));
  std::vector<cctbx::uctbx::unit_cell> unit_cells;
  for (int i_crystal=0;i_crystal<ncrystals();i_crystal++) {
    unit_cells.push_back(UnitCell(i_crystal));
  }
  double d_star_sq_min = -1;
  double d_star_sq_max = -1;
  for (int i_index=0; i_index<this->size(); i_index++) {
    for (int i_crystal=0;i_crystal<ncrystals();i_crystal++) {
      double d_star_sq = unit_cells[i_crystal].d_star_sq(
        cctbx::miller::index<>(
          H.lookup(i_index),
          K.lookup(i_index),
          L.lookup(i_index)));
      if (d_star_sq_min > d_star_sq || d_star_sq_min < 0) {
          d_star_sq_min = d_star_sq;
      }
      if (d_star_sq_max < d_star_sq) {
          d_star_sq_max = d_star_sq;
      }
    }
  }
  return af::double2(
    d_star_sq_min <= 0 ? -1 : 1/std::sqrt(d_star_sq_min),
    d_star_sq_max <= 0 ? -1 : 1/std::sqrt(d_star_sq_max));
}

af::shared< cctbx::miller::index<> > Mtz::MIx() {
  Column H(getColumn("H"));
  Column K(getColumn("K"));
  Column L(getColumn("L"));
  af::shared<cctbx::miller::index<> > miller(this->size());
  for (int j=0; j<this->size(); j++) {
    miller[j] = cctbx::miller::index<>(H.lookup(j),K.lookup(j),L.lookup(j));
  }
  return miller;
}

void Mtz::printHeader(int detail = 1){
  CMtz::ccp4_lhprt(mtz,detail);
}

void Mtz::printHeaderAdv(int detail = 1){
  CMtz::ccp4_lhprt_adv(mtz,detail);
}

cctbx::uctbx::unit_cell Mtz::UnitCell(const int& xtal){
  af::tiny<double,6> answer;
  for (int j=0; j<6; j++) {answer[j] = mtz->xtal[xtal]->cell[j];}
  return cctbx::uctbx::unit_cell(answer);
}

af::shared<cctbx::miller::index<> >
Mtz::valid_indices(
  std::string const& column_label)
{
  Column h(getColumn("H"));
  Column k(getColumn("K"));
  Column l(getColumn("L"));
  Column v(getColumn(column_label));
  af::shared<cctbx::miller::index<> > result((af::reserve(size())));
  for (int j=0; j<size(); j++) {
    if (!v.ccp4_isnan(j)) {
      result.push_back(cctbx::miller::index<>(
        static_cast<int>(h.lookup(j)),
        static_cast<int>(k.lookup(j)),
        static_cast<int>(l.lookup(j))));
    }
  }
  return result;
}

af::shared<double>
Mtz::valid_values(
  std::string const& column_label)
{
  Column v(getColumn(column_label));
  af::shared<double> result((af::reserve(size())));
  for (int j=0; j<size(); j++) {
    if (!v.ccp4_isnan(j)) result.push_back(v.lookup(j));
  }
  return result;
}

af::shared<cctbx::miller::index<> >
Mtz::valid_indices_anomalous(
  std::string const& column_label_plus,
  std::string const& column_label_minus)
{
  Column h(getColumn("H"));
  Column k(getColumn("K"));
  Column l(getColumn("L"));
  Column v_plus(getColumn(column_label_plus));
  Column v_minus(getColumn(column_label_minus));
  af::shared<cctbx::miller::index<> > result((af::reserve(2*size())));
  for (int j=0; j<size(); j++) {
    if (!v_plus.ccp4_isnan(j)) {
      result.push_back(cctbx::miller::index<>(
        h.lookup(j), k.lookup(j), l.lookup(j)));
    }
    if (!v_minus.ccp4_isnan(j)) {
      result.push_back(cctbx::miller::index<>(
        -h.lookup(j), -k.lookup(j), -l.lookup(j)));
    }
  }
  return result;
}

af::shared<double>
Mtz::valid_values_anomalous(
  std::string const& column_label_plus,
  std::string const& column_label_minus)
{
  Column v_plus(getColumn(column_label_plus));
  Column v_minus(getColumn(column_label_minus));
  af::shared<double> result((af::reserve(2*size())));
  for (int j=0; j<size(); j++) {
    if (!v_plus.ccp4_isnan(j)) result.push_back(v_plus.lookup(j));
    if (!v_minus.ccp4_isnan(j)) result.push_back(v_minus.lookup(j));
  }
  return result;
}

/* http://www.ccp4.ac.uk/dist/html/mtzMADmod.html
     F(+) = F + 0.5*D
     F(-) = F - 0.5*D
     SIGF(+) = sqrt( SIGF**2 + 0.25*SIGD**2 )
     SIGF(-) = SIGF(+)
 */
observation_arrays
Mtz::valid_delta_anomalous(
  std::string const& column_label_f_data,
  std::string const& column_label_f_sigmas,
  std::string const& column_label_d_data,
  std::string const& column_label_d_sigmas)
{
  Column h(getColumn("H"));
  Column k(getColumn("K"));
  Column l(getColumn("L"));
  Column v_fd(getColumn(column_label_f_data));
  Column v_fs(getColumn(column_label_f_sigmas));
  Column v_dd(getColumn(column_label_d_data));
  Column v_ds(getColumn(column_label_d_sigmas));
  observation_arrays result(2*size());
  for (int j=0; j<size(); j++) {
    if (   !v_fd.ccp4_isnan(j)
        && !v_fs.ccp4_isnan(j)
        && !v_dd.ccp4_isnan(j)
        && !v_ds.ccp4_isnan(j)) {
      result.indices.push_back(cctbx::miller::index<>(
        h.lookup(j), k.lookup(j), l.lookup(j)));
      result.indices.push_back(cctbx::miller::index<>(
        -h.lookup(j), -k.lookup(j), -l.lookup(j)));
      double fd = v_fd.lookup(j);
      double fs = v_fs.lookup(j);
      double ddh = v_dd.lookup(j) * .5;
      double dsh = v_ds.lookup(j) * .5;
      result.data.push_back(fd + ddh);
      result.data.push_back(fd - ddh);
      double s = std::sqrt(fs*fs + dsh*dsh);
      result.sigmas.push_back(std::sqrt(s));
      result.sigmas.push_back(std::sqrt(s));
    }
  }
  return result;
}

namespace {
  inline std::complex<double>
  polar_deg(double ampl, double phi)
  {
    return std::polar(ampl, phi*scitbx::constants::pi_180);
  }
}

af::shared<std::complex<double> >
Mtz::valid_complex(
  std::string const& column_label_ampl,
  std::string const& column_label_phi)
{
  Column v_ampl(getColumn(column_label_ampl));
  Column v_phi(getColumn(column_label_phi));
  af::shared<std::complex<double> > result((af::reserve(size())));
  for (int j=0; j<size(); j++) {
    if (v_ampl.ccp4_isnan(j) != v_phi.ccp4_isnan(j)) {
      throw Error("Unexpected NAN while extracting complex array.");
    }
    if (!v_ampl.ccp4_isnan(j)) {
      result.push_back(polar_deg(
        v_ampl.lookup(j), v_phi.lookup(j)));
    }
  }
  return result;
}

af::shared<std::complex<double> >
Mtz::valid_complex_anomalous(
  std::string const& column_label_ampl_plus,
  std::string const& column_label_phi_plus,
  std::string const& column_label_ampl_minus,
  std::string const& column_label_phi_minus)
{
  Column v_ampl_plus(getColumn(column_label_ampl_plus));
  Column v_phi_plus(getColumn(column_label_phi_plus));
  Column v_ampl_minus(getColumn(column_label_ampl_minus));
  Column v_phi_minus(getColumn(column_label_phi_minus));
  af::shared<std::complex<double> > result((af::reserve(2*size())));
  for (int j=0; j<size(); j++) {
    if (v_ampl_plus.ccp4_isnan(j) != v_phi_plus.ccp4_isnan(j)) {
      throw Error("Unexpected NAN while extracting complex array.");
    }
    if (!v_ampl_plus.ccp4_isnan(j)) {
      result.push_back(polar_deg(
        v_ampl_plus.lookup(j), v_phi_plus.lookup(j)));
    }
    if (v_ampl_minus.ccp4_isnan(j) != v_phi_minus.ccp4_isnan(j)) {
      throw Error("Unexpected NAN while extracting complex array.");
    }
    if (!v_ampl_minus.ccp4_isnan(j)) {
      result.push_back(polar_deg(
        v_ampl_minus.lookup(j), v_phi_minus.lookup(j)));
    }
  }
  return result;
}

af::shared<cctbx::hendrickson_lattman<> >
Mtz::valid_hl(
  std::string const& column_label_a,
  std::string const& column_label_b,
  std::string const& column_label_c,
  std::string const& column_label_d)
{
  Column v_a(getColumn(column_label_a));
  Column v_b(getColumn(column_label_b));
  Column v_c(getColumn(column_label_c));
  Column v_d(getColumn(column_label_d));
  af::shared<cctbx::hendrickson_lattman<> > result((af::reserve(size())));
  for (int j=0; j<size(); j++) {
    if (   v_a.ccp4_isnan(j) != v_b.ccp4_isnan(j)
        || v_a.ccp4_isnan(j) != v_c.ccp4_isnan(j)
        || v_a.ccp4_isnan(j) != v_d.ccp4_isnan(j)) {
      throw Error(
        "Unexpected NAN while extracting Hendrickson-Lattman array.");
    }
    if (!v_a.ccp4_isnan(j)) {
      result.push_back(cctbx::hendrickson_lattman<>(
        v_a.lookup(j), v_b.lookup(j), v_c.lookup(j), v_d.lookup(j)));
    }
  }
  return result;
}

af::shared<int>
Mtz::valid_integers(
  std::string const& column_label)
{
  Column v(getColumn(column_label));
  if (   v.type() != "H"
      && v.type() != "B"
      && v.type() != "Y"
      && v.type() != "I") {
    throw Error("Not an integer column.");
  }
  af::shared<int> result((af::reserve(size())));
  for (int j=0; j<size(); j++) {
    if (!v.ccp4_isnan(j)) result.push_back(static_cast<int>(v.lookup(j)));
  }
  return result;
}

std::string Crystal::crystal_name() {
  return std::string(p_xtal->xname);
}

std::string Crystal::project_name() {
  return std::string(p_xtal->pname);
}

cctbx::uctbx::unit_cell Crystal::UnitCell(){
  af::tiny<double,6> answer;
  for (int j=0; j<6; j++) {answer[j] = p_xtal->cell[j];}
  return cctbx::uctbx::unit_cell(answer);
}

int Crystal::ndatasets() {return p_xtal->nset;}

Dataset Crystal::getDataset(const int& setid){
  return Dataset(CMtz::MtzIsetInXtal(p_xtal, setid));
}

std::string Dataset::dataset_name() {
  return std::string(p_set->dname);
}

double Dataset::wavelength() {
  return p_set->wavelength;
}

int Dataset::ncolumns() {
  return p_set->ncol;
}

Column Dataset::getColumn(const int& colid){
  return Column(CMtz::MtzIcolInSet(p_set, colid));
}

std::string Column::label() {
  return std::string(p_col->label);
}

std::string Column::type() {
  return std::string(p_col->type);
}

}} // namespace iotbx::mtz
