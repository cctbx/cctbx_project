#include "iotbx/cppmtz.h"

namespace bpmtz = iotbx::mtz;

bpmtz::Mtz::Mtz(std::string s):mtz(CMtz::MtzGet(s.c_str(),1)){
    if (!mtz) {throw bpmtz::Error("Mtz read failed"); }
}
bpmtz::Mtz::~Mtz(){
  CMtz::MtzFree(mtz);}

std::string bpmtz::Mtz::title(){return std::string(mtz->title);}
std::string bpmtz::Mtz::SpaceGroup(){return mtz->mtzsymm.spcgrpname;}
int& bpmtz::Mtz::size(){return mtz->nref;}
int& bpmtz::Mtz::ncrystals() {return mtz->nxtal;}
int bpmtz::Mtz::ndatasets(const int& xtal) {return mtz->xtal[xtal]->nset;}
int& bpmtz::Mtz::ncolumns(const int& xtal, const int& set) {return mtz->xtal[xtal]->set[set]->ncol;}

af::shared<std::string> bpmtz::Mtz::history(){
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

af::shared<std::string> bpmtz::Mtz::columns(){
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

bpmtz::Column bpmtz::Mtz::getColumn(std::string s){
  return bpmtz::Column(CMtz::MtzColLookup(mtz,s.c_str()));
}

bpmtz::Crystal bpmtz::Mtz::getCrystal(const int& xtalid){
  return bpmtz::Crystal(CMtz::MtzIxtal(mtz, xtalid));
}

af::shared<double> bpmtz::Mtz::getShared(std::string s){
  bpmtz::Column col(this->getColumn(s));
  af::shared<double> shared_copy(this->size());
  for (int j=0; j<this->size(); j++) {
    shared_copy[j] = col.lookup(j);
  }
  return shared_copy;
}

af::shared< cctbx::miller::index<> > bpmtz::Mtz::MIx() {
  bpmtz::Column H(getColumn("H"));
  bpmtz::Column K(getColumn("K"));
  bpmtz::Column L(getColumn("L"));
  af::shared<cctbx::miller::index<> > miller(this->size());
  for (int j=0; j<this->size(); j++) {
    miller[j] = cctbx::miller::index<>(H.lookup(j),K.lookup(j),L.lookup(j));
  }
  return miller;
}

af::shared< cctbx::hendrickson_lattman<> > bpmtz::Mtz::HL(
  bpmtz::Column A,
  bpmtz::Column B,
  bpmtz::Column C,
  bpmtz::Column D)
{
  af::shared<cctbx::hendrickson_lattman<> > hl(this->size());
  for (int j=0; j<this->size(); j++) {
    hl[j] = cctbx::hendrickson_lattman<>(A.lookup(j),B.lookup(j),
                                         C.lookup(j),D.lookup(j));
  }
  return hl;
}

void bpmtz::Mtz::printHeader(int detail = 1){
  CMtz::ccp4_lhprt(mtz,detail);
}

void bpmtz::Mtz::printHeaderAdv(int detail = 1){
  CMtz::ccp4_lhprt_adv(mtz,detail);
}

af::shared<double> bpmtz::Mtz::UnitCell(const int& xtal){
  af::shared<double> answer;
  for (int j=0; j<6; j++) {answer.push_back(mtz->xtal[xtal]->cell[j]);}
  return answer;
}

std::string bpmtz::Crystal::crystal_name() {
  return std::string(p_xtal->xname);
}

std::string bpmtz::Crystal::project_name() {
  return std::string(p_xtal->pname);
}

af::shared<double> bpmtz::Crystal::UnitCell(){
  af::shared<double> answer;
  for (int j=0; j<6; j++) {answer.push_back(p_xtal->cell[j]);}
  return answer;
}

int bpmtz::Crystal::ndatasets() {return p_xtal->nset;}

bpmtz::Dataset bpmtz::Crystal::getDataset(const int& setid){
  return bpmtz::Dataset(CMtz::MtzIsetInXtal(p_xtal, setid));
}

std::string bpmtz::Dataset::dataset_name() {
  return std::string(p_set->dname);
}

double bpmtz::Dataset::wavelength() {
  return p_set->wavelength;
}

int bpmtz::Dataset::ncolumns() {
  return p_set->ncol;
}

bpmtz::Column bpmtz::Dataset::getColumn(const int& colid){
  return bpmtz::Column(CMtz::MtzIcolInSet(p_set, colid));
}

std::string bpmtz::Column::label() {
  return std::string(p_col->label);
}

std::string bpmtz::Column::type() {
  return std::string(p_col->type);
}

