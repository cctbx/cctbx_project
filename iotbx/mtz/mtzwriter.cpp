#include <iotbx/mtzwriter.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <scitbx/mat3.h>

using namespace CMtz;

iotbx::mtz::MtzWriter::MtzWriter():mtz(MtzMalloc(0,0)) { mtz->nref = 0;}

iotbx::mtz::MtzWriter::~MtzWriter() { MtzFree(mtz);}

void iotbx::mtz::MtzWriter::setTitle(const std::string& title){
  sprintf(mtz->title,"%s",title.c_str());
}

void iotbx::mtz::MtzWriter::setSpaceGroup(const cctbx::sgtbx::space_group& sg){
  cctbx::sgtbx::space_group_type sgtype(sg.type());
  mtz->mtzsymm.spcgrp = sgtype.number();
  std::string whole_symbol = sgtype.lookup_symbol();
  std::size_t ind;
  while ((ind=whole_symbol.find(' '))!=std::string::npos) {
    whole_symbol = whole_symbol.erase(ind,1);
  }
  sprintf(mtz->mtzsymm.spcgrpname,"%s",whole_symbol.c_str());
  mtz->mtzsymm.nsym = sg.order_z();         
  mtz->mtzsymm.nsymp = sg.order_p();            
  mtz->mtzsymm.symtyp=sg.conventional_centring_type_symbol();
  cctbx::sgtbx::matrix_group::code mgcode = sg.point_group_type();
  sprintf(mtz->mtzsymm.pgname,"%c",mgcode.label());
  for(std::size_t i=0;i<sg.order_z();i++) {
    cctbx::sgtbx::rt_mx s = sg(i);
    cctbx::sgtbx::rot_mx r = s.r();
    cctbx::sgtbx::tr_vec t = s.t();
    scitbx::mat3<int> r_num = r.num();
    int r_den = r.den();
    scitbx::vec3<int> t_num = t.num();
    int t_den = t.den();
    for (std::size_t p=0;p<3;p++) {
      for (std::size_t q=0;q<3;q++) {
        mtz->mtzsymm.sym[i][p][q] = r_num(p,q)/r_den; }
      mtz->mtzsymm.sym[i][p][3] = t_num[p];
      mtz->mtzsymm.sym[i][3][p] = 0.0;
    }
    mtz->mtzsymm.sym[i][3][3] = 1.0;
  }
}

void iotbx::mtz::MtzWriter::oneCrystal(const std::string& crystal,
  const std::string& project,const cctbx::uctbx::unit_cell& uc){

  std::size_t ind;
  std::string strip_crystal = crystal;
  std::string strip_project = project;
  while ((ind=strip_crystal.find(' '))!=std::string::npos) {
    strip_crystal.erase(ind,1);
  }
  while ((ind=strip_project.find(' '))!=std::string::npos) {
    strip_project.erase(ind,1);
  }
  float cell[6];
  for (ind=0;ind<6;++ind) {cell[ind]=uc.parameters()[ind];}
  onextal = MtzAddXtal(mtz, strip_crystal.c_str(), strip_project.c_str(),
            cell);
}

void iotbx::mtz::MtzWriter::oneDataset(const std::string& dataset,
                                       const double& w){
  float wavelength=w;
  std::size_t ind;
  std::string strip_dataset = dataset;
  while ((ind=strip_dataset.find(' '))!=std::string::npos) {
    strip_dataset.erase(ind,1);
  }
  oneset = MtzAddDataset(mtz, onextal, strip_dataset.c_str(), wavelength);
}

void iotbx::mtz::MtzWriter::write(const std::string& filename){
  MtzPut(mtz, const_cast<char*>(filename.c_str()));
}
/*

  int hmin=-1;
  int hmax=1;
  int kmin=-1;
  int kmax=1;
  int lmin=-1;
  int lmax=1;
  
  MTZCOL *colout[4];
  
  colout[0]=MtzAddColumn (P1mtz,set,"H","H");
  colout[1]=MtzAddColumn (P1mtz,set,"K","H");
  colout[2]=MtzAddColumn (P1mtz,set,"L","H");
  colout[3]=MtzAddColumn (P1mtz,set,"Fobs","F");
  
  int ncol=1;
  float adata[4];
  
  int iref=0;
  for (int h = hmin; h <= hmax; h++)
    for (int k = kmin; k <= kmax; k++)
      for (int l = lmin; l <= lmax; l++)
        if (!(l == 0 && h < 0) && 
            !(h == 0 && l == 0 && k < 0)&&
            !(h==0 && k==0 && l==0)) {
            iref+=1;
            adata[0]=(float)h;
            adata[1]=(float)k;
            adata[2]=(float)l;
            adata[3]=5728.02;
            ccp4_lwrefl(P1mtz, adata, colout, 4, iref);
        }
*/
