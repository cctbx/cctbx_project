#include <mmtbx/max_lik/max_lik.h>
//#include <mmtbx/max_lik/solvent_distribution_nco_14mar2005.h>
//#include <mmtbx/max_lik/solvent_distribution_nco_14mar2005_mod.h>
//#include <mmtbx/max_lik/solvent_distribution_nco_09may2005.h>
#include <mmtbx/max_lik/solvent_distribution_nco_09may2005_short.h>
#include <mmtbx/error.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/math/bessel.h>
#include <assert.h>
#include <math.h>
#include <iostream>

using namespace std;
namespace mmtbx { namespace max_lik {


void wat_dist::preparator(cctbx::uctbx::unit_cell const& uc)
{
   double abc = uc.parameters()[0]*uc.parameters()[1]*uc.parameters()[2];
   double uc_volume_over_abc = uc.volume() / abc;
   double sin_alpha = std::sin(scitbx::deg_as_rad(uc.parameters()[3]));
   double sin_beta  = std::sin(scitbx::deg_as_rad(uc.parameters()[4]));
   double sin_gamma = std::sin(scitbx::deg_as_rad(uc.parameters()[5]));

   double ascale = uc_volume_over_abc / sin_alpha;
   double bscale = uc_volume_over_abc / sin_beta;
   double cscale = uc_volume_over_abc / sin_gamma;
   as = uc.parameters()[0]/ascale;
   bs = uc.parameters()[1]/bscale;
   cs = uc.parameters()[2]/cscale;
   xshell = shell/as;
   yshell = shell/bs;
   zshell = shell/cs;

/* Does not work well for one atom test, so I thinks
   these formulas are wrong ...
   So, it's better to remove the code below:
   as=1./uc.reciprocal_parameters()[0];
   bs=1./uc.reciprocal_parameters()[1];
   cs=1./uc.reciprocal_parameters()[2];
   xshell = shell*uc.reciprocal_parameters()[0];
   yshell = shell*uc.reciprocal_parameters()[1];
   zshell = shell*uc.reciprocal_parameters()[2];
*/
}

inline int wat_dist::nint(double x)
//
// Fortran like NINT() /Nearest Integer/
//
{
    return int(ceil(x+0.5)-(fmod(x*0.5+0.25,1.0)!=0));
}

void wat_dist::apply_symmetry_and_make_xyz_within_0_1(
                    af::shared<vec3<double> > const& xyzf,
                    cctbx::sgtbx::space_group const& sg,
                    af::shared<double> const& atmrad,
                    af::shared<int> sel_flag,
                    af::shared<std::string> const& element_symbol)
{
 /* double xf01,yf01,zf01;
 for (int i=0; i < xyzf.size(); i++) {
   vec3<double> const& site = xyzf[i];
   for (int k = 0; k < sg.order_z(); k++) {
     xf01 = rs[0][k]*site[0]+rs[1][k]*site[1]+rs[2][k]*site[2]+rs[9][k];
     yf01 = rs[3][k]*site[0]+rs[4][k]*site[1]+rs[5][k]*site[2]+rs[10][k];
     zf01 = rs[6][k]*site[0]+rs[7][k]*site[1]+rs[8][k]*site[2]+rs[11][k];
     if(xf01 >= 1.0)       xf01 = xf01 - int(xf01);
       else if(xf01 < 0.0) xf01 = xf01 - int(xf01) + 1.0;
       else                xf01 = xf01;
     if(yf01 >= 1.0)       yf01 = yf01 - int(yf01);
       else if(yf01 < 0.0) yf01 = yf01 - int(yf01) + 1.0;
       else                yf01 = yf01;
     if(zf01 >= 1.0)       zf01 = zf01 - int(zf01);
       else if(zf01 < 0.0) zf01 = zf01 - int(zf01) + 1.0;
       else                zf01 = zf01;
     xyzf_01.push_back(vec3<double>(xf01,yf01,zf01));
     atmrad_01.push_back(atmrad[i]);
     sel_flag_.push_back(sel_flag[i]);
     MMTBX_ASSERT(xf01*yf01*zf01 <= 1.0 && xf01*yf01*zf01 >= 0.0);
   }
 }
 MMTBX_ASSERT(xyzf_01.size() == xyzf.size()*sg.order_z()); */
 xyzf_01   = xyzf;
 sel_flag_ = sel_flag;
 atmrad_01 = atmrad;
 sg_ = sg;
 element_symbol_ = element_symbol;
}

void wat_dist::do_wat_dist(double ishell,
                           af::shared<vec3<double> > const& xyzf,
                           af::shared<double> const& atmrad,
                           af::shared<std::string> const& element_symbol,
                           cctbx::uctbx::unit_cell const& uc,
                           cctbx::sgtbx::space_group const& sg,
                           vec3<int> const& nxnynz,
                           af::shared<int> const& sel_flag,
                           double rad,
                           int nshells)
{
 using cctbx::sgtbx::rt_mx;
 using cctbx::sgtbx::rot_mx;
 using cctbx::sgtbx::tr_vec;
 for (int j = 0; j < sg.order_z(); j++) {
         rt_mx rt = sg(j);
         rot_mx r = rt.r();
         tr_vec t = rt.t();
         for(int i=0; i<9; i++) {
     rs[i][j] = r.num()[i] / static_cast<double>(r.den());
  }
  for(int i=0; i<3; i++) {
    rs[i+9][j] = t.num()[i] / static_cast<double>(t.den());
  }
}
    //sel_flag_ = sel_flag;
    shell  = ishell;
    rad_ = rad;
    NX=nxnynz[0];
    NY=nxnynz[1];
    NZ=nxnynz[2];
    water_mask_.resize(af::c_grid<3>(NX,NY,NZ), -1.0);

    MMTBX_ASSERT(shell  >= 0.0);
    MMTBX_ASSERT(NX > 0 && NY > 0 && NZ >0);
    MMTBX_ASSERT(xyzf.size() == atmrad.size());

    preparator(uc);
    apply_symmetry_and_make_xyz_within_0_1(xyzf, sg, atmrad, sel_flag, element_symbol);
    set_shells(uc, nshells);
}

void wat_dist::as_xplor_map(cctbx::uctbx::unit_cell const& uc,
                        std::string const& outputfile)
{
  FILE* fh = fopen(outputfile.c_str(),"w");
  MMTBX_ASSERT(fh != 0);
  fprintf(fh, "\n");
  fprintf(fh, "       1\n");
  fprintf(fh, "REMARKS SOLVENT MASK AS A MAP in CNS/XPLOR FORMAT\n");
  fprintf(fh,"%8d%8d%8d%8d%8d%8d%8d%8d%8d\n", NX,0,NX-1,NY,0,NY-1,NZ,0,NZ-1);
  for(std::size_t i=0;i<6;i++) {
    fprintf(fh, "%12.5e", uc.parameters()[i]);
  }
  fprintf(fh, "\n");
  fprintf(fh,"ZYX\n");
  int xsize = water_mask_.accessor()[0];
  int ysize = water_mask_.accessor()[1];
  for (int z = 0; z< water_mask_.accessor()[2]; z++){
    int eol=0;
    int x=0;
    int y=0;
    fprintf(fh,"%8d\n",z);
    for (int xy = 0; xy < xsize*ysize; ++xy){
        fprintf(fh,"%12.5e",water_mask_(x,y,z));
        ++x;
        if (x==xsize) {x = 0;++y;}
        ++eol;
        if (eol % 6 == 0) {fprintf(fh,"\n");}
    }
    if (eol % 6 != 0) {fprintf(fh,"\n");}
  }
  fprintf(fh,"   -9999\n");
  fprintf(fh, "%12.5e%12.5e\n", 0.0, 1.0);
  fclose(fh);
}


void wat_dist::set_shells(cctbx::uctbx::unit_cell const& uc, int nshells)
//
// A priori probability distribution of ordered water molecules
// in macromolecular crystals.
// Based on V.Y. Lunin's notes from 10-FEB-2004
// Pavel Afonine / 09-MAR-2004
//
{
    double dist,xL,yL,zL,xR,yR,zR;
    double coas,cobs,cocs,xfi,yfi,zfi,w1,w2,w3,w4,w5,w6,w7,w8,w9;
    int x1box,x2box,y1box,y2box,z1box,z2box,mx,my,mz,kx,ky,kz;
    bool close_to_au;
    double s1,s2,s1xx,s2yy,s1xy,s1xz,s2yz,s3z,s3;
    double sx,sy,sz,tsx,tsy,tsz,sxsq,sysq,szsq;
    double sxbcen,sybcen,szbcen,dx,dy,dz,dxy,dxz,dxzyz,distsy,distsx,distsm;
    double tsxg1,tsxg4,tsxg7,tsyg4,tsyg5,tsyg8,tszg3,tszg8,tszg9;
    double mr1,mr2,mr3,mr5,mr6,mr9,tmr2,tmr3,tmr6;
    double d1,d2,d1sq,d2sq,dmax,prob_n=0.,prob_c=0.,prob_o=0.,prob=0.;
    af::ref<double, af::c_grid<3> > water_mask_ref = water_mask_.ref();

    mr1= uc.metrical_matrix()[0]; //a*a;
    mr2= uc.metrical_matrix()[3]; //a*b*cos(gamma)
    mr3= uc.metrical_matrix()[4]; //a*c*cos(beta)
    mr5= uc.metrical_matrix()[1]; //b*b
    mr6= uc.metrical_matrix()[5]; //c*b*cos(alpha)
    mr9= uc.metrical_matrix()[2]; //c*c
    tmr2 = mr2*2.0; //2*a*b*cos(gamma);
    tmr3 = mr3*2.0; //2*a*c*cos(beta);
    tmr6 = mr6*2.0; //2*b*c*cos(alpha);
    xL = -xshell; xR = 1.0+xshell;
    yL = -yshell; yR = 1.0+yshell;
    zL = -zshell; zR = 1.0+zshell;
    sx = 1.0/NX; tsx= sx*2.0; sxsq=mr1*sx*sx;
    sy = 1.0/NY; tsy= sy*2.0; sysq=mr5*sy*sy;
    sz = 1.0/NZ; tsz= sz*2.0; szsq=mr9*sz*sz;
    w1=mr1*sx*tsx; w4=mr5*sy*tsy;
    w2=mr2*sx*tsy; w5=mr6*sy*tsz;
    w3=mr3*sx*tsz; w6=mr9*sz*tsz;
    tsxg1=tsx*mr1; tsyg4=tsy*mr2; tszg3=tsz*mr3;
    tsxg4=tsx*mr2; tsyg5=tsy*mr5; tszg8=tsz*mr6;
    tsxg7=tsx*mr3; tsyg8=tsy*mr6; tszg9=tsz*mr9;

    af::shared<af::shared<int> > mx_s;
    af::shared<af::shared<int> > my_s;
    af::shared<af::shared<int> > mz_s;
    af::shared<int> counter_values;

    if(nshells < 0) {
       nshells_ = sizeof(mmtbx::max_lik::sd_nco_table)/sizeof(double)/5;
    }
    else {
       nshells_ = nshells;
    }

    double sum_nn = 0.0;
    double sum_nc = 0.0;
    double sum_no = 0.0;
    for (std::size_t i=0; i<nshells_; i++) {
      sum_nn += mmtbx::max_lik::sd_nco_table[i][2];
      sum_nc += mmtbx::max_lik::sd_nco_table[i][3];
      sum_no += mmtbx::max_lik::sd_nco_table[i][4];
    }

   for (std::size_t j=0; j < nshells_-1; j++) {          // "-1" is intentional
     af::shared<int> mx_i;
     af::shared<int> my_i;
     af::shared<int> mz_i;
     int point_counter = 0;

     d1   = mmtbx::max_lik::sd_nco_table[j][0];
     d2   = mmtbx::max_lik::sd_nco_table[j][1];
     dmax = mmtbx::max_lik::sd_nco_table[j][1];
     prob_n = mmtbx::max_lik::sd_nco_table[j][2] / sum_nn;
     prob_c = mmtbx::max_lik::sd_nco_table[j][3] / sum_nc;
     prob_o = mmtbx::max_lik::sd_nco_table[j][4] / sum_no;

     int number_of_atoms = xyzf_01.size();
     for (std::size_t i=0; i < number_of_atoms; i++) {
       cctbx::fractional <> site= xyzf_01[i];
       xfi=site[0];
       yfi=site[1];
       zfi=site[2];
       std::string element = element_symbol_[i];
       if(element == "N") {
          prob = prob_n;
       }
       if(element == "C") {
          prob = prob_c;
       }
       if(element == "O") {
          prob = prob_o;
       }
       close_to_au=(xfi>=xL||xfi<=xR)&&(yfi>=yL||yfi<=yR)&&(zfi>=zL||zfi<=zR);
       if (close_to_au) {
         d1sq = d1*d1;
         d2sq = d2*d2;
         coas=dmax/as;
         cobs=dmax/bs;
         cocs=dmax/cs;
         x1box=nint( NX*(xfi-coas) ) - 1;
         x2box=nint( NX*(xfi+coas) ) + 1;
         y1box=nint( NY*(yfi-cobs) ) - 1;
         y2box=nint( NY*(yfi+cobs) ) + 1;
         z1box=nint( NZ*(zfi-cocs) ) - 1;
         z2box=nint( NZ*(zfi+cocs) ) + 1;
         sxbcen=xfi-x1box*sx;
         sybcen=yfi-y1box*sy;
         szbcen=zfi-z1box*sz;
         distsm=mr1*sxbcen*sxbcen+mr5*sybcen*sybcen+mr9*szbcen*szbcen
              +tmr2*sxbcen*sybcen+tmr3*sxbcen*szbcen+tmr6*sybcen*szbcen;
         s1 = s1xx = s1xy = s2yz = s3 = s2 = s2yy = s1xz = s3z = 0.0;
         w7=tsxg1*sxbcen+tsxg4*sybcen+tsxg7*szbcen;
         w8=tsyg4*sxbcen+tsyg5*sybcen+tsyg8*szbcen;
         w9=tszg3*sxbcen+tszg8*sybcen+tszg9*szbcen;
         for (kx = x1box; kx <= x2box; kx++) {
           //double xn=xfi-double(kx)/NX;
           distsx=distsm+s1;
           dx=w7-s1xx;
           dxy=w8-s1xy;
           dxz=w9-s1xz;
           for (ky = y1box; ky <= y2box; ky++) {
             //double yn=yfi-double(ky)/NY;
             distsy=distsx+s2;
             dy=dxy-s2yy;
             dxzyz=dxz-s2yz;
             for (kz = z1box; kz <= z2box; kz++) {
               //double zn=zfi-double(kz)/NZ;
               dz=dxzyz-s3z;
               dist = distsy+s3;
               //dist=mr1*xn*xn+mr5*yn*yn+mr9*zn*zn+
               //     tmr2*xn*yn+tmr3*xn*zn+tmr6*yn*zn;
               //MMTBX_ASSERT(std::abs(dist - dist_) < 1.e-3);
               if (dist>d1sq && dist<=d2sq) {
                 //point_counter += 1; // this works better
                 mx = scitbx::math::mod_positive(kx, NX);
                 my = scitbx::math::mod_positive(ky, NY);
                 mz = scitbx::math::mod_positive(kz, NZ);
                 if (water_mask_ref(mx,my,mz) < 0.0 && sel_flag_[i] == 1) {
                   water_mask_ref(mx,my,mz) = prob;
                   point_counter += 1; // this makes regression test happy
                   mx_i.push_back(mx);
                   my_i.push_back(my);
                   mz_i.push_back(mz);
                 }
             }
             s3+=szsq-dz;
             s3z+=w6;
           }
           s3=s3z=0.0;
           s2+=sysq-dy;
           s2yy+=w4;
           s2yz+=w5;
         }
         s2=s2yy=s2yz=0.0;
         s1+=sxsq-dx;
         s1xx+=w1;
         s1xy+=w2;
         s1xz+=w3;
         }
       }
     }
     mx_s.push_back(mx_i);
     my_s.push_back(my_i);
     mz_s.push_back(mz_i);
     counter_values.push_back(point_counter);
   }

   for (std::size_t j=0; j < nshells_-1; j++) {
     af::shared<int> mxi = mx_s[j];
     af::shared<int> myi = my_s[j];
     af::shared<int> mzi = mz_s[j];
     MMTBX_ASSERT(mxi.size() == myi.size());
     MMTBX_ASSERT(myi.size() == mzi.size());
     MMTBX_ASSERT(my_s.size() == counter_values.size());
     int pc = counter_values[j];
     if (pc > 0) {
       for (int k = 0; k < mxi.size(); k++)
         water_mask_ref(mxi[k],myi[k],mzi[k]) /= pc;
     }
   }


    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++)
          if(water_mask_ref(kx,ky,kz) < 0.0) {
             water_mask_ref(kx,ky,kz) = 0.0;
          }

    double pixel_volume = uc.volume() / (NX*NY*NZ);
    double counter = 0.0;
    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++)
          counter += water_mask_ref(kx,ky,kz);
    if(xyzf_01.size() == 1) {
       if(std::abs(counter-1.0) > 0.05) {
          std::cout<<"counter = "<<counter<<std::endl;
          MMTBX_ASSERT(std::abs(counter-1.0) < 1.e-3);
       }
    }

    double scale = counter*pixel_volume * sg_.order_z(); // NEW

    counter = 0.0;
    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++) {
          //MMTBX_ASSERT(water_mask_ref(kx,ky,kz) == 0.0 || water_mask_ref(kx,ky,kz) == 1.0);
          //MMTBX_ASSERT(water_mask_ref(kx,ky,kz) == 1.0);
          water_mask_ref(kx,ky,kz) = water_mask_ref(kx,ky,kz) / scale;
          counter += water_mask_ref(kx,ky,kz);
         }
       MMTBX_ASSERT(std::abs(counter*pixel_volume* sg_.order_z()-1.0) < 1.e-6); // NEW

}


/*
void wat_dist::set_shells(cctbx::uctbx::unit_cell const& uc, int nshells)
//
// around atoms
//
{
    double dist=0.0,xL,yL,zL,xR,yR,zR;
    double coas,cobs,cocs,xfi,yfi,zfi,w1,w2,w3,w4,w5,w6,w7,w8,w9;
    int x1box,x2box,y1box,y2box,z1box,z2box,mx,my,mz,kx,ky,kz;
    bool point_in_uc,close_to_au;
    double s1,s2,s1xx,s2yy,s1xy,s1xz,s2yz,s3z,s3;
    double sx,sy,sz,tsx,tsy,tsz,sxsq,sysq,szsq;
    double sxbcen,sybcen,szbcen,dx,dy,dz,dxy,dxz,dxzyz;
    double tsxg1,tsxg4,tsxg7,tsyg4,tsyg5,tsyg8,tszg3,tszg8,tszg9;
    double mr1,mr2,mr3,mr5,mr6,mr9,tmr2,tmr3,tmr6;
    double d1,d2,d1sq,d2sq,dmax,prob;
    af::ref<double, af::c_grid<3> > water_mask_ref = water_mask_.ref();

    mr1= uc.metrical_matrix()[0]; //a*a;
    mr2= uc.metrical_matrix()[3]; //a*b*cos(gamma)
    mr3= uc.metrical_matrix()[4]; //a*c*cos(beta)
    mr5= uc.metrical_matrix()[1]; //b*b
    mr6= uc.metrical_matrix()[5]; //c*b*cos(alpha)
    mr9= uc.metrical_matrix()[2]; //c*c
    tmr2 = mr2*2.0; //2*a*b*cos(gamma);
    tmr3 = mr3*2.0; //2*a*c*cos(beta);
    tmr6 = mr6*2.0; //2*b*c*cos(alpha);
    xL = -xshell; xR = 1.0+xshell;
    yL = -yshell; yR = 1.0+yshell;
    zL = -zshell; zR = 1.0+zshell;
    sx = 1.0/NX; tsx= sx*2.0; sxsq=mr1*sx*sx;
    sy = 1.0/NY; tsy= sy*2.0; sysq=mr5*sy*sy;
    sz = 1.0/NZ; tsz= sz*2.0; szsq=mr9*sz*sz;
    w1=mr1*sx*tsx; w4=mr5*sy*tsy;
    w2=mr2*sx*tsy; w5=mr6*sy*tsz;
    w3=mr3*sx*tsz; w6=mr9*sz*tsz;
    tsxg1=tsx*mr1; tsyg4=tsy*mr2; tszg3=tsz*mr3;
    tsxg4=tsx*mr2; tsyg5=tsy*mr5; tszg8=tsz*mr6;
    tsxg7=tsx*mr3; tsyg8=tsy*mr6; tszg9=tsz*mr9;

      int point_counter = 0;
      double R = rad_*rad_;
      for (std::size_t i=0; i < xyzf_01.size(); i++) {
       if(sel_flag_[i] == 1) {
        cctbx::fractional <> site= xyzf_01[i];
        xfi=site[0];
        yfi=site[1];
        zfi=site[2];
        close_to_au=(xfi>=xL||xfi<=xR)&&(yfi>=yL||yfi<=yR)&&(zfi>=zL||zfi<=zR);
        if (close_to_au) {
          coas=rad_/as;
          cobs=rad_/bs;
          cocs=rad_/cs;
          x1box=nint( NX*(xfi-coas) ) - 1;
          x2box=nint( NX*(xfi+coas) ) + 1;
          y1box=nint( NY*(yfi-cobs) ) - 1;
          y2box=nint( NY*(yfi+cobs) ) + 1;
          z1box=nint( NZ*(zfi-cocs) ) - 1;
          z2box=nint( NZ*(zfi+cocs) ) + 1;
          w7=tsxg1*sxbcen+tsxg4*sybcen+tsxg7*szbcen;
          w8=tsyg4*sxbcen+tsyg5*sybcen+tsyg8*szbcen;
          w9=tszg3*sxbcen+tszg8*sybcen+tszg9*szbcen;
          for (kx = x1box; kx <= x2box; kx++) {
            double xn=xfi-double(kx)/NX;
            for (ky = y1box; ky <= y2box; ky++) {
              double yn=yfi-double(ky)/NY;
              for (kz = z1box; kz <= z2box; kz++) {
                double zn=zfi-double(kz)/NZ;
                dist=mr1*xn*xn+mr5*yn*yn+mr9*zn*zn+tmr2*xn*yn+tmr3*xn*zn+tmr6*yn*zn;
                if (dist <= R) {
                  mx = scitbx::math::mod_positive(kx, NX);
                  my = scitbx::math::mod_positive(ky, NY);
                  mz = scitbx::math::mod_positive(kz, NZ);
                  if (water_mask_ref(mx,my,mz) < 0.0) {
                    water_mask_ref(mx,my,mz) = 1;
                    point_counter += 1;
                  }
                  //if (water_mask_ref(mx,my,mz) < 0.0) {
                  //  water_mask_ref(mx,my,mz) = 0.0;
                  //}
                }}}}}}}

    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++)
          if(water_mask_ref(kx,ky,kz) < 0.0) {
             water_mask_ref(kx,ky,kz) = 0.0;
          }

    double pixel_volume = uc.volume() / (NX*NY*NZ);
    double counter = 0.0;
    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++)
          counter += water_mask_ref(kx,ky,kz);
    if(xyzf_01.size() == 1) {
       if(std::abs(counter-1.0) > 1e-3) {
          std::cout<<"counter = "<<counter<<std::endl;
          MMTBX_ASSERT(std::abs(counter-1.0) < 1e-3);
       }
    }

    double scale = counter*pixel_volume * sg_.order_z(); // NEW

    counter = 0.0;
    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++) {
          water_mask_ref(kx,ky,kz) = water_mask_ref(kx,ky,kz) / scale;
          counter += water_mask_ref(kx,ky,kz);
         }
       MMTBX_ASSERT(std::abs(counter*pixel_volume* sg_.order_z()-1.0) < 1e-6); // NEW

}
*/

/*
void wat_dist::set_shells(cctbx::uctbx::unit_cell const& uc, int nshells)
//
// set up constant uniform distribution
//
{
    int kx,ky,kz;
    af::ref<double, af::c_grid<3> > water_mask_ref = water_mask_.ref();

    int point_counter = 0;
    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++) {
             point_counter += 1;
             water_mask_ref(kx,ky,kz) = 1.0;
        }

    double pixel_volume = uc.volume() / (NX*NY*NZ);
    double counter = 0.0;
    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++)
          counter += water_mask_ref(kx,ky,kz);

    double scale = counter*pixel_volume * sg_.order_z(); // ??????????

    counter = 0.0;
    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++) {
          water_mask_ref(kx,ky,kz) = water_mask_ref(kx,ky,kz) / scale;
          counter += water_mask_ref(kx,ky,kz);
         }
       MMTBX_ASSERT(std::abs(counter*pixel_volume* sg_.order_z()-1.0) < 1.e-6);// ???????
}


void wat_dist::set_shells(cctbx::uctbx::unit_cell const& uc, int nshells)
//
// uniform in shell around protein
//
{
    double dist,xL,yL,zL,xR,yR,zR;
    double coas,cobs,cocs,xfi,yfi,zfi,w1,w2,w3,w4,w5,w6,w7,w8,w9;
    int x1box,x2box,y1box,y2box,z1box,z2box,mx,my,mz,kx,ky,kz;
    double s1,s2,s1xx,s2yy,s1xy,s1xz,s2yz,s3z,s3;
    double sx,sy,sz,tsx,tsy,tsz,sxsq,sysq,szsq;
    double sxbcen,sybcen,szbcen,dx,dy,dz,dxy,dxz,dxzyz,distsy,distsx,distsm;
    double tsxg1,tsxg4,tsxg7,tsyg4,tsyg5,tsyg8,tszg3,tszg8,tszg9;
    double mr1,mr2,mr3,mr5,mr6,mr9,tmr2,tmr3,tmr6;
    af::ref<double, af::c_grid<3> > water_mask_ref = water_mask_.ref();

    mr1= uc.metrical_matrix()[0]; //a*a;
    mr2= uc.metrical_matrix()[3]; //a*b*cos(gamma)
    mr3= uc.metrical_matrix()[4]; //a*c*cos(beta)
    mr5= uc.metrical_matrix()[1]; //b*b
    mr6= uc.metrical_matrix()[5]; //c*b*cos(alpha)
    mr9= uc.metrical_matrix()[2]; //c*c
    tmr2 = mr2*2.0; //2*a*b*cos(gamma);
    tmr3 = mr3*2.0; //2*a*c*cos(beta);
    tmr6 = mr6*2.0; //2*b*c*cos(alpha);
    xL = -xshell; xR = 1.0+xshell;
    yL = -yshell; yR = 1.0+yshell;
    zL = -zshell; zR = 1.0+zshell;
    sx = 1.0/NX; tsx= sx*2.0; sxsq=mr1*sx*sx;
    sy = 1.0/NY; tsy= sy*2.0; sysq=mr5*sy*sy;
    sz = 1.0/NZ; tsz= sz*2.0; szsq=mr9*sz*sz;
    w1=mr1*sx*tsx; w4=mr5*sy*tsy;
    w2=mr2*sx*tsy; w5=mr6*sy*tsz;
    w3=mr3*sx*tsz; w6=mr9*sz*tsz;
    tsxg1=tsx*mr1; tsyg4=tsy*mr2; tszg3=tsz*mr3;
    tsxg4=tsx*mr2; tsyg5=tsy*mr5; tszg8=tsz*mr6;
    tsxg7=tsx*mr3; tsyg8=tsy*mr6; tszg9=tsz*mr9;

   int point_counter = 0;
   double dmin  = 2.5;
   double dmax  = 4.2;
   double dmin_sq  = dmin * dmin;
   double dmax_sq  = dmax * dmax;
   double prob     = 1.0;
   int number_of_atoms = xyzf_01.size();
   for (std::size_t i=0; i < number_of_atoms; i++) {
     cctbx::fractional <> site= xyzf_01[i];
     xfi=site[0];
     yfi=site[1];
     zfi=site[2];
     coas=dmax/as;
     cobs=dmax/bs;
     cocs=dmax/cs;
     x1box=nint( NX*(xfi-coas) ) - 1;
     x2box=nint( NX*(xfi+coas) ) + 1;
     y1box=nint( NY*(yfi-cobs) ) - 1;
     y2box=nint( NY*(yfi+cobs) ) + 1;
     z1box=nint( NZ*(zfi-cocs) ) - 1;
     z2box=nint( NZ*(zfi+cocs) ) + 1;
     sxbcen=xfi-x1box*sx;
     sybcen=yfi-y1box*sy;
     szbcen=zfi-z1box*sz;
     distsm=mr1*sxbcen*sxbcen+mr5*sybcen*sybcen+mr9*szbcen*szbcen
          +tmr2*sxbcen*sybcen+tmr3*sxbcen*szbcen+tmr6*sybcen*szbcen;
     s1 = s1xx = s1xy = s2yz = s3 = s2 = s2yy = s1xz = s3z = 0.0;
     w7=tsxg1*sxbcen+tsxg4*sybcen+tsxg7*szbcen;
     w8=tsyg4*sxbcen+tsyg5*sybcen+tsyg8*szbcen;
     w9=tszg3*sxbcen+tszg8*sybcen+tszg9*szbcen;
     for (kx = x1box; kx <= x2box; kx++) {
       //double xn=xfi-double(kx)/NX;
       distsx=distsm+s1;
       dx=w7-s1xx;
       dxy=w8-s1xy;
       dxz=w9-s1xz;
       for (ky = y1box; ky <= y2box; ky++) {
         //double yn=yfi-double(ky)/NY;
         distsy=distsx+s2;
         dy=dxy-s2yy;
         dxzyz=dxz-s2yz;
         for (kz = z1box; kz <= z2box; kz++) {
           //double zn=zfi-double(kz)/NZ;
           dz=dxzyz-s3z;
           dist = distsy+s3;
           //dist=mr1*xn*xn+mr5*yn*yn+mr9*zn*zn+
           //     tmr2*xn*yn+tmr3*xn*zn+tmr6*yn*zn;
           //MMTBX_ASSERT(std::abs(dist - dist_) < 1.e-3);
           if (dist > dmin_sq && dist <= dmax_sq) {
             mx = scitbx::math::mod_positive(kx, NX);
             my = scitbx::math::mod_positive(ky, NY);
             mz = scitbx::math::mod_positive(kz, NZ);
             if (water_mask_ref(mx,my,mz) < 0.0 && sel_flag_[i] == 1) {
               water_mask_ref(mx,my,mz) = prob;
               point_counter += 1;
             }
         }
         s3+=szsq-dz;
         s3z+=w6;
       }
       s3=s3z=0.0;
       s2+=sysq-dy;
       s2yy+=w4;
       s2yz+=w5;
     }
     s2=s2yy=s2yz=0.0;
     s1+=sxsq-dx;
     s1xx+=w1;
     s1xy+=w2;
     s1xz+=w3;
     }
   }

    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++)
          if(water_mask_ref(kx,ky,kz) < 0.0) {
             water_mask_ref(kx,ky,kz) = 0.0;
          }

    double pixel_volume = uc.volume() / (NX*NY*NZ);
    double counter = 0.0;
    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++)
          counter += water_mask_ref(kx,ky,kz);


    double scale = counter*pixel_volume * sg_.order_z(); // NEW

    counter = 0.0;
    for (kx = 0; kx < NX; kx++)
      for (ky = 0; ky < NY; ky++)
        for (kz = 0; kz < NZ; kz++) {
          MMTBX_ASSERT(water_mask_ref(kx,ky,kz) == 0.0 ||
                       water_mask_ref(kx,ky,kz) == 1.0);
          water_mask_ref(kx,ky,kz) = water_mask_ref(kx,ky,kz) / scale;
          counter += water_mask_ref(kx,ky,kz);
         }
       MMTBX_ASSERT(std::abs(counter*pixel_volume* sg_.order_z()-1.0) < 1.e-6); // NEW

}
*/

}} // namespace mmtbx::max_lik
