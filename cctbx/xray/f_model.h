#ifndef CCTBX_XRAY_F_MODEL_H
#define CCTBX_XRAY_F_MODEL_H

#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>
#include <complex>
#include <cmath>
#include <scitbx/math/bessel.h>
#include <cctbx/hendrickson_lattman.h>
#include <cctbx/miller/lookup_utils.h>
#include <cctbx/uctbx.h>

namespace cctbx { namespace xray { namespace f_model_core_data {

 template<typename FloatType>
 class f_model_core_data_derivative_holder{
 protected:
   FloatType koverall_;
   FloatType ksol_;
   FloatType usol_;
   FloatType kpart_;
   FloatType upart_;
   scitbx::sym_mat3<FloatType> ustar_;
 public:
   f_model_core_data_derivative_holder():
   koverall_(0),
   ksol_(0),
   usol_(0),
   kpart_(0),
   upart_(0),
   ustar_( 0,0,0,0,0,0 )
   {}

   FloatType ksol(){ return (ksol_); }
   FloatType usol(){ return (usol_); }
   FloatType kpart(){ return (kpart_); }
   FloatType upart(){ return (upart_); }
   FloatType koverall(){ return (koverall_); }
   scitbx::sym_mat3<FloatType> ustar(){ return (ustar_); }

   void ksol(FloatType tmp ){ ksol_=tmp;}
   void usol(FloatType tmp ){ usol_=tmp; }
   void kpart(FloatType tmp ){ kpart_=tmp; }
   void upart(FloatType tmp ){ upart_=tmp; }
   void koverall(FloatType tmp ){ koverall_=tmp; }
   void ustar( scitbx::sym_mat3<FloatType> tmp ){ ustar_=tmp; }
   void accumulate( f_model_core_data_derivative_holder<FloatType> other )
   {
     koverall_ += other.koverall();
     ksol_  += other.ksol();
     usol_  += other.usol();
     kpart_ += other.kpart();
     upart_ += other.upart();
     ustar_ += other.ustar();
   }
 };


 template<typename FloatType>
 class f_model_core_data{
 public:
   // You want to use this constructor
   f_model_core_data(
      scitbx::af::const_ref<cctbx::miller::index<> >  const& hkl,           //1 indices for calculated data
      scitbx::af::const_ref<std::complex<FloatType> > const& f_atoms,       //2 f calc
      scitbx::af::const_ref<std::complex<FloatType> > const& f_mask,        //3 f bulk solvent
      cctbx::uctbx::unit_cell                         const& unit_cell,     //4 unit cell
      FloatType                                       const& k_overall,     //7 overall scale factor
      scitbx::sym_mat3<FloatType>                     const& u_star,        //8 overall aniso scale
      FloatType                                       const& k_sol,         //9 bulk solvent scale
      FloatType                                       const& u_sol,         //0 bulk solvent B
      scitbx::af::const_ref<std::complex<FloatType> > const& f_part,        //1 f partial structure
      FloatType                                       const& k_part,        //2 partial structure scale
      FloatType                                       const& u_part         //3 partial structure B
      ):
     k_overall_(k_overall),
     u_star_(u_star),
     k_sol_(k_sol),
     u_sol_(u_sol),
     k_part_(k_part),
     u_part_(u_part),
     recompute_aniso_(true),
     recompute_bulk_(true),
     recompute_part_(true),
     recompute_total_(true)
     {
       CCTBX_ASSERT( hkl.size() > 0 );
       CCTBX_ASSERT( hkl.size() == f_atoms.size() );
       CCTBX_ASSERT( (f_mask.size() == 0) || (hkl.size() == f_mask.size()) );
       CCTBX_ASSERT( (f_part.size() == 0) || (hkl.size() == f_part.size()) );
       CCTBX_ASSERT( k_overall > 0 );
       CCTBX_ASSERT( k_sol  >=0 );
       CCTBX_ASSERT( k_part >=0 );
       CCTBX_ASSERT( u_sol >= 0);
       CCTBX_ASSERT( u_part >= 0 );

       if (k_overall_==0){
         recompute_aniso_=false;
       }
       if (k_part_==0){
         recompute_part_=false;
       }
       if (k_sol_==0){
         recompute_bulk_=false;
       }


       for (std::size_t ii=0; ii< hkl.size(); ii++){
         hkl_.push_back( hkl[ii] );
         f_atoms_.push_back( f_atoms[ii] );

         if (f_mask.size()>0){
           f_mask_.push_back( f_mask[ii] );
         } else {
           f_mask_.push_back( std::complex<FloatType>(0,0) );
         }

         if (f_part.size()>0){
           f_part_.push_back( f_part[ii] );
         } else {
           f_part_.push_back( std::complex<FloatType>(0,0) );
         }
         // we want to compute d_star_sq for each hkl
         d_star_sq_.push_back( unit_cell.d_star_sq(hkl[ii]) ) ;
         // init the scale and other arrays
         aniso_scale_.push_back( 1.0 );
         //log_part_aniso_scale_.push_back( 1.0 );
         bulk_scale_.push_back( 0.0 );
         part_scale_.push_back( 0.0 );
         f_model_core_data_.push_back( std::complex<FloatType>(0,0) );
         // update all
         update_all( ii );
       }
     }
     /* or this one of course */
     f_model_core_data(
      scitbx::af::const_ref<cctbx::miller::index<> >  const& hkl,           //1 indices for calculated data
      scitbx::af::const_ref<std::complex<FloatType> > const& f_atoms,       //2 f calc
      scitbx::af::const_ref<std::complex<FloatType> > const& f_mask,        //3 f bulk solvent
      scitbx::af::const_ref<FloatType>                const& d_star_sq,     //4 unit cell
      FloatType                                       const& k_overall,     //7 overall scale factor
      scitbx::sym_mat3<FloatType>                     const& u_star,        //8 overall aniso scale
      FloatType                                       const& k_sol,         //9 bulk solvent scale
      FloatType                                       const& u_sol,         //0 bulk solvent B
      scitbx::af::const_ref<std::complex<FloatType> > const& f_part,        //1 f partial structure
      FloatType                                       const& k_part,        //2 partial structure scale
      FloatType                                       const& u_part         //3 partial structure B
      ):
     k_overall_(k_overall),
     u_star_(u_star),
     k_sol_(k_sol),
     u_sol_(u_sol),
     k_part_(k_part),
     u_part_(u_part),
     recompute_aniso_(true),
     recompute_bulk_(true),
     recompute_part_(true),
     recompute_total_(true)
     {
       CCTBX_ASSERT( hkl.size() > 0 );
       CCTBX_ASSERT( hkl.size() == f_atoms.size() );
       CCTBX_ASSERT( hkl.size() == d_star_sq.size() );
       CCTBX_ASSERT( (f_mask.size() == 0) || (hkl.size() == f_mask.size()) );
       CCTBX_ASSERT( (f_part.size() == 0) || (hkl.size() == f_part.size()) );
       CCTBX_ASSERT( k_overall > 0 );
       CCTBX_ASSERT( k_sol  >=0 );
       CCTBX_ASSERT( k_part >=0 );
       CCTBX_ASSERT( u_sol >= 0);
       CCTBX_ASSERT( u_part >= 0 );

       if (k_overall_==0){
         recompute_aniso_=false;
       }
       if (k_part_==0){
         recompute_part_=false;
       }
       if (k_sol_==0){
         recompute_bulk_=false;
       }



       for (std::size_t ii=0; ii< hkl.size(); ii++){
         hkl_.push_back( hkl[ii] );
         f_atoms_.push_back( f_atoms[ii] );
         if (f_mask.size()>0){
           f_mask_.push_back( f_mask[ii] );
         } else {
           f_mask_.push_back( std::complex<FloatType>(0,0) );
         }
         if (f_part.size()>0){
           f_part_.push_back( f_part[ii] );
         } else {
           f_part_.push_back( std::complex<FloatType>(0,0) );
         }
         // we want to compute d_star_sq for each hkl
         d_star_sq_.push_back( d_star_sq[ii] ) ;
         // init the scale and other arrays
         aniso_scale_.push_back( 1.0 );
         //log_part_aniso_scale_.push_back( 1.0 );
         bulk_scale_.push_back( 0.0 );
         part_scale_.push_back( 0.0 );
         f_model_core_data_.push_back( std::complex<FloatType>(0,0) );
         // update all
         update_all( ii );
       }
     }




       //-----------------------------


     FloatType d_target_d_koverall(FloatType const& d_target_d_a,
                                   FloatType const& d_target_d_b,
                                   std::size_t const& ii )
     {
       if (recompute_aniso_){
         compute_aniso_scale( ii );
       }
       if (recompute_bulk_){
         compute_bulk_scale( ii );
       }
       if (recompute_part_){
         compute_part_scale( ii );
       }

       FloatType va,vb,ua,ub,wa,wb,g1,g2,fa,fb;
       ua=f_atoms_[ii].real();
       ub=f_atoms_[ii].imag();
       va=f_mask_[ii].real();
       vb=f_mask_[ii].imag();
       wa=f_part_[ii].real();
       wb=f_part_[ii].imag();

       g1=k_sol_*bulk_scale_[ii];
       g2=k_part_*part_scale_[ii];

       fa=(ua+g1*va+g2*wa);

       fb=(ub+g1*vb+g2*wb);

       FloatType tmp_a, tmp_b, result;

       tmp_a = fa*aniso_scale_[ii];
       tmp_b = fb*aniso_scale_[ii];
       result = tmp_a*d_target_d_a + tmp_b*d_target_d_b;
       return (result);
     }


     FloatType d_target_d_ksol(FloatType const& d_target_d_a,
                               FloatType const& d_target_d_b,
                               std::size_t const& ii )
     {
       if (recompute_aniso_){
         compute_aniso_scale( ii );
       }
       if (recompute_bulk_){
         compute_bulk_scale( ii );
       }
       if (recompute_part_){
         compute_part_scale( ii );
       }

       FloatType va,vb,ua,ub,wa,wb,g1,g2,fa,fb;
       ua=f_atoms_[ii].real();
       ub=f_atoms_[ii].imag();
       va=f_mask_[ii].real();
       vb=f_mask_[ii].imag();
       wa=f_part_[ii].real();
       wb=f_part_[ii].imag();
       g1=k_sol_*bulk_scale_[ii];
       g2=k_part_*part_scale_[ii];
       fa=(ua+g1*va+g2*wa);
       fb=(ub+g1*vb+g2*wb);

       FloatType tmp_a, tmp_b, result;
       tmp_a = k_overall_*aniso_scale_[ii]*bulk_scale_[ii]*va;
       tmp_b = k_overall_*aniso_scale_[ii]*bulk_scale_[ii]*vb;
       result = tmp_a*d_target_d_a + tmp_b*d_target_d_b ;
       return (result);
     }

     FloatType d_target_d_kpart(FloatType const& d_target_d_a,
                                FloatType const& d_target_d_b,
                                std::size_t const& ii )
     {
       if (recompute_aniso_){
         compute_aniso_scale( ii );
       }
       if (recompute_bulk_){
         compute_bulk_scale( ii );
       }
       if (recompute_part_){
         compute_part_scale( ii );
       }

       FloatType va,vb,ua,ub,wa,wb,g1,g2,fa,fb;
       ua=f_atoms_[ii].real();
       ub=f_atoms_[ii].imag();
       va=f_mask_[ii].real();
       vb=f_mask_[ii].imag();
       wa=f_part_[ii].real();
       wb=f_part_[ii].imag();
       g1=k_sol_*bulk_scale_[ii];
       g2=k_part_*part_scale_[ii];
       fa=(ua+g1*va+g2*wa);
       fb=(ub+g1*vb+g2*wb);

       FloatType tmp_a, tmp_b, result;
       tmp_a = k_overall_*aniso_scale_[ii]*part_scale_[ii]*wa;
       tmp_b = k_overall_*aniso_scale_[ii]*part_scale_[ii]*wb;
       result = tmp_a*d_target_d_a + tmp_b*d_target_d_b ;
       return (result);
     }





     // This function compute the gradients for a single HKL
     f_model_core_data_derivative_holder<FloatType>
     d_target_d_all( FloatType const& d_target_d_a,
                     FloatType const& d_target_d_b,
                     std::size_t const& ii,
                     scitbx::af::const_ref<bool> const& gradient_flags)
     {

       CCTBX_ASSERT( gradient_flags.size() == 6 );
       if (recompute_aniso_){
         compute_aniso_scale( ii );
       }
       if (recompute_bulk_){
         compute_bulk_scale( ii );
       }
       if (recompute_part_){
         compute_part_scale( ii );
       }

       f_model_core_data_derivative_holder<FloatType> tmp_derivs;
       FloatType va,vb,ua,ub,wa,wb,g1,g2,fa,fb;
       ua=f_atoms_[ii].real();
       ub=f_atoms_[ii].imag();
       va=f_mask_[ii].real();
       vb=f_mask_[ii].imag();
       wa=f_part_[ii].real();
       wb=f_part_[ii].imag();
       g1=k_sol_*bulk_scale_[ii];
       g2=k_part_*part_scale_[ii];
       fa=(ua+g1*va+g2*wa);
       fb=(ub+g1*vb+g2*wb);

       FloatType tmp_a, tmp_b;

       // k_overall
       if (gradient_flags[0]){
         tmp_a = fa*aniso_scale_[ii];
         tmp_b = fb*aniso_scale_[ii];
         tmp_derivs.koverall( tmp_a*d_target_d_a + tmp_b*d_target_d_b );
       }
       // u_star
       FloatType tps=scitbx::constants::pi;
       tps=tps*tps*2.0;
       if (gradient_flags[1]){
          tmp_a = fa*aniso_scale_[ii]*k_overall_;
          tmp_b = fb*aniso_scale_[ii]*k_overall_;
          scitbx::sym_mat3<FloatType> tmp_ustar(
              -2.0*tps*hkl_[ii][0]*hkl_[ii][0],
              -2.0*tps*hkl_[ii][1]*hkl_[ii][1],
              -2.0*tps*hkl_[ii][2]*hkl_[ii][2],
              -4.0*tps*hkl_[ii][0]*hkl_[ii][1],
              -4.0*tps*hkl_[ii][0]*hkl_[ii][2],
              -4.0*tps*hkl_[ii][1]*hkl_[ii][2]);
          tmp_derivs.ustar( tmp_ustar*(tmp_a*d_target_d_a + tmp_b*d_target_d_b) );
       }

       // k_sol
       if (gradient_flags[2]){
         tmp_a = k_overall_*aniso_scale_[ii]*bulk_scale_[ii]*va;
         tmp_b = k_overall_*aniso_scale_[ii]*bulk_scale_[ii]*vb;
         tmp_derivs.ksol( tmp_a*d_target_d_a + tmp_b*d_target_d_b );
       }
       // k_part
       if (gradient_flags[3]){
         tmp_a = k_overall_*aniso_scale_[ii]*part_scale_[ii]*wa;
         tmp_b = k_overall_*aniso_scale_[ii]*part_scale_[ii]*wb;
         tmp_derivs.kpart( tmp_a*d_target_d_a + tmp_b*d_target_d_b );
       }
       // b_sol
       if (gradient_flags[4]){
         tmp_a =  -k_overall_*aniso_scale_[ii]*bulk_scale_[ii]*k_sol_*va*d_star_sq_[ii]*tps;
         tmp_b =  -k_overall_*aniso_scale_[ii]*bulk_scale_[ii]*k_sol_*vb*d_star_sq_[ii]*tps;
         tmp_derivs.usol( tmp_a*d_target_d_a + tmp_b*d_target_d_b );
       }
       // b_part
       if (gradient_flags[5]){
         tmp_a =  -k_overall_*aniso_scale_[ii]*part_scale_[ii]*k_part_*wa*d_star_sq_[ii]*tps;
         tmp_b =  -k_overall_*aniso_scale_[ii]*part_scale_[ii]*k_part_*wb*d_star_sq_[ii]*tps;
         tmp_derivs.upart( tmp_a*d_target_d_a + tmp_b*d_target_d_b );
       }
       return ( tmp_derivs );
     }

     f_model_core_data_derivative_holder<FloatType>
     d_target_d_all( scitbx::af::const_ref<FloatType> const& d_target_d_a,
                     scitbx::af::const_ref<FloatType> const& d_target_d_b,
                     scitbx::af::const_ref<bool> const& gradient_flags )
     {
       f_model_core_data_derivative_holder<FloatType> result;
       CCTBX_ASSERT( d_target_d_b.size()==hkl_.size() );
       CCTBX_ASSERT( d_target_d_a.size()==hkl_.size() );
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         result.accumulate( d_target_d_all(d_target_d_a[ii],
                                           d_target_d_b[ii],
                                           ii,
                                           gradient_flags) );
       }
       return (result);
     }


     f_model_core_data_derivative_holder<FloatType>
       d_target_d_all( scitbx::af::const_ref<std::complex<FloatType> > const& d_target_d_ab,
                       scitbx::af::const_ref<bool> const& gradient_flags )
     {
       f_model_core_data_derivative_holder<FloatType> result;
       CCTBX_ASSERT( d_target_d_ab.size()==hkl_.size() );
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         result.accumulate( d_target_d_all(d_target_d_ab[ii].real(),
                                           d_target_d_ab[ii].imag(),
                                           ii,
                                           gradient_flags) );
       }
       return (result);
     }


     f_model_core_data_derivative_holder<FloatType>
       d_target_d_all( scitbx::af::const_ref<FloatType> const& d_target_d_fmodel,
                       scitbx::af::const_ref<bool> const& gradient_flags )
     {
       f_model_core_data_derivative_holder<FloatType> result;
       FloatType dfda, dfdb, fm, a, b;
       CCTBX_ASSERT( d_target_d_fmodel.size()==hkl_.size() );
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         a=f_model_core_data_[ii].real();
         b=f_model_core_data_[ii].imag();
         fm=std::sqrt(a*a+b*b);
         if ( fm>0 ){
           dfda=d_target_d_fmodel[ii]*a/fm;
           dfdb=d_target_d_fmodel[ii]*b/fm;
           result.accumulate( d_target_d_all(dfda,
                                             dfdb,
                                             ii,
                                             gradient_flags) );
         }
       }
       return (result);
     }




     //! update entries please.
     void refresh()
     {
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         if (recompute_aniso_){
           compute_aniso_scale( ii );
         }
         if (recompute_bulk_){
           compute_bulk_scale( ii );
         }
         if (recompute_part_){
           compute_part_scale( ii );
         }
         if (recompute_total_){
           compute_f_model_core_data( ii );
         }
       }
       recompute_aniso_=false;
       recompute_part_=false;
       recompute_bulk_=false;
       recompute_total_=false;
     }

     // reload sf's
     void renew_fatoms( scitbx::af::const_ref< std::complex< FloatType> > const& new_f_atoms )
     {
       CCTBX_ASSERT( new_f_atoms.size() == hkl_.size() );
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         f_atoms_[ii]=new_f_atoms[ii];
       }
       recompute_total_=true;

     }
     void renew_fmask( scitbx::af::const_ref< std::complex< FloatType> > const& new_f_mask)
     {
       CCTBX_ASSERT( new_f_mask.size() == hkl_.size() );
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         f_mask_[ii]=new_f_mask[ii];
       }
       recompute_total_=true;
     }
     void renew_fpart( scitbx::af::const_ref< std::complex< FloatType> > const& new_f_part)
     {
       CCTBX_ASSERT( new_f_part.size() == hkl_.size() );
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         f_part_[ii]=new_f_part[ii];
       }
       recompute_total_=true;
     }

     // reload bulk solvent params
     void renew_overall_scale_parameters( FloatType const& new_k_overall,
                                          scitbx::sym_mat3<FloatType> const& new_u_star)
     {
       CCTBX_ASSERT( new_k_overall > 0 );
       k_overall_= new_k_overall;
       u_star_ = new_u_star;
       recompute_total_=true;
       recompute_aniso_=true;
       refresh();
     }

     void renew_bulk_solvent_scale_parameters(  FloatType const& new_k_sol,
                                                FloatType const& new_u_sol )
     {
       CCTBX_ASSERT( new_u_sol>=0 );
       CCTBX_ASSERT( new_k_sol>=0 );
       k_sol_=new_k_sol;
       u_sol_=new_u_sol;
       recompute_total_=true;
       recompute_bulk_=true;
       refresh();

     }

     void renew_partial_structure_scale_parameters(  FloatType const& new_k_part,
                                                     FloatType const& new_u_part )
     {
       CCTBX_ASSERT( new_u_part>=0 );
       CCTBX_ASSERT( new_k_part>=0 );
       k_part_=new_k_part;
       u_part_=new_u_part;
       recompute_total_=true;
       recompute_part_=true;
       refresh();
     }

     // get some scale and sf's
     scitbx::af::shared< std::complex<FloatType> > f_bulk()
     {
       scitbx::af::shared< std::complex<FloatType> > result;
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         if (recompute_bulk_){
           compute_bulk_scale(ii);
         }
         result.push_back( f_mask_[ii]*k_sol_*bulk_scale_[ii] );
       }
       return ( result );
     }

     scitbx::af::shared< FloatType > bulk_scale()
     {
       scitbx::af::shared< FloatType > result;
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         if (recompute_bulk_){
           compute_bulk_scale(ii);
         }
         result.push_back( k_sol_*bulk_scale_[ii] );
       }
       return ( result );
     }

     scitbx::af::shared< std::complex<FloatType> > f_part()
     {
       scitbx::af::shared< std::complex<FloatType> > result;
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         if (recompute_part_){
           compute_part_scale(ii);
         }
         result.push_back( f_part_[ii]*k_part_*part_scale_[ii] );
       }
       return ( result );
     }

     scitbx::af::shared< FloatType > part_scale()
     {
       scitbx::af::shared< FloatType > result;
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         if (recompute_part_){
           compute_part_scale(ii);
         }
         result.push_back( k_part_*part_scale_[ii] );
       }
       return ( result );
     }

     scitbx::af::shared<std::complex<FloatType> > get_f_model_core_data()
     {
       if (recompute_total_){
         refresh();
       }
       return ( f_model_core_data_ ) ;
     }


     scitbx::af::shared< FloatType > overall_scale()
     {
       scitbx::af::shared< FloatType > result;
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         if (recompute_aniso_){
           compute_aniso_scale(ii);
         }
         result.push_back( k_overall_*aniso_scale_[ii] );
       }
       return ( result );
     }

     scitbx::af::shared< FloatType > fu_star()
     {
       scitbx::af::shared< FloatType > result;
       for (std::size_t ii=0;ii<hkl_.size();ii++){
         if (recompute_aniso_){
           compute_aniso_scale(ii);
         }
         result.push_back( aniso_scale_[ii] );
       }
       return ( result );
     }


     scitbx::af::shared<std::complex<FloatType> > f_atoms()
     {
       return ( f_atoms_ ) ;
     }



     // model parameters
     FloatType koverall(){
       return k_overall_;
     }
     FloatType ksol()
     {
       return( k_sol_ );
     }
     FloatType kpart()
     {
       return( k_part_ );
     }
     FloatType usol()
     {
       return u_sol_;
     }
     FloatType upart()
     {
       return (u_part_ );
     }
     scitbx::sym_mat3<FloatType> ustar()
     {
       return ( u_star_ );
     }
     scitbx::af::shared< cctbx::miller::index<> > hkl()
     {
       return(hkl_);
     }

     //--------------
     void koverall(FloatType new_koverall){
       k_overall_=new_koverall;
       recompute_total_=true;
     }
     void ksol(FloatType new_ksol)
     {
       k_sol_=new_ksol ;
       recompute_total_=true;
     }
     void kpart(FloatType new_kpart)
     {
       k_part_=new_kpart;
       recompute_total_=true;
     }
     void usol(FloatType new_usol)
     {
       u_sol_=new_usol;
       recompute_bulk_=true;
       recompute_total_=true;

     }
     void upart(FloatType new_upart)
     {
       u_part_=new_upart;
       recompute_part_=true;
       recompute_total_=true;

     }
     void ustar( scitbx::sym_mat3<FloatType> new_ustar)
     {
       u_star_=new_ustar;
       recompute_aniso_=true;
       recompute_total_=true;

     }


     //--------------------------------
     // d_f_model_core_data_d_f_atoms
     scitbx::af::shared<FloatType> d_f_model_core_data_d_f_atoms()
     {
       return( overall_scale() );
     }



     // selection methods, takes in an iselection
     f_model_core_data<FloatType> select( scitbx::af::const_ref<int> const& iselection) const
     {
       scitbx::af::shared<std::complex<FloatType> > new_f_atoms;
       scitbx::af::shared<std::complex<FloatType> > new_f_mask;
       scitbx::af::shared<std::complex<FloatType> > new_f_part;
       scitbx::af::shared< cctbx::miller::index<> > new_hkl;
       scitbx::af::shared< FloatType > new_d_star_sq;


       CCTBX_ASSERT( iselection.size() <= hkl_.size() );
       for (std::size_t ii=0;ii<iselection.size();ii++){
         CCTBX_ASSERT( iselection[ii]<hkl_.size() );
         CCTBX_ASSERT( iselection[ii]>=0 );
         new_hkl.push_back( hkl_[iselection[ii]] );
         new_f_atoms.push_back( f_atoms_[iselection[ii]] );
         new_f_mask.push_back( f_mask_[iselection[ii]] );
         new_f_part.push_back( f_part_[iselection[ii]] );
         new_d_star_sq.push_back( d_star_sq_[iselection[ii]] );
       }

       f_model_core_data<FloatType> new_f_model_core_data( new_hkl.const_ref(),
                                       new_f_atoms.const_ref(),
                                       new_f_mask.const_ref(),
                                       new_d_star_sq.const_ref(), // pass in a d_star_sq in stead of unit_cell
                                       k_overall_, u_star_,
                                       k_sol_, u_sol_,
                                       new_f_part.const_ref(), k_part_, u_part_);
       return (new_f_model_core_data);
     }

     scitbx::af::shared<FloatType> d_star_sq()
     {
       return( d_star_sq_ );
     }


 protected:
    //----------------
    // Scale factors
    void compute_aniso_scale( std::size_t ii )
    {
      FloatType result = 0;
      result  = hkl_[ii][0]*(hkl_[ii][0]*u_star_[0]+hkl_[ii][1]*u_star_[3]+hkl_[ii][2]*u_star_[4]) +
                hkl_[ii][1]*(hkl_[ii][0]*u_star_[3]+hkl_[ii][1]*u_star_[1]+hkl_[ii][2]*u_star_[5]) +
                hkl_[ii][2]*(hkl_[ii][0]*u_star_[4]+hkl_[ii][1]*u_star_[5]+hkl_[ii][2]*u_star_[2]) ;
      result *= -2.0*scitbx::constants::pi*scitbx::constants::pi;
      result= std::exp(result);
      aniso_scale_[ii]=result;

    }
    void compute_bulk_scale( std::size_t ii)
    {
      FloatType result = 0;
      result = std::exp( - u_sol_*2.0*scitbx::constants::pi*scitbx::constants::pi*d_star_sq_[ii] );
      bulk_scale_[ii]=result;
    }

    void compute_part_scale( std::size_t ii )
    {
      FloatType result=0;
      result =  std::exp( - u_part_*2.0*scitbx::constants::pi*scitbx::constants::pi*d_star_sq_[ii] );
      part_scale_[ii]=result;
    }

    void compute_all_scales( std::size_t ii )
    {
      // aniso scale
      FloatType result = 0;
      if (recompute_aniso_){
        result  = hkl_[ii][0]*(hkl_[ii][0]*u_star_[0]+hkl_[ii][1]*u_star_[3]+hkl_[ii][2]*u_star_[4]) +
                  hkl_[ii][1]*(hkl_[ii][0]*u_star_[3]+hkl_[ii][1]*u_star_[1]+hkl_[ii][2]*u_star_[5]) +
                  hkl_[ii][2]*(hkl_[ii][0]*u_star_[4]+hkl_[ii][1]*u_star_[5]+hkl_[ii][2]*u_star_[2]) ;
        result *= -2.0*scitbx::constants::pi*scitbx::constants::pi;
        result= std::exp(result);
        aniso_scale_[ii]=result;
      }
      // bulk scale
      if (recompute_bulk_){
        bulk_scale_[ii] = std::exp(-u_sol_*d_star_sq_[ii]*2.0*scitbx::constants::pi*scitbx::constants::pi);
      }
      // partial scale
      if (recompute_part_){
        part_scale_[ii] = std::exp(-u_part_*d_star_sq_[ii]*2.0*scitbx::constants::pi*scitbx::constants::pi);
      }

    }

    //-----------------
    // compute f_model_core_data
    void compute_f_model_core_data( std::size_t ii )
    {
      f_model_core_data_[ii] = k_overall_*aniso_scale_[ii]*
        (                         f_atoms_[ii]
         + k_sol_*bulk_scale_[ii]*f_mask_[ii]
         + k_part_*part_scale_[ii]*f_part_[ii] );
    }

    // do it all
    void update_all( std::size_t ii )
    {
      compute_all_scales( ii );
      compute_f_model_core_data( ii );
    }


    mutable scitbx::af::shared< cctbx::miller::index<> > hkl_;
    mutable scitbx::af::shared<std::complex<FloatType> > f_atoms_;
    mutable scitbx::af::shared<std::complex<FloatType> > f_mask_;
    mutable scitbx::af::shared<std::complex<FloatType> > f_part_;

    mutable scitbx::af::shared< FloatType > d_star_sq_;

    mutable scitbx::af::shared< FloatType > aniso_scale_;
    mutable scitbx::af::shared< FloatType > log_part_aniso_scale_;
    mutable scitbx::af::shared< FloatType > bulk_scale_;
    mutable scitbx::af::shared< FloatType > part_scale_;

    mutable scitbx::af::shared< std::complex<FloatType> > f_model_core_data_;

    mutable FloatType k_overall_;
    mutable scitbx::sym_mat3<FloatType> u_star_;
    mutable FloatType k_sol_;
    mutable FloatType u_sol_;
    mutable FloatType k_part_;
    mutable FloatType u_part_;

    mutable bool recompute_aniso_;
    mutable bool recompute_bulk_;
    mutable bool recompute_part_;
    mutable bool recompute_total_;
 };

}}}

#endif // CCTBX_XRAY_TARGETS_H
