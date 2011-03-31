//! Peter Zwart April 19, 2005
#ifndef CCTBX_MILLER_LOOKUP_UTILS_H
#define CCTBX_MILLER_LOOKUP_UTILS_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/match_indices.h>
#include <cctbx/miller/asu.h>
#include <cctbx/sgtbx/reciprocal_space_asu.h>
#include <cctbx/sgtbx/space_group_type.h>

#include <map>
#include <vector>

namespace cctbx { namespace miller{
namespace lookup_utils{


  template <typename FloatType=double>
  class lookup_tensor
  {
    protected:
      typedef std::map<
        cctbx::miller::index<>,
        std::size_t,
        cctbx::miller::fast_less_than<> > lookup_map_type;

      int n_duplicates_;
      int n_indices_;
      cctbx::sgtbx::space_group space_group_;
      cctbx::sgtbx::space_group_type sgtype_;
      cctbx::sgtbx::reciprocal_space::asu asu_choice_;
      lookup_map_type hkl_lookup_;
      bool anomalous_flag_;
      public:

    /*! Default constructor */
    lookup_tensor() {}

    /*! Constructor with a given list of of HKL's */
    lookup_tensor(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
      cctbx::sgtbx::space_group const& space_group,
      bool const& anomalous_flag)
    :
      n_duplicates_( 0 ),
      n_indices_( hkl.size() ),
      space_group_(space_group),
      sgtype_(space_group_),
      asu_choice_(sgtype_),
      hkl_lookup_(),
      anomalous_flag_(anomalous_flag)
    {
      for (unsigned ii=0;ii<hkl.size();ii++){
        cctbx::miller::asym_index asumap(space_group_,
                                         asu_choice_,
                                         hkl[ii]);
        cctbx::miller::index_table_layout_adaptor asu_target_hkl;
        asu_target_hkl = asumap.one_column(anomalous_flag_);

        lookup_map_type::const_iterator l = hkl_lookup_.find( asu_target_hkl.h() );
        if (l == hkl_lookup_.end()) { // not in list
          hkl_lookup_[ asu_target_hkl.h()  ] = ii;
        } else {
          n_duplicates_++;
        }
        // checking if we can find it back would be shear paranoia
        //l = hkl_lookup_.find( asu_target_hkl.h() );
        //CCTBX_ASSERT( l  != hkl_lookup_.end() );
      }
    }

    long
    n_duplicates() const
    {
       return( n_duplicates_ );
    }

    long
    find_hkl( cctbx::miller::index<> const& target_hkl) const
    {
      long hkl_location;
      // move the index to the ASU
      cctbx::miller::asym_index asumap(space_group_,
                                       asu_choice_,
                                       target_hkl);
      cctbx::miller::index_table_layout_adaptor asu_target_hkl;
      asu_target_hkl = asumap.one_column(anomalous_flag_);
      lookup_map_type::const_iterator
        l = hkl_lookup_.find( asu_target_hkl.h() );
      if (l == hkl_lookup_.end()) {
        hkl_location = -1; // !!! negative if not found !!!
      }
      else {
        hkl_location = l->second;
      }
      if (hkl_location >= n_indices_){
        hkl_location = -1;
      }
      return (hkl_location);
    }

    scitbx::af::shared<long>
    find_hkl( scitbx::af::const_ref<cctbx::miller::index<> >
              const& target_hkl ) const
    {
      scitbx::af::shared< long > permutation_table( target_hkl.size(), -1 );
      for (unsigned ii=0;ii<target_hkl.size();ii++){
        permutation_table[ii] = find_hkl(target_hkl[ii]);
      }
      return(permutation_table);
    }
  };


  //-------------------
  // This class returns a local neighbourhood of a specific reflection
  //
  // Note that duplicates are NOT removed!
  //
  template <typename FloatType=double>
  class local_neighbourhood
  {
    public:
    /*! Default consructor */
    local_neighbourhood(){}
    /*! Constructor with all info that is needed */
    local_neighbourhood(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
      cctbx::sgtbx::space_group const& space_group,
      bool const& anomalous_flag,
      long const& radius
    )
    :
    lookup_tensor_(hkl, space_group, anomalous_flag),
    radius_(radius)
    {
      SCITBX_ASSERT( hkl.size() > 0 );
      for (unsigned ii=0;ii<hkl.size();ii++){
        hkl_.push_back( hkl[ii] );
      }
    }

    //! Access via specific miller nidex
    std::vector< unsigned >
    construct_neighbourhood( cctbx::miller::index<> const& center_hkl )
    {
      std::vector<unsigned> neighbourhood_lookup_tensor;
      long manhattan_distance;
      long ht,kt,lt,index;
      for (int hh = -radius_ ; hh< radius_ +1; hh++){
        for (int kk = -radius_ ; kk< radius_ +1; kk++){
          for (int ll = -radius_ ; ll< radius_ +1; ll++){

            manhattan_distance = std::abs(hh) +
                                 std::abs(kk) +
                                 std::abs(ll);
            if (manhattan_distance<=radius_){
              if (manhattan_distance>0){
                ht = center_hkl[0] + hh;
                kt = center_hkl[1] + kk;
                lt = center_hkl[2] + ll;
                cctbx::miller::index<> picked_neighbour(ht,kt,lt);
                index = lookup_tensor_.find_hkl(picked_neighbour);
                if (index >=0){
                  neighbourhood_lookup_tensor.push_back( unsigned(index) );
                }
              }
            }


          }
        }
      }

      return(neighbourhood_lookup_tensor);
    }

    // Acces via position in list
    std::vector< unsigned >
    construct_neighbourhood( unsigned const& center_hkl )
    {
      SCITBX_ASSERT( hkl_.size() > center_hkl);
      return( construct_neighbourhood(hkl_[center_hkl]) );
    }



    // Array version, give routnie list of miller indices please
    scitbx::af::shared< std::vector< unsigned >  >
    construct_neighbourhood(
      scitbx::af::shared< cctbx::miller::index<> > const& hkl )
    {
      scitbx::af::shared< std::vector< unsigned > > result;
      for(unsigned ii=0;ii<hkl.size();ii++){
        result.push_back( construct_neighbourhood( hkl[ii] ) );
      }
      return(result);
    }

    // Array version, give routine list of access numbers please
    scitbx::af::shared< std::vector< unsigned >  >
    construct_neighbourhood(
      scitbx::af::shared< long > const& hkl )
    {
      scitbx::af::shared< std::vector<unsigned>  > result;
      for (unsigned ii=0;ii<hkl.size();ii++){
        if ( hkl[ii]>=0 ){
          result.push_back( construct_neighbourhood( hkl[ii] ) );
        }
        else{
          std::vector<unsigned> tmp;
          result.push_back(  tmp  );
        }
      }
      return( result );

    }

    // Use the stored hkl_ array to generate a list of neighbours
    // Do not give anything
    scitbx::af::shared< std::vector< unsigned >  >
    construct_neighbourhood()
    {
      scitbx::af::shared< std::vector<unsigned>  > result;
      for (unsigned ii=0;ii<hkl_.size();ii++){
        result.push_back( construct_neighbourhood( hkl_[ii] ) );
      }
      return( result );

    }

    protected:
    lookup_tensor<FloatType> lookup_tensor_;
    scitbx::af::shared< cctbx::miller::index<> > hkl_;

    long radius_;

  };

  //--------------------------------------------------------------
  // this class allows the fast computation of
  // a neighbourhood of any radius

  template <typename FloatType=double>
  class neighbourhood_list
  {
    public:
    /*! Default constructor */
    neighbourhood_list(){}
    /*! Constructor with all info */
    neighbourhood_list(
      scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
      cctbx::sgtbx::space_group const& space_group,
      bool const& anomalous_flag,
      std::size_t const& radius)
    :
    generator_(hkl,space_group,anomalous_flag,radius)
    {

      for (std::size_t ii=0;ii<hkl.size();ii++){
        // for each miller index, we have to find the neighbour set
        nb_list_.push_back(
          generator_.construct_neighbourhood(  hkl[ii] )
        );

      }

    }

    scitbx::af::shared< std::vector<unsigned> >
    neighbour_list(){ return(nb_list_); }

    std::vector<unsigned>
    neighbour_list( unsigned const& this_hkl )
    {
      return( nb_list_[ this_hkl ] );
    }

    protected:
    local_neighbourhood<FloatType> generator_;
    scitbx::af::shared< std::vector<unsigned>  > nb_list_;

  };


  //----------------------------------------------------------
  // get a local area.
  // Duplicate removal is performed.
  // An extra property array need to be provided
  // indices for which property = True, will be listed as possible neighbours
  // missing neighbours however are used in building up the
  // local area

  template <typename FloatType = double >
  class local_area
  {
    public:
    local_area(){}
    local_area(
               scitbx::af::const_ref< cctbx::miller::index<> > const& hkl,
               scitbx::af::const_ref< bool > const& property,
               cctbx::sgtbx::space_group const& space_group,
               bool const& anomalous_flag,
               std::size_t const& radius,
               std::size_t const& depth,
               std::size_t const& at_least_this_number_of_neighbours)
    :
    depth_(depth),
    local_area_locator_( hkl, space_group, anomalous_flag, radius),
    used_( hkl.size(), 0 ),
    average_number_of_neighbours_( 0 )
    {
      // make sure the 'property' array has the same size as the HKL array
      SCITBX_ASSERT( property.size() == hkl.size() );

      // store a neighbour list please
      neighbour_list_ = local_area_locator_.construct_neighbourhood();



      unsigned last_size=0, tmp_size;
      for (unsigned ii=0;ii<hkl.size();ii++){

        bool done_searching=true;
        unsigned iter=0;
        unsigned count=0;
        std::vector<unsigned> tmp_vector;
        area_.push_back( tmp_vector );
        //place the pivot point if allowed by the property array
        if (property[ii]){
          tmp_vector.push_back(ii);
          // flag used please
          used_[ii]=1;
          done_searching=false;
        }
        last_size =  tmp_vector.size();
        // now do depth iterations


        while ( !done_searching ){
          tmp_size = tmp_vector.size();
          // loop over all elements please that we have not looked at before
          for (unsigned item=last_size-1;item<tmp_size;item++){
            std::vector<unsigned> tmp_nbl;
            tmp_nbl = neighbour_list_[ tmp_vector[item] ];
            for (unsigned kk=0;kk<tmp_nbl.size();kk++){
              // a potential new neighbour has been found
              // please check if it has been used before
              if (used_[ tmp_nbl[kk] ] == 0){
                // this has not been used before
                // add it to the list.
                tmp_vector.push_back( tmp_nbl[kk] );
                used_[ tmp_nbl[kk] ] = 1;
                if ( property[ tmp_nbl[kk] ] ){
                  count++;
                }
              }
            }
          }
          // set the last_size to what it should be
          last_size = tmp_size;
          // check if we are done
          if (iter>=depth-1){
            done_searching = true;
          }
          if (count>= at_least_this_number_of_neighbours){
            done_searching = true;
          }
          // update iteration number
          iter++;
        } // end of iterations/depth search

        average_number_of_neighbours_+=count;
        // we now have to set the used_ array back to zero again
        // as well as copy the elements into the area_ array
        for ( unsigned jj=0;jj<tmp_vector.size();jj++){
          // take care of used flags
          used_[ tmp_vector[jj] ] = 0;
          // append if property flag allows
          if (property[  tmp_vector[jj] ] ){
            area_[ii].push_back( tmp_vector[jj] );
          }
        }

      }
      average_number_of_neighbours_/=FloatType( hkl.size() );
    }

    // make the area publicly accessable.


    scitbx::af::shared< std::vector<unsigned> >
    area()
    {
      return( area_ );
    }

    scitbx::af::shared< std::vector<unsigned> > area_;

    protected:
    unsigned depth_;
    local_neighbourhood<FloatType> local_area_locator_;
    scitbx::af::shared< std::vector<unsigned> > neighbour_list_;
    scitbx::af::shared<int> used_;
    FloatType average_number_of_neighbours_;

  };








}}} // namespace cctbx::miller::lookup_utils
#endif // CCTBX_MILLER_LOOKUP_UTILS_H
