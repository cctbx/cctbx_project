#ifndef CCTBX_ELECTRON_SCATTERING_H
#define CCTBX_ELECTRON_SCATTERING_H

#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/xray_scattering.h>
#include <boost/optional.hpp>
#include <boost/config.hpp>
#include <stdexcept>
#include <ctype.h>

namespace cctbx { namespace eltbx { namespace electron_scattering {
  
  //! Coefficients for the Electron Scattering Factor.
  /*
   <p>
	 2008-06-30, translated by L. Lutterotti from
   L.-M. PENG, G. REN, S. L. DUDAREV & M. J. WHELAN 
   Robust Parameterization of Elastic and Absorptive Electron Atomic 
   Scattering Factors
   J. Appl. Cryst. A52, 257-276
   1996
   and
   L.-M. PENG 
   Electron Scattering Factors of Ions and their Parameterization
   J. Appl. Cryst. A54, 481-485 
   1998
   Electron scattering factors for atoms up to 6.0 Angstrom^-1
   <p>
   See also:
   electron_scattering::peng19965_iterator,
   electron_scattering::peng1996,
   xray_scattering::base
   */
  class peng1996: public xray_scattering::base<5>
  {
  public:
    //! Facilitates meta-programming.
    typedef xray_scattering::base<5> base_type;
    
    //! Default constructor. Calling certain methods may cause crashes!
    peng1996() {}
    
    //! Looks up coefficients for the given scattering factor label.
    /*! If exact == true, the scattering factor label must exactly
     match the tabulated label. However, the lookup is not
     case-sensitive.<br>
     E.g., "SI4" will be matched with "Si".<br>
     "Si4+" and "Si+4" will be matched with "Si4+".<br>
     See also: eltbx::basic::strip_label()
     <p>
     Note that the other methods of this class are inherited from
     class base.
     */
    peng1996(std::string const& label, bool exact=false);
    
  protected:
    friend class peng1996_iterator;
  };
  
  /*! \brief Iterator over table of Coefficients for the Analytical
   Approximation to the Scattering Factor, Waasmaier & Kirfel 1995.
   */
  class peng1996_iterator
  {
  public:
    //! Initialization of the iterator.
    peng1996_iterator();
    
    //! Retrieves the next entry from the internal table.
    /*! Use peng1996::is_valid() to detect end-of-iteration.
     */
    peng1996
    next();
    
  private:
    peng1996 current_;
  };
  
}}} // cctbx::eltbx::electron_scattering

#endif
