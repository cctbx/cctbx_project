// $Id$

#ifndef CCTBX_ELTBX_ICSD_RADII_H
#define CCTBX_ELTBX_ICSD_RADII_H

#include <string>

namespace eltbx {

  namespace detail {
    struct Label_Radius {
      const char* Label;
      float       Radius;
    };
  }

  //! Access to table of ionic radii.
  /*! Reference:<pre>
                          U s e r ' s  M a n u a l
                         I C S D  -  C R Y S T I N
                         =========================
                    Inorganic Crystal Structure Database
                            in conjunction with
                   Crystal Structure Information System
                           and its application to
                CCDF - Cambridge Crystallographic Data File
                                    and
                           MDF - Metal Data File
                    G.Bergerhoff, B.Kilger, C.Witthauer,
                            R.Hundt, R.Sievers
                                    Bonn
                                    1986
             Institut fuer Anorganische Chemie der Universitaet
            -------------------------------------------------------------
            ICSD/CRYSTIN User's Manual.
            English Version. Translated by Ruth Schubert.
            Updated Dec. 1986.</pre>

      These radii are also the basis of the distance tests,  which are
      routinely  carried  out when collecting data in the ICSD system.
      In  this  connection,   negative increments are obtained for the
      specially small atoms C+4, D+1, H+1 and N+5.
   */
  class ICSD_Radius {
    public:
      //! Lookup ionic radius for the given ion label.
      /*! If Exact == true, the ion label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.<br>
          E.g., "SI4" will be matched with "Si".<br>
          "Si4+" and "Si+4" will be matched with "Si4+".<br>
          See also: eltbx::StripLabel()
       */
      ICSD_Radius(const std::string& Label, bool Exact = false);
      //! Return label from table.
      inline const char* Label() const { return m_Label_Radius->Label; }
      //! Return radius (Angstrom) from table.
      inline float Radius() const { return m_Label_Radius->Radius; }
    private:
      const detail::Label_Radius* m_Label_Radius;
  };

} // eltbx

#endif // CCTBX_ELTBX_ICSD_RADII_H
