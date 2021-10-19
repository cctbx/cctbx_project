// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissionsand
// limitations under the License.

#pragma once

#include <string>
#include <map>
#include "Common.h"

namespace molprobity {
  namespace probe {

    /// @brief Structure that holds a vector of dots surrounding a sphere
    ///
    /// This structure stores a vector of dots surrounding a sphere with a specified
    /// radius and density.  Once constructed, its dots and parameters cannot be changed.
    class DotSphere {
    public:
      /// @brief Normally-used constructor that fills in dots based on parameters
      /// @param [in] radius The radius of the sphere in angstroms, distance from the origin to dots.
      ///         If this is zero or negative, the object will have no dots.
      /// @param [in] density The density of the dots in square angstroms.
      ///         If this is zero or negative, the object will have no dots.
      ///         This specifies a desired density; the actual density may be slightly different.
      DotSphere(double radius, double density);

      /// @brief Returns the vector of dots on the sphere
      scitbx::af::shared<Point> const &dots() const { return m_vec; }
      scitbx::af::shared<Point> dotsCopyForPythonWrapping() const { return m_vec; }

      /// @brief Accessor method for the radius used to construct the sphere.
      double radius() const { return m_rad; }

      /// @brief Accessor method for the density used to construct the sphere.
      ///
      /// This is used to normalize the scores to make them comparable for calculations
      /// done with different densities.
      double density() const { return m_dens; }

      /// @brief Report whether two spheres have been constructed with identical parameters.
      bool operator==(const DotSphere& d) const { return m_rad == d.m_rad && m_dens == d.m_dens; }

      /// @brief Report whether two spheres have been constructed with different parameters.
      bool operator!=(const DotSphere& d) const { return !(*this == d); }

      //===========================================================================
      // Seldom-used methods below here.

      /// @brief Default constructor makes an empty vector of dots
      DotSphere() : m_rad(0), m_dens(0) {};

      /// @brief Test method to verify that the class is behaving as intended.
      /// @return Empty string on success, string telling what went wrong on failure.
      static std::string test();

    protected:
      double  m_rad;            ///< Radius of the sphere in Angstroms
      double  m_dens;           ///< Dot density
      scitbx::af::shared<Point> m_vec;   ///< Stores the vector of dots
    };

    /// @brief Constructs new DotSphere objects as needed, reusing existing ones when it can.
    ///
    /// This structure minimizes the number of DotSphere objects required for a number of
    /// atoms by returning references to already-constructed spheres when asked for another
    /// with a radius that it has already constructed.  All the spheres generated have the
    /// same density.
    class DotSphereCache {
    public:
      /// @brief Constructor for spheres of a given density.
      /// @param [in] Density value for all spheres constructed by this map
      DotSphereCache(double density) : m_dens(density) {}

      /// @brief Return a reference to a sphere of specified radius, building if needed.
      /// @param [in] radius The radius of the sphere
      /// @return Reference to a sphere of the specified radius
      const DotSphere& get_sphere(double radius);

      /// @brief Return the total number of unique spheres created.
      /// @return The number of unique spheres in the cache
      size_t size() const { return m_spheres.size(); }

      //===========================================================================
      // Seldom-used methods below here.

      /// @brief Test method to verify that the class is behaving as intended.
      /// @return Empty string on success, string telling what went wrong on failure.
      static std::string test();

    protected:
      double  m_dens;                         ///< Density of the constructed spheres
      std::map<double, DotSphere> m_spheres;  ///< Already-constructed spheres by radius
    };

    /// @brief Test function to verify that all classes are behaving as intended.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string DotSpheres_test();
  }
}
