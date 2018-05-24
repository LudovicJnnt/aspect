/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_postprocess_satellite_gravity_h
#define _aspect_postprocess_satellite_gravity_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes gravity and potential for a set of points (e.g. satellites)
     * above the model surface for a user-defined range of latitudes, longitudes and height. 
     * A gravity depth-profile for a user-defined longitude and latitude is also possible. 
 
     * @ingroup Postprocessing
     */
    template <int dim>
    class GravityCalculation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Output gravity [m] to file
         */
        virtual
        std::pair<std::string,std::string> execute (TableHandler &);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

      private:
        /**
         * Gravity is calculated on a spherical surface (i.e. map at a height) at which points 
         * are spaced according to the grid_spacing value in degree.
         */
        double mapping_spherical_spacing;
        double profile_height_spacing;

        /**
         * Height above the model surface (e.g., 'satellite height') at which the gravity is calculated.
         */
        double mapping_satellite_height;
        double profile_satellite_height;

        /**
         * Longitude at which gravity depth-profile is calculated 
         */
        double profile_longitude;

        /**
         * Latitude at which gravity depth-profile is calculated 
         */
        double profile_latitude;
        /**
         * Minimum longitude required for regional gravity mapping 
         */
        double minimum_longitude;

        /**
         * Minimum longitude required for regional gravity mapping 
         */
        double maximum_longitude;
        
        /**
         * Minimum longitude required for regional gravity mapping 
         */
        double minimum_latitude;
    
        /**
         * Minimum longitude required for regional gravity mapping 
         */
        double maximum_latitude;

        /**
         * The reference density used in the gravity calculation 
         */
        double reference_density;
    };
  }
}


#endif
