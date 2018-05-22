/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_gravity_calculation_h
#define _aspect_postprocess_gravity_calculation_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that outputs the gravity to file.

     * @ingroup Postprocessing
     */
    template <int dim>
    class GravityCalculation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Output topography [m] to file
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
         * Gravity is calculated on a spherical surface at which points are spaced
         * according to the spherical_spacing values.
         */
        double spherical_spacing;

        /**
         * Gravity is calculated at a certain and constant height.
         */
        double satellite_height;

        /**
         * Gravity difference calculation - theoric. Theoric gravity requires a 
         * reference density of the all domain.
         */
        double reference_density;

        std::vector<double> long_min_max;
        std::vector<double> lat_min_max;
    };
  }
}


#endif
