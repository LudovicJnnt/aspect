/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/material_model/crust1.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>

#include <numeric>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    Crust1<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          //const double temperature = in.temperature[i];
          const std::vector<double> composition = in.composition[i];

          out.specific_heat[i] = specific_heats[0];
          out.thermal_conductivities[i] = thermal_conductivities[0];

          // Determine compositional index (e.g., rock type)
          const double c = (composition.size()>0)
                           ?
                           std::max(0.0, std::round(composition[0]))
                           :
                           0.0;
          unsigned int com = 0;
          for (unsigned int j=0; j < densities.size(); ++j)
            {
              if ((unsigned int)c == j)
                {
                  com = j;
                }
            }

          const double thermal_alpha = thermal_expansivities[0];

          const double depth = this->get_geometry_model().depth(in.position[i]);

          // Density is the second compositional field
          if ( depth < crust1_depth )
            {
              out.densities[i] = ( (in.current_cell.state() == IteratorState::valid)
                                 ?
                                 std::max(composition[1],10.0) 
                                 :
                                 densities[com]);
            }
          else
            {
              out.densities[i] = reference_rho * (1 - thermal_alpha * (in.temperature[i] - reference_T));
            }

  
          out.viscosities[i] = viscosities[com];
 
          out.thermal_expansion_coefficients[i] = thermal_alpha;

          // Compressibility at the given positions.
          // The compressibility is given as
          // $\frac 1\rho \frac{\partial\rho}{\partial p}$.
          // (here we use an incompressible medium)
          out.compressibilities[i] = 0.0;
          // Pressure derivative of entropy at the given positions.
          out.entropy_derivative_pressure[i] = 0.0;
          // Temperature derivative of entropy at the given positions.
          out.entropy_derivative_temperature[i] = 0.0;
          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

        }
    }

    template <int dim>
    double
    Crust1<dim>::
    reference_viscosity () const
    {
      return viscosities[0]; //background
    }

    template <int dim>
    double
    Crust1<dim>::
    reference_density () const
    {
      return densities[0];  //background
    }

    template <int dim>
    double
    Crust1<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_expansivities[0]; //background
    }

    template <int dim>
    double
    Crust1<dim>::
    reference_cp () const
    {
      return specific_heats[0]; //background
    }

    template <int dim>
    double
    Crust1<dim>::
    reference_thermal_diffusivity () const
    {
      return thermal_conductivities[0] /( densities[0]* specific_heats[0] ); //background
    }

    template <int dim>
    bool
    Crust1<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    Crust1<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Crust1");
        {
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Crust1 base depth", "120.e3",
                             Patterns::Double (0),
                             "Depth from model surface to base of defined crust1 density . Units: m.");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Viscosities", "1.e21",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $Pa s$");
          prm.declare_entry ("Thermal expansivities", "4.e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $1/K$");
          prm.declare_entry ("Specific heats", "1250.",
                             Patterns::List(Patterns::Double(0)),
                             "List of specific heats $C_p$ for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $J /kg /K$");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value. Units: $W/m/K$ ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Crust1<dim>::parse_parameters (ParameterHandler &prm)
    {
      //not pretty, but we need to get the number of compositional fields before
      //simulatoraccess has been initialized here...
      unsigned int n_foreground_fields;
      prm.enter_subsection ("Compositional fields");
      {
        n_foreground_fields = prm.get_integer ("Number of fields");
      }
      prm.leave_subsection();

      const unsigned int n_fields= n_foreground_fields + 1;
      const unsigned int d_fields= 10;


      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Crust1");
        {
          reference_rho = prm.get_double ("Reference density");
          reference_T   = prm.get_double ("Reference temperature");
          crust1_depth  = prm.get_double ("Crust1 base depth");
          
          // Parse multicomponent properties
          densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                              d_fields,
                                                              "Densities");
          viscosities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosities"))),
                                                                d_fields,
                                                                "Viscosities");
          thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                                                                           n_fields,
                                                                           "Thermal conductivities");
          thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),
                                                                          n_fields,
                                                                          "Thermal expansivities");
          specific_heats = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Specific heats"))),
                                                                   n_fields,
                                                                   "Specific heats");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::compositional_fields;
      this->model_dependence.thermal_conductivity = NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Crust1,
                                   "crust1",
                                   "This model is desgined specifically to be used when CRUST1.0 defines"
                                   " the compositional initial conditions.  Currently, the model only uses"
                                   " one compositional field where the value at each point represents a"
                                   " specific rock type in the CRUST1.0 model.  The user should specify 10"
                                   " density values in the input file, which represent air, water, ice,"  
                                   " sed layer 1, sed layer 2, sed layer 3, upper crust, middle crust,"
                                   " lower crust and mantle.  At some point soon, density will be determined"
                                   " from a separate compositional field that contains the CRUST1.0 densities" 
                                   " read in from a file.")
  }
}
