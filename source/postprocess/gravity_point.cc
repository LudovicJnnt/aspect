/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPEC is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/postprocess/gravity_calculation.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/global.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

/*
  This postprocessor computes gravity and potential regionally or globally for 
  a set of points located at ecoord[r,phi,theta] in the model (map).
  
  Because of the parallel nature of this postprocessor, first it is required to 
  store density values from MaterialModel in a vector per MPI. Doing so avoids to 
  use MaterialModel to get the density at quadrature points within the loops. 
  
  The triple integral goes first over longitude and latitude at previously specified
  height and then the spherical coordinates are shifted into cartesian to allow
  simplification in the mathematical equation.
  
  For each point (i.e. satellite), the third integral goes over cells and quadrature 
  points to get the unique distance between those, indispensable to calculate 
  gravity vector components x,y,z, and potential. 
*/

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    GravityCalculation<dim>::execute (TableHandler &)
    {
      Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != 0 ||
              dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()) != 0,
              ExcMessage ("This postprocessor can only be used if the geometry "
                          "is a sphere or spherical shell."));

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);

      // Get the value of the outer radius and inner radius
      const double Ro = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
                                (this->get_geometry_model()).outer_radius();
      const double Ri = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
                                (this->get_geometry_model()).inner_radius();

      // Get the value of the universal gravitational constant
      const double G = aspect::constants::big_g;
      const double density_ref = reference_density;
      std_cxx11::array<double,dim> ecoord;       // create spherical coordinate array "ecoord"
      ecoord[0]= satellite_height;               // satellite height ecoord[radius, , ]
      const double spacing  = spherical_spacing; // satellite cover resolution
      const double long_min = long_min_max[0];
      const double long_max = long_min_max[1];
      const double lat_min  = lat_min_max[0];
      const double lat_max  = lat_min_max[1];

      // now write all data to the file of choice. start with a pre-amble that
      // explains the meaning of the various fields
      const std::string filename = (this->get_output_directory() +
                                    "gravity_arrays.txt");
      std::ofstream f (filename.c_str());
      f << "#1 position_satellite_r" << '\n'
        << "#2 position_satellite_phi" << '\n'
        << "#3 position_satellite_theta" << '\n' 
        << "#4 position_satellite_x" << '\n'
        << "#5 position_satellite_y" << '\n'
        << "#6 position_satellite_z" << '\n' 
        << "#7 gravity_x" << '\n'
        << "#8 gravity_y" << '\n'
        << "#9 gravity_z" << '\n'
        << "#10 gravity_norm" << '\n'
        << "#11 gravity_theory" << '\n'
        << "#12 gravity_diff (theoric-norm)" << '\n'
        << "#13 potential" << '\n'
        << '\n';

      // Temporary output during the build of this postprocessor
      const std::string filename2 = (this->get_output_directory() +
                                    "checks.txt");
      std::ofstream f2 (filename2.c_str());
      f2 << long_min << ' ' << long_max << '\n'
         << lat_min << ' ' << lat_max << '\n'
         << spacing << '\n'
         << (long_max-long_min)/spacing << '\n'
         << ecoord[0] << ' ' << ecoord[1] << ' ' << ecoord [2] << '\n'
         << '\n';

      // Total number of element per MPI
      // use for storing density values in a vector per mpi
      int c = 0;
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
            c += 1 ;

      // Storage of density values in a vector density_all per mpi 
      int cell_n_q = c * n_q_points;
      std::vector<double> density_all;
      density_all.reserve(cell_n_q); 

      c = 0;
      cell = this->get_dof_handler().begin_active();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            MaterialModel::MaterialModelInputs<dim> in(fe_values, cell, this->introspection(), this->get_solution());
            MaterialModel::MaterialModelOutputs<dim> out(quadrature_formula.size(),this->n_compositional_fields());
            this->get_material_model().evaluate(in, out);
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
               density_all[c*n_q_points+q]=out.densities[q];
              }
            c += 1;
          }
 
      // Two cases for regional calculation are possible: 
      //   case 1  - long 000-360: normal case working as globally
      //   case 2  - long 360-000: requires separating case 2 into two case 1
      //                           OR
      //                           requires longitude shift?
      // if ((long_max-long_min) < 0)
      
      // Note that the order of spherical coordinates is r,phi,theta and not r,theta,phi since this allows
      // for dimension independent expressions. 

      // loop on phi -  satellite position [ , phi ,]
      for (int i=0; i<(long_max-long_min)/spacing; i++)
        {
          ecoord[1]=long_min+i*spacing;

          // loop on theta - satllite position [ , , theta]
          for (int j=0; j<(lat_max-lat_min)/spacing; j++)
            { 
              ecoord[2]=lat_min+j*spacing;
              // f2 << ecoord[0] << ' ' << ecoord[1] << ' ' << ecoord [2] << ' ' << i << ' ' << j << '\n';
              const Point<dim> position_satellite = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(ecoord);
              cell = this->get_dof_handler().begin_active();
     
              double local_gx = 0;
              double local_gy = 0;
              double local_gz = 0;
              double gnorm = 0;
              double gtheory = 0;      
              double g_diff = 0;      
              double local_U = 0;
              c = 0;

              for (; cell!=endc; ++cell)
                 if (cell->is_locally_owned())
                    {
                      fe_values.reinit (cell);
                      const std::vector<Point<dim> > &position_point = fe_values.get_quadrature_points();

                      for (unsigned int q = 0; q < n_q_points; ++q)
                        {
                          const double dist=(position_satellite[0]-position_point[q][0])*(position_satellite[0]-position_point[q][0]) +  
                                         (position_satellite[1]-position_point[q][1])*(position_satellite[1]-position_point[q][1]) +  
                                         (position_satellite[2]-position_point[q][2])*(position_satellite[2]-position_point[q][2]);  

                          double KK = G * density_all[c*n_q_points+q] * fe_values.JxW(q) / pow(std::sqrt(dist),3);
                          local_gx += KK * (position_satellite[0]-position_point[q][0]);  
                          local_gy += KK * (position_satellite[1]-position_point[q][1]);  
                          local_gz += KK * (position_satellite[2]-position_point[q][2]);
                          local_U -= G * density_all[c*n_q_points+q] * fe_values.JxW(q) / std::sqrt(dist);
                        }
                      c += 1;
                    }

              // assemble gravity component results from the mpi
              const double global_gx
                = Utilities::MPI::sum (local_gx, this->get_mpi_communicator());
              const double global_gy
                = Utilities::MPI::sum (local_gy, this->get_mpi_communicator());
              const double global_gz
                = Utilities::MPI::sum (local_gz, this->get_mpi_communicator());
              const double global_U
                = Utilities::MPI::sum (local_U, this->get_mpi_communicator());

              // calculate the gravity norm for each satellite position
              gnorm = std::sqrt((global_gx*global_gx)+  
                                (global_gy*global_gy)+  
                                (global_gz*global_gz));

              // analytical solution to estimate a mistfit according to satellite height 
              // can only be used if concentric density profile
              if (ecoord[0] <= Ri)
                {
                 gtheory = 0;
                 g_diff = -gnorm;
                }
              else if ((ecoord[0] > Ri) && (ecoord[0] < Ro))
                {
                 gtheory = G * numbers::PI * 4/3 * density_ref * (ecoord[0] - (pow(Ri,3) / pow(ecoord[0],2)));
                 g_diff = (gtheory - gnorm);
                }
              else
                {
                 gtheory = G * numbers::PI * 4/3 * density_ref * (pow(Ro,3)-pow(Ri,3)) / pow(ecoord[0],2);
                 g_diff = (gtheory - gnorm);
                }

              // write output
              f << ecoord[0] << ' '
                << ecoord[1] << ' '
                << ecoord[2] << ' '
                << position_satellite[0] << ' '
                << position_satellite[1] << ' '
                << position_satellite[2] << ' '
                << global_gx << ' '
                << global_gy << ' '
                << global_gz << ' '
                << gnorm << ' '
                << gtheory << ' '
                << g_diff << ' '
                << global_U
                << '\n';

            }
        }
      return std::pair<std::string, std::string> ("gravity calculation file:",
                                                  filename);
    }

    template <int dim>
    void
    GravityCalculation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Postprocess");
      {
        prm.enter_subsection ("Gravity calculation");
        {
          prm.declare_entry ("Reference density", "1",
                             Patterns::Double ());
          prm.declare_entry ("Spherical spacing", "1",
                             Patterns::Double ());
          prm.declare_entry ("Satellite height", "0",
                             Patterns::Double ());
          prm.declare_entry ("Boundary longitude", "0, 360",
                             Patterns::List(Patterns::Double ()));
          prm.declare_entry ("Boundary latitude", "0, 180",
                             Patterns::List(Patterns::Double ()));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    GravityCalculation<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Postprocess");
      {
        prm.enter_subsection ("Gravity calculation");
        {
	  reference_density = prm.get_double ("Reference density");
	  spherical_spacing = prm.get_double ("Spherical spacing");
          satellite_height  = prm.get_double ("Satellite height");
          long_min_max      = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double (Utilities::split_string_list(prm.get("Boundary longitude"))),
                                                                      2, "Boundary longitude"); 
          lat_min_max       = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double (Utilities::split_string_list(prm.get("Boundary latitude"))),
                                                                      2, "Boundary latitude"); 
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(GravityCalculation,
                                  "gravity calculation",
                                  "A postprocessor that computes gravity "
                                  "at a point.")
  }
}
