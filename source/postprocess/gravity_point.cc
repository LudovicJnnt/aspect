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
      std_cxx11::array<double,dim> ecoord;       // create spherical coordinate array
      ecoord[0]= satellite_height;               // satellite height [radius, , ]
      const double spacing = spherical_spacing;  // satellite cover resolution

      // now write all data to the file of choice. start with a pre-amble that
      // explains the meaning of the various fields
      const std::string filename = (this->get_output_directory() +
                                    "gravity_arrays.txt");
      std::ofstream f (filename.c_str());
      f << "#1 position_satellite_r" << '\n'
        << "#2 position_satellite_theta" << '\n'
        << "#3 position_satellite_phi" << '\n' 
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

      // number of element per MPI
      int c = 0;
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
            c += 1 ;

      // allocate density vector in MPI         // both I am not sure they work
      int cell_n_q = c * n_q_points;
      std::vector<double> density_all;          // defining vector density_all to store density values
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
 
      // loop on phi -  satellite position [ , ,phi]
      for (int i=0; i<360/spacing; i++)
        {
          ecoord[2]=i*spacing;

          // loop on theta - satllite position [ , theta, ]
          for (int j=0; j<180/spacing; j++)
            { 
              ecoord[1]=j*spacing;
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
              // POINTLESS WITH THE CURRENT VERSION USING S40RTS for initial condition.
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
          satellite_height = prm.get_double ("Satellite height");
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
