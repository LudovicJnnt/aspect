/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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

#include <aspect/postprocess/satellite_gravity.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/adiabatic_conditions/interface.h>
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
    GravityPointValues<dim>::execute (TableHandler &)
    {
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+1);
      const unsigned int number_quadrature_points_cell = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);

      // Get the value of the outer radius and inner radius
      // WHAT IF FULL SPHERE ???
      double model_outer_radius;
      double model_inner_radius;
      if (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != 0)
        {
          model_outer_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
                               (this->get_geometry_model()).outer_radius();
          model_inner_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
                               (this->get_geometry_model()).inner_radius();
        }
      else if (dynamic_cast<const GeometryModel::Chunk<dim>*> (&this->get_geometry_model()) != 0)
        {
          model_outer_radius = dynamic_cast<const GeometryModel::Chunk<dim>&>
                               (this->get_geometry_model()).outer_radius();
          model_inner_radius = dynamic_cast<const GeometryModel::Chunk<dim>&>
                               (this->get_geometry_model()).inner_radius();
        }
      else if (dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()) != 0)
        {
          model_outer_radius = dynamic_cast<const GeometryModel::Sphere<dim>&>
                               (this->get_geometry_model()).radius();
          model_inner_radius = 0; 
        }
      else
        {
          Assert (false, ExcMessage ("This initial condition can only be used if the geometry "
                                     "is a sphere, a spherical shell, a chunk or an "
                                     "ellipsoidal chunk."));
          model_outer_radius = 1;
          model_inner_radius = 0;
        }

      // Get the value of the universal gravitational constant
      const double G = aspect::constants::big_g;
      
      // now write all data to the file of choice. start with a pre-amble that
      // explains the meaning of the various fields
      const std::string filename = (this->get_output_directory() +
                                    "output_gravity.txt");
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

      // Temporary output during the build of this postprocessor WIP
      const std::string filename2 = (this->get_output_directory() +
                                    "checks_wip.txt");
      std::ofstream f2 (filename2.c_str());
      f2 << minimum_radius << ' ' << maximum_radius << '\n'
         << minimum_longitude << ' ' << maximum_longitude << '\n'
         << minimum_latitude << ' ' << maximum_latitude << '\n';

      // Storing cartesian coordinate, density and JzW at quadrature points in a vector 
      // avoids to use MaterialModel and fe_values within the loops.
      unsigned int local_cell_number = (this->get_triangulation().n_locally_owned_active_cells());
      const unsigned int number_quadrature_points_mpi = local_cell_number * number_quadrature_points_cell;
      std::vector<double> density_JxW (number_quadrature_points_mpi);
      std::vector<Point<dim> > position_point (number_quadrature_points_mpi);

      // Store density and  from MaterialModel and allocate to density_all array: 
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      MaterialModel::MaterialModelInputs<dim> in(quadrature_formula.size(),this->n_compositional_fields());      
      MaterialModel::MaterialModelOutputs<dim> out(quadrature_formula.size(),this->n_compositional_fields());
      local_cell_number = 0;
      for (; cell!=endc; ++cell)
        {
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              const std::vector<Point<dim> > &position_point_cell = fe_values.get_quadrature_points();
              in.reinit(fe_values, cell, this->introspection(), this->get_solution(), false);
              this->get_material_model().evaluate(in, out);
              for (unsigned int q = 0; q < number_quadrature_points_cell; ++q)
                {
                 density_JxW[local_cell_number * number_quadrature_points_cell + q] = out.densities[q] * fe_values.JxW(q);
                 position_point[local_cell_number * number_quadrature_points_cell + q] = position_point_cell[q];
                 
                }
              ++local_cell_number;
            }
         }

      // loop on r - satellite position [r, ,  ]
      for (unsigned int h=0; h < number_points_radius; ++h)
        {
          std_cxx11::array<double,dim> satellite_coordinate;
          satellite_coordinate[0] = minimum_radius + ((maximum_radius - minimum_radius) / number_points_radius) * h;
          
          // loop on phi - satellite position [ , phi ,]
          for (unsigned int i=0; i < number_points_longitude; ++i)
            {
              satellite_coordinate[1] = minimum_longitude + ((maximum_longitude - minimum_longitude) / number_points_longitude) * i;

              // loop on theta - satllite position [ , , theta]
              for (unsigned int j=0; j < number_points_latitude; ++j)
                { 
                  satellite_coordinate[2] = minimum_latitude + ((maximum_latitude - minimum_latitude) / number_points_latitude) * j;
                  Tensor<1,dim> local_g;
                  double local_U = 0;

                  // The spherical coordinates are shifted into cartesian to allow simplification in the mathematical equation.
                  const Point<dim> position_satellite = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(satellite_coordinate);

                  // For each point (i.e. satellite), the fourth integral goes over cells and quadrature 
                  // points to get the unique distance between those, indispensable to calculate 
                  // gravity vector components x,y,z, and potential. 
                  cell = this->get_dof_handler().begin_active();
                  local_cell_number = 0;
                  for (; cell!=endc; ++cell)
                    {
                      if (cell->is_locally_owned())
                        {
                          //fe_values.reinit (cell);  // required if position_point out? because of fe_values.JxW(q)?
                          for (unsigned int q = 0; q < number_quadrature_points_cell; ++q)
                            {
                              double dist = (position_satellite - position_point[local_cell_number * number_quadrature_points_cell + q]).norm();
                              double KK = G * density_JxW[local_cell_number * number_quadrature_points_cell + q] / pow(dist,3);
                              local_g += KK * (position_satellite - position_point[local_cell_number * number_quadrature_points_cell + q]);
                              local_U -= G * density_JxW[local_cell_number * number_quadrature_points_cell + q] / dist;
                            }
                          ++local_cell_number;
                        }
                    }

                  // assemble gravity component results from the mpi
                  const Tensor<1,dim> g
                    = Utilities::MPI::sum (local_g, this->get_mpi_communicator());
                  const double U
                    = Utilities::MPI::sum (local_U, this->get_mpi_communicator());

                  // analytical solution to calculate the theoritical gravity from a uniform desnity model.
                  // can only be used if concentric density profile
                  const double reference_density = (this->get_adiabatic_conditions().density(in.position[0])); //this->get_geometry_model.representative_point(0)));
                  double g_theory = 0;      
                  double g_diff = 0;      
                  if (satellite_coordinate[0] <= model_inner_radius)
                    {
                      g_theory = 0;
                      g_diff = -g.norm();
                    }
                  else if ((satellite_coordinate[0] > model_inner_radius) && (satellite_coordinate[0] < model_outer_radius))
                    {
                      g_theory = (G * numbers::PI * 4/3 * reference_density * (satellite_coordinate[0] - 
                                 (pow(model_inner_radius,3) / pow(satellite_coordinate[0],2))));
                      g_diff = (g_theory - g.norm());
                    }
                  else
                    {
                      g_theory = (G * numbers::PI * 4/3 * reference_density * (pow(model_outer_radius,3) - 
                                  pow(model_inner_radius,3)) / pow(satellite_coordinate[0],2));
                      g_diff = (g_theory - g.norm());
                    }

                  // write output
                  f << satellite_coordinate[0] << ' '
                    << satellite_coordinate[1] << ' '
                    << satellite_coordinate[2] << ' '
                    << position_satellite[0] << ' '
                    << position_satellite[1] << ' '
                    << position_satellite[2] << ' '
                    << g << ' '
                    << g.norm() << ' '
                    << g_theory << ' '
                    << g_diff << ' '
                    << U
                    << '\n';

                }
            }
        }
      return std::pair<std::string, std::string> ("gravity computation file:",
                                                  filename);
    }

    template <int dim>
    void
    GravityPointValues<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Postprocess");
      {
        prm.enter_subsection ("Gravity calculation");
        {
          prm.declare_entry ("Number points radius", "1",
                             Patterns::Double (0.0));
          prm.declare_entry ("Number points longitude", "1",
                             Patterns::Double (0.0));
          prm.declare_entry ("Number points latitude", "1",
                           Patterns::Double (0.0));
          prm.declare_entry ("Minimum radius", "0",
                             Patterns::Double (0.0));
          prm.declare_entry ("Maximum radius", "0",
                             Patterns::Double (0.0));
          prm.declare_entry ("Minimum longitude", "0",
                             Patterns::Double (0.0,360.0));
          prm.declare_entry ("Minimum latitude", "0",
                             Patterns::Double (0.0,180.0));
          prm.declare_entry ("Maximum longitude", "360",
                             Patterns::Double (0.0,360.0));
          prm.declare_entry ("Maximum latitude", "180",
                             Patterns::Double (0.0,180.0));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    GravityPointValues<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Postprocess");
      {
        prm.enter_subsection ("Gravity calculation");
        {
          number_points_radius    = prm.get_double ("Number points radius");
          number_points_longitude = prm.get_double ("Number points longitude");
	  number_points_latitude  = prm.get_double ("Number points latitude");
	  minimum_radius    = prm.get_double ("Minimum radius");
	  maximum_radius    = prm.get_double ("Maximum radius");
          minimum_longitude = prm.get_double ("Minimum longitude");
          maximum_longitude = prm.get_double ("Maximum longitude");
          minimum_latitude  = prm.get_double ("Minimum latitude");
          maximum_latitude  = prm.get_double ("Maximum latitude");
          AssertThrow (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != 0 ||
                       dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()) != 0 ||
                       dynamic_cast<const GeometryModel::Chunk<dim>*> (&this->get_geometry_model()) != 0,
                       ExcMessage ("This postprocessor can only be used if the geometry "
                                   "is a sphere or spherical shell."));
          AssertThrow (! this->get_material_model().is_compressible(), 
                       ExcMessage("This postprocessor was only tested for incompressible models so far."));
          AssertThrow (dim==3, 
                       ExcMessage("This postprocessor was only tested for 3D models."));
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
    ASPECT_REGISTER_POSTPROCESSOR(GravityPointValues,
                                  "gravity calculation",
                                  "A postprocessor that computes gravity and gravity potential "
                                  "for a set of points (e.g. satellites) in or above the model "
                                  "surface for a user-defined range of latitudes, longitudes "
                                  "and radius. ")
  }
}
