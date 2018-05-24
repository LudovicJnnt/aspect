/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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

      std_cxx11::array<double,dim> ecoord;

      int hh=1; int ii=1; int jj=1;
      double height=0;
      double h_spacing=0;
      double longitude=0;
      double latitude=0; 

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
      f2 << minimum_longitude << ' ' << maximum_longitude << '\n'
         << minimum_latitude << ' ' << maximum_latitude << '\n'
         << mapping_spherical_spacing << '\n'
         << (maximum_longitude-minimum_longitude)/mapping_spherical_spacing << '\n';

      // Storing density values from MaterialModel in a vector per MPI avoids to 
      // use MaterialModel to get the density at quadrature points within the loops. 
      // Total number of element used per MPI:
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
      f2 << profile_satellite_height << ' ' << mapping_satellite_height << '\n';

      // Decide if gravity depth-profile or mapping is activated
      if (profile_satellite_height > 0) // gravity depth-profile activated
        { 
          longitude = profile_longitude;
          latitude  = profile_latitude;
          h_spacing = profile_height_spacing;
          hh = profile_satellite_height/profile_height_spacing;
        }
      else if (mapping_satellite_height > 0) // gravity mapping activated - first integral on [r, , ] called only once
        {
          height    = mapping_satellite_height;
          h_spacing = profile_satellite_height;
          longitude = minimum_longitude;
          latitude  = minimum_latitude; 
          ii = (maximum_longitude-minimum_longitude)/mapping_spherical_spacing;
          jj = (maximum_latitude-minimum_latitude)/mapping_spherical_spacing;
        }
       else
        {
          f2 << '\n' << '\n';
        }       
 
      // loop on r - satellite position [r, ,  ]
      for (int h=0; h<hh; h++)
        {
          ecoord[0]=height+h*h_spacing;
          
          // loop on phi - satellite position [ , phi ,]
          for (int i=0; i<ii; i++)
            {
              ecoord[1]=longitude+i*mapping_spherical_spacing;

              // loop on theta - satllite position [ , , theta]
              for (int j=0; j<jj; j++)
                { 
                  ecoord[2]=latitude+j*mapping_spherical_spacing;
                  f2 << longitude  << ' ' << latitude  << ' ' << height << ' ' << h_spacing << '\n';
                  f2 << ecoord[0] << ' ' << ecoord[1] << ' ' << ecoord [2] << ' ' << j << ' ' << i << ' ' << j << ' ' << hh << ' ' << ii << ' ' << jj << '\n';
                  
                  // The spherical coordinates are shifted into cartesian to allow simplification in the mathematical equation.
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
              
                  // For each point (i.e. satellite), the forth integral goes over cells and quadrature 
                  // points to get the unique distance between those, indispensable to calculate 
                  // gravity vector components x,y,z, and potential. 
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
                      gtheory = G * numbers::PI * 4/3 * reference_density * (ecoord[0] - (pow(Ri,3) / pow(ecoord[0],2)));
                      g_diff = (gtheory - gnorm);
                    }
                  else
                    {
                      gtheory = G * numbers::PI * 4/3 * reference_density * (pow(Ro,3)-pow(Ri,3)) / pow(ecoord[0],2);
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
        }
      return std::pair<std::string, std::string> ("gravity computation file:",
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
          prm.declare_entry ("Mapping spherical spacing", "1",
                             Patterns::Double ());
          prm.declare_entry ("Profile height spacing", "1",
                             Patterns::Double ());
          prm.declare_entry ("Mapping satellite height", "0",
                             Patterns::Double ());
          prm.declare_entry ("Profile satellite height", "0",
                             Patterns::Double ());
          prm.declare_entry ("Profile longitude", "0",
                             Patterns::Double ());
          prm.declare_entry ("Profile latitude", "0",
                             Patterns::Double ());
          prm.declare_entry ("Minimum longitude", "0",
                             Patterns::Double ());
          prm.declare_entry ("Minimum latitude", "0",
                             Patterns::Double ());
          prm.declare_entry ("Maximum longitude", "360",
                             Patterns::Double ());
          prm.declare_entry ("Maximum latitude", "180",
                             Patterns::Double ());
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
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
          mapping_satellite_height  = prm.get_double ("Mapping satellite height");
          profile_satellite_height  = prm.get_double ("Profile satellite height");
	  mapping_spherical_spacing = prm.get_double ("Mapping spherical spacing");
	  profile_height_spacing    = prm.get_double ("Profile height spacing");
          profile_longitude = prm.get_double ("Profile longitude");
          profile_latitude  = prm.get_double ("Profile latitude");
          minimum_longitude = prm.get_double ("Minimum longitude");
          maximum_longitude = prm.get_double ("Maximum longitude");
          minimum_latitude  = prm.get_double ("Minimum latitude");
          maximum_latitude  = prm.get_double ("Maximum latitude");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
      prm.enter_subsection ("Material model");
      {
        prm.enter_subsection ("Simple model");
        {
	  reference_density = prm.get_double ("Reference density");
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
                                  "A postprocessor that computes gravity and potential for "
                                  "a set of points (e.g. satellites) above the model surface "
                                  "for a user-defined range of latitudes, longitudes and height. "
                                  "A gravity depth-profile for a user-defined longitude and "
                                  "latitude is also possible.")
  }
}
