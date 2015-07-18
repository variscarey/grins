//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-


#include "grins_config.h"

// This class
#include "grins/laminar_1D_flame_base.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/cached_quantities_enum.h"
#include "grins/cantera_mixture.h"
#include "grins/grins_enums.h"
#include "grins/antioch_mixture.h"

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  LaminarFlameBase::LaminarFlameBase(const std::string& physics_name, 
									    const GetPot& input)
    : Physics(physics_name, input),
      _fixed_density( input("Physics/"+reacting_low_mach_navier_stokes+"/fixed_density", false ) ),
      _fixed_rho_value( input("Physics/"+reacting_low_mach_navier_stokes+"/fixed_rho_value", 0.0 ) )
  {
    this->read_input_options(input);
    
    return;
  }

  Laminar1DFlameBase::~LaminarFlameBase()
  {
    return;
  }
  
  void Laminar1DFlameBase::read_input_options( const GetPot& input )
  {
    // Read FE family info
    this->_species_FE_family = libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>( input("Physics/"+reacting_low_mach_navier_stokes+"/species_FE_family", "LAGRANGE") );

    this->_T_FE_family = libMesh::Utility::string_to_enum<GRINSEnums::FEFamily>( input("Physics/"+reacting_low_mach_navier_stokes+"/T_FE_family", "LAGRANGE") );

    // Read FE family info
    this->_species_order = libMesh::Utility::string_to_enum<GRINSEnums::Order>( input("Physics/"+reacting_low_mach_navier_stokes+"/species_order", "SECOND") );

    this->_T_order = libMesh::Utility::string_to_enum<GRINSEnums::Order>( input("Physics/"+reacting_low_mach_navier_stokes+"/T_order", "SECOND") );

    // Read variable naming info
    this->_n_species = input.vector_variable_size("Physics/Chemistry/species");

    _species_var_names.reserve(this->_n_species);
    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	/*! \todo Make this prefix string an input option */
	std::string var_name = "w_"+std::string(input( "Physics/Chemistry/species", "DIE!", i ));
	_species_var_names.push_back( var_name );
      }

    this->_T_var_name = input("Physics/VariableNames/temperature", GRINS::T_var_name_default );

    // Read thermodynamic state info
    _p0 = input("Physics/"+reacting_low_mach_navier_stokes+"/p0", 0.0 ); /* thermodynamic pressure */

    /*_enable_thermo_press_calc = input("Physics/"+reacting_low_mach_navier_stokes+"/enable_thermo_press_calc", false );

    if( _enable_thermo_press_calc )
      {
	_p0_var_name = input("Physics/VariableNames/thermo_presure", "p0" );
	} */

    return;
  }

  void LaminarFlameBase::init_variables( libMesh::FEMSystem* system )
  {
    // Get libMesh to assign an index for each variable
    this->_dim = system->get_mesh().mesh_dimension();
    
    _species_vars.reserve(this->_n_species);
    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	_species_vars.push_back( system->add_variable( _species_var_names[i], 
						       this->_species_order, _species_FE_family) );
      }

    _T_var = system->add_variable( _T_var_name, this->_T_order, _T_FE_family);

    /* If we need to compute the thermodynamic pressure, we force this to be a first
       order scalar variable. */
    if( _propagating_flame )
      _Mdot_var = system->add_variable( _p0_var_name, libMesh::FIRST, libMesh::SCALAR);

    return;
  }

  void LaminarFlameBase::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    // const unsigned int dim = system->get_mesh().mesh_dimension();

    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	system->time_evolving( _species_vars[i] );
      }

    system->time_evolving(_T_var);
    system->time_evolving(_p_var);

    if( _propagating_flame )
      system->time_evolving(_Mdot_var);

    return;
  }

  void LaminarFlameBase::init_context( AssemblyContext& context )
  {
    // We should prerequest all the data
    // we will need to build the linear system
    // or evaluate a quantity of interest.
    context.get_element_fe(_species_vars[0])->get_JxW();
    context.get_element_fe(_species_vars[0])->get_phi();
    context.get_element_fe(_species_vars[0])->get_dphi();
    context.get_element_fe(_species_vars[0])->get_xyz();

    context.get_element_fe(_T_var)->get_JxW();
    context.get_element_fe(_T_var)->get_phi();
    context.get_element_fe(_T_var)->get_dphi();
    context.get_element_fe(_T_var)->get_xyz();

    context.get_element_fe(_Mdot_var)->get_phi();
    context.get_element_fe(_Mdot_var)->get_xyz();

    return;
  }

} // end namespace GRINS









