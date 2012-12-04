//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id:$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef AXISYM_MAGNETOSTATICS_H
#define AXISYM_MAGNETOSTATICS_H

//libMesh
#include "libmesh.h"
#include "boundary_info.h"
#include "fe_base.h"
#include "fe_interface.h"
#include "mesh.h"
#include "quadrature.h"
#include "parameters.h"
#include "string_to_enum.h"
#include "fem_system.h"
#include "fem_context.h"

//GRINS
#include "grins_config.h"
#include "physics.h"
#include "axisym_magnetostatics_bc_handling.h"

namespace GRINS
{
  //! Physics class for Axisymmetric Magnetostatics
  /*
    This physics class implements magentostatics using the potential form of the
    equations.
   */
  class AxisymmetricMagnetostatics : public Physics
  {
  public:
    AxisymmetricMagnetostatics( const std::string& physics_name, const GetPot& input );
    ~AxisymmetricMagnetostatics();

    //! Read options from GetPot input file.
    virtual void read_input_options( const GetPot& input );

    //! Initialization  AxisymmetricHeatTransfer variables
    /*!
      Add velocity and pressure variables to system.
     */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( libMesh::DiffContext &context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual bool element_time_derivative( bool request_jacobian,
					  libMesh::DiffContext& context,
					  libMesh::FEMSystem* system );

    virtual bool side_time_derivative( bool request_jacobian,
				       libMesh::DiffContext& context,
				       libMesh::FEMSystem* system );

    // Constraint part(s)
    virtual bool element_constraint( bool request_jacobian,
				     libMesh::DiffContext& context,
				     libMesh::FEMSystem* system );

    virtual bool side_constraint( bool request_jacobian,
				  libMesh::DiffContext& context,
				  libMesh::FEMSystem* system );

    // Mass matrix part(s)
    virtual bool mass_residual( bool request_jacobian,
				libMesh::DiffContext& context,
				libMesh::FEMSystem* system );

  protected:

    //! Physical dimension of problem
    /*! \todo Make this static member of base class? */
    unsigned int _dim;

    // Indices for each variable;
    //! Index for magnetic potential
    VariableIndex _A_var;

    //! Index for electric potential
    VariableIndex _V_var;

    // Names of each variable in the system
    //! Name for temperature variable
    std::string _A_var_name;

    //! Name of r-velocity
    std::string _V_var_name;

    //! Element type, read from input
    libMeshEnums::FEFamily _A_FE_family, _V_FE_family;

    //! Temperature element order, read from input
    libMeshEnums::Order _A_order, _V_order;

    Real _sigma, _mu;

  private:
    AxisymmetricMagnetostatics();

  }; // class AxisymmetricMagnetostatics

} // namespace GRINS
#endif //AXISYM_MAGNETOSTATICS_H