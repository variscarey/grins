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
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef STAB_HELPER_H
#define STAB_HELPER_H

// libMesh
#include "fe_base.h"
#include "vector_value.h"
#include "tensor_value.h"
#include "fem_context.h"

namespace GRINS
{
  class StabilizationHelper
  {
  public:
    StabilizationHelper();
    ~StabilizationHelper();

    /*! \todo Should we inline this? */
    libMesh::RealGradient compute_g( libMesh::FEBase* fe,
				     libMesh::FEMContext& c,
				     unsigned int qp ) const;
    
    /*! \todo Should we inline this? */
    libMesh::RealTensor compute_G( libMesh::FEBase* fe,
				   libMesh::FEMContext& c,
				   unsigned int qp ) const;

  }; // class StabilizationHelper

} // namespace GRINS
#endif // STAB_HELPER_H