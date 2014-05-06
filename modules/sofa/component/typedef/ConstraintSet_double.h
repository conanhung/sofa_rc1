/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/



// File automatically generated by "generateTypedef"


#ifndef SOFA_TYPEDEF_ConstraintSet_double_H
#define SOFA_TYPEDEF_ConstraintSet_double_H

//Default files containing the declaration of the vector type
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/Mat.h>


#ifdef SOFA_GPU_CUDA
#include <sofa/gpu/cuda/CudaTypesBase.h>
#include <sofa/gpu/cuda/CudaTypes.h>
#endif
#ifdef SOFA_GPU_OPENCL
#include <sofa/gpu/opencl/OpenCLTypes.h>
#endif


#include <sofa/component/constraintset/DOFBlockerLMConstraint.h>
#include <sofa/component/constraintset/DistanceLMConstraint.h>
#include <sofa/component/constraintset/DistanceLMContactConstraint.h>
#include <sofa/component/constraintset/FixedLMConstraint.h>
#include <sofa/component/constraintset/UnilateralInteractionConstraint.h>



//---------------------------------------------------------------------------------------------
//Typedef for DOFBlockerLMConstraint
typedef sofa::component::constraintset::DOFBlockerLMConstraint<sofa::defaulttype::StdRigidTypes<3, double> > DOFBlockerLMConstraintRigid3d;
typedef sofa::component::constraintset::DOFBlockerLMConstraint<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > DOFBlockerLMConstraint3d;



//---------------------------------------------------------------------------------------------
//Typedef for DistanceLMConstraint
typedef sofa::component::constraintset::DistanceLMConstraint<sofa::defaulttype::StdRigidTypes<3, double> > DistanceLMConstraintRigid3d;
typedef sofa::component::constraintset::DistanceLMConstraint<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > DistanceLMConstraint3d;



//---------------------------------------------------------------------------------------------
//Typedef for DistanceLMContactConstraint
typedef sofa::component::constraintset::DistanceLMContactConstraint<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > DistanceLMContactConstraint3d;



//---------------------------------------------------------------------------------------------
//Typedef for FixedLMConstraint
typedef sofa::component::constraintset::FixedLMConstraint<sofa::defaulttype::StdRigidTypes<3, double> > FixedLMConstraintRigid3d;
typedef sofa::component::constraintset::FixedLMConstraint<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > FixedLMConstraint3d;



//---------------------------------------------------------------------------------------------
//Typedef for UnilateralInteractionConstraint
typedef sofa::component::constraintset::UnilateralInteractionConstraint<sofa::defaulttype::StdVectorTypes<sofa::defaulttype::Vec<3, double>, sofa::defaulttype::Vec<3, double>, double> > UnilateralInteractionConstraint3d;





#ifndef SOFA_FLOAT
typedef DOFBlockerLMConstraintRigid3d DOFBlockerLMConstraintRigid3;
typedef DOFBlockerLMConstraint3d DOFBlockerLMConstraint3;
typedef DistanceLMConstraintRigid3d DistanceLMConstraintRigid3;
typedef DistanceLMConstraint3d DistanceLMConstraint3;
typedef DistanceLMContactConstraint3d DistanceLMContactConstraint3;
typedef FixedLMConstraintRigid3d FixedLMConstraintRigid3;
typedef FixedLMConstraint3d FixedLMConstraint3;
typedef UnilateralInteractionConstraint3d UnilateralInteractionConstraint3;
#endif

#endif
