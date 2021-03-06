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
*                              SOFA :: Framework                              *
*                                                                             *
* Authors: The SOFA Team (see Authors.txt)                                    *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_CORE_BEHAVIOR_BASEVECTOROPERATION_H
#define SOFA_CORE_BEHAVIOR_BASEVECTOROPERATION_H

#include <sofa/core/MultiVecId.h>
#include <sofa/core/ExecParams.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/behavior/BaseMechanicalState.h>

namespace sofa
{

namespace core
{

namespace behavior
{

class BaseVectorOperations
{
protected:
    const core::ExecParams* params;
    core::objectmodel::BaseContext* ctx;
    double result;

public:
    BaseVectorOperations(const core::ExecParams* params, core::objectmodel::BaseContext* ctx)
    : params(params),ctx(ctx)
    {}

    virtual ~BaseVectorOperations()
    {}

    /// Allocate a temporary vector
    virtual void v_alloc(sofa::core::MultiVecCoordId& id) = 0;
    virtual void v_alloc(sofa::core::MultiVecDerivId& id) = 0;
    /// Free a previously allocated temporary vector
    virtual void v_free(sofa::core::MultiVecCoordId& id) = 0;
    virtual void v_free(sofa::core::MultiVecDerivId& id) = 0;
    
    virtual void v_clear(core::MultiVecId v) = 0; ///< v=0
    virtual void v_eq(core::MultiVecId v, core::MultiVecId a) = 0; ///< v=a
    virtual void v_peq(core::MultiVecId v, core::MultiVecId a, double f=1.0) = 0; ///< v+=f*a
#ifdef SOFA_SMP
    virtual void v_peq(core::MultiVecId v, core::MultiVecId a, Shared<double> &fSh, double f=1.0) = 0; ///< v+=f*a
    virtual void v_meq(core::MultiVecId v, core::MultiVecId a, Shared<double> &fSh) = 0; ///< v+=f*a
#endif
        virtual void v_teq(core::MultiVecId v, double f) = 0; ///< v*=f
    virtual void v_op(core::MultiVecId v, core::ConstMultiVecId a, core::ConstMultiVecId b, double f=1.0) = 0; ///< v=a+b*f
#ifdef SOFA_SMP
    virtual void v_op(core::MultiVecId v, core::MultiVecId a, core::MultiVecId b, Shared<double> &f) = 0; ///< v=a+b*f
#endif
    virtual void v_multiop(const core::behavior::BaseMechanicalState::VMultiOp& o) = 0;
    virtual void v_dot(core::ConstMultiVecId a, core::ConstMultiVecId b) = 0; ///< a dot b ( get result using finish )
#ifdef SOFA_SMP
    virtual void v_dot(Shared<double> &result,core::MultiVecId a, core::MultiVecId b) = 0; ///< a dot b
#endif
    virtual void v_threshold(core::MultiVecId a, double threshold) = 0; ///< nullify the values below the given threshold 
    
    
    virtual double finish() = 0;
    
    virtual void print( core::MultiVecId v, std::ostream& out ) = 0;
    
};

} // namespace behavior

} // namespace core

} // namespace sofa

#endif //SOFA_CORE_BEHAVIOR_BASEVECTOROPERATION_H

