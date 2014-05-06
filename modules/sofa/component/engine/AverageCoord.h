/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef SOFA_COMPONENT_ENGINE_AverageCoord_H
#define SOFA_COMPONENT_ENGINE_AverageCoord_H

#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/VecId.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace engine
{

using namespace core::behavior;
using namespace core::topology;
using namespace core::objectmodel;

/** 
 * This class computes the average of a set of Coordinates
 */
template <class DataTypes>
class AverageCoord : public core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AverageCoord,DataTypes),core::DataEngine);
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef unsigned int Index;
    typedef sofa::helper::vector<Index> VecIndex;

protected:

    AverageCoord();

    virtual ~AverageCoord(){}
public:
    void init();

    void reinit();

    void update();

    Data<VecIndex> f_indices;    ///< indices of the coordinates to average
    Data<unsigned> f_vecId;  ///< index of the vector (default value corresponds to core::VecCoordId::position() )
    Data<Coord> f_average;       ///< result

    virtual std::string getTemplateName() const
    {
      return templateName(this);
    }


    static std::string templateName(const AverageCoord<DataTypes>* = NULL)
    {
      return DataTypes::Name();
    }



protected:
    MechanicalState<DataTypes> *mstate;
};

#if defined(WIN32) && !defined(SOFA_COMPONENT_ENGINE_AverageCoord_CPP)
#pragma warning(disable : 4231)
#ifndef SOFA_FLOAT
template class SOFA_ENGINE_API AverageCoord<defaulttype::Vec2dTypes>;
template class SOFA_ENGINE_API AverageCoord<defaulttype::Vec3dTypes>;
template class SOFA_ENGINE_API AverageCoord<defaulttype::Rigid2dTypes>;
template class SOFA_ENGINE_API AverageCoord<defaulttype::Rigid3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_ENGINE_API AverageCoord<defaulttype::Vec2fTypes>;
template class SOFA_ENGINE_API AverageCoord<defaulttype::Vec3fTypes>;
template class SOFA_ENGINE_API AverageCoord<defaulttype::Rigid2fTypes>;
template class SOFA_ENGINE_API AverageCoord<defaulttype::Rigid3fTypes>;
#endif //SOFA_DOUBLE
#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif