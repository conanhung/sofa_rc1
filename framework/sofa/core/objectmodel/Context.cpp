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
#include <sofa/core/objectmodel/Context.h>
// #include <sofa/simulation/common/Visitor.h>


namespace sofa
{

namespace core
{

namespace objectmodel
{

Context::Context()
: is_activated(initData(&is_activated, true, "activated", "To Activate a node"))
  , worldGravity_(initData(&worldGravity_, Vec3((SReal)0,(SReal)-9.81,(SReal)0),"gravity","Gravity in the world coordinate system"))
  , dt_(initData(&dt_,0.01,"dt","Time step"))
  , time_(initData(&time_,0.,"time","Current time"))
  , animate_(initData(&animate_,false,"animate","Animate the Simulation(applied at initialization only)"))
#ifdef SOFA_SMP
  ,  processor(initData(&processor,(int )-1,"processor","assigned processor"))
  ,  gpuPrioritary(initData(&gpuPrioritary,false,"gpuPrioritary","node should be executed on GPU")),
//  is_partition_(initData(&is_partition_,false,"partition","is a parallel partition"))
  partition_(0)
#endif
{

#ifdef SOFA_SMP
  is_partition_.setValue(false);
#endif
}

/// The Context is active
bool Context::isActive() const {return is_activated.getValue();}

/// State of the context
void Context::setActive(bool val)
{
	is_activated.setValue(val);
}



/// Simulation timestep
double Context::getDt() const
{
  return dt_.getValue();
}

/// Simulation time
double Context::getTime() const
{
  return time_.getValue();
}

/// Gravity vector in world coordinates
const Context::Vec3& Context::getGravity() const
{
  return worldGravity_.getValue();
}

/// Animation flag
bool Context::getAnimate() const
{
  return animate_.getValue();
}



//===============================================================================

/// Simulation timestep
void Context::setDt(double val)
{
  dt_.setValue(val);
}

/// Simulation time
void Context::setTime(double val)
{
  time_.setValue(val);
}

/// Gravity vector
// void Context::setGravity(const Vec3& g)
// {
// 	gravity_ = g;
// }

/// Gravity vector
void Context::setGravity(const Vec3& g)
{
  worldGravity_ .setValue(g);
}

/// Animation flag
void Context::setAnimate(bool val)
{
  animate_.setValue(val);
}



//======================


void Context::copyContext(const Context& c)
{
  // BUGFIX 12/01/06 (Jeremie A.): Can't use operator= on the class as it will copy other data in the BaseContext class (such as name)...
  // *this = c;

  copySimulationContext(c);

}
#ifdef SOFA_SMP
int Context::getProcessor() const
{
  return processor.getValue();
}
void Context::setProcessor(int p) 
{
  processor.setValue(p);
}
#endif


void Context::copySimulationContext(const Context& c)
{
  worldGravity_.setValue(c.getGravity());  ///< Gravity IN THE WORLD COORDINATE SYSTEM.
  setDt(c.getDt());
  setTime(c.getTime());
  setAnimate(c.getAnimate());
#ifdef SOFA_SMP
  if(c.gpuPrioritary.getValue())
    gpuPrioritary.setValue(true);
#endif



#ifdef SOFA_SMP
  if(!partition_){
    if(processor.getValue()!=-1)
      is_partition_.setValue(true);
    if(is_partition()){
    
      partition_= new Iterative::IterativePartition();
//          partition_->setCPU(processor.getValue());
    }
  }
  if(processor.getValue()==-1&&c.processor.getValue()!=-1){
  	processor.setValue(c.processor.getValue());
    is_partition_.setValue(true);
  }
  if(c.is_partition()&&!partition_){
    partition_=c.getPartition();
    is_partition_.setValue(true);
  }
  if((gpuPrioritary.getValue())&&partition_)
  {
    partition_->setGPUPrioritary();
  }

#endif

}




} // namespace objectmodel

} // namespace core

} // namespace sofa

