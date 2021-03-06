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
#include <sofa/component/visualmodel/VisualStyle.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/objectmodel/Context.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/simulation/common/UpdateContextVisitor.h>

namespace sofa
{
namespace component
{
namespace visualmodel
{

using namespace sofa::core::visual;
using namespace sofa::core::objectmodel;
using namespace sofa::simulation;

int VisualStyleClass = core::RegisterObject("Edit the visual style.").add<VisualStyle>();

VisualStyle::VisualStyle()
  :displayFlags(initData(&displayFlags,"displayFlags","Display Flags"))
{
  displayFlags.setWidget("widget_displayFlags");
  displayFlags.setGroup("Display Flags");
}

void VisualStyle::fwdDraw(VisualParams* vparams)
{
  backupFlags = vparams->displayFlags();
  vparams->displayFlags() = sofa::core::visual::merge_displayFlags(backupFlags, displayFlags.getValue(vparams));
}

void VisualStyle::bwdDraw(VisualParams* vparams)
{
  vparams->displayFlags() = backupFlags;

}


}
}
}

