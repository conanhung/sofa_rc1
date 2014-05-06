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
#include <sofa/core/ExecParams.h>
#include <sofa/helper/system/thread/thread_specific_ptr.h>
#include <iostream>

namespace sofa
{

namespace core
{

sofa::helper::system::atomic<int> ExecParams::g_nbThreads = 0;

ExecParams::ExecParamsThreadStorage::ExecParamsThreadStorage(int tid)
: execMode(EXEC_DEFAULT)
, threadID(tid)
, aspectID(0)
{
}

/// Get the default ExecParams, to be used to provide a default values for method parameters
ExecParams* ExecParams::defaultInstance()
{
    SOFA_THREAD_SPECIFIC_PTR(ExecParams, threadParams);
    ExecParams* ptr = threadParams;
    if (!ptr)
    {
        ptr = new ExecParams(new ExecParamsThreadStorage(g_nbThreads.exchange_and_add(1)));
        threadParams = ptr;
        if (ptr->threadID())
            std::cout << "[THREAD " << ptr->threadID() << "]: local ExecParams storage created." << std::endl;
    }
    return ptr;
}

ExecParams::ExecParamsThreadStorage* ExecParams::threadStorage()
{
    return defaultInstance()->storage;
}

/// Make sure this instance is up-to-date relative to the current thread
void ExecParams::update()
{
    storage = threadStorage();
}

} // namespace core

} // namespace sofa
