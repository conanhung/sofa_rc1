/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_GUI_QT_VIEWERFACTORY_H
#define SOFA_GUI_QT_VIEWERFACTORY_H

#include <sofa/gui/qt/SofaGUIQt.h>
#include <sofa/helper/Factory.h>
#include <sofa/helper/Factory.inl>
#include <sofa/gui/qt/viewer/SofaViewer.h>

namespace sofa{
  namespace gui{
    namespace qt{
      namespace viewer{

  struct SofaViewerArgument
  {
    QWidget* parent;
    std::string name;
  };
      }
    }
  }
}

namespace sofa{
    namespace helper{

        template < >
        class BaseCreator< sofa::gui::qt::viewer::SofaViewer, sofa::gui::qt::viewer::SofaViewerArgument>
        {
        public:
            virtual ~BaseCreator() { }
            virtual sofa::gui::qt::viewer::SofaViewer *createInstance(sofa::gui::qt::viewer::SofaViewerArgument arg) = 0;
            virtual const std::type_info& type() = 0;
            virtual const char* viewerName() = 0;
            virtual const char* acceleratedName() = 0;
        };


        class SOFA_SOFAGUIQT_API SofaViewerFactory : public sofa::helper::Factory< std::string, sofa::gui::qt::viewer::SofaViewer, sofa::gui::qt::viewer::SofaViewerArgument >
        {
        public:
            typedef sofa::helper::Factory< std::string, sofa::gui::qt::viewer::SofaViewer, sofa::gui::qt::viewer::SofaViewerArgument > Inherited;
            typedef Inherited::Key Key;
            typedef Inherited::Argument Argument;
            typedef Inherited::Object Object;
            typedef Inherited::Creator Creator;


            static SofaViewerFactory*  getInstance()
            {
                static SofaViewerFactory instance;
                return &instance;
            }

            static Object* CreateObject(Key key, Argument arg)
            {
                    return getInstance()->createObject(key, arg);
            }

            static Object* CreateAnyObject(Argument arg)
            {
                    return getInstance()->createAnyObject(arg);
            }

            static bool HasKey(Key key)
            {
                    return getInstance()->hasKey(key);
            }

            const char* getViewerName(Key key)
            {

                Creator* creator;
                std::multimap<Key, Creator*>::iterator it = this->registry.lower_bound(key);
                std::multimap<Key, Creator*>::iterator end = this->registry.upper_bound(key);
                while (it != end)
                {
                        creator = (*it).second;
                        const char* viewerName = creator->viewerName();
                        if(viewerName != NULL )
                        {
                            return viewerName;
                        }
                        ++it;
                }
        //	std::cerr<<"Object type "<<key<<" creation failed."<<std::endl;
                return NULL;
            }

            const char* getAcceleratedViewerName(Key key)
            {

                Creator* creator;
                std::multimap<Key, Creator*>::iterator it = this->registry.lower_bound(key);
                std::multimap<Key, Creator*>::iterator end = this->registry.upper_bound(key);
                while (it != end)
                {
                        creator = (*it).second;
                        const char* acceleratedName = creator->acceleratedName();
                        if(acceleratedName != NULL )
                        {
                            return acceleratedName;
                        }
                        ++it;
                }
        //	std::cerr<<"Object type "<<key<<" creation failed."<<std::endl;
                return NULL;

            }

            static const char* ViewerName( Key key)
            {
               return getInstance()->getViewerName(key);
            }
            static const char* AcceleratedName( Key key )
            {
               return getInstance()->getAcceleratedViewerName(key);
            }
        };

        template < class RealObject >
        class SofaViewerCreator : public Creator< SofaViewerFactory, RealObject >
        {
        public:
            typedef Creator< SofaViewerFactory, RealObject > Inherited;
            typedef SofaViewerFactory::Object Object;
            typedef SofaViewerFactory::Argument Argument;
            typedef SofaViewerFactory::Key Key;
            SofaViewerCreator(Key key, bool multi=false):Inherited(key,multi)
            {
            }
             const char* viewerName()
             {
                return RealObject::viewerName();
             }
             const char* acceleratedName()
             {
                return RealObject::acceleratedName();
             }

        };


    }
}

#endif //SOFA_GUI_QT_VIEWERFACTORY_H

