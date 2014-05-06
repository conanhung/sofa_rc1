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
#ifndef SOFA_VIEWER_H
#define SOFA_VIEWER_H

#include <sofa/gui/qt/SofaGUIQt.h>
#include <qstring.h>
#include <qwidget.h>
#include <sofa/helper/Factory.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/CollisionModel.h>
#ifdef SOFA_QT4 
#include <QEvent>
#include <QMouseEvent> 
#include <QKeyEvent>
#include <QTabWidget>
#include <QTimer>
#else
#include <qevent.h>
#include <qtabwidget.h>
#include <qtimer.h>
#endif
#include <sofa/helper/system/FileRepository.h>
#include <sofa/helper/system/SetDirectory.h> 
#include <sofa/gui/PickHandler.h>
#include <sofa/gui/qt/PickHandlerCallBacks.h>
#include <qcursor.h>

#include <sofa/helper/gl/Capture.h>
#include <sofa/helper/gl/Texture.h>
#include "../SofaVideoRecorderManager.h"
#ifdef SOFA_HAVE_FFMPEG
#include <sofa/helper/gl/VideoRecorder.h>
#endif //SOFA_HAVE_FFMPEG

#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/objectmodel/KeyreleasedEvent.h>
#include <sofa/core/objectmodel/MouseEvent.h>
#include <sofa/core/collision/Pipeline.h>

#include <sofa/component/configurationsetting/ViewerSetting.h>

//instruments handling
#include <sofa/component/controller/Controller.h>
#include <sofa/defaulttype/LaparoscopicRigidTypes.h>
//#include <sofa/simulation/common/GrabVisitor.h>
#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/simulation/common/UpdateMappingVisitor.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/component/visualmodel/InteractiveCamera.h>

#ifdef SOFA_QT4
#include <QEvent>
#include <QMouseEvent>
#include <QKeyEvent>
#else
#include <qevent.h>
#endif

#include <sofa/gui/qt/viewer/VisualModelPolicy.h>
#include <sofa/helper/io/Image.h>

namespace sofa
{

namespace gui
{

namespace qt
{

namespace viewer
{

using namespace sofa::defaulttype;
enum
{
  BTLEFT_MODE = 101, BTRIGHT_MODE = 102, BTMIDDLE_MODE = 103,
};




class SOFA_SOFAGUIQT_API SofaViewer
{

public:
  SofaViewer();
  virtual ~SofaViewer();

  virtual void drawColourPicking (ColourPickingVisitor::ColourCode /*code*/) {}

  /// Optional QTabWidget GUI for a concreate viewer.
  virtual void removeViewerTab(QTabWidget *){}
  /// Optional QTabWidget GUI for a concreate viewer.
  virtual void configureViewerTab(QTabWidget *){}

  /// Overload this method in instanciated viewer to notify if the viewer is embeded (like a Qt widget) or not (standalone)
  virtual bool isEmbedded(){return true;}
  virtual QWidget* getQWidget()=0;
  virtual QString helpString()=0;

  /// Overload this method in instanciated viewer to notify if the viewer is managed in another render thread or not
  virtual bool isThreaded(void);
  virtual bool startViewerThread(void);
  virtual bool stopViewerThread(void);

  virtual sofa::simulation::Node* getScene();
  virtual const std::string& getSceneFileName();
  virtual void setSceneFileName(const std::string &f);
  virtual void setScene(sofa::simulation::Node::SPtr scene, const char* filename =
                        NULL, bool /*keepParams*/= false);
  virtual void setCameraMode(core::visual::VisualParams::CameraType);

  virtual bool ready();
  virtual void wait(){}

  /// Call when set a new scene
  virtual bool loadSceneView(void);
  /// Call before delete a scene
  virtual bool unloadSceneView(void);

  /// Recompute viewer's home position so it encompass the whole scene and apply it
  virtual void viewAll(void) = 0;

  //Allow to configure your viewer using the Sofa Component, ViewerSetting
  virtual void configure(sofa::component::configurationsetting::ViewerSetting* viewerConf);

  //Fonctions needed to take a screenshot
  virtual const std::string screenshotName();
  virtual void setPrefix(const std::string filename);
  virtual void screenshot(const std::string filename, int compression_level =
                          -1);

  virtual void getView(Vec3d& pos, Quat& ori) const;
  virtual void setView(const Vec3d& pos, const Quat &ori);
  virtual void moveView(const Vec3d& pos, const Quat &ori);
  virtual void newView();   
  virtual void resetView();

  virtual void setBackgroundColour(float r, float g, float b);
  virtual void setBackgroundImage(std::string imageFileName = std::string("textures/SOFA_logo.bmp"));
  std::string getBackgroundImage();

  PickHandler* getPickHandler();

  //*************************************************************
  // QT
  //*************************************************************
  //SLOTS
  virtual void saveView()=0;
  virtual void setSizeW(int)=0;
  virtual void setSizeH(int)=0;
  virtual void captureEvent();
  virtual void fitObjectBBox(sofa::core::objectmodel::BaseObject* );
  virtual void fitNodeBBox(sofa::core::objectmodel::BaseNode*);
	

  //SIGNALS
  virtual void redrawn()=0;
  virtual void resizeW(int)=0;
  virtual void resizeH(int)=0;

protected:

  // ---------------------- Here are the Keyboard controls   ----------------------
  void keyPressEvent(QKeyEvent * e);
  void keyReleaseEvent(QKeyEvent * e);
  bool isControlPressed() const;
  // ---------------------- Here are the Mouse controls   ----------------------
  void wheelEvent(QWheelEvent *e);
  void mouseMoveEvent ( QMouseEvent *e );
  void mousePressEvent ( QMouseEvent * e);
  void mouseReleaseEvent ( QMouseEvent * e);
  void mouseEvent(QMouseEvent *e);
	
  // ---------------------- RayCasting PickHandler  ----------------------
  virtual void moveRayPickInteractor(int, int)
  {
  }
  ;

  sofa::simulation::Node::SPtr groot;
  sofa::component::visualmodel::BaseCamera::SPtr currentCamera;  
  std::string sceneFileName;
  sofa::helper::gl::Capture capture;
  
#ifdef SOFA_HAVE_FFMPEG
  sofa::helper::gl::VideoRecorder videoRecorder;
#endif //SOFA_HAVE_FFMPEG
  QTimer captureTimer;
	
  bool m_isControlPressed;
  bool _video;
  bool _axis;
  int _background;
  bool initTexturesDone;

  Vector3 backgroundColour;
  sofa::helper::gl::Texture* texLogo;
  std::string backgroundImageFile;

  Vector3 ambientColour;

  PickHandler pick;
  ColourPickingRenderCallBack colourPickingRenderCallBack;

  //instruments handling
  int _navigationMode;
  bool _mouseInteractorMoving;
  int _mouseInteractorSavedPosX;
  int _mouseInteractorSavedPosY;

  //Stereo parameters
  bool _stereoEnabled;
  bool _binocularModeEnabled;
  double _stereoShift;

};

template < typename VisualModelPolicyType >
class CustomPolicySofaViewer : public VisualModelPolicyType, public sofa::gui::qt::viewer::SofaViewer
{
public:
  using VisualModelPolicyType::load;
  using VisualModelPolicyType::unload;
  CustomPolicySofaViewer() { load(); }
  virtual ~CustomPolicySofaViewer() { unload(); }
protected:
};

typedef CustomPolicySofaViewer< OglModelPolicy > OglModelSofaViewer;


}
}
}
}

#endif
