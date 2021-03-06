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
#include <sofa/component/visualmodel/OglModel.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/system/gl.h>
#include <sofa/helper/system/glut.h>
#include <sofa/helper/gl/RAII.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Quat.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sstream>
#include <string.h>

//#define NO_VBO

namespace sofa
{

namespace component
{

namespace visualmodel
{

using namespace sofa::defaulttype;
using namespace sofa::core::loader;

SOFA_DECL_CLASS(OglModel)

int OglModelClass = core::RegisterObject("Generic visual model for OpenGL display")
.add< OglModel >()
;


OglModel::OglModel()
: premultipliedAlpha(initData(&premultipliedAlpha, (bool) false, "premultipliedAlpha", "is alpha premultiplied ?"))
#ifndef SOFA_HAVE_GLEW
, useVBO(initData(&useVBO, (bool) false, "useVBO", "Use VBO for rendering"))
#else
, useVBO(initData(&useVBO, (bool) true, "useVBO", "Use VBO for rendering"))
#endif
, writeZTransparent(initData(&writeZTransparent, (bool) false, "writeZTransparent", "Write into Z Buffer for Transparent Object"))
, alphaBlend(initData(&alphaBlend, (bool) false, "alphaBlend", "Enable alpha blending"))
, depthTest(initData(&depthTest, (bool) true, "depthTest", "Enable depth testing"))
, cullFace(initData(&cullFace, (int) 0, "cullFace", "Face culling (0 = no culling, 1 = cull back faces, 2 = cull front faces)"))
, blendEquation( initData(&blendEquation, "blendEquation", "if alpha blending is enabled this specifies how source and destination colors are combined") )
, sourceFactor( initData(&sourceFactor, "sfactor", "if alpha blending is enabled this specifies how the red, green, blue, and alpha source blending factors are computed") )
, destFactor( initData(&destFactor, "dfactor", "if alpha blending is enabled this specifies how the red, green, blue, and alpha destination blending factors are computed") )
, tex(NULL), canUseVBO(false), VBOGenDone(false), initDone(false), useTriangles(false), useQuads(false)
, oldTrianglesSize(0), oldQuadsSize(0)
{
    textures.clear();

	sofa::helper::OptionsGroup* blendEquationOptions = blendEquation.beginEdit();
    blendEquationOptions->setNames(4,"GL_FUNC_ADD", "GL_FUNC_SUBTRACT", "GL_MIN", "GL_MAX"); // .. add other options
    blendEquationOptions->setSelectedItem(0);
    //this->f_printLog.setValue(true);
    blendEquation.endEdit();

	// alpha blend values
	sofa::helper::OptionsGroup* sourceFactorOptions = sourceFactor.beginEdit();
    sourceFactorOptions->setNames(4,"GL_ZERO", "GL_ONE", "GL_SRC_ALPHA", "GL_ONE_MINUS_SRC_ALPHA"); // .. add other options
    sourceFactorOptions->setSelectedItem(2);
    //this->f_printLog.setValue(true);
    sourceFactor.endEdit();

	sofa::helper::OptionsGroup* destFactorOptions = destFactor.beginEdit();
    destFactorOptions->setNames(4,"GL_ZERO", "GL_ONE", "GL_SRC_ALPHA", "GL_ONE_MINUS_SRC_ALPHA"); // .. add other options
    destFactorOptions->setSelectedItem(3);
    //this->f_printLog.setValue(true);
    destFactor.endEdit();
}

OglModel::~OglModel()
{
    if (tex!=NULL) delete tex;

    for (unsigned int i = 0 ; i < textures.size() ; i++)
    {
        delete textures[i];
    }
}

void OglModel::drawGroup(int ig, bool transparent)
{
   const ResizableExtVector<Triangle>& triangles = this->getTriangles();
   const ResizableExtVector<Quad>& quads = this->getQuads();
   const ResizableExtVector<Coord>& vertices = this->getVertices();
   const ResizableExtVector<Deriv>& vnormals = this->getVnormals();

    FaceGroup g;
    if (ig < 0)
    {
        g.materialId = -1;
        g.tri0 = 0;
        g.nbt = triangles.size();
        g.quad0 = 0;
        g.nbq = quads.size();
    }
    else
    {
        g = this->groups.getValue()[ig];
    }
    Material m;
    if (g.materialId < 0)
        m = this->material.getValue();
    else
        m = this->materials.getValue()[g.materialId];

    if (!tex && m.useTexture && m.activated)
    {
        //get the texture id corresponding to the current material
        int indexInTextureArray = materialTextureIdMap[g.materialId];
        if (textures[indexInTextureArray])
        {
            textures[indexInTextureArray]->bind();
        }

        glEnable(GL_TEXTURE_2D);

        if(VBOGenDone && useVBO.getValue())
        {
            glBindBufferARB(GL_ARRAY_BUFFER, vbo);
            glTexCoordPointer(2, GL_FLOAT, 0, (char*)NULL + (vertices.size()*sizeof(vertices[0]))
                    + (vnormals.size()*sizeof(vnormals[0]))
                    );
            glBindBufferARB(GL_ARRAY_BUFFER, 0);
        }
        else
        {
            //get the texture coordinates
            const ResizableExtVector<TexCoord>& vtexcoords = this->getVtexcoords();
            glTexCoordPointer(2, GL_FLOAT, 0, vtexcoords.getData());
        }
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
//
//        if (hasTangents)
//        {
//            glClientActiveTexture(GL_TEXTURE1);
//            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
//            if(VBOGenDone && useVBO.getValue())
//            {
//                glBindBufferARB(GL_ARRAY_BUFFER, vbo);
//                glTexCoordPointer(3, GL_FLOAT, 0,
//                                  (char*)NULL + (vertices.size()*sizeof(vertices[0])) +
//                                  (vnormals.size()*sizeof(vnormals[0])) +
//                                  (vtexcoords.size()*sizeof(vtexcoords[0])));
//                glBindBufferARB(GL_ARRAY_BUFFER, 0);
//            }
//            else
//                glTexCoordPointer(3, GL_FLOAT, 0, vtangents.getData());
//
//            glClientActiveTexture(GL_TEXTURE2);
//            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
//            if(VBOGenDone && useVBO.getValue())
//            {
//                glBindBufferARB(GL_ARRAY_BUFFER, vbo);
//                glTexCoordPointer(3, GL_FLOAT, 0,
//                                  (char*)NULL + (vertices.size()*sizeof(vertices[0])) +
//                                  (vnormals.size()*sizeof(vnormals[0])) +
//                                  (vtexcoords.size()*sizeof(vtexcoords[0])) +
//                                  (vtangents.size()*sizeof(vtangents[0])));
//                glBindBufferARB(GL_ARRAY_BUFFER, 0);
//            }
//            else
//                glTexCoordPointer(3, GL_FLOAT, 0, vbitangents.getData());
//
//            glClientActiveTexture(GL_TEXTURE0);
//        }
    }


    bool isTransparent = (m.useDiffuse && m.diffuse[3] < 1.0);
    if (transparent ^ isTransparent) return;
    
    Vec4f ambient = m.useAmbient?m.ambient:Vec4f();
    Vec4f diffuse = m.useDiffuse?m.diffuse:Vec4f();
    Vec4f specular = m.useSpecular?m.specular:Vec4f();
    Vec4f emissive = m.useEmissive?m.emissive:Vec4f();
    float shininess = m.useShininess?m.shininess:45;

    if (shininess == 0.0f)
    {
        specular.clear();
        shininess = 1;
    }

    if (isTransparent) {
        emissive[3] = 0; //diffuse[3];
        ambient[3] = 0; //diffuse[3];
        //diffuse[3] = 0;
        specular[3] = 0;
    }
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, ambient.ptr());
    glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse.ptr());
    glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, specular.ptr());
    glMaterialfv (GL_FRONT_AND_BACK, GL_EMISSION, emissive.ptr());
    glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, shininess);
    if(VBOGenDone && useVBO.getValue())
    {
#ifdef SOFA_HAVE_GLEW
        if (g.nbt > 0)
        {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, iboTriangles);
            glDrawElements(GL_TRIANGLES, g.nbt * 3, GL_UNSIGNED_INT, (unsigned int*)NULL + (g.tri0 * 3));
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);
        }
        if (g.nbq > 0)
        {
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, iboQuads);
            glDrawElements(GL_QUADS, g.nbq * 4, GL_UNSIGNED_INT, (unsigned int*)NULL + (g.quad0 * 4));
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);
        }
#endif
    }
    else
    {
        if (g.nbt > 0)
            glDrawElements(GL_TRIANGLES, g.nbt * 3, GL_UNSIGNED_INT, triangles.getData() + g.tri0);
        if (g.nbq > 0)
            glDrawElements(GL_QUADS, g.nbq * 4, GL_UNSIGNED_INT, quads.getData() + g.quad0);
    }

    if (!tex && m.useTexture && m.activated)
    {
        int indexInTextureArray = materialTextureIdMap[g.materialId];
        if (textures[indexInTextureArray])
        {
            textures[indexInTextureArray]->unbind();
        }
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisable(GL_TEXTURE_2D);
    }
}

void OglModel::drawGroups(bool transparent)
{
    helper::ReadAccessor< Data< helper::vector<FaceGroup> > > groups = this->groups;

    //for (unsigned int i=0; i<xforms.size(); i++)
    {
    //    float matrix[16];
    //    xforms[i].writeOpenGlMatrix(matrix);
    //    pushTransformMatrix(matrix);

        if (groups.empty())
            drawGroup(-1, transparent);
        else
        {
            for (unsigned int i=0;i<groups.size();++i)
                drawGroup(i, transparent);
        }

    //    popTransformMatrix();
    }
}

void OglModel::internalDraw(const core::visual::VisualParams* vparams, bool transparent)
{
    m_vtexcoords.updateIfDirty();
//    serr<<" OglModel::internalDraw()"<<sendl;
    if (!vparams->displayFlags().getShowVisualModels()) return;

    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    const ResizableExtVector<Coord>& vertices = this->getVertices();
    const ResizableExtVector<Deriv>& vnormals = this->getVnormals();
    const ResizableExtVector<TexCoord>& vtexcoords= this->getVtexcoords();
    const ResizableExtVector<Coord>& vtangents= this->getVtangents();
    const ResizableExtVector<Coord>& vbitangents= this->getVbitangents();
    bool hasTangents = vtangents.size() && vbitangents.size();

    glEnable(GL_LIGHTING);

    //Enable<GL_BLEND> blending;
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glColor3f(1.0 , 1.0, 1.0);

    if(VBOGenDone && useVBO.getValue())
    {
//#ifdef SOFA_HAVE_GLEW
        glBindBufferARB(GL_ARRAY_BUFFER, vbo);

        glVertexPointer(3, GL_FLOAT, 0, (char*)NULL + 0);
        glNormalPointer(GL_FLOAT, 0, (char*)NULL + (vertices.size()*sizeof(vertices[0])));

        glBindBufferARB(GL_ARRAY_BUFFER, 0);
//#endif
    }
    else
    {
        glVertexPointer (3, GL_FLOAT, 0, vertices.getData());
        glNormalPointer (GL_FLOAT, 0, vnormals.getData());
    }

    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);

    if ((tex || putOnlyTexCoords.getValue()) )//&& !numberOfTextures)
    {
        glEnable(GL_TEXTURE_2D);
        if(tex)
            tex->bind();

        if(VBOGenDone && useVBO.getValue())
        {
            glBindBufferARB(GL_ARRAY_BUFFER, vbo);
            glTexCoordPointer(2, GL_FLOAT, 0, (char*)NULL + (vertices.size()*sizeof(vertices[0])) + (vnormals.size()*sizeof(vnormals[0])) );
            glBindBufferARB(GL_ARRAY_BUFFER, 0);
        }
        else
        {
            glTexCoordPointer(2, GL_FLOAT, 0, vtexcoords.getData());
        }
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);

        if (hasTangents)
        {
            glClientActiveTexture(GL_TEXTURE1);
            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
            if(VBOGenDone && useVBO.getValue())
            {
                glBindBufferARB(GL_ARRAY_BUFFER, vbo);
                glTexCoordPointer(3, GL_FLOAT, 0,
                                  (char*)NULL + (vertices.size()*sizeof(vertices[0])) +
                                  (vnormals.size()*sizeof(vnormals[0])) +
                                  (vtexcoords.size()*sizeof(vtexcoords[0])));
                glBindBufferARB(GL_ARRAY_BUFFER, 0);
            }
            else
                glTexCoordPointer(3, GL_FLOAT, 0, vtangents.getData());

            glClientActiveTexture(GL_TEXTURE2);
            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
            if(VBOGenDone && useVBO.getValue())
            {
                glBindBufferARB(GL_ARRAY_BUFFER, vbo);
                glTexCoordPointer(3, GL_FLOAT, 0,
                                  (char*)NULL + (vertices.size()*sizeof(vertices[0])) +
                                  (vnormals.size()*sizeof(vnormals[0])) +
                                  (vtexcoords.size()*sizeof(vtexcoords[0])) +
                                  (vtangents.size()*sizeof(vtangents[0])));
                glBindBufferARB(GL_ARRAY_BUFFER, 0);
            }
            else
                glTexCoordPointer(3, GL_FLOAT, 0, vbitangents.getData());

            glClientActiveTexture(GL_TEXTURE0);
        }
    }

    if (transparent)
    {
        glEnable(GL_BLEND);
        if (writeZTransparent.getValue())
            glDepthMask(GL_TRUE);
        else glDepthMask(GL_FALSE);

        glBlendFunc(GL_ZERO, GL_ONE_MINUS_SRC_ALPHA);

        drawGroups(transparent);

        if (premultipliedAlpha.getValue())
            glBlendFunc(GL_ONE, GL_ONE);
        else
            glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    }

    if (alphaBlend.getValue())
    {
        glDepthMask(GL_FALSE);
        glBlendEquation( blendEq );
        glBlendFunc( sfactor, dfactor );
        glEnable(GL_BLEND);
    }

    if (!depthTest.getValue())
        glDisable(GL_DEPTH_TEST);

    switch (cullFace.getValue())
    {
    case 1:
        glCullFace(GL_BACK);
        glEnable(GL_CULL_FACE);
        break;
    case 2:
        glCullFace(GL_FRONT);
        glEnable(GL_CULL_FACE);
        break;
    }

    drawGroups(transparent);

    switch (cullFace.getValue())
    {
    case 1:
    case 2:
        glDisable(GL_CULL_FACE);
        break;
    }

    if (!depthTest.getValue())
        glEnable(GL_DEPTH_TEST);

    if (alphaBlend.getValue())
    {
        // restore Default value
		glBlendEquation( GL_FUNC_ADD );
        glBlendFunc( GL_ONE, GL_ONE );
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
    }

    if ( (tex || putOnlyTexCoords.getValue()) )//&& !numberOfTextures)
    {
        if (tex)
            tex->unbind();
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisable(GL_TEXTURE_2D);

        if (hasTangents)
        {
            glClientActiveTexture(GL_TEXTURE1);
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
            glClientActiveTexture(GL_TEXTURE2);
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
            glClientActiveTexture(GL_TEXTURE0);
        }
    }
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisable(GL_LIGHTING);

    if (transparent) {
        glDisable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        //glBlendFunc(GL_ONE, GL_ZERO);
        glDepthMask(GL_TRUE);
    }

    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    if (vparams->displayFlags().getShowNormals())
    {
        glColor3f (1.0, 1.0, 1.0);
        for (unsigned int i=0; i<xforms.size(); i++)
        {
            float matrix[16];
            xforms[i].writeOpenGlMatrix(matrix);
            glPushMatrix();
            glMultMatrixf(matrix);

        glBegin(GL_LINES);
            for (unsigned int i = 0; i < vertices.size(); i++)
            {
                glVertex3fv (vertices[i].ptr());
                Coord p = vertices[i] + vnormals[i];
                glVertex3fv (p.ptr());
            }
        glEnd();

            glPopMatrix();
        }
    }
    m_vtexcoords.updateIfDirty();
}

bool OglModel::loadTexture(const std::string& filename)
{
    helper::io::Image *img = helper::io::Image::Create(filename);
    if (!img)
        return false;
    tex = new helper::gl::Texture(img, true, true, false, srgbTexturing.getValue());
    return true;
}

// a false result indicated problems during textures loading
bool OglModel::loadTextures()
{
    bool result = true;
    textures.clear();

    //for output
    std::string str;

    //count the total number of activated textures
    std::vector<unsigned int> activatedTextures;
    for (unsigned int i = 0 ; i < this->materials.getValue().size() ; i++)
        if (this->materials.getValue()[i].useTexture && this->materials.getValue()[i].activated)
            activatedTextures.push_back(i);

    for (std::vector< unsigned int>::iterator i = activatedTextures.begin() ; i < activatedTextures.end(); i++)
    {
        std::string textureFile(this->materials.getValue()[*i].textureFilename);

        if (!sofa::helper::system::DataRepository.findFile(textureFile))
        {
            textureFile = this->fileMesh.getFullPath();
            unsigned int position = textureFile.rfind("/");
            textureFile.replace (position+1,textureFile.length() - position, this->materials.getValue()[*i].textureFilename);

            if (!sofa::helper::system::DataRepository.findFile(textureFile))
            {
                if (i!=activatedTextures.begin()) std::cout << "\n" << std::flush;
                std::cerr   << "ERROR(OglModel(" << this->getName() << "):Texture \"" << this->materials.getValue()[*i].textureFilename << "\" not found"
                        << " in material " << this->materials.getValue()[*i].name <<  std::endl;
                result = false;
                continue;
            }
        }

        helper::io::Image *img = helper::io::Image::Create(textureFile);
        if (!img)
        {
           std::cerr << "ERROR(OglModel(" << this->getName() << "): couldn't create an image from file " << this->materials.getValue()[*i].textureFilename << std::endl;
           result = false;
           continue;
        }
        helper::gl::Texture * text = new helper::gl::Texture(img, true, true, false, srgbTexturing.getValue());
        materialTextureIdMap.insert(std::pair<int, int>(*i,textures.size()));
        textures.push_back( text );

        /****OUTPUT****/
        //clear the previous line
        std::cout << std::string(str.length(), '\b');

        str.clear();

        std::stringstream out; out << textures.size() << "/" << activatedTextures.size() << " textures loaded for OglModel " << this->getName();
        str += out.str();

        std::cout << str << std::flush;
    }
    std::cout << std::endl;

    if (textures.size() != activatedTextures.size())
        std::cerr << "ERROR(OglModel(" << this->getName()<< ")) " << (activatedTextures.size() - textures.size()) << " textures couldn't be loaded" <<  std::endl;



    /**********************************************
     * Load textures for bump mapping
     *********************************************/
//
//    for (unsigned int i = 0 ; i < this->materials.getValue().size() ; i++)
//    {
//       //we count only the bump texture with an activated material
//       if (this->materials.getValue()[i].useBumpMapping && this->materials.getValue()[i].activated)
//       {
//            std::string textureFile(this->materials.getValue()[i].bumpTextureFilename);
//
//            if (!sofa::helper::system::DataRepository.findFile(textureFile))
//            {
//                textureFile = this->fileMesh.getFullPath();
//                unsigned int position = textureFile.rfind("/");
//                textureFile.replace (position+1,textureFile.length() - position, this->materials.getValue()[i].bumpTextureFilename);
////                std::cout << "Loading texture: " << textureFile << std::endl;
//
//                if (!sofa::helper::system::DataRepository.findFile(textureFile))
//                {
//                    std::cout <<  std::endl;
//                    serr << "Texture \"" << this->materials.getValue()[i].bumpTextureFilename << "\" not found"
//                            << " in material " << this->materials.getValue()[i].name << " for OglModel " << this->name
//                            << "(\""<< this->fileMesh.getFullPath() << "\")" << sendl;
//                    break;
//                }
//            }
//
//            helper::io::Image *img = helper::io::Image::Create(textureFile);
//            if (!img)
//            {
//                std::cout <<  std::endl;
//               std::cerr << "Error:OglModel:loadTextures: couldn't create an image from file " << this->materials.getValue()[i].bumpTextureFilename << std::endl;
//               return false;
//            }
//            helper::gl::Texture * text = new helper::gl::Texture(img, true, true, false, srgbTexturing.getValue());
//            materialTextureIdMap.insert(std::pair<int, int>(i,textures.size()));
//            textures.push_back( text );
//
//            std::cout << "\r\033[K" << i+1 << "/" << this->materials.getValue().size() << " textures loaded for bump mapping for OglModel " << this->getName()
//                    << "(loading "<<textureFile << ")"<< std::flush;
//       }
//    }
//    std::cout << "\r\033[K" << std::flush;
    return result;
}

void OglModel::initVisual()
{
    initTextures();

    initDone = true;
#ifdef NO_VBO
    canUseVBO = false;
#else
#ifdef SOFA_HAVE_GLEW
    //This test is not enough to detect if we can enable the VBO.
    canUseVBO = (GLEW_ARB_vertex_buffer_object!=0);
#endif

    if (useVBO.getValue() && !canUseVBO)
    {
                serr << "OglModel : VBO is not supported by your GPU" << sendl;
    }
#endif

    updateBuffers();

	if ( alphaBlend.getValue() )
	{
		blendEq = getGLenum( blendEquation.getValue().getSelectedItem().c_str() );
		sfactor = getGLenum( sourceFactor.getValue().getSelectedItem().c_str() );
		dfactor = getGLenum( destFactor.getValue().getSelectedItem().c_str() );
	}

}

void OglModel::initTextures()
{
    if (tex)
    {
        tex->init();
    }
    else
    {
        if (!textures.empty())
        {
            for (unsigned int i = 0 ; i < textures.size() ; i++)
            {
                textures[i]->init();
            }
        }
    }
}
#ifdef SOFA_HAVE_GLEW
void OglModel::createVertexBuffer()
{


    glGenBuffersARB(1, &vbo);
    initVertexBuffer();
    VBOGenDone = true;
}

void OglModel::createTrianglesIndicesBuffer()
{
    glGenBuffersARB(1, &iboTriangles);
    initTrianglesIndicesBuffer();
    useTriangles = true;
}


void OglModel::createQuadsIndicesBuffer()
{
    glGenBuffersARB(1, &iboQuads);
    initQuadsIndicesBuffer();
    useQuads = true;
}


void OglModel::initVertexBuffer()
{
    unsigned positionsBufferSize, normalsBufferSize;
    unsigned textureCoordsBufferSize = 0, tangentsBufferSize = 0, bitangentsBufferSize = 0;
   const ResizableExtVector<Coord>& vertices = this->getVertices();
   const ResizableExtVector<Coord>& vnormals = this->getVnormals();
   const ResizableExtVector<TexCoord>& vtexcoords= this->getVtexcoords();
    const ResizableExtVector<Coord>& vtangents= this->getVtangents();
    const ResizableExtVector<Coord>& vbitangents= this->getVbitangents();
    bool hasTangents = vtangents.size() && vbitangents.size();

    positionsBufferSize = (vertices.size()*sizeof(vertices[0]));
    normalsBufferSize = (vnormals.size()*sizeof(vnormals[0]));
    if (tex || putOnlyTexCoords.getValue() || !textures.empty())
    {
        textureCoordsBufferSize = vtexcoords.size() * sizeof(vtexcoords[0]);

        if (hasTangents)
        {
            tangentsBufferSize = vtangents.size() * sizeof(vtangents[0]);
            bitangentsBufferSize = vbitangents.size() * sizeof(vbitangents[0]);
        }
    }

    unsigned int totalSize = positionsBufferSize + normalsBufferSize + textureCoordsBufferSize +
                             tangentsBufferSize + bitangentsBufferSize;

    glBindBufferARB(GL_ARRAY_BUFFER, vbo);
    //Vertex Buffer creation
    glBufferDataARB(GL_ARRAY_BUFFER,
                 totalSize,
                 NULL,
                 GL_DYNAMIC_DRAW);


    updateVertexBuffer();

    glBindBufferARB(GL_ARRAY_BUFFER, 0);
}


void OglModel::initTrianglesIndicesBuffer()
{
   const ResizableExtVector<Triangle>& triangles = this->getTriangles();

    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, iboTriangles);

    glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER, triangles.size()*sizeof(triangles[0]), NULL, GL_DYNAMIC_DRAW);
    updateTrianglesIndicesBuffer();

    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void OglModel::initQuadsIndicesBuffer()
{
   const ResizableExtVector<Quad>& quads = this->getQuads();

    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, iboQuads);
    glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER, quads.size()*sizeof(quads[0]), NULL, GL_DYNAMIC_DRAW);
    updateQuadsIndicesBuffer();
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void OglModel::updateVertexBuffer()
{
   const ResizableExtVector<Coord>& vertices = this->getVertices();
   const ResizableExtVector<Coord>& vnormals = this->getVnormals();
   const ResizableExtVector<TexCoord>& vtexcoords= this->getVtexcoords();
    const ResizableExtVector<Coord>& vtangents= this->getVtangents();
    const ResizableExtVector<Coord>& vbitangents= this->getVbitangents();
    bool hasTangents = vtangents.size() && vbitangents.size();

    unsigned positionsBufferSize, normalsBufferSize;
    unsigned textureCoordsBufferSize = 0, tangentsBufferSize = 0, bitangentsBufferSize = 0;

    positionsBufferSize = (vertices.size()*sizeof(vertices[0]));
    normalsBufferSize = (vnormals.size()*sizeof(vnormals[0]));
    if (tex || putOnlyTexCoords.getValue() || !textures.empty())
    {
        textureCoordsBufferSize = vtexcoords.size() * sizeof(vtexcoords[0]);

        if (hasTangents)
        {
            tangentsBufferSize = vtangents.size() * sizeof(vtangents[0]);
            bitangentsBufferSize = vbitangents.size() * sizeof(vbitangents[0]);
        }
    }

    glBindBufferARB(GL_ARRAY_BUFFER, vbo);
    //Positions
    glBufferSubDataARB(GL_ARRAY_BUFFER,
                    0,
                    positionsBufferSize,
                    vertices.getData());

    //Normals
    glBufferSubDataARB(GL_ARRAY_BUFFER,
                    positionsBufferSize,
                    normalsBufferSize,
                    vnormals.getData());

    //Texture coords
    if(tex || putOnlyTexCoords.getValue() ||!textures.empty())
    {
        glBufferSubDataARB(GL_ARRAY_BUFFER,
                        positionsBufferSize + normalsBufferSize,
                        textureCoordsBufferSize,
                        vtexcoords.getData());

        if (hasTangents)
        {
            glBufferSubDataARB(GL_ARRAY_BUFFER,
                            positionsBufferSize + normalsBufferSize + textureCoordsBufferSize,
                            tangentsBufferSize,
                            vtangents.getData());

            glBufferSubDataARB(GL_ARRAY_BUFFER,
                            positionsBufferSize + normalsBufferSize + textureCoordsBufferSize + tangentsBufferSize,
                            bitangentsBufferSize,
                            vbitangents.getData());
        }
    }

    glBindBufferARB(GL_ARRAY_BUFFER, 0);

}

void OglModel::updateTrianglesIndicesBuffer()
{
   const ResizableExtVector<Triangle>& triangles = this->getTriangles();
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, iboTriangles);
    glBufferSubDataARB(GL_ELEMENT_ARRAY_BUFFER, 0, triangles.size()*sizeof(triangles[0]), &triangles[0]);
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void OglModel::updateQuadsIndicesBuffer()
{
   const ResizableExtVector<Quad>& quads = this->getQuads();
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, iboQuads);
    glBufferSubDataARB(GL_ELEMENT_ARRAY_BUFFER, 0, quads.size()*sizeof(quads[0]), &quads[0]);
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER, 0);
}
#endif
void OglModel::updateBuffers()
{
   const ResizableExtVector<Triangle>& triangles = this->getTriangles();
   const ResizableExtVector<Quad>& quads = this->getQuads();
   const ResizableExtVector<Coord>& vertices = this->getVertices();

    if (initDone)
    {
#ifdef SOFA_HAVE_GLEW
        if (useVBO.getValue() && canUseVBO)
        {
            if(!VBOGenDone)
            {
                createVertexBuffer();
                //Index Buffer Object
                //Triangles indices
                if(triangles.size() > 0)
                    createTrianglesIndicesBuffer();
                //Quads indices
                if(quads.size() > 0)
                    createQuadsIndicesBuffer();
            }
            //Update VBO & IBO
            else
            {
                if(oldVerticesSize != vertices.size())
                    initVertexBuffer();
                else
                    updateVertexBuffer();
                //Indices
                //Triangles
                if(useTriangles)
                    if(oldTrianglesSize != triangles.size())
                        initTrianglesIndicesBuffer();
                    else
                        updateTrianglesIndicesBuffer();
                else if (triangles.size() > 0)
                    createTrianglesIndicesBuffer();

                //Quads
                if (useQuads)
                    if(oldQuadsSize != quads.size())
                        initQuadsIndicesBuffer();
                    else
                        updateQuadsIndicesBuffer();
                else if (quads.size() > 0)
                    createQuadsIndicesBuffer();
            }
            oldVerticesSize = vertices.size();
			oldTrianglesSize = triangles.size();
			oldQuadsSize = quads.size();
		}
#endif
	}

}


GLenum OglModel::getGLenum(const char* c ) const
{

	if ( strcmp( c, "GL_ZERO") == 0)
	{
		return GL_ZERO;
	}
	else if  ( strcmp( c, "GL_ONE") == 0)
	{
		return GL_ONE;
	}
	else if (strcmp( c, "GL_SRC_ALPHA") == 0 )
	{
		return GL_SRC_ALPHA;
	}
	else if (strcmp( c, "GL_ONE_MINUS_SRC_ALPHA") == 0 )
	{
		return GL_ONE_MINUS_SRC_ALPHA;
	}
	// .... add ohter OGL symbolic constants
	// glBlendEquation Value
	else if  ( strcmp( c, "GL_FUNC_ADD") == 0)
	{
		return GL_FUNC_ADD;
	}
	else if (strcmp( c, "GL_FUNC_SUBTRACT") == 0 )
	{
		return GL_FUNC_SUBTRACT;
	}
	else if (strcmp( c, "GL_MAX") == 0 )
	{
		return GL_MAX;
	}
	else if (strcmp( c, "GL_MIN") == 0 )
	{
		return GL_MIN;
	}

	else
	{
		// error: not valid
		std::cerr   << " OglModel - not valid or not supported openGL enum value: " << c ;
		return GL_ZERO;
	}


}


} // namespace visualmodel

} // namespace component

} // namespace sofa

