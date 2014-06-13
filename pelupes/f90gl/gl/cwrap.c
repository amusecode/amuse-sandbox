/*
!----------------------------------------------------------------------------
!
!                    f90gl version 1.2.12, September 2006
!
!     William F. Mitchell
!     Mathematical and Computational Sciences Division
!     National Institute of Standards and Technology
!     william.mitchell@nist.gov
!     http://math.nist.gov/f90gl
!
! f90gl is public domain software.  It was produced as part of work done by
! the U.S. Government, and is not subject to copyright in the United States.
!
!----------------------------------------------------------------------------
*/

#ifdef LF90
#define LF9X
#endif
#ifdef LF95
#define LF9X
#endif

#ifdef WIN32
#include <windows.h>
#else
#define APIENTRY
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(MACOSX)
#include <OpenGL/gl.h>
#include <OpenGL/glext.h>
#else
#include <GL/gl.h>
#endif

/* Make sure the extensions and version number are defined the same
   way in both fwrap and cwrap.  Not yet supported for Windows. */
#ifndef WIN32
#include "fppr.inc"
#endif

#ifdef ABS32
#define APIENTRY
#endif

#ifdef PGI32
#define APIENTRY
#endif

#ifdef LF9X
#define APIENTRY
#endif

#ifdef IFORT
#undef APIENTRY
#define APIENTRY
#endif

/* Fortran name mangling choices */
#define LOWERCASE 1
#define UNDERSCORE 2
#define UPPERCASE 3

#ifdef F_LANG
#define PURE_FUNCTIONS
#endif

/* used when Fortran must pass an assumed shape array */
#ifdef F_LANG
#define ARRAYREF [0]
#define ARRAY *
#else
#define ARRAYREF
#define ARRAY
#endif

/* Replacements for "strlen", "malloc", "free" and "ceil" which */
/* use the Windows API insead of standard C library   */
/* Thank you to Lawson Wakefield (I.S.S., Ltd.) for supplying these. */

#ifdef REPLACE_CEIL
#define strlen lstrlen

double ceil(double x)
{
double y;

    y = (double) (int) x;
    if (y < x) y++;
    return y;
}
#endif

#ifdef REPLACE_MALLOC
void* malloc(size_t nBytes)
{
    return (void*) GlobalAlloc(GPTR,(DWORD) nBytes);
}

void free(void* buffer)
{
    GlobalFree((HGLOBAL) buffer);
}

void exit(int code)
{
    ExitProcess(code);
}
#endif

/*
 * Mapping of Fortran types to C types:
 */

#define f90glboolean  GLboolean *
#define f90glbyte     GLbyte *
#define f90glubyte    GLubyte *
#define f90glshort    GLshort *
#define f90glushort   GLushort *
#define f90glint      GLint *
#define f90gluint     GLuint *
#define f90glenum     GLenum *
#define f90glbitfield GLbitfield *
#define f90glsizei    GLsizei *
#define f90glfloat    GLfloat *
#define f90glclampf   GLclampf *
#define f90gldouble   GLdouble *
#define f90glclampd   GLclampd *
#define f90glvoid     int *
#define f90glcint     int *
#define f90glcuint    unsigned int *
#define f90glcchar    void *
#define f90glculong   unsigned long *

#define end_of_string '\0'

/* procedure name mangling */

/* OpenGL 1.0 */

#if FNAME==LOWERCASE
#define fglAccum                       fglaccum
#define fglAlphaFunc                   fglalphafunc
#define fglBegin                       fglbegin
#define fglBitmap                      fglbitmap
#define f9bglBitmap                    f9bglbitmap
#define f9cglBitmap                    f9cglbitmap
#define fglBlendFunc                   fglblendfunc
#define fglCallList                    fglcalllist
#define f9y0glCallLists                f9y0glcalllists
#define f9y1glCallLists                f9y1glcalllists
#define f9y4glCallLists                f9y4glcalllists
#define f9y5glCallLists                f9y5glcalllists
#define f9y6glCallLists                f9y6glcalllists
#define fglClear                       fglclear
#define fglClearAccum                  fglclearaccum
#define fglClearColor                  fglclearcolor
#define fglClearDepth                  fglcleardepth
#define fglClearIndex                  fglclearindex
#define fglClearStencil                fglclearstencil
#define fglClipPlane                   fglclipplane
#define fglColor3b                     fglcolor3b
#define f9bglColor3b                   f9bglcolor3b
#define f9cglColor3b                   f9cglcolor3b
#define fglColor3bv                    fglcolor3bv
#define f9bglColor3bv                  f9bglcolor3bv
#define f9cglColor3bv                  f9cglcolor3bv
#define fglColor3d                     fglcolor3d
#define fglColor3dv                    fglcolor3dv
#define fglColor3f                     fglcolor3f
#define fglColor3fv                    fglcolor3fv
#define fglColor3i                     fglcolor3i
#define fglColor3iv                    fglcolor3iv
#define fglColor3s                     fglcolor3s
#define f9aglColor3s                   f9aglcolor3s
#define fglColor3sv                    fglcolor3sv
#define f9aglColor3sv                  f9aglcolor3sv
#define fglColor3ub                    fglcolor3ub
#define f9bglColor3ub                  f9bglcolor3ub
#define f9cglColor3ub                  f9cglcolor3ub
#define fglColor3ubv                   fglcolor3ubv
#define f9bglColor3ubv                 f9bglcolor3ubv
#define f9cglColor3ubv                 f9cglcolor3ubv
#define fglColor3ui                    fglcolor3ui
#define fglColor3uiv                   fglcolor3uiv
#define fglColor3us                    fglcolor3us
#define f9aglColor3us                  f9aglcolor3us
#define fglColor3usv                   fglcolor3usv
#define f9aglColor3usv                 f9aglcolor3usv
#define fglColor4b                     fglcolor4b
#define f9bglColor4b                   f9bglcolor4b
#define f9cglColor4b                   f9cglcolor4b
#define fglColor4bv                    fglcolor4bv
#define f9bglColor4bv                  f9bglcolor4bv
#define f9cglColor4bv                  f9cglcolor4bv
#define fglColor4d                     fglcolor4d
#define fglColor4dv                    fglcolor4dv
#define fglColor4f                     fglcolor4f
#define fglColor4fv                    fglcolor4fv
#define fglColor4i                     fglcolor4i
#define fglColor4iv                    fglcolor4iv
#define fglColor4s                     fglcolor4s
#define f9aglColor4s                   f9aglcolor4s
#define fglColor4sv                    fglcolor4sv
#define f9aglColor4sv                  f9aglcolor4sv
#define fglColor4ub                    fglcolor4ub
#define f9bglColor4ub                  f9bglcolor4ub
#define f9cglColor4ub                  f9cglcolor4ub
#define fglColor4ubv                   fglcolor4ubv
#define f9bglColor4ubv                 f9bglcolor4ubv
#define f9cglColor4ubv                 f9cglcolor4ubv
#define fglColor4ui                    fglcolor4ui
#define fglColor4uiv                   fglcolor4uiv
#define fglColor4us                    fglcolor4us
#define f9aglColor4us                  f9aglcolor4us
#define fglColor4usv                   fglcolor4usv
#define f9aglColor4usv                 f9aglcolor4usv
#define fglColorMask                   fglcolormask
#define f9eglColorMask                 f9eglcolormask
#define fglColorMaterial               fglcolormaterial
#define fglCopyPixels                  fglcopypixels
#define fglCullFace                    fglcullface
#define fglDeleteLists                 fgldeletelists
#define fglDepthFunc                   fgldepthfunc
#define fglDepthMask                   fgldepthmask
#define f9eglDepthMask                 f9egldepthmask
#define fglDepthRange                  fgldepthrange
#define fglDisable                     fgldisable
#define fglDrawBuffer                  fgldrawbuffer
#define f9y0glDrawPixels               f9y0gldrawpixels
#define f9y1glDrawPixels               f9y1gldrawpixels
#define f9y4glDrawPixels               f9y4gldrawpixels
#define f9y5glDrawPixels               f9y5gldrawpixels
#define f9y6glDrawPixels               f9y6gldrawpixels
#define fglEdgeFlag                    fgledgeflag
#define f9eglEdgeFlag                  f9egledgeflag
#define fglEdgeFlagv                   fgledgeflagv
#define f9eglEdgeFlagv                 f9egledgeflagv
#define fglEnable                      fglenable
#define fglEnd                         fglend
#define fglEndList                     fglendlist
#define fglEvalCoord1d                 fglevalcoord1d
#define fglEvalCoord1dv                fglevalcoord1dv
#define fglEvalCoord1f                 fglevalcoord1f
#define fglEvalCoord1fv                fglevalcoord1fv
#define fglEvalCoord2d                 fglevalcoord2d
#define fglEvalCoord2dv                fglevalcoord2dv
#define fglEvalCoord2f                 fglevalcoord2f
#define fglEvalCoord2fv                fglevalcoord2fv
#define fglEvalMesh1                   fglevalmesh1
#define fglEvalMesh2                   fglevalmesh2
#define fglEvalPoint1                  fglevalpoint1
#define fglEvalPoint2                  fglevalpoint2
#define fglFeedbackBuffer              fglfeedbackbuffer
#define fglFinish                      fglfinish
#define fglFlush                       fglflush
#define fglFogf                        fglfogf
#define fglFogfv                       fglfogfv
#define fglFogi                        fglfogi
#define fglFogiv                       fglfogiv
#define fglFrontFace                   fglfrontface
#define fglFrustum                     fglfrustum
#define fglGenLists                    fglgenlists
#define fglGetBooleanv                 fglgetbooleanv
#define f9eglGetBooleanv               f9eglgetbooleanv
#define fglGetClipPlane                fglgetclipplane
#define fglGetDoublev                  fglgetdoublev
#define fglGetError                    fglgeterror
#define fglGetFloatv                   fglgetfloatv
#define fglGetIntegerv                 fglgetintegerv
#define fglGetLightfv                  fglgetlightfv
#define fglGetLightiv                  fglgetlightiv
#define fglGetMapdv                    fglgetmapdv
#define fglGetMapfv                    fglgetmapfv
#define fglGetMapiv                    fglgetmapiv
#define fglGetMaterialfv               fglgetmaterialfv
#define fglGetMaterialiv               fglgetmaterialiv
#define fglGetPixelMapfv               fglgetpixelmapfv
#define fglGetPixelMapuiv              fglgetpixelmapuiv
#define fglGetPixelMapusv              fglgetpixelmapusv
#define f9aglGetPixelMapusv            f9aglgetpixelmapusv
#define fglGetPolygonStipple           fglgetpolygonstipple
#define f9bglGetPolygonStipple         f9bglgetpolygonstipple
#define f9cglGetPolygonStipple         f9cglgetpolygonstipple
#define f9yglGetStringLen              f9yglgetstringlen
#define f9yglGetString                 f9yglgetstring
#define f9yglGetString1                f9yglgetstring1
#define f9yglGetStringFree             f9yglgetstringfree
#define fglGetTexEnvfv                 fglgettexenvfv
#define fglGetTexEnviv                 fglgettexenviv
#define fglGetTexGendv                 fglgettexgendv
#define fglGetTexGenfv                 fglgettexgenfv
#define fglGetTexGeniv                 fglgettexgeniv
#define f9y0glGetTexImage              f9y0glgetteximage
#define f9y1glGetTexImage              f9y1glgetteximage
#define f9y4glGetTexImage              f9y4glgetteximage
#define f9y5glGetTexImage              f9y5glgetteximage
#define f9y6glGetTexImage              f9y6glgetteximage
#define fglGetTexLevelParameterfv      fglgettexlevelparameterfv
#define fglGetTexLevelParameteriv      fglgettexlevelparameteriv
#define fglGetTexParameterfv           fglgettexparameterfv
#define fglGetTexParameteriv           fglgettexparameteriv
#define fglHint                        fglhint
#define fglIndexMask                   fglindexmask
#define fglIndexd                      fglindexd
#define fglIndexdv                     fglindexdv
#define fglIndexf                      fglindexf
#define fglIndexfv                     fglindexfv
#define fglIndexi                      fglindexi
#define fglIndexiv                     fglindexiv
#define fglIndexs                      fglindexs
#define f9aglIndexs                    f9aglindexs
#define fglIndexsv                     fglindexsv
#define f9aglIndexsv                   f9aglindexsv
#define fglInitNames                   fglinitnames
#define fglIsEnabled                   fglisenabled
#define f9eglIsEnabled                 f9eglisenabled
#define fglIsList                      fglislist
#define f9eglIsList                    f9eglislist
#define fglLightModelf                 fgllightmodelf
#define fglLightModelfv                fgllightmodelfv
#define fglLightModeli                 fgllightmodeli
#define fglLightModeliv                fgllightmodeliv
#define fglLightf                      fgllightf
#define fglLightfv                     fgllightfv
#define fglLighti                      fgllighti
#define fglLightiv                     fgllightiv
#define fglLineStipple                 fgllinestipple
#define f9aglLineStipple               f9agllinestipple
#define fglLineWidth                   fgllinewidth
#define fglListBase                    fgllistbase
#define fglLoadIdentity                fglloadidentity
#define fglLoadMatrixd                 fglloadmatrixd
#define fglLoadMatrixf                 fglloadmatrixf
#define fglLoadName                    fglloadname
#define fglLogicOp                     fgllogicop
#define fglMap1d                       fglmap1d
#define fglMap1f                       fglmap1f
#define fglMap2d                       fglmap2d
#define fglMap2f                       fglmap2f
#define fglMapGrid1d                   fglmapgrid1d
#define fglMapGrid1f                   fglmapgrid1f
#define fglMapGrid2d                   fglmapgrid2d
#define fglMapGrid2f                   fglmapgrid2f
#define fglMaterialf                   fglmaterialf
#define fglMaterialfv                  fglmaterialfv
#define fglMateriali                   fglmateriali
#define fglMaterialiv                  fglmaterialiv
#define fglMatrixMode                  fglmatrixmode
#define fglMultMatrixd                 fglmultmatrixd
#define fglMultMatrixf                 fglmultmatrixf
#define fglNewList                     fglnewlist
#define fglNormal3b                    fglnormal3b
#define f9bglNormal3b                  f9bglnormal3b
#define f9cglNormal3b                  f9cglnormal3b
#define fglNormal3bv                   fglnormal3bv
#define f9bglNormal3bv                 f9bglnormal3bv
#define f9cglNormal3bv                 f9cglnormal3bv
#define fglNormal3d                    fglnormal3d
#define fglNormal3dv                   fglnormal3dv
#define fglNormal3f                    fglnormal3f
#define fglNormal3fv                   fglnormal3fv
#define fglNormal3i                    fglnormal3i
#define fglNormal3iv                   fglnormal3iv
#define fglNormal3s                    fglnormal3s
#define f9aglNormal3s                  f9aglnormal3s
#define fglNormal3sv                   fglnormal3sv
#define f9aglNormal3sv                 f9aglnormal3sv
#define fglOrtho                       fglortho
#define fglPassThrough                 fglpassthrough
#define fglPixelMapfv                  fglpixelmapfv
#define fglPixelMapuiv                 fglpixelmapuiv
#define fglPixelMapusv                 fglpixelmapusv
#define f9aglPixelMapusv               f9aglpixelmapusv
#define fglPixelStoref                 fglpixelstoref
#define fglPixelStorei                 fglpixelstorei
#define fglPixelTransferf              fglpixeltransferf
#define fglPixelTransferi              fglpixeltransferi
#define fglPixelZoom                   fglpixelzoom
#define fglPointSize                   fglpointsize
#define fglPolygonMode                 fglpolygonmode
#define fglPolygonStipple              fglpolygonstipple
#define f9bglPolygonStipple            f9bglpolygonstipple
#define f9cglPolygonStipple            f9cglpolygonstipple
#define fglPopAttrib                   fglpopattrib
#define fglPopMatrix                   fglpopmatrix
#define fglPopName                     fglpopname
#define fglPushAttrib                  fglpushattrib
#define fglPushMatrix                  fglpushmatrix
#define fglPushName                    fglpushname
#define fglRasterPos2d                 fglrasterpos2d
#define fglRasterPos2dv                fglrasterpos2dv
#define fglRasterPos2f                 fglrasterpos2f
#define fglRasterPos2fv                fglrasterpos2fv
#define fglRasterPos2i                 fglrasterpos2i
#define fglRasterPos2iv                fglrasterpos2iv
#define fglRasterPos2s                 fglrasterpos2s
#define f9aglRasterPos2s               f9aglrasterpos2s
#define fglRasterPos2sv                fglrasterpos2sv
#define f9aglRasterPos2sv              f9aglrasterpos2sv
#define fglRasterPos3d                 fglrasterpos3d
#define fglRasterPos3dv                fglrasterpos3dv
#define fglRasterPos3f                 fglrasterpos3f
#define fglRasterPos3fv                fglrasterpos3fv
#define fglRasterPos3i                 fglrasterpos3i
#define fglRasterPos3iv                fglrasterpos3iv
#define fglRasterPos3s                 fglrasterpos3s
#define f9aglRasterPos3s               f9aglrasterpos3s
#define fglRasterPos3sv                fglrasterpos3sv
#define f9aglRasterPos3sv              f9aglrasterpos3sv
#define fglRasterPos4d                 fglrasterpos4d
#define fglRasterPos4dv                fglrasterpos4dv
#define fglRasterPos4f                 fglrasterpos4f
#define fglRasterPos4fv                fglrasterpos4fv
#define fglRasterPos4i                 fglrasterpos4i
#define fglRasterPos4iv                fglrasterpos4iv
#define fglRasterPos4s                 fglrasterpos4s
#define f9aglRasterPos4s               f9aglrasterpos4s
#define fglRasterPos4sv                fglrasterpos4sv
#define f9aglRasterPos4sv              f9aglrasterpos4sv
#define fglReadBuffer                  fglreadbuffer
#define f9y0glReadPixels               f9y0glreadpixels
#define f9y1glReadPixels               f9y1glreadpixels
#define f9y4glReadPixels               f9y4glreadpixels
#define f9y5glReadPixels               f9y5glreadpixels
#define f9y6glReadPixels               f9y6glreadpixels
#define fglRectd                       fglrectd
#define fglRectdv                      fglrectdv
#define fglRectf                       fglrectf
#define fglRectfv                      fglrectfv
#define fglRecti                       fglrecti
#define fglRectiv                      fglrectiv
#define fglRects                       fglrects
#define f9aglRects                     f9aglrects
#define fglRectsv                      fglrectsv
#define f9aglRectsv                    f9aglrectsv
#define fglRenderMode                  fglrendermode
#define fglRotated                     fglrotated
#define fglRotatef                     fglrotatef
#define fglScaled                      fglscaled
#define fglScalef                      fglscalef
#define fglScissor                     fglscissor
#define fglSelectBuffer                fglselectbuffer
#define fglShadeModel                  fglshademodel
#define fglStencilFunc                 fglstencilfunc
#define fglStencilMask                 fglstencilmask
#define fglStencilOp                   fglstencilop
#define fglTexCoord1d                  fgltexcoord1d
#define fglTexCoord1dv                 fgltexcoord1dv
#define fglTexCoord1f                  fgltexcoord1f
#define fglTexCoord1fv                 fgltexcoord1fv
#define fglTexCoord1i                  fgltexcoord1i
#define fglTexCoord1iv                 fgltexcoord1iv
#define fglTexCoord1s                  fgltexcoord1s
#define f9aglTexCoord1s                f9agltexcoord1s
#define fglTexCoord1sv                 fgltexcoord1sv
#define f9aglTexCoord1sv               f9agltexcoord1sv
#define fglTexCoord2d                  fgltexcoord2d
#define fglTexCoord2dv                 fgltexcoord2dv
#define fglTexCoord2f                  fgltexcoord2f
#define fglTexCoord2fv                 fgltexcoord2fv
#define fglTexCoord2i                  fgltexcoord2i
#define fglTexCoord2iv                 fgltexcoord2iv
#define fglTexCoord2s                  fgltexcoord2s
#define f9aglTexCoord2s                f9agltexcoord2s
#define fglTexCoord2sv                 fgltexcoord2sv
#define f9aglTexCoord2sv               f9agltexcoord2sv
#define fglTexCoord3d                  fgltexcoord3d
#define fglTexCoord3dv                 fgltexcoord3dv
#define fglTexCoord3f                  fgltexcoord3f
#define fglTexCoord3fv                 fgltexcoord3fv
#define fglTexCoord3i                  fgltexcoord3i
#define fglTexCoord3iv                 fgltexcoord3iv
#define fglTexCoord3s                  fgltexcoord3s
#define f9aglTexCoord3s                f9agltexcoord3s
#define fglTexCoord3sv                 fgltexcoord3sv
#define f9aglTexCoord3sv               f9agltexcoord3sv
#define fglTexCoord4d                  fgltexcoord4d
#define fglTexCoord4dv                 fgltexcoord4dv
#define fglTexCoord4f                  fgltexcoord4f
#define fglTexCoord4fv                 fgltexcoord4fv
#define fglTexCoord4i                  fgltexcoord4i
#define fglTexCoord4iv                 fgltexcoord4iv
#define fglTexCoord4s                  fgltexcoord4s
#define f9aglTexCoord4s                f9agltexcoord4s
#define fglTexCoord4sv                 fgltexcoord4sv
#define f9aglTexCoord4sv               f9agltexcoord4sv
#define fglTexEnvf                     fgltexenvf
#define fglTexEnvfv                    fgltexenvfv
#define fglTexEnvi                     fgltexenvi
#define fglTexEnviv                    fgltexenviv
#define fglTexGend                     fgltexgend
#define fglTexGendv                    fgltexgendv
#define fglTexGenf                     fgltexgenf
#define fglTexGenfv                    fgltexgenfv
#define fglTexGeni                     fgltexgeni
#define fglTexGeniv                    fgltexgeniv
#define f9y0glTexImage1D               f9y0glteximage1d
#define f9y1glTexImage1D               f9y1glteximage1d
#define f9y4glTexImage1D               f9y4glteximage1d
#define f9y5glTexImage1D               f9y5glteximage1d
#define f9y6glTexImage1D               f9y6glteximage1d
#define f9y0glTexImage2D               f9y0glteximage2d
#define f9y1glTexImage2D               f9y1glteximage2d
#define f9y4glTexImage2D               f9y4glteximage2d
#define f9y5glTexImage2D               f9y5glteximage2d
#define f9y6glTexImage2D               f9y6glteximage2d
#define fglTexParameterf               fgltexparameterf
#define fglTexParameterfv              fgltexparameterfv
#define fglTexParameteri               fgltexparameteri
#define fglTexParameteriv              fgltexparameteriv
#define fglTranslated                  fgltranslated
#define fglTranslatef                  fgltranslatef
#define fglVertex2d                    fglvertex2d
#define fglVertex2dv                   fglvertex2dv
#define fglVertex2f                    fglvertex2f
#define fglVertex2fv                   fglvertex2fv
#define fglVertex2i                    fglvertex2i
#define fglVertex2iv                   fglvertex2iv
#define fglVertex2s                    fglvertex2s
#define f9aglVertex2s                  f9aglvertex2s
#define fglVertex2sv                   fglvertex2sv
#define f9aglVertex2sv                 f9aglvertex2sv
#define fglVertex3d                    fglvertex3d
#define fglVertex3dv                   fglvertex3dv
#define fglVertex3f                    fglvertex3f
#define fglVertex3fv                   fglvertex3fv
#define fglVertex3i                    fglvertex3i
#define fglVertex3iv                   fglvertex3iv
#define fglVertex3s                    fglvertex3s
#define f9aglVertex3s                  f9aglvertex3s
#define fglVertex3sv                   fglvertex3sv
#define f9aglVertex3sv                 f9aglvertex3sv
#define fglVertex4d                    fglvertex4d
#define fglVertex4dv                   fglvertex4dv
#define fglVertex4f                    fglvertex4f
#define fglVertex4fv                   fglvertex4fv
#define fglVertex4i                    fglvertex4i
#define fglVertex4iv                   fglvertex4iv
#define fglVertex4s                    fglvertex4s
#define f9aglVertex4s                  f9aglvertex4s
#define fglVertex4sv                   fglvertex4sv
#define f9aglVertex4sv                 f9aglvertex4sv
#define fglViewport                    fglviewport
#elif FNAME==UNDERSCORE
#define fglAccum                       fglaccum_
#define fglAlphaFunc                   fglalphafunc_
#define fglBegin                       fglbegin_
#define fglBitmap                      fglbitmap_
#define f9bglBitmap                    f9bglbitmap_
#define f9cglBitmap                    f9cglbitmap_
#define fglBlendFunc                   fglblendfunc_
#define fglCallList                    fglcalllist_
#define f9y0glCallLists                f9y0glcalllists_
#define f9y1glCallLists                f9y1glcalllists_
#define f9y4glCallLists                f9y4glcalllists_
#define f9y5glCallLists                f9y5glcalllists_
#define f9y6glCallLists                f9y6glcalllists_
#define fglClear                       fglclear_
#define fglClearAccum                  fglclearaccum_
#define fglClearColor                  fglclearcolor_
#define fglClearDepth                  fglcleardepth_
#define fglClearIndex                  fglclearindex_
#define fglClearStencil                fglclearstencil_
#define fglClipPlane                   fglclipplane_
#define fglColor3b                     fglcolor3b_
#define f9bglColor3b                   f9bglcolor3b_
#define f9cglColor3b                   f9cglcolor3b_
#define fglColor3bv                    fglcolor3bv_
#define f9bglColor3bv                  f9bglcolor3bv_
#define f9cglColor3bv                  f9cglcolor3bv_
#define fglColor3d                     fglcolor3d_
#define fglColor3dv                    fglcolor3dv_
#define fglColor3f                     fglcolor3f_
#define fglColor3fv                    fglcolor3fv_
#define fglColor3i                     fglcolor3i_
#define fglColor3iv                    fglcolor3iv_
#define fglColor3s                     fglcolor3s_
#define f9aglColor3s                   f9aglcolor3s_
#define fglColor3sv                    fglcolor3sv_
#define f9aglColor3sv                  f9aglcolor3sv_
#define fglColor3ub                    fglcolor3ub_
#define f9bglColor3ub                  f9bglcolor3ub_
#define f9cglColor3ub                  f9cglcolor3ub_
#define fglColor3ubv                   fglcolor3ubv_
#define f9bglColor3ubv                 f9bglcolor3ubv_
#define f9cglColor3ubv                 f9cglcolor3ubv_
#define fglColor3ui                    fglcolor3ui_
#define fglColor3uiv                   fglcolor3uiv_
#define fglColor3us                    fglcolor3us_
#define f9aglColor3us                  f9aglcolor3us_
#define fglColor3usv                   fglcolor3usv_
#define f9aglColor3usv                 f9aglcolor3usv_
#define fglColor4b                     fglcolor4b_
#define f9bglColor4b                   f9bglcolor4b_
#define f9cglColor4b                   f9cglcolor4b_
#define fglColor4bv                    fglcolor4bv_
#define f9bglColor4bv                  f9bglcolor4bv_
#define f9cglColor4bv                  f9cglcolor4bv_
#define fglColor4d                     fglcolor4d_
#define fglColor4dv                    fglcolor4dv_
#define fglColor4f                     fglcolor4f_
#define fglColor4fv                    fglcolor4fv_
#define fglColor4i                     fglcolor4i_
#define fglColor4iv                    fglcolor4iv_
#define fglColor4s                     fglcolor4s_
#define f9aglColor4s                   f9aglcolor4s_
#define fglColor4sv                    fglcolor4sv_
#define f9aglColor4sv                  f9aglcolor4sv_
#define fglColor4ub                    fglcolor4ub_
#define f9bglColor4ub                  f9bglcolor4ub_
#define f9cglColor4ub                  f9cglcolor4ub_
#define fglColor4ubv                   fglcolor4ubv_
#define f9bglColor4ubv                 f9bglcolor4ubv_
#define f9cglColor4ubv                 f9cglcolor4ubv_
#define fglColor4ui                    fglcolor4ui_
#define fglColor4uiv                   fglcolor4uiv_
#define fglColor4us                    fglcolor4us_
#define f9aglColor4us                  f9aglcolor4us_
#define fglColor4usv                   fglcolor4usv_
#define f9aglColor4usv                 f9aglcolor4usv_
#define fglColorMask                   fglcolormask_
#define f9eglColorMask                 f9eglcolormask_
#define fglColorMaterial               fglcolormaterial_
#define fglCopyPixels                  fglcopypixels_
#define fglCullFace                    fglcullface_
#define fglDeleteLists                 fgldeletelists_
#define fglDepthFunc                   fgldepthfunc_
#define fglDepthMask                   fgldepthmask_
#define f9eglDepthMask                 f9egldepthmask_
#define fglDepthRange                  fgldepthrange_
#define fglDisable                     fgldisable_
#define fglDrawBuffer                  fgldrawbuffer_
#define f9y0glDrawPixels               f9y0gldrawpixels_
#define f9y1glDrawPixels               f9y1gldrawpixels_
#define f9y4glDrawPixels               f9y4gldrawpixels_
#define f9y5glDrawPixels               f9y5gldrawpixels_
#define f9y6glDrawPixels               f9y6gldrawpixels_
#define fglEdgeFlag                    fgledgeflag_
#define f9eglEdgeFlag                  f9egledgeflag_
#define fglEdgeFlagv                   fgledgeflagv_
#define f9eglEdgeFlagv                 f9egledgeflagv_
#define fglEnable                      fglenable_
#define fglEnd                         fglend_
#define fglEndList                     fglendlist_
#define fglEvalCoord1d                 fglevalcoord1d_
#define fglEvalCoord1dv                fglevalcoord1dv_
#define fglEvalCoord1f                 fglevalcoord1f_
#define fglEvalCoord1fv                fglevalcoord1fv_
#define fglEvalCoord2d                 fglevalcoord2d_
#define fglEvalCoord2dv                fglevalcoord2dv_
#define fglEvalCoord2f                 fglevalcoord2f_
#define fglEvalCoord2fv                fglevalcoord2fv_
#define fglEvalMesh1                   fglevalmesh1_
#define fglEvalMesh2                   fglevalmesh2_
#define fglEvalPoint1                  fglevalpoint1_
#define fglEvalPoint2                  fglevalpoint2_
#define fglFeedbackBuffer              fglfeedbackbuffer_
#define fglFinish                      fglfinish_
#define fglFlush                       fglflush_
#define fglFogf                        fglfogf_
#define fglFogfv                       fglfogfv_
#define fglFogi                        fglfogi_
#define fglFogiv                       fglfogiv_
#define fglFrontFace                   fglfrontface_
#define fglFrustum                     fglfrustum_
#define fglGenLists                    fglgenlists_
#define fglGetBooleanv                 fglgetbooleanv_
#define f9eglGetBooleanv               f9eglgetbooleanv_
#define fglGetClipPlane                fglgetclipplane_
#define fglGetDoublev                  fglgetdoublev_
#define fglGetError                    fglgeterror_
#define fglGetFloatv                   fglgetfloatv_
#define fglGetIntegerv                 fglgetintegerv_
#define fglGetLightfv                  fglgetlightfv_
#define fglGetLightiv                  fglgetlightiv_
#define fglGetMapdv                    fglgetmapdv_
#define fglGetMapfv                    fglgetmapfv_
#define fglGetMapiv                    fglgetmapiv_
#define fglGetMaterialfv               fglgetmaterialfv_
#define fglGetMaterialiv               fglgetmaterialiv_
#define fglGetPixelMapfv               fglgetpixelmapfv_
#define fglGetPixelMapuiv              fglgetpixelmapuiv_
#define fglGetPixelMapusv              fglgetpixelmapusv_
#define f9aglGetPixelMapusv            f9aglgetpixelmapusv_
#define fglGetPolygonStipple           fglgetpolygonstipple_
#define f9bglGetPolygonStipple         f9bglgetpolygonstipple_
#define f9cglGetPolygonStipple         f9cglgetpolygonstipple_
#define f9yglGetStringLen              f9yglgetstringlen_
#define f9yglGetString                 f9yglgetstring_
#define f9yglGetString1                f9yglgetstring1_
#define f9yglGetStringFree             f9yglgetstringfree_
#define fglGetTexEnvfv                 fglgettexenvfv_
#define fglGetTexEnviv                 fglgettexenviv_
#define fglGetTexGendv                 fglgettexgendv_
#define fglGetTexGenfv                 fglgettexgenfv_
#define fglGetTexGeniv                 fglgettexgeniv_
#define f9y0glGetTexImage              f9y0glgetteximage_
#define f9y1glGetTexImage              f9y1glgetteximage_
#define f9y4glGetTexImage              f9y4glgetteximage_
#define f9y5glGetTexImage              f9y5glgetteximage_
#define f9y6glGetTexImage              f9y6glgetteximage_
#define fglGetTexLevelParameterfv      fglgettexlevelparameterfv_
#define fglGetTexLevelParameteriv      fglgettexlevelparameteriv_
#define fglGetTexParameterfv           fglgettexparameterfv_
#define fglGetTexParameteriv           fglgettexparameteriv_
#define fglHint                        fglhint_
#define fglIndexMask                   fglindexmask_
#define fglIndexd                      fglindexd_
#define fglIndexdv                     fglindexdv_
#define fglIndexf                      fglindexf_
#define fglIndexfv                     fglindexfv_
#define fglIndexi                      fglindexi_
#define fglIndexiv                     fglindexiv_
#define fglIndexs                      fglindexs_
#define f9aglIndexs                    f9aglindexs_
#define fglIndexsv                     fglindexsv_
#define f9aglIndexsv                   f9aglindexsv_
#define fglInitNames                   fglinitnames_
#define fglIsEnabled                   fglisenabled_
#define f9eglIsEnabled                 f9eglisenabled_
#define fglIsList                      fglislist_
#define f9eglIsList                    f9eglislist_
#define fglLightModelf                 fgllightmodelf_
#define fglLightModelfv                fgllightmodelfv_
#define fglLightModeli                 fgllightmodeli_
#define fglLightModeliv                fgllightmodeliv_
#define fglLightf                      fgllightf_
#define fglLightfv                     fgllightfv_
#define fglLighti                      fgllighti_
#define fglLightiv                     fgllightiv_
#define fglLineStipple                 fgllinestipple_
#define f9aglLineStipple               f9agllinestipple_
#define fglLineWidth                   fgllinewidth_
#define fglListBase                    fgllistbase_
#define fglLoadIdentity                fglloadidentity_
#define fglLoadMatrixd                 fglloadmatrixd_
#define fglLoadMatrixf                 fglloadmatrixf_
#define fglLoadName                    fglloadname_
#define fglLogicOp                     fgllogicop_
#define fglMap1d                       fglmap1d_
#define fglMap1f                       fglmap1f_
#define fglMap2d                       fglmap2d_
#define fglMap2f                       fglmap2f_
#define fglMapGrid1d                   fglmapgrid1d_
#define fglMapGrid1f                   fglmapgrid1f_
#define fglMapGrid2d                   fglmapgrid2d_
#define fglMapGrid2f                   fglmapgrid2f_
#define fglMaterialf                   fglmaterialf_
#define fglMaterialfv                  fglmaterialfv_
#define fglMateriali                   fglmateriali_
#define fglMaterialiv                  fglmaterialiv_
#define fglMatrixMode                  fglmatrixmode_
#define fglMultMatrixd                 fglmultmatrixd_
#define fglMultMatrixf                 fglmultmatrixf_
#define fglNewList                     fglnewlist_
#define fglNormal3b                    fglnormal3b_
#define f9bglNormal3b                  f9bglnormal3b_
#define f9cglNormal3b                  f9cglnormal3b_
#define fglNormal3bv                   fglnormal3bv_
#define f9bglNormal3bv                 f9bglnormal3bv_
#define f9cglNormal3bv                 f9cglnormal3bv_
#define fglNormal3d                    fglnormal3d_
#define fglNormal3dv                   fglnormal3dv_
#define fglNormal3f                    fglnormal3f_
#define fglNormal3fv                   fglnormal3fv_
#define fglNormal3i                    fglnormal3i_
#define fglNormal3iv                   fglnormal3iv_
#define fglNormal3s                    fglnormal3s_
#define f9aglNormal3s                  f9aglnormal3s_
#define fglNormal3sv                   fglnormal3sv_
#define f9aglNormal3sv                 f9aglnormal3sv_
#define fglOrtho                       fglortho_
#define fglPassThrough                 fglpassthrough_
#define fglPixelMapfv                  fglpixelmapfv_
#define fglPixelMapuiv                 fglpixelmapuiv_
#define fglPixelMapusv                 fglpixelmapusv_
#define f9aglPixelMapusv               f9aglpixelmapusv_
#define fglPixelStoref                 fglpixelstoref_
#define fglPixelStorei                 fglpixelstorei_
#define fglPixelTransferf              fglpixeltransferf_
#define fglPixelTransferi              fglpixeltransferi_
#define fglPixelZoom                   fglpixelzoom_
#define fglPointSize                   fglpointsize_
#define fglPolygonMode                 fglpolygonmode_
#define fglPolygonStipple              fglpolygonstipple_
#define f9bglPolygonStipple            f9bglpolygonstipple_
#define f9cglPolygonStipple            f9cglpolygonstipple_
#define fglPopAttrib                   fglpopattrib_
#define fglPopMatrix                   fglpopmatrix_
#define fglPopName                     fglpopname_
#define fglPushAttrib                  fglpushattrib_
#define fglPushMatrix                  fglpushmatrix_
#define fglPushName                    fglpushname_
#define fglRasterPos2d                 fglrasterpos2d_
#define fglRasterPos2dv                fglrasterpos2dv_
#define fglRasterPos2f                 fglrasterpos2f_
#define fglRasterPos2fv                fglrasterpos2fv_
#define fglRasterPos2i                 fglrasterpos2i_
#define fglRasterPos2iv                fglrasterpos2iv_
#define fglRasterPos2s                 fglrasterpos2s_
#define f9aglRasterPos2s               f9aglrasterpos2s_
#define fglRasterPos2sv                fglrasterpos2sv_
#define f9aglRasterPos2sv              f9aglrasterpos2sv_
#define fglRasterPos3d                 fglrasterpos3d_
#define fglRasterPos3dv                fglrasterpos3dv_
#define fglRasterPos3f                 fglrasterpos3f_
#define fglRasterPos3fv                fglrasterpos3fv_
#define fglRasterPos3i                 fglrasterpos3i_
#define fglRasterPos3iv                fglrasterpos3iv_
#define fglRasterPos3s                 fglrasterpos3s_
#define f9aglRasterPos3s               f9aglrasterpos3s_
#define fglRasterPos3sv                fglrasterpos3sv_
#define f9aglRasterPos3sv              f9aglrasterpos3sv_
#define fglRasterPos4d                 fglrasterpos4d_
#define fglRasterPos4dv                fglrasterpos4dv_
#define fglRasterPos4f                 fglrasterpos4f_
#define fglRasterPos4fv                fglrasterpos4fv_
#define fglRasterPos4i                 fglrasterpos4i_
#define fglRasterPos4iv                fglrasterpos4iv_
#define fglRasterPos4s                 fglrasterpos4s_
#define f9aglRasterPos4s               f9aglrasterpos4s_
#define fglRasterPos4sv                fglrasterpos4sv_
#define f9aglRasterPos4sv              f9aglrasterpos4sv_
#define fglReadBuffer                  fglreadbuffer_
#define f9y0glReadPixels               f9y0glreadpixels_
#define f9y1glReadPixels               f9y1glreadpixels_
#define f9y4glReadPixels               f9y4glreadpixels_
#define f9y5glReadPixels               f9y5glreadpixels_
#define f9y6glReadPixels               f9y6glreadpixels_
#define fglRectd                       fglrectd_
#define fglRectdv                      fglrectdv_
#define fglRectf                       fglrectf_
#define fglRectfv                      fglrectfv_
#define fglRecti                       fglrecti_
#define fglRectiv                      fglrectiv_
#define fglRects                       fglrects_
#define f9aglRects                     f9aglrects_
#define fglRectsv                      fglrectsv_
#define f9aglRectsv                    f9aglrectsv_
#define fglRenderMode                  fglrendermode_
#define fglRotated                     fglrotated_
#define fglRotatef                     fglrotatef_
#define fglScaled                      fglscaled_
#define fglScalef                      fglscalef_
#define fglScissor                     fglscissor_
#define fglSelectBuffer                fglselectbuffer_
#define fglShadeModel                  fglshademodel_
#define fglStencilFunc                 fglstencilfunc_
#define fglStencilMask                 fglstencilmask_
#define fglStencilOp                   fglstencilop_
#define fglTexCoord1d                  fgltexcoord1d_
#define fglTexCoord1dv                 fgltexcoord1dv_
#define fglTexCoord1f                  fgltexcoord1f_
#define fglTexCoord1fv                 fgltexcoord1fv_
#define fglTexCoord1i                  fgltexcoord1i_
#define fglTexCoord1iv                 fgltexcoord1iv_
#define fglTexCoord1s                  fgltexcoord1s_
#define f9aglTexCoord1s                f9agltexcoord1s_
#define fglTexCoord1sv                 fgltexcoord1sv_
#define f9aglTexCoord1sv               f9agltexcoord1sv_
#define fglTexCoord2d                  fgltexcoord2d_
#define fglTexCoord2dv                 fgltexcoord2dv_
#define fglTexCoord2f                  fgltexcoord2f_
#define fglTexCoord2fv                 fgltexcoord2fv_
#define fglTexCoord2i                  fgltexcoord2i_
#define fglTexCoord2iv                 fgltexcoord2iv_
#define fglTexCoord2s                  fgltexcoord2s_
#define f9aglTexCoord2s                f9agltexcoord2s_
#define fglTexCoord2sv                 fgltexcoord2sv_
#define f9aglTexCoord2sv               f9agltexcoord2sv_
#define fglTexCoord3d                  fgltexcoord3d_
#define fglTexCoord3dv                 fgltexcoord3dv_
#define fglTexCoord3f                  fgltexcoord3f_
#define fglTexCoord3fv                 fgltexcoord3fv_
#define fglTexCoord3i                  fgltexcoord3i_
#define fglTexCoord3iv                 fgltexcoord3iv_
#define fglTexCoord3s                  fgltexcoord3s_
#define f9aglTexCoord3s                f9agltexcoord3s_
#define fglTexCoord3sv                 fgltexcoord3sv_
#define f9aglTexCoord3sv               f9agltexcoord3sv_
#define fglTexCoord4d                  fgltexcoord4d_
#define fglTexCoord4dv                 fgltexcoord4dv_
#define fglTexCoord4f                  fgltexcoord4f_
#define fglTexCoord4fv                 fgltexcoord4fv_
#define fglTexCoord4i                  fgltexcoord4i_
#define fglTexCoord4iv                 fgltexcoord4iv_
#define fglTexCoord4s                  fgltexcoord4s_
#define f9aglTexCoord4s                f9agltexcoord4s_
#define fglTexCoord4sv                 fgltexcoord4sv_
#define f9aglTexCoord4sv               f9agltexcoord4sv_
#define fglTexEnvf                     fgltexenvf_
#define fglTexEnvfv                    fgltexenvfv_
#define fglTexEnvi                     fgltexenvi_
#define fglTexEnviv                    fgltexenviv_
#define fglTexGend                     fgltexgend_
#define fglTexGendv                    fgltexgendv_
#define fglTexGenf                     fgltexgenf_
#define fglTexGenfv                    fgltexgenfv_
#define fglTexGeni                     fgltexgeni_
#define fglTexGeniv                    fgltexgeniv_
#define f9y0glTexImage1D               f9y0glteximage1d_
#define f9y1glTexImage1D               f9y1glteximage1d_
#define f9y4glTexImage1D               f9y4glteximage1d_
#define f9y5glTexImage1D               f9y5glteximage1d_
#define f9y6glTexImage1D               f9y6glteximage1d_
#define f9y0glTexImage2D               f9y0glteximage2d_
#define f9y1glTexImage2D               f9y1glteximage2d_
#define f9y4glTexImage2D               f9y4glteximage2d_
#define f9y5glTexImage2D               f9y5glteximage2d_
#define f9y6glTexImage2D               f9y6glteximage2d_
#define fglTexParameterf               fgltexparameterf_
#define fglTexParameterfv              fgltexparameterfv_
#define fglTexParameteri               fgltexparameteri_
#define fglTexParameteriv              fgltexparameteriv_
#define fglTranslated                  fgltranslated_
#define fglTranslatef                  fgltranslatef_
#define fglVertex2d                    fglvertex2d_
#define fglVertex2dv                   fglvertex2dv_
#define fglVertex2f                    fglvertex2f_
#define fglVertex2fv                   fglvertex2fv_
#define fglVertex2i                    fglvertex2i_
#define fglVertex2iv                   fglvertex2iv_
#define fglVertex2s                    fglvertex2s_
#define f9aglVertex2s                  f9aglvertex2s_
#define fglVertex2sv                   fglvertex2sv_
#define f9aglVertex2sv                 f9aglvertex2sv_
#define fglVertex3d                    fglvertex3d_
#define fglVertex3dv                   fglvertex3dv_
#define fglVertex3f                    fglvertex3f_
#define fglVertex3fv                   fglvertex3fv_
#define fglVertex3i                    fglvertex3i_
#define fglVertex3iv                   fglvertex3iv_
#define fglVertex3s                    fglvertex3s_
#define f9aglVertex3s                  f9aglvertex3s_
#define fglVertex3sv                   fglvertex3sv_
#define f9aglVertex3sv                 f9aglvertex3sv_
#define fglVertex4d                    fglvertex4d_
#define fglVertex4dv                   fglvertex4dv_
#define fglVertex4f                    fglvertex4f_
#define fglVertex4fv                   fglvertex4fv_
#define fglVertex4i                    fglvertex4i_
#define fglVertex4iv                   fglvertex4iv_
#define fglVertex4s                    fglvertex4s_
#define f9aglVertex4s                  f9aglvertex4s_
#define fglVertex4sv                   fglvertex4sv_
#define f9aglVertex4sv                 f9aglvertex4sv_
#define fglViewport                    fglviewport_
#elif FNAME==UPPERCASE
#define fglAccum                       FGLACCUM
#define fglAlphaFunc                   FGLALPHAFUNC
#define fglBegin                       FGLBEGIN
#define fglBitmap                      FGLBITMAP
#define f9bglBitmap                    F9BGLBITMAP
#define f9cglBitmap                    F9CGLBITMAP
#define fglBlendFunc                   FGLBLENDFUNC
#define fglCallList                    FGLCALLLIST
#define f9y0glCallLists                F9Y0GLCALLLISTS
#define f9y1glCallLists                F9Y1GLCALLLISTS
#define f9y4glCallLists                F9Y4GLCALLLISTS
#define f9y5glCallLists                F9Y5GLCALLLISTS
#define f9y6glCallLists                F9Y6GLCALLLISTS
#define fglClear                       FGLCLEAR
#define fglClearAccum                  FGLCLEARACCUM
#define fglClearColor                  FGLCLEARCOLOR
#define fglClearDepth                  FGLCLEARDEPTH
#define fglClearIndex                  FGLCLEARINDEX
#define fglClearStencil                FGLCLEARSTENCIL
#define fglClipPlane                   FGLCLIPPLANE
#define fglColor3b                     FGLCOLOR3B
#define f9bglColor3b                   F9BGLCOLOR3B
#define f9cglColor3b                   F9CGLCOLOR3B
#define fglColor3bv                    FGLCOLOR3BV
#define f9bglColor3bv                  F9BGLCOLOR3BV
#define f9cglColor3bv                  F9CGLCOLOR3BV
#define fglColor3d                     FGLCOLOR3D
#define fglColor3dv                    FGLCOLOR3DV
#define fglColor3f                     FGLCOLOR3F
#define fglColor3fv                    FGLCOLOR3FV
#define fglColor3i                     FGLCOLOR3I
#define fglColor3iv                    FGLCOLOR3IV
#define fglColor3s                     FGLCOLOR3S
#define f9aglColor3s                   F9AGLCOLOR3S
#define fglColor3sv                    FGLCOLOR3SV
#define f9aglColor3sv                  F9AGLCOLOR3SV
#define fglColor3ub                    FGLCOLOR3UB
#define f9bglColor3ub                  F9BGLCOLOR3UB
#define f9cglColor3ub                  F9CGLCOLOR3UB
#define fglColor3ubv                   FGLCOLOR3UBV
#define f9bglColor3ubv                 F9BGLCOLOR3UBV
#define f9cglColor3ubv                 F9CGLCOLOR3UBV
#define fglColor3ui                    FGLCOLOR3UI
#define fglColor3uiv                   FGLCOLOR3UIV
#define fglColor3us                    FGLCOLOR3US
#define f9aglColor3us                  F9AGLCOLOR3US
#define fglColor3usv                   FGLCOLOR3USV
#define f9aglColor3usv                 F9AGLCOLOR3USV
#define fglColor4b                     FGLCOLOR4B
#define f9bglColor4b                   F9BGLCOLOR4B
#define f9cglColor4b                   F9CGLCOLOR4B
#define fglColor4bv                    FGLCOLOR4BV
#define f9bglColor4bv                  F9BGLCOLOR4BV
#define f9cglColor4bv                  F9CGLCOLOR4BV
#define fglColor4d                     FGLCOLOR4D
#define fglColor4dv                    FGLCOLOR4DV
#define fglColor4f                     FGLCOLOR4F
#define fglColor4fv                    FGLCOLOR4FV
#define fglColor4i                     FGLCOLOR4I
#define fglColor4iv                    FGLCOLOR4IV
#define fglColor4s                     FGLCOLOR4S
#define f9aglColor4s                   F9AGLCOLOR4S
#define fglColor4sv                    FGLCOLOR4SV
#define f9aglColor4sv                  F9AGLCOLOR4SV
#define fglColor4ub                    FGLCOLOR4UB
#define f9bglColor4ub                  F9BGLCOLOR4UB
#define f9cglColor4ub                  F9CGLCOLOR4UB
#define fglColor4ubv                   FGLCOLOR4UBV
#define f9bglColor4ubv                 F9BGLCOLOR4UBV
#define f9cglColor4ubv                 F9CGLCOLOR4UBV
#define fglColor4ui                    FGLCOLOR4UI
#define fglColor4uiv                   FGLCOLOR4UIV
#define fglColor4us                    FGLCOLOR4US
#define f9aglColor4us                  F9AGLCOLOR4US
#define fglColor4usv                   FGLCOLOR4USV
#define f9aglColor4usv                 F9AGLCOLOR4USV
#define fglColorMask                   FGLCOLORMASK
#define f9eglColorMask                 F9EGLCOLORMASK
#define fglColorMaterial               FGLCOLORMATERIAL
#define fglCopyPixels                  FGLCOPYPIXELS
#define fglCullFace                    FGLCULLFACE
#define fglDeleteLists                 FGLDELETELISTS
#define fglDepthFunc                   FGLDEPTHFUNC
#define fglDepthMask                   FGLDEPTHMASK
#define f9eglDepthMask                 F9EGLDEPTHMASK
#define fglDepthRange                  FGLDEPTHRANGE
#define fglDisable                     FGLDISABLE
#define fglDrawBuffer                  FGLDRAWBUFFER
#define f9y0glDrawPixels               F9Y0GLDRAWPIXELS
#define f9y1glDrawPixels               F9Y1GLDRAWPIXELS
#define f9y4glDrawPixels               F9Y4GLDRAWPIXELS
#define f9y5glDrawPixels               F9Y5GLDRAWPIXELS
#define f9y6glDrawPixels               F9Y6GLDRAWPIXELS
#define fglEdgeFlag                    FGLEDGEFLAG
#define f9eglEdgeFlag                  F9EGLEDGEFLAG
#define fglEdgeFlagv                   FGLEDGEFLAGV
#define f9eglEdgeFlagv                 F9EGLEDGEFLAGV
#define fglEnable                      FGLENABLE
#define fglEnd                         FGLEND
#define fglEndList                     FGLENDLIST
#define fglEvalCoord1d                 FGLEVALCOORD1D
#define fglEvalCoord1dv                FGLEVALCOORD1DV
#define fglEvalCoord1f                 FGLEVALCOORD1F
#define fglEvalCoord1fv                FGLEVALCOORD1FV
#define fglEvalCoord2d                 FGLEVALCOORD2D
#define fglEvalCoord2dv                FGLEVALCOORD2DV
#define fglEvalCoord2f                 FGLEVALCOORD2F
#define fglEvalCoord2fv                FGLEVALCOORD2FV
#define fglEvalMesh1                   FGLEVALMESH1
#define fglEvalMesh2                   FGLEVALMESH2
#define fglEvalPoint1                  FGLEVALPOINT1
#define fglEvalPoint2                  FGLEVALPOINT2
#define fglFeedbackBuffer              FGLFEEDBACKBUFFER
#define fglFinish                      FGLFINISH
#define fglFlush                       FGLFLUSH
#define fglFogf                        FGLFOGF
#define fglFogfv                       FGLFOGFV
#define fglFogi                        FGLFOGI
#define fglFogiv                       FGLFOGIV
#define fglFrontFace                   FGLFRONTFACE
#define fglFrustum                     FGLFRUSTUM
#define fglGenLists                    FGLGENLISTS
#define fglGetBooleanv                 FGLGETBOOLEANV
#define f9eglGetBooleanv               F9EGLGETBOOLEANV
#define fglGetClipPlane                FGLGETCLIPPLANE
#define fglGetDoublev                  FGLGETDOUBLEV
#define fglGetError                    FGLGETERROR
#define fglGetFloatv                   FGLGETFLOATV
#define fglGetIntegerv                 FGLGETINTEGERV
#define fglGetLightfv                  FGLGETLIGHTFV
#define fglGetLightiv                  FGLGETLIGHTIV
#define fglGetMapdv                    FGLGETMAPDV
#define fglGetMapfv                    FGLGETMAPFV
#define fglGetMapiv                    FGLGETMAPIV
#define fglGetMaterialfv               FGLGETMATERIALFV
#define fglGetMaterialiv               FGLGETMATERIALIV
#define fglGetPixelMapfv               FGLGETPIXELMAPFV
#define fglGetPixelMapuiv              FGLGETPIXELMAPUIV
#define fglGetPixelMapusv              FGLGETPIXELMAPUSV
#define f9aglGetPixelMapusv            F9AGLGETPIXELMAPUSV
#define fglGetPolygonStipple           FGLGETPOLYGONSTIPPLE
#define f9bglGetPolygonStipple         F9BGLGETPOLYGONSTIPPLE
#define f9cglGetPolygonStipple         F9CGLGETPOLYGONSTIPPLE
#define f9yglGetStringLen              F9YGLGETSTRINGLEN
#define f9yglGetString                 F9YGLGETSTRING
#define f9yglGetString1                F9YGLGETSTRING1
#define f9yglGetStringFree             F9YGLGETSTRINGFREE
#define fglGetTexEnvfv                 FGLGETTEXENVFV
#define fglGetTexEnviv                 FGLGETTEXENVIV
#define fglGetTexGendv                 FGLGETTEXGENDV
#define fglGetTexGenfv                 FGLGETTEXGENFV
#define fglGetTexGeniv                 FGLGETTEXGENIV
#define f9y0glGetTexImage              F9Y0GLGETTEXIMAGE
#define f9y1glGetTexImage              F9Y1GLGETTEXIMAGE
#define f9y4glGetTexImage              F9Y4GLGETTEXIMAGE
#define f9y5glGetTexImage              F9Y5GLGETTEXIMAGE
#define f9y6glGetTexImage              F9Y6GLGETTEXIMAGE
#define fglGetTexLevelParameterfv      FGLGETTEXLEVELPARAMETERFV
#define fglGetTexLevelParameteriv      FGLGETTEXLEVELPARAMETERIV
#define fglGetTexParameterfv           FGLGETTEXPARAMETERFV
#define fglGetTexParameteriv           FGLGETTEXPARAMETERIV
#define fglHint                        FGLHINT
#define fglIndexMask                   FGLINDEXMASK
#define fglIndexd                      FGLINDEXD
#define fglIndexdv                     FGLINDEXDV
#define fglIndexf                      FGLINDEXF
#define fglIndexfv                     FGLINDEXFV
#define fglIndexi                      FGLINDEXI
#define fglIndexiv                     FGLINDEXIV
#define fglIndexs                      FGLINDEXS
#define f9aglIndexs                    F9AGLINDEXS
#define fglIndexsv                     FGLINDEXSV
#define f9aglIndexsv                   F9AGLINDEXSV
#define fglInitNames                   FGLINITNAMES
#define fglIsEnabled                   FGLISENABLED
#define f9eglIsEnabled                 F9EGLISENABLED
#define fglIsList                      FGLISLIST
#define f9eglIsList                    F9EGLISLIST
#define fglLightModelf                 FGLLIGHTMODELF
#define fglLightModelfv                FGLLIGHTMODELFV
#define fglLightModeli                 FGLLIGHTMODELI
#define fglLightModeliv                FGLLIGHTMODELIV
#define fglLightf                      FGLLIGHTF
#define fglLightfv                     FGLLIGHTFV
#define fglLighti                      FGLLIGHTI
#define fglLightiv                     FGLLIGHTIV
#define fglLineStipple                 FGLLINESTIPPLE
#define f9aglLineStipple               F9AGLLINESTIPPLE
#define fglLineWidth                   FGLLINEWIDTH
#define fglListBase                    FGLLISTBASE
#define fglLoadIdentity                FGLLOADIDENTITY
#define fglLoadMatrixd                 FGLLOADMATRIXD
#define fglLoadMatrixf                 FGLLOADMATRIXF
#define fglLoadName                    FGLLOADNAME
#define fglLogicOp                     FGLLOGICOP
#define fglMap1d                       FGLMAP1D
#define fglMap1f                       FGLMAP1F
#define fglMap2d                       FGLMAP2D
#define fglMap2f                       FGLMAP2F
#define fglMapGrid1d                   FGLMAPGRID1D
#define fglMapGrid1f                   FGLMAPGRID1F
#define fglMapGrid2d                   FGLMAPGRID2D
#define fglMapGrid2f                   FGLMAPGRID2F
#define fglMaterialf                   FGLMATERIALF
#define fglMaterialfv                  FGLMATERIALFV
#define fglMateriali                   FGLMATERIALI
#define fglMaterialiv                  FGLMATERIALIV
#define fglMatrixMode                  FGLMATRIXMODE
#define fglMultMatrixd                 FGLMULTMATRIXD
#define fglMultMatrixf                 FGLMULTMATRIXF
#define fglNewList                     FGLNEWLIST
#define fglNormal3b                    FGLNORMAL3B
#define f9bglNormal3b                  F9BGLNORMAL3B
#define f9cglNormal3b                  F9CGLNORMAL3B
#define fglNormal3bv                   FGLNORMAL3BV
#define f9bglNormal3bv                 F9BGLNORMAL3BV
#define f9cglNormal3bv                 F9CGLNORMAL3BV
#define fglNormal3d                    FGLNORMAL3D
#define fglNormal3dv                   FGLNORMAL3DV
#define fglNormal3f                    FGLNORMAL3F
#define fglNormal3fv                   FGLNORMAL3FV
#define fglNormal3i                    FGLNORMAL3I
#define fglNormal3iv                   FGLNORMAL3IV
#define fglNormal3s                    FGLNORMAL3S
#define f9aglNormal3s                  F9AGLNORMAL3S
#define fglNormal3sv                   FGLNORMAL3SV
#define f9aglNormal3sv                 F9AGLNORMAL3SV
#define fglOrtho                       FGLORTHO
#define fglPassThrough                 FGLPASSTHROUGH
#define fglPixelMapfv                  FGLPIXELMAPFV
#define fglPixelMapuiv                 FGLPIXELMAPUIV
#define fglPixelMapusv                 FGLPIXELMAPUSV
#define f9aglPixelMapusv               F9AGLPIXELMAPUSV
#define fglPixelStoref                 FGLPIXELSTOREF
#define fglPixelStorei                 FGLPIXELSTOREI
#define fglPixelTransferf              FGLPIXELTRANSFERF
#define fglPixelTransferi              FGLPIXELTRANSFERI
#define fglPixelZoom                   FGLPIXELZOOM
#define fglPointSize                   FGLPOINTSIZE
#define fglPolygonMode                 FGLPOLYGONMODE
#define fglPolygonStipple              FGLPOLYGONSTIPPLE
#define f9bglPolygonStipple            F9BGLPOLYGONSTIPPLE
#define f9cglPolygonStipple            F9CGLPOLYGONSTIPPLE
#define fglPopAttrib                   FGLPOPATTRIB
#define fglPopMatrix                   FGLPOPMATRIX
#define fglPopName                     FGLPOPNAME
#define fglPushAttrib                  FGLPUSHATTRIB
#define fglPushMatrix                  FGLPUSHMATRIX
#define fglPushName                    FGLPUSHNAME
#define fglRasterPos2d                 FGLRASTERPOS2D
#define fglRasterPos2dv                FGLRASTERPOS2DV
#define fglRasterPos2f                 FGLRASTERPOS2F
#define fglRasterPos2fv                FGLRASTERPOS2FV
#define fglRasterPos2i                 FGLRASTERPOS2I
#define fglRasterPos2iv                FGLRASTERPOS2IV
#define fglRasterPos2s                 FGLRASTERPOS2S
#define f9aglRasterPos2s               F9AGLRASTERPOS2S
#define fglRasterPos2sv                FGLRASTERPOS2SV
#define f9aglRasterPos2sv              F9AGLRASTERPOS2SV
#define fglRasterPos3d                 FGLRASTERPOS3D
#define fglRasterPos3dv                FGLRASTERPOS3DV
#define fglRasterPos3f                 FGLRASTERPOS3F
#define fglRasterPos3fv                FGLRASTERPOS3FV
#define fglRasterPos3i                 FGLRASTERPOS3I
#define fglRasterPos3iv                FGLRASTERPOS3IV
#define fglRasterPos3s                 FGLRASTERPOS3S
#define f9aglRasterPos3s               F9AGLRASTERPOS3S
#define fglRasterPos3sv                FGLRASTERPOS3SV
#define f9aglRasterPos3sv              F9AGLRASTERPOS3SV
#define fglRasterPos4d                 FGLRASTERPOS4D
#define fglRasterPos4dv                FGLRASTERPOS4DV
#define fglRasterPos4f                 FGLRASTERPOS4F
#define fglRasterPos4fv                FGLRASTERPOS4FV
#define fglRasterPos4i                 FGLRASTERPOS4I
#define fglRasterPos4iv                FGLRASTERPOS4IV
#define fglRasterPos4s                 FGLRASTERPOS4S
#define f9aglRasterPos4s               F9AGLRASTERPOS4S
#define fglRasterPos4sv                FGLRASTERPOS4SV
#define f9aglRasterPos4sv              F9AGLRASTERPOS4SV
#define fglReadBuffer                  FGLREADBUFFER
#define f9y0glReadPixels               F9Y0GLREADPIXELS
#define f9y1glReadPixels               F9Y1GLREADPIXELS
#define f9y4glReadPixels               F9Y4GLREADPIXELS
#define f9y5glReadPixels               F9Y5GLREADPIXELS
#define f9y6glReadPixels               F9Y6GLREADPIXELS
#define fglRectd                       FGLRECTD
#define fglRectdv                      FGLRECTDV
#define fglRectf                       FGLRECTF
#define fglRectfv                      FGLRECTFV
#define fglRecti                       FGLRECTI
#define fglRectiv                      FGLRECTIV
#define fglRects                       FGLRECTS
#define f9aglRects                     F9AGLRECTS
#define fglRectsv                      FGLRECTSV
#define f9aglRectsv                    F9AGLRECTSV
#define fglRenderMode                  FGLRENDERMODE
#define fglRotated                     FGLROTATED
#define fglRotatef                     FGLROTATEF
#define fglScaled                      FGLSCALED
#define fglScalef                      FGLSCALEF
#define fglScissor                     FGLSCISSOR
#define fglSelectBuffer                FGLSELECTBUFFER
#define fglShadeModel                  FGLSHADEMODEL
#define fglStencilFunc                 FGLSTENCILFUNC
#define fglStencilMask                 FGLSTENCILMASK
#define fglStencilOp                   FGLSTENCILOP
#define fglTexCoord1d                  FGLTEXCOORD1D
#define fglTexCoord1dv                 FGLTEXCOORD1DV
#define fglTexCoord1f                  FGLTEXCOORD1F
#define fglTexCoord1fv                 FGLTEXCOORD1FV
#define fglTexCoord1i                  FGLTEXCOORD1I
#define fglTexCoord1iv                 FGLTEXCOORD1IV
#define fglTexCoord1s                  FGLTEXCOORD1S
#define f9aglTexCoord1s                F9AGLTEXCOORD1S
#define fglTexCoord1sv                 FGLTEXCOORD1SV
#define f9aglTexCoord1sv               F9AGLTEXCOORD1SV
#define fglTexCoord2d                  FGLTEXCOORD2D
#define fglTexCoord2dv                 FGLTEXCOORD2DV
#define fglTexCoord2f                  FGLTEXCOORD2F
#define fglTexCoord2fv                 FGLTEXCOORD2FV
#define fglTexCoord2i                  FGLTEXCOORD2I
#define fglTexCoord2iv                 FGLTEXCOORD2IV
#define fglTexCoord2s                  FGLTEXCOORD2S
#define f9aglTexCoord2s                F9AGLTEXCOORD2S
#define fglTexCoord2sv                 FGLTEXCOORD2SV
#define f9aglTexCoord2sv               F9AGLTEXCOORD2SV
#define fglTexCoord3d                  FGLTEXCOORD3D
#define fglTexCoord3dv                 FGLTEXCOORD3DV
#define fglTexCoord3f                  FGLTEXCOORD3F
#define fglTexCoord3fv                 FGLTEXCOORD3FV
#define fglTexCoord3i                  FGLTEXCOORD3I
#define fglTexCoord3iv                 FGLTEXCOORD3IV
#define fglTexCoord3s                  FGLTEXCOORD3S
#define f9aglTexCoord3s                F9AGLTEXCOORD3S
#define fglTexCoord3sv                 FGLTEXCOORD3SV
#define f9aglTexCoord3sv               F9AGLTEXCOORD3SV
#define fglTexCoord4d                  FGLTEXCOORD4D
#define fglTexCoord4dv                 FGLTEXCOORD4DV
#define fglTexCoord4f                  FGLTEXCOORD4F
#define fglTexCoord4fv                 FGLTEXCOORD4FV
#define fglTexCoord4i                  FGLTEXCOORD4I
#define fglTexCoord4iv                 FGLTEXCOORD4IV
#define fglTexCoord4s                  FGLTEXCOORD4S
#define f9aglTexCoord4s                F9AGLTEXCOORD4S
#define fglTexCoord4sv                 FGLTEXCOORD4SV
#define f9aglTexCoord4sv               F9AGLTEXCOORD4SV
#define fglTexEnvf                     FGLTEXENVF
#define fglTexEnvfv                    FGLTEXENVFV
#define fglTexEnvi                     FGLTEXENVI
#define fglTexEnviv                    FGLTEXENVIV
#define fglTexGend                     FGLTEXGEND
#define fglTexGendv                    FGLTEXGENDV
#define fglTexGenf                     FGLTEXGENF
#define fglTexGenfv                    FGLTEXGENFV
#define fglTexGeni                     FGLTEXGENI
#define fglTexGeniv                    FGLTEXGENIV
#define f9y0glTexImage1D               F9Y0GLTEXIMAGE1D
#define f9y1glTexImage1D               F9Y1GLTEXIMAGE1D
#define f9y4glTexImage1D               F9Y4GLTEXIMAGE1D
#define f9y5glTexImage1D               F9Y5GLTEXIMAGE1D
#define f9y6glTexImage1D               F9Y6GLTEXIMAGE1D
#define f9y0glTexImage2D               F9Y0GLTEXIMAGE2D
#define f9y1glTexImage2D               F9Y1GLTEXIMAGE2D
#define f9y4glTexImage2D               F9Y4GLTEXIMAGE2D
#define f9y5glTexImage2D               F9Y5GLTEXIMAGE2D
#define f9y6glTexImage2D               F9Y6GLTEXIMAGE2D
#define fglTexParameterf               FGLTEXPARAMETERF
#define fglTexParameterfv              FGLTEXPARAMETERFV
#define fglTexParameteri               FGLTEXPARAMETERI
#define fglTexParameteriv              FGLTEXPARAMETERIV
#define fglTranslated                  FGLTRANSLATED
#define fglTranslatef                  FGLTRANSLATEF
#define fglVertex2d                    FGLVERTEX2D
#define fglVertex2dv                   FGLVERTEX2DV
#define fglVertex2f                    FGLVERTEX2F
#define fglVertex2fv                   FGLVERTEX2FV
#define fglVertex2i                    FGLVERTEX2I
#define fglVertex2iv                   FGLVERTEX2IV
#define fglVertex2s                    FGLVERTEX2S
#define f9aglVertex2s                  F9AGLVERTEX2S
#define fglVertex2sv                   FGLVERTEX2SV
#define f9aglVertex2sv                 F9AGLVERTEX2SV
#define fglVertex3d                    FGLVERTEX3D
#define fglVertex3dv                   FGLVERTEX3DV
#define fglVertex3f                    FGLVERTEX3F
#define fglVertex3fv                   FGLVERTEX3FV
#define fglVertex3i                    FGLVERTEX3I
#define fglVertex3iv                   FGLVERTEX3IV
#define fglVertex3s                    FGLVERTEX3S
#define f9aglVertex3s                  F9AGLVERTEX3S
#define fglVertex3sv                   FGLVERTEX3SV
#define f9aglVertex3sv                 F9AGLVERTEX3SV
#define fglVertex4d                    FGLVERTEX4D
#define fglVertex4dv                   FGLVERTEX4DV
#define fglVertex4f                    FGLVERTEX4F
#define fglVertex4fv                   FGLVERTEX4FV
#define fglVertex4i                    FGLVERTEX4I
#define fglVertex4iv                   FGLVERTEX4IV
#define fglVertex4s                    FGLVERTEX4S
#define f9aglVertex4s                  F9AGLVERTEX4S
#define fglVertex4sv                   FGLVERTEX4SV
#define f9aglVertex4sv                 F9AGLVERTEX4SV
#define fglViewport                    FGLVIEWPORT
#endif

/* OpenGL 1.1 */

#ifdef OPENGL_1_1
#if FNAME==LOWERCASE
#define fglAreTexturesResident         fglaretexturesresident
#define f9eglAreTexturesResident       f9eglaretexturesresident
#define fglArrayElement                fglarrayelement
#define fglBindTexture                 fglbindtexture
#define f9y0glColorPointer             f9y0glcolorpointer
#define f9y1glColorPointer             f9y1glcolorpointer
#define f9y4glColorPointer             f9y4glcolorpointer
#define f9y5glColorPointer             f9y5glcolorpointer
#define f9y6glColorPointer             f9y6glcolorpointer
#define f9y7glColorPointer             f9y7glcolorpointer
#define fglCopyTexImage1D              fglcopyteximage1d
#define fglCopyTexImage2D              fglcopyteximage2d
#define fglCopyTexSubImage1D           fglcopytexsubimage1d
#define fglCopyTexSubImage2D           fglcopytexsubimage2d
#define fglDeleteTextures              fgldeletetextures
#define fglDisableClientState          fgldisableclientstate
#define fglDrawArrays                  fgldrawarrays
#define f9y0glDrawElements             f9y0gldrawelements
#define f9y4glDrawElements             f9y4gldrawelements
#define f9y5glDrawElements             f9y5gldrawelements
#define f9y6glDrawElements             f9y6gldrawelements
#define fglEdgeFlagPointer             fgledgeflagpointer
#define f9eglEdgeFlagPointer           f9egledgeflagpointer
#define fglEnableClientState           fglenableclientstate
#define fglGenTextures                 fglgentextures
#define fglGetPointerv                 fglgetpointerv
#define f9y0glIndexPointer             f9y0glindexpointer
#define f9y1glIndexPointer             f9y1glindexpointer
#define f9y4glIndexPointer             f9y4glindexpointer
#define f9y6glIndexPointer             f9y6glindexpointer
#define f9y7glIndexPointer             f9y7glindexpointer
#define fglIndexub                     fglindexub
#define f9bglIndexub                   f9bglindexub
#define f9cglIndexub                   f9cglindexub
#define fglIndexubv                    fglindexubv
#define f9bglIndexubv                  f9bglindexubv
#define f9cglIndexubv                  f9cglindexubv
#define f9y0glInterleavedArrays        f9y0glinterleavedarrays
#define f9y1glInterleavedArrays        f9y1glinterleavedarrays
#define f9y4glInterleavedArrays        f9y4glinterleavedarrays
#define f9y5glInterleavedArrays        f9y5glinterleavedarrays
#define f9y6glInterleavedArrays        f9y6glinterleavedarrays
#define fglIsTexture                   fglistexture
#define f9eglIsTexture                 f9eglistexture
#define f9y0glNormalPointer            f9y0glnormalpointer
#define f9y1glNormalPointer            f9y1glnormalpointer
#define f9y4glNormalPointer            f9y4glnormalpointer
#define f9y5glNormalPointer            f9y5glnormalpointer
#define f9y6glNormalPointer            f9y6glnormalpointer
#define f9y7glNormalPointer            f9y7glnormalpointer
#define fglPolygonOffset               fglpolygonoffset
#define fglPopClientAttrib             fglpopclientattrib
#define fglPrioritizeTextures          fglprioritizetextures
#define fglPushClientAttrib            fglpushclientattrib
#define f9y0glTexCoordPointer          f9y0gltexcoordpointer
#define f9y1glTexCoordPointer          f9y1gltexcoordpointer
#define f9y4glTexCoordPointer          f9y4gltexcoordpointer
#define f9y6glTexCoordPointer          f9y6gltexcoordpointer
#define f9y7glTexCoordPointer          f9y7gltexcoordpointer
#define f9y0glTexSubImage1D            f9y0gltexsubimage1d
#define f9y1glTexSubImage1D            f9y1gltexsubimage1d
#define f9y4glTexSubImage1D            f9y4gltexsubimage1d
#define f9y5glTexSubImage1D            f9y5gltexsubimage1d
#define f9y6glTexSubImage1D            f9y6gltexsubimage1d
#define f9y0glTexSubImage2D            f9y0gltexsubimage2d
#define f9y1glTexSubImage2D            f9y1gltexsubimage2d
#define f9y4glTexSubImage2D            f9y4gltexsubimage2d
#define f9y5glTexSubImage2D            f9y5gltexsubimage2d
#define f9y6glTexSubImage2D            f9y6gltexsubimage2d
#define f9y0glVertexPointer            f9y0glvertexpointer
#define f9y1glVertexPointer            f9y1glvertexpointer
#define f9y4glVertexPointer            f9y4glvertexpointer
#define f9y6glVertexPointer            f9y6glvertexpointer
#define f9y7glVertexPointer            f9y7glvertexpointer
#elif FNAME==UNDERSCORE
#define fglAreTexturesResident         fglaretexturesresident_
#define f9eglAreTexturesResident       f9eglaretexturesresident_
#define fglArrayElement                fglarrayelement_
#define fglBindTexture                 fglbindtexture_
#define f9y0glColorPointer             f9y0glcolorpointer_
#define f9y1glColorPointer             f9y1glcolorpointer_
#define f9y4glColorPointer             f9y4glcolorpointer_
#define f9y5glColorPointer             f9y5glcolorpointer_
#define f9y6glColorPointer             f9y6glcolorpointer_
#define f9y7glColorPointer             f9y7glcolorpointer_
#define fglCopyTexImage1D              fglcopyteximage1d_
#define fglCopyTexImage2D              fglcopyteximage2d_
#define fglCopyTexSubImage1D           fglcopytexsubimage1d_
#define fglCopyTexSubImage2D           fglcopytexsubimage2d_
#define fglDeleteTextures              fgldeletetextures_
#define fglDisableClientState          fgldisableclientstate_
#define fglDrawArrays                  fgldrawarrays_
#define f9y0glDrawElements             f9y0gldrawelements_
#define f9y4glDrawElements             f9y4gldrawelements_
#define f9y5glDrawElements             f9y5gldrawelements_
#define f9y6glDrawElements             f9y6gldrawelements_
#define fglEdgeFlagPointer             fgledgeflagpointer_
#define f9eglEdgeFlagPointer           f9egledgeflagpointer_
#define fglEnableClientState           fglenableclientstate_
#define fglGenTextures                 fglgentextures_
#define fglGetPointerv                 fglgetpointerv_
#define f9y0glIndexPointer             f9y0glindexpointer_
#define f9y1glIndexPointer             f9y1glindexpointer_
#define f9y4glIndexPointer             f9y4glindexpointer_
#define f9y6glIndexPointer             f9y6glindexpointer_
#define f9y7glIndexPointer             f9y7glindexpointer_
#define fglIndexub                     fglindexub_
#define f9bglIndexub                   f9bglindexub_
#define f9cglIndexub                   f9cglindexub_
#define fglIndexubv                    fglindexubv_
#define f9bglIndexubv                  f9bglindexubv_
#define f9cglIndexubv                  f9cglindexubv_
#define f9y0glInterleavedArrays        f9y0glinterleavedarrays_
#define f9y1glInterleavedArrays        f9y1glinterleavedarrays_
#define f9y4glInterleavedArrays        f9y4glinterleavedarrays_
#define f9y5glInterleavedArrays        f9y5glinterleavedarrays_
#define f9y6glInterleavedArrays        f9y6glinterleavedarrays_
#define fglIsTexture                   fglistexture_
#define f9eglIsTexture                 f9eglistexture_
#define f9y0glNormalPointer            f9y0glnormalpointer_
#define f9y1glNormalPointer            f9y1glnormalpointer_
#define f9y4glNormalPointer            f9y4glnormalpointer_
#define f9y5glNormalPointer            f9y5glnormalpointer_
#define f9y6glNormalPointer            f9y6glnormalpointer_
#define f9y7glNormalPointer            f9y7glnormalpointer_
#define fglPolygonOffset               fglpolygonoffset_
#define fglPopClientAttrib             fglpopclientattrib_
#define fglPrioritizeTextures          fglprioritizetextures_
#define fglPushClientAttrib            fglpushclientattrib_
#define f9y0glTexCoordPointer          f9y0gltexcoordpointer_
#define f9y1glTexCoordPointer          f9y1gltexcoordpointer_
#define f9y4glTexCoordPointer          f9y4gltexcoordpointer_
#define f9y6glTexCoordPointer          f9y6gltexcoordpointer_
#define f9y7glTexCoordPointer          f9y7gltexcoordpointer_
#define f9y0glTexSubImage1D            f9y0gltexsubimage1d_
#define f9y1glTexSubImage1D            f9y1gltexsubimage1d_
#define f9y4glTexSubImage1D            f9y4gltexsubimage1d_
#define f9y5glTexSubImage1D            f9y5gltexsubimage1d_
#define f9y6glTexSubImage1D            f9y6gltexsubimage1d_
#define f9y0glTexSubImage2D            f9y0gltexsubimage2d_
#define f9y1glTexSubImage2D            f9y1gltexsubimage2d_
#define f9y4glTexSubImage2D            f9y4gltexsubimage2d_
#define f9y5glTexSubImage2D            f9y5gltexsubimage2d_
#define f9y6glTexSubImage2D            f9y6gltexsubimage2d_
#define f9y0glVertexPointer            f9y0glvertexpointer_
#define f9y1glVertexPointer            f9y1glvertexpointer_
#define f9y4glVertexPointer            f9y4glvertexpointer_
#define f9y6glVertexPointer            f9y6glvertexpointer_
#define f9y7glVertexPointer            f9y7glvertexpointer_
#elif FNAME==UPPERCASE
#define fglAreTexturesResident         FGLARETEXTURESRESIDENT
#define f9eglAreTexturesResident       F9EGLARETEXTURESRESIDENT
#define fglArrayElement                FGLARRAYELEMENT
#define fglBindTexture                 FGLBINDTEXTURE
#define f9y0glColorPointer             F9Y0GLCOLORPOINTER
#define f9y1glColorPointer             F9Y1GLCOLORPOINTER
#define f9y4glColorPointer             F9Y4GLCOLORPOINTER
#define f9y5glColorPointer             F9Y5GLCOLORPOINTER
#define f9y6glColorPointer             F9Y6GLCOLORPOINTER
#define f9y7glColorPointer             F9Y7GLCOLORPOINTER
#define fglCopyTexImage1D              FGLCOPYTEXIMAGE1D
#define fglCopyTexImage2D              FGLCOPYTEXIMAGE2D
#define fglCopyTexSubImage1D           FGLCOPYTEXSUBIMAGE1D
#define fglCopyTexSubImage2D           FGLCOPYTEXSUBIMAGE2D
#define fglDeleteTextures              FGLDELETETEXTURES
#define fglDisableClientState          FGLDISABLECLIENTSTATE
#define fglDrawArrays                  FGLDRAWARRAYS
#define f9y0glDrawElements             F9Y0GLDRAWELEMENTS
#define f9y4glDrawElements             F9Y4GLDRAWELEMENTS
#define f9y5glDrawElements             F9Y5GLDRAWELEMENTS
#define f9y6glDrawElements             F9Y6GLDRAWELEMENTS
#define fglEdgeFlagPointer             FGLEDGEFLAGPOINTER
#define f9eglEdgeFlagPointer           F9EGLEDGEFLAGPOINTER
#define fglEnableClientState           FGLENABLECLIENTSTATE
#define fglGenTextures                 FGLGENTEXTURES
#define fglGetPointerv                 FGLGETPOINTERV
#define f9y0glIndexPointer             F9Y0GLINDEXPOINTER
#define f9y1glIndexPointer             F9Y1GLINDEXPOINTER
#define f9y4glIndexPointer             F9Y4GLINDEXPOINTER
#define f9y6glIndexPointer             F9Y6GLINDEXPOINTER
#define f9y7glIndexPointer             F9Y7GLINDEXPOINTER
#define fglIndexub                     FGLINDEXUB
#define f9bglIndexub                   F9BGLINDEXUB
#define f9cglIndexub                   F9CGLINDEXUB
#define fglIndexubv                    FGLINDEXUBV
#define f9bglIndexubv                  F9BGLINDEXUBV
#define f9cglIndexubv                  F9CGLINDEXUBV
#define f9y0glInterleavedArrays        F9Y0GLINTERLEAVEDARRAYS
#define f9y1glInterleavedArrays        F9Y1GLINTERLEAVEDARRAYS
#define f9y4glInterleavedArrays        F9Y4GLINTERLEAVEDARRAYS
#define f9y5glInterleavedArrays        F9Y5GLINTERLEAVEDARRAYS
#define f9y6glInterleavedArrays        F9Y6GLINTERLEAVEDARRAYS
#define fglIsTexture                   FGLISTEXTURE
#define f9eglIsTexture                 F9EGLISTEXTURE
#define f9y0glNormalPointer            F9Y0GLNORMALPOINTER
#define f9y1glNormalPointer            F9Y1GLNORMALPOINTER
#define f9y4glNormalPointer            F9Y4GLNORMALPOINTER
#define f9y5glNormalPointer            F9Y5GLNORMALPOINTER
#define f9y6glNormalPointer            F9Y6GLNORMALPOINTER
#define f9y7glNormalPointer            F9Y7GLNORMALPOINTER
#define fglPolygonOffset               FGLPOLYGONOFFSET
#define fglPopClientAttrib             FGLPOPCLIENTATTRIB
#define fglPrioritizeTextures          FGLPRIORITIZETEXTURES
#define fglPushClientAttrib            FGLPUSHCLIENTATTRIB
#define f9y0glTexCoordPointer          F9Y0GLTEXCOORDPOINTER
#define f9y1glTexCoordPointer          F9Y1GLTEXCOORDPOINTER
#define f9y4glTexCoordPointer          F9Y4GLTEXCOORDPOINTER
#define f9y6glTexCoordPointer          F9Y6GLTEXCOORDPOINTER
#define f9y7glTexCoordPointer          F9Y7GLTEXCOORDPOINTER
#define f9y0glTexSubImage1D            F9Y0GLTEXSUBIMAGE1D
#define f9y1glTexSubImage1D            F9Y1GLTEXSUBIMAGE1D
#define f9y4glTexSubImage1D            F9Y4GLTEXSUBIMAGE1D
#define f9y5glTexSubImage1D            F9Y5GLTEXSUBIMAGE1D
#define f9y6glTexSubImage1D            F9Y6GLTEXSUBIMAGE1D
#define f9y0glTexSubImage2D            F9Y0GLTEXSUBIMAGE2D
#define f9y1glTexSubImage2D            F9Y1GLTEXSUBIMAGE2D
#define f9y4glTexSubImage2D            F9Y4GLTEXSUBIMAGE2D
#define f9y5glTexSubImage2D            F9Y5GLTEXSUBIMAGE2D
#define f9y6glTexSubImage2D            F9Y6GLTEXSUBIMAGE2D
#define f9y0glVertexPointer            F9Y0GLVERTEXPOINTER
#define f9y1glVertexPointer            F9Y1GLVERTEXPOINTER
#define f9y4glVertexPointer            F9Y4GLVERTEXPOINTER
#define f9y6glVertexPointer            F9Y6GLVERTEXPOINTER
#define f9y7glVertexPointer            F9Y7GLVERTEXPOINTER
#endif
#endif

/* OpenGL 1.2 */

#ifdef OPENGL_1_2
#if FNAME==LOWERCASE
#define f9y0glDrawRangeElements        f9y0gldrawrangeelements
#define f9y4glDrawRangeElements        f9y4gldrawrangeelements
#define f9y5glDrawRangeElements        f9y5gldrawrangeelements
#define f9y6glDrawRangeElements        f9y6gldrawrangeelements
#define f9y0glTexImage3D               f9y0glteximage3d
#define f9y1glTexImage3D               f9y1glteximage3d
#define f9y4glTexImage3D               f9y4glteximage3d
#define f9y5glTexImage3D               f9y5glteximage3d
#define f9y6glTexImage3D               f9y6glteximage3d
#define f9y0glTexSubImage3D            f9y0gltexsubimage3d
#define f9y1glTexSubImage3D            f9y1gltexsubimage3d
#define f9y4glTexSubImage3D            f9y4gltexsubimage3d
#define f9y5glTexSubImage3D            f9y5gltexsubimage3d
#define f9y6glTexSubImage3D            f9y6gltexsubimage3d
#define fglCopyTexSubImage3D           fglcopytexsubimage3d
#elif FNAME==UNDERSCORE
#define f9y0glDrawRangeElements        f9y0gldrawrangeelements_
#define f9y4glDrawRangeElements        f9y4gldrawrangeelements_
#define f9y5glDrawRangeElements        f9y5gldrawrangeelements_
#define f9y6glDrawRangeElements        f9y6gldrawrangeelements_
#define f9y0glTexImage3D               f9y0glteximage3d_
#define f9y1glTexImage3D               f9y1glteximage3d_
#define f9y4glTexImage3D               f9y4glteximage3d_
#define f9y5glTexImage3D               f9y5glteximage3d_
#define f9y6glTexImage3D               f9y6glteximage3d_
#define f9y0glTexSubImage3D            f9y0gltexsubimage3d_
#define f9y1glTexSubImage3D            f9y1gltexsubimage3d_
#define f9y4glTexSubImage3D            f9y4gltexsubimage3d_
#define f9y5glTexSubImage3D            f9y5gltexsubimage3d_
#define f9y6glTexSubImage3D            f9y6gltexsubimage3d_
#define fglCopyTexSubImage3D           fglcopytexsubimage3d_
#elif FNAME==UPPERCASE
#define f9y0glDrawRangeElements        F9Y0GLDRAWRANGEELEMENTS
#define f9y4glDrawRangeElements        F9Y4GLDRAWRANGEELEMENTS
#define f9y5glDrawRangeElements        F9Y5GLDRAWRANGEELEMENTS
#define f9y6glDrawRangeElements        F9Y6GLDRAWRANGEELEMENTS
#define f9y0glTexImage3D               F9Y0GLTEXIMAGE3D
#define f9y1glTexImage3D               F9Y1GLTEXIMAGE3D
#define f9y4glTexImage3D               F9Y4GLTEXIMAGE3D
#define f9y5glTexImage3D               F9Y5GLTEXIMAGE3D
#define f9y6glTexImage3D               F9Y6GLTEXIMAGE3D
#define f9y0glTexSubImage3D            F9Y0GLTEXSUBIMAGE3D
#define f9y1glTexSubImage3D            F9Y1GLTEXSUBIMAGE3D
#define f9y4glTexSubImage3D            F9Y4GLTEXSUBIMAGE3D
#define f9y5glTexSubImage3D            F9Y5GLTEXSUBIMAGE3D
#define f9y6glTexSubImage3D            F9Y6GLTEXSUBIMAGE3D
#define fglCopyTexSubImage3D           FGLCOPYTEXSUBIMAGE3D
#endif
#endif

/* EXT_blend_color */

#ifdef GL_EXT_blend_color
#if FNAME==LOWERCASE
#define fglBlendColorEXT               fglblendcolorext
#elif FNAME==UNDERSCORE
#define fglBlendColorEXT               fglblendcolorext_
#elif FNAME==UPPERCASE
#define fglBlendColorEXT               FGLBLENDCOLOREXT
#endif
#endif

/* EXT_blend_minmax */

#ifdef GL_EXT_blend_minmax
#if FNAME==LOWERCASE
#define fglBlendEquationEXT            fglblendequationext
#elif FNAME==UNDERSCORE
#define fglBlendEquationEXT            fglblendequationext_
#elif FNAME==UPPERCASE
#define fglBlendEquationEXT            FGLBLENDEQUATIONEXT
#endif
#endif

/* EXT_color_table */

#ifdef GL_EXT_color_table
#if FNAME==LOWERCASE
#define f9y0glColorTableEXT            f9y0glcolortableext
#define f9y1glColorTableEXT            f9y1glcolortableext
#define f9y4glColorTableEXT            f9y4glcolortableext
#define f9y5glColorTableEXT            f9y5glcolortableext
#define f9y6glColorTableEXT            f9y6glcolortableext
#define f9y0glColorSubTableEXT         f9y0glcolorsubtableext
#define f9y1glColorSubTableEXT         f9y1glcolorsubtableext
#define f9y4glColorSubTableEXT         f9y4glcolorsubtableext
#define f9y5glColorSubTableEXT         f9y5glcolorsubtableext
#define f9y6glColorSubTableEXT         f9y6glcolorsubtableext
#define f9y0glGetColorTableEXT         f9y0glgetcolortableext
#define f9y1glGetColorTableEXT         f9y1glgetcolortableext
#define f9y4glGetColorTableEXT         f9y4glgetcolortableext
#define f9y5glGetColorTableEXT         f9y5glgetcolortableext
#define f9y6glGetColorTableEXT         f9y6glgetcolortableext
#define fglGetColorTableParameterfvEXT fglgetcolortableparameterfvext
#define fglGetColorTableParameterivEXT fglgetcolortableparameterivext
#elif FNAME==UNDERSCORE
#define f9y0glColorTableEXT            f9y0glcolortableext_
#define f9y1glColorTableEXT            f9y1glcolortableext_
#define f9y4glColorTableEXT            f9y4glcolortableext_
#define f9y5glColorTableEXT            f9y5glcolortableext_
#define f9y6glColorTableEXT            f9y6glcolortableext_
#define f9y0glColorSubTableEXT         f9y0glcolorsubtableext_
#define f9y1glColorSubTableEXT         f9y1glcolorsubtableext_
#define f9y4glColorSubTableEXT         f9y4glcolorsubtableext_
#define f9y5glColorSubTableEXT         f9y5glcolorsubtableext_
#define f9y6glColorSubTableEXT         f9y6glcolorsubtableext_
#define f9y0glGetColorTableEXT         f9y0glgetcolortableext_
#define f9y1glGetColorTableEXT         f9y1glgetcolortableext_
#define f9y4glGetColorTableEXT         f9y4glgetcolortableext_
#define f9y5glGetColorTableEXT         f9y5glgetcolortableext_
#define f9y6glGetColorTableEXT         f9y6glgetcolortableext_
#define fglGetColorTableParameterfvEXT fglgetcolortableparameterfvext_
#define fglGetColorTableParameterivEXT fglgetcolortableparameterivext_
#elif FNAME==UPPERCASE
#define f9y0glColorTableEXT            F9Y0GLCOLORTABLEEXT
#define f9y1glColorTableEXT            F9Y1GLCOLORTABLEEXT
#define f9y4glColorTableEXT            F9Y4GLCOLORTABLEEXT
#define f9y5glColorTableEXT            F9Y5GLCOLORTABLEEXT
#define f9y6glColorTableEXT            F9Y6GLCOLORTABLEEXT
#define f9y0glColorSubTableEXT         F9Y0GLCOLORSUBTABLEEXT
#define f9y1glColorSubTableEXT         F9Y1GLCOLORSUBTABLEEXT
#define f9y4glColorSubTableEXT         F9Y4GLCOLORSUBTABLEEXT
#define f9y5glColorSubTableEXT         F9Y5GLCOLORSUBTABLEEXT
#define f9y6glColorSubTableEXT         F9Y6GLCOLORSUBTABLEEXT
#define f9y0glGetColorTableEXT         F9Y0GLGETCOLORTABLEEXT
#define f9y1glGetColorTableEXT         F9Y1GLGETCOLORTABLEEXT
#define f9y4glGetColorTableEXT         F9Y4GLGETCOLORTABLEEXT
#define f9y5glGetColorTableEXT         F9Y5GLGETCOLORTABLEEXT
#define f9y6glGetColorTableEXT         F9Y6GLGETCOLORTABLEEXT
#define fglGetColorTableParameterfvEXT FGLGETCOLORTABLEPARAMETERFVEXT
#define fglGetColorTableParameterivEXT FGLGETCOLORTABLEPARAMETERIVEXT
#endif
#endif

/* EXT_point_parameters */

#ifdef GL_EXT_point_parameters
#if FNAME==LOWERCASE
#define fglPointParameterfEXT          fglpointparameterfext
#define fglPointParameterfvEXT         fglpointparameterfvext
#elif FNAME==UNDERSCORE
#define fglPointParameterfEXT          fglpointparameterfext_
#define fglPointParameterfvEXT         fglpointparameterfvext_
#elif FNAME==UPPERCASE
#define fglPointParameterfEXT          FGLPOINTPARAMETERFEXT
#define fglPointParameterfvEXT         FGLPOINTPARAMETERFVEXT
#endif
#endif

/* ARB_point_parameters (experimental)*/

#ifdef GL_ARB_point_parameters
#if FNAME==LOWERCASE
#define fglPointParameterfARB          fglpointparameterfarb
#define fglPointParameterfvARB         fglpointparameterfvarb
#elif FNAME==UNDERSCORE
#define fglPointParameterfARB          fglpointparameterfarb_
#define fglPointParameterfvARB         fglpointparameterfvarb_
#elif FNAME==UPPERCASE
#define fglPointParameterfARB          FGLPOINTPARAMETERFARB
#define fglPointParameterfvARB         FGLPOINTPARAMETERFVARB
#endif
#endif

/* EXT_polygon_offset */

#ifdef GL_EXT_polygon_offset
#if FNAME==LOWERCASE
#define fglPolygonOffsetEXT            fglpolygonoffsetext
#elif FNAME==UNDERSCORE
#define fglPolygonOffsetEXT            fglpolygonoffsetext_
#elif FNAME==UPPERCASE
#define fglPolygonOffsetEXT            FGLPOLYGONOFFSETEXT
#endif
#endif

/* EXT_subtexture */

#ifdef GL_EXT_subtexture
#if FNAME==LOWERCASE
#define f9y0glTexSubImage1DEXT         f9y0gltexsubimage1dext
#define f9y1glTexSubImage1DEXT         f9y1gltexsubimage1dext
#define f9y4glTexSubImage1DEXT         f9y4gltexsubimage1dext
#define f9y5glTexSubImage1DEXT         f9y5gltexsubimage1dext
#define f9y6glTexSubImage1DEXT         f9y6gltexsubimage1dext
#define f9y0glTexSubImage2DEXT         f9y0gltexsubimage2dext
#define f9y1glTexSubImage2DEXT         f9y1gltexsubimage2dext
#define f9y4glTexSubImage2DEXT         f9y4gltexsubimage2dext
#define f9y5glTexSubImage2DEXT         f9y5gltexsubimage2dext
#define f9y6glTexSubImage2DEXT         f9y6gltexsubimage2dext
#elif FNAME==UNDERSCORE
#define f9y0glTexSubImage1DEXT         f9y0gltexsubimage1dext_
#define f9y1glTexSubImage1DEXT         f9y1gltexsubimage1dext_
#define f9y4glTexSubImage1DEXT         f9y4gltexsubimage1dext_
#define f9y5glTexSubImage1DEXT         f9y5gltexsubimage1dext_
#define f9y6glTexSubImage1DEXT         f9y6gltexsubimage1dext_
#define f9y0glTexSubImage2DEXT         f9y0gltexsubimage2dext_
#define f9y1glTexSubImage2DEXT         f9y1gltexsubimage2dext_
#define f9y4glTexSubImage2DEXT         f9y4gltexsubimage2dext_
#define f9y5glTexSubImage2DEXT         f9y5gltexsubimage2dext_
#define f9y6glTexSubImage2DEXT         f9y6gltexsubimage2dext_
#elif FNAME==UPPERCASE
#define f9y0glTexSubImage1DEXT         F9Y0GLTEXSUBIMAGE1DEXT
#define f9y1glTexSubImage1DEXT         F9Y1GLTEXSUBIMAGE1DEXT
#define f9y4glTexSubImage1DEXT         F9Y4GLTEXSUBIMAGE1DEXT
#define f9y5glTexSubImage1DEXT         F9Y5GLTEXSUBIMAGE1DEXT
#define f9y6glTexSubImage1DEXT         F9Y6GLTEXSUBIMAGE1DEXT
#define f9y0glTexSubImage2DEXT         F9Y0GLTEXSUBIMAGE2DEXT
#define f9y1glTexSubImage2DEXT         F9Y1GLTEXSUBIMAGE2DEXT
#define f9y4glTexSubImage2DEXT         F9Y4GLTEXSUBIMAGE2DEXT
#define f9y5glTexSubImage2DEXT         F9Y5GLTEXSUBIMAGE2DEXT
#define f9y6glTexSubImage2DEXT         F9Y6GLTEXSUBIMAGE2DEXT
#endif
#endif

/* EXT_texture3D */

#ifdef GL_EXT_texture3D
#if FNAME==LOWERCASE
#define f9y0glTexImage3DEXT            f9y0glteximage3dext
#define f9y1glTexImage3DEXT            f9y1glteximage3dext
#define f9y4glTexImage3DEXT            f9y4glteximage3dext
#define f9y5glTexImage3DEXT            f9y5glteximage3dext
#define f9y6glTexImage3DEXT            f9y6glteximage3dext
#define f9y0glTexSubImage3DEXT         f9y0gltexsubimage3dext
#define f9y1glTexSubImage3DEXT         f9y1gltexsubimage3dext
#define f9y4glTexSubImage3DEXT         f9y4gltexsubimage3dext
#define f9y5glTexSubImage3DEXT         f9y5gltexsubimage3dext
#define f9y6glTexSubImage3DEXT         f9y6gltexsubimage3dext
#define fglCopyTexSubImage3DEXT        fglcopytexsubimage3dext
#elif FNAME==UNDERSCORE
#define f9y0glTexImage3DEXT            f9y0glteximage3dext_
#define f9y1glTexImage3DEXT            f9y1glteximage3dext_
#define f9y4glTexImage3DEXT            f9y4glteximage3dext_
#define f9y5glTexImage3DEXT            f9y5glteximage3dext_
#define f9y6glTexImage3DEXT            f9y6glteximage3dext_
#define f9y0glTexSubImage3DEXT         f9y0gltexsubimage3dext_
#define f9y1glTexSubImage3DEXT         f9y1gltexsubimage3dext_
#define f9y4glTexSubImage3DEXT         f9y4gltexsubimage3dext_
#define f9y5glTexSubImage3DEXT         f9y5gltexsubimage3dext_
#define f9y6glTexSubImage3DEXT         f9y6gltexsubimage3dext_
#define fglCopyTexSubImage3DEXT        fglcopytexsubimage3dext_
#elif FNAME==UPPERCASE
#define f9y0glTexImage3DEXT            F9Y0GLTEXIMAGE3DEXT
#define f9y1glTexImage3DEXT            F9Y1GLTEXIMAGE3DEXT
#define f9y4glTexImage3DEXT            F9Y4GLTEXIMAGE3DEXT
#define f9y5glTexImage3DEXT            F9Y5GLTEXIMAGE3DEXT
#define f9y6glTexImage3DEXT            F9Y6GLTEXIMAGE3DEXT
#define f9y0glTexSubImage3DEXT         F9Y0GLTEXSUBIMAGE3DEXT
#define f9y1glTexSubImage3DEXT         F9Y1GLTEXSUBIMAGE3DEXT
#define f9y4glTexSubImage3DEXT         F9Y4GLTEXSUBIMAGE3DEXT
#define f9y5glTexSubImage3DEXT         F9Y5GLTEXSUBIMAGE3DEXT
#define f9y6glTexSubImage3DEXT         F9Y6GLTEXSUBIMAGE3DEXT
#define fglCopyTexSubImage3DEXT        FGLCOPYTEXSUBIMAGE3DEXT
#endif
#endif

/* EXT_texture_object */

#ifdef GL_EXT_texture_object
#if FNAME==LOWERCASE
#define fglGenTexturesEXT              fglgentexturesext
#define fglDeleteTexturesEXT           fgldeletetexturesext
#define fglBindTextureEXT              fglbindtextureext
#define fglPrioritizeTexturesEXT       fglprioritizetexturesext
#define fglAreTexturesResidentEXT      fglaretexturesresidentext
#define f9eglAreTexturesResidentEXT    f9eglaretexturesresidentext
#define fglIsTextureEXT                fglistextureext
#define f9eglIsTextureEXT              f9eglistextureext
#elif FNAME==UNDERSCORE
#define fglGenTexturesEXT              fglgentexturesext_
#define fglDeleteTexturesEXT           fgldeletetexturesext_
#define fglBindTextureEXT              fglbindtextureext_
#define fglPrioritizeTexturesEXT       fglprioritizetexturesext_
#define fglAreTexturesResidentEXT      fglaretexturesresidentext_
#define f9eglAreTexturesResidentEXT    f9eglaretexturesresidentext_
#define fglIsTextureEXT                fglistextureext_
#define f9eglIsTextureEXT              f9eglistextureext_
#elif FNAME==UPPERCASE
#define fglGenTexturesEXT              FGLGENTEXTURESEXT
#define fglDeleteTexturesEXT           FGLDELETETEXTURESEXT
#define fglBindTextureEXT              FGLBINDTEXTUREEXT
#define fglPrioritizeTexturesEXT       FGLPRIORITIZETEXTURESEXT
#define fglAreTexturesResidentEXT      FGLARETEXTURESRESIDENTEXT
#define f9eglAreTexturesResidentEXT    F9EGLARETEXTURESRESIDENTEXT
#define fglIsTextureEXT                FGLISTEXTUREEXT
#define f9eglIsTextureEXT              F9EGLISTEXTUREEXT
#endif
#endif

/* EXT_vertex_array */

#ifdef GL_EXT_vertex_array
#if FNAME==LOWERCASE
#define f9y0glVertexPointerEXT         f9y0glvertexpointerext
#define f9y1glVertexPointerEXT         f9y1glvertexpointerext
#define f9y4glVertexPointerEXT         f9y4glvertexpointerext
#define f9y6glVertexPointerEXT         f9y6glvertexpointerext
#define f9y7glVertexPointerEXT         f9y7glvertexpointerext
#define f9y0glNormalPointerEXT         f9y0glnormalpointerext
#define f9y1glNormalPointerEXT         f9y1glnormalpointerext
#define f9y4glNormalPointerEXT         f9y4glnormalpointerext
#define f9y5glNormalPointerEXT         f9y5glnormalpointerext
#define f9y6glNormalPointerEXT         f9y6glnormalpointerext
#define f9y7glNormalPointerEXT         f9y7glnormalpointerext
#define f9y0glColorPointerEXT          f9y0glcolorpointerext
#define f9y1glColorPointerEXT          f9y1glcolorpointerext
#define f9y4glColorPointerEXT          f9y4glcolorpointerext
#define f9y5glColorPointerEXT          f9y5glcolorpointerext
#define f9y6glColorPointerEXT          f9y6glcolorpointerext
#define f9y7glColorPointerEXT          f9y7glcolorpointerext
#define f9y0glIndexPointerEXT          f9y0glindexpointerext
#define f9y1glIndexPointerEXT          f9y1glindexpointerext
#define f9y4glIndexPointerEXT          f9y4glindexpointerext
#define f9y6glIndexPointerEXT          f9y6glindexpointerext
#define f9y7glIndexPointerEXT          f9y7glindexpointerext
#define f9y0glTexCoordPointerEXT       f9y0gltexcoordpointerext
#define f9y1glTexCoordPointerEXT       f9y1gltexcoordpointerext
#define f9y4glTexCoordPointerEXT       f9y4gltexcoordpointerext
#define f9y6glTexCoordPointerEXT       f9y6gltexcoordpointerext
#define f9y7glTexCoordPointerEXT       f9y7gltexcoordpointerext
#define fglEdgeFlagPointerEXT          fgledgeflagpointerext
#define f9eglEdgeFlagPointerEXT        f9egledgeflagpointerext
#define fglGetPointervEXT              fglgetpointervext
#define fglArrayElementEXT             fglarrayelementext
#define fglDrawArraysEXT               fgldrawarraysext
#elif FNAME==UNDERSCORE
#define f9y0glVertexPointerEXT         f9y0glvertexpointerext_
#define f9y1glVertexPointerEXT         f9y1glvertexpointerext_
#define f9y4glVertexPointerEXT         f9y4glvertexpointerext_
#define f9y6glVertexPointerEXT         f9y6glvertexpointerext_
#define f9y7glVertexPointerEXT         f9y7glvertexpointerext_
#define f9y0glNormalPointerEXT         f9y0glnormalpointerext_
#define f9y1glNormalPointerEXT         f9y1glnormalpointerext_
#define f9y4glNormalPointerEXT         f9y4glnormalpointerext_
#define f9y5glNormalPointerEXT         f9y5glnormalpointerext_
#define f9y6glNormalPointerEXT         f9y6glnormalpointerext_
#define f9y7glNormalPointerEXT         f9y7glnormalpointerext_
#define f9y0glColorPointerEXT          f9y0glcolorpointerext_
#define f9y1glColorPointerEXT          f9y1glcolorpointerext_
#define f9y4glColorPointerEXT          f9y4glcolorpointerext_
#define f9y5glColorPointerEXT          f9y5glcolorpointerext_
#define f9y6glColorPointerEXT          f9y6glcolorpointerext_
#define f9y7glColorPointerEXT          f9y7glcolorpointerext_
#define f9y0glIndexPointerEXT          f9y0glindexpointerext_
#define f9y1glIndexPointerEXT          f9y1glindexpointerext_
#define f9y4glIndexPointerEXT          f9y4glindexpointerext_
#define f9y6glIndexPointerEXT          f9y6glindexpointerext_
#define f9y7glIndexPointerEXT          f9y7glindexpointerext_
#define f9y0glTexCoordPointerEXT       f9y0gltexcoordpointerext_
#define f9y1glTexCoordPointerEXT       f9y1gltexcoordpointerext_
#define f9y4glTexCoordPointerEXT       f9y4gltexcoordpointerext_
#define f9y6glTexCoordPointerEXT       f9y6gltexcoordpointerext_
#define f9y7glTexCoordPointerEXT       f9y7gltexcoordpointerext_
#define fglEdgeFlagPointerEXT          fgledgeflagpointerext_
#define f9eglEdgeFlagPointerEXT        f9egledgeflagpointerext_
#define fglGetPointervEXT              fglgetpointervext_
#define fglArrayElementEXT             fglarrayelementext_
#define fglDrawArraysEXT               fgldrawarraysext_
#elif FNAME==UPPERCASE
#define f9y0glVertexPointerEXT         F9Y0GLVERTEXPOINTEREXT
#define f9y1glVertexPointerEXT         F9Y1GLVERTEXPOINTEREXT
#define f9y4glVertexPointerEXT         F9Y4GLVERTEXPOINTEREXT
#define f9y6glVertexPointerEXT         F9Y6GLVERTEXPOINTEREXT
#define f9y7glVertexPointerEXT         F9Y7GLVERTEXPOINTEREXT
#define f9y0glNormalPointerEXT         F9Y0GLNORMALPOINTEREXT
#define f9y1glNormalPointerEXT         F9Y1GLNORMALPOINTEREXT
#define f9y4glNormalPointerEXT         F9Y4GLNORMALPOINTEREXT
#define f9y5glNormalPointerEXT         F9Y5GLNORMALPOINTEREXT
#define f9y6glNormalPointerEXT         F9Y6GLNORMALPOINTEREXT
#define f9y7glNormalPointerEXT         F9Y7GLNORMALPOINTEREXT
#define f9y0glColorPointerEXT          F9Y0GLCOLORPOINTEREXT
#define f9y1glColorPointerEXT          F9Y1GLCOLORPOINTEREXT
#define f9y4glColorPointerEXT          F9Y4GLCOLORPOINTEREXT
#define f9y5glColorPointerEXT          F9Y5GLCOLORPOINTEREXT
#define f9y6glColorPointerEXT          F9Y6GLCOLORPOINTEREXT
#define f9y7glColorPointerEXT          F9Y7GLCOLORPOINTEREXT
#define f9y0glIndexPointerEXT          F9Y0GLINDEXPOINTEREXT
#define f9y1glIndexPointerEXT          F9Y1GLINDEXPOINTEREXT
#define f9y4glIndexPointerEXT          F9Y4GLINDEXPOINTEREXT
#define f9y6glIndexPointerEXT          F9Y6GLINDEXPOINTEREXT
#define f9y7glIndexPointerEXT          F9Y7GLINDEXPOINTEREXT
#define f9y0glTexCoordPointerEXT       F9Y0GLTEXCOORDPOINTEREXT
#define f9y1glTexCoordPointerEXT       F9Y1GLTEXCOORDPOINTEREXT
#define f9y4glTexCoordPointerEXT       F9Y4GLTEXCOORDPOINTEREXT
#define f9y6glTexCoordPointerEXT       F9Y6GLTEXCOORDPOINTEREXT
#define f9y7glTexCoordPointerEXT       F9Y7GLTEXCOORDPOINTEREXT
#define fglEdgeFlagPointerEXT          FGLEDGEFLAGPOINTEREXT
#define f9eglEdgeFlagPointerEXT        F9EGLEDGEFLAGPOINTEREXT
#define fglGetPointervEXT              FGLGETPOINTERVEXT
#define fglArrayElementEXT             FGLARRAYELEMENTEXT
#define fglDrawArraysEXT               FGLDRAWARRAYSEXT
#endif
#endif

/* SGIS_multisample */

#ifdef GL_SGIS_multisample
#if FNAME==LOWERCASE
#define fglSampleMaskSGIS              fglsamplemasksgis
#define f9eglSampleMaskSGIS            f9eglsamplemasksgis
#define fglSamplePatternSGIS           fglsamplepatternsgis
#elif FNAME==UNDERSCORE
#define fglSampleMaskSGIS              fglsamplemasksgis_
#define f9eglSampleMaskSGIS            f9eglsamplemasksgis_
#define fglSamplePatternSGIS           fglsamplepatternsgis_
#elif FNAME==UPPERCASE
#define fglSampleMaskSGIS              FGLSAMPLEMASKSGIS
#define f9eglSampleMaskSGIS            F9EGLSAMPLEMASKSGIS
#define fglSamplePatternSGIS           FGLSAMPLEPATTERNSGIS
#endif
#endif

/* MESA_resize_buffers */

#ifdef GL_MESA_resize_buffers
#if FNAME==LOWERCASE
#define fglResizeBuffersMESA           fglresizebuffersmesa
#elif FNAME==UNDERSCORE
#define fglResizeBuffersMESA           fglresizebuffersmesa_
#elif FNAME==UPPERCASE
#define fglResizeBuffersMESA           FGLRESIZEBUFFERSMESA
#endif
#endif

/* MESA_window_pos */

#ifdef GL_MESA_window_pos
#if FNAME==LOWERCASE
#define fglWindowPos2iMESA             fglwindowpos2imesa
#define fglWindowPos2sMESA             fglwindowpos2smesa
#define f9aglWindowPos2sMESA           f9aglwindowpos2smesa
#define fglWindowPos2fMESA             fglwindowpos2fmesa
#define fglWindowPos2dMESA             fglwindowpos2dmesa
#define fglWindowPos2ivMESA            fglwindowpos2ivmesa
#define fglWindowPos2svMESA            fglwindowpos2svmesa
#define f9aglWindowPos2svMESA          f9aglwindowpos2svmesa
#define fglWindowPos2fvMESA            fglwindowpos2fvmesa
#define fglWindowPos2dvMESA            fglwindowpos2dvmesa
#define fglWindowPos3iMESA             fglwindowpos3imesa
#define fglWindowPos3sMESA             fglwindowpos3smesa
#define f9aglWindowPos3sMESA           f9aglwindowpos3smesa
#define fglWindowPos3fMESA             fglwindowpos3fmesa
#define fglWindowPos3dMESA             fglwindowpos3dmesa
#define fglWindowPos3ivMESA            fglwindowpos3ivmesa
#define fglWindowPos3svMESA            fglwindowpos3svmesa
#define f9aglWindowPos3svMESA          f9aglwindowpos3svmesa
#define fglWindowPos3fvMESA            fglwindowpos3fvmesa
#define fglWindowPos3dvMESA            fglwindowpos3dvmesa
#define fglWindowPos4iMESA             fglwindowpos4imesa
#define fglWindowPos4sMESA             fglwindowpos4smesa
#define f9aglWindowPos4sMESA           f9aglwindowpos4smesa
#define fglWindowPos4fMESA             fglwindowpos4fmesa
#define fglWindowPos4dMESA             fglwindowpos4dmesa
#define fglWindowPos4ivMESA            fglwindowpos4ivmesa
#define fglWindowPos4svMESA            fglwindowpos4svmesa
#define f9aglWindowPos4svMESA          f9aglwindowpos4svmesa
#define fglWindowPos4fvMESA            fglwindowpos4fvmesa
#define fglWindowPos4dvMESA            fglwindowpos4dvmesa
#elif FNAME==UNDERSCORE
#define fglWindowPos2iMESA             fglwindowpos2imesa_
#define fglWindowPos2sMESA             fglwindowpos2smesa_
#define f9aglWindowPos2sMESA           f9aglwindowpos2smesa_
#define fglWindowPos2fMESA             fglwindowpos2fmesa_
#define fglWindowPos2dMESA             fglwindowpos2dmesa_
#define fglWindowPos2ivMESA            fglwindowpos2ivmesa_
#define fglWindowPos2svMESA            fglwindowpos2svmesa_
#define f9aglWindowPos2svMESA          f9aglwindowpos2svmesa_
#define fglWindowPos2fvMESA            fglwindowpos2fvmesa_
#define fglWindowPos2dvMESA            fglwindowpos2dvmesa_
#define fglWindowPos3iMESA             fglwindowpos3imesa_
#define fglWindowPos3sMESA             fglwindowpos3smesa_
#define f9aglWindowPos3sMESA           f9aglwindowpos3smesa_
#define fglWindowPos3fMESA             fglwindowpos3fmesa_
#define fglWindowPos3dMESA             fglwindowpos3dmesa_
#define fglWindowPos3ivMESA            fglwindowpos3ivmesa_
#define fglWindowPos3svMESA            fglwindowpos3svmesa_
#define f9aglWindowPos3svMESA          f9aglwindowpos3svmesa_
#define fglWindowPos3fvMESA            fglwindowpos3fvmesa_
#define fglWindowPos3dvMESA            fglwindowpos3dvmesa_
#define fglWindowPos4iMESA             fglwindowpos4imesa_
#define fglWindowPos4sMESA             fglwindowpos4smesa_
#define f9aglWindowPos4sMESA           f9aglwindowpos4smesa_
#define fglWindowPos4fMESA             fglwindowpos4fmesa_
#define fglWindowPos4dMESA             fglwindowpos4dmesa_
#define fglWindowPos4ivMESA            fglwindowpos4ivmesa_
#define fglWindowPos4svMESA            fglwindowpos4svmesa_
#define f9aglWindowPos4svMESA          f9aglwindowpos4svmesa_
#define fglWindowPos4fvMESA            fglwindowpos4fvmesa_
#define fglWindowPos4dvMESA            fglwindowpos4dvmesa_
#elif FNAME==UPPERCASE
#define fglWindowPos2iMESA             FGLWINDOWPOS2IMESA
#define fglWindowPos2sMESA             FGLWINDOWPOS2SMESA
#define f9aglWindowPos2sMESA           F9AGLWINDOWPOS2SMESA
#define fglWindowPos2fMESA             FGLWINDOWPOS2FMESA
#define fglWindowPos2dMESA             FGLWINDOWPOS2DMESA
#define fglWindowPos2ivMESA            FGLWINDOWPOS2IVMESA
#define fglWindowPos2svMESA            FGLWINDOWPOS2SVMESA
#define f9aglWindowPos2svMESA          F9AGLWINDOWPOS2SVMESA
#define fglWindowPos2fvMESA            FGLWINDOWPOS2FVMESA
#define fglWindowPos2dvMESA            FGLWINDOWPOS2DVMESA
#define fglWindowPos3iMESA             FGLWINDOWPOS3IMESA
#define fglWindowPos3sMESA             FGLWINDOWPOS3SMESA
#define f9aglWindowPos3sMESA           F9AGLWINDOWPOS3SMESA
#define fglWindowPos3fMESA             FGLWINDOWPOS3FMESA
#define fglWindowPos3dMESA             FGLWINDOWPOS3DMESA
#define fglWindowPos3ivMESA            FGLWINDOWPOS3IVMESA
#define fglWindowPos3svMESA            FGLWINDOWPOS3SVMESA
#define f9aglWindowPos3svMESA          F9AGLWINDOWPOS3SVMESA
#define fglWindowPos3fvMESA            FGLWINDOWPOS3FVMESA
#define fglWindowPos3dvMESA            FGLWINDOWPOS3DVMESA
#define fglWindowPos4iMESA             FGLWINDOWPOS4IMESA
#define fglWindowPos4sMESA             FGLWINDOWPOS4SMESA
#define f9aglWindowPos4sMESA           F9AGLWINDOWPOS4SMESA
#define fglWindowPos4fMESA             FGLWINDOWPOS4FMESA
#define fglWindowPos4dMESA             FGLWINDOWPOS4DMESA
#define fglWindowPos4ivMESA            FGLWINDOWPOS4IVMESA
#define fglWindowPos4svMESA            FGLWINDOWPOS4SVMESA
#define f9aglWindowPos4svMESA          F9AGLWINDOWPOS4SVMESA
#define fglWindowPos4fvMESA            FGLWINDOWPOS4FVMESA
#define fglWindowPos4dvMESA            FGLWINDOWPOS4DVMESA
#endif
#endif

/* local variables */

/* places to put the result of glGetString */
GLint getstring_length;
GLuint *getstring_intstr;

/* C wrappers */

/* OpenGL 1.0 */

void APIENTRY fglAccum( f90glenum op, f90glfloat value )
{
   glAccum( (GLenum) *op, (GLfloat) *value );
}

void APIENTRY fglAlphaFunc( f90glenum func, f90glclampf ref )
{
   glAlphaFunc( (GLenum) *func, (GLclampf) *ref );
}

void APIENTRY fglBegin( f90glenum mode )
{
   glBegin( (GLenum) *mode );
}

void APIENTRY fglBitmap( f90glsizei width, f90glsizei height,
    f90glfloat xorig, f90glfloat yorig, f90glfloat xmove, f90glfloat ymove,
    f90glubyte ARRAY bitmap )
{
   glBitmap( (GLsizei) *width, (GLsizei) *height, (GLfloat) *xorig,
    (GLfloat) *yorig, (GLfloat) *xmove, (GLfloat) *ymove,
    (GLubyte *) bitmap ARRAYREF);
}

void APIENTRY f9bglBitmap(f90glsizei width,f90glsizei height,f90glfloat xorig,
   f90glfloat yorig,f90glfloat xmove,f90glfloat ymove,f90gluint ARRAY bitmap,
   f90glint size )
{
   GLubyte *byte_bitmap;
   int i;
   byte_bitmap = (GLubyte *) malloc((*size)*sizeof(GLubyte));
   for (i=0;i<*size;i++) byte_bitmap[i]=(GLubyte) bitmap ARRAYREF [i];
   glBitmap( (GLsizei) *width, (GLsizei) *height, (GLfloat) *xorig,
    (GLfloat) *yorig, (GLfloat) *xmove, (GLfloat) *ymove, byte_bitmap);
   free(byte_bitmap);
}

void APIENTRY f9cglBitmap(f90glsizei width,f90glsizei height,f90glfloat xorig,
   f90glfloat yorig,f90glfloat xmove,f90glfloat ymove,f90glushort ARRAY bitmap,
   f90glint size )
{
   GLubyte *byte_bitmap;
   int i;
   byte_bitmap = (GLubyte *) malloc((*size)*sizeof(GLubyte));
   for (i=0;i<*size;i++) byte_bitmap[i]=(GLubyte) bitmap ARRAYREF [i];
   glBitmap( (GLsizei) *width, (GLsizei) *height, (GLfloat) *xorig,
    (GLfloat) *yorig, (GLfloat) *xmove, (GLfloat) *ymove, byte_bitmap);
   free(byte_bitmap);
}

void APIENTRY fglBlendFunc( f90glenum sfactor, f90glenum dfactor )
{
   glBlendFunc( (GLenum) *sfactor, (GLenum) *dfactor );
}

void APIENTRY fglCallList( f90gluint list )
{
   glCallList( (GLuint) *list );
}

void APIENTRY f9y0glCallLists( f90glsizei n, f90glenum type, f90glvoid ARRAY lists)
{

switch (*type) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*n)*sizeof(GLubyte));
   for (i=0;i<(*n);i++) newlists[i]=(GLubyte)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*n)*sizeof(GLbyte));
   for (i=0;i<(*n);i++) newlists[i]=(GLbyte)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*n)*sizeof(GLushort));
   for (i=0;i<(*n);i++) newlists[i]=(GLushort)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*n)*sizeof(GLshort));
   for (i=0;i<(*n);i++) newlists[i]=(GLshort)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_2_BYTES:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((2*(*n))*sizeof(GLubyte));
   for (i=0;i<2*(*n);i++) newlists[i]=(GLubyte)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_3_BYTES:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((3*(*n))*sizeof(GLubyte));
   for (i=0;i<3*(*n);i++) newlists[i]=(GLubyte)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_4_BYTES:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((4*(*n))*sizeof(GLubyte));
   for (i=0;i<4*(*n);i++) newlists[i]=(GLubyte)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) lists ARRAYREF  );
   break;
}
}

void APIENTRY f9y1glCallLists( f90glsizei n, f90glenum type, f90glvoid ARRAY lists )
{
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) lists ARRAYREF);
}

void APIENTRY f9y4glCallLists( f90glsizei n, f90glenum type, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *lists;
   unsigned char *p;
   int i;
   p = (unsigned char *) &lists;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) lists );
}

void APIENTRY f9y5glCallLists( f90glsizei n, f90glenum type, f90glvoid ARRAY lists )
{
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) lists ARRAYREF);
}

void APIENTRY f9y6glCallLists( f90glsizei n, f90glenum type, f90glvoid ARRAY lists)
{

switch (*type) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*n)*sizeof(GLubyte));
   for (i=0;i<(*n);i++) newlists[i]=(GLubyte)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*n)*sizeof(GLbyte));
   for (i=0;i<(*n);i++) newlists[i]=(GLbyte)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_2_BYTES:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((2*(*n))*sizeof(GLubyte));
   for (i=0;i<2*(*n);i++) newlists[i]=(GLubyte)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_3_BYTES:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((3*(*n))*sizeof(GLubyte));
   for (i=0;i<3*(*n);i++) newlists[i]=(GLubyte)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_4_BYTES:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((4*(*n))*sizeof(GLubyte));
   for (i=0;i<4*(*n);i++) newlists[i]=(GLubyte)lists ARRAYREF [i];
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glCallLists( (GLsizei) *n, (GLenum) *type, (GLvoid *) lists  ARRAYREF );
   break;
}
}

void APIENTRY fglClear( f90glbitfield mask )
{
   glClear( (GLbitfield) *mask );
}

void APIENTRY fglClearAccum( f90glfloat red, f90glfloat green,
    f90glfloat blue, f90glfloat alpha )
{
   glClearAccum( (GLfloat) *red, (GLfloat) *green, (GLfloat) *blue,
    (GLfloat) *alpha );
}

void APIENTRY fglClearColor( f90glclampf red, f90glclampf green,
    f90glclampf blue, f90glclampf alpha )
{
   glClearColor( (GLclampf) *red, (GLclampf) *green, (GLclampf) *blue,
    (GLclampf) *alpha );
}

void APIENTRY fglClearDepth( f90glclampd depth )
{
   glClearDepth( (GLclampd) *depth );
}

void APIENTRY fglClearIndex( f90glfloat c )
{
   glClearIndex( (GLfloat) *c );
}

void APIENTRY fglClearStencil( f90glint s )
{
   glClearStencil( (GLint) *s );
}

void APIENTRY fglClipPlane( f90glenum plane, f90gldouble ARRAY equation )
{
   glClipPlane( (GLenum) *plane, (GLdouble *) equation ARRAYREF);
}

void APIENTRY fglColor3b( f90glbyte red, f90glbyte green, f90glbyte blue )
{
   glColor3b( (GLbyte) *red, (GLbyte) *green, (GLbyte) *blue );
}

void APIENTRY f9bglColor3b(f90glint red,f90glint green,f90glint blue )
{
   glColor3b( (GLbyte) *red, (GLbyte) *green, (GLbyte) *blue);
}

void APIENTRY f9cglColor3b(f90glshort red,f90glshort green,f90glshort blue )
{
   glColor3b( (GLbyte) *red, (GLbyte) *green, (GLbyte) *blue);
}

void APIENTRY fglColor3bv( f90glbyte ARRAY v )
{
   glColor3bv( (GLbyte *) v ARRAYREF);
}

void APIENTRY f9bglColor3bv(f90glint ARRAY v )
{
   GLbyte byte_v[3];
   int i;
   for (i=0;i<3;i++) byte_v[i]=(GLbyte) v ARRAYREF [i];
   glColor3bv( byte_v);
}

void APIENTRY f9cglColor3bv(f90glshort ARRAY v )
{
   GLbyte byte_v[3];
   int i;
   for (i=0;i<3;i++) byte_v[i]=(GLbyte) v ARRAYREF [i];
   glColor3bv( byte_v);
}

void APIENTRY fglColor3d( f90gldouble red, f90gldouble green,
    f90gldouble blue )
{
   glColor3d( (GLdouble) *red, (GLdouble) *green, (GLdouble) *blue );
}

void APIENTRY fglColor3dv( f90gldouble ARRAY v )
{
   glColor3dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglColor3f( f90glfloat red, f90glfloat green, f90glfloat blue )
{
   glColor3f( (GLfloat) *red, (GLfloat) *green, (GLfloat) *blue );
}

void APIENTRY fglColor3fv( f90glfloat ARRAY v )
{
   glColor3fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglColor3i( f90glint red, f90glint green, f90glint blue )
{
   glColor3i( (GLint) *red, (GLint) *green, (GLint) *blue );
}

void APIENTRY fglColor3iv( f90glint ARRAY v )
{
   glColor3iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglColor3s( f90glshort red, f90glshort green, f90glshort blue )
{
   glColor3s( (GLshort) *red, (GLshort) *green, (GLshort) *blue );
}

void APIENTRY f9aglColor3s(f90glint red,f90glint green,f90glint blue )
{
   glColor3s( (GLshort) *red, (GLshort) *green, (GLshort) *blue);
}

void APIENTRY fglColor3sv( f90glshort ARRAY v )
{
   glColor3sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglColor3sv(f90glint ARRAY v )
{
   GLshort short_v[3];
   int i;
   for (i=0;i<3;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glColor3sv( short_v);
}

void APIENTRY fglColor3ub( f90glubyte red, f90glubyte green,
    f90glubyte blue )
{
   glColor3ub( (GLubyte) *red, (GLubyte) *green, (GLubyte) *blue );
}

void APIENTRY f9bglColor3ub(f90gluint red,f90gluint green,f90gluint blue )
{
   glColor3ub( (GLubyte) *red, (GLubyte) *green, (GLubyte) *blue);
}

void APIENTRY f9cglColor3ub(f90glushort red,f90glushort green,
   f90glushort blue )
{
   glColor3ub( (GLubyte) *red, (GLubyte) *green, (GLubyte) *blue);
}

void APIENTRY fglColor3ubv( f90glubyte ARRAY v )
{
   glColor3ubv( (GLubyte *) v ARRAYREF);
}

void APIENTRY f9bglColor3ubv(f90gluint ARRAY v )
{
   GLubyte byte_v[3];
   int i;
   for (i=0;i<3;i++) byte_v[i]=(GLubyte) v ARRAYREF [i];
   glColor3ubv( byte_v);
}

void APIENTRY f9cglColor3ubv(f90glushort ARRAY v )
{
   GLubyte byte_v[3];
   int i;
   for (i=0;i<3;i++) byte_v[i]=(GLubyte) v ARRAYREF [i];
   glColor3ubv( byte_v);
}

void APIENTRY fglColor3ui( f90gluint red, f90gluint green, f90gluint blue )
{
   glColor3ui( (GLuint) *red, (GLuint) *green, (GLuint) *blue );
}

void APIENTRY fglColor3uiv( f90gluint ARRAY v )
{
   glColor3uiv( (GLuint *) v ARRAYREF);
}

void APIENTRY fglColor3us( f90glushort red, f90glushort green,
    f90glushort blue )
{
   glColor3us( (GLushort) *red, (GLushort) *green, (GLushort) *blue );
}

void APIENTRY f9aglColor3us(f90gluint red,f90gluint green,f90gluint blue )
{
   glColor3us( (GLushort) *red, (GLushort) *green, (GLushort) *blue);
}

void APIENTRY fglColor3usv( f90glushort ARRAY v )
{
   glColor3usv( (GLushort *) v ARRAYREF);
}

void APIENTRY f9aglColor3usv(f90gluint ARRAY v )
{
   GLushort short_v[3];
   int i;
   for (i=0;i<3;i++) short_v[i]=(GLushort) v ARRAYREF [i];
   glColor3usv( short_v);
}

void APIENTRY fglColor4b( f90glbyte red, f90glbyte green, f90glbyte blue,
    f90glbyte alpha )
{
   glColor4b( (GLbyte) *red, (GLbyte) *green, (GLbyte) *blue,
    (GLbyte) *alpha );
}

void APIENTRY f9bglColor4b(f90glint red,f90glint green,f90glint blue,
   f90glint alpha )
{
   glColor4b( (GLbyte) *red, (GLbyte) *green, (GLbyte) *blue,
    (GLbyte) *alpha);
}

void APIENTRY f9cglColor4b(f90glshort red,f90glshort green,f90glshort blue,
   f90glshort alpha )
{
   glColor4b( (GLbyte) *red, (GLbyte) *green, (GLbyte) *blue,
    (GLbyte) *alpha);
}

void APIENTRY fglColor4bv( f90glbyte ARRAY v )
{
   glColor4bv( (GLbyte *) v ARRAYREF);
}

void APIENTRY f9bglColor4bv(f90glint ARRAY v )
{
   GLbyte byte_v[4];
   int i;
   for (i=0;i<4;i++) byte_v[i]=(GLbyte) v ARRAYREF [i];
   glColor4bv( byte_v);
}

void APIENTRY f9cglColor4bv(f90glshort ARRAY v )
{
   GLbyte byte_v[4];
   int i;
   for (i=0;i<4;i++) byte_v[i]=(GLbyte) v ARRAYREF [i];
   glColor4bv( byte_v);
}

void APIENTRY fglColor4d( f90gldouble red, f90gldouble green,
    f90gldouble blue, f90gldouble alpha )
{
   glColor4d( (GLdouble) *red, (GLdouble) *green, (GLdouble) *blue,
    (GLdouble) *alpha );
}

void APIENTRY fglColor4dv( f90gldouble ARRAY v )
{
   glColor4dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglColor4f( f90glfloat red, f90glfloat green, f90glfloat blue,
    f90glfloat alpha )
{
   glColor4f( (GLfloat) *red, (GLfloat) *green, (GLfloat) *blue,
    (GLfloat) *alpha );
}

void APIENTRY fglColor4fv( f90glfloat ARRAY v )
{
   glColor4fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglColor4i( f90glint red, f90glint green, f90glint blue,
    f90glint alpha )
{
   glColor4i( (GLint) *red, (GLint) *green, (GLint) *blue, (GLint) *alpha );
}

void APIENTRY fglColor4iv( f90glint ARRAY v )
{
   glColor4iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglColor4s( f90glshort red, f90glshort green, f90glshort blue,
    f90glshort alpha )
{
   glColor4s( (GLshort) *red, (GLshort) *green, (GLshort) *blue,
    (GLshort) *alpha );
}

void APIENTRY f9aglColor4s(f90glint red,f90glint green,f90glint blue,
   f90glint alpha )
{
   glColor4s( (GLshort) *red, (GLshort) *green, (GLshort) *blue,
    (GLshort) *alpha);
}

void APIENTRY fglColor4sv( f90glshort ARRAY v )
{
   glColor4sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglColor4sv(f90glint ARRAY v )
{
   GLshort short_v[4];
   int i;
   for (i=0;i<4;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glColor4sv( short_v);
}

void APIENTRY fglColor4ub( f90glubyte red, f90glubyte green, f90glubyte blue,
    f90glubyte alpha )
{
   glColor4ub( (GLubyte) *red, (GLubyte) *green, (GLubyte) *blue,
    (GLubyte) *alpha );
}

void APIENTRY f9bglColor4ub(f90gluint red,f90gluint green,f90gluint blue,
   f90gluint alpha )
{
   glColor4ub( (GLubyte) *red, (GLubyte) *green, (GLubyte) *blue,
    (GLubyte) *alpha);
}

void APIENTRY f9cglColor4ub(f90glushort red,f90glushort green,
   f90glushort blue,f90glushort alpha )
{
   glColor4ub( (GLubyte) *red, (GLubyte) *green, (GLubyte) *blue,
    (GLubyte) *alpha);
}

void APIENTRY fglColor4ubv( f90glubyte ARRAY v )
{
   glColor4ubv( (GLubyte *) v ARRAYREF);
}

void APIENTRY f9bglColor4ubv(f90gluint ARRAY v )
{
   GLubyte byte_v[4];
   int i;
   for (i=0;i<4;i++) byte_v[i]=(GLubyte) v ARRAYREF [i];
   glColor4ubv( byte_v);
}

void APIENTRY f9cglColor4ubv(f90glushort ARRAY v )
{
   GLubyte byte_v[4];
   int i;
   for (i=0;i<4;i++) byte_v[i]=(GLubyte) v ARRAYREF [i];
   glColor4ubv( byte_v);
}

void APIENTRY fglColor4ui( f90gluint red, f90gluint green, f90gluint blue,
    f90gluint alpha )
{
   glColor4ui( (GLuint) *red, (GLuint) *green, (GLuint) *blue,
    (GLuint) *alpha );
}

void APIENTRY fglColor4uiv( f90gluint ARRAY v )
{
   glColor4uiv( (GLuint *) v ARRAYREF);
}

void APIENTRY fglColor4us( f90glushort red, f90glushort green,
    f90glushort blue, f90glushort alpha )
{
   glColor4us( (GLushort) *red, (GLushort) *green, (GLushort) *blue,
    (GLushort) *alpha );
}

void APIENTRY f9aglColor4us(f90gluint red,f90gluint green,f90gluint blue,
   f90gluint alpha )
{
   glColor4us( (GLushort) *red, (GLushort) *green, (GLushort) *blue,
    (GLushort) *alpha);
}

void APIENTRY fglColor4usv( f90glushort ARRAY v )
{
   glColor4usv( (GLushort *) v ARRAYREF);
}

void APIENTRY f9aglColor4usv(f90gluint ARRAY v )
{
   GLushort short_v[4];
   int i;
   for (i=0;i<4;i++) short_v[i]=(GLushort) v ARRAYREF [i];
   glColor4usv( short_v);
}

void APIENTRY fglColorMask( f90glboolean red, f90glboolean green,
    f90glboolean blue, f90glboolean alpha )
{
   glColorMask( (GLboolean) *red, (GLboolean) *green, (GLboolean) *blue,
    (GLboolean) *alpha );
}

void APIENTRY f9eglColorMask(f90glint red,f90glint green,f90glint blue,
   f90glint alpha )
{
   glColorMask( (GLboolean) *red, (GLboolean) *green, (GLboolean) *blue,
    (GLboolean) *alpha);
}

void APIENTRY fglColorMaterial( f90glenum face, f90glenum mode )
{
   glColorMaterial( (GLenum) *face, (GLenum) *mode );
}

void APIENTRY fglCopyPixels( f90glint x, f90glint y, f90glsizei width,
    f90glsizei height, f90glenum type )
{
   glCopyPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
    (GLenum) *type );
}

void APIENTRY fglCullFace( f90glenum mode )
{
   glCullFace( (GLenum) *mode );
}

void APIENTRY fglDeleteLists( f90gluint list, f90glsizei range )
{
   glDeleteLists( (GLuint) *list, (GLsizei) *range );
}

void APIENTRY fglDepthFunc( f90glenum func )
{
   glDepthFunc( (GLenum) *func );
}

void APIENTRY fglDepthMask( f90glboolean flag )
{
   glDepthMask( (GLboolean) *flag );
}

void APIENTRY f9eglDepthMask(f90glint flag )
{
   glDepthMask( (GLboolean) *flag);
}

void APIENTRY fglDepthRange( f90glclampd neer, f90glclampd faar )
{
   glDepthRange( (GLclampd) *neer, (GLclampd) *faar );
}

void APIENTRY fglDisable( f90glenum cap )
{
   glDisable( (GLenum) *cap );
}

void APIENTRY fglDrawBuffer( f90glenum mode )
{
   glDrawBuffer( (GLenum) *mode );
}

void APIENTRY f9y0glDrawPixels( f90glsizei width, f90glsizei height,
   f90glenum format, f90glenum type, f90glvoid ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newpixels;
   int i;
   newpixels = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLubyte)pixels ARRAYREF [i];
   glDrawPixels( (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
      (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newpixels;
   int i;
   newpixels = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLbyte)pixels ARRAYREF [i];
   glDrawPixels( (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
      (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newpixels;
   int i;
   newpixels = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newpixels[i]=(GLushort)pixels ARRAYREF [i];
   glDrawPixels( (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
      (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newpixels;
   int i;
   newpixels = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newpixels[i]=(GLshort)pixels ARRAYREF [i];
   glDrawPixels( (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
      (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 default:
   glDrawPixels( (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
      (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}
}

void APIENTRY f9y1glDrawPixels( f90glsizei width, f90glsizei height,
   f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glDrawPixels( (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
      (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glDrawPixels( f90glsizei width, f90glsizei height,
   f90glenum format, f90glenum type, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glDrawPixels( (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
      (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glDrawPixels( f90glsizei width, f90glsizei height,
   f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glDrawPixels( (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
      (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glDrawPixels( f90glsizei width, f90glsizei height,
   f90glenum format, f90glenum type, f90glshort ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newpixels;
   int i;
   newpixels = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLubyte)pixels ARRAYREF [i];
   glDrawPixels( (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
      (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newpixels;
   int i;
   newpixels = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLbyte)pixels ARRAYREF [i];
   glDrawPixels( (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
      (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 default:
   glDrawPixels( (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
      (GLenum) *type, (GLvoid *) pixels ARRAYREF);
   break;
}
}

void APIENTRY fglEdgeFlag( f90glboolean flag )
{
   glEdgeFlag( (GLboolean) *flag );
}

void APIENTRY f9eglEdgeFlag(f90glint flag )
{
   glEdgeFlag( (GLboolean) *flag);
}

void APIENTRY fglEdgeFlagv( f90glboolean ARRAY flag )
{
   glEdgeFlagv( (GLboolean *) flag ARRAYREF);
}

void APIENTRY f9eglEdgeFlagv(f90glint ARRAY flag, f90glint length )
{
   GLboolean *boolean_flag;
   int i;
   boolean_flag = (GLboolean *) malloc((*length)*sizeof(GLboolean));
   for (i=0;i<(*length);i++) boolean_flag[i]=(GLboolean) flag ARRAYREF [i];
   glEdgeFlagv( boolean_flag);
   free(boolean_flag);
}

void APIENTRY fglEnable( f90glenum cap )
{
   glEnable( (GLenum) *cap );
}

void APIENTRY fglEnd( void )
{
   glEnd();
}

void APIENTRY fglEndList( void )
{
   glEndList();
}

void APIENTRY fglEvalCoord1d( f90gldouble u )
{
   glEvalCoord1d( (GLdouble) *u );
}

void APIENTRY fglEvalCoord1dv( f90gldouble ARRAY u )
{
   glEvalCoord1dv( (GLdouble *) u ARRAYREF);
}

void APIENTRY fglEvalCoord1f( f90glfloat u )
{
   glEvalCoord1f( (GLfloat) *u );
}

void APIENTRY fglEvalCoord1fv( f90glfloat ARRAY u )
{
   glEvalCoord1fv( (GLfloat *) u ARRAYREF);
}

void APIENTRY fglEvalCoord2d( f90gldouble u, f90gldouble v )
{
   glEvalCoord2d( (GLdouble) *u, (GLdouble) *v );
}

void APIENTRY fglEvalCoord2dv( f90gldouble ARRAY u )
{
   glEvalCoord2dv( (GLdouble *) u ARRAYREF);
}

void APIENTRY fglEvalCoord2f( f90glfloat u, f90glfloat v )
{
   glEvalCoord2f( (GLfloat) *u, (GLfloat) *v );
}

void APIENTRY fglEvalCoord2fv( f90glfloat ARRAY u )
{
   glEvalCoord2fv( (GLfloat *) u ARRAYREF);
}

void APIENTRY fglEvalMesh1( f90glenum mode, f90glint i1, f90glint i2 )
{
   glEvalMesh1( (GLenum) *mode, (GLint) *i1, (GLint) *i2 );
}

void APIENTRY fglEvalMesh2( f90glenum mode, f90glint i1, f90glint i2,
    f90glint j1, f90glint j2 )
{
   glEvalMesh2( (GLenum) *mode, (GLint) *i1, (GLint) *i2, (GLint) *j1,
    (GLint) *j2 );
}

void APIENTRY fglEvalPoint1( f90glint i )
{
   glEvalPoint1( (GLint) *i );
}

void APIENTRY fglEvalPoint2( f90glint i, f90glint j )
{
   glEvalPoint2( (GLint) *i, (GLint) *j );
}

void APIENTRY fglFeedbackBuffer( f90glsizei size, f90glenum type,
    f90glfloat ARRAY buffer )
{
   glFeedbackBuffer( (GLsizei) *size, (GLenum) *type, (GLfloat *) buffer ARRAYREF);
}

void APIENTRY fglFinish( void )
{
   glFinish();
}

void APIENTRY fglFlush( void )
{
   glFlush();
}

void APIENTRY fglFogf( f90glenum pname, f90glfloat param )
{
   glFogf( (GLenum) *pname, (GLfloat) *param );
}

void APIENTRY fglFogfv( f90glenum pname, f90glfloat ARRAY params )
{
   glFogfv( (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglFogi( f90glenum pname, f90glint param )
{
   glFogi( (GLenum) *pname, (GLint) *param );
}

void APIENTRY fglFogiv( f90glenum pname, f90glint ARRAY params )
{
   glFogiv( (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY fglFrontFace( f90glenum mode )
{
   glFrontFace( (GLenum) *mode );
}

void APIENTRY fglFrustum( f90gldouble left, f90gldouble right,
    f90gldouble bottom, f90gldouble top, f90gldouble neer, f90gldouble faar )
{
   glFrustum( (GLdouble) *left, (GLdouble) *right, (GLdouble) *bottom,
    (GLdouble) *top, (GLdouble) *neer, (GLdouble) *faar );
}

GLuint APIENTRY fglGenLists( f90glsizei range )
{
   return (GLuint) glGenLists( (GLsizei) *range );
}

void APIENTRY fglGetBooleanv( f90glenum pname, f90glboolean ARRAY params )
{
   glGetBooleanv( (GLenum) *pname, (GLboolean *) params ARRAYREF);
}

void APIENTRY f9eglGetBooleanv(f90glenum pname,f90glint ARRAY params,f90glint length )
{
   GLboolean *boolean_params;
   int i;
   boolean_params = (GLboolean *)malloc((*length)*sizeof(GLboolean));
   glGetBooleanv( (GLenum) *pname, boolean_params);
   for (i=0;i<(*length);i++) params ARRAYREF [i]=(GLint) boolean_params [i];
   free(boolean_params);
}

void APIENTRY fglGetClipPlane( f90glenum plane, f90gldouble ARRAY equation )
{
   glGetClipPlane( (GLenum) *plane, (GLdouble *) equation ARRAYREF);
}

void APIENTRY fglGetDoublev( f90glenum pname, f90gldouble ARRAY params )
{
   glGetDoublev( (GLenum) *pname, (GLdouble *) params ARRAYREF);
}

GLenum APIENTRY fglGetError( void )
{
   return (GLenum) glGetError();
}

void APIENTRY fglGetFloatv( f90glenum pname, f90glfloat ARRAY params )
{
   glGetFloatv( (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglGetIntegerv( f90glenum pname, f90glint ARRAY params )
{
   glGetIntegerv( (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY fglGetLightfv( f90glenum light, f90glenum pname,
    f90glfloat ARRAY params )
{
   glGetLightfv( (GLenum) *light, (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglGetLightiv( f90glenum light, f90glenum pname,
    f90glint ARRAY params )
{
   glGetLightiv( (GLenum) *light, (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY fglGetMapdv( f90glenum target, f90glenum query, f90gldouble ARRAY v )
{
   glGetMapdv( (GLenum) *target, (GLenum) *query, (GLdouble *) v ARRAYREF);
}

void APIENTRY fglGetMapfv( f90glenum target, f90glenum query, f90glfloat ARRAY v )
{
   glGetMapfv( (GLenum) *target, (GLenum) *query, (GLfloat *) v ARRAYREF);
}

void APIENTRY fglGetMapiv( f90glenum target, f90glenum query, f90glint ARRAY v )
{
   glGetMapiv( (GLenum) *target, (GLenum) *query, (GLint *) v ARRAYREF);
}

void APIENTRY fglGetMaterialfv( f90glenum face, f90glenum pname,
    f90glfloat ARRAY params )
{
   glGetMaterialfv( (GLenum) *face, (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglGetMaterialiv( f90glenum face, f90glenum pname,
    f90glint ARRAY params )
{
   glGetMaterialiv( (GLenum) *face, (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY fglGetPixelMapfv( f90glenum map, f90glfloat ARRAY values )
{
   glGetPixelMapfv( (GLenum) *map, (GLfloat *) values ARRAYREF);
}

void APIENTRY fglGetPixelMapuiv( f90glenum map, f90gluint ARRAY values )
{
   glGetPixelMapuiv( (GLenum) *map, (GLuint *) values ARRAYREF);
}

void APIENTRY fglGetPixelMapusv( f90glenum map, f90glushort ARRAY values )
{
   glGetPixelMapusv( (GLenum) *map, (GLushort *) values ARRAYREF);
}

void APIENTRY f9aglGetPixelMapusv(f90glenum map,f90gluint ARRAY values )
{
   GLushort short_values[4];
   int i;
   glGetPixelMapusv( (GLenum) *map, short_values);
   for (i=0;i<4;i++) values ARRAYREF [i]=(GLuint) short_values[i];
}

void APIENTRY fglGetPolygonStipple( f90glubyte ARRAY mask , f90glint size)
{
   glGetPolygonStipple( (GLubyte *) mask ARRAYREF);
}

void APIENTRY f9bglGetPolygonStipple(f90gluint ARRAY mask , f90glint size)
{
   GLubyte *byte_mask;
   int i;
   byte_mask = (GLubyte *)malloc((*size)*sizeof(GLubyte));
   glGetPolygonStipple( byte_mask);
   for (i=0;i<*size;i++) mask ARRAYREF [i]=(GLuint) byte_mask[i];
   free(byte_mask);
}

void APIENTRY f9cglGetPolygonStipple(f90glushort ARRAY mask , f90glint size)
{
   GLubyte *byte_mask;
   int i;
   byte_mask = (GLubyte *)malloc((*size)*sizeof(GLubyte));
   glGetPolygonStipple( byte_mask);
   for (i=0;i<*size;i++) mask ARRAYREF [i]=(GLushort) byte_mask[i];
   free(byte_mask);
}

GLint APIENTRY f9yglGetStringLen(f90glenum name)
{
   int i;
   char *str;
   str = (char *) glGetString( (GLenum) *name);
   if (str) {
      getstring_length = strlen(str);
      getstring_intstr = (GLuint *)malloc((getstring_length)*sizeof(GLuint));
      for (i=0;i<getstring_length;i++) getstring_intstr[i] = (GLuint) str[i];
   }else{
      getstring_length = 0;
      getstring_intstr = (GLuint *)malloc(sizeof(GLuint));
   }
   return (GLint) getstring_length;
}

void APIENTRY f9yglGetString(f90gluint ARRAY int_str)
{
   int i;
   for (i=0;i<getstring_length;i++) int_str ARRAYREF [i] = getstring_intstr[i];
}

GLuint APIENTRY f9yglGetString1(f90glint pos)
{
   return (GLuint) getstring_intstr[(*pos)-1];
}

#ifdef PURE_FUNCTIONS
GLint APIENTRY f9yglGetStringFree( void )
#else
void APIENTRY f9yglGetStringFree( void )
#endif
{
   getstring_length = 0;
   free(getstring_intstr);
#ifdef PURE_FUNCTIONS
   return 1;
#endif
}

void APIENTRY fglGetTexEnvfv( f90glenum target, f90glenum pname,
    f90glfloat ARRAY params )
{
   glGetTexEnvfv( (GLenum) *target, (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglGetTexEnviv( f90glenum target, f90glenum pname,
    f90glint ARRAY params )
{
   glGetTexEnviv( (GLenum) *target, (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY fglGetTexGendv( f90glenum coord, f90glenum pname,
    f90gldouble ARRAY params )
{
   glGetTexGendv( (GLenum) *coord, (GLenum) *pname, (GLdouble *) params ARRAYREF);
}

void APIENTRY fglGetTexGenfv( f90glenum coord, f90glenum pname,
    f90glfloat ARRAY params )
{
   glGetTexGenfv( (GLenum) *coord, (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglGetTexGeniv( f90glenum coord, f90glenum pname,
    f90glint ARRAY params )
{
   glGetTexGeniv( (GLenum) *coord, (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY f9y0glGetTexImage( f90glenum target, f90glint level,
   f90glenum format, f90glenum type, f90glvoid ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newpixels;
   int i;
   newpixels = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   glGetTexImage( (GLenum) *target, (GLint) *level,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLint)newpixels[i];
   free(newpixels);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newpixels;
   int i;
   newpixels = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   glGetTexImage( (GLenum) *target, (GLint) *level,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLint)newpixels[i];
   free(newpixels);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newpixels;
   int i;
   newpixels = (GLushort *) malloc((*length)*sizeof(GLushort));
   glGetTexImage( (GLenum) *target, (GLint) *level,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLint)newpixels[i];
   free(newpixels);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newpixels;
   int i;
   newpixels = (GLshort *) malloc((*length)*sizeof(GLshort));
   glGetTexImage( (GLenum) *target, (GLint) *level,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLint)newpixels[i];
   free(newpixels);
   break;
  }
 default:
   glGetTexImage( (GLenum) *target, (GLint) *level,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
   break;
}
}

void APIENTRY f9y1glGetTexImage( f90glenum target, f90glint level,
   f90glenum format, f90glenum type, f90glvoid ARRAY pixels)
{
   glGetTexImage( (GLenum) *target, (GLint) *level,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glGetTexImage( f90glenum target, f90glint level,
   f90glenum format, f90glenum type, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glGetTexImage( (GLenum) *target, (GLint) *level,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glGetTexImage( f90glenum target, f90glint level,
   f90glenum format, f90glenum type, f90glvoid ARRAY pixels)
{
   glGetTexImage( (GLenum) *target, (GLint) *level,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glGetTexImage( f90glenum target, f90glint level,
   f90glenum format, f90glenum type, f90glshort ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newpixels;
   int i;
   newpixels = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   glGetTexImage( (GLenum) *target, (GLint) *level,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLshort)newpixels[i];
   free(newpixels);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newpixels;
   int i;
   newpixels = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   glGetTexImage( (GLenum) *target, (GLint) *level,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLshort)newpixels[i];
   free(newpixels);
   break;
  }
 default:
   glGetTexImage( (GLenum) *target, (GLint) *level,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}
}

void APIENTRY fglGetTexLevelParameterfv( f90glenum target, f90glint level,
    f90glenum pname, f90glfloat ARRAY params )
{
   glGetTexLevelParameterfv( (GLenum) *target, (GLint) *level,
    (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglGetTexLevelParameteriv( f90glenum target, f90glint level,
    f90glenum pname, f90glint ARRAY params )
{
   glGetTexLevelParameteriv( (GLenum) *target, (GLint) *level,
    (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY fglGetTexParameterfv( f90glenum target, f90glenum pname,
    f90glfloat ARRAY params )
{
   glGetTexParameterfv( (GLenum) *target, (GLenum) *pname,
    (GLfloat *) params ARRAYREF);
}

void APIENTRY fglGetTexParameteriv( f90glenum target, f90glenum pname,
    f90glint ARRAY params )
{
   glGetTexParameteriv( (GLenum) *target, (GLenum) *pname,
    (GLint *) params ARRAYREF);
}

void APIENTRY fglHint( f90glenum target, f90glenum mode )
{
   glHint( (GLenum) *target, (GLenum) *mode );
}

void APIENTRY fglIndexMask( f90gluint mask )
{
   glIndexMask( (GLuint) *mask );
}

void APIENTRY fglIndexd( f90gldouble c )
{
   glIndexd( (GLdouble) *c );
}

void APIENTRY fglIndexdv( f90gldouble ARRAY c )
{
   glIndexdv( (GLdouble *) c ARRAYREF);
}

void APIENTRY fglIndexf( f90glfloat c )
{
   glIndexf( (GLfloat) *c );
}

void APIENTRY fglIndexfv( f90glfloat ARRAY c )
{
   glIndexfv( (GLfloat *) c ARRAYREF);
}

void APIENTRY fglIndexi( f90glint c )
{
   glIndexi( (GLint) *c );
}

void APIENTRY fglIndexiv( f90glint ARRAY c )
{
   glIndexiv( (GLint *) c ARRAYREF);
}

void APIENTRY fglIndexs( f90glshort c )
{
   glIndexs( (GLshort) *c );
}

void APIENTRY f9aglIndexs(f90glint c )
{
   glIndexs( (GLshort) *c);
}

void APIENTRY fglIndexsv( f90glshort ARRAY c )
{
   glIndexsv( (GLshort *) c ARRAYREF);
}

void APIENTRY f9aglIndexsv(f90glint ARRAY c )
{
   GLshort short_c[1];
   int i;
   for (i=0;i<1;i++) short_c[i]=(GLshort) c ARRAYREF [i];
   glIndexsv( short_c);
}

void APIENTRY fglInitNames( void )
{
   glInitNames();
}

GLboolean APIENTRY fglIsEnabled( f90glenum cap )
{
   return (GLboolean) glIsEnabled( (GLenum) *cap );
}

GLint APIENTRY f9eglIsEnabled(f90gluint cap)
{
   GLboolean temp;
   temp = glIsEnabled( (GLenum) *cap);
   return (GLint) temp;
}

GLboolean APIENTRY fglIsList( f90gluint list )
{
   return (GLboolean) glIsList( (GLuint) *list );
}

GLint APIENTRY f9eglIsList(f90gluint list)
{
   GLboolean temp;
   temp = glIsList( (GLuint) *list);
   return (GLint) temp;
}

void APIENTRY fglLightModelf( f90glenum pname, f90glfloat param )
{
   glLightModelf( (GLenum) *pname, (GLfloat) *param );
}

void APIENTRY fglLightModelfv( f90glenum pname, f90glfloat ARRAY params )
{
   glLightModelfv( (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglLightModeli( f90glenum pname, f90glint param )
{
   glLightModeli( (GLenum) *pname, (GLint) *param );
}

void APIENTRY fglLightModeliv( f90glenum pname, f90glint ARRAY params )
{
   glLightModeliv( (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY fglLightf( f90glenum light, f90glenum pname, f90glfloat param )
{
   glLightf( (GLenum) *light, (GLenum) *pname, (GLfloat) *param );
}

void APIENTRY fglLightfv( f90glenum light, f90glenum pname,
    f90glfloat ARRAY params )
{
   glLightfv( (GLenum) *light, (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglLighti( f90glenum light, f90glenum pname, f90glint param )
{
   glLighti( (GLenum) *light, (GLenum) *pname, (GLint) *param );
}

void APIENTRY fglLightiv( f90glenum light, f90glenum pname, f90glint ARRAY params )
{
   glLightiv( (GLenum) *light, (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY fglLineStipple( f90glint factor, f90glushort pattern )
{
   glLineStipple( (GLint) *factor, (GLushort) *pattern );
}

void APIENTRY f9aglLineStipple(f90glint factor,f90gluint pattern )
{
   glLineStipple( (GLint) *factor, (GLushort) *pattern);
}

void APIENTRY fglLineWidth( f90glfloat width )
{
   glLineWidth( (GLfloat) *width );
}

void APIENTRY fglListBase( f90gluint base )
{
   glListBase( (GLuint) *base );
}

void APIENTRY fglLoadIdentity( void )
{
   glLoadIdentity();
}

void APIENTRY fglLoadMatrixd( f90gldouble ARRAY m )
{
   glLoadMatrixd( (GLdouble *) m ARRAYREF);
}

void APIENTRY fglLoadMatrixf( f90glfloat ARRAY m )
{
   glLoadMatrixf( (GLfloat *) m ARRAYREF);
}

void APIENTRY fglLoadName( f90gluint name )
{
   glLoadName( (GLuint) *name );
}

void APIENTRY fglLogicOp( f90glenum opcode )
{
   glLogicOp( (GLenum) *opcode );
}

void APIENTRY fglMap1d( f90glenum target, f90gldouble u1, f90gldouble u2,
    f90glint stride, f90glint order, f90gldouble ARRAY points )
{
   glMap1d( (GLenum) *target, (GLdouble) *u1, (GLdouble) *u2, (GLint) *stride,
    (GLint) *order, (GLdouble *) points ARRAYREF);
}

void APIENTRY fglMap1f( f90glenum target, f90glfloat u1, f90glfloat u2,
    f90glint stride, f90glint order, f90glfloat ARRAY points )
{
   glMap1f( (GLenum) *target, (GLfloat) *u1, (GLfloat) *u2, (GLint) *stride,
    (GLint) *order, (GLfloat *) points ARRAYREF);
}

void APIENTRY fglMap2d( f90glenum target, f90gldouble u1, f90gldouble u2,
    f90glint ustride, f90glint uorder, f90gldouble v1, f90gldouble v2,
    f90glint vstride, f90glint vorder, f90gldouble ARRAY points )
{
   glMap2d( (GLenum) *target, (GLdouble) *u1, (GLdouble) *u2,
    (GLint) *ustride, (GLint) *uorder, (GLdouble) *v1, (GLdouble) *v2,
    (GLint) *vstride, (GLint) *vorder, (GLdouble *) points ARRAYREF);
}

void APIENTRY fglMap2f( f90glenum target, f90glfloat u1, f90glfloat u2,
    f90glint ustride, f90glint uorder, f90glfloat v1, f90glfloat v2,
    f90glint vstride, f90glint vorder, f90glfloat ARRAY points )
{
   glMap2f( (GLenum) *target, (GLfloat) *u1, (GLfloat) *u2, (GLint) *ustride,
    (GLint) *uorder, (GLfloat) *v1, (GLfloat) *v2, (GLint) *vstride,
    (GLint) *vorder, (GLfloat *) points ARRAYREF);
}

void APIENTRY fglMapGrid1d( f90glint un, f90gldouble u1, f90gldouble u2 )
{
   glMapGrid1d( (GLint) *un, (GLdouble) *u1, (GLdouble) *u2 );
}

void APIENTRY fglMapGrid1f( f90glint un, f90glfloat u1, f90glfloat u2 )
{
   glMapGrid1f( (GLint) *un, (GLfloat) *u1, (GLfloat) *u2 );
}

void APIENTRY fglMapGrid2d( f90glint un, f90gldouble u1, f90gldouble u2,
    f90glint vn, f90gldouble v1, f90gldouble v2 )
{
   glMapGrid2d( (GLint) *un, (GLdouble) *u1, (GLdouble) *u2, (GLint) *vn,
    (GLdouble) *v1, (GLdouble) *v2 );
}

void APIENTRY fglMapGrid2f( f90glint un, f90glfloat u1, f90glfloat u2,
    f90glint vn, f90glfloat v1, f90glfloat v2 )
{
   glMapGrid2f( (GLint) *un, (GLfloat) *u1, (GLfloat) *u2, (GLint) *vn,
    (GLfloat) *v1, (GLfloat) *v2 );
}

void APIENTRY fglMaterialf( f90glenum face, f90glenum pname,
    f90glfloat param )
{
   glMaterialf( (GLenum) *face, (GLenum) *pname, (GLfloat) *param );
}

void APIENTRY fglMaterialfv( f90glenum face, f90glenum pname,
    f90glfloat ARRAY params )
{
   glMaterialfv( (GLenum) *face, (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglMateriali( f90glenum face, f90glenum pname, f90glint param )
{
   glMateriali( (GLenum) *face, (GLenum) *pname, (GLint) *param );
}

void APIENTRY fglMaterialiv( f90glenum face, f90glenum pname,
    f90glint ARRAY params )
{
   glMaterialiv( (GLenum) *face, (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY fglMatrixMode( f90glenum mode )
{
   glMatrixMode( (GLenum) *mode );
}

void APIENTRY fglMultMatrixd( f90gldouble ARRAY m )
{
   glMultMatrixd( (GLdouble *) m ARRAYREF);
}

void APIENTRY fglMultMatrixf( f90glfloat ARRAY m )
{
   glMultMatrixf( (GLfloat *) m ARRAYREF);
}

void APIENTRY fglNewList( f90gluint list, f90glenum mode )
{
   glNewList( (GLuint) *list, (GLenum) *mode );
}

void APIENTRY fglNormal3b( f90glbyte nx, f90glbyte ny, f90glbyte nz )
{
   glNormal3b( (GLbyte) *nx, (GLbyte) *ny, (GLbyte) *nz );
}

void APIENTRY f9bglNormal3b(f90glint nx,f90glint ny,f90glint nz )
{
   glNormal3b( (GLbyte) *nx, (GLbyte) *ny, (GLbyte) *nz);
}

void APIENTRY f9cglNormal3b(f90glshort nx,f90glshort ny,f90glshort nz )
{
   glNormal3b( (GLbyte) *nx, (GLbyte) *ny, (GLbyte) *nz);
}

void APIENTRY fglNormal3bv( f90glbyte ARRAY v )
{
   glNormal3bv( (GLbyte *) v ARRAYREF);
}

void APIENTRY f9bglNormal3bv(f90glint ARRAY v )
{
   GLbyte byte_v[3];
   int i;
   for (i=0;i<3;i++) byte_v[i]=(GLbyte) v ARRAYREF [i];
   glNormal3bv( byte_v);
}

void APIENTRY f9cglNormal3bv(f90glshort ARRAY v )
{
   GLbyte byte_v[3];
   int i;
   for (i=0;i<3;i++) byte_v[i]=(GLbyte) v ARRAYREF [i];
   glNormal3bv( byte_v);
}

void APIENTRY fglNormal3d( f90gldouble nx, f90gldouble ny, f90gldouble nz )
{
   glNormal3d( (GLdouble) *nx, (GLdouble) *ny, (GLdouble) *nz );
}

void APIENTRY fglNormal3dv( f90gldouble ARRAY v )
{
   glNormal3dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglNormal3f( f90glfloat nx, f90glfloat ny, f90glfloat nz )
{
   glNormal3f( (GLfloat) *nx, (GLfloat) *ny, (GLfloat) *nz );
}

void APIENTRY fglNormal3fv( f90glfloat ARRAY v )
{
   glNormal3fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglNormal3i( f90glint nx, f90glint ny, f90glint nz )
{
   glNormal3i( (GLint) *nx, (GLint) *ny, (GLint) *nz );
}

void APIENTRY fglNormal3iv( f90glint ARRAY v )
{
   glNormal3iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglNormal3s( f90glshort nx, f90glshort ny, f90glshort nz )
{
   glNormal3s( (GLshort) *nx, (GLshort) *ny, (GLshort) *nz );
}

void APIENTRY f9aglNormal3s(f90glint nx,f90glint ny,f90glint nz )
{
   glNormal3s( (GLshort) *nx, (GLshort) *ny, (GLshort) *nz);
}

void APIENTRY fglNormal3sv( f90glshort ARRAY v )
{
   glNormal3sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglNormal3sv(f90glint ARRAY v )
{
   GLshort short_v[3];
   int i;
   for (i=0;i<3;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glNormal3sv( short_v);
}

void APIENTRY fglOrtho( f90gldouble left, f90gldouble right,
    f90gldouble bottom, f90gldouble top, f90gldouble neer, f90gldouble faar )
{
   glOrtho( (GLdouble) *left, (GLdouble) *right, (GLdouble) *bottom,
    (GLdouble) *top, (GLdouble) *neer, (GLdouble) *faar );
}

void APIENTRY fglPassThrough( f90glfloat token )
{
   glPassThrough( (GLfloat) *token );
}

void APIENTRY fglPixelMapfv( f90glenum map, f90glint mapsize,
    f90glfloat ARRAY values )
{
   glPixelMapfv( (GLenum) *map, (GLint) *mapsize, (GLfloat *) values ARRAYREF);
}

void APIENTRY fglPixelMapuiv( f90glenum map, f90glint mapsize,
    f90gluint ARRAY values )
{
   glPixelMapuiv( (GLenum) *map, (GLint) *mapsize, (GLuint *) values ARRAYREF);
}

void APIENTRY fglPixelMapusv( f90glenum map, f90glint mapsize,
    f90glushort ARRAY values )
{
   glPixelMapusv( (GLenum) *map, (GLint) *mapsize, (GLushort *) values ARRAYREF);
}

void APIENTRY f9aglPixelMapusv(f90glenum map,f90glint mapsize,
   f90gluint ARRAY values )
{
   GLushort *short_values;
   int i;
   short_values = (GLushort *)malloc((*mapsize)*sizeof(GLushort));
   for (i=0;i<(*mapsize);i++) short_values[i]=(GLushort) values ARRAYREF [i];
   glPixelMapusv( (GLenum) *map, (GLint) *mapsize, short_values);
   free(short_values);
}

void APIENTRY fglPixelStoref( f90glenum pname, f90glfloat param )
{
   glPixelStoref( (GLenum) *pname, (GLfloat) *param );
}

void APIENTRY fglPixelStorei( f90glenum pname, f90glint param )
{
   glPixelStorei( (GLenum) *pname, (GLint) *param );
}

void APIENTRY fglPixelTransferf( f90glenum pname, f90glfloat param )
{
   glPixelTransferf( (GLenum) *pname, (GLfloat) *param );
}

void APIENTRY fglPixelTransferi( f90glenum pname, f90glint param )
{
   glPixelTransferi( (GLenum) *pname, (GLint) *param );
}

void APIENTRY fglPixelZoom( f90glfloat xfactor, f90glfloat yfactor )
{
   glPixelZoom( (GLfloat) *xfactor, (GLfloat) *yfactor );
}

void APIENTRY fglPointSize( f90glfloat size )
{
   glPointSize( (GLfloat) *size );
}

void APIENTRY fglPolygonMode( f90glenum face, f90glenum mode )
{
   glPolygonMode( (GLenum) *face, (GLenum) *mode );
}

void APIENTRY fglPolygonStipple( f90glubyte ARRAY mask )
{
   glPolygonStipple( (GLubyte *) mask ARRAYREF);
}

void APIENTRY f9bglPolygonStipple(f90gluint ARRAY mask )
{
   GLubyte byte_mask[128];
   int i;
   for (i=0;i<128;i++) byte_mask[i]=(GLubyte) mask ARRAYREF [i];
   glPolygonStipple( byte_mask);
}

void APIENTRY f9cglPolygonStipple(f90glushort ARRAY mask )
{
   GLubyte byte_mask[128];
   int i;
   for (i=0;i<128;i++) byte_mask[i]=(GLubyte) mask ARRAYREF [i];
   glPolygonStipple( byte_mask);
}

void APIENTRY fglPopAttrib( void )
{
   glPopAttrib();
}

void APIENTRY fglPopMatrix( void )
{
   glPopMatrix();
}

void APIENTRY fglPopName( void )
{
   glPopName();
}

void APIENTRY fglPushAttrib( f90glbitfield mask )
{
   glPushAttrib( (GLbitfield) *mask );
}

void APIENTRY fglPushMatrix( void )
{
   glPushMatrix();
}

void APIENTRY fglPushName( f90gluint name )
{
   glPushName( (GLuint) *name );
}

void APIENTRY fglRasterPos2d( f90gldouble x, f90gldouble y )
{
   glRasterPos2d( (GLdouble) *x, (GLdouble) *y );
}

void APIENTRY fglRasterPos2dv( f90gldouble ARRAY v )
{
   glRasterPos2dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglRasterPos2f( f90glfloat x, f90glfloat y )
{
   glRasterPos2f( (GLfloat) *x, (GLfloat) *y );
}

void APIENTRY fglRasterPos2fv( f90glfloat ARRAY v )
{
   glRasterPos2fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglRasterPos2i( f90glint x, f90glint y )
{
   glRasterPos2i( (GLint) *x, (GLint) *y );
}

void APIENTRY fglRasterPos2iv( f90glint ARRAY v )
{
   glRasterPos2iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglRasterPos2s( f90glshort x, f90glshort y )
{
   glRasterPos2s( (GLshort) *x, (GLshort) *y );
}

void APIENTRY f9aglRasterPos2s(f90glint x,f90glint y )
{
   glRasterPos2s( (GLshort) *x, (GLshort) *y);
}

void APIENTRY fglRasterPos2sv( f90glshort ARRAY v )
{
   glRasterPos2sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglRasterPos2sv(f90glint ARRAY v )
{
   GLshort short_v[2];
   int i;
   for (i=0;i<2;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glRasterPos2sv( short_v);
}

void APIENTRY fglRasterPos3d( f90gldouble x, f90gldouble y, f90gldouble z )
{
   glRasterPos3d( (GLdouble) *x, (GLdouble) *y, (GLdouble) *z );
}

void APIENTRY fglRasterPos3dv( f90gldouble ARRAY v )
{
   glRasterPos3dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglRasterPos3f( f90glfloat x, f90glfloat y, f90glfloat z )
{
   glRasterPos3f( (GLfloat) *x, (GLfloat) *y, (GLfloat) *z );
}

void APIENTRY fglRasterPos3fv( f90glfloat ARRAY v )
{
   glRasterPos3fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglRasterPos3i( f90glint x, f90glint y, f90glint z )
{
   glRasterPos3i( (GLint) *x, (GLint) *y, (GLint) *z );
}

void APIENTRY fglRasterPos3iv( f90glint ARRAY v )
{
   glRasterPos3iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglRasterPos3s( f90glshort x, f90glshort y, f90glshort z )
{
   glRasterPos3s( (GLshort) *x, (GLshort) *y, (GLshort) *z );
}

void APIENTRY f9aglRasterPos3s(f90glint x,f90glint y,f90glint z )
{
   glRasterPos3s( (GLshort) *x, (GLshort) *y, (GLshort) *z);
}

void APIENTRY fglRasterPos3sv( f90glshort ARRAY v )
{
   glRasterPos3sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglRasterPos3sv(f90glint ARRAY v )
{
   GLshort short_v[3];
   int i;
   for (i=0;i<3;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glRasterPos3sv( short_v);
}

void APIENTRY fglRasterPos4d( f90gldouble x, f90gldouble y, f90gldouble z,
    f90gldouble w )
{
   glRasterPos4d( (GLdouble) *x, (GLdouble) *y, (GLdouble) *z,
    (GLdouble) *w );
}

void APIENTRY fglRasterPos4dv( f90gldouble ARRAY v )
{
   glRasterPos4dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglRasterPos4f( f90glfloat x, f90glfloat y, f90glfloat z,
    f90glfloat w )
{
   glRasterPos4f( (GLfloat) *x, (GLfloat) *y, (GLfloat) *z, (GLfloat) *w );
}

void APIENTRY fglRasterPos4fv( f90glfloat ARRAY v )
{
   glRasterPos4fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglRasterPos4i( f90glint x, f90glint y, f90glint z,
    f90glint w )
{
   glRasterPos4i( (GLint) *x, (GLint) *y, (GLint) *z, (GLint) *w );
}

void APIENTRY fglRasterPos4iv( f90glint ARRAY v )
{
   glRasterPos4iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglRasterPos4s( f90glshort x, f90glshort y, f90glshort z,
    f90glshort w )
{
   glRasterPos4s( (GLshort) *x, (GLshort) *y, (GLshort) *z, (GLshort) *w );
}

void APIENTRY f9aglRasterPos4s(f90glint x,f90glint y,f90glint z,f90glint w )
{
   glRasterPos4s( (GLshort) *x, (GLshort) *y, (GLshort) *z, (GLshort) *w);
}

void APIENTRY fglRasterPos4sv( f90glshort ARRAY v )
{
   glRasterPos4sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglRasterPos4sv(f90glint ARRAY v )
{
   GLshort short_v[4];
   int i;
   for (i=0;i<4;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glRasterPos4sv( short_v);
}

void APIENTRY fglReadBuffer( f90glenum mode )
{
   glReadBuffer( (GLenum) *mode );
}

void APIENTRY f9y0glReadPixels( f90glint x, f90glint y, f90glsizei width,
   f90glsizei height,
   f90glenum format, f90glenum type, f90glvoid ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newpixels;
   int i;
   newpixels = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   glReadPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLint)newpixels[i];
   free(newpixels);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newpixels;
   int i;
   newpixels = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   glReadPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLint)newpixels[i];
   free(newpixels);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newpixels;
   int i;
   newpixels = (GLushort *) malloc((*length)*sizeof(GLushort));
   glReadPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLint)newpixels[i];
   free(newpixels);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newpixels;
   int i;
   newpixels = (GLshort *) malloc((*length)*sizeof(GLshort));
   glReadPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLint)newpixels[i];
   free(newpixels);
   break;
  }
 default:
   glReadPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}
}

void APIENTRY f9y1glReadPixels( f90glint x, f90glint y, f90glsizei width,
   f90glsizei height,
   f90glenum format, f90glenum type, f90glvoid ARRAY pixels)
{
   glReadPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glReadPixels( f90glint x, f90glint y, f90glsizei width,
   f90glsizei height,
   f90glenum format, f90glenum type, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glReadPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glReadPixels( f90glint x, f90glint y, f90glsizei width,
   f90glsizei height,
   f90glenum format, f90glenum type, f90glvoid ARRAY pixels)
{
   glReadPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glReadPixels( f90glint x, f90glint y, f90glsizei width,
   f90glsizei height,
   f90glenum format, f90glenum type, f90glshort ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newpixels;
   int i;
   newpixels = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   glReadPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLshort)newpixels[i];
   free(newpixels);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newpixels;
   int i;
   newpixels = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   glReadPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   for (i=0;i<(*length);i++) pixels ARRAYREF [i]=(GLshort)newpixels[i];
   free(newpixels);
   break;
  }
 default:
   glReadPixels( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}
}

void APIENTRY fglRectd( f90gldouble x1, f90gldouble y1, f90gldouble x2,
    f90gldouble y2 )
{
   glRectd( (GLdouble) *x1, (GLdouble) *y1, (GLdouble) *x2, (GLdouble) *y2 );
}

void APIENTRY fglRectdv( f90gldouble ARRAY v1, f90gldouble ARRAY v2 )
{
   glRectdv( (GLdouble *) v1 ARRAYREF, (GLdouble *) v2 ARRAYREF );
}

void APIENTRY fglRectf( f90glfloat x1, f90glfloat y1, f90glfloat x2,
    f90glfloat y2 )
{
   glRectf( (GLfloat) *x1, (GLfloat) *y1, (GLfloat) *x2, (GLfloat) *y2 );
}

void APIENTRY fglRectfv( f90glfloat ARRAY v1, f90glfloat ARRAY v2 )
{
   glRectfv( (GLfloat *) v1 ARRAYREF, (GLfloat *) v2 ARRAYREF );
}

void APIENTRY fglRecti( f90glint x1, f90glint y1, f90glint x2, f90glint y2 )
{
   glRecti( (GLint) *x1, (GLint) *y1, (GLint) *x2, (GLint) *y2 );
}

void APIENTRY fglRectiv( f90glint ARRAY v1, f90glint ARRAY v2 )
{
   glRectiv( (GLint *) v1 ARRAYREF, (GLint *) v2 ARRAYREF );
}

void APIENTRY fglRects( f90glshort x1, f90glshort y1, f90glshort x2,
    f90glshort y2 )
{
   glRects( (GLshort) *x1, (GLshort) *y1, (GLshort) *x2, (GLshort) *y2 );
}

void APIENTRY f9aglRects(f90glint x1,f90glint y1,f90glint x2,f90glint y2 )
{
   glRects( (GLshort) *x1, (GLshort) *y1, (GLshort) *x2, (GLshort) *y2);
}

void APIENTRY fglRectsv( f90glshort ARRAY v1, f90glshort ARRAY v2 )
{
   glRectsv( (GLshort *) v1 ARRAYREF, (GLshort *) v2 ARRAYREF );
}

void APIENTRY f9aglRectsv(f90glint ARRAY v1,f90glint ARRAY v2 )
{
   GLshort short_v1[2];
   GLshort short_v2[2];
   int i;
   for (i=0;i<2;i++) short_v1[i]=(GLshort) v1 ARRAYREF [i];
   for (i=0;i<2;i++) short_v2[i]=(GLshort) v2 ARRAYREF [i];
   glRectsv( short_v1, short_v2);
}

GLint APIENTRY fglRenderMode( f90glenum mode )
{
   return (GLint) glRenderMode( (GLenum) *mode );
}

void APIENTRY fglRotated( f90gldouble angle, f90gldouble x, f90gldouble y,
    f90gldouble z )
{
   glRotated( (GLdouble) *angle, (GLdouble) *x, (GLdouble) *y,
    (GLdouble) *z );
}

void APIENTRY fglRotatef( f90glfloat angle, f90glfloat x, f90glfloat y,
    f90glfloat z )
{
   glRotatef( (GLfloat) *angle, (GLfloat) *x, (GLfloat) *y, (GLfloat) *z );
}

void APIENTRY fglScaled( f90gldouble x, f90gldouble y, f90gldouble z )
{
   glScaled( (GLdouble) *x, (GLdouble) *y, (GLdouble) *z );
}

void APIENTRY fglScalef( f90glfloat x, f90glfloat y, f90glfloat z )
{
   glScalef( (GLfloat) *x, (GLfloat) *y, (GLfloat) *z );
}

void APIENTRY fglScissor( f90glint x, f90glint y, f90glsizei width,
    f90glsizei height )
{
   glScissor( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height );
}

void APIENTRY fglSelectBuffer( f90glsizei size, f90gluint ARRAY buffer )
{
   glSelectBuffer( (GLsizei) *size, (GLuint *) buffer ARRAYREF);
}

void APIENTRY fglShadeModel( f90glenum mode )
{
   glShadeModel( (GLenum) *mode );
}

void APIENTRY fglStencilFunc( f90glenum func, f90glint ref, f90gluint mask )
{
   glStencilFunc( (GLenum) *func, (GLint) *ref, (GLuint) *mask );
}

void APIENTRY fglStencilMask( f90gluint mask )
{
   glStencilMask( (GLuint) *mask );
}

void APIENTRY fglStencilOp( f90glenum fail, f90glenum zfail,
    f90glenum zpass )
{
   glStencilOp( (GLenum) *fail, (GLenum) *zfail, (GLenum) *zpass );
}

void APIENTRY fglTexCoord1d( f90gldouble s )
{
   glTexCoord1d( (GLdouble) *s );
}

void APIENTRY fglTexCoord1dv( f90gldouble ARRAY v )
{
   glTexCoord1dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglTexCoord1f( f90glfloat s )
{
   glTexCoord1f( (GLfloat) *s );
}

void APIENTRY fglTexCoord1fv( f90glfloat ARRAY v )
{
   glTexCoord1fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglTexCoord1i( f90glint s )
{
   glTexCoord1i( (GLint) *s );
}

void APIENTRY fglTexCoord1iv( f90glint ARRAY v )
{
   glTexCoord1iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglTexCoord1s( f90glshort s )
{
   glTexCoord1s( (GLshort) *s );
}

void APIENTRY f9aglTexCoord1s(f90glint s )
{
   glTexCoord1s( (GLshort) *s);
}

void APIENTRY fglTexCoord1sv( f90glshort ARRAY v )
{
   glTexCoord1sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglTexCoord1sv(f90glint ARRAY v )
{
   GLshort short_v[1];
   int i;
   for (i=0;i<1;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glTexCoord1sv( short_v);
}

void APIENTRY fglTexCoord2d( f90gldouble s, f90gldouble t )
{
   glTexCoord2d( (GLdouble) *s, (GLdouble) *t );
}

void APIENTRY fglTexCoord2dv( f90gldouble ARRAY v )
{
   glTexCoord2dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglTexCoord2f( f90glfloat s, f90glfloat t )
{
   glTexCoord2f( (GLfloat) *s, (GLfloat) *t );
}

void APIENTRY fglTexCoord2fv( f90glfloat ARRAY v )
{
   glTexCoord2fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglTexCoord2i( f90glint s, f90glint t )
{
   glTexCoord2i( (GLint) *s, (GLint) *t );
}

void APIENTRY fglTexCoord2iv( f90glint ARRAY v )
{
   glTexCoord2iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglTexCoord2s( f90glshort s, f90glshort t )
{
   glTexCoord2s( (GLshort) *s, (GLshort) *t );
}

void APIENTRY f9aglTexCoord2s(f90glint s,f90glint t )
{
   glTexCoord2s( (GLshort) *s, (GLshort) *t);
}

void APIENTRY fglTexCoord2sv( f90glshort ARRAY v )
{
   glTexCoord2sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglTexCoord2sv(f90glint ARRAY v )
{
   GLshort short_v[2];
   int i;
   for (i=0;i<2;i++) short_v[i]=(GLshort) v ARRAYREF[i];
   glTexCoord2sv( short_v);
}

void APIENTRY fglTexCoord3d( f90gldouble s, f90gldouble t, f90gldouble r )
{
   glTexCoord3d( (GLdouble) *s, (GLdouble) *t, (GLdouble) *r );
}

void APIENTRY fglTexCoord3dv( f90gldouble ARRAY v )
{
   glTexCoord3dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglTexCoord3f( f90glfloat s, f90glfloat t, f90glfloat r )
{
   glTexCoord3f( (GLfloat) *s, (GLfloat) *t, (GLfloat) *r );
}

void APIENTRY fglTexCoord3fv( f90glfloat ARRAY v )
{
   glTexCoord3fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglTexCoord3i( f90glint s, f90glint t, f90glint r )
{
   glTexCoord3i( (GLint) *s, (GLint) *t, (GLint) *r );
}

void APIENTRY fglTexCoord3iv( f90glint ARRAY v )
{
   glTexCoord3iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglTexCoord3s( f90glshort s, f90glshort t, f90glshort r )
{
   glTexCoord3s( (GLshort) *s, (GLshort) *t, (GLshort) *r );
}

void APIENTRY f9aglTexCoord3s(f90glint s,f90glint t,f90glint r )
{
   glTexCoord3s( (GLshort) *s, (GLshort) *t, (GLshort) *r);
}

void APIENTRY fglTexCoord3sv( f90glshort ARRAY v )
{
   glTexCoord3sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglTexCoord3sv(f90glint ARRAY v )
{
   GLshort short_v[3];
   int i;
   for (i=0;i<3;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glTexCoord3sv( short_v);
}

void APIENTRY fglTexCoord4d( f90gldouble s, f90gldouble t, f90gldouble r,
    f90gldouble q )
{
   glTexCoord4d( (GLdouble) *s, (GLdouble) *t, (GLdouble) *r,
    (GLdouble) *q );
}

void APIENTRY fglTexCoord4dv( f90gldouble ARRAY v )
{
   glTexCoord4dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglTexCoord4f( f90glfloat s, f90glfloat t, f90glfloat r,
    f90glfloat q )
{
   glTexCoord4f( (GLfloat) *s, (GLfloat) *t, (GLfloat) *r, (GLfloat) *q );
}

void APIENTRY fglTexCoord4fv( f90glfloat ARRAY v )
{
   glTexCoord4fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglTexCoord4i( f90glint s, f90glint t, f90glint r, f90glint q )
{
   glTexCoord4i( (GLint) *s, (GLint) *t, (GLint) *r, (GLint) *q );
}

void APIENTRY fglTexCoord4iv( f90glint ARRAY v )
{
   glTexCoord4iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglTexCoord4s( f90glshort s, f90glshort t, f90glshort r,
    f90glshort q )
{
   glTexCoord4s( (GLshort) *s, (GLshort) *t, (GLshort) *r, (GLshort) *q );
}

void APIENTRY f9aglTexCoord4s(f90glint s,f90glint t,f90glint r,f90glint q )
{
   glTexCoord4s( (GLshort) *s, (GLshort) *t, (GLshort) *r, (GLshort) *q);
}

void APIENTRY fglTexCoord4sv( f90glshort ARRAY v )
{
   glTexCoord4sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglTexCoord4sv(f90glint ARRAY v )
{
   GLshort short_v[4];
   int i;
   for (i=0;i<4;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glTexCoord4sv( short_v);
}

void APIENTRY fglTexEnvf( f90glenum target, f90glenum pname,
    f90glfloat param )
{
   glTexEnvf( (GLenum) *target, (GLenum) *pname, (GLfloat) *param );
}

void APIENTRY fglTexEnvfv( f90glenum target, f90glenum pname,
    f90glfloat ARRAY params )
{
   glTexEnvfv( (GLenum) *target, (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglTexEnvi( f90glenum target, f90glenum pname, f90glint param )
{
   glTexEnvi( (GLenum) *target, (GLenum) *pname, (GLint) *param );
}

void APIENTRY fglTexEnviv( f90glenum target, f90glenum pname,
    f90glint ARRAY params )
{
   glTexEnviv( (GLenum) *target, (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY fglTexGend( f90glenum coord, f90glenum pname,
    f90gldouble param )
{
   glTexGend( (GLenum) *coord, (GLenum) *pname, (GLdouble) *param );
}

void APIENTRY fglTexGendv( f90glenum coord, f90glenum pname,
    f90gldouble ARRAY params )
{
   glTexGendv( (GLenum) *coord, (GLenum) *pname, (GLdouble *) params ARRAYREF);
}

void APIENTRY fglTexGenf( f90glenum coord, f90glenum pname,
    f90glfloat param )
{
   glTexGenf( (GLenum) *coord, (GLenum) *pname, (GLfloat) *param );
}

void APIENTRY fglTexGenfv( f90glenum coord, f90glenum pname,
    f90glfloat ARRAY params )
{
   glTexGenfv( (GLenum) *coord, (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglTexGeni( f90glenum coord, f90glenum pname, f90glint param )
{
   glTexGeni( (GLenum) *coord, (GLenum) *pname, (GLint) *param );
}

void APIENTRY fglTexGeniv( f90glenum coord, f90glenum pname,
    f90glint ARRAY params )
{
   glTexGeniv( (GLenum) *coord, (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY f9y0glTexImage1D( f90glenum target, f90glint level,
   f90glint components, f90glsizei width,
   f90glint border, f90glenum format, f90glenum type, f90glvoid ARRAY pixels,
   f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newpixels;
   int i;
   newpixels = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLubyte)pixels ARRAYREF [i];
   glTexImage1D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newpixels;
   int i;
   newpixels = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLbyte)pixels ARRAYREF [i];
   glTexImage1D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newpixels;
   int i;
   newpixels = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newpixels[i]=(GLushort)pixels ARRAYREF [i];
   glTexImage1D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newpixels;
   int i;
   newpixels = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newpixels[i]=(GLshort)pixels ARRAYREF [i];
   glTexImage1D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 default:
   glTexImage1D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}
}

void APIENTRY f9y1glTexImage1D( f90glenum target, f90glint level,
   f90glint components, f90glsizei width,
   f90glint border, f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexImage1D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glTexImage1D( f90glenum target, f90glint level,
   f90glint components, f90glsizei width,
   f90glint border, f90glenum format, f90glenum type, f90gluint ARRAY temp,
   f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexImage1D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glTexImage1D( f90glenum target, f90glint level,
   f90glint components, f90glsizei width,
   f90glint border, f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexImage1D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glTexImage1D( f90glenum target, f90glint level,
   f90glint components, f90glsizei width,
   f90glint border, f90glenum format, f90glenum type, f90glshort ARRAY pixels,
   f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newpixels;
   int i;
   newpixels = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLubyte)pixels ARRAYREF [i];
   glTexImage1D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newpixels;
   int i;
   newpixels = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLbyte)pixels ARRAYREF [i];
   glTexImage1D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 default:
   glTexImage1D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}
}

void APIENTRY f9y0glTexImage2D( f90glenum target, f90glint level,
   f90glint components, f90glsizei width, f90glsizei height,
   f90glint border, f90glenum format, f90glenum type, f90glvoid ARRAY pixels,
   f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newpixels;
   int i;
   newpixels = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLubyte)pixels ARRAYREF [i];
   glTexImage2D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newpixels;
   int i;
   newpixels = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLbyte)pixels ARRAYREF [i];
   glTexImage2D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newpixels;
   int i;
   newpixels = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newpixels[i]=(GLushort)pixels ARRAYREF [i];
   glTexImage2D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newpixels;
   int i;
   newpixels = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newpixels[i]=(GLshort)pixels ARRAYREF [i];
   glTexImage2D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 default:
   glTexImage2D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}
}

void APIENTRY f9y1glTexImage2D( f90glenum target, f90glint level,
   f90glint components, f90glsizei width, f90glsizei height,
   f90glint border, f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexImage2D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glTexImage2D( f90glenum target, f90glint level,
   f90glint components, f90glsizei width, f90glsizei height,
   f90glint border, f90glenum format, f90glenum type, f90gluint ARRAY temp,
   f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexImage2D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glTexImage2D( f90glenum target, f90glint level,
   f90glint components, f90glsizei width, f90glsizei height,
   f90glint border, f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexImage2D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glTexImage2D( f90glenum target, f90glint level,
   f90glint components, f90glsizei width, f90glsizei height,
   f90glint border, f90glenum format, f90glenum type, f90glshort ARRAY pixels,
   f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newpixels;
   int i;
   newpixels = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLubyte)pixels ARRAYREF [i];
   glTexImage2D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newpixels;
   int i;
   newpixels = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newpixels[i]=(GLbyte)pixels ARRAYREF [i];
   glTexImage2D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) newpixels );
   free(newpixels);
   break;
  }
 default:
   glTexImage2D( (GLenum) *target, (GLint) *level, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, (GLint) *border,
      (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}
}

void APIENTRY fglTexParameterf( f90glenum target, f90glenum pname,
    f90glfloat param )
{
   glTexParameterf( (GLenum) *target, (GLenum) *pname, (GLfloat) *param );
}

void APIENTRY fglTexParameterfv( f90glenum target, f90glenum pname,
    f90glfloat ARRAY params )
{
   glTexParameterfv( (GLenum) *target, (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

void APIENTRY fglTexParameteri( f90glenum target, f90glenum pname,
    f90glint param )
{
   glTexParameteri( (GLenum) *target, (GLenum) *pname, (GLint) *param );
}

void APIENTRY fglTexParameteriv( f90glenum target, f90glenum pname,
    f90glint ARRAY params )
{
   glTexParameteriv( (GLenum) *target, (GLenum) *pname, (GLint *) params ARRAYREF);
}

void APIENTRY fglTranslated( f90gldouble x, f90gldouble y, f90gldouble z )
{
   glTranslated( (GLdouble) *x, (GLdouble) *y, (GLdouble) *z );
}

void APIENTRY fglTranslatef( f90glfloat x, f90glfloat y, f90glfloat z )
{
   glTranslatef( (GLfloat) *x, (GLfloat) *y, (GLfloat) *z );
}

void APIENTRY fglVertex2d( f90gldouble x, f90gldouble y )
{
   glVertex2d( (GLdouble) *x, (GLdouble) *y );
}

void APIENTRY fglVertex2dv( f90gldouble ARRAY v )
{
   glVertex2dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglVertex2f( f90glfloat x, f90glfloat y )
{
   glVertex2f( (GLfloat) *x, (GLfloat) *y );
}

void APIENTRY fglVertex2fv( f90glfloat ARRAY v )
{
   glVertex2fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglVertex2i( f90glint x, f90glint y )
{
   glVertex2i( (GLint) *x, (GLint) *y );
}

void APIENTRY fglVertex2iv( f90glint ARRAY v )
{
   glVertex2iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglVertex2s( f90glshort x, f90glshort y )
{
   glVertex2s( (GLshort) *x, (GLshort) *y );
}

void APIENTRY f9aglVertex2s(f90glint x,f90glint y )
{
   glVertex2s( (GLshort) *x, (GLshort) *y);
}

void APIENTRY fglVertex2sv( f90glshort ARRAY v )
{
   glVertex2sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglVertex2sv(f90glint ARRAY v )
{
   GLshort short_v[2];
   int i;
   for (i=0;i<2;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glVertex2sv( short_v);
}

void APIENTRY fglVertex3d( f90gldouble x, f90gldouble y, f90gldouble z )
{
   glVertex3d( (GLdouble) *x, (GLdouble) *y, (GLdouble) *z );
}

void APIENTRY fglVertex3dv( f90gldouble ARRAY v )
{
   glVertex3dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglVertex3f( f90glfloat x, f90glfloat y, f90glfloat z )
{
   glVertex3f( (GLfloat) *x, (GLfloat) *y, (GLfloat) *z );
}

void APIENTRY fglVertex3fv( f90glfloat ARRAY v )
{
   glVertex3fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglVertex3i( f90glint x, f90glint y, f90glint z )
{
   glVertex3i( (GLint) *x, (GLint) *y, (GLint) *z );
}

void APIENTRY fglVertex3iv( f90glint ARRAY v )
{
   glVertex3iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglVertex3s( f90glshort x, f90glshort y, f90glshort z )
{
   glVertex3s( (GLshort) *x, (GLshort) *y, (GLshort) *z );
}

void APIENTRY f9aglVertex3s(f90glint x,f90glint y,f90glint z )
{
   glVertex3s( (GLshort) *x, (GLshort) *y, (GLshort) *z);
}

void APIENTRY fglVertex3sv( f90glshort ARRAY v )
{
   glVertex3sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglVertex3sv(f90glint ARRAY v )
{
   GLshort short_v[3];
   int i;
   for (i=0;i<3;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glVertex3sv( short_v);
}

void APIENTRY fglVertex4d( f90gldouble x, f90gldouble y, f90gldouble z,
    f90gldouble w )
{
   glVertex4d( (GLdouble) *x, (GLdouble) *y, (GLdouble) *z, (GLdouble) *w );
}

void APIENTRY fglVertex4dv( f90gldouble ARRAY v )
{
   glVertex4dv( (GLdouble *) v ARRAYREF);
}

void APIENTRY fglVertex4f( f90glfloat x, f90glfloat y, f90glfloat z,
    f90glfloat w )
{
   glVertex4f( (GLfloat) *x, (GLfloat) *y, (GLfloat) *z, (GLfloat) *w );
}

void APIENTRY fglVertex4fv( f90glfloat ARRAY v )
{
   glVertex4fv( (GLfloat *) v ARRAYREF);
}

void APIENTRY fglVertex4i( f90glint x, f90glint y, f90glint z, f90glint w )
{
   glVertex4i( (GLint) *x, (GLint) *y, (GLint) *z, (GLint) *w );
}

void APIENTRY fglVertex4iv( f90glint ARRAY v )
{
   glVertex4iv( (GLint *) v ARRAYREF);
}

void APIENTRY fglVertex4s( f90glshort x, f90glshort y, f90glshort z,
    f90glshort w )
{
   glVertex4s( (GLshort) *x, (GLshort) *y, (GLshort) *z, (GLshort) *w );
}

void APIENTRY f9aglVertex4s(f90glint x,f90glint y,f90glint z,f90glint w )
{
   glVertex4s( (GLshort) *x, (GLshort) *y, (GLshort) *z, (GLshort) *w);
}

void APIENTRY fglVertex4sv( f90glshort ARRAY v )
{
   glVertex4sv( (GLshort *) v ARRAYREF);
}

void APIENTRY f9aglVertex4sv(f90glint ARRAY v )
{
   GLshort short_v[4];
   int i;
   for (i=0;i<4;i++) short_v[i]=(GLshort) v ARRAYREF [i];
   glVertex4sv( short_v);
}

void APIENTRY fglViewport( f90glint x, f90glint y, f90glsizei width,
    f90glsizei height )
{
   glViewport( (GLint) *x, (GLint) *y, (GLsizei) *width, (GLsizei) *height );
}


/* OpenGL 1.1 */

#ifdef OPENGL_1_1

#ifdef PURE_FUNCTIONS
void APIENTRY fglAreTexturesResident( f90glsizei n,
    f90gluint ARRAY textures, f90glboolean ARRAY residences, f90glboolean areresident )
#else
GLboolean APIENTRY fglAreTexturesResident( f90glsizei n,
    f90gluint ARRAY textures, f90glboolean ARRAY residences )
#endif
{
GLboolean retval;
   retval = (GLboolean) glAreTexturesResident( (GLsizei) *n,
    (GLuint *) textures ARRAYREF, (GLboolean *) residences  ARRAYREF);
#ifdef PURE_FUNCTIONS
   *areresident = retval;
#else
   return retval;
#endif
}
 
#ifdef PURE_FUNCTIONS
void APIENTRY f9eglAreTexturesResident(f90glsizei n,
   f90gluint ARRAY textures,f90glint ARRAY residences,f90glint areresident )
#else
GLint APIENTRY f9eglAreTexturesResident(f90glsizei n,
   f90gluint ARRAY textures,f90glint ARRAY residences )
#endif
{
   GLboolean *boolean_residences;
   GLint i, iresult;
   boolean_residences = (GLboolean *)malloc((*n)*sizeof(GLboolean));
   iresult = (GLint) glAreTexturesResident( (GLsizei) *n,
      (GLuint *) textures ARRAYREF, boolean_residences);
   for (i=0;i<*n;i++) residences ARRAYREF [i]=(GLint) boolean_residences[i];
   free(boolean_residences);
#ifdef PURE_FUNCTIONS
   *areresident = iresult;
#else
   return iresult;
#endif
}

void APIENTRY fglArrayElement( f90glint i )
{
   glArrayElement( (GLint) *i );
}

void APIENTRY fglBindTexture( f90glenum target, f90gluint texture )
{
   glBindTexture( (GLenum) *target, (GLuint) *texture );
}

void APIENTRY f9y0glColorPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{

   glColorPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y1glColorPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glColorPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y4glColorPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glColorPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr );
}

void APIENTRY f9y5glColorPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glColorPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y6glColorPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glColorPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y7glColorPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glColorPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY fglCopyTexImage1D( f90glenum target, f90glint level,
    f90glenum internalformat, f90glint x, f90glint y, f90glsizei width,
    f90glint border )
{
   glCopyTexImage1D( (GLenum) *target, (GLint) *level,
    (GLenum) *internalformat, (GLint) *x, (GLint) *y, (GLsizei) *width,
    (GLint) *border );
}

void APIENTRY fglCopyTexImage2D( f90glenum target, f90glint level,
    f90glenum internalformat, f90glint x, f90glint y, f90glsizei width,
    f90glsizei height, f90glint border )
{
   glCopyTexImage2D( (GLenum) *target, (GLint) *level,
    (GLenum) *internalformat, (GLint) *x, (GLint) *y, (GLsizei) *width,
    (GLsizei) *height, (GLint) *border );
}

void APIENTRY fglCopyTexSubImage1D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint x, f90glint y, f90glsizei width )
{
   glCopyTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *x, (GLint) *y, (GLsizei) *width );
}

void APIENTRY fglCopyTexSubImage2D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint x, f90glint y,
    f90glsizei width, f90glsizei height )
{
   glCopyTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *x, (GLint) *y, (GLsizei) *width,
    (GLsizei) *height );
}

void APIENTRY fglDeleteTextures( f90glsizei n, f90gluint ARRAY textures )
{
   glDeleteTextures( (GLsizei) *n, (GLuint *) textures ARRAYREF);
}

void APIENTRY fglDisableClientState( f90glenum cap )
{
   glDisableClientState( (GLenum) *cap );
}

void APIENTRY fglDrawArrays( f90glenum mode, f90glint first,
    f90glsizei count )
{
   glDrawArrays( (GLenum) *mode, (GLint) *first, (GLsizei) *count );
}

void APIENTRY f9y0glDrawElements( f90glenum mode,
    f90glsizei count, f90glenum type, f90glvoid ARRAY indices )
{

switch (*mode) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*count)*sizeof(GLubyte));
   for (i=0;i<(*count);i++) newlists[i]=(GLubyte)indices ARRAYREF [i];
   glDrawElements( (GLenum) *mode, 
    (GLsizei) *count, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*count)*sizeof(GLushort));
   for (i=0;i<(*count);i++) newlists[i]=(GLushort)indices ARRAYREF [i];
   glDrawElements( (GLenum) *mode, 
    (GLsizei) *count, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glDrawElements( (GLenum) *mode, 
    (GLsizei) *count, (GLenum) *type, (GLvoid *) indices ARRAYREF);
   break;
}     
}

void APIENTRY f9y4glDrawElements( f90glenum mode,
    f90glsizei count, f90glenum type, f90glvoid ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glDrawElements( (GLenum) *mode, (GLsizei) *count, (GLenum) *type,
    (GLvoid *) ptr );
}

void APIENTRY f9y5glDrawElements( f90glenum mode,
    f90glsizei count, f90glenum type, f90glvoid ARRAY indices )
{
   glDrawElements( (GLenum) *mode, 
    (GLsizei) *count, (GLenum) *type, (GLvoid *) indices ARRAYREF);
}

void APIENTRY f9y6glDrawElements( f90glenum mode,
    f90glsizei count, f90glenum type, f90glshort ARRAY indices )
{

switch (*mode) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*count)*sizeof(GLubyte));
   for (i=0;i<(*count);i++) newlists[i]=(GLubyte)indices ARRAYREF [i];
   glDrawElements( (GLenum) *mode, 
    (GLsizei) *count, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glDrawElements( (GLenum) *mode, 
    (GLsizei) *count, (GLenum) *type, (GLvoid *) indices  ARRAYREF );
   break;
}     
}

void APIENTRY fglEdgeFlagPointer( f90glsizei stride, f90glboolean ARRAY ptr )
{
   glEdgeFlagPointer( (GLsizei) *stride, (GLboolean *) ptr ARRAYREF);
}

void APIENTRY f9eglEdgeFlagPointer(f90glsizei stride, f90glbyte ARRAY ptr )
{
   glEdgeFlagPointer( (GLsizei) *stride, (GLboolean *) ptr ARRAYREF);
}

void APIENTRY fglEnableClientState( f90glenum cap )
{
   glEnableClientState( (GLenum) *cap );
}

void APIENTRY fglGenTextures( f90glsizei n, f90gluint ARRAY textures )
{
   glGenTextures( (GLsizei) *n, (GLuint *) textures ARRAYREF);
}

void APIENTRY fglGetPointerv( f90glenum pname, f90gluint ARRAY params,
   f90glint length )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   glGetPointerv( (GLenum) *pname, (GLvoid **)&ptr );
   p = (unsigned char *) ptr;
   for (i=0; i<(*length); i++) {params ARRAYREF [i]=(GLuint)*p; p++;}
}

void APIENTRY f9y0glIndexPointer( f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glIndexPointer( (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y1glIndexPointer( f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glIndexPointer( (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y4glIndexPointer( f90glenum type,
    f90glsizei stride, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glIndexPointer( (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr );
}

void APIENTRY f9y6glIndexPointer( f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{

   glIndexPointer( (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y7glIndexPointer( f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glIndexPointer( (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY fglIndexub( f90glubyte c )
{
   glIndexub( (GLubyte) *c );
}

void APIENTRY f9bglIndexub(f90gluint c )
{
   glIndexub( (GLubyte) *c);
}

void APIENTRY f9cglIndexub(f90glushort c )
{
   glIndexub( (GLubyte) *c);
}

void APIENTRY fglIndexubv( f90glubyte ARRAY c )
{
   glIndexubv( (GLubyte *) c ARRAYREF);
}

void APIENTRY f9bglIndexubv(f90gluint ARRAY c )
{
   GLubyte byte_c[1];
   int i;
   for (i=0;i<1;i++) byte_c[i]=(GLubyte) c ARRAYREF [i];
   glIndexubv( byte_c);
}

void APIENTRY f9cglIndexubv(f90glushort ARRAY c )
{
   GLubyte byte_c[1];
   int i;
   for (i=0;i<1;i++) byte_c[i]=(GLubyte) c ARRAYREF [i];
   glIndexubv( byte_c);
}

void APIENTRY f9y0glInterleavedArrays( f90glenum format,
    f90glsizei stride, f90glvoid ARRAY ptr, f90glint length )
{
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)ptr ARRAYREF [i];
   glInterleavedArrays( (GLenum) *format, (GLsizei) *stride,
    (GLvoid *) newlists );
   free(newlists);
}

void APIENTRY f9y1glInterleavedArrays( f90glenum format,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glInterleavedArrays( (GLenum) *format, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y4glInterleavedArrays( f90glenum format,
    f90glsizei stride, f90glvoid ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glInterleavedArrays( (GLenum) *format, (GLsizei) *stride,
    (GLvoid *) ptr );
}

void APIENTRY f9y5glInterleavedArrays( f90glenum format,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glInterleavedArrays( (GLenum) *format, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y6glInterleavedArrays( f90glenum format,
    f90glsizei stride, f90glshort ARRAY ptr, f90glint length )
{
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)ptr ARRAYREF [i];
   glInterleavedArrays( (GLenum) *format, (GLsizei) *stride,
    (GLvoid *) newlists );
   free(newlists);
}

GLboolean APIENTRY fglIsTexture( f90gluint texture )
{
   return (GLboolean) glIsTexture( (GLuint) *texture );
}

GLint APIENTRY f9eglIsTexture( f90gluint texture )
{
   GLboolean temp;
   temp = glIsTexture( (GLuint) *texture );
   return (GLint) temp;
}

void APIENTRY f9y0glNormalPointer( f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glNormalPointer( (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y1glNormalPointer( f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glNormalPointer( (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y4glNormalPointer( f90glenum type,
    f90glsizei stride, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glNormalPointer( (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr );
}

void APIENTRY f9y5glNormalPointer( f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glNormalPointer( (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y6glNormalPointer( f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{

   glNormalPointer( (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y7glNormalPointer( f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glNormalPointer( (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY fglPolygonOffset( f90glfloat factor, f90glfloat units )
{
   glPolygonOffset( (GLfloat) *factor, (GLfloat) *units );
}

void APIENTRY fglPopClientAttrib( void )
{
   glPopClientAttrib();
}

void APIENTRY fglPrioritizeTextures( f90glsizei n, f90gluint ARRAY textures,
    f90glclampf ARRAY priorities )
{
   glPrioritizeTextures( (GLsizei) *n, (GLuint *) textures ARRAYREF,
    (GLclampf *) priorities  ARRAYREF);
}

void APIENTRY fglPushClientAttrib( f90glbitfield mask )
{
   glPushClientAttrib( (GLbitfield) *mask );
}

void APIENTRY f9y0glTexCoordPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glTexCoordPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y1glTexCoordPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glTexCoordPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y4glTexCoordPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexCoordPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr );
}

void APIENTRY f9y6glTexCoordPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{

   glTexCoordPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y7glTexCoordPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glTexCoordPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y0glTexSubImage1D( f90glenum target, f90glint level,
    f90glint xoffset, f90glsizei width, f90glenum format, f90glenum type,
    f90glvoid ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newlists[i]=(GLushort)pixels ARRAYREF [i];
   glTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newlists[i]=(GLshort)pixels ARRAYREF [i];
   glTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y1glTexSubImage1D( f90glenum target, f90glint level,
    f90glint xoffset, f90glsizei width, f90glenum format, f90glenum type,
    f90glvoid ARRAY pixels )
{
   glTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glTexSubImage1D( f90glenum target, f90glint level,
    f90glint xoffset, f90glsizei width, f90glenum format, f90glenum type,
    f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glTexSubImage1D( f90glenum target, f90glint level,
    f90glint xoffset, f90glsizei width, f90glenum format, f90glenum type,
    f90glvoid ARRAY pixels )
{
   glTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glTexSubImage1D( f90glenum target, f90glint level,
    f90glint xoffset, f90glsizei width, f90glenum format, f90glenum type,
    f90glshort ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage1D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y0glTexSubImage2D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glsizei width, f90glsizei height,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newlists[i]=(GLushort)pixels ARRAYREF [i];
   glTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newlists[i]=(GLshort)pixels ARRAYREF [i];
   glTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y1glTexSubImage2D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glsizei width, f90glsizei height,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glTexSubImage2D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glsizei width, f90glsizei height,
    f90glenum format, f90glenum type, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glTexSubImage2D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glsizei width, f90glsizei height,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glTexSubImage2D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glsizei width, f90glsizei height,
    f90glenum format, f90glenum type, f90glshort ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage2D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y0glVertexPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glVertexPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y1glVertexPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glVertexPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y4glVertexPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glVertexPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr );
}

void APIENTRY f9y6glVertexPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glVertexPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y7glVertexPointer( f90glint size, f90glenum type,
    f90glsizei stride, f90glvoid ARRAY ptr )
{
   glVertexPointer( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLvoid *) ptr ARRAYREF);
}

#endif

/* OpenGL 1.2 */

#ifdef OPENGL_1_2

void APIENTRY f9y0glDrawRangeElements( f90glenum mode, f90gluint start,
    f90gluint end, f90glsizei count, f90glenum type, f90glvoid ARRAY indices )
{

switch (*mode) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*count)*sizeof(GLubyte));
   for (i=0;i<(*count);i++) newlists[i]=(GLubyte)indices ARRAYREF [i];
   glDrawRangeElements( (GLenum) *mode, (GLuint) *start, (GLuint) *end,
    (GLsizei) *count, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*count)*sizeof(GLushort));
   for (i=0;i<(*count);i++) newlists[i]=(GLushort)indices ARRAYREF [i];
   glDrawRangeElements( (GLenum) *mode, (GLuint) *start, (GLuint) *end,
    (GLsizei) *count, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glDrawRangeElements( (GLenum) *mode, (GLuint) *start, (GLuint) *end,
    (GLsizei) *count, (GLenum) *type, (GLvoid *) indices  ARRAYREF );
   break;
}     
}

void APIENTRY f9y4glDrawRangeElements( f90glenum mode, f90gluint start,
    f90gluint end, f90glsizei count, f90glenum type, f90glvoid ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glDrawRangeElements( (GLenum) *mode, (GLuint) *start, (GLuint) *end,
    (GLsizei) *count, (GLenum) *type, (GLvoid *) ptr );
}

void APIENTRY f9y5glDrawRangeElements( f90glenum mode, f90gluint start,
    f90gluint end, f90glsizei count, f90glenum type, f90glvoid ARRAY indices )
{
   glDrawRangeElements( (GLenum) *mode, (GLuint) *start, (GLuint) *end,
    (GLsizei) *count, (GLenum) *type, (GLvoid *) indices ARRAYREF);
}

void APIENTRY f9y6glDrawRangeElements( f90glenum mode, f90gluint start,
    f90gluint end, f90glsizei count, f90glenum type, f90glshort ARRAY indices )
{

switch (*mode) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*count)*sizeof(GLubyte));
   for (i=0;i<(*count);i++) newlists[i]=(GLubyte)indices ARRAYREF [i];
   glDrawRangeElements( (GLenum) *mode, (GLuint) *start, (GLuint) *end,
    (GLsizei) *count, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glDrawRangeElements( (GLenum) *mode, (GLuint) *start, (GLuint) *end,
    (GLsizei) *count, (GLenum) *type, (GLvoid *) indices  ARRAYREF );
   break;
}     
}

void APIENTRY f9y0glTexImage3D( f90glenum target, f90glint level,
    f90glenum internalFormat, f90glsizei width,
    f90glsizei height, f90glsizei depth, f90glint border,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexImage3D( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexImage3D( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newlists[i]=(GLushort)pixels ARRAYREF [i];
   glTexImage3D( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newlists[i]=(GLshort)pixels ARRAYREF [i];
   glTexImage3D( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexImage3D( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y1glTexImage3D( f90glenum target, f90glint level,
    f90glenum internalFormat, f90glsizei width,
    f90glsizei height, f90glsizei depth, f90glint border,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexImage3D( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glTexImage3D( f90glenum target, f90glint level,
    f90glenum internalFormat, f90glsizei width,
    f90glsizei height, f90glsizei depth, f90glint border,
    f90glenum format, f90glenum type, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexImage3D( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glTexImage3D( f90glenum target, f90glint level,
    f90glenum internalFormat, f90glsizei width,
    f90glsizei height, f90glsizei depth, f90glint border,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexImage3D( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glTexImage3D( f90glenum target, f90glint level,
    f90glenum internalFormat, f90glsizei width,
    f90glsizei height, f90glsizei depth, f90glint border,
    f90glenum format, f90glenum type, f90glshort ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexImage3D( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexImage3D( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexImage3D( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y0glTexSubImage3D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glsizei width,
    f90glsizei height, f90glsizei depth,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newlists[i]=(GLushort)pixels ARRAYREF [i];
   glTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newlists[i]=(GLshort)pixels ARRAYREF [i];
   glTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y1glTexSubImage3D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glsizei width,
    f90glsizei height, f90glsizei depth,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glTexSubImage3D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glsizei width,
    f90glsizei height, f90glsizei depth,
    f90glenum format, f90glenum type, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glTexSubImage3D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glsizei width,
    f90glsizei height, f90glsizei depth,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glTexSubImage3D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glsizei width,
    f90glsizei height, f90glsizei depth,
    f90glenum format, f90glenum type, f90glshort ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY fglCopyTexSubImage3D( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glint x,
    f90glint y, f90glsizei width, f90glsizei height )
{
   glCopyTexSubImage3D( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLint) *x, (GLint) *y,
    (GLsizei) *width, (GLsizei) *height );
}

#endif

/* EXT_blend_color */

#ifdef GL_EXT_blend_color

void APIENTRY fglBlendColorEXT( f90glclampf red, f90glclampf green,
    f90glclampf blue, f90glclampf alpha )
{
   glBlendColorEXT( (GLclampf) *red, (GLclampf) *green, (GLclampf) *blue,
    (GLclampf) *alpha );
}

#endif

/* EXT_blend_minmax */

#ifdef GL_EXT_blend_minmax

void APIENTRY fglBlendEquationEXT( f90glenum mode )
{
   glBlendEquationEXT( (GLenum) *mode );
}

#endif

/* EXT_color_table */

#ifdef GL_EXT_color_table

void APIENTRY f9y0glColorTableEXT( f90glenum target, f90glenum internalformat,
    f90glsizei width, f90glenum format, f90glenum type,
    f90glvoid ARRAY table, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)table ARRAYREF [i];
   glColorTableEXT( (GLenum) *target, (GLenum) *internalformat,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)table ARRAYREF [i];
   glColorTableEXT( (GLenum) *target, (GLenum) *internalformat,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newlists[i]=(GLushort)table ARRAYREF [i];
   glColorTableEXT( (GLenum) *target, (GLenum) *internalformat,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newlists[i]=(GLshort)table ARRAYREF [i];
   glColorTableEXT( (GLenum) *target, (GLenum) *internalformat,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glColorTableEXT( (GLenum) *target, (GLenum) *internalformat,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) table  ARRAYREF );
   break;
}     
}

void APIENTRY f9y1glColorTableEXT( f90glenum target, f90glenum internalformat,
    f90glsizei width, f90glenum format, f90glenum type,
    f90glvoid ARRAY table )
{
   glColorTableEXT( (GLenum) *target, (GLenum) *internalformat,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) table ARRAYREF);
}

void APIENTRY f9y4glColorTableEXT( f90glenum target, f90glenum internalformat,
    f90glsizei width, f90glenum format, f90glenum type,
    f90glvoid ARRAY temp, f90glint nbytes )
{
   GLvoid *table;
   unsigned char *p;
   int i;
   p = (unsigned char *) &table;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glColorTableEXT( (GLenum) *target, (GLenum) *internalformat,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) table );
}

void APIENTRY f9y5glColorTableEXT( f90glenum target, f90glenum internalformat,
    f90glsizei width, f90glenum format, f90glenum type,
    f90glvoid ARRAY table )
{
   glColorTableEXT( (GLenum) *target, (GLenum) *internalformat,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) table ARRAYREF);
}

void APIENTRY f9y6glColorTableEXT( f90glenum target, f90glenum internalformat,
    f90glsizei width, f90glenum format, f90glenum type,
    f90glshort ARRAY table, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)table ARRAYREF [i];
   glColorTableEXT( (GLenum) *target, (GLenum) *internalformat,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   glColorTableEXT( (GLenum) *target, (GLenum) *internalformat,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glColorTableEXT( (GLenum) *target, (GLenum) *internalformat,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) table  ARRAYREF );
   break;
}     
}

void APIENTRY f9y0glColorSubTableEXT( f90glenum target, f90glsizei start,
    f90glsizei count, f90glenum format, f90glenum type,
    f90glvoid ARRAY data, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)data ARRAYREF [i];
   glColorSubTableEXT( (GLenum) *target, (GLsizei) *start,
    (GLsizei) *count, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)data ARRAYREF [i];
   glColorSubTableEXT( (GLenum) *target, (GLsizei) *start,
    (GLsizei) *count, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newlists[i]=(GLushort)data ARRAYREF [i];
   glColorSubTableEXT( (GLenum) *target, (GLsizei) *start,
    (GLsizei) *count, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newlists[i]=(GLshort)data ARRAYREF [i];
   glColorSubTableEXT( (GLenum) *target, (GLsizei) *start,
    (GLsizei) *count, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glColorSubTableEXT( (GLenum) *target, (GLsizei) *start,
    (GLsizei) *count, (GLenum) *format, (GLenum) *type, (GLvoid *) data  ARRAYREF );
   break;
}     
}

void APIENTRY f9y1glColorSubTableEXT( f90glenum target, f90glsizei start,
    f90glsizei count, f90glenum format, f90glenum type,
    f90glvoid ARRAY data )
{
   glColorSubTableEXT( (GLenum) *target, (GLsizei) *start,
    (GLsizei) *count, (GLenum) *format, (GLenum) *type, (GLvoid *) data ARRAYREF);
}

void APIENTRY f9y4glColorSubTableEXT( f90glenum target, f90glsizei start,
    f90glsizei count, f90glenum format, f90glenum type,
    f90glvoid ARRAY temp, f90glint nbytes )
{
   GLvoid *data;
   unsigned char *p;
   int i;
   p = (unsigned char *) &data;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glColorSubTableEXT( (GLenum) *target, (GLsizei) *start,
    (GLsizei) *count, (GLenum) *format, (GLenum) *type, (GLvoid *) data );
}

void APIENTRY f9y5glColorSubTableEXT( f90glenum target, f90glsizei start,
    f90glsizei count, f90glenum format, f90glenum type,
    f90glvoid ARRAY data )
{
   glColorSubTableEXT( (GLenum) *target, (GLsizei) *start,
    (GLsizei) *count, (GLenum) *format, (GLenum) *type, (GLvoid *) data ARRAYREF);
}

void APIENTRY f9y6glColorSubTableEXT( f90glenum target, f90glsizei start,
    f90glsizei count, f90glenum format, f90glenum type,
    f90glshort ARRAY data, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)data ARRAYREF [i];
   glColorSubTableEXT( (GLenum) *target, (GLsizei) *start,
    (GLsizei) *count, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)data ARRAYREF [i];
   glColorSubTableEXT( (GLenum) *target, (GLsizei) *start,
    (GLsizei) *count, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glColorSubTableEXT( (GLenum) *target, (GLsizei) *start,
    (GLsizei) *count, (GLenum) *format, (GLenum) *type, (GLvoid *) data  ARRAYREF );
   break;
}     
}

void APIENTRY f9y0glGetColorTableEXT( f90glenum target, f90glenum format,
    f90glenum type, f90glvoid ARRAY table, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   glGetColorTableEXT( (GLenum) *target, (GLenum) *format, (GLenum) *type,
      (GLvoid *) newlists );
   for (i=0;i<(*length);i++) table ARRAYREF [i]=(GLint)newlists[i];
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   glGetColorTableEXT( (GLenum) *target, (GLenum) *format, (GLenum) *type,
      (GLvoid *) newlists );
   for (i=0;i<(*length);i++) table ARRAYREF [i]=(GLint)newlists[i];
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*length)*sizeof(GLushort));
   glGetColorTableEXT( (GLenum) *target, (GLenum) *format, (GLenum) *type,
      (GLvoid *) newlists );
   for (i=0;i<(*length);i++) table ARRAYREF [i]=(GLint)newlists[i];
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*length)*sizeof(GLshort));
   glGetColorTableEXT( (GLenum) *target, (GLenum) *format, (GLenum) *type,
      (GLvoid *) newlists );
   for (i=0;i<(*length);i++) table ARRAYREF [i]=(GLint)newlists[i];
   free(newlists);
   break;
  }
 default:
   glGetColorTableEXT( (GLenum) *target, (GLenum) *format, (GLenum) *type,
      (GLvoid *) table ARRAYREF);
   break;
}     
}

void APIENTRY f9y1glGetColorTableEXT( f90glenum target, f90glenum format,
   f90glenum type, f90glvoid ARRAY table, f90glint length )
{
   glGetColorTableEXT( (GLenum) *target, (GLenum) *format, (GLenum) *type,
    (GLvoid *) table ARRAYREF);
}

void APIENTRY f9y4glGetColorTableEXT( f90glenum target, f90glenum format,
   f90glenum type, f90glvoid ARRAY table, f90glint length )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   glGetColorTableEXT( (GLenum) *target, (GLenum) *format, (GLenum) *type,
    (GLvoid **)&ptr );
   p = (unsigned char *) ptr;
   for (i=0; i<(*length); i++) {table ARRAYREF [i]=(GLuint)*p; p++;}
}

void APIENTRY f9y5glGetColorTableEXT( f90glenum target, f90glenum format,
   f90glenum type, f90glvoid ARRAY table, f90glint length )
{
   glGetColorTableEXT( (GLenum) *target, (GLenum) *format, (GLenum) *type,
     (GLvoid *) table ARRAYREF);
}

void APIENTRY f9y6glGetColorTableEXT( f90glenum target, f90glenum format,
   f90glenum type, f90glshort ARRAY table, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   glGetColorTableEXT( (GLenum) *target, (GLenum) *format, (GLenum) *type,
      (GLvoid *) newlists );
   for (i=0;i<(*length);i++) table ARRAYREF [i]=(GLshort)newlists[i];
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   glGetColorTableEXT( (GLenum) *target, (GLenum) *format, (GLenum) *type,
      (GLvoid *) newlists );
   for (i=0;i<(*length);i++) table ARRAYREF [i]=(GLshort)newlists[i];
   free(newlists);
   break;
  }
 default:
   glGetColorTableEXT( (GLenum) *target, (GLenum) *format, (GLenum) *type,
      (GLvoid *) table  ARRAYREF );
   break;
}     
}

void APIENTRY fglGetColorTableParameterfvEXT( f90glenum target,
    f90glenum pname, f90glfloat ARRAY params )
{
   glGetColorTableParameterfvEXT( (GLenum) *target, (GLenum) *pname,
    (GLfloat *) params ARRAYREF);
}

void APIENTRY fglGetColorTableParameterivEXT( f90glenum target,
    f90glenum pname, f90glint ARRAY params )
{
   glGetColorTableParameterivEXT( (GLenum) *target, (GLenum) *pname,
    (GLint *) params ARRAYREF);
}

#endif

/* EXT_point_parameters */

#ifdef GL_EXT_point_parameters

void APIENTRY fglPointParameterfEXT( f90glenum pname, f90glfloat param )
{
   glPointParameterfEXT( (GLenum) *pname, (GLfloat) *param );
}

void APIENTRY fglPointParameterfvEXT( f90glenum pname, f90glfloat ARRAY params )
{
   glPointParameterfvEXT( (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

#endif

/* ARB_point_parameters */

#ifdef GL_ARB_point_parameters

void APIENTRY fglPointParameterfARB( f90glenum pname, f90glfloat param )
{
   glPointParameterfARB( (GLenum) *pname, (GLfloat) *param );
}

void APIENTRY fglPointParameterfvARB( f90glenum pname, f90glfloat ARRAY params )
{
   glPointParameterfvARB( (GLenum) *pname, (GLfloat *) params ARRAYREF);
}

#endif

/* EXT_polygon_offset */

#ifdef GL_EXT_polygon_offset

void APIENTRY fglPolygonOffsetEXT( f90glfloat factor, f90glfloat bias )
{
   glPolygonOffsetEXT( (GLfloat) *factor, (GLfloat) *bias );
}

#endif

/* EXT_subtexture */

#ifdef GL_EXT_subtexture

void APIENTRY f9y0glTexSubImage1DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glsizei width, f90glenum format, f90glenum type,
    f90glvoid ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage1DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage1DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newlists[i]=(GLushort)pixels ARRAYREF [i];
   glTexSubImage1DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newlists[i]=(GLshort)pixels ARRAYREF [i];
   glTexSubImage1DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage1DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y1glTexSubImage1DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glsizei width, f90glenum format, f90glenum type,
    f90glvoid ARRAY pixels )
{
   glTexSubImage1DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glTexSubImage1DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glsizei width, f90glenum format, f90glenum type,
    f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexSubImage1DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glTexSubImage1DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glsizei width, f90glenum format, f90glenum type,
    f90glvoid ARRAY pixels )
{
   glTexSubImage1DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glTexSubImage1DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glsizei width, f90glenum format, f90glenum type,
    f90glshort ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage1DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage1DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage1DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLsizei) *width, (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y0glTexSubImage2DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glsizei width, f90glsizei height,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage2DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage2DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newlists[i]=(GLushort)pixels ARRAYREF [i];
   glTexSubImage2DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newlists[i]=(GLshort)pixels ARRAYREF [i];
   glTexSubImage2DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage2DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y1glTexSubImage2DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glsizei width, f90glsizei height,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexSubImage2DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glTexSubImage2DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glsizei width, f90glsizei height,
    f90glenum format, f90glenum type, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexSubImage2DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glTexSubImage2DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glsizei width, f90glsizei height,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexSubImage2DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glTexSubImage2DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glsizei width, f90glsizei height,
    f90glenum format, f90glenum type, f90glshort ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage2DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage2DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage2DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLsizei) *width, (GLsizei) *height, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

#endif

/* EXT_texture3D */

#ifdef GL_EXT_texture3D

void APIENTRY f9y0glTexImage3DEXT( f90glenum target, f90glint level,
    f90glenum internalFormat, f90glsizei width,
    f90glsizei height, f90glsizei depth, f90glint border,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexImage3DEXT( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexImage3DEXT( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newlists[i]=(GLushort)pixels ARRAYREF [i];
   glTexImage3DEXT( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newlists[i]=(GLshort)pixels ARRAYREF [i];
   glTexImage3DEXT( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexImage3DEXT( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y1glTexImage3DEXT( f90glenum target, f90glint level,
    f90glenum internalFormat, f90glsizei width,
    f90glsizei height, f90glsizei depth, f90glint border,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexImage3DEXT( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glTexImage3DEXT( f90glenum target, f90glint level,
    f90glenum internalFormat, f90glsizei width,
    f90glsizei height, f90glsizei depth, f90glint border,
    f90glenum format, f90glenum type, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexImage3DEXT( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glTexImage3DEXT( f90glenum target, f90glint level,
    f90glenum internalFormat, f90glsizei width,
    f90glsizei height, f90glsizei depth, f90glint border,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexImage3DEXT( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glTexImage3DEXT( f90glenum target, f90glint level,
    f90glenum internalFormat, f90glsizei width,
    f90glsizei height, f90glsizei depth, f90glint border,
    f90glenum format, f90glenum type, f90glshort ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexImage3DEXT( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexImage3DEXT( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexImage3DEXT( (GLenum) *target, (GLint) *level, (GLenum) *internalFormat,
    (GLsizei) *width, (GLsizei) *height, (GLsizei) *depth, (GLint) *border,
    (GLenum) *format, (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y0glTexSubImage3DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glsizei width,
    f90glsizei height, f90glsizei depth,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newlists;
   int i;
   newlists = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newlists[i]=(GLushort)pixels ARRAYREF [i];
   glTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newlists;
   int i;
   newlists = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newlists[i]=(GLshort)pixels ARRAYREF [i];
   glTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY f9y1glTexSubImage3DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glsizei width,
    f90glsizei height, f90glsizei depth,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y4glTexSubImage3DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glsizei width,
    f90glsizei height, f90glsizei depth,
    f90glenum format, f90glenum type, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *pixels;
   unsigned char *p;
   int i;
   p = (unsigned char *) &pixels;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels );
}

void APIENTRY f9y5glTexSubImage3DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glsizei width,
    f90glsizei height, f90glsizei depth,
    f90glenum format, f90glenum type, f90glvoid ARRAY pixels )
{
   glTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels ARRAYREF);
}

void APIENTRY f9y6glTexSubImage3DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glsizei width,
    f90glsizei height, f90glsizei depth,
    f90glenum format, f90glenum type, f90glshort ARRAY pixels, f90glint length )
{

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newlists;
   int i;
   newlists = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLubyte)pixels ARRAYREF [i];
   glTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newlists;
   int i;
   newlists = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newlists[i]=(GLbyte)pixels ARRAYREF [i];
   glTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) newlists );
   free(newlists);
   break;
  }
 default:
   glTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLsizei) *width, (GLsizei) *height,
    (GLsizei) *depth, (GLenum) *format,
    (GLenum) *type, (GLvoid *) pixels  ARRAYREF );
   break;
}     
}

void APIENTRY fglCopyTexSubImage3DEXT( f90glenum target, f90glint level,
    f90glint xoffset, f90glint yoffset, f90glint zoffset, f90glint x,
    f90glint y, f90glsizei width, f90glsizei height )
{
   glCopyTexSubImage3DEXT( (GLenum) *target, (GLint) *level, (GLint) *xoffset,
    (GLint) *yoffset, (GLint) *zoffset, (GLint) *x, (GLint) *y,
    (GLsizei) *width, (GLsizei) *height );
}

#endif

/* EXT_texture_object */

#ifdef GL_EXT_texture_object

void APIENTRY fglGenTexturesEXT( f90glsizei n, f90gluint ARRAY textures)
{
   glGenTexturesEXT( (GLsizei) *n, (GLuint *) textures ARRAYREF);
}

void APIENTRY fglDeleteTexturesEXT( f90glsizei n, f90gluint ARRAY textures )
{
   glDeleteTexturesEXT( (GLsizei) *n, (GLuint *) textures ARRAYREF);
}

void APIENTRY fglBindTextureEXT( f90glenum target, f90gluint texture )
{
   glBindTextureEXT( (GLenum) *target, (GLuint) *texture );
}

void APIENTRY fglPrioritizeTexturesEXT( f90glsizei n, f90gluint ARRAY textures,
    f90glclampf ARRAY priorities )
{
   glPrioritizeTexturesEXT( (GLsizei) *n, (GLuint *) textures ARRAYREF,
    (GLclampf *) priorities  ARRAYREF);
}

#ifdef PURE_FUNCTIONS
void APIENTRY fglAreTexturesResidentEXT( f90glsizei n,
    f90gluint ARRAY textures, f90glboolean ARRAY residences, f90glboolean areresident )
#else
GLboolean APIENTRY fglAreTexturesResidentEXT( f90glsizei n,
    f90gluint ARRAY textures, f90glboolean ARRAY residences )
#endif
{
#ifdef PURE_FUNCTIONS
   *areresident = glAreTexturesResidentEXT( (GLsizei) *n,
    (GLuint *) textures ARRAYREF, (GLboolean *) residences  ARRAYREF);
#else
   return (GLboolean) glAreTexturesResidentEXT( (GLsizei) *n,
    (GLuint *) textures ARRAYREF, (GLboolean *) residences  ARRAYREF);
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9eglAreTexturesResidentEXT(f90glsizei n,
   f90gluint ARRAY textures,f90glint ARRAY residences, f90glint areresident )
#else
GLint APIENTRY f9eglAreTexturesResidentEXT(f90glsizei n,
   f90gluint ARRAY textures,f90glint ARRAY residences )
#endif
{
   GLboolean *boolean_residences;
   GLint i, iresult;
   boolean_residences = (GLboolean *)malloc((*n)*sizeof(GLboolean));
   iresult = (GLint) glAreTexturesResidentEXT( (GLsizei) *n,
      (GLuint *) textures ARRAYREF, boolean_residences);
   for (i=0;i<*n;i++) residences ARRAYREF [i]=(GLint) boolean_residences[i];
   free(boolean_residences);
#ifdef PURE_FUNCTIONS
   *areresident = iresult;
#else
   return iresult;
#endif
}

GLboolean APIENTRY fglIsTextureEXT( f90gluint texture )
{
   return (GLboolean) glIsTextureEXT( (GLuint) *texture );
}

GLint APIENTRY f9eglIsTextureEXT( f90gluint texture )
{
   GLboolean temp;
   temp = glIsTextureEXT( (GLuint) *texture );
   return (GLint) temp;
}

#endif

/* EXT_vertex_array */

/* Apparently the Microsoft OpenGL library does not contain this extension
   but their gl.h file defines GL_EXT_vertex_array.  So don't include
   this if we're on Win32 */

#ifndef WIN32

#ifdef GL_EXT_vertex_array

void APIENTRY f9y0glVertexPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glVertexPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y1glVertexPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glVertexPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y4glVertexPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glVertexPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr );
}

void APIENTRY f9y6glVertexPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glVertexPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y7glVertexPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glVertexPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y0glNormalPointerEXT( f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glNormalPointerEXT( (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y1glNormalPointerEXT( f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glNormalPointerEXT( (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y4glNormalPointerEXT( f90glenum type,
    f90glsizei stride, f90glsizei count, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glNormalPointerEXT( (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr );
}

void APIENTRY f9y5glNormalPointerEXT( f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glNormalPointerEXT( (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y6glNormalPointerEXT( f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glNormalPointerEXT( (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y7glNormalPointerEXT( f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glNormalPointerEXT( (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y0glColorPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glColorPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y1glColorPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glColorPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y4glColorPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glColorPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr );
}

void APIENTRY f9y5glColorPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glColorPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y6glColorPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glColorPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y7glColorPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glColorPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y0glIndexPointerEXT( f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glIndexPointerEXT( (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y1glIndexPointerEXT( f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glIndexPointerEXT( (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y4glIndexPointerEXT( f90glenum type,
    f90glsizei stride, f90glsizei count, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glIndexPointerEXT( (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr );
}

void APIENTRY f9y6glIndexPointerEXT( f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glIndexPointerEXT( (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y7glIndexPointerEXT( f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glIndexPointerEXT( (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y0glTexCoordPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glTexCoordPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y1glTexCoordPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glTexCoordPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y4glTexCoordPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90gluint ARRAY temp, f90glint nbytes )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   p = (unsigned char *) &ptr;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   glTexCoordPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr );
}

void APIENTRY f9y6glTexCoordPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glTexCoordPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY f9y7glTexCoordPointerEXT( f90glint size, f90glenum type,
    f90glsizei stride, f90glsizei count, f90glvoid ARRAY ptr )
{
   glTexCoordPointerEXT( (GLint) *size, (GLenum) *type, (GLsizei) *stride,
    (GLsizei) *count, (GLvoid *) ptr ARRAYREF);
}

void APIENTRY fglEdgeFlagPointerEXT( f90glsizei stride, f90glsizei count,
    f90glboolean ARRAY ptr )
{
   glEdgeFlagPointerEXT( (GLsizei) *stride, (GLsizei) *count,
    (GLboolean *) ptr ARRAYREF);
}

void APIENTRY f9eglEdgeFlagPointerEXT( f90glsizei stride, f90glsizei count,
    f90glbyte ARRAY ptr )
{
   glEdgeFlagPointerEXT( (GLsizei) *stride, (GLsizei) *count,
    (GLboolean *) ptr ARRAYREF);
}

void APIENTRY fglGetPointervEXT( f90glenum pname, f90gluint ARRAY params,
   f90glint length )
{
   GLvoid *ptr;
   unsigned char *p;
   int i;
   glGetPointervEXT( (GLenum) *pname, (GLvoid **)&ptr );
   p = (unsigned char *) ptr;
   for (i=0; i<(*length); i++) {params ARRAYREF [i]=(GLuint)*p; p++;}
}

void APIENTRY fglArrayElementEXT( f90glint i )
{
   glArrayElementEXT( (GLint) *i );
}

void APIENTRY fglDrawArraysEXT( f90glenum mode, f90glint first,
    f90glsizei count )
{
   glDrawArraysEXT( (GLenum) *mode, (GLint) *first, (GLsizei) *count );
}

#endif
#endif /* WIN32 */

/* SGIS_multisample */

#ifdef GL_SGIS_multisample

void APIENTRY fglSampleMaskSGIS( f90glclampf value, f90glboolean invert )
{
   glSampleMaskSGIS( (GLclampf) *value, (GLboolean) *invert );
}

void APIENTRY f9eglSampleMaskSGIS( f90glclampf value, f90glint invert )
{
   glSampleMaskSGIS( (GLclampf) *value, (GLboolean) *invert );
}

void APIENTRY fglSamplePatternSGIS( f90glenum pattern )
{
   glSamplePatternSGIS( (GLenum) *pattern );
}

#endif

/* MESA_resize_buffers */

#ifdef GL_MESA_resize_buffers

void APIENTRY fglResizeBuffersMESA( void )
{
   glResizeBuffersMESA();
}

#endif

/* MESA_window_pos */

#ifdef GL_MESA_window_pos

void APIENTRY fglWindowPos2iMESA( f90glint x, f90glint y )
{
   glWindowPos2iMESA( (GLint) *x, (GLint) *y );
}

void APIENTRY fglWindowPos2sMESA( f90glshort x, f90glshort y )
{
   glWindowPos2sMESA( (GLshort) *x, (GLshort) *y );
}

void APIENTRY f9aglWindowPos2sMESA(f90glint x,f90glint y )
{
   glWindowPos2sMESA( (GLshort) *x, (GLshort) *y);
}

void APIENTRY fglWindowPos2fMESA( f90glfloat x, f90glfloat y )
{
   glWindowPos2fMESA( (GLfloat) *x, (GLfloat) *y );
}

void APIENTRY fglWindowPos2dMESA( f90gldouble x, f90gldouble y )
{
   glWindowPos2dMESA( (GLdouble) *x, (GLdouble) *y );
}

void APIENTRY fglWindowPos2ivMESA( f90glint ARRAY p )
{
   glWindowPos2ivMESA( (GLint *) p ARRAYREF);
}

void APIENTRY fglWindowPos2svMESA( f90glshort ARRAY p )
{
   glWindowPos2svMESA( (GLshort *) p ARRAYREF);
}

void APIENTRY f9aglWindowPos2svMESA(f90glint ARRAY p )
{
   GLshort short_p[2];
   int i;
   for (i=0;i<2;i++) short_p[i]=(GLshort) p ARRAYREF [i];
   glWindowPos2svMESA( short_p);
}

void APIENTRY fglWindowPos2fvMESA( f90glfloat ARRAY p )
{
   glWindowPos2fvMESA( (GLfloat *) p ARRAYREF);
}

void APIENTRY fglWindowPos2dvMESA( f90gldouble ARRAY p )
{
   glWindowPos2dvMESA( (GLdouble *) p ARRAYREF);
}

void APIENTRY fglWindowPos3iMESA( f90glint x, f90glint y, f90glint z )
{
   glWindowPos3iMESA( (GLint) *x, (GLint) *y, (GLint) *z );
}

void APIENTRY fglWindowPos3sMESA( f90glshort x, f90glshort y, f90glshort z )
{
   glWindowPos3sMESA( (GLshort) *x, (GLshort) *y, (GLshort) *z );
}

void APIENTRY f9aglWindowPos3sMESA(f90glint x,f90glint y,f90glint z )
{
   glWindowPos3sMESA( (GLshort) *x, (GLshort) *y, (GLshort) *z);
}

void APIENTRY fglWindowPos3fMESA( f90glfloat x, f90glfloat y, f90glfloat z )
{
   glWindowPos3fMESA( (GLfloat) *x, (GLfloat) *y, (GLfloat) *z );
}

void APIENTRY fglWindowPos3dMESA( f90gldouble x, f90gldouble y,
    f90gldouble z )
{
   glWindowPos3dMESA( (GLdouble) *x, (GLdouble) *y, (GLdouble) *z );
}

void APIENTRY fglWindowPos3ivMESA( f90glint ARRAY p )
{
   glWindowPos3ivMESA( (GLint *) p ARRAYREF);
}

void APIENTRY fglWindowPos3svMESA( f90glshort ARRAY p )
{
   glWindowPos3svMESA( (GLshort *) p ARRAYREF);
}

void APIENTRY f9aglWindowPos3svMESA(f90glint ARRAY p )
{
   GLshort short_p[3];
   int i;
   for (i=0;i<3;i++) short_p[i]=(GLshort) p ARRAYREF [i];
   glWindowPos3svMESA( short_p);
}

void APIENTRY fglWindowPos3fvMESA( f90glfloat ARRAY p )
{
   glWindowPos3fvMESA( (GLfloat *) p ARRAYREF);
}

void APIENTRY fglWindowPos3dvMESA( f90gldouble ARRAY p )
{
   glWindowPos3dvMESA( (GLdouble *) p ARRAYREF);
}

void APIENTRY fglWindowPos4iMESA( f90glint x, f90glint y, f90glint z,
    f90glint w )
{
   glWindowPos4iMESA( (GLint) *x, (GLint) *y, (GLint) *z, (GLint) *w );
}

void APIENTRY fglWindowPos4sMESA( f90glshort x, f90glshort y, f90glshort z,
    f90glshort w )
{
   glWindowPos4sMESA( (GLshort) *x, (GLshort) *y, (GLshort) *z,
    (GLshort) *w );
}

void APIENTRY f9aglWindowPos4sMESA(f90glint x,f90glint y,f90glint z,
   f90glint w )
{
   glWindowPos4sMESA( (GLshort) *x, (GLshort) *y, (GLshort) *z,
    (GLshort) *w);
}

void APIENTRY fglWindowPos4fMESA( f90glfloat x, f90glfloat y, f90glfloat z,
    f90glfloat w )
{
   glWindowPos4fMESA( (GLfloat) *x, (GLfloat) *y, (GLfloat) *z,
    (GLfloat) *w );
}

void APIENTRY fglWindowPos4dMESA( f90gldouble x, f90gldouble y, f90gldouble z,
    f90gldouble w )
{
   glWindowPos4dMESA( (GLdouble) *x, (GLdouble) *y, (GLdouble) *z,
    (GLdouble) *w );
}

void APIENTRY fglWindowPos4ivMESA( f90glint ARRAY p )
{
   glWindowPos4ivMESA( (GLint *) p ARRAYREF);
}

void APIENTRY fglWindowPos4svMESA( f90glshort ARRAY p )
{
   glWindowPos4svMESA( (GLshort *) p ARRAYREF);
}

void APIENTRY f9aglWindowPos4svMESA(f90glint ARRAY p )
{
   GLshort short_p[4];
   int i;
   for (i=0;i<4;i++) short_p[i]=(GLshort) p ARRAYREF [i];
   glWindowPos4svMESA( short_p);
}

void APIENTRY fglWindowPos4fvMESA( f90glfloat ARRAY p )
{
   glWindowPos4fvMESA( (GLfloat *) p ARRAYREF);
}

void APIENTRY fglWindowPos4dvMESA( f90gldouble ARRAY p )
{
   glWindowPos4dvMESA( (GLdouble *) p ARRAYREF);
}

#endif
