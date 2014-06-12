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


#include 'fppr.inc'

module opengl_fwrap

use opengl_kinds
use opengl_glinterfaces
implicit none
private

#ifdef NASOFTWARE
#define INTENT_IN ::
#else
#define INTENT_IN ,intent(in) ::
#endif

#ifdef F_LANG
! F does not support assumed size arrays
#define DEFER_DIM :
#else
#define DEFER_DIM *
#endif

#ifdef F_LANG
#define PURE_FUNCTIONS 
#endif

!--------------------------
! Public generic interfaces
!--------------------------

! OpenGL 1.0

public :: &
   glaccum, &
   glalphafunc, &
   glbegin, &
   glbitmap, &
   glblendfunc, &
   glcalllist, &
   glcalllists, &
   glclear, &
   glclearaccum, &
   glclearcolor, &
   glcleardepth, &
   glclearindex, &
   glclearstencil, &
   glclipplane, &
   glcolor3b, &
   glcolor3bv, &
   glcolor3d, &
   glcolor3dv, &
   glcolor3f, &
   glcolor3fv, &
   glcolor3i, &
   glcolor3iv, &
   glcolor3s, &
   glcolor3sv, &
   glcolor3ub, &
   glcolor3ubv, &
   glcolor3ui, &
   glcolor3uiv, &
   glcolor3us, &
   glcolor3usv, &
   glcolor4b, &
   glcolor4bv, &
   glcolor4d, &
   glcolor4dv, &
   glcolor4f, &
   glcolor4fv
public :: &
   glcolor4i, &
   glcolor4iv, &
   glcolor4s, &
   glcolor4sv, &
   glcolor4ub, &
   glcolor4ubv, &
   glcolor4ui, &
   glcolor4uiv, &
   glcolor4us, &
   glcolor4usv, &
   glcolormask, &
   glcolormaterial, &
   glcopypixels, &
   glcullface, &
   gldeletelists, &
   gldepthfunc, &
   gldepthmask, &
   gldepthrange, &
   gldisable, &
   gldrawbuffer, &
   gldrawpixels, &
   gledgeflag, &
   gledgeflagv, &
   glenable, &
   glend, &
   glendlist, &
   glevalcoord1d, &
   glevalcoord1dv, &
   glevalcoord1f, &
   glevalcoord1fv, &
   glevalcoord2d, &
   glevalcoord2dv, &
   glevalcoord2f, &
   glevalcoord2fv, &
   glevalmesh1, &
   glevalmesh2
public :: &
   glevalpoint1, &
   glevalpoint2, &
   glfeedbackbuffer, &
   glfinish, &
   glflush, &
   glfogf, &
   glfogfv, &
   glfogi, &
   glfogiv, &
   glfrontface, &
   glfrustum, &
   glgenlists, &
   glgetbooleanv, &
   glgetclipplane, &
   glgetdoublev, &
   glgeterror, &
   glgetfloatv, &
   glgetintegerv, &
   glgetlightfv, &
   glgetlightiv, &
   glgetmapdv, &
   glgetmapfv, &
   glgetmapiv, &
   glgetmaterialfv, &
   glgetmaterialiv, &
   glgetpixelmapfv, &
   glgetpixelmapuiv, &
   glgetpixelmapusv, &
   glgetpolygonstipple, &
   glgetstring, &
   glgettexenvfv, &
   glgettexenviv, &
   glgettexgendv, &
   glgettexgenfv, &
   glgettexgeniv, &
   glgetteximage
public :: &
   glgettexlevelparameterfv, &
   glgettexlevelparameteriv, &
   glgettexparameterfv, &
   glgettexparameteriv, &
   glhint, &
   glindexmask, &
   glindexd, &
   glindexdv, &
   glindexf, &
   glindexfv, &
   glindexi, &
   glindexiv, &
   glindexs, &
   glindexsv, &
   glinitnames, &
   glisenabled, &
   glislist, &
   gllightmodelf, &
   gllightmodelfv, &
   gllightmodeli, &
   gllightmodeliv, &
   gllightf, &
   gllightfv, &
   gllighti, &
   gllightiv, &
   gllinestipple, &
   gllinewidth, &
   gllistbase, &
   glloadidentity, &
   glloadmatrixd, &
   glloadmatrixf, &
   glloadname, &
   gllogicop, &
   glmap1d, &
   glmap1f, &
   glmap2d
public :: &
   glmap2f, &
   glmapgrid1d, &
   glmapgrid1f, &
   glmapgrid2d, &
   glmapgrid2f, &
   glmaterialf, &
   glmaterialfv, &
   glmateriali, &
   glmaterialiv, &
   glmatrixmode, &
   glmultmatrixd, &
   glmultmatrixf
public :: &
   glnewlist, &
   glnormal3b, &
   glnormal3bv, &
   glnormal3d, &
   glnormal3dv, &
   glnormal3f, &
   glnormal3fv, &
   glnormal3i, &
   glnormal3iv, &
   glnormal3s, &
   glnormal3sv, &
   glortho, &
   glpassthrough, &
   glpixelmapfv, &
   glpixelmapuiv, &
   glpixelmapusv, &
   glpixelstoref, &
   glpixelstorei, &
   glpixeltransferf, &
   glpixeltransferi, &
   glpixelzoom, &
   glpointsize, &
   glpolygonmode, &
   glpolygonstipple
public :: &
   glpopattrib, &
   glpopmatrix, &
   glpopname, &
   glpushattrib, &
   glpushmatrix, &
   glpushname, &
   glrasterpos2d, &
   glrasterpos2dv, &
   glrasterpos2f, &
   glrasterpos2fv, &
   glrasterpos2i, &
   glrasterpos2iv, &
   glrasterpos2s, &
   glrasterpos2sv, &
   glrasterpos3d, &
   glrasterpos3dv, &
   glrasterpos3f, &
   glrasterpos3fv, &
   glrasterpos3i, &
   glrasterpos3iv, &
   glrasterpos3s, &
   glrasterpos3sv, &
   glrasterpos4d, &
   glrasterpos4dv, &
   glrasterpos4f, &
   glrasterpos4fv, &
   glrasterpos4i, &
   glrasterpos4iv, &
   glrasterpos4s, &
   glrasterpos4sv, &
   glreadbuffer, &
   glreadpixels, &
   glrectd, &
   glrectdv, &
   glrectf, &
   glrectfv
public :: &
   glrecti, &
   glrectiv, &
   glrects, &
   glrectsv, &
   glrendermode, &
   glrotated, &
   glrotatef, &
   glscaled, &
   glscalef, &
   glscissor, &
   glselectbuffer, &
   glshademodel, &
   glstencilfunc, &
   glstencilmask, &
   glstencilop, &
   gltexcoord1d, &
   gltexcoord1dv, &
   gltexcoord1f, &
   gltexcoord1fv, &
   gltexcoord1i, &
   gltexcoord1iv, &
   gltexcoord1s, &
   gltexcoord1sv, &
   gltexcoord2d, &
   gltexcoord2dv, &
   gltexcoord2f, &
   gltexcoord2fv, &
   gltexcoord2i, &
   gltexcoord2iv, &
   gltexcoord2s, &
   gltexcoord2sv, &
   gltexcoord3d, &
   gltexcoord3dv, &
   gltexcoord3f, &
   gltexcoord3fv, &
   gltexcoord3i
public :: &
   gltexcoord3iv, &
   gltexcoord3s, &
   gltexcoord3sv, &
   gltexcoord4d, &
   gltexcoord4dv, &
   gltexcoord4f, &
   gltexcoord4fv, &
   gltexcoord4i, &
   gltexcoord4iv, &
   gltexcoord4s, &
   gltexcoord4sv, &
   gltexenvf, &
   gltexenvfv, &
   gltexenvi, &
   gltexenviv, &
   gltexgend, &
   gltexgendv, &
   gltexgenf, &
   gltexgenfv, &
   gltexgeni, &
   gltexgeniv, &
   glteximage1d, &
   glteximage2d, &
   gltexparameterf, &
   gltexparameterfv, &
   gltexparameteri, &
   gltexparameteriv, &
   gltranslated, &
   gltranslatef, &
   glvertex2d, &
   glvertex2dv, &
   glvertex2f, &
   glvertex2fv, &
   glvertex2i, &
   glvertex2iv, &
   glvertex2s
public :: &
   glvertex2sv, &
   glvertex3d, &
   glvertex3dv, &
   glvertex3f, &
   glvertex3fv, &
   glvertex3i, &
   glvertex3iv, &
   glvertex3s, &
   glvertex3sv, &
   glvertex4d, &
   glvertex4dv, &
   glvertex4f, &
   glvertex4fv, &
   glvertex4i, &
   glvertex4iv, &
   glvertex4s, &
   glvertex4sv, &
   glviewport

! OpenGL 1.1

#ifdef OPENGL_1_1
public :: &
   glaretexturesresident, &
   glarrayelement, &
   glbindtexture, &
   glcolorpointer, &
   glcopyteximage1d, &
   glcopyteximage2d, &
   glcopytexsubimage1d, &
   glcopytexsubimage2d, &
   gldeletetextures, &
   gldisableclientstate, &
   gldrawarrays, &
   gldrawelements, &
   gledgeflagpointer, &
   glenableclientstate, &
   glgentextures, &
   glgetpointerv, &
   glindexpointer, &
   glindexub, &
   glindexubv, &
   glinterleavedarrays, &
   glistexture, &
   glnormalpointer, &
   glpolygonoffset, &
   glpopclientattrib, &
   glprioritizetextures, &
   glpushclientattrib, &
   gltexcoordpointer, &
   gltexsubimage1d, &
   gltexsubimage2d, &
   glvertexpointer
#endif

! OpenGL 1.2

#ifdef OPENGL_1_2
public :: &
   gldrawrangeelements, &
   glteximage3d, &
   gltexsubimage3d, &
   glcopytexsubimage3d
#endif

! EXT_blend_color

#ifdef GL_EXT_blend_color
public :: &
   glblendcolorext
#endif

! EXT_blend_minmax

#ifdef GL_EXT_blend_minmax
public :: &
   glblendequationext
#endif

! EXT_color_table

#ifdef GL_EXT_color_table
public :: &
   glcolortableext, &
   glcolorsubtableext, &
   glgetcolortableext, &
   glgetcolortableparameterfvext, &
   glgetcolortableparameterivext
#endif

! EXT_point_parameters

#ifdef GL_EXT_point_parameters
public :: &
   glpointparameterfext, &
   glpointparameterfvext
#endif

! EXT_polygon_offset

#ifdef GL_EXT_polygon_offset
public :: &
   glpolygonoffsetext
#endif

! EXT_subtexture

#ifdef GL_EXT_subtexture
public :: &
   gltexsubimage1dext, &
   gltexsubimage2dext
#endif

! EXT_texture3D

#ifdef GL_EXT_texture3D
public :: &
   glteximage3dext, &
   gltexsubimage3dext, &
   glcopytexsubimage3dext
#endif

! EXT_texture_object

#ifdef GL_EXT_texture_object
public :: &
   glgentexturesext, &
   gldeletetexturesext, &
   glbindtextureext, &
   glprioritizetexturesext, &
   glaretexturesresidentext, &
   glistextureext
#endif

! EXT_vertex_array

#ifdef GL_EXT_vertex_array
public :: &
   glvertexpointerext, &
   glnormalpointerext, &
   glcolorpointerext, &
   glindexpointerext, &
   gltexcoordpointerext, &
   gledgeflagpointerext, &
   glgetpointervext, &
   glarrayelementext, &
   gldrawarraysext
#endif

! SGIS_multisample

#ifdef GL_SGIS_multisample
public :: &
   glsamplemasksgis, &
   glsamplepatternsgis
#endif

! MESA_resize_buffers

#ifdef GL_MESA_resize_buffers
public :: &
   glresizebuffersmesa
#endif

! MESA_window_pos

#ifdef GL_MESA_window_pos
public :: &
   glwindowpos2imesa, &
   glwindowpos2smesa, &
   glwindowpos2fmesa, &
   glwindowpos2dmesa, &
   glwindowpos2ivmesa, &
   glwindowpos2svmesa, &
   glwindowpos2fvmesa, &
   glwindowpos2dvmesa, &
   glwindowpos3imesa, &
   glwindowpos3smesa, &
   glwindowpos3fmesa, &
   glwindowpos3dmesa, &
   glwindowpos3ivmesa, &
   glwindowpos3svmesa, &
   glwindowpos3fvmesa, &
   glwindowpos3dvmesa, &
   glwindowpos4imesa, &
   glwindowpos4smesa, &
   glwindowpos4fmesa, &
   glwindowpos4dmesa, &
   glwindowpos4ivmesa, &
   glwindowpos4svmesa, &
   glwindowpos4fvmesa, &
   glwindowpos4dvmesa
#endif

!-----------------------------------
! Public entities from other modules
!-----------------------------------

! Some older compilers require that derived xtypes from other modules be made
! public if there are public routines with arguments of that xtype.
! F forbids this, and some other compilers issue warnings about double access.
! Here we make them public except on those compilers that complain.

#ifdef F_LANG
#define DONT_PASS_TYPES
#endif

#ifndef DONT_PASS_TYPES
public :: glcptr
#endif

!------------------------
! Private module routines
!------------------------

! OpenGL 1.0

private :: &
   f9xglaccum, &
   f9xglalphafunc, &
   f9xglbegin, &
   f9xglbitmap, &
   f9xglblendfunc, &
   f9xglcalllist, &
   f9x0glcalllists, &
   f9x1glcalllists, &
   f9x4glcalllists, &
   f9x5glcalllists, &
   f9x6glcalllists, &
   f9xglclear, &
   f9xglclearaccum, &
   f9xglclearcolor, &
   f9xglcleardepth, &
   f9xglclearindex, &
   f9xglclearstencil, &
   f9xglclipplane, &
   f9xglcolor3b, &
   f9xglcolor3bv, &
   f9xglcolor3d, &
   f9xglcolor3dv, &
   f9xglcolor3f, &
   f9xglcolor3fv, &
   f9xglcolor3i, &
   f9xglcolor3iv, &
   f9xglcolor3s, &
   f9xglcolor3sv, &
   f9xglcolor3ub, &
   f9xglcolor3ubv, &
   f9xglcolor3ui, &
   f9xglcolor3uiv, &
   f9xglcolor3us, &
   f9xglcolor3usv, &
   f9xglcolor4b, &
   f9xglcolor4bv
private :: &
   f9xglcolor4d, &
   f9xglcolor4dv, &
   f9xglcolor4f, &
   f9xglcolor4fv, &
   f9xglcolor4i, &
   f9xglcolor4iv, &
   f9xglcolor4s, &
   f9xglcolor4sv, &
   f9xglcolor4ub, &
   f9xglcolor4ubv, &
   f9xglcolor4ui, &
   f9xglcolor4uiv, &
   f9xglcolor4us, &
   f9xglcolor4usv, &
   f9xglcolormask, &
   f9xglcolormaterial, &
   f9xglcopypixels, &
   f9xglcullface, &
   f9xgldeletelists, &
   f9xgldepthfunc, &
   f9xgldepthmask, &
   f9xgldepthrange, &
   f9xgldisable, &
   f9xgldrawbuffer, &
   f9x0gldrawpixels, &
   f9x1gldrawpixels, &
   f9x4gldrawpixels, &
   f9x5gldrawpixels, &
   f9x6gldrawpixels, &
   f9xgledgeflag, &
   f9xgledgeflagv, &
   f9xglenable, &
   f9xglend, &
   f9xglendlist, &
   f9xglevalcoord1d
private :: &
   f9xglevalcoord1dv, &
   f9xglevalcoord1f, &
   f9xglevalcoord1fv, &
   f9xglevalcoord2d, &
   f9xglevalcoord2dv, &
   f9xglevalcoord2f, &
   f9xglevalcoord2fv, &
   f9xglevalmesh1, &
   f9xglevalmesh2, &
   f9xglevalpoint1, &
   f9xglevalpoint2, &
   f9xglfeedbackbuffer, &
   f9xglfinish, &
   f9xglflush, &
   f9xglfogf, &
   f9xglfogfv, &
   f9xglfogi, &
   f9xglfogiv, &
   f9xglfrontface, &
   f9xglfrustum, &
   f9xglgenlists, &
   f9xglgetbooleanv, &
   f9xglgetclipplane, &
   f9xglgetdoublev, &
   f9xglgeterror, &
   f9xglgetfloatv, &
   f9xglgetintegerv, &
   f9xglgetlightfv, &
   f9xglgetlightiv, &
   f9xglgetmapdv, &
   f9xglgetmapfv, &
   f9xglgetmapiv, &
   f9xglgetmaterialfv, &
   f9xglgetmaterialiv
private :: &
   f9xglgetpixelmapfv, &
   f9xglgetpixelmapuiv, &
   f9xglgetpixelmapusv, &
   f9xglgetpolygonstipple, &
   f9xglgetstring, &
   f9xglgettexenvfv, &
   f9xglgettexenviv, &
   f9xglgettexgendv, &
   f9xglgettexgenfv, &
   f9xglgettexgeniv, &
   f9x0glgetteximage, &
   f9x1glgetteximage, &
   f9x4glgetteximage, &
   f9x5glgetteximage, &
   f9x6glgetteximage, &
   f9xglgettexlevelparameterfv, &
   f9xglgettexlevelparameteriv, &
   f9xglgettexparameterfv, &
   f9xglgettexparameteriv, &
   f9xglhint, &
   f9xglindexmask, &
   f9xglindexd, &
   f9xglindexdv, &
   f9xglindexf, &
   f9xglindexfv, &
   f9xglindexi, &
   f9xglindexiv, &
   f9xglindexs, &
   f9xglindexsv, &
   f9xglinitnames, &
   f9xglisenabled, &
   f9xglislist, &
   f9xgllightmodelf, &
   f9xgllightmodelfv, &
   f9xgllightmodeli, &
   f9xgllightmodeliv
private :: &
   f9xgllightf, &
   f9xgllightfv, &
   f9xgllighti, &
   f9xgllightiv, &
   f9xgllinestipple, &
   f9xgllinewidth, &
   f9xgllistbase, &
   f9xglloadidentity, &
   f9xglloadmatrixd, &
   f9xglloadmatrixf, &
   f9xglloadname, &
   f9xgllogicop, &
   f9xglmap1d, &
   f9xglmap1f, &
   f9xglmap2d, &
   f9xglmap2f, &
   f9xglmapgrid1d, &
   f9xglmapgrid1f, &
   f9xglmapgrid2d, &
   f9xglmapgrid2f, &
   f9xglmaterialf, &
   f9xglmaterialfv, &
   f9xglmateriali, &
   f9xglmaterialiv, &
   f9xglmatrixmode, &
   f9xglmultmatrixd, &
   f9xglmultmatrixf
private :: &
   f9xglnewlist, &
   f9xglnormal3b, &
   f9xglnormal3bv, &
   f9xglnormal3d, &
   f9xglnormal3dv, &
   f9xglnormal3f, &
   f9xglnormal3fv, &
   f9xglnormal3i, &
   f9xglnormal3iv
private :: &
   f9xglnormal3s, &
   f9xglnormal3sv, &
   f9xglortho, &
   f9xglpassthrough, &
   f9xglpixelmapfv, &
   f9xglpixelmapuiv, &
   f9xglpixelmapusv, &
   f9xglpixelstoref, &
   f9xglpixelstorei, &
   f9xglpixeltransferf, &
   f9xglpixeltransferi, &
   f9xglpixelzoom, &
   f9xglpointsize, &
   f9xglpolygonmode, &
   f9xglpolygonstipple, &
   f9xglpopattrib, &
   f9xglpopmatrix, &
   f9xglpopname, &
   f9xglpushattrib, &
   f9xglpushmatrix, &
   f9xglpushname, &
   f9xglrasterpos2d, &
   f9xglrasterpos2dv, &
   f9xglrasterpos2f, &
   f9xglrasterpos2fv, &
   f9xglrasterpos2i, &
   f9xglrasterpos2iv, &
   f9xglrasterpos2s, &
   f9xglrasterpos2sv, &
   f9xglrasterpos3d, &
   f9xglrasterpos3dv, &
   f9xglrasterpos3f, &
   f9xglrasterpos3fv, &
   f9xglrasterpos3i, &
   f9xglrasterpos3iv
private :: &
   f9xglrasterpos3s, &
   f9xglrasterpos3sv, &
   f9xglrasterpos4d, &
   f9xglrasterpos4dv, &
   f9xglrasterpos4f, &
   f9xglrasterpos4fv, &
   f9xglrasterpos4i, &
   f9xglrasterpos4iv, &
   f9xglrasterpos4s, &
   f9xglrasterpos4sv, &
   f9xglreadbuffer, &
   f9x0glreadpixels, &
   f9x1glreadpixels, &
   f9x4glreadpixels, &
   f9x5glreadpixels, &
   f9x6glreadpixels, &
   f9xglrectd, &
   f9xglrectdv, &
   f9xglrectf, &
   f9xglrectfv, &
   f9xglrecti, &
   f9xglrectiv, &
   f9xglrects, &
   f9xglrectsv, &
   f9xglrendermode, &
   f9xglrotated, &
   f9xglrotatef, &
   f9xglscaled, &
   f9xglscalef, &
   f9xglscissor, &
   f9xglselectbuffer, &
   f9xglshademodel, &
   f9xglstencilfunc, &
   f9xglstencilmask, &
   f9xglstencilop
private :: &
   f9xgltexcoord1d, &
   f9xgltexcoord1dv, &
   f9xgltexcoord1f, &
   f9xgltexcoord1fv, &
   f9xgltexcoord1i, &
   f9xgltexcoord1iv, &
   f9xgltexcoord1s, &
   f9xgltexcoord1sv, &
   f9xgltexcoord2d, &
   f9xgltexcoord2dv, &
   f9xgltexcoord2f, &
   f9xgltexcoord2fv, &
   f9xgltexcoord2i, &
   f9xgltexcoord2iv, &
   f9xgltexcoord2s, &
   f9xgltexcoord2sv, &
   f9xgltexcoord3d, &
   f9xgltexcoord3dv, &
   f9xgltexcoord3f, &
   f9xgltexcoord3fv, &
   f9xgltexcoord3i, &
   f9xgltexcoord3iv, &
   f9xgltexcoord3s, &
   f9xgltexcoord3sv, &
   f9xgltexcoord4d, &
   f9xgltexcoord4dv, &
   f9xgltexcoord4f, &
   f9xgltexcoord4fv, &
   f9xgltexcoord4i, &
   f9xgltexcoord4iv, &
   f9xgltexcoord4s, &
   f9xgltexcoord4sv
private :: &
   f9xgltexenvf, &
   f9xgltexenvfv, &
   f9xgltexenvi, &
   f9xgltexenviv, &
   f9xgltexgend, &
   f9xgltexgendv, &
   f9xgltexgenf, &
   f9xgltexgenfv, &
   f9xgltexgeni, &
   f9xgltexgeniv, &
   f9x0glteximage1d, &
   f9x1glteximage1d, &
   f9x4glteximage1d, &
   f9x5glteximage1d, &
   f9x6glteximage1d, &
   f9x0glteximage2d, &
   f9x1glteximage2d, &
   f9x4glteximage2d, &
   f9x5glteximage2d, &
   f9x6glteximage2d, &
   f9xgltexparameterf, &
   f9xgltexparameterfv, &
   f9xgltexparameteri, &
   f9xgltexparameteriv, &
   f9xgltranslated, &
   f9xgltranslatef, &
   f9xglvertex2d, &
   f9xglvertex2dv, &
   f9xglvertex2f, &
   f9xglvertex2fv, &
   f9xglvertex2i, &
   f9xglvertex2iv, &
   f9xglvertex2s
private :: &
   f9xglvertex2sv, &
   f9xglvertex3d, &
   f9xglvertex3dv, &
   f9xglvertex3f, &
   f9xglvertex3fv, &
   f9xglvertex3i, &
   f9xglvertex3iv, &
   f9xglvertex3s, &
   f9xglvertex3sv, &
   f9xglvertex4d, &
   f9xglvertex4dv, &
   f9xglvertex4f, &
   f9xglvertex4fv, &
   f9xglvertex4i, &
   f9xglvertex4iv, &
   f9xglvertex4s, &
   f9xglvertex4sv, &
   f9xglviewport

! OpenGL 1.1

#ifdef OPENGL_1_1
private :: &
   f9xglaretexturesresident, &
   f9xglarrayelement, &
   f9xglbindtexture, &
   f9x0glcolorpointer, &
   f9x1glcolorpointer, &
   f9x5glcolorpointer, &
   f9x6glcolorpointer, &
   f9x7glcolorpointer, &
   f9x4glcolorpointer, &
   f9xglcopyteximage1d, &
   f9xglcopyteximage2d, &
   f9xglcopytexsubimage1d, &
   f9xglcopytexsubimage2d, &
   f9xgldeletetextures, &
   f9xgldisableclientstate, &
   f9xgldrawarrays, &
   f9x0gldrawelements, &
   f9x4gldrawelements, &
   f9x5gldrawelements, &
   f9x6gldrawelements, &
   f9xgledgeflagpointer, &
   f9xegledgeflagpointer, &
   f9xglenableclientstate, &
   f9xglgentextures, &
   f9xglgetpointerv, &
   f9x0glindexpointer, &
   f9x1glindexpointer, &
   f9x6glindexpointer, &
   f9x7glindexpointer, &
   f9x4glindexpointer, &
   f9xglindexub, &
   f9xglindexubv, &
   f9x0glinterleavedarrays, &
   f9x4glinterleavedarrays, &
   f9x1glinterleavedarrays, &
   f9x5glinterleavedarrays, &
   f9x6glinterleavedarrays, &
   f9xglistexture
private :: &
   f9x0glnormalpointer, &
   f9x1glnormalpointer, &
   f9x5glnormalpointer, &
   f9x6glnormalpointer, &
   f9x7glnormalpointer, &
   f9x4glnormalpointer, &
   f9xglpolygonoffset, &
   f9xglpopclientattrib, &
   f9xglprioritizetextures, &
   f9xglpushclientattrib, &
   f9x0gltexcoordpointer, &
   f9x1gltexcoordpointer, &
   f9x6gltexcoordpointer, &
   f9x7gltexcoordpointer, &
   f9x4gltexcoordpointer, &
   f9x0gltexsubimage1d, &
   f9x1gltexsubimage1d, &
   f9x4gltexsubimage1d, &
   f9x5gltexsubimage1d, &
   f9x6gltexsubimage1d, &
   f9x0gltexsubimage2d, &
   f9x1gltexsubimage2d, &
   f9x4gltexsubimage2d, &
   f9x5gltexsubimage2d, &
   f9x6gltexsubimage2d, &
   f9x0glvertexpointer, &
   f9x1glvertexpointer, &
   f9x6glvertexpointer, &
   f9x7glvertexpointer, &
   f9x4glvertexpointer
#endif

! OpenGL 1.2

#ifdef OPENGL_1_2
private :: &
   f9x0gldrawrangeelements, &
   f9x4gldrawrangeelements, &
   f9x5gldrawrangeelements, &
   f9x6gldrawrangeelements, &
   f9x0glteximage3d, &
   f9x1glteximage3d, &
   f9x4glteximage3d, &
   f9x5glteximage3d, &
   f9x6glteximage3d, &
   f9x0gltexsubimage3d, &
   f9x1gltexsubimage3d, &
   f9x4gltexsubimage3d, &
   f9x5gltexsubimage3d, &
   f9x6gltexsubimage3d, &
   f9xglcopytexsubimage3d
#endif

! EXT_blend_color

#ifdef GL_EXT_blend_color
private :: &
   f9xglblendcolorext
#endif

! EXT_blend_minmax

#ifdef GL_EXT_blend_minmax
private :: &
   f9xglblendequationext
#endif

! EXT_color_table

#ifdef GL_EXT_color_table
private :: &
   f9x0glcolortableext, &
   f9x1glcolortableext, &
   f9x4glcolortableext, &
   f9x5glcolortableext, &
   f9x6glcolortableext, &
   f9x0glcolorsubtableext, &
   f9x1glcolorsubtableext, &
   f9x4glcolorsubtableext, &
   f9x5glcolorsubtableext, &
   f9x6glcolorsubtableext, &
   f9x0glgetcolortableext, &
   f9x1glgetcolortableext, &
   f9x4glgetcolortableext, &
   f9x5glgetcolortableext, &
   f9x6glgetcolortableext, &
   f9xglgetcolortableparamfvext, &
   f9xglgetcolortableparamivext
#endif

! EXT_point_parameters

#ifdef GL_EXT_point_parameters
private :: &
   f9xglpointparameterfext, &
   f9xglpointparameterfvext
#endif

! EXT_polygon_offset

#ifdef GL_EXT_polygon_offset
private :: &
   f9xglpolygonoffsetext
#endif

! EXT_subtexture

#ifdef GL_EXT_subtexture
private :: &
   f9x0gltexsubimage1dext, &
   f9x1gltexsubimage1dext, &
   f9x4gltexsubimage1dext, &
   f9x5gltexsubimage1dext, &
   f9x6gltexsubimage1dext, &
   f9x0gltexsubimage2dext, &
   f9x1gltexsubimage2dext, &
   f9x4gltexsubimage2dext, &
   f9x5gltexsubimage2dext, &
   f9x6gltexsubimage2dext
#endif

! EXT_texture3D

#ifdef GL_EXT_texture3D
private :: &
   f9x0glteximage3dext, &
   f9x1glteximage3dext, &
   f9x4glteximage3dext, &
   f9x5glteximage3dext, &
   f9x6glteximage3dext, &
   f9x0gltexsubimage3dext, &
   f9x1gltexsubimage3dext, &
   f9x4gltexsubimage3dext, &
   f9x5gltexsubimage3dext, &
   f9x6gltexsubimage3dext, &
   f9xglcopytexsubimage3dext
#endif

! EXT_texture_object

#ifdef GL_EXT_texture_object
private :: &
   f9xglgentexturesext, &
   f9xgldeletetexturesext, &
   f9xglbindtextureext, &
   f9xglprioritizetexturesext, &
   f9xglaretexturesresidentext, &
   f9xglistextureext
#endif

! EXT_vertex_array

#ifdef GL_EXT_vertex_array
private :: &
   f9x0glvertexpointerext, &
   f9x1glvertexpointerext, &
   f9x6glvertexpointerext, &
   f9x7glvertexpointerext, &
   f9x0glnormalpointerext, &
   f9x1glnormalpointerext, &
   f9x5glnormalpointerext, &
   f9x6glnormalpointerext, &
   f9x7glnormalpointerext, &
   f9x0glcolorpointerext, &
   f9x1glcolorpointerext, &
   f9x5glcolorpointerext, &
   f9x6glcolorpointerext, &
   f9x7glcolorpointerext, &
   f9x0glindexpointerext, &
   f9x1glindexpointerext, &
   f9x6glindexpointerext, &
   f9x7glindexpointerext, &
   f9x0gltexcoordpointerext, &
   f9x1gltexcoordpointerext, &
   f9x6gltexcoordpointerext, &
   f9x7gltexcoordpointerext, &
   f9xgledgeflagpointerext, &
   f9xegledgeflagpointerext, &
   f9x4glvertexpointerext, &
   f9x4glnormalpointerext, &
   f9x4glcolorpointerext, &
   f9x4glindexpointerext, &
   f9x4gltexcoordpointerext, &
   f9xglgetpointervext, &
   f9xglarrayelementext, &
   f9xgldrawarraysext
#endif

! SGIS_multisample

#ifdef GL_SGIS_multisample
private :: &
   f9xglsamplemasksgis, &
   f9xglsamplepatternsgis
#endif

! MESA_resize_buffers

#ifdef GL_MESA_resize_buffers
private :: &
   f9xglresizebuffersmesa
#endif

! MESA_window_pos

#ifdef GL_MESA_window_pos
private :: &
   f9xglwindowpos2imesa, &
   f9xglwindowpos2smesa, &
   f9xglwindowpos2fmesa, &
   f9xglwindowpos2dmesa, &
   f9xglwindowpos2ivmesa, &
   f9xglwindowpos2svmesa, &
   f9xglwindowpos2fvmesa, &
   f9xglwindowpos2dvmesa, &
   f9xglwindowpos3imesa, &
   f9xglwindowpos3smesa, &
   f9xglwindowpos3fmesa, &
   f9xglwindowpos3dmesa, &
   f9xglwindowpos3ivmesa, &
   f9xglwindowpos3svmesa, &
   f9xglwindowpos3fvmesa, &
   f9xglwindowpos3dvmesa, &
   f9xglwindowpos4imesa, &
   f9xglwindowpos4smesa, &
   f9xglwindowpos4fmesa, &
   f9xglwindowpos4dmesa, &
   f9xglwindowpos4ivmesa, &
   f9xglwindowpos4svmesa, &
   f9xglwindowpos4fvmesa, &
   f9xglwindowpos4dvmesa
#endif

!-------------------
! symbolic constants
!-------------------

! OpenGL 1.0

integer(kind=glenum), parameter, public :: GL_CURRENT_BIT = 1
integer(kind=glenum), parameter, public :: GL_POINT_BIT = 2
integer(kind=glenum), parameter, public :: GL_LINE_BIT = 4
integer(kind=glenum), parameter, public :: GL_POLYGON_BIT = 8
integer(kind=glenum), parameter, public :: GL_POLYGON_STIPPLE_BIT = 16
integer(kind=glenum), parameter, public :: GL_PIXEL_MODE_BIT = 32
integer(kind=glenum), parameter, public :: GL_LIGHTING_BIT = 64
integer(kind=glenum), parameter, public :: GL_FOG_BIT = 128
integer(kind=glenum), parameter, public :: GL_DEPTH_BUFFER_BIT = 256
integer(kind=glenum), parameter, public :: GL_ACCUM_BUFFER_BIT = 512
integer(kind=glenum), parameter, public :: GL_STENCIL_BUFFER_BIT = 1024
integer(kind=glenum), parameter, public :: GL_VIEWPORT_BIT = 2048
integer(kind=glenum), parameter, public :: GL_TRANSFORM_BIT = 4096
integer(kind=glenum), parameter, public :: GL_ENABLE_BIT = 8192
integer(kind=glenum), parameter, public :: GL_COLOR_BUFFER_BIT = 16384
integer(kind=glenum), parameter, public :: GL_HINT_BIT = 32768
integer(kind=glenum), parameter, public :: GL_EVAL_BIT = 65536
integer(kind=glenum), parameter, public :: GL_LIST_BIT = 131072
integer(kind=glenum), parameter, public :: GL_TEXTURE_BIT = 262144
integer(kind=glenum), parameter, public :: GL_SCISSOR_BIT = 524288
integer(kind=glenum), parameter, public :: GL_ALL_ATTRIB_BITS = 1048575
logical(kind=glboolean), parameter, public :: GL_FALSE = .false.
logical(kind=glboolean), parameter, public :: GL_TRUE = .true.
integer(kind=glenum), parameter, public :: GL_POINTS = 0
integer(kind=glenum), parameter, public :: GL_LINES = 1
integer(kind=glenum), parameter, public :: GL_LINE_LOOP = 2
integer(kind=glenum), parameter, public :: GL_LINE_STRIP = 3
integer(kind=glenum), parameter, public :: GL_TRIANGLES = 4
integer(kind=glenum), parameter, public :: GL_TRIANGLE_STRIP = 5
integer(kind=glenum), parameter, public :: GL_TRIANGLE_FAN = 6
integer(kind=glenum), parameter, public :: GL_QUADS = 7
integer(kind=glenum), parameter, public :: GL_QUAD_STRIP = 8
integer(kind=glenum), parameter, public :: GL_POLYGON = 9
integer(kind=glenum), parameter, public :: GL_ACCUM = 256
integer(kind=glenum), parameter, public :: GL_LOAD = 257
integer(kind=glenum), parameter, public :: GL_RETURN = 258
integer(kind=glenum), parameter, public :: GL_MULT = 259
integer(kind=glenum), parameter, public :: GL_ADD = 260
integer(kind=glenum), parameter, public :: GL_NEVER = 512
integer(kind=glenum), parameter, public :: GL_LESS = 513
integer(kind=glenum), parameter, public :: GL_EQUAL = 514
integer(kind=glenum), parameter, public :: GL_LEQUAL = 515
integer(kind=glenum), parameter, public :: GL_GREATER = 516
integer(kind=glenum), parameter, public :: GL_NOTEQUAL = 517
integer(kind=glenum), parameter, public :: GL_GEQUAL = 518
integer(kind=glenum), parameter, public :: GL_ALWAYS = 519
integer(kind=glenum), parameter, public :: GL_ZERO = 0
integer(kind=glenum), parameter, public :: GL_ONE = 1
integer(kind=glenum), parameter, public :: GL_SRC_COLOR = 768
integer(kind=glenum), parameter, public :: GL_ONE_MINUS_SRC_COLOR = 769
integer(kind=glenum), parameter, public :: GL_SRC_ALPHA = 770
integer(kind=glenum), parameter, public :: GL_ONE_MINUS_SRC_ALPHA = 771
integer(kind=glenum), parameter, public :: GL_DST_ALPHA = 772
integer(kind=glenum), parameter, public :: GL_ONE_MINUS_DST_ALPHA = 773
integer(kind=glenum), parameter, public :: GL_DST_COLOR = 774
integer(kind=glenum), parameter, public :: GL_ONE_MINUS_DST_COLOR = 775
integer(kind=glenum), parameter, public :: GL_SRC_ALPHA_SATURATE = 776
integer(kind=glenum), parameter, public :: GL_NONE = 0
integer(kind=glenum), parameter, public :: GL_FRONT_LEFT = 1024
integer(kind=glenum), parameter, public :: GL_FRONT_RIGHT = 1025
integer(kind=glenum), parameter, public :: GL_BACK_LEFT = 1026
integer(kind=glenum), parameter, public :: GL_BACK_RIGHT = 1027
integer(kind=glenum), parameter, public :: GL_FRONT = 1028
integer(kind=glenum), parameter, public :: GL_BACK = 1029
integer(kind=glenum), parameter, public :: GL_LEFT = 1030
integer(kind=glenum), parameter, public :: GL_RIGHT = 1031
integer(kind=glenum), parameter, public :: GL_FRONT_AND_BACK = 1032
integer(kind=glenum), parameter, public :: GL_AUX0 = 1033
integer(kind=glenum), parameter, public :: GL_AUX1 = 1034
integer(kind=glenum), parameter, public :: GL_AUX2 = 1035
integer(kind=glenum), parameter, public :: GL_AUX3 = 1036
integer(kind=glenum), parameter, public :: GL_NO_ERROR = 0
integer(kind=glenum), parameter, public :: GL_INVALID_ENUM = 1280
integer(kind=glenum), parameter, public :: GL_INVALID_VALUE = 1281
integer(kind=glenum), parameter, public :: GL_INVALID_OPERATION = 1282
integer(kind=glenum), parameter, public :: GL_STACK_OVERFLOW = 1283
integer(kind=glenum), parameter, public :: GL_STACK_UNDERFLOW = 1284
integer(kind=glenum), parameter, public :: GL_OUT_OF_MEMORY = 1285
integer(kind=glenum), parameter, public :: GL_2D = 1536
integer(kind=glenum), parameter, public :: GL_3D = 1537
integer(kind=glenum), parameter, public :: GL_3D_COLOR = 1538
integer(kind=glenum), parameter, public :: GL_3D_COLOR_TEXTURE = 1539
integer(kind=glenum), parameter, public :: GL_4D_COLOR_TEXTURE = 1540
integer(kind=glenum), parameter, public :: GL_PASS_THROUGH_TOKEN = 1792
integer(kind=glenum), parameter, public :: GL_POINT_TOKEN = 1793
integer(kind=glenum), parameter, public :: GL_LINE_TOKEN = 1794
integer(kind=glenum), parameter, public :: GL_POLYGON_TOKEN = 1795
integer(kind=glenum), parameter, public :: GL_BITMAP_TOKEN = 1796
integer(kind=glenum), parameter, public :: GL_DRAW_PIXEL_TOKEN = 1797
integer(kind=glenum), parameter, public :: GL_COPY_PIXEL_TOKEN = 1798
integer(kind=glenum), parameter, public :: GL_LINE_RESET_TOKEN = 1799
integer(kind=glenum), parameter, public :: GL_EXP = 2048
integer(kind=glenum), parameter, public :: GL_EXP2 = 2049
integer(kind=glenum), parameter, public :: GL_CW = 2304
integer(kind=glenum), parameter, public :: GL_CCW = 2305
integer(kind=glenum), parameter, public :: GL_COEFF = 2560
integer(kind=glenum), parameter, public :: GL_ORDER = 2561
integer(kind=glenum), parameter, public :: GL_DOMAIN = 2562
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_I_TO_I = 3184
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_S_TO_S = 3185
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_I_TO_R = 3186
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_I_TO_G = 3187
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_I_TO_B = 3188
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_I_TO_A = 3189
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_R_TO_R = 3190
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_G_TO_G = 3191
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_B_TO_B = 3192
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_A_TO_A = 3193
integer(kind=glenum), parameter, public :: GL_CURRENT_COLOR = 2816
integer(kind=glenum), parameter, public :: GL_CURRENT_INDEX = 2817
integer(kind=glenum), parameter, public :: GL_CURRENT_NORMAL = 2818
integer(kind=glenum), parameter, public :: GL_CURRENT_TEXTURE_COORDS = 2819
integer(kind=glenum), parameter, public :: GL_CURRENT_RASTER_COLOR = 2820
integer(kind=glenum), parameter, public :: GL_CURRENT_RASTER_INDEX = 2821
integer(kind=glenum), parameter, public :: GL_CURRENT_RASTER_TEXTURE_COORD =2822
integer(kind=glenum), parameter, public :: GL_CURRENT_RASTER_POSITION = 2823
integer(kind=glenum), parameter, public :: GL_CURRENT_RASTER_POSITION_VALI =2824
integer(kind=glenum), parameter, public :: GL_CURRENT_RASTER_DISTANCE = 2825
integer(kind=glenum), parameter, public :: GL_POINT_SMOOTH = 2832
integer(kind=glenum), parameter, public :: GL_POINT_SIZE = 2833
integer(kind=glenum), parameter, public :: GL_POINT_SIZE_RANGE = 2834
integer(kind=glenum), parameter, public :: GL_POINT_SIZE_GRANULARITY = 2835
integer(kind=glenum), parameter, public :: GL_LINE_SMOOTH = 2848
integer(kind=glenum), parameter, public :: GL_LINE_WIDTH = 2849
integer(kind=glenum), parameter, public :: GL_LINE_WIDTH_RANGE = 2850
integer(kind=glenum), parameter, public :: GL_LINE_WIDTH_GRANULARITY = 2851
integer(kind=glenum), parameter, public :: GL_LINE_STIPPLE = 2852
integer(kind=glenum), parameter, public :: GL_LINE_STIPPLE_PATTERN = 2853
integer(kind=glenum), parameter, public :: GL_LINE_STIPPLE_REPEAT = 2854
integer(kind=glenum), parameter, public :: GL_LIST_MODE = 2864
integer(kind=glenum), parameter, public :: GL_MAX_LIST_NESTING = 2865
integer(kind=glenum), parameter, public :: GL_LIST_BASE = 2866
integer(kind=glenum), parameter, public :: GL_LIST_INDEX = 2867
integer(kind=glenum), parameter, public :: GL_POLYGON_MODE = 2880
integer(kind=glenum), parameter, public :: GL_POLYGON_SMOOTH = 2881
integer(kind=glenum), parameter, public :: GL_POLYGON_STIPPLE = 2882
integer(kind=glenum), parameter, public :: GL_EDGE_FLAG = 2883
integer(kind=glenum), parameter, public :: GL_CULL_FACE = 2884
integer(kind=glenum), parameter, public :: GL_CULL_FACE_MODE = 2885
integer(kind=glenum), parameter, public :: GL_FRONT_FACE = 2886
integer(kind=glenum), parameter, public :: GL_LIGHTING = 2896
integer(kind=glenum), parameter, public :: GL_LIGHT_MODEL_LOCAL_VIEWER = 2897
integer(kind=glenum), parameter, public :: GL_LIGHT_MODEL_TWO_SIDE = 2898
integer(kind=glenum), parameter, public :: GL_LIGHT_MODEL_AMBIENT = 2899
integer(kind=glenum), parameter, public :: GL_SHADE_MODEL = 2900
integer(kind=glenum), parameter, public :: GL_COLOR_MATERIAL_FACE = 2901
integer(kind=glenum), parameter, public :: GL_COLOR_MATERIAL_PARAMETER = 2902
integer(kind=glenum), parameter, public :: GL_COLOR_MATERIAL = 2903
integer(kind=glenum), parameter, public :: GL_FOG = 2912
integer(kind=glenum), parameter, public :: GL_FOG_INDEX = 2913
integer(kind=glenum), parameter, public :: GL_FOG_DENSITY = 2914
integer(kind=glenum), parameter, public :: GL_FOG_START = 2915
integer(kind=glenum), parameter, public :: GL_FOG_END = 2916
integer(kind=glenum), parameter, public :: GL_FOG_MODE = 2917
integer(kind=glenum), parameter, public :: GL_FOG_COLOR = 2918
integer(kind=glenum), parameter, public :: GL_DEPTH_RANGE = 2928
integer(kind=glenum), parameter, public :: GL_DEPTH_TEST = 2929
integer(kind=glenum), parameter, public :: GL_DEPTH_WRITEMASK = 2930
integer(kind=glenum), parameter, public :: GL_DEPTH_CLEAR_VALUE = 2931
integer(kind=glenum), parameter, public :: GL_DEPTH_FUNC = 2932
integer(kind=glenum), parameter, public :: GL_ACCUM_CLEAR_VALUE = 2944
integer(kind=glenum), parameter, public :: GL_STENCIL_TEST = 2960
integer(kind=glenum), parameter, public :: GL_STENCIL_CLEAR_VALUE = 2961
integer(kind=glenum), parameter, public :: GL_STENCIL_FUNC = 2962
integer(kind=glenum), parameter, public :: GL_STENCIL_VALUE_MASK = 2963
integer(kind=glenum), parameter, public :: GL_STENCIL_FAIL = 2964
integer(kind=glenum), parameter, public :: GL_STENCIL_PASS_DEPTH_FAIL = 2965
integer(kind=glenum), parameter, public :: GL_STENCIL_PASS_DEPTH_PASS = 2966
integer(kind=glenum), parameter, public :: GL_STENCIL_REF = 2967
integer(kind=glenum), parameter, public :: GL_STENCIL_WRITEMASK = 2968
integer(kind=glenum), parameter, public :: GL_MATRIX_MODE = 2976
integer(kind=glenum), parameter, public :: GL_NORMALIZE = 2977
integer(kind=glenum), parameter, public :: GL_VIEWPORT = 2978
integer(kind=glenum), parameter, public :: GL_MODELVIEW_STACK_DEPTH = 2979
integer(kind=glenum), parameter, public :: GL_PROJECTION_STACK_DEPTH = 2980
integer(kind=glenum), parameter, public :: GL_TEXTURE_STACK_DEPTH = 2981
integer(kind=glenum), parameter, public :: GL_MODELVIEW_MATRIX = 2982
integer(kind=glenum), parameter, public :: GL_PROJECTION_MATRIX = 2983
integer(kind=glenum), parameter, public :: GL_TEXTURE_MATRIX = 2984
integer(kind=glenum), parameter, public :: GL_ATTRIB_STACK_DEPTH = 2992
integer(kind=glenum), parameter, public :: GL_ALPHA_TEST = 3008
integer(kind=glenum), parameter, public :: GL_ALPHA_TEST_FUNC = 3009
integer(kind=glenum), parameter, public :: GL_ALPHA_TEST_REF = 3010
integer(kind=glenum), parameter, public :: GL_DITHER = 3024
integer(kind=glenum), parameter, public :: GL_BLEND_DST = 3040
integer(kind=glenum), parameter, public :: GL_BLEND_SRC = 3041
integer(kind=glenum), parameter, public :: GL_BLEND = 3042
integer(kind=glenum), parameter, public :: GL_LOGIC_OP_MODE = 3056
integer(kind=glenum), parameter, public :: GL_LOGIC_OP = 3057
integer(kind=glenum), parameter, public :: GL_AUX_BUFFERS = 3072
integer(kind=glenum), parameter, public :: GL_DRAW_BUFFER = 3073
integer(kind=glenum), parameter, public :: GL_READ_BUFFER = 3074
integer(kind=glenum), parameter, public :: GL_SCISSOR_BOX = 3088
integer(kind=glenum), parameter, public :: GL_SCISSOR_TEST = 3089
integer(kind=glenum), parameter, public :: GL_INDEX_CLEAR_VALUE = 3104
integer(kind=glenum), parameter, public :: GL_INDEX_WRITEMASK = 3105
integer(kind=glenum), parameter, public :: GL_COLOR_CLEAR_VALUE = 3106
integer(kind=glenum), parameter, public :: GL_COLOR_WRITEMASK = 3107
integer(kind=glenum), parameter, public :: GL_INDEX_MODE = 3120
integer(kind=glenum), parameter, public :: GL_RGBA_MODE = 3121
integer(kind=glenum), parameter, public :: GL_DOUBLEBUFFER = 3122
integer(kind=glenum), parameter, public :: GL_STEREO = 3123
integer(kind=glenum), parameter, public :: GL_RENDER_MODE = 3136
integer(kind=glenum), parameter, public :: GL_PERSPECTIVE_CORRECTION_HINT = 3152
integer(kind=glenum), parameter, public :: GL_POINT_SMOOTH_HINT = 3153
integer(kind=glenum), parameter, public :: GL_LINE_SMOOTH_HINT = 3154
integer(kind=glenum), parameter, public :: GL_POLYGON_SMOOTH_HINT = 3155
integer(kind=glenum), parameter, public :: GL_FOG_HINT = 3156
integer(kind=glenum), parameter, public :: GL_TEXTURE_GEN_S = 3168
integer(kind=glenum), parameter, public :: GL_TEXTURE_GEN_T = 3169
integer(kind=glenum), parameter, public :: GL_TEXTURE_GEN_R = 3170
integer(kind=glenum), parameter, public :: GL_TEXTURE_GEN_Q = 3171
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_I_TO_I_SIZE = 3248
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_S_TO_S_SIZE = 3249
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_I_TO_R_SIZE = 3250
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_I_TO_G_SIZE = 3251
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_I_TO_B_SIZE = 3252
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_I_TO_A_SIZE = 3253
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_R_TO_R_SIZE = 3254
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_G_TO_G_SIZE = 3255
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_B_TO_B_SIZE = 3256
integer(kind=glenum), parameter, public :: GL_PIXEL_MAP_A_TO_A_SIZE = 3257
integer(kind=glenum), parameter, public :: GL_UNPACK_SWAP_BYTES = 3312
integer(kind=glenum), parameter, public :: GL_UNPACK_LSB_FIRST = 3313
integer(kind=glenum), parameter, public :: GL_UNPACK_ROW_LENGTH = 3314
integer(kind=glenum), parameter, public :: GL_UNPACK_SKIP_ROWS = 3315
integer(kind=glenum), parameter, public :: GL_UNPACK_SKIP_PIXELS = 3316
integer(kind=glenum), parameter, public :: GL_UNPACK_ALIGNMENT = 3317
integer(kind=glenum), parameter, public :: GL_PACK_SWAP_BYTES = 3328
integer(kind=glenum), parameter, public :: GL_PACK_LSB_FIRST = 3329
integer(kind=glenum), parameter, public :: GL_PACK_ROW_LENGTH = 3330
integer(kind=glenum), parameter, public :: GL_PACK_SKIP_ROWS = 3331
integer(kind=glenum), parameter, public :: GL_PACK_SKIP_PIXELS = 3332
integer(kind=glenum), parameter, public :: GL_PACK_ALIGNMENT = 3333
integer(kind=glenum), parameter, public :: GL_MAP_COLOR = 3344
integer(kind=glenum), parameter, public :: GL_MAP_STENCIL = 3345
integer(kind=glenum), parameter, public :: GL_INDEX_SHIFT = 3346
integer(kind=glenum), parameter, public :: GL_INDEX_OFFSET = 3347
integer(kind=glenum), parameter, public :: GL_RED_SCALE = 3348
integer(kind=glenum), parameter, public :: GL_RED_BIAS = 3349
integer(kind=glenum), parameter, public :: GL_ZOOM_X = 3350
integer(kind=glenum), parameter, public :: GL_ZOOM_Y = 3351
integer(kind=glenum), parameter, public :: GL_GREEN_SCALE = 3352
integer(kind=glenum), parameter, public :: GL_GREEN_BIAS = 3353
integer(kind=glenum), parameter, public :: GL_BLUE_SCALE = 3354
integer(kind=glenum), parameter, public :: GL_BLUE_BIAS = 3355
integer(kind=glenum), parameter, public :: GL_ALPHA_SCALE = 3356
integer(kind=glenum), parameter, public :: GL_ALPHA_BIAS = 3357
integer(kind=glenum), parameter, public :: GL_DEPTH_SCALE = 3358
integer(kind=glenum), parameter, public :: GL_DEPTH_BIAS = 3359
integer(kind=glenum), parameter, public :: GL_MAX_EVAL_ORDER = 3376
integer(kind=glenum), parameter, public :: GL_MAX_LIGHTS = 3377
integer(kind=glenum), parameter, public :: GL_MAX_CLIP_PLANES = 3378
integer(kind=glenum), parameter, public :: GL_MAX_TEXTURE_SIZE = 3379
integer(kind=glenum), parameter, public :: GL_MAX_PIXEL_MAP_TABLE = 3380
integer(kind=glenum), parameter, public :: GL_MAX_ATTRIB_STACK_DEPTH = 3381
integer(kind=glenum), parameter, public :: GL_MAX_MODELVIEW_STACK_DEPTH = 3382
integer(kind=glenum), parameter, public :: GL_MAX_NAME_STACK_DEPTH = 3383
integer(kind=glenum), parameter, public :: GL_MAX_PROJECTION_STACK_DEPTH = 3384
integer(kind=glenum), parameter, public :: GL_MAX_TEXTURE_STACK_DEPTH = 3385
integer(kind=glenum), parameter, public :: GL_MAX_VIEWPORT_DIMS = 3386
integer(kind=glenum), parameter, public :: GL_SUBPIXEL_BITS = 3408
integer(kind=glenum), parameter, public :: GL_INDEX_BITS = 3409
integer(kind=glenum), parameter, public :: GL_RED_BITS = 3410
integer(kind=glenum), parameter, public :: GL_GREEN_BITS = 3411
integer(kind=glenum), parameter, public :: GL_BLUE_BITS = 3412
integer(kind=glenum), parameter, public :: GL_ALPHA_BITS = 3413
integer(kind=glenum), parameter, public :: GL_DEPTH_BITS = 3414
integer(kind=glenum), parameter, public :: GL_STENCIL_BITS = 3415
integer(kind=glenum), parameter, public :: GL_ACCUM_RED_BITS = 3416
integer(kind=glenum), parameter, public :: GL_ACCUM_GREEN_BITS = 3417
integer(kind=glenum), parameter, public :: GL_ACCUM_BLUE_BITS = 3418
integer(kind=glenum), parameter, public :: GL_ACCUM_ALPHA_BITS = 3419
integer(kind=glenum), parameter, public :: GL_NAME_STACK_DEPTH = 3440
integer(kind=glenum), parameter, public :: GL_AUTO_NORMAL = 3456
integer(kind=glenum), parameter, public :: GL_MAP1_COLOR_4 = 3472
integer(kind=glenum), parameter, public :: GL_MAP1_INDEX = 3473
integer(kind=glenum), parameter, public :: GL_MAP1_NORMAL = 3474
integer(kind=glenum), parameter, public :: GL_MAP1_TEXTURE_COORD_1 = 3475
integer(kind=glenum), parameter, public :: GL_MAP1_TEXTURE_COORD_2 = 3476
integer(kind=glenum), parameter, public :: GL_MAP1_TEXTURE_COORD_3 = 3477
integer(kind=glenum), parameter, public :: GL_MAP1_TEXTURE_COORD_4 = 3478
integer(kind=glenum), parameter, public :: GL_MAP1_VERTEX_3 = 3479
integer(kind=glenum), parameter, public :: GL_MAP1_VERTEX_4 = 3480
integer(kind=glenum), parameter, public :: GL_MAP2_COLOR_4 = 3504
integer(kind=glenum), parameter, public :: GL_MAP2_INDEX = 3505
integer(kind=glenum), parameter, public :: GL_MAP2_NORMAL = 3506
integer(kind=glenum), parameter, public :: GL_MAP2_TEXTURE_COORD_1 = 3507
integer(kind=glenum), parameter, public :: GL_MAP2_TEXTURE_COORD_2 = 3508
integer(kind=glenum), parameter, public :: GL_MAP2_TEXTURE_COORD_3 = 3509
integer(kind=glenum), parameter, public :: GL_MAP2_TEXTURE_COORD_4 = 3510
integer(kind=glenum), parameter, public :: GL_MAP2_VERTEX_3 = 3511
integer(kind=glenum), parameter, public :: GL_MAP2_VERTEX_4 = 3512
integer(kind=glenum), parameter, public :: GL_MAP1_GRID_DOMAIN = 3536
integer(kind=glenum), parameter, public :: GL_MAP1_GRID_SEGMENTS = 3537
integer(kind=glenum), parameter, public :: GL_MAP2_GRID_DOMAIN = 3538
integer(kind=glenum), parameter, public :: GL_MAP2_GRID_SEGMENTS = 3539
integer(kind=glenum), parameter, public :: GL_TEXTURE_1D = 3552
integer(kind=glenum), parameter, public :: GL_TEXTURE_2D = 3553
integer(kind=glenum), parameter, public :: GL_TEXTURE_WIDTH = 4096
integer(kind=glenum), parameter, public :: GL_TEXTURE_HEIGHT = 4097
integer(kind=glenum), parameter, public :: GL_TEXTURE_COMPONENTS = 4099
integer(kind=glenum), parameter, public :: GL_TEXTURE_BORDER_COLOR = 4100
integer(kind=glenum), parameter, public :: GL_TEXTURE_BORDER = 4101
integer(kind=glenum), parameter, public :: GL_DONT_CARE = 4352
integer(kind=glenum), parameter, public :: GL_FASTEST = 4353
integer(kind=glenum), parameter, public :: GL_NICEST = 4354
integer(kind=glenum), parameter, public :: GL_AMBIENT = 4608
integer(kind=glenum), parameter, public :: GL_DIFFUSE = 4609
integer(kind=glenum), parameter, public :: GL_SPECULAR = 4610
integer(kind=glenum), parameter, public :: GL_POSITION = 4611
integer(kind=glenum), parameter, public :: GL_SPOT_DIRECTION = 4612
integer(kind=glenum), parameter, public :: GL_SPOT_EXPONENT = 4613
integer(kind=glenum), parameter, public :: GL_SPOT_CUTOFF = 4614
integer(kind=glenum), parameter, public :: GL_CONSTANT_ATTENUATION = 4615
integer(kind=glenum), parameter, public :: GL_LINEAR_ATTENUATION = 4616
integer(kind=glenum), parameter, public :: GL_QUADRATIC_ATTENUATION = 4617
integer(kind=glenum), parameter, public :: GL_COMPILE = 4864
integer(kind=glenum), parameter, public :: GL_COMPILE_AND_EXECUTE = 4865
integer(kind=glenum), parameter, public :: GL_BYTE = 5120
integer(kind=glenum), parameter, public :: GL_UNSIGNED_BYTE = 5121
integer(kind=glenum), parameter, public :: GL_SHORT = 5122
integer(kind=glenum), parameter, public :: GL_UNSIGNED_SHORT = 5123
integer(kind=glenum), parameter, public :: GL_INT = 5124
integer(kind=glenum), parameter, public :: GL_UNSIGNED_INT = 5125
integer(kind=glenum), parameter, public :: GL_FLOAT = 5126
integer(kind=glenum), parameter, public :: GL_2_BYTES = 5127
integer(kind=glenum), parameter, public :: GL_3_BYTES = 5128
integer(kind=glenum), parameter, public :: GL_4_BYTES = 5129
integer(kind=glenum), parameter, public :: GL_DOUBLE_EXT = 5130
integer(kind=glenum), parameter, public :: GL_CLEAR = 5376
integer(kind=glenum), parameter, public :: GL_AND = 5377
integer(kind=glenum), parameter, public :: GL_AND_REVERSE = 5378
integer(kind=glenum), parameter, public :: GL_COPY = 5379
integer(kind=glenum), parameter, public :: GL_AND_INVERTED = 5380
integer(kind=glenum), parameter, public :: GL_NOOP = 5381
integer(kind=glenum), parameter, public :: GL_XOR = 5382
integer(kind=glenum), parameter, public :: GL_OR = 5383
integer(kind=glenum), parameter, public :: GL_NOR = 5384
integer(kind=glenum), parameter, public :: GL_EQUIV = 5385
integer(kind=glenum), parameter, public :: GL_INVERT = 5386
integer(kind=glenum), parameter, public :: GL_OR_REVERSE = 5387
integer(kind=glenum), parameter, public :: GL_COPY_INVERTED = 5388
integer(kind=glenum), parameter, public :: GL_OR_INVERTED = 5389
integer(kind=glenum), parameter, public :: GL_NAND = 5390
integer(kind=glenum), parameter, public :: GL_SET = 5391
integer(kind=glenum), parameter, public :: GL_EMISSION = 5632
integer(kind=glenum), parameter, public :: GL_SHININESS = 5633
integer(kind=glenum), parameter, public :: GL_AMBIENT_AND_DIFFUSE = 5634
integer(kind=glenum), parameter, public :: GL_COLOR_INDEXES = 5635
integer(kind=glenum), parameter, public :: GL_MODELVIEW = 5888
integer(kind=glenum), parameter, public :: GL_PROJECTION = 5889
integer(kind=glenum), parameter, public :: GL_TEXTURE = 5890
integer(kind=glenum), parameter, public :: GL_COLOR = 6144
integer(kind=glenum), parameter, public :: GL_DEPTH = 6145
integer(kind=glenum), parameter, public :: GL_STENCIL = 6146
integer(kind=glenum), parameter, public :: GL_COLOR_INDEX = 6400
integer(kind=glenum), parameter, public :: GL_STENCIL_INDEX = 6401
integer(kind=glenum), parameter, public :: GL_DEPTH_COMPONENT = 6402
integer(kind=glenum), parameter, public :: GL_RED = 6403
integer(kind=glenum), parameter, public :: GL_GREEN = 6404
integer(kind=glenum), parameter, public :: GL_BLUE = 6405
integer(kind=glenum), parameter, public :: GL_ALPHA = 6406
integer(kind=glenum), parameter, public :: GL_RGB = 6407
integer(kind=glenum), parameter, public :: GL_RGBA = 6408
integer(kind=glenum), parameter, public :: GL_LUMINANCE = 6409
integer(kind=glenum), parameter, public :: GL_LUMINANCE_ALPHA = 6410
integer(kind=glenum), parameter, public :: GL_BITMAP = 6656
integer(kind=glenum), parameter, public :: GL_POINT = 6912
integer(kind=glenum), parameter, public :: GL_LINE = 6913
integer(kind=glenum), parameter, public :: GL_FILL = 6914
integer(kind=glenum), parameter, public :: GL_RENDER = 7168
integer(kind=glenum), parameter, public :: GL_FEEDBACK = 7169
integer(kind=glenum), parameter, public :: GL_SELECT = 7170
integer(kind=glenum), parameter, public :: GL_FLAT = 7424
integer(kind=glenum), parameter, public :: GL_SMOOTH = 7425
integer(kind=glenum), parameter, public :: GL_KEEP = 7680
integer(kind=glenum), parameter, public :: GL_REPLACE = 7681
integer(kind=glenum), parameter, public :: GL_INCR = 7682
integer(kind=glenum), parameter, public :: GL_DECR = 7683
integer(kind=glenum), parameter, public :: GL_VENDOR = 7936
integer(kind=glenum), parameter, public :: GL_RENDERER = 7937
integer(kind=glenum), parameter, public :: GL_VERSION = 7938
integer(kind=glenum), parameter, public :: GL_EXTENSIONS = 7939
integer(kind=glenum), parameter, public :: GL_S = 8192
integer(kind=glenum), parameter, public :: GL_T = 8193
integer(kind=glenum), parameter, public :: GL_R = 8194
integer(kind=glenum), parameter, public :: GL_Q = 8195
integer(kind=glenum), parameter, public :: GL_MODULATE = 8448
integer(kind=glenum), parameter, public :: GL_DECAL = 8449
integer(kind=glenum), parameter, public :: GL_TEXTURE_ENV_MODE = 8704
integer(kind=glenum), parameter, public :: GL_TEXTURE_ENV_COLOR = 8705
integer(kind=glenum), parameter, public :: GL_TEXTURE_ENV = 8960
integer(kind=glenum), parameter, public :: GL_EYE_LINEAR = 9216
integer(kind=glenum), parameter, public :: GL_OBJECT_LINEAR = 9217
integer(kind=glenum), parameter, public :: GL_SPHERE_MAP = 9218
integer(kind=glenum), parameter, public :: GL_TEXTURE_GEN_MODE = 9472
integer(kind=glenum), parameter, public :: GL_OBJECT_PLANE = 9473
integer(kind=glenum), parameter, public :: GL_EYE_PLANE = 9474
integer(kind=glenum), parameter, public :: GL_NEAREST = 9728
integer(kind=glenum), parameter, public :: GL_LINEAR = 9729
integer(kind=glenum), parameter, public :: GL_NEAREST_MIPMAP_NEAREST = 9984
integer(kind=glenum), parameter, public :: GL_LINEAR_MIPMAP_NEAREST = 9985
integer(kind=glenum), parameter, public :: GL_NEAREST_MIPMAP_LINEAR = 9986
integer(kind=glenum), parameter, public :: GL_LINEAR_MIPMAP_LINEAR = 9987
integer(kind=glenum), parameter, public :: GL_TEXTURE_MAG_FILTER = 10240
integer(kind=glenum), parameter, public :: GL_TEXTURE_MIN_FILTER = 10241
integer(kind=glenum), parameter, public :: GL_TEXTURE_WRAP_S = 10242
integer(kind=glenum), parameter, public :: GL_TEXTURE_WRAP_T = 10243
integer(kind=glenum), parameter, public :: GL_CLAMP = 10496
integer(kind=glenum), parameter, public :: GL_REPEAT = 10497
integer(kind=glenum), parameter, public :: GL_CLIP_PLANE0 = 12288
integer(kind=glenum), parameter, public :: GL_CLIP_PLANE1 = 12289
integer(kind=glenum), parameter, public :: GL_CLIP_PLANE2 = 12290
integer(kind=glenum), parameter, public :: GL_CLIP_PLANE3 = 12291
integer(kind=glenum), parameter, public :: GL_CLIP_PLANE4 = 12292
integer(kind=glenum), parameter, public :: GL_CLIP_PLANE5 = 12293
integer(kind=glenum), parameter, public :: GL_LIGHT0 = 16384
integer(kind=glenum), parameter, public :: GL_LIGHT1 = 16385
integer(kind=glenum), parameter, public :: GL_LIGHT2 = 16386
integer(kind=glenum), parameter, public :: GL_LIGHT3 = 16387
integer(kind=glenum), parameter, public :: GL_LIGHT4 = 16388
integer(kind=glenum), parameter, public :: GL_LIGHT5 = 16389
integer(kind=glenum), parameter, public :: GL_LIGHT6 = 16390
integer(kind=glenum), parameter, public :: GL_LIGHT7 = 16391

! OpenGL 1.1

#ifdef OPENGL_1_1
integer(kind=glenum), parameter, public :: GL_CLIENT_PIXEL_STORE_BIT = 1
integer(kind=glenum), parameter, public :: GL_CLIENT_VERTEX_ARRAY_BIT = 2
integer(kind=glenum), parameter, public :: GL_CLIENT_ALL_ATTRIB_BITS = 65535
integer(kind=glenum), parameter, public :: GL_CONSTANT_COLOR = 32769
integer(kind=glenum), parameter, public :: GL_ONE_MINUS_CONSTANT_COLOR = 32770
integer(kind=glenum), parameter, public :: GL_CONSTANT_ALPHA = 32771
integer(kind=glenum), parameter, public :: GL_ONE_MINUS_CONSTANT_ALPHA = 32772
integer(kind=glenum), parameter, public :: GL_VERTEX_ARRAY_POINTER = 32910
integer(kind=glenum), parameter, public :: GL_NORMAL_ARRAY_POINTER = 32911
integer(kind=glenum), parameter, public :: GL_COLOR_ARRAY_POINTER = 32912
integer(kind=glenum), parameter, public :: GL_INDEX_ARRAY_POINTER = 32913
integer(kind=glenum), parameter, public :: GL_TEXTURE_COORD_ARRAY_POINTER =32914
integer(kind=glenum), parameter, public :: GL_EDGE_FLAG_ARRAY_POINTER = 32915
integer(kind=glenum), parameter, public :: GL_CLIENT_ATTRIB_STACK_DEPTH = 2993
integer(kind=glenum), parameter, public :: GL_INDEX_LOGIC_OP = 3057
integer(kind=glenum), parameter, public :: GL_COLOR_LOGIC_OP = 3058
integer(kind=glenum), parameter, public :: GL_MAX_CLIENT_ATTRIB_STACK_DEPT =3387
integer(kind=glenum), parameter, public :: GL_FEEDBACK_BUFFER_POINTER = 3568
integer(kind=glenum), parameter, public :: GL_FEEDBACK_BUFFER_SIZE = 3569
integer(kind=glenum), parameter, public :: GL_FEEDBACK_BUFFER_TYPE = 3570
integer(kind=glenum), parameter, public :: GL_SELECTION_BUFFER_POINTER = 3571
integer(kind=glenum), parameter, public :: GL_SELECTION_BUFFER_SIZE = 3572
integer(kind=glenum), parameter, public :: GL_POLYGON_OFFSET_UNITS = 10752
integer(kind=glenum), parameter, public :: GL_POLYGON_OFFSET_POINT = 10753
integer(kind=glenum), parameter, public :: GL_POLYGON_OFFSET_LINE = 10754
integer(kind=glenum), parameter, public :: GL_POLYGON_OFFSET_FILL = 32823
integer(kind=glenum), parameter, public :: GL_POLYGON_OFFSET_FACTOR = 32824
integer(kind=glenum), parameter, public :: GL_TEXTURE_BINDING_1D = 32872
integer(kind=glenum), parameter, public :: GL_TEXTURE_BINDING_2D = 32873
integer(kind=glenum), parameter, public :: GL_VERTEX_ARRAY = 32884
integer(kind=glenum), parameter, public :: GL_NORMAL_ARRAY = 32885
integer(kind=glenum), parameter, public :: GL_COLOR_ARRAY = 32886
integer(kind=glenum), parameter, public :: GL_INDEX_ARRAY = 32887
integer(kind=glenum), parameter, public :: GL_TEXTURE_COORD_ARRAY = 32888
integer(kind=glenum), parameter, public :: GL_EDGE_FLAG_ARRAY = 32889
integer(kind=glenum), parameter, public :: GL_VERTEX_ARRAY_SIZE = 32890
integer(kind=glenum), parameter, public :: GL_VERTEX_ARRAY_TYPE = 32891
integer(kind=glenum), parameter, public :: GL_VERTEX_ARRAY_STRIDE = 32892
integer(kind=glenum), parameter, public :: GL_NORMAL_ARRAY_TYPE = 32894
integer(kind=glenum), parameter, public :: GL_NORMAL_ARRAY_STRIDE = 32895
integer(kind=glenum), parameter, public :: GL_COLOR_ARRAY_SIZE = 32897
integer(kind=glenum), parameter, public :: GL_COLOR_ARRAY_TYPE = 32898
integer(kind=glenum), parameter, public :: GL_COLOR_ARRAY_STRIDE = 32899
integer(kind=glenum), parameter, public :: GL_INDEX_ARRAY_TYPE = 32901
integer(kind=glenum), parameter, public :: GL_INDEX_ARRAY_STRIDE = 32902
integer(kind=glenum), parameter, public :: GL_TEXTURE_COORD_ARRAY_SIZE = 32904
integer(kind=glenum), parameter, public :: GL_TEXTURE_COORD_ARRAY_TYPE = 32905
integer(kind=glenum), parameter, public :: GL_TEXTURE_COORD_ARRAY_STRIDE = 32906
integer(kind=glenum), parameter, public :: GL_EDGE_FLAG_ARRAY_STRIDE = 32908
integer(kind=glenum), parameter, public :: GL_TEXTURE_INTERNAL_FORMAT = 4099
integer(kind=glenum), parameter, public :: GL_TEXTURE_RED_SIZE = 32860
integer(kind=glenum), parameter, public :: GL_TEXTURE_GREEN_SIZE = 32861
integer(kind=glenum), parameter, public :: GL_TEXTURE_BLUE_SIZE = 32862
integer(kind=glenum), parameter, public :: GL_TEXTURE_ALPHA_SIZE = 32863
integer(kind=glenum), parameter, public :: GL_TEXTURE_LUMINANCE_SIZE = 32864
integer(kind=glenum), parameter, public :: GL_TEXTURE_INTENSITY_SIZE = 32865
integer(kind=glenum), parameter, public :: GL_TEXTURE_PRIORITY = 32870
integer(kind=glenum), parameter, public :: GL_TEXTURE_RESIDENT = 32871
integer(kind=glenum), parameter, public :: GL_DOUBLE = 5130
integer(kind=glenum), parameter, public :: GL_PROXY_TEXTURE_1D = 32867
integer(kind=glenum), parameter, public :: GL_PROXY_TEXTURE_2D = 32868
integer(kind=glenum), parameter, public :: GL_R3_G3_B2 = 10768
integer(kind=glenum), parameter, public :: GL_ALPHA4 = 32827
integer(kind=glenum), parameter, public :: GL_ALPHA8 = 32828
integer(kind=glenum), parameter, public :: GL_ALPHA12 = 32829
integer(kind=glenum), parameter, public :: GL_ALPHA16 = 32830
integer(kind=glenum), parameter, public :: GL_LUMINANCE4 = 32831
integer(kind=glenum), parameter, public :: GL_LUMINANCE8 = 32832
integer(kind=glenum), parameter, public :: GL_LUMINANCE12 = 32833
integer(kind=glenum), parameter, public :: GL_LUMINANCE16 = 32834
integer(kind=glenum), parameter, public :: GL_LUMINANCE4_ALPHA4 = 32835
integer(kind=glenum), parameter, public :: GL_LUMINANCE6_ALPHA2 = 32836
integer(kind=glenum), parameter, public :: GL_LUMINANCE8_ALPHA8 = 32837
integer(kind=glenum), parameter, public :: GL_LUMINANCE12_ALPHA4 = 32838
integer(kind=glenum), parameter, public :: GL_LUMINANCE12_ALPHA12 = 32839
integer(kind=glenum), parameter, public :: GL_LUMINANCE16_ALPHA16 = 32840
integer(kind=glenum), parameter, public :: GL_INTENSITY = 32841
integer(kind=glenum), parameter, public :: GL_INTENSITY4 = 32842
integer(kind=glenum), parameter, public :: GL_INTENSITY8 = 32843
integer(kind=glenum), parameter, public :: GL_INTENSITY12 = 32844
integer(kind=glenum), parameter, public :: GL_INTENSITY16 = 32845
integer(kind=glenum), parameter, public :: GL_RGB4 = 32847
integer(kind=glenum), parameter, public :: GL_RGB5 = 32848
integer(kind=glenum), parameter, public :: GL_RGB8 = 32849
integer(kind=glenum), parameter, public :: GL_RGB10 = 32850
integer(kind=glenum), parameter, public :: GL_RGB12 = 32851
integer(kind=glenum), parameter, public :: GL_RGB16 = 32852
integer(kind=glenum), parameter, public :: GL_RGBA2 = 32853
integer(kind=glenum), parameter, public :: GL_RGBA4 = 32854
integer(kind=glenum), parameter, public :: GL_RGB5_A1 = 32855
integer(kind=glenum), parameter, public :: GL_RGBA8 = 32856
integer(kind=glenum), parameter, public :: GL_RGB10_A2 = 32857
integer(kind=glenum), parameter, public :: GL_RGBA12 = 32858
integer(kind=glenum), parameter, public :: GL_RGBA16 = 32859
integer(kind=glenum), parameter, public :: GL_V2F = 10784
integer(kind=glenum), parameter, public :: GL_V3F = 10785
integer(kind=glenum), parameter, public :: GL_C4UB_V2F = 10786
integer(kind=glenum), parameter, public :: GL_C4UB_V3F = 10787
integer(kind=glenum), parameter, public :: GL_C3F_V3F = 10788
integer(kind=glenum), parameter, public :: GL_N3F_V3F = 10789
integer(kind=glenum), parameter, public :: GL_C4F_N3F_V3F = 10790
integer(kind=glenum), parameter, public :: GL_T2F_V3F = 10791
integer(kind=glenum), parameter, public :: GL_T4F_V4F = 10792
integer(kind=glenum), parameter, public :: GL_T2F_C4UB_V3F = 10793
integer(kind=glenum), parameter, public :: GL_T2F_C3F_V3F = 10794
integer(kind=glenum), parameter, public :: GL_T2F_N3F_V3F = 10795
integer(kind=glenum), parameter, public :: GL_T2F_C4F_N3F_V3F = 10796
integer(kind=glenum), parameter, public :: GL_T4F_C4F_N3F_V4F = 10797
#endif

! OpenGL 1.2

#ifdef OPENGL_1_2
integer(kind=glenum), parameter, public :: GL_PACK_SKIP_IMAGES = 32875
integer(kind=glenum), parameter, public :: GL_PACK_IMAGE_HEIGHT = 32876
integer(kind=glenum), parameter, public :: GL_UNPACK_SKIP_IMAGES = 32877
integer(kind=glenum), parameter, public :: GL_UNPACK_IMAGE_HEIGHT = 32878
integer(kind=glenum), parameter, public :: GL_TEXTURE_3D = 32879
integer(kind=glenum), parameter, public :: GL_PROXY_TEXTURE_3D = 32880
integer(kind=glenum), parameter, public :: GL_TEXTURE_DEPTH = 32881
integer(kind=glenum), parameter, public :: GL_TEXTURE_WRAP_R = 32882
integer(kind=glenum), parameter, public :: GL_MAX_3D_TEXTURE_SIZE = 32883
integer(kind=glenum), parameter, public :: GL_TEXTURE_BINDING_3D = 32874
integer(kind=glenum), parameter, public :: GL_RESCALE_NORMAL = 32826
integer(kind=glenum), parameter, public :: GL_CLAMP_TO_EDGE = 33071
integer(kind=glenum), parameter, public :: GL_MAX_ELEMENTS_VERTICES = 61672
integer(kind=glenum), parameter, public :: GL_MAX_ELEMENTS_INDICES = 61673
integer(kind=glenum), parameter, public :: GL_BGR = 32992
integer(kind=glenum), parameter, public :: GL_BGRA = 32993
integer(kind=glenum), parameter, public :: GL_UNSIGNED_BYTE_3_3_2 = 32818
integer(kind=glenum), parameter, public :: GL_UNSIGNED_BYTE_2_3_3_REV = 33634
integer(kind=glenum), parameter, public :: GL_UNSIGNED_SHORT_5_6_5 = 33635
integer(kind=glenum), parameter, public :: GL_UNSIGNED_SHORT_5_6_5_REV = 33636
integer(kind=glenum), parameter, public :: GL_UNSIGNED_SHORT_4_4_4_4 = 32819
integer(kind=glenum), parameter, public :: GL_UNSIGNED_SHORT_4_4_4_4_REV = 33637
integer(kind=glenum), parameter, public :: GL_UNSIGNED_SHORT_5_5_5_1 = 32820
integer(kind=glenum), parameter, public :: GL_UNSIGNED_SHORT_1_5_5_5_REV = 33638
integer(kind=glenum), parameter, public :: GL_UNSIGNED_INT_8_8_8_8 = 32821
integer(kind=glenum), parameter, public :: GL_UNSIGNED_INT_8_8_8_8_REV = 33639
integer(kind=glenum), parameter, public :: GL_UNSIGNED_INT_10_10_10_2 = 32822
integer(kind=glenum), parameter, public :: GL_UNSIGNED_INT_2_10_10_10_REV = 33640
integer(kind=glenum), parameter, public :: GL_LIGHT_MODEL_COLOR_CONTROL = 33272
integer(kind=glenum), parameter, public :: GL_SINGLE_COLOR = 33273
integer(kind=glenum), parameter, public :: GL_SEPARATE_SPECULAR_COLOR = 33274
integer(kind=glenum), parameter, public :: GL_TEXTURE_MIN_LOD = 33082
integer(kind=glenum), parameter, public :: GL_TEXTURE_MAX_LOD = 33083
integer(kind=glenum), parameter, public :: GL_TEXTURE_BASE_LEVEL = 33084
integer(kind=glenum), parameter, public :: GL_TEXTURE_MAX_LEVEL = 33085
#endif

! EXT_abgr

#ifdef GL_EXT_abgr
integer(kind=glenum), parameter, public :: GL_ABGR_EXT = 32768
#endif

! EXT_blend_color

#ifdef GL_EXT_blend_color
integer(kind=glenum), parameter, public :: GL_CONSTANT_COLOR_EXT = 32769
integer(kind=glenum), parameter, public :: GL_ONE_MINUS_CONSTANT_COLOR_EXT=32770
integer(kind=glenum), parameter, public :: GL_CONSTANT_ALPHA_EXT = 32771
integer(kind=glenum), parameter, public :: GL_ONE_MINUS_CONSTANT_ALPHA_EXT=32772
integer(kind=glenum), parameter, public :: GL_BLEND_COLOR_EXT = 32773
#endif

! EXT_blend_minmax

#ifdef GL_EXT_blend_minmax
integer(kind=glenum), parameter, public :: GL_BLEND_EQUATION_EXT = 32777
integer(kind=glenum), parameter, public :: GL_MIN_EXT = 32775
integer(kind=glenum), parameter, public :: GL_MAX_EXT = 32776
integer(kind=glenum), parameter, public :: GL_FUNC_ADD_EXT = 32774
#endif

! EXT_blend_subtract

#ifdef GL_EXT_blend_subtract
integer(kind=glenum), parameter, public :: GL_FUNC_SUBTRACT_EXT = 32778
integer(kind=glenum), parameter, public :: GL_FUNC_REVERSE_SUBTRACT_EXT = 32779
#endif

! EXT_paletted_texture

#ifdef GL_EXT_paletted_texture
integer(kind=glenum), parameter, public :: GL_TABLE_TOO_LARGE_EXT = 32817
integer(kind=glenum), parameter, public :: GL_COLOR_TABLE_FORMAT_EXT = 32984
integer(kind=glenum), parameter, public :: GL_COLOR_TABLE_WIDTH_EXT = 32985
integer(kind=glenum), parameter, public :: GL_COLOR_TABLE_RED_SIZE_EXT = 32986
integer(kind=glenum), parameter, public :: GL_COLOR_TABLE_GREEN_SIZE_EXT = 32987
integer(kind=glenum), parameter, public :: GL_COLOR_TABLE_BLUE_SIZE_EXT = 32988
integer(kind=glenum), parameter, public :: GL_COLOR_TABLE_ALPHA_SIZE_EXT = 32989
integer(kind=glenum), parameter, public :: GL_COLOR_TABLE_LUMINANCE_SIZE_E=32990
integer(kind=glenum), parameter, public :: GL_COLOR_TABLE_INTENSITY_SIZE_E=32991
integer(kind=glenum), parameter, public :: GL_TEXTURE_INDEX_SIZE_EXT = 32993
integer(kind=glenum), parameter, public :: GL_COLOR_INDEX1_EXT = 32994
integer(kind=glenum), parameter, public :: GL_COLOR_INDEX2_EXT = 32995
integer(kind=glenum), parameter, public :: GL_COLOR_INDEX4_EXT = 32996
integer(kind=glenum), parameter, public :: GL_COLOR_INDEX8_EXT = 32997
integer(kind=glenum), parameter, public :: GL_COLOR_INDEX12_EXT = 32998
integer(kind=glenum), parameter, public :: GL_COLOR_INDEX16_EXT = 32999
#endif

! EXT_point_parameters

#ifdef GL_EXT_point_parameters
integer(kind=glenum), parameter, public :: GL_POINT_SIZE_MIN_EXT = 33062
integer(kind=glenum), parameter, public :: GL_POINT_SIZE_MAX_EXT = 33063
integer(kind=glenum), parameter, public :: GL_POINT_FADE_THRESHOLD_SIZE_EX=33064
integer(kind=glenum), parameter, public :: GL_DISTANCE_ATTENUATION_EXT = 33065
#endif

! EXT_polygon_offset

#ifdef GL_EXT_polygon_offset
integer(kind=glenum), parameter, public :: GL_POLYGON_OFFSET_EXT = 32823
integer(kind=glenum), parameter, public :: GL_POLYGON_OFFSET_FACTOR_EXT = 32824
integer(kind=glenum), parameter, public :: GL_POLYGON_OFFSET_BIAS_EXT = 32825
#endif

! EXT_rescale_normal

#ifdef GL_EXT_rescale_normal
integer(kind=glenum), parameter, public :: GL_RESCALE_NORMAL_EXT = 32826
#endif

! EXT_shared_texture_palette

#ifdef GL_EXT_shared_texture_palette
integer(kind=glenum), parameter, public :: GL_SHARED_TEXTURE_PALETTE_EXT = 33275
#endif

! EXT_texture3D

#ifdef GL_EXT_texture3D
integer(kind=glenum), parameter, public :: GL_PACK_SKIP_IMAGES_EXT = 32875
integer(kind=glenum), parameter, public :: GL_PACK_IMAGE_HEIGHT_EXT = 32876
integer(kind=glenum), parameter, public :: GL_UNPACK_SKIP_IMAGES_EXT = 32877
integer(kind=glenum), parameter, public :: GL_UNPACK_IMAGE_HEIGHT_EXT = 32878
integer(kind=glenum), parameter, public :: GL_TEXTURE_3D_EXT = 32879
integer(kind=glenum), parameter, public :: GL_PROXY_TEXTURE_3D_EXT = 32880
integer(kind=glenum), parameter, public :: GL_TEXTURE_DEPTH_EXT = 32881
integer(kind=glenum), parameter, public :: GL_TEXTURE_WRAP_R_EXT = 32882
integer(kind=glenum), parameter, public :: GL_MAX_3D_TEXTURE_SIZE_EXT = 32883
integer(kind=glenum), parameter, public :: GL_TEXTURE_3D_BINDING_EXT = 32874
#endif

! EXT_texture_object

#ifdef GL_EXT_texture_object
integer(kind=glenum), parameter, public :: GL_TEXTURE_PRIORITY_EXT = 32870
integer(kind=glenum), parameter, public :: GL_TEXTURE_RESIDENT_EXT = 32871
integer(kind=glenum), parameter, public :: GL_TEXTURE_1D_BINDING_EXT = 32872
integer(kind=glenum), parameter, public :: GL_TEXTURE_2D_BINDING_EXT = 32873
#endif

! EXT_vertex_array

#ifdef GL_EXT_vertex_array
integer(kind=glenum), parameter, public :: GL_VERTEX_ARRAY_EXT = 32884
integer(kind=glenum), parameter, public :: GL_NORMAL_ARRAY_EXT = 32885
integer(kind=glenum), parameter, public :: GL_COLOR_ARRAY_EXT = 32886
integer(kind=glenum), parameter, public :: GL_INDEX_ARRAY_EXT = 32887
integer(kind=glenum), parameter, public :: GL_TEXTURE_COORD_ARRAY_EXT = 32888
integer(kind=glenum), parameter, public :: GL_EDGE_FLAG_ARRAY_EXT = 32889
integer(kind=glenum), parameter, public :: GL_VERTEX_ARRAY_SIZE_EXT = 32890
integer(kind=glenum), parameter, public :: GL_VERTEX_ARRAY_TYPE_EXT = 32891
integer(kind=glenum), parameter, public :: GL_VERTEX_ARRAY_STRIDE_EXT = 32892
integer(kind=glenum), parameter, public :: GL_VERTEX_ARRAY_COUNT_EXT = 32893
integer(kind=glenum), parameter, public :: GL_NORMAL_ARRAY_TYPE_EXT = 32894
integer(kind=glenum), parameter, public :: GL_NORMAL_ARRAY_STRIDE_EXT = 32895
integer(kind=glenum), parameter, public :: GL_NORMAL_ARRAY_COUNT_EXT = 32896
integer(kind=glenum), parameter, public :: GL_COLOR_ARRAY_SIZE_EXT = 32897
integer(kind=glenum), parameter, public :: GL_COLOR_ARRAY_TYPE_EXT = 32898
integer(kind=glenum), parameter, public :: GL_COLOR_ARRAY_STRIDE_EXT = 32899
integer(kind=glenum), parameter, public :: GL_COLOR_ARRAY_COUNT_EXT = 32900
integer(kind=glenum), parameter, public :: GL_INDEX_ARRAY_TYPE_EXT = 32901
integer(kind=glenum), parameter, public :: GL_INDEX_ARRAY_STRIDE_EXT = 32902
integer(kind=glenum), parameter, public :: GL_INDEX_ARRAY_COUNT_EXT = 32903
integer(kind=glenum), parameter, public :: GL_TEXTURE_COORD_ARRAY_SIZE_EXT=32904
integer(kind=glenum), parameter, public :: GL_TEXTURE_COORD_ARRAY_TYPE_EXT=32905
integer(kind=glenum), parameter, public :: GL_TEXTURE_COORD_ARRAY_STRIDE_E=32906
integer(kind=glenum), parameter, public :: GL_TEXTURE_COORD_ARRAY_COUNT_EX=32907
integer(kind=glenum), parameter, public :: GL_EDGE_FLAG_ARRAY_STRIDE_EXT = 32908
integer(kind=glenum), parameter, public :: GL_EDGE_FLAG_ARRAY_COUNT_EXT = 32909
integer(kind=glenum), parameter, public :: GL_VERTEX_ARRAY_POINTER_EXT = 32910
integer(kind=glenum), parameter, public :: GL_NORMAL_ARRAY_POINTER_EXT = 32911
integer(kind=glenum), parameter, public :: GL_COLOR_ARRAY_POINTER_EXT = 32912
integer(kind=glenum), parameter, public :: GL_INDEX_ARRAY_POINTER_EXT = 32913
integer(kind=glenum), parameter, public :: GL_TEXTURE_COORD_ARRAY_POINTER_=32914
integer(kind=glenum), parameter, public :: GL_EDGE_FLAG_ARRAY_POINTER_EXT =32915
#endif

! SGIS_multisample

#ifdef GL_SGIS_multisample
integer(kind=glenum), parameter, public :: GL_MULTISAMPLE_SGIS = 32925
integer(kind=glenum), parameter, public :: GL_SAMPLE_ALPHA_TO_MASK_SGIS = 32926
integer(kind=glenum), parameter, public :: GL_SAMPLE_ALPHA_TO_ONE_SGIS = 32927
integer(kind=glenum), parameter, public :: GL_SAMPLE_MASK_SGIS = 32928
integer(kind=glenum), parameter, public :: GL_1PASS_SGIS = 32929
integer(kind=glenum), parameter, public :: GL_2PASS_0_SGIS = 32930
integer(kind=glenum), parameter, public :: GL_2PASS_1_SGIS = 32931
integer(kind=glenum), parameter, public :: GL_4PASS_0_SGIS = 32932
integer(kind=glenum), parameter, public :: GL_4PASS_1_SGIS = 32933
integer(kind=glenum), parameter, public :: GL_4PASS_2_SGIS = 32934
integer(kind=glenum), parameter, public :: GL_4PASS_3_SGIS = 32935
integer(kind=glenum), parameter, public :: GL_SAMPLE_BUFFERS_SGIS = 32936
integer(kind=glenum), parameter, public :: GL_SAMPLES_SGIS = 32937
integer(kind=glenum), parameter, public :: GL_SAMPLE_MASK_VALUE_SGIS = 32938
integer(kind=glenum), parameter, public :: GL_SAMPLE_MASK_INVERT_SGIS = 32939
integer(kind=glenum), parameter, public :: GL_SAMPLE_PATTERN_SGIS = 32940
#endif

!-----------------
! interface blocks
!-----------------

! OpenGL 1.0

interface glaccum
   module procedure f9xglaccum
end interface

interface glalphafunc
   module procedure f9xglalphafunc
end interface

interface glbegin
   module procedure f9xglbegin
end interface

interface glbitmap
   module procedure f9xglbitmap
end interface

interface glblendfunc
   module procedure f9xglblendfunc
end interface

interface glcalllist
   module procedure f9xglcalllist
end interface

interface glcalllists
   module procedure f9x0glcalllists
   module procedure f9x1glcalllists
   module procedure f9x4glcalllists
#ifndef NOBYTE
   module procedure f9x5glcalllists
#endif
#ifndef NOSHORT
   module procedure f9x6glcalllists
#endif
end interface

interface glclear
   module procedure f9xglclear
end interface

interface glclearaccum
   module procedure f9xglclearaccum
end interface

interface glclearcolor
   module procedure f9xglclearcolor
end interface

interface glcleardepth
   module procedure f9xglcleardepth
end interface

interface glclearindex
   module procedure f9xglclearindex
end interface

interface glclearstencil
   module procedure f9xglclearstencil
end interface

interface glclipplane
   module procedure f9xglclipplane
end interface

interface glcolor3b
   module procedure f9xglcolor3b
end interface

interface glcolor3bv
   module procedure f9xglcolor3bv
end interface

interface glcolor3d
   module procedure f9xglcolor3d
end interface

interface glcolor3dv
   module procedure f9xglcolor3dv
end interface

interface glcolor3f
   module procedure f9xglcolor3f
end interface

interface glcolor3fv
   module procedure f9xglcolor3fv
end interface

interface glcolor3i
   module procedure f9xglcolor3i
end interface

interface glcolor3iv
   module procedure f9xglcolor3iv
end interface

interface glcolor3s
   module procedure f9xglcolor3s
end interface

interface glcolor3sv
   module procedure f9xglcolor3sv
end interface

interface glcolor3ub
   module procedure f9xglcolor3ub
end interface

interface glcolor3ubv
   module procedure f9xglcolor3ubv
end interface

interface glcolor3ui
   module procedure f9xglcolor3ui
end interface

interface glcolor3uiv
   module procedure f9xglcolor3uiv
end interface

interface glcolor3us
   module procedure f9xglcolor3us
end interface

interface glcolor3usv
   module procedure f9xglcolor3usv
end interface

interface glcolor4b
   module procedure f9xglcolor4b
end interface

interface glcolor4bv
   module procedure f9xglcolor4bv
end interface

interface glcolor4d
   module procedure f9xglcolor4d
end interface

interface glcolor4dv
   module procedure f9xglcolor4dv
end interface

interface glcolor4f
   module procedure f9xglcolor4f
end interface

interface glcolor4fv
   module procedure f9xglcolor4fv
end interface

interface glcolor4i
   module procedure f9xglcolor4i
end interface

interface glcolor4iv
   module procedure f9xglcolor4iv
end interface

interface glcolor4s
   module procedure f9xglcolor4s
end interface

interface glcolor4sv
   module procedure f9xglcolor4sv
end interface

interface glcolor4ub
   module procedure f9xglcolor4ub
end interface

interface glcolor4ubv
   module procedure f9xglcolor4ubv
end interface

interface glcolor4ui
   module procedure f9xglcolor4ui
end interface

interface glcolor4uiv
   module procedure f9xglcolor4uiv
end interface

interface glcolor4us
   module procedure f9xglcolor4us
end interface

interface glcolor4usv
   module procedure f9xglcolor4usv
end interface

interface glcolormask
   module procedure f9xglcolormask
end interface

interface glcolormaterial
   module procedure f9xglcolormaterial
end interface

interface glcopypixels
   module procedure f9xglcopypixels
end interface

interface glcullface
   module procedure f9xglcullface
end interface

interface gldeletelists
   module procedure f9xgldeletelists
end interface

interface gldepthfunc
   module procedure f9xgldepthfunc
end interface

interface gldepthmask
   module procedure f9xgldepthmask
end interface

interface gldepthrange
   module procedure f9xgldepthrange
end interface

interface gldisable
   module procedure f9xgldisable
end interface

interface gldrawbuffer
   module procedure f9xgldrawbuffer
end interface

interface gldrawpixels
   module procedure f9x0gldrawpixels
   module procedure f9x1gldrawpixels
   module procedure f9x4gldrawpixels
#ifndef NOBYTE
   module procedure f9x5gldrawpixels
#endif
#ifndef NOSHORT
   module procedure f9x6gldrawpixels
#endif
end interface

interface gledgeflag
   module procedure f9xgledgeflag
end interface

interface gledgeflagv
   module procedure f9xgledgeflagv
end interface

interface glenable
   module procedure f9xglenable
end interface

interface glend
   module procedure f9xglend
end interface

interface glendlist
   module procedure f9xglendlist
end interface

interface glevalcoord1d
   module procedure f9xglevalcoord1d
end interface

interface glevalcoord1dv
   module procedure f9xglevalcoord1dv
end interface

interface glevalcoord1f
   module procedure f9xglevalcoord1f
end interface

interface glevalcoord1fv
   module procedure f9xglevalcoord1fv
end interface

interface glevalcoord2d
   module procedure f9xglevalcoord2d
end interface

interface glevalcoord2dv
   module procedure f9xglevalcoord2dv
end interface

interface glevalcoord2f
   module procedure f9xglevalcoord2f
end interface

interface glevalcoord2fv
   module procedure f9xglevalcoord2fv
end interface

interface glevalmesh1
   module procedure f9xglevalmesh1
end interface

interface glevalmesh2
   module procedure f9xglevalmesh2
end interface

interface glevalpoint1
   module procedure f9xglevalpoint1
end interface

interface glevalpoint2
   module procedure f9xglevalpoint2
end interface

interface glfeedbackbuffer
   module procedure f9xglfeedbackbuffer
end interface

interface glfinish
   module procedure f9xglfinish
end interface

interface glflush
   module procedure f9xglflush
end interface

interface glfogf
   module procedure f9xglfogf
end interface

interface glfogfv
   module procedure f9xglfogfv
end interface

interface glfogi
   module procedure f9xglfogi
end interface

interface glfogiv
   module procedure f9xglfogiv
end interface

interface glfrontface
   module procedure f9xglfrontface
end interface

interface glfrustum
   module procedure f9xglfrustum
end interface

interface glgenlists
   module procedure f9xglgenlists
end interface

interface glgetbooleanv
   module procedure f9xglgetbooleanv
end interface

interface glgetclipplane
   module procedure f9xglgetclipplane
end interface

interface glgetdoublev
   module procedure f9xglgetdoublev
end interface

interface glgeterror
   module procedure f9xglgeterror
end interface

interface glgetfloatv
   module procedure f9xglgetfloatv
end interface

interface glgetintegerv
   module procedure f9xglgetintegerv
end interface

interface glgetlightfv
   module procedure f9xglgetlightfv
end interface

interface glgetlightiv
   module procedure f9xglgetlightiv
end interface

interface glgetmapdv
   module procedure f9xglgetmapdv
end interface

interface glgetmapfv
   module procedure f9xglgetmapfv
end interface

interface glgetmapiv
   module procedure f9xglgetmapiv
end interface

interface glgetmaterialfv
   module procedure f9xglgetmaterialfv
end interface

interface glgetmaterialiv
   module procedure f9xglgetmaterialiv
end interface

interface glgetpixelmapfv
   module procedure f9xglgetpixelmapfv
end interface

interface glgetpixelmapuiv
   module procedure f9xglgetpixelmapuiv
end interface

interface glgetpixelmapusv
   module procedure f9xglgetpixelmapusv
end interface

interface glgetpolygonstipple
   module procedure f9xglgetpolygonstipple
end interface

interface glgetstring
   module procedure f9xglgetstring
end interface

interface glgettexenvfv
   module procedure f9xglgettexenvfv
end interface

interface glgettexenviv
   module procedure f9xglgettexenviv
end interface

interface glgettexgendv
   module procedure f9xglgettexgendv
end interface

interface glgettexgenfv
   module procedure f9xglgettexgenfv
end interface

interface glgettexgeniv
   module procedure f9xglgettexgeniv
end interface

interface glgetteximage
   module procedure f9x0glgetteximage
   module procedure f9x1glgetteximage
   module procedure f9x4glgetteximage
#ifndef NOBYTE
   module procedure f9x5glgetteximage
#endif
#ifndef NOSHORT
   module procedure f9x6glgetteximage
#endif
end interface

interface glgettexlevelparameterfv
   module procedure f9xglgettexlevelparameterfv
end interface

interface glgettexlevelparameteriv
   module procedure f9xglgettexlevelparameteriv
end interface

interface glgettexparameterfv
   module procedure f9xglgettexparameterfv
end interface

interface glgettexparameteriv
   module procedure f9xglgettexparameteriv
end interface

interface glhint
   module procedure f9xglhint
end interface

interface glindexmask
   module procedure f9xglindexmask
end interface

interface glindexd
   module procedure f9xglindexd
end interface

interface glindexdv
   module procedure f9xglindexdv
end interface

interface glindexf
   module procedure f9xglindexf
end interface

interface glindexfv
   module procedure f9xglindexfv
end interface

interface glindexi
   module procedure f9xglindexi
end interface

interface glindexiv
   module procedure f9xglindexiv
end interface

interface glindexs
   module procedure f9xglindexs
end interface

interface glindexsv
   module procedure f9xglindexsv
end interface

interface glinitnames
   module procedure f9xglinitnames
end interface

interface glisenabled
   module procedure f9xglisenabled
end interface

interface glislist
   module procedure f9xglislist
end interface

interface gllightmodelf
   module procedure f9xgllightmodelf
end interface

interface gllightmodelfv
   module procedure f9xgllightmodelfv
end interface

interface gllightmodeli
   module procedure f9xgllightmodeli
end interface

interface gllightmodeliv
   module procedure f9xgllightmodeliv
end interface

interface gllightf
   module procedure f9xgllightf
end interface

interface gllightfv
   module procedure f9xgllightfv
end interface

interface gllighti
   module procedure f9xgllighti
end interface

interface gllightiv
   module procedure f9xgllightiv
end interface

interface gllinestipple
   module procedure f9xgllinestipple
end interface

interface gllinewidth
   module procedure f9xgllinewidth
end interface

interface gllistbase
   module procedure f9xgllistbase
end interface

interface glloadidentity
   module procedure f9xglloadidentity
end interface

interface glloadmatrixd
   module procedure f9xglloadmatrixd
end interface

interface glloadmatrixf
   module procedure f9xglloadmatrixf
end interface

interface glloadname
   module procedure f9xglloadname
end interface

interface gllogicop
   module procedure f9xgllogicop
end interface

interface glmap1d
   module procedure f9xglmap1d
end interface

interface glmap1f
   module procedure f9xglmap1f
end interface

interface glmap2d
   module procedure f9xglmap2d
end interface

interface glmap2f
   module procedure f9xglmap2f
end interface

interface glmapgrid1d
   module procedure f9xglmapgrid1d
end interface

interface glmapgrid1f
   module procedure f9xglmapgrid1f
end interface

interface glmapgrid2d
   module procedure f9xglmapgrid2d
end interface

interface glmapgrid2f
   module procedure f9xglmapgrid2f
end interface

interface glmaterialf
   module procedure f9xglmaterialf
end interface

interface glmaterialfv
   module procedure f9xglmaterialfv
end interface

interface glmateriali
   module procedure f9xglmateriali
end interface

interface glmaterialiv
   module procedure f9xglmaterialiv
end interface

interface glmatrixmode
   module procedure f9xglmatrixmode
end interface

interface glmultmatrixd
   module procedure f9xglmultmatrixd
end interface

interface glmultmatrixf
   module procedure f9xglmultmatrixf
end interface

interface glnewlist
   module procedure f9xglnewlist
end interface

interface glnormal3b
   module procedure f9xglnormal3b
end interface

interface glnormal3bv
   module procedure f9xglnormal3bv
end interface

interface glnormal3d
   module procedure f9xglnormal3d
end interface

interface glnormal3dv
   module procedure f9xglnormal3dv
end interface

interface glnormal3f
   module procedure f9xglnormal3f
end interface

interface glnormal3fv
   module procedure f9xglnormal3fv
end interface

interface glnormal3i
   module procedure f9xglnormal3i
end interface

interface glnormal3iv
   module procedure f9xglnormal3iv
end interface

interface glnormal3s
   module procedure f9xglnormal3s
end interface

interface glnormal3sv
   module procedure f9xglnormal3sv
end interface

interface glortho
   module procedure f9xglortho
end interface

interface glpassthrough
   module procedure f9xglpassthrough
end interface

interface glpixelmapfv
   module procedure f9xglpixelmapfv
end interface

interface glpixelmapuiv
   module procedure f9xglpixelmapuiv
end interface

interface glpixelmapusv
   module procedure f9xglpixelmapusv
end interface

interface glpixelstoref
   module procedure f9xglpixelstoref
end interface

interface glpixelstorei
   module procedure f9xglpixelstorei
end interface

interface glpixeltransferf
   module procedure f9xglpixeltransferf
end interface

interface glpixeltransferi
   module procedure f9xglpixeltransferi
end interface

interface glpixelzoom
   module procedure f9xglpixelzoom
end interface

interface glpointsize
   module procedure f9xglpointsize
end interface

interface glpolygonmode
   module procedure f9xglpolygonmode
end interface

interface glpolygonstipple
   module procedure f9xglpolygonstipple
end interface

interface glpopattrib
   module procedure f9xglpopattrib
end interface

interface glpopmatrix
   module procedure f9xglpopmatrix
end interface

interface glpopname
   module procedure f9xglpopname
end interface

interface glpushattrib
   module procedure f9xglpushattrib
end interface

interface glpushmatrix
   module procedure f9xglpushmatrix
end interface

interface glpushname
   module procedure f9xglpushname
end interface

interface glrasterpos2d
   module procedure f9xglrasterpos2d
end interface

interface glrasterpos2dv
   module procedure f9xglrasterpos2dv
end interface

interface glrasterpos2f
   module procedure f9xglrasterpos2f
end interface

interface glrasterpos2fv
   module procedure f9xglrasterpos2fv
end interface

interface glrasterpos2i
   module procedure f9xglrasterpos2i
end interface

interface glrasterpos2iv
   module procedure f9xglrasterpos2iv
end interface

interface glrasterpos2s
   module procedure f9xglrasterpos2s
end interface

interface glrasterpos2sv
   module procedure f9xglrasterpos2sv
end interface

interface glrasterpos3d
   module procedure f9xglrasterpos3d
end interface

interface glrasterpos3dv
   module procedure f9xglrasterpos3dv
end interface

interface glrasterpos3f
   module procedure f9xglrasterpos3f
end interface

interface glrasterpos3fv
   module procedure f9xglrasterpos3fv
end interface

interface glrasterpos3i
   module procedure f9xglrasterpos3i
end interface

interface glrasterpos3iv
   module procedure f9xglrasterpos3iv
end interface

interface glrasterpos3s
   module procedure f9xglrasterpos3s
end interface

interface glrasterpos3sv
   module procedure f9xglrasterpos3sv
end interface

interface glrasterpos4d
   module procedure f9xglrasterpos4d
end interface

interface glrasterpos4dv
   module procedure f9xglrasterpos4dv
end interface

interface glrasterpos4f
   module procedure f9xglrasterpos4f
end interface

interface glrasterpos4fv
   module procedure f9xglrasterpos4fv
end interface

interface glrasterpos4i
   module procedure f9xglrasterpos4i
end interface

interface glrasterpos4iv
   module procedure f9xglrasterpos4iv
end interface

interface glrasterpos4s
   module procedure f9xglrasterpos4s
end interface

interface glrasterpos4sv
   module procedure f9xglrasterpos4sv
end interface

interface glreadbuffer
   module procedure f9xglreadbuffer
end interface

interface glreadpixels
   module procedure f9x0glreadpixels
   module procedure f9x1glreadpixels
   module procedure f9x4glreadpixels
#ifndef NOBYTE
   module procedure f9x5glreadpixels
#endif
#ifndef NOSHORT
   module procedure f9x6glreadpixels
#endif
end interface

interface glrectd
   module procedure f9xglrectd
end interface

interface glrectdv
   module procedure f9xglrectdv
end interface

interface glrectf
   module procedure f9xglrectf
end interface

interface glrectfv
   module procedure f9xglrectfv
end interface

interface glrecti
   module procedure f9xglrecti
end interface

interface glrectiv
   module procedure f9xglrectiv
end interface

interface glrects
   module procedure f9xglrects
end interface

interface glrectsv
   module procedure f9xglrectsv
end interface

interface glrendermode
   module procedure f9xglrendermode
end interface

interface glrotated
   module procedure f9xglrotated
end interface

interface glrotatef
   module procedure f9xglrotatef
end interface

interface glscaled
   module procedure f9xglscaled
end interface

interface glscalef
   module procedure f9xglscalef
end interface

interface glscissor
   module procedure f9xglscissor
end interface

interface glselectbuffer
   module procedure f9xglselectbuffer
end interface

interface glshademodel
   module procedure f9xglshademodel
end interface

interface glstencilfunc
   module procedure f9xglstencilfunc
end interface

interface glstencilmask
   module procedure f9xglstencilmask
end interface

interface glstencilop
   module procedure f9xglstencilop
end interface

interface gltexcoord1d
   module procedure f9xgltexcoord1d
end interface

interface gltexcoord1dv
   module procedure f9xgltexcoord1dv
end interface

interface gltexcoord1f
   module procedure f9xgltexcoord1f
end interface

interface gltexcoord1fv
   module procedure f9xgltexcoord1fv
end interface

interface gltexcoord1i
   module procedure f9xgltexcoord1i
end interface

interface gltexcoord1iv
   module procedure f9xgltexcoord1iv
end interface

interface gltexcoord1s
   module procedure f9xgltexcoord1s
end interface

interface gltexcoord1sv
   module procedure f9xgltexcoord1sv
end interface

interface gltexcoord2d
   module procedure f9xgltexcoord2d
end interface

interface gltexcoord2dv
   module procedure f9xgltexcoord2dv
end interface

interface gltexcoord2f
   module procedure f9xgltexcoord2f
end interface

interface gltexcoord2fv
   module procedure f9xgltexcoord2fv
end interface

interface gltexcoord2i
   module procedure f9xgltexcoord2i
end interface

interface gltexcoord2iv
   module procedure f9xgltexcoord2iv
end interface

interface gltexcoord2s
   module procedure f9xgltexcoord2s
end interface

interface gltexcoord2sv
   module procedure f9xgltexcoord2sv
end interface

interface gltexcoord3d
   module procedure f9xgltexcoord3d
end interface

interface gltexcoord3dv
   module procedure f9xgltexcoord3dv
end interface

interface gltexcoord3f
   module procedure f9xgltexcoord3f
end interface

interface gltexcoord3fv
   module procedure f9xgltexcoord3fv
end interface

interface gltexcoord3i
   module procedure f9xgltexcoord3i
end interface

interface gltexcoord3iv
   module procedure f9xgltexcoord3iv
end interface

interface gltexcoord3s
   module procedure f9xgltexcoord3s
end interface

interface gltexcoord3sv
   module procedure f9xgltexcoord3sv
end interface

interface gltexcoord4d
   module procedure f9xgltexcoord4d
end interface

interface gltexcoord4dv
   module procedure f9xgltexcoord4dv
end interface

interface gltexcoord4f
   module procedure f9xgltexcoord4f
end interface

interface gltexcoord4fv
   module procedure f9xgltexcoord4fv
end interface

interface gltexcoord4i
   module procedure f9xgltexcoord4i
end interface

interface gltexcoord4iv
   module procedure f9xgltexcoord4iv
end interface

interface gltexcoord4s
   module procedure f9xgltexcoord4s
end interface

interface gltexcoord4sv
   module procedure f9xgltexcoord4sv
end interface

interface gltexenvf
   module procedure f9xgltexenvf
end interface

interface gltexenvfv
   module procedure f9xgltexenvfv
end interface

interface gltexenvi
   module procedure f9xgltexenvi
end interface

interface gltexenviv
   module procedure f9xgltexenviv
end interface

interface gltexgend
   module procedure f9xgltexgend
end interface

interface gltexgendv
   module procedure f9xgltexgendv
end interface

interface gltexgenf
   module procedure f9xgltexgenf
end interface

interface gltexgenfv
   module procedure f9xgltexgenfv
end interface

interface gltexgeni
   module procedure f9xgltexgeni
end interface

interface gltexgeniv
   module procedure f9xgltexgeniv
end interface

interface glteximage1d
   module procedure f9x0glteximage1d
   module procedure f9x1glteximage1d
   module procedure f9x4glteximage1d
#ifndef NOBYTE
   module procedure f9x5glteximage1d
#endif
#ifndef NOSHORT
   module procedure f9x6glteximage1d
#endif
end interface

interface glteximage2d
   module procedure f9x0glteximage2d
   module procedure f9x1glteximage2d
   module procedure f9x4glteximage2d
#ifndef NOBYTE
   module procedure f9x5glteximage2d
#endif
#ifndef NOSHORT
   module procedure f9x6glteximage2d
#endif
end interface

interface gltexparameterf
   module procedure f9xgltexparameterf
end interface

interface gltexparameterfv
   module procedure f9xgltexparameterfv
end interface

interface gltexparameteri
   module procedure f9xgltexparameteri
end interface

interface gltexparameteriv
   module procedure f9xgltexparameteriv
end interface

interface gltranslated
   module procedure f9xgltranslated
end interface

interface gltranslatef
   module procedure f9xgltranslatef
end interface

interface glvertex2d
   module procedure f9xglvertex2d
end interface

interface glvertex2dv
   module procedure f9xglvertex2dv
end interface

interface glvertex2f
   module procedure f9xglvertex2f
end interface

interface glvertex2fv
   module procedure f9xglvertex2fv
end interface

interface glvertex2i
   module procedure f9xglvertex2i
end interface

interface glvertex2iv
   module procedure f9xglvertex2iv
end interface

interface glvertex2s
   module procedure f9xglvertex2s
end interface

interface glvertex2sv
   module procedure f9xglvertex2sv
end interface

interface glvertex3d
   module procedure f9xglvertex3d
end interface

interface glvertex3dv
   module procedure f9xglvertex3dv
end interface

interface glvertex3f
   module procedure f9xglvertex3f
end interface

interface glvertex3fv
   module procedure f9xglvertex3fv
end interface

interface glvertex3i
   module procedure f9xglvertex3i
end interface

interface glvertex3iv
   module procedure f9xglvertex3iv
end interface

interface glvertex3s
   module procedure f9xglvertex3s
end interface

interface glvertex3sv
   module procedure f9xglvertex3sv
end interface

interface glvertex4d
   module procedure f9xglvertex4d
end interface

interface glvertex4dv
   module procedure f9xglvertex4dv
end interface

interface glvertex4f
   module procedure f9xglvertex4f
end interface

interface glvertex4fv
   module procedure f9xglvertex4fv
end interface

interface glvertex4i
   module procedure f9xglvertex4i
end interface

interface glvertex4iv
   module procedure f9xglvertex4iv
end interface

interface glvertex4s
   module procedure f9xglvertex4s
end interface

interface glvertex4sv
   module procedure f9xglvertex4sv
end interface

interface glviewport
   module procedure f9xglviewport
end interface

! OpenGL 1.1

#ifdef OPENGL_1_1

interface glaretexturesresident
   module procedure f9xglaretexturesresident
end interface

interface glarrayelement
   module procedure f9xglarrayelement
end interface

interface glbindtexture
   module procedure f9xglbindtexture
end interface

interface glcolorpointer
   module procedure f9x0glcolorpointer
   module procedure f9x1glcolorpointer
#ifndef NOBYTE
   module procedure f9x5glcolorpointer
#endif
#ifndef NOSHORT
   module procedure f9x6glcolorpointer
#endif
   module procedure f9x7glcolorpointer
   module procedure f9x4glcolorpointer
end interface

interface glcopyteximage1d
   module procedure f9xglcopyteximage1d
end interface

interface glcopyteximage2d
   module procedure f9xglcopyteximage2d
end interface

interface glcopytexsubimage1d
   module procedure f9xglcopytexsubimage1d
end interface

interface glcopytexsubimage2d
   module procedure f9xglcopytexsubimage2d
end interface

interface gldeletetextures
   module procedure f9xgldeletetextures
end interface

interface gldisableclientstate
   module procedure f9xgldisableclientstate
end interface

interface gldrawarrays
   module procedure f9xgldrawarrays
end interface

interface gldrawelements
   module procedure f9x0gldrawelements
   module procedure f9x4gldrawelements
#ifndef NOBYTE
   module procedure f9x5gldrawelements
#endif
#ifndef NOSHORT
   module procedure f9x6gldrawelements
#endif
end interface

interface gledgeflagpointer
   module procedure f9xgledgeflagpointer
   module procedure f9xegledgeflagpointer
end interface

interface glenableclientstate
   module procedure f9xglenableclientstate
end interface

interface glgentextures
   module procedure f9xglgentextures
end interface

interface glgetpointerv
   module procedure f9xglgetpointerv
end interface

interface glindexpointer
   module procedure f9x0glindexpointer
   module procedure f9x1glindexpointer
#ifndef NOSHORT
   module procedure f9x6glindexpointer
#endif
   module procedure f9x7glindexpointer
   module procedure f9x4glindexpointer
end interface

interface glindexub
   module procedure f9xglindexub
end interface

interface glindexubv
   module procedure f9xglindexubv
end interface

interface glinterleavedarrays
   module procedure f9x0glinterleavedarrays
   module procedure f9x1glinterleavedarrays
   module procedure f9x4glinterleavedarrays
#ifndef NOBYTE
   module procedure f9x5glinterleavedarrays
#endif
#ifndef NOSHORT
   module procedure f9x6glinterleavedarrays
#endif
end interface

interface glistexture
   module procedure f9xglistexture
end interface

interface glnormalpointer
   module procedure f9x0glnormalpointer
   module procedure f9x1glnormalpointer
#ifndef NOBYTE
   module procedure f9x5glnormalpointer
#endif
#ifndef NOSHORT
   module procedure f9x6glnormalpointer
#endif
   module procedure f9x7glnormalpointer
   module procedure f9x4glnormalpointer
end interface

interface glpolygonoffset
   module procedure f9xglpolygonoffset
end interface

interface glpopclientattrib
   module procedure f9xglpopclientattrib
end interface

interface glprioritizetextures
   module procedure f9xglprioritizetextures
end interface

interface glpushclientattrib
   module procedure f9xglpushclientattrib
end interface

interface gltexcoordpointer
   module procedure f9x0gltexcoordpointer
   module procedure f9x1gltexcoordpointer
#ifndef NOSHORT
   module procedure f9x6gltexcoordpointer
#endif
   module procedure f9x7gltexcoordpointer
   module procedure f9x4gltexcoordpointer
end interface

interface gltexsubimage1d
   module procedure f9x0gltexsubimage1d
   module procedure f9x1gltexsubimage1d
   module procedure f9x4gltexsubimage1d
#ifndef NOBYTE
   module procedure f9x5gltexsubimage1d
#endif
#ifndef NOSHORT
   module procedure f9x6gltexsubimage1d
#endif
end interface

interface gltexsubimage2d
   module procedure f9x0gltexsubimage2d
   module procedure f9x1gltexsubimage2d
   module procedure f9x4gltexsubimage2d
#ifndef NOBYTE
   module procedure f9x5gltexsubimage2d
#endif
#ifndef NOSHORT
   module procedure f9x6gltexsubimage2d
#endif
end interface

interface glvertexpointer
   module procedure f9x0glvertexpointer
   module procedure f9x1glvertexpointer
#ifndef NOSHORT
   module procedure f9x6glvertexpointer
#endif
   module procedure f9x7glvertexpointer
   module procedure f9x4glvertexpointer
end interface

#endif

! OpenGL 1.2

#ifdef OPENGL_1_2

interface gldrawrangeelements
   module procedure f9x0gldrawrangeelements
   module procedure f9x4gldrawrangeelements
#ifndef NOBYTE
   module procedure f9x5gldrawrangeelements
#endif
#ifndef NOSHORT
   module procedure f9x6gldrawrangeelements
#endif
end interface

interface glteximage3d
   module procedure f9x0glteximage3d
   module procedure f9x1glteximage3d
   module procedure f9x4glteximage3d
#ifndef NOBYTE
   module procedure f9x5glteximage3d
#endif
#ifndef NOSHORT
   module procedure f9x6glteximage3d
#endif
end interface

interface gltexsubimage3d
   module procedure f9x0gltexsubimage3d
   module procedure f9x1gltexsubimage3d
   module procedure f9x4gltexsubimage3d
#ifndef NOBYTE
   module procedure f9x5gltexsubimage3d
#endif
#ifndef NOSHORT
   module procedure f9x6gltexsubimage3d
#endif
end interface

interface glcopytexsubimage3d
   module procedure f9xglcopytexsubimage3d
end interface

#endif

! EXT_blend_color

#ifdef GL_EXT_blend_color

interface glblendcolorext
   module procedure f9xglblendcolorext
end interface

#endif

! EXT_blend_minmax

#ifdef GL_EXT_blend_minmax

interface glblendequationext
   module procedure f9xglblendequationext
end interface

#endif

! EXT_color_table

#ifdef GL_EXT_color_table

interface glcolortableext
   module procedure f9x0glcolortableext
   module procedure f9x1glcolortableext
   module procedure f9x4glcolortableext
#ifndef NOBYTE
   module procedure f9x5glcolortableext
#endif
#ifndef NOSHORT
   module procedure f9x6glcolortableext
#endif
end interface

interface glcolorsubtableext
   module procedure f9x0glcolorsubtableext
   module procedure f9x1glcolorsubtableext
   module procedure f9x4glcolorsubtableext
#ifndef NOBYTE
   module procedure f9x5glcolorsubtableext
#endif
#ifndef NOSHORT
   module procedure f9x6glcolorsubtableext
#endif
end interface

interface glgetcolortableext
   module procedure f9x0glgetcolortableext
   module procedure f9x1glgetcolortableext
   module procedure f9x4glgetcolortableext
#ifndef NOBYTE
   module procedure f9x5glgetcolortableext
#endif
#ifndef NOSHORT
   module procedure f9x6glgetcolortableext
#endif
end interface

interface glgetcolortableparameterfvext
   module procedure f9xglgetcolortableparamfvext
end interface

interface glgetcolortableparameterivext
   module procedure f9xglgetcolortableparamivext
end interface

#endif

! EXT_point_parameters

#ifdef GL_EXT_point_parameters

interface glpointparameterfext
   module procedure f9xglpointparameterfext
end interface

interface glpointparameterfvext
   module procedure f9xglpointparameterfvext
end interface

#endif

! EXT_polygon_offset

#ifdef GL_EXT_polygon_offset

interface glpolygonoffsetext
   module procedure f9xglpolygonoffsetext
end interface

#endif

! EXT_subtexture

#ifdef GL_EXT_subtexture

interface gltexsubimage1dext
   module procedure f9x0gltexsubimage1dext
   module procedure f9x1gltexsubimage1dext
   module procedure f9x4gltexsubimage1dext
#ifndef NOBYTE
   module procedure f9x5gltexsubimage1dext
#endif
#ifndef NOSHORT
   module procedure f9x6gltexsubimage1dext
#endif
end interface

interface gltexsubimage2dext
   module procedure f9x0gltexsubimage2dext
   module procedure f9x1gltexsubimage2dext
   module procedure f9x4gltexsubimage2dext
#ifndef NOBYTE
   module procedure f9x5gltexsubimage2dext
#endif
#ifndef NOSHORT
   module procedure f9x6gltexsubimage2dext
#endif
end interface

#endif

! EXT_texture3D

#ifdef GL_EXT_texture3D

interface glteximage3dext
   module procedure f9x0glteximage3dext
   module procedure f9x1glteximage3dext
   module procedure f9x4glteximage3dext
#ifndef NOBYTE
   module procedure f9x5glteximage3dext
#endif
#ifndef NOSHORT
   module procedure f9x6glteximage3dext
#endif
end interface

interface gltexsubimage3dext
   module procedure f9x0gltexsubimage3dext
   module procedure f9x1gltexsubimage3dext
   module procedure f9x4gltexsubimage3dext
#ifndef NOBYTE
   module procedure f9x5gltexsubimage3dext
#endif
#ifndef NOSHORT
   module procedure f9x6gltexsubimage3dext
#endif
end interface

interface glcopytexsubimage3dext
   module procedure f9xglcopytexsubimage3dext
end interface

#endif

! EXT_texture_object

#ifdef GL_EXT_texture_object

interface glgentexturesext
   module procedure f9xglgentexturesext
end interface

interface gldeletetexturesext
   module procedure f9xgldeletetexturesext
end interface

interface glbindtextureext
   module procedure f9xglbindtextureext
end interface

interface glprioritizetexturesext
   module procedure f9xglprioritizetexturesext
end interface

interface glaretexturesresidentext
   module procedure f9xglaretexturesresidentext
end interface

interface glistextureext
   module procedure f9xglistextureext
end interface

#endif

! EXT_vertex_array

#ifdef GL_EXT_vertex_array

interface glvertexpointerext
   module procedure f9x0glvertexpointerext
   module procedure f9x1glvertexpointerext
#ifndef NOSHORT
   module procedure f9x6glvertexpointerext
#endif
   module procedure f9x7glvertexpointerext
   module procedure f9x4glvertexpointerext
end interface

interface glnormalpointerext
   module procedure f9x0glnormalpointerext
   module procedure f9x1glnormalpointerext
#ifndef NOBYTE
   module procedure f9x5glnormalpointerext
#endif
#ifndef NOSHORT
   module procedure f9x6glnormalpointerext
#endif
   module procedure f9x7glnormalpointerext
   module procedure f9x4glnormalpointerext
end interface

interface glcolorpointerext
   module procedure f9x0glcolorpointerext
   module procedure f9x1glcolorpointerext
#ifndef NOBYTE
   module procedure f9x5glcolorpointerext
#endif
#ifndef NOSHORT
   module procedure f9x6glcolorpointerext
#endif
   module procedure f9x7glcolorpointerext
   module procedure f9x4glcolorpointerext
end interface

interface glindexpointerext
   module procedure f9x0glindexpointerext
   module procedure f9x1glindexpointerext
#ifndef NOSHORT
   module procedure f9x6glindexpointerext
#endif
   module procedure f9x7glindexpointerext
   module procedure f9x4glindexpointerext
end interface

interface gltexcoordpointerext
   module procedure f9x0gltexcoordpointerext
   module procedure f9x1gltexcoordpointerext
#ifndef NOSHORT
   module procedure f9x6gltexcoordpointerext
#endif
   module procedure f9x7gltexcoordpointerext
   module procedure f9x4gltexcoordpointerext
end interface

interface gledgeflagpointerext
   module procedure f9xgledgeflagpointerext
   module procedure f9xegledgeflagpointerext
end interface

interface glgetpointervext
   module procedure f9xglgetpointervext
end interface

interface glarrayelementext
   module procedure f9xglarrayelementext
end interface

interface gldrawarraysext
   module procedure f9xgldrawarraysext
end interface

#endif

! SGIS_multisample

#ifdef GL_SGIS_multisample

interface glsamplemasksgis
   module procedure f9xglsamplemasksgis
end interface

interface glsamplepatternsgis
   module procedure f9xglsamplepatternsgis
end interface

#endif

! MESA_resize_buffers

#ifdef GL_MESA_resize_buffers

interface glresizebuffersmesa
   module procedure f9xglresizebuffersmesa
end interface

#endif

! MESA_window_pos

#ifdef GL_MESA_window_pos

interface glwindowpos2imesa
   module procedure f9xglwindowpos2imesa
end interface

interface glwindowpos2smesa
   module procedure f9xglwindowpos2smesa
end interface

interface glwindowpos2fmesa
   module procedure f9xglwindowpos2fmesa
end interface

interface glwindowpos2dmesa
   module procedure f9xglwindowpos2dmesa
end interface

interface glwindowpos2ivmesa
   module procedure f9xglwindowpos2ivmesa
end interface

interface glwindowpos2svmesa
   module procedure f9xglwindowpos2svmesa
end interface

interface glwindowpos2fvmesa
   module procedure f9xglwindowpos2fvmesa
end interface

interface glwindowpos2dvmesa
   module procedure f9xglwindowpos2dvmesa
end interface

interface glwindowpos3imesa
   module procedure f9xglwindowpos3imesa
end interface

interface glwindowpos3smesa
   module procedure f9xglwindowpos3smesa
end interface

interface glwindowpos3fmesa
   module procedure f9xglwindowpos3fmesa
end interface

interface glwindowpos3dmesa
   module procedure f9xglwindowpos3dmesa
end interface

interface glwindowpos3ivmesa
   module procedure f9xglwindowpos3ivmesa
end interface

interface glwindowpos3svmesa
   module procedure f9xglwindowpos3svmesa
end interface

interface glwindowpos3fvmesa
   module procedure f9xglwindowpos3fvmesa
end interface

interface glwindowpos3dvmesa
   module procedure f9xglwindowpos3dvmesa
end interface

interface glwindowpos4imesa
   module procedure f9xglwindowpos4imesa
end interface

interface glwindowpos4smesa
   module procedure f9xglwindowpos4smesa
end interface

interface glwindowpos4fmesa
   module procedure f9xglwindowpos4fmesa
end interface

interface glwindowpos4dmesa
   module procedure f9xglwindowpos4dmesa
end interface

interface glwindowpos4ivmesa
   module procedure f9xglwindowpos4ivmesa
end interface

interface glwindowpos4svmesa
   module procedure f9xglwindowpos4svmesa
end interface

interface glwindowpos4fvmesa
   module procedure f9xglwindowpos4fvmesa
end interface

interface glwindowpos4dvmesa
   module procedure f9xglwindowpos4dvmesa
end interface

#endif

interface size_glint
   module procedure size_glint_glint
#ifndef NOSHORT
   module procedure size_glint_glshort
#endif
#ifndef NOBYTE
   module procedure size_glint_glbyte
#endif
   module procedure size_glint_glfloat
end interface

!------------------
! Module procedures
!------------------

contains

function size_glint_glint(a) result(ressize_glint)
integer(kind=glint), dimension(:) INTENT_IN a
integer(kind=glint) :: ressize_glint
ressize_glint = size(a)
end function size_glint_glint

function size_glint_glshort(a) result(ressize_glint)
integer(kind=glshort), dimension(:) INTENT_IN a
integer(kind=glint) :: ressize_glint
ressize_glint = size(a)
end function size_glint_glshort

function size_glint_glbyte(a) result(ressize_glint)
integer(kind=glbyte), dimension(:) INTENT_IN a
integer(kind=glint) :: ressize_glint
ressize_glint = size(a)
end function size_glint_glbyte

function size_glint_glfloat(a) result(ressize_glint)
real(kind=glfloat), dimension(:) INTENT_IN a
integer(kind=glint) :: ressize_glint
ressize_glint = size(a)
end function size_glint_glfloat

! OpenGL 1.0

subroutine f9xglaccum(op,value)
integer(kind=glenum) INTENT_IN op
real(kind=glfloat) INTENT_IN value
call fglaccum(op,value)
return
end subroutine f9xglaccum

subroutine f9xglalphafunc(func,ref)
integer(kind=glenum) INTENT_IN func
real(kind=glclampf) INTENT_IN ref
call fglalphafunc(func,ref)
return
end subroutine f9xglalphafunc

subroutine f9xglbegin(mode)
integer(kind=glenum) INTENT_IN mode
call fglbegin(mode)
return
end subroutine f9xglbegin

subroutine f9xglbitmap(width,height,xorig,yorig,xmove,ymove,bitmap)
integer(kind=glsizei) INTENT_IN width, height
real(kind=glfloat) INTENT_IN xorig, yorig, xmove, ymove
! need assumed shape to determine size of array
integer(kind=glubyte), dimension(:) INTENT_IN bitmap
integer(kind=glint) :: xsize
#ifdef NOBYTE
#ifdef NOSHORT
   xsize = size(bitmap)
   call f9bglbitmap(width,height,xorig,yorig,xmove,ymove,bitmap,xsize)
#else
   xsize = size(bitmap)
   call f9cglbitmap(width,height,xorig,yorig,xmove,ymove,bitmap,xsize)
#endif
#else
   call fglbitmap(width,height,xorig,yorig,xmove,ymove,bitmap)
#endif
return
end subroutine f9xglbitmap

subroutine f9xglblendfunc(sfactor,dfactor)
integer(kind=glenum) INTENT_IN sfactor, dfactor
call fglblendfunc(sfactor,dfactor)
return
end subroutine f9xglblendfunc

subroutine f9xglcalllist(list)
integer(kind=gluint) INTENT_IN list
call fglcalllist(list)
return
end subroutine f9xglcalllist

subroutine f9x0glcalllists(n,ixtype,lists)
integer(kind=glsizei) INTENT_IN n
integer(kind=glenum) INTENT_IN ixtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN lists
call f9y0glcalllists(n,ixtype,lists)
return
end subroutine f9x0glcalllists

subroutine f9x1glcalllists(n,ixtype,lists)
integer(kind=glsizei) INTENT_IN n
integer(kind=glenum) INTENT_IN ixtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN lists
call f9y1glcalllists(n,ixtype,lists)
return
end subroutine f9x1glcalllists

subroutine f9x4glcalllists(n,ixtype,lists)
integer(kind=glsizei) INTENT_IN n
integer(kind=glenum) INTENT_IN ixtype
type(glcptr) INTENT_IN lists
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(lists%addr(i:i))
end do
i = glcptr_length
call f9y4glcalllists(n,ixtype,temp,i)
return
end subroutine f9x4glcalllists

subroutine f9x5glcalllists(n,ixtype,lists)
integer(kind=glsizei) INTENT_IN n
integer(kind=glenum) INTENT_IN ixtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN lists
call f9y5glcalllists(n,ixtype,lists)
return
end subroutine f9x5glcalllists

subroutine f9x6glcalllists(n,ixtype,lists)
integer(kind=glsizei) INTENT_IN n
integer(kind=glenum) INTENT_IN ixtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN lists
call f9y6glcalllists(n,ixtype,lists)
return
end subroutine f9x6glcalllists

subroutine f9xglclear(mask)
integer(kind=glbitfield) INTENT_IN mask
call fglclear(mask)
return
end subroutine f9xglclear

subroutine f9xglclearaccum(red,green,blue,alpha)
real(kind=glfloat) INTENT_IN red, green, blue, alpha
call fglclearaccum(red,green,blue,alpha)
return
end subroutine f9xglclearaccum

subroutine f9xglclearcolor(red,green,blue,alpha)
real(kind=glclampf) INTENT_IN red, green, blue, alpha
call fglclearcolor(red,green,blue,alpha)
return
end subroutine f9xglclearcolor

subroutine f9xglcleardepth(depth)
real(kind=glclampd) INTENT_IN depth
call fglcleardepth(depth)
return
end subroutine f9xglcleardepth

subroutine f9xglclearindex(c)
real(kind=glfloat) INTENT_IN c
call fglclearindex(c)
return
end subroutine f9xglclearindex

subroutine f9xglclearstencil(s)
integer(kind=glint) INTENT_IN s
call fglclearstencil(s)
return
end subroutine f9xglclearstencil

subroutine f9xglclipplane(plane,equation)
integer(kind=glenum) INTENT_IN plane
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN equation
call fglclipplane(plane,equation)
return
end subroutine f9xglclipplane

subroutine f9xglcolor3b(red,green,blue)
integer(kind=glbyte) INTENT_IN red, green, blue
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglcolor3b(red,green,blue)
#else
   call f9cglcolor3b(red,green,blue)
#endif
#else
   call fglcolor3b(red,green,blue)
#endif
return
end subroutine f9xglcolor3b

subroutine f9xglcolor3bv(v)
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglcolor3bv(v)
#else
   call f9cglcolor3bv(v)
#endif
#else
   call fglcolor3bv(v)
#endif
return
end subroutine f9xglcolor3bv

subroutine f9xglcolor3d(red,green,blue)
real(kind=gldouble) INTENT_IN red, green, blue
call fglcolor3d(red,green,blue)
return
end subroutine f9xglcolor3d

subroutine f9xglcolor3dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fglcolor3dv(v)
return
end subroutine f9xglcolor3dv

subroutine f9xglcolor3f(red,green,blue)
real(kind=glfloat) INTENT_IN red, green, blue
call fglcolor3f(red,green,blue)
return
end subroutine f9xglcolor3f

subroutine f9xglcolor3fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fglcolor3fv(v)
return
end subroutine f9xglcolor3fv

subroutine f9xglcolor3i(red,green,blue)
integer(kind=glint) INTENT_IN red, green, blue
call fglcolor3i(red,green,blue)
return
end subroutine f9xglcolor3i

subroutine f9xglcolor3iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fglcolor3iv(v)
return
end subroutine f9xglcolor3iv

subroutine f9xglcolor3s(red,green,blue)
integer(kind=glshort) INTENT_IN red, green, blue
#ifdef NOSHORT
   call f9aglcolor3s(red,green,blue)
#else
   call fglcolor3s(red,green,blue)
#endif
return
end subroutine f9xglcolor3s

subroutine f9xglcolor3sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9aglcolor3sv(v)
#else
   call fglcolor3sv(v)
#endif
return
end subroutine f9xglcolor3sv

subroutine f9xglcolor3ub(red,green,blue)
integer(kind=glubyte) INTENT_IN red, green, blue
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglcolor3ub(red,green,blue)
#else
   call f9cglcolor3ub(red,green,blue)
#endif
#else
   call fglcolor3ub(red,green,blue)
#endif
return
end subroutine f9xglcolor3ub

subroutine f9xglcolor3ubv(v)
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglcolor3ubv(v)
#else
   call f9cglcolor3ubv(v)
#endif
#else
   call fglcolor3ubv(v)
#endif
return
end subroutine f9xglcolor3ubv

subroutine f9xglcolor3ui(red,green,blue)
integer(kind=gluint) INTENT_IN red, green, blue
call fglcolor3ui(red,green,blue)
return
end subroutine f9xglcolor3ui

subroutine f9xglcolor3uiv(v)
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN v
call fglcolor3uiv(v)
return
end subroutine f9xglcolor3uiv

subroutine f9xglcolor3us(red,green,blue)
integer(kind=glushort) INTENT_IN red, green, blue
#ifdef NOSHORT
   call f9aglcolor3us(red,green,blue)
#else
   call fglcolor3us(red,green,blue)
#endif
return
end subroutine f9xglcolor3us

subroutine f9xglcolor3usv(v)
integer(kind=glushort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9aglcolor3usv(v)
#else
   call fglcolor3usv(v)
#endif
return
end subroutine f9xglcolor3usv

subroutine f9xglcolor4b(red,green,blue,alpha)
integer(kind=glbyte) INTENT_IN red, green, blue, alpha
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglcolor4b(red,green,blue,alpha)
#else
   call f9cglcolor4b(red,green,blue,alpha)
#endif
#else
   call fglcolor4b(red,green,blue,alpha)
#endif
return
end subroutine f9xglcolor4b

subroutine f9xglcolor4bv(v)
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglcolor4bv(v)
#else
   call f9cglcolor4bv(v)
#endif
#else
   call fglcolor4bv(v)
#endif
return
end subroutine f9xglcolor4bv

subroutine f9xglcolor4d(red,green,blue,alpha)
real(kind=gldouble) INTENT_IN red, green, blue, alpha
call fglcolor4d(red,green,blue,alpha)
return
end subroutine f9xglcolor4d

subroutine f9xglcolor4dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fglcolor4dv(v)
return
end subroutine f9xglcolor4dv

subroutine f9xglcolor4f(red,green,blue,alpha)
real(kind=glfloat) INTENT_IN red, green, blue, alpha
call fglcolor4f(red,green,blue,alpha)
return
end subroutine f9xglcolor4f

subroutine f9xglcolor4fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fglcolor4fv(v)
return
end subroutine f9xglcolor4fv

subroutine f9xglcolor4i(red,green,blue,alpha)
integer(kind=glint) INTENT_IN red, green, blue, alpha
call fglcolor4i(red,green,blue,alpha)
return
end subroutine f9xglcolor4i

subroutine f9xglcolor4iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fglcolor4iv(v)
return
end subroutine f9xglcolor4iv

subroutine f9xglcolor4s(red,green,blue,alpha)
integer(kind=glshort) INTENT_IN red, green, blue, alpha
#ifdef NOSHORT
   call f9aglcolor4s(red,green,blue,alpha)
#else
   call fglcolor4s(red,green,blue,alpha)
#endif
return
end subroutine f9xglcolor4s

subroutine f9xglcolor4sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9aglcolor4sv(v)
#else
   call fglcolor4sv(v)
#endif
return
end subroutine f9xglcolor4sv

subroutine f9xglcolor4ub(red,green,blue,alpha)
integer(kind=glubyte) INTENT_IN red, green, blue, alpha
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglcolor4ub(red,green,blue,alpha)
#else
   call f9cglcolor4ub(red,green,blue,alpha)
#endif
#else
   call fglcolor4ub(red,green,blue,alpha)
#endif
return
end subroutine f9xglcolor4ub

subroutine f9xglcolor4ubv(v)
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglcolor4ubv(v)
#else
   call f9cglcolor4ubv(v)
#endif
#else
   call fglcolor4ubv(v)
#endif
return
end subroutine f9xglcolor4ubv

subroutine f9xglcolor4ui(red,green,blue,alpha)
integer(kind=gluint) INTENT_IN red, green, blue, alpha
call fglcolor4ui(red,green,blue,alpha)
return
end subroutine f9xglcolor4ui

subroutine f9xglcolor4uiv(v)
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN v
call fglcolor4uiv(v)
return
end subroutine f9xglcolor4uiv

subroutine f9xglcolor4us(red,green,blue,alpha)
integer(kind=glushort) INTENT_IN red, green, blue, alpha
#ifdef NOSHORT
   call f9aglcolor4us(red,green,blue,alpha)
#else
   call fglcolor4us(red,green,blue,alpha)
#endif
return
end subroutine f9xglcolor4us

subroutine f9xglcolor4usv(v)
integer(kind=glushort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9aglcolor4usv(v)
#else
   call fglcolor4usv(v)
#endif
return
end subroutine f9xglcolor4usv

subroutine f9xglcolormask(red,green,blue,alpha)
logical(kind=glboolean) INTENT_IN red, green, blue, alpha
integer(kind=glint) :: int_red, int_green, int_blue, int_alpha
if (glone_byte_logical) then
   call fglcolormask(red,green,blue,alpha)
else
   if (red) then
      int_red = 1
   else
      int_red = 0
   end if
   if (green) then
      int_green = 1
   else
      int_green = 0
   end if
   if (blue) then
      int_blue = 1
   else
      int_blue = 0
   end if
   if (alpha) then
      int_alpha = 1
   else
      int_alpha = 0
   end if
   call f9eglcolormask(int_red,int_green,int_blue,int_alpha)
end if
return
end subroutine f9xglcolormask

subroutine f9xglcolormaterial(face,mode)
integer(kind=glenum) INTENT_IN face, mode
call fglcolormaterial(face,mode)
return
end subroutine f9xglcolormaterial

subroutine f9xglcopypixels(x,y,width,height,xtype)
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN xtype
call fglcopypixels(x,y,width,height,xtype)
return
end subroutine f9xglcopypixels

subroutine f9xglcullface(mode)
integer(kind=glenum) INTENT_IN mode
call fglcullface(mode)
return
end subroutine f9xglcullface

subroutine f9xgldeletelists(list,xrange)
integer(kind=gluint) INTENT_IN list
integer(kind=glsizei) INTENT_IN xrange
call fgldeletelists(list,xrange)
return
end subroutine f9xgldeletelists

subroutine f9xgldepthfunc(func)
integer(kind=glenum) INTENT_IN func
call fgldepthfunc(func)
return
end subroutine f9xgldepthfunc

subroutine f9xgldepthmask(flag)
logical(kind=glboolean) INTENT_IN flag
integer(kind=glint) :: int_flag
if (glone_byte_logical) then
   call fgldepthmask(flag)
else
   if (flag) then
      int_flag = 1
   else
      int_flag = 0
   end if
   call f9egldepthmask(int_flag)
end if
return
end subroutine f9xgldepthmask

subroutine f9xgldepthrange(near,far)
real(kind=glclampd) INTENT_IN near, far
call fgldepthrange(near,far)
return
end subroutine f9xgldepthrange

subroutine f9xgldisable(cap)
integer(kind=glenum) INTENT_IN cap
call fgldisable(cap)
return
end subroutine f9xgldisable

subroutine f9xgldrawbuffer(mode)
integer(kind=glenum) INTENT_IN mode
call fgldrawbuffer(mode)
return
end subroutine f9xgldrawbuffer

subroutine f9x0gldrawpixels(width,height,format,xtype,pixels)
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN pixels
call f9y0gldrawpixels(width,height,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x0gldrawpixels

subroutine f9x1gldrawpixels(width,height,format,xtype,pixels)
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
call f9y1gldrawpixels(width,height,format,xtype,pixels)
return
end subroutine f9x1gldrawpixels

subroutine f9x4gldrawpixels(width,height,format,xtype,pixels)
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i = glcptr_length
call f9y4gldrawpixels(width,height,format,xtype,temp,i)
return
end subroutine f9x4gldrawpixels

subroutine f9x5gldrawpixels(width,height,format,xtype,pixels)
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
call f9y5gldrawpixels(width,height,format,xtype,pixels)
return
end subroutine f9x5gldrawpixels

subroutine f9x6gldrawpixels(width,height,format,xtype,pixels)
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN pixels
call f9y6gldrawpixels(width,height,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x6gldrawpixels

subroutine f9xgledgeflag(flag)
logical(kind=glboolean) INTENT_IN flag
integer(kind=glint) :: int_flag
if (glone_byte_logical) then
   call fgledgeflag(flag)
else
   if (flag) then
      int_flag = 1
   else
      int_flag = 0
   end if
   call f9egledgeflag(int_flag)
end if
return
end subroutine f9xgledgeflag

subroutine f9xgledgeflagv(flag)
! need assumed shape to determine size of array
logical(kind=glboolean), dimension(:) INTENT_IN flag
integer(kind=glint), dimension(size(flag)) :: int_flag
integer(kind=glint) :: pos
if (glone_byte_logical) then
   call fgledgeflagv(flag)
else
   do pos=1,size(flag)
      if (flag(pos)) then
         int_flag(pos) = 1
      else
         int_flag(pos) = 0
      end if
   end do
   call f9egledgeflagv(int_flag,size_glint(int_flag))
end if
return
end subroutine f9xgledgeflagv

subroutine f9xglenable(cap)
integer(kind=glenum) INTENT_IN cap
call fglenable(cap)
return
end subroutine f9xglenable

subroutine f9xglend()
call fglend()
return
end subroutine f9xglend

subroutine f9xglendlist()
call fglendlist()
return
end subroutine f9xglendlist

subroutine f9xglevalcoord1d(u)
real(kind=gldouble) INTENT_IN u
call fglevalcoord1d(u)
return
end subroutine f9xglevalcoord1d

subroutine f9xglevalcoord1dv(u)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN u
call fglevalcoord1dv(u)
return
end subroutine f9xglevalcoord1dv

subroutine f9xglevalcoord1f(u)
real(kind=glfloat) INTENT_IN u
call fglevalcoord1f(u)
return
end subroutine f9xglevalcoord1f

subroutine f9xglevalcoord1fv(u)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN u
call fglevalcoord1fv(u)
return
end subroutine f9xglevalcoord1fv

subroutine f9xglevalcoord2d(u,v)
real(kind=gldouble) INTENT_IN u, v
call fglevalcoord2d(u,v)
return
end subroutine f9xglevalcoord2d

subroutine f9xglevalcoord2dv(u)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN u
call fglevalcoord2dv(u)
return
end subroutine f9xglevalcoord2dv

subroutine f9xglevalcoord2f(u,v)
real(kind=glfloat) INTENT_IN u, v
call fglevalcoord2f(u,v)
return
end subroutine f9xglevalcoord2f

subroutine f9xglevalcoord2fv(u)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN u
call fglevalcoord2fv(u)
return
end subroutine f9xglevalcoord2fv

subroutine f9xglevalmesh1(mode,i1,i2)
integer(kind=glenum) INTENT_IN mode
integer(kind=glint) INTENT_IN i1, i2
call fglevalmesh1(mode,i1,i2)
return
end subroutine f9xglevalmesh1

subroutine f9xglevalmesh2(mode,i1,i2,j1,j2)
integer(kind=glenum) INTENT_IN mode
integer(kind=glint) INTENT_IN i1, i2, j1, j2
call fglevalmesh2(mode,i1,i2,j1,j2)
return
end subroutine f9xglevalmesh2

subroutine f9xglevalpoint1(i)
integer(kind=glint) INTENT_IN i
call fglevalpoint1(i)
return
end subroutine f9xglevalpoint1

subroutine f9xglevalpoint2(i,j)
integer(kind=glint) INTENT_IN i, j
call fglevalpoint2(i,j)
return
end subroutine f9xglevalpoint2

subroutine f9xglfeedbackbuffer(xsize,xtype,buffer)
integer(kind=glsizei) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: buffer
call fglfeedbackbuffer(xsize,xtype,buffer)
return
end subroutine f9xglfeedbackbuffer

subroutine f9xglfinish()
call fglfinish()
return
end subroutine f9xglfinish

subroutine f9xglflush()
call fglflush()
return
end subroutine f9xglflush

subroutine f9xglfogf(pname,param)
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat) INTENT_IN param
call fglfogf(pname,param)
return
end subroutine f9xglfogf

subroutine f9xglfogfv(pname,params)
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
call fglfogfv(pname,params)
return
end subroutine f9xglfogfv

subroutine f9xglfogi(pname,param)
integer(kind=glenum) INTENT_IN pname
integer(kind=glint) INTENT_IN param
call fglfogi(pname,param)
return
end subroutine f9xglfogi

subroutine f9xglfogiv(pname,params)
integer(kind=glenum) INTENT_IN pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
call fglfogiv(pname,params)
return
end subroutine f9xglfogiv

subroutine f9xglfrontface(mode)
integer(kind=glenum) INTENT_IN mode
call fglfrontface(mode)
return
end subroutine f9xglfrontface

subroutine f9xglfrustum(left,right,bottom,top,near,far)
real(kind=gldouble) INTENT_IN left, right, bottom, top, near, far
call fglfrustum(left,right,bottom,top,near,far)
return
end subroutine f9xglfrustum

function f9xglgenlists(xrange) result(resf9xglgenlists)
integer(kind=glsizei) INTENT_IN xrange
integer(kind=gluint) :: resf9xglgenlists
resf9xglgenlists = fglgenlists(xrange)
return
end function f9xglgenlists

subroutine f9xglgetbooleanv(pname,params)
integer(kind=glenum) INTENT_IN pname
! need assumed shape to determine size of array
logical(kind=glboolean), dimension(:), intent(out) :: params
integer(kind=glint), dimension(size(params)) :: int_params
integer(kind=glint) :: i
if (glone_byte_logical) then
   call fglgetbooleanv(pname,params)
else
   call f9eglgetbooleanv(pname,int_params,size_glint(int_params))
   do i=1,size(params)
      params(i) = (int_params(i) /= 0)
   end do
end if
return
end subroutine f9xglgetbooleanv

subroutine f9xglgetclipplane(plane,equation)
integer(kind=glenum) INTENT_IN plane
real(kind=gldouble), dimension(DEFER_DIM), intent(out) :: equation
call fglgetclipplane(plane,equation)
return
end subroutine f9xglgetclipplane

subroutine f9xglgetdoublev(pname,params)
integer(kind=glenum) INTENT_IN pname
real(kind=gldouble), dimension(DEFER_DIM), intent(out) :: params
call fglgetdoublev(pname,params)
return
end subroutine f9xglgetdoublev

function f9xglgeterror() result(resf9xglgeterror)
integer(kind=glenum) :: resf9xglgeterror
resf9xglgeterror = fglgeterror()
return
end function f9xglgeterror

subroutine f9xglgetfloatv(pname,params)
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
call fglgetfloatv(pname,params)
return
end subroutine f9xglgetfloatv

subroutine f9xglgetintegerv(pname,params)
integer(kind=glenum) INTENT_IN pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
call fglgetintegerv(pname,params)
return
end subroutine f9xglgetintegerv

subroutine f9xglgetlightfv(light,pname,params)
integer(kind=glenum) INTENT_IN light, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
call fglgetlightfv(light,pname,params)
return
end subroutine f9xglgetlightfv

subroutine f9xglgetlightiv(light,pname,params)
integer(kind=glenum) INTENT_IN light, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
call fglgetlightiv(light,pname,params)
return
end subroutine f9xglgetlightiv

subroutine f9xglgetmapdv(xtarget,query,v)
integer(kind=glenum) INTENT_IN xtarget, query
real(kind=gldouble), dimension(DEFER_DIM), intent(out) :: v
call fglgetmapdv(xtarget,query,v)
return
end subroutine f9xglgetmapdv

subroutine f9xglgetmapfv(xtarget,query,v)
integer(kind=glenum) INTENT_IN xtarget, query
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: v
call fglgetmapfv(xtarget,query,v)
return
end subroutine f9xglgetmapfv

subroutine f9xglgetmapiv(xtarget,query,v)
integer(kind=glenum) INTENT_IN xtarget, query
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: v
call fglgetmapiv(xtarget,query,v)
return
end subroutine f9xglgetmapiv

subroutine f9xglgetmaterialfv(face,pname,params)
integer(kind=glenum) INTENT_IN face, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
call fglgetmaterialfv(face,pname,params)
return
end subroutine f9xglgetmaterialfv

subroutine f9xglgetmaterialiv(face,pname,params)
integer(kind=glenum) INTENT_IN face, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
call fglgetmaterialiv(face,pname,params)
return
end subroutine f9xglgetmaterialiv

subroutine f9xglgetpixelmapfv(map,values)
integer(kind=glenum) INTENT_IN map
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: values
call fglgetpixelmapfv(map,values)
return
end subroutine f9xglgetpixelmapfv

subroutine f9xglgetpixelmapuiv(map,values)
integer(kind=glenum) INTENT_IN map
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: values
call fglgetpixelmapuiv(map,values)
return
end subroutine f9xglgetpixelmapuiv

subroutine f9xglgetpixelmapusv(map,values)
integer(kind=glenum) INTENT_IN map
integer(kind=glushort), dimension(DEFER_DIM), intent(out) :: values
#ifdef NOSHORT
   call f9aglgetpixelmapusv(map,values)
#else
   call fglgetpixelmapusv(map,values)
#endif
return
end subroutine f9xglgetpixelmapusv

subroutine f9xglgetpolygonstipple(mask)
! need assumed shape to determine size of array
integer(kind=glubyte), dimension(:), intent(out) :: mask
integer(kind=glint) :: xsize
xsize = size(mask)
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglgetpolygonstipple(mask,xsize)
#else
   call f9cglgetpolygonstipple(mask,xsize)
#endif
#else
   call fglgetpolygonstipple(mask,xsize)
#endif
return
end subroutine f9xglgetpolygonstipple

function f9xglgetstring(name) result(resf9xglgetstring)
integer(kind=glenum) INTENT_IN name
character(len=1), pointer, dimension(:) :: resf9xglgetstring
integer(kind=gluint), allocatable, dimension(:) :: temp
integer(kind=glint) :: length, pos
#ifdef PURE_FUNCTIONS
integer(kind=glint) :: dummy
#endif
length = f9yglgetstringlen(name)
allocate(temp(length))
allocate(resf9xglgetstring(length))
#ifdef PURE_FUNCTIONS
do pos=1,length
   temp(pos) = f9yglgetstring1(pos)
end do
#else
call f9yglgetstring(temp)
#endif
do pos=1,length
   resf9xglgetstring(pos) = char(temp(pos))
end do
#ifdef PURE_FUNCTIONS
dummy = f9yglgetstringfree()
#else
call f9yglgetstringfree()
#endif
deallocate(temp)
return
end function f9xglgetstring

subroutine f9xglgettexenvfv(xtarget,pname,params)
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
call fglgettexenvfv(xtarget,pname,params)
return
end subroutine f9xglgettexenvfv

subroutine f9xglgettexenviv(xtarget,pname,params)
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
call fglgettexenviv(xtarget,pname,params)
return
end subroutine f9xglgettexenviv

subroutine f9xglgettexgendv(coord,pname,params)
integer(kind=glenum) INTENT_IN coord, pname
real(kind=gldouble), dimension(DEFER_DIM), intent(out) :: params
call fglgettexgendv(coord,pname,params)
return
end subroutine f9xglgettexgendv

subroutine f9xglgettexgenfv(coord,pname,params)
integer(kind=glenum) INTENT_IN coord, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
call fglgettexgenfv(coord,pname,params)
return
end subroutine f9xglgettexgenfv

subroutine f9xglgettexgeniv(coord,pname,params)
integer(kind=glenum) INTENT_IN coord, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
call fglgettexgeniv(coord,pname,params)
return
end subroutine f9xglgettexgeniv

subroutine f9x0glgetteximage(xtarget,level,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:), intent(out) :: pixels
call f9y0glgetteximage(xtarget,level,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x0glgetteximage

subroutine f9x1glgetteximage(xtarget,level,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: pixels
call f9y1glgetteximage(xtarget,level,format,xtype,pixels)
return
end subroutine f9x1glgetteximage

subroutine f9x4glgetteximage(xtarget,level,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr), intent(in) :: pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i = glcptr_length
call f9y4glgetteximage(xtarget,level,format,xtype,temp,i)
return
end subroutine f9x4glgetteximage

subroutine f9x5glgetteximage(xtarget,level,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: pixels
call f9y5glgetteximage(xtarget,level,format,xtype,pixels)
return
end subroutine f9x5glgetteximage

subroutine f9x6glgetteximage(xtarget,level,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:), intent(out) :: pixels
call f9y6glgetteximage(xtarget,level,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x6glgetteximage

subroutine f9xglgettexlevelparameterfv(xtarget,level,pname,params)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
call fglgettexlevelparameterfv(xtarget,level,pname,params)
return
end subroutine f9xglgettexlevelparameterfv

subroutine f9xglgettexlevelparameteriv(xtarget,level,pname,params)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
call fglgettexlevelparameteriv(xtarget,level,pname,params)
return
end subroutine f9xglgettexlevelparameteriv

subroutine f9xglgettexparameterfv(xtarget,pname,params)
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
call fglgettexparameterfv(xtarget,pname,params)
return
end subroutine f9xglgettexparameterfv

subroutine f9xglgettexparameteriv(xtarget,pname,params)
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
call fglgettexparameteriv(xtarget,pname,params)
return
end subroutine f9xglgettexparameteriv

subroutine f9xglhint(xtarget,mode)
integer(kind=glenum) INTENT_IN xtarget, mode
call fglhint(xtarget,mode)
return
end subroutine f9xglhint

subroutine f9xglindexmask(mask)
integer(kind=gluint) INTENT_IN mask
call fglindexmask(mask)
return
end subroutine f9xglindexmask

subroutine f9xglindexd(c)
real(kind=gldouble) INTENT_IN c
call fglindexd(c)
return
end subroutine f9xglindexd

subroutine f9xglindexdv(c)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN c
call fglindexdv(c)
return
end subroutine f9xglindexdv

subroutine f9xglindexf(c)
real(kind=glfloat) INTENT_IN c
call fglindexf(c)
return
end subroutine f9xglindexf

subroutine f9xglindexfv(c)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN c
call fglindexfv(c)
return
end subroutine f9xglindexfv

subroutine f9xglindexi(c)
integer(kind=glint) INTENT_IN c
call fglindexi(c)
return
end subroutine f9xglindexi

subroutine f9xglindexiv(c)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN c
call fglindexiv(c)
return
end subroutine f9xglindexiv

subroutine f9xglindexs(c)
integer(kind=glshort) INTENT_IN c
#ifdef NOSHORT
   call f9aglindexs(c)
#else
   call fglindexs(c)
#endif
return
end subroutine f9xglindexs

subroutine f9xglindexsv(c)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN c
#ifdef NOSHORT
   call f9aglindexsv(c)
#else
   call fglindexsv(c)
#endif
return
end subroutine f9xglindexsv

subroutine f9xglinitnames()
call fglinitnames()
return
end subroutine f9xglinitnames

function f9xglisenabled(cap) result(resf9xglisenabled)
integer(kind=glenum) INTENT_IN cap
logical(kind=glboolean) :: resf9xglisenabled
integer(kind=glint) :: temp
if (glone_byte_logical) then
   resf9xglisenabled = fglisenabled(cap)
else
   temp = f9eglisenabled(cap)
   if (temp==0) then
      resf9xglisenabled = .false.
   else
      resf9xglisenabled = .true.
   end if
end if
return
end function f9xglisenabled

function f9xglislist(list) result(resf9xglislist)
integer(kind=gluint) INTENT_IN list
logical(kind=glboolean) :: resf9xglislist
integer(kind=glint) :: temp
if (glone_byte_logical) then
   resf9xglislist = fglislist(list)
else
   temp = f9eglislist(list)
   if (temp==0) then
      resf9xglislist = .false.
   else
      resf9xglislist = .true.
   end if
end if
return
end function f9xglislist

subroutine f9xgllightmodelf(pname,param)
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat) INTENT_IN param
call fgllightmodelf(pname,param)
return
end subroutine f9xgllightmodelf

subroutine f9xgllightmodelfv(pname,params)
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
call fgllightmodelfv(pname,params)
return
end subroutine f9xgllightmodelfv

subroutine f9xgllightmodeli(pname,param)
integer(kind=glenum) INTENT_IN pname
integer(kind=glint) INTENT_IN param
call fgllightmodeli(pname,param)
return
end subroutine f9xgllightmodeli

subroutine f9xgllightmodeliv(pname,params)
integer(kind=glenum) INTENT_IN pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
call fgllightmodeliv(pname,params)
return
end subroutine f9xgllightmodeliv

subroutine f9xgllightf(light,pname,param)
integer(kind=glenum) INTENT_IN light, pname
real(kind=glfloat) INTENT_IN param
call fgllightf(light,pname,param)
return
end subroutine f9xgllightf

subroutine f9xgllightfv(light,pname,params)
integer(kind=glenum) INTENT_IN light, pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
call fgllightfv(light,pname,params)
return
end subroutine f9xgllightfv

subroutine f9xgllighti(light,pname,param)
integer(kind=glenum) INTENT_IN light, pname
integer(kind=glint) INTENT_IN param
call fgllighti(light,pname,param)
return
end subroutine f9xgllighti

subroutine f9xgllightiv(light,pname,params)
integer(kind=glenum) INTENT_IN light, pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
call fgllightiv(light,pname,params)
return
end subroutine f9xgllightiv

subroutine f9xgllinestipple(factor,pattern)
integer(kind=glint) INTENT_IN factor
integer(kind=glushort) INTENT_IN pattern
#ifdef NOSHORT
   call f9agllinestipple(factor,pattern)
#else
   call fgllinestipple(factor,pattern)
#endif
return
end subroutine f9xgllinestipple

subroutine f9xgllinewidth(width)
real(kind=glfloat) INTENT_IN width
call fgllinewidth(width)
return
end subroutine f9xgllinewidth

subroutine f9xgllistbase(base)
integer(kind=gluint) INTENT_IN base
call fgllistbase(base)
return
end subroutine f9xgllistbase

subroutine f9xglloadidentity()
call fglloadidentity()
return
end subroutine f9xglloadidentity

subroutine f9xglloadmatrixd(m)
#ifdef F_LANG
real(kind=gldouble), dimension(:,:) INTENT_IN m
#else
real(kind=gldouble), dimension(4,DEFER_DIM) INTENT_IN m
#endif
call fglloadmatrixd(m)
return
end subroutine f9xglloadmatrixd

subroutine f9xglloadmatrixf(m)
#ifdef F_LANG
real(kind=glfloat), dimension(:,:) INTENT_IN m
#else
real(kind=glfloat), dimension(4,DEFER_DIM) INTENT_IN m
#endif
call fglloadmatrixf(m)
return
end subroutine f9xglloadmatrixf

subroutine f9xglloadname(name)
integer(kind=gluint) INTENT_IN name
call fglloadname(name)
return
end subroutine f9xglloadname

subroutine f9xgllogicop(opcode)
integer(kind=glenum) INTENT_IN opcode
call fgllogicop(opcode)
return
end subroutine f9xgllogicop

subroutine f9xglmap1d(xtarget,u1,u2,stride,order,points)
integer(kind=glenum) INTENT_IN xtarget
real(kind=gldouble) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN stride, order
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN points
call fglmap1d(xtarget,u1,u2,stride,order,points)
return
end subroutine f9xglmap1d

subroutine f9xglmap1f(xtarget,u1,u2,stride,order,points)
integer(kind=glenum) INTENT_IN xtarget
real(kind=glfloat) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN stride, order
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN points
call fglmap1f(xtarget,u1,u2,stride,order,points)
return
end subroutine f9xglmap1f

subroutine f9xglmap2d(xtarget,u1,u2,ustride,uorder,v1,v2,vstride,vorder, &
   points)
integer(kind=glenum) INTENT_IN xtarget
real(kind=gldouble) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN ustride, uorder
real(kind=gldouble) INTENT_IN v1, v2
integer(kind=glint) INTENT_IN vstride, vorder
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN points
call fglmap2d(xtarget,u1,u2,ustride,uorder,v1,v2,vstride,vorder,points)
return
end subroutine f9xglmap2d

subroutine f9xglmap2f(xtarget,u1,u2,ustride,uorder,v1,v2,vstride,vorder, &
   points)
integer(kind=glenum) INTENT_IN xtarget
real(kind=glfloat) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN ustride, uorder
real(kind=glfloat) INTENT_IN v1, v2
integer(kind=glint) INTENT_IN vstride, vorder
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN points
call fglmap2f(xtarget,u1,u2,ustride,uorder,v1,v2,vstride,vorder,points)
return
end subroutine f9xglmap2f

subroutine f9xglmapgrid1d(un,u1,u2)
integer(kind=glint) INTENT_IN un
real(kind=gldouble) INTENT_IN u1, u2
call fglmapgrid1d(un,u1,u2)
return
end subroutine f9xglmapgrid1d

subroutine f9xglmapgrid1f(un,u1,u2)
integer(kind=glint) INTENT_IN un
real(kind=glfloat) INTENT_IN u1, u2
call fglmapgrid1f(un,u1,u2)
return
end subroutine f9xglmapgrid1f

subroutine f9xglmapgrid2d(un,u1,u2,vn,v1,v2)
integer(kind=glint) INTENT_IN un
real(kind=gldouble) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN vn
real(kind=gldouble) INTENT_IN v1, v2
call fglmapgrid2d(un,u1,u2,vn,v1,v2)
return
end subroutine f9xglmapgrid2d

subroutine f9xglmapgrid2f(un,u1,u2,vn,v1,v2)
integer(kind=glint) INTENT_IN un
real(kind=glfloat) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN vn
real(kind=glfloat) INTENT_IN v1, v2
call fglmapgrid2f(un,u1,u2,vn,v1,v2)
return
end subroutine f9xglmapgrid2f

subroutine f9xglmaterialf(face,pname,param)
integer(kind=glenum) INTENT_IN face, pname
real(kind=glfloat) INTENT_IN param
call fglmaterialf(face,pname,param)
return
end subroutine f9xglmaterialf

subroutine f9xglmaterialfv(face,pname,params)
integer(kind=glenum) INTENT_IN face, pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
call fglmaterialfv(face,pname,params)
return
end subroutine f9xglmaterialfv

subroutine f9xglmateriali(face,pname,param)
integer(kind=glenum) INTENT_IN face, pname
integer(kind=glint) INTENT_IN param
call fglmateriali(face,pname,param)
return
end subroutine f9xglmateriali

subroutine f9xglmaterialiv(face,pname,params)
integer(kind=glenum) INTENT_IN face, pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
call fglmaterialiv(face,pname,params)
return
end subroutine f9xglmaterialiv

subroutine f9xglmatrixmode(mode)
integer(kind=glenum) INTENT_IN mode
call fglmatrixmode(mode)
return
end subroutine f9xglmatrixmode

subroutine f9xglmultmatrixd(m)
#ifdef F_LANG
real(kind=gldouble), dimension(:,:) INTENT_IN m
#else
real(kind=gldouble), dimension(4,DEFER_DIM) INTENT_IN m
#endif
call fglmultmatrixd(m)
return
end subroutine f9xglmultmatrixd

subroutine f9xglmultmatrixf(m)
#ifdef F_LANG
real(kind=glfloat), dimension(:,:) INTENT_IN m
#else
real(kind=glfloat), dimension(4,DEFER_DIM) INTENT_IN m
#endif
call fglmultmatrixf(m)
return
end subroutine f9xglmultmatrixf

subroutine f9xglnewlist(list,mode)
integer(kind=gluint) INTENT_IN list
integer(kind=glenum) INTENT_IN mode
call fglnewlist(list,mode)
return
end subroutine f9xglnewlist

subroutine f9xglnormal3b(nx,ny,nz)
integer(kind=glbyte) INTENT_IN nx, ny, nz
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglnormal3b(nx,ny,nz)
#else
   call f9cglnormal3b(nx,ny,nz)
#endif
#else
   call fglnormal3b(nx,ny,nz)
#endif
return
end subroutine f9xglnormal3b

subroutine f9xglnormal3bv(v)
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglnormal3bv(v)
#else
   call f9cglnormal3bv(v)
#endif
#else
   call fglnormal3bv(v)
#endif
return
end subroutine f9xglnormal3bv

subroutine f9xglnormal3d(nx,ny,nz)
real(kind=gldouble) INTENT_IN nx, ny, nz
call fglnormal3d(nx,ny,nz)
return
end subroutine f9xglnormal3d

subroutine f9xglnormal3dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fglnormal3dv(v)
return
end subroutine f9xglnormal3dv

subroutine f9xglnormal3f(nx,ny,nz)
real(kind=glfloat) INTENT_IN nx, ny, nz
call fglnormal3f(nx,ny,nz)
return
end subroutine f9xglnormal3f

subroutine f9xglnormal3fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fglnormal3fv(v)
return
end subroutine f9xglnormal3fv

subroutine f9xglnormal3i(nx,ny,nz)
integer(kind=glint) INTENT_IN nx, ny, nz
call fglnormal3i(nx,ny,nz)
return
end subroutine f9xglnormal3i

subroutine f9xglnormal3iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fglnormal3iv(v)
return
end subroutine f9xglnormal3iv

subroutine f9xglnormal3s(nx,ny,nz)
integer(kind=glshort) INTENT_IN nx, ny, nz
#ifdef NOSHORT
   call f9aglnormal3s(nx,ny,nz)
#else
   call fglnormal3s(nx,ny,nz)
#endif
return
end subroutine f9xglnormal3s

subroutine f9xglnormal3sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9aglnormal3sv(v)
#else
   call fglnormal3sv(v)
#endif
return
end subroutine f9xglnormal3sv

subroutine f9xglortho(left,right,bottom,top,near,far)
real(kind=gldouble) INTENT_IN left, right, bottom, top, near, far
call fglortho(left,right,bottom,top,near,far)
return
end subroutine f9xglortho

subroutine f9xglpassthrough(token)
real(kind=glfloat) INTENT_IN token
call fglpassthrough(token)
return
end subroutine f9xglpassthrough

subroutine f9xglpixelmapfv(map,mapsize,values)
integer(kind=glenum) INTENT_IN map
integer(kind=glint) INTENT_IN mapsize
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN values
call fglpixelmapfv(map,mapsize,values)
return
end subroutine f9xglpixelmapfv

subroutine f9xglpixelmapuiv(map,mapsize,values)
integer(kind=glenum) INTENT_IN map
integer(kind=glint) INTENT_IN mapsize
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN values
call fglpixelmapuiv(map,mapsize,values)
return
end subroutine f9xglpixelmapuiv

subroutine f9xglpixelmapusv(map,mapsize,values)
integer(kind=glenum) INTENT_IN map
integer(kind=glint) INTENT_IN mapsize
integer(kind=glushort), dimension(DEFER_DIM) INTENT_IN values
#ifdef NOSHORT
   call f9aglpixelmapusv(map,mapsize,values)
#else
   call fglpixelmapusv(map,mapsize,values)
#endif
return
end subroutine f9xglpixelmapusv

subroutine f9xglpixelstoref(pname,param)
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat) INTENT_IN param
call fglpixelstoref(pname,param)
return
end subroutine f9xglpixelstoref

subroutine f9xglpixelstorei(pname,param)
integer(kind=glenum) INTENT_IN pname
integer(kind=glint) INTENT_IN param
call fglpixelstorei(pname,param)
return
end subroutine f9xglpixelstorei

subroutine f9xglpixeltransferf(pname,param)
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat) INTENT_IN param
call fglpixeltransferf(pname,param)
return
end subroutine f9xglpixeltransferf

subroutine f9xglpixeltransferi(pname,param)
integer(kind=glenum) INTENT_IN pname
integer(kind=glint) INTENT_IN param
call fglpixeltransferi(pname,param)
return
end subroutine f9xglpixeltransferi

subroutine f9xglpixelzoom(xfactor,yfactor)
real(kind=glfloat) INTENT_IN xfactor, yfactor
call fglpixelzoom(xfactor,yfactor)
return
end subroutine f9xglpixelzoom

subroutine f9xglpointsize(xsize)
real(kind=glfloat) INTENT_IN xsize
call fglpointsize(xsize)
return
end subroutine f9xglpointsize

subroutine f9xglpolygonmode(face,mode)
integer(kind=glenum) INTENT_IN face, mode
call fglpolygonmode(face,mode)
return
end subroutine f9xglpolygonmode

subroutine f9xglpolygonstipple(mask)
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN mask
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglpolygonstipple(mask)
#else
   call f9cglpolygonstipple(mask)
#endif
#else
   call fglpolygonstipple(mask)
#endif
return
end subroutine f9xglpolygonstipple

subroutine f9xglpopattrib()
call fglpopattrib()
return
end subroutine f9xglpopattrib

subroutine f9xglpopmatrix()
call fglpopmatrix()
return
end subroutine f9xglpopmatrix

subroutine f9xglpopname()
call fglpopname()
return
end subroutine f9xglpopname

subroutine f9xglpushattrib(mask)
integer(kind=glbitfield) INTENT_IN mask
call fglpushattrib(mask)
return
end subroutine f9xglpushattrib

subroutine f9xglpushmatrix()
call fglpushmatrix()
return
end subroutine f9xglpushmatrix

subroutine f9xglpushname(name)
integer(kind=gluint) INTENT_IN name
call fglpushname(name)
return
end subroutine f9xglpushname

subroutine f9xglrasterpos2d(x,y)
real(kind=gldouble) INTENT_IN x, y
call fglrasterpos2d(x,y)
return
end subroutine f9xglrasterpos2d

subroutine f9xglrasterpos2dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fglrasterpos2dv(v)
return
end subroutine f9xglrasterpos2dv

subroutine f9xglrasterpos2f(x,y)
real(kind=glfloat) INTENT_IN x, y
call fglrasterpos2f(x,y)
return
end subroutine f9xglrasterpos2f

subroutine f9xglrasterpos2fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fglrasterpos2fv(v)
return
end subroutine f9xglrasterpos2fv

subroutine f9xglrasterpos2i(x,y)
integer(kind=glint) INTENT_IN x, y
call fglrasterpos2i(x,y)
return
end subroutine f9xglrasterpos2i

subroutine f9xglrasterpos2iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fglrasterpos2iv(v)
return
end subroutine f9xglrasterpos2iv

subroutine f9xglrasterpos2s(x,y)
integer(kind=glshort) INTENT_IN x, y
#ifdef NOSHORT
   call f9aglrasterpos2s(x,y)
#else
   call fglrasterpos2s(x,y)
#endif
return
end subroutine f9xglrasterpos2s

subroutine f9xglrasterpos2sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9aglrasterpos2sv(v)
#else
   call fglrasterpos2sv(v)
#endif
return
end subroutine f9xglrasterpos2sv

subroutine f9xglrasterpos3d(x,y,z)
real(kind=gldouble) INTENT_IN x, y, z
call fglrasterpos3d(x,y,z)
return
end subroutine f9xglrasterpos3d

subroutine f9xglrasterpos3dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fglrasterpos3dv(v)
return
end subroutine f9xglrasterpos3dv

subroutine f9xglrasterpos3f(x,y,z)
real(kind=glfloat) INTENT_IN x, y, z
call fglrasterpos3f(x,y,z)
return
end subroutine f9xglrasterpos3f

subroutine f9xglrasterpos3fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fglrasterpos3fv(v)
return
end subroutine f9xglrasterpos3fv

subroutine f9xglrasterpos3i(x,y,z)
integer(kind=glint) INTENT_IN x, y, z
call fglrasterpos3i(x,y,z)
return
end subroutine f9xglrasterpos3i

subroutine f9xglrasterpos3iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fglrasterpos3iv(v)
return
end subroutine f9xglrasterpos3iv

subroutine f9xglrasterpos3s(x,y,z)
integer(kind=glshort) INTENT_IN x, y, z
#ifdef NOSHORT
   call f9aglrasterpos3s(x,y,z)
#else
   call fglrasterpos3s(x,y,z)
#endif
return
end subroutine f9xglrasterpos3s

subroutine f9xglrasterpos3sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9aglrasterpos3sv(v)
#else
   call fglrasterpos3sv(v)
#endif
return
end subroutine f9xglrasterpos3sv

subroutine f9xglrasterpos4d(x,y,z,w)
real(kind=gldouble) INTENT_IN x, y, z, w
call fglrasterpos4d(x,y,z,w)
return
end subroutine f9xglrasterpos4d

subroutine f9xglrasterpos4dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fglrasterpos4dv(v)
return
end subroutine f9xglrasterpos4dv

subroutine f9xglrasterpos4f(x,y,z,w)
real(kind=glfloat) INTENT_IN x, y, z, w
call fglrasterpos4f(x,y,z,w)
return
end subroutine f9xglrasterpos4f

subroutine f9xglrasterpos4fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fglrasterpos4fv(v)
return
end subroutine f9xglrasterpos4fv

subroutine f9xglrasterpos4i(x,y,z,w)
integer(kind=glint) INTENT_IN x, y, z, w
call fglrasterpos4i(x,y,z,w)
return
end subroutine f9xglrasterpos4i

subroutine f9xglrasterpos4iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fglrasterpos4iv(v)
return
end subroutine f9xglrasterpos4iv

subroutine f9xglrasterpos4s(x,y,z,w)
integer(kind=glshort) INTENT_IN x, y, z, w
#ifdef NOSHORT
   call f9aglrasterpos4s(x,y,z,w)
#else
   call fglrasterpos4s(x,y,z,w)
#endif
return
end subroutine f9xglrasterpos4s

subroutine f9xglrasterpos4sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9aglrasterpos4sv(v)
#else
   call fglrasterpos4sv(v)
#endif
return
end subroutine f9xglrasterpos4sv

subroutine f9xglreadbuffer(mode)
integer(kind=glenum) INTENT_IN mode
call fglreadbuffer(mode)
return
end subroutine f9xglreadbuffer

subroutine f9x0glreadpixels(x,y,width,height,format,xtype,pixels)
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:), intent(out) :: pixels
call f9y0glreadpixels(x,y,width,height,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x0glreadpixels

subroutine f9x1glreadpixels(x,y,width,height,format,xtype,pixels)
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: pixels
call f9y1glreadpixels(x,y,width,height,format,xtype,pixels)
return
end subroutine f9x1glreadpixels

subroutine f9x4glreadpixels(x,y,width,height,format,xtype,pixels)
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr), intent(in) :: pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i = glcptr_length
call f9y4glreadpixels(x,y,width,height,format,xtype,temp,i)
return
end subroutine f9x4glreadpixels

subroutine f9x5glreadpixels(x,y,width,height,format,xtype,pixels)
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: pixels
call f9y5glreadpixels(x,y,width,height,format,xtype,pixels)
return
end subroutine f9x5glreadpixels

subroutine f9x6glreadpixels(x,y,width,height,format,xtype,pixels)
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:), intent(out) :: pixels
call f9y6glreadpixels(x,y,width,height,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x6glreadpixels

subroutine f9xglrectd(x1,y1,x2,y2)
real(kind=gldouble) INTENT_IN x1, y1, x2, y2
call fglrectd(x1,y1,x2,y2)
return
end subroutine f9xglrectd

subroutine f9xglrectdv(v1,v2)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v1, v2
call fglrectdv(v1,v2)
return
end subroutine f9xglrectdv

subroutine f9xglrectf(x1,y1,x2,y2)
real(kind=glfloat) INTENT_IN x1, y1, x2, y2
call fglrectf(x1,y1,x2,y2)
return
end subroutine f9xglrectf

subroutine f9xglrectfv(v1,v2)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v1, v2
call fglrectfv(v1,v2)
return
end subroutine f9xglrectfv

subroutine f9xglrecti(x1,y1,x2,y2)
integer(kind=glint) INTENT_IN x1, y1, x2, y2
call fglrecti(x1,y1,x2,y2)
return
end subroutine f9xglrecti

subroutine f9xglrectiv(v1,v2)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v1, v2
call fglrectiv(v1,v2)
return
end subroutine f9xglrectiv

subroutine f9xglrects(x1,y1,x2,y2)
integer(kind=glshort) INTENT_IN x1, y1, x2, y2
#ifdef NOSHORT
   call f9aglrects(x1,y1,x2,y2)
#else
   call fglrects(x1,y1,x2,y2)
#endif
return
end subroutine f9xglrects

subroutine f9xglrectsv(v1,v2)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v1, v2
#ifdef NOSHORT
   call f9aglrectsv(v1,v2)
#else
   call fglrectsv(v1,v2)
#endif
return
end subroutine f9xglrectsv

function f9xglrendermode(mode) result(resf9xglrendermode)
integer(kind=glenum) INTENT_IN mode
integer(kind=glint) :: resf9xglrendermode
resf9xglrendermode = fglrendermode(mode)
return
end function f9xglrendermode

subroutine f9xglrotated(angle,x,y,z)
real(kind=gldouble) INTENT_IN angle, x, y, z
call fglrotated(angle,x,y,z)
return
end subroutine f9xglrotated

subroutine f9xglrotatef(angle,x,y,z)
real(kind=glfloat) INTENT_IN angle, x, y, z
call fglrotatef(angle,x,y,z)
return
end subroutine f9xglrotatef

subroutine f9xglscaled(x,y,z)
real(kind=gldouble) INTENT_IN x, y, z
call fglscaled(x,y,z)
return
end subroutine f9xglscaled

subroutine f9xglscalef(x,y,z)
real(kind=glfloat) INTENT_IN x, y, z
call fglscalef(x,y,z)
return
end subroutine f9xglscalef

subroutine f9xglscissor(x,y,width,height)
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
call fglscissor(x,y,width,height)
return
end subroutine f9xglscissor

subroutine f9xglselectbuffer(xsize,buffer)
integer(kind=glsizei) INTENT_IN xsize
integer(kind=gluint), dimension(DEFER_DIM), target, intent(in out) :: buffer
call fglselectbuffer(xsize,buffer)
return
end subroutine f9xglselectbuffer

subroutine f9xglshademodel(mode)
integer(kind=glenum) INTENT_IN mode
call fglshademodel(mode)
return
end subroutine f9xglshademodel

subroutine f9xglstencilfunc(func,ref,mask)
integer(kind=glenum) INTENT_IN func
integer(kind=glint) INTENT_IN ref
integer(kind=gluint) INTENT_IN mask
call fglstencilfunc(func,ref,mask)
return
end subroutine f9xglstencilfunc

subroutine f9xglstencilmask(mask)
integer(kind=gluint) INTENT_IN mask
call fglstencilmask(mask)
return
end subroutine f9xglstencilmask

subroutine f9xglstencilop(fail,zfail,zpass)
integer(kind=glenum) INTENT_IN fail, zfail, zpass
call fglstencilop(fail,zfail,zpass)
return
end subroutine f9xglstencilop

subroutine f9xgltexcoord1d(s)
real(kind=gldouble) INTENT_IN s
call fgltexcoord1d(s)
return
end subroutine f9xgltexcoord1d

subroutine f9xgltexcoord1dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord1dv(v)
return
end subroutine f9xgltexcoord1dv

subroutine f9xgltexcoord1f(s)
real(kind=glfloat) INTENT_IN s
call fgltexcoord1f(s)
return
end subroutine f9xgltexcoord1f

subroutine f9xgltexcoord1fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord1fv(v)
return
end subroutine f9xgltexcoord1fv

subroutine f9xgltexcoord1i(s)
integer(kind=glint) INTENT_IN s
call fgltexcoord1i(s)
return
end subroutine f9xgltexcoord1i

subroutine f9xgltexcoord1iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord1iv(v)
return
end subroutine f9xgltexcoord1iv

subroutine f9xgltexcoord1s(s)
integer(kind=glshort) INTENT_IN s
#ifdef NOSHORT
   call f9agltexcoord1s(s)
#else
   call fgltexcoord1s(s)
#endif
return
end subroutine f9xgltexcoord1s

subroutine f9xgltexcoord1sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9agltexcoord1sv(v)
#else
   call fgltexcoord1sv(v)
#endif
return
end subroutine f9xgltexcoord1sv

subroutine f9xgltexcoord2d(s,t)
real(kind=gldouble) INTENT_IN s, t
call fgltexcoord2d(s,t)
return
end subroutine f9xgltexcoord2d

subroutine f9xgltexcoord2dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord2dv(v)
return
end subroutine f9xgltexcoord2dv

subroutine f9xgltexcoord2f(s,t)
real(kind=glfloat) INTENT_IN s, t
call fgltexcoord2f(s,t)
return
end subroutine f9xgltexcoord2f

subroutine f9xgltexcoord2fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord2fv(v)
return
end subroutine f9xgltexcoord2fv

subroutine f9xgltexcoord2i(s,t)
integer(kind=glint) INTENT_IN s, t
call fgltexcoord2i(s,t)
return
end subroutine f9xgltexcoord2i

subroutine f9xgltexcoord2iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord2iv(v)
return
end subroutine f9xgltexcoord2iv

subroutine f9xgltexcoord2s(s,t)
integer(kind=glshort) INTENT_IN s, t
#ifdef NOSHORT
   call f9agltexcoord2s(s,t)
#else
   call fgltexcoord2s(s,t)
#endif
return
end subroutine f9xgltexcoord2s

subroutine f9xgltexcoord2sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9agltexcoord2sv(v)
#else
   call fgltexcoord2sv(v)
#endif
return
end subroutine f9xgltexcoord2sv

subroutine f9xgltexcoord3d(s,t,r)
real(kind=gldouble) INTENT_IN s, t, r
call fgltexcoord3d(s,t,r)
return
end subroutine f9xgltexcoord3d

subroutine f9xgltexcoord3dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord3dv(v)
return
end subroutine f9xgltexcoord3dv

subroutine f9xgltexcoord3f(s,t,r)
real(kind=glfloat) INTENT_IN s, t, r
call fgltexcoord3f(s,t,r)
return
end subroutine f9xgltexcoord3f

subroutine f9xgltexcoord3fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord3fv(v)
return
end subroutine f9xgltexcoord3fv

subroutine f9xgltexcoord3i(s,t,r)
integer(kind=glint) INTENT_IN s, t, r
call fgltexcoord3i(s,t,r)
return
end subroutine f9xgltexcoord3i

subroutine f9xgltexcoord3iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord3iv(v)
return
end subroutine f9xgltexcoord3iv

subroutine f9xgltexcoord3s(s,t,r)
integer(kind=glshort) INTENT_IN s, t, r
#ifdef NOSHORT
   call f9agltexcoord3s(s,t,r)
#else
   call fgltexcoord3s(s,t,r)
#endif
return
end subroutine f9xgltexcoord3s

subroutine f9xgltexcoord3sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9agltexcoord3sv(v)
#else
   call fgltexcoord3sv(v)
#endif
return
end subroutine f9xgltexcoord3sv

subroutine f9xgltexcoord4d(s,t,r,q)
real(kind=gldouble) INTENT_IN s, t, r, q
call fgltexcoord4d(s,t,r,q)
return
end subroutine f9xgltexcoord4d

subroutine f9xgltexcoord4dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord4dv(v)
return
end subroutine f9xgltexcoord4dv

subroutine f9xgltexcoord4f(s,t,r,q)
real(kind=glfloat) INTENT_IN s, t, r, q
call fgltexcoord4f(s,t,r,q)
return
end subroutine f9xgltexcoord4f

subroutine f9xgltexcoord4fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord4fv(v)
return
end subroutine f9xgltexcoord4fv

subroutine f9xgltexcoord4i(s,t,r,q)
integer(kind=glint) INTENT_IN s, t, r, q
call fgltexcoord4i(s,t,r,q)
return
end subroutine f9xgltexcoord4i

subroutine f9xgltexcoord4iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fgltexcoord4iv(v)
return
end subroutine f9xgltexcoord4iv

subroutine f9xgltexcoord4s(s,t,r,q)
integer(kind=glshort) INTENT_IN s, t, r, q
#ifdef NOSHORT
   call f9agltexcoord4s(s,t,r,q)
#else
   call fgltexcoord4s(s,t,r,q)
#endif
return
end subroutine f9xgltexcoord4s

subroutine f9xgltexcoord4sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9agltexcoord4sv(v)
#else
   call fgltexcoord4sv(v)
#endif
return
end subroutine f9xgltexcoord4sv

subroutine f9xgltexenvf(xtarget,pname,param)
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat) INTENT_IN param
call fgltexenvf(xtarget,pname,param)
return
end subroutine f9xgltexenvf

subroutine f9xgltexenvfv(xtarget,pname,params)
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
call fgltexenvfv(xtarget,pname,params)
return
end subroutine f9xgltexenvfv

subroutine f9xgltexenvi(xtarget,pname,param)
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint) INTENT_IN param
call fgltexenvi(xtarget,pname,param)
return
end subroutine f9xgltexenvi

subroutine f9xgltexenviv(xtarget,pname,params)
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
call fgltexenviv(xtarget,pname,params)
return
end subroutine f9xgltexenviv

subroutine f9xgltexgend(coord,pname,param)
integer(kind=glenum) INTENT_IN coord, pname
real(kind=gldouble) INTENT_IN param
call fgltexgend(coord,pname,param)
return
end subroutine f9xgltexgend

subroutine f9xgltexgendv(coord,pname,params)
integer(kind=glenum) INTENT_IN coord, pname
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN params
call fgltexgendv(coord,pname,params)
return
end subroutine f9xgltexgendv

subroutine f9xgltexgenf(coord,pname,param)
integer(kind=glenum) INTENT_IN coord, pname
real(kind=glfloat) INTENT_IN param
call fgltexgenf(coord,pname,param)
return
end subroutine f9xgltexgenf

subroutine f9xgltexgenfv(coord,pname,params)
integer(kind=glenum) INTENT_IN coord, pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
call fgltexgenfv(coord,pname,params)
return
end subroutine f9xgltexgenfv

subroutine f9xgltexgeni(coord,pname,param)
integer(kind=glenum) INTENT_IN coord, pname
integer(kind=glint) INTENT_IN param
call fgltexgeni(coord,pname,param)
return
end subroutine f9xgltexgeni

subroutine f9xgltexgeniv(coord,pname,params)
integer(kind=glenum) INTENT_IN coord, pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
call fgltexgeniv(coord,pname,params)
return
end subroutine f9xgltexgeniv

subroutine f9x0glteximage1d(xtarget,level,components,width,border,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN pixels
call f9y0glteximage1d(xtarget,level,components,width,border,format,xtype,pixels, &
   size_glint(pixels))
return
end subroutine f9x0glteximage1d

subroutine f9x1glteximage1d(xtarget,level,components,width,border,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
call f9y1glteximage1d(xtarget,level,components,width,border,format,xtype,pixels)
return
end subroutine f9x1glteximage1d

subroutine f9x4glteximage1d(xtarget,level,components,width,border,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i = glcptr_length
call f9y4glteximage1d(xtarget,level,components,width,border,format,xtype,temp,i)
return
end subroutine f9x4glteximage1d

subroutine f9x5glteximage1d(xtarget,level,components,width,border,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
call f9y5glteximage1d(xtarget,level,components,width,border,format,xtype,pixels)
return
end subroutine f9x5glteximage1d

subroutine f9x6glteximage1d(xtarget,level,components,width,border,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN pixels
call f9y6glteximage1d(xtarget,level,components,width,border,format,xtype,pixels, &
   size_glint(pixels))
return
end subroutine f9x6glteximage1d

subroutine f9x0glteximage2d(xtarget,level,components,width,height,border, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN pixels
call f9y0glteximage2d(xtarget,level,components,width,height,border,format,xtype, &
   pixels,size_glint(pixels))
return
end subroutine f9x0glteximage2d

subroutine f9x1glteximage2d(xtarget,level,components,width,height,border, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
call f9y1glteximage2d(xtarget,level,components,width,height,border,format,xtype, &
   pixels)
return
end subroutine f9x1glteximage2d

subroutine f9x4glteximage2d(xtarget,level,components,width,height,border, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i = glcptr_length
call f9y4glteximage2d(xtarget,level,components,width,height,border,format,xtype, &
   temp,i)
return
end subroutine f9x4glteximage2d

subroutine f9x5glteximage2d(xtarget,level,components,width,height,border, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
call f9y5glteximage2d(xtarget,level,components,width,height,border,format,xtype, &
   pixels)
return
end subroutine f9x5glteximage2d

subroutine f9x6glteximage2d(xtarget,level,components,width,height,border, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN pixels
call f9y6glteximage2d(xtarget,level,components,width,height,border,format,xtype, &
   pixels,size_glint(pixels))
return
end subroutine f9x6glteximage2d

subroutine f9xgltexparameterf(xtarget,pname,param)
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat) INTENT_IN param
call fgltexparameterf(xtarget,pname,param)
return
end subroutine f9xgltexparameterf

subroutine f9xgltexparameterfv(xtarget,pname,params)
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
call fgltexparameterfv(xtarget,pname,params)
return
end subroutine f9xgltexparameterfv

subroutine f9xgltexparameteri(xtarget,pname,param)
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint) INTENT_IN param
call fgltexparameteri(xtarget,pname,param)
return
end subroutine f9xgltexparameteri

subroutine f9xgltexparameteriv(xtarget,pname,params)
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
call fgltexparameteriv(xtarget,pname,params)
return
end subroutine f9xgltexparameteriv

subroutine f9xgltranslated(x,y,z)
real(kind=gldouble) INTENT_IN x, y, z
call fgltranslated(x,y,z)
return
end subroutine f9xgltranslated

subroutine f9xgltranslatef(x,y,z)
real(kind=glfloat) INTENT_IN x, y, z
call fgltranslatef(x,y,z)
return
end subroutine f9xgltranslatef

subroutine f9xglvertex2d(x,y)
real(kind=gldouble) INTENT_IN x, y
call fglvertex2d(x,y)
return
end subroutine f9xglvertex2d

subroutine f9xglvertex2dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fglvertex2dv(v)
return
end subroutine f9xglvertex2dv

subroutine f9xglvertex2f(x,y)
real(kind=glfloat) INTENT_IN x, y
call fglvertex2f(x,y)
return
end subroutine f9xglvertex2f

subroutine f9xglvertex2fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fglvertex2fv(v)
return
end subroutine f9xglvertex2fv

subroutine f9xglvertex2i(x,y)
integer(kind=glint) INTENT_IN x, y
call fglvertex2i(x,y)
return
end subroutine f9xglvertex2i

subroutine f9xglvertex2iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fglvertex2iv(v)
return
end subroutine f9xglvertex2iv

subroutine f9xglvertex2s(x,y)
integer(kind=glshort) INTENT_IN x, y
#ifdef NOSHORT
   call f9aglvertex2s(x,y)
#else
   call fglvertex2s(x,y)
#endif
return
end subroutine f9xglvertex2s

subroutine f9xglvertex2sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9aglvertex2sv(v)
#else
   call fglvertex2sv(v)
#endif
return
end subroutine f9xglvertex2sv

subroutine f9xglvertex3d(x,y,z)
real(kind=gldouble) INTENT_IN x, y, z
call fglvertex3d(x,y,z)
return
end subroutine f9xglvertex3d

subroutine f9xglvertex3dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fglvertex3dv(v)
return
end subroutine f9xglvertex3dv

subroutine f9xglvertex3f(x,y,z)
real(kind=glfloat) INTENT_IN x, y, z
call fglvertex3f(x,y,z)
return
end subroutine f9xglvertex3f

subroutine f9xglvertex3fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fglvertex3fv(v)
return
end subroutine f9xglvertex3fv

subroutine f9xglvertex3i(x,y,z)
integer(kind=glint) INTENT_IN x, y, z
call fglvertex3i(x,y,z)
return
end subroutine f9xglvertex3i

subroutine f9xglvertex3iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fglvertex3iv(v)
return
end subroutine f9xglvertex3iv

subroutine f9xglvertex3s(x,y,z)
integer(kind=glshort) INTENT_IN x, y, z
#ifdef NOSHORT
   call f9aglvertex3s(x,y,z)
#else
   call fglvertex3s(x,y,z)
#endif
return
end subroutine f9xglvertex3s

subroutine f9xglvertex3sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9aglvertex3sv(v)
#else
   call fglvertex3sv(v)
#endif
return
end subroutine f9xglvertex3sv

subroutine f9xglvertex4d(x,y,z,w)
real(kind=gldouble) INTENT_IN x, y, z, w
call fglvertex4d(x,y,z,w)
return
end subroutine f9xglvertex4d

subroutine f9xglvertex4dv(v)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
call fglvertex4dv(v)
return
end subroutine f9xglvertex4dv

subroutine f9xglvertex4f(x,y,z,w)
real(kind=glfloat) INTENT_IN x, y, z, w
call fglvertex4f(x,y,z,w)
return
end subroutine f9xglvertex4f

subroutine f9xglvertex4fv(v)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
call fglvertex4fv(v)
return
end subroutine f9xglvertex4fv

subroutine f9xglvertex4i(x,y,z,w)
integer(kind=glint) INTENT_IN x, y, z, w
call fglvertex4i(x,y,z,w)
return
end subroutine f9xglvertex4i

subroutine f9xglvertex4iv(v)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
call fglvertex4iv(v)
return
end subroutine f9xglvertex4iv

subroutine f9xglvertex4s(x,y,z,w)
integer(kind=glshort) INTENT_IN x, y, z, w
#ifdef NOSHORT
   call f9aglvertex4s(x,y,z,w)
#else
   call fglvertex4s(x,y,z,w)
#endif
return
end subroutine f9xglvertex4s

subroutine f9xglvertex4sv(v)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
#ifdef NOSHORT
   call f9aglvertex4sv(v)
#else
   call fglvertex4sv(v)
#endif
return
end subroutine f9xglvertex4sv

subroutine f9xglviewport(x,y,width,height)
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
call fglviewport(x,y,width,height)
return
end subroutine f9xglviewport

! OpenGL 1.1

#ifdef OPENGL_1_1

#ifdef PURE_FUNCTIONS
subroutine f9xglaretexturesresident(n,textures, &
   residences,areresident)
#else
function f9xglaretexturesresident(n,textures, &
   residences) result(resf9xglaretexturesresident)
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
logical(kind=glboolean), dimension(DEFER_DIM), intent(out) :: residences
#ifdef PURE_FUNCTIONS
logical(kind=glboolean), intent(out) :: areresident
#else
logical(kind=glboolean) :: resf9xglaretexturesresident
#endif
integer(kind=glint), dimension(n) :: int_residences
integer(kind=glint) :: int_result,i
if (glone_byte_logical) then
#ifdef PURE_FUNCTIONS
   call fglaretexturesresident(n,textures,residences,areresident)
#else
   resf9xglaretexturesresident = fglaretexturesresident(n,textures, &
      residences)
#endif
else
#ifdef PURE_FUNCTIONS
   call f9eglaretexturesresident(n,textures,int_residences,int_result)
#else
   int_result = f9eglaretexturesresident(n,textures,int_residences)
#endif
   do i=1,n
      residences(i) = (int_residences(i) /= 0)
   end do
#ifdef PURE_FUNCTIONS
   areresident = (int_result /= 0)
#else
   resf9xglaretexturesresident = (int_result /= 0)
#endif
end if
return
#ifdef PURE_FUNCTIONS
end subroutine f9xglaretexturesresident
#else
end function f9xglaretexturesresident
#endif

subroutine f9xglarrayelement(i)
integer(kind=glint) INTENT_IN i
call fglarrayelement(i)
return
end subroutine f9xglarrayelement

subroutine f9xglbindtexture(xtarget,texture)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=gluint) INTENT_IN texture
call fglbindtexture(xtarget,texture)
return
end subroutine f9xglbindtexture

subroutine f9x0glcolorpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y0glcolorpointer(xsize,xtype,stride,ptr)
return
end subroutine f9x0glcolorpointer

subroutine f9x1glcolorpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y1glcolorpointer(xsize,xtype,stride,ptr)
return
end subroutine f9x1glcolorpointer

subroutine f9x4glcolorpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
type(glcptr) INTENT_IN ptr
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(ptr%addr(i:i))
end do
i=glcptr_length
call f9y4glcolorpointer(xsize,xtype,stride,temp,i)
return
end subroutine f9x4glcolorpointer

subroutine f9x5glcolorpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glbyte), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOBYTE
#ifdef NOSHORT
call f9y0glcolorpointer(xsize,xtype,stride,ptr)
#else
call f9y6glcolorpointer(xsize,xtype,stride,ptr)
#endif
#else
call f9y5glcolorpointer(xsize,xtype,stride,ptr)
#endif
return
end subroutine f9x5glcolorpointer

subroutine f9x6glcolorpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glshort), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOSHORT
call f9y0glcolorpointer(xsize,xtype,stride,ptr)
#else
call f9y6glcolorpointer(xsize,xtype,stride,ptr)
#endif
return
end subroutine f9x6glcolorpointer

subroutine f9x7glcolorpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=gldouble), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y7glcolorpointer(xsize,xtype,stride,ptr)
return
end subroutine f9x7glcolorpointer

subroutine f9xglcopyteximage1d(xtarget,level,internalformat,x,y,width,border)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
call fglcopyteximage1d(xtarget,level,internalformat,x,y,width,border)
return
end subroutine f9xglcopyteximage1d

subroutine f9xglcopyteximage2d(xtarget,level,internalformat,x,y,width,height, &
   border)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
call fglcopyteximage2d(xtarget,level,internalformat,x,y,width,height,border)
return
end subroutine f9xglcopyteximage2d

subroutine f9xglcopytexsubimage1d(xtarget,level,xoffset,x,y,width)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, x, y
integer(kind=glsizei) INTENT_IN width
call fglcopytexsubimage1d(xtarget,level,xoffset,x,y,width)
return
end subroutine f9xglcopytexsubimage1d

subroutine f9xglcopytexsubimage2d(xtarget,level,xoffset,yoffset,x,y,width, &
   height)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, x, y
integer(kind=glsizei) INTENT_IN width, height
call fglcopytexsubimage2d(xtarget,level,xoffset,yoffset,x,y,width,height)
return
end subroutine f9xglcopytexsubimage2d

subroutine f9xgldeletetextures(n,textures)
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
call fgldeletetextures(n,textures)
return
end subroutine f9xgldeletetextures

subroutine f9xgldisableclientstate(cap)
integer(kind=glenum) INTENT_IN cap
call fgldisableclientstate(cap)
return
end subroutine f9xgldisableclientstate

subroutine f9xgldrawarrays(mode,first,xcount)
integer(kind=glenum) INTENT_IN mode
integer(kind=glint) INTENT_IN first
integer(kind=glsizei) INTENT_IN xcount
call fgldrawarrays(mode,first,xcount)
return
end subroutine f9xgldrawarrays

subroutine f9x0gldrawelements(mode,xcount,xtype,indices)
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN indices
call f9y0gldrawelements(mode,xcount,xtype,indices)
return
end subroutine f9x0gldrawelements

subroutine f9x4gldrawelements(mode,xcount,xtype,indices)
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=glsizei) INTENT_IN xcount
type(glcptr) INTENT_IN indices
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(indices%addr(i:i))
end do
i = glcptr_length
call f9y4gldrawelements(mode,xcount,xtype,temp,i)
return
end subroutine f9x4gldrawelements

subroutine f9x5gldrawelements(mode,xcount,xtype,indices)
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN indices
call f9y5gldrawelements(mode,xcount,xtype,indices)
return
end subroutine f9x5gldrawelements

subroutine f9x6gldrawelements(mode,xcount,xtype,indices)
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN indices
call f9y6gldrawelements(mode,xcount,xtype,indices)
return
end subroutine f9x6gldrawelements

subroutine f9xgledgeflagpointer(stride,ptr)
integer(kind=glsizei) INTENT_IN stride
logical(kind=glboolean), dimension(DEFER_DIM), target, intent(in out) :: ptr
call fgledgeflagpointer(stride,ptr)
return
end subroutine f9xgledgeflagpointer

subroutine f9xegledgeflagpointer(stride,ptr)
integer(kind=glsizei) INTENT_IN stride
integer(kind=glbyte), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9egledgeflagpointer(stride,ptr)
return
end subroutine f9xegledgeflagpointer

subroutine f9xglenableclientstate(cap)
integer(kind=glenum) INTENT_IN cap
call fglenableclientstate(cap)
return
end subroutine f9xglenableclientstate

subroutine f9xglgentextures(n,textures)
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: textures
call fglgentextures(n,textures)
return
end subroutine f9xglgentextures

subroutine f9xglgetpointerv(pname,params)
integer(kind=glenum) INTENT_IN pname
type(glcptr), intent(out) :: params
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
call fglgetpointerv(pname,temp,glcptr_length)
do i=1,glcptr_length
   params%addr(i:i) = char(temp(i))
end do
return
end subroutine f9xglgetpointerv

subroutine f9x0glindexpointer(xtype,stride,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y0glindexpointer(xtype,stride,ptr)
return
end subroutine f9x0glindexpointer

subroutine f9x1glindexpointer(xtype,stride,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y1glindexpointer(xtype,stride,ptr)
return
end subroutine f9x1glindexpointer

subroutine f9x4glindexpointer(xtype,stride,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
type(glcptr) INTENT_IN ptr
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(ptr%addr(i:i))
end do
i=glcptr_length
call f9y4glindexpointer(xtype,stride,temp,i)
return
end subroutine f9x4glindexpointer

subroutine f9x6glindexpointer(xtype,stride,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glshort), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOSHORT
call f9y0glindexpointer(xtype,stride,ptr)
#else
call f9y6glindexpointer(xtype,stride,ptr)
#endif
return
end subroutine f9x6glindexpointer

subroutine f9x7glindexpointer(xtype,stride,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=gldouble), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y7glindexpointer(xtype,stride,ptr)
return
end subroutine f9x7glindexpointer

subroutine f9xglindexub(c)
integer(kind=glubyte) INTENT_IN c
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglindexub(c)
#else
   call f9cglindexub(c)
#endif
#else
   call fglindexub(c)
#endif
return
end subroutine f9xglindexub

subroutine f9xglindexubv(c)
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN c
#ifdef NOBYTE
#ifdef NOSHORT
   call f9bglindexubv(c)
#else
   call f9cglindexubv(c)
#endif
#else
   call fglindexubv(c)
#endif
return
end subroutine f9xglindexubv

subroutine f9x0glinterleavedarrays(format,stride,ptr)
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN stride
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN ptr
call f9y0glinterleavedarrays(format,stride,ptr,size_glint(ptr))
return
end subroutine f9x0glinterleavedarrays

subroutine f9x1glinterleavedarrays(format,stride,ptr)
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN ptr
call f9y1glinterleavedarrays(format,stride,ptr)
return
end subroutine f9x1glinterleavedarrays

subroutine f9x4glinterleavedarrays(format,stride,ptr)
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN stride
type(glcptr) INTENT_IN ptr
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(ptr%addr(i:i))
end do
i = glcptr_length
call f9y4glinterleavedarrays(format,stride,temp,i)
return
end subroutine f9x4glinterleavedarrays

subroutine f9x5glinterleavedarrays(format,stride,ptr)
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN stride
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN ptr
call f9y5glinterleavedarrays(format,stride,ptr)
return
end subroutine f9x5glinterleavedarrays

subroutine f9x6glinterleavedarrays(format,stride,ptr)
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN stride
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN ptr
call f9y6glinterleavedarrays(format,stride,ptr,size_glint(ptr))
return
end subroutine f9x6glinterleavedarrays

function f9xglistexture(texture) result(resf9xglistexture)
integer(kind=gluint) INTENT_IN texture
logical(kind=glboolean) :: resf9xglistexture
integer(kind=glint) :: itemp
if (glone_byte_logical) then
   resf9xglistexture = fglistexture(texture)
else
   itemp = f9eglistexture(texture)
   resf9xglistexture = (itemp /= 0)
end if
return
end function f9xglistexture

subroutine f9x0glnormalpointer(xtype,stride,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y0glnormalpointer(xtype,stride,ptr)
return
end subroutine f9x0glnormalpointer

subroutine f9x1glnormalpointer(xtype,stride,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y1glnormalpointer(xtype,stride,ptr)
return
end subroutine f9x1glnormalpointer

subroutine f9x4glnormalpointer(xtype,stride,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
type(glcptr) INTENT_IN ptr
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(ptr%addr(i:i))
end do
i=glcptr_length
call f9y4glnormalpointer(xtype,stride,temp,i)
return
end subroutine f9x4glnormalpointer

subroutine f9x5glnormalpointer(xtype,stride,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glbyte), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOBYTE
#ifdef NOSHORT
call f9y0glnormalpointer(xtype,stride,ptr)
#else
call f9y6glnormalpointer(xtype,stride,ptr)
#endif
#else
call f9y5glnormalpointer(xtype,stride,ptr)
#endif
return
end subroutine f9x5glnormalpointer

subroutine f9x6glnormalpointer(xtype,stride,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glshort), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOSHORT
call f9y0glnormalpointer(xtype,stride,ptr)
#else
call f9y6glnormalpointer(xtype,stride,ptr)
#endif
return
end subroutine f9x6glnormalpointer

subroutine f9x7glnormalpointer(xtype,stride,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=gldouble), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y7glnormalpointer(xtype,stride,ptr)
return
end subroutine f9x7glnormalpointer

subroutine f9xglpolygonoffset(factor,units)
real(kind=glfloat) INTENT_IN factor, units
call fglpolygonoffset(factor,units)
return
end subroutine f9xglpolygonoffset

subroutine f9xglpopclientattrib()
call fglpopclientattrib()
return
end subroutine f9xglpopclientattrib

subroutine f9xglprioritizetextures(n,textures,priorities)
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
real(kind=glclampf), dimension(DEFER_DIM) INTENT_IN priorities
call fglprioritizetextures(n,textures,priorities)
return
end subroutine f9xglprioritizetextures

subroutine f9xglpushclientattrib(mask)
integer(kind=glbitfield) INTENT_IN mask
call fglpushclientattrib(mask)
return
end subroutine f9xglpushclientattrib

subroutine f9x0gltexcoordpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y0gltexcoordpointer(xsize,xtype,stride,ptr)
return
end subroutine f9x0gltexcoordpointer

subroutine f9x1gltexcoordpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y1gltexcoordpointer(xsize,xtype,stride,ptr)
return
end subroutine f9x1gltexcoordpointer

subroutine f9x4gltexcoordpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
type(glcptr) INTENT_IN ptr
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(ptr%addr(i:i))
end do
i=glcptr_length
call f9y4gltexcoordpointer(xsize,xtype,stride,temp,i)
return
end subroutine f9x4gltexcoordpointer

subroutine f9x6gltexcoordpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glshort), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOSHORT
call f9y0gltexcoordpointer(xsize,xtype,stride,ptr)
#else
call f9y6gltexcoordpointer(xsize,xtype,stride,ptr)
#endif
return
end subroutine f9x6gltexcoordpointer

subroutine f9x7gltexcoordpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=gldouble), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y7gltexcoordpointer(xsize,xtype,stride,ptr)
return
end subroutine f9x7gltexcoordpointer

subroutine f9x0gltexsubimage1d(xtarget,level,xoffset,width,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN pixels
call f9y0gltexsubimage1d(xtarget,level,xoffset,width,format,xtype,pixels, &
   size_glint(pixels))
return
end subroutine f9x0gltexsubimage1d

subroutine f9x1gltexsubimage1d(xtarget,level,xoffset,width,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
call f9y1gltexsubimage1d(xtarget,level,xoffset,width,format,xtype,pixels)
return
end subroutine f9x1gltexsubimage1d

subroutine f9x4gltexsubimage1d(xtarget,level,xoffset,width,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i=glcptr_length
call f9y4gltexsubimage1d(xtarget,level,xoffset,width,format,xtype,temp,i)
return
end subroutine f9x4gltexsubimage1d

subroutine f9x5gltexsubimage1d(xtarget,level,xoffset,width,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
call f9y5gltexsubimage1d(xtarget,level,xoffset,width,format,xtype,pixels)
return
end subroutine f9x5gltexsubimage1d

subroutine f9x6gltexsubimage1d(xtarget,level,xoffset,width,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN pixels
call f9y6gltexsubimage1d(xtarget,level,xoffset,width,format,xtype,pixels, &
   size_glint(pixels))
return
end subroutine f9x6gltexsubimage1d

subroutine f9x0gltexsubimage2d(xtarget,level,xoffset,yoffset,width,height, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN pixels
call f9y0gltexsubimage2d(xtarget,level,xoffset,yoffset,width,height,format, &
   xtype,pixels,size_glint(pixels))
return
end subroutine f9x0gltexsubimage2d

subroutine f9x1gltexsubimage2d(xtarget,level,xoffset,yoffset,width,height, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
call f9y1gltexsubimage2d(xtarget,level,xoffset,yoffset,width,height,format, &
   xtype,pixels)
return
end subroutine f9x1gltexsubimage2d

subroutine f9x4gltexsubimage2d(xtarget,level,xoffset,yoffset,width,height, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i=glcptr_length
call f9y4gltexsubimage2d(xtarget,level,xoffset,yoffset,width,height,format, &
   xtype,temp,i)
return
end subroutine f9x4gltexsubimage2d

subroutine f9x5gltexsubimage2d(xtarget,level,xoffset,yoffset,width,height, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
call f9y5gltexsubimage2d(xtarget,level,xoffset,yoffset,width,height,format, &
   xtype,pixels)
return
end subroutine f9x5gltexsubimage2d

subroutine f9x6gltexsubimage2d(xtarget,level,xoffset,yoffset,width,height, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN pixels
call f9y6gltexsubimage2d(xtarget,level,xoffset,yoffset,width,height,format, &
   xtype,pixels,size_glint(pixels))
return
end subroutine f9x6gltexsubimage2d

subroutine f9x0glvertexpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y0glvertexpointer(xsize,xtype,stride,ptr)
return
end subroutine f9x0glvertexpointer

subroutine f9x1glvertexpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y1glvertexpointer(xsize,xtype,stride,ptr)
return
end subroutine f9x1glvertexpointer

subroutine f9x4glvertexpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
type(glcptr) INTENT_IN ptr
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(ptr%addr(i:i))
end do
i=glcptr_length
call f9y4glvertexpointer(xsize,xtype,stride,temp,i)
return
end subroutine f9x4glvertexpointer

subroutine f9x6glvertexpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glshort), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOSHORT
call f9y0glvertexpointer(xsize,xtype,stride,ptr)
#else
call f9y6glvertexpointer(xsize,xtype,stride,ptr)
#endif
return
end subroutine f9x6glvertexpointer

subroutine f9x7glvertexpointer(xsize,xtype,stride,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=gldouble), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y7glvertexpointer(xsize,xtype,stride,ptr)
return
end subroutine f9x7glvertexpointer

#endif

! OpenGL 1.2

#ifdef OPENGL_1_2

subroutine f9x0gldrawrangeelements(mode,start,xend,xcount,xtype,indices)
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=gluint) INTENT_IN start,xend
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN indices
call f9y0gldrawrangeelements(mode,start,xend,xcount,xtype,indices)
return
end subroutine f9x0gldrawrangeelements

subroutine f9x4gldrawrangeelements(mode,start,xend,xcount,xtype,indices)
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=gluint) INTENT_IN start,xend
integer(kind=glsizei) INTENT_IN xcount
type(glcptr) INTENT_IN indices
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(indices%addr(i:i))
end do
i = glcptr_length
call f9y4gldrawrangeelements(mode,start,xend,xcount,xtype,temp,i)
return
end subroutine f9x4gldrawrangeelements

subroutine f9x5gldrawrangeelements(mode,start,xend,xcount,xtype,indices)
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=gluint) INTENT_IN start,xend
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN indices
call f9y5gldrawrangeelements(mode,start,xend,xcount,xtype,indices)
return
end subroutine f9x5gldrawrangeelements

subroutine f9x6gldrawrangeelements(mode,start,xend,xcount,xtype,indices)
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=gluint) INTENT_IN start,xend
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN indices
call f9y6gldrawrangeelements(mode,start,xend,xcount,xtype,indices)
return
end subroutine f9x6gldrawrangeelements

subroutine f9x0glteximage3d(xtarget,level,internalformat,width, &
   height,depth,border,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN pixels
call f9y0glteximage3d(xtarget,level,internalformat,width,height, &
   depth,border,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x0glteximage3d

subroutine f9x1glteximage3d(xtarget,level,internalformat,width, &
   height,depth,border,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
call f9y1glteximage3d(xtarget,level,internalformat,width,height, &
   depth,border,format,xtype,pixels)
return
end subroutine f9x1glteximage3d

subroutine f9x4glteximage3d(xtarget,level,internalformat,width, &
  height,depth,border,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i=glcptr_length
call f9y4glteximage3d(xtarget,level,internalformat,width,height, &
   depth,border,format,xtype,temp,i)
return
end subroutine f9x4glteximage3d

subroutine f9x5glteximage3d(xtarget,level,internalformat,width, &
   height,depth,border,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
call f9y5glteximage3d(xtarget,level,internalformat,width,height, &
   depth,border,format,xtype,pixels)
return
end subroutine f9x5glteximage3d

subroutine f9x6glteximage3d(xtarget,level,internalformat,width, &
   height,depth,border,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN pixels
call f9y6glteximage3d(xtarget,level,internalformat,width,height, &
   depth,border,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x6glteximage3d

subroutine f9x0gltexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,width, &
   height,depth,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN pixels
call f9y0gltexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,width,height, &
   depth,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x0gltexsubimage3d

subroutine f9x1gltexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,width, &
   height,depth,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
call f9y1gltexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,width,height, &
   depth,format,xtype,pixels)
return
end subroutine f9x1gltexsubimage3d

subroutine f9x4gltexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,width, &
  height,depth,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i=glcptr_length
call f9y4gltexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,width,height, &
   depth,format,xtype,temp,i)
return
end subroutine f9x4gltexsubimage3d

subroutine f9x5gltexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,width, &
   height,depth,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
call f9y5gltexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,width,height, &
   depth,format,xtype,pixels)
return
end subroutine f9x5gltexsubimage3d

subroutine f9x6gltexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,width, &
   height,depth,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN pixels
call f9y6gltexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,width,height, &
   depth,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x6gltexsubimage3d

subroutine f9xglcopytexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,x,y, &
   width,height)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset, x, y
integer(kind=glsizei) INTENT_IN width, height
call fglcopytexsubimage3d(xtarget,level,xoffset,yoffset,zoffset,x,y,width, &
   height)
return
end subroutine f9xglcopytexsubimage3d

#endif

! EXT_blend_color

#ifdef GL_EXT_blend_color

subroutine f9xglblendcolorext(red,green,blue,alpha)
real(kind=glclampf) INTENT_IN red, green, blue, alpha
call fglblendcolorext(red,green,blue,alpha)
return
end subroutine f9xglblendcolorext

#endif

! EXT_blend_minmax

#ifdef GL_EXT_blend_minmax

subroutine f9xglblendequationext(mode)
integer(kind=glenum) INTENT_IN mode
call fglblendequationext(mode)
return
end subroutine f9xglblendequationext

#endif

! EXT_color_table

#ifdef GL_EXT_color_table

subroutine f9x0glcolortableext(xtarget,internalformat,width,format,xtype,table)
integer(kind=glenum) INTENT_IN xtarget, internalformat
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN table
call f9y0glcolortableext(xtarget,internalformat,width,format,xtype,table, &
   size_glint(table))
return
end subroutine f9x0glcolortableext

subroutine f9x1glcolortableext(xtarget,internalformat,width,format,xtype,table)
integer(kind=glenum) INTENT_IN xtarget, internalformat
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN table
call f9y1glcolortableext(xtarget,internalformat,width,format,xtype,table)
return
end subroutine f9x1glcolortableext

subroutine f9x4glcolortableext(xtarget,internalformat,width,format,xtype,table)
integer(kind=glenum) INTENT_IN xtarget, internalformat
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN table
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(table%addr(i:i))
end do
i=glcptr_length
call f9y4glcolortableext(xtarget,internalformat,width,format,xtype,temp,i)
return
end subroutine f9x4glcolortableext

subroutine f9x5glcolortableext(xtarget,internalformat,width,format,xtype,table)
integer(kind=glenum) INTENT_IN xtarget, internalformat
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN table
call f9y5glcolortableext(xtarget,internalformat,width,format,xtype,table)
return
end subroutine f9x5glcolortableext

subroutine f9x6glcolortableext(xtarget,internalformat,width,format,xtype,table)
integer(kind=glenum) INTENT_IN xtarget, internalformat
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN table
call f9y6glcolortableext(xtarget,internalformat,width,format,xtype,table, &
   size_glint(table))
return
end subroutine f9x6glcolortableext

subroutine f9x0glcolorsubtableext(xtarget,start,xcount,format,xtype,data)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glsizei) INTENT_IN start, xcount
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN data
call f9y0glcolorsubtableext(xtarget,start,xcount,format,xtype,data, &
   size_glint(data))
return
end subroutine f9x0glcolorsubtableext

subroutine f9x1glcolorsubtableext(xtarget,start,xcount,format,xtype,data)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glsizei) INTENT_IN start, xcount
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN data
call f9y1glcolorsubtableext(xtarget,start,xcount,format,xtype,data)
return
end subroutine f9x1glcolorsubtableext

subroutine f9x4glcolorsubtableext(xtarget,start,xcount,format,xtype,data)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glsizei) INTENT_IN start, xcount
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN data
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(data%addr(i:i))
end do
i=glcptr_length
call f9y4glcolorsubtableext(xtarget,start,xcount,format,xtype,temp,i)
return
end subroutine f9x4glcolorsubtableext

subroutine f9x5glcolorsubtableext(xtarget,start,xcount,format,xtype,data)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glsizei) INTENT_IN start, xcount
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN data
call f9y5glcolorsubtableext(xtarget,start,xcount,format,xtype,data)
return
end subroutine f9x5glcolorsubtableext

subroutine f9x6glcolorsubtableext(xtarget,start,xcount,format,xtype,data)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glsizei) INTENT_IN start, xcount
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN data
call f9y6glcolorsubtableext(xtarget,start,xcount,format,xtype,data, &
   size_glint(data))
return
end subroutine f9x6glcolorsubtableext

subroutine f9x0glgetcolortableext(xtarget,format,xtype,table)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:), intent(out) :: table
call f9y0glgetcolortableext(xtarget,format,xtype,table,size_glint(table))
return
end subroutine f9x0glgetcolortableext

subroutine f9x1glgetcolortableext(xtarget,format,xtype,table)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
real(kind=glfloat), dimension(:), intent(out) :: table
call f9y1glgetcolortableext(xtarget,format,xtype,table,size_glint(table))
return
end subroutine f9x1glgetcolortableext

subroutine f9x4glgetcolortableext(xtarget,format,xtype,table)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr), intent(out) :: table
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
call f9y4glgetcolortableext(xtarget,format,xtype,temp,glcptr_length)
do i=1,glcptr_length
   table%addr(i:i) = char(temp(i))
end do
return
end subroutine f9x4glgetcolortableext

subroutine f9x5glgetcolortableext(xtarget,format,xtype,table)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glbyte), dimension(:), intent(out) :: table
call f9y5glgetcolortableext(xtarget,format,xtype,table,size_glint(table))
return
end subroutine f9x5glgetcolortableext

subroutine f9x6glgetcolortableext(xtarget,format,xtype,table)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:), intent(out) :: table
call f9y6glgetcolortableext(xtarget,format,xtype,table,size_glint(table))
return
end subroutine f9x6glgetcolortableext

subroutine f9xglgetcolortableparamfvext(xtarget,pname,params)
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
call fglgetcolortableparameterfvext(xtarget,pname,params)
return
end subroutine f9xglgetcolortableparamfvext

subroutine f9xglgetcolortableparamivext(xtarget,pname,params)
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
call fglgetcolortableparameterivext(xtarget,pname,params)
return
end subroutine f9xglgetcolortableparamivext

#endif

! EXT_point_parameters

#ifdef GL_EXT_point_parameters

subroutine f9xglpointparameterfext(pname,param)
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat) INTENT_IN param
call fglpointparameterfext(pname,param)
return
end subroutine f9xglpointparameterfext

subroutine f9xglpointparameterfvext(pname,params)
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
call fglpointparameterfvext(pname,params)
return
end subroutine f9xglpointparameterfvext

#endif

! EXT_polygon_offset

#ifdef GL_EXT_polygon_offset

subroutine f9xglpolygonoffsetext(factor,bias)
real(kind=glfloat) INTENT_IN factor, bias
call fglpolygonoffsetext(factor,bias)
return
end subroutine f9xglpolygonoffsetext

#endif

! EXT_subtexture

#ifdef GL_EXT_subtexture

subroutine f9x0gltexsubimage1dext(xtarget,level,xoffset,width,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN pixels
call f9y0gltexsubimage1dext(xtarget,level,xoffset,width,format,xtype,pixels, &
   size_glint(pixels))
return
end subroutine f9x0gltexsubimage1dext

subroutine f9x1gltexsubimage1dext(xtarget,level,xoffset,width,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
call f9y1gltexsubimage1dext(xtarget,level,xoffset,width,format,xtype,pixels)
return
end subroutine f9x1gltexsubimage1dext

subroutine f9x4gltexsubimage1dext(xtarget,level,xoffset,width,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i=glcptr_length
call f9y4gltexsubimage1dext(xtarget,level,xoffset,width,format,xtype,temp,i)
return
end subroutine f9x4gltexsubimage1dext

subroutine f9x5gltexsubimage1dext(xtarget,level,xoffset,width,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
call f9y5gltexsubimage1dext(xtarget,level,xoffset,width,format,xtype,pixels)
return
end subroutine f9x5gltexsubimage1dext

subroutine f9x6gltexsubimage1dext(xtarget,level,xoffset,width,format,xtype, &
   pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN pixels
call f9y6gltexsubimage1dext(xtarget,level,xoffset,width,format,xtype,pixels, &
   size_glint(pixels))
return
end subroutine f9x6gltexsubimage1dext

subroutine f9x0gltexsubimage2dext(xtarget,level,xoffset,yoffset,width,height, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN pixels
call f9y0gltexsubimage2dext(xtarget,level,xoffset,yoffset,width,height,format, &
   xtype,pixels,size_glint(pixels))
return
end subroutine f9x0gltexsubimage2dext

subroutine f9x1gltexsubimage2dext(xtarget,level,xoffset,yoffset,width,height, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
call f9y1gltexsubimage2dext(xtarget,level,xoffset,yoffset,width,height,format, &
   xtype,pixels)
return
end subroutine f9x1gltexsubimage2dext

subroutine f9x4gltexsubimage2dext(xtarget,level,xoffset,yoffset,width,height, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i=glcptr_length
call f9y4gltexsubimage2dext(xtarget,level,xoffset,yoffset,width,height,format, &
   xtype,temp,i)
return
end subroutine f9x4gltexsubimage2dext

subroutine f9x5gltexsubimage2dext(xtarget,level,xoffset,yoffset,width,height, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
call f9y5gltexsubimage2dext(xtarget,level,xoffset,yoffset,width,height,format, &
   xtype,pixels)
return
end subroutine f9x5gltexsubimage2dext

subroutine f9x6gltexsubimage2dext(xtarget,level,xoffset,yoffset,width,height, &
   format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN pixels
call f9y6gltexsubimage2dext(xtarget,level,xoffset,yoffset,width,height,format, &
   xtype,pixels,size_glint(pixels))
return
end subroutine f9x6gltexsubimage2dext

#endif

! EXT_texture3D

#ifdef GL_EXT_texture3D

subroutine f9x0glteximage3dext(xtarget,level,internalformat,width, &
   height,depth,border,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN pixels
call f9y0glteximage3dext(xtarget,level,internalformat,width,height, &
   depth,border,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x0glteximage3dext

subroutine f9x1glteximage3dext(xtarget,level,internalformat,width, &
   height,depth,border,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
call f9y1glteximage3dext(xtarget,level,internalformat,width,height, &
   depth,border,format,xtype,pixels)
return
end subroutine f9x1glteximage3dext

subroutine f9x4glteximage3dext(xtarget,level,internalformat,width, &
  height,depth,border,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i=glcptr_length
call f9y4glteximage3dext(xtarget,level,internalformat,width,height, &
   depth,border,format,xtype,temp,i)
return
end subroutine f9x4glteximage3dext

subroutine f9x5glteximage3dext(xtarget,level,internalformat,width, &
   height,depth,border,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
call f9y5glteximage3dext(xtarget,level,internalformat,width,height, &
   depth,border,format,xtype,pixels)
return
end subroutine f9x5glteximage3dext

subroutine f9x6glteximage3dext(xtarget,level,internalformat,width, &
   height,depth,border,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN pixels
call f9y6glteximage3dext(xtarget,level,internalformat,width,height, &
   depth,border,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x6glteximage3dext

subroutine f9x0gltexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,width, &
   height,depth,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN pixels
call f9y0gltexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,width,height, &
   depth,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x0gltexsubimage3dext

subroutine f9x1gltexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,width, &
   height,depth,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
call f9y1gltexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,width,height, &
   depth,format,xtype,pixels)
return
end subroutine f9x1gltexsubimage3dext

subroutine f9x4gltexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,width, &
  height,depth,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN pixels
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(pixels%addr(i:i))
end do
i=glcptr_length
call f9y4gltexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,width,height, &
   depth,format,xtype,temp,i)
return
end subroutine f9x4gltexsubimage3dext

subroutine f9x5gltexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,width, &
   height,depth,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
call f9y5gltexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,width,height, &
   depth,format,xtype,pixels)
return
end subroutine f9x5gltexsubimage3dext

subroutine f9x6gltexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,width, &
   height,depth,format,xtype,pixels)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN pixels
call f9y6gltexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,width,height, &
   depth,format,xtype,pixels,size_glint(pixels))
return
end subroutine f9x6gltexsubimage3dext

subroutine f9xglcopytexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,x,y, &
   width,height)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset, x, y
integer(kind=glsizei) INTENT_IN width, height
call fglcopytexsubimage3dext(xtarget,level,xoffset,yoffset,zoffset,x,y,width, &
   height)
return
end subroutine f9xglcopytexsubimage3dext

#endif

! EXT_texture_object

#ifdef GL_EXT_texture_object

subroutine f9xglgentexturesext(n,textures)
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: textures
call fglgentexturesext(n,textures)
return
end subroutine f9xglgentexturesext

subroutine f9xgldeletetexturesext(n,textures)
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
call fgldeletetexturesext(n,textures)
return
end subroutine f9xgldeletetexturesext

subroutine f9xglbindtextureext(xtarget,texture)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=gluint) INTENT_IN texture
call fglbindtextureext(xtarget,texture)
return
end subroutine f9xglbindtextureext

subroutine f9xglprioritizetexturesext(n,textures,priorities)
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
real(kind=glclampf), dimension(DEFER_DIM) INTENT_IN priorities
call fglprioritizetexturesext(n,textures,priorities)
return
end subroutine f9xglprioritizetexturesext

#ifdef PURE_FUNCTIONS
subroutine f9xglaretexturesresidentext(n,textures, &
   residences,areresident)
#else
function f9xglaretexturesresidentext(n,textures, &
   residences) result(resf9xglaretexturesresidentext)
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
logical(kind=glboolean), dimension(DEFER_DIM), intent(out) :: residences
#ifdef PURE_FUNCTIONS
logical(kind=glboolean), intent(out) :: areresident
#else
logical(kind=glboolean) :: resf9xglaretexturesresidentext
#endif
integer(kind=glint), dimension(n) :: int_residences
integer(kind=glint) :: int_result,i
if (glone_byte_logical) then
#ifdef PURE_FUNCTIONS
   call fglaretexturesresidentext(n,textures,residences,areresident)
#else
   resf9xglaretexturesresidentext = fglaretexturesresidentext(n,textures, &
      residences)
#endif
else
#ifdef PURE_FUNCTIONS
   call f9eglaretexturesresidentext(n,textures,int_residences,int_result)
#else
   int_result = f9eglaretexturesresidentext(n,textures,int_residences)
#endif
   do i=1,n
      residences(i) = (int_residences(i) /= 0)
   end do
#ifdef PURE_FUNCTIONS
   areresident = (int_result /= 0)
#else
   resf9xglaretexturesresidentext = (int_result /= 0)
#endif
end if
return
#ifdef PURE_FUNCTIONS
end subroutine f9xglaretexturesresidentext
#else
end function f9xglaretexturesresidentext
#endif

function f9xglistextureext(texture) result(resf9xglistextureext)
integer(kind=gluint) INTENT_IN texture
logical(kind=glboolean) :: resf9xglistextureext
integer(kind=glint) :: itemp
if (glone_byte_logical) then
   resf9xglistextureext = fglistextureext(texture)
else
   itemp = f9eglistextureext(texture)
   resf9xglistextureext = (itemp /= 0)
end if
return
end function f9xglistextureext

#endif

! EXT_vertex_array

#ifdef GL_EXT_vertex_array

subroutine f9x0glvertexpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y0glvertexpointerext(xsize,xtype,stride,xcount,ptr)
return
end subroutine f9x0glvertexpointerext

subroutine f9x1glvertexpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y1glvertexpointerext(xsize,xtype,stride,xcount,ptr)
return
end subroutine f9x1glvertexpointerext

subroutine f9x4glvertexpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
type(glcptr) INTENT_IN ptr
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(ptr%addr(i:i))
end do
i=glcptr_length
call f9y4glvertexpointerext(xsize,xtype,stride,xcount,temp,i)
return
end subroutine f9x4glvertexpointerext

subroutine f9x6glvertexpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glshort), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOSHORT
call f9y0glvertexpointerext(xsize,xtype,stride,xcount,ptr)
#else
call f9y6glvertexpointerext(xsize,xtype,stride,xcount,ptr)
#endif
return
end subroutine f9x6glvertexpointerext

subroutine f9x7glvertexpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=gldouble), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y7glvertexpointerext(xsize,xtype,stride,xcount,ptr)
return
end subroutine f9x7glvertexpointerext

subroutine f9x0glnormalpointerext(xtype,stride,xcount,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y0glnormalpointerext(xtype,stride,xcount,ptr)
return
end subroutine f9x0glnormalpointerext

subroutine f9x1glnormalpointerext(xtype,stride,xcount,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y1glnormalpointerext(xtype,stride,xcount,ptr)
return
end subroutine f9x1glnormalpointerext

subroutine f9x4glnormalpointerext(xtype,stride,xcount,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
type(glcptr) INTENT_IN ptr
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(ptr%addr(i:i))
end do
i=glcptr_length
call f9y4glnormalpointerext(xtype,stride,xcount,temp,i)
return
end subroutine f9x4glnormalpointerext

subroutine f9x5glnormalpointerext(xtype,stride,xcount,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glbyte), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOBYTE
#ifdef NOSHORT
call f9y0glnormalpointerext(xtype,stride,xcount,ptr)
#else
call f9y6glnormalpointerext(xtype,stride,xcount,ptr)
#endif
#else
call f9y5glnormalpointerext(xtype,stride,xcount,ptr)
#endif
return
end subroutine f9x5glnormalpointerext

subroutine f9x6glnormalpointerext(xtype,stride,xcount,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glshort), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOSHORT
call f9y0glnormalpointerext(xtype,stride,xcount,ptr)
#else
call f9y6glnormalpointerext(xtype,stride,xcount,ptr)
#endif
return
end subroutine f9x6glnormalpointerext

subroutine f9x7glnormalpointerext(xtype,stride,xcount,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=gldouble), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y7glnormalpointerext(xtype,stride,xcount,ptr)
return
end subroutine f9x7glnormalpointerext

subroutine f9x0glcolorpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y0glcolorpointerext(xsize,xtype,stride,xcount,ptr)
return
end subroutine f9x0glcolorpointerext

subroutine f9x1glcolorpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y1glcolorpointerext(xsize,xtype,stride,xcount,ptr)
return
end subroutine f9x1glcolorpointerext

subroutine f9x4glcolorpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
type(glcptr) INTENT_IN ptr
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(ptr%addr(i:i))
end do
i=glcptr_length
call f9y4glcolorpointerext(xsize,xtype,stride,xcount,temp,i)
return
end subroutine f9x4glcolorpointerext

subroutine f9x5glcolorpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glbyte), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOBYTE
#ifdef NOSHORT
call f9y0glcolorpointerext(xsize,xtype,stride,xcount,ptr)
#else
call f9y6glcolorpointerext(xsize,xtype,stride,xcount,ptr)
#endif
#else
call f9y5glcolorpointerext(xsize,xtype,stride,xcount,ptr)
#endif
return
end subroutine f9x5glcolorpointerext

subroutine f9x6glcolorpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glshort), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOSHORT
call f9y0glcolorpointerext(xsize,xtype,stride,xcount,ptr)
#else
call f9y6glcolorpointerext(xsize,xtype,stride,xcount,ptr)
#endif
return
end subroutine f9x6glcolorpointerext

subroutine f9x7glcolorpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=gldouble), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y7glcolorpointerext(xsize,xtype,stride,xcount,ptr)
return
end subroutine f9x7glcolorpointerext

subroutine f9x0glindexpointerext(xtype,stride,xcount,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y0glindexpointerext(xtype,stride,xcount,ptr)
return
end subroutine f9x0glindexpointerext

subroutine f9x1glindexpointerext(xtype,stride,xcount,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y1glindexpointerext(xtype,stride,xcount,ptr)
return
end subroutine f9x1glindexpointerext

subroutine f9x4glindexpointerext(xtype,stride,xcount,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
type(glcptr) INTENT_IN ptr
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(ptr%addr(i:i))
end do
i=glcptr_length
call f9y4glindexpointerext(xtype,stride,xcount,temp,i)
return
end subroutine f9x4glindexpointerext

subroutine f9x6glindexpointerext(xtype,stride,xcount,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glshort), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOSHORT
call f9y0glindexpointerext(xtype,stride,xcount,ptr)
#else
call f9y6glindexpointerext(xtype,stride,xcount,ptr)
#endif
return
end subroutine f9x6glindexpointerext

subroutine f9x7glindexpointerext(xtype,stride,xcount,ptr)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=gldouble), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y7glindexpointerext(xtype,stride,xcount,ptr)
return
end subroutine f9x7glindexpointerext

subroutine f9x0gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y0gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
return
end subroutine f9x0gltexcoordpointerext

subroutine f9x1gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y1gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
return
end subroutine f9x1gltexcoordpointerext

subroutine f9x4gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
type(glcptr) INTENT_IN ptr
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(ptr%addr(i:i))
end do
i=glcptr_length
call f9y4gltexcoordpointerext(xsize,xtype,stride,xcount,temp,i)
return
end subroutine f9x4gltexcoordpointerext

subroutine f9x6gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glshort), dimension(DEFER_DIM), target, intent(in out) :: ptr
#ifdef NOSHORT
call f9y0gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
#else
call f9y6gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
#endif
return
end subroutine f9x6gltexcoordpointerext

subroutine f9x7gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=gldouble), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9y7gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
return
end subroutine f9x7gltexcoordpointerext

subroutine f9xgledgeflagpointerext(stride,xcount,ptr)
integer(kind=glsizei) INTENT_IN stride, xcount
logical(kind=glboolean), dimension(DEFER_DIM), target, intent(in out) :: ptr
call fgledgeflagpointerext(stride,xcount,ptr)
return
end subroutine f9xgledgeflagpointerext

subroutine f9xegledgeflagpointerext(stride,xcount,ptr)
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glbyte), dimension(DEFER_DIM), target, intent(in out) :: ptr
call f9egledgeflagpointerext(stride,xcount,ptr)
return
end subroutine f9xegledgeflagpointerext

subroutine f9xglgetpointervext(pname,params)
integer(kind=glenum) INTENT_IN pname
type(glcptr), intent(out) :: params
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
call fglgetpointervext(pname,temp,glcptr_length)
do i=1,glcptr_length
   params%addr(i:i) = char(temp(i))
end do
return
end subroutine f9xglgetpointervext

subroutine f9xglarrayelementext(i)
integer(kind=glint) INTENT_IN i
call fglarrayelementext(i)
return
end subroutine f9xglarrayelementext

subroutine f9xgldrawarraysext(mode,first,xcount)
integer(kind=glenum) INTENT_IN mode
integer(kind=glint) INTENT_IN first
integer(kind=glsizei) INTENT_IN xcount
call fgldrawarraysext(mode,first,xcount)
return
end subroutine f9xgldrawarraysext

#endif

! SGIS_multisample

#ifdef GL_SGIS_multisample

subroutine f9xglsamplemasksgis(value,invert)
real(kind=glclampf) INTENT_IN value
logical(kind=glboolean) INTENT_IN invert
integer(kind=glint) :: int_invert
if (glone_byte_logical) then
   call fglsamplemasksgis(value,invert)
else
   if (invert) then
      int_invert = 1
   else
      int_invert = 0
   endif
   call f9eglsamplemasksgis(value,int_invert)
endif
return
end subroutine f9xglsamplemasksgis

subroutine f9xglsamplepatternsgis(pattern)
integer(kind=glenum) INTENT_IN pattern
call fglsamplepatternsgis(pattern)
return
end subroutine f9xglsamplepatternsgis

#endif

! MESA_resize_buffers

#ifdef GL_MESA_resize_buffers

subroutine f9xglresizebuffersmesa()
call fglresizebuffersmesa()
return
end subroutine f9xglresizebuffersmesa

#endif

! MESA_window_pos

#ifdef GL_MESA_window_pos

subroutine f9xglwindowpos2imesa(x,y)
integer(kind=glint) INTENT_IN x, y
call fglwindowpos2imesa(x,y)
return
end subroutine f9xglwindowpos2imesa

subroutine f9xglwindowpos2smesa(x,y)
integer(kind=glshort) INTENT_IN x, y
#ifdef NOSHORT
   call f9aglwindowpos2smesa(x,y)
#else
   call fglwindowpos2smesa(x,y)
#endif
return
end subroutine f9xglwindowpos2smesa

subroutine f9xglwindowpos2fmesa(x,y)
real(kind=glfloat) INTENT_IN x, y
call fglwindowpos2fmesa(x,y)
return
end subroutine f9xglwindowpos2fmesa

subroutine f9xglwindowpos2dmesa(x,y)
real(kind=gldouble) INTENT_IN x, y
call fglwindowpos2dmesa(x,y)
return
end subroutine f9xglwindowpos2dmesa

subroutine f9xglwindowpos2ivmesa(p)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN p
call fglwindowpos2ivmesa(p)
return
end subroutine f9xglwindowpos2ivmesa

subroutine f9xglwindowpos2svmesa(p)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN p
#ifdef NOSHORT
   call f9aglwindowpos2svmesa(p)
#else
   call fglwindowpos2svmesa(p)
#endif
return
end subroutine f9xglwindowpos2svmesa

subroutine f9xglwindowpos2fvmesa(p)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN p
call fglwindowpos2fvmesa(p)
return
end subroutine f9xglwindowpos2fvmesa

subroutine f9xglwindowpos2dvmesa(p)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN p
call fglwindowpos2dvmesa(p)
return
end subroutine f9xglwindowpos2dvmesa

subroutine f9xglwindowpos3imesa(x,y,z)
integer(kind=glint) INTENT_IN x, y, z
call fglwindowpos3imesa(x,y,z)
return
end subroutine f9xglwindowpos3imesa

subroutine f9xglwindowpos3smesa(x,y,z)
integer(kind=glshort) INTENT_IN x, y, z
#ifdef NOSHORT
   call f9aglwindowpos3smesa(x,y,z)
#else
   call fglwindowpos3smesa(x,y,z)
#endif
return
end subroutine f9xglwindowpos3smesa

subroutine f9xglwindowpos3fmesa(x,y,z)
real(kind=glfloat) INTENT_IN x, y, z
call fglwindowpos3fmesa(x,y,z)
return
end subroutine f9xglwindowpos3fmesa

subroutine f9xglwindowpos3dmesa(x,y,z)
real(kind=gldouble) INTENT_IN x, y, z
call fglwindowpos3dmesa(x,y,z)
return
end subroutine f9xglwindowpos3dmesa

subroutine f9xglwindowpos3ivmesa(p)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN p
call fglwindowpos3ivmesa(p)
return
end subroutine f9xglwindowpos3ivmesa

subroutine f9xglwindowpos3svmesa(p)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN p
#ifdef NOSHORT
   call f9aglwindowpos3svmesa(p)
#else
   call fglwindowpos3svmesa(p)
#endif
return
end subroutine f9xglwindowpos3svmesa

subroutine f9xglwindowpos3fvmesa(p)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN p
call fglwindowpos3fvmesa(p)
return
end subroutine f9xglwindowpos3fvmesa

subroutine f9xglwindowpos3dvmesa(p)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN p
call fglwindowpos3dvmesa(p)
return
end subroutine f9xglwindowpos3dvmesa

subroutine f9xglwindowpos4imesa(x,y,z,w)
integer(kind=glint) INTENT_IN x, y, z, w
call fglwindowpos4imesa(x,y,z,w)
return
end subroutine f9xglwindowpos4imesa

subroutine f9xglwindowpos4smesa(x,y,z,w)
integer(kind=glshort) INTENT_IN x, y, z, w
#ifdef NOSHORT
   call f9aglwindowpos4smesa(x,y,z,w)
#else
   call fglwindowpos4smesa(x,y,z,w)
#endif
return
end subroutine f9xglwindowpos4smesa

subroutine f9xglwindowpos4fmesa(x,y,z,w)
real(kind=glfloat) INTENT_IN x, y, z, w
call fglwindowpos4fmesa(x,y,z,w)
return
end subroutine f9xglwindowpos4fmesa

subroutine f9xglwindowpos4dmesa(x,y,z,w)
real(kind=gldouble) INTENT_IN x, y, z, w
call fglwindowpos4dmesa(x,y,z,w)
return
end subroutine f9xglwindowpos4dmesa

subroutine f9xglwindowpos4ivmesa(p)
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN p
call fglwindowpos4ivmesa(p)
return
end subroutine f9xglwindowpos4ivmesa

subroutine f9xglwindowpos4svmesa(p)
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN p
#ifdef NOSHORT
   call f9aglwindowpos4svmesa(p)
#else
   call fglwindowpos4svmesa(p)
#endif
return
end subroutine f9xglwindowpos4svmesa

subroutine f9xglwindowpos4fvmesa(p)
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN p
call fglwindowpos4fvmesa(p)
return
end subroutine f9xglwindowpos4fvmesa

subroutine f9xglwindowpos4dvmesa(p)
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN p
call fglwindowpos4dvmesa(p)
return
end subroutine f9xglwindowpos4dvmesa

#endif

end module opengl_fwrap

!--------------------------
! The primary OpenGL module
!--------------------------

module opengl_gl
use opengl_kinds
use opengl_fwrap
public
end module opengl_gl
