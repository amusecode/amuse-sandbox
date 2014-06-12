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


#include "fppr.inc"

module opengl_glinterfaces

! interface blocks for the C wrappers.

#ifdef F_LANG
! F does not allow multiple instances of implicit none
#define IMPLICIT_NONE
#else
#ifdef DVF
! Digital visual fortran does not allow implicit none in interface blocks
#define IMPLICIT_NONE
#else
#define IMPLICIT_NONE implicit none
#endif
#endif


#ifdef NASOFTWARE
#define INTENT_IN ::
#else
#define INTENT_IN ,intent(in) ::
#endif

#ifdef OLDXLF
! xlf before 3.2.5.2 and 4.1.0.1 doesn't allow assumed size target array (in
! fwrap) to be passed to an assumed size target array
#define TARGET
#else
#define TARGET target,
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

use opengl_kinds

#ifdef SALFORD

! define STDCALL name mangling for Salford compiler

 STDCALL fglaccum 'fglaccum' (REF,REF)
 STDCALL fglalphafunc 'fglalphafunc' (REF,REF)
 STDCALL fglbegin 'fglbegin' (REF)
 STDCALL fglbitmap 'fglbitmap' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9bglbitmap 'f9bglbitmap' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9cglbitmap 'f9cglbitmap' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglblendfunc 'fglblendfunc' (REF,REF)
 STDCALL fglcalllist 'fglcalllist' (REF)
 STDCALL f9y0glcalllists 'f9y0glcalllists' (REF,REF,REF)
 STDCALL f9y1glcalllists 'f9y1glcalllists' (REF,REF,REF)
 STDCALL f9y4glcalllists 'f9y4glcalllists' (REF,REF,REF,REF)
 STDCALL f9y5glcalllists 'f9y5glcalllists' (REF,REF,REF)
 STDCALL f9y6glcalllists 'f9y6glcalllists' (REF,REF,REF)
 STDCALL fglclear 'fglclear' (REF)
 STDCALL fglclearaccum 'fglclearaccum' (REF,REF,REF,REF)
 STDCALL fglclearcolor 'fglclearcolor' (REF,REF,REF,REF)
 STDCALL fglcleardepth 'fglcleardepth' (REF)
 STDCALL fglclearindex 'fglclearindex' (REF)
 STDCALL fglclearstencil 'fglclearstencil' (REF)
 STDCALL fglclipplane 'fglclipplane' (REF,REF)
 STDCALL fglcolor3b 'fglcolor3b' (REF,REF,REF)
 STDCALL f9bglcolor3b 'f9bglcolor3b' (REF,REF,REF)
 STDCALL f9cglcolor3b 'f9cglcolor3b' (REF,REF,REF)
 STDCALL fglcolor3bv 'fglcolor3bv' (REF)
 STDCALL f9bglcolor3bv 'f9bglcolor3bv' (REF)
 STDCALL f9cglcolor3bv 'f9cglcolor3bv' (REF)
 STDCALL fglcolor3d 'fglcolor3d' (REF,REF,REF)
 STDCALL fglcolor3dv 'fglcolor3dv' (REF)
 STDCALL fglcolor3f 'fglcolor3f' (REF,REF,REF)
 STDCALL fglcolor3fv 'fglcolor3fv' (REF)
 STDCALL fglcolor3i 'fglcolor3i' (REF,REF,REF)
 STDCALL fglcolor3iv 'fglcolor3iv' (REF)
 STDCALL fglcolor3s 'fglcolor3s' (REF,REF,REF)
 STDCALL f9aglcolor3s 'f9aglcolor3s' (REF,REF,REF)
 STDCALL fglcolor3sv 'fglcolor3sv' (REF)
 STDCALL f9aglcolor3sv 'f9aglcolor3sv' (REF)
 STDCALL fglcolor3ub 'fglcolor3ub' (REF,REF,REF)
 STDCALL f9bglcolor3ub 'f9bglcolor3ub' (REF,REF,REF)
 STDCALL f9cglcolor3ub 'f9cglcolor3ub' (REF,REF,REF)
 STDCALL fglcolor3ubv 'fglcolor3ubv' (REF)
 STDCALL f9bglcolor3ubv 'f9bglcolor3ubv' (REF)
 STDCALL f9cglcolor3ubv 'f9cglcolor3ubv' (REF)
 STDCALL fglcolor3ui 'fglcolor3ui' (REF,REF,REF)
 STDCALL fglcolor3uiv 'fglcolor3uiv' (REF)
 STDCALL fglcolor3us 'fglcolor3us' (REF,REF,REF)
 STDCALL f9aglcolor3us 'f9aglcolor3us' (REF,REF,REF)
 STDCALL fglcolor3usv 'fglcolor3usv' (REF)
 STDCALL f9aglcolor3usv 'f9aglcolor3usv' (REF)
 STDCALL fglcolor4b 'fglcolor4b' (REF,REF,REF,REF)
 STDCALL f9bglcolor4b 'f9bglcolor4b' (REF,REF,REF,REF)
 STDCALL f9cglcolor4b 'f9cglcolor4b' (REF,REF,REF,REF)
 STDCALL fglcolor4bv 'fglcolor4bv' (REF)
 STDCALL f9bglcolor4bv 'f9bglcolor4bv' (REF)
 STDCALL f9cglcolor4bv 'f9cglcolor4bv' (REF)
 STDCALL fglcolor4d 'fglcolor4d' (REF,REF,REF,REF)
 STDCALL fglcolor4dv 'fglcolor4dv' (REF)
 STDCALL fglcolor4f 'fglcolor4f' (REF,REF,REF,REF)
 STDCALL fglcolor4fv 'fglcolor4fv' (REF)
 STDCALL fglcolor4i 'fglcolor4i' (REF,REF,REF,REF)
 STDCALL fglcolor4iv 'fglcolor4iv' (REF)
 STDCALL fglcolor4s 'fglcolor4s' (REF,REF,REF,REF)
 STDCALL f9aglcolor4s 'f9aglcolor4s' (REF,REF,REF,REF)
 STDCALL fglcolor4sv 'fglcolor4sv' (REF)
 STDCALL f9aglcolor4sv 'f9aglcolor4sv' (REF)
 STDCALL fglcolor4ub 'fglcolor4ub' (REF,REF,REF,REF)
 STDCALL f9bglcolor4ub 'f9bglcolor4ub' (REF,REF,REF,REF)
 STDCALL f9cglcolor4ub 'f9cglcolor4ub' (REF,REF,REF,REF)
 STDCALL fglcolor4ubv 'fglcolor4ubv' (REF)
 STDCALL f9bglcolor4ubv 'f9bglcolor4ubv' (REF)
 STDCALL f9cglcolor4ubv 'f9cglcolor4ubv' (REF)
 STDCALL fglcolor4ui 'fglcolor4ui' (REF,REF,REF,REF)
 STDCALL fglcolor4uiv 'fglcolor4uiv' (REF)
 STDCALL fglcolor4us 'fglcolor4us' (REF,REF,REF,REF)
 STDCALL f9aglcolor4us 'f9aglcolor4us' (REF,REF,REF,REF)
 STDCALL fglcolor4usv 'fglcolor4usv' (REF)
 STDCALL f9aglcolor4usv 'f9aglcolor4usv' (REF)
 STDCALL fglcolormask 'fglcolormask' (REF,REF,REF,REF)
 STDCALL f9eglcolormask 'f9eglcolormask' (REF,REF,REF,REF)
 STDCALL fglcolormaterial 'fglcolormaterial' (REF,REF)
 STDCALL fglcopypixels 'fglcopypixels' (REF,REF,REF,REF,REF)
 STDCALL fglcullface 'fglcullface' (REF)
 STDCALL fgldeletelists 'fgldeletelists' (REF,REF)
 STDCALL fgldepthfunc 'fgldepthfunc' (REF)
 STDCALL fgldepthmask 'fgldepthmask' (REF)
 STDCALL f9egldepthmask 'f9egldepthmask' (REF)
 STDCALL fgldepthrange 'fgldepthrange' (REF,REF)
 STDCALL fgldisable 'fgldisable' (REF)
 STDCALL fgldrawbuffer 'fgldrawbuffer' (REF)
 STDCALL f9y0gldrawpixels 'f9y0gldrawpixels' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y1gldrawpixels 'f9y1gldrawpixels' (REF,REF,REF,REF,REF)
 STDCALL f9y4gldrawpixels 'f9y4gldrawpixels' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y5gldrawpixels 'f9y5gldrawpixels' (REF,REF,REF,REF,REF)
 STDCALL f9y6gldrawpixels 'f9y6gldrawpixels' (REF,REF,REF,REF,REF,REF)
 STDCALL fgledgeflag 'fgledgeflag' (REF)
 STDCALL f9egledgeflag 'f9egledgeflag' (REF)
 STDCALL fgledgeflagv 'fgledgeflagv' (REF)
 STDCALL f9egledgeflagv 'f9egledgeflagv' (REF,REF)
 STDCALL fglenable 'fglenable' (REF)
 STDCALL fglend 'fglend'
 STDCALL fglendlist 'fglendlist'
 STDCALL fglevalcoord1d 'fglevalcoord1d' (REF)
 STDCALL fglevalcoord1dv 'fglevalcoord1dv' (REF)
 STDCALL fglevalcoord1f 'fglevalcoord1f' (REF)
 STDCALL fglevalcoord1fv 'fglevalcoord1fv' (REF)
 STDCALL fglevalcoord2d 'fglevalcoord2d' (REF,REF)
 STDCALL fglevalcoord2dv 'fglevalcoord2dv' (REF)
 STDCALL fglevalcoord2f 'fglevalcoord2f' (REF,REF)
 STDCALL fglevalcoord2fv 'fglevalcoord2fv' (REF)
 STDCALL fglevalmesh1 'fglevalmesh1' (REF,REF,REF)
 STDCALL fglevalmesh2 'fglevalmesh2' (REF,REF,REF,REF,REF)
 STDCALL fglevalpoint1 'fglevalpoint1' (REF)
 STDCALL fglevalpoint2 'fglevalpoint2' (REF,REF)
 STDCALL fglfeedbackbuffer 'fglfeedbackbuffer' (REF,REF,REF)
 STDCALL fglfinish 'fglfinish'
 STDCALL fglflush 'fglflush'
 STDCALL fglfogf 'fglfogf' (REF,REF)
 STDCALL fglfogfv 'fglfogfv' (REF,REF)
 STDCALL fglfogi 'fglfogi' (REF,REF)
 STDCALL fglfogiv 'fglfogiv' (REF,REF)
 STDCALL fglfrontface 'fglfrontface' (REF)
 STDCALL fglfrustum 'fglfrustum' (REF,REF,REF,REF,REF,REF)
 STDCALL fglgenlists 'fglgenlists' (REF)
 STDCALL fglgetbooleanv 'fglgetbooleanv' (REF,REF)
 STDCALL f9eglgetbooleanv 'f9eglgetbooleanv' (REF,REF,REF)
 STDCALL fglgetclipplane 'fglgetclipplane' (REF,REF)
 STDCALL fglgetdoublev 'fglgetdoublev' (REF,REF)
 STDCALL fglgeterror 'fglgeterror'
 STDCALL fglgetfloatv 'fglgetfloatv' (REF,REF)
 STDCALL fglgetintegerv 'fglgetintegerv' (REF,REF)
 STDCALL fglgetlightfv 'fglgetlightfv' (REF,REF,REF)
 STDCALL fglgetlightiv 'fglgetlightiv' (REF,REF,REF)
 STDCALL fglgetmapdv 'fglgetmapdv' (REF,REF,REF)
 STDCALL fglgetmapfv 'fglgetmapfv' (REF,REF,REF)
 STDCALL fglgetmapiv 'fglgetmapiv' (REF,REF,REF)
 STDCALL fglgetmaterialfv 'fglgetmaterialfv' (REF,REF,REF)
 STDCALL fglgetmaterialiv 'fglgetmaterialiv' (REF,REF,REF)
 STDCALL fglgetpixelmapfv 'fglgetpixelmapfv' (REF,REF)
 STDCALL fglgetpixelmapuiv 'fglgetpixelmapuiv' (REF,REF)
 STDCALL fglgetpixelmapusv 'fglgetpixelmapusv' (REF,REF)
 STDCALL f9aglgetpixelmapusv 'f9aglgetpixelmapusv' (REF,REF)
 STDCALL fglgetpolygonstipple 'fglgetpolygonstipple' (REF,REF)
 STDCALL f9bglgetpolygonstipple 'f9bglgetpolygonstipple' (REF,REF)
 STDCALL f9cglgetpolygonstipple 'f9cglgetpolygonstipple' (REF,REF)
 STDCALL f9yglgetstringlen 'f9yglgetstringlen' (REF)
 STDCALL f9yglgetstring 'f9yglgetstring' (REF)
 STDCALL f9yglgetstring1 'f9yglgetstring1' (REF)
 STDCALL f9yglgetstringfree 'f9yglgetstringfree'
 STDCALL fglgettexenvfv 'fglgettexenvfv' (REF,REF,REF)
 STDCALL fglgettexenviv 'fglgettexenviv' (REF,REF,REF)
 STDCALL fglgettexgendv 'fglgettexgendv' (REF,REF,REF)
 STDCALL fglgettexgenfv 'fglgettexgenfv' (REF,REF,REF)
 STDCALL fglgettexgeniv 'fglgettexgeniv' (REF,REF,REF)
 STDCALL f9y0glgetteximage 'f9y0glgetteximage' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glgetteximage 'f9y1glgetteximage' (REF,REF,REF,REF,REF)
 STDCALL f9y4glgetteximage 'f9y4glgetteximage' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y5glgetteximage 'f9y5glgetteximage' (REF,REF,REF,REF,REF)
 STDCALL f9y6glgetteximage 'f9y6glgetteximage' (REF,REF,REF,REF,REF,REF)
 STDCALL fglgettexlevelparameterfv 'fglgettexlevelparameterfv' (REF,REF,REF,REF)
 STDCALL fglgettexlevelparameteriv 'fglgettexlevelparameteriv' (REF,REF,REF,REF)
 STDCALL fglgettexparameterfv 'fglgettexparameterfv' (REF,REF,REF)
 STDCALL fglgettexparameteriv 'fglgettexparameteriv' (REF,REF,REF)
 STDCALL fglhint 'fglhint' (REF,REF)
 STDCALL fglindexmask 'fglindexmask' (REF)
 STDCALL fglindexd 'fglindexd' (REF)
 STDCALL fglindexdv 'fglindexdv' (REF)
 STDCALL fglindexf 'fglindexf' (REF)
 STDCALL fglindexfv 'fglindexfv' (REF)
 STDCALL fglindexi 'fglindexi' (REF)
 STDCALL fglindexiv 'fglindexiv' (REF)
 STDCALL fglindexs 'fglindexs' (REF)
 STDCALL f9aglindexs 'f9aglindexs' (REF)
 STDCALL fglindexsv 'fglindexsv' (REF)
 STDCALL f9aglindexsv 'f9aglindexsv' (REF)
 STDCALL fglinitnames 'fglinitnames'
 STDCALL fglisenabled 'fglisenabled' (REF)
 STDCALL f9eglisenabled 'f9eglisenabled' (REF)
 STDCALL fglislist 'fglislist' (REF)
 STDCALL f9eglislist 'f9eglislist' (REF)
 STDCALL fgllightmodelf 'fgllightmodelf' (REF,REF)
 STDCALL fgllightmodelfv 'fgllightmodelfv' (REF,REF)
 STDCALL fgllightmodeli 'fgllightmodeli' (REF,REF)
 STDCALL fgllightmodeliv 'fgllightmodeliv' (REF,REF)
 STDCALL fgllightf 'fgllightf' (REF,REF,REF)
 STDCALL fgllightfv 'fgllightfv' (REF,REF,REF)
 STDCALL fgllighti 'fgllighti' (REF,REF,REF)
 STDCALL fgllightiv 'fgllightiv' (REF,REF,REF)
 STDCALL fgllinestipple 'fgllinestipple' (REF,REF)
 STDCALL f9agllinestipple 'f9agllinestipple' (REF,REF)
 STDCALL fgllinewidth 'fgllinewidth' (REF)
 STDCALL fgllistbase 'fgllistbase' (REF)
 STDCALL fglloadidentity 'fglloadidentity'
 STDCALL fglloadmatrixd 'fglloadmatrixd' (REF)
 STDCALL fglloadmatrixf 'fglloadmatrixf' (REF)
 STDCALL fglloadname 'fglloadname' (REF)
 STDCALL fgllogicop 'fgllogicop' (REF)
 STDCALL fglmap1d 'fglmap1d' (REF,REF,REF,REF,REF,REF)
 STDCALL fglmap1f 'fglmap1f' (REF,REF,REF,REF,REF,REF)
 STDCALL fglmap2d 'fglmap2d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglmap2f 'fglmap2f' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglmapgrid1d 'fglmapgrid1d' (REF,REF,REF)
 STDCALL fglmapgrid1f 'fglmapgrid1f' (REF,REF,REF)
 STDCALL fglmapgrid2d 'fglmapgrid2d' (REF,REF,REF,REF,REF,REF)
 STDCALL fglmapgrid2f 'fglmapgrid2f' (REF,REF,REF,REF,REF,REF)
 STDCALL fglmaterialf 'fglmaterialf' (REF,REF,REF)
 STDCALL fglmaterialfv 'fglmaterialfv' (REF,REF,REF)
 STDCALL fglmateriali 'fglmateriali' (REF,REF,REF)
 STDCALL fglmaterialiv 'fglmaterialiv' (REF,REF,REF)
 STDCALL fglmatrixmode 'fglmatrixmode' (REF)
 STDCALL fglmultmatrixd 'fglmultmatrixd' (REF)
 STDCALL fglmultmatrixf 'fglmultmatrixf' (REF)
 STDCALL fglnewlist 'fglnewlist' (REF,REF)
 STDCALL fglnormal3b 'fglnormal3b' (REF,REF,REF)
 STDCALL f9bglnormal3b 'f9bglnormal3b' (REF,REF,REF)
 STDCALL f9cglnormal3b 'f9cglnormal3b' (REF,REF,REF)
 STDCALL fglnormal3bv 'fglnormal3bv' (REF)
 STDCALL f9bglnormal3bv 'f9bglnormal3bv' (REF)
 STDCALL f9cglnormal3bv 'f9cglnormal3bv' (REF)
 STDCALL fglnormal3d 'fglnormal3d' (REF,REF,REF)
 STDCALL fglnormal3dv 'fglnormal3dv' (REF)
 STDCALL fglnormal3f 'fglnormal3f' (REF,REF,REF)
 STDCALL fglnormal3fv 'fglnormal3fv' (REF)
 STDCALL fglnormal3i 'fglnormal3i' (REF,REF,REF)
 STDCALL fglnormal3iv 'fglnormal3iv' (REF)
 STDCALL fglnormal3s 'fglnormal3s' (REF,REF,REF)
 STDCALL f9aglnormal3s 'f9aglnormal3s' (REF,REF,REF)
 STDCALL fglnormal3sv 'fglnormal3sv' (REF)
 STDCALL f9aglnormal3sv 'f9aglnormal3sv' (REF)
 STDCALL fglortho 'fglortho' (REF,REF,REF,REF,REF,REF)
 STDCALL fglpassthrough 'fglpassthrough' (REF)
 STDCALL fglpixelmapfv 'fglpixelmapfv' (REF,REF,REF)
 STDCALL fglpixelmapuiv 'fglpixelmapuiv' (REF,REF,REF)
 STDCALL fglpixelmapusv 'fglpixelmapusv' (REF,REF,REF)
 STDCALL f9aglpixelmapusv 'f9aglpixelmapusv' (REF,REF,REF)
 STDCALL fglpixelstoref 'fglpixelstoref' (REF,REF)
 STDCALL fglpixelstorei 'fglpixelstorei' (REF,REF)
 STDCALL fglpixeltransferf 'fglpixeltransferf' (REF,REF)
 STDCALL fglpixeltransferi 'fglpixeltransferi' (REF,REF)
 STDCALL fglpixelzoom 'fglpixelzoom' (REF,REF)
 STDCALL fglpointsize 'fglpointsize' (REF)
 STDCALL fglpolygonmode 'fglpolygonmode' (REF,REF)
 STDCALL fglpolygonstipple 'fglpolygonstipple' (REF)
 STDCALL f9bglpolygonstipple 'f9bglpolygonstipple' (REF)
 STDCALL f9cglpolygonstipple 'f9cglpolygonstipple' (REF)
 STDCALL fglpopattrib 'fglpopattrib'
 STDCALL fglpopmatrix 'fglpopmatrix'
 STDCALL fglpopname 'fglpopname'
 STDCALL fglpushattrib 'fglpushattrib' (REF)
 STDCALL fglpushmatrix 'fglpushmatrix'
 STDCALL fglpushname 'fglpushname' (REF)
 STDCALL fglrasterpos2d 'fglrasterpos2d' (REF,REF)
 STDCALL fglrasterpos2dv 'fglrasterpos2dv' (REF)
 STDCALL fglrasterpos2f 'fglrasterpos2f' (REF,REF)
 STDCALL fglrasterpos2fv 'fglrasterpos2fv' (REF)
 STDCALL fglrasterpos2i 'fglrasterpos2i' (REF,REF)
 STDCALL fglrasterpos2iv 'fglrasterpos2iv' (REF)
 STDCALL fglrasterpos2s 'fglrasterpos2s' (REF,REF)
 STDCALL f9aglrasterpos2s 'f9aglrasterpos2s' (REF,REF)
 STDCALL fglrasterpos2sv 'fglrasterpos2sv' (REF)
 STDCALL f9aglrasterpos2sv 'f9aglrasterpos2sv' (REF)
 STDCALL fglrasterpos3d 'fglrasterpos3d' (REF,REF,REF)
 STDCALL fglrasterpos3dv 'fglrasterpos3dv' (REF)
 STDCALL fglrasterpos3f 'fglrasterpos3f' (REF,REF,REF)
 STDCALL fglrasterpos3fv 'fglrasterpos3fv' (REF)
 STDCALL fglrasterpos3i 'fglrasterpos3i' (REF,REF,REF)
 STDCALL fglrasterpos3iv 'fglrasterpos3iv' (REF)
 STDCALL fglrasterpos3s 'fglrasterpos3s' (REF,REF,REF)
 STDCALL f9aglrasterpos3s 'f9aglrasterpos3s' (REF,REF,REF)
 STDCALL fglrasterpos3sv 'fglrasterpos3sv' (REF)
 STDCALL f9aglrasterpos3sv 'f9aglrasterpos3sv' (REF)
 STDCALL fglrasterpos4d 'fglrasterpos4d' (REF,REF,REF,REF)
 STDCALL fglrasterpos4dv 'fglrasterpos4dv' (REF)
 STDCALL fglrasterpos4f 'fglrasterpos4f' (REF,REF,REF,REF)
 STDCALL fglrasterpos4fv 'fglrasterpos4fv' (REF)
 STDCALL fglrasterpos4i 'fglrasterpos4i' (REF,REF,REF,REF)
 STDCALL fglrasterpos4iv 'fglrasterpos4iv' (REF)
 STDCALL fglrasterpos4s 'fglrasterpos4s' (REF,REF,REF,REF)
 STDCALL f9aglrasterpos4s 'f9aglrasterpos4s' (REF,REF,REF,REF)
 STDCALL fglrasterpos4sv 'fglrasterpos4sv' (REF)
 STDCALL f9aglrasterpos4sv 'f9aglrasterpos4sv' (REF)
 STDCALL fglreadbuffer 'fglreadbuffer' (REF)
 STDCALL f9y0glreadpixels 'f9y0glreadpixels' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glreadpixels 'f9y1glreadpixels' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4glreadpixels 'f9y4glreadpixels' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5glreadpixels 'f9y5glreadpixels' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6glreadpixels 'f9y6glreadpixels' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglrectd 'fglrectd' (REF,REF,REF,REF)
 STDCALL fglrectdv 'fglrectdv' (REF,REF)
 STDCALL fglrectf 'fglrectf' (REF,REF,REF,REF)
 STDCALL fglrectfv 'fglrectfv' (REF,REF)
 STDCALL fglrecti 'fglrecti' (REF,REF,REF,REF)
 STDCALL fglrectiv 'fglrectiv' (REF,REF)
 STDCALL fglrects 'fglrects' (REF,REF,REF,REF)
 STDCALL f9aglrects 'f9aglrects' (REF,REF,REF,REF)
 STDCALL fglrectsv 'fglrectsv' (REF,REF)
 STDCALL f9aglrectsv 'f9aglrectsv' (REF,REF)
 STDCALL fglrendermode 'fglrendermode' (REF)
 STDCALL fglrotated 'fglrotated' (REF,REF,REF,REF)
 STDCALL fglrotatef 'fglrotatef' (REF,REF,REF,REF)
 STDCALL fglscaled 'fglscaled' (REF,REF,REF)
 STDCALL fglscalef 'fglscalef' (REF,REF,REF)
 STDCALL fglscissor 'fglscissor' (REF,REF,REF,REF)
 STDCALL fglselectbuffer 'fglselectbuffer' (REF,REF)
 STDCALL fglshademodel 'fglshademodel' (REF)
 STDCALL fglstencilfunc 'fglstencilfunc' (REF,REF,REF)
 STDCALL fglstencilmask 'fglstencilmask' (REF)
 STDCALL fglstencilop 'fglstencilop' (REF,REF,REF)
 STDCALL fgltexcoord1d 'fgltexcoord1d' (REF)
 STDCALL fgltexcoord1dv 'fgltexcoord1dv' (REF)
 STDCALL fgltexcoord1f 'fgltexcoord1f' (REF)
 STDCALL fgltexcoord1fv 'fgltexcoord1fv' (REF)
 STDCALL fgltexcoord1i 'fgltexcoord1i' (REF)
 STDCALL fgltexcoord1iv 'fgltexcoord1iv' (REF)
 STDCALL fgltexcoord1s 'fgltexcoord1s' (REF)
 STDCALL f9agltexcoord1s 'f9agltexcoord1s' (REF)
 STDCALL fgltexcoord1sv 'fgltexcoord1sv' (REF)
 STDCALL f9agltexcoord1sv 'f9agltexcoord1sv' (REF)
 STDCALL fgltexcoord2d 'fgltexcoord2d' (REF,REF)
 STDCALL fgltexcoord2dv 'fgltexcoord2dv' (REF)
 STDCALL fgltexcoord2f 'fgltexcoord2f' (REF,REF)
 STDCALL fgltexcoord2fv 'fgltexcoord2fv' (REF)
 STDCALL fgltexcoord2i 'fgltexcoord2i' (REF,REF)
 STDCALL fgltexcoord2iv 'fgltexcoord2iv' (REF)
 STDCALL fgltexcoord2s 'fgltexcoord2s' (REF,REF)
 STDCALL f9agltexcoord2s 'f9agltexcoord2s' (REF,REF)
 STDCALL fgltexcoord2sv 'fgltexcoord2sv' (REF)
 STDCALL f9agltexcoord2sv 'f9agltexcoord2sv' (REF)
 STDCALL fgltexcoord3d 'fgltexcoord3d' (REF,REF,REF)
 STDCALL fgltexcoord3dv 'fgltexcoord3dv' (REF)
 STDCALL fgltexcoord3f 'fgltexcoord3f' (REF,REF,REF)
 STDCALL fgltexcoord3fv 'fgltexcoord3fv' (REF)
 STDCALL fgltexcoord3i 'fgltexcoord3i' (REF,REF,REF)
 STDCALL fgltexcoord3iv 'fgltexcoord3iv' (REF)
 STDCALL fgltexcoord3s 'fgltexcoord3s' (REF,REF,REF)
 STDCALL f9agltexcoord3s 'f9agltexcoord3s' (REF,REF,REF)
 STDCALL fgltexcoord3sv 'fgltexcoord3sv' (REF)
 STDCALL f9agltexcoord3sv 'f9agltexcoord3sv' (REF)
 STDCALL fgltexcoord4d 'fgltexcoord4d' (REF,REF,REF,REF)
 STDCALL fgltexcoord4dv 'fgltexcoord4dv' (REF)
 STDCALL fgltexcoord4f 'fgltexcoord4f' (REF,REF,REF,REF)
 STDCALL fgltexcoord4fv 'fgltexcoord4fv' (REF)
 STDCALL fgltexcoord4i 'fgltexcoord4i' (REF,REF,REF,REF)
 STDCALL fgltexcoord4iv 'fgltexcoord4iv' (REF)
 STDCALL fgltexcoord4s 'fgltexcoord4s' (REF,REF,REF,REF)
 STDCALL f9agltexcoord4s 'f9agltexcoord4s' (REF,REF,REF,REF)
 STDCALL fgltexcoord4sv 'fgltexcoord4sv' (REF)
 STDCALL f9agltexcoord4sv 'f9agltexcoord4sv' (REF)
 STDCALL fgltexenvf 'fgltexenvf' (REF,REF,REF)
 STDCALL fgltexenvfv 'fgltexenvfv' (REF,REF,REF)
 STDCALL fgltexenvi 'fgltexenvi' (REF,REF,REF)
 STDCALL fgltexenviv 'fgltexenviv' (REF,REF,REF)
 STDCALL fgltexgend 'fgltexgend' (REF,REF,REF)
 STDCALL fgltexgendv 'fgltexgendv' (REF,REF,REF)
 STDCALL fgltexgenf 'fgltexgenf' (REF,REF,REF)
 STDCALL fgltexgenfv 'fgltexgenfv' (REF,REF,REF)
 STDCALL fgltexgeni 'fgltexgeni' (REF,REF,REF)
 STDCALL fgltexgeniv 'fgltexgeniv' (REF,REF,REF)
 STDCALL f9y0glteximage1d 'f9y0glteximage1d' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glteximage1d 'f9y1glteximage1d' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4glteximage1d 'f9y4glteximage1d' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5glteximage1d 'f9y5glteximage1d' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6glteximage1d 'f9y6glteximage1d' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y0glteximage2d 'f9y0glteximage2d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glteximage2d 'f9y1glteximage2d' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4glteximage2d 'f9y4glteximage2d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5glteximage2d 'f9y5glteximage2d' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6glteximage2d 'f9y6glteximage2d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fgltexparameterf 'fgltexparameterf' (REF,REF,REF)
 STDCALL fgltexparameterfv 'fgltexparameterfv' (REF,REF,REF)
 STDCALL fgltexparameteri 'fgltexparameteri' (REF,REF,REF)
 STDCALL fgltexparameteriv 'fgltexparameteriv' (REF,REF,REF)
 STDCALL fgltranslated 'fgltranslated' (REF,REF,REF)
 STDCALL fgltranslatef 'fgltranslatef' (REF,REF,REF)
 STDCALL fglvertex2d 'fglvertex2d' (REF,REF)
 STDCALL fglvertex2dv 'fglvertex2dv' (REF)
 STDCALL fglvertex2f 'fglvertex2f' (REF,REF)
 STDCALL fglvertex2fv 'fglvertex2fv' (REF)
 STDCALL fglvertex2i 'fglvertex2i' (REF,REF)
 STDCALL fglvertex2iv 'fglvertex2iv' (REF)
 STDCALL fglvertex2s 'fglvertex2s' (REF,REF)
 STDCALL f9aglvertex2s 'f9aglvertex2s' (REF,REF)
 STDCALL fglvertex2sv 'fglvertex2sv' (REF)
 STDCALL f9aglvertex2sv 'f9aglvertex2sv' (REF)
 STDCALL fglvertex3d 'fglvertex3d' (REF,REF,REF)
 STDCALL fglvertex3dv 'fglvertex3dv' (REF)
 STDCALL fglvertex3f 'fglvertex3f' (REF,REF,REF)
 STDCALL fglvertex3fv 'fglvertex3fv' (REF)
 STDCALL fglvertex3i 'fglvertex3i' (REF,REF,REF)
 STDCALL fglvertex3iv 'fglvertex3iv' (REF)
 STDCALL fglvertex3s 'fglvertex3s' (REF,REF,REF)
 STDCALL f9aglvertex3s 'f9aglvertex3s' (REF,REF,REF)
 STDCALL fglvertex3sv 'fglvertex3sv' (REF)
 STDCALL f9aglvertex3sv 'f9aglvertex3sv' (REF)
 STDCALL fglvertex4d 'fglvertex4d' (REF,REF,REF,REF)
 STDCALL fglvertex4dv 'fglvertex4dv' (REF)
 STDCALL fglvertex4f 'fglvertex4f' (REF,REF,REF,REF)
 STDCALL fglvertex4fv 'fglvertex4fv' (REF)
 STDCALL fglvertex4i 'fglvertex4i' (REF,REF,REF,REF)
 STDCALL fglvertex4iv 'fglvertex4iv' (REF)
 STDCALL fglvertex4s 'fglvertex4s' (REF,REF,REF,REF)
 STDCALL f9aglvertex4s 'f9aglvertex4s' (REF,REF,REF,REF)
 STDCALL fglvertex4sv 'fglvertex4sv' (REF)
 STDCALL f9aglvertex4sv 'f9aglvertex4sv' (REF)
 STDCALL fglviewport 'fglviewport' (REF,REF,REF,REF)
#ifdef OPENGL_1_1
 STDCALL fglaretexturesresident 'fglaretexturesresident' (REF,REF,REF)
 STDCALL f9eglaretexturesresident 'f9eglaretexturesresident' (REF,REF,REF)
 STDCALL fglarrayelement 'fglarrayelement' (REF)
 STDCALL fglbindtexture 'fglbindtexture' (REF,REF)
 STDCALL f9y0glcolorpointer 'f9y0glcolorpointer' (REF,REF,REF,REF,REF)
 STDCALL f9y1glcolorpointer 'f9y1glcolorpointer' (REF,REF,REF,REF)
 STDCALL f9y4glcolorpointer 'f9y4glcolorpointer' (REF,REF,REF,REF,REF)
 STDCALL f9y5glcolorpointer 'f9y5glcolorpointer' (REF,REF,REF,REF)
 STDCALL f9y6glcolorpointer 'f9y6glcolorpointer' (REF,REF,REF,REF,REF)
 STDCALL f9y7glcolorpointer 'f9y7glcolorpointer' (REF,REF,REF,REF)
 STDCALL fglcopyteximage1d 'fglcopyteximage1d' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglcopyteximage2d 'fglcopyteximage2d' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglcopytexsubimage1d 'fglcopytexsubimage1d' (REF,REF,REF,REF,REF,REF)
 STDCALL fglcopytexsubimage2d 'fglcopytexsubimage2d' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fgldeletetextures 'fgldeletetextures' (REF,REF)
 STDCALL fgldisableclientstate 'fgldisableclientstate' (REF)
 STDCALL fgldrawarrays 'fgldrawarrays' (REF,REF,REF)
 STDCALL f9y0gldrawelements 'f9y0gldrawelements' (REF,REF,REF,REF)
 STDCALL f9y4gldrawelements 'f9y4gldrawelements' (REF,REF,REF,REF,REF)
 STDCALL f9y5gldrawelements 'f9y5gldrawelements' (REF,REF,REF,REF)
 STDCALL f9y6gldrawelements 'f9y6gldrawelements' (REF,REF,REF,REF)
 STDCALL fgledgeflagpointer 'fgledgeflagpointer' (REF,REF)
 STDCALL f9egledgeflagpointer 'f9egledgeflagpointer' (REF,REF,REF)
 STDCALL fglenableclientstate 'fglenableclientstate' (REF)
 STDCALL fglgentextures 'fglgentextures' (REF,REF)
 STDCALL fglgetpointerv 'fglgetpointerv' (REF,REF,REF)
 STDCALL f9y0glindexpointer 'f9y0glindexpointer' (REF,REF,REF,REF)
 STDCALL f9y1glindexpointer 'f9y1glindexpointer' (REF,REF,REF)
 STDCALL f9y4glindexpointer 'f9y4glindexpointer' (REF,REF,REF,REF)
 STDCALL f9y6glindexpointer 'f9y6glindexpointer' (REF,REF,REF)
 STDCALL f9y7glindexpointer 'f9y7glindexpointer' (REF,REF,REF)
 STDCALL fglindexub 'fglindexub' (REF)
 STDCALL f9bglindexub 'f9bglindexub' (REF)
 STDCALL f9cglindexub 'f9cglindexub'
 STDCALL fglindexubv 'fglindexubv' (REF)
 STDCALL f9bglindexubv 'f9bglindexubv' (REF)
 STDCALL f9cglindexubv 'f9cglindexubv' (REF)
 STDCALL f9y0glinterleavedarrays 'f9y0glinterleavedarrays' (REF,REF,REF,REF)
 STDCALL f9y4glinterleavedarrays 'f9y4glinterleavedarrays' (REF,REF,REF,REF,REF)
 STDCALL f9y1glinterleavedarrays 'f9y1glinterleavedarrays' (REF,REF,REF)
 STDCALL f9y5glinterleavedarrays 'f9y5glinterleavedarrays' (REF,REF,REF)
 STDCALL f9y6glinterleavedarrays 'f9y6glinterleavedarrays' (REF,REF,REF,REF)
 STDCALL fglistexture 'fglistexture' (REF)
 STDCALL f9eglistexture 'f9eglistexture' (REF)
 STDCALL f9y0glnormalpointer 'f9y0glnormalpointer' (REF,REF,REF,REF)
 STDCALL f9y1glnormalpointer 'f9y1glnormalpointer' (REF,REF,REF)
 STDCALL f9y4glnormalpointer 'f9y4glnormalpointer' (REF,REF,REF,REF)
 STDCALL f9y5glnormalpointer 'f9y5glnormalpointer' (REF,REF,REF)
 STDCALL f9y6glnormalpointer 'f9y6glnormalpointer' (REF,REF,REF,REF)
 STDCALL f9y7glnormalpointer 'f9y7glnormalpointer' (REF,REF,REF)
 STDCALL fglpolygonoffset 'fglpolygonoffset' (REF,REF)
 STDCALL fglpopclientattrib 'fglpopclientattrib'
 STDCALL fglprioritizetextures 'fglprioritizetextures' (REF,REF,REF)
 STDCALL fglpushclientattrib 'fglpushclientattrib' (REF)
 STDCALL f9y0gltexcoordpointer 'f9y0gltexcoordpointer' (REF,REF,REF,REF,REF)
 STDCALL f9y1gltexcoordpointer 'f9y1gltexcoordpointer' (REF,REF,REF,REF)
 STDCALL f9y4gltexcoordpointer 'f9y4gltexcoordpointer' (REF,REF,REF,REF,REF)
 STDCALL f9y6gltexcoordpointer 'f9y6gltexcoordpointer' (REF,REF,REF,REF)
 STDCALL f9y7gltexcoordpointer 'f9y7gltexcoordpointer' (REF,REF,REF,REF)
 STDCALL f9y0gltexsubimage1d 'f9y0gltexsubimage1d' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1gltexsubimage1d 'f9y1gltexsubimage1d' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4gltexsubimage1d 'f9y4gltexsubimage1d' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5gltexsubimage1d 'f9y5gltexsubimage1d' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6gltexsubimage1d 'f9y6gltexsubimage1d' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y0gltexsubimage2d 'f9y0gltexsubimage2d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1gltexsubimage2d 'f9y1gltexsubimage2d' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4gltexsubimage2d 'f9y4gltexsubimage2d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5gltexsubimage2d 'f9y5gltexsubimage2d' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6gltexsubimage2d 'f9y6gltexsubimage2d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y0glvertexpointer 'f9y0glvertexpointer' (REF,REF,REF,REF,REF)
 STDCALL f9y1glvertexpointer 'f9y1glvertexpointer' (REF,REF,REF,REF)
 STDCALL f9y4glvertexpointer 'f9y4glvertexpointer' (REF,REF,REF,REF,REF)
 STDCALL f9y6glvertexpointer 'f9y6glvertexpointer' (REF,REF,REF,REF)
 STDCALL f9y7glvertexpointer 'f9y7glvertexpointer' (REF,REF,REF,REF)
#endif
#ifdef OPENGL_1_2
 STDCALL f9y0gldrawrangeelements 'f9y0gldrawrangeelements' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y4gldrawrangeelements 'f9y4gldrawrangeelements' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5gldrawrangeelements 'f9y5gldrawrangeelements' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y6gldrawrangeelements 'f9y6gldrawrangeelements' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y0glteximage3d 'f9y0glteximage3d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glteximage3d 'f9y1glteximage3d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4glteximage3d 'f9y4glteximage3d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5glteximage3d 'f9y5glteximage3d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6glteximage3d 'f9y6glteximage3d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y0gltexsubimage3d 'f9y0gltexsubimage3d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1gltexsubimage3d 'f9y1gltexsubimage3d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4gltexsubimage3d 'f9y4gltexsubimage3d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5gltexsubimage3d 'f9y5gltexsubimage3d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6gltexsubimage3d 'f9y6gltexsubimage3d' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglcopytexsubimage3d 'fglcopytexsubimage3d' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
#endif
#ifdef GL_EXT_blend_color
 STDCALL fglblendcolorext 'fglblendcolorext' (REF,REF,REF,REF)
#endif
#ifdef GL_EXT_blend_minmax
 STDCALL fglblendequationext 'fglblendequationext' (REF)
#endif
#ifdef GL_EXT_color_table
 STDCALL f9y0glcolortableext 'f9y0glcolortableext' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glcolortableext 'f9y1glcolortableext' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y4glcolortableext 'f9y4glcolortableext' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5glcolortableext 'f9y5glcolortableext' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y6glcolortableext 'f9y6glcolortableext' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y0glcolorsubtableext 'f9y0glcolorsubtableext' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glcolorsubtableext 'f9y1glcolorsubtableext' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y4glcolorsubtableext 'f9y4glcolorsubtableext' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5glcolorsubtableext 'f9y5glcolorsubtableext' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y6glcolorsubtableext 'f9y6glcolorsubtableext' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y0glgetcolortableext 'f9y0glgetcolortableext' (REF,REF,REF,REF,REF)
 STDCALL f9y4glgetcolortableext 'f9y4glgetcolortableext' (REF,REF,REF,REF,REF)
 STDCALL f9y1glgetcolortableext 'f9y1glgetcolortableext' (REF,REF,REF,REF)
 STDCALL f9y5glgetcolortableext 'f9y5glgetcolortableext' (REF,REF,REF,REF)
 STDCALL f9y6glgetcolortableext 'f9y6glgetcolortableext' (REF,REF,REF,REF,REF)
 STDCALL fglgetcolortableparameterfvext 'fglgetcolortableparameterfvext' (REF,REF,REF)
 STDCALL fglgetcolortableparameterivext 'fglgetcolortableparameterivext' (REF,REF,REF)
#endif
#ifdef GL_EXT_point_parameters
 STDCALL fglpointparameterfext 'fglpointparameterfext' (REF,REF)
 STDCALL fglpointparameterfvext 'fglpointparameterfvext' (REF,REF)
#endif
#ifdef GL_EXT_polygon_offset
 STDCALL fglpolygonoffsetext 'fglpolygonoffsetext' (REF,REF)
#endif
#ifdef GL_EXT_subtexture
 STDCALL f9y0gltexsubimage1dext 'f9y0gltexsubimage1dext' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1gltexsubimage1dext 'f9y1gltexsubimage1dext' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4gltexsubimage1dext 'f9y4gltexsubimage1dext' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5gltexsubimage1dext 'f9y5gltexsubimage1dext' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6gltexsubimage1dext 'f9y6gltexsubimage1dext' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y0gltexsubimage2dext 'f9y0gltexsubimage2dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1gltexsubimage2dext 'f9y1gltexsubimage2dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4gltexsubimage2dext 'f9y4gltexsubimage2dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5gltexsubimage2dext 'f9y5gltexsubimage2dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6gltexsubimage2dext 'f9y6gltexsubimage2dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
#endif
#ifdef GL_EXT_texture3D
 STDCALL f9y0glteximage3dext 'f9y0glteximage3dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glteximage3dext 'f9y1glteximage3dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4glteximage3dext 'f9y4glteximage3dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5glteximage3dext 'f9y5glteximage3dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6glteximage3dext 'f9y6glteximage3dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y0gltexsubimage3dext 'f9y0gltexsubimage3dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1gltexsubimage3dext 'f9y1gltexsubimage3dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4gltexsubimage3dext 'f9y4gltexsubimage3dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5gltexsubimage3dext 'f9y5gltexsubimage3dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6gltexsubimage3dext 'f9y6gltexsubimage3dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglcopytexsubimage3dext 'fglcopytexsubimage3dext' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
#endif
#ifdef GL_EXT_texture_object
 STDCALL fglgentexturesext 'fglgentexturesext' (REF,REF)
 STDCALL fgldeletetexturesext 'fgldeletetexturesext' (REF,REF)
 STDCALL fglbindtextureext 'fglbindtextureext' (REF,REF)
 STDCALL fglprioritizetexturesext 'fglprioritizetexturesext' (REF,REF,REF)
 STDCALL fglaretexturesresidentext 'fglaretexturesresidentext' (REF,REF,REF)
 STDCALL f9eglaretexturesresidentext 'f9eglaretexturesresidentext' (REF,REF,REF)
 STDCALL fglistextureext 'fglistextureext' (REF)
 STDCALL f9eglistextureext 'f9eglistextureext' (REF)
#endif
#ifdef GL_EXT_vertex_array
 STDCALL f9y0glvertexpointerext 'f9y0glvertexpointerext' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glvertexpointerext 'f9y1glvertexpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y4glvertexpointerext 'f9y4glvertexpointerext' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y6glvertexpointerext 'f9y6glvertexpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y7glvertexpointerext 'f9y7glvertexpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y0glnormalpointerext 'f9y0glnormalpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y1glnormalpointerext 'f9y1glnormalpointerext' (REF,REF,REF,REF)
 STDCALL f9y4glnormalpointerext 'f9y4glnormalpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y5glnormalpointerext 'f9y5glnormalpointerext' (REF,REF,REF,REF)
 STDCALL f9y6glnormalpointerext 'f9y6glnormalpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y7glnormalpointerext 'f9y7glnormalpointerext' (REF,REF,REF,REF)
 STDCALL f9y0glcolorpointerext 'f9y0glcolorpointerext' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glcolorpointerext 'f9y1glcolorpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y4glcolorpointerext 'f9y4glcolorpointerext' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y5glcolorpointerext 'f9y5glcolorpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y6glcolorpointerext 'f9y6glcolorpointerext' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y7glcolorpointerext 'f9y7glcolorpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y0glindexpointerext 'f9y0glindexpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y1glindexpointerext 'f9y1glindexpointerext' (REF,REF,REF,REF)
 STDCALL f9y4glindexpointerext 'f9y4glindexpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y6glindexpointerext 'f9y6glindexpointerext' (REF,REF,REF,REF)
 STDCALL f9y7glindexpointerext 'f9y7glindexpointerext' (REF,REF,REF,REF)
 STDCALL f9y0gltexcoordpointerext 'f9y0gltexcoordpointerext' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y1gltexcoordpointerext 'f9y1gltexcoordpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y4gltexcoordpointerext 'f9y4gltexcoordpointerext' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y6gltexcoordpointerext 'f9y6gltexcoordpointerext' (REF,REF,REF,REF,REF)
 STDCALL f9y7gltexcoordpointerext 'f9y7gltexcoordpointerext' (REF,REF,REF,REF,REF)
 STDCALL fgledgeflagpointerext 'fgledgeflagpointerext' (REF,REF,REF)
 STDCALL f9egledgeflagpointerext 'f9egledgeflagpointerext' (REF,REF,REF,REF)
 STDCALL fglgetpointervext 'fglgetpointervext' (REF,REF,REF)
 STDCALL fglarrayelementext 'fglarrayelementext' (REF)
 STDCALL fgldrawarraysext 'fgldrawarraysext' (REF,REF,REF)
#endif
#ifdef GL_SGIS_multisample
 STDCALL fglsamplemasksgis 'fglsamplemasksgis' (REF,REF)
 STDCALL f9eglsamplemasksgis 'f9eglsamplemasksgis' (REF,REF)
 STDCALL fglsamplepatternsgis 'fglsamplepatternsgis' (REF,REF)
#endif
#ifdef GL_MESA_resize_buffers
 STDCALL fglresizebuffersmesa 'fglresizebuffersmesa'
#endif
#ifdef GL_MESA_window_pos
 STDCALL fglwindowpos2imesa 'fglwindowpos2imesa' (REF,REF)
 STDCALL fglwindowpos2smesa 'fglwindowpos2smesa' (REF,REF)
 STDCALL f9aglwindowpos2smesa 'f9aglwindowpos2smesa' (REF,REF)
 STDCALL fglwindowpos2fmesa 'fglwindowpos2fmesa' (REF,REF)
 STDCALL fglwindowpos2dmesa 'fglwindowpos2dmesa' (REF,REF)
 STDCALL fglwindowpos2ivmesa 'fglwindowpos2ivmesa' (REF)
 STDCALL fglwindowpos2svmesa 'fglwindowpos2svmesa' (REF)
 STDCALL f9aglwindowpos2svmesa 'f9aglwindowpos2svmesa' (REF)
 STDCALL fglwindowpos2fvmesa 'fglwindowpos2fvmesa' (REF)
 STDCALL fglwindowpos2dvmesa 'fglwindowpos2dvmesa' (REF)
 STDCALL fglwindowpos3imesa 'fglwindowpos3imesa' (REF,REF,REF)
 STDCALL fglwindowpos3smesa 'fglwindowpos3smesa' (REF,REF,REF)
 STDCALL f9aglwindowpos3smesa 'f9aglwindowpos3smesa' (REF,REF,REF)
 STDCALL fglwindowpos3fmesa 'fglwindowpos3fmesa' (REF,REF,REF)
 STDCALL fglwindowpos3dmesa 'fglwindowpos3dmesa' (REF,REF,REF)
 STDCALL fglwindowpos3ivmesa 'fglwindowpos3ivmesa' (REF)
 STDCALL fglwindowpos3svmesa 'fglwindowpos3svmesa' (REF)
 STDCALL f9aglwindowpos3svmesa 'f9aglwindowpos3svmesa' (REF)
 STDCALL fglwindowpos3fvmesa 'fglwindowpos3fvmesa' (REF)
 STDCALL fglwindowpos3dvmesa 'fglwindowpos3dvmesa' (REF)
 STDCALL fglwindowpos4imesa 'fglwindowpos4imesa' (REF,REF,REF,REF)
 STDCALL fglwindowpos4smesa 'fglwindowpos4smesa' (REF,REF,REF,REF)
 STDCALL f9aglwindowpos4smesa 'f9aglwindowpos4smesa' (REF,REF,REF,REF)
 STDCALL fglwindowpos4fmesa 'fglwindowpos4fmesa' (REF,REF,REF,REF)
 STDCALL fglwindowpos4dmesa 'fglwindowpos4dmesa' (REF,REF,REF,REF)
 STDCALL fglwindowpos4ivmesa 'fglwindowpos4ivmesa' (REF)
 STDCALL fglwindowpos4svmesa 'fglwindowpos4svmesa' (REF)
 STDCALL f9aglwindowpos4svmesa 'f9aglwindowpos4svmesa' (REF)
 STDCALL fglwindowpos4fvmesa 'fglwindowpos4fvmesa' (REF)
 STDCALL fglwindowpos4dvmesa 'fglwindowpos4dvmesa' (REF)
#endif
#endif

private

!--------------------------
! C wrapper procedure names
!--------------------------

! OpenGL 1.0

public :: &
   fglaccum, &
   fglalphafunc, &
   fglbegin, &
   fglbitmap, &
   f9bglbitmap, &
   f9cglbitmap, &
   fglblendfunc, &
   fglcalllist, &
   f9y0glcalllists, &
   f9y1glcalllists, &
   f9y4glcalllists, &
   f9y5glcalllists, &
   f9y6glcalllists, &
   fglclear, &
   fglclearaccum, &
   fglclearcolor, &
   fglcleardepth, &
   fglclearindex, &
   fglclearstencil, &
   fglclipplane, &
   fglcolor3b, &
   f9bglcolor3b, &
   f9cglcolor3b, &
   fglcolor3bv, &
   f9bglcolor3bv, &
   f9cglcolor3bv, &
   fglcolor3d, &
   fglcolor3dv, &
   fglcolor3f, &
   fglcolor3fv, &
   fglcolor3i, &
   fglcolor3iv, &
   fglcolor3s, &
   f9aglcolor3s, &
   fglcolor3sv, &
   f9aglcolor3sv
public :: &
   fglcolor3ub, &
   f9bglcolor3ub, &
   f9cglcolor3ub, &
   fglcolor3ubv, &
   f9bglcolor3ubv, &
   f9cglcolor3ubv, &
   fglcolor3ui, &
   fglcolor3uiv, &
   fglcolor3us, &
   f9aglcolor3us, &
   fglcolor3usv, &
   f9aglcolor3usv, &
   fglcolor4b, &
   f9bglcolor4b, &
   f9cglcolor4b, &
   fglcolor4bv, &
   f9bglcolor4bv, &
   f9cglcolor4bv, &
   fglcolor4d, &
   fglcolor4dv, &
   fglcolor4f, &
   fglcolor4fv, &
   fglcolor4i, &
   fglcolor4iv, &
   fglcolor4s, &
   f9aglcolor4s, &
   fglcolor4sv, &
   f9aglcolor4sv, &
   fglcolor4ub, &
   f9bglcolor4ub, &
   f9cglcolor4ub, &
   fglcolor4ubv, &
   f9bglcolor4ubv, &
   f9cglcolor4ubv, &
   fglcolor4ui, &
   fglcolor4uiv
public :: &
   fglcolor4us, &
   f9aglcolor4us, &
   fglcolor4usv, &
   f9aglcolor4usv, &
   fglcolormask, &
   f9eglcolormask, &
   fglcolormaterial, &
   fglcopypixels, &
   fglcullface, &
   fgldeletelists, &
   fgldepthfunc, &
   fgldepthmask, &
   f9egldepthmask, &
   fgldepthrange, &
   fgldisable, &
   fgldrawbuffer, &
   f9y0gldrawpixels, &
   f9y1gldrawpixels, &
   f9y4gldrawpixels, &
   f9y5gldrawpixels, &
   f9y6gldrawpixels, &
   fgledgeflag, &
   f9egledgeflag, &
   fgledgeflagv, &
   f9egledgeflagv, &
   fglenable, &
   fglend, &
   fglendlist, &
   fglevalcoord1d, &
   fglevalcoord1dv, &
   fglevalcoord1f, &
   fglevalcoord1fv, &
   fglevalcoord2d, &
   fglevalcoord2dv
public :: &
   fglevalcoord2f, &
   fglevalcoord2fv, &
   fglevalmesh1, &
   fglevalmesh2, &
   fglevalpoint1, &
   fglevalpoint2, &
   fglfeedbackbuffer, &
   fglfinish, &
   fglflush, &
   fglfogf, &
   fglfogfv, &
   fglfogi, &
   fglfogiv, &
   fglfrontface, &
   fglfrustum, &
   fglgenlists, &
   fglgetbooleanv, &
   f9eglgetbooleanv, &
   fglgetclipplane, &
   fglgetdoublev, &
   fglgeterror, &
   fglgetfloatv, &
   fglgetintegerv, &
   fglgetlightfv, &
   fglgetlightiv, &
   fglgetmapdv, &
   fglgetmapfv, &
   fglgetmapiv, &
   fglgetmaterialfv, &
   fglgetmaterialiv, &
   fglgetpixelmapfv, &
   fglgetpixelmapuiv, &
   fglgetpixelmapusv, &
   f9aglgetpixelmapusv
public :: &
   fglgetpolygonstipple, &
   f9bglgetpolygonstipple, &
   f9cglgetpolygonstipple, &
   f9yglgetstringlen, &
   f9yglgetstring, &
   f9yglgetstring1, &
   f9yglgetstringfree, &
   fglgettexenvfv, &
   fglgettexenviv, &
   fglgettexgendv, &
   fglgettexgenfv, &
   fglgettexgeniv, &
   f9y0glgetteximage, &
   f9y1glgetteximage, &
   f9y4glgetteximage, &
   f9y5glgetteximage, &
   f9y6glgetteximage, &
   fglgettexlevelparameterfv, &
   fglgettexlevelparameteriv, &
   fglgettexparameterfv, &
   fglgettexparameteriv, &
   fglhint, &
   fglindexmask, &
   fglindexd, &
   fglindexdv, &
   fglindexf, &
   fglindexfv, &
   fglindexi, &
   fglindexiv, &
   fglindexs, &
   f9aglindexs, &
   fglindexsv, &
   f9aglindexsv, &
   fglinitnames, &
   fglisenabled, &
   f9eglisenabled
public :: &
   fglislist, &
   f9eglislist, &
   fgllightmodelf, &
   fgllightmodelfv, &
   fgllightmodeli, &
   fgllightmodeliv, &
   fgllightf, &
   fgllightfv, &
   fgllighti, &
   fgllightiv, &
   fgllinestipple, &
   f9agllinestipple, &
   fgllinewidth, &
   fgllistbase, &
   fglloadidentity, &
   fglloadmatrixd, &
   fglloadmatrixf, &
   fglloadname, &
   fgllogicop, &
   fglmap1d, &
   fglmap1f, &
   fglmap2d, &
   fglmap2f, &
   fglmapgrid1d, &
   fglmapgrid1f, &
   fglmapgrid2d, &
   fglmapgrid2f, &
   fglmaterialf, &
   fglmaterialfv, &
   fglmateriali, &
   fglmaterialiv, &
   fglmatrixmode, &
   fglmultmatrixd, &
   fglmultmatrixf
public :: &
   fglnewlist, &
   fglnormal3b, &
   f9bglnormal3b, &
   f9cglnormal3b, &
   fglnormal3bv, &
   f9bglnormal3bv, &
   f9cglnormal3bv, &
   fglnormal3d, &
   fglnormal3dv, &
   fglnormal3f, &
   fglnormal3fv, &
   fglnormal3i, &
   fglnormal3iv, &
   fglnormal3s, &
   f9aglnormal3s, &
   fglnormal3sv, &
   f9aglnormal3sv, &
   fglortho, &
   fglpassthrough, &
   fglpixelmapfv, &
   fglpixelmapuiv, &
   fglpixelmapusv, &
   f9aglpixelmapusv, &
   fglpixelstoref, &
   fglpixelstorei, &
   fglpixeltransferf, &
   fglpixeltransferi, &
   fglpixelzoom, &
   fglpointsize, &
   fglpolygonmode, &
   fglpolygonstipple, &
   f9bglpolygonstipple, &
   f9cglpolygonstipple, &
   fglpopattrib
public :: &
   fglpopmatrix, &
   fglpopname, &
   fglpushattrib, &
   fglpushmatrix, &
   fglpushname, &
   fglrasterpos2d, &
   fglrasterpos2dv, &
   fglrasterpos2f, &
   fglrasterpos2fv, &
   fglrasterpos2i, &
   fglrasterpos2iv, &
   fglrasterpos2s, &
   f9aglrasterpos2s, &
   fglrasterpos2sv, &
   f9aglrasterpos2sv, &
   fglrasterpos3d, &
   fglrasterpos3dv, &
   fglrasterpos3f, &
   fglrasterpos3fv, &
   fglrasterpos3i, &
   fglrasterpos3iv, &
   fglrasterpos3s, &
   f9aglrasterpos3s, &
   fglrasterpos3sv, &
   f9aglrasterpos3sv, &
   fglrasterpos4d, &
   fglrasterpos4dv, &
   fglrasterpos4f, &
   fglrasterpos4fv, &
   fglrasterpos4i, &
   fglrasterpos4iv, &
   fglrasterpos4s, &
   f9aglrasterpos4s, &
   fglrasterpos4sv, &
   f9aglrasterpos4sv
public :: &
   fglreadbuffer, &
   f9y0glreadpixels, &
   f9y1glreadpixels, &
   f9y4glreadpixels, &
   f9y5glreadpixels, &
   f9y6glreadpixels, &
   fglrectd, &
   fglrectdv, &
   fglrectf, &
   fglrectfv, &
   fglrecti, &
   fglrectiv, &
   fglrects, &
   f9aglrects, &
   fglrectsv, &
   f9aglrectsv, &
   fglrendermode, &
   fglrotated, &
   fglrotatef, &
   fglscaled, &
   fglscalef, &
   fglscissor, &
   fglselectbuffer, &
   fglshademodel, &
   fglstencilfunc, &
   fglstencilmask, &
   fglstencilop, &
   fgltexcoord1d, &
   fgltexcoord1dv, &
   fgltexcoord1f
public :: &
   fgltexcoord1fv, &
   fgltexcoord1i, &
   fgltexcoord1iv, &
   fgltexcoord1s, &
   f9agltexcoord1s, &
   fgltexcoord1sv, &
   f9agltexcoord1sv, &
   fgltexcoord2d, &
   fgltexcoord2dv, &
   fgltexcoord2f, &
   fgltexcoord2fv, &
   fgltexcoord2i, &
   fgltexcoord2iv, &
   fgltexcoord2s, &
   f9agltexcoord2s, &
   fgltexcoord2sv, &
   f9agltexcoord2sv, &
   fgltexcoord3d, &
   fgltexcoord3dv, &
   fgltexcoord3f, &
   fgltexcoord3fv, &
   fgltexcoord3i, &
   fgltexcoord3iv, &
   fgltexcoord3s, &
   f9agltexcoord3s, &
   fgltexcoord3sv, &
   f9agltexcoord3sv, &
   fgltexcoord4d, &
   fgltexcoord4dv, &
   fgltexcoord4f, &
   fgltexcoord4fv, &
   fgltexcoord4i, &
   fgltexcoord4iv, &
   fgltexcoord4s, &
   f9agltexcoord4s, &
   fgltexcoord4sv
public :: &
   f9agltexcoord4sv, &
   fgltexenvf, &
   fgltexenvfv, &
   fgltexenvi, &
   fgltexenviv, &
   fgltexgend, &
   fgltexgendv, &
   fgltexgenf, &
   fgltexgenfv, &
   fgltexgeni, &
   fgltexgeniv, &
   f9y0glteximage1d, &
   f9y1glteximage1d, &
   f9y4glteximage1d, &
   f9y5glteximage1d, &
   f9y6glteximage1d, &
   f9y0glteximage2d, &
   f9y1glteximage2d, &
   f9y4glteximage2d, &
   f9y5glteximage2d, &
   f9y6glteximage2d, &
   fgltexparameterf, &
   fgltexparameterfv, &
   fgltexparameteri, &
   fgltexparameteriv, &
   fgltranslated, &
   fgltranslatef, &
   fglvertex2d, &
   fglvertex2dv, &
   fglvertex2f, &
   fglvertex2fv, &
   fglvertex2i, &
   fglvertex2iv, &
   fglvertex2s, &
   f9aglvertex2s, &
   fglvertex2sv, &
   f9aglvertex2sv
public :: &
   fglvertex3d, &
   fglvertex3dv, &
   fglvertex3f, &
   fglvertex3fv, &
   fglvertex3i, &
   fglvertex3iv, &
   fglvertex3s, &
   f9aglvertex3s, &
   fglvertex3sv, &
   f9aglvertex3sv, &
   fglvertex4d, &
   fglvertex4dv, &
   fglvertex4f, &
   fglvertex4fv, &
   fglvertex4i, &
   fglvertex4iv, &
   fglvertex4s, &
   f9aglvertex4s, &
   fglvertex4sv, &
   f9aglvertex4sv, &
   fglviewport

! OpenGL 1.1

#ifdef OPENGL_1_1
public :: &
   fglaretexturesresident, &
   f9eglaretexturesresident, &
   fglarrayelement, &
   fglbindtexture, &
   f9y0glcolorpointer, &
   f9y1glcolorpointer, &
   f9y4glcolorpointer, &
   f9y5glcolorpointer, &
   f9y6glcolorpointer, &
   f9y7glcolorpointer, &
   fglcopyteximage1d, &
   fglcopyteximage2d, &
   fglcopytexsubimage1d, &
   fglcopytexsubimage2d, &
   fgldeletetextures, &
   fgldisableclientstate, &
   fgldrawarrays, &
   f9y0gldrawelements, &
   f9y4gldrawelements, &
   f9y5gldrawelements, &
   f9y6gldrawelements, &
   fgledgeflagpointer, &
   f9egledgeflagpointer, &
   fglenableclientstate, &
   fglgentextures, &
   fglgetpointerv, &
   f9y0glindexpointer, &
   f9y1glindexpointer, &
   f9y4glindexpointer, &
   f9y6glindexpointer, &
   f9y7glindexpointer, &
   fglindexub, &
   f9bglindexub, &
   f9cglindexub, &
   fglindexubv, &
   f9bglindexubv
public :: &
   f9cglindexubv, &
   f9y0glinterleavedarrays, &
   f9y4glinterleavedarrays, &
   f9y1glinterleavedarrays, &
   f9y5glinterleavedarrays, &
   f9y6glinterleavedarrays, &
   fglistexture, &
   f9eglistexture, &
   f9y0glnormalpointer, &
   f9y1glnormalpointer, &
   f9y4glnormalpointer, &
   f9y5glnormalpointer, &
   f9y6glnormalpointer, &
   f9y7glnormalpointer, &
   fglpolygonoffset, &
   fglpopclientattrib, &
   fglprioritizetextures, &
   fglpushclientattrib, &
   f9y0gltexcoordpointer, &
   f9y1gltexcoordpointer, &
   f9y4gltexcoordpointer, &
   f9y6gltexcoordpointer, &
   f9y7gltexcoordpointer, &
   f9y0gltexsubimage1d, &
   f9y1gltexsubimage1d, &
   f9y4gltexsubimage1d, &
   f9y5gltexsubimage1d, &
   f9y6gltexsubimage1d, &
   f9y0gltexsubimage2d, &
   f9y1gltexsubimage2d, &
   f9y4gltexsubimage2d
public :: &
   f9y5gltexsubimage2d, &
   f9y6gltexsubimage2d, &
   f9y0glvertexpointer, &
   f9y1glvertexpointer, &
   f9y4glvertexpointer, &
   f9y6glvertexpointer, &
   f9y7glvertexpointer
#endif

! OpenGL 1.2

#ifdef OPENGL_1_2
public :: &
   f9y0gldrawrangeelements, &
   f9y4gldrawrangeelements, &
   f9y5gldrawrangeelements, &
   f9y6gldrawrangeelements, &
   f9y0glteximage3d, &
   f9y1glteximage3d, &
   f9y4glteximage3d, &
   f9y5glteximage3d, &
   f9y6glteximage3d, &
   f9y0gltexsubimage3d, &
   f9y1gltexsubimage3d, &
   f9y4gltexsubimage3d, &
   f9y5gltexsubimage3d, &
   f9y6gltexsubimage3d, &
   fglcopytexsubimage3d
#endif

! EXT_blend_color

#ifdef GL_EXT_blend_color
public :: &
   fglblendcolorext
#endif

! EXT_blend_minmax

#ifdef GL_EXT_blend_minmax
public :: &
   fglblendequationext
#endif

! EXT_color_table

#ifdef GL_EXT_color_table
public :: &
   f9y0glcolortableext, &
   f9y1glcolortableext, &
   f9y4glcolortableext, &
   f9y5glcolortableext, &
   f9y6glcolortableext, &
   f9y0glcolorsubtableext, &
   f9y1glcolorsubtableext, &
   f9y4glcolorsubtableext, &
   f9y5glcolorsubtableext, &
   f9y6glcolorsubtableext, &
   f9y0glgetcolortableext, &
   f9y4glgetcolortableext, &
   f9y1glgetcolortableext, &
   f9y5glgetcolortableext, &
   f9y6glgetcolortableext, &
   fglgetcolortableparameterfvext, &
   fglgetcolortableparameterivext
#endif

! EXT_point_parameters

#ifdef GL_EXT_point_parameters
public :: &
   fglpointparameterfext, &
   fglpointparameterfvext
#endif

! EXT_polygon_offset

#ifdef GL_EXT_polygon_offset
public :: &
   fglpolygonoffsetext
#endif

! EXT_subtexture

#ifdef GL_EXT_subtexture
public :: &
   f9y0gltexsubimage1dext, &
   f9y1gltexsubimage1dext, &
   f9y4gltexsubimage1dext, &
   f9y5gltexsubimage1dext, &
   f9y6gltexsubimage1dext, &
   f9y0gltexsubimage2dext, &
   f9y1gltexsubimage2dext, &
   f9y4gltexsubimage2dext, &
   f9y5gltexsubimage2dext, &
   f9y6gltexsubimage2dext
#endif

! EXT_texture3D

#ifdef GL_EXT_texture3D
public :: &
   f9y0glteximage3dext, &
   f9y1glteximage3dext, &
   f9y4glteximage3dext, &
   f9y5glteximage3dext, &
   f9y6glteximage3dext, &
   f9y0gltexsubimage3dext, &
   f9y1gltexsubimage3dext, &
   f9y4gltexsubimage3dext, &
   f9y5gltexsubimage3dext, &
   f9y6gltexsubimage3dext, &
   fglcopytexsubimage3dext
#endif

! EXT_texture_object

#ifdef GL_EXT_texture_object
public :: &
   fglgentexturesext, &
   fgldeletetexturesext, &
   fglbindtextureext, &
   fglprioritizetexturesext, &
   fglaretexturesresidentext, &
   f9eglaretexturesresidentext, &
   fglistextureext, &
   f9eglistextureext
#endif

! EXT_vertex_array

#ifdef GL_EXT_vertex_array
public :: &
   f9y0glvertexpointerext, &
   f9y1glvertexpointerext, &
   f9y4glvertexpointerext, &
   f9y6glvertexpointerext, &
   f9y7glvertexpointerext, &
   f9y0glnormalpointerext, &
   f9y1glnormalpointerext, &
   f9y4glnormalpointerext, &
   f9y5glnormalpointerext, &
   f9y6glnormalpointerext, &
   f9y7glnormalpointerext, &
   f9y0glcolorpointerext, &
   f9y1glcolorpointerext, &
   f9y4glcolorpointerext, &
   f9y5glcolorpointerext, &
   f9y6glcolorpointerext, &
   f9y7glcolorpointerext, &
   f9y0glindexpointerext, &
   f9y1glindexpointerext, &
   f9y4glindexpointerext, &
   f9y6glindexpointerext, &
   f9y7glindexpointerext, &
   f9y0gltexcoordpointerext, &
   f9y1gltexcoordpointerext, &
   f9y4gltexcoordpointerext, &
   f9y6gltexcoordpointerext, &
   f9y7gltexcoordpointerext, &
   fgledgeflagpointerext, &
   f9egledgeflagpointerext, &
   fglgetpointervext, &
   fglarrayelementext, &
   fgldrawarraysext
#endif

! SGIS_multisample

#ifdef GL_SGIS_multisample
public :: &
   fglsamplemasksgis, &
   f9eglsamplemasksgis, &
   fglsamplepatternsgis
#endif

! MESA_resize_buffers

#ifdef GL_MESA_resize_buffers
public :: &
   fglresizebuffersmesa
#endif

! MESA_window_pos

#ifdef GL_MESA_window_pos
public :: &
   fglwindowpos2imesa, &
   fglwindowpos2smesa, &
   f9aglwindowpos2smesa, &
   fglwindowpos2fmesa, &
   fglwindowpos2dmesa, &
   fglwindowpos2ivmesa, &
   fglwindowpos2svmesa, &
   f9aglwindowpos2svmesa, &
   fglwindowpos2fvmesa, &
   fglwindowpos2dvmesa, &
   fglwindowpos3imesa, &
   fglwindowpos3smesa, &
   f9aglwindowpos3smesa, &
   fglwindowpos3fmesa, &
   fglwindowpos3dmesa, &
   fglwindowpos3ivmesa, &
   fglwindowpos3svmesa, &
   f9aglwindowpos3svmesa, &
   fglwindowpos3fvmesa, &
   fglwindowpos3dvmesa, &
   fglwindowpos4imesa, &
   fglwindowpos4smesa, &
   f9aglwindowpos4smesa, &
   fglwindowpos4fmesa, &
   fglwindowpos4dmesa, &
   fglwindowpos4ivmesa, &
   fglwindowpos4svmesa, &
   f9aglwindowpos4svmesa, &
   fglwindowpos4fvmesa, &
   fglwindowpos4dvmesa
#endif

!-----------------
! Interface blocks
!-----------------

! OpenGL 1.0

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglaccum"
#endif
subroutine fglaccum(op,value)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: op,value
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglaccum
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglaccum@8' :: fglaccum
#endif
integer(kind=glenum) INTENT_IN op
real(kind=glfloat) INTENT_IN value
end subroutine fglaccum
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglalphafunc"
#endif
subroutine fglalphafunc(func,ref)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: func,ref
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglalphafunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglalphafunc@8' :: fglalphafunc
#endif
integer(kind=glenum) INTENT_IN func
real(kind=glclampf) INTENT_IN ref
end subroutine fglalphafunc
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglbegin"
#endif
subroutine fglbegin(mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglbegin
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglbegin@4' :: fglbegin
#endif
integer(kind=glenum) INTENT_IN mode
end subroutine fglbegin
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglbitmap"
#endif
subroutine fglbitmap(width,height,xorig,yorig,xmove,ymove, &
   bitmap)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: width,height,xorig,yorig,xmove,ymove
!MS$ATTRIBUTES REFERENCE :: bitmap
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglbitmap
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglbitmap@28' :: fglbitmap
#endif
integer(kind=glsizei) INTENT_IN width, height
real(kind=glfloat) INTENT_IN xorig, yorig, xmove, ymove
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN bitmap
end subroutine fglbitmap
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglbitmap"
#endif
subroutine f9bglbitmap(width,height,xorig,yorig,xmove,ymove, &
   bitmap,xsize)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: width,height,xorig,yorig,xmove,ymove
!MS$ATTRIBUTES REFERENCE :: bitmap,xsize
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglbitmap
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglbitmap@32' :: f9bglbitmap
#endif
integer(kind=glsizei) INTENT_IN width, height
real(kind=glfloat) INTENT_IN xorig, yorig, xmove, ymove
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN bitmap
integer(kind=glint) INTENT_IN xsize
end subroutine f9bglbitmap
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglbitmap"
#endif
subroutine f9cglbitmap(width,height,xorig,yorig,xmove,ymove, &
   bitmap,xsize)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: width,height,xorig,yorig,xmove,ymove
!MS$ATTRIBUTES REFERENCE :: bitmap,xsize
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglbitmap
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglbitmap@32' :: f9cglbitmap
#endif
integer(kind=glsizei) INTENT_IN width, height
real(kind=glfloat) INTENT_IN xorig, yorig, xmove, ymove
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN bitmap
integer(kind=glint) INTENT_IN xsize
end subroutine f9cglbitmap
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglblendfunc"
#endif
subroutine fglblendfunc(sfactor,dfactor)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: sfactor,dfactor
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglblendfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglblendfunc@8' :: fglblendfunc
#endif
integer(kind=glenum) INTENT_IN sfactor, dfactor
end subroutine fglblendfunc
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcalllist"
#endif
subroutine fglcalllist(list)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: list
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcalllist
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcalllist@4' :: fglcalllist
#endif
integer(kind=gluint) INTENT_IN list
end subroutine fglcalllist
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glcalllists"
#endif
subroutine f9y0glcalllists(n,itype,lists)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,itype,lists
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glcalllists
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glcalllists@12' :: f9y0glcalllists
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=glenum) INTENT_IN itype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN lists
end subroutine f9y0glcalllists
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glcalllists"
#endif
subroutine f9y1glcalllists(n,itype,lists)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,itype,lists
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glcalllists
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glcalllists@12' :: f9y1glcalllists
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=glenum) INTENT_IN itype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN lists
end subroutine f9y1glcalllists
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glcalllists"
#endif
subroutine f9y4glcalllists(n,itype,lists,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,itype,lists,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glcalllists
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glcalllists@16' :: f9y4glcalllists
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=glenum) INTENT_IN itype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN lists
integer(kind=glint) INTENT_IN length
end subroutine f9y4glcalllists
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glcalllists"
#endif
subroutine f9y5glcalllists(n,itype,lists)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,itype,lists
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glcalllists
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glcalllists@12' :: f9y5glcalllists
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=glenum) INTENT_IN itype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN lists
end subroutine f9y5glcalllists
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glcalllists"
#endif
subroutine f9y6glcalllists(n,itype,lists)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,itype,lists
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glcalllists
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glcalllists@12' :: f9y6glcalllists
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=glenum) INTENT_IN itype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN lists
end subroutine f9y6glcalllists
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglclear"
#endif
subroutine fglclear(mask)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mask
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglclear
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglclear@4' :: fglclear
#endif
integer(kind=glbitfield) INTENT_IN mask
end subroutine fglclear
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglclearaccum"
#endif
subroutine fglclearaccum(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglclearaccum
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglclearaccum@16' :: fglclearaccum
#endif
real(kind=glfloat) INTENT_IN red, green, blue, alpha
end subroutine fglclearaccum
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglclearcolor"
#endif
subroutine fglclearcolor(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglclearcolor
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglclearcolor@16' :: fglclearcolor
#endif
real(kind=glclampf) INTENT_IN red, green, blue, alpha
end subroutine fglclearcolor
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcleardepth"
#endif
subroutine fglcleardepth(depth)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: depth
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcleardepth
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcleardepth@4' :: fglcleardepth
#endif
real(kind=glclampd) INTENT_IN depth
end subroutine fglcleardepth
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglclearindex"
#endif
subroutine fglclearindex(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglclearindex
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglclearindex@4' :: fglclearindex
#endif
real(kind=glfloat) INTENT_IN c
end subroutine fglclearindex
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglclearstencil"
#endif
subroutine fglclearstencil(s)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglclearstencil
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglclearstencil@4' :: fglclearstencil
#endif
integer(kind=glint) INTENT_IN s
end subroutine fglclearstencil
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglclipplane"
#endif
subroutine fglclipplane(plane,equation)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: plane,equation
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglclipplane
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglclipplane@8' :: fglclipplane
#endif
integer(kind=glenum) INTENT_IN plane
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN equation
end subroutine fglclipplane
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3b"
#endif
subroutine fglcolor3b(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3b
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3b@12' :: fglcolor3b
#endif
integer(kind=glbyte) INTENT_IN red, green, blue
end subroutine fglcolor3b
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglcolor3b"
#endif
subroutine f9bglcolor3b(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglcolor3b
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglcolor3b@12' :: f9bglcolor3b
#endif
integer(kind=glbyte) INTENT_IN red, green, blue
end subroutine f9bglcolor3b
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglcolor3b"
#endif
subroutine f9cglcolor3b(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglcolor3b
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglcolor3b@12' :: f9cglcolor3b
#endif
integer(kind=glbyte) INTENT_IN red, green, blue
end subroutine f9cglcolor3b
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3bv"
#endif
subroutine fglcolor3bv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3bv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3bv@4' :: fglcolor3bv
#endif
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor3bv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglcolor3bv"
#endif
subroutine f9bglcolor3bv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglcolor3bv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglcolor3bv@4' :: f9bglcolor3bv
#endif
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9bglcolor3bv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglcolor3bv"
#endif
subroutine f9cglcolor3bv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglcolor3bv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglcolor3bv@4' :: f9cglcolor3bv
#endif
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9cglcolor3bv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3d"
#endif
subroutine fglcolor3d(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3d@12' :: fglcolor3d
#endif
real(kind=gldouble) INTENT_IN red, green, blue
end subroutine fglcolor3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3dv"
#endif
subroutine fglcolor3dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3dv@4' :: fglcolor3dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor3dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3f"
#endif
subroutine fglcolor3f(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3f@12' :: fglcolor3f
#endif
real(kind=glfloat) INTENT_IN red, green, blue
end subroutine fglcolor3f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3fv"
#endif
subroutine fglcolor3fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3fv@4' :: fglcolor3fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor3fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3i"
#endif
subroutine fglcolor3i(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3i@12' :: fglcolor3i
#endif
integer(kind=glint) INTENT_IN red, green, blue
end subroutine fglcolor3i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3iv"
#endif
subroutine fglcolor3iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3iv@4' :: fglcolor3iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor3iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3s"
#endif
subroutine fglcolor3s(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3s@12' :: fglcolor3s
#endif
integer(kind=glshort) INTENT_IN red, green, blue
end subroutine fglcolor3s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglcolor3s"
#endif
subroutine f9aglcolor3s(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglcolor3s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglcolor3s@12' :: f9aglcolor3s
#endif
integer(kind=glshort) INTENT_IN red, green, blue
end subroutine f9aglcolor3s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3sv"
#endif
subroutine fglcolor3sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3sv@4' :: fglcolor3sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor3sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglcolor3sv"
#endif
subroutine f9aglcolor3sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglcolor3sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglcolor3sv@4' :: f9aglcolor3sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9aglcolor3sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3ub"
#endif
subroutine fglcolor3ub(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3ub
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3ub@12' :: fglcolor3ub
#endif
integer(kind=glubyte) INTENT_IN red, green, blue
end subroutine fglcolor3ub
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglcolor3ub"
#endif
subroutine f9bglcolor3ub(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglcolor3ub
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglcolor3ub@12' :: f9bglcolor3ub
#endif
integer(kind=glubyte) INTENT_IN red, green, blue
end subroutine f9bglcolor3ub
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglcolor3ub"
#endif
subroutine f9cglcolor3ub(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglcolor3ub
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglcolor3ub@12' :: f9cglcolor3ub
#endif
integer(kind=glubyte) INTENT_IN red, green, blue
end subroutine f9cglcolor3ub
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3ubv"
#endif
subroutine fglcolor3ubv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3ubv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3ubv@4' :: fglcolor3ubv
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor3ubv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglcolor3ubv"
#endif
subroutine f9bglcolor3ubv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglcolor3ubv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglcolor3ubv@4' :: f9bglcolor3ubv
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9bglcolor3ubv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglcolor3ubv"
#endif
subroutine f9cglcolor3ubv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglcolor3ubv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglcolor3ubv@4' :: f9cglcolor3ubv
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9cglcolor3ubv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3ui"
#endif
subroutine fglcolor3ui(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3ui
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3ui@12' :: fglcolor3ui
#endif
integer(kind=gluint) INTENT_IN red, green, blue
end subroutine fglcolor3ui
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3uiv"
#endif
subroutine fglcolor3uiv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3uiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3uiv@4' :: fglcolor3uiv
#endif
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor3uiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3us"
#endif
subroutine fglcolor3us(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3us
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3us@12' :: fglcolor3us
#endif
integer(kind=glushort) INTENT_IN red, green, blue
end subroutine fglcolor3us
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglcolor3us"
#endif
subroutine f9aglcolor3us(red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglcolor3us
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglcolor3us@12' :: f9aglcolor3us
#endif
integer(kind=glushort) INTENT_IN red, green, blue
end subroutine f9aglcolor3us
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor3usv"
#endif
subroutine fglcolor3usv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor3usv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor3usv@4' :: fglcolor3usv
#endif
integer(kind=glushort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor3usv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglcolor3usv"
#endif
subroutine f9aglcolor3usv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglcolor3usv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglcolor3usv@4' :: f9aglcolor3usv
#endif
integer(kind=glushort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9aglcolor3usv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4b"
#endif
subroutine fglcolor4b(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4b
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4b@16' :: fglcolor4b
#endif
integer(kind=glbyte) INTENT_IN red, green, blue, alpha
end subroutine fglcolor4b
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglcolor4b"
#endif
subroutine f9bglcolor4b(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglcolor4b
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglcolor4b@16' :: f9bglcolor4b
#endif
integer(kind=glbyte) INTENT_IN red, green, blue, alpha
end subroutine f9bglcolor4b
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglcolor4b"
#endif
subroutine f9cglcolor4b(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglcolor4b
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglcolor4b@16' :: f9cglcolor4b
#endif
integer(kind=glbyte) INTENT_IN red, green, blue, alpha
end subroutine f9cglcolor4b
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4bv"
#endif
subroutine fglcolor4bv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4bv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4bv@4' :: fglcolor4bv
#endif
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor4bv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglcolor4bv"
#endif
subroutine f9bglcolor4bv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglcolor4bv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglcolor4bv@4' :: f9bglcolor4bv
#endif
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9bglcolor4bv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglcolor4bv"
#endif
subroutine f9cglcolor4bv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglcolor4bv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglcolor4bv@4' :: f9cglcolor4bv
#endif
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9cglcolor4bv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4d"
#endif
subroutine fglcolor4d(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4d@16' :: fglcolor4d
#endif
real(kind=gldouble) INTENT_IN red, green, blue, alpha
end subroutine fglcolor4d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4dv"
#endif
subroutine fglcolor4dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4dv@4' :: fglcolor4dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor4dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4f"
#endif
subroutine fglcolor4f(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4f@16' :: fglcolor4f
#endif
real(kind=glfloat) INTENT_IN red, green, blue, alpha
end subroutine fglcolor4f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4fv"
#endif
subroutine fglcolor4fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4fv@4' :: fglcolor4fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor4fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4i"
#endif
subroutine fglcolor4i(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4i@16' :: fglcolor4i
#endif
integer(kind=glint) INTENT_IN red, green, blue, alpha
end subroutine fglcolor4i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4iv"
#endif
subroutine fglcolor4iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4iv@4' :: fglcolor4iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor4iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4s"
#endif
subroutine fglcolor4s(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4s@16' :: fglcolor4s
#endif
integer(kind=glshort) INTENT_IN red, green, blue, alpha
end subroutine fglcolor4s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglcolor4s"
#endif
subroutine f9aglcolor4s(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglcolor4s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglcolor4s@16' :: f9aglcolor4s
#endif
integer(kind=glshort) INTENT_IN red, green, blue, alpha
end subroutine f9aglcolor4s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4sv"
#endif
subroutine fglcolor4sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4sv@4' :: fglcolor4sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor4sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglcolor4sv"
#endif
subroutine f9aglcolor4sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglcolor4sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglcolor4sv@4' :: f9aglcolor4sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9aglcolor4sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4ub"
#endif
subroutine fglcolor4ub(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4ub
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4ub@16' :: fglcolor4ub
#endif
integer(kind=glubyte) INTENT_IN red, green, blue, alpha
end subroutine fglcolor4ub
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglcolor4ub"
#endif
subroutine f9bglcolor4ub(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglcolor4ub
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglcolor4ub@16' :: f9bglcolor4ub
#endif
integer(kind=glubyte) INTENT_IN red, green, blue, alpha
end subroutine f9bglcolor4ub
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglcolor4ub"
#endif
subroutine f9cglcolor4ub(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglcolor4ub
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglcolor4ub@16' :: f9cglcolor4ub
#endif
integer(kind=glubyte) INTENT_IN red, green, blue, alpha
end subroutine f9cglcolor4ub
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4ubv"
#endif
subroutine fglcolor4ubv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4ubv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4ubv@4' :: fglcolor4ubv
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor4ubv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglcolor4ubv"
#endif
subroutine f9bglcolor4ubv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglcolor4ubv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglcolor4ubv@4' :: f9bglcolor4ubv
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9bglcolor4ubv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglcolor4ubv"
#endif
subroutine f9cglcolor4ubv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglcolor4ubv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglcolor4ubv@4' :: f9cglcolor4ubv
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9cglcolor4ubv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4ui"
#endif
subroutine fglcolor4ui(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4ui
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4ui@16' :: fglcolor4ui
#endif
integer(kind=gluint) INTENT_IN red, green, blue, alpha
end subroutine fglcolor4ui
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4uiv"
#endif
subroutine fglcolor4uiv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4uiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4uiv@4' :: fglcolor4uiv
#endif
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor4uiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4us"
#endif
subroutine fglcolor4us(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4us
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4us@16' :: fglcolor4us
#endif
integer(kind=glushort) INTENT_IN red, green, blue, alpha
end subroutine fglcolor4us
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglcolor4us"
#endif
subroutine f9aglcolor4us(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglcolor4us
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglcolor4us@16' :: f9aglcolor4us
#endif
integer(kind=glushort) INTENT_IN red, green, blue, alpha
end subroutine f9aglcolor4us
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolor4usv"
#endif
subroutine fglcolor4usv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolor4usv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolor4usv@4' :: fglcolor4usv
#endif
integer(kind=glushort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglcolor4usv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglcolor4usv"
#endif
subroutine f9aglcolor4usv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglcolor4usv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglcolor4usv@4' :: f9aglcolor4usv
#endif
integer(kind=glushort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9aglcolor4usv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolormask"
#endif
subroutine fglcolormask(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolormask
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolormask@16' :: fglcolormask
#endif
logical(kind=glboolean) INTENT_IN red, green, blue, alpha
end subroutine fglcolormask
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9eglcolormask"
#endif
subroutine f9eglcolormask(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9eglcolormask
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9eglcolormask@16' :: f9eglcolormask
#endif
integer(kind=glint) INTENT_IN red, green, blue, alpha
end subroutine f9eglcolormask
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcolormaterial"
#endif
subroutine fglcolormaterial(face,mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: face,mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcolormaterial
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcolormaterial@8' :: fglcolormaterial
#endif
integer(kind=glenum) INTENT_IN face, mode
end subroutine fglcolormaterial
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcopypixels"
#endif
subroutine fglcopypixels(x,y,width,height,xtype)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,width,height,xtype
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcopypixels
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcopypixels@20' :: fglcopypixels
#endif
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN xtype
end subroutine fglcopypixels
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcullface"
#endif
subroutine fglcullface(mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcullface
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcullface@4' :: fglcullface
#endif
integer(kind=glenum) INTENT_IN mode
end subroutine fglcullface
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgldeletelists"
#endif
subroutine fgldeletelists(list,xrange)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: list,xrange
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgldeletelists
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgldeletelists@8' :: fgldeletelists
#endif
integer(kind=gluint) INTENT_IN list
integer(kind=glsizei) INTENT_IN xrange
end subroutine fgldeletelists
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgldepthfunc"
#endif
subroutine fgldepthfunc(func)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: func
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgldepthfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgldepthfunc@4' :: fgldepthfunc
#endif
integer(kind=glenum) INTENT_IN func
end subroutine fgldepthfunc
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgldepthmask"
#endif
subroutine fgldepthmask(flag)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: flag
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgldepthmask
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgldepthmask@4' :: fgldepthmask
#endif
logical(kind=glboolean) INTENT_IN flag
end subroutine fgldepthmask
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9egldepthmask"
#endif
subroutine f9egldepthmask(flag)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: flag
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9egldepthmask
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9egldepthmask@4' :: f9egldepthmask
#endif
integer(kind=glint) INTENT_IN flag
end subroutine f9egldepthmask
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgldepthrange"
#endif
subroutine fgldepthrange(near,far)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: near,far
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgldepthrange
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgldepthrange@8' :: fgldepthrange
#endif
real(kind=glclampd) INTENT_IN near, far
end subroutine fgldepthrange
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgldisable"
#endif
subroutine fgldisable(cap)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: cap
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgldisable
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgldisable@4' :: fgldisable
#endif
integer(kind=glenum) INTENT_IN cap
end subroutine fgldisable
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgldrawbuffer"
#endif
subroutine fgldrawbuffer(mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgldrawbuffer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgldrawbuffer@4' :: fgldrawbuffer
#endif
integer(kind=glenum) INTENT_IN mode
end subroutine fgldrawbuffer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0gldrawpixels"
#endif
subroutine f9y0gldrawpixels(width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gldrawpixels
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gldrawpixels@24' :: f9y0gldrawpixels
#endif
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0gldrawpixels
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1gldrawpixels"
#endif
subroutine f9y1gldrawpixels(width,height,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gldrawpixels
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gldrawpixels@20' :: f9y1gldrawpixels
#endif
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y1gldrawpixels
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4gldrawpixels"
#endif
subroutine f9y4gldrawpixels(width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gldrawpixels
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gldrawpixels@24' :: f9y4gldrawpixels
#endif
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4gldrawpixels
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5gldrawpixels"
#endif
subroutine f9y5gldrawpixels(width,height,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5gldrawpixels
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5gldrawpixels@20' :: f9y5gldrawpixels
#endif
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y5gldrawpixels
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6gldrawpixels"
#endif
subroutine f9y6gldrawpixels(width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gldrawpixels
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gldrawpixels@24' :: f9y6gldrawpixels
#endif
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6gldrawpixels
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgledgeflag"
#endif
subroutine fgledgeflag(flag)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: flag
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgledgeflag
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgledgeflag@4' :: fgledgeflag
#endif
logical(kind=glboolean) INTENT_IN flag
end subroutine fgledgeflag
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9egledgeflag"
#endif
subroutine f9egledgeflag(flag)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: flag
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9egledgeflag
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9egledgeflag@4' :: f9egledgeflag
#endif
integer(kind=glint) INTENT_IN flag
end subroutine f9egledgeflag
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgledgeflagv"
#endif
subroutine fgledgeflagv(flag)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: flag
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgledgeflagv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgledgeflagv@4' :: fgledgeflagv
#endif
logical(kind=glboolean), dimension(DEFER_DIM) INTENT_IN flag
end subroutine fgledgeflagv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9egledgeflagv"
#endif
subroutine f9egledgeflagv(flag,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: flag,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9egledgeflagv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9egledgeflagv@8' :: f9egledgeflagv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN flag
integer(kind=glint) INTENT_IN length
end subroutine f9egledgeflagv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglenable"
#endif
subroutine fglenable(cap)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: cap
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglenable
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglenable@4' :: fglenable
#endif
integer(kind=glenum) INTENT_IN cap
end subroutine fglenable
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglend"
#endif
subroutine fglend()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglend
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglend@0' :: fglend
#endif
end subroutine fglend
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglendlist"
#endif
subroutine fglendlist()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglendlist
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglendlist@0' :: fglendlist
#endif
end subroutine fglendlist
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalcoord1d"
#endif
subroutine fglevalcoord1d(u)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: u
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalcoord1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalcoord1d@4' :: fglevalcoord1d
#endif
real(kind=gldouble) INTENT_IN u
end subroutine fglevalcoord1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalcoord1dv"
#endif
subroutine fglevalcoord1dv(u)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: u
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalcoord1dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalcoord1dv@4' :: fglevalcoord1dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN u
end subroutine fglevalcoord1dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalcoord1f"
#endif
subroutine fglevalcoord1f(u)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: u
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalcoord1f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalcoord1f@4' :: fglevalcoord1f
#endif
real(kind=glfloat) INTENT_IN u
end subroutine fglevalcoord1f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalcoord1fv"
#endif
subroutine fglevalcoord1fv(u)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: u
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalcoord1fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalcoord1fv@4' :: fglevalcoord1fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN u
end subroutine fglevalcoord1fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalcoord2d"
#endif
subroutine fglevalcoord2d(u,v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: u,v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalcoord2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalcoord2d@8' :: fglevalcoord2d
#endif
real(kind=gldouble) INTENT_IN u, v
end subroutine fglevalcoord2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalcoord2dv"
#endif
subroutine fglevalcoord2dv(u)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: u
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalcoord2dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalcoord2dv@4' :: fglevalcoord2dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN u
end subroutine fglevalcoord2dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalcoord2f"
#endif
subroutine fglevalcoord2f(u,v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: u,v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalcoord2f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalcoord2f@8' :: fglevalcoord2f
#endif
real(kind=glfloat) INTENT_IN u, v
end subroutine fglevalcoord2f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalcoord2fv"
#endif
subroutine fglevalcoord2fv(u)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: u
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalcoord2fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalcoord2fv@4' :: fglevalcoord2fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN u
end subroutine fglevalcoord2fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalmesh1"
#endif
subroutine fglevalmesh1(mode,i1,i2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,i1,i2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalmesh1
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalmesh1@12' :: fglevalmesh1
#endif
integer(kind=glenum) INTENT_IN mode
integer(kind=glint) INTENT_IN i1, i2
end subroutine fglevalmesh1
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalmesh2"
#endif
subroutine fglevalmesh2(mode,i1,i2,j1,j2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,i1,i2,j1,j2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalmesh2
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalmesh2@20' :: fglevalmesh2
#endif
integer(kind=glenum) INTENT_IN mode
integer(kind=glint) INTENT_IN i1, i2, j1, j2
end subroutine fglevalmesh2
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalpoint1"
#endif
subroutine fglevalpoint1(i)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: i
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalpoint1
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalpoint1@4' :: fglevalpoint1
#endif
integer(kind=glint) INTENT_IN i
end subroutine fglevalpoint1
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglevalpoint2"
#endif
subroutine fglevalpoint2(i,j)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: i,j
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglevalpoint2
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglevalpoint2@8' :: fglevalpoint2
#endif
integer(kind=glint) INTENT_IN i, j
end subroutine fglevalpoint2
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglfeedbackbuffer"
#endif
subroutine fglfeedbackbuffer(xsize,xtype,buffer)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,buffer
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglfeedbackbuffer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglfeedbackbuffer@12' :: fglfeedbackbuffer
#endif
integer(kind=glsizei) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: buffer
end subroutine fglfeedbackbuffer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglfinish"
#endif
subroutine fglfinish()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglfinish
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglfinish@0' :: fglfinish
#endif
end subroutine fglfinish
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglflush"
#endif
subroutine fglflush()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglflush
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglflush@0' :: fglflush
#endif
end subroutine fglflush
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglfogf"
#endif
subroutine fglfogf(pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglfogf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglfogf@8' :: fglfogf
#endif
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat) INTENT_IN param
end subroutine fglfogf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglfogfv"
#endif
subroutine fglfogfv(pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglfogfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglfogfv@8' :: fglfogfv
#endif
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
end subroutine fglfogfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglfogi"
#endif
subroutine fglfogi(pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglfogi
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglfogi@8' :: fglfogi
#endif
integer(kind=glenum) INTENT_IN pname
integer(kind=glint) INTENT_IN param
end subroutine fglfogi
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglfogiv"
#endif
subroutine fglfogiv(pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglfogiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglfogiv@8' :: fglfogiv
#endif
integer(kind=glenum) INTENT_IN pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
end subroutine fglfogiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglfrontface"
#endif
subroutine fglfrontface(mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglfrontface
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglfrontface@4' :: fglfrontface
#endif
integer(kind=glenum) INTENT_IN mode
end subroutine fglfrontface
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglfrustum"
#endif
subroutine fglfrustum(left,right,bottom,top,near,far)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: left,right,bottom,top,near,far
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglfrustum
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglfrustum@24' :: fglfrustum
#endif
real(kind=gldouble) INTENT_IN left, right, bottom, top, near, far
end subroutine fglfrustum
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgenlists"
#endif
function fglgenlists(xrange) result(resfglgenlists)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xrange
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgenlists
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgenlists@4' :: fglgenlists
#endif
integer(kind=glsizei) INTENT_IN xrange
integer(kind=gluint) :: resfglgenlists
end function fglgenlists
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetbooleanv"
#endif
subroutine fglgetbooleanv(pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetbooleanv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetbooleanv@8' :: fglgetbooleanv
#endif
integer(kind=glenum) INTENT_IN pname
logical(kind=glboolean), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgetbooleanv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9eglgetbooleanv"
#endif
subroutine f9eglgetbooleanv(pname,params,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9eglgetbooleanv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9eglgetbooleanv@12' :: f9eglgetbooleanv
#endif
integer(kind=glenum) INTENT_IN pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
integer(kind=glint) INTENT_IN length
end subroutine f9eglgetbooleanv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetclipplane"
#endif
subroutine fglgetclipplane(plane,equation)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: plane,equation
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetclipplane
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetclipplane@8' :: fglgetclipplane
#endif
integer(kind=glenum) INTENT_IN plane
real(kind=gldouble), dimension(DEFER_DIM), intent(out) :: equation
end subroutine fglgetclipplane
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetdoublev"
#endif
subroutine fglgetdoublev(pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetdoublev
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetdoublev@8' :: fglgetdoublev
#endif
integer(kind=glenum) INTENT_IN pname
real(kind=gldouble), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgetdoublev
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgeterror"
#endif
function fglgeterror() result(resfglgeterror)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgeterror
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgeterror@0' :: fglgeterror
#endif
integer(kind=glenum) :: resfglgeterror
end function fglgeterror
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetfloatv"
#endif
subroutine fglgetfloatv(pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetfloatv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetfloatv@8' :: fglgetfloatv
#endif
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgetfloatv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetintegerv"
#endif
subroutine fglgetintegerv(pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetintegerv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetintegerv@8' :: fglgetintegerv
#endif
integer(kind=glenum) INTENT_IN pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgetintegerv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetlightfv"
#endif
subroutine fglgetlightfv(light,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: light,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetlightfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetlightfv@12' :: fglgetlightfv
#endif
integer(kind=glenum) INTENT_IN light, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgetlightfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetlightiv"
#endif
subroutine fglgetlightiv(light,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: light,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetlightiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetlightiv@12' :: fglgetlightiv
#endif
integer(kind=glenum) INTENT_IN light, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgetlightiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetmapdv"
#endif
subroutine fglgetmapdv(xtarget,query,v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,query,v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetmapdv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetmapdv@12' :: fglgetmapdv
#endif
integer(kind=glenum) INTENT_IN xtarget, query
real(kind=gldouble), dimension(DEFER_DIM), intent(out) :: v
end subroutine fglgetmapdv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetmapfv"
#endif
subroutine fglgetmapfv(xtarget,query,v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,query,v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetmapfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetmapfv@12' :: fglgetmapfv
#endif
integer(kind=glenum) INTENT_IN xtarget, query
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: v
end subroutine fglgetmapfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetmapiv"
#endif
subroutine fglgetmapiv(xtarget,query,v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,query,v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetmapiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetmapiv@12' :: fglgetmapiv
#endif
integer(kind=glenum) INTENT_IN xtarget, query
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: v
end subroutine fglgetmapiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetmaterialfv"
#endif
subroutine fglgetmaterialfv(face,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: face,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetmaterialfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetmaterialfv@12' :: fglgetmaterialfv
#endif
integer(kind=glenum) INTENT_IN face, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgetmaterialfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetmaterialiv"
#endif
subroutine fglgetmaterialiv(face,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: face,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetmaterialiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetmaterialiv@12' :: fglgetmaterialiv
#endif
integer(kind=glenum) INTENT_IN face, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgetmaterialiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetpixelmapfv"
#endif
subroutine fglgetpixelmapfv(map,values)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: map,values
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetpixelmapfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetpixelmapfv@8' :: fglgetpixelmapfv
#endif
integer(kind=glenum) INTENT_IN map
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: values
end subroutine fglgetpixelmapfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetpixelmapuiv"
#endif
subroutine fglgetpixelmapuiv(map,values)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: map,values
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetpixelmapuiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetpixelmapuiv@8' :: fglgetpixelmapuiv
#endif
integer(kind=glenum) INTENT_IN map
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: values
end subroutine fglgetpixelmapuiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetpixelmapusv"
#endif
subroutine fglgetpixelmapusv(map,values)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: map,values
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetpixelmapusv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetpixelmapusv@8' :: fglgetpixelmapusv
#endif
integer(kind=glenum) INTENT_IN map
integer(kind=glushort), dimension(DEFER_DIM), intent(out) :: values
end subroutine fglgetpixelmapusv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglgetpixelmapusv"
#endif
subroutine f9aglgetpixelmapusv(map,values)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: map,values
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglgetpixelmapusv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglgetpixelmapusv@8' :: f9aglgetpixelmapusv
#endif
integer(kind=glenum) INTENT_IN map
integer(kind=glushort), dimension(DEFER_DIM), intent(out) :: values
end subroutine f9aglgetpixelmapusv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetpolygonstipple"
#endif
subroutine fglgetpolygonstipple(mask,xsize)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mask,xsize
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetpolygonstipple
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetpolygonstipple@8' :: fglgetpolygonstipple
#endif
integer(kind=glubyte), dimension(DEFER_DIM), intent(out) :: mask
integer(kind=glint) INTENT_IN xsize
end subroutine fglgetpolygonstipple
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglgetpolygonstipple"
#endif
subroutine f9bglgetpolygonstipple(mask,xsize)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mask,xsize
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglgetpolygonstipple
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglgetpolygonstipple@8' :: f9bglgetpolygonstipple
#endif
integer(kind=glubyte), dimension(DEFER_DIM), intent(out) :: mask
integer(kind=glint) INTENT_IN xsize
end subroutine f9bglgetpolygonstipple
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglgetpolygonstipple"
#endif
subroutine f9cglgetpolygonstipple(mask,xsize)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mask,xsize
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglgetpolygonstipple
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglgetpolygonstipple@8' :: f9cglgetpolygonstipple
#endif
integer(kind=glubyte), dimension(DEFER_DIM), intent(out) :: mask
integer(kind=glint) INTENT_IN xsize
end subroutine f9cglgetpolygonstipple
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9yglgetstringlen"
#endif
function f9yglgetstringlen(name) result(resf9yglgetstringlen)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: name
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglgetstringlen
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglgetstringlen@4' :: f9yglgetstringlen
#endif
integer(kind=glenum) INTENT_IN name
integer(kind=glint) :: resf9yglgetstringlen
end function f9yglgetstringlen
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9yglgetstring"
#endif
subroutine f9yglgetstring(resf9yglgetstring)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: resf9yglgetstring
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglgetstring
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglgetstring@4' :: f9yglgetstring
#endif
integer(kind=gluint), intent(out), dimension(DEFER_DIM) :: resf9yglgetstring
end subroutine f9yglgetstring
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9yglgetstring1"
#endif
function f9yglgetstring1(pos) result(resf9yglgetstring1)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pos
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglgetstring1
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglgetstring1@4' :: f9yglgetstring1
#endif
integer(kind=glint) INTENT_IN pos
integer(kind=gluint) :: resf9yglgetstring1
end function f9yglgetstring1
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9yglgetstringfree"
#endif
#ifdef PURE_FUNCTIONS
function f9yglgetstringfree() result(resf9yglgetstringfree)
#else
subroutine f9yglgetstringfree()
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglgetstringfree
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglgetstringfree@0' :: f9yglgetstringfree
#endif
#ifdef PURE_FUNCTIONS
integer(kind=glint) :: resf9yglgetstringfree
#endif
#ifdef PURE_FUNCTIONS
end function f9yglgetstringfree
#else
end subroutine f9yglgetstringfree
#endif
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgettexenvfv"
#endif
subroutine fglgettexenvfv(xtarget,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgettexenvfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgettexenvfv@12' :: fglgettexenvfv
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgettexenvfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgettexenviv"
#endif
subroutine fglgettexenviv(xtarget,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgettexenviv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgettexenviv@12' :: fglgettexenviv
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgettexenviv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgettexgendv"
#endif
subroutine fglgettexgendv(coord,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: coord,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgettexgendv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgettexgendv@12' :: fglgettexgendv
#endif
integer(kind=glenum) INTENT_IN coord, pname
real(kind=gldouble), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgettexgendv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgettexgenfv"
#endif
subroutine fglgettexgenfv(coord,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: coord,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgettexgenfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgettexgenfv@12' :: fglgettexgenfv
#endif
integer(kind=glenum) INTENT_IN coord, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgettexgenfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgettexgeniv"
#endif
subroutine fglgettexgeniv(coord,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: coord,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgettexgeniv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgettexgeniv@12' :: fglgettexgeniv
#endif
integer(kind=glenum) INTENT_IN coord, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgettexgeniv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glgetteximage"
#endif
subroutine f9y0glgetteximage(xtarget,level,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glgetteximage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glgetteximage@24' :: f9y0glgetteximage
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0glgetteximage
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glgetteximage"
#endif
subroutine f9y1glgetteximage(xtarget,level,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glgetteximage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glgetteximage@20' :: f9y1glgetteximage
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: pixels
end subroutine f9y1glgetteximage
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glgetteximage"
#endif
subroutine f9y4glgetteximage(xtarget,level,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glgetteximage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glgetteximage@24' :: f9y4glgetteximage
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4glgetteximage
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glgetteximage"
#endif
subroutine f9y5glgetteximage(xtarget,level,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glgetteximage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glgetteximage@20' :: f9y5glgetteximage
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: pixels
end subroutine f9y5glgetteximage
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glgetteximage"
#endif
subroutine f9y6glgetteximage(xtarget,level,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glgetteximage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glgetteximage@24' :: f9y6glgetteximage
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM), intent(out) :: pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6glgetteximage
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgettexlevelparameterfv"
#endif
subroutine fglgettexlevelparameterfv(xtarget,level,pname, &
   params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,pname
!MS$ATTRIBUTES REFERENCE :: params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgettexlevelparameterfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgettexlevelparameterfv@16' :: fglgettexlevelparameterfv
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgettexlevelparameterfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgettexlevelparameteriv"
#endif
subroutine fglgettexlevelparameteriv(xtarget,level,pname, &
   params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,pname
!MS$ATTRIBUTES REFERENCE :: params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgettexlevelparameteriv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgettexlevelparameteriv@16' :: fglgettexlevelparameteriv
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgettexlevelparameteriv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgettexparameterfv"
#endif
subroutine fglgettexparameterfv(xtarget,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgettexparameterfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgettexparameterfv@12' :: fglgettexparameterfv
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgettexparameterfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgettexparameteriv"
#endif
subroutine fglgettexparameteriv(xtarget,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgettexparameteriv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgettexparameteriv@12' :: fglgettexparameteriv
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgettexparameteriv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglhint"
#endif
subroutine fglhint(xtarget,mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglhint
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglhint@8' :: fglhint
#endif
integer(kind=glenum) INTENT_IN xtarget, mode
end subroutine fglhint
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglindexmask"
#endif
subroutine fglindexmask(mask)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mask
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglindexmask
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglindexmask@4' :: fglindexmask
#endif
integer(kind=gluint) INTENT_IN mask
end subroutine fglindexmask
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglindexd"
#endif
subroutine fglindexd(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglindexd
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglindexd@4' :: fglindexd
#endif
real(kind=gldouble) INTENT_IN c
end subroutine fglindexd
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglindexdv"
#endif
subroutine fglindexdv(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglindexdv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglindexdv@4' :: fglindexdv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN c
end subroutine fglindexdv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglindexf"
#endif
subroutine fglindexf(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglindexf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglindexf@4' :: fglindexf
#endif
real(kind=glfloat) INTENT_IN c
end subroutine fglindexf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglindexfv"
#endif
subroutine fglindexfv(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglindexfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglindexfv@4' :: fglindexfv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN c
end subroutine fglindexfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglindexi"
#endif
subroutine fglindexi(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglindexi
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglindexi@4' :: fglindexi
#endif
integer(kind=glint) INTENT_IN c
end subroutine fglindexi
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglindexiv"
#endif
subroutine fglindexiv(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglindexiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglindexiv@4' :: fglindexiv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN c
end subroutine fglindexiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglindexs"
#endif
subroutine fglindexs(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglindexs
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglindexs@4' :: fglindexs
#endif
integer(kind=glshort) INTENT_IN c
end subroutine fglindexs
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglindexs"
#endif
subroutine f9aglindexs(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglindexs
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglindexs@4' :: f9aglindexs
#endif
integer(kind=glshort) INTENT_IN c
end subroutine f9aglindexs
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglindexsv"
#endif
subroutine fglindexsv(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglindexsv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglindexsv@4' :: fglindexsv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN c
end subroutine fglindexsv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglindexsv"
#endif
subroutine f9aglindexsv(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglindexsv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglindexsv@4' :: f9aglindexsv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN c
end subroutine f9aglindexsv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglinitnames"
#endif
subroutine fglinitnames()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglinitnames
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglinitnames@0' :: fglinitnames
#endif
end subroutine fglinitnames
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglisenabled"
#endif
function fglisenabled(cap) result(resfglisenabled)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: cap
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglisenabled
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglisenabled@4' :: fglisenabled
#endif
integer(kind=gluint) INTENT_IN cap
logical(kind=glboolean) :: resfglisenabled
end function fglisenabled
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9eglisenabled"
#endif
function f9eglisenabled(cap) result(resf9eglisenabled)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: cap
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9eglisenabled
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9eglisenabled@4' :: f9eglisenabled
#endif
integer(kind=gluint) INTENT_IN cap
integer(kind=glint) :: resf9eglisenabled
end function f9eglisenabled
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglislist"
#endif
function fglislist(list) result(resfglislist)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: list
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglislist
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglislist@4' :: fglislist
#endif
integer(kind=gluint) INTENT_IN list
logical(kind=glboolean) :: resfglislist
end function fglislist
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9eglislist"
#endif
function f9eglislist(list) result(resf9eglislist)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: list
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9eglislist
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9eglislist@4' :: f9eglislist
#endif
integer(kind=gluint) INTENT_IN list
integer(kind=glint) :: resf9eglislist
end function f9eglislist
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllightmodelf"
#endif
subroutine fgllightmodelf(pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllightmodelf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllightmodelf@8' :: fgllightmodelf
#endif
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat) INTENT_IN param
end subroutine fgllightmodelf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllightmodelfv"
#endif
subroutine fgllightmodelfv(pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllightmodelfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllightmodelfv@8' :: fgllightmodelfv
#endif
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
end subroutine fgllightmodelfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllightmodeli"
#endif
subroutine fgllightmodeli(pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllightmodeli
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllightmodeli@8' :: fgllightmodeli
#endif
integer(kind=glenum) INTENT_IN pname
integer(kind=glint) INTENT_IN param
end subroutine fgllightmodeli
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllightmodeliv"
#endif
subroutine fgllightmodeliv(pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllightmodeliv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllightmodeliv@8' :: fgllightmodeliv
#endif
integer(kind=glenum) INTENT_IN pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
end subroutine fgllightmodeliv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllightf"
#endif
subroutine fgllightf(light,pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: light,pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllightf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllightf@12' :: fgllightf
#endif
integer(kind=glenum) INTENT_IN light, pname
real(kind=glfloat) INTENT_IN param
end subroutine fgllightf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllightfv"
#endif
subroutine fgllightfv(light,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: light,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllightfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllightfv@12' :: fgllightfv
#endif
integer(kind=glenum) INTENT_IN light, pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
end subroutine fgllightfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllighti"
#endif
subroutine fgllighti(light,pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: light,pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllighti
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllighti@12' :: fgllighti
#endif
integer(kind=glenum) INTENT_IN light, pname
integer(kind=glint) INTENT_IN param
end subroutine fgllighti
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllightiv"
#endif
subroutine fgllightiv(light,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: light,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllightiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllightiv@12' :: fgllightiv
#endif
integer(kind=glenum) INTENT_IN light, pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
end subroutine fgllightiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllinestipple"
#endif
subroutine fgllinestipple(factor,pattern)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: factor,pattern
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllinestipple
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllinestipple@8' :: fgllinestipple
#endif
integer(kind=glint) INTENT_IN factor
integer(kind=glushort) INTENT_IN pattern
end subroutine fgllinestipple
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9agllinestipple"
#endif
subroutine f9agllinestipple(factor,pattern)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: factor,pattern
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9agllinestipple
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9agllinestipple@8' :: f9agllinestipple
#endif
integer(kind=glint) INTENT_IN factor
integer(kind=glushort) INTENT_IN pattern
end subroutine f9agllinestipple
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllinewidth"
#endif
subroutine fgllinewidth(width)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: width
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllinewidth
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllinewidth@4' :: fgllinewidth
#endif
real(kind=glfloat) INTENT_IN width
end subroutine fgllinewidth
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllistbase"
#endif
subroutine fgllistbase(base)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: base
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllistbase
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllistbase@4' :: fgllistbase
#endif
integer(kind=gluint) INTENT_IN base
end subroutine fgllistbase
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglloadidentity"
#endif
subroutine fglloadidentity()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglloadidentity
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglloadidentity@0' :: fglloadidentity
#endif
end subroutine fglloadidentity
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglloadmatrixd"
#endif
subroutine fglloadmatrixd(m)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: m
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglloadmatrixd
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglloadmatrixd@4' :: fglloadmatrixd
#endif
#ifdef F_LANG
real(kind=gldouble), dimension(:,:) INTENT_IN m
#else
real(kind=gldouble), dimension(4,DEFER_DIM) INTENT_IN m
#endif
end subroutine fglloadmatrixd
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglloadmatrixf"
#endif
subroutine fglloadmatrixf(m)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: m
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglloadmatrixf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglloadmatrixf@4' :: fglloadmatrixf
#endif
#ifdef F_LANG
real(kind=glfloat), dimension(:,:) INTENT_IN m
#else
real(kind=glfloat), dimension(4,DEFER_DIM) INTENT_IN m
#endif
end subroutine fglloadmatrixf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglloadname"
#endif
subroutine fglloadname(name)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: name
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglloadname
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglloadname@4' :: fglloadname
#endif
integer(kind=gluint) INTENT_IN name
end subroutine fglloadname
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgllogicop"
#endif
subroutine fgllogicop(opcode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: opcode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgllogicop
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgllogicop@4' :: fgllogicop
#endif
integer(kind=glenum) INTENT_IN opcode
end subroutine fgllogicop
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmap1d"
#endif
subroutine fglmap1d(xtarget,u1,u2,stride,order,points)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,u1,u2,stride,order,points
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmap1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmap1d@24' :: fglmap1d
#endif
integer(kind=glenum) INTENT_IN xtarget
real(kind=gldouble) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN stride, order
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN points
end subroutine fglmap1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmap1f"
#endif
subroutine fglmap1f(xtarget,u1,u2,stride,order,points)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,u1,u2,stride,order,points
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmap1f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmap1f@24' :: fglmap1f
#endif
integer(kind=glenum) INTENT_IN xtarget
real(kind=glfloat) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN stride, order
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN points
end subroutine fglmap1f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmap2d"
#endif
subroutine fglmap2d(xtarget,u1,u2,ustride,uorder,v1,v2,vstride, &
   vorder,points)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,u1,u2,ustride,uorder,v1,v2,vstride
!MS$ATTRIBUTES REFERENCE :: vorder,points
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmap2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmap2d@40' :: fglmap2d
#endif
integer(kind=glenum) INTENT_IN xtarget
real(kind=gldouble) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN ustride, uorder
real(kind=gldouble) INTENT_IN v1, v2
integer(kind=glint) INTENT_IN vstride, vorder
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN points
end subroutine fglmap2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmap2f"
#endif
subroutine fglmap2f(xtarget,u1,u2,ustride,uorder,v1,v2,vstride, &
   vorder,points)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,u1,u2,ustride,uorder,v1,v2,vstride
!MS$ATTRIBUTES REFERENCE :: vorder,points
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmap2f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmap2f@40' :: fglmap2f
#endif
integer(kind=glenum) INTENT_IN xtarget
real(kind=glfloat) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN ustride, uorder
real(kind=glfloat) INTENT_IN v1, v2
integer(kind=glint) INTENT_IN vstride, vorder
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN points
end subroutine fglmap2f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmapgrid1d"
#endif
subroutine fglmapgrid1d(un,u1,u2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: un,u1,u2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmapgrid1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmapgrid1d@12' :: fglmapgrid1d
#endif
integer(kind=glint) INTENT_IN un
real(kind=gldouble) INTENT_IN u1, u2
end subroutine fglmapgrid1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmapgrid1f"
#endif
subroutine fglmapgrid1f(un,u1,u2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: un,u1,u2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmapgrid1f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmapgrid1f@12' :: fglmapgrid1f
#endif
integer(kind=glint) INTENT_IN un
real(kind=glfloat) INTENT_IN u1, u2
end subroutine fglmapgrid1f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmapgrid2d"
#endif
subroutine fglmapgrid2d(un,u1,u2,vn,v1,v2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: un,u1,u2,vn,v1,v2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmapgrid2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmapgrid2d@24' :: fglmapgrid2d
#endif
integer(kind=glint) INTENT_IN un
real(kind=gldouble) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN vn
real(kind=gldouble) INTENT_IN v1, v2
end subroutine fglmapgrid2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmapgrid2f"
#endif
subroutine fglmapgrid2f(un,u1,u2,vn,v1,v2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: un,u1,u2,vn,v1,v2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmapgrid2f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmapgrid2f@24' :: fglmapgrid2f
#endif
integer(kind=glint) INTENT_IN un
real(kind=glfloat) INTENT_IN u1, u2
integer(kind=glint) INTENT_IN vn
real(kind=glfloat) INTENT_IN v1, v2
end subroutine fglmapgrid2f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmaterialf"
#endif
subroutine fglmaterialf(face,pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: face,pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmaterialf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmaterialf@12' :: fglmaterialf
#endif
integer(kind=glenum) INTENT_IN face, pname
real(kind=glfloat) INTENT_IN param
end subroutine fglmaterialf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmaterialfv"
#endif
subroutine fglmaterialfv(face,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: face,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmaterialfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmaterialfv@12' :: fglmaterialfv
#endif
integer(kind=glenum) INTENT_IN face, pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
end subroutine fglmaterialfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmateriali"
#endif
subroutine fglmateriali(face,pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: face,pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmateriali
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmateriali@12' :: fglmateriali
#endif
integer(kind=glenum) INTENT_IN face, pname
integer(kind=glint) INTENT_IN param
end subroutine fglmateriali
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmaterialiv"
#endif
subroutine fglmaterialiv(face,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: face,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmaterialiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmaterialiv@12' :: fglmaterialiv
#endif
integer(kind=glenum) INTENT_IN face, pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
end subroutine fglmaterialiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmatrixmode"
#endif
subroutine fglmatrixmode(mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmatrixmode
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmatrixmode@4' :: fglmatrixmode
#endif
integer(kind=glenum) INTENT_IN mode
end subroutine fglmatrixmode
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmultmatrixd"
#endif
subroutine fglmultmatrixd(m)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: m
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmultmatrixd
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmultmatrixd@4' :: fglmultmatrixd
#endif
#ifdef F_LANG
real(kind=gldouble), dimension(:,:) INTENT_IN m
#else
real(kind=gldouble), dimension(4,DEFER_DIM) INTENT_IN m
#endif
end subroutine fglmultmatrixd
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglmultmatrixf"
#endif
subroutine fglmultmatrixf(m)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: m
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglmultmatrixf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglmultmatrixf@4' :: fglmultmatrixf
#endif
#ifdef F_LANG
real(kind=glfloat), dimension(:,:) INTENT_IN m
#else
real(kind=glfloat), dimension(4,DEFER_DIM) INTENT_IN m
#endif
end subroutine fglmultmatrixf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglnewlist"
#endif
subroutine fglnewlist(list,mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: list,mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglnewlist
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglnewlist@8' :: fglnewlist
#endif
integer(kind=gluint) INTENT_IN list
integer(kind=glenum) INTENT_IN mode
end subroutine fglnewlist
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglnormal3b"
#endif
subroutine fglnormal3b(nx,ny,nz)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nx,ny,nz
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglnormal3b
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglnormal3b@12' :: fglnormal3b
#endif
integer(kind=glbyte) INTENT_IN nx, ny, nz
end subroutine fglnormal3b
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglnormal3b"
#endif
subroutine f9bglnormal3b(nx,ny,nz)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nx,ny,nz
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglnormal3b
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglnormal3b@12' :: f9bglnormal3b
#endif
integer(kind=glbyte) INTENT_IN nx, ny, nz
end subroutine f9bglnormal3b
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglnormal3b"
#endif
subroutine f9cglnormal3b(nx,ny,nz)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nx,ny,nz
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglnormal3b
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglnormal3b@12' :: f9cglnormal3b
#endif
integer(kind=glbyte) INTENT_IN nx, ny, nz
end subroutine f9cglnormal3b
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglnormal3bv"
#endif
subroutine fglnormal3bv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglnormal3bv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglnormal3bv@4' :: fglnormal3bv
#endif
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglnormal3bv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglnormal3bv"
#endif
subroutine f9bglnormal3bv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglnormal3bv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglnormal3bv@4' :: f9bglnormal3bv
#endif
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9bglnormal3bv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglnormal3bv"
#endif
subroutine f9cglnormal3bv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglnormal3bv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglnormal3bv@4' :: f9cglnormal3bv
#endif
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9cglnormal3bv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglnormal3d"
#endif
subroutine fglnormal3d(nx,ny,nz)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nx,ny,nz
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglnormal3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglnormal3d@12' :: fglnormal3d
#endif
real(kind=gldouble) INTENT_IN nx, ny, nz
end subroutine fglnormal3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglnormal3dv"
#endif
subroutine fglnormal3dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglnormal3dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglnormal3dv@4' :: fglnormal3dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglnormal3dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglnormal3f"
#endif
subroutine fglnormal3f(nx,ny,nz)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nx,ny,nz
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglnormal3f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglnormal3f@12' :: fglnormal3f
#endif
real(kind=glfloat) INTENT_IN nx, ny, nz
end subroutine fglnormal3f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglnormal3fv"
#endif
subroutine fglnormal3fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglnormal3fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglnormal3fv@4' :: fglnormal3fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglnormal3fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglnormal3i"
#endif
subroutine fglnormal3i(nx,ny,nz)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nx,ny,nz
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglnormal3i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglnormal3i@12' :: fglnormal3i
#endif
integer(kind=glint) INTENT_IN nx, ny, nz
end subroutine fglnormal3i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglnormal3iv"
#endif
subroutine fglnormal3iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglnormal3iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglnormal3iv@4' :: fglnormal3iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglnormal3iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglnormal3s"
#endif
subroutine fglnormal3s(nx,ny,nz)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nx,ny,nz
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglnormal3s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglnormal3s@12' :: fglnormal3s
#endif
integer(kind=glshort) INTENT_IN nx, ny, nz
end subroutine fglnormal3s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglnormal3s"
#endif
subroutine f9aglnormal3s(nx,ny,nz)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nx,ny,nz
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglnormal3s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglnormal3s@12' :: f9aglnormal3s
#endif
integer(kind=glshort) INTENT_IN nx, ny, nz
end subroutine f9aglnormal3s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglnormal3sv"
#endif
subroutine fglnormal3sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglnormal3sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglnormal3sv@4' :: fglnormal3sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglnormal3sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglnormal3sv"
#endif
subroutine f9aglnormal3sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglnormal3sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglnormal3sv@4' :: f9aglnormal3sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9aglnormal3sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglortho"
#endif
subroutine fglortho(left,right,bottom,top,near,far)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: left,right,bottom,top,near,far
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglortho
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglortho@24' :: fglortho
#endif
real(kind=gldouble) INTENT_IN left, right, bottom, top, near, far
end subroutine fglortho
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpassthrough"
#endif
subroutine fglpassthrough(token)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: token
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpassthrough
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpassthrough@4' :: fglpassthrough
#endif
real(kind=glfloat) INTENT_IN token
end subroutine fglpassthrough
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpixelmapfv"
#endif
subroutine fglpixelmapfv(map,mapsize,values)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: map,mapsize,values
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpixelmapfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpixelmapfv@12' :: fglpixelmapfv
#endif
integer(kind=glenum) INTENT_IN map
integer(kind=glint) INTENT_IN mapsize
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN values
end subroutine fglpixelmapfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpixelmapuiv"
#endif
subroutine fglpixelmapuiv(map,mapsize,values)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: map,mapsize,values
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpixelmapuiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpixelmapuiv@12' :: fglpixelmapuiv
#endif
integer(kind=glenum) INTENT_IN map
integer(kind=glint) INTENT_IN mapsize
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN values
end subroutine fglpixelmapuiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpixelmapusv"
#endif
subroutine fglpixelmapusv(map,mapsize,values)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: map,mapsize,values
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpixelmapusv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpixelmapusv@12' :: fglpixelmapusv
#endif
integer(kind=glenum) INTENT_IN map
integer(kind=glint) INTENT_IN mapsize
integer(kind=glushort), dimension(DEFER_DIM) INTENT_IN values
end subroutine fglpixelmapusv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglpixelmapusv"
#endif
subroutine f9aglpixelmapusv(map,mapsize,values)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: map,mapsize,values
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglpixelmapusv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglpixelmapusv@12' :: f9aglpixelmapusv
#endif
integer(kind=glenum) INTENT_IN map
integer(kind=glint) INTENT_IN mapsize
integer(kind=glushort), dimension(DEFER_DIM) INTENT_IN values
end subroutine f9aglpixelmapusv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpixelstoref"
#endif
subroutine fglpixelstoref(pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpixelstoref
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpixelstoref@8' :: fglpixelstoref
#endif
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat) INTENT_IN param
end subroutine fglpixelstoref
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpixelstorei"
#endif
subroutine fglpixelstorei(pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpixelstorei
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpixelstorei@8' :: fglpixelstorei
#endif
integer(kind=glenum) INTENT_IN pname
integer(kind=glint) INTENT_IN param
end subroutine fglpixelstorei
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpixeltransferf"
#endif
subroutine fglpixeltransferf(pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpixeltransferf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpixeltransferf@8' :: fglpixeltransferf
#endif
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat) INTENT_IN param
end subroutine fglpixeltransferf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpixeltransferi"
#endif
subroutine fglpixeltransferi(pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpixeltransferi
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpixeltransferi@8' :: fglpixeltransferi
#endif
integer(kind=glenum) INTENT_IN pname
integer(kind=glint) INTENT_IN param
end subroutine fglpixeltransferi
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpixelzoom"
#endif
subroutine fglpixelzoom(xfactor,yfactor)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xfactor,yfactor
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpixelzoom
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpixelzoom@8' :: fglpixelzoom
#endif
real(kind=glfloat) INTENT_IN xfactor, yfactor
end subroutine fglpixelzoom
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpointsize"
#endif
subroutine fglpointsize(xsize)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpointsize
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpointsize@4' :: fglpointsize
#endif
real(kind=glfloat) INTENT_IN xsize
end subroutine fglpointsize
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpolygonmode"
#endif
subroutine fglpolygonmode(face,mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: face,mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpolygonmode
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpolygonmode@8' :: fglpolygonmode
#endif
integer(kind=glenum) INTENT_IN face, mode
end subroutine fglpolygonmode
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpolygonstipple"
#endif
subroutine fglpolygonstipple(mask)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mask
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpolygonstipple
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpolygonstipple@4' :: fglpolygonstipple
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN mask
end subroutine fglpolygonstipple
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglpolygonstipple"
#endif
subroutine f9bglpolygonstipple(mask)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mask
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglpolygonstipple
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglpolygonstipple@4' :: f9bglpolygonstipple
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN mask
end subroutine f9bglpolygonstipple
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglpolygonstipple"
#endif
subroutine f9cglpolygonstipple(mask)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mask
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglpolygonstipple
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglpolygonstipple@4' :: f9cglpolygonstipple
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN mask
end subroutine f9cglpolygonstipple
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpopattrib"
#endif
subroutine fglpopattrib()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpopattrib
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpopattrib@0' :: fglpopattrib
#endif
end subroutine fglpopattrib
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpopmatrix"
#endif
subroutine fglpopmatrix()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpopmatrix
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpopmatrix@0' :: fglpopmatrix
#endif
end subroutine fglpopmatrix
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpopname"
#endif
subroutine fglpopname()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpopname
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpopname@0' :: fglpopname
#endif
end subroutine fglpopname
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpushattrib"
#endif
subroutine fglpushattrib(mask)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mask
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpushattrib
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpushattrib@4' :: fglpushattrib
#endif
integer(kind=glbitfield) INTENT_IN mask
end subroutine fglpushattrib
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpushmatrix"
#endif
subroutine fglpushmatrix()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpushmatrix
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpushmatrix@0' :: fglpushmatrix
#endif
end subroutine fglpushmatrix
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpushname"
#endif
subroutine fglpushname(name)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: name
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpushname
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpushname@4' :: fglpushname
#endif
integer(kind=gluint) INTENT_IN name
end subroutine fglpushname
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos2d"
#endif
subroutine fglrasterpos2d(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos2d@8' :: fglrasterpos2d
#endif
real(kind=gldouble) INTENT_IN x, y
end subroutine fglrasterpos2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos2dv"
#endif
subroutine fglrasterpos2dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos2dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos2dv@4' :: fglrasterpos2dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos2dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos2f"
#endif
subroutine fglrasterpos2f(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos2f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos2f@8' :: fglrasterpos2f
#endif
real(kind=glfloat) INTENT_IN x, y
end subroutine fglrasterpos2f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos2fv"
#endif
subroutine fglrasterpos2fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos2fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos2fv@4' :: fglrasterpos2fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos2fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos2i"
#endif
subroutine fglrasterpos2i(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos2i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos2i@8' :: fglrasterpos2i
#endif
integer(kind=glint) INTENT_IN x, y
end subroutine fglrasterpos2i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos2iv"
#endif
subroutine fglrasterpos2iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos2iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos2iv@4' :: fglrasterpos2iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos2iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos2s"
#endif
subroutine fglrasterpos2s(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos2s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos2s@8' :: fglrasterpos2s
#endif
integer(kind=glshort) INTENT_IN x, y
end subroutine fglrasterpos2s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglrasterpos2s"
#endif
subroutine f9aglrasterpos2s(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglrasterpos2s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglrasterpos2s@8' :: f9aglrasterpos2s
#endif
integer(kind=glshort) INTENT_IN x, y
end subroutine f9aglrasterpos2s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos2sv"
#endif
subroutine fglrasterpos2sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos2sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos2sv@4' :: fglrasterpos2sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos2sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglrasterpos2sv"
#endif
subroutine f9aglrasterpos2sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglrasterpos2sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglrasterpos2sv@4' :: f9aglrasterpos2sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9aglrasterpos2sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos3d"
#endif
subroutine fglrasterpos3d(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos3d@12' :: fglrasterpos3d
#endif
real(kind=gldouble) INTENT_IN x, y, z
end subroutine fglrasterpos3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos3dv"
#endif
subroutine fglrasterpos3dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos3dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos3dv@4' :: fglrasterpos3dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos3dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos3f"
#endif
subroutine fglrasterpos3f(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos3f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos3f@12' :: fglrasterpos3f
#endif
real(kind=glfloat) INTENT_IN x, y, z
end subroutine fglrasterpos3f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos3fv"
#endif
subroutine fglrasterpos3fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos3fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos3fv@4' :: fglrasterpos3fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos3fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos3i"
#endif
subroutine fglrasterpos3i(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos3i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos3i@12' :: fglrasterpos3i
#endif
integer(kind=glint) INTENT_IN x, y, z
end subroutine fglrasterpos3i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos3iv"
#endif
subroutine fglrasterpos3iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos3iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos3iv@4' :: fglrasterpos3iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos3iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos3s"
#endif
subroutine fglrasterpos3s(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos3s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos3s@12' :: fglrasterpos3s
#endif
integer(kind=glshort) INTENT_IN x, y, z
end subroutine fglrasterpos3s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglrasterpos3s"
#endif
subroutine f9aglrasterpos3s(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglrasterpos3s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglrasterpos3s@12' :: f9aglrasterpos3s
#endif
integer(kind=glshort) INTENT_IN x, y, z
end subroutine f9aglrasterpos3s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos3sv"
#endif
subroutine fglrasterpos3sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos3sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos3sv@4' :: fglrasterpos3sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos3sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglrasterpos3sv"
#endif
subroutine f9aglrasterpos3sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglrasterpos3sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglrasterpos3sv@4' :: f9aglrasterpos3sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9aglrasterpos3sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos4d"
#endif
subroutine fglrasterpos4d(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos4d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos4d@16' :: fglrasterpos4d
#endif
real(kind=gldouble) INTENT_IN x, y, z, w
end subroutine fglrasterpos4d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos4dv"
#endif
subroutine fglrasterpos4dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos4dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos4dv@4' :: fglrasterpos4dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos4dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos4f"
#endif
subroutine fglrasterpos4f(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos4f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos4f@16' :: fglrasterpos4f
#endif
real(kind=glfloat) INTENT_IN x, y, z, w
end subroutine fglrasterpos4f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos4fv"
#endif
subroutine fglrasterpos4fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos4fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos4fv@4' :: fglrasterpos4fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos4fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos4i"
#endif
subroutine fglrasterpos4i(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos4i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos4i@16' :: fglrasterpos4i
#endif
integer(kind=glint) INTENT_IN x, y, z, w
end subroutine fglrasterpos4i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos4iv"
#endif
subroutine fglrasterpos4iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos4iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos4iv@4' :: fglrasterpos4iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos4iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos4s"
#endif
subroutine fglrasterpos4s(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos4s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos4s@16' :: fglrasterpos4s
#endif
integer(kind=glshort) INTENT_IN x, y, z, w
end subroutine fglrasterpos4s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglrasterpos4s"
#endif
subroutine f9aglrasterpos4s(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglrasterpos4s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglrasterpos4s@16' :: f9aglrasterpos4s
#endif
integer(kind=glshort) INTENT_IN x, y, z, w
end subroutine f9aglrasterpos4s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrasterpos4sv"
#endif
subroutine fglrasterpos4sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrasterpos4sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrasterpos4sv@4' :: fglrasterpos4sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglrasterpos4sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglrasterpos4sv"
#endif
subroutine f9aglrasterpos4sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglrasterpos4sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglrasterpos4sv@4' :: f9aglrasterpos4sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9aglrasterpos4sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglreadbuffer"
#endif
subroutine fglreadbuffer(mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglreadbuffer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglreadbuffer@4' :: fglreadbuffer
#endif
integer(kind=glenum) INTENT_IN mode
end subroutine fglreadbuffer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glreadpixels"
#endif
subroutine f9y0glreadpixels(x,y,width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glreadpixels
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glreadpixels@32' :: f9y0glreadpixels
#endif
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0glreadpixels
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glreadpixels"
#endif
subroutine f9y1glreadpixels(x,y,width,height,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,width,height,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glreadpixels
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glreadpixels@28' :: f9y1glreadpixels
#endif
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: pixels
end subroutine f9y1glreadpixels
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glreadpixels"
#endif
subroutine f9y4glreadpixels(x,y,width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glreadpixels
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glreadpixels@32' :: f9y4glreadpixels
#endif
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4glreadpixels
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glreadpixels"
#endif
subroutine f9y5glreadpixels(x,y,width,height,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,width,height,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glreadpixels
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glreadpixels@28' :: f9y5glreadpixels
#endif
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: pixels
end subroutine f9y5glreadpixels
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glreadpixels"
#endif
subroutine f9y6glreadpixels(x,y,width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glreadpixels
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glreadpixels@32' :: f9y6glreadpixels
#endif
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM), intent(out) :: pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6glreadpixels
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrectd"
#endif
subroutine fglrectd(x1,y1,x2,y2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x1,y1,x2,y2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrectd
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrectd@16' :: fglrectd
#endif
real(kind=gldouble) INTENT_IN x1, y1, x2, y2
end subroutine fglrectd
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrectdv"
#endif
subroutine fglrectdv(v1,v2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v1,v2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrectdv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrectdv@8' :: fglrectdv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v1, v2
end subroutine fglrectdv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrectf"
#endif
subroutine fglrectf(x1,y1,x2,y2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x1,y1,x2,y2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrectf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrectf@16' :: fglrectf
#endif
real(kind=glfloat) INTENT_IN x1, y1, x2, y2
end subroutine fglrectf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrectfv"
#endif
subroutine fglrectfv(v1,v2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v1,v2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrectfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrectfv@8' :: fglrectfv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v1, v2
end subroutine fglrectfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrecti"
#endif
subroutine fglrecti(x1,y1,x2,y2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x1,y1,x2,y2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrecti
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrecti@16' :: fglrecti
#endif
integer(kind=glint) INTENT_IN x1, y1, x2, y2
end subroutine fglrecti
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrectiv"
#endif
subroutine fglrectiv(v1,v2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v1,v2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrectiv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrectiv@8' :: fglrectiv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v1, v2
end subroutine fglrectiv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrects"
#endif
subroutine fglrects(x1,y1,x2,y2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x1,y1,x2,y2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrects
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrects@16' :: fglrects
#endif
integer(kind=glshort) INTENT_IN x1, y1, x2, y2
end subroutine fglrects
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglrects"
#endif
subroutine f9aglrects(x1,y1,x2,y2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x1,y1,x2,y2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglrects
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglrects@16' :: f9aglrects
#endif
integer(kind=glshort) INTENT_IN x1, y1, x2, y2
end subroutine f9aglrects
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrectsv"
#endif
subroutine fglrectsv(v1,v2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v1,v2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrectsv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrectsv@8' :: fglrectsv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v1, v2
end subroutine fglrectsv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglrectsv"
#endif
subroutine f9aglrectsv(v1,v2)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v1,v2
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglrectsv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglrectsv@8' :: f9aglrectsv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v1, v2
end subroutine f9aglrectsv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrendermode"
#endif
function fglrendermode(mode) result(resfglrendermode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrendermode
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrendermode@4' :: fglrendermode
#endif
integer(kind=glenum) INTENT_IN mode
integer(kind=glint) :: resfglrendermode
end function fglrendermode
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrotated"
#endif
subroutine fglrotated(angle,x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: angle,x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrotated
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrotated@16' :: fglrotated
#endif
real(kind=gldouble) INTENT_IN angle, x, y, z
end subroutine fglrotated
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglrotatef"
#endif
subroutine fglrotatef(angle,x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: angle,x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglrotatef
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglrotatef@16' :: fglrotatef
#endif
real(kind=glfloat) INTENT_IN angle, x, y, z
end subroutine fglrotatef
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglscaled"
#endif
subroutine fglscaled(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglscaled
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglscaled@12' :: fglscaled
#endif
real(kind=gldouble) INTENT_IN x, y, z
end subroutine fglscaled
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglscalef"
#endif
subroutine fglscalef(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglscalef
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglscalef@12' :: fglscalef
#endif
real(kind=glfloat) INTENT_IN x, y, z
end subroutine fglscalef
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglscissor"
#endif
subroutine fglscissor(x,y,width,height)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,width,height
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglscissor
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglscissor@16' :: fglscissor
#endif
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
end subroutine fglscissor
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglselectbuffer"
#endif
subroutine fglselectbuffer(xsize,buffer)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,buffer
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglselectbuffer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglselectbuffer@8' :: fglselectbuffer
#endif
integer(kind=glsizei) INTENT_IN xsize
integer(kind=gluint), dimension(DEFER_DIM), TARGET intent(in out) :: buffer
end subroutine fglselectbuffer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglshademodel"
#endif
subroutine fglshademodel(mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglshademodel
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglshademodel@4' :: fglshademodel
#endif
integer(kind=glenum) INTENT_IN mode
end subroutine fglshademodel
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglstencilfunc"
#endif
subroutine fglstencilfunc(func,ref,mask)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: func,ref,mask
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglstencilfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglstencilfunc@12' :: fglstencilfunc
#endif
integer(kind=glenum) INTENT_IN func
integer(kind=glint) INTENT_IN ref
integer(kind=gluint) INTENT_IN mask
end subroutine fglstencilfunc
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglstencilmask"
#endif
subroutine fglstencilmask(mask)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mask
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglstencilmask
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglstencilmask@4' :: fglstencilmask
#endif
integer(kind=gluint) INTENT_IN mask
end subroutine fglstencilmask
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglstencilop"
#endif
subroutine fglstencilop(fail,zfail,zpass)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: fail,zfail,zpass
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglstencilop
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglstencilop@12' :: fglstencilop
#endif
integer(kind=glenum) INTENT_IN fail, zfail, zpass
end subroutine fglstencilop
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord1d"
#endif
subroutine fgltexcoord1d(s)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord1d@4' :: fgltexcoord1d
#endif
real(kind=gldouble) INTENT_IN s
end subroutine fgltexcoord1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord1dv"
#endif
subroutine fgltexcoord1dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord1dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord1dv@4' :: fgltexcoord1dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord1dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord1f"
#endif
subroutine fgltexcoord1f(s)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord1f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord1f@4' :: fgltexcoord1f
#endif
real(kind=glfloat) INTENT_IN s
end subroutine fgltexcoord1f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord1fv"
#endif
subroutine fgltexcoord1fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord1fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord1fv@4' :: fgltexcoord1fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord1fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord1i"
#endif
subroutine fgltexcoord1i(s)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord1i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord1i@4' :: fgltexcoord1i
#endif
integer(kind=glint) INTENT_IN s
end subroutine fgltexcoord1i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord1iv"
#endif
subroutine fgltexcoord1iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord1iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord1iv@4' :: fgltexcoord1iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord1iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord1s"
#endif
subroutine fgltexcoord1s(s)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord1s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord1s@4' :: fgltexcoord1s
#endif
integer(kind=glshort) INTENT_IN s
end subroutine fgltexcoord1s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9agltexcoord1s"
#endif
subroutine f9agltexcoord1s(s)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9agltexcoord1s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9agltexcoord1s@4' :: f9agltexcoord1s
#endif
integer(kind=glshort) INTENT_IN s
end subroutine f9agltexcoord1s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord1sv"
#endif
subroutine fgltexcoord1sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord1sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord1sv@4' :: fgltexcoord1sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord1sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9agltexcoord1sv"
#endif
subroutine f9agltexcoord1sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9agltexcoord1sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9agltexcoord1sv@4' :: f9agltexcoord1sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9agltexcoord1sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord2d"
#endif
subroutine fgltexcoord2d(s,t)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord2d@8' :: fgltexcoord2d
#endif
real(kind=gldouble) INTENT_IN s, t
end subroutine fgltexcoord2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord2dv"
#endif
subroutine fgltexcoord2dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord2dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord2dv@4' :: fgltexcoord2dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord2dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord2f"
#endif
subroutine fgltexcoord2f(s,t)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord2f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord2f@8' :: fgltexcoord2f
#endif
real(kind=glfloat) INTENT_IN s, t
end subroutine fgltexcoord2f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord2fv"
#endif
subroutine fgltexcoord2fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord2fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord2fv@4' :: fgltexcoord2fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord2fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord2i"
#endif
subroutine fgltexcoord2i(s,t)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord2i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord2i@8' :: fgltexcoord2i
#endif
integer(kind=glint) INTENT_IN s, t
end subroutine fgltexcoord2i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord2iv"
#endif
subroutine fgltexcoord2iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord2iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord2iv@4' :: fgltexcoord2iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord2iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord2s"
#endif
subroutine fgltexcoord2s(s,t)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord2s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord2s@8' :: fgltexcoord2s
#endif
integer(kind=glshort) INTENT_IN s, t
end subroutine fgltexcoord2s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9agltexcoord2s"
#endif
subroutine f9agltexcoord2s(s,t)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9agltexcoord2s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9agltexcoord2s@8' :: f9agltexcoord2s
#endif
integer(kind=glshort) INTENT_IN s, t
end subroutine f9agltexcoord2s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord2sv"
#endif
subroutine fgltexcoord2sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord2sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord2sv@4' :: fgltexcoord2sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord2sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9agltexcoord2sv"
#endif
subroutine f9agltexcoord2sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9agltexcoord2sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9agltexcoord2sv@4' :: f9agltexcoord2sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9agltexcoord2sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord3d"
#endif
subroutine fgltexcoord3d(s,t,r)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t,r
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord3d@12' :: fgltexcoord3d
#endif
real(kind=gldouble) INTENT_IN s, t, r
end subroutine fgltexcoord3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord3dv"
#endif
subroutine fgltexcoord3dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord3dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord3dv@4' :: fgltexcoord3dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord3dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord3f"
#endif
subroutine fgltexcoord3f(s,t,r)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t,r
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord3f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord3f@12' :: fgltexcoord3f
#endif
real(kind=glfloat) INTENT_IN s, t, r
end subroutine fgltexcoord3f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord3fv"
#endif
subroutine fgltexcoord3fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord3fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord3fv@4' :: fgltexcoord3fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord3fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord3i"
#endif
subroutine fgltexcoord3i(s,t,r)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t,r
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord3i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord3i@12' :: fgltexcoord3i
#endif
integer(kind=glint) INTENT_IN s, t, r
end subroutine fgltexcoord3i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord3iv"
#endif
subroutine fgltexcoord3iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord3iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord3iv@4' :: fgltexcoord3iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord3iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord3s"
#endif
subroutine fgltexcoord3s(s,t,r)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t,r
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord3s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord3s@12' :: fgltexcoord3s
#endif
integer(kind=glshort) INTENT_IN s, t, r
end subroutine fgltexcoord3s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9agltexcoord3s"
#endif
subroutine f9agltexcoord3s(s,t,r)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t,r
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9agltexcoord3s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9agltexcoord3s@12' :: f9agltexcoord3s
#endif
integer(kind=glshort) INTENT_IN s, t, r
end subroutine f9agltexcoord3s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord3sv"
#endif
subroutine fgltexcoord3sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord3sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord3sv@4' :: fgltexcoord3sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord3sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9agltexcoord3sv"
#endif
subroutine f9agltexcoord3sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9agltexcoord3sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9agltexcoord3sv@4' :: f9agltexcoord3sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9agltexcoord3sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord4d"
#endif
subroutine fgltexcoord4d(s,t,r,q)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t,r,q
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord4d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord4d@16' :: fgltexcoord4d
#endif
real(kind=gldouble) INTENT_IN s, t, r, q
end subroutine fgltexcoord4d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord4dv"
#endif
subroutine fgltexcoord4dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord4dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord4dv@4' :: fgltexcoord4dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord4dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord4f"
#endif
subroutine fgltexcoord4f(s,t,r,q)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t,r,q
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord4f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord4f@16' :: fgltexcoord4f
#endif
real(kind=glfloat) INTENT_IN s, t, r, q
end subroutine fgltexcoord4f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord4fv"
#endif
subroutine fgltexcoord4fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord4fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord4fv@4' :: fgltexcoord4fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord4fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord4i"
#endif
subroutine fgltexcoord4i(s,t,r,q)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t,r,q
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord4i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord4i@16' :: fgltexcoord4i
#endif
integer(kind=glint) INTENT_IN s, t, r, q
end subroutine fgltexcoord4i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord4iv"
#endif
subroutine fgltexcoord4iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord4iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord4iv@4' :: fgltexcoord4iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord4iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord4s"
#endif
subroutine fgltexcoord4s(s,t,r,q)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t,r,q
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord4s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord4s@16' :: fgltexcoord4s
#endif
integer(kind=glshort) INTENT_IN s, t, r, q
end subroutine fgltexcoord4s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9agltexcoord4s"
#endif
subroutine f9agltexcoord4s(s,t,r,q)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: s,t,r,q
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9agltexcoord4s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9agltexcoord4s@16' :: f9agltexcoord4s
#endif
integer(kind=glshort) INTENT_IN s, t, r, q
end subroutine f9agltexcoord4s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexcoord4sv"
#endif
subroutine fgltexcoord4sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexcoord4sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexcoord4sv@4' :: fgltexcoord4sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fgltexcoord4sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9agltexcoord4sv"
#endif
subroutine f9agltexcoord4sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9agltexcoord4sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9agltexcoord4sv@4' :: f9agltexcoord4sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9agltexcoord4sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexenvf"
#endif
subroutine fgltexenvf(xtarget,pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexenvf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexenvf@12' :: fgltexenvf
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat) INTENT_IN param
end subroutine fgltexenvf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexenvfv"
#endif
subroutine fgltexenvfv(xtarget,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexenvfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexenvfv@12' :: fgltexenvfv
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
end subroutine fgltexenvfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexenvi"
#endif
subroutine fgltexenvi(xtarget,pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexenvi
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexenvi@12' :: fgltexenvi
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint) INTENT_IN param
end subroutine fgltexenvi
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexenviv"
#endif
subroutine fgltexenviv(xtarget,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexenviv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexenviv@12' :: fgltexenviv
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
end subroutine fgltexenviv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexgend"
#endif
subroutine fgltexgend(coord,pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: coord,pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexgend
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexgend@12' :: fgltexgend
#endif
integer(kind=glenum) INTENT_IN coord, pname
real(kind=gldouble) INTENT_IN param
end subroutine fgltexgend
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexgendv"
#endif
subroutine fgltexgendv(coord,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: coord,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexgendv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexgendv@12' :: fgltexgendv
#endif
integer(kind=glenum) INTENT_IN coord, pname
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN params
end subroutine fgltexgendv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexgenf"
#endif
subroutine fgltexgenf(coord,pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: coord,pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexgenf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexgenf@12' :: fgltexgenf
#endif
integer(kind=glenum) INTENT_IN coord, pname
real(kind=glfloat) INTENT_IN param
end subroutine fgltexgenf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexgenfv"
#endif
subroutine fgltexgenfv(coord,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: coord,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexgenfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexgenfv@12' :: fgltexgenfv
#endif
integer(kind=glenum) INTENT_IN coord, pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
end subroutine fgltexgenfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexgeni"
#endif
subroutine fgltexgeni(coord,pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: coord,pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexgeni
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexgeni@12' :: fgltexgeni
#endif
integer(kind=glenum) INTENT_IN coord, pname
integer(kind=glint) INTENT_IN param
end subroutine fgltexgeni
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexgeniv"
#endif
subroutine fgltexgeniv(coord,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: coord,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexgeniv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexgeniv@12' :: fgltexgeniv
#endif
integer(kind=glenum) INTENT_IN coord, pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
end subroutine fgltexgeniv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glteximage1d"
#endif
subroutine f9y0glteximage1d(xtarget,level,components,width,border,format,xtype, &
   pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,components,width,border
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glteximage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glteximage1d@36' :: f9y0glteximage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0glteximage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glteximage1d"
#endif
subroutine f9y1glteximage1d(xtarget,level,components,width,border,format,xtype, &
   pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,components,width,border
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glteximage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glteximage1d@32' :: f9y1glteximage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y1glteximage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glteximage1d"
#endif
subroutine f9y4glteximage1d(xtarget,level,components,width,border,format,xtype, &
   pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,components,width,border
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glteximage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glteximage1d@36' :: f9y4glteximage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4glteximage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glteximage1d"
#endif
subroutine f9y5glteximage1d(xtarget,level,components,width,border,format,xtype, &
   pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,components,width,border
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glteximage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glteximage1d@32' :: f9y5glteximage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y5glteximage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glteximage1d"
#endif
subroutine f9y6glteximage1d(xtarget,level,components,width,border,format,xtype, &
   pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,components,width,border
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glteximage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glteximage1d@36' :: f9y6glteximage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6glteximage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glteximage2d"
#endif
subroutine f9y0glteximage2d(xtarget,level,components,width,height,border, &
   format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,components,width,height,border
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glteximage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glteximage2d@40' :: f9y0glteximage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0glteximage2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glteximage2d"
#endif
subroutine f9y1glteximage2d(xtarget,level,components,width,height,border, &
   format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,components,width,height,border
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glteximage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glteximage2d@36' :: f9y1glteximage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y1glteximage2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glteximage2d"
#endif
subroutine f9y4glteximage2d(xtarget,level,components,width,height,border, &
   format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,components,width,height,border
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glteximage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glteximage2d@40' :: f9y4glteximage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4glteximage2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glteximage2d"
#endif
subroutine f9y5glteximage2d(xtarget,level,components,width,height,border, &
   format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,components,width,height,border
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glteximage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glteximage2d@36' :: f9y5glteximage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y5glteximage2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glteximage2d"
#endif
subroutine f9y6glteximage2d(xtarget,level,components,width,height,border, &
   format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,components,width,height,border
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glteximage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glteximage2d@40' :: f9y6glteximage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, components
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6glteximage2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexparameterf"
#endif
subroutine fgltexparameterf(xtarget,pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexparameterf
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexparameterf@12' :: fgltexparameterf
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat) INTENT_IN param
end subroutine fgltexparameterf
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexparameterfv"
#endif
subroutine fgltexparameterfv(xtarget,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexparameterfv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexparameterfv@12' :: fgltexparameterfv
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
end subroutine fgltexparameterfv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexparameteri"
#endif
subroutine fgltexparameteri(xtarget,pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexparameteri
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexparameteri@12' :: fgltexparameteri
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint) INTENT_IN param
end subroutine fgltexparameteri
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltexparameteriv"
#endif
subroutine fgltexparameteriv(xtarget,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltexparameteriv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltexparameteriv@12' :: fgltexparameteriv
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN params
end subroutine fgltexparameteriv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltranslated"
#endif
subroutine fgltranslated(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltranslated
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltranslated@12' :: fgltranslated
#endif
real(kind=gldouble) INTENT_IN x, y, z
end subroutine fgltranslated
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgltranslatef"
#endif
subroutine fgltranslatef(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgltranslatef
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgltranslatef@12' :: fgltranslatef
#endif
real(kind=glfloat) INTENT_IN x, y, z
end subroutine fgltranslatef
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex2d"
#endif
subroutine fglvertex2d(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex2d@8' :: fglvertex2d
#endif
real(kind=gldouble) INTENT_IN x, y
end subroutine fglvertex2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex2dv"
#endif
subroutine fglvertex2dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex2dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex2dv@4' :: fglvertex2dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex2dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex2f"
#endif
subroutine fglvertex2f(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex2f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex2f@8' :: fglvertex2f
#endif
real(kind=glfloat) INTENT_IN x, y
end subroutine fglvertex2f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex2fv"
#endif
subroutine fglvertex2fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex2fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex2fv@4' :: fglvertex2fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex2fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex2i"
#endif
subroutine fglvertex2i(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex2i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex2i@8' :: fglvertex2i
#endif
integer(kind=glint) INTENT_IN x, y
end subroutine fglvertex2i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex2iv"
#endif
subroutine fglvertex2iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex2iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex2iv@4' :: fglvertex2iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex2iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex2s"
#endif
subroutine fglvertex2s(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex2s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex2s@8' :: fglvertex2s
#endif
integer(kind=glshort) INTENT_IN x, y
end subroutine fglvertex2s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglvertex2s"
#endif
subroutine f9aglvertex2s(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglvertex2s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglvertex2s@8' :: f9aglvertex2s
#endif
integer(kind=glshort) INTENT_IN x, y
end subroutine f9aglvertex2s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex2sv"
#endif
subroutine fglvertex2sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex2sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex2sv@4' :: fglvertex2sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex2sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglvertex2sv"
#endif
subroutine f9aglvertex2sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglvertex2sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglvertex2sv@4' :: f9aglvertex2sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9aglvertex2sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex3d"
#endif
subroutine fglvertex3d(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex3d@12' :: fglvertex3d
#endif
real(kind=gldouble) INTENT_IN x, y, z
end subroutine fglvertex3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex3dv"
#endif
subroutine fglvertex3dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex3dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex3dv@4' :: fglvertex3dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex3dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex3f"
#endif
subroutine fglvertex3f(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex3f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex3f@12' :: fglvertex3f
#endif
real(kind=glfloat) INTENT_IN x, y, z
end subroutine fglvertex3f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex3fv"
#endif
subroutine fglvertex3fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex3fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex3fv@4' :: fglvertex3fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex3fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex3i"
#endif
subroutine fglvertex3i(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex3i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex3i@12' :: fglvertex3i
#endif
integer(kind=glint) INTENT_IN x, y, z
end subroutine fglvertex3i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex3iv"
#endif
subroutine fglvertex3iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex3iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex3iv@4' :: fglvertex3iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex3iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex3s"
#endif
subroutine fglvertex3s(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex3s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex3s@12' :: fglvertex3s
#endif
integer(kind=glshort) INTENT_IN x, y, z
end subroutine fglvertex3s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglvertex3s"
#endif
subroutine f9aglvertex3s(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglvertex3s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglvertex3s@12' :: f9aglvertex3s
#endif
integer(kind=glshort) INTENT_IN x, y, z
end subroutine f9aglvertex3s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex3sv"
#endif
subroutine fglvertex3sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex3sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex3sv@4' :: fglvertex3sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex3sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglvertex3sv"
#endif
subroutine f9aglvertex3sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglvertex3sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglvertex3sv@4' :: f9aglvertex3sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9aglvertex3sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex4d"
#endif
subroutine fglvertex4d(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex4d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex4d@16' :: fglvertex4d
#endif
real(kind=gldouble) INTENT_IN x, y, z, w
end subroutine fglvertex4d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex4dv"
#endif
subroutine fglvertex4dv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex4dv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex4dv@4' :: fglvertex4dv
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex4dv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex4f"
#endif
subroutine fglvertex4f(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex4f
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex4f@16' :: fglvertex4f
#endif
real(kind=glfloat) INTENT_IN x, y, z, w
end subroutine fglvertex4f
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex4fv"
#endif
subroutine fglvertex4fv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex4fv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex4fv@4' :: fglvertex4fv
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex4fv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex4i"
#endif
subroutine fglvertex4i(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex4i
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex4i@16' :: fglvertex4i
#endif
integer(kind=glint) INTENT_IN x, y, z, w
end subroutine fglvertex4i
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex4iv"
#endif
subroutine fglvertex4iv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex4iv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex4iv@4' :: fglvertex4iv
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex4iv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex4s"
#endif
subroutine fglvertex4s(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex4s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex4s@16' :: fglvertex4s
#endif
integer(kind=glshort) INTENT_IN x, y, z, w
end subroutine fglvertex4s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglvertex4s"
#endif
subroutine f9aglvertex4s(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglvertex4s
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglvertex4s@16' :: f9aglvertex4s
#endif
integer(kind=glshort) INTENT_IN x, y, z, w
end subroutine f9aglvertex4s
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglvertex4sv"
#endif
subroutine fglvertex4sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglvertex4sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglvertex4sv@4' :: fglvertex4sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine fglvertex4sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglvertex4sv"
#endif
subroutine f9aglvertex4sv(v)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: v
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglvertex4sv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglvertex4sv@4' :: f9aglvertex4sv
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN v
end subroutine f9aglvertex4sv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglviewport"
#endif
subroutine fglviewport(x,y,width,height)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,width,height
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglviewport
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglviewport@16' :: fglviewport
#endif
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
end subroutine fglviewport
end interface

! OpenGL 1.1

#ifdef OPENGL_1_1

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglaretexturesresident"
#endif
#ifdef PURE_FUNCTIONS
subroutine fglaretexturesresident(n,textures, &
   residences,areresident)
#else
function fglaretexturesresident(n,textures, &
   residences) result(resfglaretexturesresident)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,textures
!MS$ATTRIBUTES REFERENCE :: residences
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglaretexturesresident
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglaretexturesresident@12' :: fglaretexturesresident
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
logical(kind=glboolean), dimension(DEFER_DIM), intent(out) :: residences
#ifdef PURE_FUNCTIONS
logical(kind=glboolean), intent(out) :: areresident
#else
logical(kind=glboolean) :: resfglaretexturesresident
#endif
#ifdef PURE_FUNCTIONS
end subroutine fglaretexturesresident
#else
end function fglaretexturesresident
#endif
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9eglaretexturesresident"
#endif
#ifdef PURE_FUNCTIONS
subroutine f9eglaretexturesresident(n,textures, &
   residences,areresident)
#else
function f9eglaretexturesresident(n,textures, &
   residences) result(resf9eglaretexturesresident)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,textures
!MS$ATTRIBUTES REFERENCE :: residences
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9eglaretexturesresident
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9eglaretexturesresident@12' :: f9eglaretexturesresident
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: residences
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: areresident
#else
integer(kind=glint) :: resf9eglaretexturesresident
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9eglaretexturesresident
#else
end function f9eglaretexturesresident
#endif
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglarrayelement"
#endif
subroutine fglarrayelement(i)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: i
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglarrayelement
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglarrayelement@4' :: fglarrayelement
#endif
integer(kind=glint) INTENT_IN i
end subroutine fglarrayelement
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglbindtexture"
#endif
subroutine fglbindtexture(xtarget,texture)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,texture
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglbindtexture
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglbindtexture@8' :: fglbindtexture
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=gluint) INTENT_IN texture
end subroutine fglbindtexture
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glcolorpointer"
#endif
subroutine f9y0glcolorpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glcolorpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glcolorpointer@16' :: f9y0glcolorpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glint), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y0glcolorpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glcolorpointer"
#endif
subroutine f9y1glcolorpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glcolorpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glcolorpointer@16' :: f9y1glcolorpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y1glcolorpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glcolorpointer"
#endif
subroutine f9y4glcolorpointer(xsize,xtype,stride,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glcolorpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glcolorpointer@20' :: f9y4glcolorpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y4glcolorpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glcolorpointer"
#endif
subroutine f9y5glcolorpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glcolorpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glcolorpointer@16' :: f9y5glcolorpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glbyte), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y5glcolorpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glcolorpointer"
#endif
subroutine f9y6glcolorpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glcolorpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glcolorpointer@16' :: f9y6glcolorpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glshort), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y6glcolorpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y7glcolorpointer"
#endif
subroutine f9y7glcolorpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y7glcolorpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y7glcolorpointer@16' :: f9y7glcolorpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=gldouble), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y7glcolorpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcopyteximage1d"
#endif
subroutine fglcopyteximage1d(xtarget,level,internalformat,x,y, &
   width,border)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,x,y
!MS$ATTRIBUTES REFERENCE :: width,border
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcopyteximage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcopyteximage1d@28' :: fglcopyteximage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width
integer(kind=glint) INTENT_IN border
end subroutine fglcopyteximage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcopyteximage2d"
#endif
subroutine fglcopyteximage2d(xtarget,level,internalformat,x,y, &
   width,height,border)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,x,y
!MS$ATTRIBUTES REFERENCE :: width,height,border
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcopyteximage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcopyteximage2d@32' :: fglcopyteximage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glint) INTENT_IN x, y
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glint) INTENT_IN border
end subroutine fglcopyteximage2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcopytexsubimage1d"
#endif
subroutine fglcopytexsubimage1d(xtarget,level,xoffset,x,y, &
   width)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,x,y
!MS$ATTRIBUTES REFERENCE :: width
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcopytexsubimage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcopytexsubimage1d@24' :: fglcopytexsubimage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, x, y
integer(kind=glsizei) INTENT_IN width
end subroutine fglcopytexsubimage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcopytexsubimage2d"
#endif
subroutine fglcopytexsubimage2d(xtarget,level,xoffset,yoffset,x, &
   y,width,height)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset,x
!MS$ATTRIBUTES REFERENCE :: y,width,height
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcopytexsubimage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcopytexsubimage2d@32' :: fglcopytexsubimage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, x, y
integer(kind=glsizei) INTENT_IN width, height
end subroutine fglcopytexsubimage2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgldeletetextures"
#endif
subroutine fgldeletetextures(n,textures)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,textures
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgldeletetextures
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgldeletetextures@8' :: fgldeletetextures
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
end subroutine fgldeletetextures
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgldisableclientstate"
#endif
subroutine fgldisableclientstate(cap)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: cap
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgldisableclientstate
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgldisableclientstate@4' :: fgldisableclientstate
#endif
integer(kind=glenum) INTENT_IN cap
end subroutine fgldisableclientstate
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgldrawarrays"
#endif
subroutine fgldrawarrays(mode,first,xcount)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,first,xcount
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgldrawarrays
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgldrawarrays@12' :: fgldrawarrays
#endif
integer(kind=glenum) INTENT_IN mode
integer(kind=glint) INTENT_IN first
integer(kind=glsizei) INTENT_IN xcount
end subroutine fgldrawarrays
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0gldrawelements"
#endif
subroutine f9y0gldrawelements(mode,xcount,xtype,indices)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,xcount,xtype,indices
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gldrawelements
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gldrawelements@16' :: f9y0gldrawelements
#endif
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN indices
end subroutine f9y0gldrawelements
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4gldrawelements"
#endif
subroutine f9y4gldrawelements(mode,xcount,xtype,indices,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,xcount,xtype,indices,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gldrawelements
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gldrawelements@20' :: f9y4gldrawelements
#endif
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=glsizei) INTENT_IN xcount
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN indices
integer(kind=glint) INTENT_IN length
end subroutine f9y4gldrawelements
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5gldrawelements"
#endif
subroutine f9y5gldrawelements(mode,xcount,xtype,indices)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,xcount,xtype,indices
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5gldrawelements
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5gldrawelements@16' :: f9y5gldrawelements
#endif
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN indices
end subroutine f9y5gldrawelements
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6gldrawelements"
#endif
subroutine f9y6gldrawelements(mode,xcount,xtype,indices)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,xcount,xtype,indices
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gldrawelements
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gldrawelements@16' :: f9y6gldrawelements
#endif
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN indices
end subroutine f9y6gldrawelements
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgledgeflagpointer"
#endif
subroutine fgledgeflagpointer(stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgledgeflagpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgledgeflagpointer@8' :: fgledgeflagpointer
#endif
integer(kind=glsizei) INTENT_IN stride
logical(kind=glboolean), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine fgledgeflagpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9egledgeflagpointer"
#endif
subroutine f9egledgeflagpointer(stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9egledgeflagpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9egledgeflagpointer@8' :: f9egledgeflagpointer
#endif
integer(kind=glsizei) INTENT_IN stride
integer(kind=glbyte), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9egledgeflagpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglenableclientstate"
#endif
subroutine fglenableclientstate(cap)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: cap
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglenableclientstate
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglenableclientstate@4' :: fglenableclientstate
#endif
integer(kind=glenum) INTENT_IN cap
end subroutine fglenableclientstate
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgentextures"
#endif
subroutine fglgentextures(n,textures)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,textures
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgentextures
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgentextures@8' :: fglgentextures
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: textures
end subroutine fglgentextures
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetpointerv"
#endif
subroutine fglgetpointerv(pname,params,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetpointerv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetpointerv@12' :: fglgetpointerv
#endif
integer(kind=glenum) INTENT_IN pname
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: params
integer(kind=glint) INTENT_IN length
end subroutine fglgetpointerv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glindexpointer"
#endif
subroutine f9y0glindexpointer(xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glindexpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glindexpointer@12' :: f9y0glindexpointer
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glint), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y0glindexpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glindexpointer"
#endif
subroutine f9y1glindexpointer(xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glindexpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glindexpointer@12' :: f9y1glindexpointer
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y1glindexpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glindexpointer"
#endif
subroutine f9y4glindexpointer(xtype,stride,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glindexpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glindexpointer@16' :: f9y4glindexpointer
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y4glindexpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glindexpointer"
#endif
subroutine f9y6glindexpointer(xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glindexpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glindexpointer@12' :: f9y6glindexpointer
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glshort), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y6glindexpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y7glindexpointer"
#endif
subroutine f9y7glindexpointer(xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y7glindexpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y7glindexpointer@12' :: f9y7glindexpointer
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=gldouble), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y7glindexpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglindexub"
#endif
subroutine fglindexub(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglindexub
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglindexub@4' :: fglindexub
#endif
integer(kind=glubyte) INTENT_IN c
end subroutine fglindexub
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglindexub"
#endif
subroutine f9bglindexub(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglindexub
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglindexub@4' :: f9bglindexub
#endif
integer(kind=glubyte) INTENT_IN c
end subroutine f9bglindexub
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglindexub"
#endif
subroutine f9cglindexub(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglindexub
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglindexub@2' :: f9cglindexub
#endif
integer(kind=glubyte) INTENT_IN c
end subroutine f9cglindexub
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglindexubv"
#endif
subroutine fglindexubv(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglindexubv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglindexubv@4' :: fglindexubv
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN c
end subroutine fglindexubv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9bglindexubv"
#endif
subroutine f9bglindexubv(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9bglindexubv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9bglindexubv@4' :: f9bglindexubv
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN c
end subroutine f9bglindexubv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9cglindexubv"
#endif
subroutine f9cglindexubv(c)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: c
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9cglindexubv
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9cglindexubv@4' :: f9cglindexubv
#endif
integer(kind=glubyte), dimension(DEFER_DIM) INTENT_IN c
end subroutine f9cglindexubv
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glinterleavedarrays"
#endif
subroutine f9y0glinterleavedarrays(format,stride,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,stride,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glinterleavedarrays
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glinterleavedarrays@16' :: f9y0glinterleavedarrays
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN stride
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y0glinterleavedarrays
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glinterleavedarrays"
#endif
subroutine f9y1glinterleavedarrays(format,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glinterleavedarrays
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glinterleavedarrays@12' :: f9y1glinterleavedarrays
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y1glinterleavedarrays
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glinterleavedarrays"
#endif
subroutine f9y4glinterleavedarrays(format,stride,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,stride,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glinterleavedarrays
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glinterleavedarrays@16' :: f9y4glinterleavedarrays
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN stride
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y4glinterleavedarrays
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glinterleavedarrays"
#endif
subroutine f9y5glinterleavedarrays(format,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glinterleavedarrays
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glinterleavedarrays@12' :: f9y5glinterleavedarrays
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN stride
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y5glinterleavedarrays
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glinterleavedarrays"
#endif
subroutine f9y6glinterleavedarrays(format,stride,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,stride,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glinterleavedarrays
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glinterleavedarrays@16' :: f9y6glinterleavedarrays
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN stride
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y6glinterleavedarrays
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglistexture"
#endif
function fglistexture(texture) result(resfglistexture)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: texture
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglistexture
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglistexture@4' :: fglistexture
#endif
integer(kind=gluint) INTENT_IN texture
logical(kind=glboolean) :: resfglistexture
end function fglistexture
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9eglistexture"
#endif
function f9eglistexture(texture) result(resf9eglistexture)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: texture
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9eglistexture
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9eglistexture@4' :: f9eglistexture
#endif
integer(kind=gluint) INTENT_IN texture
integer(kind=glint) :: resf9eglistexture
end function f9eglistexture
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glnormalpointer"
#endif
subroutine f9y0glnormalpointer(xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glnormalpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glnormalpointer@12' :: f9y0glnormalpointer
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y0glnormalpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glnormalpointer"
#endif
subroutine f9y1glnormalpointer(xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glnormalpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glnormalpointer@12' :: f9y1glnormalpointer
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y1glnormalpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glnormalpointer"
#endif
subroutine f9y4glnormalpointer(xtype,stride,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glnormalpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glnormalpointer@16' :: f9y4glnormalpointer
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y4glnormalpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glnormalpointer"
#endif
subroutine f9y5glnormalpointer(xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glnormalpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glnormalpointer@12' :: f9y5glnormalpointer
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y5glnormalpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glnormalpointer"
#endif
subroutine f9y6glnormalpointer(xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glnormalpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glnormalpointer@12' :: f9y6glnormalpointer
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y6glnormalpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y7glnormalpointer"
#endif
subroutine f9y7glnormalpointer(xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y7glnormalpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y7glnormalpointer@12' :: f9y7glnormalpointer
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y7glnormalpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpolygonoffset"
#endif
subroutine fglpolygonoffset(factor,units)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: factor,units
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpolygonoffset
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpolygonoffset@8' :: fglpolygonoffset
#endif
real(kind=glfloat) INTENT_IN factor, units
end subroutine fglpolygonoffset
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpopclientattrib"
#endif
subroutine fglpopclientattrib()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpopclientattrib
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpopclientattrib@0' :: fglpopclientattrib
#endif
end subroutine fglpopclientattrib
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglprioritizetextures"
#endif
subroutine fglprioritizetextures(n,textures,priorities)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,textures,priorities
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglprioritizetextures
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglprioritizetextures@12' :: fglprioritizetextures
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
real(kind=glclampf), dimension(DEFER_DIM) INTENT_IN priorities
end subroutine fglprioritizetextures
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpushclientattrib"
#endif
subroutine fglpushclientattrib(mask)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mask
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpushclientattrib
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpushclientattrib@4' :: fglpushclientattrib
#endif
integer(kind=glbitfield) INTENT_IN mask
end subroutine fglpushclientattrib
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0gltexcoordpointer"
#endif
subroutine f9y0gltexcoordpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gltexcoordpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gltexcoordpointer@16' :: f9y0gltexcoordpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y0gltexcoordpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1gltexcoordpointer"
#endif
subroutine f9y1gltexcoordpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gltexcoordpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gltexcoordpointer@16' :: f9y1gltexcoordpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y1gltexcoordpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4gltexcoordpointer"
#endif
subroutine f9y4gltexcoordpointer(xsize,xtype,stride,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gltexcoordpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gltexcoordpointer@20' :: f9y4gltexcoordpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y4gltexcoordpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6gltexcoordpointer"
#endif
subroutine f9y6gltexcoordpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gltexcoordpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gltexcoordpointer@16' :: f9y6gltexcoordpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y6gltexcoordpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y7gltexcoordpointer"
#endif
subroutine f9y7gltexcoordpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y7gltexcoordpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y7gltexcoordpointer@16' :: f9y7gltexcoordpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y7gltexcoordpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0gltexsubimage1d"
#endif
subroutine f9y0gltexsubimage1d(xtarget,level,xoffset,width, &
   format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,width
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gltexsubimage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gltexsubimage1d@32' :: f9y0gltexsubimage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0gltexsubimage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1gltexsubimage1d"
#endif
subroutine f9y1gltexsubimage1d(xtarget,level,xoffset,width, &
   format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,width
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gltexsubimage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gltexsubimage1d@28' :: f9y1gltexsubimage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y1gltexsubimage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4gltexsubimage1d"
#endif
subroutine f9y4gltexsubimage1d(xtarget,level,xoffset,width, &
   format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,width
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gltexsubimage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gltexsubimage1d@32' :: f9y4gltexsubimage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4gltexsubimage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5gltexsubimage1d"
#endif
subroutine f9y5gltexsubimage1d(xtarget,level,xoffset,width, &
   format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,width
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5gltexsubimage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5gltexsubimage1d@28' :: f9y5gltexsubimage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y5gltexsubimage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6gltexsubimage1d"
#endif
subroutine f9y6gltexsubimage1d(xtarget,level,xoffset,width, &
   format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,width
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gltexsubimage1d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gltexsubimage1d@32' :: f9y6gltexsubimage1d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6gltexsubimage1d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0gltexsubimage2d"
#endif
subroutine f9y0gltexsubimage2d(xtarget,level,xoffset,yoffset, &
   width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gltexsubimage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gltexsubimage2d@40' :: f9y0gltexsubimage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0gltexsubimage2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1gltexsubimage2d"
#endif
subroutine f9y1gltexsubimage2d(xtarget,level,xoffset,yoffset, &
   width,height,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gltexsubimage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gltexsubimage2d@36' :: f9y1gltexsubimage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y1gltexsubimage2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4gltexsubimage2d"
#endif
subroutine f9y4gltexsubimage2d(xtarget,level,xoffset,yoffset, &
   width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gltexsubimage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gltexsubimage2d@40' :: f9y4gltexsubimage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4gltexsubimage2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5gltexsubimage2d"
#endif
subroutine f9y5gltexsubimage2d(xtarget,level,xoffset,yoffset, &
   width,height,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5gltexsubimage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5gltexsubimage2d@36' :: f9y5gltexsubimage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y5gltexsubimage2d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6gltexsubimage2d"
#endif
subroutine f9y6gltexsubimage2d(xtarget,level,xoffset,yoffset, &
   width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gltexsubimage2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gltexsubimage2d@40' :: f9y6gltexsubimage2d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6gltexsubimage2d
end interface
interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glvertexpointer"
#endif
subroutine f9y0glvertexpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glvertexpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glvertexpointer@16' :: f9y0glvertexpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y0glvertexpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glvertexpointer"
#endif
subroutine f9y1glvertexpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glvertexpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glvertexpointer@16' :: f9y1glvertexpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y1glvertexpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glvertexpointer"
#endif
subroutine f9y4glvertexpointer(xsize,xtype,stride,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glvertexpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glvertexpointer@20' :: f9y4glvertexpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y4glvertexpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glvertexpointer"
#endif
subroutine f9y6glvertexpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glvertexpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glvertexpointer@16' :: f9y6glvertexpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y6glvertexpointer
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y7glvertexpointer"
#endif
subroutine f9y7glvertexpointer(xsize,xtype,stride,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
! !MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y7glvertexpointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y7glvertexpointer@16' :: f9y7glvertexpointer
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN ptr
end subroutine f9y7glvertexpointer
end interface

#endif

! OpenGL 1.2

#ifdef OPENGL_1_2

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0gldrawrangeelements"
#endif
subroutine f9y0gldrawrangeelements(mode,start,xend,xcount,xtype,indices)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,start,xend,xcount,xtype,indices
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gldrawrangeelements
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gldrawrangeelements@24' :: f9y0gldrawrangeelements
#endif
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=gluint) INTENT_IN start, xend
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN indices
end subroutine f9y0gldrawrangeelements
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4gldrawrangeelements"
#endif
subroutine f9y4gldrawrangeelements(mode,start,xend,xcount,xtype,indices,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,start,xend,xcount,xtype,indices,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gldrawrangeelements
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gldrawrangeelements@28' :: f9y4gldrawrangeelements
#endif
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=gluint) INTENT_IN start, xend
integer(kind=glsizei) INTENT_IN xcount
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN indices
integer(kind=glint) INTENT_IN length
end subroutine f9y4gldrawrangeelements
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5gldrawrangeelements"
#endif
subroutine f9y5gldrawrangeelements(mode,start,xend,xcount,xtype,indices)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,start,xend,xcount,xtype,indices
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5gldrawrangeelements
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5gldrawrangeelements@24' :: f9y5gldrawrangeelements
#endif
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=gluint) INTENT_IN start, xend
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN indices
end subroutine f9y5gldrawrangeelements
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6gldrawrangeelements"
#endif
subroutine f9y6gldrawrangeelements(mode,start,xend,xcount,xtype,indices)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,start,xend,xcount,xtype,indices
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gldrawrangeelements
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gldrawrangeelements@24' :: f9y6gldrawrangeelements
#endif
integer(kind=glenum) INTENT_IN mode, xtype
integer(kind=gluint) INTENT_IN start, xend
integer(kind=glsizei) INTENT_IN xcount
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN indices
end subroutine f9y6gldrawrangeelements
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glteximage3d"
#endif
subroutine f9y0glteximage3d(xtarget,level,internalformat, &
   width,height,depth,border,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,width,height
!MS$ATTRIBUTES REFERENCE :: depth,border,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glteximage3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glteximage3d@44' :: f9y0glteximage3d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0glteximage3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glteximage3d"
#endif
subroutine f9y1glteximage3d(xtarget,level,internalformat, &
   width,height,depth,border,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,width
!MS$ATTRIBUTES REFERENCE :: height,depth,border,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glteximage3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glteximage3d@40' :: f9y1glteximage3d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y1glteximage3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glteximage3d"
#endif
subroutine f9y4glteximage3d(xtarget,level,internalformat, &
   width,height,depth,border,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,width,height
!MS$ATTRIBUTES REFERENCE :: depth,border,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glteximage3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glteximage3d@44' :: f9y4glteximage3d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4glteximage3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glteximage3d"
#endif
subroutine f9y5glteximage3d(xtarget,level,internalformat, &
   width,height,depth,border,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,width
!MS$ATTRIBUTES REFERENCE :: height,depth,border,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glteximage3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glteximage3d@40' :: f9y5glteximage3d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y5glteximage3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glteximage3d"
#endif
subroutine f9y6glteximage3d(xtarget,level,internalformat, &
   width,height,depth,border,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,width,height
!MS$ATTRIBUTES REFERENCE :: depth,border,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glteximage3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glteximage3d@44' :: f9y6glteximage3d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6glteximage3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0gltexsubimage3d"
#endif
subroutine f9y0gltexsubimage3d(xtarget,level,xoffset,yoffset, &
   zoffset,width,height,depth,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset,zoffset,width
!MS$ATTRIBUTES REFERENCE :: height,depth,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gltexsubimage3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gltexsubimage3d@48' :: f9y0gltexsubimage3d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0gltexsubimage3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1gltexsubimage3d"
#endif
subroutine f9y1gltexsubimage3d(xtarget,level,xoffset,yoffset, &
   zoffset,width,height,depth,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset,zoffset
!MS$ATTRIBUTES REFERENCE :: width,height,depth,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gltexsubimage3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gltexsubimage3d@44' :: f9y1gltexsubimage3d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y1gltexsubimage3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4gltexsubimage3d"
#endif
subroutine f9y4gltexsubimage3d(xtarget,level,xoffset,yoffset, &
   zoffset,width,height,depth,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset,zoffset,width
!MS$ATTRIBUTES REFERENCE :: height,depth,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gltexsubimage3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gltexsubimage3d@48' :: f9y4gltexsubimage3d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4gltexsubimage3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5gltexsubimage3d"
#endif
subroutine f9y5gltexsubimage3d(xtarget,level,xoffset,yoffset, &
   zoffset,width,height,depth,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset,zoffset
!MS$ATTRIBUTES REFERENCE :: width,height,depth,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5gltexsubimage3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5gltexsubimage3d@44' :: f9y5gltexsubimage3d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y5gltexsubimage3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6gltexsubimage3d"
#endif
subroutine f9y6gltexsubimage3d(xtarget,level,xoffset,yoffset, &
   zoffset,width,height,depth,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset,zoffset,width
!MS$ATTRIBUTES REFERENCE :: height,depth,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gltexsubimage3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gltexsubimage3d@48' :: f9y6gltexsubimage3d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6gltexsubimage3d
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcopytexsubimage3d"
#endif
subroutine fglcopytexsubimage3d(xtarget,level,xoffset, &
   yoffset,zoffset,x,y,width,height)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset
!MS$ATTRIBUTES REFERENCE :: yoffset,zoffset,x,y,width,height
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcopytexsubimage3d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcopytexsubimage3d@36' :: fglcopytexsubimage3d
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset, x, y
integer(kind=glsizei) INTENT_IN width, height
end subroutine fglcopytexsubimage3d
end interface
#endif

! EXT_blend_color

#ifdef GL_EXT_blend_color

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglblendcolorext"
#endif
subroutine fglblendcolorext(red,green,blue,alpha)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: red,green,blue,alpha
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglblendcolorext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglblendcolorext@16' :: fglblendcolorext
#endif
real(kind=glclampf) INTENT_IN red, green, blue, alpha
end subroutine fglblendcolorext
end interface

#endif

! EXT_blend_minmax

#ifdef GL_EXT_blend_minmax

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglblendequationext"
#endif
subroutine fglblendequationext(mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglblendequationext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglblendequationext@4' :: fglblendequationext
#endif
integer(kind=glenum) INTENT_IN mode
end subroutine fglblendequationext
end interface

#endif

! EXT_color_table

#ifdef GL_EXT_color_table

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glcolortableext"
#endif
subroutine f9y0glcolortableext(xtarget,internalformat,width,format,xtype,table, &
   length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,internalformat,width,format,xtype
!MS$ATTRIBUTES REFERENCE :: table,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glcolortableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glcolortableext@28' :: f9y0glcolortableext
#endif
integer(kind=glenum) INTENT_IN xtarget, internalformat
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN table
integer(kind=glint) INTENT_IN length
end subroutine f9y0glcolortableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glcolortableext"
#endif
subroutine f9y1glcolortableext(xtarget,internalformat,width,format,xtype,table)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,internalformat,width,format,xtype,table
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glcolortableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glcolortableext@24' :: f9y1glcolortableext
#endif
integer(kind=glenum) INTENT_IN xtarget, internalformat
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN table
end subroutine f9y1glcolortableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glcolortableext"
#endif
subroutine f9y4glcolortableext(xtarget,internalformat,width,format,xtype,table, &
   length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,internalformat,width,format,xtype
!MS$ATTRIBUTES REFERENCE :: table,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glcolortableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glcolortableext@28' :: f9y4glcolortableext
#endif
integer(kind=glenum) INTENT_IN xtarget, internalformat
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN table
integer(kind=glint) INTENT_IN length
end subroutine f9y4glcolortableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glcolortableext"
#endif
subroutine f9y5glcolortableext(xtarget,internalformat,width,format,xtype,table)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,internalformat,width,format,xtype,table
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glcolortableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glcolortableext@24' :: f9y5glcolortableext
#endif
integer(kind=glenum) INTENT_IN xtarget, internalformat
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN table
end subroutine f9y5glcolortableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glcolortableext"
#endif
subroutine f9y6glcolortableext(xtarget,internalformat,width,format,xtype,table, &
   length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,internalformat,width,format,xtype
!MS$ATTRIBUTES REFERENCE :: table,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glcolortableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glcolortableext@28' :: f9y6glcolortableext
#endif
integer(kind=glenum) INTENT_IN xtarget, internalformat
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN table
integer(kind=glint) INTENT_IN length
end subroutine f9y6glcolortableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glcolorsubtableext"
#endif
subroutine f9y0glcolorsubtableext(xtarget,start,xcount,format,xtype,data, &
   length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,start,xcount,format,xtype,data
!MS$ATTRIBUTES REFERENCE :: length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glcolorsubtableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glcolorsubtableext@28' :: f9y0glcolorsubtableext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glsizei) INTENT_IN start, xcount
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) INTENT_IN length
end subroutine f9y0glcolorsubtableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glcolorsubtableext"
#endif
subroutine f9y1glcolorsubtableext(xtarget,start,xcount,format,xtype,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,start,xcount,format,xtype,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glcolorsubtableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glcolorsubtableext@24' :: f9y1glcolorsubtableext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glsizei) INTENT_IN start, xcount
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN data
end subroutine f9y1glcolorsubtableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glcolorsubtableext"
#endif
subroutine f9y4glcolorsubtableext(xtarget,start,xcount,format,xtype,data, &
   length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,start,xcount,format,xtype,data
!MS$ATTRIBUTES REFERENCE :: length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glcolorsubtableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glcolorsubtableext@28' :: f9y4glcolorsubtableext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glsizei) INTENT_IN start, xcount
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) INTENT_IN length
end subroutine f9y4glcolorsubtableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glcolorsubtableext"
#endif
subroutine f9y5glcolorsubtableext(xtarget,start,xcount,format,xtype,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,start,xcount,format,xtype,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glcolorsubtableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glcolorsubtableext@24' :: f9y5glcolorsubtableext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glsizei) INTENT_IN start, xcount
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN data
end subroutine f9y5glcolorsubtableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glcolorsubtableext"
#endif
subroutine f9y6glcolorsubtableext(xtarget,start,xcount,format,xtype,data, &
   length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,start,xcount,format,xtype,data
!MS$ATTRIBUTES REFERENCE :: length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glcolorsubtableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glcolorsubtableext@28' :: f9y6glcolorsubtableext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glsizei) INTENT_IN start, xcount
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) INTENT_IN length
end subroutine f9y6glcolorsubtableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glgetcolortableext"
#endif
subroutine f9y0glgetcolortableext(xtarget,format,xtype,table,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,format,xtype,table,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glgetcolortableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glgetcolortableext@20' :: f9y0glgetcolortableext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: table
integer(kind=glint) INTENT_IN length
end subroutine f9y0glgetcolortableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glgetcolortableext"
#endif
subroutine f9y1glgetcolortableext(xtarget,format,xtype,table,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,format,xtype,table,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glgetcolortableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glgetcolortableext@16' :: f9y1glgetcolortableext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: table
integer(kind=glint) INTENT_IN length
end subroutine f9y1glgetcolortableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glgetcolortableext"
#endif
subroutine f9y4glgetcolortableext(xtarget,format,xtype,table,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,format,xtype,table,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glgetcolortableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glgetcolortableext@16' :: f9y4glgetcolortableext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: table
integer(kind=glint) INTENT_IN length
end subroutine f9y4glgetcolortableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glgetcolortableext"
#endif
subroutine f9y5glgetcolortableext(xtarget,format,xtype,table,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,format,xtype,table,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glgetcolortableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glgetcolortableext@16' :: f9y5glgetcolortableext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: table
integer(kind=glint) INTENT_IN length
end subroutine f9y5glgetcolortableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glgetcolortableext"
#endif
subroutine f9y6glgetcolortableext(xtarget,format,xtype,table,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,format,xtype,table,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glgetcolortableext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glgetcolortableext@20' :: f9y6glgetcolortableext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM), intent(out) :: table
integer(kind=glint) INTENT_IN length
end subroutine f9y6glgetcolortableext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetcolortableparameterfvext"
#endif
subroutine fglgetcolortableparameterfvext(xtarget,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetcolortableparameterfvext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetcolortableparameterfvext@12' :: fglgetcolortableparameterfvext
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgetcolortableparameterfvext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetcolortableparameterivext"
#endif
subroutine fglgetcolortableparameterivext(xtarget,pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetcolortableparameterivext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetcolortableparameterivext@12' :: fglgetcolortableparameterivext
#endif
integer(kind=glenum) INTENT_IN xtarget, pname
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: params
end subroutine fglgetcolortableparameterivext
end interface

#endif

! EXT_point_parameters

#ifdef GL_EXT_point_parameters

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpointparameterfext"
#endif
subroutine fglpointparameterfext(pname,param)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpointparameterfext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpointparameterfext@8' :: fglpointparameterfext
#endif
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat) INTENT_IN param
end subroutine fglpointparameterfext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpointparameterfvext"
#endif
subroutine fglpointparameterfvext(pname,params)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpointparameterfvext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpointparameterfvext@8' :: fglpointparameterfvext
#endif
integer(kind=glenum) INTENT_IN pname
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN params
end subroutine fglpointparameterfvext
end interface

#endif

! EXT_polygon_offset

#ifdef GL_EXT_polygon_offset

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglpolygonoffsetext"
#endif
subroutine fglpolygonoffsetext(factor,bias)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: factor,bias
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglpolygonoffsetext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglpolygonoffsetext@8' :: fglpolygonoffsetext
#endif
real(kind=glfloat) INTENT_IN factor, bias
end subroutine fglpolygonoffsetext
end interface

#endif

! EXT_subtexture

#ifdef GL_EXT_subtexture

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0gltexsubimage1dext"
#endif
subroutine f9y0gltexsubimage1dext(xtarget,level,xoffset,width, &
   format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,width
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gltexsubimage1dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gltexsubimage1dext@32' :: f9y0gltexsubimage1dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0gltexsubimage1dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1gltexsubimage1dext"
#endif
subroutine f9y1gltexsubimage1dext(xtarget,level,xoffset,width, &
   format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,width
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gltexsubimage1dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gltexsubimage1dext@28' :: f9y1gltexsubimage1dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y1gltexsubimage1dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4gltexsubimage1dext"
#endif
subroutine f9y4gltexsubimage1dext(xtarget,level,xoffset,width, &
   format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,width
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gltexsubimage1dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gltexsubimage1dext@32' :: f9y4gltexsubimage1dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4gltexsubimage1dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5gltexsubimage1dext"
#endif
subroutine f9y5gltexsubimage1dext(xtarget,level,xoffset,width, &
   format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,width
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5gltexsubimage1dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5gltexsubimage1dext@28' :: f9y5gltexsubimage1dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y5gltexsubimage1dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6gltexsubimage1dext"
#endif
subroutine f9y6gltexsubimage1dext(xtarget,level,xoffset,width, &
   format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,width
!MS$ATTRIBUTES REFERENCE :: format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gltexsubimage1dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gltexsubimage1dext@32' :: f9y6gltexsubimage1dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6gltexsubimage1dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0gltexsubimage2dext"
#endif
subroutine f9y0gltexsubimage2dext(xtarget,level,xoffset,yoffset, &
   width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gltexsubimage2dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gltexsubimage2dext@40' :: f9y0gltexsubimage2dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0gltexsubimage2dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1gltexsubimage2dext"
#endif
subroutine f9y1gltexsubimage2dext(xtarget,level,xoffset,yoffset, &
   width,height,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gltexsubimage2dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gltexsubimage2dext@36' :: f9y1gltexsubimage2dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y1gltexsubimage2dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4gltexsubimage2dext"
#endif
subroutine f9y4gltexsubimage2dext(xtarget,level,xoffset,yoffset, &
   width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gltexsubimage2dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gltexsubimage2dext@40' :: f9y4gltexsubimage2dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4gltexsubimage2dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5gltexsubimage2dext"
#endif
subroutine f9y5gltexsubimage2dext(xtarget,level,xoffset,yoffset, &
   width,height,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5gltexsubimage2dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5gltexsubimage2dext@36' :: f9y5gltexsubimage2dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y5gltexsubimage2dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6gltexsubimage2dext"
#endif
subroutine f9y6gltexsubimage2dext(xtarget,level,xoffset,yoffset, &
   width,height,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset
!MS$ATTRIBUTES REFERENCE :: width,height,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gltexsubimage2dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gltexsubimage2dext@40' :: f9y6gltexsubimage2dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6gltexsubimage2dext
end interface

#endif

! EXT_texture3D

#ifdef GL_EXT_texture3D

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glteximage3dext"
#endif
subroutine f9y0glteximage3dext(xtarget,level,internalformat, &
   width,height,depth,border,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,width,height
!MS$ATTRIBUTES REFERENCE :: depth,border,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glteximage3dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glteximage3dext@44' :: f9y0glteximage3dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0glteximage3dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glteximage3dext"
#endif
subroutine f9y1glteximage3dext(xtarget,level,internalformat, &
   width,height,depth,border,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,width
!MS$ATTRIBUTES REFERENCE :: height,depth,border,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glteximage3dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glteximage3dext@40' :: f9y1glteximage3dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y1glteximage3dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glteximage3dext"
#endif
subroutine f9y4glteximage3dext(xtarget,level,internalformat, &
   width,height,depth,border,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,width,height
!MS$ATTRIBUTES REFERENCE :: depth,border,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glteximage3dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glteximage3dext@44' :: f9y4glteximage3dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4glteximage3dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glteximage3dext"
#endif
subroutine f9y5glteximage3dext(xtarget,level,internalformat, &
   width,height,depth,border,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,width
!MS$ATTRIBUTES REFERENCE :: height,depth,border,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glteximage3dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glteximage3dext@40' :: f9y5glteximage3dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y5glteximage3dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glteximage3dext"
#endif
subroutine f9y6glteximage3dext(xtarget,level,internalformat, &
   width,height,depth,border,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,internalformat,width,height
!MS$ATTRIBUTES REFERENCE :: depth,border,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glteximage3dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glteximage3dext@44' :: f9y6glteximage3dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level
integer(kind=glenum) INTENT_IN internalformat
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glint) INTENT_IN border
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6glteximage3dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0gltexsubimage3dext"
#endif
subroutine f9y0gltexsubimage3dext(xtarget,level,xoffset,yoffset, &
   zoffset,width,height,depth,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset,zoffset,width
!MS$ATTRIBUTES REFERENCE :: height,depth,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gltexsubimage3dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gltexsubimage3dext@48' :: f9y0gltexsubimage3dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y0gltexsubimage3dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1gltexsubimage3dext"
#endif
subroutine f9y1gltexsubimage3dext(xtarget,level,xoffset,yoffset, &
   zoffset,width,height,depth,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset,zoffset
!MS$ATTRIBUTES REFERENCE :: width,height,depth,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gltexsubimage3dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gltexsubimage3dext@44' :: f9y1gltexsubimage3dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y1gltexsubimage3dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4gltexsubimage3dext"
#endif
subroutine f9y4gltexsubimage3dext(xtarget,level,xoffset,yoffset, &
   zoffset,width,height,depth,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset,zoffset,width
!MS$ATTRIBUTES REFERENCE :: height,depth,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gltexsubimage3dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gltexsubimage3dext@48' :: f9y4gltexsubimage3dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y4gltexsubimage3dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5gltexsubimage3dext"
#endif
subroutine f9y5gltexsubimage3dext(xtarget,level,xoffset,yoffset, &
   zoffset,width,height,depth,format,xtype,pixels)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset,zoffset
!MS$ATTRIBUTES REFERENCE :: width,height,depth,format,xtype,pixels
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5gltexsubimage3dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5gltexsubimage3dext@44' :: f9y5gltexsubimage3dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN pixels
end subroutine f9y5gltexsubimage3dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6gltexsubimage3dext"
#endif
subroutine f9y6gltexsubimage3dext(xtarget,level,xoffset,yoffset, &
   zoffset,width,height,depth,format,xtype,pixels,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset,yoffset,zoffset,width
!MS$ATTRIBUTES REFERENCE :: height,depth,format,xtype,pixels,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gltexsubimage3dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gltexsubimage3dext@48' :: f9y6gltexsubimage3dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset
integer(kind=glsizei) INTENT_IN width, height, depth
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN pixels
integer(kind=glint) INTENT_IN length
end subroutine f9y6gltexsubimage3dext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglcopytexsubimage3dext"
#endif
subroutine fglcopytexsubimage3dext(xtarget,level,xoffset, &
   yoffset,zoffset,x,y,width,height)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,level,xoffset
!MS$ATTRIBUTES REFERENCE :: yoffset,zoffset,x,y,width,height
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglcopytexsubimage3dext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglcopytexsubimage3dext@36' :: fglcopytexsubimage3dext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN level, xoffset, yoffset, zoffset, x, y
integer(kind=glsizei) INTENT_IN width, height
end subroutine fglcopytexsubimage3dext
end interface

#endif

! EXT_texture_object

#ifdef GL_EXT_texture_object

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgentexturesext"
#endif
subroutine fglgentexturesext(n,textures)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,textures
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgentexturesext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgentexturesext@8' :: fglgentexturesext
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: textures
end subroutine fglgentexturesext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgldeletetexturesext"
#endif
subroutine fgldeletetexturesext(n,textures)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,textures
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgldeletetexturesext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgldeletetexturesext@8' :: fgldeletetexturesext
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
end subroutine fgldeletetexturesext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglbindtextureext"
#endif
subroutine fglbindtextureext(xtarget,texture)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,texture
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglbindtextureext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglbindtextureext@8' :: fglbindtextureext
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=gluint) INTENT_IN texture
end subroutine fglbindtextureext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglprioritizetexturesext"
#endif
subroutine fglprioritizetexturesext(n,textures,priorities)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,textures,priorities
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglprioritizetexturesext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglprioritizetexturesext@12' :: fglprioritizetexturesext
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
real(kind=glclampf), dimension(DEFER_DIM) INTENT_IN priorities
end subroutine fglprioritizetexturesext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglaretexturesresidentext"
#endif
#ifdef PURE_FUNCTIONS
subroutine fglaretexturesresidentext(n,textures, &
   residences,areresident)
#else
function fglaretexturesresidentext(n,textures, &
   residences) result(resfglaretexturesresidentext)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,textures
!MS$ATTRIBUTES REFERENCE :: residences
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglaretexturesresidentext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglaretexturesresidentext@12' :: fglaretexturesresidentext
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
logical(kind=glboolean), dimension(DEFER_DIM), intent(out) :: residences
#ifdef PURE_FUNCTIONS
logical(kind=glboolean), intent(out) :: areresident
#else
logical(kind=glboolean) :: resfglaretexturesresidentext
#endif
#ifdef PURE_FUNCTIONS
end subroutine fglaretexturesresidentext
#else
end function fglaretexturesresidentext
#endif
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9eglaretexturesresidentext"
#endif
#ifdef PURE_FUNCTIONS
subroutine f9eglaretexturesresidentext(n,textures, &
   residences,areresident)
#else
function f9eglaretexturesresidentext(n,textures, &
   residences) result(resf9eglaretexturesresidentext)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: n,textures
!MS$ATTRIBUTES REFERENCE :: residences
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9eglaretexturesresidentext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9eglaretexturesresidentext@12' :: f9eglaretexturesresidentext
#endif
integer(kind=glsizei) INTENT_IN n
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN textures
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: residences
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: areresident
#else
integer(kind=glint) :: resf9eglaretexturesresidentext
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9eglaretexturesresidentext
#else
end function f9eglaretexturesresidentext
#endif
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglistextureext"
#endif
function fglistextureext(texture) result(resfglistextureext)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: texture
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglistextureext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglistextureext@4' :: fglistextureext
#endif
integer(kind=gluint) INTENT_IN texture
logical(kind=glboolean) :: resfglistextureext
end function fglistextureext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9eglistextureext"
#endif
function f9eglistextureext(texture) result(resf9eglistextureext)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: texture
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9eglistextureext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9eglistextureext@4' :: f9eglistextureext
#endif
integer(kind=gluint) INTENT_IN texture
integer(kind=glint) :: resf9eglistextureext
end function f9eglistextureext
end interface

#endif

! EXT_vertex_array

#ifdef GL_EXT_vertex_array

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glvertexpointerext"
#endif
subroutine f9y0glvertexpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glvertexpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glvertexpointerext@24' :: f9y0glvertexpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glint), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y0glvertexpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glvertexpointerext"
#endif
subroutine f9y1glvertexpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glvertexpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glvertexpointerext@20' :: f9y1glvertexpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y1glvertexpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glvertexpointerext"
#endif
subroutine f9y4glvertexpointerext(xsize,xtype,stride,xcount,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glvertexpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glvertexpointerext@24' :: f9y4glvertexpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y4glvertexpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glvertexpointerext"
#endif
subroutine f9y6glvertexpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glvertexpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glvertexpointerext@20' :: f9y6glvertexpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glshort), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y6glvertexpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y7glvertexpointerext"
#endif
subroutine f9y7glvertexpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y7glvertexpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y7glvertexpointerext@20' :: f9y7glvertexpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=gldouble), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y7glvertexpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glnormalpointerext"
#endif
subroutine f9y0glnormalpointerext(xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glnormalpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glnormalpointerext@20' :: f9y0glnormalpointerext
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glint), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y0glnormalpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glnormalpointerext"
#endif
subroutine f9y1glnormalpointerext(xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glnormalpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glnormalpointerext@16' :: f9y1glnormalpointerext
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y1glnormalpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glnormalpointerext"
#endif
subroutine f9y4glnormalpointerext(xtype,stride,xcount,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,xcount,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glnormalpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glnormalpointerext@20' :: f9y4glnormalpointerext
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y4glnormalpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glnormalpointerext"
#endif
subroutine f9y5glnormalpointerext(xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glnormalpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glnormalpointerext@16' :: f9y5glnormalpointerext
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glbyte), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y5glnormalpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glnormalpointerext"
#endif
subroutine f9y6glnormalpointerext(xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glnormalpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glnormalpointerext@20' :: f9y6glnormalpointerext
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glshort), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y6glnormalpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y7glnormalpointerext"
#endif
subroutine f9y7glnormalpointerext(xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y7glnormalpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y7glnormalpointerext@16' :: f9y7glnormalpointerext
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=gldouble), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y7glnormalpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glcolorpointerext"
#endif
subroutine f9y0glcolorpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glcolorpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glcolorpointerext@24' :: f9y0glcolorpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glint), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y0glcolorpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glcolorpointerext"
#endif
subroutine f9y1glcolorpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glcolorpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glcolorpointerext@20' :: f9y1glcolorpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y1glcolorpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glcolorpointerext"
#endif
subroutine f9y4glcolorpointerext(xsize,xtype,stride,xcount,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glcolorpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glcolorpointerext@24' :: f9y4glcolorpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y4glcolorpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y5glcolorpointerext"
#endif
subroutine f9y5glcolorpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glcolorpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glcolorpointerext@20' :: f9y5glcolorpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glbyte), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y5glcolorpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glcolorpointerext"
#endif
subroutine f9y6glcolorpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glcolorpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glcolorpointerext@24' :: f9y6glcolorpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glshort), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y6glcolorpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y7glcolorpointerext"
#endif
subroutine f9y7glcolorpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y7glcolorpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y7glcolorpointerext@20' :: f9y7glcolorpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=gldouble), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y7glcolorpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0glindexpointerext"
#endif
subroutine f9y0glindexpointerext(xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glindexpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glindexpointerext@16' :: f9y0glindexpointerext
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glint), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y0glindexpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1glindexpointerext"
#endif
subroutine f9y1glindexpointerext(xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glindexpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glindexpointerext@16' :: f9y1glindexpointerext
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y1glindexpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4glindexpointerext"
#endif
subroutine f9y4glindexpointerext(xtype,stride,xcount,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,xcount,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glindexpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glindexpointerext@20' :: f9y4glindexpointerext
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y4glindexpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6glindexpointerext"
#endif
subroutine f9y6glindexpointerext(xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glindexpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glindexpointerext@16' :: f9y6glindexpointerext
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glshort), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y6glindexpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y7glindexpointerext"
#endif
subroutine f9y7glindexpointerext(xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y7glindexpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y7glindexpointerext@16' :: f9y7glindexpointerext
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=gldouble), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y7glindexpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y0gltexcoordpointerext"
#endif
subroutine f9y0gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gltexcoordpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gltexcoordpointerext@24' :: f9y0gltexcoordpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glint), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y0gltexcoordpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y1gltexcoordpointerext"
#endif
subroutine f9y1gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gltexcoordpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gltexcoordpointerext@20' :: f9y1gltexcoordpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y1gltexcoordpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y4gltexcoordpointerext"
#endif
subroutine f9y4gltexcoordpointerext(xsize,xtype,stride,xcount,ptr,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gltexcoordpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gltexcoordpointerext@24' :: f9y4gltexcoordpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN ptr
integer(kind=glint) INTENT_IN length
end subroutine f9y4gltexcoordpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y6gltexcoordpointerext"
#endif
subroutine f9y6gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gltexcoordpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gltexcoordpointerext@20' :: f9y6gltexcoordpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glshort), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y6gltexcoordpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9y7gltexcoordpointerext"
#endif
subroutine f9y7gltexcoordpointerext(xsize,xtype,stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize,xtype,stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y7gltexcoordpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y7gltexcoordpointerext@20' :: f9y7gltexcoordpointerext
#endif
integer(kind=glint) INTENT_IN xsize
integer(kind=glenum) INTENT_IN xtype
integer(kind=glsizei) INTENT_IN stride, xcount
real(kind=gldouble), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9y7gltexcoordpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgledgeflagpointerext"
#endif
subroutine fgledgeflagpointerext(stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgledgeflagpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgledgeflagpointerext@12' :: fgledgeflagpointerext
#endif
integer(kind=glsizei) INTENT_IN stride, xcount
logical(kind=glboolean), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine fgledgeflagpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9egledgeflagpointerext"
#endif
subroutine f9egledgeflagpointerext(stride,xcount,ptr)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: stride,xcount,ptr
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9egledgeflagpointerext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9egledgeflagpointerext@12' :: f9egledgeflagpointerext
#endif
integer(kind=glsizei) INTENT_IN stride, xcount
integer(kind=glbyte), dimension(DEFER_DIM), TARGET intent(in out) :: ptr
end subroutine f9egledgeflagpointerext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglgetpointervext"
#endif
subroutine fglgetpointervext(pname,params,length)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pname,params,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglgetpointervext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglgetpointervext@12' :: fglgetpointervext
#endif
integer(kind=glenum) INTENT_IN pname
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: params
integer(kind=glint) INTENT_IN length
end subroutine fglgetpointervext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglarrayelementext"
#endif
subroutine fglarrayelementext(i)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: i
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglarrayelementext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglarrayelementext@4' :: fglarrayelementext
#endif
integer(kind=glint) INTENT_IN i
end subroutine fglarrayelementext
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fgldrawarraysext"
#endif
subroutine fgldrawarraysext(mode,first,xcount)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode,first,xcount
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgldrawarraysext
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgldrawarraysext@12' :: fgldrawarraysext
#endif
integer(kind=glenum) INTENT_IN mode
integer(kind=glint) INTENT_IN first
integer(kind=glsizei) INTENT_IN xcount
end subroutine fgldrawarraysext
end interface

#endif

! SGIS_multisample

#ifdef GL_SGIS_multisample

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglsamplemasksgis"
#endif
subroutine fglsamplemasksgis(value,invert)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: value,invert
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglsamplemasksgis
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglsamplemasksgis@8' :: fglsamplemasksgis
#endif
real(kind=glclampf) INTENT_IN value
logical(kind=glboolean) INTENT_IN invert
end subroutine fglsamplemasksgis
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9eglsamplemasksgis"
#endif
subroutine f9eglsamplemasksgis(value,invert)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: value,invert
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9eglsamplemasksgis
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9eglsamplemasksgis@8' :: f9eglsamplemasksgis
#endif
real(kind=glclampf) INTENT_IN value
integer(kind=glint) INTENT_IN invert
end subroutine f9eglsamplemasksgis
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglsamplepatternsgis"
#endif
subroutine fglsamplepatternsgis(pattern)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pattern
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglsamplepatternsgis
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglsamplepatternsgis@4' :: fglsamplepatternsgis
#endif
integer(kind=glenum) INTENT_IN pattern
end subroutine fglsamplepatternsgis
end interface

#endif

! MESA_resize_buffers

#ifdef GL_MESA_resize_buffers

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglresizebuffersmesa"
#endif
subroutine fglresizebuffersmesa()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglresizebuffersmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglresizebuffersmesa@0' :: fglresizebuffersmesa
#endif
end subroutine fglresizebuffersmesa
end interface

#endif

! MESA_window_pos

#ifdef GL_MESA_window_pos

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos2imesa"
#endif
subroutine fglwindowpos2imesa(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos2imesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos2imesa@8' :: fglwindowpos2imesa
#endif
integer(kind=glint) INTENT_IN x, y
end subroutine fglwindowpos2imesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos2smesa"
#endif
subroutine fglwindowpos2smesa(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos2smesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos2smesa@8' :: fglwindowpos2smesa
#endif
integer(kind=glshort) INTENT_IN x, y
end subroutine fglwindowpos2smesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglwindowpos2smesa"
#endif
subroutine f9aglwindowpos2smesa(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglwindowpos2smesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglwindowpos2smesa@8' :: f9aglwindowpos2smesa
#endif
integer(kind=glshort) INTENT_IN x, y
end subroutine f9aglwindowpos2smesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos2fmesa"
#endif
subroutine fglwindowpos2fmesa(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos2fmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos2fmesa@8' :: fglwindowpos2fmesa
#endif
real(kind=glfloat) INTENT_IN x, y
end subroutine fglwindowpos2fmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos2dmesa"
#endif
subroutine fglwindowpos2dmesa(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos2dmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos2dmesa@8' :: fglwindowpos2dmesa
#endif
real(kind=gldouble) INTENT_IN x, y
end subroutine fglwindowpos2dmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos2ivmesa"
#endif
subroutine fglwindowpos2ivmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos2ivmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos2ivmesa@4' :: fglwindowpos2ivmesa
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos2ivmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos2svmesa"
#endif
subroutine fglwindowpos2svmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos2svmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos2svmesa@4' :: fglwindowpos2svmesa
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos2svmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglwindowpos2svmesa"
#endif
subroutine f9aglwindowpos2svmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglwindowpos2svmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglwindowpos2svmesa@4' :: f9aglwindowpos2svmesa
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN p
end subroutine f9aglwindowpos2svmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos2fvmesa"
#endif
subroutine fglwindowpos2fvmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos2fvmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos2fvmesa@4' :: fglwindowpos2fvmesa
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos2fvmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos2dvmesa"
#endif
subroutine fglwindowpos2dvmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos2dvmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos2dvmesa@4' :: fglwindowpos2dvmesa
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos2dvmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos3imesa"
#endif
subroutine fglwindowpos3imesa(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos3imesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos3imesa@12' :: fglwindowpos3imesa
#endif
integer(kind=glint) INTENT_IN x, y, z
end subroutine fglwindowpos3imesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos3smesa"
#endif
subroutine fglwindowpos3smesa(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos3smesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos3smesa@12' :: fglwindowpos3smesa
#endif
integer(kind=glshort) INTENT_IN x, y, z
end subroutine fglwindowpos3smesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglwindowpos3smesa"
#endif
subroutine f9aglwindowpos3smesa(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglwindowpos3smesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglwindowpos3smesa@12' :: f9aglwindowpos3smesa
#endif
integer(kind=glshort) INTENT_IN x, y, z
end subroutine f9aglwindowpos3smesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos3fmesa"
#endif
subroutine fglwindowpos3fmesa(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos3fmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos3fmesa@12' :: fglwindowpos3fmesa
#endif
real(kind=glfloat) INTENT_IN x, y, z
end subroutine fglwindowpos3fmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos3dmesa"
#endif
subroutine fglwindowpos3dmesa(x,y,z)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos3dmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos3dmesa@12' :: fglwindowpos3dmesa
#endif
real(kind=gldouble) INTENT_IN x, y, z
end subroutine fglwindowpos3dmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos3ivmesa"
#endif
subroutine fglwindowpos3ivmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos3ivmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos3ivmesa@4' :: fglwindowpos3ivmesa
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos3ivmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos3svmesa"
#endif
subroutine fglwindowpos3svmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos3svmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos3svmesa@4' :: fglwindowpos3svmesa
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos3svmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglwindowpos3svmesa"
#endif
subroutine f9aglwindowpos3svmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglwindowpos3svmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglwindowpos3svmesa@4' :: f9aglwindowpos3svmesa
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN p
end subroutine f9aglwindowpos3svmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos3fvmesa"
#endif
subroutine fglwindowpos3fvmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos3fvmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos3fvmesa@4' :: fglwindowpos3fvmesa
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos3fvmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos3dvmesa"
#endif
subroutine fglwindowpos3dvmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos3dvmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos3dvmesa@4' :: fglwindowpos3dvmesa
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos3dvmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos4imesa"
#endif
subroutine fglwindowpos4imesa(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos4imesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos4imesa@16' :: fglwindowpos4imesa
#endif
integer(kind=glint) INTENT_IN x, y, z, w
end subroutine fglwindowpos4imesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos4smesa"
#endif
subroutine fglwindowpos4smesa(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos4smesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos4smesa@16' :: fglwindowpos4smesa
#endif
integer(kind=glshort) INTENT_IN x, y, z, w
end subroutine fglwindowpos4smesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglwindowpos4smesa"
#endif
subroutine f9aglwindowpos4smesa(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglwindowpos4smesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglwindowpos4smesa@16' :: f9aglwindowpos4smesa
#endif
integer(kind=glshort) INTENT_IN x, y, z, w
end subroutine f9aglwindowpos4smesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos4fmesa"
#endif
subroutine fglwindowpos4fmesa(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos4fmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos4fmesa@16' :: fglwindowpos4fmesa
#endif
real(kind=glfloat) INTENT_IN x, y, z, w
end subroutine fglwindowpos4fmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos4dmesa"
#endif
subroutine fglwindowpos4dmesa(x,y,z,w)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,z,w
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos4dmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos4dmesa@16' :: fglwindowpos4dmesa
#endif
real(kind=gldouble) INTENT_IN x, y, z, w
end subroutine fglwindowpos4dmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos4ivmesa"
#endif
subroutine fglwindowpos4ivmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos4ivmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos4ivmesa@4' :: fglwindowpos4ivmesa
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos4ivmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos4svmesa"
#endif
subroutine fglwindowpos4svmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos4svmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos4svmesa@4' :: fglwindowpos4svmesa
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos4svmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 f9aglwindowpos4svmesa"
#endif
subroutine f9aglwindowpos4svmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9aglwindowpos4svmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9aglwindowpos4svmesa@4' :: f9aglwindowpos4svmesa
#endif
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN p
end subroutine f9aglwindowpos4svmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos4fvmesa"
#endif
subroutine fglwindowpos4fvmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos4fvmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos4fvmesa@4' :: fglwindowpos4fvmesa
#endif
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos4fvmesa
end interface

interface
#ifdef NASOFTWARE
!NAS$ ALIEN "F77 fglwindowpos4dvmesa"
#endif
subroutine fglwindowpos4dvmesa(p)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: p
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglwindowpos4dvmesa
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglwindowpos4dvmesa@4' :: fglwindowpos4dvmesa
#endif
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN p
end subroutine fglwindowpos4dvmesa
end interface

#endif

end module opengl_glinterfaces
