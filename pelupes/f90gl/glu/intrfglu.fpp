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

module opengl_gluinterfaces
use opengl_kinds
#ifdef SALFORD
 STDCALL fglubegincurve 'fglubegincurve' (REF,REF)
 STDCALL fglubeginpolygon 'fglubeginpolygon' (REF,REF)
 STDCALL fglubeginsurface 'fglubeginsurface' (REF,REF)
 STDCALL fglubegintrim 'fglubegintrim' (REF,REF)
 STDCALL f9y0glubuild1dmipmaps 'f9y0glubuild1dmipmaps' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glubuild1dmipmaps 'f9y1glubuild1dmipmaps' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y4glubuild1dmipmaps 'f9y4glubuild1dmipmaps' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5glubuild1dmipmaps 'f9y5glubuild1dmipmaps' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y6glubuild1dmipmaps 'f9y6glubuild1dmipmaps' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y0glubuild2dmipmaps 'f9y0glubuild2dmipmaps' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1glubuild2dmipmaps 'f9y1glubuild2dmipmaps' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4glubuild2dmipmaps 'f9y4glubuild2dmipmaps' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5glubuild2dmipmaps 'f9y5glubuild2dmipmaps' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6glubuild2dmipmaps 'f9y6glubuild2dmipmaps' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglucylinder 'fglucylinder' (REF,REF,REF,REF,REF,REF,REF)
 STDCALL fgludeletenurbsrenderer 'fgludeletenurbsrenderer' (REF,REF)
 STDCALL fgludeletequadric 'fgludeletequadric' (REF,REF)
 STDCALL fgludeletetess 'fgludeletetess' (REF,REF)
 STDCALL fgludisk 'fgludisk' (REF,REF,REF,REF,REF,REF)
 STDCALL fgluendcurve 'fgluendcurve' (REF,REF)
 STDCALL fgluendpolygon 'fgluendpolygon' (REF,REF)
 STDCALL fgluendsurface 'fgluendsurface' (REF,REF)
 STDCALL fgluendtrim 'fgluendtrim' (REF,REF)
 STDCALL f9ygluerrorstringlen 'f9ygluerrorstringlen' (REF)
 STDCALL fglugetnurbsproperty 'fglugetnurbsproperty' (REF,REF,REF,REF)
 STDCALL f9yglugetstringlen 'f9yglugetstringlen' (REF)
 STDCALL f9yglugetstring 'f9yglugetstring' (REF)
 STDCALL f9yglugetstring1 'f9yglugetstring1' (REF)
 STDCALL f9yglugetstringfree 'f9yglugetstringfree'
 STDCALL fglugettessproperty 'fglugettessproperty' (REF,REF,REF,REF)
 STDCALL fgluloadsamplingmatrices 'fgluloadsamplingmatrices' (REF,REF,REF,REF,REF)
 STDCALL fglulookat 'fglulookat' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9yglunewnurbsrenderer 'f9yglunewnurbsrenderer' (REF,REF)
 STDCALL f9yglunewquadric 'f9yglunewquadric' (REF,REF)
 STDCALL f9yglunewtess 'f9yglunewtess' (REF,REF)
 STDCALL fglunextcontour 'fglunextcontour' (REF,REF,REF)
 STDCALL f9y1glunurbscallback 'f9y1glunurbscallback' (REF,REF,REF)
 STDCALL fglunurbscurve 'fglunurbscurve' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglunurbsproperty 'fglunurbsproperty' (REF,REF,REF,REF)
 STDCALL fglunurbssurface 'fglunurbssurface' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fgluortho2d 'fgluortho2d' (REF,REF,REF,REF)
 STDCALL fglupartialdisk 'fglupartialdisk' (REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fgluperspective 'fgluperspective' (REF,REF,REF,REF)
 STDCALL fglupickmatrix 'fglupickmatrix' (REF,REF,REF,REF,REF)
 STDCALL fgluproject 'fgluproject' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglupwlcurve 'fglupwlcurve' (REF,REF,REF,REF,REF,REF)
 STDCALL f9y1gluquadriccallback 'f9y1gluquadriccallback' (REF,REF,REF)
 STDCALL fgluquadricdrawstyle 'fgluquadricdrawstyle' (REF,REF,REF)
 STDCALL fgluquadricnormals 'fgluquadricnormals' (REF,REF,REF)
 STDCALL fgluquadricorientation 'fgluquadricorientation' (REF,REF,REF)
 STDCALL fgluquadrictexture 'fgluquadrictexture' (REF,REF,REF)
 STDCALL f9egluquadrictexture 'f9egluquadrictexture' (REF,REF,REF)
 STDCALL f9y0gluscaleimage 'f9y0gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y1gluscaleimage 'f9y1gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y4gluscaleimage 'f9y4gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y5gluscaleimage 'f9y5gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y6gluscaleimage 'f9y6gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y01gluscaleimage 'f9y01gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y04gluscaleimage 'f9y04gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y05gluscaleimage 'f9y05gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y06gluscaleimage 'f9y06gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y00gluscaleimage 'f9y10gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y04gluscaleimage 'f9y14gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y05gluscaleimage 'f9y15gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y06gluscaleimage 'f9y16gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y00gluscaleimage 'f9y40gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y01gluscaleimage 'f9y41gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y05gluscaleimage 'f9y45gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y06gluscaleimage 'f9y46gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y00gluscaleimage 'f9y50gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y01gluscaleimage 'f9y51gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y04gluscaleimage 'f9y54gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y06gluscaleimage 'f9y56gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y00gluscaleimage 'f9y60gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y01gluscaleimage 'f9y61gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y04gluscaleimage 'f9y64gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL f9y05gluscaleimage 'f9y65gluscaleimage' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglusphere 'fglusphere' (REF,REF,REF,REF,REF)
 STDCALL fglutessbegincontour 'fglutessbegincontour' (REF,REF)
 STDCALL f9y0glutessbeginpolygon 'f9y0glutessbeginpolygon' (REF,REF,REF)
 STDCALL f9y1glutessbeginpolygon 'f9y1glutessbeginpolygon' (REF,REF,REF)
 STDCALL f9y1glutesscallback 'f9y1glutesscallback' (REF,REF,REF)
 STDCALL fglutessendcontour 'fglutessendcontour' (REF,REF)
 STDCALL fglutessendpolygon 'fglutessendpolygon' (REF,REF)
 STDCALL fglutessnormal 'fglutessnormal' (REF,REF,REF,REF,REF)
 STDCALL fglutessproperty 'fglutessproperty' (REF,REF,REF,REF)
 STDCALL f9y0glutessvertex 'f9y0glutessvertex' (REF,REF,REF,REF)
 STDCALL f9y1glutessvertex 'f9y1glutessvertex' (REF,REF,REF,REF)
 STDCALL f9y4glutessvertex 'f9y4glutessvertex' (REF,REF,REF)
 STDCALL f9y5glutessvertex 'f9y5glutessvertex' (REF,REF,REF,REF)
 STDCALL f9y6glutessvertex 'f9y6glutessvertex' (REF,REF,REF,REF)
 STDCALL f9y7glutessvertex 'f9y7glutessvertex' (REF,REF,REF,REF)
 STDCALL fgluunproject 'fgluunproject' (REF,REF,REF,REF,REF,REF,REF,REF,REF)
 STDCALL fglunurbscallback 'fglunurbscallback' (REF,REF,REF,REF)
 STDCALL fgluquadriccallback 'fgluquadriccallback' (REF,REF,REF,REF)
 STDCALL f9y0glutesscallback 'f9y0glutesscallback' (REF,REF,REF,REF)
#endif
private

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
! workaround for FPlus beta test bug report 15
#define INTENT_IN ::
#define C_BINDING bind(c)
#define PASS_BY_REF ,pass_by("*")
#else
#define INTENT_IN ,intent(in) ::
#define C_BINDING
#define PASS_BY_REF
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

#ifdef GLU_1_1_TESS
#define glutesselator glutriangulatorobj
#endif

#ifdef NASOFTWARE
#define FULL_TESSCALLBACK
#endif
#ifdef F_LANG
#define FULL_TESSCALLBACK
#endif
#ifdef DIGITAL_NO_ECO_02
#define FULL_TESSCALLBACK
#endif

#ifdef F_LANG
#define PURE_FUNCTIONS
#endif

!--------------------------
! C wrapper procedure names
!--------------------------

public :: &
   fglubegincurve, &
   fglubeginpolygon, &
   fglubeginsurface, &
   fglubegintrim, &
   f9y0glubuild1dmipmaps, &
   f9y1glubuild1dmipmaps, &
   f9y4glubuild1dmipmaps, &
   f9y5glubuild1dmipmaps, &
   f9y6glubuild1dmipmaps, &
   f9y0glubuild2dmipmaps, &
   f9y1glubuild2dmipmaps, &
   f9y4glubuild2dmipmaps, &
   f9y5glubuild2dmipmaps, &
   f9y6glubuild2dmipmaps, &
   fglucylinder, &
   fgludeletenurbsrenderer, &
   fgludeletequadric, &
   fgludeletetess, &
   fgludisk, &
   fgluendcurve, &
   fgluendpolygon, &
   fgluendsurface, &
   fgluendtrim, &
   f9ygluerrorstringlen, &
   fglugetnurbsproperty, &
   f9yglugetstring, &
   f9yglugetstring1, &
   f9yglugetstringlen, &
   f9yglugetstringfree, &
   fglugettessproperty, &
   fgluloadsamplingmatrices, &
   fglulookat, &
   f9yglunewnurbsrenderer, &
   f9yglunewquadric, &
   f9yglunewtess, &
   fglunextcontour
public :: &
   fglunurbscallback, &
   f9y1glunurbscallback, &
   fglunurbscurve, &
   fglunurbsproperty, &
   fglunurbssurface, &
   fgluortho2d, &
   fglupartialdisk, &
   fgluperspective, &
   fglupickmatrix, &
   fgluproject, &
   fglupwlcurve, &
   fgluquadriccallback, &
   f9y1gluquadriccallback, &
   fgluquadricdrawstyle, &
   fgluquadricnormals, &
   fgluquadricorientation, &
   fgluquadrictexture, &
   f9egluquadrictexture
public :: &
   f9y0gluscaleimage, &
   f9y1gluscaleimage, &
   f9y4gluscaleimage, &
   f9y5gluscaleimage, &
   f9y6gluscaleimage, &
   f9y01gluscaleimage, &
   f9y04gluscaleimage, &
   f9y05gluscaleimage, &
   f9y06gluscaleimage, &
   f9y10gluscaleimage, &
   f9y14gluscaleimage, &
   f9y15gluscaleimage, &
   f9y16gluscaleimage, &
   f9y40gluscaleimage, &
   f9y41gluscaleimage, &
   f9y45gluscaleimage, &
   f9y46gluscaleimage, &
   f9y50gluscaleimage, &
   f9y51gluscaleimage, &
   f9y54gluscaleimage, &
   f9y56gluscaleimage, &
   f9y60gluscaleimage, &
   f9y61gluscaleimage, &
   f9y64gluscaleimage, &
   f9y65gluscaleimage, &
   fglusphere, &
   fglutessbegincontour, &
   f9y0glutessbeginpolygon, &
   f9y1glutessbeginpolygon, &
   f9y1glutesscallback
public :: &
#ifdef FULL_TESSCALLBACK
   f9y01glutesscallback, &
   f9y02glutesscallback, &
   f9y03iglutesscallback, &
   f9y03rglutesscallback, &
   f9y04glutesscallback, &
   f9y05glutesscallback, &
   f9y06iglutesscallback, &
   f9y06rglutesscallback, &
   f9y11iglutesscallback, &
   f9y12iglutesscallback, &
   f9y13iglutesscallback, &
   f9y14iglutesscallback, &
   f9y15iglutesscallback, &
   f9y16iglutesscallback, &
   f9y11rglutesscallback, &
   f9y12rglutesscallback, &
   f9y13rglutesscallback, &
   f9y14rglutesscallback, &
   f9y15rglutesscallback, &
   f9y16rglutesscallback, &
#else
   f9y0glutesscallback, &
#endif
   fglutessendcontour, &
   fglutessendpolygon, &
   fglutessnormal, &
   fglutessproperty, &
   f9y0glutessvertex, &
   f9y1glutessvertex, &
   f9y4glutessvertex, &
   f9y5glutessvertex, &
   f9y6glutessvertex, &
   f9y7glutessvertex, &
   fgluunproject

#ifdef F_LANG
! special routines for F passing array args one element at a time

public :: &
   glu_f_uint_array_free, &
   glu_f_uint_array_recv_beg, &
   glu_f_uint_array_recv
#endif

interface
!NAS$ ALIEN "F77 fglubegincurve"
subroutine fglubegincurve(nurb,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglubegincurve
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglubegincurve@8' :: fglubegincurve
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
end subroutine fglubegincurve
end interface

interface
!NAS$ ALIEN "F77 fglubeginpolygon"
subroutine fglubeginpolygon(tess,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglubeginpolygon
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglubeginpolygon@8' :: fglubeginpolygon
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
end subroutine fglubeginpolygon
end interface

interface
!NAS$ ALIEN "F77 fglubeginsurface"
subroutine fglubeginsurface(nurb,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglubeginsurface
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglubeginsurface@8' :: fglubeginsurface
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
end subroutine fglubeginsurface
end interface

interface
!NAS$ ALIEN "F77 fglubegintrim"
subroutine fglubegintrim(nurb,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglubegintrim
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglubegintrim@8' :: fglubegintrim
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
end subroutine fglubegintrim
end interface

interface
!NAS$ ALIEN "F77 f9y0glubuild1dmipmaps"
function f9y0glubuild1dmipmaps(xtarget,component,width,format,xtype, &
   data,length) result(resf9y0glubuild1dmipmaps)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,component,width,format,xtype
!MS$ATTRIBUTES REFERENCE :: data,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glubuild1dmipmaps
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glubuild1dmipmaps@28' :: f9y0glubuild1dmipmaps
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) INTENT_IN length
integer(kind=glint) :: resf9y0glubuild1dmipmaps
end function f9y0glubuild1dmipmaps
end interface

interface
!NAS$ ALIEN "F77 f9y1glubuild1dmipmaps"
function f9y1glubuild1dmipmaps(xtarget,component,width,format,xtype, &
   data) result(resf9y1glubuild1dmipmaps)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,component,width,format,xtype
!MS$ATTRIBUTES REFERENCE :: data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glubuild1dmipmaps
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glubuild1dmipmaps@24' :: f9y1glubuild1dmipmaps
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) :: resf9y1glubuild1dmipmaps
end function f9y1glubuild1dmipmaps
end interface

interface
!NAS$ ALIEN "F77 f9y4glubuild1dmipmaps"
function f9y4glubuild1dmipmaps(xtarget,component,width,format,xtype, &
   temp,nbytes) result(resf9y4glubuild1dmipmaps)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,component,width,format,xtype
!MS$ATTRIBUTES REFERENCE :: temp,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glubuild1dmipmaps
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glubuild1dmipmaps@28' :: f9y4glubuild1dmipmaps
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN temp
integer(kind=glint) INTENT_IN nbytes
integer(kind=glint) :: resf9y4glubuild1dmipmaps
end function f9y4glubuild1dmipmaps
end interface

interface
!NAS$ ALIEN "F77 f9y5glubuild1dmipmaps"
function f9y5glubuild1dmipmaps(xtarget,component,width,format,xtype, &
   data) result(resf9y5glubuild1dmipmaps)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,component,width,format,xtype
!MS$ATTRIBUTES REFERENCE :: data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glubuild1dmipmaps
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glubuild1dmipmaps@24' :: f9y5glubuild1dmipmaps
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) :: resf9y5glubuild1dmipmaps
end function f9y5glubuild1dmipmaps
end interface

interface
!NAS$ ALIEN "F77 f9y6glubuild1dmipmaps"
function f9y6glubuild1dmipmaps(xtarget,component,width,format,xtype, &
   data,length) result(resf9y6glubuild1dmipmaps)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,component,width,format,xtype
!MS$ATTRIBUTES REFERENCE :: data,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glubuild1dmipmaps
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glubuild1dmipmaps@28' :: f9y6glubuild1dmipmaps
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) INTENT_IN length
integer(kind=glint) :: resf9y6glubuild1dmipmaps
end function f9y6glubuild1dmipmaps
end interface

interface
!NAS$ ALIEN "F77 f9y0glubuild2dmipmaps"
function f9y0glubuild2dmipmaps(xtarget,component,width,height,format,xtype, &
   data,length) result(resf9y0glubuild2dmipmaps)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,component,width,height,format,xtype
!MS$ATTRIBUTES REFERENCE :: data,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glubuild2dmipmaps
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glubuild2dmipmaps@32' :: f9y0glubuild2dmipmaps
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width,height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) INTENT_IN length
integer(kind=glint) :: resf9y0glubuild2dmipmaps
end function f9y0glubuild2dmipmaps
end interface

interface
!NAS$ ALIEN "F77 f9y1glubuild2dmipmaps"
function f9y1glubuild2dmipmaps(xtarget,component,width,height,format,xtype, &
   data) result(resf9y1glubuild2dmipmaps)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,component,width,height,format,xtype
!MS$ATTRIBUTES REFERENCE :: data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glubuild2dmipmaps
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glubuild2dmipmaps@28' :: f9y1glubuild2dmipmaps
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width,height
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) :: resf9y1glubuild2dmipmaps
end function f9y1glubuild2dmipmaps
end interface

interface
!NAS$ ALIEN "F77 f9y4glubuild2dmipmaps"
function f9y4glubuild2dmipmaps(xtarget,component,width,height,format,xtype, &
   temp,nbytes) result(resf9y4glubuild2dmipmaps)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,component,width,height,format,xtype
!MS$ATTRIBUTES REFERENCE :: temp,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glubuild2dmipmaps
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glubuild2dmipmaps@32' :: f9y4glubuild2dmipmaps
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width,height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN temp
integer(kind=glint) INTENT_IN nbytes
integer(kind=glint) :: resf9y4glubuild2dmipmaps
end function f9y4glubuild2dmipmaps
end interface

interface
!NAS$ ALIEN "F77 f9y5glubuild2dmipmaps"
function f9y5glubuild2dmipmaps(xtarget,component,width,height,format,xtype, &
   data) result(resf9y5glubuild2dmipmaps)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,component,width,height,format,xtype
!MS$ATTRIBUTES REFERENCE :: data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glubuild2dmipmaps
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glubuild2dmipmaps@28' :: f9y5glubuild2dmipmaps
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width,height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) :: resf9y5glubuild2dmipmaps
end function f9y5glubuild2dmipmaps
end interface

interface
!NAS$ ALIEN "F77 f9y6glubuild2dmipmaps"
function f9y6glubuild2dmipmaps(xtarget,component,width,height,format,xtype, &
   data,length) result(resf9y6glubuild2dmipmaps)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtarget,component,width,height,format,xtype
!MS$ATTRIBUTES REFERENCE :: data,length
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glubuild2dmipmaps
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glubuild2dmipmaps@32' :: f9y6glubuild2dmipmaps
#endif
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width,height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) INTENT_IN length
integer(kind=glint) :: resf9y6glubuild2dmipmaps
end function f9y6glubuild2dmipmaps
end interface

interface
!NAS$ ALIEN "F77 fglucylinder"
subroutine fglucylinder(quad,nbytes,base,top,height,slices,stacks)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes,base,top,height,slices,stacks
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglucylinder
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglucylinder@28' :: fglucylinder
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN quad
integer(kind=glint) INTENT_IN nbytes
real(kind=gldouble) INTENT_IN base, top, height
integer(kind=glint) INTENT_IN slices, stacks
end subroutine fglucylinder
end interface

interface
!NAS$ ALIEN "F77 fgludeletenurbsrenderer"
subroutine fgludeletenurbsrenderer(nurb,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgludeletenurbsrenderer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgludeletenurbsrenderer@8' :: fgludeletenurbsrenderer
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
end subroutine fgludeletenurbsrenderer
end interface

interface
!NAS$ ALIEN "F77 fgludeletequadric"
subroutine fgludeletequadric(quad,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgludeletequadric
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgludeletequadric@8' :: fgludeletequadric
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN quad
integer(kind=glint) INTENT_IN nbytes
end subroutine fgludeletequadric
end interface

interface
!NAS$ ALIEN "F77 fgludeletetess"
subroutine fgludeletetess(tess,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgludeletetess
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgludeletetess@8' :: fgludeletetess
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
end subroutine fgludeletetess
end interface

interface
!NAS$ ALIEN "F77 fgludisk"
subroutine fgludisk(quad,nbytes,inner,outer,slices,loops)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes,inner,outer,slices,loops
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgludisk
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgludisk@24' :: fgludisk
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN quad
integer(kind=glint) INTENT_IN nbytes
real(kind=gldouble) INTENT_IN inner, outer
integer(kind=glint) INTENT_IN slices, loops
end subroutine fgludisk
end interface

interface
!NAS$ ALIEN "F77 fgluendcurve"
subroutine fgluendcurve(nurb,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluendcurve
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluendcurve@8' :: fgluendcurve
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
end subroutine fgluendcurve
end interface

interface
!NAS$ ALIEN "F77 fgluendpolygon"
subroutine fgluendpolygon(tess,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluendpolygon
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluendpolygon@8' :: fgluendpolygon
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
end subroutine fgluendpolygon
end interface

interface
!NAS$ ALIEN "F77 fgluendsurface"
subroutine fgluendsurface(nurb,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluendsurface
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluendsurface@8' :: fgluendsurface
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
end subroutine fgluendsurface
end interface

interface
!NAS$ ALIEN "F77 fgluendtrim"
subroutine fgluendtrim(nurb,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluendtrim
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluendtrim@8' :: fgluendtrim
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
end subroutine fgluendtrim
end interface

interface
!NAS$ ALIEN "F77 f9ygluerrorstringlen"
function f9ygluerrorstringlen(error) result(resf9ygluerrorstringlen)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: error
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9ygluerrorstringlen
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9ygluerrorstringlen@4' :: f9ygluerrorstringlen
#endif
integer(kind=glenum) INTENT_IN error
integer(kind=glint) :: resf9ygluerrorstringlen
end function f9ygluerrorstringlen
end interface

interface
!NAS$ ALIEN "F77 fglugetnurbsproperty"
subroutine fglugetnurbsproperty(nurb,nbytes,property,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes,property,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglugetnurbsproperty
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglugetnurbsproperty@16' :: fglugetnurbsproperty
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN property
real(kind=glfloat), intent(out) :: data
end subroutine fglugetnurbsproperty
end interface

interface
!NAS$ ALIEN "F77 f9yglugetstringlen"
function f9yglugetstringlen(error) result(resf9yglugetstringlen)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: error
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglugetstringlen
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglugetstringlen@4' :: f9yglugetstringlen
#endif
integer(kind=glenum) INTENT_IN error
integer(kind=glint) :: resf9yglugetstringlen
end function f9yglugetstringlen
end interface

interface
!NAS$ ALIEN "F77 f9yglugetstring"
subroutine f9yglugetstring(int_str)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: int_str
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglugetstring
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglugetstring@4' :: f9yglugetstring
#endif
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: int_str
end subroutine f9yglugetstring
end interface

interface
!NAS$ ALIEN "F77 f9yglugetstring1"
function f9yglugetstring1(pos) result(resf9yglugetstring1)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pos
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglugetstring1
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglugetstring1@4' :: f9yglugetstring1
#endif
integer(kind=glenum) INTENT_IN pos
integer(kind=glint) :: resf9yglugetstring1
end function f9yglugetstring1
end interface

interface
!NAS$ ALIEN "F77 f9yglugetstringfree"
#ifdef PURE_FUNCTIONS
function f9yglugetstringfree() result(resf9yglugetstringfree)
#else
subroutine f9yglugetstringfree()
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglugetstringfree
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglugetstringfree@0' :: f9yglugetstringfree
#endif
#ifdef PURE_FUNCTIONS
integer(kind=glint) :: resf9yglugetstringfree
#endif
#ifdef PURE_FUNCTIONS
end function f9yglugetstringfree
#else
end subroutine f9yglugetstringfree
#endif
end interface

interface
!NAS$ ALIEN "F77 fglugettessproperty"
subroutine fglugettessproperty(tess,nbytes,which,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglugettessproperty
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglugettessproperty@16' :: fglugettessproperty
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN which
real(kind=gldouble), intent(out) :: data
end subroutine fglugettessproperty
end interface

interface
!NAS$ ALIEN "F77 fgluloadsamplingmatrices"
subroutine fgluloadsamplingmatrices(nurb,nbytes,model,perspective, &
   view)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes,model,perspective,view
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluloadsamplingmatrices
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluloadsamplingmatrices@20' :: fgluloadsamplingmatrices
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN model, perspective
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN view
end subroutine fgluloadsamplingmatrices
end interface

interface
!NAS$ ALIEN "F77 fglulookat"
subroutine fglulookat(eyex,eyey,eyez,centerx,centery,centerz, &
   upx,upy,upz)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: eyex,eyey,eyez,centerx,centery,centerz
!MS$ATTRIBUTES REFERENCE :: upx,upy,upz
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglulookat
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglulookat@36' :: fglulookat
#endif
real(kind=gldouble) INTENT_IN eyex, eyey, eyez, centerx, centery, centerz, &
    upx, upy, upz
end subroutine fglulookat
end interface

interface
!NAS$ ALIEN "F77 f9yglunewnurbsrenderer"
#ifdef F_LANG
function f9yglunewnurbsrenderer(temp,nbytes) result(resf9yglunewnurbsrenderer)
#else
subroutine f9yglunewnurbsrenderer(temp,nbytes)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglunewnurbsrenderer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglunewnurbsrenderer@8' :: f9yglunewnurbsrenderer
#endif
#ifdef F_LANG
integer(kind=gluint), dimension(DEFER_DIM), intent(in) :: temp
#else
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: temp
#endif
integer(kind=glint) INTENT_IN nbytes
#ifdef F_LANG
integer(kind=glint) :: resf9yglunewnurbsrenderer
#endif
#ifdef F_LANG
end function f9yglunewnurbsrenderer
#else
end subroutine f9yglunewnurbsrenderer
#endif
end interface

interface
!NAS$ ALIEN "F77 f9yglunewquadric"
#ifdef F_LANG
function f9yglunewquadric(temp,nbytes) result(resf9yglunewquadric)
#else
subroutine f9yglunewquadric(temp,nbytes)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglunewquadric
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglunewquadric@8' :: f9yglunewquadric
#endif
#ifdef F_LANG
integer(kind=gluint), dimension(DEFER_DIM), intent(in) :: temp
#else
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: temp
#endif
integer(kind=glint) INTENT_IN nbytes
#ifdef F_LANG
integer(kind=glint) :: resf9yglunewquadric
#endif
#ifdef F_LANG
end function f9yglunewquadric
#else
end subroutine f9yglunewquadric
#endif
end interface

interface
!NAS$ ALIEN "F77 f9yglunewtess"
#ifdef F_LANG
function f9yglunewtess(temp,nbytes) result(resf9yglunewtess)
#else
subroutine f9yglunewtess(temp,nbytes)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglunewtess
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglunewtess@8' :: f9yglunewtess
#endif
#ifdef F_LANG
integer(kind=gluint), dimension(DEFER_DIM), intent(in) :: temp
#else
integer(kind=gluint), dimension(DEFER_DIM), intent(out) :: temp
#endif
integer(kind=glint) INTENT_IN nbytes
#ifdef F_LANG
integer(kind=glint) :: resf9yglunewtess
#endif
#ifdef F_LANG
end function f9yglunewtess
#else
end subroutine f9yglunewtess
#endif
end interface

interface
!NAS$ ALIEN "F77 fglunextcontour"
subroutine fglunextcontour(tess,nbytes,xtype)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,xtype
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglunextcontour
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglunextcontour@12' :: fglunextcontour
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN xtype
end subroutine fglunextcontour
end interface

interface
!NAS$ ALIEN "F77 f9y1glunurbscallback"
subroutine f9y1glunurbscallback(nurb,nbytes,which)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes,which
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glunurbscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glunurbscallback@12' :: f9y1glunurbscallback
#endif
integer(kind=gluint), dimension(DEFER_DIM), intent(in out) :: nurb
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN which
end subroutine f9y1glunurbscallback
end interface

interface
!NAS$ ALIEN "F77 fglunurbscurve"
subroutine fglunurbscurve(nurb,nbytes,knotcount,knots,stride,control, &
   order,xtype)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes,knotcount,knots,stride,control
!MS$ATTRIBUTES REFERENCE :: order,xtype
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglunurbscurve
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglunurbscurve@32' :: fglunurbscurve
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
integer(kind=glint) INTENT_IN knotcount
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: knots
integer(kind=glint) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: control
integer(kind=glint) INTENT_IN order
integer(kind=glenum) INTENT_IN xtype
end subroutine fglunurbscurve
end interface

interface
!NAS$ ALIEN "F77 fglunurbsproperty"
subroutine fglunurbsproperty(nurb,nbytes,property,value)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes,property,value
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglunurbsproperty
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglunurbsproperty@16' :: fglunurbsproperty
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN property
real(kind=glfloat) INTENT_IN value
end subroutine fglunurbsproperty
end interface

interface
!NAS$ ALIEN "F77 fglunurbssurface"
subroutine fglunurbssurface(nurb,nbytes,sknotcount,sknots,tknotcount, &
   tknots,sstride,tstride,control,sorder,torder,xtype)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes,sknotcount,sknots,tknotcount,tknots
!MS$ATTRIBUTES REFERENCE :: sstride,tstride,control,sorder,torder,xtype
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglunurbssurface
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglunurbssurface@48' :: fglunurbssurface
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
integer(kind=glint) INTENT_IN sknotcount
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: sknots
integer(kind=glint) INTENT_IN tknotcount
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: tknots
integer(kind=glint) INTENT_IN sstride, tstride
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: control
integer(kind=glint) INTENT_IN sorder, torder
integer(kind=glenum) INTENT_IN xtype
end subroutine fglunurbssurface
end interface

interface
!NAS$ ALIEN "F77 fgluortho2d"
subroutine fgluortho2d(left,right,bottom,top)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: left,right,bottom,top
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluortho2d
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluortho2d@16' :: fgluortho2d
#endif
real(kind=gldouble) INTENT_IN left, right, bottom, top
end subroutine fgluortho2d
end interface

interface
!NAS$ ALIEN "F77 fglupartialdisk"
subroutine fglupartialdisk(quad,nbytes,inner,outer,slices,loops,start, &
   sweep)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes,inner,outer,slices,loops,start
!MS$ATTRIBUTES REFERENCE :: sweep
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglupartialdisk
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglupartialdisk@32' :: fglupartialdisk
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN quad
integer(kind=glint) INTENT_IN nbytes
real(kind=gldouble) INTENT_IN inner, outer
integer(kind=glint) INTENT_IN slices, loops
real(kind=gldouble) INTENT_IN start, sweep
end subroutine fglupartialdisk
end interface

interface
!NAS$ ALIEN "F77 fgluperspective"
subroutine fgluperspective(fovy,aspect,znear,zfar)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: fovy,aspect,znear,zfar
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluperspective
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluperspective@16' :: fgluperspective
#endif
real(kind=gldouble) INTENT_IN fovy, aspect, znear, zfar
end subroutine fgluperspective
end interface

interface
!NAS$ ALIEN "F77 fglupickmatrix"
subroutine fglupickmatrix(x,y,delx,dely,viewport)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,delx,dely,viewport
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglupickmatrix
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglupickmatrix@20' :: fglupickmatrix
#endif
real(kind=gldouble) INTENT_IN x, y, delx, dely
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN viewport
end subroutine fglupickmatrix
end interface

interface
!NAS$ ALIEN "F77 fgluproject"
#ifdef PURE_FUNCTIONS
subroutine fgluproject(objx,objy,objz,model,proj,view,winx,winy, &
   winz,errcode)
#else
function fgluproject(objx,objy,objz,model,proj,view,winx,winy, &
   winz) result(resfgluproject)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: objx,objy,objz,model,proj,view,winx,winy
!MS$ATTRIBUTES REFERENCE :: winz
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluproject
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluproject@36' :: fgluproject
#endif
real(kind=gldouble) INTENT_IN objx, objy, objz
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN model, proj
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN view
real(kind=gldouble), intent(out) :: winx, winy, winz
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resfgluproject
#endif
#ifdef PURE_FUNCTIONS
end subroutine fgluproject
#else
end function fgluproject
#endif
end interface

interface
!NAS$ ALIEN "F77 fglupwlcurve"
subroutine fglupwlcurve(nurb,nbytes,xcount,data,stride,xtype)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes,xcount,data,stride,xtype
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglupwlcurve
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglupwlcurve@24' :: fglupwlcurve
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN nurb
integer(kind=glint) INTENT_IN nbytes
integer(kind=glint) INTENT_IN xcount
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: data
integer(kind=glint) INTENT_IN stride
integer(kind=glenum) INTENT_IN xtype
end subroutine fglupwlcurve
end interface

interface
!NAS$ ALIEN "F77 f9y1gluquadriccallback"
subroutine f9y1gluquadriccallback(quad,nbytes,which)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes,which
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gluquadriccallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gluquadriccallback@12' :: f9y1gluquadriccallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: quad
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN which
end subroutine f9y1gluquadriccallback
end interface

interface
!NAS$ ALIEN "F77 fgluquadricdrawstyle"
subroutine fgluquadricdrawstyle(quad,nbytes,draw)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes,draw
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluquadricdrawstyle
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluquadricdrawstyle@12' :: fgluquadricdrawstyle
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN quad
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN draw
end subroutine fgluquadricdrawstyle
end interface

interface
!NAS$ ALIEN "F77 fgluquadricnormals"
subroutine fgluquadricnormals(quad,nbytes,normal)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes,normal
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluquadricnormals
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluquadricnormals@12' :: fgluquadricnormals
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN quad
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN normal
end subroutine fgluquadricnormals
end interface

interface
!NAS$ ALIEN "F77 fgluquadricorientation"
subroutine fgluquadricorientation(quad,nbytes,orientation)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes,orientation
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluquadricorientation
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluquadricorientation@12' :: fgluquadricorientation
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN quad
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN orientation
end subroutine fgluquadricorientation
end interface

interface
!NAS$ ALIEN "F77 fgluquadrictexture"
subroutine fgluquadrictexture(quad,nbytes,texture)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes,texture
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluquadrictexture
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluquadrictexture@12' :: fgluquadrictexture
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN quad
integer(kind=glint) INTENT_IN nbytes
logical(kind=glboolean) INTENT_IN texture
end subroutine fgluquadrictexture
end interface

interface
!NAS$ ALIEN "F77 f9egluquadrictexture"
subroutine f9egluquadrictexture(quad,nbytes,texture)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes,texture
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9egluquadrictexture
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9egluquadrictexture@12' :: f9egluquadrictexture
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN quad
integer(kind=glint) INTENT_IN nbytes
integer(kind=glint) INTENT_IN texture
end subroutine f9egluquadrictexture
end interface

interface
!NAS$ ALIEN "F77 f9y0gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y0gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,lengthout,errcode)
#else
function f9y0gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,lengthout) result(resf9y0gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout,lengthout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0gluscaleimage@44' :: f9y0gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: dataout
integer(kind=glint) INTENT_IN lengthout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y0gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y0gluscaleimage
#else
end function f9y0gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y1gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y1gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,errcode)
#else
function f9y1gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout) result(resf9y1gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gluscaleimage@36' :: f9y1gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y1gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y1gluscaleimage
#else
end function f9y1gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y4gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y4gluscaleimage(format,win,hin,typein,wout,hout, &
   typeout,temp,nbytes,errcode)
#else
function f9y4gluscaleimage(format,win,hin,typein,wout,hout, &
   typeout,temp,nbytes) result(resf9y4gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,temp,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4gluscaleimage@36' :: f9y4gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN temp
integer(kind=glint) INTENT_IN nbytes
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y4gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y4gluscaleimage
#else
end function f9y4gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y5gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y5gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,errcode)
#else
function f9y5gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout) result(resf9y5gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5gluscaleimage@36' :: f9y5gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y5gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y5gluscaleimage
#else
end function f9y5gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y6gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y6gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,lengthout,errcode)
#else
function f9y6gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,lengthout) result(resf9y6gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout,lengthout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6gluscaleimage@44' :: f9y6gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glshort), dimension(DEFER_DIM), intent(out) :: dataout
integer(kind=glint) INTENT_IN lengthout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y6gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y6gluscaleimage
#else
end function f9y6gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y01gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y01gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,errcode)
#else
function f9y01gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout) result(resf9y01gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y01gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y01gluscaleimage@40' :: f9y01gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y01gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y01gluscaleimage
#else
end function f9y01gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y04gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y04gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,temp,nbytes,errcode)
#else
function f9y04gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,temp,nbytes) result(resf9y04gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,temp,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y04gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y04gluscaleimage@44' :: f9y04gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN temp
integer(kind=glint) INTENT_IN nbytes
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y04gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y04gluscaleimage
#else
end function f9y04gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y05gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y05gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,errcode)
#else
function f9y05gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout) result(resf9y05gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y05gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y05gluscaleimage@40' :: f9y05gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y05gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y05gluscaleimage
#else
end function f9y05gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y06gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y06gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,lengthout,errcode)
#else
function f9y06gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,lengthout) result(resf9y06gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout,lengthout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y06gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y06gluscaleimage@44' :: f9y06gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glshort), dimension(DEFER_DIM), intent(out) :: dataout
integer(kind=glint) INTENT_IN lengthout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y06gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y06gluscaleimage
#else
end function f9y06gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y10gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y10gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,lengthout,errcode)
#else
function f9y10gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,lengthout) result(resf9y10gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout,lengthout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y10gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y10gluscaleimage@40' :: f9y10gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: dataout
integer(kind=glint) INTENT_IN lengthout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y10gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y10gluscaleimage
#else
end function f9y10gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y14gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y14gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,temp,nbytes,errcode)
#else
function f9y14gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,temp,nbytes) result(resf9y14gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,temp,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y14gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y14gluscaleimage@40' :: f9y14gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN temp
integer(kind=glint) INTENT_IN nbytes
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y14gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y14gluscaleimage
#else
end function f9y14gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y15gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y15gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,errcode)
#else
function f9y15gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout) result(resf9y15gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y15gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y15gluscaleimage@36' :: f9y15gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y15gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y15gluscaleimage
#else
end function f9y15gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y16gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y16gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,lengthout,errcode)
#else
function f9y16gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,lengthout) result(resf9y16gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout,lengthout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y16gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y16gluscaleimage@40' :: f9y16gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glshort), dimension(DEFER_DIM), intent(out) :: dataout
integer(kind=glint) INTENT_IN lengthout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y16gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y16gluscaleimage
#else
end function f9y16gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y40gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y40gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,lengthout,errcode)
#else
function f9y40gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,lengthout) result(resf9y40gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout,lengthout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y40gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y40gluscaleimage@44' :: f9y40gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: dataout
integer(kind=glint) INTENT_IN lengthout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y40gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y40gluscaleimage
#else
end function f9y40gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y41gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y41gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,errcode)
#else
function f9y41gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout) result(resf9y41gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y41gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y41gluscaleimage@40' :: f9y01gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y41gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y41gluscaleimage
#else
end function f9y41gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y45gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y45gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,errcode)
#else
function f9y45gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout) result(resf9y45gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y45gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y45gluscaleimage@40' :: f9y45gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y45gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y45gluscaleimage
#else
end function f9y45gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y46gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y46gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,lengthout,errcode)
#else
function f9y46gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,lengthout) result(resf9y46gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout,lengthout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y46gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y46gluscaleimage@44' :: f9y46gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glshort), dimension(DEFER_DIM), intent(out) :: dataout
integer(kind=glint) INTENT_IN lengthout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y46gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y46gluscaleimage
#else
end function f9y46gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y50gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y50gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,lengthout,errcode)
#else
function f9y50gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,lengthout) result(resf9y50gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout,lengthout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y50gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y50gluscaleimage@40' :: f9y50gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: dataout
integer(kind=glint) INTENT_IN lengthout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y50gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y50gluscaleimage
#else
end function f9y50gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y51gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y51gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,errcode)
#else
function f9y51gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout) result(resf9y51gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y51gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y51gluscaleimage@36' :: f9y51gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y51gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y51gluscaleimage
#else
end function f9y51gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y54gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y54gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,temp,nbytes,errcode)
#else
function f9y54gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,temp,nbytes) result(resf9y54gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,temp,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y54gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y54gluscaleimage@40' :: f9y54gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN temp
integer(kind=glint) INTENT_IN nbytes
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y54gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y54gluscaleimage
#else
end function f9y54gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y56gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y56gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,lengthout,errcode)
#else
function f9y56gluscaleimage(format,win,hin,typein,datain,wout,hout, &
   typeout,dataout,lengthout) result(resf9y56gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout,lengthout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y56gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y56gluscaleimage@40' :: f9y56gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glshort), dimension(DEFER_DIM), intent(out) :: dataout
integer(kind=glint) INTENT_IN lengthout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y56gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y56gluscaleimage
#else
end function f9y56gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y60gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y60gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,errcode)
#else
function f9y60gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout) result(resf9y60gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y60gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y60gluscaleimage@40' :: f9y60gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glint), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y60gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y60gluscaleimage
#else
end function f9y60gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y61gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y61gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,errcode)
#else
function f9y61gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout) result(resf9y61gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y61gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y61gluscaleimage@40' :: f9y61gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y61gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y61gluscaleimage
#else
end function f9y61gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y64gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y64gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,temp,nbytes,errcode)
#else
function f9y64gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,temp,nbytes) result(resf9y64gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,temp,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y64gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y64gluscaleimage@44' :: f9y64gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=gluint), dimension(DEFER_DIM) INTENT_IN temp
integer(kind=glint) INTENT_IN nbytes
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y64gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y64gluscaleimage
#else
end function f9y64gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 f9y65gluscaleimage"
#ifdef PURE_FUNCTIONS
subroutine f9y65gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout,errcode)
#else
function f9y65gluscaleimage(format,win,hin,typein,datain,lengthin,wout,hout, &
   typeout,dataout) result(resf9y65gluscaleimage)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: format,win,hin,typein,datain,lengthin
!MS$ATTRIBUTES REFERENCE :: wout,hout,typeout,dataout
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y65gluscaleimage
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y65gluscaleimage@40' :: f9y65gluscaleimage
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glshort), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glint) INTENT_IN lengthin
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9y65gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
end subroutine f9y65gluscaleimage
#else
end function f9y65gluscaleimage
#endif
end interface

interface
!NAS$ ALIEN "F77 fglusphere"
subroutine fglusphere(quad,nbytes,radius,slices,stacks)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes,radius,slices,stacks
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglusphere
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglusphere@20' :: fglusphere
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN quad
integer(kind=glint) INTENT_IN nbytes
real(kind=gldouble) INTENT_IN radius
integer(kind=glint) INTENT_IN slices, stacks
end subroutine fglusphere
end interface

interface
!NAS$ ALIEN "F77 fglutessbegincontour"
subroutine fglutessbegincontour(tess,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutessbegincontour
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutessbegincontour@8' :: fglutessbegincontour
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
end subroutine fglutessbegincontour
end interface

interface
!NAS$ ALIEN "F77 f9y0glutessbeginpolygon"
subroutine f9y0glutessbeginpolygon(tess,nbytes,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glutessbeginpolygon
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glutessbeginpolygon@12' :: f9y0glutessbeginpolygon
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
integer(kind=glint), dimension(DEFER_DIM), TARGET intent(in out) :: data
end subroutine f9y0glutessbeginpolygon
end interface

interface
!NAS$ ALIEN "F77 f9y1glutessbeginpolygon"
subroutine f9y1glutessbeginpolygon(tess,nbytes,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutessbeginpolygon
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutessbeginpolygon@12' :: f9y1glutessbeginpolygon
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: data
end subroutine f9y1glutessbeginpolygon
end interface

interface
!NAS$ ALIEN "F77 f9y1glutesscallback"
subroutine f9y1glutesscallback(tess,nbytes,which)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutesscallback@12' :: f9y1glutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN which
end subroutine f9y1glutesscallback
end interface

interface
!NAS$ ALIEN "F77 fglutessendcontour"
subroutine fglutessendcontour(tess,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutessendcontour
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutessendcontour@8' :: fglutessendcontour
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
end subroutine fglutessendcontour
end interface

interface
!NAS$ ALIEN "F77 fglutessendpolygon"
subroutine fglutessendpolygon(tess,nbytes)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutessendpolygon
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutessendpolygon@8' :: fglutessendpolygon
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
end subroutine fglutessendpolygon
end interface

interface
!NAS$ ALIEN "F77 fglutessnormal"
subroutine fglutessnormal(tess,nbytes,valuex,valuey,valuez)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,valuex,valuey,valuez
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutessnormal
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutessnormal@20' :: fglutessnormal
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
real(kind=gldouble) INTENT_IN valuex, valuey, valuez
end subroutine fglutessnormal
end interface

interface
!NAS$ ALIEN "F77 fglutessproperty"
subroutine fglutessproperty(tess,nbytes,which,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutessproperty
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutessproperty@16' :: fglutessproperty
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN which
real(kind=gldouble) INTENT_IN data
end subroutine fglutessproperty
end interface

interface
!NAS$ ALIEN "F77 f9y0glutessvertex"
subroutine f9y0glutessvertex(tess,nbytes,location,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,location,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glutessvertex
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glutessvertex@16' :: f9y0glutessvertex
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
integer(kind=glint), dimension(DEFER_DIM), TARGET intent(in out) :: data
end subroutine f9y0glutessvertex
end interface

interface
!NAS$ ALIEN "F77 f9y1glutessvertex"
subroutine f9y1glutessvertex(tess,nbytes,location,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,location,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutessvertex
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutessvertex@16' :: f9y1glutessvertex
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
real(kind=glfloat), dimension(DEFER_DIM), TARGET intent(in out) :: data
end subroutine f9y1glutessvertex
end interface

interface
!NAS$ ALIEN "F77 f9y4glutessvertex"
subroutine f9y4glutessvertex(tess,nbytes,location)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,location
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y4glutessvertex
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y4glutessvertex@12' :: f9y4glutessvertex
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
end subroutine f9y4glutessvertex
end interface

interface
!NAS$ ALIEN "F77 f9y5glutessvertex"
subroutine f9y5glutessvertex(tess,nbytes,location,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,location,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y5glutessvertex
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y5glutessvertex@16' :: f9y5glutessvertex
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
integer(kind=glbyte), dimension(DEFER_DIM), TARGET intent(in out) :: data
end subroutine f9y5glutessvertex
end interface

interface
!NAS$ ALIEN "F77 f9y6glutessvertex"
subroutine f9y6glutessvertex(tess,nbytes,location,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,location,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y6glutessvertex
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y6glutessvertex@16' :: f9y6glutessvertex
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
integer(kind=glshort), dimension(DEFER_DIM), TARGET intent(in out) :: data
end subroutine f9y6glutessvertex
end interface

interface
!NAS$ ALIEN "F77 f9y7glutessvertex"
subroutine f9y7glutessvertex(tess,nbytes,location,data)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,location,data
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y7glutessvertex
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y7glutessvertex@16' :: f9y7glutessvertex
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN tess
integer(kind=glint) INTENT_IN nbytes
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
real(kind=gldouble), dimension(DEFER_DIM), TARGET intent(in out) :: data
end subroutine f9y7glutessvertex
end interface

interface
!NAS$ ALIEN "F77 fgluunproject"
#ifdef PURE_FUNCTIONS
subroutine fgluunproject(winx,winy,winz,model,proj,view,objx, &
   objy,objz,errcode)
#else
function fgluunproject(winx,winy,winz,model,proj,view,objx, &
   objy,objz) result(resfgluunproject)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: winx,winy,winz,model,proj,view,objx
!MS$ATTRIBUTES REFERENCE :: objy,objz
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluunproject
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluunproject@36' :: fgluunproject
#endif
real(kind=gldouble) INTENT_IN winx, winy, winz
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN model, proj
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN view
real(kind=gldouble), intent(out) :: objx, objy, objz
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resfgluunproject
#endif
#ifdef PURE_FUNCTIONS
end subroutine fgluunproject
#else
end function fgluunproject
#endif
end interface

#ifdef F_LANG
! special routines for F for passing array arguments one element at a time

interface
function glu_f_uint_array_free() result(fres)
use opengl_kinds
IMPLICIT_NONE
integer(kind=glint) :: fres
end function glu_f_uint_array_free
end interface

interface
function glu_f_uint_array_recv_beg() result(fres)
use opengl_kinds
IMPLICIT_NONE
integer(kind=glint) :: fres
end function glu_f_uint_array_recv_beg
end interface

interface
function glu_f_uint_array_recv(pos) result(fres)
use opengl_kinds
IMPLICIT_NONE
integer(kind=glint), intent(in) :: pos
integer(kind=gluint) :: fres
end function glu_f_uint_array_recv
end interface

#endif

interface
#ifdef NASOFTWARE
bind(c,name="fglunurbscallback_") &
#endif
subroutine fglunurbscallback(nurb,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: nurb,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglunurbscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglunurbscallback@16' :: fglunurbscallback
#endif
integer(kind=gluint), dimension(DEFER_DIM), intent(in out) :: nurb
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(errno)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   end subroutine callbackfunc
end interface
end subroutine fglunurbscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,name="fgluquadriccallback_") &
#endif
subroutine fgluquadriccallback(quad,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: quad,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluquadriccallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluquadriccallback@16' :: fgluquadriccallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: quad
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(errno)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   end subroutine callbackfunc
end interface
end subroutine fgluquadriccallback
end interface

#ifdef FULL_TESSCALLBACK

interface
#ifdef NASOFTWARE
bind(c,"f9y01glutesscallback_") &
#endif
subroutine f9y01glutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y01glutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y01glutesscallback@16' :: f9y01glutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(xtype)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: xtype
   end subroutine callbackfunc
end interface
end subroutine f9y01glutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y02glutesscallback_") &
#endif
subroutine f9y02glutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y02glutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y02glutesscallback@16' :: f9y02glutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(flag)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint) PASS_BY_REF, intent(in out) :: flag
   end subroutine callbackfunc
end interface
end subroutine f9y02glutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y03iglutesscallback_") &
#endif
subroutine f9y03iglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y03iglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y03iglutesscallback@16' :: f9y03iglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(vertex_data)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: vertex_data
   end subroutine callbackfunc
end interface
end subroutine f9y03iglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y03rglutesscallback_") &
#endif
subroutine f9y03rglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y03rglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y03rglutesscallback@16' :: f9y03rglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(vertex_data)
   use opengl_kinds
   IMPLICIT_NONE
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: vertex_data
   end subroutine callbackfunc
end interface
end subroutine f9y03rglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y04glutesscallback_") &
#endif
subroutine f9y04glutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y04glutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y04glutesscallback@16' :: f9y04glutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc()
   use opengl_kinds
   IMPLICIT_NONE
   end subroutine callbackfunc
end interface
end subroutine f9y04glutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y05glutesscallback_") &
#endif
subroutine f9y05glutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y05glutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y05glutesscallback@16' :: f9y05glutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(errno)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   end subroutine callbackfunc
end interface
end subroutine f9y05glutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y06iglutesscallback_") &
#endif
subroutine f9y06iglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y06iglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y06iglutesscallback@16' :: f9y06iglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(coords, vertex_data, weight, outdata)
   use opengl_kinds
   IMPLICIT_NONE
   real(kind=gldouble), dimension(DEFER_DIM), intent(in out) :: coords
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: vertex_data, outdata
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: weight
   end subroutine callbackfunc
end interface
end subroutine f9y06iglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y06rglutesscallback_") &
#endif
subroutine f9y06rglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y06rglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y06rglutesscallback@16' :: f9y06rglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(coords, vertex_data, weight, outdata)
   use opengl_kinds
   IMPLICIT_NONE
   real(kind=gldouble), dimension(DEFER_DIM), intent(in out) :: coords
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: vertex_data, weight, &
                                                       outdata
   end subroutine callbackfunc
end interface
end subroutine f9y06rglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y11iglutesscallback_") &
#endif
subroutine f9y11iglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y11iglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y11iglutesscallback@16' :: f9y11iglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(xtype,polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: xtype
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y11iglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y12iglutesscallback_") &
#endif
subroutine f9y12iglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y12iglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y12iglutesscallback@16' :: f9y12iglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(flag,polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint) PASS_BY_REF, intent(in out) :: flag
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y12iglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y13iglutesscallback_") &
#endif
subroutine f9y13iglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y13iglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y13iglutesscallback@16' :: f9y13iglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(vertex_data,polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: vertex_data
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y13iglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y14iglutesscallback_") &
#endif
subroutine f9y14iglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y14iglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y14iglutesscallback@16' :: f9y14iglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y14iglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y15iglutesscallback_") &
#endif
subroutine f9y15iglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y15iglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y15iglutesscallback@16' :: f9y15iglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(errno,polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y15iglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y16iglutesscallback_") &
#endif
subroutine f9y16iglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y16iglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y16iglutesscallback@16' :: f9y16iglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(coords, vertex_data, weight, outdata, polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   real(kind=gldouble), dimension(DEFER_DIM), intent(in out) :: coords
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: vertex_data, &
                                                                outdata
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: weight
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y16iglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y11rglutesscallback_") &
#endif
subroutine f9y11rglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y11rglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y11rglutesscallback@16' :: f9y11rglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(xtype,polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: xtype
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y11rglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y12rglutesscallback_") &
#endif
subroutine f9y12rglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y12rglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y12rglutesscallback@16' :: f9y12rglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(flag,polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint) PASS_BY_REF, intent(in out) :: flag
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y12rglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y13rglutesscallback_") &
#endif
subroutine f9y13rglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y13rglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y13rglutesscallback@16' :: f9y13rglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(vertex_data,polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: vertex_data
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y13rglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y14rglutesscallback_") &
#endif
subroutine f9y14rglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y14rglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y14rglutesscallback@16' :: f9y14rglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y14rglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y15rglutesscallback_") &
#endif
subroutine f9y15rglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y15rglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y15rglutesscallback@16' :: f9y15rglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(errno,polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y15rglutesscallback
end interface

interface
#ifdef NASOFTWARE
bind(c,"f9y16rglutesscallback_") &
#endif
subroutine f9y16rglutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y16rglutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y16rglutesscallback@16' :: f9y16rglutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) PASS_BY_REF INTENT_IN nbytes
integer(kind=glenum) PASS_BY_REF INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(coords, vertex_data, weight, outdata, polygon_data)
   use opengl_kinds
   IMPLICIT_NONE
   real(kind=gldouble), dimension(DEFER_DIM), intent(in out) :: coords
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: vertex_data, weight, &
                                                       outdata
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
end subroutine f9y16rglutesscallback
end interface

#else

interface
#ifdef NASOFTWARE
bind(c,"f9y0glutesscallback_") &
#endif
subroutine f9y0glutesscallback(tess,nbytes,which,callbackfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tess,nbytes,which,callbackfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y0glutesscallback
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y0glutesscallback@16' :: f9y0glutesscallback
#endif
integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: tess
integer(kind=glint) INTENT_IN nbytes
integer(kind=glenum) INTENT_IN which
external callbackfunc
end subroutine f9y0glutesscallback
end interface

#endif

#endif

end module opengl_gluinterfaces
