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

module opengl_glu

use opengl_kinds
use opengl_gluinterfaces
implicit none
private

#ifdef NASOFTWARE
#define INTENT_IN ::
#define C_BINDING bind(c)
#define PASS_BY_REF ,pass_by("*")
#else 
#define INTENT_IN ,intent(in) :: 
#define C_BINDING
#define PASS_BY_REF
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

#ifdef DONT_PASS_COMMENTS
! GLUNULLFUNC used to be an integer, and a generic interface that
! distinguished between integer and procedure arguments determined
! if NULL should be passed as the callback function.  This is not
! standard conforming, so the approach has been changed.  The old
! approach has been kept in the code in case the new approach doesn't
! work with some compilers.  OLDNULLFUNC determines which approach to
! use, and gets #defined here for any that need the old approach.
#endif
#ifdef NASOFTWARE
#define OLDNULLFUNC
#endif
#ifdef EPC
#define OLDNULLFUNC
#endif
#ifdef DVF
#define OLDNULLFUNC
#endif
#ifdef COMPAQ
#define OLDNULLFUNC
#endif
#ifdef MIPSPRO
#define OLDNULLFUNC
#endif
#ifdef F_LANG
#define OLDNULLFUNC
#endif
#ifdef GNU
#define OLDNULLFUNC
#endif

#ifdef DONT_PASS_COMMENTS
! If you can't find a way to get the compiler to match the interface
! blocks for the callback functions, declare them external.
#endif
#ifdef NAGWARE
#define EXTERN_CALLBACK
#endif
#ifdef INTEL
#define EXTERN_CALLBACK
#endif

#ifdef DONT_PASS_COMMENTS
! Compilers that can't distinguish integer from dummy procedure in a
! generic interface define NO_NULLFUNC (only needed with OLDNULLFUNC)
#ifdef F_LANG
#define NO_NULLFUNC
#endif
#endif

!--------------------------
! Public generic interfaces
!--------------------------

public :: &
#ifndef OLDNULLFUNC
   glunullfunc, &
#endif
   glubegincurve, &
   glubeginpolygon, &
   glubeginsurface, &
   glubegintrim, &
   glubuild1dmipmaps, &
   glubuild2dmipmaps, &
   glucylinder, &
   gludeletenurbsrenderer, &
   gludeletequadric, &
   gludeletetess, &
   gludisk, &
   gluendcurve, &
   gluendpolygon, &
   gluendsurface, &
   gluendtrim, &
   gluerrorstring, &
   glugetnurbsproperty, &
   glugetstring, &
#ifdef GLU_1_1_TESS
#else
   glugettessproperty, &
#endif
   gluloadsamplingmatrices, &
   glulookat, &
   glunewnurbsrenderer, &
   glunewquadric, &
   glunewtess, &
   glunextcontour
public :: &
   glunurbscallback, &
   glunurbscurve, &
   glunurbsproperty, &
   glunurbssurface, &
   gluortho2d, &
   glupartialdisk, &
   gluperspective, &
   glupickmatrix, &
   gluproject, &
   glupwlcurve, &
   gluquadriccallback, &
   gluquadricdrawstyle, &
   gluquadricnormals, &
   gluquadricorientation, &
   gluquadrictexture, &
   gluscaleimage, &
   glusphere, &
#ifdef GLU_1_1_TESS
#else
   glutessbegincontour, &
   glutessbeginpolygon, &
#endif
   glutesscallback, &
#ifdef GLU_1_1_TESS
#else
   glutessendcontour, &
   glutessendpolygon, &
   glutessnormal, &
   glutessproperty, &
#endif
   glutessvertex, &
   gluunproject

!-----------------------------------
! Public entities from other modules
!-----------------------------------

! Some older compilers require that derived types from other modules be made
! public if there are public routines with arguments of that type.
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

private :: &
   f9xglubegincurve, &
   f9xglubeginpolygon, &
   f9xglubeginsurface, &
   f9xglubegintrim, &
   f9x0glubuild1dmipmaps, &
   f9x1glubuild1dmipmaps, &
   f9x4glubuild1dmipmaps, &
   f9x5glubuild1dmipmaps, &
   f9x6glubuild1dmipmaps, &
   f9x0glubuild2dmipmaps, &
   f9x1glubuild2dmipmaps, &
   f9x4glubuild2dmipmaps, &
   f9x5glubuild2dmipmaps, &
   f9x6glubuild2dmipmaps, &
   f9xglucylinder, &
   f9xgludeletenurbsrenderer, &
   f9xgludeletequadric, &
   f9xgludeletetess, &
   f9xgludisk, &
   f9xgluendcurve, &
   f9xgluendpolygon, &
   f9xgluendsurface, &
   f9xgluendtrim, &
   f9xgluerrorstring, &
   f9xglugetnurbsproperty, &
   f9xglugetstring, &
#ifdef GLU_1_1_TESS
#else
   f9xglugettessproperty, &
#endif
   f9xgluloadsamplingmatrices, &
   f9xglulookat, &
   f9xglunewnurbsrenderer, &
   f9xglunewquadric, &
   f9xglunewtess, &
   f9xglunextcontour
private :: &
   f9xglunurbscallback, &
#ifdef OLDNULLFUNC
   f9x1glunurbscallback, &
#endif
   f9xglunurbscurve, &
   f9xglunurbsproperty, &
   f9xglunurbssurface, &
   f9xgluortho2d, &
   f9xglupartialdisk, &
   f9xgluperspective, &
   f9xglupickmatrix, &
   f9xgluproject, &
   f9xglupwlcurve, &
   f9xgluquadriccallback, &
#ifdef OLDNULLFUNC
   f9x1gluquadriccallback, &
#endif
   f9xgluquadricdrawstyle, &
   f9xgluquadricnormals, &
   f9xgluquadricorientation, &
   f9xgluquadrictexture, &
   f9x0gluscaleimage, &
   f9x1gluscaleimage, &
   f9x4gluscaleimage, &
   f9x5gluscaleimage, &
   f9x6gluscaleimage, &
   f9xglusphere
private :: &
   f9x01gluscaleimage, &
   f9x04gluscaleimage, &
   f9x05gluscaleimage, &
   f9x06gluscaleimage, &
   f9x10gluscaleimage, &
   f9x14gluscaleimage, &
   f9x15gluscaleimage, &
   f9x16gluscaleimage, &
   f9x40gluscaleimage, &
   f9x41gluscaleimage, &
   f9x45gluscaleimage, &
   f9x46gluscaleimage, &
   f9x50gluscaleimage, &
   f9x51gluscaleimage, &
   f9x54gluscaleimage, &
   f9x56gluscaleimage, &
   f9x60gluscaleimage, &
   f9x61gluscaleimage, &
   f9x64gluscaleimage, &
   f9x65gluscaleimage
private :: &
#ifdef GLU_1_1_TESS
#else
   f9xglutessbegincontour, &
   f9x0glutessbeginpolygon, &
   f9x1glutessbeginpolygon, &
#endif
#ifdef FULL_TESSCALLBACK
   f9x01glutesscallback, &
   f9x02glutesscallback, &
   f9x03iglutesscallback, &
   f9x03rglutesscallback, &
   f9x04glutesscallback, &
   f9x05glutesscallback, &
   f9x06iglutesscallback, &
   f9x06rglutesscallback, &
   f9x11iglutesscallback, &
   f9x12iglutesscallback, &
   f9x13iglutesscallback, &
   f9x14iglutesscallback, &
   f9x15iglutesscallback, &
   f9x16iglutesscallback, &
   f9x11rglutesscallback, &
   f9x12rglutesscallback, &
   f9x13rglutesscallback, &
   f9x14rglutesscallback, &
   f9x15rglutesscallback, &
   f9x16rglutesscallback
#ifdef OLDNULLFUNC
private :: &
   f9x01aglutesscallback, &
   f9x02aglutesscallback, &
   f9x03iaglutesscallback, &
   f9x03raglutesscallback, &
   f9x04aglutesscallback, &
   f9x05aglutesscallback, &
   f9x06iaglutesscallback, &
   f9x06raglutesscallback, &
   f9x11iaglutesscallback, &
   f9x12iaglutesscallback, &
   f9x13iaglutesscallback, &
   f9x14iaglutesscallback, &
   f9x15iaglutesscallback, &
   f9x16iaglutesscallback, &
   f9x11raglutesscallback, &
   f9x12raglutesscallback, &
   f9x13raglutesscallback, &
   f9x14raglutesscallback, &
   f9x15raglutesscallback, &
   f9x16raglutesscallback
#endif
#else
#ifdef OLDNULLFUNC
   f9x1glutesscallback, &
#endif
   f9x0glutesscallback
#endif
private :: &
#ifdef GLU_1_1_TESS
#else
   f9xglutessendcontour, &
   f9xglutessendpolygon, &
   f9xglutessnormal, &
   f9xglutessproperty, &
#endif
   f9x0glutessvertex, &
   f9x1glutessvertex, &
   f9x4glutessvertex, &
   f9x5glutessvertex, &
   f9x6glutessvertex, &
   f9x7glutessvertex, &
   f9xgluunproject

!------
! Types
!------

type, public :: gluquadricobj
   type(glcptr), dimension(2) :: addr
end type gluquadricobj

type, public :: glunurbsobj
   type(glcptr), dimension(2) :: addr
end type glunurbsobj

integer(kind=glint), parameter, private :: num_tess_cb = 12
type, public :: glutesselator
   type(glcptr), dimension(num_tess_cb+1) :: addr
end type glutesselator

#ifdef FULL_TESSCALLBACK

type, public :: f90glu_tess_begin_type
   integer :: val
end type f90glu_tess_begin_type

type, public :: f90glu_tess_edgeflag_type
   integer :: val
end type f90glu_tess_edgeflag_type

type, public :: f90glu_tess_vertex_type
   integer :: val
end type f90glu_tess_vertex_type

type, public :: f90glu_tess_vertex_int_type
   integer :: val
end type f90glu_tess_vertex_int_type

type, public :: f90glu_tess_end_type
   integer :: val
end type f90glu_tess_end_type

type, public :: f90glu_tess_error_type
   integer :: val
end type f90glu_tess_error_type

type, public :: f90glu_tess_combine_type
   integer :: val
end type f90glu_tess_combine_type

type, public :: f90glu_tess_combine_int_type
   integer :: val
end type f90glu_tess_combine_int_type

type, public :: f90glu_tess_begindata_type
   integer :: val
end type f90glu_tess_begindata_type

type, public :: f90glu_tess_edgeflagdata_type
   integer :: val
end type f90glu_tess_edgeflagdata_type

type, public :: f90glu_tess_vertexdata_type
   integer :: val
end type f90glu_tess_vertexdata_type

type, public :: f90glu_tess_enddata_type
   integer :: val
end type f90glu_tess_enddata_type

type, public :: f90glu_tess_errordata_type
   integer :: val
end type f90glu_tess_errordata_type

type, public :: f90glu_tess_combinedata_type
   integer :: val
end type f90glu_tess_combinedata_type

type, public :: f90glu_tess_begindata_int_type
   integer :: val
end type f90glu_tess_begindata_int_type

type, public :: f90glu_tess_efdata_int_type
   integer :: val
end type f90glu_tess_efdata_int_type

type, public :: f90glu_tess_vertexdata_int_type
   integer :: val
end type f90glu_tess_vertexdata_int_type

type, public :: f90glu_tess_enddata_int_type
   integer :: val
end type f90glu_tess_enddata_int_type

type, public :: f90glu_tess_errordata_int_type
   integer :: val
end type f90glu_tess_errordata_int_type

type, public :: f90glu_tess_combdata_int_type
   integer :: val
end type f90glu_tess_combdata_int_type

#endif

!-------------------
! symbolic constants
!-------------------

integer(kind=glenum), parameter, public :: GLU_FALSE = 0
integer(kind=glenum), parameter, public :: GLU_TRUE = 1
integer(kind=glenum), parameter, public :: GLU_VERSION_1_1 = 1
integer(kind=glenum), parameter, public :: GLU_VERSION_1_2 = 1
integer(kind=glenum), parameter, public :: GLU_VERSION = 100800
integer(kind=glenum), parameter, public :: GLU_EXTENSIONS = 100801
integer(kind=glenum), parameter, public :: GLU_INVALID_ENUM = 100900
integer(kind=glenum), parameter, public :: GLU_INVALID_VALUE = 100901
integer(kind=glenum), parameter, public :: GLU_OUT_OF_MEMORY = 100902
integer(kind=glenum), parameter, public :: GLU_INCOMPATIBLE_GL_VERSION = 100903
integer(kind=glenum), parameter, public :: GLU_OUTLINE_POLYGON = 100240
integer(kind=glenum), parameter, public :: GLU_OUTLINE_PATCH = 100241
integer(kind=glenum), parameter, public :: GLU_ERROR = 100103
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR1 = 100251
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR2 = 100252
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR3 = 100253
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR4 = 100254
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR5 = 100255
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR6 = 100256
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR7 = 100257
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR8 = 100258
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR9 = 100259
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR10 = 100260
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR11 = 100261
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR12 = 100262
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR13 = 100263
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR14 = 100264
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR15 = 100265
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR16 = 100266
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR17 = 100267
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR18 = 100268
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR19 = 100269
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR20 = 100270
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR21 = 100271
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR22 = 100272
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR23 = 100273
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR24 = 100274
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR25 = 100275
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR26 = 100276
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR27 = 100277
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR28 = 100278
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR29 = 100279
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR30 = 100280
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR31 = 100281
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR32 = 100282
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR33 = 100283
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR34 = 100284
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR35 = 100285
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR36 = 100286
integer(kind=glenum), parameter, public :: GLU_NURBS_ERROR37 = 100287
integer(kind=glenum), parameter, public :: GLU_AUTO_LOAD_MATRIX = 100200
integer(kind=glenum), parameter, public :: GLU_CULLING = 100201
integer(kind=glenum), parameter, public :: GLU_SAMPLING_TOLERANCE = 100203
integer(kind=glenum), parameter, public :: GLU_DISPLAY_MODE = 100204
integer(kind=glenum), parameter, public :: GLU_PARAMETRIC_TOLERANCE = 100202
integer(kind=glenum), parameter, public :: GLU_SAMPLING_METHOD = 100205
integer(kind=glenum), parameter, public :: GLU_U_STEP = 100206
integer(kind=glenum), parameter, public :: GLU_V_STEP = 100207
integer(kind=glenum), parameter, public :: GLU_PATH_LENGTH = 100215
integer(kind=glenum), parameter, public :: GLU_PARAMETRIC_ERROR = 100216
integer(kind=glenum), parameter, public :: GLU_DOMAIN_DISTANCE = 100217
integer(kind=glenum), parameter, public :: GLU_MAP1_TRIM_2 = 100210
integer(kind=glenum), parameter, public :: GLU_MAP1_TRIM_3 = 100211
integer(kind=glenum), parameter, public :: GLU_POINT = 100010
integer(kind=glenum), parameter, public :: GLU_LINE = 100011
integer(kind=glenum), parameter, public :: GLU_FILL = 100012
integer(kind=glenum), parameter, public :: GLU_SILHOUETTE = 100013
integer(kind=glenum), parameter, public :: GLU_SMOOTH = 100000
integer(kind=glenum), parameter, public :: GLU_FLAT = 100001
integer(kind=glenum), parameter, public :: GLU_NONE = 100002
integer(kind=glenum), parameter, public :: GLU_OUTSIDE = 100020
integer(kind=glenum), parameter, public :: GLU_INSIDE = 100021
integer(kind=glenum), parameter, public :: GLU_BEGIN = 100100
integer(kind=glenum), parameter, public :: GLU_VERTEX = 100101
integer(kind=glenum), parameter, public :: GLU_END = 100102
integer(kind=glenum), parameter, public :: GLU_EDGE_FLAG = 100104
#ifdef FULL_TESSCALLBACK
type(f90glu_tess_begin_type), parameter, public :: GLU_TESS_BEGIN = f90glu_tess_begin_type(100100)
type(f90glu_tess_vertex_type), parameter, public :: GLU_TESS_VERTEX = f90glu_tess_vertex_type(100101)
type(f90glu_tess_vertex_int_type), parameter, public :: GLU_TESS_VERTEX_INT = f90glu_tess_vertex_int_type(100101)
type(f90glu_tess_end_type), parameter, public :: GLU_TESS_END = f90glu_tess_end_type(100102)
type(f90glu_tess_error_type), parameter, public :: GLU_TESS_ERROR = f90glu_tess_error_type(100103)
type(f90glu_tess_edgeflag_type), parameter, public :: GLU_TESS_EDGE_FLAG = f90glu_tess_edgeflag_type(100104)
type(f90glu_tess_combine_type), parameter, public :: GLU_TESS_COMBINE = f90glu_tess_combine_type(100105)
type(f90glu_tess_combine_int_type), parameter, public :: GLU_TESS_COMBINE_INT = f90glu_tess_combine_int_type(100105)
type(f90glu_tess_begindata_type), parameter, public :: GLU_TESS_BEGIN_DATA = f90glu_tess_begindata_type(100106)
type(f90glu_tess_vertexdata_type), parameter, public :: GLU_TESS_VERTEX_DATA = f90glu_tess_vertexdata_type(100107)
type(f90glu_tess_enddata_type), parameter, public :: GLU_TESS_END_DATA = f90glu_tess_enddata_type(100108)
type(f90glu_tess_errordata_type), parameter, public :: GLU_TESS_ERROR_DATA = f90glu_tess_errordata_type(100109)
type(f90glu_tess_edgeflagdata_type), parameter, public :: GLU_TESS_EDGE_FLAG_DATA = f90glu_tess_edgeflagdata_type(100110)
type(f90glu_tess_combinedata_type), parameter, public :: GLU_TESS_COMBINE_DATA = f90glu_tess_combinedata_type(100111)
type(f90glu_tess_begindata_int_type), parameter, public :: GLU_TESS_BEGIN_DATA_INT = f90glu_tess_begindata_int_type(100106)
type(f90glu_tess_vertexdata_int_type), parameter, public :: GLU_TESS_VERTEX_DATA_INT = f90glu_tess_vertexdata_int_type(100107)
type(f90glu_tess_enddata_int_type), parameter, public :: GLU_TESS_END_DATA_INT = f90glu_tess_enddata_int_type(100108)
type(f90glu_tess_errordata_int_type), parameter, public :: GLU_TESS_ERROR_DATA_INT = f90glu_tess_errordata_int_type(100109)
type(f90glu_tess_efdata_int_type), parameter, public :: GLU_TESS_EDGE_FLAG_DATA_INT = f90glu_tess_efdata_int_type(100110)
type(f90glu_tess_combdata_int_type), parameter, public :: GLU_TESS_COMBINE_DATA_INT = f90glu_tess_combdata_int_type(100111)
#else
integer(kind=glenum), parameter, public :: GLU_TESS_BEGIN = 100100
integer(kind=glenum), parameter, public :: GLU_TESS_VERTEX = 100101
integer(kind=glenum), parameter, public :: GLU_TESS_END = 100102
integer(kind=glenum), parameter, public :: GLU_TESS_ERROR = 100103
integer(kind=glenum), parameter, public :: GLU_TESS_EDGE_FLAG = 100104
integer(kind=glenum), parameter, public :: GLU_TESS_COMBINE = 100105
integer(kind=glenum), parameter, public :: GLU_TESS_BEGIN_DATA = 100106
integer(kind=glenum), parameter, public :: GLU_TESS_VERTEX_DATA = 100107
integer(kind=glenum), parameter, public :: GLU_TESS_END_DATA = 100108
integer(kind=glenum), parameter, public :: GLU_TESS_ERROR_DATA = 100109
integer(kind=glenum), parameter, public :: GLU_TESS_EDGE_FLAG_DATA = 100110
integer(kind=glenum), parameter, public :: GLU_TESS_COMBINE_DATA = 100111
#endif
integer(kind=glenum), parameter, public :: GLU_CW = 100120
integer(kind=glenum), parameter, public :: GLU_CCW = 100121
integer(kind=glenum), parameter, public :: GLU_INTERIOR = 100122
integer(kind=glenum), parameter, public :: GLU_EXTERIOR = 100123
integer(kind=glenum), parameter, public :: GLU_UNKNOWN = 100124
integer(kind=glenum), parameter, public :: GLU_TESS_WINDING_RULE = 100140
integer(kind=glenum), parameter, public :: GLU_TESS_BOUNDARY_ONLY = 100141
integer(kind=glenum), parameter, public :: GLU_TESS_TOLERANCE = 100142
integer(kind=glenum), parameter, public :: GLU_TESS_ERROR1 = 100151
integer(kind=glenum), parameter, public :: GLU_TESS_ERROR2 = 100152
integer(kind=glenum), parameter, public :: GLU_TESS_ERROR3 = 100153
integer(kind=glenum), parameter, public :: GLU_TESS_ERROR4 = 100154
integer(kind=glenum), parameter, public :: GLU_TESS_ERROR5 = 100155
integer(kind=glenum), parameter, public :: GLU_TESS_ERROR6 = 100156
integer(kind=glenum), parameter, public :: GLU_TESS_ERROR7 = 100157
integer(kind=glenum), parameter, public :: GLU_TESS_ERROR8 = 100158
integer(kind=glenum), parameter, public :: GLU_TESS_ERROR9 = 100159
integer(kind=glenum), parameter, public :: GLU_TESS_MISSING_BEGIN_POLYGON=100151
integer(kind=glenum), parameter, public :: GLU_TESS_MISSING_BEGIN_CONTOUR=100152
integer(kind=glenum), parameter, public :: GLU_TESS_MISSING_END_POLYGON = 100153
integer(kind=glenum), parameter, public :: GLU_TESS_MISSING_END_CONTOUR = 100154
integer(kind=glenum), parameter, public :: GLU_TESS_COORD_TOO_LARGE = 100155
integer(kind=glenum), parameter, public :: GLU_TESS_NEED_COMBINE_CALLBACK=100156
integer(kind=glenum), parameter, public :: GLU_TESS_WINDING_ODD = 100130
integer(kind=glenum), parameter, public :: GLU_TESS_WINDING_NONZERO = 100131
integer(kind=glenum), parameter, public :: GLU_TESS_WINDING_POSITIVE = 100132
integer(kind=glenum), parameter, public :: GLU_TESS_WINDING_NEGATIVE = 100133
integer(kind=glenum), parameter, public :: GLU_TESS_WINDING_ABS_GEQ_TWO = 100134
#ifdef OLDNULLFUNC
integer(kind=glint), parameter, public :: GLUNULLFUNC = 1
#endif

! interface blocks

#ifdef OLDNULLFUNC
#else
integer(kind=glint), external :: myglucomparefuncs
#ifdef SALFORD
STDCALL myglucomparefuncs 'myglucomparefuncs' (REF)
#endif
#endif

interface glubegincurve
   module procedure f9xglubegincurve
end interface

interface glubeginpolygon
   module procedure f9xglubeginpolygon
end interface

interface glubeginsurface
   module procedure f9xglubeginsurface
end interface

interface glubegintrim
   module procedure f9xglubegintrim
end interface

interface glubuild1dmipmaps
   module procedure f9x0glubuild1dmipmaps
   module procedure f9x1glubuild1dmipmaps
   module procedure f9x4glubuild1dmipmaps
#ifndef NOBYTE
   module procedure f9x5glubuild1dmipmaps
#endif
#ifndef NOSHORT
   module procedure f9x6glubuild1dmipmaps
#endif
end interface

interface glubuild2dmipmaps
   module procedure f9x0glubuild2dmipmaps
   module procedure f9x1glubuild2dmipmaps
   module procedure f9x4glubuild2dmipmaps
#ifndef NOBYTE
   module procedure f9x5glubuild2dmipmaps
#endif
#ifndef NOSHORT
   module procedure f9x6glubuild2dmipmaps
#endif
end interface

interface glucylinder
   module procedure f9xglucylinder
end interface

interface gludeletenurbsrenderer
   module procedure f9xgludeletenurbsrenderer
end interface

interface gludeletequadric
   module procedure f9xgludeletequadric
end interface

interface gludeletetess
   module procedure f9xgludeletetess
end interface

interface gludisk
   module procedure f9xgludisk
end interface

interface gluendcurve
   module procedure f9xgluendcurve
end interface

interface gluendpolygon
   module procedure f9xgluendpolygon
end interface

interface gluendsurface
   module procedure f9xgluendsurface
end interface

interface gluendtrim
   module procedure f9xgluendtrim
end interface

interface gluerrorstring
   module procedure f9xgluerrorstring
end interface

interface glugetnurbsproperty
   module procedure f9xglugetnurbsproperty
end interface

interface glugetstring
   module procedure f9xglugetstring
end interface

#ifdef GLU_1_1_TESS
#else
interface glugettessproperty
   module procedure f9xglugettessproperty
end interface
#endif

interface gluloadsamplingmatrices
   module procedure f9xgluloadsamplingmatrices
end interface

interface glulookat
   module procedure f9xglulookat
end interface

interface glunewnurbsrenderer
   module procedure f9xglunewnurbsrenderer
end interface

interface glunewquadric
   module procedure f9xglunewquadric
end interface

interface glunewtess
   module procedure f9xglunewtess
end interface

interface glunextcontour
   module procedure f9xglunextcontour
end interface

interface glunurbscallback
   module procedure f9xglunurbscallback
#ifdef OLDNULLFUNC
   module procedure f9x1glunurbscallback
#endif
end interface

interface glunurbscurve
   module procedure f9xglunurbscurve
end interface

interface glunurbsproperty
   module procedure f9xglunurbsproperty
end interface

interface glunurbssurface
   module procedure f9xglunurbssurface
end interface

interface gluortho2d
   module procedure f9xgluortho2d
end interface

interface glupartialdisk
   module procedure f9xglupartialdisk
end interface

interface gluperspective
   module procedure f9xgluperspective
end interface

interface glupickmatrix
   module procedure f9xglupickmatrix
end interface

interface gluproject
   module procedure f9xgluproject
end interface

interface glupwlcurve
   module procedure f9xglupwlcurve
end interface

interface gluquadriccallback
   module procedure f9xgluquadriccallback
#ifdef OLDNULLFUNC
   module procedure f9x1gluquadriccallback
#endif
end interface

interface gluquadricdrawstyle
   module procedure f9xgluquadricdrawstyle
end interface

interface gluquadricnormals
   module procedure f9xgluquadricnormals
end interface

interface gluquadricorientation
   module procedure f9xgluquadricorientation
end interface

interface gluquadrictexture
   module procedure f9xgluquadrictexture
end interface

interface gluscaleimage
   module procedure f9x0gluscaleimage
   module procedure f9x1gluscaleimage
   module procedure f9x4gluscaleimage
#ifndef NOBYTE
   module procedure f9x5gluscaleimage
#endif
#ifndef NOSHORT
   module procedure f9x6gluscaleimage
#endif
   module procedure f9x01gluscaleimage
   module procedure f9x04gluscaleimage
#ifndef NOBYTE
   module procedure f9x05gluscaleimage
#endif
#ifndef NOSHORT
   module procedure f9x06gluscaleimage
#endif
   module procedure f9x10gluscaleimage
   module procedure f9x14gluscaleimage
#ifndef NOBYTE
   module procedure f9x15gluscaleimage
#endif
#ifndef NOSHORT
   module procedure f9x16gluscaleimage
#endif
   module procedure f9x40gluscaleimage
   module procedure f9x41gluscaleimage
#ifndef NOBYTE
   module procedure f9x45gluscaleimage
#endif
#ifndef NOSHORT
   module procedure f9x46gluscaleimage
#endif
#ifndef NOBYTE
   module procedure f9x50gluscaleimage
   module procedure f9x51gluscaleimage
   module procedure f9x54gluscaleimage
#ifndef NOSHORT
   module procedure f9x56gluscaleimage
#endif
#endif
#ifndef NOSHORT
   module procedure f9x60gluscaleimage
   module procedure f9x61gluscaleimage
   module procedure f9x64gluscaleimage
#ifndef NOBYTE
   module procedure f9x65gluscaleimage
#endif
#endif
end interface

interface glusphere
   module procedure f9xglusphere
end interface

#ifdef GLU_1_1_TESS
#else
interface glutessbegincontour
   module procedure f9xglutessbegincontour
end interface

interface glutessbeginpolygon
   module procedure f9x0glutessbeginpolygon
   module procedure f9x1glutessbeginpolygon
end interface
#endif

interface glutesscallback
#ifdef FULL_TESSCALLBACK
   module procedure f9x01glutesscallback
   module procedure f9x02glutesscallback
   module procedure f9x03iglutesscallback
   module procedure f9x03rglutesscallback
   module procedure f9x04glutesscallback
   module procedure f9x05glutesscallback
   module procedure f9x06iglutesscallback
   module procedure f9x06rglutesscallback
   module procedure f9x11iglutesscallback
   module procedure f9x12iglutesscallback
   module procedure f9x13iglutesscallback
   module procedure f9x14iglutesscallback
   module procedure f9x15iglutesscallback
   module procedure f9x16iglutesscallback
   module procedure f9x11rglutesscallback
   module procedure f9x12rglutesscallback
   module procedure f9x13rglutesscallback
   module procedure f9x14rglutesscallback
   module procedure f9x15rglutesscallback
   module procedure f9x16rglutesscallback
#ifdef OLDNULLFUNC
   module procedure f9x01aglutesscallback
   module procedure f9x02aglutesscallback
   module procedure f9x03iaglutesscallback
   module procedure f9x03raglutesscallback
   module procedure f9x04aglutesscallback
   module procedure f9x05aglutesscallback
   module procedure f9x06iaglutesscallback
   module procedure f9x06raglutesscallback
   module procedure f9x11iaglutesscallback
   module procedure f9x12iaglutesscallback
   module procedure f9x13iaglutesscallback
   module procedure f9x14iaglutesscallback
   module procedure f9x15iaglutesscallback
   module procedure f9x16iaglutesscallback
   module procedure f9x11raglutesscallback
   module procedure f9x12raglutesscallback
   module procedure f9x13raglutesscallback
   module procedure f9x14raglutesscallback
   module procedure f9x15raglutesscallback
   module procedure f9x16raglutesscallback
#endif
#else
   module procedure f9x0glutesscallback
#ifdef OLDNULLFUNC
   module procedure f9x1glutesscallback
#endif
#endif
end interface

#ifdef GLU_1_1_TESS
#else
interface glutessendcontour
   module procedure f9xglutessendcontour
end interface

interface glutessendpolygon
   module procedure f9xglutessendpolygon
end interface

interface glutessnormal
   module procedure f9xglutessnormal
end interface

interface glutessproperty
   module procedure f9xglutessproperty
end interface
#endif

interface glutessvertex
   module procedure f9x0glutessvertex
   module procedure f9x1glutessvertex
   module procedure f9x4glutessvertex
#ifndef NOBYTE
   module procedure f9x5glutessvertex
#endif
#ifndef NOSHORT
   module procedure f9x6glutessvertex
#endif
   module procedure f9x7glutessvertex
end interface

interface gluunproject
   module procedure f9xgluunproject
end interface

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

#ifndef OLDNULLFUNC
subroutine glunullfunc()
! this is the GLUNULLFUNC symbol
end subroutine glunullfunc
#endif

subroutine f9xglubegincurve(nurb)
type(glunurbsobj), pointer :: nurb
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglubegincurve(temp,i)
return
end subroutine f9xglubegincurve

subroutine f9xglubeginpolygon(tess)
type(glutesselator), pointer :: tess
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call fglubeginpolygon(temp,i)
return
end subroutine f9xglubeginpolygon

subroutine f9xglubeginsurface(nurb)
type(glunurbsobj), pointer :: nurb
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglubeginsurface(temp,i)
return
end subroutine f9xglubeginsurface

subroutine f9xglubegintrim(nurb)
type(glunurbsobj), pointer :: nurb
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglubegintrim(temp,i)
return
end subroutine f9xglubegintrim

function f9x0glubuild1dmipmaps(xtarget,component,width,format,xtype, &
   data) result(resf9x0glubuild1dmipmaps)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN data
integer(kind=glint) :: resf9x0glubuild1dmipmaps
resf9x0glubuild1dmipmaps = f9y0glubuild1dmipmaps(xtarget,component,width, &
   format,xtype,data,size_glint(data))
return
end function f9x0glubuild1dmipmaps

function f9x1glubuild1dmipmaps(xtarget,component,width,format,xtype, &
   data) result(resf9x1glubuild1dmipmaps)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) :: resf9x1glubuild1dmipmaps
resf9x1glubuild1dmipmaps = f9y1glubuild1dmipmaps(xtarget,component,width, &
   format,xtype,data)
return
end function f9x1glubuild1dmipmaps

function f9x4glubuild1dmipmaps(xtarget,component,width,format,xtype, &
   data) result(resf9x4glubuild1dmipmaps)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN data
integer(kind=glint) :: resf9x4glubuild1dmipmaps
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(data%addr(i:i))
end do
i = glcptr_length
resf9x4glubuild1dmipmaps = f9y4glubuild1dmipmaps(xtarget,component,width, &
   format,xtype,temp,i)
return
end function f9x4glubuild1dmipmaps

function f9x5glubuild1dmipmaps(xtarget,component,width,format,xtype, &
   data) result(resf9x5glubuild1dmipmaps)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) :: resf9x5glubuild1dmipmaps
resf9x5glubuild1dmipmaps = f9y5glubuild1dmipmaps(xtarget,component,width, &
   format,xtype,data)
return
end function f9x5glubuild1dmipmaps

function f9x6glubuild1dmipmaps(xtarget,component,width,format,xtype, &
   data) result(resf9x6glubuild1dmipmaps)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN data
integer(kind=glint) :: resf9x6glubuild1dmipmaps
resf9x6glubuild1dmipmaps = f9y6glubuild1dmipmaps(xtarget,component,width, &
   format,xtype,data,size_glint(data))
return
end function f9x6glubuild1dmipmaps

function f9x0glubuild2dmipmaps(xtarget,component,width,height, &
   format,xtype,data) result(resf9x0glubuild2dmipmaps)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN data
integer(kind=glint) :: resf9x0glubuild2dmipmaps
resf9x0glubuild2dmipmaps = f9y0glubuild2dmipmaps(xtarget,component,width, &
   height,format,xtype,data,size_glint(data))
return
end function f9x0glubuild2dmipmaps

function f9x1glubuild2dmipmaps(xtarget,component,width,height, &
   format,xtype,data) result(resf9x1glubuild2dmipmaps)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) :: resf9x1glubuild2dmipmaps
resf9x1glubuild2dmipmaps = f9y1glubuild2dmipmaps(xtarget,component,width, &
   height,format,xtype,data)
return
end function f9x1glubuild2dmipmaps

function f9x4glubuild2dmipmaps(xtarget,component,width,height, &
   format,xtype,data) result(resf9x4glubuild2dmipmaps)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
type(glcptr) INTENT_IN data
integer(kind=glint) :: resf9x4glubuild2dmipmaps
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(data%addr(i:i))
end do
i = glcptr_length
resf9x4glubuild2dmipmaps = f9y4glubuild2dmipmaps(xtarget,component,width, &
   height,format,xtype,temp,i)
return
end function f9x4glubuild2dmipmaps

function f9x5glubuild2dmipmaps(xtarget,component,width,height, &
   format,xtype,data) result(resf9x5glubuild2dmipmaps)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN data
integer(kind=glint) :: resf9x5glubuild2dmipmaps
resf9x5glubuild2dmipmaps = f9y5glubuild2dmipmaps(xtarget,component,width, &
   height,format,xtype,data)
return
end function f9x5glubuild2dmipmaps

function f9x6glubuild2dmipmaps(xtarget,component,width,height, &
   format,xtype,data) result(resf9x6glubuild2dmipmaps)
integer(kind=glenum) INTENT_IN xtarget
integer(kind=glint) INTENT_IN component
integer(kind=glsizei) INTENT_IN width, height
integer(kind=glenum) INTENT_IN format, xtype
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN data
integer(kind=glint) :: resf9x6glubuild2dmipmaps
resf9x6glubuild2dmipmaps = f9y6glubuild2dmipmaps(xtarget,component,width, &
   height,format,xtype,data,size_glint(data))
return
end function f9x6glubuild2dmipmaps

subroutine f9xglucylinder(quad,base,top,height,slices,stacks)
type(gluquadricobj), pointer :: quad
real(kind=gldouble) INTENT_IN base, top, height
integer(kind=glint) INTENT_IN slices, stacks
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglucylinder(temp,i,base,top,height,slices,stacks)
return
end subroutine f9xglucylinder

subroutine f9xgludeletenurbsrenderer(nurb)
type(glunurbsobj), pointer :: nurb
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fgludeletenurbsrenderer(temp,i)
deallocate(nurb)
return
end subroutine f9xgludeletenurbsrenderer

subroutine f9xgludeletequadric(quad)
type(gluquadricobj), pointer :: quad
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
call fgludeletequadric(temp,i)
deallocate(quad)
return
end subroutine f9xgludeletequadric

subroutine f9xgludeletetess(tess)
type(glutesselator), pointer :: tess
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call fgludeletetess(temp,i)
deallocate(tess)
return
end subroutine f9xgludeletetess

subroutine f9xgludisk(quad,inner,outer,slices,loops)
type(gluquadricobj), pointer :: quad
real(kind=gldouble) INTENT_IN inner, outer
integer(kind=glint) INTENT_IN slices, loops
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
call fgludisk(temp,i,inner,outer,slices,loops)
return
end subroutine f9xgludisk

subroutine f9xgluendcurve(nurb)
type(glunurbsobj), pointer :: nurb
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fgluendcurve(temp,i)
return
end subroutine f9xgluendcurve

subroutine f9xgluendpolygon(tess)
type(glutesselator), pointer :: tess
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call fgluendpolygon(temp,i)
return
end subroutine f9xgluendpolygon

subroutine f9xgluendsurface(nurb)
type(glunurbsobj), pointer :: nurb
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fgluendsurface(temp,i)
return
end subroutine f9xgluendsurface

subroutine f9xgluendtrim(nurb)
type(glunurbsobj), pointer :: nurb
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fgluendtrim(temp,i)
return
end subroutine f9xgluendtrim

function f9xgluerrorstring(error) result(resf9xgluerrorstring)
integer(kind=glenum) INTENT_IN error
character(len=1), pointer, dimension(:) :: resf9xgluerrorstring
integer(kind=gluint), allocatable, dimension(:) :: temp
integer(kind=glint) :: length, pos
#ifdef PURE_FUNCTIONS
integer(kind=glint) :: dummy
#endif
length = f9ygluerrorstringlen(error)
allocate(temp(length))
allocate(resf9xgluerrorstring(length))
#ifdef PURE_FUNCTIONS
do pos=1,length
   temp(pos) = f9yglugetstring1(pos)
end do
#else
call f9yglugetstring(temp)
#endif
do pos=1,length
   resf9xgluerrorstring(pos) = char(temp(pos))
end do
#ifdef PURE_FUNCTIONS
dummy = f9yglugetstringfree()
#else
call f9yglugetstringfree()
#endif
deallocate(temp)
return
end function f9xgluerrorstring

subroutine f9xglugetnurbsproperty(nurb,property,data)
type(glunurbsobj), pointer :: nurb
integer(kind=glenum) INTENT_IN property
real(kind=glfloat), intent(out) :: data
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglugetnurbsproperty(temp,i,property,data)
return
end subroutine f9xglugetnurbsproperty

function f9xglugetstring(name) result(resf9xglugetstring)
integer(kind=glenum) INTENT_IN name
character(len=1), pointer, dimension(:) :: resf9xglugetstring
integer(kind=gluint), allocatable, dimension(:) :: temp
integer(kind=glint) :: length, pos
#ifdef PURE_FUNCTIONS
integer(kind=glint) :: dummy
#endif
length = f9yglugetstringlen(name)
allocate(temp(length))
allocate(resf9xglugetstring(length))
#ifdef PURE_FUNCTIONS
do pos=1,length
   temp(pos) = f9yglugetstring1(pos)
end do
#else
call f9yglugetstring(temp)
#endif
do pos=1,length
   resf9xglugetstring(pos) = char(temp(pos))
end do
#ifdef PURE_FUNCTIONS
dummy = f9yglugetstringfree()
#else
call f9yglugetstringfree()
#endif
deallocate(temp)
return
end function f9xglugetstring

#ifdef GLU_1_1_TESS
#else
subroutine f9xglugettessproperty(tess,which,data)
type(glutesselator), pointer :: tess
integer(kind=glenum) INTENT_IN which
real(kind=gldouble), intent(out) :: data
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call fglugettessproperty(temp,i,which,data)
return
end subroutine f9xglugettessproperty
#endif

subroutine f9xgluloadsamplingmatrices(nurb,model,perspective,view)
type(glunurbsobj), pointer :: nurb
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN model, perspective
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN view
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fgluloadsamplingmatrices(temp,i,model,perspective,view)
return
end subroutine f9xgluloadsamplingmatrices

subroutine f9xglulookat(eyex,eyey,eyez,centerx,centery,centerz,upx,upy,upz)
real(kind=gldouble) INTENT_IN eyex, eyey, eyez, centerx, centery, centerz, &
    upx, upy, upz
call fglulookat(eyex,eyey,eyez,centerx,centery,centerz,upx,upy,upz)
return
end subroutine f9xglulookat

function f9xglunewnurbsrenderer() result(resf9xglunewnurbsrenderer)
type(glunurbsobj), pointer :: resf9xglunewnurbsrenderer
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glcint) :: nbytes
#ifdef F_LANG
integer(kind=glint) :: length,pos
#endif
integer :: i
logical :: isnull
allocate(resf9xglunewnurbsrenderer)
nbytes = glcptr_length
#ifdef F_LANG
pos = f9yglunewnurbsrenderer(temp,nbytes)
#else
call f9yglunewnurbsrenderer(temp,nbytes)
#endif
#ifdef F_LANG
length = glu_f_uint_array_recv_beg()
do pos=1,min(length,size(temp))
   temp(pos) = glu_f_uint_array_recv(pos)
end do
#endif
#ifdef F_LANG
pos = glu_f_uint_array_free()
#endif
do i=1,glcptr_length
   resf9xglunewnurbsrenderer%addr(1)%addr(i:i) = char(temp(i))
   resf9xglunewnurbsrenderer%addr(2)%addr(i:i) = char(temp(i+glcptr_length))
end do
isnull = (resf9xglunewnurbsrenderer%addr(1) == glnullptr)
if (isnull) then
   deallocate(resf9xglunewnurbsrenderer)
#ifdef F_LANG
   resf9xglunewnurbsrenderer => NULL()
#else
   nullify(resf9xglunewnurbsrenderer)
#endif
end if
return
end function f9xglunewnurbsrenderer

function f9xglunewquadric() result(resf9xglunewquadric)
type(gluquadricobj), pointer :: resf9xglunewquadric
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glcint) :: nbytes
#ifdef F_LANG
integer(kind=glint) :: length,pos
#endif
integer :: i
logical :: isnull
allocate(resf9xglunewquadric)
nbytes = glcptr_length
#ifdef F_LANG
pos = f9yglunewquadric(temp,nbytes)
#else
call f9yglunewquadric(temp,nbytes)
#endif
#ifdef F_LANG
length = glu_f_uint_array_recv_beg()
do pos=1,min(length,size(temp))
   temp(pos) = glu_f_uint_array_recv(pos)
end do
#endif
#ifdef F_LANG
pos = glu_f_uint_array_free()
#endif
do i=1,glcptr_length
   resf9xglunewquadric%addr(1)%addr(i:i) = char(temp(i))
   resf9xglunewquadric%addr(2)%addr(i:i) = char(temp(i+glcptr_length))
end do
isnull = (resf9xglunewquadric%addr(1) == glnullptr)
if (isnull) then
   deallocate(resf9xglunewquadric)
#ifdef F_LANG
   resf9xglunewquadric => NULL()
#else
   nullify(resf9xglunewquadric)
#endif
end if
return
end function f9xglunewquadric

function f9xglunewtess() result(resf9xglunewtess)
type(glutesselator), pointer :: resf9xglunewtess
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glcint) :: nbytes
#ifdef F_LANG
integer(kind=glint) :: length,pos
#endif
integer :: i,j
logical :: isnull
allocate(resf9xglunewtess)
nbytes = glcptr_length
#ifdef F_LANG
pos = f9yglunewtess(temp,nbytes)
#else
call f9yglunewtess(temp,nbytes)
#endif
#ifdef F_LANG
length = glu_f_uint_array_recv_beg()
do pos=1,min(length,size(temp))
   temp(pos) = glu_f_uint_array_recv(pos)
end do
#endif
#ifdef F_LANG
pos = glu_f_uint_array_free()
#endif
do i=1,glcptr_length
   resf9xglunewtess%addr(1)%addr(i:i) = char(temp(i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      resf9xglunewtess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
isnull = (resf9xglunewtess%addr(1) == glnullptr)
if (isnull) then
   deallocate(resf9xglunewtess)
#ifdef F_LANG
   resf9xglunewtess => NULL()
#else
   nullify(resf9xglunewtess)
#endif
end if
return
end function f9xglunewtess

subroutine f9xglunextcontour(tess,xtype)
type(glutesselator), pointer :: tess
integer(kind=glenum) INTENT_IN xtype
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call fglunextcontour(temp,i,xtype)
return
end subroutine f9xglunextcontour

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glunurbscallback(nurb,which)
#else
subroutine f9x1glunurbscallback(nurb,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glunurbsobj), pointer :: nurb
integer(kind=glenum) INTENT_IN which
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then
#endif
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call f9y1glunurbscallback(temp,i,which)
do i=1,glcptr_length
   nurb%addr(2)%addr(i:i) = char(temp(i+glcptr_length))
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x1glunurbscallback

subroutine f9xglunurbscallback(nurb,which,callbackfunc)
type(glunurbsobj), pointer :: nurb
integer(kind=glenum) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(errno)
   use opengl_kinds
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglunurbscallback(temp,i,which,callbackfunc)
do i=1,glcptr_length
   nurb%addr(2)%addr(i:i) = char(temp(i+glcptr_length))
end do
return
end subroutine f9xglunurbscallback

#else

subroutine f9xglunurbscallback(nurb,which,callbackfunc)
type(glunurbsobj), pointer :: nurb
integer(kind=glenum) INTENT_IN which
#ifdef EXTERN_CALLBACK
external callbackfunc
#else
interface
   C_BINDING subroutine callbackfunc(errno)
   use opengl_kinds
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   end subroutine callbackfunc
end interface
#endif
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glunurbscallback(temp,i,which)
else
   call fglunurbscallback(temp,i,which,callbackfunc)
endif
do i=1,glcptr_length
   nurb%addr(2)%addr(i:i) = char(temp(i+glcptr_length))
end do
return
end subroutine f9xglunurbscallback

#endif

subroutine f9xglunurbscurve(nurb,knotcount,knots,stride,control,order,xtype)
type(glunurbsobj), pointer :: nurb
integer(kind=glint) INTENT_IN knotcount
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: knots
integer(kind=glint) INTENT_IN stride
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: control
integer(kind=glint) INTENT_IN order
integer(kind=glenum) INTENT_IN xtype
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglunurbscurve(temp,i,knotcount,knots,stride,control,order,xtype)
return
end subroutine f9xglunurbscurve

subroutine f9xglunurbsproperty(nurb,property,value)
type(glunurbsobj), pointer :: nurb
integer(kind=glenum) INTENT_IN property
real(kind=glfloat) INTENT_IN value
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglunurbsproperty(temp,i,property,value)
return
end subroutine f9xglunurbsproperty

subroutine f9xglunurbssurface(nurb,sknotcount,sknots,tknotcount,tknots, &
   sstride,tstride,control,sorder,torder,xtype)
type(glunurbsobj), pointer :: nurb
integer(kind=glint) INTENT_IN sknotcount
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: sknots
integer(kind=glint) INTENT_IN tknotcount
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: tknots
integer(kind=glint) INTENT_IN sstride, tstride
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: control
integer(kind=glint) INTENT_IN sorder, torder
integer(kind=glenum) INTENT_IN xtype
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglunurbssurface(temp,i,sknotcount,sknots,tknotcount,tknots,sstride, &
   tstride,control,sorder,torder,xtype)
return
end subroutine f9xglunurbssurface

subroutine f9xgluortho2d(left,right,bottom,top)
real(kind=gldouble) INTENT_IN left, right, bottom, top
call fgluortho2d(left,right,bottom,top)
return
end subroutine f9xgluortho2d

subroutine f9xglupartialdisk(quad,inner,outer,slices,loops,start,sweep)
type(gluquadricobj), pointer :: quad
real(kind=gldouble) INTENT_IN inner, outer
integer(kind=glint) INTENT_IN slices, loops
real(kind=gldouble) INTENT_IN start, sweep
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglupartialdisk(temp,i,inner,outer,slices,loops,start,sweep)
return
end subroutine f9xglupartialdisk

subroutine f9xgluperspective(fovy,aspect,znear,zfar)
real(kind=gldouble) INTENT_IN fovy, aspect, znear, zfar
call fgluperspective(fovy,aspect,znear,zfar)
return
end subroutine f9xgluperspective

subroutine f9xglupickmatrix(x,y,delx,dely,viewport)
real(kind=gldouble) INTENT_IN x, y, delx, dely
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN viewport
call fglupickmatrix(x,y,delx,dely,viewport)
return
end subroutine f9xglupickmatrix

#ifdef PURE_FUNCTIONS
subroutine f9xgluproject(objx,objy,objz,model,proj,view,winx,winy,winz,errcode)
#else
function f9xgluproject(objx,objy,objz,model,proj,view,winx,winy, &
   winz) result(resf9xgluproject)
#endif
real(kind=gldouble) INTENT_IN objx, objy, objz
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN model, proj
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN view
real(kind=gldouble), intent(out) :: winx, winy, winz
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9xgluproject
#endif
#ifdef PURE_FUNCTIONS
call fgluproject(objx,objy,objz,model,proj,view,winx,winy,winz,errcode)
#else
resf9xgluproject = fgluproject(objx,objy,objz,model,proj,view,winx,winy,winz)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9xgluproject
#else
end function f9xgluproject
#endif

subroutine f9xglupwlcurve(nurb,xcount,data,stride,xtype)
type(glunurbsobj), pointer :: nurb
integer(kind=glint) INTENT_IN xcount
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: data
integer(kind=glint) INTENT_IN stride
integer(kind=glenum) INTENT_IN xtype
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(nurb%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(nurb%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglupwlcurve(temp,i,xcount,data,stride,xtype)
return
end subroutine f9xglupwlcurve

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1gluquadriccallback(quad,which)
#else
subroutine f9x1gluquadriccallback(quad,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(gluquadricobj), pointer :: quad
integer(kind=glenum) INTENT_IN which
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then
#endif
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
call f9y1gluquadriccallback(temp,i,which)
do i=1,glcptr_length
   quad%addr(2)%addr(i:i) = char(temp(i+glcptr_length))
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x1gluquadriccallback

subroutine f9xgluquadriccallback(quad,which,callbackfunc)
type(gluquadricobj), pointer :: quad
integer(kind=glenum) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(errno)
   use opengl_kinds
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
call fgluquadriccallback(temp,i,which,callbackfunc)
do i=1,glcptr_length
   quad%addr(2)%addr(i:i) = char(temp(i+glcptr_length))
end do
return
end subroutine f9xgluquadriccallback

#else

subroutine f9xgluquadriccallback(quad,which,callbackfunc)
type(gluquadricobj), pointer :: quad
integer(kind=glenum) INTENT_IN which
#ifdef EXTERN_CALLBACK
external callbackfunc
#else
interface
   C_BINDING subroutine callbackfunc(errno)
   use opengl_kinds
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   end subroutine callbackfunc
end interface
#endif
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1gluquadriccallback(temp,i,which)
else
   call fgluquadriccallback(temp,i,which,callbackfunc)
endif
do i=1,glcptr_length
   quad%addr(2)%addr(i:i) = char(temp(i+glcptr_length))
end do
return
end subroutine f9xgluquadriccallback

#endif

subroutine f9xgluquadricdrawstyle(quad,draw)
type(gluquadricobj), pointer :: quad
integer(kind=glenum) INTENT_IN draw
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
call fgluquadricdrawstyle(temp,i,draw)
return
end subroutine f9xgluquadricdrawstyle

subroutine f9xgluquadricnormals(quad,normal)
type(gluquadricobj), pointer :: quad
integer(kind=glenum) INTENT_IN normal
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
call fgluquadricnormals(temp,i,normal)
return
end subroutine f9xgluquadricnormals

subroutine f9xgluquadricorientation(quad,orientation)
type(gluquadricobj), pointer :: quad
integer(kind=glenum) INTENT_IN orientation
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
call fgluquadricorientation(temp,i,orientation)
return
end subroutine f9xgluquadricorientation

subroutine f9xgluquadrictexture(quad,texture)
type(gluquadricobj), pointer :: quad
logical(kind=glboolean) INTENT_IN texture
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: int_texture
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
if (glone_byte_logical) then
   call fgluquadrictexture(temp,i,texture)
else
   if (texture) then
      int_texture = 1
   else
      int_texture = 0
   end if
   call f9egluquadrictexture(temp,i,int_texture)
end if
return
end subroutine f9xgluquadrictexture

#ifdef PURE_FUNCTIONS
subroutine f9x0gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x0gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x0gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
! need assumed shape to determine size of array
integer(kind=glint), dimension(:), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x0gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y0gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout,size_glint(dataout),errcode)
#else
resf9x0gluscaleimage = f9y0gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout,size_glint(dataout))
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x0gluscaleimage
#else
end function f9x0gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x1gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x1gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x1gluscaleimage)
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
integer(kind=glint) :: resf9x1gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y1gluscaleimage(format,win,hin,typein,datain,wout, &
   hout,typeout,dataout,errcode)
#else
resf9x1gluscaleimage = f9y1gluscaleimage(format,win,hin,typein,datain,wout, &
   hout,typeout,dataout)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x1gluscaleimage
#else
end function f9x1gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x4gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x4gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x4gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
type(glcptr) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
type(glcptr) INTENT_IN dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x4gluscaleimage
#endif
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(datain%addr(i:i))
end do
do i=1,glcptr_length
   temp(i+glcptr_length) = ichar(dataout%addr(i:i))
end do
i = glcptr_length
#ifdef PURE_FUNCTIONS
call f9y4gluscaleimage(format,win,hin,typein,wout, &
   hout,typeout,temp,i,errcode)
#else
resf9x4gluscaleimage = f9y4gluscaleimage(format,win,hin,typein,wout, &
   hout,typeout,temp,i)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x4gluscaleimage
#else
end function f9x4gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x5gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x5gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x5gluscaleimage)
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
integer(kind=glint) :: resf9x5gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y5gluscaleimage(format,win,hin,typein,datain,wout, &
   hout,typeout,dataout,errcode)
#else
resf9x5gluscaleimage = f9y5gluscaleimage(format,win,hin,typein,datain,wout, &
   hout,typeout,dataout)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x5gluscaleimage
#else
end function f9x5gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x6gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x6gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x6gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x6gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y6gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout,size_glint(dataout),errcode)
#else
resf9x6gluscaleimage = f9y6gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout,size_glint(dataout))
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x6gluscaleimage
#else
end function f9x6gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x01gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x01gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x01gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x01gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y01gluscaleimage(format,win,hin,typein,datain,size_glint(datain),wout, &
   hout,typeout,dataout,errcode)
#else
resf9x01gluscaleimage = f9y01gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x01gluscaleimage
#else
end function f9x01gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x04gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x04gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x04gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
type(glcptr) INTENT_IN dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x04gluscaleimage
#endif
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
i = glcptr_length
#ifdef PURE_FUNCTIONS
call f9y04gluscaleimage(format,win,hin,typein,datain,size_glint(datain),wout, &
   hout,typeout,temp,i,errcode)
#else
resf9x04gluscaleimage = f9y04gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,temp,i)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x04gluscaleimage
#else
end function f9x04gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x05gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x05gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x05gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x05gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y05gluscaleimage(format,win,hin,typein,datain,size_glint(datain),wout, &
   hout,typeout,dataout,errcode)
#else
resf9x05gluscaleimage = f9y05gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x05gluscaleimage
#else
end function f9x05gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x06gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x06gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x06gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
! need assumed shape to determine size of array
integer(kind=glint), dimension(:) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x06gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y06gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout,size_glint(dataout),errcode)
#else
resf9x06gluscaleimage = f9y06gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout,size_glint(dataout))
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x06gluscaleimage
#else
end function f9x06gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x10gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x10gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x10gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
! need assumed shape to determine size of array
integer(kind=glint), dimension(:), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x10gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y10gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,dataout,size_glint(dataout),errcode)
#else
resf9x10gluscaleimage = f9y10gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,dataout,size_glint(dataout))
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x10gluscaleimage
#else
end function f9x10gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x14gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x14gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x14gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
type(glcptr) INTENT_IN dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x14gluscaleimage
#endif
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
i = glcptr_length
#ifdef PURE_FUNCTIONS
call f9y14gluscaleimage(format,win,hin,typein,datain,wout, &
   hout,typeout,temp,i,errcode)
#else
resf9x14gluscaleimage = f9y14gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,temp,i)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x14gluscaleimage
#else
end function f9x14gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x15gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x15gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x15gluscaleimage)
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
integer(kind=glint) :: resf9x15gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y15gluscaleimage(format,win,hin,typein,datain,wout, &
   hout,typeout,dataout,errcode)
#else
resf9x15gluscaleimage = f9y15gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,dataout)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x15gluscaleimage
#else
end function f9x15gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x16gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x16gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x16gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
real(kind=glfloat), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x16gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y16gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,dataout,size_glint(dataout),errcode)
#else
resf9x16gluscaleimage = f9y16gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,dataout,size_glint(dataout))
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x16gluscaleimage
#else
end function f9x16gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x40gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x40gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x40gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
type(glcptr) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
! need assumed shape to determine size of array
integer(kind=glint), dimension(:), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x40gluscaleimage
#endif
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(datain%addr(i:i))
end do
i = glcptr_length
#ifdef PURE_FUNCTIONS
call f9y40gluscaleimage(format,win,hin,typein,temp,i, &
   wout,hout,typeout,dataout,size_glint(dataout),errcode)
#else
resf9x40gluscaleimage = f9y40gluscaleimage(format,win,hin,typein,temp,i, &
   wout,hout,typeout,dataout,size_glint(dataout))
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x40gluscaleimage
#else
end function f9x40gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x41gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x41gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x41gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
type(glcptr) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x41gluscaleimage
#endif
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(datain%addr(i:i))
end do
#ifdef PURE_FUNCTIONS
call f9y41gluscaleimage(format,win,hin,typein,temp,i,wout, &
   hout,typeout,dataout,errcode)
#else
resf9x41gluscaleimage = f9y41gluscaleimage(format,win,hin,typein,temp,i, &
   wout,hout,typeout,dataout)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x41gluscaleimage
#else
end function f9x41gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x45gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x45gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x45gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
type(glcptr) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x45gluscaleimage
#endif
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(datain%addr(i:i))
end do
#ifdef PURE_FUNCTIONS
call f9y45gluscaleimage(format,win,hin,typein,temp,i,wout, &
   hout,typeout,dataout,errcode)
#else
resf9x45gluscaleimage = f9y45gluscaleimage(format,win,hin,typein,temp,i, &
   wout,hout,typeout,dataout)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x45gluscaleimage
#else
end function f9x45gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x46gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x46gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x46gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
type(glcptr) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x46gluscaleimage
#endif
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(datain%addr(i:i))
end do
#ifdef PURE_FUNCTIONS
call f9y46gluscaleimage(format,win,hin,typein,temp,i, &
   wout,hout,typeout,dataout,size_glint(dataout),errcode)
#else
resf9x46gluscaleimage = f9y46gluscaleimage(format,win,hin,typein,temp,i, &
   wout,hout,typeout,dataout,size_glint(dataout))
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x46gluscaleimage
#else
end function f9x46gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x50gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x50gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x50gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
! need assumed shape to determine size of array
integer(kind=glint), dimension(:), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x50gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y50gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,dataout,size_glint(dataout),errcode)
#else
resf9x50gluscaleimage = f9y50gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,dataout,size_glint(dataout))
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x50gluscaleimage
#else
end function f9x50gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x51gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x51gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x51gluscaleimage)
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
integer(kind=glint) :: resf9x51gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y51gluscaleimage(format,win,hin,typein,datain,wout, &
   hout,typeout,dataout,errcode)
#else
resf9x51gluscaleimage = f9y51gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,dataout)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x51gluscaleimage
#else
end function f9x51gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x54gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x54gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x54gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
type(glcptr) INTENT_IN dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x54gluscaleimage
#endif
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
i = glcptr_length
#ifdef PURE_FUNCTIONS
call f9y54gluscaleimage(format,win,hin,typein,datain,wout, &
   hout,typeout,temp,i,errcode)
#else
resf9x54gluscaleimage = f9y54gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,temp,i)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x54gluscaleimage
#else
end function f9x54gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x56gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x56gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x56gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
integer(kind=glbyte), dimension(DEFER_DIM) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x56gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y56gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,dataout,size_glint(dataout),errcode)
#else
resf9x56gluscaleimage = f9y56gluscaleimage(format,win,hin,typein,datain, &
   wout,hout,typeout,dataout,size_glint(dataout))
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x56gluscaleimage
#else
end function f9x56gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x60gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x60gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x60gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
! need assumed shape to determine size of array
integer(kind=glint), dimension(:), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x60gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y60gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout,errcode)
#else
resf9x60gluscaleimage = f9y60gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x60gluscaleimage
#else
end function f9x60gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x61gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x61gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x61gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
real(kind=glfloat), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x61gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y61gluscaleimage(format,win,hin,typein,datain,size_glint(datain),wout, &
   hout,typeout,dataout,errcode)
#else
resf9x61gluscaleimage = f9y61gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x61gluscaleimage
#else
end function f9x61gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x64gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x64gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x64gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
type(glcptr) INTENT_IN dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x64gluscaleimage
#endif
integer(kind=gluint), dimension(glcptr_length) :: temp
integer(kind=glint) :: i
i = glcptr_length
#ifdef PURE_FUNCTIONS
call f9y64gluscaleimage(format,win,hin,typein,datain,size_glint(datain),wout, &
   hout,typeout,temp,i,errcode)
#else
resf9x64gluscaleimage = f9y64gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,temp,i)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x64gluscaleimage
#else
end function f9x64gluscaleimage
#endif

#ifdef PURE_FUNCTIONS
subroutine f9x65gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout,errcode)
#else
function f9x65gluscaleimage(format,win,hin,typein,datain,wout,hout,typeout, &
   dataout) result(resf9x65gluscaleimage)
#endif
integer(kind=glenum) INTENT_IN format
integer(kind=glsizei) INTENT_IN win, hin
integer(kind=glenum) INTENT_IN typein
! need assumed shape to determine size of array
integer(kind=glshort), dimension(:) INTENT_IN datain
integer(kind=glsizei) INTENT_IN wout, hout
integer(kind=glenum) INTENT_IN typeout
integer(kind=glbyte), dimension(DEFER_DIM), intent(out) :: dataout
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9x65gluscaleimage
#endif
#ifdef PURE_FUNCTIONS
call f9y65gluscaleimage(format,win,hin,typein,datain,size_glint(datain),wout, &
   hout,typeout,dataout,errcode)
#else
resf9x65gluscaleimage = f9y65gluscaleimage(format,win,hin,typein,datain, &
   size_glint(datain),wout,hout,typeout,dataout)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9x65gluscaleimage
#else
end function f9x65gluscaleimage
#endif

subroutine f9xglusphere(quad,radius,slices,stacks)
type(gluquadricobj), pointer :: quad
real(kind=gldouble) INTENT_IN radius
integer(kind=glint) INTENT_IN slices, stacks
integer(kind=gluint), dimension(2*glcptr_length) :: temp
integer(kind=glint) :: i
do i=1,glcptr_length
   temp(i) = ichar(quad%addr(1)%addr(i:i))
   temp(i+glcptr_length) = ichar(quad%addr(2)%addr(i:i))
end do
i = glcptr_length
call fglusphere(temp,i,radius,slices,stacks)
return
end subroutine f9xglusphere

#ifdef GLU_1_1_TESS
#else
subroutine f9xglutessbegincontour(tess)
type(glutesselator), pointer :: tess
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call fglutessbegincontour(temp,i)
return
end subroutine f9xglutessbegincontour

! I'm not sure if this will work, but since polygon_data is just passed on
! to a callback function, with F we'll try just passing it like with any
! other compiler.
subroutine f9x0glutessbeginpolygon(tess,polygon_data)
type(glutesselator), pointer :: tess
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: polygon_data
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y0glutessbeginpolygon(temp,i,polygon_data)
return
end subroutine f9x0glutessbeginpolygon

subroutine f9x1glutessbeginpolygon(tess,polygon_data)
type(glutesselator), pointer :: tess
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: polygon_data
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutessbeginpolygon(temp,i,polygon_data)
return
end subroutine f9x1glutessbeginpolygon
#endif

#ifdef FULL_TESSCALLBACK
#else

subroutine f9x0glutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
integer(kind=glenum) INTENT_IN which
external callbackfunc
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y0glutesscallback(temp,i,which,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which)
else
   call f9y0glutesscallback(temp,i,which,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x0glutesscallback

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutesscallback(tess,which)
#else
subroutine f9x1glutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
integer(kind=glenum) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x1glutesscallback

#endif

#endif

#ifdef FULL_TESSCALLBACK

subroutine f9x01glutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_begin_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(xtype)
   use opengl_kinds
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: xtype
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y01glutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y01glutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x01glutesscallback

subroutine f9x02glutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_edgeflag_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(flag)
   use opengl_kinds
   integer(kind=glint) PASS_BY_REF, intent(in out) :: flag
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y02glutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y02glutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x02glutesscallback

subroutine f9x03rglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_vertex_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(vertex_data)
   use opengl_kinds
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: vertex_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y03rglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y03rglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x03rglutesscallback

subroutine f9x03iglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_vertex_int_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(vertex_data)
   use opengl_kinds
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: vertex_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y03iglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y03iglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x03iglutesscallback

subroutine f9x04glutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_end_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc()
   use opengl_kinds
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y04glutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y04glutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x04glutesscallback

subroutine f9x05glutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_error_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(errno)
   use opengl_kinds
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y05glutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y05glutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x05glutesscallback

subroutine f9x06rglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_combine_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(coords, vertex_data, weight, outdata)
   use opengl_kinds
   real(kind=gldouble), dimension(DEFER_DIM), intent(in out) :: coords
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: vertex_data, weight, &
                                                       outdata
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y06rglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y06rglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x06rglutesscallback

subroutine f9x06iglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_combine_int_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(coords, vertex_data, weight, outdata)
   use opengl_kinds
   real(kind=gldouble), dimension(DEFER_DIM), intent(in out) :: coords
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: vertex_data, outdata
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: weight
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y06iglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y06iglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x06iglutesscallback

subroutine f9x11rglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_begindata_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(xtype,polygon_data)
   use opengl_kinds
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: xtype
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y11rglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y11rglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x11rglutesscallback

subroutine f9x12rglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_edgeflagdata_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(flag,polygon_data)
   use opengl_kinds
   integer(kind=glint) PASS_BY_REF, intent(in out) :: flag
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y12rglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y12rglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x12rglutesscallback

subroutine f9x13rglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_vertexdata_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(vertex_data,polygon_data)
   use opengl_kinds
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: vertex_data
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y13rglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y13rglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x13rglutesscallback

subroutine f9x14rglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_enddata_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(polygon_data)
   use opengl_kinds
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y14rglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y14rglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x14rglutesscallback

subroutine f9x15rglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_errordata_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(errno,polygon_data)
   use opengl_kinds
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y15rglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y15rglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x15rglutesscallback

subroutine f9x16rglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_combinedata_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(coords, vertex_data, weight, outdata, polygon_data)
   use opengl_kinds
   real(kind=gldouble), dimension(DEFER_DIM), intent(in out) :: coords
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: vertex_data, weight, &
                                                       outdata
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y16rglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y16rglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x16rglutesscallback

subroutine f9x11iglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_begindata_int_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(xtype,polygon_data)
   use opengl_kinds
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: xtype
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y11iglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y11iglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x11iglutesscallback

subroutine f9x12iglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_efdata_int_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(flag,polygon_data)
   use opengl_kinds
   integer(kind=glint) PASS_BY_REF, intent(in out) :: flag
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y12iglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y12iglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x12iglutesscallback

subroutine f9x13iglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_vertexdata_int_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(vertex_data,polygon_data)
   use opengl_kinds
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: vertex_data
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y13iglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y13iglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x13iglutesscallback

subroutine f9x14iglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_enddata_int_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(polygon_data)
   use opengl_kinds
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y14iglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y14iglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x14iglutesscallback

subroutine f9x15iglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_errordata_int_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(errno,polygon_data)
   use opengl_kinds
   integer(kind=glenum) PASS_BY_REF, intent(in out) :: errno
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y15iglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y15iglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x15iglutesscallback

subroutine f9x16iglutesscallback(tess,which,callbackfunc)
type(glutesselator), pointer :: tess
type(f90glu_tess_combdata_int_type) INTENT_IN which
interface
   C_BINDING subroutine callbackfunc(coords, vertex_data, weight, outdata, polygon_data)
   use opengl_kinds
   real(kind=gldouble), dimension(DEFER_DIM), intent(in out) :: coords
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: vertex_data, &
                                                                outdata
   real(kind=glfloat), dimension(DEFER_DIM), intent(in out) :: weight
   integer(kind=glint), dimension(DEFER_DIM), intent(in out) :: polygon_data
   end subroutine callbackfunc
end interface
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
#ifdef OLDNULLFUNC
call f9y16iglutesscallback(temp,i,which%val,callbackfunc)
#else
if (myglucomparefuncs(glunullfunc,callbackfunc)==1) then
   call f9y1glutesscallback(temp,i,which%val)
else
   call f9y16iglutesscallback(temp,i,which%val,callbackfunc)
endif
#endif
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
return
end subroutine f9x16iglutesscallback

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x01aglutesscallback(tess,which)
#else
subroutine f9x01aglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_begin_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x01aglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x02aglutesscallback(tess,which)
#else
subroutine f9x02aglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_edgeflag_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x02aglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x03raglutesscallback(tess,which)
#else
subroutine f9x03raglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_vertex_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x03raglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x03iaglutesscallback(tess,which)
#else
subroutine f9x03iaglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_vertex_int_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x03iaglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x04aglutesscallback(tess,which)
#else
subroutine f9x04aglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_end_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x04aglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x05aglutesscallback(tess,which)
#else
subroutine f9x05aglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_error_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x05aglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x06raglutesscallback(tess,which)
#else
subroutine f9x06raglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_combine_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x06raglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x06iaglutesscallback(tess,which)
#else
subroutine f9x06iaglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_combine_int_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x06iaglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x11raglutesscallback(tess,which)
#else
subroutine f9x11raglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_begindata_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x11raglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x12raglutesscallback(tess,which)
#else
subroutine f9x12raglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_edgeflagdata_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x12raglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x13raglutesscallback(tess,which)
#else
subroutine f9x13raglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_vertexdata_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x13raglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x14raglutesscallback(tess,which)
#else
subroutine f9x14raglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_enddata_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x14raglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x15raglutesscallback(tess,which)
#else
subroutine f9x15raglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_errordata_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x15raglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x16raglutesscallback(tess,which)
#else
subroutine f9x16raglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_combinedata_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x16raglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x11iaglutesscallback(tess,which)
#else
subroutine f9x11iaglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_begindata_int_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x11iaglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x12iaglutesscallback(tess,which)
#else
subroutine f9x12iaglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_efdata_int_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x12iaglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x13iaglutesscallback(tess,which)
#else
subroutine f9x13iaglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_vertexdata_int_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x13iaglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x14iaglutesscallback(tess,which)
#else
subroutine f9x14iaglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_enddata_int_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x14iaglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x15iaglutesscallback(tess,which)
#else
subroutine f9x15iaglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_errordata_int_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x15iaglutesscallback

#ifdef NO_NULLFUNC
subroutine f9x16iaglutesscallback(tess,which)
#else
subroutine f9x16iaglutesscallback(tess,which,callbackfunc)
integer(kind=glint) INTENT_IN callbackfunc
#endif
type(glutesselator), pointer :: tess
type(f90glu_tess_combdata_int_type) INTENT_IN which
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
#ifdef NO_NULLFUNC
#else
if (callbackfunc == glunullfunc) then 
#endif
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutesscallback(temp,i,which%val)
do j=1,num_tess_cb
   do i=1,glcptr_length
      tess%addr(j+1)%addr(i:i) = char(temp(i+j*glcptr_length))
   end do
end do
#ifdef NO_NULLFUNC
#else
end if
#endif
return
end subroutine f9x16iaglutesscallback

#endif ! OLDNULLFUNC

#endif ! FULL_TESSCALLBACK

#ifdef GLU_1_1_TESS
#else
subroutine f9xglutessendcontour(tess)
type(glutesselator), pointer :: tess
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call fglutessendcontour(temp,i)
return
end subroutine f9xglutessendcontour

subroutine f9xglutessendpolygon(tess)
type(glutesselator), pointer :: tess
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call fglutessendpolygon(temp,i)
return
end subroutine f9xglutessendpolygon

subroutine f9xglutessnormal(tess,valuex,valuey,valuez)
type(glutesselator), pointer :: tess
real(kind=gldouble) INTENT_IN valuex, valuey, valuez
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call fglutessnormal(temp,i,valuex,valuey,valuez)
return
end subroutine f9xglutessnormal

subroutine f9xglutessproperty(tess,which,data)
type(glutesselator), pointer :: tess
integer(kind=glenum) INTENT_IN which
real(kind=gldouble) INTENT_IN data
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call fglutessproperty(temp,i,which,data)
return
end subroutine f9xglutessproperty
#endif

! I'm not sure if this will work, but since vertex_data is just passed on
! to a callback function, with F we'll try just passing it like with any
! other compiler.
subroutine f9x0glutessvertex(tess,location,vertex_data)
type(glutesselator), pointer :: tess
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
integer(kind=glint), dimension(DEFER_DIM), target, intent(in out) :: vertex_data
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y0glutessvertex(temp,i,location,vertex_data)
return
end subroutine f9x0glutessvertex

subroutine f9x1glutessvertex(tess,location,vertex_data)
type(glutesselator), pointer :: tess
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
real(kind=glfloat), dimension(DEFER_DIM), target, intent(in out) :: vertex_data
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y1glutessvertex(temp,i,location,vertex_data)
return
end subroutine f9x1glutessvertex

subroutine f9x4glutessvertex(tess,location,vertex_data)
type(glutesselator), pointer :: tess
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
type(glcptr) INTENT_IN vertex_data
integer(kind=gluint), dimension((num_tess_cb+2)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
do i=1,glcptr_length
   temp(i+(num_tess_cb+1)*glcptr_length) = ichar(vertex_data%addr(i:i))
end do
i = glcptr_length
call f9y4glutessvertex(temp,i,location)
return
end subroutine f9x4glutessvertex

subroutine f9x5glutessvertex(tess,location,vertex_data)
type(glutesselator), pointer :: tess
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
integer(kind=glbyte), dimension(DEFER_DIM), target, intent(in out) :: vertex_data
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y5glutessvertex(temp,i,location,vertex_data)
return
end subroutine f9x5glutessvertex

subroutine f9x6glutessvertex(tess,location,vertex_data)
type(glutesselator), pointer :: tess
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
integer(kind=glshort), dimension(DEFER_DIM), target, intent(in out) :: vertex_data
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y6glutessvertex(temp,i,location,vertex_data)
return
end subroutine f9x6glutessvertex

subroutine f9x7glutessvertex(tess,location,vertex_data)
type(glutesselator), pointer :: tess
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN location
real(kind=gldouble), dimension(DEFER_DIM), target, intent(in out) :: vertex_data
integer(kind=gluint), dimension((num_tess_cb+1)*glcptr_length) :: temp
integer(kind=glint) :: i,j
do i=1,glcptr_length
   temp(i) = ichar(tess%addr(1)%addr(i:i))
end do
do j=1,num_tess_cb
   do i=1,glcptr_length
      temp(i+j*glcptr_length) = ichar(tess%addr(j+1)%addr(i:i))
   end do
end do
i = glcptr_length
call f9y7glutessvertex(temp,i,location,vertex_data)
return
end subroutine f9x7glutessvertex

#ifdef PURE_FUNCTIONS
subroutine f9xgluunproject(winx,winy,winz,model,proj,view,objx,objy, &
   objz,errcode)
#else
function f9xgluunproject(winx,winy,winz,model,proj,view,objx,objy, &
   objz) result(resf9xgluunproject)
#endif
real(kind=gldouble) INTENT_IN winx, winy, winz
real(kind=gldouble), dimension(DEFER_DIM) INTENT_IN model, proj
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN view
real(kind=gldouble), intent(out) :: objx, objy, objz
#ifdef PURE_FUNCTIONS
integer(kind=glint), intent(out) :: errcode
#else
integer(kind=glint) :: resf9xgluunproject
#endif
#ifdef PURE_FUNCTIONS
call fgluunproject(winx,winy,winz,model,proj,view,objx,objy,objz,errcode)
#else
resf9xgluunproject = fgluunproject(winx,winy,winz,model,proj,view,objx,objy, &
   objz)
#endif
return
#ifdef PURE_FUNCTIONS
end subroutine f9xgluunproject
#else
end function f9xgluunproject
#endif

end module opengl_glu
