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
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

/* Make sure the extensions and version numbers are defined the same
   way in both fwrapglu and cwrapglu.  Not yet supported for Windows. */
#ifndef WIN32
#include "fppr.inc"
#endif

#define LOWERCASE 1
#define UNDERSCORE 2
#define UPPERCASE 3

/* used when Fortran must pass an assumed shape array */
#ifdef F_LANG
#define ARRAYREF [0]
#define ARRAY *
#define INVARRAY &
#else
#define ARRAYREF
#define ARRAY
#define INVARRAY
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

/* For OpenGL libraries that do not have the GLU 1.2 version of tesselation */

#ifndef GLU_VERSION_1_2
#define GLU_1_1_TESS
#endif

#ifdef GLU_1_1_TESS
#define GLUtesselator GLUtriangulatorObj
#endif

#ifdef F_LANG
#define PURE_FUNCTIONS
#endif

#ifdef REPLACE_CEIL
#define strlen lstrlen
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

#define end_of_string '\0'

#define f90gluquadricobj      GLUquadricObj **
#define f90glunurbsobj        GLUnurbsObj **
#define f90glutesselator      GLUtesselator **

/* procedure name mangling */

#if FNAME==LOWERCASE
#define mygluCompareFuncs              myglucomparefuncs
#define fgluBeginCurve                 fglubegincurve
#define fgluBeginPolygon               fglubeginpolygon
#define fgluBeginSurface               fglubeginsurface
#define fgluBeginTrim                  fglubegintrim
#define f9y0gluBuild1DMipmaps          f9y0glubuild1dmipmaps
#define f9y1gluBuild1DMipmaps          f9y1glubuild1dmipmaps
#define f9y4gluBuild1DMipmaps          f9y4glubuild1dmipmaps
#define f9y5gluBuild1DMipmaps          f9y5glubuild1dmipmaps
#define f9y6gluBuild1DMipmaps          f9y6glubuild1dmipmaps
#define f9y0gluBuild2DMipmaps          f9y0glubuild2dmipmaps
#define f9y1gluBuild2DMipmaps          f9y1glubuild2dmipmaps
#define f9y4gluBuild2DMipmaps          f9y4glubuild2dmipmaps
#define f9y5gluBuild2DMipmaps          f9y5glubuild2dmipmaps
#define f9y6gluBuild2DMipmaps          f9y6glubuild2dmipmaps
#define fgluCylinder                   fglucylinder
#define fgluDeleteNurbsRenderer        fgludeletenurbsrenderer
#define fgluDeleteQuadric              fgludeletequadric
#define fgluDeleteTess                 fgludeletetess
#define fgluDisk                       fgludisk
#define fgluEndCurve                   fgluendcurve
#define fgluEndPolygon                 fgluendpolygon
#define fgluEndSurface                 fgluendsurface
#define fgluEndTrim                    fgluendtrim
#define f9ygluErrorStringLen           f9ygluerrorstringlen
#define fgluGetNurbsProperty           fglugetnurbsproperty
#define f9ygluGetString                f9yglugetstring
#define f9ygluGetString1               f9yglugetstring1
#define f9ygluGetStringLen             f9yglugetstringlen
#define f9ygluGetStringFree            f9yglugetstringfree
#define fgluGetTessProperty            fglugettessproperty
#define fgluLoadSamplingMatrices       fgluloadsamplingmatrices
#define fgluLookAt                     fglulookat
#define f9ygluNewNurbsRenderer         f9yglunewnurbsrenderer
#define f9ygluNewQuadric               f9yglunewquadric
#define f9ygluNewTess                  f9yglunewtess
#define fgluNextContour                fglunextcontour
#define fgluNurbsCallback              fglunurbscallback
#define f9y1gluNurbsCallback           f9y1glunurbscallback
#define fgluNurbsCurve                 fglunurbscurve
#define fgluNurbsProperty              fglunurbsproperty
#define fgluNurbsSurface               fglunurbssurface
#define fgluOrtho2D                    fgluortho2d
#define fgluPartialDisk                fglupartialdisk
#define fgluPerspective                fgluperspective
#define fgluPickMatrix                 fglupickmatrix
#define fgluProject                    fgluproject
#define fgluPwlCurve                   fglupwlcurve
#define fgluQuadricCallback            fgluquadriccallback
#define f9y1gluQuadricCallback         f9y1gluquadriccallback
#define fgluQuadricDrawStyle           fgluquadricdrawstyle
#define fgluQuadricNormals             fgluquadricnormals
#define fgluQuadricOrientation         fgluquadricorientation
#define fgluQuadricTexture             fgluquadrictexture
#define f9egluQuadricTexture           f9egluquadrictexture
#define f9y0gluScaleImage              f9y0gluscaleimage
#define f9y1gluScaleImage              f9y1gluscaleimage
#define f9y4gluScaleImage              f9y4gluscaleimage
#define f9y5gluScaleImage              f9y5gluscaleimage
#define f9y6gluScaleImage              f9y6gluscaleimage
#define f9y01gluScaleImage             f9y01gluscaleimage
#define f9y04gluScaleImage             f9y04gluscaleimage
#define f9y05gluScaleImage             f9y05gluscaleimage
#define f9y06gluScaleImage             f9y06gluscaleimage
#define f9y10gluScaleImage             f9y10gluscaleimage
#define f9y14gluScaleImage             f9y14gluscaleimage
#define f9y15gluScaleImage             f9y15gluscaleimage
#define f9y16gluScaleImage             f9y16gluscaleimage
#define f9y40gluScaleImage             f9y40gluscaleimage
#define f9y41gluScaleImage             f9y41gluscaleimage
#define f9y45gluScaleImage             f9y45gluscaleimage
#define f9y46gluScaleImage             f9y46gluscaleimage
#define f9y50gluScaleImage             f9y50gluscaleimage
#define f9y51gluScaleImage             f9y51gluscaleimage
#define f9y54gluScaleImage             f9y54gluscaleimage
#define f9y56gluScaleImage             f9y56gluscaleimage
#define f9y60gluScaleImage             f9y60gluscaleimage
#define f9y61gluScaleImage             f9y61gluscaleimage
#define f9y64gluScaleImage             f9y64gluscaleimage
#define f9y65gluScaleImage             f9y65gluscaleimage
#define fgluSphere                     fglusphere
#define fgluTessBeginContour           fglutessbegincontour
#define f9y0gluTessBeginPolygon        f9y0glutessbeginpolygon
#define f9y1gluTessBeginPolygon        f9y1glutessbeginpolygon
#define f9y0gluTessCallback            f9y0glutesscallback
#define f9y1gluTessCallback            f9y1glutesscallback
#define f9y01gluTessCallback           f9y01glutesscallback
#define f9y02gluTessCallback           f9y02glutesscallback
#define f9y03igluTessCallback          f9y03iglutesscallback
#define f9y03rgluTessCallback          f9y03rglutesscallback
#define f9y04gluTessCallback           f9y04glutesscallback
#define f9y05gluTessCallback           f9y05glutesscallback
#define f9y06igluTessCallback          f9y06iglutesscallback
#define f9y06rgluTessCallback          f9y06rglutesscallback
#define f9y11igluTessCallback          f9y11iglutesscallback
#define f9y12igluTessCallback          f9y12iglutesscallback
#define f9y13igluTessCallback          f9y13iglutesscallback
#define f9y14igluTessCallback          f9y14iglutesscallback
#define f9y15igluTessCallback          f9y15iglutesscallback
#define f9y16igluTessCallback          f9y16iglutesscallback
#define f9y11rgluTessCallback          f9y11rglutesscallback
#define f9y12rgluTessCallback          f9y12rglutesscallback
#define f9y13rgluTessCallback          f9y13rglutesscallback
#define f9y14rgluTessCallback          f9y14rglutesscallback
#define f9y15rgluTessCallback          f9y15rglutesscallback
#define f9y16rgluTessCallback          f9y16rglutesscallback
#define fgluTessEndContour             fglutessendcontour
#define fgluTessEndPolygon             fglutessendpolygon
#define fgluTessNormal                 fglutessnormal
#define fgluTessProperty               fglutessproperty
#define f9y0gluTessVertex              f9y0glutessvertex
#define f9y1gluTessVertex              f9y1glutessvertex
#define f9y4gluTessVertex              f9y4glutessvertex
#define f9y5gluTessVertex              f9y5glutessvertex
#define f9y6gluTessVertex              f9y6glutessvertex
#define f9y7gluTessVertex              f9y7glutessvertex
#define fgluUnProject                  fgluunproject
#elif FNAME==UNDERSCORE
#define mygluCompareFuncs              myglucomparefuncs_
#define fgluBeginCurve                 fglubegincurve_
#define fgluBeginPolygon               fglubeginpolygon_
#define fgluBeginSurface               fglubeginsurface_
#define fgluBeginTrim                  fglubegintrim_
#define f9y0gluBuild1DMipmaps          f9y0glubuild1dmipmaps_
#define f9y1gluBuild1DMipmaps          f9y1glubuild1dmipmaps_
#define f9y4gluBuild1DMipmaps          f9y4glubuild1dmipmaps_
#define f9y5gluBuild1DMipmaps          f9y5glubuild1dmipmaps_
#define f9y6gluBuild1DMipmaps          f9y6glubuild1dmipmaps_
#define f9y0gluBuild2DMipmaps          f9y0glubuild2dmipmaps_
#define f9y1gluBuild2DMipmaps          f9y1glubuild2dmipmaps_
#define f9y4gluBuild2DMipmaps          f9y4glubuild2dmipmaps_
#define f9y5gluBuild2DMipmaps          f9y5glubuild2dmipmaps_
#define f9y6gluBuild2DMipmaps          f9y6glubuild2dmipmaps_
#define fgluCylinder                   fglucylinder_
#define fgluDeleteNurbsRenderer        fgludeletenurbsrenderer_
#define fgluDeleteQuadric              fgludeletequadric_
#define fgluDeleteTess                 fgludeletetess_
#define fgluDisk                       fgludisk_
#define fgluEndCurve                   fgluendcurve_
#define fgluEndPolygon                 fgluendpolygon_
#define fgluEndSurface                 fgluendsurface_
#define fgluEndTrim                    fgluendtrim_
#define f9ygluErrorStringLen           f9ygluerrorstringlen_
#define fgluGetNurbsProperty           fglugetnurbsproperty_
#define f9ygluGetString                f9yglugetstring_
#define f9ygluGetString1               f9yglugetstring1_
#define f9ygluGetStringLen             f9yglugetstringlen_
#define f9ygluGetStringFree            f9yglugetstringfree_
#define fgluGetTessProperty            fglugettessproperty_
#define fgluLoadSamplingMatrices       fgluloadsamplingmatrices_
#define fgluLookAt                     fglulookat_
#define f9ygluNewNurbsRenderer         f9yglunewnurbsrenderer_
#define f9ygluNewQuadric               f9yglunewquadric_
#define f9ygluNewTess                  f9yglunewtess_
#define fgluNextContour                fglunextcontour_
#define fgluNurbsCallback              fglunurbscallback_
#define f9y1gluNurbsCallback           f9y1glunurbscallback_
#define fgluNurbsCurve                 fglunurbscurve_
#define fgluNurbsProperty              fglunurbsproperty_
#define fgluNurbsSurface               fglunurbssurface_
#define fgluOrtho2D                    fgluortho2d_
#define fgluPartialDisk                fglupartialdisk_
#define fgluPerspective                fgluperspective_
#define fgluPickMatrix                 fglupickmatrix_
#define fgluProject                    fgluproject_
#define fgluPwlCurve                   fglupwlcurve_
#define fgluQuadricCallback            fgluquadriccallback_
#define f9y1gluQuadricCallback         f9y1gluquadriccallback_
#define fgluQuadricDrawStyle           fgluquadricdrawstyle_
#define fgluQuadricNormals             fgluquadricnormals_
#define fgluQuadricOrientation         fgluquadricorientation_
#define fgluQuadricTexture             fgluquadrictexture_
#define f9egluQuadricTexture           f9egluquadrictexture_
#define f9y0gluScaleImage              f9y0gluscaleimage_
#define f9y1gluScaleImage              f9y1gluscaleimage_
#define f9y4gluScaleImage              f9y4gluscaleimage_
#define f9y5gluScaleImage              f9y5gluscaleimage_
#define f9y6gluScaleImage              f9y6gluscaleimage_
#define f9y01gluScaleImage             f9y01gluscaleimage_
#define f9y04gluScaleImage             f9y04gluscaleimage_
#define f9y05gluScaleImage             f9y05gluscaleimage_
#define f9y06gluScaleImage             f9y06gluscaleimage_
#define f9y10gluScaleImage             f9y10gluscaleimage_
#define f9y14gluScaleImage             f9y14gluscaleimage_
#define f9y15gluScaleImage             f9y15gluscaleimage_
#define f9y16gluScaleImage             f9y16gluscaleimage_
#define f9y40gluScaleImage             f9y40gluscaleimage_
#define f9y41gluScaleImage             f9y41gluscaleimage_
#define f9y45gluScaleImage             f9y45gluscaleimage_
#define f9y46gluScaleImage             f9y46gluscaleimage_
#define f9y50gluScaleImage             f9y50gluscaleimage_
#define f9y51gluScaleImage             f9y51gluscaleimage_
#define f9y54gluScaleImage             f9y54gluscaleimage_
#define f9y56gluScaleImage             f9y56gluscaleimage_
#define f9y60gluScaleImage             f9y60gluscaleimage_
#define f9y61gluScaleImage             f9y61gluscaleimage_
#define f9y64gluScaleImage             f9y64gluscaleimage_
#define f9y65gluScaleImage             f9y65gluscaleimage_
#define fgluSphere                     fglusphere_
#define fgluTessBeginContour           fglutessbegincontour_
#define f9y0gluTessBeginPolygon        f9y0glutessbeginpolygon_
#define f9y1gluTessBeginPolygon        f9y1glutessbeginpolygon_
#define f9y0gluTessCallback            f9y0glutesscallback_
#define f9y1gluTessCallback            f9y1glutesscallback_
#define f9y01gluTessCallback           f9y01glutesscallback_
#define f9y02gluTessCallback           f9y02glutesscallback_
#define f9y03igluTessCallback          f9y03iglutesscallback_
#define f9y03rgluTessCallback          f9y03rglutesscallback_
#define f9y04gluTessCallback           f9y04glutesscallback_
#define f9y05gluTessCallback           f9y05glutesscallback_
#define f9y06igluTessCallback          f9y06iglutesscallback_
#define f9y06rgluTessCallback          f9y06rglutesscallback_
#define f9y11igluTessCallback          f9y11iglutesscallback_
#define f9y12igluTessCallback          f9y12iglutesscallback_
#define f9y13igluTessCallback          f9y13iglutesscallback_
#define f9y14igluTessCallback          f9y14iglutesscallback_
#define f9y15igluTessCallback          f9y15iglutesscallback_
#define f9y16igluTessCallback          f9y16iglutesscallback_
#define f9y11rgluTessCallback          f9y11rglutesscallback_
#define f9y12rgluTessCallback          f9y12rglutesscallback_
#define f9y13rgluTessCallback          f9y13rglutesscallback_
#define f9y14rgluTessCallback          f9y14rglutesscallback_
#define f9y15rgluTessCallback          f9y15rglutesscallback_
#define f9y16rgluTessCallback          f9y16rglutesscallback_
#define fgluTessEndContour             fglutessendcontour_
#define fgluTessEndPolygon             fglutessendpolygon_
#define fgluTessNormal                 fglutessnormal_
#define fgluTessProperty               fglutessproperty_
#define f9y0gluTessVertex              f9y0glutessvertex_
#define f9y1gluTessVertex              f9y1glutessvertex_
#define f9y4gluTessVertex              f9y4glutessvertex_
#define f9y5gluTessVertex              f9y5glutessvertex_
#define f9y6gluTessVertex              f9y6glutessvertex_
#define f9y7gluTessVertex              f9y7glutessvertex_
#define fgluUnProject                  fgluunproject_
#elif FNAME==UPPERCASE
#define mygluCompareFuncs              MYGLUCOMPAREFUNCS
#define fgluBeginCurve                 FGLUBEGINCURVE
#define fgluBeginPolygon               FGLUBEGINPOLYGON
#define fgluBeginSurface               FGLUBEGINSURFACE
#define fgluBeginTrim                  FGLUBEGINTRIM
#define f9y0gluBuild1DMipmaps          F9Y0GLUBUILD1DMIPMAPS
#define f9y1gluBuild1DMipmaps          F9Y1GLUBUILD1DMIPMAPS
#define f9y4gluBuild1DMipmaps          F9Y4GLUBUILD1DMIPMAPS
#define f9y5gluBuild1DMipmaps          F9Y5GLUBUILD1DMIPMAPS
#define f9y6gluBuild1DMipmaps          F9Y6GLUBUILD1DMIPMAPS
#define f9y0gluBuild2DMipmaps          F9Y0GLUBUILD2DMIPMAPS
#define f9y1gluBuild2DMipmaps          F9Y1GLUBUILD2DMIPMAPS
#define f9y4gluBuild2DMipmaps          F9Y4GLUBUILD2DMIPMAPS
#define f9y5gluBuild2DMipmaps          F9Y5GLUBUILD2DMIPMAPS
#define f9y6gluBuild2DMipmaps          F9Y6GLUBUILD2DMIPMAPS
#define fgluCylinder                   FGLUCYLINDER
#define fgluDeleteNurbsRenderer        FGLUDELETENURBSRENDERER
#define fgluDeleteQuadric              FGLUDELETEQUADRIC
#define fgluDeleteTess                 FGLUDELETETESS
#define fgluDisk                       FGLUDISK
#define fgluEndCurve                   FGLUENDCURVE
#define fgluEndPolygon                 FGLUENDPOLYGON
#define fgluEndSurface                 FGLUENDSURFACE
#define fgluEndTrim                    FGLUENDTRIM
#define f9ygluErrorStringLen           F9YGLUERRORSTRINGLEN
#define fgluGetNurbsProperty           FGLUGETNURBSPROPERTY
#define fgluGetString                  FGLUGETSTRING
#define fgluGetString1                 FGLUGETSTRING1
#define f9ygluGetString                F9YGLUGETSTRING
#define f9ygluGetStringLen             F9YGLUGETSTRINGLEN
#define f9ygluGetStringFree            F9YGLUGETSTRINGFREE
#define fgluGetTessProperty            FGLUGETTESSPROPERTY
#define fgluLoadSamplingMatrices       FGLULOADSAMPLINGMATRICES
#define fgluLookAt                     FGLULOOKAT
#define f9ygluNewNurbsRenderer         F9YGLUNEWNURBSRENDERER
#define f9ygluNewQuadric               F9YGLUNEWQUADRIC
#define f9ygluNewTess                  F9YGLUNEWTESS
#define fgluNextContour                FGLUNEXTCONTOUR
#define fgluNurbsCallback              FGLUNURBSCALLBACK
#define f9y1gluNurbsCallback           F9Y1GLUNURBSCALLBACK
#define fgluNurbsCurve                 FGLUNURBSCURVE
#define fgluNurbsProperty              FGLUNURBSPROPERTY
#define fgluNurbsSurface               FGLUNURBSSURFACE
#define fgluOrtho2D                    FGLUORTHO2D
#define fgluPartialDisk                FGLUPARTIALDISK
#define fgluPerspective                FGLUPERSPECTIVE
#define fgluPickMatrix                 FGLUPICKMATRIX
#define fgluProject                    FGLUPROJECT
#define fgluPwlCurve                   FGLUPWLCURVE
#define fgluQuadricCallback            FGLUQUADRICCALLBACK
#define f9y1gluQuadricCallback         F9Y1GLUQUADRICCALLBACK
#define fgluQuadricDrawStyle           FGLUQUADRICDRAWSTYLE
#define fgluQuadricNormals             FGLUQUADRICNORMALS
#define fgluQuadricOrientation         FGLUQUADRICORIENTATION
#define fgluQuadricTexture             FGLUQUADRICTEXTURE
#define f9egluQuadricTexture           F9EGLUQUADRICTEXTURE
#define f9y0gluScaleImage              F9Y0GLUSCALEIMAGE
#define f9y1gluScaleImage              F9Y1GLUSCALEIMAGE
#define f9y4gluScaleImage              F9Y4GLUSCALEIMAGE
#define f9y5gluScaleImage              F9Y5GLUSCALEIMAGE
#define f9y6gluScaleImage              F9Y6GLUSCALEIMAGE
#define f9y01gluScaleImage             F9Y01GLUSCALEIMAGE
#define f9y04gluScaleImage             F9Y04GLUSCALEIMAGE
#define f9y05gluScaleImage             F9Y05GLUSCALEIMAGE
#define f9y06gluScaleImage             F9Y06GLUSCALEIMAGE
#define f9y10gluScaleImage             F9Y10GLUSCALEIMAGE
#define f9y14gluScaleImage             F9Y14GLUSCALEIMAGE
#define f9y15gluScaleImage             F9Y15GLUSCALEIMAGE
#define f9y16gluScaleImage             F9Y16GLUSCALEIMAGE
#define f9y40gluScaleImage             F9Y40GLUSCALEIMAGE
#define f9y41gluScaleImage             F9Y41GLUSCALEIMAGE
#define f9y45gluScaleImage             F9Y45GLUSCALEIMAGE
#define f9y46gluScaleImage             F9Y46GLUSCALEIMAGE
#define f9y50gluScaleImage             F9Y50GLUSCALEIMAGE
#define f9y51gluScaleImage             F9Y51GLUSCALEIMAGE
#define f9y54gluScaleImage             F9Y54GLUSCALEIMAGE
#define f9y56gluScaleImage             F9Y56GLUSCALEIMAGE
#define f9y60gluScaleImage             F9Y60GLUSCALEIMAGE
#define f9y61gluScaleImage             F9Y61GLUSCALEIMAGE
#define f9y64gluScaleImage             F9Y64GLUSCALEIMAGE
#define f9y65gluScaleImage             F9Y65GLUSCALEIMAGE
#define fgluSphere                     FGLUSPHERE
#define fgluTessBeginContour           FGLUTESSBEGINCONTOUR
#define f9y0gluTessBeginPolygon        F9Y0GLUTESSBEGINPOLYGON
#define f9y1gluTessBeginPolygon        F9Y1GLUTESSBEGINPOLYGON
#define f9y0gluTessCallback            F9Y0GLUTESSCALLBACK
#define f9y1gluTessCallback            F9Y1GLUTESSCALLBACK
#define f9y01gluTessCallback           F9Y01GLUTESSCALLBACK
#define f9y02gluTessCallback           F9Y02GLUTESSCALLBACK
#define f9y03igluTessCallback          F9Y03IGLUTESSCALLBACK
#define f9y03rgluTessCallback          F9Y03RGLUTESSCALLBACK
#define f9y04gluTessCallback           F9Y04GLUTESSCALLBACK
#define f9y05gluTessCallback           F9Y05GLUTESSCALLBACK
#define f9y06igluTessCallback          F9Y06IGLUTESSCALLBACK
#define f9y06rgluTessCallback          F9Y06RGLUTESSCALLBACK
#define f9y11igluTessCallback          F9Y11IGLUTESSCALLBACK
#define f9y12igluTessCallback          F9Y12IGLUTESSCALLBACK
#define f9y13igluTessCallback          F9Y13IGLUTESSCALLBACK
#define f9y14igluTessCallback          F9Y14IGLUTESSCALLBACK
#define f9y15igluTessCallback          F9Y15IGLUTESSCALLBACK
#define f9y16igluTessCallback          F9Y16IGLUTESSCALLBACK
#define f9y11rgluTessCallback          F9Y11RGLUTESSCALLBACK
#define f9y12rgluTessCallback          F9Y12RGLUTESSCALLBACK
#define f9y13rgluTessCallback          F9Y13RGLUTESSCALLBACK
#define f9y14rgluTessCallback          F9Y14RGLUTESSCALLBACK
#define f9y15rgluTessCallback          F9Y15RGLUTESSCALLBACK
#define f9y16rgluTessCallback          F9Y16RGLUTESSCALLBACK
#define fgluTessEndContour             FGLUTESSENDCONTOUR
#define fgluTessEndPolygon             FGLUTESSENDPOLYGON
#define fgluTessNormal                 FGLUTESSNORMAL
#define fgluTessProperty               FGLUTESSPROPERTY
#define f9y0gluTessVertex              F9Y0GLUTESSVERTEX
#define f9y1gluTessVertex              F9Y1GLUTESSVERTEX
#define f9y4gluTessVertex              F9Y4GLUTESSVERTEX
#define f9y5gluTessVertex              F9Y5GLUTESSVERTEX
#define f9y6gluTessVertex              F9Y6GLUTESSVERTEX
#define f9y7gluTessVertex              F9Y7GLUTESSVERTEX
#define fgluUnProject                  FGLUUNPROJECT
#endif

/* name mangling for the special routines for passing arrays in F */

#ifdef F_LANG
#if FNAME==UNDERSCORE

#define glu_F_uint_array_free         glu_f_uint_array_free_
#define glu_F_uint_array_recv_beg     glu_f_uint_array_recv_beg_
#define glu_F_uint_array_recv         glu_f_uint_array_recv_

#endif
#endif

/* type definitions for callback functions */

typedef void (APIENTRY *f90glNurbErrorFunc) (GLenum *);
typedef void (APIENTRY *f90glQuadErrorFunc) (GLenum *);
typedef void (APIENTRY *f90glTessBeginFunc) (GLenum *);
typedef void (APIENTRY *f90glTessEdgeFlagFunc) (GLint *);
typedef void (APIENTRY *f90glTessVertexFunc) (void *);
typedef void (APIENTRY *f90glTessEndFunc) (void);
typedef void (APIENTRY *f90glTessErrorFunc) (GLenum *);
typedef void (APIENTRY *f90glTessCombineFunc) (GLdouble *, void *, GLfloat *, void **);
typedef void (APIENTRY *f90glTessBeginDataFunc) (GLenum *, void *);
typedef void (APIENTRY *f90glTessEdgeFlagDataFunc) (GLint *, void *);
typedef void (APIENTRY *f90glTessEndDataFunc) (void *);
typedef void (APIENTRY *f90glTessVertexDataFunc) (void *, void *);
typedef void (APIENTRY *f90glTessErrorDataFunc) (GLenum *, void *);
typedef void (APIENTRY *f90glTessCombineDataFunc) (GLdouble *, void *, GLfloat *, void **, void *);

/* places to store the pointers for the current nurb, quadric and tess */

struct {
   GLUnurbsObj* nurb;
   f90glNurbErrorFunc error;
} CurrentNurb;

struct {
   GLUquadricObj* quad;
   f90glQuadErrorFunc error;
} CurrentQuadric;

struct {
   GLUtesselator* tess;
   f90glTessBeginFunc begin;
   f90glTessEdgeFlagFunc edgeflag;
   f90glTessVertexFunc vertex;
   f90glTessEndFunc end;
   f90glTessErrorFunc error;
   f90glTessCombineFunc combine;
   f90glTessBeginDataFunc begindata;
   f90glTessEdgeFlagDataFunc edgeflagdata;
   f90glTessEndDataFunc enddata;
   f90glTessVertexDataFunc vertexdata;
   f90glTessErrorDataFunc errordata;
   f90glTessCombineDataFunc combinedata;
} CurrentTess;

/* number of entries in the above structures */

#define NURB_SIZE 2
#define QUAD_SIZE 2
#define TESS_SIZE 13

/* places to put the result of glGetString */
GLint glugetstring_length;
GLuint *glugetstring_intstr;

#ifdef F_LANG
/* places to put arrays passed one value at a time by F */
GLuint   *F_uint_array;
GLint     F_array_size;
#endif

#ifdef F_LANG
/* special procedures for F passing arrays */

GLint glu_F_uint_array_free(void)
{
   free(F_uint_array);
   return 1;
}

GLint glu_F_uint_array_recv_beg(void)
{
   return F_array_size;
}

GLuint glu_F_uint_array_recv(f90glint pos)
{
   return F_uint_array[(*pos)-1];
}

#endif /* F_LANG */

/* callback wrappers called by the C routines to call the Fortran callbacks */

void APIENTRY f90fcbNurbError(GLenum errno)
{
   CurrentNurb.error(&errno);
}

void APIENTRY f90fcbQuadError(GLenum errno)
{
   CurrentQuadric.error(&errno);
}

void APIENTRY f90fcbTessBegin(GLenum type)
{
   CurrentTess.begin(&type);
}

void APIENTRY f90fcbTessEdgeFlag(GLboolean flag)
{
   GLint int_flag;
   int_flag = flag;
   CurrentTess.edgeflag(&int_flag);
}

void APIENTRY f90fcbTessVertex(void *vertex_data)
{
   CurrentTess.vertex(INVARRAY vertex_data);
}

void APIENTRY f90fcbTessEnd(void)
{
   CurrentTess.end();
}

void APIENTRY f90fcbTessError(GLenum errno)
{
   CurrentTess.error(&errno);
}

void APIENTRY f90fcbTessCombine(GLdouble *coords, void **vertex_data, GLfloat *weight,
   void **outData)
{
   CurrentTess.combine(INVARRAY coords,INVARRAY vertex_data,INVARRAY weight,INVARRAY outData);
}

void APIENTRY f90fcbTessBeginData(GLenum type, void *polygon_data)
{
   CurrentTess.begindata(&type,INVARRAY polygon_data);
}

void APIENTRY f90fcbTessEdgeFlagData(GLboolean flag, void *polygon_data)
{
   GLint int_flag;
   int_flag = flag;
   CurrentTess.edgeflagdata(&int_flag,INVARRAY polygon_data);
}

void APIENTRY f90fcbTessEndData(void *polygon_data)
{
   CurrentTess.enddata(INVARRAY polygon_data);
}

void APIENTRY f90fcbTessVertexData(void *vertex_data, void *polygon_data)
{
   CurrentTess.vertexdata(INVARRAY vertex_data,INVARRAY polygon_data);
}

void APIENTRY f90fcbTessErrorData(GLenum errno, void *polygon_data)
{
   CurrentTess.errordata(&errno,INVARRAY polygon_data);
}

void APIENTRY f90fcbTessCombineData(GLdouble *coords, void **vertex_data,
   GLfloat *weight, void **outData, void *polygon_data)
{
   CurrentTess.combinedata(INVARRAY coords,INVARRAY vertex_data,INVARRAY weight,INVARRAY outData,INVARRAY polygon_data);
}

/* set the pointers in the current nurb, quadric and tesselator */

void APIENTRY f90gluiSetCurrentNurb(f90glunurbsobj nurbObject)
{
   CurrentNurb.nurb = (GLUnurbsObj *) nurbObject[0];
   CurrentNurb.error = (f90glNurbErrorFunc) nurbObject[1];
}

void APIENTRY f90gluiSetCurrentQuadric(f90gluquadricobj quadObject)
{
   CurrentQuadric.quad = (GLUquadricObj *) quadObject[0];
   CurrentQuadric.error = (f90glQuadErrorFunc) quadObject[1];
}

void APIENTRY f90gluiSetCurrentTess(f90glutesselator tessObject)
{
   CurrentTess.tess = (GLUtesselator*) tessObject[0];
   CurrentTess.begin = (f90glTessBeginFunc) tessObject[1];
   CurrentTess.edgeflag = (f90glTessEdgeFlagFunc) tessObject[2];
   CurrentTess.vertex = (f90glTessVertexFunc) tessObject[3];
   CurrentTess.end = (f90glTessEndFunc) tessObject[4];
   CurrentTess.error = (f90glTessErrorFunc) tessObject[5];
   CurrentTess.combine = (f90glTessCombineFunc) tessObject[6];
   CurrentTess.begindata = (f90glTessBeginDataFunc) tessObject[7];
   CurrentTess.edgeflagdata = (f90glTessEdgeFlagDataFunc) tessObject[8];
   CurrentTess.enddata = (f90glTessEndDataFunc) tessObject[9];
   CurrentTess.vertexdata = (f90glTessVertexDataFunc) tessObject[10];
   CurrentTess.errordata = (f90glTessErrorDataFunc) tessObject[11];
   CurrentTess.combinedata = (f90glTessCombineDataFunc) tessObject[12];
}

/* compare two pointers for equality -- used to check for GLUNULLFUNC */
GLint APIENTRY mygluCompareFuncs(void * ptr1, void * ptr2)
{
   if (ptr1==ptr2) {
      return 1;
   }
   else {
      return 0;
   }
}

/* GLU 1.2 */

void APIENTRY fgluBeginCurve( f90gluint ARRAY nurbaddr, f90glint nbytes )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluBeginCurve( (GLUnurbsObj *) nurbObject[0] );
}

void APIENTRY fgluBeginPolygon( f90gluint ARRAY tessaddr, f90glint nbytes )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluBeginPolygon( (GLUtesselator *) tessObject[0] );
}

void APIENTRY fgluBeginSurface( f90gluint ARRAY nurbaddr, f90glint nbytes )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluBeginSurface( (GLUnurbsObj *) nurbObject[0] );
}

void APIENTRY fgluBeginTrim( f90gluint ARRAY nurbaddr, f90glint nbytes )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluBeginTrim( (GLUnurbsObj *) nurbObject[0] );
}

GLint APIENTRY f9y0gluBuild1DMipmaps( f90glenum target, 
   f90glint components, f90glsizei width,
   f90glenum format, f90glenum type, f90glvoid ARRAY data,
   f90glint length )
{

GLint retval;

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdata;
   int i;
   newdata = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newdata[i]=(GLubyte)data ARRAYREF [i];
   retval = gluBuild1DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdata;
   int i;
   newdata = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newdata[i]=(GLbyte)data ARRAYREF [i];
   retval = gluBuild1DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newdata;
   int i;
   newdata = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newdata[i]=(GLushort)data ARRAYREF [i];
   retval = gluBuild1DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newdata;
   int i;
   newdata = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newdata[i]=(GLshort)data ARRAYREF [i];
   retval = gluBuild1DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 default:
   retval = gluBuild1DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) data  ARRAYREF );
   break;
}
return retval;
}

GLint APIENTRY f9y1gluBuild1DMipmaps( f90glenum target, 
   f90glint components, f90glsizei width,
   f90glenum format, f90glenum type, f90glvoid ARRAY data )
{
   return gluBuild1DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) data ARRAYREF );
}

GLint APIENTRY f9y4gluBuild1DMipmaps( f90glenum target, 
   f90glint components, f90glsizei width,
   f90glenum format, f90glenum type, f90gluint ARRAY temp,
   f90glint nbytes )
{
   GLvoid *data;
   unsigned char *p;
   int i;
   p = (unsigned char *) &data;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   return gluBuild1DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) data );
}

GLint APIENTRY f9y5gluBuild1DMipmaps( f90glenum target, 
   f90glint components, f90glsizei width,
   f90glenum format, f90glenum type, f90glvoid ARRAY data )
{
   return gluBuild1DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) data ARRAYREF);
}

GLint APIENTRY f9y6gluBuild1DMipmaps( f90glenum target, 
   f90glint components, f90glsizei width,
   f90glenum format, f90glenum type, f90glshort ARRAY data,
   f90glint length )
{

GLint retval;

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdata;
   int i;
   newdata = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newdata[i]=(GLubyte)data ARRAYREF [i];
   retval = gluBuild1DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdata;
   int i;
   newdata = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newdata[i]=(GLbyte)data ARRAYREF [i];
   retval = gluBuild1DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 default:
   retval = gluBuild1DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) data  ARRAYREF );
   break;
}
return retval;
}

GLint APIENTRY f9y0gluBuild2DMipmaps( f90glenum target, 
   f90glint components, f90glsizei width, f90glsizei height,
   f90glenum format, f90glenum type, f90glvoid ARRAY data,
   f90glint length )
{

GLint retval;

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdata;
   int i;
   newdata = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newdata[i]=(GLubyte)data ARRAYREF [i];
   retval = gluBuild2DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdata;
   int i;
   newdata = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newdata[i]=(GLbyte)data ARRAYREF [i];
   retval = gluBuild2DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newdata;
   int i;
   newdata = (GLushort *) malloc((*length)*sizeof(GLushort));
   for (i=0;i<(*length);i++) newdata[i]=(GLushort)data ARRAYREF [i];
   retval = gluBuild2DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newdata;
   int i;
   newdata = (GLshort *) malloc((*length)*sizeof(GLshort));
   for (i=0;i<(*length);i++) newdata[i]=(GLshort)data ARRAYREF [i];
   retval = gluBuild2DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 default:
   retval = gluBuild2DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) data  ARRAYREF );
   break;
}
return retval;
}

GLint APIENTRY f9y1gluBuild2DMipmaps( f90glenum target, 
   f90glint components, f90glsizei width, f90glsizei height,
   f90glenum format, f90glenum type, f90glvoid ARRAY data )
{
   return gluBuild2DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) data ARRAYREF);
}

GLint APIENTRY f9y4gluBuild2DMipmaps( f90glenum target, 
   f90glint components, f90glsizei width, f90glsizei height,
   f90glenum format, f90glenum type, f90gluint ARRAY temp,
   f90glint nbytes )
{
   GLvoid *data;
   unsigned char *p;
   int i;
   p = (unsigned char *) &data;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   return gluBuild2DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) data );
}

GLint APIENTRY f9y5gluBuild2DMipmaps( f90glenum target, 
   f90glint components, f90glsizei width, f90glsizei height,
   f90glenum format, f90glenum type, f90glvoid ARRAY data )
{
   return gluBuild2DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) data ARRAYREF);
}

GLint APIENTRY f9y6gluBuild2DMipmaps( f90glenum target, 
   f90glint components, f90glsizei width, f90glsizei height,
   f90glenum format, f90glenum type, f90glshort ARRAY data,
   f90glint length )
{

GLint retval;

switch (*type) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdata;
   int i;
   newdata = (GLubyte *) malloc((*length)*sizeof(GLubyte));
   for (i=0;i<(*length);i++) newdata[i]=(GLubyte)data ARRAYREF [i];
   retval = gluBuild2DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdata;
   int i;
   newdata = (GLbyte *) malloc((*length)*sizeof(GLbyte));
   for (i=0;i<(*length);i++) newdata[i]=(GLbyte)data ARRAYREF [i];
   retval = gluBuild2DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) newdata );
   free(newdata);
   break;
  }
 default:
   retval = gluBuild2DMipmaps( (GLenum) *target, (GLint) *components,
      (GLsizei) *width, (GLsizei) *height, 
      (GLenum) *format, (GLenum) *type, (GLvoid *) data ARRAYREF);
   break;
}
return retval;
}

void APIENTRY fgluCylinder( f90gluint ARRAY quadaddr, f90glint nbytes, f90gldouble base,
    f90gldouble top, f90gldouble height, f90glint slices, f90glint stacks )
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   gluCylinder( (GLUquadricObj *) quadObject[0], (GLdouble) *base, (GLdouble) *top,
    (GLdouble) *height, (GLint) *slices, (GLint) *stacks );
}

void APIENTRY fgluDeleteNurbsRenderer( f90gluint ARRAY nurbaddr, f90glint nbytes )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluDeleteNurbsRenderer( (GLUnurbsObj *) nurbObject[0] );
}

void APIENTRY fgluDeleteQuadric( f90gluint ARRAY quadaddr, f90glint nbytes )
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   gluDeleteQuadric( (GLUquadricObj *) quadObject[0] );
}

void APIENTRY fgluDeleteTess( f90gluint ARRAY tessaddr, f90glint nbytes )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluDeleteTess( (GLUtesselator *) tessObject[0] );
}

void APIENTRY fgluDisk( f90gluint ARRAY quadaddr, f90glint nbytes, f90gldouble inner,
    f90gldouble outer, f90glint slices, f90glint loops )
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   gluDisk( (GLUquadricObj *) quadObject[0], (GLdouble) *inner, (GLdouble) *outer,
    (GLint) *slices, (GLint) *loops );
}

void APIENTRY fgluEndCurve( f90gluint ARRAY nurbaddr, f90glint nbytes )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluEndCurve( (GLUnurbsObj *) nurbObject[0] );
}

void APIENTRY fgluEndPolygon( f90gluint ARRAY tessaddr, f90glint nbytes )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluEndPolygon( (GLUtesselator *) tessObject[0] );
}

void APIENTRY fgluEndSurface( f90gluint ARRAY nurbaddr, f90glint nbytes )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluEndSurface( (GLUnurbsObj *) nurbObject[0] );
}

void APIENTRY fgluEndTrim( f90gluint ARRAY nurbaddr, f90glint nbytes )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluEndTrim( (GLUnurbsObj *) nurbObject[0] );
}

GLint APIENTRY f9ygluErrorStringLen(f90glenum error)
{
   int i;
   char *str;
   str = (char *) gluErrorString( (GLenum) *error);
   if (str) {
      glugetstring_length = strlen(str);
      glugetstring_intstr = (GLuint *)malloc((glugetstring_length)*sizeof(GLuint));
      for (i=0;i<glugetstring_length;i++) glugetstring_intstr[i] = (GLuint) str[i];
   }else{
      glugetstring_length = 0;
      glugetstring_intstr = (GLuint *)malloc(sizeof(GLuint));
   }
   return (GLint) glugetstring_length;
}

void APIENTRY fgluGetNurbsProperty( f90gluint ARRAY nurbaddr, f90glint nbytes, f90glenum property,
    f90glfloat data )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluGetNurbsProperty( (GLUnurbsObj *) nurbObject[0], (GLenum) *property,
    (GLfloat *) data );
}

GLint APIENTRY f9ygluGetStringLen(f90glenum name)
{
   int i;
   char *str;
   str = (char *) gluGetString( (GLenum) *name);
   if (str) {
      glugetstring_length = strlen(str);
      glugetstring_intstr = (GLuint *)malloc((glugetstring_length)*sizeof(GLuint));
      for (i=0;i<glugetstring_length;i++) glugetstring_intstr[i] = (GLuint) str[i];
   }else{
      glugetstring_length = 0;
      glugetstring_intstr = (GLuint *)malloc(sizeof(GLuint));
   }
   return (GLint) glugetstring_length;
}
 
void APIENTRY f9ygluGetString(f90gluint ARRAY int_str)
{
   int i;
   for (i=0;i<glugetstring_length;i++) int_str ARRAYREF [i] = glugetstring_intstr[i];
}

GLuint APIENTRY f9ygluGetString1(f90glint pos)
{
   return (GLuint) glugetstring_intstr[(*pos)-1];
}

#ifdef PURE_FUNCTIONS
GLint APIENTRY f9ygluGetStringFree( void )
#else
void APIENTRY f9ygluGetStringFree( void )
#endif
{
   glugetstring_length = 0;
   free(glugetstring_intstr);
#ifdef PURE_FUNCTIONS
   return 1;
#endif
}

#ifndef GLU_1_1_TESS

void APIENTRY fgluGetTessProperty( f90gluint ARRAY tessaddr, f90glint nbytes,
    f90glenum which, f90gldouble data )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluGetTessProperty( (GLUtesselator *) tessObject[0], (GLenum) *which,
    (GLdouble *) data );
}

#endif

void APIENTRY fgluLoadSamplingMatrices( f90gluint ARRAY nurbaddr, f90glint nbytes, f90glfloat ARRAY model,
    f90glfloat ARRAY perspective, f90glint ARRAY view )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   gluLoadSamplingMatrices( (GLUnurbsObj *) nurbObject[0], (GLfloat *) model ARRAYREF,
    (GLfloat *) perspective ARRAYREF, (GLint *) view  ARRAYREF);
}

void APIENTRY fgluLookAt( f90gldouble eyeX, f90gldouble eyeY,
    f90gldouble eyeZ, f90gldouble centerX, f90gldouble centerY,
    f90gldouble centerZ, f90gldouble upX, f90gldouble upY, f90gldouble upZ )
{
   gluLookAt( (GLdouble) *eyeX, (GLdouble) *eyeY, (GLdouble) *eyeZ,
    (GLdouble) *centerX, (GLdouble) *centerY, (GLdouble) *centerZ,
    (GLdouble) *upX, (GLdouble) *upY, (GLdouble) *upZ );
}

#ifdef F_LANG
GLint APIENTRY f9ygluNewNurbsRenderer(f90gluint ARRAY temp, f90glint nbytes)
#else
void APIENTRY f9ygluNewNurbsRenderer(f90gluint ARRAY temp, f90glint nbytes)
#endif
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   nurbObject[0] = (GLUnurbsObj *) gluNewNurbsRenderer();
   nurbObject[1] = (GLUnurbsObj *) NULL;
   f90gluiSetCurrentNurb((f90glunurbsobj)nurbObject);
   p = (unsigned char *) nurbObject;
#ifdef F_LANG
   F_uint_array = (GLuint *) malloc(NURB_SIZE*(*nbytes)*sizeof(GLuint));
   for (i=0; i<NURB_SIZE*(*nbytes); i++) {F_uint_array[i]=(GLuint)*p; p++;}
   F_array_size = NURB_SIZE*(*nbytes);
#else
   for (i=0; i<NURB_SIZE*(*nbytes); i++) {temp ARRAYREF [i]=(GLuint)*p; p++;}
#endif
#ifdef F_LANG
   return 1;
#endif
}

#ifdef F_LANG
GLint APIENTRY f9ygluNewQuadric(f90gluint ARRAY temp, f90glint nbytes)
#else
void APIENTRY f9ygluNewQuadric(f90gluint ARRAY temp, f90glint nbytes)
#endif
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   quadObject[0] = (GLUquadricObj *) gluNewQuadric();
   quadObject[1] = (GLUquadricObj *) NULL;
   f90gluiSetCurrentQuadric((f90gluquadricobj)quadObject);
   p = (unsigned char *) quadObject;
#ifdef F_LANG
   F_uint_array = (GLuint *) malloc(QUAD_SIZE*(*nbytes)*sizeof(GLuint));
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) {F_uint_array[i]=(GLuint)*p; p++;}
   F_array_size = QUAD_SIZE*(*nbytes);
#else
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) {temp ARRAYREF [i]=(GLuint)*p; p++;}
#endif
#ifdef F_LANG
   return 1;
#endif
}

#ifdef F_LANG
GLint APIENTRY f9ygluNewTess(f90gluint ARRAY temp, f90glint nbytes)
#else
void APIENTRY f9ygluNewTess(f90gluint ARRAY temp, f90glint nbytes)
#endif
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   tessObject[0] = (GLUtesselator *) gluNewTess();
   for (i=1; i<TESS_SIZE; i++) tessObject[i] = (GLUtesselator *) NULL;
   f90gluiSetCurrentTess(tessObject);
   p = (unsigned char *) tessObject;
#ifdef F_LANG
   F_uint_array = (GLuint *) malloc(TESS_SIZE*(*nbytes)*sizeof(GLuint));
   for (i=0; i<TESS_SIZE*(*nbytes); i++) {F_uint_array[i]=(GLuint)*p; p++;}
   F_array_size = TESS_SIZE*(*nbytes);
#else
   for (i=0; i<TESS_SIZE*(*nbytes); i++) {temp ARRAYREF [i]=(GLuint)*p; p++;}
#endif
#ifdef F_LANG
   return 1;
#endif
}

void APIENTRY fgluNextContour( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum type )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluNextContour( (GLUtesselator *) tessObject[0], (GLenum) *type );
}

void APIENTRY f9y1gluNurbsCallback( f90gluint ARRAY nurbaddr, f90glint nbytes, f90glenum which)
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   nurbObject[1] = (GLUnurbsObj *) NULL;
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluNurbsCallback( (GLUnurbsObj *) nurbObject[0], (GLenum) *which, NULL);
   f90gluiSetCurrentNurb(nurbObject);
   for (i=0; i<NURB_SIZE*(*nbytes); i++) {nurbaddr ARRAYREF [i]=(GLuint)*p; p++;}
}

void APIENTRY fgluNurbsCallback( f90gluint ARRAY nurbaddr, f90glint nbytes, f90glenum which,
    GLvoid (*CallBackFunc)() )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   switch (*which) {
   case GLU_ERROR:
      nurbObject[1] = (GLUnurbsObj *) CallBackFunc;
      gluNurbsCallback( (GLUnurbsObj *) nurbObject[0], (GLenum) *which,
#ifdef WIN32
         (GLUnurbsErrorProc)
#endif
         f90fcbNurbError );
      break;
   }
   f90gluiSetCurrentNurb(nurbObject);
   for (i=0; i<NURB_SIZE*(*nbytes); i++) {nurbaddr [i]=(GLuint)*p; p++;}
}

void APIENTRY fgluNurbsCurve( f90gluint ARRAY nurbaddr, f90glint nbytes, f90glint knotCount,
    f90glfloat ARRAY knots, f90glint stride, f90glfloat ARRAY control, f90glint order,
    f90glenum type )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluNurbsCurve( (GLUnurbsObj *) nurbObject[0], (GLint) *knotCount,
    (GLfloat *) knots ARRAYREF, (GLint) *stride, (GLfloat *) control ARRAYREF, (GLint) *order,
    (GLenum) *type );
}

void APIENTRY fgluNurbsProperty( f90gluint ARRAY nurbaddr, f90glint nbytes, f90glenum property,
    f90glfloat value )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluNurbsProperty( (GLUnurbsObj *) nurbObject[0], (GLenum) *property,
    (GLfloat) *value );
}

void APIENTRY fgluNurbsSurface( f90gluint ARRAY nurbaddr, f90glint nbytes, f90glint sKnotCount,
    f90glfloat ARRAY sKnots, f90glint tKnotCount, f90glfloat ARRAY tKnots, f90glint sStride,
    f90glint tStride, f90glfloat ARRAY control, f90glint sOrder, f90glint tOrder,
    f90glenum type )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluNurbsSurface( (GLUnurbsObj *) nurbObject[0], (GLint) *sKnotCount,
    (GLfloat *) sKnots ARRAYREF, (GLint) *tKnotCount, (GLfloat *) tKnots ARRAYREF, (GLint) *sStride,
    (GLint) *tStride, (GLfloat *) control ARRAYREF, (GLint) *sOrder, (GLint) *tOrder,
    (GLenum) *type );
}

void APIENTRY fgluOrtho2D( f90gldouble left, f90gldouble right,
    f90gldouble bottom, f90gldouble top )
{
   gluOrtho2D( (GLdouble) *left, (GLdouble) *right, (GLdouble) *bottom,
    (GLdouble) *top );
}

void APIENTRY fgluPartialDisk( f90gluint ARRAY quadaddr, f90glint nbytes, f90gldouble inner,
    f90gldouble outer, f90glint slices, f90glint loops, f90gldouble start,
    f90gldouble sweep )
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   gluPartialDisk( (GLUquadricObj *) quadObject[0], (GLdouble) *inner,
    (GLdouble) *outer, (GLint) *slices, (GLint) *loops, (GLdouble) *start,
    (GLdouble) *sweep );
}

void APIENTRY fgluPerspective( f90gldouble fovy, f90gldouble aspect,
    f90gldouble zNear, f90gldouble zFar )
{
   gluPerspective( (GLdouble) *fovy, (GLdouble) *aspect, (GLdouble) *zNear,
    (GLdouble) *zFar );
}

void APIENTRY fgluPickMatrix( f90gldouble x, f90gldouble y, f90gldouble delX,
    f90gldouble delY, f90glint ARRAY viewport )
{
   gluPickMatrix( (GLdouble) *x, (GLdouble) *y, (GLdouble) *delX,
    (GLdouble) *delY, (GLint *) viewport ARRAYREF);
}

#ifdef PURE_FUNCTIONS
void APIENTRY fgluProject( f90gldouble objX, f90gldouble objY,
    f90gldouble objZ, f90gldouble ARRAY model, f90gldouble ARRAY proj, f90glint ARRAY view,
    f90gldouble winX, f90gldouble winY, f90gldouble winZ, f90glint errcode)
#else
GLint APIENTRY fgluProject( f90gldouble objX, f90gldouble objY,
    f90gldouble objZ, f90gldouble ARRAY model, f90gldouble ARRAY proj, f90glint ARRAY view,
    f90gldouble winX, f90gldouble winY, f90gldouble winZ )
#endif
{
#ifdef PURE_FUNCTIONS
   *errcode = gluProject( (GLdouble) *objX, (GLdouble) *objY,
    (GLdouble) *objZ, (GLdouble *) model ARRAYREF, (GLdouble *) proj ARRAYREF, (GLint *) view ARRAYREF,
    (GLdouble *) winX, (GLdouble *) winY, (GLdouble *) winZ );
#else
   return (GLint) gluProject( (GLdouble) *objX, (GLdouble) *objY,
    (GLdouble) *objZ, (GLdouble *) model ARRAYREF, (GLdouble *) proj ARRAYREF, (GLint *) view ARRAYREF,
    (GLdouble *) winX, (GLdouble *) winY, (GLdouble *) winZ );
#endif
}

void APIENTRY fgluPwlCurve( f90gluint ARRAY nurbaddr, f90glint nbytes, f90glint count,
    f90glfloat ARRAY data, f90glint stride, f90glenum type )
{
   GLUnurbsObj *nurbObject[NURB_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) nurbObject;
   for (i=0; i<NURB_SIZE*(*nbytes); i++) { *p = (unsigned char)nurbaddr ARRAYREF [i]; p++;}
   if (CurrentNurb.nurb != (GLUnurbsObj *) nurbObject[0]) {
      f90gluiSetCurrentNurb(nurbObject);
   }
   gluPwlCurve( (GLUnurbsObj *) nurbObject[0], (GLint) *count, (GLfloat *) data ARRAYREF,
    (GLint) *stride, (GLenum) *type );
}

void APIENTRY f9y1gluQuadricCallback( f90gluint ARRAY quadaddr, f90glint nbytes,
    f90glenum which)
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   quadObject[1] = (GLUquadricObj *) NULL;
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   gluQuadricCallback( (GLUquadricObj *) quadObject[0], (GLenum) *which,
      NULL );
   f90gluiSetCurrentQuadric(quadObject);
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) {quadaddr ARRAYREF [i]=(GLuint)*p; p++;}
}

void APIENTRY fgluQuadricCallback( f90gluint ARRAY quadaddr, f90glint nbytes, f90glenum which,
    GLvoid (*CallBackFunc)() )
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   switch (*which) {
   case GLU_ERROR:
      quadObject[1] = (GLUquadricObj *) CallBackFunc;
      gluQuadricCallback( (GLUquadricObj *) quadObject[0], (GLenum) *which,
#ifdef WIN32
         (GLUquadricErrorProc)
#endif
         f90fcbQuadError );
      break;
   }
   f90gluiSetCurrentQuadric(quadObject);
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) {quadaddr ARRAYREF [i]=(GLuint)*p; p++;}
}

void APIENTRY fgluQuadricDrawStyle( f90gluint ARRAY quadaddr, f90glint nbytes, f90glenum draw )
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   gluQuadricDrawStyle( (GLUquadricObj *) quadObject[0], (GLenum) *draw );
}

void APIENTRY fgluQuadricNormals( f90gluint ARRAY quadaddr, f90glint nbytes, f90glenum normal )
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   gluQuadricNormals( (GLUquadricObj *) quadObject[0], (GLenum) *normal );
}

void APIENTRY fgluQuadricOrientation( f90gluint ARRAY quadaddr, f90glint nbytes,
    f90glenum orientation )
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   gluQuadricOrientation( (GLUquadricObj *) quadObject[0], (GLenum) *orientation );
}

void APIENTRY fgluQuadricTexture( f90gluint ARRAY quadaddr, f90glint nbytes,
    f90glboolean texture )
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   gluQuadricTexture( (GLUquadricObj *) quadObject[0], (GLboolean) *texture );
}

void APIENTRY f9egluQuadricTexture(f90gluint ARRAY quadaddr, f90glint nbytes,f90glint texture )
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   gluQuadricTexture( (GLUquadricObj *) quadObject[0], (GLboolean) *texture);
}

/* TEMP f9y0gluScaleImage and f9y6gluScaleImage require that typeIn == typeOut
   This should only affect compilers that do not support short integers */

#ifdef PURE_FUNCTIONS
void APIENTRY f9y0gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout, f90glint errcode )
#else
GLint APIENTRY f9y0gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout )
#endif
{
GLint retval;

switch (*typeIn) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdatain, *newdataout;
   int i;
   newdatain = (GLubyte *) malloc((*lengthin)*sizeof(GLubyte));
   newdataout = (GLubyte *) malloc((*lengthout)*sizeof(GLubyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLubyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLint)newdataout[i];
   free(newdatain);
   free(newdataout);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdatain, *newdataout;
   int i;
   newdatain = (GLbyte *) malloc((*lengthin)*sizeof(GLbyte));
   newdataout = (GLbyte *) malloc((*lengthout)*sizeof(GLbyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLbyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLint)newdataout[i];
   free(newdatain);
   free(newdataout);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newdatain, *newdataout;
   int i;
   newdatain = (GLushort *) malloc((*lengthin)*sizeof(GLushort));
   newdataout = (GLushort *) malloc((*lengthout)*sizeof(GLushort));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLushort)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLint)newdataout[i];
   free(newdatain);
   free(newdataout);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newdatain, *newdataout;
   int i;
   newdatain = (GLshort *) malloc((*lengthin)*sizeof(GLshort));
   newdataout = (GLshort *) malloc((*lengthout)*sizeof(GLshort));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLshort)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLint)newdataout[i];
   free(newdatain);
   free(newdataout);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y1gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut, f90glint errcode )
#else
GLint APIENTRY f9y1gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut )
#endif
{
   GLint retval;
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF);
#ifdef PURE_FUNCTIONS
   *errcode = retval;
#else
   return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y4gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90gluint ARRAY temp,
   f90glint nbytes, f90glint errcode )
#else
GLint APIENTRY f9y4gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90gluint ARRAY temp,
   f90glint nbytes )
#endif
{
   GLint retval;
   GLvoid *dataIn, *dataOut;
   unsigned char *p;
   int i;
   p = (unsigned char *) &dataIn;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   p = (unsigned char *) &dataOut;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i+(*nbytes)]; p++;}
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut );
#ifdef PURE_FUNCTIONS
   *errcode = retval;
#else
   return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y5gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, const void * ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut, f90glint errcode )
#else
GLint APIENTRY f9y5gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, const void * ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut )
#endif
{
   GLint retval;
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF);
#ifdef PURE_FUNCTIONS
   *errcode = retval;
#else
   return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y6gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout, f90glint errcode )
#else
GLint APIENTRY f9y6gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glshort ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glshort ARRAY dataOut,
   f90glint lengthout )
#endif
{
   GLint retval;

switch (*typeIn) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdatain, *newdataout;
   int i;
   newdatain = (GLubyte *) malloc((*lengthin)*sizeof(GLubyte));
   newdataout = (GLubyte *) malloc((*lengthout)*sizeof(GLubyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLubyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLshort)newdataout[i];
   free(newdatain);
   free(newdataout);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdatain, *newdataout;
   int i;
   newdatain = (GLbyte *) malloc((*lengthin)*sizeof(GLbyte));
   newdataout = (GLbyte *) malloc((*lengthout)*sizeof(GLbyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLbyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLshort)newdataout[i];
   free(newdatain);
   free(newdataout);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y01gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint errcode )
#else
GLint APIENTRY f9y01gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut)
#endif
{
GLint retval;

switch (*typeIn) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdatain;
   int i;
   newdatain = (GLubyte *) malloc((*lengthin)*sizeof(GLubyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLubyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdatain;
   int i;
   newdatain = (GLbyte *) malloc((*lengthin)*sizeof(GLbyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLbyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newdatain;
   int i;
   newdatain = (GLushort *) malloc((*lengthin)*sizeof(GLushort));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLushort)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newdatain;
   int i;
   newdatain = (GLshort *) malloc((*lengthin)*sizeof(GLshort));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLshort)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y04gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90gluint ARRAY temp,
   f90glint nbytes, f90glint errcode )
#else
GLint APIENTRY f9y04gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90gluint ARRAY temp,
   f90glint nbytes )
#endif
{
   GLint retval;
   GLvoid *dataOut;
   unsigned char *p;
   int i;
   p = (unsigned char *) &dataOut;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i+(*nbytes)]; p++;}

switch (*typeIn) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdatain;
   int i;
   newdatain = (GLubyte *) malloc((*lengthin)*sizeof(GLubyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLubyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut);
   free(newdatain);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdatain;
   int i;
   newdatain = (GLbyte *) malloc((*lengthin)*sizeof(GLbyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLbyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut);
   free(newdatain);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newdatain;
   int i;
   newdatain = (GLushort *) malloc((*lengthin)*sizeof(GLushort));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLushort)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut);
   free(newdatain);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newdatain;
   int i;
   newdatain = (GLshort *) malloc((*lengthin)*sizeof(GLshort));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLshort)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut);
   free(newdatain);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut);
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y05gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint errcode )
#else
GLint APIENTRY f9y05gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut)
#endif
{
GLint retval;

switch (*typeIn) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdatain;
   int i;
   newdatain = (GLubyte *) malloc((*lengthin)*sizeof(GLubyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLubyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdatain;
   int i;
   newdatain = (GLbyte *) malloc((*lengthin)*sizeof(GLbyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLbyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newdatain;
   int i;
   newdatain = (GLushort *) malloc((*lengthin)*sizeof(GLushort));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLushort)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newdatain;
   int i;
   newdatain = (GLshort *) malloc((*lengthin)*sizeof(GLshort));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLshort)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y06gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout, f90glint errcode )
#else
GLint APIENTRY f9y06gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glshort ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glshort ARRAY dataOut,
   f90glint lengthout )
#endif
{
   GLint retval;

switch (*typeOut) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdataout;
   int i;
   newdataout = (GLubyte *) malloc((*lengthout)*sizeof(GLubyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLshort)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdataout;
   int i;
   newdataout = (GLbyte *) malloc((*lengthout)*sizeof(GLbyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLshort)newdataout[i];
   free(newdataout);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF);
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y10gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout, f90glint errcode )
#else
GLint APIENTRY f9y10gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout )
#endif
{
GLint retval;

switch (*typeOut) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdataout;
   int i;
   newdataout = (GLubyte *) malloc((*lengthout)*sizeof(GLubyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdataout;
   int i;
   newdataout = (GLbyte *) malloc((*lengthout)*sizeof(GLbyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newdataout;
   int i;
   newdataout = (GLushort *) malloc((*lengthout)*sizeof(GLushort));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newdataout;
   int i;
   newdataout = (GLshort *) malloc((*lengthout)*sizeof(GLshort));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y14gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90gluint ARRAY temp,
   f90glint nbytes, f90glint errcode )
#else
GLint APIENTRY f9y14gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90gluint ARRAY temp,
   f90glint nbytes )
#endif
{
   GLint retval;
   GLvoid *dataOut;
   unsigned char *p;
   int i;
   p = (unsigned char *) &dataOut;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i+(*nbytes)]; p++;}
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut );
#ifdef PURE_FUNCTIONS
   *errcode = retval;
#else
   return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y15gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut, f90glint errcode )
#else
GLint APIENTRY f9y15gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut )
#endif
{
   GLint retval;
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF);
#ifdef PURE_FUNCTIONS
   *errcode = retval;
#else
   return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y16gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout, f90glint errcode )
#else
GLint APIENTRY f9y16gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glshort ARRAY dataIn,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glshort ARRAY dataOut,
   f90glint lengthout )
#endif
{
   GLint retval;

switch (*typeOut) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdataout;
   int i;
   newdataout = (GLubyte *) malloc((*lengthout)*sizeof(GLubyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLshort)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdataout;
   int i;
   newdataout = (GLbyte *) malloc((*lengthout)*sizeof(GLbyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLshort)newdataout[i];
   free(newdataout);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF);
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y40gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90gluint ARRAY temp, f90glint nbytes,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout, f90glint errcode )
#else
GLint APIENTRY f9y40gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90gluint ARRAY temp, f90glint nbytes,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout )
#endif
{
   GLint retval;
   GLvoid *dataIn;
   unsigned char *p;
   int i;
   p = (unsigned char *) &dataIn;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}

switch (*typeOut) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdataout;
   int i;
   newdataout = (GLubyte *) malloc((*lengthout)*sizeof(GLubyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn , (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdataout;
   int i;
   newdataout = (GLbyte *) malloc((*lengthout)*sizeof(GLbyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn , (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newdataout;
   int i;
   newdataout = (GLushort *) malloc((*lengthout)*sizeof(GLushort));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn , (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newdataout;
   int i;
   newdataout = (GLshort *) malloc((*lengthout)*sizeof(GLshort));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn , (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn , (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF);
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y41gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90gluint ARRAY temp, f90glint nbytes,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint errcode )
#else
GLint APIENTRY f9y41gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90gluint ARRAY temp, f90glint nbytes,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut)
#endif
{
   GLint retval;
   GLvoid *dataIn;
   unsigned char *p;
   int i;
   p = (unsigned char *) &dataIn;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
#ifdef PURE_FUNCTIONS
   *errcode = retval;
#else
   return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y45gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90gluint ARRAY temp, f90glint nbytes,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint errcode )
#else
GLint APIENTRY f9y45gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90gluint ARRAY temp, f90glint nbytes,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut)
#endif
{
   GLint retval;
   GLvoid *dataIn;
   unsigned char *p;
   int i;
   p = (unsigned char *) &dataIn;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
#ifdef PURE_FUNCTIONS
   *errcode = retval;
#else
   return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y46gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90gluint ARRAY temp, f90glint nbytes,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout, f90glint errcode )
#else
GLint APIENTRY f9y46gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90gluint ARRAY temp, f90glint nbytes,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glshort ARRAY dataOut,
   f90glint lengthout )
#endif
{
   GLint retval;
   GLvoid *dataIn;
   unsigned char *p;
   int i;
   p = (unsigned char *) &dataIn;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i]; p++;}

switch (*typeIn) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdataout;
   int i;
   newdataout = (GLubyte *) malloc((*lengthout)*sizeof(GLubyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLshort)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdataout;
   int i;
   newdataout = (GLbyte *) malloc((*lengthout)*sizeof(GLbyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF [i] = (GLshort)newdataout[i];
   free(newdataout);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn, (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y50gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout, f90glint errcode )
#else
GLint APIENTRY f9y50gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout )
#endif
{
GLint retval;

switch (*typeOut) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdataout;
   int i;
   newdataout = (GLubyte *) malloc((*lengthout)*sizeof(GLubyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdataout;
   int i;
   newdataout = (GLbyte *) malloc((*lengthout)*sizeof(GLbyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_UNSIGNED_SHORT:
  {
   GLushort *newdataout;
   int i;
   newdataout = (GLushort *) malloc((*lengthout)*sizeof(GLushort));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_SHORT:
  {
   GLshort *newdataout;
   int i;
   newdataout = (GLshort *) malloc((*lengthout)*sizeof(GLshort));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLint)newdataout[i];
   free(newdataout);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF);
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y51gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut, f90glint errcode )
#else
GLint APIENTRY f9y51gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut )
#endif
{
   GLint retval;
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF);
#ifdef PURE_FUNCTIONS
   *errcode = retval;
#else
   return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y54gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90gluint ARRAY temp,
   f90glint nbytes, f90glint errcode )
#else
GLint APIENTRY f9y54gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glsizei wOut,
   f90glsizei hOut, f90glenum typeOut, f90gluint ARRAY temp,
   f90glint nbytes )
#endif
{
   GLint retval;
   GLvoid *dataOut;
   unsigned char *p;
   int i;
   p = (unsigned char *) &dataOut;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i+(*nbytes)]; p++;}
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut );
#ifdef PURE_FUNCTIONS
   *errcode = retval;
#else
   return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y56gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint lengthout, f90glint errcode )
#else
GLint APIENTRY f9y56gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glshort ARRAY dataIn,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glshort ARRAY dataOut,
   f90glint lengthout )
#endif
{
   GLint retval;

switch (*typeOut) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdataout;
   int i;
   newdataout = (GLubyte *) malloc((*lengthout)*sizeof(GLubyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLshort)newdataout[i];
   free(newdataout);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdataout;
   int i;
   newdataout = (GLbyte *) malloc((*lengthout)*sizeof(GLbyte));
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) newdataout );
   for (i=0;i<(*lengthout);i++) dataOut ARRAYREF[i] = (GLshort)newdataout[i];
   free(newdataout);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF, (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF);
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y60gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint errcode )
#else
GLint APIENTRY f9y60gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glshort ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glint ARRAY dataOut)
#endif
{
   GLint retval;

switch (*typeIn) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdatain;
   int i;
   newdatain = (GLubyte *) malloc((*lengthin)*sizeof(GLubyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLubyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdatain;
   int i;
   newdatain = (GLbyte *) malloc((*lengthin)*sizeof(GLbyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLbyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y61gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint errcode )
#else
GLint APIENTRY f9y61gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glshort ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glshort ARRAY dataOut)
#endif
{
   GLint retval;

switch (*typeIn) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdatain;
   int i;
   newdatain = (GLubyte *) malloc((*lengthin)*sizeof(GLubyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLubyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdatain;
   int i;
   newdatain = (GLbyte *) malloc((*lengthin)*sizeof(GLbyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLbyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y64gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90gluint ARRAY temp,
   f90glint nbytes, f90glint errcode )
#else
GLint APIENTRY f9y64gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glshort ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90gluint ARRAY temp,
   f90glint nbytes)
#endif
{
   GLint retval;
   GLvoid *dataOut;
   unsigned char *p;
   int i;
   p = (unsigned char *) &dataOut;
   for (i=0; i<(*nbytes); i++) { *p = (unsigned char)temp ARRAYREF [i+(*nbytes)]; p++;}

switch (*typeIn) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdatain;
   int i;
   newdatain = (GLubyte *) malloc((*lengthin)*sizeof(GLubyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLubyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut );
   free(newdatain);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdatain;
   int i;
   newdatain = (GLbyte *) malloc((*lengthin)*sizeof(GLbyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLbyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut );
   free(newdatain);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut );
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}

#ifdef PURE_FUNCTIONS
void APIENTRY f9y65gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glvoid ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glvoid ARRAY dataOut,
   f90glint errcode )
#else
GLint APIENTRY f9y65gluScaleImage( f90glenum format, f90glsizei wIn, 
   f90glsizei hIn, f90glenum typeIn, f90glshort ARRAY dataIn, f90glint lengthin,
   f90glsizei wOut, f90glsizei hOut, f90glenum typeOut, f90glshort ARRAY dataOut)
#endif
{
   GLint retval;

switch (*typeIn) {

 case GL_UNSIGNED_BYTE: case GL_BITMAP:
  {
   GLubyte *newdatain;
   int i;
   newdatain = (GLubyte *) malloc((*lengthin)*sizeof(GLubyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLubyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 case GL_BYTE:
  {
   GLbyte *newdatain;
   int i;
   newdatain = (GLbyte *) malloc((*lengthin)*sizeof(GLbyte));
   for (i=0;i<(*lengthin);i++) newdatain[i]=(GLbyte)dataIn ARRAYREF [i];
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) newdatain, (GLsizei) *wOut, 
      (GLsizei) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   free(newdatain);
   break;
  }
 default:
   retval = gluScaleImage( (GLenum) *format, (GLsizei) *wIn,
      (GLsizei) *hIn, (GLenum) *typeIn, (const void *) dataIn ARRAYREF , (GLsizei) *wOut, 
      (GLenum) *hOut, (GLenum) *typeOut, (GLvoid *) dataOut  ARRAYREF );
   break;
}
#ifdef PURE_FUNCTIONS
*errcode = retval;
#else
return retval;
#endif
}
void APIENTRY fgluSphere( f90gluint ARRAY quadaddr, f90glint nbytes, f90gldouble radius,
    f90glint slices, f90glint stacks )
{
   GLUquadricObj *quadObject[QUAD_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) quadObject;
   for (i=0; i<QUAD_SIZE*(*nbytes); i++) { *p = (unsigned char)quadaddr ARRAYREF [i]; p++;}
   if (CurrentQuadric.quad != (GLUquadricObj *) quadObject[0]) {
      f90gluiSetCurrentQuadric(quadObject);
   }
   gluSphere( (GLUquadricObj *) quadObject[0], (GLdouble) *radius, (GLint) *slices,
    (GLint) *stacks );
}

#ifndef GLU_1_1_TESS

void APIENTRY fgluTessBeginContour( f90gluint ARRAY tessaddr, f90glint nbytes )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessBeginContour( (GLUtesselator *) tessObject[0] );
}

void APIENTRY f9y0gluTessBeginPolygon( f90gluint ARRAY tessaddr, f90glint nbytes,
    f90glvoid ARRAY data )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessBeginPolygon( (GLUtesselator *) tessObject[0], (GLvoid *) data ARRAYREF);
}

void APIENTRY f9y1gluTessBeginPolygon( f90gluint ARRAY tessaddr, f90glint nbytes,
    f90glvoid ARRAY data )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessBeginPolygon( (GLUtesselator *) tessObject[0], (GLvoid *) data  ARRAYREF );
}

#endif

void APIENTRY f9y1gluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which)
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   switch (*which) {
#ifdef GLU_1_1_TESS
   case GLU_BEGIN:
#else
   case GLU_TESS_BEGIN:
#endif
      tessObject[1] = (GLUtesselator *) NULL;
      break;
#ifdef GLU_1_1_TESS
   case GLU_EDGE_FLAG:
#else
   case GLU_TESS_EDGE_FLAG:
#endif
      tessObject[2] = (GLUtesselator *) NULL;
      break;
#ifdef GLU_1_1_TESS
   case GLU_VERTEX:
#else
   case GLU_TESS_VERTEX:
#endif
      tessObject[3] = (GLUtesselator *) NULL;
      break;
#ifdef GLU_1_1_TESS
   case GLU_END:
#else
   case GLU_TESS_END:
#endif
      tessObject[4] = (GLUtesselator *) NULL;
      break;
#ifdef GLU_1_1_TESS
   case GLU_ERROR:
#else
   case GLU_TESS_ERROR:
#endif
      tessObject[5] = (GLUtesselator *) NULL;
      break;
#ifndef GLU_1_1_TESS
   case GLU_TESS_COMBINE:
      tessObject[6] = (GLUtesselator *) NULL;
      break;
   case GLU_TESS_BEGIN_DATA:
      tessObject[7] = (GLUtesselator *) NULL;
      break;
   case GLU_TESS_EDGE_FLAG_DATA:
      tessObject[8] = (GLUtesselator *) NULL;
      break;
   case GLU_TESS_END_DATA:
      tessObject[9] = (GLUtesselator *) NULL;
      break;
   case GLU_TESS_VERTEX_DATA:
      tessObject[10] = (GLUtesselator *) NULL;
      break;
   case GLU_TESS_ERROR_DATA:
      tessObject[11] = (GLUtesselator *) NULL;
      break;
   case GLU_TESS_COMBINE_DATA:
      tessObject[12] = (GLUtesselator *) NULL;
      break;
#endif
   }
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
      NULL );
   f90gluiSetCurrentTess(tessObject);
   for (i=0; i<TESS_SIZE*(*nbytes); i++) {tessaddr ARRAYREF [i]=(GLuint)*p; p++;}
}

void APIENTRY f9y0gluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   switch (*which) {
#ifdef GLU_1_1_TESS
   case GLU_BEGIN:
#else
   case GLU_TESS_BEGIN:
#endif
      tessObject[1] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessBeginProc)
#endif
      f90fcbTessBegin );
      break;
#ifdef GLU_1_1_TESS
   case GLU_EDGE_FLAG:
#else
   case GLU_TESS_EDGE_FLAG:
#endif
      tessObject[2] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessEdgeFlagProc)
#endif
      f90fcbTessEdgeFlag );
      break;
#ifdef GLU_1_1_TESS
   case GLU_VERTEX:
#else
   case GLU_TESS_VERTEX:
#endif
      tessObject[3] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessVertexProc)
#endif
      f90fcbTessVertex );
      break;
#ifdef GLU_1_1_TESS
   case GLU_END:
#else
   case GLU_TESS_END:
#endif
      tessObject[4] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessEndProc)
#endif
      f90fcbTessEnd );
      break;
#ifdef GLU_1_1_TESS
   case GLU_ERROR:
#else
   case GLU_TESS_ERROR:
#endif
      tessObject[5] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessErrorProc)
#endif
      f90fcbTessError );
      break;
#ifndef GLU_1_1_TESS
   case GLU_TESS_COMBINE:
      tessObject[6] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessCombineProc)
#endif
      f90fcbTessCombine );
      break;
   case GLU_TESS_BEGIN_DATA:
      tessObject[7] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessBeginDataProc)
#endif
      f90fcbTessBeginData );
      break;
   case GLU_TESS_EDGE_FLAG_DATA:
      tessObject[8] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessEdgeFlagDataProc)
#endif
      f90fcbTessEdgeFlagData );
      break;
   case GLU_TESS_END_DATA:
      tessObject[9] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessEndDataProc)
#endif
      f90fcbTessEndData );
      break;
   case GLU_TESS_VERTEX_DATA:
      tessObject[10] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessVertexDataProc)
#endif
      f90fcbTessVertexData );
      break;
   case GLU_TESS_ERROR_DATA:
      tessObject[11] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessErrorDataProc)
#endif
      f90fcbTessErrorData );
      break;
   case GLU_TESS_COMBINE_DATA:
      tessObject[12] = (GLUtesselator *) CallBackFunc;
      gluTessCallback( (GLUtesselator *) tessObject[0], (GLenum) *which,
#ifdef WIN32
      (GLUtessCombineDataProc)
#endif
      f90fcbTessCombineData );
      break;
#endif
   }
   f90gluiSetCurrentTess(tessObject);
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) {tessaddr ARRAYREF [i]=(GLuint)*p; p++;}
}

void APIENTRY f9y01gluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y02gluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y03igluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y03rgluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y04gluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y05gluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y06igluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y06rgluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y11igluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y12igluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y13igluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y14igluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y15igluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y16igluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y11rgluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y12rgluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y13rgluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y14rgluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y15rgluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

void APIENTRY f9y16rgluTessCallback( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    GLvoid (APIENTRY *CallBackFunc)() )
{
  f9y0gluTessCallback(tessaddr, nbytes, which, CallBackFunc );
}

#ifndef GLU_1_1_TESS

void APIENTRY fgluTessEndContour( f90gluint ARRAY tessaddr, f90glint nbytes )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessEndContour( (GLUtesselator *) tessObject[0] );
}

void APIENTRY fgluTessEndPolygon( f90gluint ARRAY tessaddr, f90glint nbytes )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessEndPolygon( (GLUtesselator *) tessObject[0] );
}

void APIENTRY fgluTessNormal( f90gluint ARRAY tessaddr, f90glint nbytes, f90gldouble valueX,
    f90gldouble valueY, f90gldouble valueZ )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessNormal( (GLUtesselator *) tessObject[0], (GLdouble) *valueX,
    (GLdouble) *valueY, (GLdouble) *valueZ );
}

void APIENTRY fgluTessProperty( f90gluint ARRAY tessaddr, f90glint nbytes, f90glenum which,
    f90gldouble data )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessProperty( (GLUtesselator *) tessObject[0], (GLenum) *which,
    (GLdouble) *data );
}

#endif

void APIENTRY f9y0gluTessVertex( f90gluint ARRAY tessaddr, f90glint nbytes,
    f90gldouble ARRAY location, f90glvoid ARRAY data )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessVertex( (GLUtesselator *) tessObject[0], (GLdouble *) location ARRAYREF ,
    (GLvoid *) data  ARRAYREF );
}

void APIENTRY f9y1gluTessVertex( f90gluint ARRAY tessaddr, f90glint nbytes,
    f90gldouble ARRAY location, f90glvoid ARRAY data )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessVertex( (GLUtesselator *) tessObject[0], (GLdouble *) location ARRAYREF ,
    (GLvoid *) data  ARRAYREF );
}

void APIENTRY f9y4gluTessVertex( f90gluint ARRAY tessaddr, f90glint nbytes,
    f90gldouble ARRAY location )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   GLvoid *data;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   p = (unsigned char *) &data;
/* end of tessaddr contains the address for data */
   for (i=TESS_SIZE*(*nbytes); i<TESS_SIZE*(1+(*nbytes)); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessVertex( (GLUtesselator *) tessObject[0], (GLdouble *) location ARRAYREF,
    (GLvoid *) data );
}

void APIENTRY f9y5gluTessVertex( f90gluint ARRAY tessaddr, f90glint nbytes,
    f90gldouble ARRAY location, f90glvoid ARRAY data )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessVertex( (GLUtesselator *) tessObject[0], (GLdouble *) location ARRAYREF ,
    (GLvoid *) data  ARRAYREF );
}

void APIENTRY f9y6gluTessVertex( f90gluint ARRAY tessaddr, f90glint nbytes,
    f90gldouble ARRAY location, f90glvoid ARRAY data )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessVertex( (GLUtesselator *) tessObject[0], (GLdouble *) location ARRAYREF ,
    (GLvoid *) data  ARRAYREF );
}

void APIENTRY f9y7gluTessVertex( f90gluint ARRAY tessaddr, f90glint nbytes,
    f90gldouble ARRAY location, f90glvoid ARRAY data )
{
   GLUtesselator *tessObject[TESS_SIZE];
   unsigned char *p;
   int i;
   p = (unsigned char *) tessObject;
   for (i=0; i<TESS_SIZE*(*nbytes); i++) { *p = (unsigned char)tessaddr ARRAYREF [i]; p++;}
   if (CurrentTess.tess != (GLUtesselator *) tessObject[0]) {
      f90gluiSetCurrentTess(tessObject);
   }
   gluTessVertex( (GLUtesselator *) tessObject[0], (GLdouble *) location ARRAYREF ,
    (GLvoid *) data  ARRAYREF );
}

#ifdef PURE_FUNCTIONS
void APIENTRY fgluUnProject( f90gldouble winX, f90gldouble winY,
    f90gldouble winZ, f90gldouble ARRAY model, f90gldouble ARRAY proj, f90glint ARRAY view,
    f90gldouble objX, f90gldouble objY, f90gldouble objZ, f90glint errcode)
#else
GLint APIENTRY fgluUnProject( f90gldouble winX, f90gldouble winY,
    f90gldouble winZ, f90gldouble ARRAY model, f90gldouble ARRAY proj, f90glint ARRAY view,
    f90gldouble objX, f90gldouble objY, f90gldouble objZ )
#endif
{
#ifdef PURE_FUNCTIONS
   *errcode = gluUnProject( (GLdouble) *winX, (GLdouble) *winY,
    (GLdouble) *winZ, (GLdouble *) model ARRAYREF, (GLdouble *) proj ARRAYREF, (GLint *) view ARRAYREF,
    (GLdouble *) objX, (GLdouble *) objY, (GLdouble *) objZ );
#else
   return (GLint) gluUnProject( (GLdouble) *winX, (GLdouble) *winY,
    (GLdouble) *winZ, (GLdouble *) model ARRAYREF, (GLdouble *) proj ARRAYREF, (GLint *) view ARRAYREF,
    (GLdouble *) objX, (GLdouble *) objY, (GLdouble *) objZ );
#endif
}
