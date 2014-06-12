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

/* include support for new routines in GLUT 3.7 */
#define GLUT37

#ifdef LF90
#define LF9X
#endif
#ifdef LF95
#define LF9X
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef WIN32
#include <windows.h>
#endif

#if defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

/* make sure glut.h didn't undef APIENTRY */

#ifndef APIENTRY
#define APIENTRY
#endif

/* There are some differences between what I put in the patches for GLUT 3.6
   and what finally made it into GLUT 3.7.  This flag selects which to use. The
   old forms will probably be removed sometime well after GLUT 3.7 release */

#if (GLUT_XLIB_IMPLEMENTATION < 14 )
#define TMP_GLUT_VERSION 36
#else
#define TMP_GLUT_VERSION 37
#endif

/* glutf90.h moved from the GLUT source directory to include/GL at GLUT 3.7 */
#if (TMP_GLUT_VERSION == 36)
#include "glutf90.h"
#else
#if defined(MACOSX)
#include <GLUT/glutf90.h>
#else
#include <GL/glutf90.h>
#endif
#endif

/* environments that don't want to use APIENTRY */

#ifndef WIN32
#define APIENTRY
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

/* used when Fortran must pass an assumed shape array */
#ifdef F_LANG
#define ARRAYREF [0]
#define ARRAY *
#else
#define ARRAYREF
#define ARRAY
#endif

#ifdef F_LANG
#define PURE_FUNCTIONS
#endif

/* The patches for GLUT 3.6 used a different name for get/set FCB */
#if (TMP_GLUT_VERSION == 36 && FREEGLUT<1)
#define __glutGetFCB glut_Get_FCB
#define __glutSetFCB glut_Set_FCB
#endif

#ifdef FREEGLUT
GLUTAPI void* APIENTRY __glutGetTimerFCB(int which, int num, int *value);
GLUTAPI void APIENTRY __glutSetTimerFCB(int which, void *func, int num, int *value);
#endif

/* The patches for GLUT 3.6 used character strings instead of constants for
   the "which" argument to __glut[Get|Set]FCB */

#if (TMP_GLUT_VERSION == 36 && FREEGLUT<1)
/* for GLUT 3.6 */
#define GLUT_TMP_DISPLAY         "display"
#define GLUT_TMP_RESHAPE         "reshape"
#define GLUT_TMP_MOUSE           "mouse"
#define GLUT_TMP_MOTION          "motion"
#define GLUT_TMP_PASSIVE         "passive"
#define GLUT_TMP_ENTRY           "entry"
#define GLUT_TMP_KEYBOARD        "keyboard"
#define GLUT_TMP_KEYBOARD_UP     "keyboardup"
#define GLUT_TMP_WINDOW_STATUS   "windowStatus"
#define GLUT_TMP_VISIBILITY      "visibility"
#define GLUT_TMP_SPECIAL         "special"
#define GLUT_TMP_SPECIAL_UP      "specialup"
#define GLUT_TMP_BUTTON_BOX      "buttonBox"
#define GLUT_TMP_DIALS           "dials"
#define GLUT_TMP_SPACE_MOTION    "spaceMotion"
#define GLUT_TMP_SPACE_ROTATE    "spaceRotate"
#define GLUT_TMP_SPACE_BUTTON    "spaceButton"
#define GLUT_TMP_TABLET_MOTION   "tabletMotion"
#define GLUT_TMP_TABLET_BUTTON   "tabletButton"
#define GLUT_TMP_JOYSTICK        "joystick"
#define GLUT_TMP_OVERLAY_DISPLAY "overlay"
#define GLUT_TMP_SELECT          "select"
#define GLUT_TMP_TIMER           "timer"

#else
/* for GLUT 3.7 and later */

#define GLUT_TMP_DISPLAY         GLUT_FCB_DISPLAY
#define GLUT_TMP_RESHAPE         GLUT_FCB_RESHAPE
#define GLUT_TMP_MOUSE           GLUT_FCB_MOUSE
#define GLUT_TMP_MOTION          GLUT_FCB_MOTION
#define GLUT_TMP_PASSIVE         GLUT_FCB_PASSIVE
#define GLUT_TMP_ENTRY           GLUT_FCB_ENTRY
#define GLUT_TMP_KEYBOARD        GLUT_FCB_KEYBOARD
#define GLUT_TMP_KEYBOARD_UP     GLUT_FCB_KEYBOARD_UP
#define GLUT_TMP_WINDOW_STATUS   GLUT_FCB_WINDOW_STATUS
#define GLUT_TMP_VISIBILITY      GLUT_FCB_VISIBILITY
#define GLUT_TMP_SPECIAL         GLUT_FCB_SPECIAL
#define GLUT_TMP_SPECIAL_UP      GLUT_FCB_SPECIAL_UP
#define GLUT_TMP_BUTTON_BOX      GLUT_FCB_BUTTON_BOX
#define GLUT_TMP_DIALS           GLUT_FCB_DIALS
#define GLUT_TMP_SPACE_MOTION    GLUT_FCB_SPACE_MOTION
#define GLUT_TMP_SPACE_ROTATE    GLUT_FCB_SPACE_ROTATE
#define GLUT_TMP_SPACE_BUTTON    GLUT_FCB_SPACE_BUTTON
#define GLUT_TMP_TABLET_MOTION   GLUT_FCB_TABLET_MOTION
#define GLUT_TMP_TABLET_BUTTON   GLUT_FCB_TABLET_BUTTON
#define GLUT_TMP_JOYSTICK        GLUT_FCB_JOYSTICK
#define GLUT_TMP_OVERLAY_DISPLAY GLUT_FCB_OVERLAY_DISPLAY
#define GLUT_TMP_SELECT          GLUT_FCB_SELECT
#define GLUT_TMP_TIMER           GLUT_FCB_TIMER

#endif

/* GLUT Fortran callback function types from wrapper to user code */

#ifdef WIN32
#ifdef NAS
#define GLUTCALLBACKFCB __stdcall
#else
#ifdef DVF
#define GLUTCALLBACKFCB __stdcall
#else
#define GLUTCALLBACKFCB __cdecl
#endif
#endif
#else
#define GLUTCALLBACKFCB
#endif

typedef void (GLUTCALLBACKFCB *GLUTdisplayFCBUSR) (void);
typedef void (GLUTCALLBACKFCB *GLUTreshapeFCBUSR) (int *, int *);
/* NOTE the pressed key is int, not unsigned char for Fortran! */
typedef void (GLUTCALLBACKFCB *GLUTkeyboardFCBUSR) (int *, int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTmouseFCBUSR) (int *, int *, int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTmotionFCBUSR) (int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTpassiveFCBUSR) (int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTentryFCBUSR) (int *);
typedef void (GLUTCALLBACKFCB *GLUTwindowStatusFCBUSR) (int *);
typedef void (GLUTCALLBACKFCB *GLUTvisibilityFCBUSR) (int *);
typedef void (GLUTCALLBACKFCB *GLUTspecialFCBUSR) (int *, int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTbuttonBoxFCBUSR) (int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTdialsFCBUSR) (int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTspaceMotionFCBUSR) (int *, int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTspaceRotateFCBUSR) (int *, int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTspaceButtonFCBUSR) (int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTtabletMotionFCBUSR) (int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTtabletButtonFCBUSR) (int *, int *, int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTjoystickFCBUSR) (unsigned int *buttonMask, int *x, int *y, int *z);
typedef void (GLUTCALLBACKFCB *GLUTselectFCBUSR) (int *);
typedef void (GLUTCALLBACKFCB *GLUTtimerFCBUSR) (int *);
typedef void (GLUTCALLBACKFCB *GLUTmenuStateFCBUSR) (int *);  /* DEPRICATED. */
typedef void (GLUTCALLBACKFCB *GLUTmenuStatusFCBUSR) (int *, int *, int *);
typedef void (GLUTCALLBACKFCB *GLUTidleFCBUSR) (void);

/* Fortran name mangling choices */
#define LOWERCASE 1
#define UNDERSCORE 2
#define UPPERCASE 3

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

/* procedure name mangling */

#if FNAME==LOWERCASE
#define myglutCompareFuncs             myglutcomparefuncs
#define f9yglutAddMenuEntry            f9yglutaddmenuentry
#define f9yglutAddSubMenu              f9yglutaddsubmenu
#define fglutAttachMenu                fglutattachmenu
#define fglutBitmapCharacter           fglutbitmapcharacter
#define fglutBitmapLength              fglutbitmaplength
#define fglutBitmapWidth               fglutbitmapwidth
#define fglutButtonBoxFunc             fglutbuttonboxfunc
#define f9y1glutButtonBoxFunc          f9y1glutbuttonboxfunc
#define f9yglutChangeToMenuEntry       f9yglutchangetomenuentry
#define f9yglutChangeToSubMenu         f9yglutchangetosubmenu
#define fglutCopyColormap              fglutcopycolormap
#define fglutCreateMenu                fglutcreatemenu
#define fglutCreateSubWindow           fglutcreatesubwindow
#define f9yglutCreateWindow            f9yglutcreatewindow
#define fglutDestroyMenu               fglutdestroymenu
#define fglutDestroyWindow             fglutdestroywindow
#define fglutDetachMenu                fglutdetachmenu
#define fglutDeviceGet                 fglutdeviceget
#define fglutDialsFunc                 fglutdialsfunc
#define f9y1glutDialsFunc              f9y1glutdialsfunc
#define fglutDisplayFunc               fglutdisplayfunc
#define fglutEntryFunc                 fglutentryfunc
#define f9y1glutEntryFunc              f9y1glutentryfunc
#define fglutEstablishOverlay          fglutestablishoverlay
#define f9yglutExtensionSupported      f9yglutextensionsupported
#define fglutFullScreen                fglutfullscreen
#define fglutGet                       fglutget
#define fglutGetColor                  fglutgetcolor
#define fglutGetMenu                   fglutgetmenu
#define fglutGetModifiers              fglutgetmodifiers
#define fglutGetWindow                 fglutgetwindow
#define fglutHideOverlay               fgluthideoverlay
#define fglutHideWindow                fgluthidewindow
#define fglutIconifyWindow             fgluticonifywindow
#define fglutIdleFunc                  fglutidlefunc
#define f9y1glutIdleFunc               f9y1glutidlefunc
#define fglutInit                      fglutinit
#define fglutInitArgs                  fglutinitargs
#define fglutInitDisplayMode           fglutinitdisplaymode
#define f9yglutInitDisplayString       f9yglutinitdisplaystring
#define fglutInitWindowPosition        fglutinitwindowposition
#define fglutInitWindowSize            fglutinitwindowsize
#define fglutKeyboardFunc              fglutkeyboardfunc
#define f9y1glutKeyboardFunc           f9y1glutkeyboardfunc
#define fglutLayerGet                  fglutlayerget
#define fglutMainLoop                  fglutmainloop
#define fglutMenuStateFunc             fglutmenustatefunc
#define f9y1glutMenuStateFunc          f9y1glutmenustatefunc
#define fglutMenuStatusFunc            fglutmenustatusfunc
#define f9y1glutMenuStatusFunc         f9y1glutmenustatusfunc
#define fglutMotionFunc                fglutmotionfunc
#define f9y1glutMotionFunc             f9y1glutmotionfunc
#define fglutMouseFunc                 fglutmousefunc
#define f9y1glutMouseFunc              f9y1glutmousefunc
#define fglutOverlayDisplayFunc        fglutoverlaydisplayfunc
#define f9y1glutOverlayDisplayFunc     f9y1glutoverlaydisplayfunc
#define fglutPassiveMotionFunc         fglutpassivemotionfunc
#define f9y1glutPassiveMotionFunc      f9y1glutpassivemotionfunc
#define fglutPopWindow                 fglutpopwindow
#define fglutPositionWindow            fglutpositionwindow
#define fglutPostOverlayRedisplay      fglutpostoverlayredisplay
#define fglutPostRedisplay             fglutpostredisplay
#define fglutPostWindowOverlayRedisplay fglutpostwindowoverlayredisplay
#define fglutPostWindowRedisplay       fglutpostwindowredisplay
#define fglutPushWindow                fglutpushwindow
#define fglutRemoveMenuItem            fglutremovemenuitem
#define fglutRemoveOverlay             fglutremoveoverlay
#define fglutReportErrors              fglutreporterrors
#define fglutReshapeFunc               fglutreshapefunc
#define f9y1glutReshapeFunc            f9y1glutreshapefunc
#define fglutReshapeWindow             fglutreshapewindow
#define fglutSetColor                  fglutsetcolor
#define fglutSetCursor                 fglutsetcursor
#define f9yglutSetIconTitle            f9yglutseticontitle
#define fglutSetMenu                   fglutsetmenu
#define fglutSetWindow                 fglutsetwindow
#define f9yglutSetWindowTitle          f9yglutsetwindowtitle
#define fglutSetupVideoResizing        fglutsetupvideoresizing
#define fglutShowOverlay               fglutshowoverlay
#define fglutShowWindow                fglutshowwindow
#define fglutSolidCone                 fglutsolidcone
#define fglutSolidCube                 fglutsolidcube
#define fglutSolidDodecahedron         fglutsoliddodecahedron
#define fglutSolidIcosahedron          fglutsolidicosahedron
#define fglutSolidOctahedron           fglutsolidoctahedron
#define fglutSolidSphere               fglutsolidsphere
#define fglutSolidTeapot               fglutsolidteapot
#define fglutSolidTetrahedron          fglutsolidtetrahedron
#define fglutSolidTorus                fglutsolidtorus
#define fglutSpaceballButtonFunc       fglutspaceballbuttonfunc
#define f9y1glutSpaceballButtonFunc    f9y1glutspaceballbuttonfunc
#define fglutSpaceballMotionFunc       fglutspaceballmotionfunc
#define f9y1glutSpaceballMotionFunc    f9y1glutspaceballmotionfunc
#define fglutSpaceballRotateFunc       fglutspaceballrotatefunc
#define f9y1glutSpaceballRotateFunc    f9y1glutspaceballrotatefunc
#define fglutSpecialFunc               fglutspecialfunc
#define f9y1glutSpecialFunc            f9y1glutspecialfunc
#define fglutStopVideoResizing         fglutstopvideoresizing
#define fglutStrokeCharacter           fglutstrokecharacter
#define fglutStrokeLength              fglutstrokelength
#define fglutStrokeWidth               fglutstrokewidth
#define fglutSwapBuffers               fglutswapbuffers
#define fglutTabletButtonFunc          fgluttabletbuttonfunc
#define f9y1glutTabletButtonFunc       f9y1gluttabletbuttonfunc
#define fglutTabletMotionFunc          fgluttabletmotionfunc
#define f9y1glutTabletMotionFunc       f9y1gluttabletmotionfunc
#define fglutTimerFunc                 fgluttimerfunc
#define fglutUseLayer                  fglutuselayer
#define fglutVideoPan                  fglutvideopan
#define fglutVideoResize               fglutvideoresize
#define fglutVideoResizeGet            fglutvideoresizeget
#define fglutVisibilityFunc            fglutvisibilityfunc
#define f9y1glutVisibilityFunc         f9y1glutvisibilityfunc
#define fglutWarpPointer               fglutwarppointer
#define fglutWindowStatusFunc          fglutwindowstatusfunc
#define f9y1glutWindowStatusFunc       f9y1glutwindowstatusfunc
#define fglutWireCone                  fglutwirecone
#define fglutWireCube                  fglutwirecube
#define fglutWireDodecahedron          fglutwiredodecahedron
#define fglutWireIcosahedron           fglutwireicosahedron
#define fglutWireOctahedron            fglutwireoctahedron
#define fglutWireSphere                fglutwiresphere
#define fglutWireTeapot                fglutwireteapot
#define fglutWireTetrahedron           fglutwiretetrahedron
#define fglutWireTorus                 fglutwiretorus
#ifdef GLUT37
#define fglutKeyboardUpFunc            fglutkeyboardupfunc
#define f9y1glutKeyboardUpFunc         f9y1glutkeyboardupfunc
#define fglutSpecialUpFunc             fglutspecialupfunc
#define f9y1glutSpecialUpFunc          f9y1glutspecialupfunc
#define fglutJoystickFunc              fglutjoystickfunc
#define f9y1glutJoystickFunc           f9y1glutjoystickfunc
#define fglutIgnoreKeyRepeat           fglutignorekeyrepeat
#define fglutSetKeyRepeat              fglutsetkeyrepeat
#define fglutForceJoystickFunc         fglutforcejoystickfunc
#define f9yglutGameModeString          f9yglutgamemodestring
#define fglutEnterGameMode             fglutentergamemode
#define fglutLeaveGameMode             fglutleavegamemode
#define fglutGameModeGet               fglutgamemodeget
#ifdef FREEGLUT
#define fglutLeaveMainLoop             fglutleavemainloop
#define fglutSetOption                 fglutsetoption
#endif
#endif
#elif FNAME==UNDERSCORE
#define myglutCompareFuncs             myglutcomparefuncs_
#define f9yglutAddMenuEntry            f9yglutaddmenuentry_
#define f9yglutAddSubMenu              f9yglutaddsubmenu_
#define fglutAttachMenu                fglutattachmenu_
#define fglutBitmapCharacter           fglutbitmapcharacter_
#define fglutBitmapLength              fglutbitmaplength_
#define fglutBitmapWidth               fglutbitmapwidth_
#define fglutButtonBoxFunc             fglutbuttonboxfunc_
#define f9y1glutButtonBoxFunc          f9y1glutbuttonboxfunc_
#define f9yglutChangeToMenuEntry       f9yglutchangetomenuentry_
#define f9yglutChangeToSubMenu         f9yglutchangetosubmenu_
#define fglutCopyColormap              fglutcopycolormap_
#define fglutCreateMenu                fglutcreatemenu_
#define fglutCreateSubWindow           fglutcreatesubwindow_
#define f9yglutCreateWindow            f9yglutcreatewindow_
#define fglutDestroyMenu               fglutdestroymenu_
#define fglutDestroyWindow             fglutdestroywindow_
#define fglutDetachMenu                fglutdetachmenu_
#define fglutDeviceGet                 fglutdeviceget_
#define fglutDialsFunc                 fglutdialsfunc_
#define f9y1glutDialsFunc              f9y1glutdialsfunc_
#define fglutDisplayFunc               fglutdisplayfunc_
#define fglutEntryFunc                 fglutentryfunc_
#define f9y1glutEntryFunc              f9y1glutentryfunc_
#define fglutEstablishOverlay          fglutestablishoverlay_
#define f9yglutExtensionSupported      f9yglutextensionsupported_
#define fglutFullScreen                fglutfullscreen_
#define fglutGet                       fglutget_
#define fglutGetColor                  fglutgetcolor_
#define fglutGetMenu                   fglutgetmenu_
#define fglutGetModifiers              fglutgetmodifiers_
#define fglutGetWindow                 fglutgetwindow_
#define fglutHideOverlay               fgluthideoverlay_
#define fglutHideWindow                fgluthidewindow_
#define fglutIconifyWindow             fgluticonifywindow_
#define fglutIdleFunc                  fglutidlefunc_
#define f9y1glutIdleFunc               f9y1glutidlefunc_
#define fglutInit                      fglutinit_
#define fglutInitArgs                  fglutinitargs_
#define fglutInitDisplayMode           fglutinitdisplaymode_
#define f9yglutInitDisplayString       f9yglutinitdisplaystring_
#define fglutInitWindowPosition        fglutinitwindowposition_
#define fglutInitWindowSize            fglutinitwindowsize_
#define fglutKeyboardFunc              fglutkeyboardfunc_
#define f9y1glutKeyboardFunc           f9y1glutkeyboardfunc_
#define fglutLayerGet                  fglutlayerget_
#define fglutMainLoop                  fglutmainloop_
#define fglutMenuStateFunc             fglutmenustatefunc_
#define f9y1glutMenuStateFunc          f9y1glutmenustatefunc_
#define fglutMenuStatusFunc            fglutmenustatusfunc_
#define f9y1glutMenuStatusFunc         f9y1glutmenustatusfunc_
#define fglutMotionFunc                fglutmotionfunc_
#define f9y1glutMotionFunc             f9y1glutmotionfunc_
#define fglutMouseFunc                 fglutmousefunc_
#define f9y1glutMouseFunc              f9y1glutmousefunc_
#define fglutOverlayDisplayFunc        fglutoverlaydisplayfunc_
#define f9y1glutOverlayDisplayFunc     f9y1glutoverlaydisplayfunc_
#define fglutPassiveMotionFunc         fglutpassivemotionfunc_
#define f9y1glutPassiveMotionFunc      f9y1glutpassivemotionfunc_
#define fglutPopWindow                 fglutpopwindow_
#define fglutPositionWindow            fglutpositionwindow_
#define fglutPostOverlayRedisplay      fglutpostoverlayredisplay_
#define fglutPostRedisplay             fglutpostredisplay_
#define fglutPostWindowOverlayRedisplay fglutpostwindowoverlayredisplay_
#define fglutPostWindowRedisplay       fglutpostwindowredisplay_
#define fglutPushWindow                fglutpushwindow_
#define fglutRemoveMenuItem            fglutremovemenuitem_
#define fglutRemoveOverlay             fglutremoveoverlay_
#define fglutReportErrors              fglutreporterrors_
#define fglutReshapeFunc               fglutreshapefunc_
#define f9y1glutReshapeFunc            f9y1glutreshapefunc_
#define fglutReshapeWindow             fglutreshapewindow_
#define fglutSetColor                  fglutsetcolor_
#define fglutSetCursor                 fglutsetcursor_
#define f9yglutSetIconTitle            f9yglutseticontitle_
#define fglutSetMenu                   fglutsetmenu_
#define fglutSetWindow                 fglutsetwindow_
#define f9yglutSetWindowTitle          f9yglutsetwindowtitle_
#define fglutSetupVideoResizing        fglutsetupvideoresizing_
#define fglutShowOverlay               fglutshowoverlay_
#define fglutShowWindow                fglutshowwindow_
#define fglutSolidCone                 fglutsolidcone_
#define fglutSolidCube                 fglutsolidcube_
#define fglutSolidDodecahedron         fglutsoliddodecahedron_
#define fglutSolidIcosahedron          fglutsolidicosahedron_
#define fglutSolidOctahedron           fglutsolidoctahedron_
#define fglutSolidSphere               fglutsolidsphere_
#define fglutSolidTeapot               fglutsolidteapot_
#define fglutSolidTetrahedron          fglutsolidtetrahedron_
#define fglutSolidTorus                fglutsolidtorus_
#define fglutSpaceballButtonFunc       fglutspaceballbuttonfunc_
#define f9y1glutSpaceballButtonFunc    f9y1glutspaceballbuttonfunc_
#define fglutSpaceballMotionFunc       fglutspaceballmotionfunc_
#define f9y1glutSpaceballMotionFunc    f9y1glutspaceballmotionfunc_
#define fglutSpaceballRotateFunc       fglutspaceballrotatefunc_
#define f9y1glutSpaceballRotateFunc    f9y1glutspaceballrotatefunc_
#define fglutSpecialFunc               fglutspecialfunc_
#define f9y1glutSpecialFunc            f9y1glutspecialfunc_
#define fglutStopVideoResizing         fglutstopvideoresizing_
#define fglutStrokeCharacter           fglutstrokecharacter_
#define fglutStrokeLength              fglutstrokelength_
#define fglutStrokeWidth               fglutstrokewidth_
#define fglutSwapBuffers               fglutswapbuffers_
#define fglutTabletButtonFunc          fgluttabletbuttonfunc_
#define f9y1glutTabletButtonFunc       f9y1gluttabletbuttonfunc_
#define fglutTabletMotionFunc          fgluttabletmotionfunc_
#define f9y1glutTabletMotionFunc       f9y1gluttabletmotionfunc_
#define fglutTimerFunc                 fgluttimerfunc_
#define fglutUseLayer                  fglutuselayer_
#define fglutVideoPan                  fglutvideopan_
#define fglutVideoResize               fglutvideoresize_
#define fglutVideoResizeGet            fglutvideoresizeget_
#define fglutVisibilityFunc            fglutvisibilityfunc_
#define f9y1glutVisibilityFunc         f9y1glutvisibilityfunc_
#define fglutWarpPointer               fglutwarppointer_
#define fglutWindowStatusFunc          fglutwindowstatusfunc_
#define f9y1glutWindowStatusFunc       f9y1glutwindowstatusfunc_
#define fglutWireCone                  fglutwirecone_
#define fglutWireCube                  fglutwirecube_
#define fglutWireDodecahedron          fglutwiredodecahedron_
#define fglutWireIcosahedron           fglutwireicosahedron_
#define fglutWireOctahedron            fglutwireoctahedron_
#define fglutWireSphere                fglutwiresphere_
#define fglutWireTeapot                fglutwireteapot_
#define fglutWireTetrahedron           fglutwiretetrahedron_
#define fglutWireTorus                 fglutwiretorus_
#ifdef GLUT37
#define fglutKeyboardUpFunc            fglutkeyboardupfunc_
#define f9y1glutKeyboardUpFunc         f9y1glutkeyboardupfunc_
#define fglutSpecialUpFunc             fglutspecialupfunc_
#define f9y1glutSpecialUpFunc          f9y1glutspecialupfunc_
#define fglutJoystickFunc              fglutjoystickfunc_
#define f9y1glutJoystickFunc           f9y1glutjoystickfunc_
#define fglutIgnoreKeyRepeat           fglutignorekeyrepeat_
#define fglutSetKeyRepeat              fglutsetkeyrepeat_
#define fglutForceJoystickFunc         fglutforcejoystickfunc_
#define f9yglutGameModeString          f9yglutgamemodestring_
#define fglutEnterGameMode             fglutentergamemode_
#define fglutLeaveGameMode             fglutleavegamemode_
#define fglutGameModeGet               fglutgamemodeget_
#ifdef FREEGLUT
#define fglutLeaveMainLoop             fglutleavemainloop_
#define fglutSetOption                 fglutsetoption_
#endif
#endif
#elif FNAME==UPPERCASE
#define myglutCompareFuncs             MYGLUTCOMPAREFUNCS
#define f9yglutAddMenuEntry            F9YGLUTADDMENUENTRY
#define f9yglutAddSubMenu              F9YGLUTADDSUBMENU
#define fglutAttachMenu                FGLUTATTACHMENU
#define fglutBitmapCharacter           FGLUTBITMAPCHARACTER
#define fglutBitmapLength              FGLUTBITMAPLENGTH
#define fglutBitmapWidth               FGLUTBITMAPWIDTH
#define fglutButtonBoxFunc             FGLUTBUTTONBOXFUNC
#define f9y1glutButtonBoxFunc          F9Y1GLUTBUTTONBOXFUNC
#define f9yglutChangeToMenuEntry       F9YGLUTCHANGETOMENUENTRY
#define f9yglutChangeToSubMenu         F9YGLUTCHANGETOSUBMENU
#define fglutCopyColormap              FGLUTCOPYCOLORMAP
#define fglutCreateMenu                FGLUTCREATEMENU
#define fglutCreateSubWindow           FGLUTCREATESUBWINDOW
#define f9yglutCreateWindow            F9YGLUTCREATEWINDOW
#define fglutDestroyMenu               FGLUTDESTROYMENU
#define fglutDestroyWindow             FGLUTDESTROYWINDOW
#define fglutDetachMenu                FGLUTDETACHMENU
#define fglutDeviceGet                 FGLUTDEVICEGET
#define fglutDialsFunc                 FGLUTDIALSFUNC
#define f9y1glutDialsFunc              F9Y1GLUTDIALSFUNC
#define fglutDisplayFunc               FGLUTDISPLAYFUNC
#define fglutEntryFunc                 FGLUTENTRYFUNC
#define f9y1glutEntryFunc              F9Y1GLUTENTRYFUNC
#define fglutEstablishOverlay          FGLUTESTABLISHOVERLAY
#define f9yglutExtensionSupported      F9YGLUTEXTENSIONSUPPORTED
#define fglutFullScreen                FGLUTFULLSCREEN
#define fglutGet                       FGLUTGET
#define fglutGetColor                  FGLUTGETCOLOR
#define fglutGetMenu                   FGLUTGETMENU
#define fglutGetModifiers              FGLUTGETMODIFIERS
#define fglutGetWindow                 FGLUTGETWINDOW
#define fglutHideOverlay               FGLUTHIDEOVERLAY
#define fglutHideWindow                FGLUTHIDEWINDOW
#define fglutIconifyWindow             FGLUTICONIFYWINDOW
#define fglutIdleFunc                  FGLUTIDLEFUNC
#define f9y1glutIdleFunc               F9Y1GLUTIDLEFUNC
#define fglutInit                      FGLUTINIT
#define fglutInitArgs                  FGLUTINITARGS
#define fglutInitDisplayMode           FGLUTINITDISPLAYMODE
#define f9yglutInitDisplayString       F9YGLUTINITDISPLAYSTRING
#define fglutInitWindowPosition        FGLUTINITWINDOWPOSITION
#define fglutInitWindowSize            FGLUTINITWINDOWSIZE
#define fglutKeyboardFunc              FGLUTKEYBOARDFUNC
#define f9y1glutKeyboardFunc           F9Y1GLUTKEYBOARDFUNC
#define fglutLayerGet                  FGLUTLAYERGET
#define fglutMainLoop                  FGLUTMAINLOOP
#define fglutMenuStateFunc             FGLUTMENUSTATEFUNC
#define f9y1glutMenuStateFunc          F9Y1GLUTMENUSTATEFUNC
#define fglutMenuStatusFunc            FGLUTMENUSTATUSFUNC
#define f9y1glutMenuStatusFunc         F9Y1GLUTMENUSTATUSFUNC
#define fglutMotionFunc                FGLUTMOTIONFUNC
#define f9y1glutMotionFunc             F9Y1GLUTMOTIONFUNC
#define fglutMouseFunc                 FGLUTMOUSEFUNC
#define f9y1glutMouseFunc              F9Y1GLUTMOUSEFUNC
#define fglutOverlayDisplayFunc        FGLUTOVERLAYDISPLAYFUNC
#define f9y1glutOverlayDisplayFunc     F9Y1GLUTOVERLAYDISPLAYFUNC
#define fglutPassiveMotionFunc         FGLUTPASSIVEMOTIONFUNC
#define f9y1glutPassiveMotionFunc      F9Y1GLUTPASSIVEMOTIONFUNC
#define fglutPopWindow                 FGLUTPOPWINDOW
#define fglutPositionWindow            FGLUTPOSITIONWINDOW
#define fglutPostOverlayRedisplay      FGLUTPOSTOVERLAYREDISPLAY
#define fglutPostRedisplay             FGLUTPOSTREDISPLAY
#define fglutPostWindowOverlayRedisplay FGLUTPOSTWINDOWOVERLAYREDISPLAY
#define fglutPostWindowRedisplay       FGLUTPOSTWINDOWREDISPLAY
#define fglutPushWindow                FGLUTPUSHWINDOW
#define fglutRemoveMenuItem            FGLUTREMOVEMENUITEM
#define fglutRemoveOverlay             FGLUTREMOVEOVERLAY
#define fglutReportErrors              FGLUTREPORTERRORS
#define fglutReshapeFunc               FGLUTRESHAPEFUNC
#define f9y1glutReshapeFunc            F9Y1GLUTRESHAPEFUNC
#define fglutReshapeWindow             FGLUTRESHAPEWINDOW
#define fglutSetColor                  FGLUTSETCOLOR
#define fglutSetCursor                 FGLUTSETCURSOR
#define f9yglutSetIconTitle            F9YGLUTSETICONTITLE
#define fglutSetMenu                   FGLUTSETMENU
#define fglutSetWindow                 FGLUTSETWINDOW
#define f9yglutSetWindowTitle          F9YGLUTSETWINDOWTITLE
#define fglutSetupVideoResizing        FGLUTSETUPVIDEORESIZING
#define fglutShowOverlay               FGLUTSHOWOVERLAY
#define fglutShowWindow                FGLUTSHOWWINDOW
#define fglutSolidCone                 FGLUTSOLIDCONE
#define fglutSolidCube                 FGLUTSOLIDCUBE
#define fglutSolidDodecahedron         FGLUTSOLIDDODECAHEDRON
#define fglutSolidIcosahedron          FGLUTSOLIDICOSAHEDRON
#define fglutSolidOctahedron           FGLUTSOLIDOCTAHEDRON
#define fglutSolidSphere               FGLUTSOLIDSPHERE
#define fglutSolidTeapot               FGLUTSOLIDTEAPOT
#define fglutSolidTetrahedron          FGLUTSOLIDTETRAHEDRON
#define fglutSolidTorus                FGLUTSOLIDTORUS
#define fglutSpaceballButtonFunc       FGLUTSPACEBALLBUTTONFUNC
#define f9y1glutSpaceballButtonFunc    F9Y1GLUTSPACEBALLBUTTONFUNC
#define fglutSpaceballMotionFunc       FGLUTSPACEBALLMOTIONFUNC
#define f9y1glutSpaceballMotionFunc    F9Y1GLUTSPACEBALLMOTIONFUNC
#define fglutSpaceballRotateFunc       FGLUTSPACEBALLROTATEFUNC
#define f9y1glutSpaceballRotateFunc    F9Y1GLUTSPACEBALLROTATEFUNC
#define fglutSpecialFunc               FGLUTSPECIALFUNC
#define f9y1glutSpecialFunc            F9Y1GLUTSPECIALFUNC
#define fglutStopVideoResizing         FGLUTSTOPVIDEORESIZING
#define fglutStrokeCharacter           FGLUTSTROKECHARACTER
#define fglutStrokeLength              FGLUTSTROKELENGTH
#define fglutStrokeWidth               FGLUTSTROKEWIDTH
#define fglutSwapBuffers               FGLUTSWAPBUFFERS
#define fglutTabletButtonFunc          FGLUTTABLETBUTTONFUNC
#define f9y1glutTabletButtonFunc       F9Y1GLUTTABLETBUTTONFUNC
#define fglutTabletMotionFunc          FGLUTTABLETMOTIONFUNC
#define f9y1glutTabletMotionFunc       F9Y1GLUTTABLETMOTIONFUNC
#define fglutTimerFunc                 FGLUTTIMERFUNC
#define fglutUseLayer                  FGLUTUSELAYER
#define fglutVideoPan                  FGLUTVIDEOPAN
#define fglutVideoResize               FGLUTVIDEORESIZE
#define fglutVideoResizeGet            FGLUTVIDEORESIZEGET
#define fglutVisibilityFunc            FGLUTVISIBILITYFUNC
#define f9y1glutVisibilityFunc         F9Y1GLUTVISIBILITYFUNC
#define fglutWarpPointer               FGLUTWARPPOINTER
#define fglutWindowStatusFunc          FGLUTWINDOWSTATUSFUNC
#define f9y1glutWindowStatusFunc       F9Y1GLUTWINDOWSTATUSFUNC
#define fglutWireCone                  FGLUTWIRECONE
#define fglutWireCube                  FGLUTWIRECUBE
#define fglutWireDodecahedron          FGLUTWIREDODECAHEDRON
#define fglutWireIcosahedron           FGLUTWIREICOSAHEDRON
#define fglutWireOctahedron            FGLUTWIREOCTAHEDRON
#define fglutWireSphere                FGLUTWIRESPHERE
#define fglutWireTeapot                FGLUTWIRETEAPOT
#define fglutWireTetrahedron           FGLUTWIRETETRAHEDRON
#define fglutWireTorus                 FGLUTWIRETORUS
#ifdef GLUT37
#define fglutKeyboardUpFunc            FGLUTKEYBOARDUPFUNC
#define f9y1glutKeyboardUpFunc         F9Y1GLUTKEYBOARDUPFUNC
#define fglutSpecialUpFunc             FGLUTSPECIALUPFUNC
#define f9y1glutSpecialUpFunc          F9Y1GLUTSPECIALUPFUNC
#define fglutJoystickFunc              FGLUTJOYSTICKFUNC
#define f9y1glutJoystickFunc           F9Y1GLUTJOYSTICKFUNC
#define fglutIgnoreKeyRepeat           FGLUTIGNOREKEYREPEAT
#define fglutSetKeyRepeat              FGLUTSETKEYREPEAT
#define fglutForceJoystickFunc         FGLUTFORCEJOYSTICKFUNC
#define f9yglutGameModeString          F9YGLUTGAMEMODESTRING
#define fglutEnterGameMode             FGLUTENTERGAMEMODE
#define fglutLeaveGameMode             FGLUTLEAVEGAMEMODE
#define fglutGameModeGet               FGLUTGAMEMODEGET
#endif
#ifdef FREEGLUT
#define fglutLeaveMainLoop             FGLUTLEAVEMAINLOOP
#define fglutSetOption                 FGLUTSETOPTION
#endif
#endif

/* fonts are passed in as this integer */
#define F90GLUT_STROKE_ROMAN               1
#define F90GLUT_STROKE_MONO_ROMAN          2
#define F90GLUT_BITMAP_9_BY_15             3
#define F90GLUT_BITMAP_8_BY_13             4
#define F90GLUT_BITMAP_TIMES_ROMAN_10      5
#define F90GLUT_BITMAP_TIMES_ROMAN_24      6
#define F90GLUT_BITMAP_HELVETICA_10        7
#define F90GLUT_BITMAP_HELVETICA_12        8
#define F90GLUT_BITMAP_HELVETICA_18        9

/* places to put command line arguments, if provided */
int __myArgc;
char **__myArgv;

/* places to keep the idle, menuState and menuStatus callback functions */
GLUTidleFCBUSR fidle;
GLUTmenuStateFCBUSR fmenustate;
GLUTmenuStatusFCBUSR fmenustatus;

/* translate integer to font */
void * f90glut_get_font(int name)
{
switch (name) {
   case F90GLUT_STROKE_ROMAN:
      return GLUT_STROKE_ROMAN;
   case F90GLUT_STROKE_MONO_ROMAN:
      return GLUT_STROKE_MONO_ROMAN;
   case F90GLUT_BITMAP_9_BY_15:
      return GLUT_BITMAP_9_BY_15;
   case F90GLUT_BITMAP_8_BY_13:
      return GLUT_BITMAP_8_BY_13;
   case F90GLUT_BITMAP_TIMES_ROMAN_10:
      return GLUT_BITMAP_TIMES_ROMAN_10;
   case F90GLUT_BITMAP_TIMES_ROMAN_24:
      return GLUT_BITMAP_TIMES_ROMAN_24;
   case F90GLUT_BITMAP_HELVETICA_10:
      return GLUT_BITMAP_HELVETICA_10;
   case F90GLUT_BITMAP_HELVETICA_12:
      return GLUT_BITMAP_HELVETICA_12;
   case F90GLUT_BITMAP_HELVETICA_18:
      return GLUT_BITMAP_HELVETICA_18;
   }
return NULL;
}

/* compare two pointers for equality -- used to check for GLUTNULLFUNC */
GLint APIENTRY myglutCompareFuncs(void * ptr1, void * ptr2)
{
   if (ptr1==ptr2) {
      return 1;
   }
   else {
      return 0;
   }
}

/* GLUT 3.6 */

void APIENTRY f9yglutAddMenuEntry(f90glint ARRAY label, f90glint size_label,
   f90glcint value )
{
   char *char_label = (char *)malloc((*size_label+1)*sizeof(char));
   int i;
   for (i=0;i<*size_label;i++) char_label[i]=(char) label ARRAYREF [i];
   char_label[*size_label] = end_of_string;
   glutAddMenuEntry( char_label, (int) *value);
   free(char_label);
}

void APIENTRY f9yglutAddSubMenu(f90glint ARRAY label, f90glint size_label,
   f90glcint submenu )
{
   char *char_label = (char *)malloc((*size_label+1)*sizeof(char));
   int i;
   for (i=0;i<*size_label;i++) char_label[i]=(char) label ARRAYREF [i];
   char_label[*size_label] = end_of_string;
   glutAddSubMenu( char_label, (int) *submenu);
   free(char_label);
}

void APIENTRY fglutAttachMenu( f90glcint button )
{
   glutAttachMenu( (int) *button );
}

void APIENTRY fglutBitmapCharacter( f90glint font, f90glcint character )
{
   glutBitmapCharacter( f90glut_get_font(*font), (int) *character );
}

int APIENTRY fglutBitmapLength( f90glint font, f90glint ARRAY string,
   f90glint size_string )
{
   unsigned char *char_string = (unsigned char *)malloc((*size_string+1)*sizeof(unsigned char));
   int i;
   for (i=0;i<*size_string;i++) char_string[i]=(unsigned char) string ARRAYREF [i];
   char_string[*size_string] = end_of_string;
   i = glutBitmapLength( f90glut_get_font(*font), char_string );
   free(char_string);
   return i;
}

int APIENTRY fglutBitmapWidth( f90glint font, f90glcint character )
{
   return (int) glutBitmapWidth( f90glut_get_font(*font), (int) *character );
}

void GLUTCALLBACK fortranButtonBoxWrapper(int button, int state)
{
  ((GLUTbuttonBoxFCBUSR)__glutGetFCB(GLUT_TMP_BUTTON_BOX))(&button, &state);
}

void APIENTRY f9y1glutButtonBoxFunc()
{
    glutButtonBoxFunc(NULL);
}

void APIENTRY fglutButtonBoxFunc(GLUTbuttonBoxFCBUSR buttonbox)
{
    __glutSetFCB(GLUT_TMP_BUTTON_BOX,(void *)buttonbox);
    glutButtonBoxFunc(fortranButtonBoxWrapper);
}

void APIENTRY f9yglutChangeToMenuEntry(f90glcint item,f90glint ARRAY label,
    f90glint size_label,f90glcint value )
{
   char *char_label = (char *)malloc((*size_label+1)*sizeof(char));
   int i;
   for (i=0;i<*size_label;i++) char_label[i]=(char) label ARRAYREF [i];
   char_label[*size_label] = end_of_string;
   glutChangeToMenuEntry( (int) *item, char_label, (int) *value);
   free(char_label);
}

void APIENTRY f9yglutChangeToSubMenu(f90glcint item,f90glint ARRAY label,
    f90glint size_label,f90glcint submenu )
{
   char *char_label = (char *)malloc((*size_label+1)*sizeof(char));
   int i;
   for (i=0;i<*size_label;i++) char_label[i]=(char) label ARRAYREF [i];
   char_label[*size_label] = end_of_string;
   glutChangeToSubMenu( (int) *item, char_label, (int) *submenu);
   free(char_label);
}

void APIENTRY fglutCopyColormap( f90glcint win )
{
   glutCopyColormap( (int) *win );
}

void GLUTCALLBACK fortranSelectWrapper(int value)
{
  ((GLUTselectFCBUSR)__glutGetFCB(GLUT_TMP_SELECT))(&value);
}

#ifdef PURE_FUNCTIONS
void APIENTRY fglutCreateMenu(GLUTselectFCBUSR select, f90glint menuid)
#else
int APIENTRY fglutCreateMenu(GLUTselectFCBUSR select)
#endif
{
#ifndef PURE_FUNCTIONS
int *menuid;
int tempint;
menuid = &tempint;
#endif
  *menuid = (int) glutCreateMenu(fortranSelectWrapper);
  __glutSetFCB(GLUT_TMP_SELECT,(void *)select);
#ifndef PURE_FUNCTIONS
  return *menuid;
#endif
}

int APIENTRY fglutCreateSubWindow( f90glcint win, f90glcint x, f90glcint y,
    f90glcint width, f90glcint height )
{
   return (int) glutCreateSubWindow( (int) *win, (int) *x, (int) *y,
    (int) *width, (int) *height );
}

int APIENTRY f9yglutCreateWindow(f90glint ARRAY title, f90glint size_title )
{
   char *char_title = (char *)malloc((*size_title+1)*sizeof(char));
   int i, retval;
   for (i=0;i<*size_title;i++) char_title[i]=(char) title ARRAYREF [i];
   char_title[*size_title] = end_of_string;
   retval = glutCreateWindow( char_title);
#if !defined(MACOSX)
   free(char_title);
#endif
   return retval;
}

void APIENTRY fglutDestroyMenu( f90glcint menu )
{
   glutDestroyMenu( (int) *menu );
}

void APIENTRY fglutDestroyWindow( f90glcint win )
{
   glutDestroyWindow( (int) *win );
}

void APIENTRY fglutDetachMenu( f90glcint button )
{
   glutDetachMenu( (int) *button );
}

int APIENTRY fglutDeviceGet( f90glenum type )
{
   return (int) glutDeviceGet( (GLenum) *type );
}

void GLUTCALLBACK fortranDialsWrapper(int dial, int value)
{
  ((GLUTdialsFCBUSR)__glutGetFCB(GLUT_TMP_DIALS))(&dial, &value);
}

void APIENTRY f9y1glutDialsFunc()
{
    glutDialsFunc(NULL);
}

void APIENTRY fglutDialsFunc(GLUTdialsFCBUSR dials)
{
    __glutSetFCB(GLUT_TMP_DIALS,(void *)dials);
    glutDialsFunc(fortranDialsWrapper);
}

void GLUTCALLBACK fortranDisplayWrapper(void)
{
  ((GLUTdisplayFCBUSR)__glutGetFCB(GLUT_TMP_DISPLAY))();
}

void APIENTRY fglutDisplayFunc(GLUTdisplayFCBUSR display)
{
    __glutSetFCB(GLUT_TMP_DISPLAY,(void *)display);
    glutDisplayFunc(fortranDisplayWrapper);
}

void GLUTCALLBACK fortranEntryWrapper(int state)
{
  ((GLUTentryFCBUSR)__glutGetFCB(GLUT_TMP_ENTRY))(&state);
}

void APIENTRY f9y1glutEntryFunc()
{
    glutEntryFunc(NULL);
}

void APIENTRY fglutEntryFunc(GLUTentryFCBUSR entry)
{
    __glutSetFCB(GLUT_TMP_ENTRY,(void *)entry);
    glutEntryFunc(fortranEntryWrapper);
}

void APIENTRY fglutEstablishOverlay( void )
{
   glutEstablishOverlay();
}

int APIENTRY f9yglutExtensionSupported(f90glint ARRAY name, f90glint size_name )
{
   char *char_name = (char *)malloc((*size_name+1)*sizeof(char));
   int i, retval;
   for (i=0;i<*size_name;i++) char_name[i]=(char) name ARRAYREF [i];
   char_name[*size_name] = end_of_string;
   retval = glutExtensionSupported( char_name);
   free(char_name);
   return retval;
}

void APIENTRY fglutFullScreen( void )
{
   glutFullScreen();
}

int APIENTRY fglutGet( f90glenum type )
{
   return (int) glutGet( (GLenum) *type );
}

GLfloat APIENTRY fglutGetColor( f90glcint ndx, f90glcint component )
{
   return (GLfloat) glutGetColor( (int) *ndx, (int) *component );
}

int APIENTRY fglutGetMenu( void )
{
   return (int) glutGetMenu();
}

int APIENTRY fglutGetModifiers( void )
{
   return (int) glutGetModifiers();
}

int APIENTRY fglutGetWindow( void )
{
   return (int) glutGetWindow();
}

void APIENTRY fglutHideOverlay( void )
{
   glutHideOverlay();
}

void APIENTRY fglutHideWindow( void )
{
   glutHideWindow();
}

void APIENTRY fglutIconifyWindow( void )
{
   glutIconifyWindow();
}

void GLUTCALLBACK fortranIdleWrapper(void)
{
    (*fidle) ();
}

void APIENTRY f9y1glutIdleFunc()
{
    fidle = (GLUTidleFCBUSR) NULL;
    glutIdleFunc(NULL);
}

void APIENTRY fglutIdleFunc(GLUTidleFCBUSR idle)
{
    fidle = (GLUTidleFCBUSR) idle;
    glutIdleFunc(fortranIdleWrapper);
}

void APIENTRY fglutInit(void)
{
  __myArgc = 1;
  __myArgv = (char **) malloc((__myArgc+1)*sizeof(char));
  __myArgv[0] = "noname";
  __myArgv[1] = NULL;
  glutInit(&__myArgc, __myArgv);
}

void APIENTRY fglutInitArgs(f90glcint argc, f90glint ARRAY argv, f90glint ARRAY starts)
{
  int i, j;
  __myArgc = *argc;
  __myArgv = (char **) malloc((__myArgc+1)*sizeof(char *));
  for (i=0; i<__myArgc; i++) {
    __myArgv[i] = (char *) malloc((starts[i+1]-starts[i]+1)*sizeof(char));
    for (j=0; j<starts ARRAYREF [i+1]-starts ARRAYREF [i]; j++) {
      __myArgv[i][j] = (char) argv ARRAYREF [starts ARRAYREF [i]+j-1];
    }
    __myArgv[i][starts ARRAYREF [i+1]-starts ARRAYREF [i]] = end_of_string;
  }
  __myArgv[__myArgc] = NULL;
  
  glutInit(&__myArgc, __myArgv);
}

void APIENTRY fglutInitDisplayMode( f90glcuint mode )
{
   glutInitDisplayMode( (unsigned int) *mode );
}

void APIENTRY f9yglutInitDisplayString(f90glint ARRAY string,
    f90glint size_string )
{
   char *char_string = (char *)malloc((*size_string+1)*sizeof(char));
   int i;
   for (i=0;i<*size_string;i++) char_string[i]=(char) string ARRAYREF [i];
   char_string[*size_string] = end_of_string;
   glutInitDisplayString( char_string);
   free(char_string);
}

void APIENTRY fglutInitWindowPosition( f90glcint x, f90glcint y )
{
   glutInitWindowPosition( (int) *x, (int) *y );
}

void APIENTRY fglutInitWindowSize( f90glcint width, f90glcint height )
{
   glutInitWindowSize( (int) *width, (int) *height );
}

void GLUTCALLBACK fortranKeyboardWrapper(unsigned char key, int x, int y)
{
  int int_key;
  int_key = (int)key;
  ((GLUTkeyboardFCBUSR)__glutGetFCB(GLUT_TMP_KEYBOARD))(&int_key, &x, &y);
}

void APIENTRY f9y1glutKeyboardFunc()
{
    glutKeyboardFunc(NULL);
}

void APIENTRY fglutKeyboardFunc(GLUTkeyboardFCBUSR keyboard)
{
    __glutSetFCB(GLUT_TMP_KEYBOARD,(void *)keyboard);
    glutKeyboardFunc(fortranKeyboardWrapper);
}

int APIENTRY fglutLayerGet( f90glenum type )
{
   return (int) glutLayerGet( (GLenum) *type );
}

void APIENTRY fglutMainLoop( void )
{
   glutMainLoop();
}

#ifdef FREEGLUT
void APIENTRY fglutLeaveMainLoop( void )
{
   glutLeaveMainLoop();
}

void APIENTRY fglutSetOption( f90glenum option, f90glcint value )
{
   glutSetOption( (GLenum) *option, (int) *value);
}
#endif

void GLUTCALLBACK fortranMenuStateWrapper(int state)
{
    (*fmenustate)(&state);
}

void APIENTRY f9y1glutMenuStateFunc()
{
    fmenustate = (GLUTmenuStateFCBUSR) NULL;
    glutMenuStateFunc(NULL);
}

void APIENTRY fglutMenuStateFunc(GLUTmenuStateFCBUSR menustate)
{
    fmenustate = (GLUTmenuStateFCBUSR) menustate;
    glutMenuStateFunc(fortranMenuStateWrapper);
}

void GLUTCALLBACK fortranMenuStatusWrapper(int state, int x, int y)
{
    (*fmenustatus) (&state, &x, &y);
}

void APIENTRY f9y1glutMenuStatusFunc()
{
    fmenustatus = (GLUTmenuStatusFCBUSR) NULL;
    glutMenuStatusFunc(NULL);
}

void APIENTRY fglutMenuStatusFunc(GLUTmenuStatusFCBUSR menustatus)
{
    fmenustatus = (GLUTmenuStatusFCBUSR) menustatus;
    glutMenuStatusFunc(fortranMenuStatusWrapper);
}

void GLUTCALLBACK fortranMotionWrapper(int x, int y)
{
  ((GLUTmotionFCBUSR)__glutGetFCB(GLUT_TMP_MOTION))(&x, &y);
}

void APIENTRY f9y1glutMotionFunc()
{
    glutMotionFunc(NULL);
}

void APIENTRY fglutMotionFunc(GLUTmotionFCBUSR motion)
{
    __glutSetFCB(GLUT_TMP_MOTION,(void *)motion);
    glutMotionFunc(fortranMotionWrapper);
}

void GLUTCALLBACK fortranMouseWrapper(int btn, int state, int x, int y)
{
  ((GLUTmouseFCBUSR)__glutGetFCB(GLUT_TMP_MOUSE))(&btn, &state, &x, &y);
}

void APIENTRY f9y1glutMouseFunc()
{
    glutMouseFunc(NULL);
}

void APIENTRY fglutMouseFunc(GLUTmouseFCBUSR mouse)
{
    __glutSetFCB(GLUT_TMP_MOUSE,(void *)mouse);
    glutMouseFunc(fortranMouseWrapper);
}

void APIENTRY f9y1glutOverlayDisplayFunc()
{
    glutOverlayDisplayFunc(NULL);
}

void APIENTRY fglutOverlayDisplayFunc(GLUTdisplayFCBUSR display)
{
    __glutSetFCB(GLUT_TMP_OVERLAY_DISPLAY,(void *)display);
    glutOverlayDisplayFunc(fortranDisplayWrapper);
}

void GLUTCALLBACK fortranPassiveMotionWrapper(int x, int y)
{
  ((GLUTpassiveFCBUSR)__glutGetFCB(GLUT_TMP_PASSIVE))(&x, &y);
}

void APIENTRY f9y1glutPassiveMotionFunc()
{
    glutPassiveMotionFunc(NULL);
}

void APIENTRY fglutPassiveMotionFunc(GLUTpassiveFCBUSR passivemotion)
{
    __glutSetFCB(GLUT_TMP_PASSIVE,(void *)passivemotion);
    glutPassiveMotionFunc(fortranPassiveMotionWrapper);
}

void APIENTRY fglutPopWindow( void )
{
   glutPopWindow();
}

void APIENTRY fglutPositionWindow( f90glcint x, f90glcint y )
{
   glutPositionWindow( (int) *x, (int) *y );
}

void APIENTRY fglutPostOverlayRedisplay( void )
{
   glutPostOverlayRedisplay();
}

void APIENTRY fglutPostRedisplay( void )
{
   glutPostRedisplay();
}

void APIENTRY fglutPostWindowOverlayRedisplay( f90glcint win )
{
   glutPostWindowOverlayRedisplay( (int) *win );
}

void APIENTRY fglutPostWindowRedisplay( f90glcint win )
{
   glutPostWindowRedisplay( (int) *win );
}

void APIENTRY fglutPushWindow( void )
{
   glutPushWindow();
}

void APIENTRY fglutRemoveMenuItem( f90glcint item )
{
   glutRemoveMenuItem( (int) *item );
}

void APIENTRY fglutRemoveOverlay( void )
{
   glutRemoveOverlay();
}

void APIENTRY fglutReportErrors( void )
{
   glutReportErrors();
}

void GLUTCALLBACK fortranReshapeWrapper(int width, int height)
{
  ((GLUTreshapeFCBUSR)__glutGetFCB(GLUT_TMP_RESHAPE))(&width, &height);
}

void APIENTRY f9y1glutReshapeFunc()
{
    glutReshapeFunc(NULL);
}

void APIENTRY fglutReshapeFunc(GLUTreshapeFCBUSR reshape)
{
    __glutSetFCB(GLUT_TMP_RESHAPE,(void *)reshape);
    glutReshapeFunc(fortranReshapeWrapper);
}

void APIENTRY fglutReshapeWindow( f90glcint width, f90glcint height )
{
   glutReshapeWindow( (int) *width, (int) *height );
}

void APIENTRY fglutSetColor( f90glcint cell, f90glfloat red, f90glfloat green,
    f90glfloat blue )
{
   glutSetColor( (int) *cell, (GLfloat) *red, (GLfloat) *green,
    (GLfloat) *blue );
}

void APIENTRY fglutSetCursor( f90glcint cursor )
{
   glutSetCursor( (int) *cursor );
}

void APIENTRY f9yglutSetIconTitle( f90glint ARRAY title, f90glint size_title)
{
   char *char_title = (char *)malloc((*size_title+1)*sizeof(char));
   int i;
   for (i=0;i<*size_title;i++) char_title[i]=(char) title ARRAYREF [i];
   char_title[*size_title] = end_of_string;
   glutSetIconTitle( char_title );
   free(char_title);
}

void APIENTRY fglutSetMenu( f90glcint menu )
{
   glutSetMenu( (int) *menu );
}

void APIENTRY fglutSetWindow( f90glcint win )
{
   glutSetWindow( (int) *win );
}

void APIENTRY f9yglutSetWindowTitle(f90glint ARRAY title, f90glint size_title )
{
   char *char_title = (char *)malloc((*size_title+1)*sizeof(char));
   int i;
   for (i=0;i<*size_title;i++) char_title[i]=(char) title ARRAYREF [i];
   char_title[*size_title] = end_of_string;
   glutSetWindowTitle( char_title);
#if !defined(MACOSX)
   free(char_title);
#endif
}

void APIENTRY fglutSetupVideoResizing( void )
{
   glutSetupVideoResizing();
}

void APIENTRY fglutShowOverlay( void )
{
   glutShowOverlay();
}

void APIENTRY fglutShowWindow( void )
{
   glutShowWindow();
}

void APIENTRY fglutSolidCone( f90gldouble base, f90gldouble height,
    f90glint slices, f90glint stacks )
{
   glutSolidCone( (GLdouble) *base, (GLdouble) *height, (GLint) *slices,
    (GLint) *stacks );
}

void APIENTRY fglutSolidCube( f90gldouble size )
{
   glutSolidCube( (GLdouble) *size );
}

void APIENTRY fglutSolidDodecahedron( void )
{
   glutSolidDodecahedron();
}

void APIENTRY fglutSolidIcosahedron( void )
{
   glutSolidIcosahedron();
}

void APIENTRY fglutSolidOctahedron( void )
{
   glutSolidOctahedron();
}

void APIENTRY fglutSolidSphere( f90gldouble radius, f90glint slices,
    f90glint stacks )
{
   glutSolidSphere( (GLdouble) *radius, (GLint) *slices, (GLint) *stacks );
}

void APIENTRY fglutSolidTeapot( f90gldouble size )
{
   glutSolidTeapot( (GLdouble) *size );
}

void APIENTRY fglutSolidTetrahedron( void )
{
   glutSolidTetrahedron();
}

void APIENTRY fglutSolidTorus( f90gldouble innerRadius,
    f90gldouble outerRadius, f90glint sides, f90glint rings )
{
   glutSolidTorus( (GLdouble) *innerRadius, (GLdouble) *outerRadius,
    (GLint) *sides, (GLint) *rings );
}

void GLUTCALLBACK fortranSpaceballButtonWrapper(int button, int state)
{
  ((GLUTspaceButtonFCBUSR)__glutGetFCB(GLUT_TMP_SPACE_BUTTON))(&button, &state);
}

void APIENTRY f9y1glutSpaceballButtonFunc()
{
    glutSpaceballButtonFunc(NULL);
}

void APIENTRY fglutSpaceballButtonFunc(GLUTspaceButtonFCBUSR spaceballbutton)
{
    __glutSetFCB(GLUT_TMP_SPACE_BUTTON,(void *)spaceballbutton);
    glutSpaceballButtonFunc(fortranSpaceballButtonWrapper);
}

void GLUTCALLBACK fortranSpaceballMotionWrapper(int x, int y, int z)
{
  ((GLUTspaceMotionFCBUSR)__glutGetFCB(GLUT_TMP_SPACE_MOTION))(&x, &y, &z);
}

void APIENTRY f9y1glutSpaceballMotionFunc()
{
    glutSpaceballMotionFunc(NULL);
}

void APIENTRY fglutSpaceballMotionFunc(GLUTspaceMotionFCBUSR spaceballmotion)
{
    __glutSetFCB(GLUT_TMP_SPACE_MOTION,(void *)spaceballmotion);
    glutSpaceballMotionFunc(fortranSpaceballMotionWrapper);
}

void GLUTCALLBACK fortranSpaceballRotateWrapper(int x, int y, int z)
{
  ((GLUTspaceRotateFCBUSR)__glutGetFCB(GLUT_TMP_SPACE_ROTATE))(&x, &y, &z);
}

void APIENTRY f9y1glutSpaceballRotateFunc()
{
    glutSpaceballRotateFunc(NULL);
}

void APIENTRY fglutSpaceballRotateFunc(GLUTspaceRotateFCBUSR spaceballrotate)
{
    __glutSetFCB(GLUT_TMP_SPACE_ROTATE,(void *)spaceballrotate);
    glutSpaceballRotateFunc(fortranSpaceballRotateWrapper);
}

void GLUTCALLBACK fortranSpecialWrapper(int key, int x, int y)
{
  ((GLUTspecialFCBUSR)__glutGetFCB(GLUT_TMP_SPECIAL))(&key, &x, &y);
}

void APIENTRY f9y1glutSpecialFunc()
{
    glutSpecialFunc(NULL);
}

void APIENTRY fglutSpecialFunc(GLUTspecialFCBUSR special)
{
    __glutSetFCB(GLUT_TMP_SPECIAL,(void *)special);
    glutSpecialFunc(fortranSpecialWrapper);
}

void APIENTRY fglutStopVideoResizing( void )
{
   glutStopVideoResizing();
}

void APIENTRY fglutStrokeCharacter( f90glint font, f90glcint character )
{
   glutStrokeCharacter( f90glut_get_font(*font), (int) *character );
}

int APIENTRY fglutStrokeLength( f90glint font, f90glint ARRAY string,
   f90glint size_string )
{
   unsigned char *char_string = (unsigned char *)malloc((*size_string+1)*sizeof(unsigned char));
   int i;
   for (i=0;i<*size_string;i++) char_string[i]=(unsigned char) string ARRAYREF [i];
   char_string[*size_string] = end_of_string;
   i = glutStrokeLength( f90glut_get_font(*font), char_string );
   free(char_string);
   return i;
}

int APIENTRY fglutStrokeWidth( f90glint font, f90glcint character )
{
   return (int) glutStrokeWidth( f90glut_get_font(*font), (int) *character );
}

void APIENTRY fglutSwapBuffers( void )
{
   glutSwapBuffers();
}

void GLUTCALLBACK fortranTabletButtonWrapper(int button, int state, int x, int y)
{
  ((GLUTtabletButtonFCBUSR)__glutGetFCB(GLUT_TMP_TABLET_BUTTON))(&button, &state, &x, &y);
}

void APIENTRY f9y1glutTabletButtonFunc()
{
    glutTabletButtonFunc(NULL);
}

void APIENTRY fglutTabletButtonFunc(GLUTtabletButtonFCBUSR tabletbutton)
{
    __glutSetFCB(GLUT_TMP_TABLET_BUTTON,(void *)tabletbutton);
    glutTabletButtonFunc(fortranTabletButtonWrapper);
}

void GLUTCALLBACK fortranTabletMotionWrapper(int x, int y)
{
  ((GLUTtabletMotionFCBUSR)__glutGetFCB(GLUT_TMP_TABLET_MOTION))(&x, &y);
}

void APIENTRY f9y1glutTabletMotionFunc()
{
    glutTabletMotionFunc(NULL);
}

void APIENTRY fglutTabletMotionFunc(GLUTtabletMotionFCBUSR tabletmotion)
{
    __glutSetFCB(GLUT_TMP_TABLET_MOTION,(void *)tabletmotion);
    glutTabletMotionFunc(fortranTabletMotionWrapper);
}

void GLUTCALLBACK fortranTimerWrapper(int value)
{
#ifdef FREEGLUT
  int num;
  ((GLUTtimerFCBUSR)__glutGetTimerFCB(GLUT_TMP_TIMER,value,&num))(&num);
#else
  ((GLUTtimerFCBUSR)__glutGetFCB(GLUT_TMP_TIMER))(&value);
#endif
}

void APIENTRY fglutTimerFunc(unsigned int *msecs, GLUTtimerFCBUSR func,
                             int *value)
{
#ifdef FREEGLUT
    int num;
    __glutSetTimerFCB(GLUT_TMP_TIMER,(void *)func, *value, &num);
    glutTimerFunc(*msecs, fortranTimerWrapper, num);
#else
    glutTimerFunc(*msecs, fortranTimerWrapper, *value);
    __glutSetFCB(GLUT_TMP_TIMER,(void *)func);
#endif
}

void APIENTRY fglutUseLayer( f90glenum layer )
{
   glutUseLayer( (GLenum) *layer );
}

void APIENTRY fglutVideoPan( f90glcint x, f90glcint y, f90glcint width,
    f90glcint height )
{
   glutVideoPan( (int) *x, (int) *y, (int) *width, (int) *height );
}

void APIENTRY fglutVideoResize( f90glcint x, f90glcint y, f90glcint width,
    f90glcint height )
{
   glutVideoResize( (int) *x, (int) *y, (int) *width, (int) *height );
}

int APIENTRY fglutVideoResizeGet( f90glenum param )
{
   return (int) glutVideoResizeGet( (GLenum) *param );
}

void GLUTCALLBACK fortranVisibilityWrapper(int state)
{
  ((GLUTvisibilityFCBUSR)__glutGetFCB(GLUT_TMP_VISIBILITY))(&state);
}

void APIENTRY f9y1glutVisibilityFunc()
{
    glutVisibilityFunc(NULL);
}

void APIENTRY fglutVisibilityFunc(GLUTvisibilityFCBUSR visibility)
{
    __glutSetFCB(GLUT_TMP_VISIBILITY,(void *)visibility);
    glutVisibilityFunc(fortranVisibilityWrapper);
}

void APIENTRY fglutWarpPointer( f90glcint x, f90glcint y )
{
   glutWarpPointer( (int) *x, (int) *y );
}

void GLUTCALLBACK fortranWindowStatusWrapper(int state)
{
  ((GLUTwindowStatusFCBUSR)__glutGetFCB(GLUT_TMP_WINDOW_STATUS))(&state);
}

void APIENTRY f9y1glutWindowStatusFunc()
{
    glutWindowStatusFunc(NULL);
}

void APIENTRY fglutWindowStatusFunc(GLUTwindowStatusFCBUSR windowstatus)
{
    __glutSetFCB(GLUT_TMP_WINDOW_STATUS,(void *)windowstatus);
    glutWindowStatusFunc(fortranWindowStatusWrapper);
}

void APIENTRY fglutWireCone( f90gldouble base, f90gldouble height,
    f90glint slices, f90glint stacks )
{
   glutWireCone( (GLdouble) *base, (GLdouble) *height, (GLint) *slices,
    (GLint) *stacks );
}

void APIENTRY fglutWireCube( f90gldouble size )
{
   glutWireCube( (GLdouble) *size );
}

void APIENTRY fglutWireDodecahedron( void )
{
   glutWireDodecahedron();
}

void APIENTRY fglutWireIcosahedron( void )
{
   glutWireIcosahedron();
}

void APIENTRY fglutWireOctahedron( void )
{
   glutWireOctahedron();
}

void APIENTRY fglutWireSphere( f90gldouble radius, f90glint slices,
    f90glint stacks )
{
   glutWireSphere( (GLdouble) *radius, (GLint) *slices, (GLint) *stacks );
}

void APIENTRY fglutWireTeapot( f90gldouble size )
{
   glutWireTeapot( (GLdouble) *size );
}

void APIENTRY fglutWireTetrahedron( void )
{
   glutWireTetrahedron();
}

void APIENTRY fglutWireTorus( f90gldouble innerRadius,
    f90gldouble outerRadius, f90glint sides, f90glint rings )
{
   glutWireTorus( (GLdouble) *innerRadius, (GLdouble) *outerRadius,
    (GLint) *sides, (GLint) *rings );
}

#ifdef GLUT37
/* GLUT 3.7 */

#if (TMP_GLUT_VERSION == 36  && FREEGLUT<1)
/* for the old GLUT, provide dummies */

void APIENTRY f9y1glutKeyboardUpFunc()
{
}

void APIENTRY fglutKeyboardUpFunc(GLUTkeyboardFCBUSR keyboard)
{
}

void APIENTRY f9y1glutSpecialUpFunc()
{
}

void APIENTRY fglutSpecialUpFunc(GLUTspecialFCBUSR special)
{
}

void APIENTRY f9y1glutJoystickFunc(f90glint pollinterval)
{
}

void APIENTRY fglutJoystickFunc(GLUTjoystickFCBUSR joystick,
    f90glint pollinterval)
{
}

void APIENTRY fglutIgnoreKeyRepeat( f90glint ignore )
{
}

void APIENTRY fglutSetKeyRepeat( f90glint repeatmode )
{
}

void APIENTRY fglutForceJoystickFunc(void)
{
}

void APIENTRY f9yglutGameModeString(f90glint string,
    f90glint size_string )
{
}

int APIENTRY fglutEnterGameMode( void )
{
   return (int) 0;
}

void APIENTRY fglutLeaveGameMode( void )
{
}

int APIENTRY fglutGameModeGet( f90glenum mode )
{
   return (int) 0;
}

#else

void GLUTCALLBACK fortranKeyboardUpWrapper(unsigned char key, int x, int y)
{
  int int_key;
  int_key = (int)key;
  ((GLUTkeyboardFCBUSR)__glutGetFCB(GLUT_TMP_KEYBOARD_UP))(&int_key, &x, &y);
}

void APIENTRY f9y1glutKeyboardUpFunc()
{
    glutKeyboardUpFunc(NULL);
}

void APIENTRY fglutKeyboardUpFunc(GLUTkeyboardFCBUSR keyboard)
{
    __glutSetFCB(GLUT_TMP_KEYBOARD_UP,(void *)keyboard);
    glutKeyboardUpFunc(fortranKeyboardUpWrapper);
}

void GLUTCALLBACK fortranSpecialUpWrapper(int key, int x, int y)
{
  ((GLUTspecialFCBUSR)__glutGetFCB(GLUT_TMP_SPECIAL_UP))(&key, &x, &y);
}

void APIENTRY f9y1glutSpecialUpFunc()
{
    glutSpecialUpFunc(NULL);
}

void APIENTRY fglutSpecialUpFunc(GLUTspecialFCBUSR special)
{
    __glutSetFCB(GLUT_TMP_SPECIAL_UP,(void *)special);
    glutSpecialUpFunc(fortranSpecialUpWrapper);
}

void GLUTCALLBACK fortranJoystickWrapper(unsigned int buttonmask, int x,
    int y, int z)
{
  ((GLUTjoystickFCBUSR)__glutGetFCB(GLUT_TMP_JOYSTICK))(&buttonmask,
    &x, &y, &z);
}

void APIENTRY f9y1glutJoystickFunc(f90glint pollinterval)
{
    glutJoystickFunc(NULL,(GLint) *pollinterval);
}

void APIENTRY fglutJoystickFunc(GLUTjoystickFCBUSR joystick,
    f90glint pollinterval)
{
    __glutSetFCB(GLUT_TMP_JOYSTICK,(void *)joystick);
    glutJoystickFunc(fortranJoystickWrapper,(GLint) *pollinterval);
}

void APIENTRY fglutIgnoreKeyRepeat( f90glint ignore )
{
   glutIgnoreKeyRepeat( (GLint) *ignore );
}

void APIENTRY fglutSetKeyRepeat( f90glint repeatmode )
{
   glutSetKeyRepeat( (GLint) *repeatmode );
}

void APIENTRY fglutForceJoystickFunc(void)
{
   glutForceJoystickFunc();
}

void APIENTRY f9yglutGameModeString(f90glint ARRAY string,
    f90glint size_string )
{
   char *char_string = (char *)malloc((*size_string+1)*sizeof(char));
   int i;
   for (i=0;i<*size_string;i++) char_string[i]=(char) string ARRAYREF [i];
   char_string[*size_string] = end_of_string;
   glutGameModeString( char_string);
   free(char_string);
}

int APIENTRY fglutEnterGameMode( void )
{
   return (int) glutEnterGameMode();
}

void APIENTRY fglutLeaveGameMode( void )
{
   glutLeaveGameMode();
}

int APIENTRY fglutGameModeGet( f90glenum mode )
{
   return (int) glutGameModeGet( (GLenum) *mode );
}

#endif
#endif
