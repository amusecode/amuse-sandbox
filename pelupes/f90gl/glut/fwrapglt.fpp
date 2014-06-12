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

! include support for new routines in GLUT 3.7
#define GLUT37

module opengl_glut

use opengl_kinds
use opengl_glutinterfaces
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

#ifdef F_LANG
#define IMPLICIT_NONE
#else
#define IMPLICIT_NONE implicit none
#endif

#ifdef DONT_PASS_COMMENTS
! GLUTNULLFUNC used to be an integer, and a generic interface that
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
#ifdef IFORT
#define OLDNULLFUNC
#endif
#ifdef F_LANG
#define OLDNULLFUNC
#endif
#ifdef GNU
#define OLDNULLFUNC
#endif
#ifdef PATHSCALE
#define OLDNULLFUNC
#endif
#ifdef HP
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
#endif
#ifdef F_LANG
#define NO_NULLFUNC
#endif

!--------------------------
! Public generic interfaces
!--------------------------

public :: &
#ifndef OLDNULLFUNC
   glutnullfunc, &
#endif
   glutaddmenuentry, &
   glutaddsubmenu, &
   glutattachmenu, &
   glutbitmapcharacter, &
   glutbitmaplength, &
   glutbitmapwidth, &
   glutbuttonboxfunc, &
   glutchangetomenuentry, &
   glutchangetosubmenu, &
   glutcopycolormap, &
   glutcreatemenu, &
   glutcreatesubwindow, &
   glutcreatewindow, &
   glutdestroymenu, &
   glutdestroywindow, &
   glutdetachmenu, &
   glutdeviceget, &
   glutdialsfunc, &
   glutdisplayfunc, &
   glutentryfunc, &
   glutestablishoverlay, &
   glutextensionsupported, &
   glutfullscreen, &
   glutget, &
   glutgetcolor, &
   glutgetmenu, &
   glutgetmodifiers, &
   glutgetwindow, &
   gluthideoverlay, &
   gluthidewindow, &
   gluticonifywindow, &
   glutidlefunc, &
   glutinit, &
   glutinitdisplaymode, &
   glutinitdisplaystring, &
   glutinitwindowposition
public :: &
   glutinitwindowsize, &
   glutkeyboardfunc, &
   glutlayerget, &
   glutmainloop, &
   glutmenustatefunc, &
   glutmenustatusfunc, &
   glutmotionfunc, &
   glutmousefunc, &
   glutoverlaydisplayfunc, &
   glutpassivemotionfunc, &
   glutpopwindow, &
   glutpositionwindow, &
   glutpostoverlayredisplay, &
   glutpostredisplay, &
   glutpostwindowoverlayredisplay, &
   glutpostwindowredisplay, &
   glutpushwindow, &
   glutremovemenuitem, &
   glutremoveoverlay, &
   glutreporterrors, &
   glutreshapefunc, &
   glutreshapewindow, &
   glutsetcolor, &
   glutsetcursor, &
   glutseticontitle, &
   glutsetmenu, &
   glutsetwindow, &
   glutsetwindowtitle, &
   glutsetupvideoresizing, &
   glutshowoverlay, &
   glutshowwindow, &
   glutsolidcone, &
   glutsolidcube, &
   glutsoliddodecahedron, &
   glutsolidicosahedron, &
   glutsolidoctahedron
public :: &
   glutsolidsphere, &
   glutsolidteapot, &
   glutsolidtetrahedron, &
   glutsolidtorus, &
   glutspaceballbuttonfunc, &
   glutspaceballmotionfunc, &
   glutspaceballrotatefunc, &
   glutspecialfunc, &
   glutstopvideoresizing, &
   glutstrokecharacter, &
   glutstrokelength, &
   glutstrokewidth, &
   glutswapbuffers, &
   gluttabletbuttonfunc, &
   gluttabletmotionfunc, &
   gluttimerfunc, &
   glutuselayer, &
   glutvideopan, &
   glutvideoresize, &
   glutvideoresizeget, &
   glutvisibilityfunc, &
   glutwarppointer, &
   glutwindowstatusfunc, &
   glutwirecone, &
   glutwirecube, &
   glutwiredodecahedron, &
   glutwireicosahedron, &
   glutwireoctahedron, &
   glutwiresphere, &
   glutwireteapot, &
   glutwiretetrahedron, &
   glutwiretorus
#ifdef GLUT37
public :: &
   glutkeyboardupfunc, &
   glutspecialupfunc, &
   glutjoystickfunc, &
   glutignorekeyrepeat, &
   glutsetkeyrepeat, &
   glutforcejoystickfunc, &
   glutgamemodestring, &
   glutentergamemode, &
   glutleavegamemode, &
   glutgamemodeget
#endif
#ifdef FREEGLUT
public :: &
   glutleavemainloop, &
   glutsetoption 
#endif

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
   f9xglutaddmenuentry, &
   f9xglutaddsubmenu, &
   f9xglutattachmenu, &
   f9xglutbitmapcharacter, &
   f9xglutbitmaplength, &
   f9xglutbitmapwidth, &
   f9xglutbuttonboxfunc, &
   f9xglutchangetomenuentry, &
   f9xglutchangetosubmenu, &
   f9xglutcopycolormap, &
   f9xglutcreatemenu, &
   f9xglutcreatesubwindow, &
   f9xglutcreatewindow, &
   f9xglutdestroymenu, &
   f9xglutdestroywindow, &
   f9xglutdetachmenu, &
   f9xglutdeviceget, &
   f9xglutdialsfunc, &
   f9xglutdisplayfunc, &
   f9xglutentryfunc, &
   f9xglutestablishoverlay, &
   f9xglutextensionsupported, &
   f9xglutfullscreen, &
   f9xglutget, &
   f9xglutgetcolor, &
   f9xglutgetmenu, &
   f9xglutgetmodifiers, &
   f9xglutgetwindow, &
   f9xgluthideoverlay, &
   f9xgluthidewindow, &
   f9xgluticonifywindow, &
   f9xglutidlefunc
private :: &
   f9xglutinit, &
   f9xglutinitdisplaymode, &
   f9xglutinitdisplaystring, &
   f9xglutinitwindowposition, &
   f9xglutinitwindowsize, &
   f9xglutkeyboardfunc, &
   f9xglutlayerget, &
   f9xglutmainloop, &
   f9xglutmenustatefunc, &
   f9xglutmenustatusfunc, &
   f9xglutmotionfunc, &
   f9xglutmousefunc, &
   f9xglutoverlaydisplayfunc, &
   f9xglutpassivemotionfunc, &
   f9xglutpopwindow, &
   f9xglutpositionwindow, &
   f9xglutpostoverlayredisplay, &
   f9xglutpostredisplay, &
   f9xglutpostwinoverlayredisplay, &
   f9xglutpostwindowredisplay, &
   f9xglutpushwindow, &
   f9xglutremovemenuitem, &
   f9xglutremoveoverlay, &
   f9xglutreporterrors, &
   f9xglutreshapefunc, &
   f9xglutreshapewindow, &
   f9xglutsetcolor
private :: &
   f9xglutsetcursor, &
   f9xglutseticontitle, &
   f9xglutsetmenu, &
   f9xglutsetwindow, &
   f9xglutsetwindowtitle, &
   f9xglutsetupvideoresizing, &
   f9xglutshowoverlay, &
   f9xglutshowwindow, &
   f9xglutsolidcone, &
   f9xglutsolidcube, &
   f9xglutsoliddodecahedron, &
   f9xglutsolidicosahedron, &
   f9xglutsolidoctahedron, &
   f9xglutsolidsphere, &
   f9xglutsolidteapot, &
   f9xglutsolidtetrahedron, &
   f9xglutsolidtorus, &
   f9xglutspaceballbuttonfunc, &
   f9xglutspaceballmotionfunc, &
   f9xglutspaceballrotatefunc, &
   f9xglutspecialfunc, &
   f9xglutstopvideoresizing, &
   f9xglutstrokecharacter, &
   f9xglutstrokelength, &
   f9xglutstrokewidth, &
   f9xglutswapbuffers, &
   f9xgluttabletbuttonfunc, &
   f9xgluttabletmotionfunc, &
   f9xgluttimerfunc
private :: &
   f9xglutuselayer, &
   f9xglutvideopan, &
   f9xglutvideoresize, &
   f9xglutvideoresizeget, &
   f9xglutvisibilityfunc, &
   f9xglutwarppointer, &
   f9xglutwindowstatusfunc, &
   f9xglutwirecone, &
   f9xglutwirecube, &
   f9xglutwiredodecahedron, &
   f9xglutwireicosahedron, &
   f9xglutwireoctahedron, &
   f9xglutwiresphere, &
   f9xglutwireteapot, &
   f9xglutwiretetrahedron, &
   f9xglutwiretorus
#ifdef OLDNULLFUNC
private :: &
   f9x1glutidlefunc, &
   f9x1glutbuttonboxfunc, &
   f9x1glutdialsfunc, &
   f9x1glutentryfunc, &
   f9x1glutkeyboardfunc, &
   f9x1glutmenustatefunc, &
   f9x1glutmenustatusfunc, &
   f9x1glutmotionfunc, &
   f9x1glutmousefunc, &
   f9x1glutoverlaydisplayfunc, &
   f9x1glutpassivemotionfunc, &
   f9x1glutreshapefunc, &
   f9x1glutspaceballbuttonfunc, &
   f9x1glutspaceballmotionfunc, &
   f9x1glutspaceballrotatefunc, &
   f9x1glutspecialfunc, &
   f9x1gluttabletbuttonfunc, &
   f9x1gluttabletmotionfunc, &
   f9x1glutvisibilityfunc, &
   f9x1glutwindowstatusfunc, &
   f9x1glutkeyboardupfunc, &
   f9x1glutspecialupfunc, &
   f9x1glutjoystickfunc
#endif
#ifdef GLUT37
private :: &
   f9xglutkeyboardupfunc, &
   f9xglutspecialupfunc, &
   f9xglutjoystickfunc, &
   f9xglutignorekeyrepeat, &
   f9xglutsetkeyrepeat, &
   f9xglutforcejoystickfunc, &
   f9xglutgamemodestring, &
   f9xglutentergamemode, &
   f9xglutleavegamemode, &
   f9xglutgamemodeget
#endif
#ifdef FREEGLUT
private :: &
   f9xglutleavemainloop, &
   f9xglutsetoption
#endif

!-------------------
! symbolic constants
!-------------------

integer(kind=glenum), parameter, public :: GLUT_RGB = 0
integer(kind=glenum), parameter, public :: GLUT_RGBA = 0
integer(kind=glenum), parameter, public :: GLUT_INDEX = 1
integer(kind=glenum), parameter, public :: GLUT_SINGLE = 0
integer(kind=glenum), parameter, public :: GLUT_DOUBLE = 2
integer(kind=glenum), parameter, public :: GLUT_ACCUM = 4
integer(kind=glenum), parameter, public :: GLUT_ALPHA = 8
integer(kind=glenum), parameter, public :: GLUT_DEPTH = 16
integer(kind=glenum), parameter, public :: GLUT_STENCIL = 32
integer(kind=glenum), parameter, public :: GLUT_MULTISAMPLE = 128
integer(kind=glenum), parameter, public :: GLUT_STEREO = 256
integer(kind=glenum), parameter, public :: GLUT_LUMINANCE = 512
integer(kind=glenum), parameter, public :: GLUT_LEFT_BUTTON = 0
integer(kind=glenum), parameter, public :: GLUT_MIDDLE_BUTTON = 1
integer(kind=glenum), parameter, public :: GLUT_RIGHT_BUTTON = 2
integer(kind=glenum), parameter, public :: GLUT_DOWN = 0
integer(kind=glenum), parameter, public :: GLUT_UP = 1
integer(kind=glenum), parameter, public :: GLUT_KEY_F1 = 1
integer(kind=glenum), parameter, public :: GLUT_KEY_F2 = 2
integer(kind=glenum), parameter, public :: GLUT_KEY_F3 = 3
integer(kind=glenum), parameter, public :: GLUT_KEY_F4 = 4
integer(kind=glenum), parameter, public :: GLUT_KEY_F5 = 5
integer(kind=glenum), parameter, public :: GLUT_KEY_F6 = 6
integer(kind=glenum), parameter, public :: GLUT_KEY_F7 = 7
integer(kind=glenum), parameter, public :: GLUT_KEY_F8 = 8
integer(kind=glenum), parameter, public :: GLUT_KEY_F9 = 9
integer(kind=glenum), parameter, public :: GLUT_KEY_F10 = 10
integer(kind=glenum), parameter, public :: GLUT_KEY_F11 = 11
integer(kind=glenum), parameter, public :: GLUT_KEY_F12 = 12
integer(kind=glenum), parameter, public :: GLUT_KEY_LEFT = 100
integer(kind=glenum), parameter, public :: GLUT_KEY_UP = 101
integer(kind=glenum), parameter, public :: GLUT_KEY_RIGHT = 102
integer(kind=glenum), parameter, public :: GLUT_KEY_DOWN = 103
integer(kind=glenum), parameter, public :: GLUT_KEY_PAGE_UP = 104
integer(kind=glenum), parameter, public :: GLUT_KEY_PAGE_DOWN = 105
integer(kind=glenum), parameter, public :: GLUT_KEY_HOME = 106
integer(kind=glenum), parameter, public :: GLUT_KEY_END = 107
integer(kind=glenum), parameter, public :: GLUT_KEY_INSERT = 108
integer(kind=glenum), parameter, public :: GLUT_LEFT = 0
integer(kind=glenum), parameter, public :: GLUT_ENTERED = 1
integer(kind=glenum), parameter, public :: GLUT_MENU_NOT_IN_USE = 0
integer(kind=glenum), parameter, public :: GLUT_MENU_IN_USE = 1
integer(kind=glenum), parameter, public :: GLUT_NOT_VISIBLE = 0
integer(kind=glenum), parameter, public :: GLUT_VISIBLE = 1
integer(kind=glenum), parameter, public :: GLUT_HIDDEN = 0
integer(kind=glenum), parameter, public :: GLUT_FULLY_RETAINED = 1
integer(kind=glenum), parameter, public :: GLUT_PARTIALLY_RETAINED = 2
integer(kind=glenum), parameter, public :: GLUT_FULLY_COVERED = 3
integer(kind=glenum), parameter, public :: GLUT_RED = 0
integer(kind=glenum), parameter, public :: GLUT_GREEN = 1
integer(kind=glenum), parameter, public :: GLUT_BLUE = 2
integer(kind=glenum), parameter, public :: GLUT_NORMAL = 0
integer(kind=glenum), parameter, public :: GLUT_OVERLAY = 1
integer(kind=glenum), parameter, public :: GLUT_WINDOW_X = 100
integer(kind=glenum), parameter, public :: GLUT_WINDOW_Y = 101
integer(kind=glenum), parameter, public :: GLUT_WINDOW_WIDTH = 102
integer(kind=glenum), parameter, public :: GLUT_WINDOW_HEIGHT = 103
integer(kind=glenum), parameter, public :: GLUT_WINDOW_BUFFER_SIZE = 104
integer(kind=glenum), parameter, public :: GLUT_WINDOW_STENCIL_SIZE = 105
integer(kind=glenum), parameter, public :: GLUT_WINDOW_DEPTH_SIZE = 106
integer(kind=glenum), parameter, public :: GLUT_WINDOW_RED_SIZE = 107
integer(kind=glenum), parameter, public :: GLUT_WINDOW_GREEN_SIZE = 108
integer(kind=glenum), parameter, public :: GLUT_WINDOW_BLUE_SIZE = 109
integer(kind=glenum), parameter, public :: GLUT_WINDOW_ALPHA_SIZE = 110
integer(kind=glenum), parameter, public :: GLUT_WINDOW_ACCUM_RED_SIZE = 111
integer(kind=glenum), parameter, public :: GLUT_WINDOW_ACCUM_GREEN_SIZE = 112
integer(kind=glenum), parameter, public :: GLUT_WINDOW_ACCUM_BLUE_SIZE = 113
integer(kind=glenum), parameter, public :: GLUT_WINDOW_ACCUM_ALPHA_SIZE = 114
integer(kind=glenum), parameter, public :: GLUT_WINDOW_DOUBLEBUFFER = 115
integer(kind=glenum), parameter, public :: GLUT_WINDOW_RGBA = 116
integer(kind=glenum), parameter, public :: GLUT_WINDOW_PARENT = 117
integer(kind=glenum), parameter, public :: GLUT_WINDOW_NUM_CHILDREN = 118
integer(kind=glenum), parameter, public :: GLUT_WINDOW_COLORMAP_SIZE = 119
integer(kind=glenum), parameter, public :: GLUT_WINDOW_NUM_SAMPLES = 120
integer(kind=glenum), parameter, public :: GLUT_WINDOW_STEREO = 121
integer(kind=glenum), parameter, public :: GLUT_WINDOW_CURSOR = 122
integer(kind=glenum), parameter, public :: GLUT_SCREEN_WIDTH = 200
integer(kind=glenum), parameter, public :: GLUT_SCREEN_HEIGHT = 201
integer(kind=glenum), parameter, public :: GLUT_SCREEN_WIDTH_MM = 202
integer(kind=glenum), parameter, public :: GLUT_SCREEN_HEIGHT_MM = 203
integer(kind=glenum), parameter, public :: GLUT_MENU_NUM_ITEMS = 300
integer(kind=glenum), parameter, public :: GLUT_DISPLAY_MODE_POSSIBLE = 400
integer(kind=glenum), parameter, public :: GLUT_INIT_WINDOW_X = 500
integer(kind=glenum), parameter, public :: GLUT_INIT_WINDOW_Y = 501
integer(kind=glenum), parameter, public :: GLUT_INIT_WINDOW_WIDTH = 502
integer(kind=glenum), parameter, public :: GLUT_INIT_WINDOW_HEIGHT = 503
integer(kind=glenum), parameter, public :: GLUT_INIT_DISPLAY_MODE = 504
integer(kind=glenum), parameter, public :: GLUT_ELAPSED_TIME = 700
integer(kind=glenum), parameter, public :: GLUT_HAS_KEYBOARD = 600
integer(kind=glenum), parameter, public :: GLUT_HAS_MOUSE = 601
integer(kind=glenum), parameter, public :: GLUT_HAS_SPACEBALL = 602
integer(kind=glenum), parameter, public :: GLUT_HAS_DIAL_AND_BUTTON_BOX = 603
integer(kind=glenum), parameter, public :: GLUT_HAS_TABLET = 604
integer(kind=glenum), parameter, public :: GLUT_NUM_MOUSE_BUTTONS = 605
integer(kind=glenum), parameter, public :: GLUT_NUM_SPACEBALL_BUTTONS = 606
integer(kind=glenum), parameter, public :: GLUT_NUM_BUTTON_BOX_BUTTONS = 607
integer(kind=glenum), parameter, public :: GLUT_NUM_DIALS = 608
integer(kind=glenum), parameter, public :: GLUT_NUM_TABLET_BUTTONS = 609
integer(kind=glenum), parameter, public :: GLUT_OVERLAY_POSSIBLE = 800
integer(kind=glenum), parameter, public :: GLUT_LAYER_IN_USE = 801
integer(kind=glenum), parameter, public :: GLUT_HAS_OVERLAY = 802
integer(kind=glenum), parameter, public :: GLUT_TRANSPARENT_INDEX = 803
integer(kind=glenum), parameter, public :: GLUT_NORMAL_DAMAGED = 804
integer(kind=glenum), parameter, public :: GLUT_OVERLAY_DAMAGED = 805
integer(kind=glenum), parameter, public :: GLUT_VIDEO_RESIZE_POSSIBLE = 900
integer(kind=glenum), parameter, public :: GLUT_VIDEO_RESIZE_IN_USE = 901
integer(kind=glenum), parameter, public :: GLUT_VIDEO_RESIZE_X_DELTA = 902
integer(kind=glenum), parameter, public :: GLUT_VIDEO_RESIZE_Y_DELTA = 903
integer(kind=glenum), parameter, public :: GLUT_VIDEO_RESIZE_WIDTH_DELTA = 904
integer(kind=glenum), parameter, public :: GLUT_VIDEO_RESIZE_HEIGHT_DELTA = 905
integer(kind=glenum), parameter, public :: GLUT_VIDEO_RESIZE_X = 906
integer(kind=glenum), parameter, public :: GLUT_VIDEO_RESIZE_Y = 907
integer(kind=glenum), parameter, public :: GLUT_VIDEO_RESIZE_WIDTH = 908
integer(kind=glenum), parameter, public :: GLUT_VIDEO_RESIZE_HEIGHT = 909
integer(kind=glenum), parameter, public :: GLUT_ACTIVE_SHIFT = 1
integer(kind=glenum), parameter, public :: GLUT_ACTIVE_CTRL = 2
integer(kind=glenum), parameter, public :: GLUT_ACTIVE_ALT = 4
integer(kind=glenum), parameter, public :: GLUT_CURSOR_RIGHT_ARROW = 0
integer(kind=glenum), parameter, public :: GLUT_CURSOR_LEFT_ARROW = 1
integer(kind=glenum), parameter, public :: GLUT_CURSOR_INFO = 2
integer(kind=glenum), parameter, public :: GLUT_CURSOR_DESTROY = 3
integer(kind=glenum), parameter, public :: GLUT_CURSOR_HELP = 4
integer(kind=glenum), parameter, public :: GLUT_CURSOR_CYCLE = 5
integer(kind=glenum), parameter, public :: GLUT_CURSOR_SPRAY = 6
integer(kind=glenum), parameter, public :: GLUT_CURSOR_WAIT = 7
integer(kind=glenum), parameter, public :: GLUT_CURSOR_TEXT = 8
integer(kind=glenum), parameter, public :: GLUT_CURSOR_CROSSHAIR = 9
integer(kind=glenum), parameter, public :: GLUT_CURSOR_UP_DOWN = 10
integer(kind=glenum), parameter, public :: GLUT_CURSOR_LEFT_RIGHT = 11
integer(kind=glenum), parameter, public :: GLUT_CURSOR_TOP_SIDE = 12
integer(kind=glenum), parameter, public :: GLUT_CURSOR_BOTTOM_SIDE = 13
integer(kind=glenum), parameter, public :: GLUT_CURSOR_LEFT_SIDE = 14
integer(kind=glenum), parameter, public :: GLUT_CURSOR_RIGHT_SIDE = 15
integer(kind=glenum), parameter, public :: GLUT_CURSOR_TOP_LEFT_CORNER = 16
integer(kind=glenum), parameter, public :: GLUT_CURSOR_TOP_RIGHT_CORNER = 17
integer(kind=glenum), parameter, public :: GLUT_CURSOR_BOTTOM_RIGHT_CORNER = 18
integer(kind=glenum), parameter, public :: GLUT_CURSOR_BOTTOM_LEFT_CORNER = 19
integer(kind=glenum), parameter, public :: GLUT_CURSOR_INHERIT = 100
integer(kind=glenum), parameter, public :: GLUT_CURSOR_NONE = 101
integer(kind=glenum), parameter, public :: GLUT_CURSOR_FULL_CROSSHAIR = 102
integer(kind=glint), parameter, public :: GLUT_STROKE_ROMAN = 1
integer(kind=glint), parameter, public :: GLUT_STROKE_MONO_ROMAN = 2
integer(kind=glint), parameter, public :: GLUT_BITMAP_9_BY_15 = 3
integer(kind=glint), parameter, public :: GLUT_BITMAP_8_BY_13 = 4
integer(kind=glint), parameter, public :: GLUT_BITMAP_TIMES_ROMAN_10 = 5
integer(kind=glint), parameter, public :: GLUT_BITMAP_TIMES_ROMAN_24 = 6
integer(kind=glint), parameter, public :: GLUT_BITMAP_HELVETICA_10 = 7
integer(kind=glint), parameter, public :: GLUT_BITMAP_HELVETICA_12 = 8
integer(kind=glint), parameter, public :: GLUT_BITMAP_HELVETICA_18 = 9
#ifdef GLUT37
integer(kind=glint), parameter, public :: GLUT_WINDOW_FORMAT_ID = 123
integer(kind=glint), parameter, public :: GLUT_DEVICE_IGNORE_KEY_REPEAT = 610
integer(kind=glint), parameter, public :: GLUT_DEVICE_KEY_REPEAT = 611
integer(kind=glint), parameter, public :: GLUT_HAS_JOYSTICK = 612
integer(kind=glint), parameter, public :: GLUT_OWNS_JOYSTICK = 613
integer(kind=glint), parameter, public :: GLUT_JOYSTICK_BUTTONS = 614
integer(kind=glint), parameter, public :: GLUT_JOYSTICK_AXES = 615
integer(kind=glint), parameter, public :: GLUT_JOYSTICK_POLL_RATE = 616
integer(kind=glint), parameter, public :: GLUT_KEY_REPEAT_OFF = 0
integer(kind=glint), parameter, public :: GLUT_KEY_REPEAT_ON = 1
integer(kind=glint), parameter, public :: GLUT_KEY_REPEAT_DEFAULT = 2
integer(kind=glint), parameter, public :: GLUT_JOYSTICK_BUTTON_A = 1
integer(kind=glint), parameter, public :: GLUT_JOYSTICK_BUTTON_B = 2
integer(kind=glint), parameter, public :: GLUT_JOYSTICK_BUTTON_C = 4
integer(kind=glint), parameter, public :: GLUT_JOYSTICK_BUTTON_D = 8
integer(kind=glint), parameter, public :: GLUT_GAME_MODE_ACTIVE = 0
integer(kind=glint), parameter, public :: GLUT_GAME_MODE_POSSIBLE = 1
integer(kind=glint), parameter, public :: GLUT_GAME_MODE_WIDTH = 2
integer(kind=glint), parameter, public :: GLUT_GAME_MODE_HEIGHT = 3
integer(kind=glint), parameter, public :: GLUT_GAME_MODE_PIXEL_DEPTH = 4
integer(kind=glint), parameter, public :: GLUT_GAME_MODE_REFRESH_RATE = 5
integer(kind=glint), parameter, public :: GLUT_GAME_MODE_DISPLAY_CHANGED = 6
#endif
#ifdef OLDNULLFUNC
integer(kind=glint), parameter, public :: GLUTNULLFUNC = 1
#endif

#ifdef FREEGLUT
integer(kind=glenum), parameter, public :: GLUT_ACTION_ON_WINDOW_CLOSE = 505
integer(kind=glint), parameter, public :: GLUT_ACTION_EXIT  = 0
integer(kind=glint), parameter, public :: GLUT_ACTION_GLUTMAINLOOP_RETURNS = 1
integer(kind=glint), parameter, public :: GLUT_ACTION_CONTINUE_EXECUTION = 2
#endif

!-----------------
! interface blocks
!-----------------

#ifdef OLDNULLFUNC
#else
integer(kind=glint), external :: myglutcomparefuncs
#ifdef SALFORD
STDCALL myglutcomparefuncs 'myglutcomparefuncs' (REF)
#endif
#endif

interface glutaddmenuentry
   module procedure f9xglutaddmenuentry
end interface

interface glutaddsubmenu
   module procedure f9xglutaddsubmenu
end interface

interface glutattachmenu
   module procedure f9xglutattachmenu
end interface

interface glutbitmapcharacter
   module procedure f9xglutbitmapcharacter
end interface

interface glutbitmaplength
   module procedure f9xglutbitmaplength
end interface

interface glutbitmapwidth
   module procedure f9xglutbitmapwidth
end interface

interface glutbuttonboxfunc
   module procedure f9xglutbuttonboxfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutbuttonboxfunc
#endif
end interface

interface glutchangetomenuentry
   module procedure f9xglutchangetomenuentry
end interface

interface glutchangetosubmenu
   module procedure f9xglutchangetosubmenu
end interface

interface glutcopycolormap
   module procedure f9xglutcopycolormap
end interface

interface glutcreatemenu
   module procedure f9xglutcreatemenu
end interface

interface glutcreatesubwindow
   module procedure f9xglutcreatesubwindow
end interface

interface glutcreatewindow
   module procedure f9xglutcreatewindow
end interface

interface glutdestroymenu
   module procedure f9xglutdestroymenu
end interface

interface glutdestroywindow
   module procedure f9xglutdestroywindow
end interface

interface glutdetachmenu
   module procedure f9xglutdetachmenu
end interface

interface glutdeviceget
   module procedure f9xglutdeviceget
end interface

interface glutdialsfunc
   module procedure f9xglutdialsfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutdialsfunc
#endif
end interface

interface glutdisplayfunc
   module procedure f9xglutdisplayfunc
end interface

interface glutentryfunc
   module procedure f9xglutentryfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutentryfunc
#endif
end interface

interface glutestablishoverlay
   module procedure f9xglutestablishoverlay
end interface

interface glutextensionsupported
   module procedure f9xglutextensionsupported
end interface

interface glutfullscreen
   module procedure f9xglutfullscreen
end interface

interface glutget
   module procedure f9xglutget
end interface

interface glutgetcolor
   module procedure f9xglutgetcolor
end interface

interface glutgetmenu
   module procedure f9xglutgetmenu
end interface

interface glutgetmodifiers
   module procedure f9xglutgetmodifiers
end interface

interface glutgetwindow
   module procedure f9xglutgetwindow
end interface

interface gluthideoverlay
   module procedure f9xgluthideoverlay
end interface

interface gluthidewindow
   module procedure f9xgluthidewindow
end interface

interface gluticonifywindow
   module procedure f9xgluticonifywindow
end interface

interface glutidlefunc
   module procedure f9xglutidlefunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutidlefunc
#endif
end interface

interface glutinit
   module procedure f9xglutinit
end interface

interface glutinitdisplaymode
   module procedure f9xglutinitdisplaymode
end interface

interface glutinitdisplaystring
   module procedure f9xglutinitdisplaystring
end interface

interface glutinitwindowposition
   module procedure f9xglutinitwindowposition
end interface

interface glutinitwindowsize
   module procedure f9xglutinitwindowsize
end interface

interface glutkeyboardfunc
   module procedure f9xglutkeyboardfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutkeyboardfunc
#endif
end interface

interface glutlayerget
   module procedure f9xglutlayerget
end interface

interface glutmainloop
   module procedure f9xglutmainloop
end interface

#ifdef FREEGLUT
interface glutleavemainloop
   module procedure f9xglutleavemainloop
end interface

interface glutsetoption
   module procedure f9xglutsetoption
end interface
#endif

interface glutmenustatefunc
   module procedure f9xglutmenustatefunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutmenustatefunc
#endif
end interface

interface glutmenustatusfunc
   module procedure f9xglutmenustatusfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutmenustatusfunc
#endif
end interface

interface glutmotionfunc
   module procedure f9xglutmotionfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutmotionfunc
#endif
end interface

interface glutmousefunc
   module procedure f9xglutmousefunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutmousefunc
#endif
end interface

interface glutoverlaydisplayfunc
   module procedure f9xglutoverlaydisplayfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutoverlaydisplayfunc
#endif
end interface

interface glutpassivemotionfunc
   module procedure f9xglutpassivemotionfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutpassivemotionfunc
#endif
end interface

interface glutpopwindow
   module procedure f9xglutpopwindow
end interface

interface glutpositionwindow
   module procedure f9xglutpositionwindow
end interface

interface glutpostoverlayredisplay
   module procedure f9xglutpostoverlayredisplay
end interface

interface glutpostredisplay
   module procedure f9xglutpostredisplay
end interface

interface glutpostwindowoverlayredisplay
   module procedure f9xglutpostwinoverlayredisplay
end interface

interface glutpostwindowredisplay
   module procedure f9xglutpostwindowredisplay
end interface

interface glutpushwindow
   module procedure f9xglutpushwindow
end interface

interface glutremovemenuitem
   module procedure f9xglutremovemenuitem
end interface

interface glutremoveoverlay
   module procedure f9xglutremoveoverlay
end interface

interface glutreporterrors
   module procedure f9xglutreporterrors
end interface

interface glutreshapefunc
   module procedure f9xglutreshapefunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutreshapefunc
#endif
end interface

interface glutreshapewindow
   module procedure f9xglutreshapewindow
end interface

interface glutsetcolor
   module procedure f9xglutsetcolor
end interface

interface glutsetcursor
   module procedure f9xglutsetcursor
end interface

interface glutseticontitle
   module procedure f9xglutseticontitle
end interface

interface glutsetmenu
   module procedure f9xglutsetmenu
end interface

interface glutsetwindow
   module procedure f9xglutsetwindow
end interface

interface glutsetwindowtitle
   module procedure f9xglutsetwindowtitle
end interface

interface glutsetupvideoresizing
   module procedure f9xglutsetupvideoresizing
end interface

interface glutshowoverlay
   module procedure f9xglutshowoverlay
end interface

interface glutshowwindow
   module procedure f9xglutshowwindow
end interface

interface glutsolidcone
   module procedure f9xglutsolidcone
end interface

interface glutsolidcube
   module procedure f9xglutsolidcube
end interface

interface glutsoliddodecahedron
   module procedure f9xglutsoliddodecahedron
end interface

interface glutsolidicosahedron
   module procedure f9xglutsolidicosahedron
end interface

interface glutsolidoctahedron
   module procedure f9xglutsolidoctahedron
end interface

interface glutsolidsphere
   module procedure f9xglutsolidsphere
end interface

interface glutsolidteapot
   module procedure f9xglutsolidteapot
end interface

interface glutsolidtetrahedron
   module procedure f9xglutsolidtetrahedron
end interface

interface glutsolidtorus
   module procedure f9xglutsolidtorus
end interface

interface glutspaceballbuttonfunc
   module procedure f9xglutspaceballbuttonfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutspaceballbuttonfunc
#endif
end interface

interface glutspaceballmotionfunc
   module procedure f9xglutspaceballmotionfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutspaceballmotionfunc
#endif
end interface

interface glutspaceballrotatefunc
   module procedure f9xglutspaceballrotatefunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutspaceballrotatefunc
#endif
end interface

interface glutspecialfunc
   module procedure f9xglutspecialfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutspecialfunc
#endif
end interface

interface glutstopvideoresizing
   module procedure f9xglutstopvideoresizing
end interface

interface glutstrokecharacter
   module procedure f9xglutstrokecharacter
end interface

interface glutstrokelength
   module procedure f9xglutstrokelength
end interface

interface glutstrokewidth
   module procedure f9xglutstrokewidth
end interface

interface glutswapbuffers
   module procedure f9xglutswapbuffers
end interface

interface gluttabletbuttonfunc
   module procedure f9xgluttabletbuttonfunc
#ifdef OLDNULLFUNC
   module procedure f9x1gluttabletbuttonfunc
#endif
end interface

interface gluttabletmotionfunc
   module procedure f9xgluttabletmotionfunc
#ifdef OLDNULLFUNC
   module procedure f9x1gluttabletmotionfunc
#endif
end interface

interface gluttimerfunc
   module procedure f9xgluttimerfunc
end interface

interface glutuselayer
   module procedure f9xglutuselayer
end interface

interface glutvideopan
   module procedure f9xglutvideopan
end interface

interface glutvideoresize
   module procedure f9xglutvideoresize
end interface

interface glutvideoresizeget
   module procedure f9xglutvideoresizeget
end interface

interface glutvisibilityfunc
   module procedure f9xglutvisibilityfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutvisibilityfunc
#endif
end interface

interface glutwarppointer
   module procedure f9xglutwarppointer
end interface

interface glutwindowstatusfunc
   module procedure f9xglutwindowstatusfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutwindowstatusfunc
#endif
end interface

interface glutwirecone
   module procedure f9xglutwirecone
end interface

interface glutwirecube
   module procedure f9xglutwirecube
end interface

interface glutwiredodecahedron
   module procedure f9xglutwiredodecahedron
end interface

interface glutwireicosahedron
   module procedure f9xglutwireicosahedron
end interface

interface glutwireoctahedron
   module procedure f9xglutwireoctahedron
end interface

interface glutwiresphere
   module procedure f9xglutwiresphere
end interface

interface glutwireteapot
   module procedure f9xglutwireteapot
end interface

interface glutwiretetrahedron
   module procedure f9xglutwiretetrahedron
end interface

interface glutwiretorus
   module procedure f9xglutwiretorus
end interface

#ifdef GLUT37
interface glutkeyboardupfunc
   module procedure f9xglutkeyboardupfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutkeyboardupfunc
#endif
end interface

interface glutspecialupfunc
   module procedure f9xglutspecialupfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutspecialupfunc
#endif
end interface

interface glutjoystickfunc
   module procedure f9xglutjoystickfunc
#ifdef OLDNULLFUNC
   module procedure f9x1glutjoystickfunc
#endif
end interface

interface glutignorekeyrepeat
   module procedure f9xglutignorekeyrepeat
end interface

interface glutsetkeyrepeat
   module procedure f9xglutsetkeyrepeat
end interface

interface glutforcejoystickfunc
   module procedure f9xglutforcejoystickfunc
end interface

interface glutgamemodestring
   module procedure f9xglutgamemodestring
end interface

interface glutentergamemode
   module procedure f9xglutentergamemode
end interface

interface glutleavegamemode
   module procedure f9xglutleavegamemode
end interface

interface glutgamemodeget
   module procedure f9xglutgamemodeget
end interface

#endif

! Module procedures

contains

#ifndef OLDNULLFUNC
subroutine glutnullfunc()
! this is the GLUTNULLFUNC symbol
end subroutine glutnullfunc
#endif

subroutine f9xglutaddmenuentry(label,value)
character(len=*) INTENT_IN label
integer(kind=glcint) INTENT_IN value
integer(kind=glint), dimension(len_trim(label)) :: int_label
integer :: i
int_label = (/ (ichar(label(i:i)),i=1,len_trim(label)) /)
call f9yglutaddmenuentry(int_label,len_trim(label),value)
return
end subroutine f9xglutaddmenuentry

subroutine f9xglutaddsubmenu(label,submenu)
character(len=*) INTENT_IN label
integer(kind=glcint) INTENT_IN submenu
integer(kind=glint), dimension(len_trim(label)) :: int_label
integer :: i
int_label = (/ (ichar(label(i:i)),i=1,len_trim(label)) /)
call f9yglutaddsubmenu(int_label,len_trim(label),submenu)
return
end subroutine f9xglutaddsubmenu

subroutine f9xglutattachmenu(button)
integer(kind=glcint) INTENT_IN button
call fglutattachmenu(button)
return
end subroutine f9xglutattachmenu

subroutine f9xglutbitmapcharacter(font,xcharacter)
integer(kind=glint) INTENT_IN font
integer(kind=glcint) INTENT_IN xcharacter
call fglutbitmapcharacter(font,xcharacter)
return
end subroutine f9xglutbitmapcharacter

function f9xglutbitmaplength(font,string) result(resf9xglutbitmaplength)
integer(kind=glint) INTENT_IN font
character(len=*) INTENT_IN string
integer(kind=glcint) :: resf9xglutbitmaplength
integer(kind=gluint), dimension(len_trim(string)) :: int_string
integer :: j
integer(kind=glint) :: size_string
int_string = (/ (ichar(string(j:j)),j=1,len_trim(string)) /)
size_string = size(int_string)
resf9xglutbitmaplength = fglutbitmaplength(font,int_string,size_string)
return
end function f9xglutbitmaplength

function f9xglutbitmapwidth(font,xcharacter) result(resf9xglutbitmapwidth)
integer(kind=glint) INTENT_IN font
integer(kind=glcint) INTENT_IN xcharacter
integer(kind=glcint) :: resf9xglutbitmapwidth
resf9xglutbitmapwidth = fglutbitmapwidth(font,xcharacter)
return
end function f9xglutbitmapwidth

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutbuttonboxfunc()
call f9y1glutbuttonboxfunc()
#else
subroutine f9x1glutbuttonboxfunc(buttonboxfunc)
integer(kind=glint) INTENT_IN buttonboxfunc
if (buttonboxfunc == glutnullfunc) then
   call f9y1glutbuttonboxfunc()
end if
#endif
return
end subroutine f9x1glutbuttonboxfunc

subroutine f9xglutbuttonboxfunc(buttonboxfunc)
interface
   C_BINDING subroutine buttonboxfunc(button,state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state
   end subroutine buttonboxfunc
end interface
call fglutbuttonboxfunc(buttonboxfunc)
return
end subroutine f9xglutbuttonboxfunc

#else

subroutine f9xglutbuttonboxfunc(buttonboxfunc)
#ifdef EXTERN_CALLBACK
external buttonboxfunc
#else
interface
   C_BINDING subroutine buttonboxfunc(button,state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state
   end subroutine buttonboxfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,buttonboxfunc)==1) then
   call f9y1glutbuttonboxfunc()
else
   call fglutbuttonboxfunc(buttonboxfunc)
endif
return
end subroutine f9xglutbuttonboxfunc

#endif

subroutine f9xglutchangetomenuentry(item,label,value)
integer(kind=glcint) INTENT_IN item
character(len=*) INTENT_IN label
integer(kind=glcint) INTENT_IN value
integer(kind=glint), dimension(len_trim(label)) :: int_label
integer :: i
int_label = (/ (ichar(label(i:i)),i=1,len_trim(label)) /)
call f9yglutchangetomenuentry(item,int_label,len_trim(label),value)
return
end subroutine f9xglutchangetomenuentry

subroutine f9xglutchangetosubmenu(item,label,submenu)
integer(kind=glcint) INTENT_IN item
character(len=*) INTENT_IN label
integer(kind=glcint) INTENT_IN submenu
integer(kind=glint), dimension(len_trim(label)) :: int_label
integer :: i
int_label = (/ (ichar(label(i:i)),i=1,len_trim(label)) /)
call f9yglutchangetosubmenu(item,int_label,len_trim(label),submenu)
return
end subroutine f9xglutchangetosubmenu

subroutine f9xglutcopycolormap(win)
integer(kind=glcint) INTENT_IN win
call fglutcopycolormap(win)
return
end subroutine f9xglutcopycolormap

#ifdef F_LANG
subroutine f9xglutcreatemenu(selectfunc,menuid)
#else
function f9xglutcreatemenu(selectfunc) result(resf9xglutcreatemenu)
#endif
#ifdef EXTERN_CALLBACK
external selectfunc
#else
interface
   C_BINDING subroutine selectfunc(value)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: value
   end subroutine selectfunc
end interface
#endif
#ifdef F_LANG
integer(kind=glcint), intent(out) :: menuid
#else
integer(kind=glcint) :: resf9xglutcreatemenu
#endif
#ifdef F_LANG
call fglutcreatemenu(selectfunc,menuid)
#else
resf9xglutcreatemenu = fglutcreatemenu(selectfunc)
#endif
return
#ifdef F_LANG
end subroutine f9xglutcreatemenu
#else
end function f9xglutcreatemenu
#endif

function f9xglutcreatesubwindow(win,x,y,width, &
   height) result(resf9xglutcreatesubwindow)
integer(kind=glcint) INTENT_IN win, x, y, width, height
integer(kind=glcint) :: resf9xglutcreatesubwindow
resf9xglutcreatesubwindow = fglutcreatesubwindow(win,x,y,width,height)
return
end function f9xglutcreatesubwindow

function f9xglutcreatewindow(title) result(resf9xglutcreatewindow)
character(len=*) INTENT_IN title
integer(kind=glcint) :: resf9xglutcreatewindow
integer(kind=glint), dimension(len_trim(title)) :: int_title
integer :: i
int_title = (/ (ichar(title(i:i)),i=1,len_trim(title)) /)
resf9xglutcreatewindow = f9yglutcreatewindow(int_title,len_trim(title))
return
end function f9xglutcreatewindow

subroutine f9xglutdestroymenu(menu)
integer(kind=glcint) INTENT_IN menu
call fglutdestroymenu(menu)
return
end subroutine f9xglutdestroymenu

subroutine f9xglutdestroywindow(win)
integer(kind=glcint) INTENT_IN win
call fglutdestroywindow(win)
return
end subroutine f9xglutdestroywindow

subroutine f9xglutdetachmenu(button)
integer(kind=glcint) INTENT_IN button
call fglutdetachmenu(button)
return
end subroutine f9xglutdetachmenu

function f9xglutdeviceget(xtype) result(resf9xglutdeviceget)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glcint) :: resf9xglutdeviceget
resf9xglutdeviceget = fglutdeviceget(xtype)
return
end function f9xglutdeviceget

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutdialsfunc()
call f9y1glutdialsfunc()
#else
subroutine f9x1glutdialsfunc(dialsfunc)
integer(kind=glint) INTENT_IN dialsfunc
if (dialsfunc == glutnullfunc) then
   call f9y1glutdialsfunc()
end if
#endif
return
end subroutine f9x1glutdialsfunc

subroutine f9xglutdialsfunc(dialsfunc)
interface
   C_BINDING subroutine dialsfunc(dial,value)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: dial, value
   end subroutine dialsfunc
end interface
call fglutdialsfunc(dialsfunc)
return
end subroutine f9xglutdialsfunc

#else

subroutine f9xglutdialsfunc(dialsfunc)
#ifdef EXTERN_CALLBACK
external dialsfunc
#else
interface
   C_BINDING subroutine dialsfunc(dial,value)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: dial, value
   end subroutine dialsfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,dialsfunc)==1) then
   call f9y1glutdialsfunc()
else
   call fglutdialsfunc(dialsfunc)
endif
return
end subroutine f9xglutdialsfunc

#endif

subroutine f9xglutdisplayfunc(displayfunc)
#ifdef EXTERN_CALLBACK
external displayfunc
#else
interface
   C_BINDING subroutine displayfunc()
   use opengl_kinds
   IMPLICIT_NONE
   end subroutine displayfunc
end interface
#endif
call fglutdisplayfunc(displayfunc)
return
end subroutine f9xglutdisplayfunc

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutentryfunc()
call f9y1glutentryfunc()
#else
subroutine f9x1glutentryfunc(entryfunc)
integer(kind=glint) INTENT_IN entryfunc
if (entryfunc == glutnullfunc) then
   call f9y1glutentryfunc()
end if
#endif
return
end subroutine f9x1glutentryfunc

subroutine f9xglutentryfunc(entryfunc)
interface
   C_BINDING subroutine entryfunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine entryfunc
end interface
call fglutentryfunc(entryfunc)
return
end subroutine f9xglutentryfunc

#else

subroutine f9xglutentryfunc(entryfunc)
#ifdef EXTERN_CALLBACK
external entryfunc
#else
interface
   C_BINDING subroutine entryfunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine entryfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,entryfunc)==1) then
   call f9y1glutentryfunc()
else
   call fglutentryfunc(entryfunc)
endif
return
end subroutine f9xglutentryfunc

#endif

subroutine f9xglutestablishoverlay()
call fglutestablishoverlay()
return
end subroutine f9xglutestablishoverlay

function f9xglutextensionsupported(name) result(resf9xglutextensionsupported)
character(len=*) INTENT_IN name
integer(kind=glcint) :: resf9xglutextensionsupported
integer(kind=glint), dimension(len_trim(name)) :: int_name
integer :: i
int_name = (/ (ichar(name(i:i)),i=1,len_trim(name)) /)
resf9xglutextensionsupported = f9yglutextensionsupported(int_name,len_trim(name))
return
end function f9xglutextensionsupported

subroutine f9xglutfullscreen()
call fglutfullscreen()
return
end subroutine f9xglutfullscreen

function f9xglutget(xtype) result(resf9xglutget)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glcint) :: resf9xglutget
resf9xglutget = fglutget(xtype)
return
end function f9xglutget

function f9xglutgetcolor(ndx,component) result(resf9xglutgetcolor)
integer(kind=glcint) INTENT_IN ndx, component
real(kind=glfloat) :: resf9xglutgetcolor
resf9xglutgetcolor = fglutgetcolor(ndx,component)
return
end function f9xglutgetcolor

function f9xglutgetmenu() result(resf9xglutgetmenu)
integer(kind=glcint) :: resf9xglutgetmenu
resf9xglutgetmenu = fglutgetmenu()
return
end function f9xglutgetmenu

function f9xglutgetmodifiers() result(resf9xglutgetmodifiers)
integer(kind=glcint) :: resf9xglutgetmodifiers
resf9xglutgetmodifiers = fglutgetmodifiers()
return
end function f9xglutgetmodifiers

function f9xglutgetwindow() result(resf9xglutgetwindow)
integer(kind=glcint) :: resf9xglutgetwindow
resf9xglutgetwindow = fglutgetwindow()
return
end function f9xglutgetwindow

subroutine f9xgluthideoverlay()
call fgluthideoverlay()
return
end subroutine f9xgluthideoverlay

subroutine f9xgluthidewindow()
call fgluthidewindow()
return
end subroutine f9xgluthidewindow

subroutine f9xgluticonifywindow()
call fgluticonifywindow()
return
end subroutine f9xgluticonifywindow

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutidlefunc()
call f9y1glutidlefunc()
#else
subroutine f9x1glutidlefunc(idlefunc)
integer(kind=glint) INTENT_IN idlefunc
if (idlefunc == glutnullfunc) then
   call f9y1glutidlefunc()
end if
#endif
return
end subroutine f9x1glutidlefunc

subroutine f9xglutidlefunc(idlefunc)
interface
   C_BINDING subroutine idlefunc()
   use opengl_kinds
   IMPLICIT_NONE
   end subroutine idlefunc
end interface
call fglutidlefunc(idlefunc)
return
end subroutine f9xglutidlefunc

#else

subroutine f9xglutidlefunc(idlefunc)
#ifdef EXTERN_CALLBACK
external idlefunc
#else
interface
   C_BINDING subroutine idlefunc()
   use opengl_kinds
   IMPLICIT_NONE
   end subroutine idlefunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,idlefunc)==1) then
   call f9y1glutidlefunc()
else
   call fglutidlefunc(idlefunc)
endif
return
end subroutine f9xglutidlefunc

#endif

subroutine f9xglutinit(argc,argv)
integer(kind=glcint), optional INTENT_IN argc
character(len=*), optional, dimension(DEFER_DIM) INTENT_IN argv
integer(kind=glint), allocatable, dimension(:) :: intargv, starts
integer(kind=glint) :: i, j, leng
if (present(argc) .and. present(argv)) then
   allocate(starts(argc+1))
! TEMP workaround for gfortran bug 29941. argv should be a whole array here.
   allocate(intargv(len(argv(1))*argc))
   starts(1) = 1
   do i=1,argc
      leng = len_trim(argv(i))
      do j=1,leng
         intargv(j+starts(i)-1) = ichar(argv(i)(j:j))
      end do
      starts(i+1) = starts(i) + leng
   end do
   call fglutinitargs(argc,intargv,starts)
   deallocate(starts,intargv)
else
   call fglutinit()
end if
return
end subroutine f9xglutinit

subroutine f9xglutinitdisplaymode(mode)
integer(kind=glcuint) INTENT_IN mode
call fglutinitdisplaymode(mode)
return
end subroutine f9xglutinitdisplaymode

subroutine f9xglutinitdisplaystring(string)
character(len=*) INTENT_IN string
integer(kind=glint), dimension(len_trim(string)) :: int_string
integer :: i
int_string = (/ (ichar(string(i:i)),i=1,len_trim(string)) /)
call f9yglutinitdisplaystring(int_string,len_trim(string))
return
end subroutine f9xglutinitdisplaystring

subroutine f9xglutinitwindowposition(x,y)
integer(kind=glcint) INTENT_IN x, y
call fglutinitwindowposition(x,y)
return
end subroutine f9xglutinitwindowposition

subroutine f9xglutinitwindowsize(width,height)
integer(kind=glcint) INTENT_IN width, height
call fglutinitwindowsize(width,height)
return
end subroutine f9xglutinitwindowsize

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutkeyboardfunc()
call f9y1glutkeyboardfunc()
#else
subroutine f9x1glutkeyboardfunc(keyboardfunc)
integer(kind=glint) INTENT_IN keyboardfunc
if (keyboardfunc == glutnullfunc) then
   call f9y1glutkeyboardfunc()
end if
#endif
return
end subroutine f9x1glutkeyboardfunc

subroutine f9xglutkeyboardfunc(keyboardfunc)
interface
   C_BINDING subroutine keyboardfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint) PASS_BY_REF, intent(in out) :: key
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine keyboardfunc
end interface
call fglutkeyboardfunc(keyboardfunc)
return
end subroutine f9xglutkeyboardfunc

#else

subroutine f9xglutkeyboardfunc(keyboardfunc)
#ifdef EXTERN_CALLBACK
external keyboardfunc
#else
interface
   C_BINDING subroutine keyboardfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint) PASS_BY_REF, intent(in out) :: key
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine keyboardfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,keyboardfunc)==1) then
   call f9y1glutkeyboardfunc()
else
   call fglutkeyboardfunc(keyboardfunc)
endif
return
end subroutine f9xglutkeyboardfunc

#endif

function f9xglutlayerget(xtype) result(resf9xglutlayerget)
integer(kind=glenum) INTENT_IN xtype
integer(kind=glcint) :: resf9xglutlayerget
resf9xglutlayerget = fglutlayerget(xtype)
return
end function f9xglutlayerget

subroutine f9xglutmainloop()
call fglutmainloop()
return
end subroutine f9xglutmainloop

#ifdef FREEGLUT
subroutine f9xglutleavemainloop()
call fglutleavemainloop()
return
end subroutine f9xglutleavemainloop

subroutine f9xglutsetoption(option,value)
integer(kind=glenum) INTENT_IN option
integer(kind=glcint) INTENT_IN value
call fglutsetoption(option,value)
return
end subroutine f9xglutsetoption
#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutmenustatefunc()
call f9y1glutmenustatefunc()
#else
subroutine f9x1glutmenustatefunc(menustatefunc)
integer(kind=glint) INTENT_IN menustatefunc
if (menustatefunc == glutnullfunc) then
   call f9y1glutmenustatefunc()
end if
#endif
return
end subroutine f9x1glutmenustatefunc

subroutine f9xglutmenustatefunc(menustatefunc)
interface
   C_BINDING subroutine menustatefunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine menustatefunc
end interface
call fglutmenustatefunc(menustatefunc)
return
end subroutine f9xglutmenustatefunc

#else

subroutine f9xglutmenustatefunc(menustatefunc)
#ifdef EXTERN_CALLBACK
external menustatefunc
#else
interface
   C_BINDING subroutine menustatefunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine menustatefunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,menustatefunc)==1) then
   call f9y1glutmenustatefunc()
else
   call fglutmenustatefunc(menustatefunc)
endif
return
end subroutine f9xglutmenustatefunc

#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutmenustatusfunc()
call f9y1glutmenustatusfunc()
#else
subroutine f9x1glutmenustatusfunc(menustatusfunc)
integer(kind=glint) INTENT_IN menustatusfunc
if (menustatusfunc == glutnullfunc) then
   call f9y1glutmenustatusfunc()
end if
#endif
return
end subroutine f9x1glutmenustatusfunc

subroutine f9xglutmenustatusfunc(menustatusfunc)
interface
   C_BINDING subroutine menustatusfunc(status,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: status, x, y
   end subroutine menustatusfunc
end interface
call fglutmenustatusfunc(menustatusfunc)
return
end subroutine f9xglutmenustatusfunc

#else

subroutine f9xglutmenustatusfunc(menustatusfunc)
#ifdef EXTERN_CALLBACK
external menustatusfunc
#else
interface
   C_BINDING subroutine menustatusfunc(status,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: status, x, y
   end subroutine menustatusfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,menustatusfunc)==1) then
   call f9y1glutmenustatusfunc()
else
   call fglutmenustatusfunc(menustatusfunc)
endif
return
end subroutine f9xglutmenustatusfunc

#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutmotionfunc()
call f9y1glutmotionfunc()
#else
subroutine f9x1glutmotionfunc(motionfunc)
integer(kind=glint) INTENT_IN motionfunc
if (motionfunc == glutnullfunc) then
   call f9y1glutmotionfunc()
end if
#endif
return
end subroutine f9x1glutmotionfunc

subroutine f9xglutmotionfunc(motionfunc)
interface
   C_BINDING subroutine motionfunc(x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine motionfunc
end interface
call fglutmotionfunc(motionfunc)
return
end subroutine f9xglutmotionfunc

#else

subroutine f9xglutmotionfunc(motionfunc)
#ifdef EXTERN_CALLBACK
external motionfunc
#else
interface
   C_BINDING subroutine motionfunc(x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine motionfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,motionfunc)==1) then
   call f9y1glutmotionfunc()
else
   call fglutmotionfunc(motionfunc)
endif
return
end subroutine f9xglutmotionfunc

#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutmousefunc()
call f9y1glutmousefunc()
#else
subroutine f9x1glutmousefunc(mousefunc)
integer(kind=glint) INTENT_IN mousefunc
if (mousefunc == glutnullfunc) then
   call f9y1glutmousefunc()
end if
#endif
return
end subroutine f9x1glutmousefunc

subroutine f9xglutmousefunc(mousefunc)
interface
   C_BINDING subroutine mousefunc(button,state,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state, x, y
   end subroutine mousefunc
end interface
call fglutmousefunc(mousefunc)
return
end subroutine f9xglutmousefunc

#else

subroutine f9xglutmousefunc(mousefunc)
#ifdef EXTERN_CALLBACK
external mousefunc
#else
interface
   C_BINDING subroutine mousefunc(button,state,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state, x, y
   end subroutine mousefunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,mousefunc)==1) then
   call f9y1glutmousefunc()
else
   call fglutmousefunc(mousefunc)
endif
return
end subroutine f9xglutmousefunc

#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutoverlaydisplayfunc()
call f9y1glutoverlaydisplayfunc()
#else
subroutine f9x1glutoverlaydisplayfunc(overlaydisplayfunc)
integer(kind=glint) INTENT_IN overlaydisplayfunc
if (overlaydisplayfunc == glutnullfunc) then
   call f9y1glutoverlaydisplayfunc()
end if
#endif
return
end subroutine f9x1glutoverlaydisplayfunc

subroutine f9xglutoverlaydisplayfunc(overlaydisplayfunc)
interface
   C_BINDING subroutine overlaydisplayfunc()
   use opengl_kinds
   IMPLICIT_NONE
   end subroutine overlaydisplayfunc
end interface
call fglutoverlaydisplayfunc(overlaydisplayfunc)
return
end subroutine f9xglutoverlaydisplayfunc

#else

subroutine f9xglutoverlaydisplayfunc(overlaydisplayfunc)
#ifdef EXTERN_CALLBACK
external overlaydisplayfunc
#else
interface
   C_BINDING subroutine overlaydisplayfunc()
   use opengl_kinds
   IMPLICIT_NONE
   end subroutine overlaydisplayfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,overlaydisplayfunc)==1) then
   call f9y1glutoverlaydisplayfunc()
else
   call fglutoverlaydisplayfunc(overlaydisplayfunc)
endif
return
end subroutine f9xglutoverlaydisplayfunc

#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutpassivemotionfunc()
call f9y1glutpassivemotionfunc()
#else
subroutine f9x1glutpassivemotionfunc(passivemotionfunc)
integer(kind=glint) INTENT_IN passivemotionfunc
if (passivemotionfunc == glutnullfunc) then
   call f9y1glutpassivemotionfunc()
end if
#endif
return
end subroutine f9x1glutpassivemotionfunc

subroutine f9xglutpassivemotionfunc(passivemotionfunc)
interface
   C_BINDING subroutine passivemotionfunc(x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine passivemotionfunc
end interface
call fglutpassivemotionfunc(passivemotionfunc)
return
end subroutine f9xglutpassivemotionfunc

#else

subroutine f9xglutpassivemotionfunc(passivemotionfunc)
#ifdef EXTERN_CALLBACK
external passivemotionfunc
#else
interface
   C_BINDING subroutine passivemotionfunc(x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine passivemotionfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,passivemotionfunc)==1) then
   call f9y1glutpassivemotionfunc()
else
   call fglutpassivemotionfunc(passivemotionfunc)
endif
return
end subroutine f9xglutpassivemotionfunc

#endif

subroutine f9xglutpopwindow()
call fglutpopwindow()
return
end subroutine f9xglutpopwindow

subroutine f9xglutpositionwindow(x,y)
integer(kind=glcint) INTENT_IN x, y
call fglutpositionwindow(x,y)
return
end subroutine f9xglutpositionwindow

subroutine f9xglutpostoverlayredisplay()
call fglutpostoverlayredisplay()
return
end subroutine f9xglutpostoverlayredisplay

subroutine f9xglutpostredisplay()
call fglutpostredisplay()
return
end subroutine f9xglutpostredisplay

subroutine f9xglutpostwinoverlayredisplay(win)
integer(kind=glcint) INTENT_IN win
call fglutpostwindowoverlayredisplay(win)
return
end subroutine f9xglutpostwinoverlayredisplay

subroutine f9xglutpostwindowredisplay(win)
integer(kind=glcint) INTENT_IN win
call fglutpostwindowredisplay(win)
return
end subroutine f9xglutpostwindowredisplay

subroutine f9xglutpushwindow()
call fglutpushwindow()
return
end subroutine f9xglutpushwindow

subroutine f9xglutremovemenuitem(item)
integer(kind=glcint) INTENT_IN item
call fglutremovemenuitem(item)
return
end subroutine f9xglutremovemenuitem

subroutine f9xglutremoveoverlay()
call fglutremoveoverlay()
return
end subroutine f9xglutremoveoverlay

subroutine f9xglutreporterrors()
call fglutreporterrors()
return
end subroutine f9xglutreporterrors

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutreshapefunc()
call f9y1glutreshapefunc()
#else
subroutine f9x1glutreshapefunc(reshapefunc)
integer(kind=glint) INTENT_IN reshapefunc
if (reshapefunc == glutnullfunc) then
   call f9y1glutreshapefunc()
end if
#endif
return
end subroutine f9x1glutreshapefunc

subroutine f9xglutreshapefunc(reshapefunc)
interface
   C_BINDING subroutine reshapefunc(width,height)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: width, height
   end subroutine reshapefunc
end interface
call fglutreshapefunc(reshapefunc)
return
end subroutine f9xglutreshapefunc

#else

subroutine f9xglutreshapefunc(reshapefunc)
#ifdef EXTERN_CALLBACK
external reshapefunc
#else
interface
   C_BINDING subroutine reshapefunc(width,height)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: width, height
   end subroutine reshapefunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,reshapefunc)==1) then
   call f9y1glutreshapefunc()
else
   call fglutreshapefunc(reshapefunc)
endif
return
end subroutine f9xglutreshapefunc

#endif

subroutine f9xglutreshapewindow(width,height)
integer(kind=glcint) INTENT_IN width, height
call fglutreshapewindow(width,height)
return
end subroutine f9xglutreshapewindow

subroutine f9xglutsetcolor(cell,red,green,blue)
integer(kind=glcint) INTENT_IN cell
real(kind=glfloat) INTENT_IN red, green, blue
call fglutsetcolor(cell,red,green,blue)
return
end subroutine f9xglutsetcolor

subroutine f9xglutsetcursor(cursor)
integer(kind=glcint) INTENT_IN cursor
call fglutsetcursor(cursor)
return
end subroutine f9xglutsetcursor

subroutine f9xglutseticontitle(title)
character(len=*) INTENT_IN title
integer(kind=glint), dimension(len_trim(title)) :: int_title
integer :: i
int_title = (/ (ichar(title(i:i)),i=1,len_trim(title)) /)
call f9yglutseticontitle(int_title,len_trim(title))
return
end subroutine f9xglutseticontitle

subroutine f9xglutsetmenu(menu)
integer(kind=glcint) INTENT_IN menu
call fglutsetmenu(menu)
return
end subroutine f9xglutsetmenu

subroutine f9xglutsetwindow(win)
integer(kind=glcint) INTENT_IN win
call fglutsetwindow(win)
return
end subroutine f9xglutsetwindow

subroutine f9xglutsetwindowtitle(title)
character(len=*) INTENT_IN title
integer(kind=glint), dimension(len_trim(title)) :: int_title
integer :: i
int_title = (/ (ichar(title(i:i)),i=1,len_trim(title)) /)
call f9yglutsetwindowtitle(int_title,len_trim(title))
return
end subroutine f9xglutsetwindowtitle

subroutine f9xglutsetupvideoresizing()
call fglutsetupvideoresizing()
return
end subroutine f9xglutsetupvideoresizing

subroutine f9xglutshowoverlay()
call fglutshowoverlay()
return
end subroutine f9xglutshowoverlay

subroutine f9xglutshowwindow()
call fglutshowwindow()
return
end subroutine f9xglutshowwindow

subroutine f9xglutsolidcone(base,height,slices,stacks)
real(kind=gldouble) INTENT_IN base, height
integer(kind=glint) INTENT_IN slices, stacks
call fglutsolidcone(base,height,slices,stacks)
return
end subroutine f9xglutsolidcone

subroutine f9xglutsolidcube(xsize)
real(kind=gldouble) INTENT_IN xsize
call fglutsolidcube(xsize)
return
end subroutine f9xglutsolidcube

subroutine f9xglutsoliddodecahedron()
call fglutsoliddodecahedron()
return
end subroutine f9xglutsoliddodecahedron

subroutine f9xglutsolidicosahedron()
call fglutsolidicosahedron()
return
end subroutine f9xglutsolidicosahedron

subroutine f9xglutsolidoctahedron()
call fglutsolidoctahedron()
return
end subroutine f9xglutsolidoctahedron

subroutine f9xglutsolidsphere(radius,slices,stacks)
real(kind=gldouble) INTENT_IN radius
integer(kind=glint) INTENT_IN slices, stacks
call fglutsolidsphere(radius,slices,stacks)
return
end subroutine f9xglutsolidsphere

subroutine f9xglutsolidteapot(xsize)
real(kind=gldouble) INTENT_IN xsize
call fglutsolidteapot(xsize)
return
end subroutine f9xglutsolidteapot

subroutine f9xglutsolidtetrahedron()
call fglutsolidtetrahedron()
return
end subroutine f9xglutsolidtetrahedron

subroutine f9xglutsolidtorus(innerradius,outerradius,sides,rings)
real(kind=gldouble) INTENT_IN innerradius, outerradius
integer(kind=glint) INTENT_IN sides, rings
call fglutsolidtorus(innerradius,outerradius,sides,rings)
return
end subroutine f9xglutsolidtorus

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutspaceballbuttonfunc()
call f9y1glutspaceballbuttonfunc()
#else
subroutine f9x1glutspaceballbuttonfunc(spaceballbuttonfunc)
integer(kind=glint) INTENT_IN spaceballbuttonfunc
if (spaceballbuttonfunc == glutnullfunc) then
   call f9y1glutspaceballbuttonfunc()
end if
#endif
return
end subroutine f9x1glutspaceballbuttonfunc

subroutine f9xglutspaceballbuttonfunc(spaceballbuttonfunc)
interface
   C_BINDING subroutine spaceballbuttonfunc(button,state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state
   end subroutine spaceballbuttonfunc
end interface
call fglutspaceballbuttonfunc(spaceballbuttonfunc)
return
end subroutine f9xglutspaceballbuttonfunc

#else

subroutine f9xglutspaceballbuttonfunc(spaceballbuttonfunc)
#ifdef EXTERN_CALLBACK
external spaceballbuttonfunc
#else
interface
   C_BINDING subroutine spaceballbuttonfunc(button,state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state
   end subroutine spaceballbuttonfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,spaceballbuttonfunc)==1) then
   call f9y1glutspaceballbuttonfunc()
else
   call fglutspaceballbuttonfunc(spaceballbuttonfunc)
endif
return
end subroutine f9xglutspaceballbuttonfunc

#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutspaceballmotionfunc()
call f9y1glutspaceballmotionfunc()
#else
subroutine f9x1glutspaceballmotionfunc(spaceballmotionfunc)
integer(kind=glint) INTENT_IN spaceballmotionfunc
if (spaceballmotionfunc == glutnullfunc) then
   call f9y1glutspaceballmotionfunc()
end if
#endif
return
end subroutine f9x1glutspaceballmotionfunc

subroutine f9xglutspaceballmotionfunc(spaceballmotionfunc)
interface
   C_BINDING subroutine spaceballmotionfunc(x,y,z)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y, z
   end subroutine spaceballmotionfunc
end interface
call fglutspaceballmotionfunc(spaceballmotionfunc)
return
end subroutine f9xglutspaceballmotionfunc

#else

subroutine f9xglutspaceballmotionfunc(spaceballmotionfunc)
#ifdef EXTERN_CALLBACK
external spaceballmotionfunc
#else
interface
   C_BINDING subroutine spaceballmotionfunc(x,y,z)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y, z
   end subroutine spaceballmotionfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,spaceballmotionfunc)==1) then
   call f9y1glutspaceballmotionfunc()
else
   call fglutspaceballmotionfunc(spaceballmotionfunc)
endif
return
end subroutine f9xglutspaceballmotionfunc

#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutspaceballrotatefunc()
call f9y1glutspaceballrotatefunc()
#else
subroutine f9x1glutspaceballrotatefunc(spaceballrotatefunc)
integer(kind=glint) INTENT_IN spaceballrotatefunc
if (spaceballrotatefunc == glutnullfunc) then
   call f9y1glutspaceballrotatefunc()
end if
#endif
return
end subroutine f9x1glutspaceballrotatefunc

subroutine f9xglutspaceballrotatefunc(spaceballrotatefunc)
interface
   C_BINDING subroutine spaceballrotatefunc(x,y,z)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y, z
   end subroutine spaceballrotatefunc
end interface
call fglutspaceballrotatefunc(spaceballrotatefunc)
return
end subroutine f9xglutspaceballrotatefunc

#else

subroutine f9xglutspaceballrotatefunc(spaceballrotatefunc)
#ifdef EXTERN_CALLBACK
external spaceballrotatefunc
#else
interface
   C_BINDING subroutine spaceballrotatefunc(x,y,z)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y, z
   end subroutine spaceballrotatefunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,spaceballrotatefunc)==1) then
   call f9y1glutspaceballrotatefunc()
else
   call fglutspaceballrotatefunc(spaceballrotatefunc)
endif
return
end subroutine f9xglutspaceballrotatefunc

#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutspecialfunc()
call f9y1glutspecialfunc()
#else
subroutine f9x1glutspecialfunc(specialfunc)
integer(kind=glint) INTENT_IN specialfunc
if (specialfunc == glutnullfunc) then
   call f9y1glutspecialfunc()
end if
#endif
return
end subroutine f9x1glutspecialfunc

subroutine f9xglutspecialfunc(specialfunc)
interface
   C_BINDING subroutine specialfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: key, x, y
   end subroutine specialfunc
end interface
call fglutspecialfunc(specialfunc)
return
end subroutine f9xglutspecialfunc

#else

subroutine f9xglutspecialfunc(specialfunc)
#ifdef EXTERN_CALLBACK
external specialfunc
#else
interface
   C_BINDING subroutine specialfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: key, x, y
   end subroutine specialfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,specialfunc)==1) then
   call f9y1glutspecialfunc()
else
   call fglutspecialfunc(specialfunc)
endif
return
end subroutine f9xglutspecialfunc

#endif

subroutine f9xglutstopvideoresizing()
call fglutstopvideoresizing()
return
end subroutine f9xglutstopvideoresizing

subroutine f9xglutstrokecharacter(font,xcharacter)
integer(kind=glint) INTENT_IN font
integer(kind=glcint) INTENT_IN xcharacter
call fglutstrokecharacter(font,xcharacter)
return
end subroutine f9xglutstrokecharacter

function f9xglutstrokelength(font,string) result(resf9xglutstrokelength)
integer(kind=glint) INTENT_IN font
character(len=*) INTENT_IN string
integer(kind=glcint) :: resf9xglutstrokelength
integer(kind=gluint), dimension(len_trim(string)) :: int_string
integer :: j
int_string = (/ (ichar(string(j:j)),j=1,len_trim(string)) /)
resf9xglutstrokelength = fglutstrokelength(font,int_string,size(int_string))
return
end function f9xglutstrokelength

function f9xglutstrokewidth(font,xcharacter) result(resf9xglutstrokewidth)
integer(kind=glint) INTENT_IN font
integer(kind=glcint) INTENT_IN xcharacter
integer(kind=glcint) :: resf9xglutstrokewidth
resf9xglutstrokewidth = fglutstrokewidth(font,xcharacter)
return
end function f9xglutstrokewidth

subroutine f9xglutswapbuffers()
call fglutswapbuffers()
return
end subroutine f9xglutswapbuffers

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1gluttabletbuttonfunc()
call f9y1gluttabletbuttonfunc()
#else
subroutine f9x1gluttabletbuttonfunc(tabletbuttonfunc)
integer(kind=glint) INTENT_IN tabletbuttonfunc
if (tabletbuttonfunc == glutnullfunc) then
   call f9y1gluttabletbuttonfunc()
end if
#endif
return
end subroutine f9x1gluttabletbuttonfunc

subroutine f9xgluttabletbuttonfunc(tabletbuttonfunc)
interface
   C_BINDING subroutine tabletbuttonfunc(button,state,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state, x, y
   end subroutine tabletbuttonfunc
end interface
call fgluttabletbuttonfunc(tabletbuttonfunc)
return
end subroutine f9xgluttabletbuttonfunc

#else

subroutine f9xgluttabletbuttonfunc(tabletbuttonfunc)
#ifdef EXTERN_CALLBACK
external tabletbuttonfunc
#else
interface
   C_BINDING subroutine tabletbuttonfunc(button,state,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state, x, y
   end subroutine tabletbuttonfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,tabletbuttonfunc)==1) then
   call f9y1gluttabletbuttonfunc()
else
   call fgluttabletbuttonfunc(tabletbuttonfunc)
endif
return
end subroutine f9xgluttabletbuttonfunc

#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1gluttabletmotionfunc()
call f9y1gluttabletmotionfunc()
#else
subroutine f9x1gluttabletmotionfunc(tabletmotionfunc)
integer(kind=glint) INTENT_IN tabletmotionfunc
if (tabletmotionfunc == glutnullfunc) then
   call f9y1gluttabletmotionfunc()
end if
#endif
return
end subroutine f9x1gluttabletmotionfunc

subroutine f9xgluttabletmotionfunc(tabletmotionfunc)
interface
   C_BINDING subroutine tabletmotionfunc(x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine tabletmotionfunc
end interface
call fgluttabletmotionfunc(tabletmotionfunc)
return
end subroutine f9xgluttabletmotionfunc

#else

subroutine f9xgluttabletmotionfunc(tabletmotionfunc)
#ifdef EXTERN_CALLBACK
external tabletmotionfunc
#else
interface
   C_BINDING subroutine tabletmotionfunc(x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine tabletmotionfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,tabletmotionfunc)==1) then
   call f9y1gluttabletmotionfunc()
else
   call fgluttabletmotionfunc(tabletmotionfunc)
endif
return
end subroutine f9xgluttabletmotionfunc

#endif

subroutine f9xgluttimerfunc(millis,timerfunc,value)
integer(kind=glcuint) INTENT_IN millis
#ifdef EXTERN_CALLBACK
external timerfunc
#else
interface
   C_BINDING subroutine timerfunc(value)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: value
   end subroutine timerfunc
end interface
#endif
integer(kind=glcint) INTENT_IN value
call fgluttimerfunc(millis,timerfunc,value)
return
end subroutine f9xgluttimerfunc

subroutine f9xglutuselayer(layer)
integer(kind=glenum) INTENT_IN layer
call fglutuselayer(layer)
return
end subroutine f9xglutuselayer

subroutine f9xglutvideopan(x,y,width,height)
integer(kind=glcint) INTENT_IN x, y, width, height
call fglutvideopan(x,y,width,height)
return
end subroutine f9xglutvideopan

subroutine f9xglutvideoresize(x,y,width,height)
integer(kind=glcint) INTENT_IN x, y, width, height
call fglutvideoresize(x,y,width,height)
return
end subroutine f9xglutvideoresize

function f9xglutvideoresizeget(param) result(resf9xglutvideoresizeget)
integer(kind=glenum) INTENT_IN param
integer(kind=glcint) :: resf9xglutvideoresizeget
resf9xglutvideoresizeget = fglutvideoresizeget(param)
return
end function f9xglutvideoresizeget

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutvisibilityfunc()
call f9y1glutvisibilityfunc()
#else
subroutine f9x1glutvisibilityfunc(visibilityfunc)
integer(kind=glint) INTENT_IN visibilityfunc
if (visibilityfunc == glutnullfunc) then
   call f9y1glutvisibilityfunc()
end if
#endif
return
end subroutine f9x1glutvisibilityfunc

subroutine f9xglutvisibilityfunc(visibilityfunc)
interface
   C_BINDING subroutine visibilityfunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine visibilityfunc
end interface
call fglutvisibilityfunc(visibilityfunc)
return
end subroutine f9xglutvisibilityfunc

#else

subroutine f9xglutvisibilityfunc(visibilityfunc)
#ifdef EXTERN_CALLBACK
external visibilityfunc
#else
interface
   C_BINDING subroutine visibilityfunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine visibilityfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,visibilityfunc)==1) then
   call f9y1glutvisibilityfunc()
else
   call fglutvisibilityfunc(visibilityfunc)
endif
return
end subroutine f9xglutvisibilityfunc

#endif

subroutine f9xglutwarppointer(x,y)
integer(kind=glcint) INTENT_IN x, y
call fglutwarppointer(x,y)
return
end subroutine f9xglutwarppointer

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutwindowstatusfunc()
call f9y1glutwindowstatusfunc()
#else
subroutine f9x1glutwindowstatusfunc(windowstatusfunc)
integer(kind=glint) INTENT_IN windowstatusfunc
if (windowstatusfunc == glutnullfunc) then
   call f9y1glutwindowstatusfunc()
end if
#endif
return
end subroutine f9x1glutwindowstatusfunc

subroutine f9xglutwindowstatusfunc(windowstatusfunc)
interface
   C_BINDING subroutine windowstatusfunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine windowstatusfunc
end interface
call fglutwindowstatusfunc(windowstatusfunc)
return
end subroutine f9xglutwindowstatusfunc

#else

subroutine f9xglutwindowstatusfunc(windowstatusfunc)
#ifdef EXTERN_CALLBACK
external windowstatusfunc
#else
interface
   C_BINDING subroutine windowstatusfunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine windowstatusfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,windowstatusfunc)==1) then
   call f9y1glutwindowstatusfunc()
else
   call fglutwindowstatusfunc(windowstatusfunc)
endif
return
end subroutine f9xglutwindowstatusfunc

#endif

subroutine f9xglutwirecone(base,height,slices,stacks)
real(kind=gldouble) INTENT_IN base, height
integer(kind=glint) INTENT_IN slices, stacks
call fglutwirecone(base,height,slices,stacks)
return
end subroutine f9xglutwirecone

subroutine f9xglutwirecube(xsize)
real(kind=gldouble) INTENT_IN xsize
call fglutwirecube(xsize)
return
end subroutine f9xglutwirecube

subroutine f9xglutwiredodecahedron()
call fglutwiredodecahedron()
return
end subroutine f9xglutwiredodecahedron

subroutine f9xglutwireicosahedron()
call fglutwireicosahedron()
return
end subroutine f9xglutwireicosahedron

subroutine f9xglutwireoctahedron()
call fglutwireoctahedron()
return
end subroutine f9xglutwireoctahedron

subroutine f9xglutwiresphere(radius,slices,stacks)
real(kind=gldouble) INTENT_IN radius
integer(kind=glint) INTENT_IN slices, stacks
call fglutwiresphere(radius,slices,stacks)
return
end subroutine f9xglutwiresphere

subroutine f9xglutwireteapot(xsize)
real(kind=gldouble) INTENT_IN xsize
call fglutwireteapot(xsize)
return
end subroutine f9xglutwireteapot

subroutine f9xglutwiretetrahedron()
call fglutwiretetrahedron()
return
end subroutine f9xglutwiretetrahedron

subroutine f9xglutwiretorus(innerradius,outerradius,sides,rings)
real(kind=gldouble) INTENT_IN innerradius, outerradius
integer(kind=glint) INTENT_IN sides, rings
call fglutwiretorus(innerradius,outerradius,sides,rings)
return
end subroutine f9xglutwiretorus

#ifdef GLUT37

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutkeyboardupfunc()
call f9y1glutkeyboardupfunc()
#else
subroutine f9x1glutkeyboardupfunc(keyboardupfunc)
integer(kind=glint) INTENT_IN keyboardupfunc
if (keyboardupfunc == glutnullfunc) then
   call f9y1glutkeyboardupfunc()
end if
#endif
return
end subroutine f9x1glutkeyboardupfunc

subroutine f9xglutkeyboardupfunc(keyboardupfunc)
interface
   C_BINDING subroutine keyboardupfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint) PASS_BY_REF, intent(in out) :: key
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine keyboardupfunc
end interface
call fglutkeyboardupfunc(keyboardupfunc)
return
end subroutine f9xglutkeyboardupfunc

#else

subroutine f9xglutkeyboardupfunc(keyboardupfunc)
#ifdef EXTERN_CALLBACK
external keyboardupfunc
#else
interface
   C_BINDING subroutine keyboardupfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint) PASS_BY_REF, intent(in out) :: key
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine keyboardupfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,keyboardupfunc)==1) then
   call f9y1glutkeyboardupfunc()
else
   call fglutkeyboardupfunc(keyboardupfunc)
endif
return
end subroutine f9xglutkeyboardupfunc

#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutspecialupfunc()
call f9y1glutspecialupfunc()
#else
subroutine f9x1glutspecialupfunc(specialupfunc)
integer(kind=glint) INTENT_IN specialupfunc
if (specialupfunc == glutnullfunc) then
   call f9y1glutspecialupfunc()
end if
#endif
return
end subroutine f9x1glutspecialupfunc

subroutine f9xglutspecialupfunc(specialupfunc)
interface
   C_BINDING subroutine specialupfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: key, x, y
   end subroutine specialupfunc
end interface
call fglutspecialupfunc(specialupfunc)
return
end subroutine f9xglutspecialupfunc

#else

subroutine f9xglutspecialupfunc(specialupfunc)
#ifdef EXTERN_CALLBACK
external specialupfunc
#else
interface
   C_BINDING subroutine specialupfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: key, x, y
   end subroutine specialupfunc
end interface
#endif
if (myglutcomparefuncs(glutnullfunc,specialupfunc)==1) then
   call f9y1glutspecialupfunc()
else
   call fglutspecialupfunc(specialupfunc)
endif
return
end subroutine f9xglutspecialupfunc

#endif

#ifdef OLDNULLFUNC

#ifdef NO_NULLFUNC
subroutine f9x1glutjoystickfunc(pollinterval)
integer(kind=glint) INTENT_IN pollinterval
call f9y1glutjoystickfunc(pollinterval)
#else
subroutine f9x1glutjoystickfunc(joystickfunc,pollinterval)
integer(kind=glint) INTENT_IN joystickfunc
integer(kind=glint) INTENT_IN pollinterval
if (joystickfunc == glutnullfunc) then
   call f9y1glutjoystickfunc(pollinterval)
end if
#endif
return
end subroutine f9x1glutjoystickfunc

subroutine f9xglutjoystickfunc(joystickfunc,pollinterval)
interface
   C_BINDING subroutine joystickfunc(buttonmask,x,y,z)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcuint) PASS_BY_REF, intent(in out) :: buttonmask
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y, z
   end subroutine joystickfunc
end interface
integer(kind=glint) INTENT_IN pollinterval
call fglutjoystickfunc(joystickfunc,pollinterval)
return
end subroutine f9xglutjoystickfunc

#else

subroutine f9xglutjoystickfunc(joystickfunc,pollinterval)
#ifdef EXTERN_CALLBACK
external joystickfunc
#else
interface
   C_BINDING subroutine joystickfunc(buttonmask,x,y,z)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcuint) PASS_BY_REF, intent(in out) :: buttonmask
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y, z
   end subroutine joystickfunc
end interface
#endif
integer(kind=glint) INTENT_IN pollinterval
if (myglutcomparefuncs(glutnullfunc,joystickfunc)==1) then
   call f9y1glutjoystickfunc(pollinterval)
else
   call fglutjoystickfunc(joystickfunc,pollinterval)
endif
return
end subroutine f9xglutjoystickfunc

#endif

subroutine f9xglutignorekeyrepeat(ignore)
integer(kind=glint) INTENT_IN ignore
call fglutignorekeyrepeat(ignore)
return
end subroutine f9xglutignorekeyrepeat

subroutine f9xglutsetkeyrepeat(repeatmode)
integer(kind=glint) INTENT_IN repeatmode
call fglutsetkeyrepeat(repeatmode)
return
end subroutine f9xglutsetkeyrepeat

subroutine f9xglutforcejoystickfunc()
call fglutforcejoystickfunc()
return
end subroutine f9xglutforcejoystickfunc

subroutine f9xglutgamemodestring(string)
character(len=*) INTENT_IN string
integer(kind=glint), dimension(len_trim(string)) :: int_string
integer :: i
int_string = (/ (ichar(string(i:i)),i=1,len_trim(string)) /)
call f9yglutgamemodestring(int_string,len_trim(string))
return
end subroutine f9xglutgamemodestring

function f9xglutentergamemode() result(resf9xglutentergamemode)
integer(kind=glcint) :: resf9xglutentergamemode
resf9xglutentergamemode = fglutentergamemode()
return
end function f9xglutentergamemode

subroutine f9xglutleavegamemode()
call fglutleavegamemode()
return
end subroutine f9xglutleavegamemode

function f9xglutgamemodeget(mode) result(resf9xglutgamemodeget)
integer(kind=glenum) INTENT_IN mode
integer(kind=glcint) :: resf9xglutgamemodeget
resf9xglutgamemodeget = fglutgamemodeget(mode)
return
end function f9xglutgamemodeget

#endif

end module opengl_glut
