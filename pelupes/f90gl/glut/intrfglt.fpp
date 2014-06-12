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

! include support for new routines in GLUT 3.7
#define GLUT37

module opengl_glutinterfaces
use opengl_kinds
#ifdef SALFORD
 STDCALL f9yglutaddmenuentry 'f9yglutaddmenuentry' (REF,REF,REF)
 STDCALL f9yglutaddsubmenu 'f9yglutaddsubmenu' (REF,REF,REF)
 STDCALL fglutattachmenu 'fglutattachmenu' (REF)
 STDCALL fglutbitmapcharacter 'fglutbitmapcharacter' (REF,REF)
 STDCALL fglutbitmaplength 'fglutbitmaplength' (REF,REF,REF)
 STDCALL fglutbitmapwidth 'fglutbitmapwidth' (REF,REF)
 STDCALL f9y1glutbuttonboxfunc 'f9y1glutbuttonboxfunc'
 STDCALL f9yglutchangetomenuentry 'f9yglutchangetomenuentry' (REF,REF,REF,REF)
 STDCALL f9yglutchangetosubmenu 'f9yglutchangetosubmenu' (REF,REF,REF,REF)
 STDCALL fglutcopycolormap 'fglutcopycolormap' (REF)
 STDCALL fglutcreatesubwindow 'fglutcreatesubwindow' (REF,REF,REF,REF,REF)
 STDCALL f9yglutcreatewindow 'f9yglutcreatewindow' (REF,REF)
 STDCALL fglutdestroymenu 'fglutdestroymenu' (REF)
 STDCALL fglutdestroywindow 'fglutdestroywindow' (REF)
 STDCALL fglutdetachmenu 'fglutdetachmenu' (REF)
 STDCALL fglutdeviceget 'fglutdeviceget' (REF)
 STDCALL f9y1glutdialsfunc 'f9y1glutdialsfunc'
 STDCALL f9y1glutentryfunc 'f9y1glutentryfunc'
 STDCALL fglutestablishoverlay 'fglutestablishoverlay'
 STDCALL f9yglutextensionsupported 'f9yglutextensionsupported' (REF,REF)
 STDCALL fglutfullscreen 'fglutfullscreen'
 STDCALL fglutget 'fglutget' (REF)
 STDCALL fglutgetcolor 'fglutgetcolor' (REF,REF)
 STDCALL fglutgetmenu 'fglutgetmenu'
 STDCALL fglutgetmodifiers 'fglutgetmodifiers'
 STDCALL fglutgetwindow 'fglutgetwindow'
 STDCALL fgluthideoverlay 'fgluthideoverlay'
 STDCALL fgluthidewindow 'fgluthidewindow'
 STDCALL fgluticonifywindow 'fgluticonifywindow'
 STDCALL f9y1glutidlefunc 'f9y1glutidlefunc'
 STDCALL fglutinit 'fglutinit'
 STDCALL fglutinitargs 'fglutinitargs' (REF,REF,REF)
 STDCALL fglutinitdisplaymode 'fglutinitdisplaymode' (REF)
 STDCALL f9yglutinitdisplaystring 'f9yglutinitdisplaystring' (REF,REF)
 STDCALL fglutinitwindowposition 'fglutinitwindowposition' (REF,REF)
 STDCALL fglutinitwindowsize 'fglutinitwindowsize' (REF,REF)
 STDCALL f9y1glutkeyboardfunc 'f9y1glutkeyboardfunc'
 STDCALL fglutlayerget 'fglutlayerget' (REF)
 STDCALL fglutmainloop 'fglutmainloop'
 STDCALL f9y1glutmenustatefunc 'f9y1glutmenustatefunc'
 STDCALL f9y1glutmenustatusfunc 'f9y1glutmenustatusfunc'
 STDCALL f9y1glutmotionfunc 'f9y1glutmotionfunc'
 STDCALL f9y1glutmousefunc 'f9y1glutmousefunc'
 STDCALL f9y1glutoverlaydisplayfunc 'f9y1glutoverlaydisplayfunc'
 STDCALL f9y1glutpassivemotionfunc 'f9y1glutpassivemotionfunc'
 STDCALL fglutpopwindow 'fglutpopwindow'
 STDCALL fglutpositionwindow 'fglutpositionwindow' (REF,REF)
 STDCALL fglutpostoverlayredisplay 'fglutpostoverlayredisplay'
 STDCALL fglutpostredisplay 'fglutpostredisplay'
 STDCALL fglutpostwindowoverlayredisplay 'fglutpostwindowoverlayredisplay' (REF)
 STDCALL fglutpostwindowredisplay 'fglutpostwindowredisplay' (REF)
 STDCALL fglutpushwindow 'fglutpushwindow'
 STDCALL fglutremovemenuitem 'fglutremovemenuitem' (REF)
 STDCALL fglutremoveoverlay 'fglutremoveoverlay'
 STDCALL fglutreporterrors 'fglutreporterrors'
 STDCALL f9y1glutreshapefunc 'f9y1glutreshapefunc'
 STDCALL fglutreshapewindow 'fglutreshapewindow' (REF,REF)
 STDCALL fglutsetcolor 'fglutsetcolor' (REF,REF,REF,REF)
 STDCALL fglutsetcursor 'fglutsetcursor' (REF)
 STDCALL f9yglutseticontitle 'f9yglutseticontitle' (REF,REF)
 STDCALL fglutsetmenu 'fglutsetmenu' (REF)
 STDCALL fglutsetwindow 'fglutsetwindow' (REF)
 STDCALL f9yglutsetwindowtitle 'f9yglutsetwindowtitle' (REF,REF)
 STDCALL fglutsetupvideoresizing 'fglutsetupvideoresizing'
 STDCALL fglutshowoverlay 'fglutshowoverlay'
 STDCALL fglutshowwindow 'fglutshowwindow'
 STDCALL fglutsolidcone 'fglutsolidcone' (REF,REF,REF,REF)
 STDCALL fglutsolidcube 'fglutsolidcube' (REF)
 STDCALL fglutsoliddodecahedron 'fglutsoliddodecahedron'
 STDCALL fglutsolidicosahedron 'fglutsolidicosahedron'
 STDCALL fglutsolidoctahedron 'fglutsolidoctahedron'
 STDCALL fglutsolidsphere 'fglutsolidsphere' (REF,REF,REF)
 STDCALL fglutsolidteapot 'fglutsolidteapot' (REF)
 STDCALL fglutsolidtetrahedron 'fglutsolidtetrahedron'
 STDCALL fglutsolidtorus 'fglutsolidtorus' (REF,REF,REF,REF)
 STDCALL f9y1glutspaceballbuttonfunc 'f9y1glutspaceballbuttonfunc'
 STDCALL f9y1glutspaceballmotionfunc 'f9y1glutspaceballmotionfunc'
 STDCALL f9y1glutspaceballrotatefunc 'f9y1glutspaceballrotatefunc'
 STDCALL f9y1glutspecialfunc 'f9y1glutspecialfunc'
 STDCALL fglutstopvideoresizing 'fglutstopvideoresizing'
 STDCALL fglutstrokecharacter 'fglutstrokecharacter' (REF,REF)
 STDCALL fglutstrokelength 'fglutstrokelength' (REF,REF,REF)
 STDCALL fglutstrokewidth 'fglutstrokewidth' (REF,REF)
 STDCALL fglutswapbuffers 'fglutswapbuffers'
 STDCALL f9y1gluttabletbuttonfunc 'f9y1gluttabletbuttonfunc'
 STDCALL f9y1gluttabletmotionfunc 'f9y1gluttabletmotionfunc'
 STDCALL fglutuselayer 'fglutuselayer' (REF)
 STDCALL fglutvideopan 'fglutvideopan' (REF,REF,REF,REF)
 STDCALL fglutvideoresize 'fglutvideoresize' (REF,REF,REF,REF)
 STDCALL fglutvideoresizeget 'fglutvideoresizeget' (REF)
 STDCALL f9y1glutvisibilityfunc 'f9y1glutvisibilityfunc'
 STDCALL fglutwarppointer 'fglutwarppointer' (REF,REF)
 STDCALL f9y1glutwindowstatusfunc 'f9y1glutwindowstatusfunc'
 STDCALL fglutwirecone 'fglutwirecone' (REF,REF,REF,REF)
 STDCALL fglutwirecube 'fglutwirecube' (REF)
 STDCALL fglutwiredodecahedron 'fglutwiredodecahedron'
 STDCALL fglutwireicosahedron 'fglutwireicosahedron'
 STDCALL fglutwireoctahedron 'fglutwireoctahedron'
 STDCALL fglutwiresphere 'fglutwiresphere' (REF,REF,REF)
 STDCALL fglutwireteapot 'fglutwireteapot' (REF)
 STDCALL fglutwiretetrahedron 'fglutwiretetrahedron'
 STDCALL fglutwiretorus 'fglutwiretorus' (REF,REF,REF,REF)
 STDCALL fglutbuttonboxfunc 'fglutbuttonboxfunc' (REF)
 STDCALL fglutcreatemenu 'fglutcreatemenu' (REF)
 STDCALL fglutdialsfunc 'fglutdialsfunc' (REF)
 STDCALL fglutdisplayfunc 'fglutdisplayfunc' (REF)
 STDCALL fglutentryfunc 'fglutentryfunc' (REF)
 STDCALL fglutidlefunc 'fglutidlefunc' (REF)
 STDCALL fglutkeyboardfunc 'fglutkeyboardfunc' (REF)
 STDCALL fglutmenustatefunc 'fglutmenustatefunc' (REF)
 STDCALL fglutmenustatusfunc 'fglutmenustatusfunc' (REF)
 STDCALL fglutmotionfunc 'fglutmotionfunc' (REF)
 STDCALL fglutmousefunc 'fglutmousefunc' (REF)
 STDCALL fglutoverlaydisplayfunc 'fglutoverlaydisplayfunc' (REF)
 STDCALL fglutpassivemotionfunc 'fglutpassivemotionfunc' (REF)
 STDCALL fglutreshapefunc 'fglutreshapefunc' (REF)
 STDCALL fglutspaceballbuttonfunc 'fglutspaceballbuttonfunc' (REF)
 STDCALL fglutspaceballmotionfunc 'fglutspaceballmotionfunc' (REF)
 STDCALL fglutspaceballrotatefunc 'fglutspaceballrotatefunc' (REF)
 STDCALL fglutspecialfunc 'fglutspecialfunc' (REF)
 STDCALL fgluttabletbuttonfunc 'fgluttabletbuttonfunc' (REF)
 STDCALL fgluttabletmotionfunc 'fgluttabletmotionfunc' (REF)
 STDCALL fgluttimerfunc 'fgluttimerfunc' (REF,REF,REF)
 STDCALL fglutvisibilityfunc 'fglutvisibilityfunc' (REF)
 STDCALL fglutwindowstatusfunc 'fglutwindowstatusfunc' (REF)
#ifdef GLUT37
 STDCALL fglutkeyboardupfunc 'fglutkeyboardupfunc' (REF)
 STDCALL f9y1glutkeyboardupfunc 'f9y1glutkeyboardupfunc'
 STDCALL fglutspecialupfunc 'fglutspecialupfunc' (REF)
 STDCALL f9y1glutspecialupfunc 'f9y1glutspecialupfunc'
 STDCALL fglutjoystickfunc 'fglutjoystickfunc' (REF,REF)
 STDCALL f9y1glutjoystickfunc 'f9y1glutjoystickfunc' (REF)
 STDCALL fglutignorekeyrepeat 'fglutignorekeyrepeat' (REF)
 STDCALL fglutsetkeyrepeat 'fglutsetkeyrepeat' (REF)
 STDCALL fglutforcejoystickfunc 'fglutforcejoystickfunc'
 STDCALL f9yglutgamemodestring 'f9yglutgamemodestring' (REF)
 STDCALL fglutentergamemode 'fglutentergamemode'
 STDCALL fglutleavegamemode 'fglutleavegamemode'
 STDCALL fglutgamemodeget 'fglutgamemodeget' (REF)
#endif
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
! C wrapper procedure names
!--------------------------

public :: &
   f9yglutaddmenuentry, &
   f9yglutaddsubmenu, &
   fglutattachmenu, &
   fglutbitmapcharacter, &
   fglutbitmaplength, &
   fglutbitmapwidth, &
   fglutbuttonboxfunc, &
   f9y1glutbuttonboxfunc, &
   f9yglutchangetomenuentry, &
   f9yglutchangetosubmenu, &
   fglutcopycolormap, &
   fglutcreatemenu, &
   fglutcreatesubwindow, &
   f9yglutcreatewindow, &
   fglutdestroymenu, &
   fglutdestroywindow, &
   fglutdetachmenu, &
   fglutdeviceget, &
   fglutdialsfunc, &
   f9y1glutdialsfunc, &
   fglutdisplayfunc, &
   fglutentryfunc, &
   f9y1glutentryfunc, &
   fglutestablishoverlay, &
   f9yglutextensionsupported, &
   fglutfullscreen, &
   fglutget, &
   fglutgetcolor, &
   fglutgetmenu, &
   fglutgetmodifiers, &
   fglutgetwindow, &
   fgluthideoverlay, &
   fgluthidewindow, &
   fgluticonifywindow
public :: &
   fglutidlefunc, &
   f9y1glutidlefunc, &
   fglutinit, &
   fglutinitargs, &
   fglutinitdisplaymode, &
   f9yglutinitdisplaystring, &
   fglutinitwindowposition, &
   fglutinitwindowsize, &
   fglutkeyboardfunc, &
   f9y1glutkeyboardfunc, &
   fglutlayerget, &
   fglutmainloop, &
   fglutmenustatefunc, &
   f9y1glutmenustatefunc, &
   fglutmenustatusfunc, &
   f9y1glutmenustatusfunc, &
   fglutmotionfunc, &
   f9y1glutmotionfunc, &
   fglutmousefunc, &
   f9y1glutmousefunc, &
   fglutoverlaydisplayfunc, &
   f9y1glutoverlaydisplayfunc, &
   fglutpassivemotionfunc, &
   f9y1glutpassivemotionfunc, &
   fglutpopwindow, &
   fglutpositionwindow, &
   fglutpostoverlayredisplay, &
   fglutpostredisplay, &
   fglutpostwindowoverlayredisplay, &
   fglutpostwindowredisplay, &
   fglutpushwindow, &
   fglutremovemenuitem, &
   fglutremoveoverlay, &
   fglutreporterrors, &
   fglutreshapefunc, &
   f9y1glutreshapefunc
public :: &
   fglutreshapewindow, &
   fglutsetcolor, &
   fglutsetcursor, &
   f9yglutseticontitle, &
   fglutsetmenu, &
   fglutsetwindow, &
   f9yglutsetwindowtitle, &
   fglutsetupvideoresizing, &
   fglutshowoverlay, &
   fglutshowwindow, &
   fglutsolidcone, &
   fglutsolidcube, &
   fglutsoliddodecahedron, &
   fglutsolidicosahedron, &
   fglutsolidoctahedron, &
   fglutsolidsphere, &
   fglutsolidteapot, &
   fglutsolidtetrahedron, &
   fglutsolidtorus, &
   fglutspaceballbuttonfunc, &
   f9y1glutspaceballbuttonfunc, &
   fglutspaceballmotionfunc, &
   f9y1glutspaceballmotionfunc, &
   fglutspaceballrotatefunc, &
   f9y1glutspaceballrotatefunc, &
   fglutspecialfunc, &
   f9y1glutspecialfunc, &
   fglutstopvideoresizing, &
   fglutstrokecharacter, &
   fglutstrokelength, &
   fglutstrokewidth, &
   fglutswapbuffers, &
   fgluttabletbuttonfunc, &
   f9y1gluttabletbuttonfunc, &
   fgluttabletmotionfunc
public :: &
   f9y1gluttabletmotionfunc, &
   fgluttimerfunc, &
   fglutuselayer, &
   fglutvideopan, &
   fglutvideoresize, &
   fglutvideoresizeget, &
   fglutvisibilityfunc, &
   f9y1glutvisibilityfunc, &
   fglutwarppointer, &
   fglutwindowstatusfunc, &
   f9y1glutwindowstatusfunc, &
   fglutwirecone, &
   fglutwirecube, &
   fglutwiredodecahedron, &
   fglutwireicosahedron, &
   fglutwireoctahedron, &
   fglutwiresphere, &
   fglutwireteapot, &
   fglutwiretetrahedron, &
   fglutwiretorus
#ifdef GLUT37
public :: &
   fglutkeyboardupfunc, &
   f9y1glutkeyboardupfunc, &
   fglutspecialupfunc, &
   f9y1glutspecialupfunc, &
   fglutjoystickfunc, &
   f9y1glutjoystickfunc, &
   fglutignorekeyrepeat, &
   fglutsetkeyrepeat, &
   fglutforcejoystickfunc, &
   f9yglutgamemodestring, &
   fglutentergamemode, &
   fglutleavegamemode, &
   fglutgamemodeget
#endif

interface
!NAS$ ALIEN "F77 f9yglutaddmenuentry"
subroutine f9yglutaddmenuentry(label,size_label,value)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: label,size_label,value
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglutaddmenuentry
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglutaddmenuentry@12' :: f9yglutaddmenuentry
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN label
integer(kind=glint) INTENT_IN size_label
integer(kind=glcint) INTENT_IN value
end subroutine f9yglutaddmenuentry
end interface

interface
!NAS$ ALIEN "F77 f9yglutaddsubmenu"
subroutine f9yglutaddsubmenu(label,size_label,submenu)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: label,size_label,submenu
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglutaddsubmenu
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglutaddsubmenu@12' :: f9yglutaddsubmenu
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN label
integer(kind=glint) INTENT_IN size_label
integer(kind=glcint) INTENT_IN submenu
end subroutine f9yglutaddsubmenu
end interface

interface
!NAS$ ALIEN "F77 fglutattachmenu"
subroutine fglutattachmenu(button)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: button
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutattachmenu
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutattachmenu@4' :: fglutattachmenu
#endif
integer(kind=glcint) INTENT_IN button
end subroutine fglutattachmenu
end interface

interface
!NAS$ ALIEN "F77 fglutbitmapcharacter"
subroutine fglutbitmapcharacter(font,xcharacter)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: font,xcharacter
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutbitmapcharacter
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutbitmapcharacter@8' :: fglutbitmapcharacter
#endif
integer(kind=glint) INTENT_IN font
integer(kind=glcint) INTENT_IN xcharacter
end subroutine fglutbitmapcharacter
end interface

interface
!NAS$ ALIEN "F77 fglutbitmaplength"
function fglutbitmaplength(font,string,size_string) result(resfglutbitmaplength)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: font,string,size_string
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutbitmaplength
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutbitmaplength@12' :: fglutbitmaplength
#endif
integer(kind=glint) INTENT_IN font
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN string
integer(kind=glint) INTENT_IN size_string
integer(kind=glcint) :: resfglutbitmaplength
end function fglutbitmaplength
end interface

interface
!NAS$ ALIEN "F77 fglutbitmapwidth"
function fglutbitmapwidth(font,xcharacter) result(resfglutbitmapwidth)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: font,xcharacter
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutbitmapwidth
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutbitmapwidth@8' :: fglutbitmapwidth
#endif
integer(kind=glint) INTENT_IN font
integer(kind=glcint) INTENT_IN xcharacter
integer(kind=glcint) :: resfglutbitmapwidth
end function fglutbitmapwidth
end interface

interface
!NAS$ ALIEN "F77 f9y1glutbuttonboxfunc"
subroutine f9y1glutbuttonboxfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutbuttonboxfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutbuttonboxfunc@0' :: f9y1glutbuttonboxfunc
#endif
end subroutine f9y1glutbuttonboxfunc
end interface

interface
!NAS$ ALIEN "F77 f9yglutchangetomenuentry"
subroutine f9yglutchangetomenuentry(item,label,size_label,value)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: item,label,value
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglutchangetomenuentry
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglutchangetomenuentry@16' :: f9yglutchangetomenuentry
#endif
integer(kind=glcint) INTENT_IN item
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN label
integer(kind=glint) INTENT_IN size_label
integer(kind=glcint) INTENT_IN value
end subroutine f9yglutchangetomenuentry
end interface

interface
!NAS$ ALIEN "F77 f9yglutchangetosubmenu"
subroutine f9yglutchangetosubmenu(item,label,size_label,submenu)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: item,label,size_label,submenu
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglutchangetosubmenu
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglutchangetosubmenu@16' :: f9yglutchangetosubmenu
#endif
integer(kind=glcint) INTENT_IN item
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN label
integer(kind=glint) INTENT_IN size_label
integer(kind=glcint) INTENT_IN submenu
end subroutine f9yglutchangetosubmenu
end interface

interface
!NAS$ ALIEN "F77 fglutcopycolormap"
subroutine fglutcopycolormap(win)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: win
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutcopycolormap
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutcopycolormap@4' :: fglutcopycolormap
#endif
integer(kind=glcint) INTENT_IN win
end subroutine fglutcopycolormap
end interface

interface
!NAS$ ALIEN "F77 fglutcreatesubwindow"
function fglutcreatesubwindow(win,x,y,width, &
   height) result(resfglutcreatesubwindow)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: win,x,y,width,height
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutcreatesubwindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutcreatesubwindow@20' :: fglutcreatesubwindow
#endif
integer(kind=glcint) INTENT_IN win, x, y, width, height
integer(kind=glcint) :: resfglutcreatesubwindow
end function fglutcreatesubwindow
end interface

interface
!NAS$ ALIEN "F77 f9yglutcreatewindow"
function f9yglutcreatewindow(title,size_title) result(resf9yglutcreatewindow)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: title,size_title
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglutcreatewindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglutcreatewindow@8' :: f9yglutcreatewindow
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN title
integer(kind=glint) INTENT_IN size_title
integer(kind=glcint) :: resf9yglutcreatewindow
end function f9yglutcreatewindow
end interface

interface
!NAS$ ALIEN "F77 fglutdestroymenu"
subroutine fglutdestroymenu(menu)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: menu
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutdestroymenu
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutdestroymenu@4' :: fglutdestroymenu
#endif
integer(kind=glcint) INTENT_IN menu
end subroutine fglutdestroymenu
end interface

interface
!NAS$ ALIEN "F77 fglutdestroywindow"
subroutine fglutdestroywindow(win)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: win
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutdestroywindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutdestroywindow@4' :: fglutdestroywindow
#endif
integer(kind=glcint) INTENT_IN win
end subroutine fglutdestroywindow
end interface

interface
!NAS$ ALIEN "F77 fglutdetachmenu"
subroutine fglutdetachmenu(button)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: button
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutdetachmenu
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutdetachmenu@4' :: fglutdetachmenu
#endif
integer(kind=glcint) INTENT_IN button
end subroutine fglutdetachmenu
end interface

interface
!NAS$ ALIEN "F77 fglutdeviceget"
function fglutdeviceget(xtype) result(resfglutdeviceget)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutdeviceget
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutdeviceget@4' :: fglutdeviceget
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glcint) :: resfglutdeviceget
end function fglutdeviceget
end interface

interface
!NAS$ ALIEN "F77 f9y1glutdialsfunc"
subroutine f9y1glutdialsfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutdialsfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutdialsfunc@0' :: f9y1glutdialsfunc
#endif
end subroutine f9y1glutdialsfunc
end interface

interface
!NAS$ ALIEN "F77 f9y1glutentryfunc"
subroutine f9y1glutentryfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutentryfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutentryfunc@0' :: f9y1glutentryfunc
#endif
end subroutine f9y1glutentryfunc
end interface

interface
!NAS$ ALIEN "F77 fglutestablishoverlay"
subroutine fglutestablishoverlay()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutestablishoverlay
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutestablishoverlay@0' :: fglutestablishoverlay
#endif
end subroutine fglutestablishoverlay
end interface

interface
!NAS$ ALIEN "F77 f9yglutextensionsupported"
function f9yglutextensionsupported(name,size_name) &
   result(resf9yglutextensionsupported)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: name,size_name
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglutextensionsupported
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglutextensionsupported@8' :: f9yglutextensionsupported
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN name
integer(kind=glint) INTENT_IN size_name
integer(kind=glcint) :: resf9yglutextensionsupported
end function f9yglutextensionsupported
end interface

interface
!NAS$ ALIEN "F77 fglutfullscreen"
subroutine fglutfullscreen()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutfullscreen
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutfullscreen@0' :: fglutfullscreen
#endif
end subroutine fglutfullscreen
end interface

interface
!NAS$ ALIEN "F77 fglutget"
function fglutget(xtype) result(resfglutget)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutget
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutget@4' :: fglutget
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glcint) :: resfglutget
end function fglutget
end interface

interface
!NAS$ ALIEN "F77 fglutgetcolor"
function fglutgetcolor(ndx,component) result(resfglutgetcolor)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: ndx,component
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutgetcolor
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutgetcolor@8' :: fglutgetcolor
#endif
integer(kind=glcint) INTENT_IN ndx, component
real(kind=glfloat) :: resfglutgetcolor
end function fglutgetcolor
end interface

interface
!NAS$ ALIEN "F77 fglutgetmenu"
function fglutgetmenu() result(resfglutgetmenu)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutgetmenu
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutgetmenu@0' :: fglutgetmenu
#endif
integer(kind=glcint) :: resfglutgetmenu
end function fglutgetmenu
end interface

interface
!NAS$ ALIEN "F77 fglutgetmodifiers"
function fglutgetmodifiers() result(resfglutgetmodifiers)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutgetmodifiers
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutgetmodifiers@0' :: fglutgetmodifiers
#endif
integer(kind=glcint) :: resfglutgetmodifiers
end function fglutgetmodifiers
end interface

interface
!NAS$ ALIEN "F77 fglutgetwindow"
function fglutgetwindow() result(resfglutgetwindow)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutgetwindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutgetwindow@0' :: fglutgetwindow
#endif
integer(kind=glcint) :: resfglutgetwindow
end function fglutgetwindow
end interface

interface
!NAS$ ALIEN "F77 fgluthideoverlay"
subroutine fgluthideoverlay()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluthideoverlay
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluthideoverlay@0' :: fgluthideoverlay
#endif
end subroutine fgluthideoverlay
end interface

interface
!NAS$ ALIEN "F77 fgluthidewindow"
subroutine fgluthidewindow()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluthidewindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluthidewindow@0' :: fgluthidewindow
#endif
end subroutine fgluthidewindow
end interface

interface
!NAS$ ALIEN "F77 fgluticonifywindow"
subroutine fgluticonifywindow()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluticonifywindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluticonifywindow@0' :: fgluticonifywindow
#endif
end subroutine fgluticonifywindow
end interface

interface
!NAS$ ALIEN "F77 f9y1glutidlefunc"
subroutine f9y1glutidlefunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutidlefunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutidlefunc@0' :: f9y1glutidlefunc
#endif
end subroutine f9y1glutidlefunc
end interface

interface
!NAS$ ALIEN "F77 fglutinit"
subroutine fglutinit()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutinit
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutinit@0' :: fglutinit
#endif
end subroutine fglutinit
end interface

interface
!NAS$ ALIEN "F77 fglutinitargs"
subroutine fglutinitargs(argcp,argv,starts)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: argcp,argv,starts
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutinitargs
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutinitargs@12' :: fglutinitargs
#endif
integer(kind=glcint) INTENT_IN argcp
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN argv, starts
end subroutine fglutinitargs
end interface

interface
!NAS$ ALIEN "F77 fglutinitdisplaymode"
subroutine fglutinitdisplaymode(mode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutinitdisplaymode
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutinitdisplaymode@4' :: fglutinitdisplaymode
#endif
integer(kind=glcuint) INTENT_IN mode
end subroutine fglutinitdisplaymode
end interface

interface
!NAS$ ALIEN "F77 f9yglutinitdisplaystring"
subroutine f9yglutinitdisplaystring(string,size_string)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: string,size_string
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglutinitdisplaystring
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglutinitdisplaystring@8' :: f9yglutinitdisplaystring
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN string
integer(kind=glint) INTENT_IN size_string
end subroutine f9yglutinitdisplaystring
end interface

interface
!NAS$ ALIEN "F77 fglutinitwindowposition"
subroutine fglutinitwindowposition(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutinitwindowposition
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutinitwindowposition@8' :: fglutinitwindowposition
#endif
integer(kind=glcint) INTENT_IN x, y
end subroutine fglutinitwindowposition
end interface

interface
!NAS$ ALIEN "F77 fglutinitwindowsize"
subroutine fglutinitwindowsize(width,height)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: width,height
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutinitwindowsize
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutinitwindowsize@8' :: fglutinitwindowsize
#endif
integer(kind=glcint) INTENT_IN width, height
end subroutine fglutinitwindowsize
end interface

interface
!NAS$ ALIEN "F77 f9y1glutkeyboardfunc"
subroutine f9y1glutkeyboardfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutkeyboardfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutkeyboardfunc@0' :: f9y1glutkeyboardfunc
#endif
end subroutine f9y1glutkeyboardfunc
end interface

interface
!NAS$ ALIEN "F77 fglutlayerget"
function fglutlayerget(xtype) result(resfglutlayerget)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xtype
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutlayerget
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutlayerget@4' :: fglutlayerget
#endif
integer(kind=glenum) INTENT_IN xtype
integer(kind=glcint) :: resfglutlayerget
end function fglutlayerget
end interface

interface
!NAS$ ALIEN "F77 fglutmainloop"
subroutine fglutmainloop()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutmainloop
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutmainloop@0' :: fglutmainloop
#endif
end subroutine fglutmainloop
end interface

interface
!NAS$ ALIEN "F77 f9y1glutmenustatefunc"
subroutine f9y1glutmenustatefunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutmenustatefunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutmenustatefunc@0' :: f9y1glutmenustatefunc
#endif
end subroutine f9y1glutmenustatefunc
end interface

interface
!NAS$ ALIEN "F77 f9y1glutmenustatusfunc"
subroutine f9y1glutmenustatusfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutmenustatusfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutmenustatusfunc@0' :: f9y1glutmenustatusfunc
#endif
end subroutine f9y1glutmenustatusfunc
end interface

interface
!NAS$ ALIEN "F77 f9y1glutmotionfunc"
subroutine f9y1glutmotionfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutmotionfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutmotionfunc@0' :: f9y1glutmotionfunc
#endif
end subroutine f9y1glutmotionfunc
end interface

interface
!NAS$ ALIEN "F77 f9y1glutmousefunc"
subroutine f9y1glutmousefunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutmousefunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutmousefunc@0' :: f9y1glutmousefunc
#endif
end subroutine f9y1glutmousefunc
end interface

interface
!NAS$ ALIEN "F77 f9y1glutoverlaydisplayfunc"
subroutine f9y1glutoverlaydisplayfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutoverlaydisplayfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutoverlaydisplayfunc@0' :: f9y1glutoverlaydisplayfunc
#endif
end subroutine f9y1glutoverlaydisplayfunc
end interface

interface
!NAS$ ALIEN "F77 f9y1glutpassivemotionfunc"
subroutine f9y1glutpassivemotionfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutpassivemotionfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutpassivemotionfunc@0' :: f9y1glutpassivemotionfunc
#endif
end subroutine f9y1glutpassivemotionfunc
end interface

interface
!NAS$ ALIEN "F77 fglutpopwindow"
subroutine fglutpopwindow()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutpopwindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutpopwindow@0' :: fglutpopwindow
#endif
end subroutine fglutpopwindow
end interface

interface
!NAS$ ALIEN "F77 fglutpositionwindow"
subroutine fglutpositionwindow(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutpositionwindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutpositionwindow@8' :: fglutpositionwindow
#endif
integer(kind=glcint) INTENT_IN x, y
end subroutine fglutpositionwindow
end interface

interface
!NAS$ ALIEN "F77 fglutpostoverlayredisplay"
subroutine fglutpostoverlayredisplay()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutpostoverlayredisplay
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutpostoverlayredisplay@0' :: fglutpostoverlayredisplay
#endif
end subroutine fglutpostoverlayredisplay
end interface

interface
!NAS$ ALIEN "F77 fglutpostredisplay"
subroutine fglutpostredisplay()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutpostredisplay
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutpostredisplay@0' :: fglutpostredisplay
#endif
end subroutine fglutpostredisplay
end interface

interface
!NAS$ ALIEN "F77 fglutpostwindowoverlayredisplay"
subroutine fglutpostwindowoverlayredisplay(win)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: win
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutpostwindowoverlayredisplay
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutpostwindowoverlayredisplay@4' :: fglutpostwindowoverlayredisplay
#endif
integer(kind=glcint) INTENT_IN win
end subroutine fglutpostwindowoverlayredisplay
end interface

interface
!NAS$ ALIEN "F77 fglutpostwindowredisplay"
subroutine fglutpostwindowredisplay(win)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: win
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutpostwindowredisplay
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutpostwindowredisplay@4' :: fglutpostwindowredisplay
#endif
integer(kind=glcint) INTENT_IN win
end subroutine fglutpostwindowredisplay
end interface

interface
!NAS$ ALIEN "F77 fglutpushwindow"
subroutine fglutpushwindow()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutpushwindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutpushwindow@0' :: fglutpushwindow
#endif
end subroutine fglutpushwindow
end interface

interface
!NAS$ ALIEN "F77 fglutremovemenuitem"
subroutine fglutremovemenuitem(item)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: item
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutremovemenuitem
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutremovemenuitem@4' :: fglutremovemenuitem
#endif
integer(kind=glcint) INTENT_IN item
end subroutine fglutremovemenuitem
end interface

interface
!NAS$ ALIEN "F77 fglutremoveoverlay"
subroutine fglutremoveoverlay()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutremoveoverlay
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutremoveoverlay@0' :: fglutremoveoverlay
#endif
end subroutine fglutremoveoverlay
end interface

interface
!NAS$ ALIEN "F77 fglutreporterrors"
subroutine fglutreporterrors()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutreporterrors
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutreporterrors@0' :: fglutreporterrors
#endif
end subroutine fglutreporterrors
end interface

interface
!NAS$ ALIEN "F77 f9y1glutreshapefunc"
subroutine f9y1glutreshapefunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutreshapefunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutreshapefunc@0' :: f9y1glutreshapefunc
#endif
end subroutine f9y1glutreshapefunc
end interface

interface
!NAS$ ALIEN "F77 fglutreshapewindow"
subroutine fglutreshapewindow(width,height)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: width,height
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutreshapewindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutreshapewindow@8' :: fglutreshapewindow
#endif
integer(kind=glcint) INTENT_IN width, height
end subroutine fglutreshapewindow
end interface

interface
!NAS$ ALIEN "F77 fglutsetcolor"
subroutine fglutsetcolor(cell,red,green,blue)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: cell,red,green,blue
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsetcolor
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsetcolor@16' :: fglutsetcolor
#endif
integer(kind=glcint) INTENT_IN cell
real(kind=glfloat) INTENT_IN red, green, blue
end subroutine fglutsetcolor
end interface

interface
!NAS$ ALIEN "F77 fglutsetcursor"
subroutine fglutsetcursor(cursor)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: cursor
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsetcursor
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsetcursor@4' :: fglutsetcursor
#endif
integer(kind=glcint) INTENT_IN cursor
end subroutine fglutsetcursor
end interface

interface
!NAS$ ALIEN "F77 f9yglutseticontitle"
subroutine f9yglutseticontitle(title,size_title)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: title,size_title
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglutseticontitle
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglutseticontitle@8' :: f9yglutseticontitle
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN title
integer(kind=glint) INTENT_IN size_title
end subroutine f9yglutseticontitle
end interface

interface
!NAS$ ALIEN "F77 fglutsetmenu"
subroutine fglutsetmenu(menu)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: menu
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsetmenu
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsetmenu@4' :: fglutsetmenu
#endif
integer(kind=glcint) INTENT_IN menu
end subroutine fglutsetmenu
end interface

interface
!NAS$ ALIEN "F77 fglutsetwindow"
subroutine fglutsetwindow(win)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: win
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsetwindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsetwindow@4' :: fglutsetwindow
#endif
integer(kind=glcint) INTENT_IN win
end subroutine fglutsetwindow
end interface

interface
!NAS$ ALIEN "F77 f9yglutsetwindowtitle"
subroutine f9yglutsetwindowtitle(title,size_title)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: title,size_title
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglutsetwindowtitle
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglutsetwindowtitle@8' :: f9yglutsetwindowtitle
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN title
integer(kind=glint) INTENT_IN size_title
end subroutine f9yglutsetwindowtitle
end interface

interface
!NAS$ ALIEN "F77 fglutsetupvideoresizing"
subroutine fglutsetupvideoresizing()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsetupvideoresizing
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsetupvideoresizing@0' :: fglutsetupvideoresizing
#endif
end subroutine fglutsetupvideoresizing
end interface

interface
!NAS$ ALIEN "F77 fglutshowoverlay"
subroutine fglutshowoverlay()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutshowoverlay
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutshowoverlay@0' :: fglutshowoverlay
#endif
end subroutine fglutshowoverlay
end interface

interface
!NAS$ ALIEN "F77 fglutshowwindow"
subroutine fglutshowwindow()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutshowwindow
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutshowwindow@0' :: fglutshowwindow
#endif
end subroutine fglutshowwindow
end interface

interface
!NAS$ ALIEN "F77 fglutsolidcone"
subroutine fglutsolidcone(base,height,slices,stacks)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: base,height,slices,stacks
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsolidcone
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsolidcone@16' :: fglutsolidcone
#endif
real(kind=gldouble) INTENT_IN base, height
integer(kind=glint) INTENT_IN slices, stacks
end subroutine fglutsolidcone
end interface

interface
!NAS$ ALIEN "F77 fglutsolidcube"
subroutine fglutsolidcube(xsize)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsolidcube
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsolidcube@4' :: fglutsolidcube
#endif
real(kind=gldouble) INTENT_IN xsize
end subroutine fglutsolidcube
end interface

interface
!NAS$ ALIEN "F77 fglutsoliddodecahedron"
subroutine fglutsoliddodecahedron()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsoliddodecahedron
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsoliddodecahedron@0' :: fglutsoliddodecahedron
#endif
end subroutine fglutsoliddodecahedron
end interface

interface
!NAS$ ALIEN "F77 fglutsolidicosahedron"
subroutine fglutsolidicosahedron()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsolidicosahedron
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsolidicosahedron@0' :: fglutsolidicosahedron
#endif
end subroutine fglutsolidicosahedron
end interface

interface
!NAS$ ALIEN "F77 fglutsolidoctahedron"
subroutine fglutsolidoctahedron()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsolidoctahedron
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsolidoctahedron@0' :: fglutsolidoctahedron
#endif
end subroutine fglutsolidoctahedron
end interface

interface
!NAS$ ALIEN "F77 fglutsolidsphere"
subroutine fglutsolidsphere(radius,slices,stacks)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: radius,slices,stacks
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsolidsphere
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsolidsphere@12' :: fglutsolidsphere
#endif
real(kind=gldouble) INTENT_IN radius
integer(kind=glint) INTENT_IN slices, stacks
end subroutine fglutsolidsphere
end interface

interface
!NAS$ ALIEN "F77 fglutsolidteapot"
subroutine fglutsolidteapot(xsize)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsolidteapot
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsolidteapot@4' :: fglutsolidteapot
#endif
real(kind=gldouble) INTENT_IN xsize
end subroutine fglutsolidteapot
end interface

interface
!NAS$ ALIEN "F77 fglutsolidtetrahedron"
subroutine fglutsolidtetrahedron()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsolidtetrahedron
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsolidtetrahedron@0' :: fglutsolidtetrahedron
#endif
end subroutine fglutsolidtetrahedron
end interface

interface
!NAS$ ALIEN "F77 fglutsolidtorus"
subroutine fglutsolidtorus(innerradius,outerradius,sides, &
   rings)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: innerradius,outerradius,sides
!MS$ATTRIBUTES REFERENCE :: rings
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsolidtorus
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsolidtorus@16' :: fglutsolidtorus
#endif
real(kind=gldouble) INTENT_IN innerradius, outerradius
integer(kind=glint) INTENT_IN sides, rings
end subroutine fglutsolidtorus
end interface

interface
!NAS$ ALIEN "F77 f9y1glutspaceballbuttonfunc"
subroutine f9y1glutspaceballbuttonfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutspaceballbuttonfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutspaceballbuttonfunc@0' :: f9y1glutspaceballbuttonfunc
#endif
end subroutine f9y1glutspaceballbuttonfunc
end interface

interface
!NAS$ ALIEN "F77 f9y1glutspaceballmotionfunc"
subroutine f9y1glutspaceballmotionfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutspaceballmotionfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutspaceballmotionfunc@0' :: f9y1glutspaceballmotionfunc
#endif
end subroutine f9y1glutspaceballmotionfunc
end interface

interface
!NAS$ ALIEN "F77 f9y1glutspaceballrotatefunc"
subroutine f9y1glutspaceballrotatefunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutspaceballrotatefunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutspaceballrotatefunc@0' :: f9y1glutspaceballrotatefunc
#endif
end subroutine f9y1glutspaceballrotatefunc
end interface

interface
!NAS$ ALIEN "F77 f9y1glutspecialfunc"
subroutine f9y1glutspecialfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutspecialfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutspecialfunc@0' :: f9y1glutspecialfunc
#endif
end subroutine f9y1glutspecialfunc
end interface

interface
!NAS$ ALIEN "F77 fglutstopvideoresizing"
subroutine fglutstopvideoresizing()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutstopvideoresizing
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutstopvideoresizing@0' :: fglutstopvideoresizing
#endif
end subroutine fglutstopvideoresizing
end interface

interface
!NAS$ ALIEN "F77 fglutstrokecharacter"
subroutine fglutstrokecharacter(font,xcharacter)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: font,xcharacter
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutstrokecharacter
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutstrokecharacter@8' :: fglutstrokecharacter
#endif
integer(kind=glint) INTENT_IN font
integer(kind=glcint) INTENT_IN xcharacter
end subroutine fglutstrokecharacter
end interface

interface
!NAS$ ALIEN "F77 fglutstrokelength"
function fglutstrokelength(font,string,size_string) result(resfglutstrokelength)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: font,string,size_string
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutstrokelength
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutstrokelength@12' :: fglutstrokelength
#endif
integer(kind=glint) INTENT_IN font
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN string
integer(kind=glint) INTENT_IN size_string
integer(kind=glcint) :: resfglutstrokelength
end function fglutstrokelength
end interface

interface
!NAS$ ALIEN "F77 fglutstrokewidth"
function fglutstrokewidth(font,xcharacter) result(resfglutstrokewidth)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: font,xcharacter
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutstrokewidth
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutstrokewidth@8' :: fglutstrokewidth
#endif
integer(kind=glint) INTENT_IN font
integer(kind=glcint) INTENT_IN xcharacter
integer(kind=glcint) :: resfglutstrokewidth
end function fglutstrokewidth
end interface

interface
!NAS$ ALIEN "F77 fglutswapbuffers"
subroutine fglutswapbuffers()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutswapbuffers
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutswapbuffers@0' :: fglutswapbuffers
#endif
end subroutine fglutswapbuffers
end interface

interface
!NAS$ ALIEN "F77 f9y1gluttabletbuttonfunc"
subroutine f9y1gluttabletbuttonfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gluttabletbuttonfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gluttabletbuttonfunc@0' :: f9y1gluttabletbuttonfunc
#endif
end subroutine f9y1gluttabletbuttonfunc
end interface

interface
!NAS$ ALIEN "F77 f9y1gluttabletmotionfunc"
subroutine f9y1gluttabletmotionfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1gluttabletmotionfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1gluttabletmotionfunc@0' :: f9y1gluttabletmotionfunc
#endif
end subroutine f9y1gluttabletmotionfunc
end interface

interface
!NAS$ ALIEN "F77 fglutuselayer"
subroutine fglutuselayer(layer)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: layer
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutuselayer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutuselayer@4' :: fglutuselayer
#endif
integer(kind=glenum) INTENT_IN layer
end subroutine fglutuselayer
end interface

interface
!NAS$ ALIEN "F77 fglutvideopan"
subroutine fglutvideopan(x,y,width,height)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,width,height
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutvideopan
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutvideopan@16' :: fglutvideopan
#endif
integer(kind=glcint) INTENT_IN x, y, width, height
end subroutine fglutvideopan
end interface

interface
!NAS$ ALIEN "F77 fglutvideoresize"
subroutine fglutvideoresize(x,y,width,height)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y,width,height
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutvideoresize
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutvideoresize@16' :: fglutvideoresize
#endif
integer(kind=glcint) INTENT_IN x, y, width, height
end subroutine fglutvideoresize
end interface

interface
!NAS$ ALIEN "F77 fglutvideoresizeget"
function fglutvideoresizeget(param) result(resfglutvideoresizeget)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: param
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutvideoresizeget
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutvideoresizeget@4' :: fglutvideoresizeget
#endif
integer(kind=glenum) INTENT_IN param
integer(kind=glcint) :: resfglutvideoresizeget
end function fglutvideoresizeget
end interface

interface
!NAS$ ALIEN "F77 f9y1glutvisibilityfunc"
subroutine f9y1glutvisibilityfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutvisibilityfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutvisibilityfunc@0' :: f9y1glutvisibilityfunc
#endif
end subroutine f9y1glutvisibilityfunc
end interface

interface
!NAS$ ALIEN "F77 fglutwarppointer"
subroutine fglutwarppointer(x,y)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: x,y
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutwarppointer
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutwarppointer@8' :: fglutwarppointer
#endif
integer(kind=glcint) INTENT_IN x, y
end subroutine fglutwarppointer
end interface

interface
!NAS$ ALIEN "F77 f9y1glutwindowstatusfunc"
subroutine f9y1glutwindowstatusfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutwindowstatusfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutwindowstatusfunc@0' :: f9y1glutwindowstatusfunc
#endif
end subroutine f9y1glutwindowstatusfunc
end interface

interface
!NAS$ ALIEN "F77 fglutwirecone"
subroutine fglutwirecone(base,height,slices,stacks)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: base,height,slices,stacks
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutwirecone
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutwirecone@16' :: fglutwirecone
#endif
real(kind=gldouble) INTENT_IN base, height
integer(kind=glint) INTENT_IN slices, stacks
end subroutine fglutwirecone
end interface

interface
!NAS$ ALIEN "F77 fglutwirecube"
subroutine fglutwirecube(xsize)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutwirecube
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutwirecube@4' :: fglutwirecube
#endif
real(kind=gldouble) INTENT_IN xsize
end subroutine fglutwirecube
end interface

interface
!NAS$ ALIEN "F77 fglutwiredodecahedron"
subroutine fglutwiredodecahedron()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutwiredodecahedron
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutwiredodecahedron@0' :: fglutwiredodecahedron
#endif
end subroutine fglutwiredodecahedron
end interface

interface
!NAS$ ALIEN "F77 fglutwireicosahedron"
subroutine fglutwireicosahedron()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutwireicosahedron
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutwireicosahedron@0' :: fglutwireicosahedron
#endif
end subroutine fglutwireicosahedron
end interface

interface
!NAS$ ALIEN "F77 fglutwireoctahedron"
subroutine fglutwireoctahedron()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutwireoctahedron
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutwireoctahedron@0' :: fglutwireoctahedron
#endif
end subroutine fglutwireoctahedron
end interface

interface
!NAS$ ALIEN "F77 fglutwiresphere"
subroutine fglutwiresphere(radius,slices,stacks)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: radius,slices,stacks
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutwiresphere
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutwiresphere@12' :: fglutwiresphere
#endif
real(kind=gldouble) INTENT_IN radius
integer(kind=glint) INTENT_IN slices, stacks
end subroutine fglutwiresphere
end interface

interface
!NAS$ ALIEN "F77 fglutwireteapot"
subroutine fglutwireteapot(xsize)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: xsize
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutwireteapot
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutwireteapot@4' :: fglutwireteapot
#endif
real(kind=gldouble) INTENT_IN xsize
end subroutine fglutwireteapot
end interface

interface
!NAS$ ALIEN "F77 fglutwiretetrahedron"
subroutine fglutwiretetrahedron()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutwiretetrahedron
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutwiretetrahedron@0' :: fglutwiretetrahedron
#endif
end subroutine fglutwiretetrahedron
end interface

interface
!NAS$ ALIEN "F77 fglutwiretorus"
subroutine fglutwiretorus(innerradius,outerradius,sides,rings)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: innerradius,outerradius,sides,rings
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutwiretorus
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutwiretorus@16' :: fglutwiretorus
#endif
real(kind=gldouble) INTENT_IN innerradius, outerradius
integer(kind=glint) INTENT_IN sides, rings
end subroutine fglutwiretorus
end interface

#ifdef GLUT37

interface
!NAS$ ALIEN "F77 f9y1glutkeyboardupfunc"
subroutine f9y1glutkeyboardupfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutkeyboardupfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutkeyboardupfunc@0' :: f9y1glutkeyboardupfunc
#endif
end subroutine f9y1glutkeyboardupfunc
end interface

interface
!NAS$ ALIEN "F77 f9y1glutspecialupfunc"
subroutine f9y1glutspecialupfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutspecialupfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutspecialupfunc@0' :: f9y1glutspecialupfunc
#endif
end subroutine f9y1glutspecialupfunc
end interface

interface
!NAS$ ALIEN "F77 f9y1glutjoystickfunc"
subroutine f9y1glutjoystickfunc(pollinterval)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: pollinterval
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9y1glutjoystickfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9y1glutjoystickfunc@4' :: f9y1glutjoystickfunc
#endif
integer(kind=glint) INTENT_IN pollinterval
end subroutine f9y1glutjoystickfunc
end interface

interface
!NAS$ ALIEN "F77 fglutignorekeyrepeat"
subroutine fglutignorekeyrepeat(ignore)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: ignore
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutignorekeyrepeat
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutignorekeyrepeat@4' :: fglutignorekeyrepeat
#endif
integer(kind=glint) INTENT_IN ignore
end subroutine fglutignorekeyrepeat
end interface

interface
!NAS$ ALIEN "F77 fglutsetkeyrepeat"
subroutine fglutsetkeyrepeat(repeatmode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: repeatmode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutsetkeyrepeat
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutsetkeyrepeat@4' :: fglutsetkeyrepeat
#endif
integer(kind=glint) INTENT_IN repeatmode
end subroutine fglutsetkeyrepeat
end interface

interface
!NAS$ ALIEN "F77 fglutforcejoystickfunc"
subroutine fglutforcejoystickfunc()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutforcejoystickfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutforcejoystickfunc@0' :: fglutforcejoystickfunc
#endif
end subroutine fglutforcejoystickfunc
end interface

interface
!NAS$ ALIEN "F77 f9yglutgamemodestring"
subroutine f9yglutgamemodestring(string,size_string)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: string,size_string
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: f9yglutgamemodestring
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_f9yglutgamemodestring@8' :: f9yglutgamemodestring
#endif
integer(kind=glint), dimension(DEFER_DIM) INTENT_IN string
integer(kind=glint) INTENT_IN size_string
end subroutine f9yglutgamemodestring
end interface

interface
!NAS$ ALIEN "F77 fglutentergamemode"
function fglutentergamemode() result(resfglutentergamemode)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutentergamemode
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutentergamemode@0' :: fglutentergamemode
#endif
integer(kind=glcint) :: resfglutentergamemode
end function fglutentergamemode
end interface

interface
!NAS$ ALIEN "F77 fglutleavegamemode"
subroutine fglutleavegamemode()
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutleavegamemode
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutleavegamemode@0' :: fglutleavegamemode
#endif
end subroutine fglutleavegamemode
end interface

interface
!NAS$ ALIEN "F77 fglutgamemodeget"
function fglutgamemodeget(mode) result(resfglutgamemodeget)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mode
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutgamemodeget
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutgamemodeget@4' :: fglutgamemodeget
#endif
integer(kind=glenum) INTENT_IN mode
integer(kind=glcint) :: resfglutgamemodeget
end function fglutgamemodeget
end interface

#endif

interface
#ifdef NASOFTWARE
bind(c,"fglutbuttonboxfunc_") &
#endif
subroutine fglutbuttonboxfunc(buttonboxfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: buttonboxfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutbuttonboxfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutbuttonboxfunc@4' :: fglutbuttonboxfunc
#endif
interface
   C_BINDING subroutine buttonboxfunc(button,state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state
   end subroutine buttonboxfunc
end interface
end subroutine fglutbuttonboxfunc
end interface

interface
#ifdef PURE_FUNCTIONS
subroutine fglutcreatemenu(selectfunc,menuid)
#else
#ifdef NASOFTWARE
bind(c,"fglutcreatemenu_") &
#endif
function fglutcreatemenu(selectfunc) result(resfglutcreatemenu)
#endif
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: selectfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutcreatemenu
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutcreatemenu@4' :: fglutcreatemenu
#endif
interface
   C_BINDING subroutine selectfunc(value)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: value
   end subroutine selectfunc
end interface
#ifdef PURE_FUNCTIONS
integer(kind=glcint), intent(out) :: menuid
#else
integer(kind=glcint) :: resfglutcreatemenu
#endif
#ifdef PURE_FUNCTIONS
end subroutine fglutcreatemenu
#else
end function fglutcreatemenu
#endif
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutdialsfunc_") &
#endif
subroutine fglutdialsfunc(dialsfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: dialsfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutdialsfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutdialsfunc@4' :: fglutdialsfunc
#endif
interface
   C_BINDING subroutine dialsfunc(dial,value)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: dial, value
   end subroutine dialsfunc
end interface
end subroutine fglutdialsfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutdisplayfunc_") &
#endif
subroutine fglutdisplayfunc(displayfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: displayfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutdisplayfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutdisplayfunc@4' :: fglutdisplayfunc
#endif
interface
   C_BINDING subroutine displayfunc()
   use opengl_kinds
   IMPLICIT_NONE
   end subroutine displayfunc
end interface
end subroutine fglutdisplayfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutentryfunc_") &
#endif
subroutine fglutentryfunc(entryfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: entryfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutentryfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutentryfunc@4' :: fglutentryfunc
#endif
interface
   C_BINDING subroutine entryfunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine entryfunc
end interface
end subroutine fglutentryfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutidlefunc_") &
#endif
subroutine fglutidlefunc(idlefunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: idlefunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutidlefunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutidlefunc@4' :: fglutidlefunc
#endif
interface
   C_BINDING subroutine idlefunc()
   use opengl_kinds
   IMPLICIT_NONE
   end subroutine idlefunc
end interface
end subroutine fglutidlefunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutkeyboardfunc_") &
#endif
subroutine fglutkeyboardfunc(keyboardfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: keyboardfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutkeyboardfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutkeyboardfunc@4' :: fglutkeyboardfunc
#endif
interface
   C_BINDING subroutine keyboardfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint) PASS_BY_REF, intent(in out) :: key
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine keyboardfunc
end interface
end subroutine fglutkeyboardfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutmenustatefunc_") &
#endif
subroutine fglutmenustatefunc(menustatefunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: menustatefunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutmenustatefunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutmenustatefunc@4' :: fglutmenustatefunc
#endif
interface
   C_BINDING subroutine menustatefunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine menustatefunc
end interface
end subroutine fglutmenustatefunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutmenustatusfunc_") &
#endif
subroutine fglutmenustatusfunc(menustatusfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: menustatusfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutmenustatusfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutmenustatusfunc@4' :: fglutmenustatusfunc
#endif
interface
   C_BINDING subroutine menustatusfunc(status,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: status, x, y
   end subroutine menustatusfunc
end interface
end subroutine fglutmenustatusfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutmotionfunc_") &
#endif
subroutine fglutmotionfunc(motionfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: motionfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutmotionfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutmotionfunc@4' :: fglutmotionfunc
#endif
interface
   C_BINDING subroutine motionfunc(x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine motionfunc
end interface
end subroutine fglutmotionfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutmousefunc_") &
#endif
subroutine fglutmousefunc(mousefunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: mousefunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutmousefunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutmousefunc@4' :: fglutmousefunc
#endif
interface
   C_BINDING subroutine mousefunc(button,state,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state, x, y
   end subroutine mousefunc
end interface
end subroutine fglutmousefunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutoverlaydisplayfunc_") &
#endif
subroutine fglutoverlaydisplayfunc(overlaydisplayfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: overlaydisplayfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutoverlaydisplayfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutoverlaydisplayfunc@4' :: fglutoverlaydisplayfunc
#endif
interface
   C_BINDING subroutine overlaydisplayfunc()
   use opengl_kinds
   IMPLICIT_NONE
   end subroutine overlaydisplayfunc
end interface
end subroutine fglutoverlaydisplayfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutpassivemotionfunc_") &
#endif
subroutine fglutpassivemotionfunc(passivemotionfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: passivemotionfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutpassivemotionfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutpassivemotionfunc@4' :: fglutpassivemotionfunc
#endif
interface
   C_BINDING subroutine passivemotionfunc(x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine passivemotionfunc
end interface
end subroutine fglutpassivemotionfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutreshapefunc_") &
#endif
subroutine fglutreshapefunc(reshapefunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: reshapefunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutreshapefunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutreshapefunc@4' :: fglutreshapefunc
#endif
interface
   C_BINDING subroutine reshapefunc(width,height)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: width, height
   end subroutine reshapefunc
end interface
end subroutine fglutreshapefunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutspaceballbuttonfunc_") &
#endif
subroutine fglutspaceballbuttonfunc(spaceballbuttonfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: spaceballbuttonfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutspaceballbuttonfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutspaceballbuttonfunc@4' :: fglutspaceballbuttonfunc
#endif
interface
   C_BINDING subroutine spaceballbuttonfunc(button,state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state
   end subroutine spaceballbuttonfunc
end interface
end subroutine fglutspaceballbuttonfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutspaceballmotionfunc_") &
#endif
subroutine fglutspaceballmotionfunc(spaceballmotionfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: spaceballmotionfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutspaceballmotionfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutspaceballmotionfunc@4' :: fglutspaceballmotionfunc
#endif
interface
   C_BINDING subroutine spaceballmotionfunc(x,y,z)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y, z
   end subroutine spaceballmotionfunc
end interface
end subroutine fglutspaceballmotionfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutspaceballrotatefunc_") &
#endif
subroutine fglutspaceballrotatefunc(spaceballrotatefunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: spaceballrotatefunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutspaceballrotatefunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutspaceballrotatefunc@4' :: fglutspaceballrotatefunc
#endif
interface
   C_BINDING subroutine spaceballrotatefunc(x,y,z)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y, z
   end subroutine spaceballrotatefunc
end interface
end subroutine fglutspaceballrotatefunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutspecialfunc_") &
#endif
subroutine fglutspecialfunc(specialfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: specialfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutspecialfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutspecialfunc@4' :: fglutspecialfunc
#endif
interface
   C_BINDING subroutine specialfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: key, x, y
   end subroutine specialfunc
end interface
end subroutine fglutspecialfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fgluttabletbuttonfunc_") &
#endif
subroutine fgluttabletbuttonfunc(tabletbuttonfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tabletbuttonfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluttabletbuttonfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluttabletbuttonfunc@4' :: fgluttabletbuttonfunc
#endif
interface
   C_BINDING subroutine tabletbuttonfunc(button,state,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: button, state, x, y
   end subroutine tabletbuttonfunc
end interface
end subroutine fgluttabletbuttonfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fgluttabletmotionfunc_") &
#endif
subroutine fgluttabletmotionfunc(tabletmotionfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: tabletmotionfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluttabletmotionfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluttabletmotionfunc@4' :: fgluttabletmotionfunc
#endif
interface
   C_BINDING subroutine tabletmotionfunc(x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine tabletmotionfunc
end interface
end subroutine fgluttabletmotionfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fgluttimerfunc_") &
#endif
subroutine fgluttimerfunc(millis,timerfunc,value)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: millis,timerfunc,value
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fgluttimerfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fgluttimerfunc@12' :: fgluttimerfunc
#endif
integer(kind=glcuint) PASS_BY_REF INTENT_IN millis
interface
   C_BINDING subroutine timerfunc(value)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: value
   end subroutine timerfunc
end interface
integer(kind=glcint) PASS_BY_REF INTENT_IN value
end subroutine fgluttimerfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutvisibilityfunc_") &
#endif
subroutine fglutvisibilityfunc(visibilityfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: visibilityfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutvisibilityfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutvisibilityfunc@4' :: fglutvisibilityfunc
#endif
interface
   C_BINDING subroutine visibilityfunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine visibilityfunc
end interface
end subroutine fglutvisibilityfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutwindowstatusfunc_") &
#endif
subroutine fglutwindowstatusfunc(windowstatusfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: windowstatusfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutwindowstatusfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutwindowstatusfunc@4' :: fglutwindowstatusfunc
#endif
interface
   C_BINDING subroutine windowstatusfunc(state)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: state
   end subroutine windowstatusfunc
end interface
end subroutine fglutwindowstatusfunc
end interface

#ifdef GLUT37

interface
#ifdef NASOFTWARE
bind(c,"fglutkeyboardupfunc_") &
#endif
subroutine fglutkeyboardupfunc(keyboardupfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: keyboardupfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutkeyboardupfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutkeyboardupfunc@4' :: fglutkeyboardupfunc
#endif
interface
   C_BINDING subroutine keyboardupfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glint) PASS_BY_REF, intent(in out) :: key
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y
   end subroutine keyboardupfunc
end interface
end subroutine fglutkeyboardupfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutspecialupfunc_") &
#endif
subroutine fglutspecialupfunc(specialupfunc)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: specialupfunc
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutspecialupfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutspecialupfunc@4' :: fglutspecialupfunc
#endif
interface
   C_BINDING subroutine specialupfunc(key,x,y)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: key, x, y
   end subroutine specialupfunc
end interface
end subroutine fglutspecialupfunc
end interface

interface
#ifdef NASOFTWARE
bind(c,"fglutjoystickfunc_") &
#endif
subroutine fglutjoystickfunc(joystickfunc,pollinterval)
use opengl_kinds
IMPLICIT_NONE
#ifdef MSATTR_REF
!MS$ATTRIBUTES REFERENCE :: joystickfunc, pollinterval
#endif
#ifdef MSATTR_STD
!MS$ATTRIBUTES STDCALL :: fglutjoystickfunc
#endif
#ifdef MSATTR_ALIAS
!MS$ATTRIBUTES ALIAS : '_fglutjoystickfunc@8' :: fglutjoystickfunc
#endif
interface
   C_BINDING subroutine joystickfunc(buttonmask,x,y,z)
   use opengl_kinds
   IMPLICIT_NONE
   integer(kind=glcuint) PASS_BY_REF, intent(in out) :: buttonmask
   integer(kind=glcint) PASS_BY_REF, intent(in out) :: x, y, z
   end subroutine joystickfunc
end interface
integer(kind=glint) INTENT_IN pollinterval
end subroutine fglutjoystickfunc
end interface

#endif

#endif

end module opengl_glutinterfaces
