module view_modifier

! This module provides facilities to modify the view in an OpenGL window.
! The mouse buttons and keyboard arrow keys can be used to zoom, pan,
! rotate and change the scale.  A menu or submenu can be used to select which
! buttons perform which function and to reset the view to the initial settings.
! This is limited to one window.

! William F. Mitchell
! william.mitchell@nist.gov
! Mathematical and Computational Sciences Division
! National Institute of Standards and Technology
! April, 1998

! To use this module:
!
! 1) put a USE view_modifier statement in any program unit that calls a
!    procedure in this module
!
! 2) set the initial operation assignments, view and scale below the
!    "Initial configuration" comment below
!
! 3) call view_modifier_init after glutcreateWindow
!    This is a subroutine that returns integer(kind=glcint) menuid.  The menuid
!    is the ID returned by glutcreatemenu.  You can either use the view_modifier
!    menu as your menu by calling glutattachmenu immediately after
!    view_modifier_init, as in
!       call view_modifier_init(menuid)
!       call glutattachmenu(GLUT_RIGHT_BUTTON)
!    or by using the menuid to attach a submenu to your own menu, as in
!       call glutaddsubmenu("View Modifier",menuid)
!
! 4) in any callback functions that update the display, put
!       call reset_view
!    as the first executable statement
!
! Note that view_modifier_init sets the callback functions for glutmousefunc,
! glutmotionfunc and glutspecialFunc, so don't call these yourself
!
! The menu allows you to select what operation is attached to the left and
! middle mouse buttons and arrow keys, reset to the initial view, and quit.
! The right mouse button should be used for the menu.

use opengl_gl
use opengl_glu
use opengl_glut
implicit none
private
public :: view_modifier_init, reset_view
public :: display ! TEMP -- see comments immediately after this module
private :: cart2sphere, sphere2cart, cart3D_plus_cart3D, cart3D_minus_cart3D, &
           reset_to_init, mouse, motion, arrows, &
           menu_handler, set_left_button, set_middle_button, set_arrow_keys, &
           operator(+), operator(-)

integer(kind=glcint), parameter, private :: ZOOM = 1, PAN = 2, ROTATE = 3, SCALEX = 4, &
                      SCALEY = 5, SCALEZ = 6
integer(kind=glcint), parameter, private :: RESET = 10, QUIT = 11
real(kind=gldouble), parameter, private :: PI = 3.141592653589793_gldouble

type, private :: cart2D ! 2D cartesian coordinates
   real(kind=gldouble) :: x, y
end type cart2D

type, private :: cart3D ! 3D cartesian coordinates
   real(kind=gldouble) :: x, y, z
end type cart3D

type, private :: sphere3D ! 3D spherical coordinates
   real(kind=gldouble) :: theta, phi, rho
end type sphere3D

type(cart2D), save, private :: angle
type(cart3D), save, private :: shift
real(kind=gldouble), save, private :: xscale_factor, yscale_factor, zscale_factor
logical, save, private :: moving_left, moving_middle
type(cart2D), save, private :: begin_left, begin_middle

interface operator(+)
   module procedure cart3D_plus_cart3D
end interface
interface operator(-)
   module procedure cart3D_minus_cart3D
end interface

! ------- Initial configuration -------

! Set the initial operation performed by each button and the arrow keys.
! The operations are ZOOM, PAN, ROTATE, SCALEX, SCALEY, and SCALEZ

integer, save, private ::   left_button_func = ROTATE, &
                 middle_button_func = ZOOM, &
                     arrow_key_func = PAN

! Set the initial view as the point you are looking at, the point you are
! looking from, and the scale factors

type(cart3D), parameter, private :: &
     init_lookat = cart3D(0.0_gldouble, 0.0_gldouble, 0.0_gldouble), &
   init_lookfrom = cart3D(10.0_gldouble, -20.0_gldouble, 5.0_gldouble)

real(kind=gldouble), parameter, private :: &
   init_xscale_factor = 1.0_gldouble, &
   init_yscale_factor = 1.0_gldouble, &
   init_zscale_factor = 1.0_gldouble

! -------- end of Initial configuration ------

contains

!          -------------
subroutine reset_to_init()
!          -------------

! This resets the view to the initial configuration

type(sphere3D) :: slookfrom

slookfrom = cart2sphere(init_lookfrom-init_lookat)
angle%x = -180.0_gldouble*slookfrom%theta/PI - 90.0_gldouble
angle%y = -180.0_gldouble*slookfrom%phi/PI
shift%x = 0.0_gldouble
shift%y = 0.0_gldouble
shift%z = -slookfrom%rho
xscale_factor = init_xscale_factor
yscale_factor = init_yscale_factor
zscale_factor = init_zscale_factor

call glutpostredisplay()

return
end subroutine reset_to_init

!          ----------
subroutine reset_view()
!          ----------

! This routine resets the view to the current orientation and scale

call glmatrixmode(GL_MODELVIEW)
call glpopmatrix()
call glpushmatrix()
call gltranslated(shift%x, shift%y, shift%z)
call glrotated(angle%x, 0.0_gldouble, 0.0_gldouble, 1.0_gldouble)
call glrotated(angle%y, cos(PI*angle%x/180.0_gldouble), &
               -sin(PI*angle%x/180.0_gldouble), 0.0_gldouble)
call gltranslated(-init_lookat%x, -init_lookat%y, -init_lookat%z)
call glscaled(xscale_factor,yscale_factor,zscale_factor)

return
end subroutine reset_view

!          -----
subroutine mouse(button, state, x, y)
!          -----
integer(kind=glcint), intent(in out) :: button, state, x, y

! This gets called when a mouse button changes
 
  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_DOWN) then
    moving_left = .true.
    begin_left = cart2D(x,y)
  end if
  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_UP) then
    moving_left = .false.
  end if
  if (button == GLUT_MIDDLE_BUTTON .and. state == GLUT_DOWN) then
    moving_middle = .true.
    begin_middle = cart2D(x,y)
  end if
  if (button == GLUT_MIDDLE_BUTTON .and. state == GLUT_UP) then
    moving_middle = .false.
  end if
  return
end subroutine mouse

!          ------
subroutine motion(x, y)
!          ------
integer(kind=glcint), intent(in out) :: x, y

! This gets called when the mouse moves

integer :: button_function
type(cart2D) :: begin
real(kind=gldouble) :: factor

! Determine and apply the button function

if (moving_left) then
   button_function = left_button_func
   begin = begin_left
else if(moving_middle) then
   button_function = middle_button_func
   begin = begin_middle
end if

select case(button_function)
case (ZOOM)
   if (y < begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + 0.002_gldouble*(begin%y-y))
   else if (y > begin%y) then
      factor = 1.0_gldouble + 0.002_gldouble*(y-begin%y)
   else
      factor = 1.0_gldouble
   end if
   shift%z = factor*shift%z
case (PAN)
   shift%x = shift%x + 0.01*(x - begin%x)
   shift%y = shift%y - 0.01*(y - begin%y)
case (ROTATE)
   angle%x = angle%x + (x - begin%x)
   angle%y = angle%y + (y - begin%y)
case (SCALEX)
   if (y < begin%y) then
      factor = 1.0_gldouble + 0.002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + 0.002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   xscale_factor = xscale_factor * factor
case (SCALEY)
   if (y < begin%y) then
      factor = 1.0_gldouble + 0.002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + 0.002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   yscale_factor = yscale_factor * factor
case (SCALEZ)
   if (y < begin%y) then
      factor = 1.0_gldouble + 0.002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + 0.002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   zscale_factor = zscale_factor * factor
end select

! update private variables and redisplay

if (moving_left) then
   begin_left = cart2D(x,y)
else if(moving_middle) then
   begin_middle = cart2D(x,y)
end if

if (moving_left .or. moving_middle) then
   call glutpostredisplay()
end if

return
end subroutine motion

!          ------
subroutine arrows(key, x, y)
!          ------
integer(kind=glcint), intent(in out) :: key, x, y

! This routine handles the arrow key operations

real(kind=gldouble) :: factor

select case(arrow_key_func)
case(ZOOM)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble + 0.02_gldouble
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble/(1.0_gldouble + 0.02_gldouble)
   case default
      factor = 1.0_gldouble
   end select
   shift%z = factor*shift%z
case(PAN)
   select case(key)
   case(GLUT_KEY_LEFT)
      shift%x = shift%x - 0.02
   case(GLUT_KEY_RIGHT)
      shift%x = shift%x + 0.02
   case(GLUT_KEY_DOWN)
      shift%y = shift%y - 0.02
   case(GLUT_KEY_UP)
      shift%y = shift%y + 0.02
   end select
case(ROTATE)
   select case(key)
   case(GLUT_KEY_LEFT)
      angle%x = angle%x - 1.0_gldouble
   case(GLUT_KEY_RIGHT)
      angle%x = angle%x + 1.0_gldouble
   case(GLUT_KEY_DOWN)
      angle%y = angle%y + 1.0_gldouble
   case(GLUT_KEY_UP)
      angle%y = angle%y - 1.0_gldouble
   end select
case(SCALEX)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + 0.02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + 0.02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   xscale_factor = xscale_factor * factor
case(SCALEY)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + 0.02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + 0.02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   yscale_factor = yscale_factor * factor
case(SCALEZ)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + 0.02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + 0.02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   zscale_factor = zscale_factor * factor

end select
   
call glutpostredisplay()

return
end subroutine arrows

!          ------------
subroutine menu_handler(value)
!          ------------
integer(kind=glcint), intent(in out) :: value

! This routine handles the first level entries in the menu

select case(value)

case(RESET)
   call reset_to_init()
case(QUIT)
   stop

end select

return
end subroutine menu_handler

!          ---------------
subroutine set_left_button(value)
!          ---------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the left button as given by menu selection

left_button_func = value

return
end subroutine set_left_button

!          -----------------
subroutine set_middle_button(value)
!          -----------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the middle button as given by menu selection

middle_button_func = value

return
end subroutine set_middle_button

!          --------------
subroutine set_arrow_keys(value)
!          --------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the arrow keys as given by menu selection

arrow_key_func = value

return
end subroutine set_arrow_keys

!          ------------------
subroutine view_modifier_init(menuid)
!          ------------------
integer(kind=glcint), intent(out) :: menuid

! This initializes the view modifier variables and sets initial view.
! It should be called immediately after glutcreateWindow

integer(kind=glcint) :: button_left, button_middle, arrow_keys

! set the callback functions

call glutmousefunc(mouse)
call glutmotionfunc(motion)
call glutspecialfunc(arrows)

! create the menu

call glutcreatemenu(set_left_button,button_left)
call glutaddmenuentry("rotate",ROTATE)
call glutaddmenuentry("zoom",ZOOM)
call glutaddmenuentry("pan",PAN)
call glutaddmenuentry("scale x",SCALEX)
call glutaddmenuentry("scale y",SCALEY)
call glutaddmenuentry("scale z", SCALEZ)
call glutcreatemenu(set_middle_button,button_middle)
call glutaddmenuentry("rotate",ROTATE)
call glutaddmenuentry("zoom",ZOOM)
call glutaddmenuentry("pan",PAN)
call glutaddmenuentry("scale x",SCALEX)
call glutaddmenuentry("scale y",SCALEY)
call glutaddmenuentry("scale z", SCALEZ)
call glutcreatemenu(set_arrow_keys,arrow_keys)
call glutaddmenuentry("rotate",ROTATE)
call glutaddmenuentry("zoom",ZOOM)
call glutaddmenuentry("pan",PAN)
call glutaddmenuentry("scale x",SCALEX)
call glutaddmenuentry("scale y",SCALEY)
call glutaddmenuentry("scale z", SCALEZ)
call glutcreatemenu(menu_handler,menuid)
call glutaddsubmenu("left mouse button",button_left)
call glutaddsubmenu("middle mouse button",button_middle)
call glutaddsubmenu("arrow keys",arrow_keys)
call glutaddmenuentry("reset to initial view",RESET)
call glutaddmenuentry("quit",QUIT)

! set the perspective

call glmatrixmode(GL_PROJECTION)
call gluperspective(10.0_gldouble, 1.0_gldouble, 0.1_gldouble, 200.0_gldouble)

! set the initial view

call glpushmatrix()
call reset_to_init()

return
end subroutine view_modifier_init

!        -----------
function sphere2cart(spoint) result(cpoint)
!        -----------
type(sphere3D), intent(in) :: spoint
type(cart3D) :: cpoint

! This converts a 3D point from spherical to cartesean coordinates

real(kind=gldouble) :: t,p,r

t=spoint%theta
p=spoint%phi
r=spoint%rho

cpoint%x = r*cos(t)*sin(p)
cpoint%y = r*sin(t)*sin(p)
cpoint%z = r*cos(p)

return
end function sphere2cart

!        -----------
function cart2sphere(cpoint) result(spoint)
!        -----------
type(cart3D), intent(in) :: cpoint
type(sphere3D) :: spoint

! This converts a 3D point from cartesean to spherical coordinates

real(kind=gldouble) :: x,y,z

x=cpoint%x
y=cpoint%y
z=cpoint%z

spoint%rho = sqrt(x*x+y*y+z*z)
if (x==0.0_gldouble .and. y==0.0_gldouble) then
   spoint%theta = 0.0_gldouble
else
   spoint%theta = atan2(y,x)
end if
if (spoint%rho == 0.0_gldouble) then
   spoint%phi = 0.0_gldouble
else
   spoint%phi = acos(z/spoint%rho)
end if

return
end function cart2sphere

!        ------------------
function cart3D_plus_cart3D(cart1,cart2) result(cart3)
!        ------------------
type(cart3D), intent(in) :: cart1, cart2
type(cart3D) :: cart3

! Compute the sum of two 3D cartesean points

cart3%x = cart1%x + cart2%x
cart3%y = cart1%y + cart2%y
cart3%z = cart1%z + cart2%z

return
end function cart3D_plus_cart3D

!        -------------------
function cart3D_minus_cart3D(cart1,cart2) result(cart3)
!        -------------------
type(cart3D), intent(in) :: cart1, cart2
type(cart3D) :: cart3

! Compute the difference of two 3D cartesean points

cart3%x = cart1%x - cart2%x
cart3%y = cart1%y - cart2%y
cart3%z = cart1%z - cart2%z

return
end function cart3D_minus_cart3D

! This doesn't belong in this module -- see below.
subroutine display()

! This gets called when the display needs to be redrawn

call reset_view()

call glclear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
call glcalllist(1)
call glutswapbuffers()

return
end subroutine display

end module view_modifier

!--------------------------------------------------------------------------

! This is a simple program to demonstrate the use of the view_modifier module
! It consists of a module with the callback functions and a main program.

! The F compiler on Linux complained about not finding this module. Because
! of lack of time I'm just going to move subroutine display to the
! view_modifier module instead of tracking down the problem.
!module view_demo_callbacks
!use opengl_gl
!use opengl_glut
!use view_modifier
!private
!public :: display
!
!contains
!
!subroutine display()
!
!! This gets called when the display needs to be redrawn
!
!call reset_view()
!
!call glclear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
!call glcalllist(1)
!call glutswapbuffers()
!
!return
!end subroutine display
!
!end module view_demo_callbacks

program view_demo

use opengl_gl
use opengl_glut
use view_modifier
!use view_demo_callbacks
implicit none

integer :: winid, menuid
real(kind=glfloat), dimension(3), parameter :: & ! colors for bronze from Redbook teapots
   ambient = (/ 0.2125_glfloat, 0.1275_glfloat, 0.054_glfloat /), &
   diffuse = (/ 0.714_glfloat, 0.4284_glfloat, 0.18144_glfloat /), &
  specular = (/ 0.393548_glfloat, 0.271906_glfloat, 0.166721_glfloat /)
real(kind=glfloat), dimension(4), parameter :: &
       pos = (/ 1.0_glfloat, 1.0_glfloat, 1.0_glfloat, 0.0_glfloat /), &
     white = (/ 1.0_glfloat, 1.0_glfloat, 1.0_glfloat, 1.0_glfloat /)

! Initializations

call glutinit()
call glutinitdisplaymode(ior(GLUT_DOUBLE,ior(GLUT_RGB,GLUT_DEPTH)))

! Create a window

winid = glutcreatewindow("View Modifier Demo")
call view_modifier_init(menuid)
call glutattachmenu(GLUT_RIGHT_BUTTON)

! Set the display callback

call glutdisplayfunc(display)

! Create the image

call glnewlist(1,GL_COMPILE)

! Draw axes so we know the orientation

call glbegin(GL_LINES)
call glvertex3f(0.0_glfloat,0.0_glfloat,0.0_glfloat)
call glvertex3f(2.0_glfloat,0.0_glfloat,0.0_glfloat)
call glvertex3f(0.0_glfloat,0.0_glfloat,0.0_glfloat)
call glvertex3f(0.0_glfloat,2.0_glfloat,0.0_glfloat)
call glvertex3f(0.0_glfloat,0.0_glfloat,0.0_glfloat)
call glvertex3f(0.0_glfloat,0.0_glfloat,2.0_glfloat)

! Draw crude x, y and z to label the axes

call glvertex3f(2.1_glfloat,-0.1_glfloat,0.1_glfloat) ! X
call glvertex3f(2.1_glfloat,0.1_glfloat,-0.1_glfloat)
call glvertex3f(2.1_glfloat,-0.1_glfloat,-0.1_glfloat)
call glvertex3f(2.1_glfloat,0.1_glfloat,0.1_glfloat)
call glvertex3f(0.1_glfloat,2.1_glfloat,0.1_glfloat) ! Y
call glvertex3f(0.0_glfloat,2.1_glfloat,0.0_glfloat)
call glvertex3f(-0.1_glfloat,2.1_glfloat,0.1_glfloat)
call glvertex3f(0.1_glfloat,2.1_glfloat,-0.1_glfloat)
call glvertex3f(-0.1_glfloat,0.1_glfloat,2.1_glfloat) ! Z
call glvertex3f(0.1_glfloat,0.1_glfloat,2.1_glfloat)
call glvertex3f(0.1_glfloat,0.1_glfloat,2.1_glfloat)
call glvertex3f(-0.1_glfloat,-0.1_glfloat,2.1_glfloat)
call glvertex3f(-0.1_glfloat,-0.1_glfloat,2.1_glfloat)
call glvertex3f(0.1_glfloat,-0.1_glfloat,2.1_glfloat)
call glend()

! Draw a teapot

! rotate so the z-axis comes out the top, x-axis out the spout
call glrotated(90.0_gldouble,1.0_gldouble,0.0_gldouble,0.0_gldouble)
call glmaterialfv(GL_FRONT, GL_AMBIENT, ambient)
call glmaterialfv(GL_FRONT, GL_DIFFUSE, diffuse)
call glmaterialfv(GL_FRONT, GL_SPECULAR, specular)
call glmaterialf(GL_FRONT, GL_SHININESS, 25.6_glfloat)
call glutsolidteapot(1.0_gldouble)

call glendlist()

! Set the lighting

call glclearcolor(0.8_glclampf, 0.8_glclampf, 0.8_glclampf, 1.0_glclampf)
call gllightfv(GL_LIGHT0, GL_DIFFUSE, white)
call gllightfv(GL_LIGHT0, GL_POSITION, pos)
call glenable(GL_LIGHTING)
call glenable(GL_LIGHT0)
call glenable(GL_DEPTH_TEST)

! Let glut take over

call glutmainloop()

stop
end program view_demo
