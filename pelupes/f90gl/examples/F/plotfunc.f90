! This program plots a function of two variables.  The function, called
! func_to_plot, is in module function at the beginning of this file.
! We then have the same module used in the modview example, followed by
! another module for plotting the function, called function_plotter.
! You might want to change default initial settings in modules modview and
! function_plotter.

! William F. Mitchell
! william.mitchell@nist.gov
! Mathematical and Computational Sciences Division
! National Institute of Standards and Technology
! August, 1999

!---------------------------------------------------------------------------

module function

use opengl_gl
implicit none
private
public :: func_to_plot

contains

! The function to plot

function func_to_plot(x,y) result(resfunc)
real(kind=GLDOUBLE) :: resfunc
real(kind=GLDOUBLE), intent(in) :: x,y

resfunc = 0.5_GLDOUBLE * ( &
       cos(0.3_GLDOUBLE*sqrt(80.0_GLDOUBLE*x)-16.0_GLDOUBLE*y/3.0_GLDOUBLE)* &
       cos(16.0_GLDOUBLE*x/3.0_GLDOUBLE) + x-y )

end function func_to_plot

end module function

!---------------------------------------------------------------------------

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
     init_lookat = cart3D(0.5_gldouble, 0.5_gldouble, 0.0_gldouble), &
   init_lookfrom = cart3D(5.0_gldouble, 10.0_gldouble, 2.5_gldouble)

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

end module view_modifier

!---------------------------------------------------------------------------

module function_plotter
use opengl_gl
use opengl_glut
use view_modifier
use function
implicit none
private
public :: display, draw_func, menu_handler, make_menu
private :: surface_color_handler, contour_color_handler, get_rainbow, &
           normcrossprod, param_handler

! symbolic constants

integer, private, parameter :: surfgrid_toggle = 1, &
                      surfsolid_toggle = 2, &
                      contour_toggle = 3, &
                      quit_selected = 4

integer, private, parameter :: set_nx = 1, &
                      set_ny = 2, &
                      set_ncontour = 3, &
                      set_contour_val = 4, &
                      set_xrange = 5, &
                      set_yrange = 6, &
                      reset_params = 7

integer, private, parameter :: black_contour = 1, &
                      rainbow_contour = 2

integer, private, parameter :: white_surface = 1, &
                      red_surface = 2, &
                      rainbow_surface = 3

! Default initial settings

integer, private, parameter :: init_ngridx = 40, &
                      init_ngridy = 40, &
                      init_num_contour = 20, &
                      init_contour_color = black_contour, &
                      init_surface_color = rainbow_surface

real(kind=GLDOUBLE), private, parameter :: init_minx = 0.0_GLDOUBLE, &
                             init_maxx = 1.0_GLDOUBLE, &
                             init_miny = 0.0_GLDOUBLE, &
                             init_maxy = 1.0_GLDOUBLE

logical, private, parameter :: init_draw_surface_grid = .false., &
                      init_draw_surface_solid = .true., &
                      init_draw_contour = .true.

! Current settings

integer, private, save :: ngridx = init_ngridx, &
           ngridy = init_ngridy, &
           num_contour = init_num_contour, &
           contour_color = init_contour_color, &
           surface_color = init_surface_color

real(kind=GLDOUBLE), private, save :: minx = init_minx, &
                  maxx = init_maxx, &
                  miny = init_miny, &
                  maxy = init_maxy, &
                  minz = 0.0_GLDOUBLE, &
                  maxz = 0.0_GLDOUBLE

logical, private, save :: draw_surface_grid = init_draw_surface_grid, &
           draw_surface_solid = init_draw_surface_solid, &
           draw_contour = init_draw_contour, &
           contour_values_given = .false.

real(kind=GLDOUBLE), private, save, allocatable, dimension(:) :: actual_contours

contains

subroutine display()

! This gets called when the display needs to be redrawn

call reset_view()

call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
call glCallList(1)
call glutSwapBuffers()

return
end subroutine display

subroutine draw_func()
real(kind=GLDOUBLE), dimension(0:ngridx) :: gridx
real(kind=GLDOUBLE), dimension(0:ngridy) :: gridy
real(kind=GLDOUBLE), dimension(0:ngridy,2) :: zval
integer :: i,j,k,cont
real(kind=GLDOUBLE) :: x1,x2,x3,xt,y1,y2,y3,yt,z1,z2,z3,zt
real(kind=GLDOUBLE) :: frac,xcross1,xcross2,ycross1,ycross2
real(kind=GLDOUBLE), dimension(num_contour) :: contour_value
real(kind=GLFLOAT), dimension(4) :: color, &
                 red = (/1.0,0.0,0.0,1.0/), &
                 black = (/0.0,0.0,0.0,1.0/), &
                 white = (/1.0,1.0,1.0,1.0/)
real(kind=GLFLOAT), dimension(3) :: normal

! prepare to make a new display list

call reset_view()
call glDeleteLists(1_gluint, 1_glsizei)
call glNewList(1_gluint, gl_compile_and_execute)

! set the grid points

gridx = (/ ((minx + i*(maxx-minx)/ngridx),i=0,ngridx) /)
gridy = (/ ((miny + i*(maxy-miny)/ngridy),i=0,ngridy) /)

! if this is the first call and either rainbow coloring of a solid surface
! or contours are being drawn, minz and maxz need to be set

if (minz == 0.0_GLDOUBLE .and. maxz == 0.0_GLDOUBLE .and. &
    ( (draw_surface_solid .and. surface_color == rainbow_surface) .or. &
      draw_contour ) ) then
   do i=0,ngridx
      do j=0,ngridy
         z1 = func_to_plot(gridx(i),gridy(j))
         minz = min(z1,minz)
         maxz = max(z1,maxz)
      end do
   end do
endif

! draw the solid surface

if (draw_surface_solid) then

   call glPolygonMode(gl_front_and_back, gl_fill)
   call glBegin(gl_triangles)

! set the color for a red or white surface.  For white, set the lighting
! such that there is uniform brightness.

   if (surface_color == red_surface) then
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,red)
   elseif (surface_color == white_surface) then
      call glDisable(gl_light0)
      call glLightModelfv(gl_light_model_ambient, (/1.0,1.0,1.0,1.0/))
      call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,white)
   endif

! compute the function values for the first line of points

   do j=0,ngridy
      zval(j,2) = func_to_plot(gridx(0),gridy(j))
   end do

! for each x grid interval...

   do i=1,ngridx

! copy left side function values from the right side of the previous interval

      zval(:,1) = zval(:,2)

! compute the function values for the right side of the interval

      do j=0,ngridy
         zval(j,2) = func_to_plot(gridx(i),gridy(j))
      end do

      minz = min(minz,minval(zval))
      maxz = max(maxz,maxval(zval))

! for each y grid interval ...

      do j=1,ngridy

! add two triangles to the display list
! Before each triangle, set the normal.  Before each vertex, set the
! color if we're coloring by height

         normal = normcrossprod((/gridx(i-1),gridx(i),gridx(i)/), &
                                (/gridy(j-1),gridy(j-1),gridy(j)/), &
                                (/zval(j-1,1),zval(j-1,2),zval(j,2)/))
         call glNormal3fv(normal)
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j-1,1),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i-1),gridy(j-1),zval(j-1,1))
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j-1,2),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i  ),gridy(j-1),zval(j-1,2))
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j,2),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i  ),gridy(j  ),zval(j  ,2))
         normal = normcrossprod((/gridx(i-1),gridx(i-1),gridx(i)/), &
                                (/gridy(j),gridy(j-1),gridy(j)/), &
                                (/zval(j,1),zval(j-1,1),zval(j,2)/))
         call glNormal3fv(normal)
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j,2),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i  ),gridy(j  ),zval(j  ,2))
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j-1,1),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i-1),gridy(j-1),zval(j-1,1))
         if (surface_color == rainbow_surface) then
            call get_rainbow(zval(j,1),minz,maxz,color)
            call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,color)
         endif
         call glvertex3d(gridx(i-1),gridy(j  ),zval(j  ,1))

      end do
   end do

   call glEnd()

! if the surface is white, reset the lighting conditions

   if (surface_color == white_surface) then
      call glEnable(gl_light0)
      call glLightModelfv(gl_light_model_ambient, (/0.5,0.5,0.5,1.0/))
   endif
   

endif ! draw_surface_solid

! draw the surface grid

if (draw_surface_grid) then

   call glPolygonMode(gl_front_and_back, gl_line)
   call glBegin(gl_triangles)

! draw surface grid in black

   call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,black)

! compute the function values for the first line of points

   do j=0,ngridy
      zval(j,2) = func_to_plot(gridx(0),gridy(j))
   end do

! for each x grid interval...

   do i=1,ngridx

! copy left side function values from the right side of the previous interval

      zval(:,1) = zval(:,2)

! compute the function values for the right side of the interval

      do j=0,ngridy
         zval(j,2) = func_to_plot(gridx(i),gridy(j))
      end do

      minz = min(minz,minval(zval))
      maxz = max(maxz,maxval(zval))

! for each y grid interval ...

      do j=1,ngridy

! add two triangles to the display list

         call glvertex3d(gridx(i-1),gridy(j-1),zval(j-1,1))
         call glvertex3d(gridx(i  ),gridy(j-1),zval(j-1,2))
         call glvertex3d(gridx(i  ),gridy(j  ),zval(j  ,2))
         call glvertex3d(gridx(i  ),gridy(j  ),zval(j  ,2))
         call glvertex3d(gridx(i-1),gridy(j-1),zval(j-1,1))
         call glvertex3d(gridx(i-1),gridy(j  ),zval(j  ,1))

      end do
   end do

   call glEnd()

endif ! draw_surface_grid

! draw the contour plot

if (draw_contour) then

   call glPolygonMode(gl_front_and_back, gl_line)
   call glBegin(gl_lines)
   call glNormal3fv((/0.0_glfloat, 0.0_glfloat, 1.0_glfloat/))

! draw the domain in black, also sets color to black for black_contour

   call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,black)
   call glVertex3d(minx,miny,0.0_GLDOUBLE)
   call glVertex3d(maxx,miny,0.0_GLDOUBLE)
   call glVertex3d(maxx,miny,0.0_GLDOUBLE)
   call glVertex3d(maxx,maxy,0.0_GLDOUBLE)
   call glVertex3d(maxx,maxy,0.0_GLDOUBLE)
   call glVertex3d(minx,maxy,0.0_GLDOUBLE)
   call glVertex3d(minx,maxy,0.0_GLDOUBLE)
   call glVertex3d(minx,miny,0.0_GLDOUBLE)

! set the contour values

   if (contour_values_given) then
      contour_value = actual_contours
   else
      do i=1,num_contour
         contour_value(i) = minz+(maxz-minz)*(i-1)/real(num_contour-1,GLDOUBLE)
      end do
   endif

! compute the function values for the first line of points

   do j=0,ngridy
      zval(j,2) = func_to_plot(gridx(0),gridy(j))
   end do

! for each x grid interval...

   do i=1,ngridx

! copy left side function values from the right side of the previous interval

      zval(:,1) = zval(:,2)

! compute the function values for the right side of the interval

      do j=0,ngridy
         zval(j,2) = func_to_plot(gridx(i),gridy(j))
      end do

      minz = min(minz,minval(zval))
      maxz = max(maxz,maxval(zval))

! for each y grid interval ...

      do j=1,ngridy

! for two triangles

         do k=1,2

! set the vertices

            if (k==1) then
               x1 = gridx(i-1)
               y1 = gridy(j-1)
               z1 = zval(j-1,1)
               x2 = gridx(i  )
               y2 = gridy(j-1)
               z2 = zval(j-1,2)
               x3 = gridx(i  )
               y3 = gridy(j  )
               z3 = zval(j  ,2)
            else
               x1 = gridx(i-1)
               y1 = gridy(j-1)
               z1 = zval(j-1,1)
               x2 = gridx(i-1)
               y2 = gridy(j  )
               z2 = zval(j  ,1)
               x3 = gridx(i  )
               y3 = gridy(j  )
               z3 = zval(j  ,2)
            endif

! order the vertices by z value

            xt = x1
            yt = y1
            zt = z1
            if (z2 < z1) then
               xt = x1
               yt = y1
               zt = z1
               if (z3 < z1) then
                  if (z2 < z3) then
                     x1 = x2
                     y1 = y2
                     z1 = z2
                     x2 = x3
                     y2 = y3
                     z2 = z3
                  else
                     x1 = x3
                     y1 = y3
                     z1 = z3
                  endif
                  x3 = xt
                  y3 = yt
                  z3 = zt
               else
                  x1 = x2
                  y1 = y2
                  z1 = z2
                  x2 = xt
                  y2 = yt
                  z2 = zt
               endif
            elseif (z3 < z1) then
               x1 = x3
               y1 = y3
               z1 = z3
               x3 = x2
               y3 = y2
               z3 = z2
               x2 = xt
               y2 = yt
               z2 = zt
            elseif (z3 < z2) then
               xt = x2
               yt = y2
               zt = z2
               x2 = x3
               y2 = y3
               z2 = z3
               x3 = xt
               y3 = yt
               z3 = zt
            endif

! if z1==z3, the triangle is horizontal and no contours pass through it

            if (z1==z3) cycle

! for each contour value

            do cont = 1,num_contour

! see if it passes through this triangle

               if (contour_value(cont) < z1) cycle
               if (contour_value(cont) > z3) exit

! set the color for contours colored by solution

               if (contour_color == rainbow_contour) then
                  call get_rainbow(contour_value(cont),minz,maxz,color)
                  call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse, &
                                    color)
               endif

! see where it crosses the 1-3 edge

               frac = (contour_value(cont)-z1)/(z3-z1)
               xcross1 = (1.0_GLDOUBLE - frac)*x1 + frac*x3
               ycross1 = (1.0_GLDOUBLE - frac)*y1 + frac*y3

! see where it crosses one of the other edges

               if (contour_value(cont) == z2) then
                  xcross2 = x2
                  ycross2 = y2
               elseif (contour_value(cont) < z2) then
                  frac = (contour_value(cont)-z1)/(z2-z1)
                  xcross2 = (1.0_GLDOUBLE - frac)*x1 + frac*x2
                  ycross2 = (1.0_GLDOUBLE - frac)*y1 + frac*y2
               else
                  frac = (contour_value(cont)-z2)/(z3-z2)
                  xcross2 = (1.0_GLDOUBLE - frac)*x2 + frac*x3
                  ycross2 = (1.0_GLDOUBLE - frac)*y2 + frac*y3
               endif

! add the line segment to the display list

               call glVertex3d(xcross1,ycross1,0.0_GLDOUBLE)
               call glVertex3d(xcross2,ycross2,0.0_GLDOUBLE)

            end do
         end do
      end do
   end do

   call glEnd()

endif ! draw_contour

! finish off display list

call glEndList()
call glutPostRedisplay()

end subroutine draw_func

subroutine get_rainbow(val,minval,maxval,c)
real(kind=GLDOUBLE), intent(in) :: val,maxval,minval
real(kind=GLFLOAT), intent(out), dimension(4) :: c

real(kind=GLFLOAT) :: f

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = 0.5_glfloat
endif

if (f < 0.25) then
   c(1) = 0.0_glfloat
   c(2) = 4.0_glfloat * f
   c(3) = 1.0_glfloat
   c(4) = 1.0_glfloat
elseif (f < 0.5) then
   c(1) = 0.0_glfloat
   c(2) = 1.0_glfloat
   c(3) = 2.0_glfloat - 4.0_glfloat*f
   c(4) = 1.0_glfloat
elseif (f < 0.75) then
   c(1) = 4.0_glfloat * f - 2.0_glfloat
   c(2) = 1.0_glfloat
   c(3) = 0.0_glfloat
   c(4) = 1.0_glfloat
else
   c(1) = 1.0_glfloat
   c(2) = 4.0_glfloat - 4.0_glfloat*f
   c(3) = 0.0_glfloat
   c(4) = 1.0_glfloat
endif

end subroutine get_rainbow

function normcrossprod(x,y,z) result(resncp)
real(kind=glfloat), dimension(3) :: resncp
real(kind=gldouble), dimension(3), intent(in) :: x,y,z
real(kind=glfloat), dimension(3) :: t1,t2
real(kind=glfloat) :: norm
t1(1) = x(2) - x(1)
t1(2) = y(2) - y(1)
t1(3) = z(2) - z(1)
t2(1) = x(3) - x(1)
t2(2) = y(3) - y(1)
t2(3) = z(3) - z(1)
resncp(1) = t1(2)*t2(3) - t1(3)*t2(2)
resncp(2) = t1(3)*t2(1) - t1(1)*t2(3)
resncp(3) = t1(1)*t2(2) - t1(2)*t2(1)
norm = sqrt(dot_product(resncp,resncp))
if (norm /= 0.0_glfloat) resncp = resncp/norm
end function normcrossprod

subroutine menu_handler(selection)
integer(kind=glcint), intent(in out) :: selection

select case (selection)

case (surfgrid_toggle)
   draw_surface_grid = .not. draw_surface_grid
   call draw_func()

case (surfsolid_toggle)
   draw_surface_solid = .not. draw_surface_solid
   call draw_func()

case (contour_toggle)
   draw_contour = .not. draw_contour
   call draw_func()

case (quit_selected)
   stop

end select

return
end subroutine menu_handler

subroutine param_handler(selection)
integer(kind=glcint), intent(in out) :: selection

select case (selection)

case (set_nx)
   print *,"Enter number of x subintervals:"
   read *, ngridx
   call draw_func()

case (set_ny)
   print *,"Enter number of y subintervals:"
   read *, ngridy
   call draw_func()

case (set_ncontour)
   print *,"Enter number of contour lines:"
   read *, num_contour
   contour_values_given = .false.
   call draw_func()

case (set_contour_val)
   print *,"enter number of contours:"
   read *, num_contour
   if (allocated(actual_contours)) deallocate(actual_contours)
   allocate(actual_contours(num_contour))
   print *,"enter ",num_contour," contour values:"
   read *,actual_contours
   contour_values_given = .true.
   call draw_func()

case (set_xrange)
   print *,"Enter minimum and maximum x values:"
   read *,minx,maxx
   call draw_func()

case (set_yrange)
   print *,"Enter minimum and maximum y values:"
   read *,miny,maxy
   call draw_func()

case (reset_params)
   ngridx = init_ngridx
   ngridy = init_ngridy
   num_contour = init_num_contour
   contour_color = init_contour_color
   surface_color = init_surface_color
   minx = init_minx
   maxx = init_maxx
   miny = init_miny
   maxy = init_maxy
   draw_surface_grid = init_draw_surface_grid
   draw_surface_solid = init_draw_surface_solid
   draw_contour = init_draw_contour
   call draw_func()

end select

end subroutine param_handler

subroutine contour_color_handler(selection)
integer(kind=glcint), intent(in out) :: selection

contour_color = selection
call draw_func()

end subroutine contour_color_handler

subroutine surface_color_handler(selection)
integer(kind=glcint), intent(in out) :: selection

surface_color = selection
call draw_func()

end subroutine surface_color_handler

subroutine make_menu(submenuid)
integer, intent(in) :: submenuid
integer :: menuid, param_id, contour_color_menu, surface_color_menu

call glutCreateMenu(contour_color_handler,contour_color_menu)
call glutAddMenuEntry("black",black_contour)
call glutAddMenuEntry("contour value",rainbow_contour)

call glutCreateMenu(surface_color_handler,surface_color_menu)
call glutAddMenuEntry("red",red_surface)
call glutAddMenuEntry("white",white_surface)
call glutAddMenuEntry("surface height",rainbow_surface)

call glutCreateMenu(param_handler,param_id)
call glutAddMenuEntry("number of x grid intervals",set_nx)
call glutAddMenuEntry("number of y grid intervals",set_ny)
call glutAddMenuEntry("number of uniform contour lines",set_ncontour)
call glutAddMenuEntry("contour values",set_contour_val)
call glutAddMenuEntry("x range",set_xrange)
call glutAddMenuEntry("y range",set_yrange)
call glutAddSubMenu("contour color",contour_color_menu)
call glutAddSubMenu("surface color",surface_color_menu)
call glutAddMenuEntry("reset to initial parameters",reset_params)

call glutCreateMenu(menu_handler,menuid)
call glutAddSubMenu("View Modifier",submenuid)
call glutAddMenuEntry("toggle surface grid",surfgrid_toggle)
call glutAddMenuEntry("toggle solid surface",surfsolid_toggle)
call glutAddMenuEntry("toggle contour",contour_toggle)
call glutAddSubMenu("plotting parameters",param_id)
call glutAddMenuEntry("quit",quit_selected)
call glutAttachMenu(GLUT_RIGHT_BUTTON)
end subroutine make_menu

end module function_plotter

!---------------------------------------------------------------------------

program plot_func

use opengl_gl
use opengl_glut
use view_modifier
use function_plotter
implicit none

integer :: winid, menuid, submenuid

! Initializations

call glutInit()
call glutInitDisplayMode(ior(GLUT_DOUBLE,ior(GLUT_RGB,GLUT_DEPTH)))
call glutInitWindowSize(500_glcint,500_glcint)

! Create a window

winid = glutCreateWindow("Function plotter")

! initialize view_modifier, receiving the id for it's submenu

call view_modifier_init(submenuid)

! create the menu

call make_menu(submenuid)

! Set the display callback

call glutDisplayFunc(display)

! set the lighting conditions

call glClearColor(0.9_glclampf, 0.9_glclampf, 0.9_glclampf, 1.0_glclampf)
call glLightfv(gl_light0, gl_diffuse, (/1.0,1.0,1.0,1.0/))
call glLightfv(gl_light0, gl_position, (/1.5,-0.5,2.0,0.0/))
call glEnable(gl_lighting)
call glEnable(gl_light0)
call glLightModelfv(gl_light_model_ambient, (/0.5,0.5,0.5,1.0/))
call glDepthFunc(gl_lequal)
call glEnable(gl_depth_test)

! Create the image

call draw_func()

! Let glut take over

call glutMainLoop()

end program plot_func
