!  blender renders two spinning icosahedrons (red and green).
!  The blending factors for the two icosahedrons vary sinusoidally
!  and slightly out of phase.  blender also renders two lines of
!  text in a stroke font: one line antialiased, the other not.

module blender_data

use opengl_gl
use opengl_glu
use opengl_glut
implicit none
private
public :: output, display, idle, visible

real(kind=glfloat), dimension(4), save, public :: light0_ambient = (/0.2, 0.2, 0.2, 1.0/)
real(kind=glfloat), dimension(4), save, public :: light0_diffuse = (/0.0, 0.0, 0.0, 1.0/)
real(kind=glfloat), dimension(4), save, public :: light1_diffuse = (/1.0, 0.0, 0.0, 1.0/)
real(kind=glfloat), dimension(4), save, public :: light1_position = (/1.0, 1.0, 1.0, 0.0/)
real(kind=glfloat), dimension(4), save, public :: light2_diffuse = (/0.0, 1.0, 0.0, 1.0/)
real(kind=glfloat), dimension(4), save, public :: light2_position = (/-1.0, -1.0, 1.0, 0.0/)
real, save, private :: s = 0.0
real(kind=glfloat), save, private :: angle1 = 0.0, angle2 = 0.0

contains

subroutine output(x, y, text)
real(kind=glfloat), intent(in) :: x,y
character(len=*), intent(in) :: text
integer(kind=glcint) :: p,i

  call glpushmatrix()
  call gltranslatef(x, y, 0.0_glfloat)
  do i=1,len(text)
    p = ichar(text(i:i))
    call glutstrokecharacter(GLUT_STROKE_ROMAN, p)
  end do
  call glpopmatrix()
  return
end subroutine output

subroutine display()

  real(kind=glfloat), save, dimension(4) :: amb = (/0.4, 0.4, 0.4, 0.0/)
  real(kind=glfloat), save, dimension(4) :: dif = (/1.0, 1.0, 1.0, 0.0/)

  call glclear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
  call glenable(GL_LIGHT1)
  call gldisable(GL_LIGHT2)
  dif(4) = cos(s) / 2.0 + 0.5
  amb(4) = dif(4)
  call glmaterialfv(GL_FRONT, GL_AMBIENT, amb)
  call glmaterialfv(GL_FRONT, GL_DIFFUSE, dif)

  call glpushmatrix()
! let's take a chance that the default integer is the same kind as
! glint, and not bother with the _glint on constants
  call gltranslatef(-0.3, -0.3, 0.0)
  call glrotatef(angle1, 1.0, 5.0, 0.0)
  call glcalllist(1)        ! render ico display list
  call glpopmatrix()

  call glclear(GL_DEPTH_BUFFER_BIT)
  call glenable(GL_LIGHT2)
  call gldisable(GL_LIGHT1)
  dif(4) = 0.5 - cos(s * 0.95) / 2.0
  amb(4) = dif(4)
  call glmaterialfv(GL_FRONT, GL_AMBIENT, amb)
  call glmaterialfv(GL_FRONT, GL_DIFFUSE, dif)

  call glpushmatrix()
  call gltranslatef(0.3, 0.3, 0.0)
  call glrotatef(angle2, 1.0, 0.0, 5.0)
  call glcalllist(1)        ! render ico display list
  call glpopmatrix()

  call glpushattrib(GL_ENABLE_BIT)
  call gldisable(GL_DEPTH_TEST)
  call gldisable(GL_LIGHTING)
  call glmatrixmode(GL_PROJECTION)
  call glpushmatrix()
  call glloadidentity()
  call gluortho2d(0.0_gldouble, 1500.0_gldouble, 0.0_gldouble, 1500.0_gldouble)
  call glmatrixmode(GL_MODELVIEW)
  call glpushmatrix()
  call glloadidentity()
!  Rotate text slightly to help show jaggies.
  call glrotatef(4.0, 0.0, 0.0, 1.0)
  call output(200.0, 225.0, "This is antialiased.")
  call glscalef(0.5,0.5,0.5)
  call gldisable(GL_LINE_SMOOTH)
  call gldisable(GL_BLEND)
  call output(160.0, 100.0, "This text is not.")
  call glpopmatrix()
  call glmatrixmode(GL_PROJECTION)
  call glpopmatrix()
  call glpopattrib()
  call glmatrixmode(GL_MODELVIEW)

  call glutswapbuffers()
  return
end subroutine display

subroutine idle()
  angle1 = modulo(angle1 + 0.8, 360.0)
  angle2 = modulo(angle2 + 1.1, 360.0)
  s = s + 0.05
  call glutpostredisplay()
  return
end subroutine idle

subroutine visible(vis)
integer, intent(in out) :: vis

  if (vis == GLUT_VISIBLE) then
    call glutidlefunc(idle)
  else
    call glutidlefunc()
  end if
  return
end subroutine visible

end module blender_data

program main
use blender_data
use opengl_gl
use opengl_glu
use opengl_glut
implicit none
integer(kind=glcint) :: i

  call glutinit()
  call glutinitdisplaymode(ior(ior(GLUT_DOUBLE,GLUT_RGB),GLUT_DEPTH))
  i = glutcreatewindow("blender")
  call glutdisplayfunc(display)
  call glutvisibilityfunc(visible)

  call glnewlist(1, GL_COMPILE)  ! create ico display list
  call glutsolidicosahedron()
  call glendlist()

  call glenable(GL_LIGHTING)
  call glenable(GL_LIGHT0)
  call gllightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient)
  call gllightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse)
  call gllightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse)
  call gllightfv(GL_LIGHT1, GL_POSITION, light1_position)
  call gllightfv(GL_LIGHT2, GL_DIFFUSE, light2_diffuse)
  call gllightfv(GL_LIGHT2, GL_POSITION, light2_position)
  call glenable(GL_DEPTH_TEST)
  call glenable(GL_CULL_FACE)
  call glenable(GL_BLEND)
  call glblendfunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
  call glenable(GL_LINE_SMOOTH)
  call gllinewidth(2.0)

  call glmatrixmode(GL_PROJECTION)
  call gluperspective( 40.0_gldouble, & ! field of view in degree
                           1.0_gldouble, & ! aspect ratio
                           1.0_gldouble, & ! Z near
                          10.0_gldouble)   ! Z far
  call glmatrixmode(GL_MODELVIEW)
  call glulookat( &
     0.0_gldouble, 0.0_gldouble, 5.0_gldouble, & ! eye is at (0,0,5)
     0.0_gldouble, 0.0_gldouble, 0.0_gldouble, & ! center is at (0,0,0)
     0.0_gldouble, 1.0_gldouble, 0.0_gldouble)    ! up is in positive Y direction
  call gltranslatef(0.0, 0.6, -1.0)

  call glutmainloop()
  stop
end program main
