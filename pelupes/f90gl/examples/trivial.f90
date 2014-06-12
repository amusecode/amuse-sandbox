! This is a simple, introductory OpenGL program.  It originally comes from
! page 6 of the "red book".  This is a Fortran translation and slight
! modification of the program hello.c from the GLUT distribution.

module callbacks
implicit none
public :: display

contains

subroutine display()
use opengl_gl

! clear all pixels
call glclear(GL_COLOR_BUFFER_BIT)

! draw white polygon (rectangle) with corners at
! ((0.25, 0.25, 0.0) and (0.75, 0.75, 0.0)
call glcolor3f(1.0_glfloat, 1.0_glfloat, 1.0_glfloat)
call glbegin(GL_POLYGON)
   call glvertex3f(0.25_glfloat, 0.25_glfloat, 0.0_glfloat)
   call glvertex3f(0.75_glfloat, 0.25_glfloat, 0.0_glfloat)
   call glvertex3f(0.75_glfloat, 0.75_glfloat, 0.0_glfloat)
   call glvertex3f(0.25_glfloat, 0.75_glfloat, 0.0_glfloat)
call glend()

! don't wait!
! start processing buffered OpenGL routines
call glflush()

end subroutine display

end module callbacks

program trivial
use opengl_gl
use opengl_glut
use callbacks
implicit none

real(kind=glfloat),  parameter :: zero  = 0.0_glfloat
real(kind=gldouble), parameter :: dzero = 0.0_gldouble, &
                                  one   = 1.0_gldouble
integer(kind=glcint) :: iwin
integer(kind=glenum) :: itemp

! Declare initial window size, position, and display mode
! (single buffer and RGBA).  Open window with "trivial"
! in its title bar.

call glutinit()
call glutinitdisplaymode(ior(GLUT_SINGLE,GLUT_RGB))
iwin = glutcreatewindow("trivial")

! select clearing color
call glclearcolor(zero, zero, zero, zero)

! initialize viewing values
itemp = GL_PROJECTION
call glmatrixmode(itemp)
call glloadidentity()
call glortho(dzero, one, dzero, one, -one, one)

! Register callback function to display graphics.
! Enter main loop and process events.
call glutdisplayfunc(display)
call glutmainloop()

end program trivial
