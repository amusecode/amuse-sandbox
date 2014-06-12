module bitfont_callbacks
use opengl_gl
use opengl_glu
use opengl_glut
private
public :: output, display, freshape

contains

subroutine output(x,y,s)
  real, intent(in) :: x,y
  character(len=*), intent(in) :: s
  character(len=1) :: c
  integer :: i,lenc
  
  call glrasterpos2f(x,y)
  lenc = len(s)
  do i=1,lenc
     c = s(i:i)
     call glutbitmapcharacter(GLUT_BITMAP_TIMES_ROMAN_24, &
          ichar(c))
  end do
  return
end subroutine output

subroutine display()
  call glclear(GL_COLOR_BUFFER_BIT)
  call output(0.0,24.0, &
         "This is written in a GLUT bitmap font.")
  call output(100.0,100.0,"ABCDEFGabcdefg")
  call output(50.0,145.0, &
       "(positioned in pixels with upper-left origin)")
  return
end subroutine display

subroutine freshape(w,h)
  integer, intent(in out) :: w, h
  call glviewport(0, 0, w, h)
  call glmatrixmode(GL_PROJECTION)
  call glloadidentity()
  call gluortho2d(0.0_gldouble, real(w,kind(0.0_gldouble)), 0.0_gldouble, &
                  real(h,kind(0.0_gldouble)))
  call glscalef(1.0, -1.0, 1.0)
  call gltranslatef(0.0_glfloat, real(-h,kind(0.0_glfloat)), 0.0_glfloat)
  call glmatrixmode(GL_MODELVIEW)
  return
end subroutine freshape

end module bitfont_callbacks

program main
  use bitfont_callbacks
  use opengl_gl
  use opengl_glut
  implicit none
  integer :: win
  call glutinitdisplaymode(GLUT_RGB + GLUT_SINGLE)
  call glutinitwindowsize(500, 150)
  call glutinit()
  win = glutcreatewindow("Fortran GLUT bitmap A")
  call glclearcolor(0.0_glfloat, 0.0_glfloat, 0.0_glfloat, 1.0_glfloat)
  call glutdisplayfunc(display)
  call glutreshapefunc(freshape)
  call glutmainloop()
  stop
end program main
