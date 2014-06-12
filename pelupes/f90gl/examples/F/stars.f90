module stars_mod
use opengl_gl
use opengl_glu
use opengl_glut
implicit none

private
public :: Init, Reshape, Key, Visible, Display
private :: mySin, myCos, NewStar, RotatePoint, MoveStars, StarPoint, &
           ShowStar, UpdateStars, ShowStars, Idle, rand

real(kind=gldouble), private, parameter :: M_PI = 3.14159265358979323846_gldouble

integer(kind=glenum), private, parameter :: &
  NORMAL = 0, &
  WEIRD = 1

integer(kind=glenum), private, parameter :: &
  STREAK = 0, &
  CIRCLE = 1

integer, private, parameter :: MAXSTARS = 400
integer, private, parameter :: MAXPOS = 10000
integer, private, parameter :: MAXWARP = 10
integer, private, parameter :: MAXANGLES = 6000

type, private :: starRec
  integer(kind=GLint) :: type
  real, dimension(0:1) :: x, y, z
  real :: offsetX, offsetY, offsetR, rotation
end type starRec

integer(kind=GLenum), public :: doubleBuffer
integer(kind=GLint), public :: windW = 300, windH = 300

integer(kind=GLenum), private :: flag = NORMAL
integer(kind=GLint), private :: starCount = MAXSTARS / 2
real, private :: speed = 1.0
integer(kind=GLint), private :: nitro = 0
type(starRec), private, dimension(0:MAXSTARS) :: stars
real, private, dimension(0:MAXANGLES) :: sinTable

contains

function mySin(angle) result(resSin)
real, intent(in) :: angle
real :: resSin
  resSin = sinTable(modulo(int(angle),MAXANGLES))
return
end function mySin

function myCos(angle) result(resCos)
real, intent(in) :: angle
real :: resCos
   resCos = sinTable(modulo(int(angle) + (MAXANGLES / 4),MAXANGLES))
return
end function myCos

subroutine NewStar(n, d)
integer(kind=glint), intent(in) :: n,d
integer :: resrand
  call rand(resrand)
  if (modulo(resrand,4) == 0) then
    stars(n)%type = CIRCLE
  else
    stars(n)%type = STREAK
  end if
  call rand(resrand)
  stars(n)%x(0) = modulo(resrand, MAXPOS) - MAXPOS / 2
  call rand(resrand)
  stars(n)%y(0) = modulo(resrand, MAXPOS) - MAXPOS / 2
  call rand(resrand)
  stars(n)%z(0) = modulo(resrand, MAXPOS) + d
  stars(n)%x(1) = stars(n)%x(0)
  stars(n)%y(1) = stars(n)%y(0)
  stars(n)%z(1) = stars(n)%z(0)
  call rand(resrand)
  if (modulo(resrand, 4) == 0 .and. flag == WEIRD) then
    call rand(resrand)
    stars(n)%offsetX = modulo(resrand, 100) - 100 / 2
    call rand(resrand)
    stars(n)%offsetY = modulo(resrand, 100) - 100 / 2
    call rand(resrand)
    stars(n)%offsetR = modulo(resrand, 25) - 25 / 2
  else
    stars(n)%offsetX = 0.0
    stars(n)%offsetY = 0.0
    stars(n)%offsetR = 0.0
  end if
return
end subroutine NewStar

subroutine RotatePoint(x, y, rotation)
real, intent(in out) :: x,y
real, intent(in) :: rotation
  real :: tmpX, tmpY

  tmpX = x * myCos(rotation) - y * mySin(rotation)
  tmpY = y * myCos(rotation) + x * mySin(rotation)
  x = tmpX
  y = tmpY
return
end subroutine RotatePoint

subroutine MoveStars()
  real :: offset
  integer(kind=GLint) :: n

  offset = speed * 60.0

  do n=0,starCount
    stars(n)%x(1) = stars(n)%x(0)
    stars(n)%y(1) = stars(n)%y(0)
    stars(n)%z(1) = stars(n)%z(0)
    stars(n)%x(0) = stars(n)%x(0) + stars(n)%offsetX
    stars(n)%y(0) = stars(n)%y(0) + stars(n)%offsetY
    stars(n)%z(0) = stars(n)%z(0) - offset
    stars(n)%rotation = stars(n)%rotation + stars(n)%offsetR
    if (stars(n)%rotation >= MAXANGLES) then
      stars(n)%rotation = 0.0
    end if
  end do
return
end subroutine MoveStars

subroutine StarPoint(n,respoint)
integer(kind=glint), intent(in) :: n
integer(kind=glenum), intent(out) :: respoint
  real :: x0, y0

  x0 = stars(n)%x(0) * windW / stars(n)%z(0)
  y0 = stars(n)%y(0) * windH / stars(n)%z(0)
  call RotatePoint(x0, y0, stars(n)%rotation)
  x0 = x0 + windW / 2.0
  y0 = y0 + windH / 2.0

  if (x0 >= 0.0 .and. x0 < windW .and. y0 >= 0.0 .and. y0 < windH) then
    respoint = GL_TRUE
  else
    respoint = GL_FALSE
  end if
return
end subroutine StarPoint

subroutine ShowStar(n)
integer(kind=glint), intent(in) :: n
  real :: x0, y0, x1, y1, width, x, y
  integer(kind=GLint) :: i

  x0 = stars(n)%x(0) * windW / stars(n)%z(0)
  y0 = stars(n)%y(0) * windH / stars(n)%z(0)
  call RotatePoint(x0, y0, stars(n)%rotation)
  x0 = x0 + windW / 2.0
  y0 = y0 + windH / 2.0

  if (x0 >= 0.0 .and. x0 < windW .and. y0 >= 0.0 .and. y0 < windH) then
    if (stars(n)%type == STREAK) then
      x1 = stars(n)%x(1) * windW / stars(n)%z(1)
      y1 = stars(n)%y(1) * windH / stars(n)%z(1)
      call RotatePoint(x1, y1, stars(n)%rotation)
      x1 = x1 + windW / 2.0
      y1 = y1 + windH / 2.0

      call glLineWidth(MAXPOS / 100.0 / stars(n)%z(0) + 1.0)
      call glColor3f(1.0, (MAXWARP - speed) / MAXWARP, (MAXWARP - speed) / MAXWARP)
      if (abs(x0 - x1) < 1.0 .and. abs(y0 - y1) < 1.0) then
        call glBegin(GL_POINTS)
        call glVertex2f(x0, y0)
        call glEnd()
      else
        call glBegin(GL_LINES)
        call glVertex2f(x0, y0)
        call glVertex2f(x1, y1)
        call glEnd()
      end if
    else
      width = MAXPOS / 10.0 / stars(n)%z(0) + 1.0
      call glColor3f(1.0, 0.0, 0.0)
      call glBegin(GL_POLYGON)
      do i=0,7
        x = x0 + width * myCos(i * MAXANGLES / 8.0)
        y = y0 + width * mySin(i * MAXANGLES / 8.0)
        call glVertex2f(x, y)
      end do
      call glEnd()
    end if
  end if
return
end subroutine ShowStar

subroutine UpdateStars()
  integer(kind=GLint) :: n
  integer(kind=glenum) :: respoint

  call glClear(GL_COLOR_BUFFER_BIT)

  do n=0,starCount
    if (stars(n)%z(0) > speed .or. (stars(n)%z(0) > 0.0 .and. speed < MAXWARP)) then
      call StarPoint(n,respoint)
      if (respoint == GL_FALSE) then
        call NewStar(n, MAXPOS)
      end if
    else
      call NewStar(n, MAXPOS)
    end if
  end do
return
end subroutine UpdateStars

subroutine ShowStars()
  integer(kind=GLint) :: n

  call glClear(GL_COLOR_BUFFER_BIT)

  do n=0,starCount
    if (stars(n)%z(0) > speed .or. (stars(n)%z(0) > 0.0 .and. speed < MAXWARP)) then
      call ShowStar(n)
    end if
  end do
return
end subroutine ShowStars

subroutine Init()
  real :: angle
  integer(kind=GLint) :: n

  call random_seed()

  do n=0,MAXSTARS
    call NewStar(n, 100)
  end do

  angle = 0.0
  do n=0,MAXANGLES
    sinTable(n) = sin(angle)
    angle = angle + M_PI / (MAXANGLES / 2.0)
  end do

  call glClearColor(0.0, 0.0, 0.0, 0.0)

  call glDisable(GL_DITHER)
return
end subroutine Init

subroutine Reshape(width, height)
integer(kind=glcint), intent(in out) :: width, height
  windW = width
  windH = height

  call glViewport(0, 0, windW, windH)

  call glMatrixMode(GL_PROJECTION)
  call glLoadIdentity()
  call gluOrtho2D(-0.5_gldouble, windW + 0.5_gldouble, -0.5_gldouble, windH + 0.5_gldouble)
  call glMatrixMode(GL_MODELVIEW)
return
end subroutine Reshape

!/* ARGSUSED1 */
subroutine Key(ikey, x, y)
integer(kind=glcint), intent(in out) :: ikey, x, y
  select case(ikey)
  case (ichar(" "))
    if (flag == NORMAL) then
       flag = WEIRD
    else
       flag = NORMAL
    end if
  case (ichar("t"))
    nitro = 1
  case (27)
    stop
  end select
return
end subroutine key

subroutine Idle()
  call MoveStars()
  call UpdateStars()
  if (nitro > 0) then
    speed = (nitro / 10.0) + 1.0
    if (speed > MAXWARP) then
      speed = MAXWARP
    end if
    nitro = nitro + 1
    if (nitro > MAXWARP * 10) then
      nitro = -nitro
    end if
  elseif (nitro < 0) then
    nitro = nitro + 1
    speed = (-nitro / 10.0) + 1.0
    if (speed > MAXWARP) then
      speed = MAXWARP
    end if
  end if
  call glutPostRedisplay()
return
end subroutine Idle

subroutine Display()
  call ShowStars()
  if (doubleBuffer == GL_TRUE) then
    call glutSwapBuffers()
  else
    call glFlush()
  end if
return
end subroutine Display

subroutine Visible(state)
integer(kind=glcint), intent(in out) :: state
  if (state == GLUT_VISIBLE) then
    call glutIdleFunc(Idle)
  else
    call glutIdleFunc()
  end if
return
end subroutine visible

!static void
!Args(int argc, char **argv)
!{
!  GLint i;
!
!  doubleBuffer = GL_TRUE;
!
!  for (i = 1; i < argc; i++) {
!    if (strcmp(argv(i), "-sb") == 0) {
!      doubleBuffer = GL_FALSE;
!    } else if (strcmp(argv(i), "-db") == 0) {
!      doubleBuffer = GL_TRUE;
!    }
!  }
!}

subroutine rand(resrand)
integer, intent(out) :: resrand
real :: frand
  call random_number(frand)
  resrand = 32768*frand
return
end subroutine rand

end module stars_mod

program stars_prog
use opengl_gl
use opengl_glut
use stars_mod
  integer(kind=GLenum) :: type
  integer(kind=glcint) :: win

  call glutInitWindowSize(windW, windH)
  call glutInit()
  doubleBuffer = GL_TRUE !  Args(argc, argv);

  type = GLUT_RGB
  if (doubleBuffer == GL_TRUE) then
     type = ior(type,GLUT_DOUBLE)
  else
     type = ior(type,GLUT_SINGLE)
  end if
  call glutInitDisplayMode(type)
  win = glutCreateWindow("Stars")

  call Init()

  call glutReshapeFunc(Reshape)
  call glutKeyboardFunc(Key)
  call glutVisibilityFunc(Visible)
  call glutDisplayFunc(Display)
  call glutMainLoop()
stop
end program stars_prog
