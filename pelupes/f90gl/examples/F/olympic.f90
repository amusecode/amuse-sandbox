! In addition to modifications for the F and ELF90 subsets, rgb_colors was
! changed from glubyte to glfloat

module olympic_mod
use opengl_gl
use opengl_glu
use opengl_glut
implicit none
private

public :: fReshape, Key, visible, DrawScene, Init
private :: FillTorus, Clamp, MyRand, ReInit, Idle

integer, parameter, private :: &
   XSIZE      = 100, &
   YSIZE      = 75, &
   RINGS      = 5, &
   BLUERING   = 0, &
   BLACKRING  = 1, &
   REDRING    = 2, &
   YELLOWRING = 3, &
   GREENRING  = 4, &
   BLACK      = 0, &
   RED        = 1, &
   GREEN      = 2, &
   YELLOW     = 3, &
   BLUE       = 4, &
   MAGENTA    = 5, &
   CYAN       = 6, &
   WHITE      = 7

real(kind=glfloat), parameter, private :: BACKGROUND = 8.0

integer, parameter, private :: double = kind(0.0_gldouble)

real, parameter, private :: M_PI = 3.141592654

integer(kind=glenum), public :: rgb, doubleBuffer

!integer(kind=glubyte), private, dimension(0:RINGS-1,0:2) :: rgb_colors
real(kind=glfloat),   private, dimension(0:RINGS-1,0:2) :: rgb_colors
integer(kind=glint), private, dimension(0:RINGS-1) :: mapped_colors
real(kind=glfloat), private, dimension(0:RINGS-1,0:2) :: dests
real(kind=glfloat), private, dimension(0:RINGS-1,0:2) :: offsets
real(kind=glfloat), private, dimension(0:RINGS-1) :: angs
real(kind=glfloat), private, dimension(0:RINGS-1,0:2) :: rotAxis
integer, private, dimension(0:RINGS-1) :: iters
integer(kind=gluint), private :: theTorus

contains

subroutine FillTorus(rc, numc, rt, numt)
real, intent(in) :: rc, rt
integer, intent(in) :: numc, numt

  integer :: i, j, k
  real :: s, t
  real(kind=glfloat) :: x, y, z
  real :: pi, twopi

  pi = M_PI
  twopi = 2 * pi

  do i = 0, numc-1
    call glbegin(GL_QUAD_STRIP)
    do j = 0, numt
      do k = 1, 0, -1
        s = modulo((i + k), numc) + 0.5
        t = modulo(j, numt)

        x = cos(t * twopi / numt) * cos(s * twopi / numc)
        y = sin(t * twopi / numt) * cos(s * twopi / numc)
        z = sin(s * twopi / numc)
        call glnormal3f(x, y, z)

        x = (rt + rc * cos(s * twopi / numc)) * cos(t * twopi / numt)
        y = (rt + rc * cos(s * twopi / numc)) * sin(t * twopi / numt)
        z = rc * sin(s * twopi / numc)
        call glvertex3f(x, y, z)
      end do
    end do
    call glend()
  end do
return
end subroutine FillTorus

function Clamp(iters_left,t) result(resClamp)
integer, intent(in) :: iters_left
real, intent(in) :: t
real :: resClamp

  if (iters_left < 3) then
    resClamp = 0.0
  else
    resClamp = (iters_left - 2) * t / iters_left
  end if
return
end function Clamp

subroutine MyRand(resMyRand)
real, intent(out) :: resMyRand
real :: rval

  call random_number(rval)
  resMyRand = 10.0 * (rval - 0.5)

return
end subroutine MyRand

subroutine ReInit()

  integer :: i
  real :: deviation, temp

  call MyRand(deviation)
  deviation = deviation / 2
  deviation = deviation * deviation
  do i = 0, RINGS-1
    call MyRand(offsets(i,0))
    call MyRand(offsets(i,1))
    call MyRand(offsets(i,2))
    call MyRand(angs(i))
    angs(i) = 260.0 * angs(i)
    call MyRand(rotAxis(i,0))
    call MyRand(rotAxis(i,1))
    call MyRand(rotAxis(i,2))
    call MyRand(temp)
    iters(i) = (deviation * temp + 60.0)
  end do

return
end subroutine ReInit

subroutine Init()

  real(kind=glfloat), save :: top_y = 1.0
  real(kind=glfloat), save :: bottom_y = 0.0
  real(kind=glfloat), save :: top_z = 0.15
  real(kind=glfloat), save :: bottom_z = 0.69
  real(kind=glfloat), save :: fspacing = 2.5
  real(kind=glfloat), save, dimension(4) :: lmodel_ambient = (/0.0, 0.0, 0.0, 0.0/)
  real(kind=glfloat), save, dimension(1) :: lmodel_twoside = (/GL_FALSE/)
  real(kind=glfloat), save, dimension(1) :: lmodel_local = (/GL_FALSE/)
  real(kind=glfloat), save, dimension(4) :: light0_ambient = (/0.1, 0.1, 0.1, 1.0/)
  real(kind=glfloat), save, dimension(4) :: light0_diffuse = (/1.0, 1.0, 1.0, 0.0/)
  real(kind=glfloat), save, dimension(4) :: light0_position = (/0.8660254, 0.5, 1.0, 0.0/)
  real(kind=glfloat), save, dimension(4) :: light0_specular = (/1.0, 1.0, 1.0, 0.0/)
  real(kind=glfloat), save, dimension(4) :: bevel_mat_ambient = (/0.0, 0.0, 0.0, 1.0/)
  real(kind=glfloat), save, dimension(1) :: bevel_mat_shininess = (/40.0/)
  real(kind=glfloat), save, dimension(4) :: bevel_mat_specular = (/1.0, 1.0, 1.0, 0.0/)
  real(kind=glfloat), save, dimension(4) :: bevel_mat_diffuse = (/1.0, 0.0, 0.0, 0.0/)

  call random_seed()
  call ReInit()
! changed for glfloat rgb_colors
  rgb_colors = 0.0_glfloat
  rgb_colors(BLUERING,2) = 1.0_glfloat
  rgb_colors(REDRING,0) = 1.0_glfloat
  rgb_colors(GREENRING,1) = 1.0_glfloat
  rgb_colors(YELLOWRING,0) = 1.0_glfloat
  rgb_colors(YELLOWRING,1) = 1.0_glfloat
  mapped_colors(BLUERING) = BLUE
  mapped_colors(REDRING) = RED
  mapped_colors(GREENRING) = GREEN
  mapped_colors(YELLOWRING) = YELLOW
  mapped_colors(BLACKRING) = BLACK

  dests(BLUERING,:) = (/-fspacing, top_y, top_z/)
  dests(BLACKRING,:) = (/0.0, top_y, top_z/)
  dests(REDRING,:) = (/fspacing, top_y, top_z/)
  dests(YELLOWRING,:) = (/-fspacing / 2.0, bottom_y, bottom_z/)
  dests(GREENRING,:) = (/fspacing / 2.0, bottom_y, bottom_z/)

  theTorus = glgenlists(1)
  call glnewlist(theTorus, GL_COMPILE)
  call FillTorus(0.1, 8, 1.0, 25)
  call glendlist()

  call glenable(GL_CULL_FACE)
  call glcullface(GL_BACK)
  call glenable(GL_DEPTH_TEST)
  call glcleardepth(1.0_glclampd)

  if (rgb == GL_TRUE) then
    call glclearcolor(0.5_glclampf, 0.5_glclampf, 0.5_glclampf, 0.0_glclampf)
    call gllightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient)
    call gllightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse)
    call gllightfv(GL_LIGHT0, GL_SPECULAR, light0_specular)
    call gllightfv(GL_LIGHT0, GL_POSITION, light0_position)
    call glenable(GL_LIGHT0)

    call gllightmodelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, lmodel_local)
    call gllightmodelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside)
    call gllightmodelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient)
    call glenable(GL_LIGHTING)

    call glmaterialfv(GL_FRONT, GL_AMBIENT, bevel_mat_ambient)
    call glmaterialfv(GL_FRONT, GL_SHININESS, bevel_mat_shininess)
    call glmaterialfv(GL_FRONT, GL_SPECULAR, bevel_mat_specular)
    call glmaterialfv(GL_FRONT, GL_DIFFUSE, bevel_mat_diffuse)

    call glcolormaterial(GL_FRONT_AND_BACK, GL_DIFFUSE)
    call glenable(GL_COLOR_MATERIAL)
    call glshademodel(GL_SMOOTH)
  else
    call glclearindex(BACKGROUND)
    call glshademodel(GL_FLAT)
  end if

  call glmatrixmode(GL_PROJECTION)
  call gluperspective(45.0_gldouble, 1.33_gldouble, 0.1_gldouble, 100.0_gldouble)
  call glmatrixmode(GL_MODELVIEW)

return
end subroutine Init

subroutine Idle()

  integer :: i, j
  integer(kind=glenum), save :: more = GL_FALSE

  do i = 0, RINGS-1
    if (iters(i) /= 0) then
      do j = 0, 2
        offsets(i,j) = Clamp(iters(i), offsets(i,j))
      end do
      angs(i) = Clamp(iters(i), angs(i))
      iters(i) = iters(i) - 1
      more = GL_TRUE
    end if
  end do
  if (more == GL_TRUE) then
    call glutpostredisplay()
  else
    call glutidlefunc()
  end if

return
end subroutine Idle

subroutine fReshape(width,height)
integer(kind=glcint), intent(in out) ::  width, height

! if glcint is not the same as glsizei, width and height will
! need to be copied to variables of the later kind

  call glviewport(0_glint, 0_glint, width, height)

return
end subroutine fReshape

subroutine Key(ikey, x, y)
integer(kind=glcint), intent(in out) :: ikey, x, y

  select case(ikey)
  case (27) ! esc
    stop
  case (ichar(" "))
    call ReInit()
    call glutidlefunc(Idle)
  end select

return
end subroutine Key

! fortran handling of command line arguments is nonstandard, so
! this feature is omitted.  Here is the original C code.

!GLenum
!Args(int argc, char **argv)
!{
!  GLint i;
!
!  rgb = GL_TRUE;
!  doubleBuffer = GL_TRUE;
!
!  for (i = 1; i < argc; i++) {
!    if (strcmp(argv[i], "-ci") == 0) {
!      rgb = GL_FALSE;
!    } else if (strcmp(argv[i], "-rgb") == 0) {
!      rgb = GL_TRUE;
!    } else if (strcmp(argv[i], "-sb") == 0) {
!      doubleBuffer = GL_FALSE;
!    } else if (strcmp(argv[i], "-db") == 0) {
!      doubleBuffer = GL_TRUE;
!    } else {
!      printf("%s (Bad option).\n", argv[i]);
!      return GL_FALSE;
!    }
!  }
!  return GL_TRUE;
!}

subroutine visible(vis)
integer(kind=glcint), intent(in out) :: vis

  if (vis == GLUT_VISIBLE) then
    call glutidlefunc(Idle)
  else
    call glutidlefunc()
  end if

return
end subroutine visible

subroutine DrawScene()

  integer :: i

  call glpushmatrix()

  call glclear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
  call glulookat(0.0_gldouble, 0.0_gldouble, 10.0_gldouble, &
                    0.0_gldouble, 0.0_gldouble, 0.0_gldouble, &
                    0.0_gldouble, 1.0_gldouble, 0.0_gldouble)

  do i = 0, RINGS-1
    if (rgb == GL_TRUE) then
!      call glcolor3ubv(rgb_colors(i,:))
      call glcolor3fv(rgb_colors(i,:))
    else
      call glindexi(mapped_colors(i))
    end if
    call glpushmatrix()
    call gltranslatef(dests(i,0) + offsets(i,0), dests(i,1) + offsets(i,1), &
      dests(i,2) + offsets(i,2))
    call glrotatef(angs(i), rotAxis(i,0), rotAxis(i,1), rotAxis(i,2))
    call glcalllist(theTorus)
    call glpopmatrix()
  end do

  call glpopmatrix()
  if (doubleBuffer == GL_TRUE) then
    call glutswapbuffers()
  else
    call glflush()
  end if
return
end subroutine DrawScene

end module olympic_mod

program main

  use opengl_gl
  use opengl_glu
  use opengl_glut
  use olympic_mod
  implicit none

  integer(kind=glenum) :: ftype
  integer :: i

  call glutinitwindowsize(400_glcint, 300_glcint)

! not checking command line arguments
!  glutInit(&argc, argv);
!  if (Args(argc, argv) == GL_FALSE) {
!    exit(1);
!  }
  call glutinit()
  rgb = GL_TRUE           ! default values which could have been
  doubleBuffer = GL_TRUE  ! overwritten by command line arguments

  if (rgb == GL_TRUE) then
     ftype = GLUT_RGB
  else
     ftype = GLUT_INDEX
  end if

  if (doubleBuffer == GL_TRUE) then
     ftype = ior(ftype,GLUT_DOUBLE)
  else
     ftype = ior(ftype,GLUT_SINGLE)
  end if
  call glutinitdisplaymode(ftype)

  i = glutcreatewindow("Olympic")

  call Init()

  call glutreshapefunc(fReshape)
  call glutkeyboardfunc(Key)
  call glutdisplayfunc(DrawScene)

  call glutvisibilityfunc(visible)

  call glutmainloop()

  stop
end program main
