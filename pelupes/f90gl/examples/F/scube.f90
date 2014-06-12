module scube_mod
use opengl_gl
use opengl_glut
implicit none
private
public :: keyboard, display, visible, fog_select, menu_select, buildColormap
private :: setColor, drawCube, drawCheck, myShadowMatrix, idle


logical, save, public :: &
   useRGB      = .true., &
   useLighting = .true., &
   useFog      = .false., &
   useDB       = .true., &
   useLogo     = .true., &
   useQuads    = .true.

integer, save, private :: tick = -1
logical, save, private :: moving = .true.

integer, parameter, private :: &
   GREY    = 0, &
   RED     = 1, &
   GREEN   = 2, &
   BLUE    = 3, &
   CYAN    = 4, &
   MAGENTA = 5, &
   YELLOW  = 6, &
   BLACK   = 7

real(kind=glfloat), save, private, dimension(8,4) :: materialColor = reshape( &
  (/ 0.8, 0.8, 0.0, 0.0, 0.0, 0.8, 0.8, 0.0, &
     0.8, 0.0, 0.8, 0.0, 0.8, 0.0, 0.8, 0.0, &
     0.8, 0.0, 0.0, 0.8, 0.8, 0.8, 0.0, 0.0, &
     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.6 /), &
     (/8,4/))

real(kind=glfloat), save, public, dimension(4) :: &
 lightPos = (/2.0, 4.0, 2.0, 1.0/), &
 lightDir = (/-2.0, -4.0, -2.0, 1.0/), &
 lightAmb = (/0.2, 0.2, 0.2, 1.0/), &
 lightDiff = (/0.8, 0.8, 0.8, 1.0/), &
 lightSpec = (/0.4, 0.4, 0.4, 1.0/)

real(kind=glfloat), save, private, dimension(4) :: &
 groundPlane = (/0.0, 1.0, 0.0, 1.499/), &
 backPlane = (/0.0, 0.0, 1.0, 0.899/)

real(kind=glfloat), save, public, dimension(4) :: &
 fogColor = (/0.0, 0.0, 0.0, 0.0/)
real(kind=glfloat), save, public, dimension(1) :: &
 fogIndex = (/0.0/)

integer(kind=glubyte), private, parameter :: v1 = -42, v2 = 85
integer(kind=glubyte), save, public, dimension(128) :: shadowPattern = (/ &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2, &
  v1, v1, v1, v1, v2, v2, v2, v2    /)

character(len=30), save, public :: windowNameRGBDB = "shadow cube (OpenGL RGB DB)"
character(len=30), save, public :: windowNameRGB = "shadow cube (OpenGL RGB)"
character(len=30), save, public :: windowNameIndexDB = "shadow cube (OpenGL Index DB)"
character(len=30), save, public :: windowNameIndex = "shadow cube (OpenGL Index)"

contains

subroutine buildColormap()

integer :: mapSize,rampSize,entry,i,hue
real(kind=glfloat) :: val,r,g,b

mapSize = 2**glutget(GLUT_WINDOW_BUFFER_SIZE)
rampSize = mapSize / 8

  if (useRGB) then
    return
  else

    do entry=0,mapSize-1
      hue = entry / rampSize
      val = modulo(entry,rampSize) * (1.0 / (rampSize - 1))

      if (hue==0 .or. hue==1 .or. hue==5 .or. hue==6) then
         r = val
      else
         r = 0
      end if
      if (hue==0 .or. hue==2 .or. hue==4 .or. hue==6) then
         g = val
      else
         g = 0
      end if
      if (hue==0 .or. hue==3 .or. hue==4 .or. hue==5) then
         b = val
      else
         b = 0
      end if

      call glutsetcolor(entry, r, g, b)
    end do

    do i=1,8
      materialColor(i,1) = i * rampSize + 0.2 * (rampSize - 1)
      materialColor(i,2) = i * rampSize + 0.8 * (rampSize - 1)
      materialColor(i,3) = i * rampSize + 1.0 * (rampSize - 1)
      materialColor(i,4) = 0.0
    end do

    fogIndex(1) = -0.2 * (rampSize - 1)
  end if
return
end subroutine buildColormap

subroutine setColor(c)
integer, intent(in) :: c

  if (useLighting) then
    if (useRGB) then
      call glmaterialfv(GL_FRONT_AND_BACK, &
        GL_AMBIENT_AND_DIFFUSE, materialColor(c+1,:))
    else
      call glmaterialfv(GL_FRONT_AND_BACK, &
        GL_COLOR_INDEXES, materialColor(c+1,:))
    end if
  else
    if (useRGB) then
      call glcolor4fv(materialColor(c+1,:))
    else
      call glindexf(materialColor(c+1,1))
    end if
  end if
  return
end subroutine setColor

subroutine drawCube(color)
integer, intent(in) :: color

real(kind=glfloat), save, dimension(4,4,6) :: cube_vertexes = reshape( (/ &
    -1.0, -1.0, -1.0, 1.0, &
    -1.0, -1.0, 1.0, 1.0, &
    -1.0, 1.0, 1.0, 1.0, &
    -1.0, 1.0, -1.0, 1.0, &

    1.0, 1.0, 1.0, 1.0, &
    1.0, -1.0, 1.0, 1.0, &
    1.0, -1.0, -1.0, 1.0, &
    1.0, 1.0, -1.0, 1.0, &

    -1.0, -1.0, -1.0, 1.0, &
    1.0, -1.0, -1.0, 1.0, &
    1.0, -1.0, 1.0, 1.0, &
    -1.0, -1.0, 1.0, 1.0, &

    1.0, 1.0, 1.0, 1.0, &
    1.0, 1.0, -1.0, 1.0, &
    -1.0, 1.0, -1.0, 1.0, &
    -1.0, 1.0, 1.0, 1.0, &

    -1.0, -1.0, -1.0, 1.0, &
    -1.0, 1.0, -1.0, 1.0, &
    1.0, 1.0, -1.0, 1.0, &
    1.0, -1.0, -1.0, 1.0, &

    1.0, 1.0, 1.0, 1.0, &
    -1.0, 1.0, 1.0, 1.0, &
    -1.0, -1.0, 1.0, 1.0, &
    1.0, -1.0, 1.0, 1.0 /), &
    (/4,4,6/) )

real(kind=glfloat), save, dimension(4,6) :: cube_normals = reshape( (/ &
  -1.0, 0.0, 0.0, 0.0, &
  1.0, 0.0, 0.0, 0.0, &
  0.0, -1.0, 0.0, 0.0, &
  0.0, 1.0, 0.0, 0.0, &
  0.0, 0.0, -1.0, 0.0, &
  0.0, 0.0, 1.0, 0.0 /), &
  (/4,6/) )

  integer :: i

  call setColor(color)

  do i=1,6
    call glnormal3fv(cube_normals(:,i))
    call glbegin(GL_POLYGON)
    call glvertex4fv(cube_vertexes(:,1,i))
    call glvertex4fv(cube_vertexes(:,2,i))
    call glvertex4fv(cube_vertexes(:,3,i))
    call glvertex4fv(cube_vertexes(:,4,i))
    call glend()
  end do
  return
end subroutine drawCube

subroutine drawCheck(w,h,evenColor,oddColor)
integer, intent(in) :: w,h,evenColor,oddColor

  logical, save :: initialized = .false., &
                   usedLighting = .false.
  integer(kind=gluint), save :: checklist = 0
  real, save, dimension(4) :: square_normal = (/0.0, 0.0, 1.0, 0.0/)
  real, save, dimension(4,4) :: square
  integer :: i,j

  if (.not. initialized .or. (usedLighting .eqv. useLighting)) then

    if (checklist == 0) then
      checklist = glgenlists(1)
    end if
    call glnewlist(checklist, GL_COMPILE_AND_EXECUTE)

    if (useQuads) then
      call glnormal3fv(square_normal)
      call glbegin(GL_QUADS)
    end if
    do j=0,h-1
      do i=0,w-1
        square(:,1) = (/ -1.0 + 2.0/w * i, -1.0 + 2.0/h * (j+1), 0.0, 1.0/)
        square(:,2) = (/ -1.0 + 2.0/w * i, -1.0 + 2.0/h * j, 0.0, 1.0/)
        square(:,3) = (/ -1.0 + 2.0/w * (i+1), -1.0 + 2.0/h * j, 0.0, 1.0/)
        square(:,4) = (/ -1.0 + 2.0/w * (i+1), -1.0 + 2.0/h * (j+1), 0.0, 1.0/)

        if (ieor(iand(i,1),iand(j,1)) /= 0) then
          call setColor(oddColor)
        else
          call setColor(evenColor)
        end if

        if (.not.useQuads) then
          call glbegin(GL_POLYGON)
        end if
        call glvertex4fv(square(:,1))
        call glvertex4fv(square(:,2))
        call glvertex4fv(square(:,3))
        call glvertex4fv(square(:,4))
        if (.not.useQuads) then
          call glend()
        end if
      end do
    end do

    if (useQuads) then
      call glend()
    end if
    call glendlist()

    initialized = .true.
    usedLighting = useLighting
  else
    call glcalllist(checklist)
  end if
return
end subroutine drawCheck

subroutine myShadowMatrix(ground,light)
real, dimension(:), intent(in) :: ground, light

  real :: dot
  real(kind=glfloat), dimension(4,4) :: shadowMat
  integer :: i

  dot = dot_product(ground,light)

  do i=1,4
     shadowMat(i,:) = -light(i)*ground
     shadowMat(i,i) = shadowMat(i,i) + dot
  end do

  call glmultmatrixf(shadowMat)
return
end subroutine myShadowMatrix

subroutine idle()

  tick = tick + 1
  if (tick >= 120) then
    tick = 0
  end if
  call glutpostredisplay()
return
end subroutine idle

subroutine keyboard(ich, x, y)
integer, intent(in out) :: ich,x,y

  character(len=1) :: ch
  real(kind=glfloat) :: rGL_LINEAR, rGL_EXP, rGL_EXP2

  ch = char(ich)
  select case(ch)
  case (char(27))      ! escape
    stop
  case ("L","l")
    useLighting = .not. useLighting
    if (useLighting) then
       call glenable(GL_LIGHTING) 
    else
       call gldisable(GL_LIGHTING)
    end if
    call glutpostredisplay()
  case ("F","f")
    useFog = .not. useFog
    if (useFog) then
       call glenable(GL_FOG)
    else
       call gldisable(GL_FOG)
    end if
    call glutpostredisplay()
  case ("1")
    rGL_LINEAR = GL_LINEAR
    call glfogf(GL_FOG_MODE, rGL_LINEAR)
    call glutpostredisplay()
  case ("2")
    rGL_EXP = GL_EXP
    call glfogf(GL_FOG_MODE, rGL_EXP)
    call glutpostredisplay()
  case ("3")
    rGL_EXP2 = GL_EXP2
    call glfogf(GL_FOG_MODE, rGL_EXP2)
    call glutpostredisplay()
  case (" ")
    if (.not. moving) then
      call idle()
      call glutpostredisplay()
    end if
  end select
return
end subroutine keyboard

subroutine display()

  real(kind=glfloat), dimension(16) :: cubeXform

  call glclear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))

  call glpushmatrix()
  call gltranslatef(0.0, -1.5, 0.0) ! taking a chance that glfloat is
  call glrotatef(-90.0, 1.0, 0.0, 0.0)    ! the same as the default real
  call glscalef(2.0, 2.0, 2.0)

  call drawCheck(6, 6, BLUE, YELLOW)  ! draw ground
  call glpopmatrix()

  call glpushmatrix()
  call gltranslatef(0.0, 0.0, -0.9)
  call glscalef(2.0, 2.0, 2.0)

  call drawCheck(6, 6, BLUE, YELLOW)  ! draw back
  call glpopmatrix()

  call glpushmatrix()
  call gltranslatef(0.0, 0.2, 0.0)
  call glscalef(0.3, 0.3, 0.3)
  call glrotatef((360.0 / (30 * 1)) * tick, 1.0, 0.0, 0.0)
  call glrotatef((360.0 / (30 * 2)) * tick, 0.0, 1.0, 0.0)
  call glrotatef((360.0 / (30 * 4)) * tick, 0.0, 0.0, 1.0)
  call glscalef(1.0, 2.0, 1.0)
  call glgetfloatv(GL_MODELVIEW_MATRIX, cubeXform)

  call drawCube(RED)        ! draw cube
  call glpopmatrix()

  call gldepthmask(.false._glboolean)
  if (useRGB) then
    call glenable(GL_BLEND)
  else
    call glenable(GL_POLYGON_STIPPLE)
  end if
  if (useFog) then
    call gldisable(GL_FOG)
  end if
  call glpushmatrix()
  call myShadowMatrix(groundPlane, lightPos)
  call gltranslatef(0.0, 0.0, 2.0)
  call glmultmatrixf(reshape(cubeXform,(/4,4/)))

  call drawCube(BLACK)      ! draw ground shadow
  call glpopmatrix()

  call glpushmatrix()
  call myShadowMatrix(backPlane, lightPos)
  call gltranslatef(0.0, 0.0, 2.0)
  call glmultmatrixf(reshape(cubeXform,(/4,4/)))

  call drawCube(BLACK)      ! draw back shadow
  call glpopmatrix()

  call gldepthmask(.true._glboolean)
  if (useRGB) then
    call gldisable(GL_BLEND)
  else
    call gldisable(GL_POLYGON_STIPPLE)
  end if
  if (useFog) then
    call glenable(GL_FOG)
  end if
  if(useDB) then
    call glutswapbuffers()
  else
    call glflush()
  end if
return
end subroutine display

subroutine fog_select(fog)
integer, intent(in out) :: fog

  real(kind=glfloat) :: rfog
  rfog = fog
  call glfogf(GL_FOG_MODE, rfog)
  call glutpostredisplay()
return
end subroutine fog_select

subroutine menu_select(mode)
integer, intent(in out) :: mode

  select case (mode)
  case (1)
    moving = .true.
    call glutidlefunc(idle)
  case (2)
    moving = .false.
    call glutidlefunc()
  case (3)
    useFog = .not. useFog
    if (useFog) then
       call glenable(GL_FOG)
    else
       call gldisable(GL_FOG)
    end if
    call glutpostredisplay()
  case (4)
    useLighting = .not. useLighting
    if (useLighting) then
       call glenable(GL_LIGHTING)
    else
       call gldisable(GL_LIGHTING)
    end if
    call glutpostredisplay()
  case (5)
    stop
  end select
return
end subroutine menu_select

subroutine visible(state)
integer, intent(in out) :: state

  if (state == GLUT_VISIBLE) then
    if (moving) then
      call glutidlefunc(idle)
    end if
  else
    if (moving) then
      call glutidlefunc()
    end if
  end if
return
end subroutine visible

end module scube_mod

program main
use opengl_gl
use opengl_glut
use scube_mod
implicit none

  integer, parameter :: width = 350, height = 350
  integer :: i, win
  character(len=30) :: name
  integer :: fog_menu
  real(kind=glfloat) :: rGL_EXP

  call glutinitwindowsize(width, height)
  call glutinit()

  ! choose visual
  if (useRGB) then
    if (useDB) then
      call glutinitdisplaymode(ior(ior(GLUT_DOUBLE,GLUT_RGB),GLUT_DEPTH))
      name = windowNameRGBDB
    else
      call glutinitdisplaymode(ior(ior(GLUT_SINGLE,GLUT_RGB),GLUT_DEPTH))
      name = windowNameRGB
    end if
  else
    if (useDB) then
      call glutinitdisplaymode(ior(ior(GLUT_DOUBLE,GLUT_INDEX),GLUT_DEPTH))
      name = windowNameIndexDB
    else
      call glutinitdisplaymode(ior(ior(GLUT_SINGLE,GLUT_INDEX),GLUT_DEPTH))
      name = windowNameIndex
    end if
  end if

  win = glutcreatewindow(name)

  call buildColormap()

  call glutkeyboardfunc(keyboard)
  call glutdisplayfunc(display)
  call glutvisibilityfunc(visible)

  call glutcreatemenu(fog_select,fog_menu)
  call glutaddmenuentry("Linear fog", GL_LINEAR)
  call glutaddmenuentry("Exp fog", GL_EXP)
  call glutaddmenuentry("Exp^2 fog", GL_EXP2)

  call glutcreatemenu(menu_select,i)
  call glutaddmenuentry("Start motion", 1)
  call glutaddmenuentry("Stop motion", 2)
  call glutaddmenuentry("Toggle fog", 3)
  call glutaddmenuentry("Toggle lighting", 4)
  call glutaddsubmenu("Fog type", fog_menu)
  call glutaddmenuentry("Quit", 5)
  call glutattachmenu(GLUT_RIGHT_BUTTON)

  ! setup context
  call glmatrixmode(GL_PROJECTION)
  call glloadidentity()
  call glfrustum(-1.0_gldouble, 1.0_gldouble, -1.0_gldouble, &
                     1.0_gldouble, 1.0_gldouble, 3.0_gldouble)

  call glmatrixmode(GL_MODELVIEW)
  call glloadidentity()
  call gltranslatef(0.0, 0.0, -2.0)

  call glenable(GL_DEPTH_TEST)

  if (useLighting) then
    call glenable(GL_LIGHTING)
  end if
  call glenable(GL_LIGHT0)
  call gllightfv(GL_LIGHT0, GL_POSITION, lightPos)
  call gllightfv(GL_LIGHT0, GL_AMBIENT, lightAmb)
  call gllightfv(GL_LIGHT0, GL_DIFFUSE, lightDiff)
  call gllightfv(GL_LIGHT0, GL_SPECULAR, lightSpec)
  
   ! call gllightfv(GL_LIGHT0, GL_SPOT_DIRECTION, lightDir);
   ! call gllightf(GL_LIGHT0, GL_SPOT_EXPONENT, 80);
   ! call gllightf(GL_LIGHT0, GL_SPOT_CUTOFF, 25);

  call glenable(GL_NORMALIZE)

  if (useFog) then
    call glenable(GL_FOG)
  end if
  call glfogfv(GL_FOG_COLOR, fogColor)
  call glfogfv(GL_FOG_INDEX, fogIndex)
  rGL_EXP = GL_EXP
  call glfogf(GL_FOG_MODE, rGL_EXP)
  call glfogf(GL_FOG_DENSITY, 0.5)
  call glfogf(GL_FOG_START, 1.0)
  call glfogf(GL_FOG_END, 3.0)

  call glenable(GL_CULL_FACE)
  call glcullface(GL_BACK)

  call glshademodel(GL_SMOOTH)

  call glblendfunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
  call glpolygonstipple(shadowPattern)

  call glclearcolor(0.0, 0.0, 0.0, 1.0)
  call glclearindex(0.0)
  call glcleardepth(1.0_gldouble)

  call glutmainloop()
  stop
end program main
