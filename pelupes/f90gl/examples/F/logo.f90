module logo_mod

use opengl_gl
use opengl_glu
use opengl_glut
implicit none
private

private :: make_letter, recalcModelView, number, animate
public :: redraw, mouse, motion, freshape, controlLights, draw_logo

type, private :: coord
   real(kind=glfloat) :: x,y
end type coord

type(coord), private, save :: angle = coord(0.0,0.0) ! in degrees
type(coord), private :: begin
logical, private :: moving
logical, private, save :: newModel = .true.
logical, private, save :: show_colorbars = .false.
real(kind=glfloat), private, dimension(4), save :: &
        lightZeroPosition = (/10.0, 4.0, 10.0, 1.0/), &
        lightZeroColor = (/1.0, 1.0, 1.0, 1.0/), &
        lightOnePosition = (/-1.0, -2.0, 1.0, 0.0/), &
        lightOneColor = (/0.5, 0.5, 0.5, 1.0/), &
        logocolor = (/1.0, 0.72, 0.0, 1.0/)

logical, private, save :: lightZeroSwitch = .true., lightOneSwitch = .false.
integer(kind=glint), private, parameter :: thelist = 20
real (kind=glfloat), private, parameter :: PI = 3.1415926
integer(kind=glint), private :: window_width, window_height

contains

subroutine make_letter(centers,radii)
type (coord), dimension(:), intent(in) :: centers
real (kind=glfloat), dimension(:), intent(in) :: radii
integer, parameter :: npoints_around_circle = 10
integer :: npoints
real (kind=glfloat) :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,alpha
type (coord) :: center1, center2
real (kind=glfloat), dimension(size(radii)) :: orient
real (kind=glfloat) :: radius1, radius2, orient1, orient2
integer :: point, i

npoints = size(radii)

do point=1,npoints
   center1 = centers(max(1,point-1))
   center2 = centers(min(npoints,point+1))
   if (center1%x==center2%x .and. center1%y==center2%y) then
      orient(point) = 0.0
   else
      orient(point) = atan2(center2%y-center1%y,center1%x-center2%x)
   end if
end do

call glbegin(GL_QUADS)

do point=1,npoints-1
   center1 = centers(point)
   center2 = centers(point+1)
   radius1 = radii(point)
   radius2 = radii(point+1)
   orient1 = orient(point)
   orient2 = orient(point+1)
   x3 = radius2*sin(orient2)+center2%x
   y3 = radius2*cos(orient2)+center2%y
   z3 = 0.0_glfloat
   x4 = radius1*sin(orient1)+center1%x
   y4 = radius1*cos(orient1)+center1%y
   z4 = 0.0_glfloat
   do i=1,npoints_around_circle
      alpha = 2.0_glfloat * PI * i / real(npoints_around_circle,glfloat)
      x1 = x4
      y1 = y4
      z1 = z4
      x2 = x3
      y2 = y3
      z2 = z3
      x3 = radius2*sin(orient2)*cos(alpha)+center2%x
      y3 = radius2*cos(orient2)*cos(alpha)+center2%y
      z3 = radius2*sin(alpha)
      x4 = radius1*sin(orient1)*cos(alpha)+center1%x
      y4 = radius1*cos(orient1)*cos(alpha)+center1%y
      z4 = radius1*sin(alpha)
      call glnormal3f(x1-center1%x,y1-center1%y,z1)
      call glvertex3f(x1,y1,z1)
      call glnormal3f(x2-center2%x,y2-center2%y,z2)
      call glvertex3f(x2,y2,z2)
      call glnormal3f(x3-center2%x,y3-center2%y,z3)
      call glvertex3f(x3,y3,z3)
      call glnormal3f(x4-center1%x,y4-center1%y,z4)
      call glvertex3f(x4,y4,z4)
   end do
end do

call glend()

return
end subroutine make_letter

subroutine recalcModelView()

  call glpopmatrix()
  call glpushmatrix()
  call gltranslatef(2.5_glfloat,0.0_glfloat,0.0_glfloat)
  call glrotatef(angle%x, 0.0, 0.0, 1.0)
  call glrotatef(angle%y, cos(PI*angle%x/180.0), -sin(PI*angle%x/180.0),0.0)
  call gltranslatef(-2.5_glfloat,0.0_glfloat,0.0_glfloat)
  newModel = .false.
  return
end subroutine recalcModelView

!          ------
subroutine number(n,x,y)
!          ------
integer, intent(in) :: n
real(kind=gldouble), intent(in) :: x,y

character(len=4) :: chn
integer :: i

write(unit=chn,fmt="(i4)") n
call glmatrixmode(GL_MODELVIEW)
call glpushmatrix()
call gltranslated(x,y,0.0_gldouble)
do i=1,4
   call glutbitmapcharacter(GLUT_BITMAP_8_BY_13,ichar(chn(i:i)))
end do
call glpopmatrix()

return
end subroutine number

subroutine redraw()

  if (newModel) then
    call recalcModelView()
  end if
  call glclear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
  call glcalllist(thelist)
  call glutswapbuffers()
  return
end subroutine redraw

subroutine mouse(button, state, x, y)
integer(kind=glcint), intent(in out) :: button, state, x, y
real :: value,xfract,yfract

  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_DOWN) then
    moving = .true.
    begin = coord(x,y)
  end if
  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_UP) then
    moving = .false.
  end if
  if (button == GLUT_MIDDLE_BUTTON .and. state == GLUT_DOWN) then
    xfract = real(x)/real(window_width)
    yfract = real(y)/real(window_height)
    value = (xfract-0.1)*2.0
    if (value >= 0.0 .and. value <= 1.0) then
       if (yfract > 0.05 .and. yfract < 0.10) then
          logocolor(1) = value
       else if (yfract > 0.15 .and. yfract < 0.20) then
          logocolor(2) = value
       else if (yfract > 0.25 .and. yfract < 0.30) then
          logocolor(3) = value
       end if
!       print *,logocolor
       call draw_logo()
       call glutpostredisplay()
    end if
  end if
          
  return
end subroutine mouse

subroutine motion(x, y)
integer(kind=glcint), intent(in out) :: x, y

  if (moving) then
    angle%x = angle%x + (x - begin%x)
    angle%y = angle%y + (y - begin%y)
    begin = coord(x,y)
    newModel = .true.
    call glutpostredisplay()
  end if
  return
end subroutine motion

subroutine freshape(w,h)
integer(kind=glcint), intent(in out) :: w,h
window_width = w
window_height = h
call glviewport(0,0,w,h)
return
end subroutine freshape

subroutine controlLights(value)
integer(kind=glcint), intent(in out) :: value

  select case (value)
  case (1)
    lightZeroSwitch = .not. lightZeroSwitch
    if (lightZeroSwitch) then
      call glenable(GL_LIGHT0)
    else
      call gldisable(GL_LIGHT0)
    end if
  case (2)
    lightOneSwitch = .not. lightOneSwitch
    if (lightOneSwitch) then
      call glenable(GL_LIGHT1)
    else
      call gldisable(GL_LIGHT1)
    end if
  case (3)
    show_colorbars = .not. show_colorbars
  case (4)
    call animate()
  end select
  call draw_logo()
  call glutpostredisplay()
  return
end subroutine controlLights

subroutine animate()
real(kind=glfloat) :: transx,transy,transz,angx,angy
integer :: nstep,i

nstep = 19
do i=0,nstep

  transx = -6.0*(nstep-i)/real(nstep) + 2.5
  transy = -4.0*(nstep-i)/real(nstep)
  transz = -160.0*(nstep-i)/real(nstep)
  angx = -10.0 -80.0*(nstep-i)/real(nstep)
  angy = -50.0 + 320.0*(nstep-i)/real(nstep)
  call glpopmatrix()
  call glpushmatrix()
  call gltranslatef(transx,transy,transz)
  call glrotatef(angx, 0.0, 0.0, 1.0)
  call glrotatef(angy, cos(PI*angx/180.0), -sin(PI*angx/180.0),0.0)
  call gltranslatef(-transx,0.0,0.0)
  call glclear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
  call glcalllist(thelist)
  call glutswapbuffers()

end do
return
end subroutine animate

subroutine draw_logo()
type (coord),allocatable,dimension(:) :: centers
real (kind=glfloat),allocatable,dimension(:) :: radii
real (kind=glfloat) :: x, alpha
integer :: int_color, i

  x = 247.0/255.0
!  x = 0.
  call glclearcolor(x,x,x,1.0)
  call gldeletelists(thelist,1_glsizei)
  call glnewlist(thelist,GL_COMPILE)

! color bars

  if (show_colorbars) then
    call gldisable(GL_LIGHTING)
    call glmatrixmode(GL_MODELVIEW)
    call glpushmatrix()
    call glloadidentity()
    call glmatrixmode(GL_PROJECTION)
    call glpushmatrix()
    call glloadidentity()
    call glortho(0.0_gldouble,1.0_gldouble,0.0_gldouble,1.0_gldouble,-1.0_gldouble,1.0_gldouble)
    call glcolor3f(1.0_glfloat,0.0_glfloat,0.0_glfloat)
    call glrectf(0.1_glfloat,0.9_glfloat,0.1_glfloat+0.5_glfloat*logocolor(1),0.95_glfloat)
    call glcolor3f(0.0_glfloat,1.0_glfloat,0.0_glfloat)
    call glrectf(0.1_glfloat,0.8_glfloat,0.1_glfloat+0.5_glfloat*logocolor(2),0.85_glfloat)
    call glcolor3f(0.0_glfloat,0.0_glfloat,1.0_glfloat)
    call glrectf(0.1_glfloat,0.7_glfloat,0.1_glfloat+0.5_glfloat*logocolor(3),0.75_glfloat)
    call glcolor3f(1.0_glfloat,1.0_glfloat,1.0_glfloat)
    call glrectf(0.1_glfloat,0.9_glfloat,0.6_glfloat,0.95_glfloat)
    call glrectf(0.1_glfloat,0.8_glfloat,0.6_glfloat,0.85_glfloat)
    call glrectf(0.1_glfloat,0.7_glfloat,0.6_glfloat,0.75_glfloat)
!    int_color = int(logocolor(1)*256)
!    call number(int_color,.7_gldouble,.9_gldouble)
!    int_color = int(logocolor(2)*256)
!    call number(int_color,.7_gldouble,.8_gldouble)
!    int_color = int(logocolor(3)*256)
!    call number(int_color,.7_gldouble,.7_gldouble)
    call glmatrixmode(GL_PROJECTION)
    call glpopmatrix()
    call glmatrixmode(GL_MODELVIEW)
    call glpopmatrix()
  end if

  call gllightmodelfv(GL_LIGHT_MODEL_AMBIENT, (/0.25,0.25,0.25,1.0/))
  call gllightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition)
  call gllightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor)
  call gllightfv(GL_LIGHT0, GL_SPECULAR, lightZeroColor)
  call gllightfv(GL_LIGHT1, GL_POSITION, lightOnePosition)
  call gllightfv(GL_LIGHT1, GL_DIFFUSE, lightOneColor)
  if (lightZeroSwitch) then
    call glenable(GL_LIGHT0)
  else
    call gldisable(GL_LIGHT0)
  end if
  if (lightOneSwitch) then
    call glenable(GL_LIGHT1)
  else
    call gldisable(GL_LIGHT1)
  end if
  call glenable(GL_LIGHTING)
  call glpolygonmode(GL_FRONT_AND_BACK, GL_FILL)
! curve of f
  allocate(centers(41),radii(41))
  call glpushmatrix()
  call glrotatef(25.0_glfloat,0.0_glfloat,0.0_glfloat,1.0_glfloat)
  do i=1,41
     x = (i-21)/20.0_glfloat
     centers(i)%x = x
     centers(i)%y = -x*(x-1.3)*(x+1.3)
     if (x<=-0.5) then
        radii(i) = 0.2*x+0.2
     else if (x>=0.5) then
        radii(i) = 0.2-0.2*x
     else
        radii(i) = 0.1_glfloat
     end if
  end do
  call glmaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, logocolor)
  call glmaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, logocolor)
  call glmaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128.0_glfloat)
  call make_letter(centers,radii)
  call glpopmatrix()
  deallocate(centers,radii)
! crossbar of f
  allocate(centers(4),radii(4))
  centers(1) = coord(-0.4,0.0)
  centers(2) = coord(-0.4,0.0)
  centers(3) = coord(0.4,0.0)
  centers(4) = coord(0.4,0.0)
  radii(1) = 0.0
  radii(2) = 0.1
  radii(3) = 0.1
  radii(4) = 0.0
  call make_letter(centers,radii)
  deallocate(centers,radii)
! left part of 9
  allocate(centers(28),radii(28))
  call glpushmatrix()
  call gltranslatef(1.0_glfloat,0.0_glfloat,0.0_glfloat)
  call glrotatef(10.0_glfloat,0.0_glfloat,0.0_glfloat,-1.0_glfloat)
  do i=1,28
     alpha = (i-6)*PI/20.0 + PI/2.0 + PI/10.0
     centers(i)%x = 0.4*cos(alpha)
     centers(i)%y = 0.4*sin(alpha) + 0.4
     radii(i) = 0.1
     if (i>20) then
        radii(i) = (28-i)/80.0
     end if
  end do
  call make_letter(centers,radii)
  call glpopmatrix()
  deallocate(centers,radii)
! right part of 9
  allocate(centers(25),radii(25))
  call glpushmatrix()
  call gltranslatef(1.0_glfloat,0.0_glfloat,0.0_glfloat)
  call glrotatef(10.0_glfloat,0.0_glfloat,0.0_glfloat,-1.0_glfloat)
  do i=1,25
     alpha = (i-1)*PI/20.0 - PI/2.0 - 4.0*PI/20.0
     centers(i)%x = 0.5*cos(alpha)
     centers(i)%y = 0.8*sin(alpha)
     radii(i) = 0.1
     if (i<5) then
        radii(i) = i/50.0
     end if
  end do
  call make_letter(centers,radii)
  call glpopmatrix()
  deallocate(centers,radii)
! 0
  allocate(centers(42),radii(42))
  call glpushmatrix()
  call gltranslatef(2.4_glfloat,0.0_glfloat,0.0_glfloat)
  do i=1,42
     alpha = (i-1)*PI/20.0
     centers(i)%x = 0.6*cos(alpha)
     centers(i)%y = 0.8*sin(alpha)
     radii(i) = 0.1
  end do
  call make_letter(centers,radii)
  call glpopmatrix()
  deallocate(centers,radii)
! g
  allocate(centers(15+56),radii(15+56))
! right part of g
  call glpushmatrix()
  call gltranslatef(3.6_glfloat,-0.5_glfloat,0.0_glfloat)
!  call glrotatef(10.0_glfloat,0.0_glfloat,0.0_glfloat,-1.0_glfloat)
  do i=1,15
     alpha = (i-1)*PI/20.0 - PI/2.0 - 4.0*PI/20.0
     centers(i)%x = 0.5*cos(alpha) - 0.1
     centers(i)%y =1.2*sin(alpha)
     radii(i) = 0.1
     if (i<5) then
         radii(i) = i/50.0
     end if
  end do
! left part of g
!  call glrotatef(20.0_glfloat,0.0_glfloat,0.0_glfloat,-1.0_glfloat)
  do i=1+15,55+15
     alpha = 3.0*(i-1-15)*PI/80.0
     centers(i)%x = 0.4*cos(alpha)
     centers(i)%y = 0.4*sin(alpha)
     radii(i) = 0.1
  end do
  centers(56+15) = centers(55+15)
  radii(56+15) = 0.0
  call make_letter(centers,radii)
  call glpopmatrix()
  deallocate(centers,radii)
! l
  allocate(centers(41),radii(41))
  call glpushmatrix()
  call gltranslatef(4.6_glfloat,0.0_glfloat,0.0_glfloat)
  call glrotatef(25.0_glfloat,0.0_glfloat,0.0_glfloat,1.0_glfloat)
  do i=1,41
     x = (i-21)/20.0_glfloat
     centers(i)%x = x
     centers(i)%y = -x*(x-1.4)*(x+1.4)
     if (x<=-0.5) then
        radii(i) = 0.2*x+0.2
     else if (x>=0.5) then
        radii(i) = 0.2-0.2*x
     else
        radii(i) = 0.1_glfloat
     end if
  end do
  call make_letter(centers,radii)
  call glpopmatrix()
  deallocate(centers,radii)

  call glendlist()

  return
end subroutine draw_logo

end module logo_mod

program main

use opengl_gl
use opengl_glu
use opengl_glut
use logo_mod
implicit none

integer :: i, window_width, window_height

  call glutinit()
  call glutinitdisplaymode(ior(ior(GLUT_RGB,GLUT_DOUBLE),GLUT_DEPTH))
  i = glutcreatewindow("f90gl logo")
  call glutreshapewindow(320,320)
  window_width = 320
  window_height = 320
  call glutdisplayfunc(redraw)
  call glutmousefunc(mouse)
  call glutmotionfunc(motion)
  call glutreshapefunc(freshape)
  call glutcreatemenu(controlLights,i)
  call glutaddmenuentry("Toggle right light", 1)
  call glutaddmenuentry("Toggle left light", 2)
  call glutaddmenuentry("Toggle color bars", 3)
  call glutaddmenuentry("Animate",4)
  call glutattachmenu(GLUT_RIGHT_BUTTON)
  call glenable(GL_CULL_FACE)
  call glenable(GL_DEPTH_TEST)
  call glenable(GL_LIGHTING)
  call glenable(GL_NORMALIZE)
  call glmatrixmode(GL_PROJECTION)
  call gluperspective(10.0_gldouble, & ! field of view in degree
                          1.0_gldouble, & ! aspect ratio
                          1.0_gldouble, & ! Z near
                        400.0_gldouble)   ! Z far
  call glmatrixmode(GL_MODELVIEW)
  call glulookat( &
   2.0_gldouble, 0.0_gldouble, 40.0_gldouble, & ! eye is at (0,0,30)
   2.0_gldouble, 0.0_gldouble, 0.0_gldouble, & ! center is at (0,0,0)
   0.0_gldouble, 1.0_gldouble, 0.0_gldouble) ! up is in postivie Y direction
  call glpushmatrix()  ! dummy push so we can pop on model recalc 
  call draw_logo()

  call glutmainloop()
  stop
end program main
