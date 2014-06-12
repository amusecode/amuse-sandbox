module dino_mod

use opengl_gl
use opengl_glu
use opengl_glut

implicit none
private
public :: redraw, mouse, motion, controlLights, makeDinosaur
private :: glbegin_cb, glvertex2fv_cb, glend_cb, extrudeSolidFrompolygon, &
           recalcmodelView

integer, private, parameter :: &
  RESERVED   =  0, &
  BODY_SIDE  =  1, &
  BODY_EDGE  =  2, &
  BODY_WHOLE =  3, &
  ARM_SIDE   =  4, &
  ARM_EDGE   =  5, &
  ARM_WHOLE  =  6, &
  LEG_SIDE   =  7, &
  LEG_EDGE   =  8, &
  LEG_WHOLE  =  9, &
  EYE_SIDE   = 10, &
  EYE_EDGE   = 11, &
  EYE_WHOLE  = 12, &
  DINOSAUR   = 13

type, private :: coord
   real(kind=glfloat) :: x,y
end type coord

type, private :: tessvertex
   real(kind=glfloat), dimension(2) :: val
end type tessvertex

type(coord), private :: angle = coord(-150.0,0.0) ! in degrees
type(coord), private :: begin
logical, private :: moving
integer, private :: W = 300, H = 300
real(kind=glfloat), private :: bodyWidth = 2.0
logical, private :: newModel = .true.

real(kind=glfloat), dimension(4), private :: skinColor = (/0.1, 1.0, 0.1, 1.0/), &
        eyeColor = (/1.0, 0.2, 0.2, 1.0/)
logical, private :: lightZeroSwitch = .true., lightOneSwitch = .true.

contains

subroutine glBegin_cb(type)
integer(kind=glenum), intent(in out) :: type
call glBegin(type)
end subroutine glBegin_cb

subroutine glVertex2fv_cb(vertex_data)
real(kind=glfloat), dimension(:), intent(in out) :: vertex_data
call glVertex2fv(vertex_data)
end subroutine glVertex2fv_cb

subroutine glEnd_cb()
call glEnd()
end subroutine glEnd_cb

subroutine extrudeSolidFromPolygon(data, dataSize, &
                                   thickness, side, edge, whole)
type(tessvertex), target, intent(in out), dimension(:) :: data
integer, intent(in) :: dataSize
real(kind=glfloat), intent(in) :: thickness
integer(kind=gluint), intent(in) :: side, edge, whole

  type(GLUtesselator), save, pointer :: tobj
  logical :: tobj_null = .true.

  
  real(kind=gldouble), dimension(3) :: vertex
  real(kind=glfloat) :: dx, dy, len
  integer :: i, count

  count = dataSize

  if (tobj_null) then
    tobj => gluNewTess()  ! create and initialize a GLU
                             ! polygon * * tesselation object
    call gluTessCallback(tobj, GLU_TESS_BEGIN, glBegin_cb)
    call gluTessCallback(tobj, GLU_TESS_VERTEX, glVertex2fv_cb)
    call gluTessCallback(tobj, GLU_TESS_END, glEnd_cb)
  endif
  call glNewList(side, GL_COMPILE)
  call glShadeModel(GL_SMOOTH)  ! smooth minimizes seeing tessellation
  call gluBeginPolygon(tobj)
  do i=1,count
     vertex(1) = data(i)%val(1)
     vertex(2) = data(i)%val(2)
     vertex(3) = 0.0
     call gluTessVertex(tobj, vertex, data(i)%val)
  end do
  call gluEndPolygon(tobj)
  call glEndList()
  call glNewList(edge, GL_COMPILE)
  call glShadeModel(GL_FLAT)  ! flat shade keeps angular hands
                            ! from being * * "smoothed"
  call glBegin(GL_QUAD_STRIP)
  do i=0,count
!      mod function handles closing the edge
    call glVertex3f(data(modulo(i,count)+1)%val(1), data(modulo(i,count)+1)%val(2), 0.0)
    call glVertex3f(data(modulo(i,count)+1)%val(1), data(modulo(i,count)+1)%val(2), thickness)
!       Calculate a unit normal by dividing by Euclidean
!       distance. We * could be lazy and use
!       glEnable(GL_NORMALIZE) so we could pass in * arbitrary
!       normals for a very slight performance hit.
    dx = data(modulo((i + 1),count)+1)%val(2) - data(modulo(i,count)+1)%val(2)
    dy = data(modulo(i,count)+1)%val(1) - data(modulo((i + 1), count)+1)%val(1)
    len = sqrt(dx * dx + dy * dy)
    call glNormal3f(dx / len, dy / len, 0.0)
  end do
  call glEnd()
  call glEndList()
  call glNewList(whole, GL_COMPILE)
  call glFrontFace(GL_CW)
  call glCallList(edge)
  call glNormal3f(0.0, 0.0, -1.0)  ! constant normal for side
  call glCallList(side)
  call glPushMatrix()
  call glTranslatef(0.0, 0.0, thickness)
  call glFrontFace(GL_CCW)
  call glNormal3f(0.0, 0.0, 1.0)  ! opposite normal for other side
  call glCallList(side)
  call glPopMatrix()
  call glEndList()
return
end subroutine extrudeSolidFromPolygon

subroutine makeDinosaur()

  real(kind=glfloat) :: bodyWidth = 3.0

type(tessvertex), target, save, dimension(22) :: body = &
  (/ tessvertex((/0.0, 3.0/)), tessvertex((/1.0, 1.0/)), tessvertex((/5.0, 1.0/)), &
     tessvertex((/8.0, 4.0/)), tessvertex((/10.0, 4.0/)), tessvertex((/11.0, 5.0/)), &
tessvertex((/11.0, 11.5/)), tessvertex((/13.0, 12.0/)), tessvertex((/13.0, 13.0/)), &
 tessvertex((/10.0, 13.5/)), tessvertex((/13.0, 14.0/)), tessvertex((/13.0, 15.0/)),&
 tessvertex((/11.0, 16.0/)), tessvertex((/8.0, 16.0/)), tessvertex((/7.0, 15.0/)), &
 tessvertex((/7.0, 13.0/)), tessvertex((/8.0, 12.0/)), tessvertex((/7.0, 11.0/)), &
 tessvertex((/6.0, 6.0/)), tessvertex((/4.0, 3.0/)), tessvertex((/3.0, 2.0/)), &
 tessvertex((/1.0, 2.0/)) /)

type(tessvertex), target, save, dimension(16) :: arm = &
  (/ tessvertex((/8.0, 10.0/)), tessvertex((/9.0, 9.0/)), tessvertex((/10.0, 9.0/)), &
    tessvertex((/13.0, 8.0/)), tessvertex((/14.0, 9.0/)), tessvertex((/16.0, 9.0/)), &
  tessvertex((/15.0, 9.5/)), tessvertex((/16.0, 10.0/)), tessvertex((/15.0, 10.0/)),&
tessvertex((/15.5, 11.0/)), tessvertex((/14.5, 10.0/)), tessvertex((/14.0, 11.0/)),&
  tessvertex((/14.0, 10.0/)), tessvertex((/13.0, 9.0/)), tessvertex((/11.0, 11.0/)), &
  tessvertex((/9.0, 11.0/)) /)

type(tessvertex), target, save, dimension(14) :: leg = &
  (/ tessvertex((/8.0, 6.0/)), tessvertex((/8.0, 4.0/)), tessvertex((/9.0, 3.0/)), &
     tessvertex((/9.0, 2.0/)), tessvertex((/8.0, 1.0/)), tessvertex((/8.0, 0.5/)), &
     tessvertex((/9.0, 0.0/)), tessvertex((/12.0, 0.0/)), tessvertex((/10.0, 1.0/)), &
     tessvertex((/10.0, 2.0/)), tessvertex((/12.0, 4.0/)), tessvertex((/11.0, 6.0/)),&
     tessvertex((/10.0, 7.0/)), tessvertex((/9.0, 7.0/)) /)

type(tessvertex), target, save, dimension(6) :: eye = &
  (/ tessvertex((/8.75, 15.0/)), tessvertex((/9.0, 14.7/)), &
     tessvertex((/9.6, 14.7/)), tessvertex((/10.1, 15.0/)), &
     tessvertex((/9.6, 15.25/)), tessvertex((/9.0, 15.25/)) /)

  call extrudeSolidFromPolygon(body, size(body), bodyWidth, &
    BODY_SIDE, BODY_EDGE, BODY_WHOLE)
  call extrudeSolidFromPolygon(arm, size(arm), bodyWidth / 4, &
    ARM_SIDE, ARM_EDGE, ARM_WHOLE)
  call extrudeSolidFromPolygon(leg, size(leg), bodyWidth / 2, &
    LEG_SIDE, LEG_EDGE, LEG_WHOLE)
  call extrudeSolidFromPolygon(eye, size(eye), bodyWidth + 0.2, &
    EYE_SIDE, EYE_EDGE, EYE_WHOLE)
  call glNewList(DINOSAUR, GL_COMPILE)
  call glMaterialfv(GL_FRONT, GL_DIFFUSE, skinColor)
  call glCallList(BODY_WHOLE)
  call glPushMatrix()
  call glTranslatef(0.0, 0.0, bodyWidth)
  call glCallList(ARM_WHOLE)
  call glCallList(LEG_WHOLE)
  call glTranslatef(0.0, 0.0, -bodyWidth - bodyWidth / 4)
  call glCallList(ARM_WHOLE)
  call glTranslatef(0.0, 0.0, -bodyWidth / 4)
  call glCallList(LEG_WHOLE)
  call glTranslatef(0.0, 0.0, bodyWidth / 2 - 0.1)
  call glMaterialfv(GL_FRONT, GL_DIFFUSE, eyeColor)
  call glCallList(EYE_WHOLE)
  call glPopMatrix()
  call glEndList()
end subroutine makedinosaur

subroutine recalcModelView()

real, parameter :: pi = 3.1415926
  call glPopMatrix()
  call glPushMatrix()
  call glRotatef(angle%x, 0.0, 1.0, 0.0)
  call glRotatef(angle%y, cos(pi*angle%x/180.0), 0.0, sin(pi*angle%x/180.0))
  call glTranslatef(-8.0, -8.0, -bodyWidth / 2)
  newModel = .false.
end subroutine recalcModelView

subroutine redraw()

  if (newModel) then
    call recalcModelView()
  endif
  call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
  call glCallList(DINOSAUR)
  call glutSwapBuffers()
end subroutine redraw

subroutine mouse(button, state, x, y)
integer(kind=glcint), intent(in out) :: button, state, x, y

  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_DOWN) then
    moving = .true.
    begin = coord(x,y)
  endif
  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_UP) then
    moving = .false.
  endif
end subroutine mouse

subroutine motion(x, y)
integer(kind=glcint), intent(in out) :: x, y

  if (moving) then
    angle%x = angle%x + (x - begin%x)
    angle%y = angle%y + (y - begin%y)
    begin = coord(x,y)
    newModel = .true.
    call glutPostRedisplay()
  endif
end subroutine motion

subroutine controlLights(value)
integer(kind=glcint), intent(in out) :: value

  select case (value)
  case (1)
    lightZeroSwitch = .not. lightZeroSwitch
    if (lightZeroSwitch) then
      call glEnable(GL_LIGHT0)
    else
      call glDisable(GL_LIGHT0)
    endif
  case (2)
    lightOneSwitch = .not. lightOneSwitch
    if (lightOneSwitch) then
      call glEnable(GL_LIGHT1)
    else
      call glDisable(GL_LIGHT1)
    endif
  end select
  call glutPostRedisplay()
end subroutine controlLights

end module dino_mod

program main

use opengl_gl
use opengl_glu
use opengl_glut
use dino_mod
implicit none

integer :: i
real(kind=glfloat), dimension(4) :: lightZeroPosition = (/10.0, 4.0, 10.0, 1.0/), &
        lightZeroColor = (/0.8, 1.0, 0.8, 1.0/), & ! green-tinted
        lightOnePosition = (/-1.0, -2.0, 1.0, 0.0/), &
        lightOneColor = (/0.6, 0.3, 0.2, 1.0/) ! red-tinted

  call glutInit()
  call glutInitDisplayMode(ior(ior(GLUT_RGB,GLUT_DOUBLE),GLUT_DEPTH))
  i = glutCreateWindow("glutdino")
  call glutDisplayFunc(redraw)
  call glutMouseFunc(mouse)
  call glutMotionFunc(motion)
  call glutCreateMenu(controlLights,i)
  call glutAddMenuEntry("Toggle right light", 1)
  call glutAddMenuEntry("Toggle left light", 2)
  call glutAttachMenu(GLUT_RIGHT_BUTTON)
  call makeDinosaur()
  call glEnable(GL_CULL_FACE)
  call glEnable(GL_DEPTH_TEST)
  call glEnable(GL_LIGHTING)
  call glMatrixMode(GL_PROJECTION)
  call gluPerspective(40.0_gldouble, & ! field of view in degree
                          1.0_gldouble, & ! aspect ratio
                          1.0_gldouble, & ! Z near
                         40.0_gldouble)   ! Z far
  call glMatrixMode(GL_MODELVIEW)
  call gluLookAt( &
   0.0_gldouble, 0.0_gldouble, 30.0_gldouble, & ! eye is at (0,0,30)
   0.0_gldouble, 0.0_gldouble, 0.0_gldouble, & ! center is at (0,0,0)
   0.0_gldouble, 1.0_gldouble, 0.0_gldouble) ! up is in postivie Y direction
  call glPushMatrix()  ! dummy push so we can pop on model recalc 
  call glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1)
  call glLightfv(GL_LIGHT0, GL_POSITION, lightZeroPosition)
  call glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor)
  call glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1)
  call glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05)
  call glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition)
  call glLightfv(GL_LIGHT1, GL_DIFFUSE, lightOneColor)
  call glEnable(GL_LIGHT0)
  call glEnable(GL_LIGHT1)
  call glutMainLoop()
end program main
