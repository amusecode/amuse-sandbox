!  Example showing how to use OpenGL's feedback mode to capture
!  transformed vertices and output them as Encapsulated PostScript.
!  Handles limited hidden surface removal by sorting and does
!  smooth shading (albeit limited due to PostScript).

! Modified for use with PHAML graphics, William F. Mitchell, 5/15/97
! Translated to Fortran by William F. Mitchell, 6/16/00
! Extended to an example for f90gl by William F. Mitchell, 1/3/03
! Translated to F by William F. Mitchell, 1/23/04

! This is a translation of a C program written by Mark J. Kilgard.  The
! original C program had the following notice:

! /* Copyright (c) Mark J. Kilgard, 1997. */
! 
! /* This program is freely distributable without licensing fees 
!    and is provided without guarantee or warrantee expressed or 
!    implied. This program is -not- in the public domain. */

! (end of Mark's notice)

! The Fortran version is a contribution of NIST, not subject to copyright.

! The OpenGL image must be drawn in display list number 1, since this
! routine uses "call glCallList(1_GLint)" to render the image.

! The interface for the public routine is:

! subroutine outputEPS(bsize, doSort, filename)
! integer, intent(in) :: bsize
! logical, intent(in) :: doSort
! character(len=*), intent(in), optional :: filename

! bsize is the amount of space (in words) to allocate for the feedback buffer.
! doSort indicates whether or not to sort the rendered objects, back to front.
! With F, doSort is ignored (because my sort routine is an old FORTRAN 77
! routine from Netlib and I don't want to translate it).
! filename is the name of the file to which to write the postscript.  If it
! is not present, the feedback buffer is printed (for debugging)

module rendereps
use opengl_gl
implicit none
private
public :: outputEPS

real(kind=GLfloat), private, dimension(1) :: pointSize

character(len=74), private, dimension(25) :: gouraudtriangleEPS = (/ &
  "/bd{bind def}bind def /triangle { aload pop   setrgbcolor  aload pop 5 3  ",&
  "roll 4 2 roll 3 2 roll exch moveto lineto lineto closepath fill } bd      ",&
  "/computediff1 { 2 copy sub abs threshold ge {pop pop pop true} { exch 2   ",&
  "index sub abs threshold ge { pop pop true} { sub abs threshold ge } ifelse",&
  "} ifelse } bd /computediff3 { 3 copy 0 get 3 1 roll 0 get 3 1 roll 0 get  ",&
  "computediff1 {true} { 3 copy 1 get 3 1 roll 1 get 3 1 roll 1 get          ",&
  "computediff1 {true} { 3 copy 2 get 3 1 roll  2 get 3 1 roll 2 get         ",&
  "computediff1 } ifelse } ifelse } bd /middlecolor { aload pop 4 -1 roll    ",&
  "aload pop 4 -1 roll add 2 div 5 1 roll 3 -1 roll add 2 div 3 1 roll add 2 ",&
  "div 3 1 roll exch 3 array astore } bd /gouraudtriangle { computediff3 { 4 ",&
  "-1 roll aload 7 1 roll 6 -1 roll pop 3 -1 roll pop add 2 div 3 1 roll add ",&
  "2 div exch 3 -1 roll aload 7 1 roll exch pop 4 -1 roll pop add 2 div 3 1  ",&
  "roll add 2 div exch 3 -1 roll aload 7 1 roll pop 3 -1 roll pop add 2 div 3",&
  "1 roll add 2 div exch 7 3 roll 10 -3 roll dup 3 index middlecolor 4 1 roll",&
  "2 copy middlecolor 4 1 roll 3 copy pop middlecolor 4 1 roll 13 -1 roll    ",&
  "aload pop 17 index 6 index 15 index 19 index 6 index 17 index 6 array     ",&
  "astore 10 index 10 index 14 index gouraudtriangle 17 index 5 index 17     ",&
  "index 19 index 5 index 19 index 6 array astore 10 index 9 index 13 index  ",&
  "gouraudtriangle 13 index 16 index 5 index 15 index 18 index 5 index 6     ",&
  "array astore 12 index 12 index 9 index gouraudtriangle 17 index 16 index  ",&
  "15 index 19 index 18 index 17 index 6 array astore 10 index 12 index 14   ",&
  "index gouraudtriangle 18 {pop} repeat } { aload pop 5 3 roll aload pop 7 3",&
  "roll aload pop 9 3 roll 4 index 6 index 4 index add add 3 div 10 1 roll 7 ",&
  "index 5 index 3 index add add 3 div 10 1 roll 6 index 4 index 2 index add ",&
  "add 3 div 10 1 roll 9 {pop} repeat 3 array astore triangle } ifelse } bd  " &
  /)

! OpenGL's GL_3D_COLOR feedback vertex format. Use as offsets from base entry.

integer, private, parameter :: X=0, Y=1, Z=2, CRED=3, CGREEN=4, CBLUE=5, ALPHA=6

!type DepthIndex 
!  integer(kind=GLint) :: ptr
!  real(kind=GLfloat) :: depth
!end type DepthIndex

private :: print3DcolorVertex, spewPrimitiveEPS, printBuffer, &
           spewUnsortedFeedback, spewWireFrameEPS

contains

!          ------------------
subroutine print3DcolorVertex(count, buffer)
!          ------------------
integer(kind=GLint), intent(inout) :: count
real(kind=GLfloat), intent(in), dimension(:) :: buffer

! Write contents of one vertex to stdout.

!write(unit=*,fmt="('  ',7(f4.2,' '))") buffer(count:count+6)
write(unit=*,fmt="(a2,7(f4.2))") "  ",buffer(count:count+6)
count = count + 7

end subroutine print3DcolorVertex

!          -----------
subroutine printBuffer(bsize, buffer)
!          -----------
integer, intent(in) :: bsize
real(kind=GLfloat),intent(in),dimension(:) :: buffer

integer(kind=GLint) :: count
integer(kind=GLint) :: token
integer :: nvertices, i

  count = 1
  do
    if (count >= bsize) exit
    token = buffer(count)
    count = count + 1
    select case (token)
    case (GL_PASS_THROUGH_TOKEN)
      write(unit=*,fmt=*) "GL_PASS_THROUGH_TOKEN"
      write(unit=*,fmt="(a2,f4.2)") "  ",buffer(count)
      count = count + 1
    case (GL_POINT_TOKEN)
      write(unit=*,fmt=*) "GL_POINT_TOKEN"
      call print3DcolorVertex(count, buffer)
    case (GL_LINE_TOKEN)
      write(unit=*,fmt=*) "GL_LINE_TOKEN"
      call print3DcolorVertex(count, buffer)
      call print3DcolorVertex(count, buffer)
    case (GL_LINE_RESET_TOKEN)
      write(unit=*,fmt=*) "GL_LINE_RESET_TOKEN"
      call print3DcolorVertex(count, buffer)
      call print3DcolorVertex(count, buffer)
    case (GL_POLYGON_TOKEN)
      write(unit=*,fmt=*) "GL_POLYGON_TOKEN"
      nvertices = buffer(count)
      count = count + 1
      do i=1,nvertices
        call print3DcolorVertex(count, buffer)
      end do
    case default
      write(unit=*,fmt=*) "Incomplete implementation.  Unexpected token ",token
    end select
  end do
end subroutine printBuffer

!          ----------------
subroutine spewPrimitiveEPS(file, loc, buffer)
!          ----------------
integer, intent(in) :: file
integer, intent(inout) :: loc
real(kind=GLfloat), intent(in),dimension(:) :: buffer

  integer(kind=GLint) :: token
  integer :: nvertices, i
  real(kind=GLfloat) :: red, green, blue
  logical :: smooth
  real(kind=GLfloat) :: dx, dy, dr, dg, db, absR, absG, absB, colormax
  integer :: steps
  real(kind=GLfloat) :: xstep, ystep, rstep, gstep, bstep
  real(kind=GLfloat) :: xnext, ynext, rnext, gnext, bnext, distance
! Lower for better smooth lines.
  real(kind=GLfloat), parameter :: EPS_SMOOTH_LINE_FACTOR = 0.06

  token = buffer(loc)
  loc = loc + 1
  select case (token)
  case (GL_LINE_RESET_TOKEN, GL_LINE_TOKEN)

    dr = buffer(loc+7+CRED) - buffer(loc+CRED)
    dg = buffer(loc+7+CGREEN) - buffer(loc+CGREEN)
    db = buffer(loc+7+CBLUE) - buffer(loc+CBLUE)

    if (dr /= 0 .or. dg /= 0 .or. db /= 0) then
!        Smooth shaded line.
      dx = buffer(loc+7+X) - buffer(loc+X)
      dy = buffer(loc+7+Y) - buffer(loc+Y)
      distance = sqrt(dx * dx + dy * dy)

      absR = abs(dr)
      absG = abs(dg)
      absB = abs(db)

      colormax = max(absR, absG, absB)
      steps = max(1.0_GLfloat, colormax * distance * EPS_SMOOTH_LINE_FACTOR)

      xstep = dx / steps
      ystep = dy / steps

      rstep = dr / steps
      gstep = dg / steps
      bstep = db / steps

      xnext = buffer(loc+X)
      ynext = buffer(loc+Y)
      rnext = buffer(loc+CRED)
      gnext = buffer(loc+CGREEN)
      bnext = buffer(loc+CBLUE)

!       Back up half a step; we want the end points to be
!       exactly the their endpoint colors.
      xnext = xnext - xstep / 2.0_GLfloat
      ynext = ynext - ystep / 2.0_GLfloat
      rnext = rnext - rstep / 2.0_GLfloat
      gnext = gnext - gstep / 2.0_GLfloat
      bnext = bnext - bstep / 2.0_GLfloat
    else
!       Single color line.
      steps = 0
    endif

    write(unit=file,fmt=*) buffer(loc+CRED),buffer(loc+CGREEN),buffer(loc+CBLUE),&
                           " setrgbcolor"
    write(unit=file,fmt=*) buffer(loc+X),buffer(loc+Y)," moveto"

    do i=1,steps
      xnext = xnext + xstep
      ynext = ynext + ystep
      rnext = rnext + rstep
      gnext = gnext + gstep
      bnext = bnext + bstep
      write(unit=file,fmt=*) xnext, ynext, " lineto stroke"
      write(unit=file,fmt=*) rnext, gnext, bnext, " setrgbcolor"
      write(unit=file,fmt=*) xnext, ynext, " moveto"
    end do
    write(unit=file,fmt=*) buffer(loc+7+X), buffer(loc+7+Y), " lineto stroke"

    loc = loc + 14  ! Each vertex element in the feedback buffer is 7 GLfloats.

  case (GL_POLYGON_TOKEN)
    nvertices = buffer(loc)
    loc = loc + 1

    if (nvertices > 0) then
      red = buffer(loc+CRED)
      green = buffer(loc+CGREEN)
      blue = buffer(loc+CBLUE)
      smooth = .false.
      do i=1,nvertices-1
        if (red /= buffer(loc+7*i+CRED) .or. green /= buffer(loc+7*i+CGREEN) &
            .or. blue /= buffer(loc+7*i+CBLUE)) then
          smooth = .true.
          exit
        endif
      end do
      if (smooth) then
!         Smooth shaded polygon; varying colors at vetices.

!         Break polygon into "nvertices-2" triangle fans.
        do i=1,nvertices-2
          write(unit=file,fmt=*) "[",buffer(loc+X),buffer(loc+7*i+X), &
             buffer(loc+7*(i+1)+X),buffer(loc+Y),buffer(loc+7*i+Y), &
             buffer(loc+7*(i+1)+Y),"] [",buffer(loc+CRED),buffer(loc+CGREEN), &
             buffer(loc+CBLUE),"] [",buffer(loc+7*i+CRED),buffer(loc+7*i+CGREEN), &
             buffer(loc+7*i+CBLUE),"] [",buffer(loc+7*(i+1)+CRED), &
             buffer(loc+7*(i+1)+CGREEN),buffer(loc+7*(i+1)+CBLUE), &
             "] gouraudtriangle"
        end do
      else
!         Flat shaded polygon; all vertex colors the same.
        write(unit=file,fmt=*) "newpath"
        write(unit=file,fmt=*) red, green, blue, " setrgbcolor"

!         Draw a filled triangle.
        write(unit=file,fmt=*) buffer(loc+X),buffer(loc+Y)," moveto"
        do i=1,nvertices-1
          write(unit=file,fmt=*) buffer(loc+7*i+X),buffer(loc+7*i+Y)," lineto"
        end do
        write(unit=file,fmt=*) "closepath fill"
        write(unit=file,fmt=*) ""
      endif
    endif
    loc = loc + nvertices * 7  ! Each vertex element in the
                               !  feedback buffer is 7 GLfloats.
  case (GL_POINT_TOKEN)
    write(unit=file,fmt=*) buffer(loc+CRED),buffer(loc+CGREEN),buffer(loc+CBLUE), &
                           " setrgbcolor"
    write(unit=file,fmt=*) buffer(loc+X),buffer(loc+Y),pointSize(1)/2.0,0, &
                           360," arc fill"
    write(unit=file,fmt=*) ""
    loc = loc + 7  ! Each vertex element in the feedback buffer is 7 GLfloats.

  case default
    write(unit=*,fmt=*) "Incomplete implementation.  Unexpected token ",token

  end select
end subroutine spewPrimitiveEPS

!          --------------------
subroutine spewUnsortedFeedback(file, bsize, buffer)
!          --------------------
integer, intent(in) :: file, bsize
real(kind=GLfloat), intent(in), dimension(:) :: buffer

  integer :: loc

  loc = 1
  do
    if (loc >= bsize) exit
    call spewPrimitiveEPS(file, loc, buffer)
  end do
end subroutine spewUnsortedFeedback

!!          ------------------
!subroutine spewSortedFeedback(file, bsize, buffer)
!!          ------------------
!integer, intent(in) :: file, bsize
!real(GLfloat), intent(in) :: buffer(:)
!
!  integer(GLint) :: token
!  integer :: loc
!  real(GLfloat) :: depthSum
!  integer :: nprimitives, item
!  type(DepthIndex), allocatable :: prims(:)
!  integer :: nvertices, i
!  real :: ydum(1)
!
!!   Count how many primitives there are.
!  nprimitives = 0
!  loc = 1
!  do while (loc < bsize)
!    token = buffer(loc)
!    loc = loc + 1
!    select case (token)
!    case (GL_LINE_TOKEN, GL_LINE_RESET_TOKEN)
!      loc = loc + 14
!      nprimitives = nprimitives + 1
!    case (GL_POLYGON_TOKEN)
!      nvertices = buffer(loc)
!      loc = loc + 1
!      loc = loc + 7*nvertices
!      nprimitives = nprimitives + 1
!    case (GL_POINT_TOKEN)
!      loc = loc + 7
!      nprimitives = nprimitives + 1
!    case default
!      write(unit=*,fmt=*) "Incomplete implementation.  Unexpected token ",token
!    end select
!  end do
!
!!    Allocate an array of pointers that will point back at
!!    primitives in the feedback buffer.  There will be one
!!    entry per primitive.  This array is also where we keep the
!!    primitive's average depth.  There is one entry per
!!    primitive  in the feedback buffer.
!  allocate(prims(nprimitives))
!
!  item = 1
!  loc = 1
!  do while (loc < bsize)
!    prims(item)%ptr = loc  ! Save this primitive's location.
!    token = buffer(loc)
!    loc = loc + 1
!    select case (token)
!    case (GL_LINE_TOKEN, GL_LINE_RESET_TOKEN)
!      depthSum = buffer(loc+Z) + buffer(loc+7+Z)
!      prims(item)%depth = depthSum / 2.0
!! WFM to force triangle edges on top of filled triangle
!!      prims(item)%depth = prims(item)%depth - .0001
!      loc = loc + 14
!    case (GL_POLYGON_TOKEN)
!      nvertices = buffer(loc)
!      loc = loc + 1
!      depthSum = buffer(loc+Z)
!      do i=1,nvertices-1
!        depthSum = depthSum + buffer(loc+7*i+Z)
!      end do
!      prims(item)%depth = depthSum / nvertices
!      loc = loc + 7*nvertices
!    case (GL_POINT_TOKEN)
!      prims(item)%depth = buffer(loc+Z)
!      loc = loc + 7
!    case default
!      write(unit=*,fmt=*) "Incomplete implementation.  Unexpected token ",token
!    end select
!    item = item + 1
!  end do
!
!!  Sort the primitives back to front.
!!  call ssort(prims, ydum, nprimitives, -1)
!
!!    XXX Understand that sorting by a primitives average depth
!!    doesn't allow us to disambiguate some cases like self
!!    intersecting polygons.  Handling these cases would require
!!    breaking up the primitives.  That's too involved for this
!!    example.  Sorting by depth is good enough for lots of
!!    applications.
!
!!  Emit the Encapsulated PostScript for the primitives in
!!  back to front order.
!  do item=1,nprimitives
!    call spewPrimitiveEPS(file, prims(item)%ptr, buffer)
!  end do
!
!  deallocate(prims)
!end subroutine spewSortedFeedback

!          ----------------
subroutine spewWireFrameEPS(file, doSort, bsize, buffer, creator)
!          ----------------
integer, intent(in) :: file, bsize
logical, intent(in) :: doSort
real(kind=GLfloat), intent(in), dimension(:) :: buffer
character(len=*), intent(in) :: creator

  real(kind=GLfloat), dimension(4) :: clearColor, viewport
  real(kind=GLfloat), dimension(1) :: lineWidth
  integer :: i
! Lower for better (slower) smooth shading.
  real(kind=GLfloat), parameter :: EPS_GOURAUD_THRESHOLD=0.1

!    Read back a bunch of OpenGL state to help make the EPS
!    consistent with the OpenGL clear color, line width, point
!    bsize, and viewport.
  call glGetFloatv(GL_VIEWPORT, viewport)
  call glGetFloatv(GL_COLOR_CLEAR_VALUE, clearColor)
  call glGetFloatv(GL_LINE_WIDTH, lineWidth)
  call glGetFloatv(GL_POINT_SIZE, pointSize)

!    Emit EPS header. */
  write(unit=file,fmt="(a)") "%!PS-Adobe-2.0 EPSF-2.0"
  write(unit=file,fmt="(a)") "%%Creator: eps.f90 (using OpenGL feedback)"
  write(unit=file,fmt="(a,4i10)") "%%BoundingBox: ", int(viewport)
  write(unit=file,fmt="(a)") "%%EndComments"
  write(unit=file,fmt=*) ""
  write(unit=file,fmt=*) "gsave"
  write(unit=file,fmt=*) ""

!   Output Frederic Delhoume's "gouraudtriangle" PostScript fragment.
  write(unit=file,fmt=*) "% the gouraudtriangle PostScript fragement below is free"
  write(unit=file,fmt=*) "% written by Frederic Delhoume (delhoume@ilog.fr)"
  write(unit=file,fmt=*) "/threshold ",EPS_GOURAUD_THRESHOLD," def"
  do i=1,size(gouraudtriangleEPS)
    write(unit=file,fmt=*) gouraudtriangleEPS(i)
  end do

  write(unit=file,fmt=*) ""
  write(unit=file,fmt=*) lineWidth(1)," setlinewidth"

!   Clear the background like OpenGL had it.
  write(unit=file,fmt=*) clearColor(1:3)," setrgbcolor"
  write(unit=file,fmt=*) viewport," rectfill"
  write(unit=file,fmt=*) ""

!  if (doSort) then
!    call spewSortedFeedback(file, bsize, buffer)
!  else
    call spewUnsortedFeedback(file, bsize, buffer)
!  endif

!   Emit EPS trailer.
  write(unit=file,fmt=*) "grestore"
  write(unit=file,fmt=*) ""
! WFM I don't know why Mark put this comment instead of showpage, unless some
!     apps that take eps input don't like it (ghostview is fine with it)
!  write(unit=file,fmt=*) "%Add `showpage' to the end of this file to be able to print to a printer."

  write(unit=file,fmt=*) "% 'showpage' is needed to print to a printer, but there"
  write(unit=file,fmt=*) "% may be some apps for which it must be removed"
  write(unit=file,fmt=*) "showpage"

  close(file)
end subroutine spewWireFrameEPS

subroutine outputEPS(bsize, doSort, filename)
integer, intent(in) :: bsize
logical, intent(in) :: doSort
character(len=*), intent(in), optional :: filename

  real(kind=GLfloat), allocatable, dimension(:) :: feedbackBuffer
  integer(kind=GLint) :: returned, idum
  integer :: file, iostat
  logical :: opened

  allocate(feedbackBuffer(bsize))
  call glFeedbackBuffer(bsize, GL_3D_COLOR, feedbackBuffer)
  idum = glRenderMode(GL_FEEDBACK)
  call glCallList(1_GLint)
  returned = glRenderMode(GL_RENDER)
  if (present(filename)) then
    file = 11
    do
      inquire(unit=file,opened=opened)
      if (.not. opened) exit
      file = file + 1
    end do
    open(unit=file,file=filename,iostat=iostat,status="NEW",action="WRITE")
    if (iostat==0) then
      call spewWireFrameEPS(file, doSort, returned, feedbackBuffer, "rendereps")
    else
      write(unit=*,fmt=*) "Could not open ", filename
    endif
  else
!     Helps debugging to be able to see the decode feedback buffer as text.
    call printBuffer(returned, feedbackBuffer)
  endif
  deallocate(feedbackBuffer)
end subroutine outputEPS

end module rendereps

!-----------------------------------------------------------
!
! Use the fscene.f90 program to illustrate use of rendereps.  Contrast the
! remainder of this file to the file fscene.f90.

!  GLUT Fortran program to render simple red scene.

! In a risky move, this program assumes that the default kind of
! integer is the same as glint and glsizei and also that the
! default real is the same kind as glfloat

module callbacks

private
public ::  mainmenu, display, myinit, submenu

contains

        subroutine display()
        use opengl_gl
        use opengl_glu
        use opengl_glut
        call gldeletelists(1_gluint, 1_glsizei)
! put the display in list 1 for rendereps
        call glnewlist(1_gluint, gl_compile_and_execute)
        call glclear(GL_COLOR_BUFFER_BIT + GL_DEPTH_BUFFER_BIT)
        call glpushmatrix()
        call glscalef(1.3, 1.3, 1.3)
        call glrotatef(20.0, 1.0, 0.0, 0.0)

        call glpushmatrix()
        call gltranslatef(-0.75, 0.5, 0.0)
        call glrotatef(90.0, 1.0, 0.0, 0.0)
        call glutsolidtorus(0.275_gldouble, 0.85_gldouble, 10, 15)
        call glpopmatrix()

        call glpushmatrix()
        call gltranslatef(-0.75, -0.5, 0.0)
        call glrotatef(270.0, 1.0, 0.0, 0.0)
        call glutsolidtetrahedron()
        call glpopmatrix()

        call glpushmatrix()
        call gltranslatef(0.75, 0.0, -1.0)
        call glutsolidicosahedron()
        call glpopmatrix()

        call glpopmatrix()
        call glendlist()
        call glflush()
        end subroutine display

        subroutine reshape(w,h)
        use opengl_gl
        use opengl_glu
        use opengl_glut
        integer, intent(in out) :: w,h
        real(kind=gldouble) :: wr,hr,d
        call glviewport(0, 0, w, h)
        call glmatrixmode(GL_PROJECTION)
        call glloadidentity()
        wr = w
        hr = h
        d = 1.0_gldouble
        if ( w <= h ) then
           call glortho(-2.5_gldouble, 2.5_gldouble,      &
             -2.5_gldouble * hr/wr, 2.5_gldouble * hr/wr, &
             -10.0_gldouble, 10.0_gldouble)
        else
           call glortho(-2.5_gldouble * hr/wr,            &
              2.5_gldouble * hr/wr,                       &
             -2.5_gldouble, 2.5_gldouble, -10.0_gldouble, &
             10.0_gldouble)
        end if
        call glmatrixmode(GL_MODELVIEW)
        end subroutine reshape
        
        subroutine submenu(value)
        use opengl_gl
        use opengl_glu
        use opengl_glut
        integer, intent(in out) :: value
        if ( value == 1 ) then
          call glenable(GL_DEPTH_TEST)
          call glenable(GL_LIGHTING)
          call gldisable(GL_BLEND)
          call glpolygonmode(GL_FRONT_AND_BACK, GL_FILL)
        else
          call gldisable(GL_DEPTH_TEST)
          call gldisable(GL_LIGHTING)
          call glcolor3f(1.0, 1.0, 1.0)
          call glpolygonmode(GL_FRONT_AND_BACK, GL_LINE)
          call glenable(GL_LINE_SMOOTH)
          call glenable(GL_BLEND)
          call glblendfunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        end if
        call glutpostredisplay()
        end subroutine submenu

        subroutine mainmenu(value)
        use rendereps
        use opengl_glut
        integer, intent(in out) :: value
        if (value == 10) then
           call glutsetcursor(glut_cursor_wait)
           call outputeps(100000,.true.,"epsout.eps")
           call glutsetcursor(glut_cursor_inherit)
        endif
        if (value == 666) stop
        end subroutine mainmenu

        subroutine myinit()
        use opengl_gl
        use opengl_glu
        use opengl_glut
        real, dimension(4) :: lambient = (/ 0.0, 0.0, 0.0, 1.0 /), &
                              ldiffuse = (/ 1.0, 0.0, 0.0, 1.0 /), &
                              lspecular= (/ 1.0, 1.0, 1.0, 1.0 /), &
                              lposition= (/ 1.0, 1.0, 1.0, 0.0 /)

        call gllightfv(GL_LIGHT0, GL_AMBIENT, lambient)
        call gllightfv(GL_LIGHT0, GL_DIFFUSE, ldiffuse)
        call gllightfv(GL_LIGHT0, GL_SPECULAR, lspecular)
        call gllightfv(GL_LIGHT0, GL_POSITION, lposition)
        call glenable(GL_LIGHT0)
        call gldepthfunc(GL_LESS)
        call glenable(GL_DEPTH_TEST)
        call glenable(GL_LIGHTING)
        end subroutine myinit

end module callbacks

        program main
        use opengl_gl
        use opengl_glu
        use opengl_glut
        use callbacks
        integer :: i,j
        call glutinitwindowposition(500,500)
        call glutinitwindowsize(500,500)
        call glutinit()
        i = glutcreatewindow("Fortran GLUT program")
        call myinit()
        call glutdisplayfunc(display)
        call glutreshapefunc(reshape)
        call glutcreatemenu(submenu,i)
        call glutaddmenuentry("Filled", 1)
        call glutaddmenuentry("Outline", 2)
        call glutcreatemenu(mainmenu,j)
        call glutaddsubmenu("Polygon mode", i)
        call glutaddmenuentry("Save postscript",10)
        call glutaddmenuentry("Quit", 666)
        call glutattachmenu(GLUT_RIGHT_BUTTON)
        call glutmainloop()
        end program main
