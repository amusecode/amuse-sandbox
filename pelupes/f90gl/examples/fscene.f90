!  GLUT Fortran program to render simple red scene.

! This is a fortran 90 program in fixed source form.
! In a risky move, this program assumes that the default kind of
! integer is the same as glint and glsizei and also that the
! default real is the same kind as glfloat

        subroutine display
        use opengl_gl
        use opengl_glu
        use opengl_glut
        call glclear(GL_COLOR_BUFFER_BIT + GL_DEPTH_BUFFER_BIT)
        call glpushmatrix
        call glscalef(1.3, 1.3, 1.3)
        call glrotatef(20.0, 1.0, 0.0, 0.0)

        call glpushmatrix
        call gltranslatef(-0.75, 0.5, 0.0)
        call glrotatef(90.0, 1.0, 0.0, 0.0)
        call glutsolidtorus(0.275_gldouble, 0.85_gldouble, 10, 15)
        call glpopmatrix

        call glpushmatrix
        call gltranslatef(-0.75, -0.5, 0.0)
        call glrotatef(270.0, 1.0, 0.0, 0.0)
        call glutsolidtetrahedron
        call glpopmatrix

        call glpushmatrix
        call gltranslatef(0.75, 0.0, -1.0)
        call glutsolidicosahedron
        call glpopmatrix

        call glpopmatrix
        call glflush
        end

        subroutine reshape(w,h)
        use opengl_gl
        use opengl_glu
        use opengl_glut
        integer, intent(inout) :: w,h
        real(gldouble) wr,hr
        double precision d
        call glviewport(0, 0, w, h)
        call glmatrixmode(GL_PROJECTION)
        call glloadidentity
        wr = w
        hr = h
        d = 1.0d0
        if ( w .le. h ) then
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
        end
        
        subroutine submenu(value)
        use opengl_gl
        use opengl_glu
        use opengl_glut
        integer, intent(inout) :: value
        if ( value .eq. 1 ) then
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
        call glutpostredisplay
        end

        subroutine mainmenu(value)
        integer, intent(inout) :: value
        if (value == 666) stop
        end

        subroutine myinit
        use opengl_gl
        use opengl_glu
        use opengl_glut
        real lambient(4), ldiffuse(4), lspecular(4), lposition(4)
        data lambient /0.0, 0.0, 0.0, 1.0/
        data ldiffuse /1.0, 0.0, 0.0, 1.0/
        data lspecular /1.0, 1.0, 1.0, 1.0/
        data lposition /1.0, 1.0, 1.0, 0.0/

        call gllightfv(GL_LIGHT0, GL_AMBIENT, lambient)
        call gllightfv(GL_LIGHT0, GL_DIFFUSE, ldiffuse)
        call gllightfv(GL_LIGHT0, GL_SPECULAR, lspecular)
        call gllightfv(GL_LIGHT0, GL_POSITION, lposition)
        call glenable(GL_LIGHT0)
        call gldepthfunc(GL_LESS)
        call glenable(GL_DEPTH_TEST)
        call glenable(GL_LIGHTING)
        end

        program main
        use opengl_gl
        use opengl_glu
        use opengl_glut
        interface
          subroutine display()
          end subroutine display
          subroutine reshape(w,h)
          integer, intent(inout):: w,h
          end subroutine reshape
          subroutine submenu(value)
          integer, intent(inout):: value
          end subroutine submenu
          subroutine mainmenu(value)
          integer, intent(inout):: value
          end subroutine mainmenu
        end interface
        call glutinitwindowposition(500,500)
        call glutinitwindowsize(500,500)
        call glutinit
        i = glutcreatewindow('Fortran GLUT program')
        call myinit
        call glutdisplayfunc(display)
        call glutreshapefunc(reshape)
        i = glutcreatemenu(submenu)
        call glutaddmenuentry('Filled', 1)
        call glutaddmenuentry('Outline', 2)
        j = glutcreatemenu(mainmenu)
        call glutaddsubmenu('Polygon mode', i)
        call glutaddmenuentry('Quit', 666)
        call glutattachmenu(GLUT_RIGHT_BUTTON)
        call glutmainloop
        end

