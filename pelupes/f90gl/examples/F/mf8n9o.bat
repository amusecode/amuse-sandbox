REM This batch file has been tested on
REM   computer: Pentium II
REM         OS: Windows NT 4.0 (service pack 4)
REM   compiler: Essential Lahey F90 4.00b
REM          C: Microsoft Visual C++ 4.0
REM     OpenGL: Microsoft OpenGL, Glut 3.6 from Lahey binaries

REM ------------- User configuration parameters ---------------

REM modify these for your system

REM the directories containing the OpenGL libraries, f90gl libraries, GLUT
REM libraries, and f90gl GLUT libraries
set OGLLIBDIR=-libp c:\LaheyGL\lib_lf90

REM the fortran 90 libraries for OpenGL, including GLUT, GLU and OpenGL

set F90GLUTLIB=-l f90GL.lib,f90GLU.lib,f90GLUT.lib,glut32.lib,GLU32.lib,OPENGL32.lib

REM other libraries
set WINLIB=

REM the f90 compiler flag for specifying the location of MOD files
set MODS=-mod c:\LaheyGL\lib_lf90

REM fortran 90 compiler and compiler flags
set F90C=elf90
set F90FLAGS=-win -nomap

REM ----------- end of user configuration parameters ------------

REM select example to compile; default is all of them

if "%1"=="all" goto all
if "%1"=="blender" goto blender
if "%1"=="fbitfont" goto fbitfont
if "%1"=="fscene" goto fscene
if "%1"=="glutdino" goto glutdino
if "%1"=="logo" goto logo
if "%1"=="molehill" goto molehill
if "%1"=="olympic" goto olympic
if "%1"=="scube" goto scube
if "%1"=="sphere" goto sphere
if "%1"=="trivial" goto trivial
if "%1"=="modview" goto modview
if "%1"=="clean" goto clean

:all
call mf8n9o blender
call mf8n9o fbitfont
call mf8n9o fscene
call mf8n9o glutdino
call mf8n9o logo
call mf8n9o molehill
call mf8n9o olympic
REM call mf8n9o scube
call mf8n9o sphere
call mf8n9o trivial                                  
call mf8n9o modview                                  
goto done

REM blender

:blender
set 386link=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%F90C% blender.f90 %F90FLAGS% %MODS%
set 386link=
goto done

REM fbitfont

:fbitfont
set 386link=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%F90C% fbitfont.f90 %F90FLAGS% %MODS%
set 386link=
goto done

REM fscene

:fscene
set 386link=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%F90C% fscene.f90 %F90FLAGS% %MODS%
set 386link=
goto done

REM glutdino

:glutdino
set 386link=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%F90C% gltdino2.f90 %F90FLAGS% %MODS%
set 386link=
goto done

REM logo

:logo
set 386link=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%F90C% logo.f90 %F90FLAGS% %MODS%
set 386link=
goto done

REM molehill

:molehill
set 386link=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%F90C% molehill.f90 %F90FLAGS% %MODS%
set 386link=
goto done

REM olympic

:olympic
set 386link=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%F90C% olympic.f90 %F90FLAGS% %MODS%
set 386link=
goto done

REM scube

:scube
set 386link=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%F90C% scube.f90 %F90FLAGS% %MODS%
set 386link=
goto done

REM sphere

:sphere
set 386link=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%F90C% sphere.f90 %F90FLAGS% %MODS%
set 386link=
goto done

REM trivial

:trivial
set 386link=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%F90C% trivial.f90 %F90FLAGS% %MODS%
set 386link=
goto done

REM modview

:modview
set 386link=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%F90C% modview.f90 %F90FLAGS% %MODS%
set 386link=
goto done

REM clean

:clean
del *.obj *.mod *.exe *.map

:done
