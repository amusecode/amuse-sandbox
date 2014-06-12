@echo off
REM This batch file has been tested on
REM   computer: Pentium II
REM         OS: Windows NT 4.0 (service pack 4)
REM   compiler: Salford FTN95
REM          C: Microsoft Visual C++ 6.0
REM     OpenGL: Microsoft OpenGL, Glut 3.7.1, f90gl 1.1.4

REM NOTE f90gl does not currently work with Salford FTN95

REM ------------- User configuration parameters ---------------

REM modify these for your system

REM the directories containing the OpenGL libraries, f90gl libraries, GLUT
REM libraries, and f90gl GLUT libraries
set OGLLIBDIR=
set LIBSAV=%LIB%
set LIB=%LIB%;..\lib

REM the fortran 90 libraries for OpenGL, including GLUT, GLU and OpenGL

set F90GLUTLIB= f90GL.lib f90GLU.lib f90GLUT.lib glut32.lib GLU32.lib OPENGL32.lib

REM other libraries
set WINLIB=

REM the f90 compiler flag for specifying the location of MOD files
set MODS=/MOD_PATH ..\lib

REM fortran 90 compiler and compiler flags
set F90C=ftn95
set F90FLAGS= /FPS_EXTERNAL

REM fortran 90 compiler flag for fixed source form
set FFIXED=

REM Linker
set LINKER=SLINK
 
REM ----------- end of user configuration parameters ------------

REM select example to compile; default is all of them

if "%1"=="all" goto all
if "%1"=="blender" goto blender
if "%1"=="checker" goto checker
if "%1"=="eps" goto eps
if "%1"=="fbitfont" goto fbitfont
if "%1"=="fscene" goto fscene
if "%1"=="glutdino" goto glutdino
if "%1"=="logo" goto logo
if "%1"=="molehill" goto molehill
if "%1"=="olympic" goto olympic
if "%1"=="plotfunc" goto plotfunc
if "%1"=="polyoff" goto polyoff
if "%1"=="scube" goto scube
if "%1"=="sphere" goto sphere
if "%1"=="stars" goto stars
if "%1"=="trivial" goto trivial
if "%1"=="modview" goto modview
if "%1"=="clean" goto clean

:all
call mf8noo blender
call mf8noo checker
call mf8noo eps
call mf8noo fbitfont
call mf8noo fscene
call mf8noo glutdino
call mf8noo logo
call mf8noo molehill
call mf8noo olympic
call mf8noo plotfunc
call mf8noo polyoff
call mf8noo scube
call mf8noo sphere
call mf8noo stars
call mf8noo trivial
call mf8noo modview
goto done

REM blender

:blender
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  blender.f90 %F90FLAGS% %MODS%
%LINKER% blender.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM checker

:checker
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  checker.f90 %F90FLAGS% %MODS%
%LINKER% checker.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM eps

:eps
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  eps.f90 %F90FLAGS% %MODS%
%LINKER% eps.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM fbitfont

:fbitfont
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  fbitfont.f90 %F90FLAGS% %MODS%
%LINKER% fbitfont.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM fscene

:fscene
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  fscene.f90 %F90FLAGS% %MODS%
%LINKER% fscene.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM glutdino

:glutdino
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  glutdino.f90 %F90FLAGS% %MODS%
%LINKER% glutdino.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM logo

:logo
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  logo.f90 %F90FLAGS% %MODS%
%LINKER% logo.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM molehill

:molehill
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  molehill.f90 %F90FLAGS% %MODS%
%LINKER% molehill.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM olympic

:olympic
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  olympic.f90 %F90FLAGS% %MODS%
%LINKER% olympic.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM plotfunc

:plotfunc
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  plotfunc.f90 %F90FLAGS% %MODS%
%LINKER% plotfunc.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM polyoff

:polyoff
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  polyoff.f90 %F90FLAGS% %MODS%
%LINKER% polyoff.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM scube

:scube
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  scube.f90 %F90FLAGS% %MODS%
%LINKER% scube.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM sphere

:sphere
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  sphere.f90 %F90FLAGS% %MODS%
%LINKER% sphere.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM stars

:stars
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  stars.f90 %F90FLAGS% %MODS%
%LINKER% stars.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM trivial

:trivial
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  trivial.f90 %F90FLAGS% %MODS%
%LINKER% trivial.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM modview

:modview
set LINKLIB=%OGLLIBDIR% %F90GLUTLIB% %WINLIB%
%f90C%  modview.f90 %F90FLAGS% %MODS%
%LINKER% modview.obj %LINKLIB% >cc
rem set LINKLIB=
set 386link=
goto done

REM clean

:clean
del *.obj *.mod *.exe *.map

:done
set LIB=%LIBSAV%

