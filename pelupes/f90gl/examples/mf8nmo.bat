@echo off
REM This batch file has been tested on
REM   computer: Pentium
REM         OS: Windows NT 4.0
REM   compiler: Microsoft Fortran PowerStation 4.0
REM          C: Microsoft Visual C++ 4.0
REM     OpenGL: Microsoft OpenGL, Glut 3.6 with patch, f90gl 1.1.?
REM contributed by: Nicholas Yue yue@cwr.uwa.edu.au

REM ------------- User configuration parameters ---------------

REM modify these for your system

REM the directories containing the OpenGL libraries, f90gl libraries, GLUT
REM libraries, and f90gl GLUT libraries
set OGLLIBDIR=

REM the directories containing the Windows libraries
set WINLIBDIR=c:\MSDEV\lib

REM the fortran 90 libraries for OpenGL, including GLUT, GLU and OpenGL
set F90GLUTLIB= gdi32.lib user32.lib glut32.lib glu32.lib opengl32.lib ..\lib\f90GL.lib ..\lib\f90GLU.lib ..\lib\f90GLUT.lib 

REM the Windows libraries
set WINLIB=

REM the f90 compiler flag for specifying the location of MOD files
set MODS=-I"..\lib"

REM fortran 90 compiler and compiler flags
set F90C=fl32
set F90FLAGS=

REM fortran 90 compiler flag for fixed source form
set FFIXED=

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
call mf8nmo blender
call mf8nmo checker
call mf8nmo eps
call mf8nmo fbitfont
call mf8nmo fscene
call mf8nmo glutdino
call mf8nmo logo
call mf8nmo molehill
call mf8nmo olympic
call mf8nmo plotfunc
call mf8nmo polyoff
call mf8nmo scube
call mf8nmo sphere
call mf8nmo stars
call mf8nmo trivial
call mf8nmo modview
goto done

REM blender

:blender
%F90C% blender.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM checker

:checker
%F90C% checker.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM eps

:eps
%F90C% eps.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM fbitfont

:fbitfont
%F90C% fbitfont.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM fscene

:fscene
%F90C% fscene.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM glutdino

:glutdino
%F90C% glutdino.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM logo

:logo
%F90C% logo.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM molehill

:molehill
%F90C% molehill.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM olympic

:olympic
%F90C% olympic.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM plotfunc

:plotfunc
%F90C% plotfunc.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM polyoff

:polyoff
%F90C% polyoff.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM scube

:scube
%F90C% scube.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM sphere

:sphere
%F90C% sphere.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM stars

:stars
%F90C% stars.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM trivial

:trivial
%F90C% trivial.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM modview

:modview
%F90C% modview.f90 %F90FLAGS% %MODS% %OGLLIBDIR% %F90GLUTLIB% %WINLIB%
goto done

REM clean

:clean
del *.obj *.mod *.exe *.f9 fppr.inc

:done
