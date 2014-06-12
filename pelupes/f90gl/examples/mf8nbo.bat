@echo off
REM This batch file has been tested on
REM   computer: Pentium III
REM         OS: Windows 2000
REM   compiler: Absoft Pro Fortran 9.0
REM          C: Microsoft Visual C++ 6.0
REM     OpenGL: Microsoft OpenGL (Absoft import lib), Glut 3.7.1, f90gl 1.2.9

REM ------------- User configuration parameters ---------------

REM modify these for your system

REM the directories containing the OpenGL libraries, f90gl libraries, GLUT
REM libraries, and f90gl GLUT libraries
set OGLLIBDIR=

REM the directories containing the Windows libraries
set WINLIBDIR=

REM the fortran 90 libraries for OpenGL, including GLUT, GLU and OpenGL
set F90GLUTLIB= f90glut.lib f90glu.lib f90gl.lib glut32.lib glu32.lib opengl32.lib

REM the Windows libraries
set WINLIB= absRT0.lib fio.lib fmath.lib f90math.lib comdlg32.lib libac.lib kernel32.lib

REM the f90 compiler flag for specifying the location of MOD files
set MODS=-p ..\include\GL

REM fortran 90 compiler and compiler flags
set F90C=f90
set F90FLAGS=-O

REM fortran 90 compiler flag for fixed source form
set FFIXED=-f fixed

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
call mf8nbo blender
call mf8nbo checker
call mf8nbo eps
call mf8nbo fbitfont
call mf8nbo fscene
call mf8nbo glutdino
call mf8nbo logo
call mf8nbo molehill
call mf8nbo olympic
call mf8nbo plotfunc
call mf8nbo polyoff
call mf8nbo scube
call mf8nbo sphere
call mf8nbo stars
call mf8nbo trivial
call mf8nbo modview
goto done

REM blender

:blender
%F90C% -c blender.f90 %F90FLAGS% %MODS%
link -force:multiple blender.obj %F90GLUTLIB% %WINLIB%
goto done

REM checker

:checker
%F90C% -c checker.f90 %F90FLAGS% %MODS%
link -force:multiple checker.obj %F90GLUTLIB% %WINLIB%
goto done

REM eps

:eps
%F90C% -c eps.f90 %F90FLAGS% %MODS%
link -force:multiple eps.obj %F90GLUTLIB% %WINLIB%
goto done

REM fbitfont

:fbitfont
%F90C% -c fbitfont.f90 %F90FLAGS% %MODS%
link -force:multiple fbitfont.obj %F90GLUTLIB% %WINLIB%
goto done

REM fscene

:fscene
%F90C% -v -c fscene.f90 %F90FLAGS% %MODS%
link -force:multiple fscene.obj %F90GLUTLIB% %WINLIB%
goto done

REM glutdino

:glutdino
%F90C% -c glutdino.f90 %F90FLAGS% %MODS%
link -force:multiple glutdino.obj %F90GLUTLIB% %WINLIB%
goto done

REM logo

:logo
%F90C% -c logo.f90 %F90FLAGS% %MODS%
link -force:multiple logo.obj %F90GLUTLIB% %WINLIB%
goto done

REM molehill

:molehill
%F90C% -c molehill.f90 %F90FLAGS% %MODS%
link -force:multiple molehill.obj %F90GLUTLIB% %WINLIB%
goto done

REM olympic

:olympic
REM %F90C% olympic.f90 %F90FLAGS% %MODS%
%F90C% -c olympic.f90 %F90FLAGS% %MODS%
link -force:multiple olympic.obj %F90GLUTLIB% %WINLIB%
goto done

REM plotfunc

:plotfunc
REM %F90C% plotfunc.f90 %F90FLAGS% %MODS%
%F90C% -c plotfunc.f90 %F90FLAGS% %MODS%
link -force:multiple plotfunc.obj %F90GLUTLIB% %WINLIB%
goto done

REM polyoff

:polyoff
REM %F90C% polyoff.f90 %F90FLAGS% %MODS%
%F90C% -c polyoff.f90 %F90FLAGS% %MODS%
link -force:multiple polyoff.obj %F90GLUTLIB% %WINLIB%
goto done

REM scube

:scube
%F90C% -c scube.f90 %F90FLAGS% %MODS%
link -force:multiple scube.obj %F90GLUTLIB% %WINLIB%
goto done

REM sphere

:sphere
%F90C% -c sphere.f90 %F90FLAGS% %MODS%
link -force:multiple sphere.obj %F90GLUTLIB% %WINLIB%
goto done

REM stars

:stars
%F90C% -c stars.f90 %F90FLAGS% %MODS%
link -force:multiple stars.obj %F90GLUTLIB% %WINLIB%
goto done

REM trivial

:trivial
%F90C% -c trivial.f90 %F90FLAGS% %MODS%
REM %F90C% trivial.f90 %F90FLAGS% %MODS%
link -force:multiple trivial.obj %F90GLUTLIB% %WINLIB%
goto done

REM modview

:modview
REM %F90C% modview.f90 %F90FLAGS% %MODS%
%F90C% -c modview.f90 %F90FLAGS% %MODS%
link -force:multiple modview.obj %F90GLUTLIB% %WINLIB%
goto done

REM clean

:clean
del *.obj *.mod *.exe *.f9 fppr.inc *.exp *.lib

:done
