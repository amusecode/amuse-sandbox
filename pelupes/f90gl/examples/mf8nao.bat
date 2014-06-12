@echo off
REM This batch file has been tested on
REM   computer: Pentium II
REM         OS: Windows NT 4.0 (service pack 4)
REM   compiler: N.A.Software FortranPlus version 2.0.55r
REM          C: Microsoft Visual C++ 6.0
REM     OpenGL: Microsoft OpenGL, Glut 3.7.1, f90gl 1.2.0

REM ------------- User configuration parameters ---------------
REM ****
REM **** To run the demos, make sure that the f90gl\lib directory
REM **** is on the PATH and the LIB environments.
REM ****

REM modify these for your system

REM the directories containing the OpenGL libraries, f90gl libraries, GLUT
REM libraries, and f90gl GLUT libraries
set OGLLIBDIR=

REM the directories containing the Windows libraries
set WINLIBDIR=

REM the fortran 90 libraries for OpenGL, including GLUT, GLU and OpenGL
set F90GLUTLIB=

REM the Windows libraries
set WINLIB=

REM the f90 compiler flag for specifying the location of MOD files
set MODS=

REM Use the f90glut module
f95 -use f90glut

REM fortran 90 compiler and compiler flags
set F90C=f95
set F90FLAGS=-fast -mslink

REM fortran 90 compiler flag for fixed source form
set FFIXED=-fixed

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
call mf8nao blender
call mf8nao checker
call mf8nao eps
call mf8nao fbitfont
call mf8nao fscene
call mf8nao glutdino
call mf8nao logo
call mf8nao molehill
call mf8nao olympic
call mf8nao plotfunc
call mf8nao polyoff
call mf8nao scube
call mf8nao sphere
call mf8nao stars
call mf8nao trivial
call mf8nao modview
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
del *.obj *.mod *.exe *.f9 fppr.inc *.exp *.lib *.mds
rmdir /s/q album

:done
