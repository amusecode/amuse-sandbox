@echo off
REM This batch file has been tested on
REM   computer: Pentium III
REM         OS: Windows 2000
REM   compiler: Absoft Pro Fortran 9.0
REM          C: Microsoft Visual C++ 6.0
REM     OpenGL: Microsoft OpenGL (Absoft import lib), Glut 3.7.1, f90gl 1.2.9

REM ------------- User configuration parameters ---------------

REM modify these for your system

REM any relative paths should be relative to a subdirectory of the directory
REM containing this makefile

REM procedure name mangling approach used by your fortran compiler
REM  LOWERCASE  - convert to lower case
REM  UNDERSCORE - convert to lower case and append underscore
REM  UPPERCASE  - convert to upper case
set FNAME=UPPERCASE

REM the include directory for Windows visual C
set WININC=c:\Program Files\Microsoft Visual Studio\vc98\include

REM fortran 90 compiler and compiler flags
set F90C=f90
set F90FLAGS=

REM C compiler and compiler flags
set CC=cl
set CFLAGS=-Zl /I"%WININC%" /DWIN32 /DABS32 /DVERSION_1_1
 
REM The compiler flag to get at the module file for opengl_kinds in ../include/GL
set USEMOD=-p ..\include\GL

REM ----------- end of user configuration parameters ------------

REM The preprocessor
set F90PPR=..\util\sppr

REM The include file for f90ppr directives
set F90PPR_INC=..\fppincs\fpp8nbo

REM select action; default is install

if "%1"=="install" goto install
if "%1"=="clean" goto clean
if "%1"=="reallyclean" goto reallyclean

REM build and install the libraries

:install
cd util
call mf8nao sppr
cd ..
cd gl
call mf8nbo install
cd ..
cd glu
call mf8nbo install
cd ..
cd glut
call mf8nbo install
cd ..
goto done

REM delete files created during compilation, but keep the libraries

:clean
cd gl
call mf8nao clean
cd ..
cd glu
call mf8nao clean
cd ..
cd glut
call mf8nao clean
cd ..
cd examples
call mf8nao clean
cd ..
cd util
call mf8nao clean
cd ..
goto done

REM delete everything, including libraries and utilities

:reallyclean
cd gl
call mf8nao clean
cd ..
cd glu
call mf8nao clean
cd ..
cd glut
call mf8nao clean
cd ..
cd examples
call mf8nao clean
cd ..
cd util
call mf8nao clean
cd ..
del lib\*.lib
del include\GL\*.mod
goto done

:done
