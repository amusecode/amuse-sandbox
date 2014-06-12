@echo off
REM This batch file has been tested on
REM   computer: Pentium
REM         OS: Windows 95 OSR 2
REM   compiler: Lahey F90 4.50
REM          C: Borland C++ 4.5
REM     OpenGL: Microsoft OpenGL, Glut 3.6 with patch, f90gl 1.1.4
REM updated by: Lawson B. Wakefield (I.S.S. Ltd.)

REM ------------- User configuration parameters ---------------

REM modify these for your system

REM any relative paths should be relative to a subdirectory of the directory
REM containing this makefile

REM procedure name mangling approach used by your fortran compiler
REM  LOWERCASE  - convert to lower case
REM  UNDERSCORE - convert to lower case and append underscore
REM  UPPERCASE  - convert to upper case
set FNAME=UPPERCASE

REM the include directory for the C compiler.  The OpenGL and GLUT libraries
REM should be there, too.
set WININC=

REM fortran 90 compiler and compiler flags
set F90C=lf90
set F90FLAGS=-nomap

REM C compiler and compiler flags
set CC=bcc32
set CFLAGS=+..\bcclf90.cfg
 
REM The compiler flag to get at the module file for opengl_kinds in ../include/GL
set USEMOD=-mod ..\lib

REM ----------- end of user configuration parameters ------------

REM The preprocessor
set F90PPR=..\util\sppr

REM The include file for f90ppr directives
set F90PPR_INC=..\fppincs\fpp8nlob

REM select action; default is install

if "%1"=="install" goto install
if "%1"=="clean" goto clean
if "%1"=="reallyclean" goto reallyclean

REM build and install the libraries

:install
cd util
call mf8nlob sppr
cd ..
cd gl
call mf8nlob install
cd ..
cd glu
call mf8nlob install
cd ..
cd glut
call mf8nlob install
cd ..
goto done

REM delete files created during compilation, but keep the libraries

:clean
cd gl
call mf8nlob clean
cd ..
cd glu
call mf8nlob clean
cd ..
cd glut
call mf8nlob clean
cd ..
cd examples
call mf8nlob clean
cd ..
cd util
call mf8nlob clean
cd ..
goto done

REM delete everything, including libraries and utilities

:reallyclean
cd gl
call mf8nlob clean
cd ..
cd glu
call mf8nlob clean
cd ..
cd glut
call mf8nlob clean
cd ..
cd examples
call mf8nlob clean
cd ..
cd util
call mf8nlob clean
cd ..
del lib\*.lib
del lib\*.mod
del lib\modtable.txt
goto done

:done
