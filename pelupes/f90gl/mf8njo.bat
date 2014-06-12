@echo off
REM This batch file has been tested on
REM   computer: Pentium II
REM         OS: Windows NT 4.0 (service pack 4)
REM   compiler: Lahey F95 5.00b
REM          C: Microsoft Visual C++ 6.0
REM     OpenGL: Microsoft OpenGL (Lahey import lib), Glut 3.7.1, f90gl 1.2.0

REM ------------- User configuration parameters ---------------

REM modify these for your system

REM any relative paths should be relative to a subdirectory of the directory
REM containing this makefile

REM procedure name mangling approach used by your fortran compiler
REM  LOWERCASE  - convert to lower case
REM  UNDERSCORE - convert to lower case and append underscore
REM  UPPERCASE  - convert to upper case
set FNAME=UNDERSCORE

REM the include directory for the C compiler.  The OpenGL and GLUT libraries
REM should be there, too.
set WININC=c:\Program Files\Microsoft Visual Studio\VC98\Include

REM fortran 90 compiler and compiler flags
set F90C=lf95
set F90FLAGS=-nomap

REM C compiler and compiler flags
set CC=cl
set CFLAGS=-DWIN32 -DLF95 -I"%WININC%"
 
REM The compiler flag to get at the module file for opengl_kinds in ../include/GL
set USEMOD=-mod ..\lib

REM ----------- end of user configuration parameters ------------

REM The preprocessor
set F90PPR=..\util\sppr

REM The include file for f90ppr directives
set F90PPR_INC=..\fppincs\fpp8njo

REM select action; default is install

if "%1"=="install" goto install
if "%1"=="clean" goto clean
if "%1"=="reallyclean" goto reallyclean

REM build and install the libraries

:install
cd util
call mf8njo sppr
cd ..
cd gl
call mf8njo install
cd ..
cd glu
call mf8njo install
cd ..
cd glut
call mf8njo install
cd ..
goto done

REM delete files created during compilation, but keep the libraries

:clean
cd gl
call mf8njo clean
cd ..
cd glu
call mf8njo clean
cd ..
cd glut
call mf8njo clean
cd ..
cd examples
call mf8njo clean
cd ..
cd util
call mf8njo clean
cd ..
goto done

REM delete everything, including libraries and utilities

:reallyclean
cd gl
call mf8njo clean
cd ..
cd glu
call mf8njo clean
cd ..
cd glut
call mf8njo clean
cd ..
cd examples
call mf8njo clean
cd ..
cd util
call mf8njo clean
cd ..
del lib\*.lib
del lib\*.mod
goto done

:done
