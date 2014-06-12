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

REM any relative paths should be relative to a subdirectory of the directory
REM containing this makefile

REM procedure name mangling approach used by your fortran compiler
REM  LOWERCASE  - convert to lower case
REM  UNDERSCORE - convert to lower case and append underscore
REM  UPPERCASE  - convert to upper case
set FNAME=LOWERCASE

REM the include directory for the C compiler.  The OpenGL and GLUT libraries
REM should be there, too.
REM set WININC=%SCCINCLUDE%
set WININC=c:\Program Files\Microsoft Visual Studio\vc98\include

REM fortran 90 compiler and compiler flags
set F90C=ftn95
set F90FLAGS= /link

REM C compiler and compiler flags
REM set CC=scc
REM set CFLAGS= /ansi_c
set CC=cl
set CFLAGS=-Zl /I"%WININC%" /DWIN32
 
REM The compiler flag to get at the module file for opengl_kinds in ../include/GL
set USEMOD=/MOD_PATH ..\lib

REM ----------- end of user configuration parameters ------------

REM The preprocessor
set F90PPR=..\util\sppr

REM The include file for f90ppr directives
set F90PPR_INC=..\fppincs\fpp8noo

REM select action; default is install

if "%1"=="install" goto install
if "%1"=="clean" goto clean
if "%1"=="reallyclean" goto reallyclean

REM build and install the libraries

:install
cd util
call mf8noo sppr
cd ..
cd gl
call mf8noo install
cd ..
cd glu
call mf8noo install
cd ..
cd glut
call mf8noo install
cd ..
goto done

REM delete files created during compilation, but keep the libraries

:clean
cd gl
call mf8noo clean
cd ..
cd glu
call mf8noo clean
cd ..
cd glut
call mf8noo clean
cd ..
cd examples
call mf8noo clean
cd ..
cd util
call mf8noo clean
cd ..
goto done

REM delete everything, including libraries and utilities

:reallyclean
cd gl
call mf8noo clean
cd ..
cd glu
call mf8noo clean
cd ..
cd glut
call mf8noo clean
cd ..
cd examples
call mf8noo clean
cd ..
cd util
call mf8noo clean
cd ..
del lib\*.lib
del lib\*.mod
del lib\modtable.txt
goto done

:done
