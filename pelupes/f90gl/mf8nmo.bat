@echo off
REM This batch file has been tested on
REM   computer: Pentium 120MHz
REM         OS: Windows NT 4.0 (stock, no patch)
REM   compiler: Microsoft Fortran Powerstation 4.0
REM          C: Microsoft Visual C++ 4.0
REM     OpenGL: Microsoft OpenGL, Glut 3.6 with patch, f90gl 1.1.?
REM contributed by: Nicholas Yue yue@cwr.uwa.edu.au

REM ------------- User configuration parameters ---------------

REM modify these for your system

REM any relative paths should be relative to a subdirectory of the directory
REM containing this makefile

REM procedure name mangling approach used by your fortran compiler
REM  LOWERCASE  - convert to lower case
REM  UNDERSCORE - convert to lower case and append underscore
REM  UPPERCASE  - convert to upper case
set FNAME=LOWERCASE

REM the include directory for Windows
set WININC=c:\MSDEV\INCLUDE

REM fortran 90 compiler and compiler flags
set F90C=fl32
set F90FLAGS=

REM C compiler and compiler flags
set CC=cl
set CFLAGS=
 
REM The compiler flag to get at the module file for opengl_kinds in ../include/GL
set USEMOD=/I"..\lib"

REM ----------- end of user configuration parameters ------------

REM The preprocessor
set F90PPR=..\util\sppr

REM The include file for f90ppr directives
set F90PPR_INC=..\fppincs\fpp8nmo

REM select action; default is install

if "%1"=="install" goto install
if "%1"=="clean" goto clean
if "%1"=="reallyclean" goto reallyclean

REM build and install the libraries

:install
cd util
call mf8nmo sppr
cd ..
cd gl
call mf8nmo install
cd ..
cd glu
call mf8nmo install
cd ..
cd glut
call mf8nmo install
cd ..
goto done

REM delete files created during compilation, but keep the libraries

:clean
cd gl
call mf8nmo clean
cd ..
cd glu
call mf8nmo clean
cd ..
cd glut
call mf8nmo clean
cd ..
cd examples
call mf8nmo clean
cd ..
cd util
call mf8nmo clean
cd ..
goto done

REM delete everything, including libraries and utilities

:reallyclean
cd gl
call mf8nmo clean
cd ..
cd glu
call mf8nmo clean
cd ..
cd glut
call mf8nmo clean
cd ..
cd examples
call mf8nmo clean
cd ..
cd util
call mf8nmo clean
cd ..
del lib\*.lib
del lib\*.mod
del lib\modtable.txt
goto done

:done
