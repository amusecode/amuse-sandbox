@echo off
REM This batch file has been tested on
REM   computer: Intel Xeon-based system
REM         OS: Windows Server 2003 Enterprise x64 Edition
REM   compiler: Intel Visual Fortran 9.1
REM          C: Intel C++ Compiler 9.1
REM     OpenGL: Microsoft OpenGL, Glut 3.7.6, f90gl 1.2.12

REM ------------- User configuration parameters ---------------

REM modify these for your system

REM any relative paths should be relative to a subdirectory of the directory
REM containing this makefile

REM procedure name mangling approach used by your fortran compiler
REM  LOWERCASE  - convert to lower case
REM  UNDERSCORE - convert to lower case and append underscore
REM  UPPERCASE  - convert to upper case
REM  MODULE     - prepend module_name_MOD_ (used for faking out F)
set FNAME=UPPERCASE

REM the include directory for Windows visual C
set WININC=F:\Program Files\Microsoft Platform SDK for Windows Server 2003 R2\Include

REM the include directory for GLUT
set GLUTINC=C:\glut-3.7.1-Intel64\include

REM fortran 90 compiler and compiler flags
set F90C=ifort
set F90FLAGS=/nologo

REM C compiler and compiler flags
set CC=icl
set CFLAGS=-DWIN32 -DIFORT /I"%WININC%" /I"%GLUTINC%" /DOPENGL_1_1 /nologo
 
REM The compiler flag to get at the module file for opengl_kinds in ../include/GL
set USEMOD=/I"..\lib"

REM ----------- end of user configuration parameters ------------

REM The preprocessor
set F90PPR=..\util\sppr

REM The include file for f90ppr directives
set F90PPR_INC=..\fppincs\fpp8nio

REM select action; default is install

if "%1"=="install" goto install
if "%1"=="clean" goto clean
if "%1"=="reallyclean" goto reallyclean

REM build and install the libraries

:install
cd util
call mfxnio sppr
cd ..
cd gl
call mfxnio install
cd ..
cd glu
call mfxnio install
cd ..
cd glut
call mfxnio install
cd ..
goto done

REM delete files created during compilation, but keep the libraries

:clean
cd gl
call mfxnio clean
cd ..
cd glu
call mfxnio clean
cd ..
cd glut
call mfxnio clean
cd ..
cd examples
call mfxnio clean
cd ..
cd util
call mfxnio clean
cd ..
goto done

REM delete everything, including libraries and utilities

:reallyclean
cd gl
call mfxnio clean
cd ..
cd glu
call mfxnio clean
cd ..
cd glut
call mfxnio clean
cd ..
cd examples
call mfxnio clean
cd ..
cd util
call mfxnio clean
cd ..
del lib\*.lib
del lib\*.mod
del lib\modtable.txt
goto done

:done
