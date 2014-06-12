REM Many of the macros are inherited from the makefile in the parent directory

REM select what to make; default is library

if "%1"=="library" goto library
if "%1"=="install" goto install
if "%1"=="clean" goto clean

REM build libraries

:library

REM get the include file of preprocessor directives

copy %F90PPR_INC% fppr.inc

REM make the C wrappers

%CC% %CFLAGS% /DFNAME=%FNAME% /c cwrapglt.c

REM make the interface module, Fortran wrappers and module opengl_gl

f95 -use f90gl
f95 -use f90glu

%F90PPR% < intrfglt.fpp > intrfglt.f90
%F90C% -c intrfglt.f90 %F90FLAGS% %USEMOD%
%F90PPR% < fwrapglt.fpp > fwrapglt.f90
%F90C% -c fwrapglt.f90 %F90FLAGS% %USEMOD%

REM package them into the libraries

f95 -package f90glut intrfglt.obj fwrapglt.obj -mslink glut32.lib cwrapglt.obj f90gl.lib f90glu.lib glu32.lib msvcrt.lib

goto done

REM install the libraries

:install

REM first build them
call mf8nao library

REM then move them

if exist ..\lib\f90glut.lib del ..\lib\f90glut.lib
if exist ..\lib\f90glut.dll del ..\lib\f90glut.dll
if exist ..\lib\f90glut rmdir /s/q ..\lib\f90glut
move f90glut.lib ..\lib\
move f90glut.dll ..\lib\
move f90glut ..\lib\
goto done

REM delete everything created while building the libraries

:clean

del fppr.inc fwrapglt.f90 intrfglt.f90 *.obj *.map *.mod *.lib *.lnk modtable.txt *.exp
rmdir /s/q album
goto done

:done
