REM Many of the macros are inherited from the makefile in the parent directory

REM select what to make; default is library

if "%1"=="library" goto library
if "%1"=="install" goto install
if "%1"=="clean" goto clean

REM build libraries

:library

REM get the include file of preprocessor directives

copy %F90PPR_INC% fppr.inc
..\gl\isshort
type defshort >> fppr.inc

REM make the C wrappers

%CC% %CFLAGS% /DFNAME=%FNAME% /c cwrapglu.c

REM make the interface module, Fortran wrappers and module opengl_gl

f95 -use f90gl
%F90PPR% < intrfglu.fpp > intrfglu.f90
%F90C% -c intrfglu.f90 %F90FLAGS% %USEMOD%
%F90PPR% < fwrapglu.fpp > fwrapglu.f90
%F90C% -c fwrapglu.f90 %F90FLAGS% %USEMOD%

REM package them into the libraries

f95 -package f90glu intrfglu.obj fwrapglu.obj -mslink cwrapglu.obj f90gl.lib glu32.lib msvcrt.lib
goto done

REM install the libraries

:install

REM first build them
call mf8nao library

REM then move them

if exist ..\lib\f90glu.lib del ..\lib\f90glu.lib
if exist ..\lib\f90glu.dll del ..\lib\f90glu.dll
if exist ..\lib\f90glu rmdir /s/q ..\lib\f90glu
move f90glu.lib ..\lib\
move f90glu.dll ..\lib\
move f90glu ..\lib\

goto done

REM delete everything created while building the libraries

:clean

del fppr.inc fwrapglu.f90 intrfglu.f90 *.obj *.map *.mod *.lib *.lnk modtable.txt defshort *.exp
rmdir /s/q album
goto done

:done
