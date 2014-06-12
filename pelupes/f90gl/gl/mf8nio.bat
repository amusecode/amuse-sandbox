REM Many of the macros are inherited from the makefile in the parent directory

REM select what to make; default is library

if "%1"=="library" goto library
if "%1"=="install" goto install
if "%1"=="clean" goto clean

REM build libraries

:library

REM make the include file of preprocessor directives

%F90C% isshort.f90 %F90FLAGS%
type %F90PPR_INC% > fppr.inc
isshort
type defshort >> fppr.inc

REM make the C wrappers

%CC% %CFLAGS% /DFNAME=%FNAME% /c cwrap.c

REM make the opengl_kinds module

%F90PPR% < glkinds.fpp > glkinds.f90
%F90C% -c glkinds.f90 %F90FLAGS% %USEMOD%

REM make the interface module, Fortran wrappers and module opengl_gl

%F90PPR% < interf.fpp > interf.f90
%F90C% -c interf.f90 %F90FLAGS% %USEMOD%
%F90PPR% < fwrap.fpp > fwrap.f90
%F90C% -c fwrap.f90 %F90FLAGS% %USEMOD%

REM package them into the libraries

lib /out:f90GL.lib cwrap.obj interf.obj fwrap.obj glkinds.obj

goto done

REM install the libraries

:install

REM first build them
call mf8nio library

REM then move them

move f90GL.lib ..\lib\
move opengl_fwrap.mod ..\lib\
move opengl_gl.mod ..\lib\
move opengl_glinterfaces.mod ..\lib\
move opengl_kinds.mod ..\lib\

goto done

REM delete everything created while building the libraries

:clean

del fppr.inc fwrap.f90 interf.f90 glkinds.f90 isshort.exe *.obj *.map *.mod *.lib *.lnk modtable.txt defshort
goto done

:done
