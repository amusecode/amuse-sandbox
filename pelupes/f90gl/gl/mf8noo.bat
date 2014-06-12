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

%CC% %CFLAGS%  /DFNAME=%FNAME% /c cwrap.c

REM change the Fortran compiler flag for DLL name mangling

@rem set F90FLAGS=-win -ml msvc -dll
set F90FLAGS=/FPS_EXTERNAL

REM make the opengl_kinds module

%F90PPR% < glkinds.fpp > glkinds.f90
%F90C%  glkinds.f90 %F90FLAGS% %USEMOD%

REM make the interface module, Fortran wrappers and module opengl_gl

%F90PPR% < interf.fpp > interf.f90
%F90C%  interf.f90 %F90FLAGS% %USEMOD%
%F90PPR% < fwrap.fpp > fwrap.f90
%F90C%  fwrap.f90 %F90FLAGS% %USEMOD%

REM package them into the libraries

slink /out:f90GL.lib -addobj:cwrap.obj  -addobj:fwrap.obj -addobj:glkinds.obj
goto done                                     

REM install the libraries

:install

REM first build them
call mf8noo library

REM then move them

copy *.lib ..\lib
move *.mod ..\lib
goto done

REM delete everything created while building the libraries

:clean

del fppr.inc fwrap.f90 interf.f90 glkinds.f90 isshort.exe *.obj *.map *.mod *.lib *.lnk modtable.txt defshort
goto done

:done  
