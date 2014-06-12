REM Many of the macros are inherited from the makefile in the parent directory

REM select what to make; default is library

if "%1"=="library" goto library
if "%1"=="install" goto install
if "%1"=="clean" goto clean

REM build libraries

:library

REM get the include file of preprocessor directives

type %F90PPR_INC% > fppr.inc
..\gl\isshort
type defshort >> fppr.inc

REM make the C wrappers

%CC% %CFLAGS% -DFNAME=%FNAME% cwrapglu.c

REM make the interface module, Fortran wrappers and module opengl_gl

%F90PPR% < intrfglu.fpp > intrfglu.f90
%F90C% -c intrfglu.f90 %F90FLAGS% %USEMOD%
%F90PPR% < fwrapglu.fpp > fwrapglu.f90
%F90C% -c fwrapglu.f90 %F90FLAGS% %USEMOD%

REM package them into the libraries

lib /out:f90GLU.lib *.obj
del *.obj

goto done

REM install the libraries

:install

REM first build them
call mf8njob library

REM then move them

move *.lib ..\lib
goto done

REM delete everything created while building the libraries

:clean

del fppr.inc
del intrfglu.f90 fwrapglu.f90
del *.obj
del *.mod
del *.lib
del *.lnk
del defshort
goto done

:done
