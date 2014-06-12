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

%F90PPR% < intrfglt.fpp > intrfglt.f90
%F90C% -c intrfglt.f90 %F90FLAGS% %USEMOD%
%F90PPR% < fwrapglt.fpp > fwrapglt.f90
%F90C% -c fwrapglt.f90 %F90FLAGS% %USEMOD%

REM package them into the libraries

lib /out:f90GLUT.lib cwrapglt.obj intrfglt.obj fwrapglt.obj

goto done

REM install the libraries

:install

REM first build them
call mf8nio library

REM then move them

move f90GLUT.lib ..\lib\
move opengl_glut.mod ..\lib\
move opengl_glutinterfaces.mod ..\lib\
goto done

REM delete everything created while building the libraries

:clean

del fppr.inc fwrapglt.f90 intrfglt.f90 *.obj *.map *.mod *.lib *.lnk modtable.txt
goto done

:done
