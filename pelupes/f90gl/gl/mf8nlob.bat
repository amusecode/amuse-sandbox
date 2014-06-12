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

%CC% %CFLAGS% -DFNAME=%FNAME% cwrap.c

REM make the opengl_kinds module

%F90PPR% < glkinds.fpp > glkinds.f90
%F90C% -c glkinds.f90 %F90FLAGS% %USEMOD%

REM make the interface module, Fortran wrappers and module opengl_gl

%F90PPR% < interf.fpp > interf.f90
%F90C% -c interf.f90 %F90FLAGS% %USEMOD%
%F90PPR% < fwrap.fpp > fwrap.f90
%F90C% -c fwrap.f90 %F90FLAGS% %USEMOD%

REM package them into the libraries

ren cwrap.obj cwrap.tmp
del *.obj
move ..\lib\opengl_*.lib .
lm opengl_0.lib/extractall
ren temp*.obj ogl0*.obj
lm opengl_2.lib/extractall
ren temp*.obj ogl2*.obj
ren cwrap.tmp cwrap.obj
echo f90GL.lib > response.tmp
if a%OS% == aWindows_NT goto WNT
rem Windows 9x
for %%f in (*.obj) do echo -+%%f & >> response.tmp
goto NEXT
:WNT
rem Windows NT
for %%f in (*.obj) do echo -+%%f ^& >> response.tmp
:NEXT
echo ; >> response.tmp
lm @response.tmp
del response.tmp
del *.obj
move opengl_*.lib ..\lib

goto done

REM install the libraries

:install

REM first build them
call mf8nlob library

REM then move them

move *.lib ..\lib
goto done

REM delete everything created while building the libraries

:clean

del fppr.inc
del glkinds.f90 fwrap.f90 interf.f90
del isshort.exe
del *.obj
del *.map
del *.mod
del *.lib
del *.lnk
del modtable.txt
del defshort
goto done

:done
