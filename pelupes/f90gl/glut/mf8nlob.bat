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

%CC% %CFLAGS% -DFNAME=%FNAME% cwrapglt.c

REM make the interface module, Fortran wrappers and module opengl_gl

%F90PPR% < intrfglt.fpp > intrfglt.f90
%F90C% -c intrfglt.f90 %F90FLAGS% %USEMOD%
%F90PPR% < fwrapglt.fpp > fwrapglt.f90
%F90C% -c fwrapglt.f90 %F90FLAGS% %USEMOD%

REM package them into the libraries

ren cwrapglt.obj cwrapglt.tmp
del *.obj
move ..\lib\opengl_*.lib .
lm opengl_7.lib/extractall
ren temp*.obj ogl7*.obj
ren cwrapglt.tmp cwrapglt.obj
echo f90GLUT.lib > response.tmp
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
del intrfglt.f90 fwrapglt.f90
del *.obj
del *.map
del *.mod
del *.lib
del *.lnk
del modtable.txt
goto done

:done
