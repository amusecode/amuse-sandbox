REM Many of the macros are inherited from the batch file in the parent directory

REM select what to make

if "%1"=="sppr" goto sppr
if "%1"=="clean" goto clean
goto done

:sppr
%F90C% sppr.f90 %F90FLAGS%
cfig386 sppr.exe -nosignon
goto done

:clean

del sppr.exe
del *.obj
del *.mod
del *.lib
del modtable.txt
goto done

:done
