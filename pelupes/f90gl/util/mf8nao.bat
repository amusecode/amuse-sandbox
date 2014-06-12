REM Many of the macros are inherited from the batch file in the parent directory

REM select what to make

if "%1"=="sppr" goto sppr
if "%1"=="clean" goto clean
goto done

:sppr
%F90C% sppr.f90 %F90FLAGS%
goto done

:clean

del sppr.exe *.obj *.map *.mod *.lib modtable.txt
rmdir /s/q album
goto done

:done
