program isshort

! This program checks to see if one byte and two byte integers exist in the
! Fortran compiler that compiles it

implicit none
integer :: bytekind, shortkind, intkind

bytekind = selected_int_kind(2)
shortkind = selected_int_kind(4)
intkind = selected_int_kind(9)

open(unit=10,action="write",file="defshort",status="replace")

if (bytekind==shortkind .or. bytekind==intkind) then
   write(unit=10,fmt="(a)") "#define NOBYTE"
end if
if (shortkind==intkind) then
   write(unit=10,fmt="(a)") "#define NOSHORT"
end if

close(unit=10)

stop
end program isshort
