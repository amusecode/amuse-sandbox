!----------------------------------------------------------------------------
!
!                    f90gl version 1.2.12, September 2006
!
!     William F. Mitchell
!     Mathematical and Computational Sciences Division
!     National Institute of Standards and Technology
!     william.mitchell@nist.gov
!     http://math.nist.gov/f90gl
!
! f90gl is public domain software.  It was produced as part of work done by
! the U.S. Government, and is not subject to copyright in the United States.
!
!----------------------------------------------------------------------------


#include 'fppr.inc'

module opengl_kinds
implicit none

public :: operator(==)
private :: ptrcompare

! Define the types

! default kind of logical; probably not a one byte logical
integer, parameter, public :: glboolean = kind(.true.)
logical, parameter, public :: glone_byte_logical = .false.

! The remaining kind numbers are set automatically, but may need to be
! modified on some systems to match the C types.

integer, parameter, public :: &
#ifdef SUNSOFT
! The SunSoft compiler before version 2.0 on Solaris uses 4 bytes for the short
! integers even though they have a different kind parameter.  Don't use them.
   glbyte     = selected_int_kind(9), &
   glubyte    = selected_int_kind(9), &
   glshort    = selected_int_kind(9), &
   glushort   = selected_int_kind(9), &
#else
   glbyte     = selected_int_kind(2), &
   glubyte    = selected_int_kind(2), &
   glshort    = selected_int_kind(4), &
   glushort   = selected_int_kind(4), &
#endif
   glint      = selected_int_kind(9), &
   gluint     = selected_int_kind(9), &
   glenum     = selected_int_kind(9), &
   glbitfield = selected_int_kind(9), &
   glsizei    = selected_int_kind(9), &
   glfloat    = selected_real_kind(6), &
   glclampf   = selected_real_kind(6), &
   gldouble   = selected_real_kind(15), &
   glclampd   = selected_real_kind(15), &
   glcint     = selected_int_kind(9), &
   glcuint    = selected_int_kind(9), &
   glcfloat   = selected_real_kind(6)

! Determine if the short types exist

logical, parameter, public :: &
   glshort_is_int = (glshort == glint), &
   glbyte_is_int = (glbyte == glint), &
   glbyte_is_short = (glbyte == glshort .and. .not. glbyte_is_int)

! C pointers
#ifdef PTR8
integer(kind=glint), parameter, public :: glcptr_length = 8
#else
integer(kind=glint), parameter, public :: glcptr_length = 4
#endif
type, public :: glcptr
   character(len=glcptr_length) :: addr
end type glcptr

type (glcptr), parameter, public :: glnullptr = glcptr( &
#ifdef PTR8
   char(0)//char(0)//char(0)//char(0)// &
#endif
   char(0)//char(0)//char(0)//char(0))

interface operator(==)
   module procedure ptrcompare
end interface

contains

! comparison operator for C pointers

function ptrcompare(p1,p2) result(resptrcompare)
type (glcptr), intent(in) :: p1, p2
logical(kind=glboolean) :: resptrcompare
#ifdef NAGWARE21
! An apparent bug in NAGWare f90 version 2.1 requires that these strings
! be compared one character at a time
integer :: i
resptrcompare = .true.
do i=1,glcptr_length
   if (p1%addr(i:i) /= p2%addr(i:i)) then
      resptrcompare = .false.
   end if
end do
#else
resptrcompare = (p1%addr == p2%addr)
#endif
return
end function ptrcompare

end module opengl_kinds
