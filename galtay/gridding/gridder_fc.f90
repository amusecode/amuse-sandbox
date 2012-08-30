! I'm going to always use the Python convention for arrays and 
! start their indexing at 0 as opposed to the default Fortran 1
!=================================================================

module fortran_code_module
implicit none



real(8), parameter :: PI = 3.1415926535897932d0
real(8), parameter :: SQRT_PI = 1.7724538509055159d0
real(8), parameter :: SEVEN_OVER_EIGHTY = 7.0d0 / 80
real(8), parameter :: FORTY_OVER_SEVEN_PI = 40.0d0 / ( 7 * PI )
real(8), parameter :: SQRT_PI_OVER_TWO = SQRT_PI / 2



contains



subroutine return_erf_integrals( lims, h, nn, out )
  real(8), dimension(0:nn-1), intent(in) :: lims  ! limits of integration
  real(8), intent(in) :: h                        ! smoothing length
  integer(4), intent(in) :: nn                    ! size of lims
  real(8), dimension(0:nn-2), intent(out) :: out  ! integrals between lims


  real(8) :: var   ! variance of equivalent gaussian
  real(8) :: sqv   ! sqrt(2*var)
  real(8) :: AA    ! 1/sqv

  var =  h * h * SEVEN_OVER_EIGHTY
  sqv = sqrt(2*var)
  AA = 1.0d0 / sqv

  out = (SQRT_PI_OVER_TWO * sqv) * &
       ( erf( lims(1:nn-1) * AA ) - erf( lims(0:nn-2) * AA ) )

end subroutine return_erf_integrals



subroutine prolong_grid( ngc, coarse, ngf, fine )

  integer(4), intent(in) :: ngc  ! ng for coarse grid
  real(8), dimension(0:ngc-1,0:ngc-1), intent(in) :: coarse
  integer(4), intent(in) :: ngf  ! ng for coarse grid
  real(8), dimension(0:ngf-1,0:ngf-1), intent(out) :: fine

  integer(4) :: ic, jc, ifine, jfine


  do jfine = 0, ngf-1
     jc = jfine / 2

     do ifine = 0, ngf-1
        ic = ifine / 2
     
        fine(ifine,jfine) = coarse(ic,jc) * 0.25d0
        
     end do

  end do


end subroutine prolong_grid



subroutine restrict_grid( ngf, fine, ngc, coarse )

  integer(4), intent(in) :: ngf  ! ng for coarse grid
  real(8), dimension(0:ngf-1,0:ngf-1), intent(in) :: fine
  integer(4), intent(in) :: ngc  ! ng for coarse grid
  real(8), dimension(0:ngc-1,0:ngc-1), intent(out) :: coarse

  integer(4) :: jlo, jhi, ilo, ihi, icoarse, jcoarse

  do jcoarse = 0, ngc-1
     jlo = 2 * jcoarse
     jhi = jlo + 1

     do icoarse = 0, ngc-1
        ilo = 2 * icoarse
        ihi = ilo + 1
     
        coarse(icoarse,jcoarse) = &
             fine(ilo,jlo) + fine(ilo,jhi) + &
             fine(ihi,jlo) + fine(ihi,jhi)
        
     end do

  end do

end subroutine restrict_grid







subroutine return_sph_grid( xp, yp, hsml, w, &
     grid_x_lo, grid_y_lo, grid_len, ng, periodic, &
     verbose, npar, grid )

  real(8), dimension(0:npar-1), intent(in) :: xp
  real(8), dimension(0:npar-1), intent(in) :: yp
  real(8), dimension(0:npar-1), intent(in) :: hsml
  real(8), dimension(0:npar-1), intent(in) :: w
  real(8), intent(in) :: grid_x_lo
  real(8), intent(in) :: grid_y_lo
  real(8), intent(in) :: grid_len
  integer(4), intent(in) :: ng
  integer(4), intent(in) :: periodic
  integer(4), intent(in) :: verbose
  integer(4), intent(in) :: npar
  real(8), dimension(0:ng-1,0:ng-1), intent(out) :: grid

  real(8) :: dl
  real(8) :: xp_g
  real(8) :: yp_g

  integer(4) :: xsq_lo_g 
  integer(4) :: xsq_hi_g   
  integer(4) :: ysq_lo_g 
  integer(4) :: ysq_hi_g 
  integer(4) :: xsq_dl
  integer(4) :: ysq_dl
  integer(4) :: prd_sq_dl

  integer(4) :: nxcells, nycells
  integer(4) :: nxedges, nyedges

  real(8), dimension(:), allocatable :: xlims
  real(8), dimension(:), allocatable :: ylims
  real(8), dimension(:), allocatable :: xout
  real(8), dimension(:), allocatable :: yout
  real(8), dimension(:,:), allocatable :: mat


  integer(4) :: ip, igx, igy, imx, imy, i



  ! linear extent of a grid cell
  !----------------------------------
  dl = grid_len / ng


  ! report
  !----------------------------------
  !write(*,*) 
  !write(*,*) ' fortran called '
  !write(*,*) ' grid_x_lo: ', grid_x_lo
  !write(*,*) ' grid_y_lo: ', grid_y_lo
  !write(*,*) ' grid_len:  ', grid_len
  !write(*,*) ' ng:        ', ng
  !write(*,*) ' periodic:  ', periodic
  !write(*,*) ' npar:      ', npar
  !write(*,*) 



  ! begin loop over particles
  !----------------------------------
  do ip = 0, npar-1

     ! report progress  
     !-----------------------------------
     if (verbose==1) then
        if (npar > 1000) then
           if ( mod(ip,npar/4)==0 ) then
              write(*,*) 'done ', ip, 'of ', npar
           endif
        endif
     endif

     ! position in grid units.  
     !-----------------------------------
     xp_g = ( xp(ip) - grid_x_lo ) / dl
     yp_g = ( yp(ip) - grid_y_lo ) / dl

     ! low/high corner of square in grid units
     !-----------------------------------
     xsq_lo_g = xp_g - hsml(ip) / dl
     xsq_hi_g = xp_g + hsml(ip) / dl

     ysq_lo_g = yp_g - hsml(ip) / dl
     ysq_hi_g = yp_g + hsml(ip) / dl

     ! cycle if no part of particle is on grid
     !-----------------------------------
     if ( xsq_hi_g < 0 ) then
        cycle
     endif

     if ( ysq_hi_g < 0 ) then
        cycle
     endif

     if ( xsq_lo_g >= ng ) then
        cycle
     endif

     if ( ysq_lo_g >= ng ) then
        cycle
     endif


     ! number of cells covered in each dimension
     !-----------------------------------
     xsq_dl = ( xsq_hi_g - xsq_lo_g ) + 1
     ysq_dl = ( ysq_hi_g - ysq_lo_g ) + 1

     ! total number of cells covered
     !-----------------------------------
     prd_sq_dl = xsq_dl * ysq_dl

     !
     ! [xy]sq_dl determines the extent of the particle
     ! in each dimension in units of grid cells
     !   [xy]sq_dl     = # of cells covered
     !   [xy]sq_dl + 1 = # of cell edges
     !

     nxcells = xsq_dl
     nycells = ysq_dl

     nxedges = nxcells + 1
     nyedges = nycells + 1

     ! if particle all in one cell
     !---------------------------------
     if ( prd_sq_dl == 1 ) then
        igx = xsq_lo_g
        igy = ysq_lo_g
        grid(igx,igy) = grid(igx,igy) + w(ip)

     ! if particle spread over multiple cells
     !---------------------------------
     else
        allocate( xlims( 0:nxedges-1 ) )  ! number of x-edges
        allocate( ylims( 0:nyedges-1 ) )  ! number of y-edges
        allocate( xout(  0:nxcells-1 ) )  ! number of x-cells
        allocate( yout(  0:nycells-1 ) )  ! number of y-cells

        allocate( mat( 0:nxcells-1, 0:nycells-1 ) )

        do i = 0, nxedges-1
           xlims(i) = ( xsq_lo_g + i ) * dl
        end do

        do i = 0, nyedges-1
           ylims(i) = ( ysq_lo_g + i ) * dl
        end do

        xlims = xlims + grid_x_lo - xp(ip)
        ylims = ylims + grid_y_lo - yp(ip)

        call return_erf_integrals( xlims, hsml(ip), size(xlims), xout )
        call return_erf_integrals( ylims, hsml(ip), size(ylims), yout )

        mat = spread(xout, dim=2, ncopies=size(yout)) * &
              spread(yout, dim=1, ncopies=size(xout))

        mat = mat / sum(mat)
        mat = mat * w(ip)

        ! outer loop over x-cells
        !----------------------------------------
        over_xcells: do imx = 0, nxcells-1

           ! increment forward from low x corner
           !----------------------------------------
           igx = xsq_lo_g + imx 

           ! if were out of x-bounds wrap or skip
           !----------------------------------------
           if (igx >= ng) then
              if (periodic==1) then                 
                 igx = igx-ng
                 if ( igx >= ng ) cycle ! if one wrap doesnt bring us in SKIP
              else
                 cycle
              endif
           endif
           if (igx < 0) then
              if (periodic==1) then
                 igx = igx+ng
                 if ( igx < 0 ) cycle ! if one wrap doesnt bring us in SKIP
              else
                 cycle
              endif
           endif

           ! inner loop over y-cells
           !----------------------------------------
           over_ycells: do imy = 0, nycells-1

              ! increment forward from low y corner
              !----------------------------------------
              igy = ysq_lo_g + imy 

              ! if were out of y-bounds wrap or skip
              !----------------------------------------
              if (igy >= ng) then
                 if (periodic==1) then
                    igy = igy-ng
                    if ( igy >= ng ) cycle ! if one wrap doesnt bring us in SKIP
                 else
                    cycle
                 endif
              endif
              if (igy < 0) then
                 if (periodic==1) then
                    igy = igy+ng
                    if ( igy < 0 ) cycle ! if one wrap doesnt bring us in SKIP
                 else
                    cycle
                 endif
              endif


              !---------------------------------------------
              grid(igx,igy) = grid(igx,igy) + mat(imx,imy)
              !---------------------------------------------


           end do over_ycells

        end do over_xcells

        deallocate( xlims, ylims, xout, yout, mat )

     endif

  end do


end subroutine return_sph_grid












end module fortran_code_module



!!$program test
!!$use fortran_code_module
!!$implicit none
!!$
!!$real(8) :: xlims(5)
!!$real(8) :: h
!!$real(8) :: out(4)
!!$
!!$
!!$real(8) :: veca(4)
!!$real(8) :: vecb(3)
!!$real(8) :: mat(4,3)
!!$
!!$call return_erf_integrals(xlims, h, size(xlims), out)
!!$call return_erf_integrals(xlims, h, size(xlims), out)
!!$
!!$mat = spread(veca, dim=2, ncopies=size(vecb)) * &
!!$      spread(vecb, dim=1, ncopies=size(veca))
!!$
!!$
!!$end program test
