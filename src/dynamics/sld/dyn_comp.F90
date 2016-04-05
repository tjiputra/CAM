module dyn_comp
!----------------------------------------------------------------------- 
! 
! SLD dycore interface module
!
!-----------------------------------------------------------------------

use shr_kind_mod,    only: r8 => shr_kind_r8
use spmd_utils,      only: masterproc
use ppgrid,          only: pver, pverp
use pmgrid,          only: plev, plevp
use constituents,    only: pcnst, cnst_name, cnst_longname
use constituents,    only: tendnam, fixcnam, tottnam, hadvnam, vadvnam
use hycoef,          only: hycoef_init
use sld_control_mod, only: dif2, dif4, divdampn, eps, kmxhdc
use cam_history,     only: addfld, add_default, horiz_only
use pio,             only: file_desc_t
use cam_abortutils,  only: endrun
use cam_logfile,     only: iulog

#if (defined SPMD)
use spmd_dyn,        only: spmd_readnl, spmdinit_dyn
#endif

implicit none
private
save

public :: dyn_readnl, dyn_init, dyn_import_t, dyn_export_t

! these structures are not used in this dycore, but are included
! for source code compatibility.  
type dyn_import_t
   integer :: placeholder
end type dyn_import_t

type dyn_export_t
   integer :: placeholder
end type dyn_export_t

!=============================================================================================
contains
!=============================================================================================

subroutine dyn_readnl(nlfile)

   ! Read SLD namelist group.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_real8

   ! args
   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    
   ! local vars
   integer :: unitn, ierr

   real(r8) :: sld_dif2_coef          ! del2 horizontal diffusion coeff.
   real(r8) :: sld_dif4_coef          ! del4 horizontal diffusion coeff.
   real(r8) :: sld_divdampn      ! Number of days to invoke divergence damper
   real(r8) :: sld_tfilt_eps           ! time filter coefficient. Defaults to 0.06.
   integer  :: sld_kmxhdc        ! Number of levels to apply Courant limiter
    
   namelist /dyn_sld_inparm/ sld_dif2_coef, sld_dif4_coef, &
      sld_divdampn, sld_tfilt_eps, sld_kmxhdc

   character(len=*), parameter :: sub = 'dyn_readnl'
   !--------------------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'dyn_sld_inparm', status=ierr)
      if (ierr == 0) then
         read(unitn, dyn_sld_inparm, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub//': ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(sld_dif2_coef, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: sld_dif2_coef")

   call mpi_bcast(sld_dif4_coef, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: sld_dif4_coef")

   call mpi_bcast(sld_divdampn, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: sld_divdampn")

   call mpi_bcast(sld_tfilt_eps, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: sld_tfilt_eps")

   call mpi_bcast(sld_kmxhdc, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: sld_kmxhdc")

   dif2     = sld_dif2_coef
   dif4     = sld_dif4_coef
   divdampn = sld_divdampn
   eps      = sld_tfilt_eps
   kmxhdc   = sld_kmxhdc

   ! Write namelist variables to logfile
   if (masterproc) then

      write(iulog,*) 'SLD Dycore Parameters:'

      if (divdampn > 0._r8) then
         write(iulog,*) '  Divergence damper for spectral dycore invoked for days 0. to ',divdampn,' of this case'
      elseif (divdampn < 0._r8) then
         call endrun ('READ_NAMELIST: divdampn must be a positive number')
      else
         write(iulog,*) '  Divergence damper for spectral dycore NOT invoked'
      endif

      if (kmxhdc >= plev .or. kmxhdc < 0) then
         call endrun ('DYN_SLD_READNL:  ERROR:  KMXHDC must be between 0 and plev-1')
      end if

      write(iulog,9108) eps, dif2, dif4, kmxhdc
   end if

#if (defined SPMD)
   call spmd_readnl(nlfile)
#endif 

9108 format('   Time filter coefficient (EPS)                 ',f10.3,/,&
            '   DEL2 Horizontal diffusion coefficient (DIF2)  ',e10.3/, &
            '   DEL4 Horizontal diffusion coefficient (DIF4)  ',e10.3/, &
            '   Number of levels Courant limiter applied      ',i10)

end subroutine dyn_readnl

!=============================================================================================

subroutine dyn_init(file)

   use phys_control,        only: phys_getopts
   use dyn_grid,            only: define_cam_grids

   ! ARGUMENTS:
   type(file_desc_t), intent(in) :: file       ! PIO file handle for initial or restart file

   logical :: history_amwg       ! output for AMWG diagnostics
   ! Local workspace
   integer m

   call trunc()

#if (defined SPMD)
   call spmdinit_dyn()
#endif 

   ! Initialize hybrid coordinate arrays
   call hycoef_init(file)

   ! Define the CAM grids (must be before addfld calls)
   call define_cam_grids()

   call addfld ('ETADOT',(/ 'ilev' /),'A','1/s','Vertical (eta) velocity',gridname='gauss_grid')
   call addfld ('U&IC',  (/ 'lev' /), 'I','m/s','Zonal wind',             gridname='gauss_grid' )
   call addfld ('V&IC',  (/ 'lev' /), 'I','m/s','Meridional wind',        gridname='gauss_grid' )
   call add_default ('U&IC',0, 'I')
   call add_default ('V&IC',0, 'I')

   call addfld ('PS&IC',horiz_only, 'I','Pa','Surface pressure',          gridname='gauss_grid' )
   call addfld ('T&IC', (/ 'lev' /),'I','K', 'Temperature',               gridname='gauss_grid' )
   call add_default ('PS&IC',0, 'I')
   call add_default ('T&IC',0, 'I')

   do m = 1, pcnst
      call addfld (trim(cnst_name(m))//'&IC',(/ 'lev' /),'I', 'kg/kg',cnst_longname(m),                      gridname='gauss_grid')
      call add_default(trim(cnst_name(m))//'&IC',0, 'I')
      call addfld (hadvnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' horizontal advection tendency',gridname='gauss_grid')
      call addfld (vadvnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' vertical advection tendency',  gridname='gauss_grid')
      call addfld (tendnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' total tendency',               gridname='gauss_grid')
      call addfld (tottnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' horz + vert + fixer tendency', gridname='gauss_grid')
      call addfld (fixcnam(m),(/ 'lev' /),'A','kg/kg/s',trim(cnst_name(m))//' tendency due to slt fixer',    gridname='gauss_grid')
   end do

   call addfld ('DUH     ',(/ 'lev' /),'A','K/s     ','U horizontal diffusive heating',              gridname='gauss_grid')
   call addfld ('DVH     ',(/ 'lev' /),'A','K/s     ','V horizontal diffusive heating',              gridname='gauss_grid')
   call addfld ('DTH     ',(/ 'lev' /),'A','K/s     ','T horizontal diffusive heating',              gridname='gauss_grid')
   call addfld ('ENGYCORR',(/ 'lev' /),'A','W/m2    ','Energy correction for over-all conservation', gridname='gauss_grid')
   call addfld ('TFIX    ',horiz_only, 'A','K/s     ','T fixer (T equivalent of Energy correction)', gridname='gauss_grid')

   call phys_getopts(history_amwg_out=history_amwg)

   if (history_amwg) call add_default ('DTH     ', 1, ' ')

   call addfld ('FU      ',(/ 'lev' /),'A','m/s2    ','Zonal wind forcing term',            gridname='gauss_grid')
   call addfld ('FV      ',(/ 'lev' /),'A','m/s2    ','Meridional wind forcing term',       gridname='gauss_grid')
   call addfld ('UTEND   ',(/ 'lev' /),'A','m/s2    ','U tendency',                         gridname='gauss_grid')
   call addfld ('VTEND   ',(/ 'lev' /),'A','m/s2    ','V tendency',                         gridname='gauss_grid')
   call addfld ('TTEND   ',(/ 'lev' /),'A','K/s     ','T tendency',                         gridname='gauss_grid')
   call addfld ('LPSTEN  ',horiz_only, 'A','Pa/s    ','Surface pressure tendency',          gridname='gauss_grid')
   call addfld ('VAT     ',(/ 'lev' /),'A','K/s     ','Vertical advective tendency of T',   gridname='gauss_grid')
   call addfld ('KTOOP   ',(/ 'lev' /),'A','K/s     ','(Kappa*T)*(omega/P)',                gridname='gauss_grid')

end subroutine dyn_init

!#######################################################################
! Private routines
!#######################################################################

subroutine trunc()

! Check consistency of truncation parameters and evaluate pointers
! and displacements for spectral arrays
!
! Original version:  CCM1

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use pspect,         only: ptrm, ptrn, ptrk, pmax, pmmax
  use comspe,         only: ncoefi, nalp, nco2, nm, nstart, nlen, ncutoff, locm
  use cam_abortutils, only: endrun

  implicit none

!---------------------------Local variables-----------------------------
!
  integer n              ! Loop index over diagonals
  integer ik2            ! K+2
  integer m              ! loop index
!
!-----------------------------------------------------------------------
!
! trunc first evaluates truncation parameters for a general pentagonal 
! truncation for which the following parameter relationships are true
!
! 0 .le. |m| .le. ptrm
!
! |m| .le. n .le. |m|+ptrn for |m| .le. ptrk-ptrn
!
! |m| .le. n .le. ptrk     for (ptrk-ptrn) .le. |m| .le. ptrm
!
! Most commonly utilized truncations include:
!  1: triangular  truncation for which ptrk=ptrm=ptrn
!  2: rhomboidal  truncation for which ptrk=ptrm+ptrn
!  3: trapezoidal truncation for which ptrn=ptrk .gt. ptrm
!
! Simple sanity check
! It is necessary that ptrm .ge. ptrk-ptrn .ge. 0
!
  if (ptrm.lt.(ptrk-ptrn)) then
     call endrun ('TRUNC: Error in truncation parameters. ntrm.lt.(ptrk-ptrn)')
  end if
  if (ptrk.lt.ptrn) then
     call endrun ('TRUNC: Error in truncation parameters. ptrk.lt.ptrn')
  end if
!
! Evaluate pointers and displacement info based on truncation params
!
  ncoefi(1) = 1
  ik2 = ptrk + 2
  do n=1,pmax
     ncoefi(n+1) = ncoefi(n) + min0(pmmax,ik2-n)
     nalp(n) = ncoefi(n) - 1
     nco2(n) = ncoefi(n)*2
     nm(n) = ncoefi(n+1) - ncoefi(n)
  end do
  nstart(1) = 0
  nlen(1) = ptrn + 1
  do m=2,pmmax
     nstart(m) = nstart(m-1) + nlen(m-1)
     nlen(m) = min0(ptrn+1,ptrk+2-m)
  end do
!
! Define break-even point for vector lengths in GRCALC.  Don't implement
! for non-PVM machines
!
  ncutoff = pmax
!
! Assign wavenumbers if not SPMD.
!
#if ( ! defined SPMD )
   do m=1,pmmax
      locm(m,0) = m
   enddo
#endif

end subroutine trunc

end module dyn_comp
