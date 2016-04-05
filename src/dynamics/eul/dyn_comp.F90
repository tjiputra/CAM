module dyn_comp
!----------------------------------------------------------------------- 
! 
! Eulerian dycore interface module
!
!-----------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
use spmd_utils,   only: masterproc
use constituents, only: pcnst, cnst_name, cnst_longname
use constituents, only: sflxnam, tendnam, fixcnam, tottnam, hadvnam, vadvnam, cnst_get_ind
use pmgrid,       only: plev, plevp
use hycoef,       only: hycoef_init
use cam_history,  only: addfld, add_default, horiz_only
use phys_control, only: phys_getopts

use eul_control_mod, only: dif2, hdif_order, kmnhdn, hdif_coef, divdampn, eps, &
                           kmxhdc, eul_nsplit
use cam_logfile,  only: iulog
use cam_abortutils,  only: endrun

use pio,          only: file_desc_t

#if (defined SPMD)
use spmd_dyn,     only: spmd_readnl, spmdinit_dyn
#endif

implicit none
private

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
CONTAINS
!=============================================================================================

subroutine dyn_readnl(nlfile)

   ! Read dynamics namelist group.
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_integer, mpi_real8

   ! args
   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    
   ! local vars
   integer :: unitn, ierr
     
   real(r8) :: eul_dif2_coef     ! del2 horizontal diffusion coeff.
   integer  :: eul_hdif_order    ! Order of horizontal diffusion operator
   integer  :: eul_hdif_kmnhdn   ! Nth order horizontal diffusion operator top level.
   real(r8) :: eul_hdif_coef     ! Nth order horizontal diffusion coefficient.
   real(r8) :: eul_divdampn      ! Number of days to invoke divergence damper
   real(r8) :: eul_tfilt_eps     ! Time filter coefficient. Defaults to 0.06.
   integer  :: eul_kmxhdc        ! Number of levels to apply Courant limiter
    
   namelist /dyn_eul_inparm/ eul_dif2_coef, eul_hdif_order, eul_hdif_kmnhdn, &
      eul_hdif_coef, eul_divdampn, eul_tfilt_eps, eul_kmxhdc, eul_nsplit

   character(len=*), parameter :: sub = 'dyn_readnl'
   !-----------------------------------------------------------------------------

   ! Read namelist 
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'dyn_eul_inparm', status=ierr)
      if (ierr == 0) then
         read(unitn, dyn_eul_inparm, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub//': ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(eul_dif2_coef, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_dif2_coef")

   call mpi_bcast(eul_hdif_order, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_hdif_order")

   call mpi_bcast(eul_hdif_kmnhdn, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_hdif_kmnhdn")

   call mpi_bcast(eul_hdif_coef, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_hdif_coef")

   call mpi_bcast(eul_divdampn, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_divdampn")

   call mpi_bcast(eul_tfilt_eps, 1, mpi_real8, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_tfilt_eps")

   call mpi_bcast(eul_kmxhdc, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_kmxhdc")

   call mpi_bcast(eul_nsplit, 1, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: eul_nsplit")

   dif2       = eul_dif2_coef
   hdif_order = eul_hdif_order
   kmnhdn     = eul_hdif_kmnhdn
   hdif_coef  = eul_hdif_coef
   divdampn   = eul_divdampn
   eps        = eul_tfilt_eps
   kmxhdc     = eul_kmxhdc

   ! Write namelist variables to logfile
   if (masterproc) then

      write(iulog,*) 'Eulerian Dycore Parameters:'


      ! Order of diffusion
      if (hdif_order < 2 .or. mod(hdif_order, 2) /= 0) then
         write(iulog,*) sub//': Order of diffusion must be greater than 0 and multiple of 2'
         write(iulog,*) 'hdif_order = ', hdif_order
         call endrun(sub//': ERROR: invalid eul_hdif_order specified')
      end if

      if (divdampn > 0._r8) then
         write(iulog,*) '  Divergence damper for spectral dycore invoked for days 0. to ',divdampn,' of this case'
      elseif (divdampn < 0._r8) then
         call endrun (sub//': divdampn must be non-negative')
      else
         write(iulog,*) '  Divergence damper for spectral dycore NOT invoked'
      endif

      if (kmxhdc >= plev .or. kmxhdc < 0) then
         call endrun (sub//':  ERROR:  KMXHDC must be between 0 and plev-1')
      end if

      write(iulog,9108) eps, hdif_order, kmnhdn, hdif_coef, kmxhdc, eul_nsplit

      if (kmnhdn > 1) then
         write(iulog,9109) dif2
      end if

   end if

#if (defined SPMD)
   call spmd_readnl(nlfile)
#endif 

9108 format('   Time filter coefficient (EPS)                 ',f10.3,/,&
            '   Horizontal diffusion order (N)                ',i10/, &
            '   Top layer for Nth order horizontal diffusion  ',i10/, &
            '   Nth order horizontal diffusion coefficient    ',e10.3/, &
            '   Number of levels Courant limiter applied      ',i10/,   &
            '   Dynamics Subcycling                           ',i10)

9109 format('   DEL2 horizontal diffusion applied above Nth order diffusion',/,&
            '   DEL2 Horizontal diffusion coefficient (DIF2)  ',e10.3)


end subroutine dyn_readnl

!=============================================================================================

subroutine dyn_init(file)
   use dyn_grid,     only: define_cam_grids, initgrid
   use scamMod,      only: single_column

   ! ARGUMENTS:
   type(file_desc_t), intent(in) :: file       ! PIO file handle for initial or restart file

   ! Local workspace
   integer m                     ! Index
   integer :: ixcldice, ixcldliq ! constituent indices for cloud liquid and ice water.
   logical :: history_amwg       ! output for AMWG diagnostics
   logical :: history_budget     ! output tendencies and state variables for CAM4
                                 ! temperature, water vapor, cloud ice and cloud
                                 ! liquid budgets.
   integer :: history_budget_histfile_num  ! output history file number for budget fields
   !----------------------------------------------------------------------------

   call trunc()

#if (defined SPMD)
   call spmdinit_dyn()
#endif 

   ! Initialize hybrid coordinate arrays
   call hycoef_init(file)

   ! Run initgrid (the old initcom) which sets up coordinates and weights
   call initgrid()

   ! Define the CAM grids (must be before addfld calls)
   call define_cam_grids()

   call addfld ('ETADOT',(/ 'ilev' /),'A', '1/s','Vertical (eta) velocity',             gridname='gauss_grid')
   call addfld ('U&IC',  (/ 'lev' /), 'I', 'm/s','Zonal wind',                          gridname='gauss_grid' )
   call addfld ('V&IC',  (/ 'lev' /), 'I', 'm/s','Meridional wind',                     gridname='gauss_grid' )
   call add_default ('U&IC',0, 'I')
   call add_default ('V&IC',0, 'I')

   call addfld ('PS&IC',horiz_only,'I',    'Pa','Surface pressure',                     gridname='gauss_grid' )
   call addfld ('T&IC',(/ 'lev' /),'I', 'K','Temperature',                              gridname='gauss_grid' )
   call add_default ('PS&IC',0, 'I')
   call add_default ('T&IC',0, 'I')

   do m = 1, pcnst
      call addfld (trim(cnst_name(m))//'&IC',(/ 'lev' /),'I', 'kg/kg',cnst_longname(m), gridname='gauss_grid' )
      call add_default(trim(cnst_name(m))//'&IC',0, 'I')
      call addfld (hadvnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' horizontal advection tendency',  &
           gridname='gauss_grid')
      call addfld (vadvnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' vertical advection tendency',    &
           gridname='gauss_grid')
      call addfld (tendnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' total tendency',                 &
           gridname='gauss_grid')
      call addfld (tottnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' horz + vert + fixer tendency',   &
           gridname='gauss_grid')
      call addfld (fixcnam(m), (/ 'lev' /), 'A', 'kg/kg/s',trim(cnst_name(m))//' tendency due to slt fixer',      &
           gridname='gauss_grid')
   end do

   call addfld ('DUH     ',(/ 'lev' /),'A', 'K/s     ','U horizontal diffusive heating',              gridname='gauss_grid')
   call addfld ('DVH     ',(/ 'lev' /),'A', 'K/s     ','V horizontal diffusive heating',              gridname='gauss_grid')
   call addfld ('DTH     ',(/ 'lev' /),'A', 'K/s     ','T horizontal diffusive heating',              gridname='gauss_grid')

   call addfld ('ENGYCORR',(/ 'lev' /),'A', 'W/m2    ','Energy correction for over-all conservation', gridname='gauss_grid')
   call addfld ('TFIX    ',horiz_only ,'A', 'K/s     ','T fixer (T equivalent of Energy correction)', gridname='gauss_grid')

   call addfld ('FU      ',(/ 'lev' /),'A', 'm/s2    ','Zonal wind forcing term',                     gridname='gauss_grid')
   call addfld ('FV      ',(/ 'lev' /),'A', 'm/s2    ','Meridional wind forcing term',                gridname='gauss_grid')
   call addfld ('UTEND   ',(/ 'lev' /),'A', 'm/s2    ','U tendency',                                  gridname='gauss_grid')
   call addfld ('VTEND   ',(/ 'lev' /),'A', 'm/s2    ','V tendency',                                  gridname='gauss_grid')
   call addfld ('TTEND   ',(/ 'lev' /),'A', 'K/s     ','T tendency',                                  gridname='gauss_grid')
   call addfld ('LPSTEN  ',horiz_only ,'A', 'Pa/s    ','Surface pressure tendency',                   gridname='gauss_grid')
   call addfld ('VAT     ',(/ 'lev' /),'A', 'K/s     ','Vertical advective tendency of T',            gridname='gauss_grid')
   call addfld ('KTOOP   ',(/ 'lev' /),'A', 'K/s     ','(Kappa*T)*(omega/P)',                         gridname='gauss_grid')

   call phys_getopts(history_amwg_out=history_amwg, &
        history_budget_out = history_budget, &
        history_budget_histfile_num_out = history_budget_histfile_num)

   if (history_amwg) then
      call add_default ('DTH     ', 1, ' ')
   end if

   if ( history_budget ) then
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('CLDICE', ixcldice)
      ! The following variables are not defined for single column
      if (.not. single_column) then 
         call add_default(hadvnam(       1), history_budget_histfile_num, ' ')
         call add_default(hadvnam(ixcldliq), history_budget_histfile_num, ' ')
         call add_default(hadvnam(ixcldice), history_budget_histfile_num, ' ')
         call add_default(vadvnam(       1), history_budget_histfile_num, ' ')
         call add_default(vadvnam(ixcldliq), history_budget_histfile_num, ' ')
         call add_default(vadvnam(ixcldice), history_budget_histfile_num, ' ')
      end if
      call add_default(fixcnam(       1), history_budget_histfile_num, ' ')
      call add_default(fixcnam(ixcldliq), history_budget_histfile_num, ' ')
      call add_default(fixcnam(ixcldice), history_budget_histfile_num, ' ')
      call add_default(tottnam(       1), history_budget_histfile_num, ' ')
      call add_default(tottnam(ixcldliq), history_budget_histfile_num, ' ')
      call add_default(tottnam(ixcldice), history_budget_histfile_num, ' ')
      call add_default(tendnam(       1), history_budget_histfile_num, ' ')
      call add_default(tendnam(ixcldliq), history_budget_histfile_num, ' ')
      call add_default(tendnam(ixcldice), history_budget_histfile_num, ' ')
      call add_default('TTEND',           history_budget_histfile_num, ' ')
      call add_default('TFIX',            history_budget_histfile_num, ' ')
      call add_default('KTOOP',           history_budget_histfile_num, ' ')
      call add_default('VAT',             history_budget_histfile_num, ' ')
      call add_default('DTH',             history_budget_histfile_num, ' ')
   end if

end subroutine dyn_init

!#######################################################################
! Private routines
!#######################################################################

subroutine trunc()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check consistency of truncation parameters and evaluate pointers
! and displacements for spectral arrays
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      L. Bath, June 1992
!                    T. Acker, March 1996
! Reviewed:          J. Hack, D. Williamson, August 1992
! Reviewed:          J. Hack, D. Williamson, April 1996
!-----------------------------------------------------------------------
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use pspect,         only: ptrm, ptrn, ptrk, pmmax
   use comspe,         only: nstart, nlen, locm, lnstart
   use cam_abortutils, only: endrun

   implicit none

!---------------------------Local variables-----------------------------
!
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
      call endrun ('TRUNC: Error in truncation parameters.  ntrm < (ptrk-ptrn)')
   end if
   if (ptrk.lt.ptrn) then
      call endrun ('TRUNC: Error in truncation parameters.  ptrk < ptrn')
   end if
!
! Evaluate pointers and displacement info based on truncation params
!
   nstart(1) = 0
   nlen(1) = ptrn + 1
   do m=2,pmmax
      nstart(m) = nstart(m-1) + nlen(m-1)
      nlen(m) = min0(ptrn+1,ptrk+2-m)
   end do
!
! Assign wavenumbers  and spectral offsets if not SPMD
!
#if ( ! defined SPMD )
   do m=1,pmmax
      locm(m,0) = m
      lnstart(m) = nstart(m)
   enddo
#endif

end subroutine trunc

end module dyn_comp
