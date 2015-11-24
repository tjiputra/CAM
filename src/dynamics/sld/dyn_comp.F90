module dyn_comp
!----------------------------------------------------------------------- 
! 
! Dycore interface module for SLD
!
!-----------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8
use constituents, only: pcnst, cnst_name, cnst_longname
use constituents, only: tendnam, fixcnam, tottnam, hadvnam, vadvnam
use ppgrid,       only: pver, pverp
use pmgrid,       only: plev, plevp, dyndecomp_set
use hycoef,       only: hycoef_init
use cam_history,  only: addfld, add_default, horiz_only
use pio,          only: file_desc_t


implicit none
private

public :: dyn_init, dyn_import_t, dyn_export_t

! these structures are not used in this dycore, but are included
! for source code compatibility.  
type dyn_import_t
   integer :: placeholder
end type dyn_import_t

type dyn_export_t
   integer :: placeholder
end type dyn_export_t

!#######################################################################
CONTAINS
!#######################################################################

subroutine dyn_init(file, nlfilename)

   use spmd_utils,          only: masterproc
   use sld_control_mod,     only: dyn_sld_readnl
   use phys_control,        only: phys_getopts
#if (defined SPMD)
   use spmd_dyn,            only: spmd_readnl,spmdinit_dyn
#endif
   use dyn_grid,            only: define_cam_grids

   ! ARGUMENTS:
   type(file_desc_t), intent(in) :: file       ! PIO file handle for initial or restart file
   character(len=*),  intent(in) :: nlfilename

   logical :: history_amwg       ! output for AMWG diagnostics
   ! Local workspace
   integer m

   dyndecomp_set = .true.

   call dyn_sld_readnl(nlfilename)

#if (defined SPMD)
   call spmd_readnl(nlfilename)
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

   call addfld ('DUH     ',(/ 'lev' /),'A','K/s     ','U horizontal diffusive heating')
   call addfld ('DVH     ',(/ 'lev' /),'A','K/s     ','V horizontal diffusive heating')
   call addfld ('DTH     ',(/ 'lev' /),'A','K/s     ','T horizontal diffusive heating')
   call addfld ('ENGYCORR',(/ 'lev' /),'A','W/m2    ','Energy correction for over-all conservation')
   call addfld ('TFIX    ',horiz_only, 'A','K/s     ','T fixer (T equivalent of Energy correction)')

   call phys_getopts(history_amwg_out=history_amwg)

   if (history_amwg) call add_default ('DTH     ', 1, ' ')

   call addfld ('FU      ',(/ 'lev' /),'A','m/s2    ','Zonal wind forcing term')
   call addfld ('FV      ',(/ 'lev' /),'A','m/s2    ','Meridional wind forcing term')
   call addfld ('UTEND   ',(/ 'lev' /),'A','m/s2    ','U tendency')
   call addfld ('VTEND   ',(/ 'lev' /),'A','m/s2    ','V tendency')
   call addfld ('TTEND   ',(/ 'lev' /),'A','K/s     ','T tendency')
   call addfld ('LPSTEN  ',horiz_only, 'A','Pa/s    ','Surface pressure tendency')
   call addfld ('VAT     ',(/ 'lev' /),'A','K/s     ','Vertical advective tendency of T')
   call addfld ('KTOOP   ',(/ 'lev' /),'A','K/s     ','(Kappa*T)*(omega/P)')

end subroutine dyn_init

end module dyn_comp
