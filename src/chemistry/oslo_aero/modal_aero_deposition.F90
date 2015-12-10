module modal_aero_deposition

!------------------------------------------------------------------------------------------------
! Purpose:
!
! Partition the contributions from modal components of wet and dry 
! deposition at the surface into the fields passed to the coupler.
!
! *** N.B. *** Currently only a simple scheme for the 3-mode version
!              of MAM has been implemented.
!
! Revision history:
! Feb 2009  M. Flanner, B. Eaton   Original version for trop_mam3.
! Jul 2011  F Vitt -- made avaliable to be used in a prescribed modal aerosol mode (no prognostic MAM)
! Mar 2012  F Vitt -- made changes for to prevent abort when 7-mode aeroslol model is used
!                     some of the needed consituents do not exist in 7-mode so bin_fluxes will be false
!------------------------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use camsrfexch,       only: cam_out_t     
use constituents,     only: pcnst, cnst_get_ind
use ppgrid,           only: pcols
use abortutils,       only: endrun
use aerosoldef,       only: l_bc_n,l_bc_ax,l_bc_ni,l_bc_a,l_bc_ai,l_bc_ac
use aerosoldef,       only: l_om_ni,l_om_ai,l_om_ac,l_dst_a2,l_dst_a3

implicit none
private
save

public :: &
   modal_aero_deposition_init, &
   set_srf_drydep,             &
   set_srf_wetdep

! Private module data
integer :: idx_bc1  = -1
integer :: idx_pom1 = -1
integer :: idx_soa1 = -1
integer :: idx_soa2 = -1
integer :: idx_dst1 = -1
integer :: idx_dst3 = -1
integer :: idx_ncl3 = -1
integer :: idx_so43 = -1
logical :: bin_fluxes = .false.

logical :: initialized = .false.

!==============================================================================
contains
!==============================================================================

subroutine modal_aero_deposition_init(bc1_ndx,pom1_ndx,soa1_ndx,soa2_ndx,dst1_ndx,dst3_ndx,ncl3_ndx,so43_ndx,num3_ndx)

! set aerosol indices for re-mapping surface deposition fluxes:
! *_a1 = accumulation mode
! *_a2 = aitken mode
! *_a3 = coarse mode

   ! can be initialized with user specified indices
   ! if called from aerodep_flx module (for prescribed modal aerosol fluxes) then these indices are specified

   integer, optional, intent(in) :: bc1_ndx,pom1_ndx,soa1_ndx,soa2_ndx,dst1_ndx,dst3_ndx,ncl3_ndx,so43_ndx,num3_ndx

   ! if already initialized abort the run
   if (initialized) then
     call endrun('modal_aero_deposition_init is already initialized')
   endif

   if (present(bc1_ndx)) then
      idx_bc1  = bc1_ndx
   else
      call cnst_get_ind('bc_a1',  idx_bc1)
   endif
   if (present(pom1_ndx)) then
      idx_pom1 = pom1_ndx
   else
      call cnst_get_ind('pom_a1', idx_pom1)
   endif
   if (present(soa1_ndx)) then
      idx_soa1 = soa1_ndx
   else
      call cnst_get_ind('soa_a1', idx_soa1)
   endif
   if (present(soa2_ndx)) then
      idx_soa2 = soa2_ndx
   else
      call cnst_get_ind('soa_a2', idx_soa2)
   endif
   if (present(dst1_ndx)) then
      idx_dst1 = dst1_ndx
   else
      call cnst_get_ind('dst_a1', idx_dst1,abort=.false.)
   endif
   if (present(dst3_ndx)) then
      idx_dst3 = dst3_ndx
   else
      call cnst_get_ind('dst_a3', idx_dst3,abort=.false.)
   endif
   if (present(ncl3_ndx)) then
      idx_ncl3 = ncl3_ndx
   else
      call cnst_get_ind('ncl_a3', idx_ncl3,abort=.false.)
   endif
   if (present(so43_ndx)) then
      idx_so43 = so43_ndx
   else
      call cnst_get_ind('so4_a3', idx_so43,abort=.false.)
   endif

!  for 7 mode bin_fluxes will be false
   bin_fluxes = idx_dst1>0 .and. idx_dst3>0 .and.idx_ncl3>0 .and. idx_so43>0
   initialized = .true.

end subroutine modal_aero_deposition_init

!==============================================================================
subroutine set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)

! Set surface wet deposition fluxes passed to coupler.

   ! Arguments:
   real(r8), intent(in) :: aerdepwetis(:,:)  ! aerosol wet deposition (interstitial)
   real(r8), intent(in) :: aerdepwetcw(:,:)  ! aerosol wet deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i
   integer :: ncol                      ! number of columns
   !----------------------------------------------------------------------------

   ncol = cam_out%ncol

   do i = 1, ncol
      ! black carbon fluxes
      cam_out%bcphiwet(i) = -(aerdepwetis(i,l_bc_ni)+aerdepwetcw(i,l_bc_ni)+ &
      aerdepwetis(i,l_bc_ai)+aerdepwetcw(i,l_bc_ai)+aerdepwetis(i,l_bc_a)+aerdepwetcw(i,l_bc_a)+aerdepwetis(i,l_bc_ac)+aerdepwetcw(i,l_bc_ac))

      ! organic carbon fluxes
      cam_out%ocphiwet(i) = -(aerdepwetis(i,l_om_ni)+aerdepwetcw(i,l_om_ni)+ &
      aerdepwetis(i,l_om_ai)+aerdepwetcw(i,l_om_ai)+aerdepwetis(i,l_om_ac)+aerdepwetcw(i,l_om_ac))

      ! dust fluxes
      !
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      cam_out%dstwet1(i) = -(aerdepwetis(i,l_dst_a2)+aerdepwetcw(i,l_dst_a2))
      
      !  A. Simple: Assign all coarse-mode dust to bulk size bin 3:
      cam_out%dstwet2(i) = 0._r8
      cam_out%dstwet3(i) = -(aerdepwetis(i,l_dst_a3)+aerdepwetcw(i,l_dst_a3))
      cam_out%dstwet4(i) = 0._r8

   enddo

end subroutine set_srf_wetdep

!==============================================================================

subroutine set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)

! Set surface dry deposition fluxes passed to coupler.
   
   ! Arguments:
   real(r8), intent(in) :: aerdepdryis(:,:)  ! aerosol dry deposition (interstitial)
   real(r8), intent(in) :: aerdepdrycw(:,:)  ! aerosol dry deposition (cloud water)
   type(cam_out_t), intent(inout) :: cam_out     ! cam export state

   ! Local variables:
   integer :: i
   integer :: ncol                      ! number of columns
   !----------------------------------------------------------------------------

   ncol = cam_out%ncol

   ! derive cam_out variables from deposition fluxes
   !  note: wet deposition fluxes are negative into surface, 
   !        dry deposition fluxes are positive into surface.
   !        CLM wants positive definite fluxes.
   do i = 1, ncol
      ! black carbon fluxes
      cam_out%bcphidry(i) = aerdepdryis(i,l_bc_ni)+aerdepdrycw(i,l_bc_ni)+ &
      aerdepdryis(i,l_bc_ai)+aerdepdrycw(i,l_bc_ai)+aerdepdryis(i,l_bc_a)+aerdepdrycw(i,l_bc_a)+aerdepdryis(i,l_bc_ac)+aerdepdrycw(i,l_bc_ac)
      cam_out%bcphodry(i) = aerdepdryis(i,l_bc_n)+aerdepdrycw(i,l_bc_n)+ aerdepdryis(i,l_bc_ax)+aerdepdrycw(i,l_bc_ax)

      ! organic carbon fluxes
      cam_out%ocphidry(i) = aerdepdryis(i,l_om_ni)+aerdepdrycw(i,l_om_ni)+ &
      aerdepdryis(i,l_om_ai)+aerdepdrycw(i,l_om_ai)+aerdepdryis(i,l_bc_a)+aerdepdrycw(i,l_bc_a)+aerdepdryis(i,l_om_ac)+aerdepdrycw(i,l_om_ac)
      cam_out%ocphidry(i) = 0._r8 !aerdepdryis(i,idx_pom1)+aerdepdryis(i,idx_soa1)+aerdepdrycw(i,idx_pom1)+aerdepdrycw(i,idx_soa1)
      cam_out%ocphodry(i) = 0._r8 !aerdepdryis(i,idx_soa2)+aerdepdrycw(i,idx_soa2)

      ! dust fluxes
      !
      ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
      cam_out%dstdry1(i) = aerdepdryis(i,l_dst_a2)+aerdepdrycw(i,l_dst_a2)
      
      ! Two options for partitioning deposition into bins 2-4:
      !  A. Simple: Assign all coarse-mode dust to bulk size bin 3:
      cam_out%dstdry2(i) = 0._r8
      cam_out%dstdry3(i) = aerdepdryis(i,l_dst_a3)+aerdepdrycw(i,l_dst_a3)
      cam_out%dstdry4(i) = 0._r8
   enddo

end subroutine set_srf_drydep


!==============================================================================

end module modal_aero_deposition
