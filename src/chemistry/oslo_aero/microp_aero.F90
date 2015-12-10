module microp_aero

!---------------------------------------------------------------------------------
! Purpose:
!   CAM Interface for aerosol activation 
!
! ***N.B.*** This module is currently hardcoded to recognize only the aerosols/modes that
!            affect the climate calculation.  This is implemented by using list
!            index 0 in all the calls to rad_constituent interfaces.
!
! Author: Andrew Gettelman
! Based on code from: Hugh Morrison, Xiaohong Liu and Steve Ghan
! May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)         
! for questions contact Andrew Gettelman  (andrew@ucar.edu)
! Modifications: A. Gettelman Nov 2010  - changed to support separation of 
!                microphysics and macrophysics and concentrate aerosol information here
!
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp
use physconst,        only: rair, tmelt
use constituents,     only: cnst_get_ind, pcnst
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,     only: phys_getopts
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num, rad_cnst_get_mode_props
use shr_spfn_mod,     only: erf => shr_spfn_erf_nonintrinsic, &
                            erfc => shr_spfn_erfc_nonintrinsic
use wv_saturation,    only: qsat_water
use nucleate_ice,     only: nucleati
use ndrop,            only: ndrop_init, dropmixnuc
use ndrop_bam,        only: ndrop_bam_init, ndrop_bam_run, ndrop_bam_ccn
use cam_history,      only: addfld, phys_decomp, add_default, outfld
use cam_logfile,      only: iulog
use abortutils,       only: endrun
use ref_pres,         only: top_lev => trop_cloud_top_lev

use classnuc,         only: classnuc_in, preexisting_ice

#ifdef OSLO_AERO
use commondefinitions, only:  nmodes_oslo => nmodes
use aerosoldef, only: MODE_IDX_DST_A2, MODE_IDX_DST_A3, MODE_IDX_SO4_AC &
                      ,MODE_IDX_OMBC_INTMIX_COAT_AIT, lifeCycleNumberMedianRadius, &
                      l_dst_a2, l_dst_a3, l_bc_ai, getNumberOfTracersInMode, &
                      getTracerIndex, getCloudTracerIndex
use oslo_utils, only: CalculateNumberConcentration
use parmix_progncdnc
#endif

implicit none
private
save

public :: microp_aero_init, microp_aero_run, microp_aero_readnl, microp_aero_register

!++ MH_2015/04/10
integer, public :: frzimm_idx, frzcnt_idx, frzdep_idx
!-- MH_2015/04/10

! Private module data

character(len=16)   :: eddy_scheme  ! eddy scheme

! contact freezing due to dust
! dust number mean radius (m), Zender et al JGR 2003 assuming number mode radius of 0.6 micron, sigma=2
real(r8), parameter :: rn_dst1 = 0.258e-6_r8
real(r8), parameter :: rn_dst2 = 0.717e-6_r8
real(r8), parameter :: rn_dst3 = 1.576e-6_r8
real(r8), parameter :: rn_dst4 = 3.026e-6_r8

real(r8), public :: bulk_scale    ! prescribed aerosol bulk sulfur scale factor

! smallest mixing ratio considered in microphysics
real(r8), parameter :: qsmall = 1.e-18_r8

! minimum allowed cloud fraction
real(r8), parameter :: mincld = 0.0001_r8

! indices in state%q and pbuf structures
integer :: cldliq_idx = -1
integer :: cldice_idx = -1
integer :: numliq_idx = -1
integer :: numice_idx = -1
integer :: kvh_idx = -1
integer :: tke_idx = -1
integer :: wp2_idx = -1
integer :: ast_idx = -1
integer :: cldo_idx = -1
integer :: dgnum_idx    = -1
integer :: dgnumwet_idx = -1

! Bulk aerosols
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all      ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate
integer :: idxdst1  = -1 ! index in aerosol list for dust1
integer :: idxdst2  = -1 ! index in aerosol list for dust2
integer :: idxdst3  = -1 ! index in aerosol list for dust3
integer :: idxdst4  = -1 ! index in aerosol list for dust4
integer :: idxbcphi = -1 ! index in aerosol list for Soot (BCPHIL)

! modal aerosols
logical :: prog_modal_aero
logical :: clim_modal_aero

integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
integer :: mode_coarse_idx = -1  ! index of coarse mode
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode

integer :: naai_idx, naai_hom_idx, npccn_idx, rndst_idx, nacon_idx

real(r8) :: sigmag_aitken
logical  :: separate_dust = .false.

!===============================================================================
contains
!===============================================================================
subroutine microp_aero_register
   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Register pbuf fields for aerosols needed by microphysics
   ! 
   ! Author: Cheryl Craig October 2012
   ! 
   !-----------------------------------------------------------------------
   use ppgrid,         only: pcols
   use physics_buffer, only: pbuf_add_field, dtype_r8

   call pbuf_add_field('NAAI',       'physpkg',dtype_r8,(/pcols,pver/), naai_idx)
   call pbuf_add_field('NAAI_HOM',   'physpkg',dtype_r8,(/pcols,pver/), naai_hom_idx)
   call pbuf_add_field('NPCCN',      'physpkg',dtype_r8,(/pcols,pver/), npccn_idx)
   call pbuf_add_field('RNDST',      'physpkg',dtype_r8,(/pcols,pver,4/), rndst_idx)
   call pbuf_add_field('NACON',      'physpkg',dtype_r8,(/pcols,pver,4/), nacon_idx)
 
!++ MH_2015/04/10
   call pbuf_add_field('FRZIMM','physpkg',dtype_r8,(/pcols,pver/), frzimm_idx)
   call pbuf_add_field('FRZCNT','physpkg',dtype_r8,(/pcols,pver/), frzcnt_idx)
   call pbuf_add_field('FRZDEP','physpkg',dtype_r8,(/pcols,pver/), frzdep_idx)
!-- MH_2015/04/10

end subroutine microp_aero_register

subroutine microp_aero_init

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Initialize constants for aerosols needed by microphysics
   ! 
   ! Author: Andrew Gettelman May 2010
   ! 
   !-----------------------------------------------------------------------

   ! local variables
   integer  :: iaer
   integer  :: m, n, nmodes, nspec

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'microp_aero_init'
   logical :: history_amwg
   !-----------------------------------------------------------------------

   ! Query the PBL eddy scheme
   call phys_getopts(eddy_scheme_out          = eddy_scheme, &
                     history_amwg_out = history_amwg)

   ! Access the physical properties of the aerosols that are affecting the climate
   ! by using routines from the rad_constituents module.

   ! get indices into state and pbuf structures
   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMLIQ', numliq_idx)
   call cnst_get_ind('NUMICE', numice_idx)

   select case(trim(eddy_scheme))
   case ('diag_TKE')
      tke_idx      = pbuf_get_index('tke')
   case ('CLUBB_SGS')
      wp2_idx = pbuf_get_index('WP2')
   case default
      kvh_idx      = pbuf_get_index('kvh')
   end select

   ! prog_modal_aero determines whether prognostic modal aerosols are present in the run.
   call phys_getopts(prog_modal_aero_out=prog_modal_aero)

   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   ast_idx      = pbuf_get_index('AST')

#if (defined OSLO_AERO)
      cldo_idx     = pbuf_get_index('CLDO')
      clim_modal_aero = .true. !Needed to avoid ending up in BAM routines

      call ndrop_init()
#else
   if (clim_modal_aero) then

      cldo_idx     = pbuf_get_index('CLDO')
      dgnum_idx    = pbuf_get_index('DGNUM' )
      dgnumwet_idx = pbuf_get_index('DGNUMWET')

      call ndrop_init()

      ! Init indices for specific modes/species

      ! mode index for specified mode types
      do m = 1, nmodes
         call rad_cnst_get_info(0, m, mode_type=str32)
         select case (trim(str32))
         case ('accum')
            mode_accum_idx = m
         case ('aitken')
            mode_aitken_idx = m
         case ('coarse')
            mode_coarse_idx = m
         case ('coarse_dust')
            mode_coarse_dst_idx = m
         case ('coarse_seasalt')
            mode_coarse_slt_idx = m
         end select
      end do

      ! check if coarse dust is in separate mode
      separate_dust = mode_coarse_dst_idx > 0

      ! for 3-mode 
      if ( mode_coarse_dst_idx<0 ) mode_coarse_dst_idx = mode_coarse_idx
      if ( mode_coarse_slt_idx<0 ) mode_coarse_slt_idx = mode_coarse_idx

      ! Check that required mode types were found
      if (mode_accum_idx == -1 .or. mode_aitken_idx == -1 .or. &
          mode_coarse_dst_idx == -1.or. mode_coarse_slt_idx == -1) then
         write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
            mode_accum_idx, mode_aitken_idx, mode_coarse_dst_idx, mode_coarse_slt_idx
         call endrun(routine//': ERROR required mode type not found')
      end if

      ! species indices for specified types
      ! find indices for the dust and seasalt species in the coarse mode
      call rad_cnst_get_info(0, mode_coarse_dst_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_dst_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('dust')
            coarse_dust_idx = n
         end select
      end do
      call rad_cnst_get_info(0, mode_coarse_slt_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_slt_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('seasalt')
            coarse_nacl_idx = n
         end select
      end do

      ! Check that required mode specie types were found
      if ( coarse_dust_idx == -1 .or. coarse_nacl_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_dust_idx, coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

      ! get specific mode properties
      call rad_cnst_get_mode_props(0, mode_aitken_idx, sigmag=sigmag_aitken)

   else

      ! Props needed for BAM number concentration calcs.

      call rad_cnst_get_info(0, naero=naer_all)
      allocate( &
         aername(naer_all),        &
         num_to_mass_aer(naer_all) )

      do iaer = 1, naer_all
         call rad_cnst_get_aer_props(0, iaer, &
            aername         = aername(iaer), &
            num_to_mass_aer = num_to_mass_aer(iaer) )

         ! Look for sulfate, dust, and soot in this list (Bulk aerosol only)
         if (trim(aername(iaer)) == 'SULFATE') idxsul = iaer
         if (trim(aername(iaer)) == 'DUST1') idxdst1 = iaer
         if (trim(aername(iaer)) == 'DUST2') idxdst2 = iaer
         if (trim(aername(iaer)) == 'DUST3') idxdst3 = iaer
         if (trim(aername(iaer)) == 'DUST4') idxdst4 = iaer
         if (trim(aername(iaer)) == 'BCPHIL') idxbcphi = iaer
      end do

      call ndrop_bam_init()

   end if

#endif 

   call addfld('LCLOUD', ' ', pver, 'A', 'Liquid cloud fraction used in stratus activation', phys_decomp)

   call addfld('WSUB     ', 'm/s     ', pver, 'A', 'Diagnostic sub-grid vertical velocity'                   ,phys_decomp)
   call addfld('WSUBI    ', 'm/s     ', pver, 'A', 'Diagnostic sub-grid vertical velocity for ice'           ,phys_decomp)
   call addfld('NIHF',  '1/m3', pver, 'A', 'Activated Ice Number Concentation due to homogenous freezing',  phys_decomp)
   call addfld('NIDEP', '1/m3', pver, 'A', 'Activated Ice Number Concentation due to deposition nucleation',phys_decomp)
   call addfld('NIIMM', '1/m3', pver, 'A', 'Activated Ice Number Concentation due to immersion freezing',   phys_decomp)
   call addfld('NIMEY', '1/m3', pver, 'A', 'Activated Ice Number Concentation due to meyers deposition',    phys_decomp)

   if (history_amwg) then
      call add_default ('WSUB     ', 1, ' ')
   end if

!++ MH_2015/04/10
if(classnuc_in) then
    call addfld('FREQIMM', 'fraction', pver, 'A', 'Fractional occurance of immersion  freezing', phys_decomp)
    call addfld('FREQCNT', 'fraction', pver, 'A', 'Fractional occurance of contact    freezing', phys_decomp)
    call addfld('FREQDEP', 'fraction', pver, 'A', 'Fractional occurance of deposition freezing', phys_decomp)
    call addfld('FREQMIX', 'fraction', pver, 'A', 'Fractional occurance of mixed-phase clouds' , phys_decomp)
    call add_default('FREQIMM', 1, ' ')
    call add_default('FREQCNT', 1, ' ')
    call add_default('FREQDEP', 1, ' ')
    call add_default('FREQMIX', 1, ' ')

    call addfld('DSTFREZIMM', 'm-3s-1', pver, 'A', 'dust immersion freezing rate', phys_decomp)
    call addfld('DSTFREZCNT', 'm-3s-1', pver, 'A', 'dust contact freezing rate', phys_decomp)
    call addfld('DSTFREZDEP', 'm-3s-1', pver, 'A', 'dust deposition freezing rate', phys_decomp)
    call add_default('DSTFREZIMM', 1, ' ')
    call add_default('DSTFREZCNT', 1, ' ')
    call add_default('DSTFREZDEP', 1, ' ')

    call addfld('BCFREZIMM', 'm-3s-1', pver, 'A', 'bc immersion freezing rate', phys_decomp)
    call addfld('BCFREZCNT', 'm-3s-1', pver, 'A', 'bc contact freezing rate', phys_decomp)
    call addfld('BCFREZDEP', 'm-3s-1', pver, 'A', 'bc deposition freezing rate', phys_decomp)
    call add_default('BCFREZIMM', 1, ' ')
    call add_default('BCFREZCNT', 1, ' ')
    call add_default('BCFREZDEP', 1, ' ')

    call addfld('NIMIX_IMM', '#/m3', pver, 'A', 'Activated Ice Number Concentration due to het immersion freezing in Mixed Clouds', phys_decomp)
    call add_default('NIMIX_IMM', 1, ' ')
    call addfld('NIMIX_CNT', '#/m3', pver, 'A', 'Activated Ice Number Concentration due to het contact freezing in Mixed Clouds', phys_decomp)
    call add_default('NIMIX_CNT', 1, ' ')  
    call addfld('NIMIX_DEP', '#/m3', pver, 'A', 'Activated Ice Number Concentration due to het deposition freezing in Mixed Clouds', phys_decomp)
    call add_default('NIMIX_DEP', 1, ' ')

    call addfld('DSTNIDEP', '#/m3', pver, 'A', 'Activated Ice Number Concentration due to dst dep freezing in Mixed Clouds', phys_decomp)
    call add_default('DSTNIDEP', 1, ' ')
    call addfld('DSTNICNT', '#/m3', pver, 'A', 'Activated Ice Number Concentration due to dst cnt freezing in Mixed Clouds', phys_decomp)
    call add_default('DSTNICNT', 1, ' ')
    call addfld('DSTNIIMM', '#/m3', pver, 'A', 'Activated Ice Number Concentration due to dst imm freezing in Mixed Clouds', phys_decomp)
    call add_default('DSTNIIMM', 1, ' ')

    call addfld('BCNIDEP', '#/m3', pver, 'A', 'Activated Ice Number Concentration due to bc dep freezing in Mixed Clouds', phys_decomp)
    call add_default('BCNIDEP', 1, ' ')
    call addfld('BCNICNT', '#/m3', pver, 'A', 'Activated Ice Number Concentration due to bc cnt freezing in Mixed Clouds', phys_decomp)
    call add_default('BCNICNT', 1, ' ')
    call addfld('BCNIIMM', '#/m3', pver, 'A', 'Activated Ice Number Concentration due to bc imm freezing in Mixed Clouds', phys_decomp)
    call add_default('BCNIIMM', 1, ' ')

    call addfld('NUMICE10s', '#/m3', pver, 'A', 'Ice Number Concentration due to het freezing in Mixed Clouds during 10-s period', phys_decomp)
    call add_default('NUMICE10s', 1, ' ')
    call addfld('NUMIMM10sDST', '#/m3', pver, 'A', 'Ice Number Concentration due to imm freezing by dst in Mixed Clouds during 10-s period', phys_decomp)
    call add_default('NUMIMM10sDST', 1, ' ')
    call addfld('NUMIMM10sBC', '#/m3', pver, 'A', 'Ice Number Concentration due to imm freezing by bc in Mixed Clouds during 10-s period', phys_decomp)
    call add_default('NUMIMM10sBC', 1, ' ')
end if


   if(preexisting_ice) then
     call addfld('fhom     ', 'fraction', pver, 'A', 'Fraction of cirrus where homogeneous freezing occur'   ,phys_decomp) 
     call addfld ('WICE      ', 'm/s   ', pver, 'A','Vertical velocity Reduction caused by preexisting ice'  ,phys_decomp)
     call addfld ('WEFF      ', 'm/s   ', pver, 'A','Effective Vertical velocity for ice nucleation' ,phys_decomp)
     call add_default ('fhom     ', 1, ' ') 
     call add_default ('WICE    ', 1, ' ')
     call add_default ('WEFF    ', 1, ' ')
   endif
!-- MH_2015/04/10
   
end subroutine microp_aero_init

!===============================================================================

subroutine microp_aero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Namelist variables
   real(r8) :: microp_aero_bulk_scale = 2._r8  ! prescribed aerosol bulk sulfur scale factor
 
   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'microp_aero_readnl'

   namelist /microp_aero_nl/ microp_aero_bulk_scale
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'microp_aero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, microp_aero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variable
   call mpibcast(microp_aero_bulk_scale, 1, mpir8, 0, mpicom)
#endif

   ! set local variables
   bulk_scale = microp_aero_bulk_scale

end subroutine microp_aero_readnl

!===============================================================================

subroutine microp_aero_run ( &
   state, ptend, deltatin, pbuf)

   !++ MH_2015/04/10
   use modal_aero_data, only: qqcw_get_field
   use physconst,       only: rhoh2o 
   use wv_saturation,   only: svp_water, svp_ice
   !-- MH_2015/04/10
   
   ! input arguments
   type(physics_state), target, intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: deltatin     ! time step (s)
   type(physics_buffer_desc),   pointer       :: pbuf(:)




   ! local workspace
   ! all units mks unless otherwise stated

   integer :: i, k, m, n
   integer :: itim
   integer :: lchnk
   integer :: ncol
   integer :: nmodes
   integer :: nucboast

   real(r8), pointer :: ast(:,:)        

   real(r8)          :: icecldf(pcols,pver)    ! ice cloud fraction   
   real(r8)          :: liqcldf(pcols,pver)    ! liquid cloud fraction

   real(r8), pointer :: naai(:,:)       ! number of activated aerosol for ice nucleation 
   real(r8), pointer :: naai_hom(:,:)   ! number of activated aerosol for ice nucleation (homogeneous freezing only)
   real(r8), pointer :: npccn(:,:)      ! number of CCN (liquid activated)
   real(r8), pointer :: rndst(:,:,:)    ! radius of 4 dust bins for contact freezing
   real(r8), pointer :: nacon(:,:,:)    ! number in 4 dust bins for contact freezing

   real(r8), pointer :: t(:,:)          ! input temperature (K)
   real(r8), pointer :: qn(:,:)         ! input water vapor mixing ratio (kg/kg)
   ! note: all input cloud variables are grid-averaged
   real(r8), pointer :: qc(:,:)         ! cloud water mixing ratio (kg/kg)
   real(r8), pointer :: qi(:,:)         ! cloud ice mixing ratio (kg/kg)
   real(r8), pointer :: nc(:,:)         ! cloud water number conc (1/kg)
   real(r8), pointer :: ni(:,:)         ! cloud ice number conc (1/kg)
   real(r8), pointer :: pmid(:,:)       ! pressure at layer midpoints (pa)
   real(r8), pointer :: pdel(:,:)       ! pressure difference across level (pa)
   real(r8), pointer :: pint(:,:)       ! air pressure layer interfaces (pa)
   real(r8), pointer :: rpdel(:,:)      ! inverse pressure difference across level (pa)
   real(r8), pointer :: zm(:,:)         ! geopotential height of model levels (m)
   real(r8), pointer :: omega(:,:)      ! vertical velocity (Pa/s)
   real(r8), pointer :: num_accum(:,:)  ! number m.r. of accumulation mode
   real(r8), pointer :: num_aitken(:,:) ! number m.r. of aitken mode
   real(r8), pointer :: num_coarse(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: coarse_dust(:,:) ! mass m.r. of coarse dust
   real(r8), pointer :: coarse_nacl(:,:) ! mass m.r. of coarse nacl

   real(r8), pointer :: kvh(:,:)        ! vertical eddy diff coef (m2 s-1)
   real(r8), pointer :: tke(:,:)        ! TKE from the UW PBL scheme (m2 s-2)
   real(r8), pointer :: wp2(:,:)        ! CLUBB vertical velocity variance

   real(r8), pointer :: cldn(:,:)       ! cloud fraction
   real(r8), pointer :: cldo(:,:)       ! old cloud fraction

   real(r8), pointer :: dgnum(:,:,:)    ! aerosol mode dry diameter
   real(r8), pointer :: dgnumwet(:,:,:) ! aerosol mode diameter

   real(r8), pointer :: aer_mmr(:,:)    ! aerosol mass mixing ratio

   real(r8) :: rho(pcols,pver)     ! air density (kg m-3)
   real(r8) :: relhum(pcols,pver)  ! relative humidity
   real(r8) :: icldm(pcols,pver)   ! ice cloud fraction
   real(r8) :: lcldm(pcols,pver)   ! liq cloud fraction
   real(r8) :: nfice(pcols,pver)   ! fice variable
   real(r8) :: dumfice             ! dummy var in fice calc
   real(r8) :: lcldn(pcols,pver)   ! fractional coverage of new liquid cloud
   real(r8) :: lcldo(pcols,pver)   ! fractional coverage of old liquid cloud
   real(r8) :: qcld                ! total cloud water
   real(r8) :: nctend_mixnuc(pcols,pver)
   real(r8) :: dum, dum2           ! temporary dummy variable
   real(r8) :: dmc, ssmc           ! variables for modal scheme.

   real(r8) :: so4_num                               ! so4 aerosol number (#/cm^3)
   real(r8) :: soot_num                              ! soot (hydrophilic) aerosol number (#/cm^3)
   real(r8) :: dst1_num,dst2_num,dst3_num,dst4_num   ! dust aerosol number (#/cm^3)
   real(r8) :: dst_num                               ! total dust aerosol number (#/cm^3)

   real(r8) :: qs(pcols)            ! liquid-ice weighted sat mixing rat (kg/kg)
   real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
   real(r8) :: gammas(pcols)        ! parameter for cond/evap of cloud water

   ! bulk aerosol variables
   real(r8), allocatable :: naer2(:,:,:)    ! bulk aerosol number concentration (1/m3)
   real(r8), allocatable :: maerosol(:,:,:) ! bulk aerosol mass conc (kg/m3)

   real(r8) :: wsub(pcols,pver)    ! diagnosed sub-grid vertical velocity st. dev. (m/s)
   real(r8) :: wsubi(pcols,pver)   ! diagnosed sub-grid vertical velocity ice (m/s)

   ! history output for ice nucleation
   real(r8) :: nihf(pcols,pver)  !output number conc of ice nuclei due to heterogenous freezing (1/m3)
   real(r8) :: niimm(pcols,pver) !output number conc of ice nuclei due to immersion freezing (hetero nuc) (1/m3)
   real(r8) :: nidep(pcols,pver) !output number conc of ice nuclei due to deoposion nucleation (hetero nuc) (1/m3)
   real(r8) :: nimey(pcols,pver) !output number conc of ice nuclei due to meyers deposition (1/m3)

   real(r8) :: wght
   
   !++ MH_2015/08/17
   real(r8) :: fhom(pcols,pver)    ! how much fraction of cloud can reach Shom
   real(r8) :: wice(pcols,pver)    ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at Shom 
   real(r8) :: weff(pcols,pver)    ! effective Vertical velocity for ice nucleation (m/s); weff=wsubi-wice
   !-- MH_2015/08/17
   
   !++ MH_2015/04/10
#ifdef OSLO_AERO
   real(r8) :: fn(pcols,pver,0:nmodes_oslo)
#else
   real(r8) :: fn(pcols,pver,nmodes)
#endif
   real(r8) :: awcam(pcols,pver,3),awfacm(pcols,pver,3)
   real(r8) :: fn_in(3)
   real(r8) :: hetraer(pcols,pver,3),dstcoat(pcols,pver,3)
   real(r8) :: total_interstitial_aer_num(pcols,pver,3),total_cloudborne_aer_num(pcols,pver,3)
   real(r8) :: total_aer_num(pcols,pver,3),coated_aer_num(pcols,pver,3),uncoated_aer_num(pcols,pver,3)
   type qqcw_type
   real(r8), pointer :: fldcw(:,:) 
   end type qqcw_type
   type(qqcw_type) :: qqcw(pcnst)
   real(r8) :: qaercwpt(pcols,pver,pcnst)
   integer :: kk
   
   real(r8), parameter :: pi = 3.14159265358979323846_r8 
   real(r8) :: con1, r3lx, mi0l, supersatice
   real(r8), parameter :: rhow = rhoh2o
   real(r8) :: qcic, ncic
   
   real(r8), pointer :: frzimm(:,:), frzcnt(:,:), frzdep(:,:)
   real(r8) :: frzbcimm(pcols,pver), frzduimm(pcols,pver)
   real(r8) :: frzbccnt(pcols,pver), frzducnt(pcols,pver)
   real(r8) :: frzbcdep(pcols,pver), frzdudep(pcols,pver)
   
   real(r8) :: freqimm(pcols,pver), freqcnt(pcols,pver), freqdep(pcols,pver), freqmix(pcols,pver)
   real(r8) :: nnuccc_bc(pcols,pver), nnucct_bc(pcols,pver), nnudep_bc(pcols,pver)
   real(r8) :: nnuccc_dst(pcols,pver), nnucct_dst(pcols,pver), nnudep_dst(pcols,pver)
   real(r8) :: niimm_bc(pcols,pver), nicnt_bc(pcols,pver), nidep_bc(pcols,pver)
   real(r8) :: niimm_dst(pcols,pver), nicnt_dst(pcols,pver), nidep_dst(pcols,pver)
   real(r8) :: numice10s(pcols,pver)
   real(r8) :: numice10s_imm_dst(pcols,pver)
   real(r8) :: numice10s_imm_bc(pcols,pver)
   !++ wy4.0
   real(r8) :: na500(pcols,pver)
   real(r8) :: tot_na500(pcols,pver)
   !-- wy4.0
   !-- MH_2015/04/10
   
   !++ MH_2015/04/10
#ifdef OSLO_AERO
      logical  :: hasAerosol(pcols, pver, nmodes_oslo)
      real(r8) :: f_acm(pcols,pver, nmodes_oslo)
      real(r8) :: f_bcm(pcols,pver, nmodes_oslo)
      real(r8) :: f_aqm(pcols, pver, nmodes_oslo)
      real(r8) :: f_so4_condm(pcols, pver, nmodes_oslo)         !Needed in "get component fraction"
      real(r8) :: f_soam(pcols, pver, nmodes_oslo)              !Needed in "get component fraction"
      real(r8) :: numberConcentration(pcols,pver,0:nmodes_oslo) ![#/m3] number concentraiton
      real(r8) :: volumeConcentration(pcols,pver,nmodes_oslo)   ![m3/m3] volume concentration
      real(r8) :: hygroscopicity(pcols,pver,nmodes_oslo)        ![mol_{aer}/mol_{water}] hygroscopicity
      real(r8) :: lnsigma(pcols,pver,nmodes_oslo)              ![-] log(base e) sigma
      real(r8) :: CProcessModes(pcols,pver)
      real(r8) :: cam(pcols,pver,nmodes_oslo)
      real(r8) :: f_c(pcols, pver)
      real(r8) :: f_aq(pcols,pver)
      real(r8) :: f_bc(pcols,pver)
      real(r8) :: f_so4_cond(pcols,pver)
      real(r8) :: f_soa(pcols,pver)
      real(r8) :: volumeCore(pcols,pver,nmodes_oslo)
      real(r8) :: volumeCoat(pcols,pver,nmodes_oslo)
      real(r8) :: sigmag_amode(3)
      real(r8) :: CloudnumberConcentration(pcols,pver,0:nmodes_oslo)
      
      real(r8) :: fn_bc(pcols,pver), fn_dst1(pcols,pver), fn_dst3(pcols,pver)
      real(r8) :: hetraer_bc(pcols,pver), hetraer_dst1(pcols,pver), hetraer_dst3(pcols,pver)
      real(r8) :: dstcoat_bc(pcols,pver), dstcoat_dst1(pcols,pver), dstcoat_dst3(pcols,pver)
#endif
      !-- MH_2015/04/10

   !-------------------------------------------------------------------------------
   
   lchnk = state%lchnk
   ncol  = state%ncol
   t     => state%t
   qn    => state%q(:,:,1)
   qc    => state%q(:,:,cldliq_idx)
   qi    => state%q(:,:,cldice_idx)
   nc    => state%q(:,:,numliq_idx)
   ni    => state%q(:,:,numice_idx)
   pmid  => state%pmid
   pdel  => state%pdel
   pint  => state%pint
   rpdel => state%rpdel
   zm    => state%zm
   omega => state%omega

   itim = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx,      ast, start=(/1,1,itim/), kount=(/pcols,pver,1/))

   liqcldf(:ncol,:pver) = ast(:ncol,:pver)
   icecldf(:ncol,:pver) = ast(:ncol,:pver)

   call pbuf_get_field(pbuf, naai_idx, naai)
   call pbuf_get_field(pbuf, naai_hom_idx, naai_hom)
   call pbuf_get_field(pbuf, npccn_idx, npccn)
   call pbuf_get_field(pbuf, nacon_idx, nacon)
   call pbuf_get_field(pbuf, rndst_idx, rndst)

   if (clim_modal_aero) then
      itim = pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, ast_idx,  cldn, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
      call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
#ifndef OSLO_AERO
      call rad_cnst_get_info(0, nmodes=nmodes)
      call pbuf_get_field(pbuf, dgnum_idx,    dgnum,    start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )
#endif
   end if

   ! initialize output
   naai(1:ncol,1:pver)     = 0._r8  
   naai_hom(1:ncol,1:pver) = 0._r8  
   npccn(1:ncol,1:pver)    = 0._r8  
   nacon(1:ncol,1:pver,:)  = 0._r8

   ! set default or fixed dust bins for contact freezing
   rndst(1:ncol,1:pver,1) = rn_dst1
   rndst(1:ncol,1:pver,2) = rn_dst2
   rndst(1:ncol,1:pver,3) = rn_dst3
   rndst(1:ncol,1:pver,4) = rn_dst4

   ! initialize history output fields for ice nucleation
   nihf(1:ncol,1:pver)  = 0._r8  
   niimm(1:ncol,1:pver) = 0._r8  
   nidep(1:ncol,1:pver) = 0._r8 
   nimey(1:ncol,1:pver) = 0._r8 

   ! initialize time-varying parameters
   do k = top_lev, pver
      do i = 1, ncol
         rho(i,k) = pmid(i,k)/(rair*t(i,k))
      end do
   end do

!++ MH_2015/04/10
#ifdef OSLO_AERO
    fn(1:ncol,1:pver,0:nmodes_oslo) = 0._r8
#else
    fn(1:ncol,1:pver,1:nmodes) = 0._r8
#endif    
    hetraer(1:ncol,1:pver,1:3) = 0._r8
    total_aer_num(1:ncol,1:pver,1:3) = 0._r8
    coated_aer_num(1:ncol,1:pver,1:3) = 0._r8
    uncoated_aer_num(1:ncol,1:pver,1:3) = 0._r8
    total_interstitial_aer_num(1:ncol,1:pver,1:3) = 0._r8
    total_cloudborne_aer_num(1:ncol,1:pver,1:3) = 0._r8
    awcam(1:ncol,1:pver,1:3) = 0._r8
    awfacm(1:ncol,1:pver,1:3) = 0._r8
    dstcoat(1:ncol,1:pver,1:3) = 0._r8
    !++ wy4.0
    na500(1:ncol,1:pver) = 0._r8
    tot_na500(1:ncol,1:pver) = 0._r8
    !-- wy4.0
    
#ifdef OSLO_AERO
   qaercwpt(1:ncol,1:pver,:) = 0.0_r8
       do m=1,nmodes_oslo
               do n=1,getNumberOfTracersInMode(m)
                       kk=getTracerIndex(m,n,.false.)! This gives the tracer index used in the q-array
                       qqcw(kk)%fldcw => qqcw_get_field(pbuf,kk,lchnk)
                       qaercwpt(:,:,kk) = qqcw(kk)%fldcw
               end do
       end do
#endif
!-- MH_2015/04/10

   
#ifndef OSLO_AERO
   if (clim_modal_aero) then
      ! mode number mixing ratios
      call rad_cnst_get_mode_num(0, mode_accum_idx,  'a', state, pbuf, num_accum)
      call rad_cnst_get_mode_num(0, mode_aitken_idx, 'a', state, pbuf, num_aitken)
      call rad_cnst_get_mode_num(0, mode_coarse_dst_idx, 'a', state, pbuf, num_coarse)

      ! mode specie mass m.r.
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'a', state, pbuf, coarse_dust)
      call rad_cnst_get_aer_mmr(0, mode_coarse_slt_idx, coarse_nacl_idx, 'a', state, pbuf, coarse_nacl)

   else
      ! init number/mass arrays for bulk aerosols
      allocate( &
         naer2(pcols,pver,naer_all), &
         maerosol(pcols,pver,naer_all))

      do m = 1, naer_all
         call rad_cnst_get_aer_mmr(0, m, state, pbuf, aer_mmr)
         maerosol(:ncol,:,m) = aer_mmr(:ncol,:)*rho(:ncol,:)

         if (m .eq. idxsul) then
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)*bulk_scale
         else
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)
         end if
      end do
   end if
#endif
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! More refined computation of sub-grid vertical velocity 
   ! Set to be zero at the surface by initialization.

   select case (trim(eddy_scheme))
   case ('diag_TKE')
      call pbuf_get_field(pbuf, tke_idx, tke)
   case ('CLUBB_SGS')
      itim = pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, wp2_idx, wp2, start=(/1,1,itim/),kount=(/pcols,pverp,1/))
      allocate(tke(pcols,pverp))
      tke(:ncol,:) = (3._r8/2._r8)*wp2(:ncol,:)
   case default
      call pbuf_get_field(pbuf, kvh_idx, kvh)
   end select

   ! Set minimum values above top_lev.
   wsub(:ncol,:top_lev-1)  = 0.20_r8
   wsubi(:ncol,:top_lev-1) = 0.001_r8
   
   !++ MH_2015/08/17
   if(preexisting_ice) then
   fhom(:,:) = 0.0_r8
   wice(:,:) = 0.0_r8
   weff(:,:) = 0.0_r8
   endif
   !-- MH_2015/08/17

   do k = top_lev, pver
      do i = 1, ncol

         select case (trim(eddy_scheme))
         case ('diag_TKE', 'CLUBB_SGS')
               wsub(i,k) = sqrt(0.5_r8*(tke(i,k) + tke(i,k+1))*(2._r8/3._r8))
               wsub(i,k) = min(wsub(i,k),10._r8)
         case default 
            ! get sub-grid vertical velocity from diff coef.
            ! following morrison et al. 2005, JAS
            ! assume mixing length of 30 m
               dum = (kvh(i,k) + kvh(i,k+1))/2._r8/30._r8
            ! use maximum sub-grid vertical vel of 10 m/s
               dum = min(dum, 10._r8)
            ! set wsub to value at current vertical level
               wsub(i,k)  = dum
     end select

         wsubi(i,k) = max(0.001_r8, wsub(i,k))
         !++ MH_2015/09/09
         if(.not. preexisting_ice) then
           wsubi(i,k) = min(wsubi(i,k), 0.2_r8)
         endif
         !-- MH_2015/09/09
     
#ifdef CLUBB_SGS
     if (wsubi(i,k) .le. 0.04_r8) then
           nucboast=100._r8
       wsubi(i,k)=nucboast*wsubi(i,k)  ! boost ice SGS vertical velocity in CAM-CLUBB
                       ! to force nucleation in upper-level stratiform 
                       ! clouds.  Temporary fix until cloud-top radiative
                       ! cooling parameterization is added to CLUBB similar
                       ! to the one of appendix C of Bretherton and Park (2009).  
     endif
#endif
     
         wsub(i,k)  = max(0.20_r8, wsub(i,k))
      end do
   end do
   call outfld( 'WSUB'       , wsub,      pcols, lchnk )
   call outfld( 'WSUBI'      , wsubi,     pcols, lchnk )

   if (trim(eddy_scheme) == 'CLUBB_SGS') deallocate(tke)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !Get humidity and saturation vapor pressures

   ! find wet bulk temperature and saturation value for provisional t and q without
   ! condensation

   do k = top_lev, pver

      call qsat_water(t(:ncol,k), pmid(:ncol,k), &
           es(:ncol), qs(:ncol), gam=gammas(:ncol))

      do i = 1, ncol

         relhum(i,k) = qn(i,k)/qs(i)

         ! get cloud fraction, check for minimum
         icldm(i,k) = max(icecldf(i,k), mincld)
         lcldm(i,k) = max(liqcldf(i,k), mincld)

         ! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)
         nfice(i,k) = 0._r8
         dumfice    = qc(i,k) + qi(i,k)
         if (dumfice > qsmall .and. qi(i,k) > qsmall) then
            nfice(i,k) = qi(i,k)/dumfice
         end if
      end do
   end do

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !ICE Nucleation

#ifdef OSLO_AERO
   call calculateNumberConcentration(ncol, qaercwpt, rho, CloudnumberConcentration)
#endif

!++ MH_2015/04/10
#ifdef OSLO_AERO

!Get size distributed interstitial aerosol
        call parmix_progncdnc_sub(        &   
                    ncol                  &        !I [nbr] number of columns used
                    ,state%q              &        !I [kg/kg] mass mixing ratio of tracers
                    ,rho                  &        !I [kg/m3] air density
                    ,CProcessModes        &        !O [kg/m3] added mass (total distributed all background modes)
                    ,f_c                  &        !O 
                    ,f_bc                 &        !O 
                    ,f_aq                 &        !O 
                    ,f_so4_cond           &        !O 
                    ,f_soa                &    
                    ,cam                  &        !O
                    ,f_acm             &        !O [frc] carbon fraction in mode
                    ,f_bcm             &        !O [frc] fraction of c being bc
                    ,f_aqm             &        !O [frc] fraction of sulfate being aquous
                    ,f_so4_condm           &    !O [frc] fraction of non-aquous SO4 being condensate
                    ,f_soam                &
                    ,numberConcentration &      !O [#/m3] number concentration
                    ,volumeConcentration &      !O [m3/m3] volume concentration
                    ,hygroscopicity    &        !O [mol/mol]
                    ,lnsigma           &        !O [-] log sigma 
                    ,hasAerosol        &        !I [t/f] do we have this type of aerosol here?
                    ,volumeCore        &
                    ,volumeCoat        &
                  )
#endif
!-- MH_2015/04/10

   do k = top_lev, pver
      do i = 1, ncol

         if (t(i,k).lt.tmelt - 5._r8) then


            ! compute aerosol number for so4, soot, and dust with units #/cm^3
            so4_num  = 0._r8
            soot_num = 0._r8
            dst1_num = 0._r8
            dst2_num = 0._r8
            dst3_num = 0._r8
            dst4_num = 0._r8
            dst_num  = 0._r8

            if (clim_modal_aero) then
#ifdef OSLO_AERO
               soot_num =  numberConcentration(i,k,MODE_IDX_OMBC_INTMIX_COAT_AIT)*1.0e-6_r8

               dst_num = (numberConcentration(i,k,MODE_IDX_DST_A2) &
                          + numberConcentration(i,k,MODE_IDX_DST_A3))*1.0e-6_r8

               so4_num = (numberConcentration(i,k,MODE_IDX_SO4_AC))*1.0e-6_r8
#else
               !For modal aerosols, assume for the upper troposphere:
               ! soot = accumulation mode
               ! sulfate = aiken mode
               ! dust = coarse mode
               ! since modal has internal mixtures.
               soot_num = num_accum(i,k)*rho(i,k)*1.0e-6_r8
               dmc  = coarse_dust(i,k)*rho(i,k)
               ssmc = coarse_nacl(i,k)*rho(i,k)

               if ( separate_dust ) then
                  ! 7-mode -- has separate dust and seasalt mode types and no need for weighting 
                  wght = 1._r8
               else
                  ! 3-mode -- needs weighting for dust since dust and seasalt are combined in the "coarse" mode type
                  wght = dmc/(ssmc + dmc)
               endif

               if (dmc > 0._r8) then
                  dst_num = wght * num_coarse(i,k)*rho(i,k)*1.0e-6_r8
               else 
                  dst_num = 0.0_r8
               end if

               if (dgnum(i,k,mode_aitken_idx) > 0._r8) then
                  ! only allow so4 with D>0.1 um in ice nucleation
                  !++ MH_2015/09/09
                  if(.not. preexisting_ice) then
                  so4_num  = num_aitken(i,k)*rho(i,k)*1.0e-6_r8 &
                     * (0.5_r8 - 0.5_r8*erf(log(0.1e-6_r8/dgnum(i,k,mode_aitken_idx))/  &
                     (2._r8**0.5_r8*log(sigmag_aitken))))
                  else
                  so4_num  = num_aitken(i,k)*rho(i,k)*1.0e-6_r8          ! sxj, all so4 from aitken
                  end if
                  !-- MH_2015/09/09
               else 
                  so4_num = 0.0_r8 
               end if
               so4_num = max(0.0_r8, so4_num)

            else

               if (idxsul > 0) then 
                  so4_num = naer2(i,k,idxsul)/25._r8 *1.0e-6_r8
               end if
               if (idxbcphi > 0) then 
                  soot_num = naer2(i,k,idxbcphi)/25._r8 *1.0e-6_r8
               end if
               if (idxdst1 > 0) then 
                  dst1_num = naer2(i,k,idxdst1)/25._r8 *1.0e-6_r8
               end if
               if (idxdst2 > 0) then 
                  dst2_num = naer2(i,k,idxdst2)/25._r8 *1.0e-6_r8
               end if
               if (idxdst3 > 0) then 
                  dst3_num = naer2(i,k,idxdst3)/25._r8 *1.0e-6_r8
               end if
               if (idxdst4 > 0) then 
                  dst4_num = naer2(i,k,idxdst4)/25._r8 *1.0e-6_r8
               end if
               dst_num = dst1_num + dst2_num + dst3_num + dst4_num

#endif 
            end if
            ! *** Turn off soot nucleation ***
            soot_num = 0.0_r8

!++ MH_2015/08/17
            if(preexisting_ice) then
            call nucleati( &
               wsubi(i,k), t(i,k), relhum(i,k), icldm(i,k), qc(i,k), &
               nfice(i,k), rho(i,k), so4_num, dst_num, soot_num,     &
               naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
               qi(i,k), ni(i,k), pmid(i,k),wice(i,k), weff(i,k), fhom(i,k))      ! sxj,   input: qi,ni,pmid;   output: wice,weff,fhom
            else
            call nucleati( &
               wsubi(i,k), t(i,k), relhum(i,k), icldm(i,k), qc(i,k), &
               nfice(i,k), rho(i,k), so4_num, dst_num, soot_num,     &
               naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k))
            end if
!-- MH_2015/08/17

            naai_hom(i,k) = nihf(i,k)

            ! output activated ice (convert from #/kg -> #/m3)
            nihf(i,k)     = nihf(i,k) *rho(i,k)
            niimm(i,k)    = niimm(i,k)*rho(i,k)
            nidep(i,k)    = nidep(i,k)*rho(i,k)
            nimey(i,k)    = nimey(i,k)*rho(i,k)
         end if
      end do
   end do

   call outfld('NIHF',   nihf, pcols, lchnk)
   call outfld('NIIMM', niimm, pcols, lchnk)
   call outfld('NIDEP', nidep, pcols, lchnk)
   call outfld('NIMEY', nimey, pcols, lchnk)
   !++ MH_2015/08/17
   if(preexisting_ice) then
   call outfld( 'fhom' , fhom, pcols, lchnk)
   call outfld( 'WICE' , wice, pcols, lchnk)
   call outfld( 'WEFF' , weff, pcols, lchnk)
   endif
   !-- MH_2015/08/17


   if (clim_modal_aero) then

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !droplet activation for modal aerosol

      ! partition cloud fraction into liquid water part
      lcldn = 0._r8
      lcldo = 0._r8
      do k = top_lev, pver
         do i = 1, ncol
            qcld = qc(i,k) + qi(i,k)
            if (qcld > qsmall) then
               lcldn(i,k) = cldn(i,k)*qc(i,k)/qcld
               lcldo(i,k) = cldo(i,k)*qc(i,k)/qcld
            end if
         end do
      end do

      call outfld('LCLOUD', lcldn, pcols, lchnk)

      call dropmixnuc( &
         state, ptend, deltatin, pbuf, wsub,         &  ! Input
         lcldn, lcldo,                               &
         !++ MH_2015/09/07
         hasAerosol, &
         CProcessModes, f_c, f_bc, f_aq, f_so4_cond, &
         f_soa,                                      &
         cam, f_acm, f_bcm, f_aqm, f_so4_condm,      &
         f_soam,                                     &
         numberConcentration, volumeConcentration,   &
         hygroscopicity, lnsigma,                    &
         !-- MH_2015/09/07
         nctend_mixnuc,                              &  ! Output
         !++ MH_2015/04/10
         fn )
         !-- MH_2015/04/10

      npccn(:ncol,:) = nctend_mixnuc(:ncol,:)

   else

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !droplet activation for bulk aerosol

      ! no tendencies returned from ndrop_bam_run, so just init ptend here
      call physics_ptend_init(ptend, state%psetcols, 'none')

      do k = top_lev, pver
         do i = 1, ncol

            if (qc(i,k) >= qsmall) then

               ! get droplet activation rate

               call ndrop_bam_run( &
                  wsub(i,k), t(i,k), rho(i,k), naer2(i,k,:), naer_all, &
                  naer_all, maerosol(i,k,:),  &
                  dum2)
               dum = dum2
            else
               dum = 0._r8
            end if

            ! note: deltatin/2.  accounts for sub step in microphysics
            ! ***** This assumes two sub-steps in microphysics.  It's dangerous to 
            ! ***** make that assumption here.  Should move all coding related to 
            ! ***** microphysics substepping into the microphysics.
            npccn(i,k) = (dum - nc(i,k)/lcldm(i,k))/(deltatin/2._r8)*lcldm(i,k)
         end do
      end do

   end if


   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Contact freezing  (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
   ! estimate rndst and nanco for 4 dust bins here to pass to MG microphysics
   
   do k = top_lev, pver
      do i = 1, ncol

         if (t(i,k) < 269.15_r8) then

            if (clim_modal_aero) then
#if OSLO_AERO
               !fxm: I think model uses bins, not modes.. But to get it 
               !approximately correct, use mode radius in first version
               nacon(i,k,2) = numberConcentration(i,k,MODE_IDX_DST_A2)
               nacon(i,k,3) = numberConcentration(i,k,MODE_IDX_DST_A3) 
               rndst(i,k,2) = lifeCycleNumberMedianRadius(MODE_IDX_DST_A2)
               rndst(i,k,3) = lifeCycleNumberMedianRadius(MODE_IDX_DST_A3)
               nacon(i,k,1) = 0.0_r8 !Set to zero to make sure
               nacon(i,k,4) = 0.0_r8 !Set to zero to make sure
#else

               ! For modal aerosols:
               !  use size '3' for dust coarse mode...
               !  scale by dust fraction in coarse mode
               
               dmc  = coarse_dust(i,k)
               ssmc = coarse_nacl(i,k)

               if ( separate_dust ) then
                  ! 7-mode -- has separate dust and seasalt mode types and no need for weighting 
                  wght = 1._r8
               else
                  ! 3-mode -- needs weighting for dust since dust and seasalt are combined in the "coarse" mode type
                  wght = dmc/(ssmc + dmc)
               endif

               if (dmc > 0.0_r8) then
                  nacon(i,k,3) = wght*num_coarse(i,k)*rho(i,k)
               else
                  nacon(i,k,3) = 0._r8
               end if

               !also redefine parameters based on size...

               rndst(i,k,3) = 0.5_r8*dgnumwet(i,k,mode_coarse_dst_idx)
               if (rndst(i,k,3) <= 0._r8) then 
                  rndst(i,k,3) = rn_dst3
               end if

#endif
            else

               !For Bulk Aerosols: set equal to aerosol number for dust for bins 2-4 (bin 1=0)

               if (idxdst2 > 0) then 
                  nacon(i,k,2) = naer2(i,k,idxdst2)
               end if
               if (idxdst3 > 0) then 
                  nacon(i,k,3) = naer2(i,k,idxdst3)
               end if
               if (idxdst4 > 0) then 
                  nacon(i,k,4) = naer2(i,k,idxdst4)
               end if
            end if

         end if
      end do
   end do

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !bulk aerosol ccn concentration (modal does it in ndrop, from dropmixnuc)

   if (.not. clim_modal_aero) then

      ! ccn concentration as diagnostic
      call ndrop_bam_ccn(lchnk, ncol, maerosol, naer2)

      deallocate( &
         naer2,    &
         maerosol)

   end if


!++ MH_2015/04/10
if(classnuc_in) then
! output aerosols as reference information for heterogeneous freezing
        do i = 1, ncol
            do k = 1, pver
!            call get_aer_num(state%q(i,k,:)*rho(i,k), qaercwpt(i,k,:)*rho(i,k), rho(i,k),         &
            call get_aer_num(numberConcentration(i,k,:), CloudnumberConcentration(i,k,:), rho(i,k),         &
                        !++ MH_2015/04/10
                        f_acm(i,k,:), f_so4_condm(i,k,:), cam(i,k,:), volumeCore(i,k,:), volumeCoat(i,k,:), &
                        !-- MH_2015/04/10
                        total_aer_num(i,k,:), coated_aer_num(i,k,:), uncoated_aer_num(i,k,:),  &
                        total_interstitial_aer_num(i,k,:), total_cloudborne_aer_num(i,k,:),    &
                        hetraer(i,k,:), awcam(i,k,:), awfacm(i,k,:), dstcoat(i,k,:),           &
                        !++ wy4.0
                        na500(i,k), tot_na500(i,k))
                        !-- wy4.0
            end do
        end do
        
        call outfld('bc_num', total_aer_num(:,:,1), pcols,lchnk)
        call outfld('dst1_num', total_aer_num(:,:,2), pcols,lchnk)
        call outfld('dst3_num', total_aer_num(:,:,3), pcols,lchnk)
        call outfld('bcc_num', coated_aer_num(:,:,1), pcols,lchnk)
        call outfld('dst1c_num', coated_aer_num(:,:,2), pcols,lchnk)
        call outfld('dst3c_num', coated_aer_num(:,:,3), pcols,lchnk)
        call outfld('bcuc_num', uncoated_aer_num(:,:,1), pcols,lchnk)
        call outfld('dst1uc_num', uncoated_aer_num(:,:,2), pcols,lchnk)
        call outfld('dst3uc_num', uncoated_aer_num(:,:,3), pcols,lchnk)
        call outfld('bc_a1_num', total_interstitial_aer_num(:,:,1), pcols,lchnk)
        call outfld('dst_a1_num', total_interstitial_aer_num(:,:,2), pcols,lchnk)
        call outfld('dst_a3_num', total_interstitial_aer_num(:,:,3), pcols,lchnk)
        call outfld('bc_c1_num', total_cloudborne_aer_num(:,:,1), pcols,lchnk)
        call outfld('dst_c1_num', total_cloudborne_aer_num(:,:,2), pcols,lchnk)
        call outfld('dst_c3_num', total_cloudborne_aer_num(:,:,3), pcols,lchnk)
        
        !++ wy4.0
        call outfld('na500', na500, pcols, lchnk)
        call outfld('totna500', tot_na500, pcols, lchnk)
        !-- wy4.0   

    call pbuf_get_field(pbuf,frzimm_idx,frzimm)
    call pbuf_get_field(pbuf,frzcnt_idx,frzcnt)
    call pbuf_get_field(pbuf,frzdep_idx,frzdep)
    
    frzbcimm(1:ncol,1:pver) = 0._r8
    frzduimm(1:ncol,1:pver) = 0._r8
    frzbccnt(1:ncol,1:pver) = 0._r8
    frzducnt(1:ncol,1:pver) = 0._r8
    frzbcdep(1:ncol,1:pver) = 0._r8
    frzdudep(1:ncol,1:pver) = 0._r8

    freqimm(1:ncol,1:pver) = 0._r8
    freqcnt(1:ncol,1:pver) = 0._r8
    freqdep(1:ncol,1:pver) = 0._r8
    freqmix(1:ncol,1:pver) = 0._r8

    numice10s(1:ncol,1:pver) = 0._r8
    numice10s_imm_dst(1:ncol,1:pver) = 0._r8
    numice10s_imm_bc(1:ncol,1:pver) = 0._r8
    
    do i = 1,ncol
    do k = 1,pver
     if(t(i,k).gt.235.15_r8) then
       qcic = min(qc(i,k)/lcldm(i,k),5.e-3_r8)
       ncic = max(nc(i,k)/lcldm(i,k),0._r8)

       con1 = 1._r8/(1.333_r8*pi)**0.333_r8
       r3lx = con1*(rho(i,k)*qcic/(rhow*max(ncic*rho(i,k), 1.0e6_r8)))**0.333 ! in m
       r3lx = max(4.e-6, r3lx)
       mi0l = 4._r8/3._r8*pi*rhow*r3lx**3_r8
       supersatice = svp_water(t(i,k))/svp_ice(t(i,k))

       fn_in(1) = fn(i,k,MODE_IDX_SO4_AC)  ! bc accumulation mode
       fn_in(2) = fn(i,k,MODE_IDX_DST_A2)  ! dust_a1 accumulation mode
       fn_in(3) = fn(i,k,MODE_IDX_DST_A3) ! dust_a3 coarse mode

       call hetfrz_classnuc(deltatin, t(i,k), state%pmid(i,k), supersatice,               &
                            fn_in(:),                                          &
                            r3lx, ncic*rho(i,k)*1.0e-6_r8,                     &
                            frzbcimm(i,k), frzduimm(i,k),                               &
                            frzbccnt(i,k), frzducnt(i,k),                               &
                            frzbcdep(i,k), frzdudep(i,k),                               &
                            hetraer(i,k,:), awcam(i,k,:), awfacm(i,k,:), dstcoat(i,k,:),                  &
                            total_aer_num(i,k,:), coated_aer_num(i,k,:), uncoated_aer_num(i,k,:),  &
                            total_interstitial_aer_num(i,k,:), total_cloudborne_aer_num(i,k,:))

       frzimm(i,k) = frzbcimm(i,k) + frzduimm(i,k)
       frzcnt(i,k) = frzbccnt(i,k) + frzducnt(i,k)
       frzdep(i,k) = frzbcdep(i,k) + frzdudep(i,k)
       
       if( frzimm(i,k) .gt. 0._r8 ) freqimm(i,k) = 1._r8
       if( frzcnt(i,k) .gt. 0._r8 ) freqcnt(i,k) = 1._r8
       if( frzdep(i,k) .gt. 0._r8 ) freqdep(i,k) = 1._r8
       if( (frzimm(i,k)+frzcnt(i,k)+frzdep(i,k)) .gt. 0.) freqmix(i,k) = 1._r8
     else
       frzimm(i,k) = 0._r8
       frzcnt(i,k) = 0._r8
       frzdep(i,k) = 0._r8
     end if
     
     nnuccc_bc(i,k) = frzbcimm(i,k)*1.0e6_r8*cldn(i,k)
     nnucct_bc(i,k) = frzbccnt(i,k)*1.0e6_r8*cldn(i,k)
     nnudep_bc(i,k) = frzbcdep(i,k)*1.0e6_r8*cldn(i,k)
     nnuccc_dst(i,k) = frzduimm(i,k)*1.0e6_r8*cldn(i,k)
     nnucct_dst(i,k) = frzducnt(i,k)*1.0e6_r8*cldn(i,k)     
     nnudep_dst(i,k) = frzdudep(i,k)*1.0e6_r8*cldn(i,k)
     niimm_bc(i,k) = frzbcimm(i,k)*1.0e6_r8*deltatin
     nicnt_bc(i,k) = frzbccnt(i,k)*1.0e6_r8*deltatin
     nidep_bc(i,k) = frzbcdep(i,k)*1.0e6_r8*deltatin
     niimm_dst(i,k) = frzduimm(i,k)*1.0e6_r8*deltatin
     nicnt_dst(i,k) = frzducnt(i,k)*1.0e6_r8*deltatin
     nidep_dst(i,k) = frzdudep(i,k)*1.0e6_r8*deltatin
     fn_bc(i,k) = fn(i,k,MODE_IDX_SO4_AC)
     fn_dst1(i,k) = fn(i,k,MODE_IDX_DST_A2)
     fn_dst3(i,k) = fn(i,k,MODE_IDX_DST_A3)
     hetraer_bc(i,k) = hetraer(i,k,1)
     hetraer_dst1(i,k) = hetraer(i,k,2)
     hetraer_dst3(i,k) = hetraer(i,k,3)
     dstcoat_bc(i,k) = dstcoat(i,k,1)
     dstcoat_dst1(i,k) = dstcoat(i,k,2)
     dstcoat_dst3(i,k) = dstcoat(i,k,3)
     
     numice10s(i,k) = (frzimm(i,k)+frzcnt(i,k)+frzdep(i,k))*1.0e6_r8*deltatin*(10._r8/deltatin)
     numice10s_imm_dst(i,k) = frzduimm(i,k)*1.0e6_r8*deltatin*(10._r8/deltatin)
     numice10s_imm_bc(i,k) = frzbcimm(i,k)*1.0e6_r8*deltatin*(10._r8/deltatin)
    end do
    end do

   call outfld('FREQIMM', freqimm, pcols, lchnk)
   call outfld('FREQCNT', freqcnt, pcols, lchnk)
   call outfld('FREQDEP', freqdep, pcols, lchnk)
   call outfld('FREQMIX', freqmix, pcols, lchnk)

   call outfld('DSTFREZIMM', nnuccc_dst, pcols, lchnk)
   call outfld('DSTFREZCNT', nnucct_dst, pcols, lchnk)
   call outfld('DSTFREZDEP', nnudep_dst, pcols, lchnk)

   call outfld('BCFREZIMM', nnuccc_bc, pcols, lchnk)
   call outfld('BCFREZCNT', nnucct_bc, pcols, lchnk)
   call outfld('BCFREZDEP', nnudep_bc, pcols, lchnk)

!  call outfld('NIMIX', nimix, pcols, lchnk)
   call outfld('NIMIX_IMM', niimm_bc+niimm_dst, pcols, lchnk)
   call outfld('NIMIX_CNT', nicnt_bc+nicnt_dst, pcols, lchnk)   
   call outfld('NIMIX_DEP', nidep_bc+nidep_dst, pcols, lchnk)

   call outfld('DSTNICNT', nicnt_dst, pcols, lchnk)
   call outfld('DSTNIDEP', nidep_dst, pcols, lchnk)
   call outfld('DSTNIIMM', niimm_dst, pcols, lchnk)

   call outfld('BCNICNT', nicnt_bc, pcols, lchnk)
   call outfld('BCNIDEP', nidep_bc, pcols, lchnk)
   call outfld('BCNIIMM', niimm_bc, pcols, lchnk)

   call outfld('NUMICE10s', numice10s, pcols, lchnk)
   call outfld('NUMIMM10sDST', numice10s_imm_dst, pcols, lchnk)
   call outfld('NUMIMM10sBC', numice10s_imm_bc, pcols, lchnk)
   
   call outfld('fn_bc_c1_num', fn_bc, pcols, lchnk)
   call outfld('fn_dst_c1_num', fn_dst1, pcols, lchnk)
   call outfld('fn_dst_c3_num', fn_dst3, pcols, lchnk)
 end if
!-- MH_2015/04/10

end subroutine microp_aero_run

!===============================================================================


!===============================================================================

end module microp_aero

