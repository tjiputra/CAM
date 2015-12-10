module oslo_dust_intr
!---------------------------------------------------------------------------------

use aerosoldef,       only: l_dst_a2, l_dst_a3
use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
use camsrfexch,       only: cam_in_t
use ppgrid,           only: pcols
use abortutils,       only: endrun
use cam_logfile,      only: iulog

implicit none
private
save

   integer, parameter :: numberOfDustModes = 2  !define in aerosoldef?

   !This can be refined, but the fractions in coarse/fine mode are approx ok
   real(r8), parameter, dimension(numberOfDustModes) :: emis_fraction_in_mode = (/0.13_r8, 0.87_r8 /)
   integer, dimension(numberOfDustModes)             :: tracerMap = (/-99, -99/) !index of dust tracers in the modes

   !Related to soil erodibility
   real(r8), allocatable ::  soil_erodibility(:,:)     ! soil erodibility factor
   real(r8), allocatable ::  soil_erodibility_in(:,:)  ! temporary input array

   real(r8)      :: dust_emis_fact = -1.e36_r8   ! tuning parameter for dust emissions
   character(cl) :: soil_erod = 'soil_erod'   ! full pathname for soil erodibility dataset

public oslo_dust_emis_intr
public getNumberOfDustModes
public getDustTracerIndexInMode
public getEmissionFractionInDustMode
public isOsloDustTracer
public oslo_dust_initialize


!===============================================================================
contains
!===============================================================================

   function getEmissionFractionInDustMode(modeIndex) RESULT(fraction)
      integer, intent(in) :: modeIndex
      real(r8)            :: fraction
      fraction = emis_fraction_in_mode(modeIndex)
   end function getEmissionFractionInDustMode

   function getNumberOfDustModes() RESULT(answer)
      integer answer
      answer = numberOfDustModes
   end function getNumberOfDustModes


   subroutine oslo_dust_initialize(dust_emis_fact_in, soil_erod_in)
      implicit none

      character(cl), intent(in) :: soil_erod_in
      real(r8), intent(in)      :: dust_emis_fact_in

      !Nail the overall tuning factor
      dust_emis_fact = dust_emis_fact_in

      call read_soil_erodibility_data(soil_erod_in)

      call set_oslo_indices()

   end subroutine oslo_dust_initialize

   subroutine set_oslo_indices()
      implicit none
      tracerMap(1) = l_dst_a2
      tracerMap(2) = l_dst_a3
   end subroutine set_oslo_indices


   !****************************************************
   !This is copied from the MAM aerosols. Should not really
   !be necessary since the land model could calculate emissions
   !based on soil erodibility. 
   
   !However, the following code in dustMod.F90 (land model) makes it 
   !necessary to apply it here!
   !715 Set basin factor to 1 for now
   !716 
   !717     call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
   !718     do c = begc, endc
   !719       l = clm3%g%l%c%landunit(c)
   !720       if (.not. clm3%g%l%lakpoi(l)) then
   !721          mbl_bsn_fct(c) = 1.0_r8
   !722       end if
   !723     end do

   !For a general discussion of these factors, see: 
   !Zender et al JGR (vol 108, D16, 2003) 
   !http://onlinelibrary.wiley.com/doi/10.1029/2002JD003039/abstract

   subroutine read_soil_erodibility_data(soil_erod_in)

      use cam_history,      only: addfld, add_default, phys_decomp
      use ioFileMod,        only: getfil
      use pio
      use cam_pio_utils,    only : cam_pio_openfile
      use phys_grid,        only : get_ncols_p, get_rlat_all_p, get_rlon_all_p
      use interpolate_data, only : lininterp_init, lininterp, lininterp_finish, interp_type
      use ppgrid,           only : begchunk, endchunk
      use mo_constants,     only : pi, d2r

      implicit none

      character(cl), intent(in) :: soil_erod_in      !filename of soil erodibility data set
      character(cl)      :: soil_erod_file

      integer :: did, vid, nlat, nlon
      type(file_desc_t) :: ncid

      type(interp_type) :: lon_wgts, lat_wgts
      real(r8) :: to_lats(pcols), to_lons(pcols)
      real(r8), allocatable :: dst_lons(:)
      real(r8), allocatable :: dst_lats(:)
      integer :: c, ncols, ierr
      real(r8), parameter :: zero=0._r8, twopi=2._r8*pi

      !-----------------------------------------------------------------------

      ! set module data from namelist vars read in aerosol_intr module
      soil_erod      = soil_erod_in

      ! for soil erodibility in mobilization, apply inside CAM instead of lsm.
      ! read in soil erodibility factors, similar to Zender's boundary conditions

      ! Get file name.  
      call getfil(soil_erod, soil_erod_file, 0)
      call cam_pio_openfile (ncid, trim(soil_erod_file), PIO_NOWRITE)

      ! Get input data resolution.
      ierr = pio_inq_dimid( ncid, 'lon', did )
      ierr = pio_inq_dimlen( ncid, did, nlon )

      ierr = pio_inq_dimid( ncid, 'lat', did )
      ierr = pio_inq_dimlen( ncid, did, nlat )

      allocate(dst_lons(nlon))
      allocate(dst_lats(nlat))
      allocate(soil_erodibility_in(nlon,nlat))

      ierr = pio_inq_varid( ncid, 'lon', vid )
      ierr = pio_get_var( ncid, vid, dst_lons  )

      ierr = pio_inq_varid( ncid, 'lat', vid )
      ierr = pio_get_var( ncid, vid, dst_lats  )

      ierr = pio_inq_varid( ncid, 'mbl_bsn_fct_geo', vid )
      ierr = pio_get_var( ncid, vid, soil_erodibility_in )

   !-----------------------------------------------------------------------
   !     	... convert to radians and setup regridding
   !-----------------------------------------------------------------------
       dst_lats(:) = d2r * dst_lats(:)
       dst_lons(:) = d2r * dst_lons(:)

       allocate( soil_erodibility(pcols,begchunk:endchunk) )

   !-----------------------------------------------------------------------
   !     	... regrid ..
   !-----------------------------------------------------------------------
       do c=begchunk,endchunk
          ncols = get_ncols_p(c)
          call get_rlat_all_p(c, pcols, to_lats)
          call get_rlon_all_p(c, pcols, to_lons)
          
          call lininterp_init(dst_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
          call lininterp_init(dst_lats, nlat, to_lats, ncols, 1, lat_wgts)

          call lininterp(soil_erodibility_in(:,:), nlon,nlat , soil_erodibility(:,c), ncols, lon_wgts,lat_wgts)

          call lininterp_finish(lat_wgts)
          call lininterp_finish(lon_wgts)
       end do
       deallocate( soil_erodibility_in, stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'dust_initialize: failed to deallocate soil_erod_in, ierr = ',ierr
          call endrun
       end if

      call addfld('MBL_BSN_FCT','frac ',1, 'A','Soil erodibility factor',phys_decomp)
      call addfld('OSLO_DUST_EMIS','kg/m2/s ',1, 'A','oslo dust emissions',phys_decomp)
      

   end subroutine read_soil_erodibility_data

   function getDustTracerIndexInMode(modeIndex)RESULT(answer)
      integer, intent(in) :: modeIndex
      integer answer

      answer = tracerMap(modeIndex)
   
   end function getDustTracerIndexInMode

   function isOsloDustTracer(physTracerIndex) RESULT(answer)
      implicit none
      integer, intent(in) :: physTracerIndex
      integer             :: n
      logical             :: answer
      answer = .FALSE.
      do n = 1, numberOfDustModes
         if(tracerMap(n) .eq. physTracerIndex)then
            answer = .TRUE.
         end if
      end do
   end function isOsloDustTracer

   subroutine oslo_dust_emis_intr(state, cam_in)

      !----------------------------------------------------------------------- 
      ! 
      ! Purpose: 
      ! Interface to emission of all dusts.
      ! Notice that the mobilization is calculated in the land model (need #define BGC) and
      ! the soil erodibility factor is applied here.
      ! 
      ! see comments above in subroutine read_soil_erodibility_data 
      !-----------------------------------------------------------------------
      use cam_history,   only: outfld
      use physics_types, only: physics_state

      ! Arguments:

      type(physics_state),    intent(in)    :: state   ! Physics state variables
      type(cam_in_t), target, intent(inout) :: cam_in  ! import state

      ! Local variables

      integer :: lchnk
      integer :: ncol
      integer :: i,n
      real(r8) :: soil_erod_tmp(pcols)
      real(r8) :: totalEmissionFlux(pcols)
      real(r8), pointer :: cflx(:,:)

      lchnk = state%lchnk
      ncol = state%ncol

      !Filter away unreasonable values for soil erodibility
      !(using low values e.g. gives emissions in greenland..)
      where(soil_erodibility(:,lchnk) .lt. 0.1_r8)
         soil_erod_tmp(:)=0.0_r8
      elsewhere
         soil_erod_tmp(:)=soil_erodibility(:,lchnk)
      end where

      cflx => cam_in%cflx
   
      !Note that following CESM use of "dust_emis_fact", the emissions are 
      !scaled by the INVERSE of the factor!!
      !There is another random scale factor of 1.15 there. Adapting the exact
      !same formulation as MAM now and tune later
      !As of NE-380: Oslo dust emissions are 2/3 of CAM emissions
      do n=1, numberOfDustModes
!cak         cflx(:ncol, tracerMap(n)) = cflx(:ncol,tracerMap(n))*soil_erod_tmp(:ncol)/(dust_emis_fact)*1.15_r8*(2.0_r8) !factor 2 is CAM-Oslo specific
         cflx(:ncol, tracerMap(n)) = cflx(:ncol,tracerMap(n))*soil_erod_tmp(:ncol)/(dust_emis_fact)*1.15_r8  ! gives better AOD close to dust sources
      end do
     
      totalEmissionFlux(:)=0.0_r8
      do n=1,numberOfDustModes
         totalEmissionFlux(:ncol) = totalEmissionFlux(:ncol) + cflx(:ncol,tracerMap(n))
      end do

      call outfld('MBL_BSN_FCT',soil_erod_tmp,pcols,lchnk) 
      call outfld('OSLO_DUST_EMIS',totalEmissionFlux,pcols,lchnk)

      return
   end subroutine oslo_dust_emis_intr

end module oslo_dust_intr
