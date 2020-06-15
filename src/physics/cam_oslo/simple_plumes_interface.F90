
SUBROUTINE simple_plumes_interface (lchnk,ncol, nswbands,surf_geopot, &
!ak+                                    deltah_km, cdnc, clon, clat, year_fr, &
                                    deltah_km, clon, clat, year_fr, &
                                    sp_tau, sp_ssa, sp_asy, re_mult, xcdnc)
!ak-

!---------------------------------------------------------------------------
! Purpose:
! - interface between the MAC-SP ("simple plumes") aerosol climatology
!   and CAM/NorESM radiation scheme
!    * called from "radiation_tend" in radiation.F90
!    * calls "sp_aop_profile" in mo_simple_plumes.F90
!
! May 2016, P. Räisänen : initial version
! Sept 2019, A. Kirkevåg : modified for CAM6-Oslo / NorESM2 (+ extra output)
!---------------------------------------------------------------------------

#include <preprocessorDefinitions.h>

  USE ppgrid,          only: pcols, pver
  USE shr_kind_mod,    only: r8=>shr_kind_r8
  USE physconst,        ONLY: pi, rga
#ifdef SPAERO
  USE radconstants,     ONLY: wav_sp
#endif
  USE CAM_HISTORY,      ONLY: outfld
  USE mo_simple_plumes, ONLY: sp_aop_profile
!ak+
  use cam_logfile,      only: iulog
!ak-

  IMPLICIT NONE

!ak (NorESM1-M):  INTEGER, PARAMETER :: nextra = 10 ! Number of layers added between the 
  INTEGER, PARAMETER :: nextra = 13 ! Number of layers added between the
                                    ! sea level and model surface for computing
                                    ! the simple plumes aerosols 

! Arguments
  INTEGER, INTENT(in) :: lchnk                   ! chunk identifier
  INTEGER, INTENT(in) :: ncol ! number of columns
  INTEGER, INTENT(in) :: nswbands ! number of shortwave bands
    
  REAL(r8), INTENT(in) :: surf_geopot(pcols)  ! surface geopotential
  REAL(r8), INTENT(in) :: deltah_km(pcols,pver) ! layer thicknesses  (1=uppermost layer)
!ak  REAL(r8), INTENT(in) :: cdnc(pcols,pver)    ! in-cloud CDNC [cm-3?] (1=uppermost layer)
  REAL(r8), INTENT(in) :: clon(pcols)  ! longitude (in radians)
  REAL(r8), INTENT(in) :: clat(pcols)  ! latitude (in radians)
  REAL(r8), INTENT(in) :: year_fr ! Fractional year (1903.0 is the 0Z on 
                                  !the first of January 1903, Gregorian)

  REAL(r8), INTENT(out) :: sp_tau (pcols,pver,nswbands) ! aerosol extinction optical depth (1=uppermost layer)
  REAL(r8), INTENT(out) :: sp_ssa (pcols,pver,nswbands) ! aerosol single scattering albedo (1=uppermost layer)
  REAL(r8), INTENT(out) :: sp_asy (pcols,pver,nswbands) ! aerosol assymetry parameter (1=uppermost layer)
  REAL(r8), INTENT(out) :: re_mult(pcols,pver) ! Multiplication factor of liquid cloud effective radius 
  REAL(r8), INTENT(out) :: xcdnc(pcols)  ! CDNC modification factor

! Local variables: input for SP_AOP_PROFILE
 
  REAL(r8) :: col_lon(pcols)       ! longitudes in degrees
  REAL(r8) :: col_lat(pcols)       ! latitudes in degrees
  REAL(r8) :: oro(pcols)           ! surface orography (height above sea level in [m])
  REAL(r8) :: z_sp(pcols,pver+nextra)  ! layer mid-point height above sea level [m], (1=lowermost layer)    
  REAL(r8) :: dz_sp(pcols,pver+nextra) ! layer thickness [m] (1=lowermost layer)   

  REAL(r8) :: lambda  ! wavelength in [nm]

! Local variables: output from SP_AOP_PROFILE

  REAL(r8) :: aod_prof(pcols,pver+nextra) ! aerosol optical depth (1=lowermost layer)
  REAL(r8) :: ssa_prof(pcols,pver+nextra) ! aerosol single-scattering albedo (1=lowermost layer)
  REAL(r8) :: asy_prof(pcols,pver+nextra) ! aerosol asymmetry parameter (1=lowermost layer)

! Local variables: fields written to model output. 
  REAL(r8) :: aodvis_sp(pcols) ! AOD at 0.35-0.64 µm
  REAL(r8) :: absvis_sp(pcols) ! Absorption AOD at 0.35-0.64 µm
!ak+
  REAL(r8) :: aodv3d_sp(pcols,pver) ! 3D AOD at 0.35-0.64 µm
  REAL(r8) :: absv3d_sp(pcols,pver) ! 3D absorption AOD at 0.35-0.64 µm
!ak-  
! Other local variables
  REAL(r8) :: zhalf(pcols,pver+1)       ! layer interface height above sea level [m] (1=uppermost)
  REAL(r8) :: zhalf_sp(pcols,pver+nextra+1) ! layer interface height above sea level [m]
                                         ! grid used for simple plume climatology (1=surface)        

  REAL(r8), PARAMETER :: rad2deg = 180._r8 / pi

!ak+  REAL(r8) :: zmin, zmax, deltaz, rv_mult, cdnc_sp, epsilon, epsilon_sp, &
!ak              beta, beta_sp, cdnc_test
  REAL(r8) :: zmin, zmax, deltaz
!ak-

  INTEGER :: i,j,k,ns

!****************************************
! Prepare input data for SP_AOP_PROFILE
!****************************************

! Longitude and latitude in degrees
  col_lon(:) = rad2deg * clon(:)
  col_lat(:) = rad2deg * clat(:)
! Surface elevation 
  oro(:) = surf_geopot(:) * rga

!test
   do i=1,pcols
!    write(*,*) 'i, surf_geopot(i) = ', i,  surf_geopot(i)
     if(surf_geopot(i).gt.100000._r8) oro(i)=0._r8
   enddo
!test

! Define the vertical grids
! 1) Model half-levels, layers indices decreasing upwards

  zhalf(:,pver+1) = oro(:)
  DO k=pver,1,-1
    zhalf(:,k) = zhalf(:,k+1) +  1000._r8 * deltah_km(:,k)  
  ENDDO 

! 2) Layer mid-point and half-level elevations for calling the simple plume
! climatology. A separate grid mist be used for this, since the grid must 
! start from the sea level (0 m), not from the model orographic height. 
! Otherwise, elevated terrain does not result in reduced total AOD, as it 
! should.
!
! To accomplish this, NEXTRA equally thick layers are added between the 
! sea-level and the model orographic height. The numbering of the simple
! plume layers goes upwards, lowest level = 1. 

  zhalf_sp(:,1) = 0._r8 
  zhalf_sp(:,nextra+1) = MAX(zhalf(:,pver+1),0.001_r8)  

  DO j=2,nextra
    zhalf_sp(:,j) = REAL(j-1)/REAL(nextra) * zhalf_sp(:,nextra+1)
  ENDDO

! Above the model surface, use the model grid for computing the simple
! plumes aerosols. Indexing:
! - lowest model half-level plev+1 = half-level nextra+1 for simple plumes
! - uppermost model half-level 1 = half-level pver+nextra+1 for simple plumes

  DO j=nextra+2,pver+nextra+1
    k = pver+nextra+2-j
    zhalf_sp(:,j) = zhalf(:,k)
  ENDDO
 
  DO j=1,pver+nextra  
    z_sp(:,j) = 0.5_r8*(zhalf_sp(:,j)+zhalf_sp(:,j+1))  
    dz_sp(:,j) = zhalf_sp(:,j+1)-zhalf_sp(:,j)
  ENDDO

!***********************************************************************
! Loop over spectral bands (SW only)

  DO ns=1,nswbands
! Wavelength in nanometres
!ak    lambda = 1000._r8*wav_sp(ns)   ! wav_sp is now a wave number (cm-1) 
#ifdef SPAERO
    lambda = 1.e7_r8/wav_sp(ns)
#else   ! Hack to make the model compile without SPAERO defined
    lambda = 1._r8
#endif

! Calculate aerosol profiles from the "simple plumes" climatology

    CALL sp_aop_profile(pver+nextra, ncol, lambda, oro, col_lon, col_lat, &
                        year_fr, z_sp, dz_sp, &
                        xcdnc, aod_prof, ssa_prof, asy_prof)

! Mapping to model vertical grid
!   model layer 1    <=> simple plumes layer pver+nextra
!   model layer pver <=> simple plumes layer nextra+1
     
    DO i=1,ncol
      DO k=nextra+1,nextra+pver 
        j = pver+nextra+1-k
        sp_tau(i,j,ns) = aod_prof(i,k)
        sp_ssa(i,j,ns) = ssa_prof(i,k)
        sp_asy(i,j,ns) = asy_prof(i,k)
! Security settings (are these needed?)
        sp_ssa(i,j,ns) = MIN(MAX(sp_ssa(i,j,ns),1.E-6_r8),0.999999_r8)
        sp_asy(i,j,ns) = MIN(MAX(sp_asy(i,j,ns),0._r8),0.99_r8)
      ENDDO
    ENDDO
  ENDDO  !1,nswbands 

! Conversion from CDNC modification factor to effective radius multiplier
! This assumes that XCDNC = total CDNC / CDNC without simple-plume aerosols 

  DO i=1,ncol 
    xcdnc(i) = MIN(MAX(xcdnc(i),1._r8),100._r8)  
! First, multiplication factor for volume-mean radius (r_v is proportional to CDNC**(-1/3))
!ak+    rv_mult = xcdnc(i)**(-0.333333_r8)
!ak-
! Include dispersion effect following Rotstayn & Liu (GRL 2009, 36, L10801),
!  eq. 2 (OLDBETA), similar to "cldwat.f90". 
!    beta = effective radius / volume-mean radius

    DO j=1,pver
!ak      cdnc_sp = cdnc(i,j) * xcdnc(i)
!ak      epsilon = 1._r8-0.7_r8*EXP(-0.003*cdnc(i,j))
!ak      epsilon_sp = 1._r8-0.7_r8*EXP(-0.003*cdnc_sp)
!ak      beta =    (1._r8+2._r8*epsilon**2   )**0.666667_r8 / (1._r8+epsilon**2   )**0.333333_r8
!ak      beta_sp = (1._r8+2._r8*epsilon_sp**2)**0.666667_r8 / (1._r8+epsilon_sp**2)**0.333333_r8
!
! Multiplication factor for effective radius, including the dispersion effect 
!     re_mult(i,j) = beta_sp/beta * rv_mult 
! For BACCHUS, without the dispersion effect (P. Räisänen, 14 March 2017)
!ak      re_mult(i,:) = xcdnc(i)**(-0.333333_r8)  
      re_mult(i,j) = xcdnc(i)**(-0.333333_r8)  
    ENDDO
  ENDDO 
 
! Write fields to model output

  aodvis_sp(:) = 0._r8
  absvis_sp(:) = 0._r8
!ak+
  aodv3d_sp(:,:) = 0._r8
  absv3d_sp(:,:) = 0._r8
!ak-
    
  DO j=1,pver    
!ak    aodvis_sp(:) = aodvis_sp(:)+sp_tau(:,j,8)       
!ak    absvis_sp(:) = absvis_sp(:)+sp_tau(:,j,8)*(1._r8-sp_ssa(:,j,8))       
    aodv3d_sp(:,j) = sp_tau(:,j,10)       
    absv3d_sp(:,j) = sp_tau(:,j,10)*(1._r8-sp_ssa(:,j,10))       
    aodvis_sp(:) = aodvis_sp(:)+aodv3d_sp(:,j)
    absvis_sp(:) = absvis_sp(:)+absv3d_sp(:,j)
  ENDDO

!ak+
  call outfld('AODVISSP',aodvis_sp ,pcols,lchnk)
  call outfld('ABSVISSP',absvis_sp ,pcols,lchnk)
  call outfld('XCDNC_SP',xcdnc     ,pcols,lchnk)
  call outfld('AODV3DSP',aodv3d_sp ,pcols,lchnk)
  call outfld('ABSV3DSP',absv3d_sp ,pcols,lchnk)
!ak-

  RETURN
END SUBROUTINE simple_plumes_interface
