!*****************************************************************************
! Purpose: Calculate BC and Dust number, including total number(interstitial+
!          cloud borne), one monolayer coated number, and uncoated number
!
! Public interfaces:
!
! get_aer_num
!
!Subroutine will be called in ndrop (after parmix progncdnd_sub), 
!so it needs additional input 
!ndrop is called in mircop_aero and needs further output which is
!provided by get_aer_num
!*****************************************************************************


subroutine get_aer_num(qaerpt, qaercwpt, rhoair,           &   ! input
                       f_acm, f_condm,                     &
                       cam, volumeCore, volumeCoat,        &
                       total_aer_num,                      &   ! output
                       coated_aer_num,                     &
                       uncoated_aer_num,                   &
                       total_interstial_aer_num,           &
                       total_cloudborne_aer_num,           &
                       hetraer, awcam, awfacm, dstcoat,    &
!++ wy4.0
                       na500, tot_na500)
!-- wy4.0

    use spmd_utils, only: iam
    use shr_kind_mod,  only: r8 => shr_kind_r8
!    use ppgrid, only : pcols, pver
    use constituents,  only: pcnst
    use commondefinitions, only:  nmodes_oslo => nmodes
    use aerosoldef,    only:MODE_IDX_DST_A2, MODE_IDX_DST_A3, &
                            l_dst_a2, l_dst_a3, l_bc_ai,      &
                            MODE_IDX_OMBC_INTMIX_COAT_AIT, l_bc_ac,         &
                            lifeCycleNumberMedianRadius,      &
                            lifeCycleSigma


    implicit none

    ! input
    real(r8), intent(in) :: qaerpt(0:nmodes_oslo)           ! aerosol number and mass mixing ratios(instertitial)
    real(r8), intent(in) :: qaercwpt(0:nmodes_oslo)         ! cloud borne aerosol number and mass mixing ratios
    real(r8), intent(in) :: rhoair                  ! air density (kg/m3)
    real(r8), intent(in) :: f_acm(nmodes_oslo)
    real(r8), intent(in) :: f_condm(nmodes_oslo)
    real(r8), intent(in) :: cam(nmodes_oslo)
    real(r8), intent(in) :: volumeCoat(nmodes_oslo)
    real(r8), intent(in) :: volumeCore(nmodes_oslo)
    real(r8) :: sigmag_amode(3)
    
    
    ! output
    real(r8), intent(out) :: total_aer_num(3)            ! #/cm^3
    real(r8), intent(out) :: total_interstial_aer_num(3) ! #/cm^3
    real(r8), intent(out) :: total_cloudborne_aer_num(3) ! #/cm^3
    real(r8), intent(out) :: coated_aer_num(3)           ! #/cm^3 
    real(r8), intent(out) :: uncoated_aer_num(3)         ! #/cm^3
    real(r8), intent(out) :: hetraer(3)                  ! BC and Dust mass mean radius [m]
    real(r8), intent(out) :: awcam(3)                    ! modal added mass [mug m-3]
    real(r8), intent(out) :: awfacm(3)                   ! (OC+BC)/(OC+BC+SO4)
    real(r8), intent(out) :: dstcoat(3)                  ! coated fraction
    real(r8), intent(out) :: na500                       ! #/cm^3 interstitial aerosol number with D>500 nm (#/cm^3)
    real(r8), intent(out) :: tot_na500                   ! #/cm^3 total aerosol number with D>500 nm (#/cm^3)
    !local variables
    !------------coated variables--------------------
    real(r8), parameter :: n_so4_monolayers_dust = 1.0_r8 ! number of so4(+nh4) monolayers needed to coat a dust particle
    real(r8), parameter :: dr_so4_monolayers_dust = n_so4_monolayers_dust * 4.76e-10
    real(r8) :: tmp1, tmp2
    
    real(r8) :: bc_num                                    ! bc number in accumulation mode
    real(r8) :: dst1_num, dst3_num                        ! dust number in accumulation and corase mode
    real(r8) :: as_acm, as_condm
    real(r8) :: dst1_num_imm, dst3_num_imm, bc_num_imm
    real(r8) :: fac_volsfc_bc, fac_volsfc_dust_a1, fac_volsfc_dust_a3
    
    real(r8) :: r_bc                         ! model radii of BC modes [m]   
    real(r8) :: r_dust_a1, r_dust_a3         ! model radii of dust modes [m]  
    
    integer :: i
    
   integer  :: num_bc_idx, num_dst1_idx, num_dst3_idx    ! mode indices
    
    num_bc_idx = MODE_IDX_OMBC_INTMIX_COAT_AIT
    num_dst1_idx = MODE_IDX_DST_A2
    num_dst3_idx = MODE_IDX_DST_A3


!*****************************************************************************
!                calculate intersitial aerosol 
!*****************************************************************************

         dst1_num = qaerpt(num_dst1_idx)*1.0e-6_r8    ! #/cm3
         dst3_num = qaerpt(num_dst3_idx)*1.0e-6_r8    ! #/cm3
         bc_num = qaerpt(num_bc_idx)*1.0e-6_r8    ! #/cm3
     

!*****************************************************************************
!                calculate cloud borne aerosol 
!*****************************************************************************

     dst1_num_imm = qaercwpt(num_dst1_idx)*1.0e-6_r8    ! #/cm3
     dst3_num_imm = qaercwpt(num_dst3_idx)*1.0e-6_r8    ! #/cm3
     bc_num_imm = qaercwpt(num_bc_idx)*1.0e-6_r8    ! #/cm3
 
!   calculate mass mean radius
      r_dust_a1 = lifeCycleNumberMedianRadius(num_dst1_idx)
      r_dust_a3 = lifeCycleNumberMedianRadius(num_dst3_idx)
      r_bc = lifeCycleNumberMedianRadius(num_bc_idx)
  
    hetraer(1) = r_bc
    hetraer(2) = r_dust_a1
    hetraer(3) = r_dust_a3


!*****************************************************************************
!                calculate coated fraction 
!*****************************************************************************

! volumeCore and volumeCoat from subroutine calculateHygroscopicity in paramix_progncdnc.f90

    sigmag_amode(1) = lifeCycleSigma(num_bc_idx)
    sigmag_amode(2) = lifeCycleSigma(num_dst1_idx)
    sigmag_amode(3) = lifeCycleSigma(num_dst3_idx)

    fac_volsfc_bc = exp(2.5*(log(sigmag_amode(1))**2))
    fac_volsfc_dust_a1 = exp(2.5*(log(sigmag_amode(2))**2))
    fac_volsfc_dust_a3 = exp(2.5*(log(sigmag_amode(3))**2))

    tmp1 = volumeCoat(num_bc_idx)*(r_bc*2._r8)*fac_volsfc_bc
    tmp2 = max(6.0_r8*dr_so4_monolayers_dust*volumeCore(num_bc_idx), 0.0_r8) ! dr_so4_monolayers_dust = n_so4_monolayers_dust (=1) * 4.67e-10
    dstcoat(1) = tmp1/tmp2

    tmp1 = volumeCoat(num_dst1_idx)*(r_dust_a1*2._r8)*fac_volsfc_dust_a1
    tmp2 = max(6.0_r8*dr_so4_monolayers_dust*volumeCore(num_dst1_idx), 0.0_r8) ! dr_so4_monolayers_dust = n_so4_monolayers_dust (=1) * 4.67e-10
    dstcoat(2) = tmp1/tmp2
    
    tmp1 = volumeCoat(num_dst3_idx)*(r_dust_a3*2._r8)*fac_volsfc_dust_a3
    tmp2 = max(6.0_r8*dr_so4_monolayers_dust*volumeCore(num_dst3_idx), 0.0_r8) ! dr_so4_monolayers_dust = n_so4_monolayers_dust (=1) * 4.67e-10
    dstcoat(3) = tmp1/tmp2
    
    if (dstcoat(1) > 1._r8) dstcoat(1) = 1._r8
    if (dstcoat(1) < 0.001_r8) dstcoat(1) = 0.001_r8
    if (dstcoat(2) > 1._r8) dstcoat(2) = 1._r8
    if (dstcoat(2) < 0.001_r8) dstcoat(2) = 0.001_r8
    if (dstcoat(3) > 1._r8) dstcoat(3) = 1._r8
    if (dstcoat(3) < 0.001_r8) dstcoat(3) = 0.001_r8
 
!*****************************************************************************
!                prepare some variables for water activity 
!*****************************************************************************
! cam ([kg/m3] added mass distributed to modes) from paramix_progncdnc.f90
  
    ! accumulation mode for dust_a1 
    if (qaerpt(num_dst1_idx) > 0._r8) then 
       awcam(2) = cam(num_dst1_idx)*1.e9_r8    ! kg/m3 -> ug/m3
    else
        awcam(2) = 0._r8
    end if
    if (awcam(2) >0._r8) then
        awfacm(2) = qaerpt(num_dst1_idx)*(as_acm)/(as_acm+as_condm)
    else
        awfacm(2) = 0._r8
    end if
    
    ! accumulation mode for dust_a3
    if (qaerpt(num_dst3_idx) > 0._r8) then 
        awcam(3) = cam(num_dst3_idx)*1.e9_r8    ! kg/m3 -> ug/m3
    else
        awcam(3) = 0._r8
    end if
    awfacm(3) = 0._r8
        
        
    ! accumulation mode for bc
    if (qaerpt(num_bc_idx) > 0._r8) then 
        awcam(1) = cam(num_bc_idx)*1.e9_r8    ! kg/m3 -> ug/m3
    else
        awcam(1) = 0._r8
    end if
    awfacm(1) = awfacm(2)


!*****************************************************************************
!                prepare output
!*****************************************************************************

    total_interstial_aer_num(1) = bc_num 
    total_interstial_aer_num(2) = dst1_num 
    total_interstial_aer_num(3) = dst3_num 

    total_cloudborne_aer_num(1) = bc_num_imm 
    total_cloudborne_aer_num(2) = dst1_num_imm 
    total_cloudborne_aer_num(3) = dst3_num_imm
  
    do i = 1, 3
        total_aer_num(i) = total_interstial_aer_num(i)+total_cloudborne_aer_num(i)
        coated_aer_num(i) = total_interstial_aer_num(i)*dstcoat(i)
        uncoated_aer_num(i) = total_interstial_aer_num(i)*(1._r8-dstcoat(i))
    end do


    tot_na500 = total_aer_num(1)*0.0256_r8          & ! scaled for D>0.5 um using Clarke et al., 1997; 2004; 2007: rg=0.1um, sig=1.6
!#ifdef MODAL_AERO
!#if (defined MODAL_AERO_3MODE)
        +total_aer_num(2)*0.488_r8                  &  ! scaled for D>0.5-1 um from 0.1-1 um
!#elif (defined MODAL_AERO_7MODE)
!        +total_aer_num(2)*0.566_r8                  &  ! scaled for D>0.5-2 um from 0.1-2 um
!#endif
        +total_aer_num(3)
!#endif

    na500 = total_interstial_aer_num(1)*0.0256_r8   & ! scaled for D>0.5 um using Clarke et al., 1997; 2004; 2007: rg=0.1um, sig=1.6
!#ifdef MODAL_AERO
!#if (defined MODAL_AERO_3MODE)
        +total_interstial_aer_num(2)*0.488_r8       &  ! scaled for D>0.5-1 um from 0.1-1 um
!#elif (defined MODAL_AERO_7MODE)
!        +total_interstial_aer_num(2)*0.566_r8       &  ! scaled for D>0.5-2 um from 0.1-2 um
!#endif
        +total_interstial_aer_num(3)
!#endif

!-- wy4.0
  
end subroutine get_aer_num
