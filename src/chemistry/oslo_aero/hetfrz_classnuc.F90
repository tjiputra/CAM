!-----------------------------------------------------------------------
!
! Purpose: Calculate heterogeneous freezing rates from classical nucleation theory
!
! Public interfaces:
!
! hetfrz_classnuc 
!
! Author: Corinna Hoose, UiO, May 2009
!
! Modifications: Yong Wang - Implement classical nucleation theory
!                            with 2 different ice nucleation schemes(single theta,
!                            PDF theta) into CAM5.3
!-----------------------------------------------------------------------
      
subroutine hetfrz_classnuc(deltat, t, p, supersatice, &
                            !++wy
                            fn_in, &                                   &
                            !--wy
                            r3lx, icnlx,                               &
                            frzbcimm, frzduimm,                        &
                            frzbccnt, frzducnt,                        &
                            frzbcdep, frzdudep,                        &
                            !++wy
                            hetraer, awcam, awfacm, dstcoat,                         &
                            total_aer_num, coated_aer_num, uncoated_aer_num,  &
                            total_interstitial_aer_num, total_cloudborne_aer_num)
!                            total_interstitial_aer_num, total_cloudborne_aer_num, errstring)
                            !--wy

    use shr_kind_mod,  only: r8 => shr_kind_r8,r4 => shr_kind_r4
    use physconst,     only: rhoh2o, mwh2o
    !++wy
    use wv_saturation, only: svp_water, svp_ice
    !--wy
    use abortutils,    only: endrun
    !++wy
!    use rad_constituents, only: rad_cnst_get_info
!#ifdef MODAL_AERO
    use modal_aero_data, only: ntot_amode
!#endif
    use cam_logfile,    only: iulog
    !--wy

    implicit none
 
!   input
    real(r8), intent(in) :: deltat            ! timestep [s]
    real(r8), intent(in) :: t                 ! temperature [K]
    real(r8), intent(in) :: p                 ! pressure [Pa]
    real(r8), intent(in) :: supersatice       ! supersaturation ratio wrt ice at 100%rh over water [ ]
    real(r8), intent(in) :: r3lx              ! volume mean drop radius [m]
    real(r8), intent(in) :: icnlx             ! in-cloud droplet concentration [cm-3]
    !++wy
!    real(r8), intent(in) :: fn_in(ntot_amode)   ! fraction activated [ ] for aerosols in cloud droplets
    real(r8), intent(in) :: fn_in(3)            ! fraction activated [ ] for cloud borne aerosol number
                                                ! index values are 1:bc, 2:dust_a1, 3:dust_a3
    real(r8), intent(in) :: hetraer(3)          ! bc and dust mass mean radius [m]
    real(r8), intent(in) :: awcam(3)            ! modal added mass [mug m-3]
    real(r8), intent(in) :: awfacm(3)           ! (OC+BC)/(OC+BC+SO4)
    real(r8), intent(in) :: dstcoat(3)          ! coated fraction
    real(r8), intent(in) :: total_aer_num(3)    ! total bc and dust number concentration(interstitial+cloudborne) [#/cm^3]
    real(r8), intent(in) :: coated_aer_num(3)   ! coated bc and dust number concentration(interstitial)
    real(r8), intent(in) :: uncoated_aer_num(3) ! uncoated bc and dust number concentration(interstitial)
    real(r8), intent(in) :: total_interstitial_aer_num(3) ! total bc and dust concentration(interstitial)
    real(r8), intent(in) :: total_cloudborne_aer_num(3)   ! total bc and dust concentration(cloudborne)
    !--wy
!   output
    real(r8), intent(out) :: frzbcimm           ! het. frz by BC immersion nucleation [cm-3 s-1]
    real(r8), intent(out) :: frzduimm           ! het. frz by dust immersion nucleation [cm-3 s-1]
    real(r8), intent(out) :: frzbccnt           ! het. frz by BC contact nucleation [cm-3 s-1]
    real(r8), intent(out) :: frzducnt           ! het. frz by dust contact nucleation [cm-3 s-1]
    real(r8), intent(out) :: frzbcdep           ! het. frz by BC deposition nucleation [cm-3 s-1]
    real(r8), intent(out) :: frzdudep           ! het. frz by dust deposition nucleation [cm-3 s-1]

!    character(len=*), intent(out) :: errstring

!   local variables
!++wy
!------------water activity variables-------------
!    integer :: ntot_amode
    real(r8) :: aw(3)                           ! water activity [ ]
    real(r8) :: molal(3)                        ! molality [moles/kg]
    real(r8), parameter :: Mso4 = 96.06
!-------------------------------------------------
!--wy
!    integer  :: id_bc, id_dst1, id_dst3         ! index
    integer, parameter :: id_bc   = 1
    integer, parameter :: id_dst1 = 2
    integer, parameter :: id_dst3 = 3
    logical :: do_bc, do_dst1, do_dst3

    real(r8), parameter :: t0 = 273.15_r8       ! freezing temperature [K]
    real(r8), parameter :: n1 = 1.e19_r8        ! number of water molecules in contact with unit area of substrate [m-2]
    real(r8), parameter :: kboltz = 1.38e-23_r8    
    real(r8), parameter :: hplanck = 6.63e-34_r8
    real(r8), parameter :: rhplanck = 1._r8/hplanck
    real(r8), parameter :: amu = 1.66053886e-27_r8 
    real(r8), parameter :: nus = 1.e13_r8       ! frequ. of vibration [s-1] higher freq. (as in P&K, consistent with Anupam's data) 
    real(r8), parameter :: taufrz = 195.435_r8  ! time constant for falloff of freezing rate [s]
    real(r8), parameter :: rhwincloud = 0.98_r8 ! 98% RH in mixed-phase clouds (Korolev & Isaac, JAS 2006)
    real(r8), parameter :: limfacbc = 0.01_r8   ! max. ice nucleating fraction soot
    real(r8), parameter :: pi = 4._r8*atan(1.0_r8)   
    real(r8) :: tc   
    real(r8) :: vwice   
    real(r8) :: rhoice   
    real(r8) :: sigma_iw                        ! [J/m2]   
    real(r8) :: sigma_iv                        ! [J/m2]   
    real(r8) :: esice                           ! [Pa]   
    real(r8) :: eswtr                           ! [Pa]   
    real(r8) :: rgimm   
    real(r8) :: rgdep   
    real(r8) :: dg0dep   
    real(r8) :: Adep   
    real(r8) :: dg0cnt   
    real(r8) :: Acnt   
    real(r8) :: rgimm_bc   
    real(r8) :: rgimm_dust_a1, rgimm_dust_a3   
    real(r8) :: dg0imm_bc  
    real(r8) :: dg0imm_dust_a1, dg0imm_dust_a3  
    real(r8) :: Aimm_bc
    real(r8) :: Aimm_dust_a1, Aimm_dust_a3
    real(r8) :: q, m, phi   
    real(r8) :: r_bc                            ! model radii of BC modes [m]   
    real(r8) :: r_dust_a1, r_dust_a3            ! model radii of dust modes [m]   
    real(r8) :: f_imm_bc 
    real(r8) :: f_imm_dust_a1, f_imm_dust_a3 
    real(r8) :: Jimm_bc
    real(r8) :: Jimm_dust_a1, Jimm_dust_a3
    real(r8) :: f_dep_bc 
    real(r8) :: f_dep_dust_a1, f_dep_dust_a3  
    real(r8) :: Jdep_bc 
    real(r8) :: Jdep_dust_a1, Jdep_dust_a3 
    real(r8) :: f_cnt_bc   
    real(r8) :: f_cnt_dust_a1,f_cnt_dust_a3   
    real(r8) :: Jcnt_bc
    real(r8) :: Jcnt_dust_a1,Jcnt_dust_a3 
    integer :: i
    !++ wy4.0
    !********************************************************
    ! Hoose et al., 2010 fitting parameters
    ! Use equations labelled with H10
    !********************************************************
    !freezing parameters for immersion freezing
    !real(r8),parameter :: theta_imm_bc = 40.17         ! contact angle [deg], converted to rad later
    !real(r8),parameter :: dga_imm_bc = 14.4E-20        ! activation energy [J]
    !real(r8),parameter :: theta_imm_dust = 30.98       ! contact angle [deg], converted to rad later
    !real(r8),parameter :: dga_imm_dust = 15.7E-20      ! activation energy [J]
    !freezing parameters for deposition nucleation
    !real(r8),parameter :: theta_dep_dust = 12.7        ! contact angle [deg], converted to rad later !Zimmermann et al (2008), illite
    !real(r8),parameter :: dga_dep_dust = -6.21E-21     ! activation energy [J]
    !real(r8),parameter :: theta_dep_bc = 28.           ! contact angle [deg], converted to rad later !Moehler et al (2005), soot
    !real(r8),parameter :: dga_dep_bc = -2.E-19         ! activation energy [J]
    !********************************************************
    ! Our fitting parameters
    !********************************************************
    ! freezing parameters for immersion freezing
    real(r8),parameter :: theta_imm_bc = 48.0_r8            ! contact angle [deg], converted to rad later !DeMott et al (1990)
    real(r8),parameter :: dga_imm_bc = 14.15E-20_r8         ! activation energy [J]
    real(r8),parameter :: theta_imm_dust = 46.0_r8          ! contact angle [deg], converted to rad later !DeMott et al (2011) SD
    real(r8),parameter :: dga_imm_dust = 14.75E-20_r8       ! activation energy [J]
    ! freezing parameters for deposition nucleation
    real(r8),parameter :: theta_dep_dust = 20.0_r8          ! contact angle [deg], converted to rad later !Koehler et al (2010) SD
    real(r8),parameter :: dga_dep_dust = -8.1E-21_r8        ! activation energy [J]
    real(r8),parameter :: theta_dep_bc = 28._r8             ! contact angle [deg], converted to rad later !Moehler et al (2005), soot
    real(r8),parameter :: dga_dep_bc = -2.E-19_r8           ! activation energy [J]
    !-- wy4.0
    !real(r8),parameter :: Kcoll_bc = 3.E-6             ! collision kernel [cm3 s-1]
    !real(r8),parameter :: Kcoll_dust_a1 = 4.E-7        ! collision kernel [cm3 s-1]
    !real(r8),parameter :: Kcoll_dust_a3 = 2.E-7        ! collision kernel [cm3 s-1]
    !real(r8) :: Kcoll_bc = 3.E-6             ! collision kernel [cm3 s-1]
    !real(r8) :: Kcoll_dust_a1 = 4.E-7        ! collision kernel [cm3 s-1]
    !real(r8) :: Kcoll_dust_a3 = 2.E-7        ! collision kernel [cm3 s-1]
    !++ wy4.0
    real(r8) :: Kcoll_bc                                    ! collision kernel [cm3 s-1]
    real(r8) :: Kcoll_dust_a1                               ! collision kernel [cm3 s-1]
    real(r8) :: Kcoll_dust_a3                               ! collision kernel [cm3 s-1]
    !-- wy4.0
!++wy
    logical :: tot_in = .false.                             ! MH false when cloudborne aerosol concentration is available; uses fn
    real(r8) :: fn(3) ! 1:bc,2:dust_a1,3:dust_a3
!--wy 

    !++ wy4.0
    !*****************************************************************************
    !                PDF theta model 
    !*****************************************************************************
    ! some variables for PDF theta model
    ! immersion freezing
    real(r8),parameter :: theta_min = 1._r8/180._r8*pi
    real(r8),parameter :: theta_max = 179._r8/180._r8*pi
    real(r8) :: x1_imm
    real(r8) :: x2_imm
    real(r8) :: norm_theta_imm
    real(r8),parameter :: imm_dust_mean_theta = 46.0_r8/180.0_r8*pi 
    real(r8),parameter :: imm_dust_var_theta = 0.01_r8
    real(r8) :: pdf_d_theta
    integer,parameter :: pdf_n_theta = 501 !101   ! MH_2015/06/01
    real(r8) :: dim_theta(pdf_n_theta)
    real(r8) :: dim_f_imm_dust_a1(pdf_n_theta), dim_f_imm_dust_a3(pdf_n_theta)
    real(r8) :: dim_Jimm_dust_a1(pdf_n_theta), dim_Jimm_dust_a3(pdf_n_theta)
    real(r8) :: pdf_imm_theta(pdf_n_theta)
    real(r8) :: sum_imm_dust_a1, sum_imm_dust_a3
    logical :: pdf_imm_in = .true.
    !-- wy4.0
   
    
!++ classnuc wy
!    call rad_cnst_get_info(0, nmodes=ntot_amode)
!-- classnuc wy

!    errstring = ' '

!    id_bc = 1
!    id_dst1 = 2
!    id_dst3 = 3

    !++ wy4.0
    if (pdf_imm_in) then
        pdf_d_theta = (179._r8-1._r8)/180._r8*pi/(pdf_n_theta-1)
        ! calculate the integral in the denominator
        x1_imm = (LOG(theta_min)-LOG(imm_dust_mean_theta))/(sqrt(2.0_r8)*imm_dust_var_theta)
        x2_imm = (LOG(theta_max)-LOG(imm_dust_mean_theta))/(sqrt(2.0_r8)*imm_dust_var_theta)
        norm_theta_imm = (ERFAPP(x2_imm)-ERFAPP(x1_imm))*0.5_r8
        do i = 1, pdf_n_theta
            dim_theta(i) = 1._r8/180._r8*pi+(i-1)*pdf_d_theta
            pdf_imm_theta(i) = exp(-((LOG(dim_theta(i))-LOG(imm_dust_mean_theta))**2._r8)/(2._r8*imm_dust_var_theta**2._r8))/ &
                                (dim_theta(i)*imm_dust_var_theta*SQRT(2._r8*pi))/norm_theta_imm
        end do
    end if
    !-- wy4.0


!++wy
    fn(1) = fn_in(1) ! bc accumulation mode
    fn(2) = fn_in(1) ! dust_a1 accumulation mode
    fn(3) = fn_in(3) ! dust_a3 coarse mode

!   get saturation vapor pressures
    eswtr = svp_water(t)  ! 0 for liquid
    esice = svp_ice(t)  ! 1 for ice
!--wy

!   calculation
    tc = t-t0
    rhoice = 916.7-0.175*tc-5.e-4*tc**2
    vwice = mwh2o*amu/rhoice
    sigma_iw = (28.5+0.25*tc)*1E-3
    sigma_iv = (76.1-0.155*tc + 28.5+0.25*tc)*1E-3
!++wy
!   get mass mean radius
    r_bc = hetraer(1)    
    r_dust_a1 = hetraer(2)    
    r_dust_a3 = hetraer(3)    
!--wy


!   calculate collision kernels as a function of environmental parameters and aerosol/droplet sizes
    call collkernel(t, p, eswtr, rhwincloud, r3lx,         &
                    r_bc,                                  &  ! BC modes
                    r_dust_a1, r_dust_a3,                  &  ! dust modes
                    Kcoll_bc,                              &  ! collision kernel [cm3 s-1]
                    Kcoll_dust_a1, Kcoll_dust_a3)
        
!*****************************************************************************
!                take water activity into account 
!*****************************************************************************
!   solute effect
    aw(:) = 1._r8
    molal(:) = 0._r8

!++wy
! The heterogeneous ice freezing temperatures of all IN generally decrease with
! increasing total solute mole fraction. Therefore, the large solution concentration
! will cause the freezing point depression and the ice freezing temperatures of all
! IN will get close to the homogeneous ice freezing temperatures. Since we take into
! account water activity for three heterogeneous freezing modes(immersion, deposition, 
! and contact), we utilize interstitial aerosols(not cloudborne aerosols) to calculate 
! water activity. 
! If the index of IN is 0, it means three freezing modes of this aerosol are depressed.

    do i = 1, 3
        !calculate molality
        if ( total_interstitial_aer_num(i) > 0._r8 ) then
            molal(i) = (1.e-6*awcam(i)*(1.-awfacm(i))/(Mso4*total_interstitial_aer_num(i)*1.e6))/(4*pi/3*rhoh2o*(MAX(r3lx,4.e-6))**3)
            aw(i) = 1./(1.+2.9244948e-2*molal(i)+2.3141243e-3*molal(i)**2+7.8184854e-7*molal(i)**3)
        end if
    end do

!*****************************************************************************
!                immersion freezing begin 
!*****************************************************************************    
!--wy

!   initialization
    frzbcimm = 0.
    frzduimm = 0.
    frzbccnt = 0.
    frzducnt = 0.
    frzbcdep = 0.
    frzdudep = 0.

!   critical germ size
    rgimm = 2*vwice*sigma_iw/(kboltz*t*LOG(supersatice))
!   take solute effect into account
    rgimm_bc = rgimm !initialize
    rgimm_dust_a1 = rgimm
    rgimm_dust_a3 = rgimm

   if (aw(id_bc)*supersatice > 1._r8 ) then
      do_bc   = .true.
      rgimm_bc = 2*vwice*sigma_iw/(kboltz*t*LOG(aw(id_bc)*supersatice))
   else
      do_bc = .false.
   end if

   if (aw(id_dst1)*supersatice > 1._r8 ) then
      do_dst1 = .true.
      rgimm_dust_a1 = 2*vwice*sigma_iw/(kboltz*t*LOG(aw(id_dst1)*supersatice))
   else
      do_dst1 = .false.
   end if

   if (aw(id_dst3)*supersatice > 1._r8 ) then
      do_dst3 = .true.
      rgimm_dust_a3 = 2*vwice*sigma_iw/(kboltz*t*LOG(aw(id_dst3)*supersatice))
   else
      do_dst3 = .false.
   end if

!    if ( id_bc   .gt. 0 .and. aw(id_bc)*supersatice .gt. 1. ) rgimm_bc = 2*vwice*sigma_iw/(kboltz*t*LOG(aw(id_bc)*supersatice))
!    if ( id_dst1 .gt. 0 .and. aw(id_dst1)*supersatice .gt. 1. ) rgimm_dust_a1 = 2*vwice*sigma_iw/(kboltz*t*LOG(aw(id_dst1)*supersatice))
!    if ( id_dst3 .gt. 0 .and. aw(id_dst3)*supersatice .gt. 1. ) rgimm_dust_a3 = 2*vwice*sigma_iw/(kboltz*t*LOG(aw(id_dst3)*supersatice))
!!   if aw*Si<=1, the freezing point depression is strong enough to prevent freezing for this mode:
!    if ( id_bc   .gt. 0 .and. aw(id_bc)*supersatice .le. 1.) id_bc = 0
!    if ( id_dst1 .gt. 0 .and. aw(id_dst1)*supersatice .le. 1.) id_dst1 = 0
!    if ( id_dst3 .gt. 0 .and. aw(id_dst3)*supersatice .le. 1.) id_dst3 = 0
    
!   form factor
!++ wy4.0
!only consider flat surfaces due to uncertainty of curved surfaces
!   flat surface
    m = COS(theta_imm_bc*pi/180.)
    f_imm_bc = (2+m)*(1-m)**2/4.
    if (.not. pdf_imm_in) then
        m = COS(theta_imm_dust*pi/180.)
        f_imm_dust_a1 = (2+m)*(1-m)**2/4.

        m = COS(theta_imm_dust*pi/180.)
        f_imm_dust_a3 = (2+m)*(1-m)**2/4.
    else
        do i = 1, pdf_n_theta
            m = cos(dim_theta(i))
            dim_f_imm_dust_a1(i) = (2+m)*(1-m)**2/4.

            m = cos(dim_theta(i))
            dim_f_imm_dust_a3(i) = (2+m)*(1-m)**2/4.
        end do
    end if
!-- wy4.0

! H10  ++ MH_2015/02/23
!    m = COS(theta_imm_bc*pi/180.)
!    q =r_bc/rgimm_bc
!    phi = SQRT(1-2*m*q+q*q)
!    f_imm_bc = 0.5*(1+((1-m*q)/phi)**3+q**3*(2-3*((q-m)/phi)+((q-m)/phi)**3) &
!                +3*m*q**2*((q-m)/phi-1))
!    m = COS(theta_imm_dust*pi/180.)
!    q = r_dust_a1/rgimm_dust_a1
!    phi = SQRT(1-2*m*q+q*q)
!    f_imm_dust_a1 = 0.5*(1+((1-m*q)/phi)**3+q**3*(2-3*((q-m)/phi)+((q-m)/phi)**3) &
!                +3*m*q**2*((q-m)/phi-1))
!    m = COS(theta_imm_dust*pi/180.)
!    q = r_dust_a3/rgimm_dust_a3
!    phi = SQRT(1-2*m*q+q*q)
!    f_imm_dust_a3 = 0.5*(1+((1-m*q)/phi)**3+q**3*(2-3*((q-m)/phi)+((q-m)/phi)**3) &
!                +3*m*q**2*((q-m)/phi-1))
!-- MH_2015/02/23

!   homogeneous energy of germ formation
    dg0imm_bc = 4*pi/3.*sigma_iw*rgimm_bc**2
    dg0imm_dust_a1 = 4*pi/3.*sigma_iw*rgimm_dust_a1**2
    dg0imm_dust_a3 = 4*pi/3.*sigma_iw*rgimm_dust_a3**2

!   prefactor
    Aimm_bc = n1*((vwice*rhplanck)/(rgimm_bc**3)*SQRT(3./pi*kboltz*T*dg0imm_bc))
    Aimm_dust_a1 = n1*((vwice*rhplanck)/(rgimm_dust_a1**3)*SQRT(3./pi*kboltz*T*dg0imm_dust_a1))
    Aimm_dust_a3 = n1*((vwice*rhplanck)/(rgimm_dust_a3**3)*SQRT(3./pi*kboltz*T*dg0imm_dust_a3))

!   nucleation rate per particle
    !++ wy4.0
    Jimm_bc = Aimm_bc*r_bc**2/SQRT(f_imm_bc)*EXP((-dga_imm_bc-f_imm_bc*dg0imm_bc)/(kboltz*T))
    if (.not. pdf_imm_in) then
        ! 1/sqrt(f)
        ! the expression of Chen et al. (sqrt(f)) may however lead to unphysical
        ! behavior as it implies J->0 when f->0 (i.e. ice nucleation would be
        ! more difficult on easily wettable materials). 
        Jimm_dust_a1 = Aimm_dust_a1*r_dust_a1**2/SQRT(f_imm_dust_a1)*EXP((-dga_imm_dust-f_imm_dust_a1*dg0imm_dust_a1)/(kboltz*T))
        Jimm_dust_a3 = Aimm_dust_a3*r_dust_a3**2/SQRT(f_imm_dust_a3)*EXP((-dga_imm_dust-f_imm_dust_a3*dg0imm_dust_a3)/(kboltz*T))
    end if
    if (pdf_imm_in) then
        do i = 1, pdf_n_theta
            ! 1/sqrt(f)
            dim_Jimm_dust_a1(i) = Aimm_dust_a1*r_dust_a1**2/SQRT(dim_f_imm_dust_a1(i))*EXP((-dga_imm_dust-dim_f_imm_dust_a1(i)*dg0imm_dust_a1)/(kboltz*T))
            dim_Jimm_dust_a1(i) = max(dim_Jimm_dust_a1(i), 0._r8)
            dim_Jimm_dust_a3(i) = Aimm_dust_a3*r_dust_a3**2/SQRT(dim_f_imm_dust_a3(i))*EXP((-dga_imm_dust-dim_f_imm_dust_a3(i)*dg0imm_dust_a3)/(kboltz*T))
            dim_Jimm_dust_a3(i) = max(dim_Jimm_dust_a3(i), 0._r8)
        end do
    end if
    !-- wy4.0

! H10  ++ MH_2015/02/23
!    Jimm_bc = Aimm_bc*r_bc**2*SQRT(f_imm_bc)*EXP((-dga_imm_bc-f_imm_bc*dg0imm_bc)/(kboltz*T))
!    Jimm_dust_a1 = Aimm_dust_a1*r_dust_a1**2*SQRT(f_imm_dust_a1)*EXP((-dga_imm_dust-f_imm_dust_a1*dg0imm_dust_a1)/(kboltz*T))
!    Jimm_dust_a3 = Aimm_dust_a3*r_dust_a3**2*SQRT(f_imm_dust_a3)*EXP((-dga_imm_dust-f_imm_dust_a3*dg0imm_dust_a3)/(kboltz*T))
!-- MH_2015/02/23

!   Limit to 1% of available potential IN (for BC), no limit for dust 
!++wy 4.0
    if (pdf_imm_in) then
        sum_imm_dust_a1 = 0._r8
        sum_imm_dust_a3 = 0._r8
        do i = 1, pdf_n_theta-1
            if (sum_imm_dust_a1 > 0.99_r8) then
                sum_imm_dust_a1 = 1.0_r8
            else
                sum_imm_dust_a1 = sum_imm_dust_a1+0.5_r8*((pdf_imm_theta(i)*exp(-dim_Jimm_dust_a1(i)*deltat)+ &
                                pdf_imm_theta(i+1)*exp(-dim_Jimm_dust_a1(i+1)*deltat)))*pdf_d_theta
           end if
           if (sum_imm_dust_a3 > 0.99_r8) then
               sum_imm_dust_a3 = 1.0_r8
          else
               sum_imm_dust_a3 = sum_imm_dust_a3+0.5_r8*((pdf_imm_theta(i)*exp(-dim_Jimm_dust_a3(i)*deltat)+ &
                                pdf_imm_theta(i+1)*exp(-dim_Jimm_dust_a3(i+1)*deltat)))*pdf_d_theta
          end if
        end do
    end if
!--wy 4.0
!++wy
    if (.not.tot_in) then
        if ( do_bc )   frzbcimm = frzbcimm+MIN(limfacbc*total_cloudborne_aer_num(id_bc)/deltat, &
                                        total_cloudborne_aer_num(id_bc)/deltat*(1.-exp(-Jimm_bc*deltat))) 
        !++ wy4.0
        if (.not. pdf_imm_in) then
            if ( do_dst1 ) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst1)/deltat, & 
                                            total_cloudborne_aer_num(id_dst1)/deltat*(1.-exp(-Jimm_dust_a1*deltat)))
            if ( do_dst3 ) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst3)/deltat, &
                                            total_cloudborne_aer_num(id_dst3)/deltat*(1.-exp(-Jimm_dust_a3*deltat)))
        else
            if ( do_dst1 ) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst1)/deltat,        &
                                            total_cloudborne_aer_num(id_dst1)/deltat*(1.-sum_imm_dust_a1))
            if ( do_dst3 ) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst3)/deltat,        &
                                            total_cloudborne_aer_num(id_dst3)/deltat*(1.-sum_imm_dust_a3))
        end if
        !--wy4.0
    else
        if ( do_bc )   frzbcimm = frzbcimm+MIN(limfacbc*fn(id_bc)*total_aer_num(id_bc)/deltat, & 
                                        fn(id_bc)*total_aer_num(id_bc)/deltat*(1.-exp(-Jimm_bc*deltat))) 
        !++wy4.0
        if (.not. pdf_imm_in) then
            if ( do_dst1 ) frzduimm = frzduimm+MIN(1*fn(id_dst1)*total_aer_num(id_dst1)/deltat, &
                                            fn(id_dst1)*total_aer_num(id_dst1)/deltat*(1.-exp(-Jimm_dust_a1*deltat)))
            if ( do_dst3 ) frzduimm = frzduimm+MIN(1*fn(id_dst3)*total_aer_num(id_dst3)/deltat, &
                                            fn(id_dst3)*total_aer_num(id_dst3)/deltat*(1.-exp(-Jimm_dust_a3*deltat)))
        else
            if ( do_dst1 ) frzduimm = frzduimm+MIN(1*fn(id_dst1)*total_aer_num(id_dst1)/deltat,        &
                                            fn(id_dst1)*total_aer_num(id_dst1)/deltat*(1.-sum_imm_dust_a1))
            if ( do_dst3 ) frzduimm = frzduimm+MIN(1*fn(id_dst3)*total_aer_num(id_dst3)/deltat,        &
                                            fn(id_dst3)*total_aer_num(id_dst3)/deltat*(1.-sum_imm_dust_a3))
        end if
        !-- wy4.0
!        if ( id_bc .gt. 0 )   frzbcimm = frzbcimm+MIN(limfacbc*total_cloudborne_aer_num(id_bc)/deltat, &
!                                        total_cloudborne_aer_num(id_bc)/deltat*(1.-exp(-Jimm_bc*deltat))) 
!        if ( id_dst1 .gt. 0 ) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst1)/deltat, & 
!                                        total_cloudborne_aer_num(id_dst1)/deltat*(1.-exp(-Jimm_dust_a1*deltat)))
!        if ( id_dst3 .gt. 0 ) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst3)/deltat, &
!                                        total_cloudborne_aer_num(id_dst3)/deltat*(1.-exp(-Jimm_dust_a3*deltat)))
!    else
!        if ( id_bc .gt. 0 )   frzbcimm = frzbcimm+MIN(limfacbc*fn(id_bc)*total_aer_num(id_bc)/deltat, & 
!                                        fn(id_bc)*total_aer_num(id_bc)/deltat*(1.-exp(-Jimm_bc*deltat))) 
!        if ( id_dst1 .gt. 0 ) frzduimm = frzduimm+MIN(1*fn(id_dst1)*total_aer_num(id_dst1)/deltat, &
!                                        fn(id_dst1)*total_aer_num(id_dst1)/deltat*(1.-exp(-Jimm_dust_a1*deltat)))
!        if ( id_dst3 .gt. 0 ) frzduimm = frzduimm+MIN(1*fn(id_dst3)*total_aer_num(id_dst3)/deltat, &
!                                        fn(id_dst3)*total_aer_num(id_dst3)/deltat*(1.-exp(-Jimm_dust_a3*deltat)))
    end if
!--wy


   if (t > 263.15_r8) then
      frzduimm = 0._r8
      frzbcimm = 0._r8
   end if

!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Deposition nucleation
!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   critical germ size
!   assume 98% RH in mixed-phase clouds (Korolev & Isaac, JAS 2006)
    rgdep=2*vwice*sigma_iv/(kboltz*t*LOG(rhwincloud*supersatice)) 

!   form factor
    !++ wy4.0
    m = COS(theta_dep_bc*pi/180.)
    f_dep_bc = (2+m)*(1-m)**2/4.

    m = COS(theta_dep_dust*pi/180.)
    f_dep_dust_a1 = (2+m)*(1-m)**2/4.

    m = COS(theta_dep_dust*pi/180.)
    f_dep_dust_a3 = (2+m)*(1-m)**2/4.
    !-- wy4.0

!H10  ++ MH_2015/02/23
!    m = COS(theta_dep_bc*pi/180.)
!    q = r_bc/rgdep
!    phi = SQRT(1-2*m*q+q*q)
!    f_dep_bc = 0.5*(1+((1-m*q)/phi)**3+q**3*(2-3*((q-m)/phi)+((q-m)/phi)**3) &
!            +3*m*q**2*((q-m)/phi-1))
!    m = COS(theta_dep_dust*pi/180.)
!    q = r_dust_a1/rgdep
!    phi = SQRT(1-2*m*q+q*q)
!    f_dep_dust_a1 =0.5*(1+((1-m*q)/phi)**3+q**3*(2-3*((q-m)/phi)+((q-m)/phi)**3) &
!            +3*m*q**2*((q-m)/phi-1))
!    m = COS(theta_dep_dust*pi/180.)
!    q = r_dust_a3/rgdep
!    phi = SQRT(1-2*m*q+q*q)
!    f_dep_dust_a3 = 0.5*(1+((1-m*q)/phi)**3+q**3*(2-3*((q-m)/phi)+((q-m)/phi)**3) &
!            +3*m*q**2*((q-m)/phi-1))
!-- MH_2015/02/23

!   homogeneous energy of germ formation
    dg0dep = 4*pi/3.*sigma_iv*rgdep**2

!   prefactor
!   attention: division of small numbers
    Adep = (rhwincloud*eswtr)**2*(vwice/(mwh2o*amu))/(kboltz*T*nus)*SQRT(sigma_iv/(kboltz*T))

!   nucleation rate per particle
    if ( rgdep .gt. 0 ) then
        !++ wy4.0

        Jdep_bc = Adep*r_bc**2/SQRT(f_dep_bc)*EXP((-dga_dep_bc-f_dep_bc*dg0dep)/(kboltz*T))
        Jdep_dust_a1 = Adep*r_dust_a1**2/SQRT(f_dep_dust_a1)*EXP((-dga_dep_dust-f_dep_dust_a1*dg0dep)/(kboltz*T))
        Jdep_dust_a3 = Adep*r_dust_a3**2/SQRT(f_dep_dust_a3)*EXP((-dga_dep_dust-f_dep_dust_a3*dg0dep)/(kboltz*T))
        !-- wy4.0

!H10  ++ MH_2015/02/23
!        Jdep_bc = Adep*r_bc**2*SQRT(f_dep_bc)*EXP((-dga_dep_bc-f_dep_bc*dg0dep)/(kboltz*T))
!        Jdep_dust_a1 = Adep*r_dust_a1**2*SQRT(f_dep_dust_a1)*EXP((-dga_dep_dust-f_dep_dust_a1*dg0dep)/(kboltz*T))
!        Jdep_dust_a3 = Adep*r_dust_a3**2*SQRT(f_dep_dust_a3)*EXP((-dga_dep_dust-f_dep_dust_a3*dg0dep)/(kboltz*T))
!-- MH_2015/02/23
    else
        Jdep_bc = 0.
        Jdep_dust_a1 = 0.
        Jdep_dust_a3 = 0.
    end if

!   Limit to 1% of available potential IN (for BC), no limit for dust 
!++wy
    if (.not.tot_in) then
        if ( do_bc ) frzbcdep = frzbcdep+MIN(limfacbc*uncoated_aer_num(id_bc)/deltat, &
                                               uncoated_aer_num(id_bc)/deltat &
                                               *(1.-exp(-Jdep_bc*deltat))) 
        if ( do_dst1 ) frzdudep = frzdudep+MIN(1*uncoated_aer_num(id_dst1)/deltat, &
                                               uncoated_aer_num(id_dst1)/deltat &
                                               *(1.-exp(-Jdep_dust_a1*deltat)))
        if ( do_dst3 ) frzdudep = frzdudep+MIN(1*uncoated_aer_num(id_dst3)/deltat, &
                                               uncoated_aer_num(id_dst3)/deltat &
                                               *(1.-exp(-Jdep_dust_a3*deltat)))
    else
        if ( do_bc ) frzbcdep = frzbcdep+MIN(limfacbc*(1.-fn(id_bc)) &
                                               *(1.-dstcoat(1))*total_aer_num(id_bc)/deltat, &
                                               (1.-fn(id_bc))*(1.-dstcoat(1))*total_aer_num(id_bc)/deltat &
                                               *(1.-exp(-Jdep_bc*deltat))) 
        if ( do_dst1 ) frzdudep = frzdudep+MIN(1.*(1.-fn(id_dst1)) &
                                               *(1.-dstcoat(2))*total_aer_num(id_dst1)/deltat, &
                                               (1.-fn(id_dst1))*(1.-dstcoat(2))*total_aer_num(id_dst1)/deltat &
                                               *(1.-exp(-Jdep_dust_a1*deltat))) 
        if ( do_dst3 ) frzdudep = frzdudep+MIN(1.*(1.-fn(id_dst3)) &
                                               *(1.-dstcoat(3))*total_aer_num(id_dst3)/deltat, &
                                               (1.-fn(id_dst3))*(1.-dstcoat(3))*total_aer_num(id_dst3)/deltat &
                                               *(1.-exp(-Jdep_dust_a3*deltat))) 
    end if
!--wy

!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   contact nucleation
!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   form factor
    !++wy4.0
    m = COS(theta_dep_bc*pi/180.)
    f_cnt_bc = (2+m)*(1-m)**2/4.

    m = COS(theta_dep_dust*pi/180.)
    f_cnt_dust_a1 = (2+m)*(1-m)**2/4.

    m = COS(theta_dep_dust*pi/180.)
    f_cnt_dust_a3 = (2+m)*(1-m)**2/4.
    !--wy4.0

!H10  ++ MH_2015/02/23
!    m = COS(theta_dep_bc*pi/180.)
!    q = r_bc/rgimm
!    phi = SQRT(1-2*m*q+q*q)
!    f_cnt_bc = 0.5*(1+((1-m*q)/phi)**3+q**3*(2-3*((q-m)/phi)+((q-m)/phi)**3) &
!                  +3*m*q**2*((q-m)/phi-1))

!    m = COS(theta_dep_dust*pi/180.)
!    q = r_dust_a1/rgimm
!    phi = SQRT(1-2*m*q+q*q)
!    f_cnt_dust_a1 = 0.5*(1+((1-m*q)/phi)**3+q**3*(2-3*((q-m)/phi)+((q-m)/phi)**3) &
!                       +3*m*q**2*((q-m)/phi-1))

!    m = COS(theta_dep_dust*pi/180.)
!    q = r_dust_a3/rgimm
!    phi = SQRT(1-2*m*q+q*q)
!    f_cnt_dust_a3 = 0.5*(1+((1-m*q)/phi)**3+q**3*(2-3*((q-m)/phi)+((q-m)/phi)**3) &
!                      +3*m*q**2*((q-m)/phi-1))
!-- MH_2015/02/23

!   homogeneous energy of germ formation
    dg0cnt = 4*pi/3.*sigma_iv*rgimm**2

!   prefactor
!   attention: division of small numbers
    Acnt = rhwincloud*eswtr*4*pi/(nus*SQRT(2*pi*mwh2o*amu*kboltz*T))

!   nucleation rate per particle
    Jcnt_bc = Acnt*r_bc**2*EXP((-dga_dep_bc-f_cnt_bc*dg0cnt)/(kboltz*T))*Kcoll_bc*icnlx
    Jcnt_dust_a1 = Acnt*r_dust_a1**2*EXP((-dga_dep_dust-f_cnt_dust_a1*dg0cnt)/(kboltz*T))*Kcoll_dust_a1*icnlx
    Jcnt_dust_a3 = Acnt*r_dust_a3**2*EXP((-dga_dep_dust-f_cnt_dust_a3*dg0cnt)/(kboltz*T))*Kcoll_dust_a3*icnlx

!   Limit to 1% of available potential IN (for BC), no limit for dust 
!++wy
    if (.not.tot_in) then
        if ( do_bc ) frzbccnt = frzbccnt+MIN(limfacbc*uncoated_aer_num(id_bc)/deltat, &
                                               uncoated_aer_num(id_bc)/deltat &
                                               *(1.-exp(-Jcnt_bc*deltat))) 
        if ( do_dst1 ) frzducnt = frzducnt+MIN(1.*uncoated_aer_num(id_dst1)/deltat, &
                                               uncoated_aer_num(id_dst1)/deltat &
                                               *(1.-exp(-Jcnt_dust_a1*deltat)))
        if ( do_dst3 ) frzducnt = frzducnt+MIN(1.*uncoated_aer_num(id_dst3)/deltat, &
                                               uncoated_aer_num(id_dst3)/deltat &
                                               *(1.-exp(-Jcnt_dust_a3*deltat)))
    else
        if ( do_bc ) frzbccnt = frzbccnt+MIN(limfacbc*(1.-fn(id_bc))*(1.-dstcoat(1))*total_aer_num(id_bc)/deltat, &
                                               (1.-fn(id_bc))*(1.-dstcoat(1))*total_aer_num(id_bc)/deltat &
                                               *(1.-exp(-Jcnt_bc*deltat))) 
        if ( do_dst1 ) frzducnt = frzducnt+MIN(1*(1.-fn(id_dst1))*(1.-dstcoat(2))*total_aer_num(id_dst1)/deltat, &
                                               (1.-fn(id_dst1))*(1.-dstcoat(2))*total_aer_num(id_dst1)/deltat &
                                               *(1.-exp(-Jcnt_dust_a1*deltat)))
        if ( do_dst3 ) frzducnt = frzducnt+MIN(1*(1.-fn(id_dst3))*(1.-dstcoat(3))*total_aer_num(id_dst3)/deltat, &
                                               (1.-fn(id_dst3))*(1.-dstcoat(3))*total_aer_num(id_dst3)/deltat &
                                               *(1.-exp(-Jcnt_dust_a3*deltat)))
    end if
!--wy
    
!    if ( .not. ( frzducnt*frzducnt .ge. 0. ) ) PRINT*,'frzducnt', frzducnt, Jcnt_dust_a1, Jcnt_dust_a3, &
!                                                      Kcoll_dust_a1, Kcoll_dust_a3
!    if ( .not. ( frzducnt*frzducnt .ge. 0. ) ) STOP 'ERROR in frzducnt'
!    if ( frzducnt .le. -1. ) PRINT*,'frzducnt', frzducnt, Jcnt_dust_a1,Jcnt_dust_a3, &    
!                                    Kcoll_dust_a1, Kcoll_dust_a3
!    if ( frzducnt .le. -1. ) STOP 'ERROR in frzducnt'

   if (frzducnt <= -1._r8) then
      write(iulog,*) 'hetfrz_classnuc_calc: frzducnt', frzducnt, Jcnt_dust_a1,Jcnt_dust_a3, &    
                                    Kcoll_dust_a1, Kcoll_dust_a3
!      errstring = 'ERROR in hetfrz_classnuc_calc::frzducnt'
      return
   end if     

!++ wy4.0    
    contains
    !++3.1
    ! Approximation to the error function
    real(r8) function ERFAPP(x)
        real(r8),  intent(in) :: x
        real(r8) :: ax
        ax=x*x*(1.27324+(0.147*x*x))/(1.0+(0.147*x*x))
        if (x >= 0._r8) then
            ERFAPP=SQRT(1.0-exp(-ax))
        else
            ERFAPP=-1.0_r8*SQRT(1.0-exp(-ax))
        end if
    end function ERFAPP
    !--3.1
!-- wy4.0
    
end subroutine  hetfrz_classnuc

