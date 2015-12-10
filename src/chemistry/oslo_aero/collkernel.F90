!-----------------------------------------------------------------------
!
! Purpose: calculate collision kernels as a function of environmental parameters and aerosol/droplet sizes
!
! Public interfaces:
!
! collkernel
!
! Author: Corinna Hoose, UiO, October 2009
!
! Modifications: Yong Wang
!-----------------------------------------------------------------------

      subroutine collkernel(t, pres, eswtr, rhwincloud, r3lx, &
                      !++wy
                      r_bc,                                   &  ! BC modes
                      r_dust_a1, r_dust_a3,                   &  ! dust modes
                      Kcoll_bc,                               &  ! collision kernel [cm3 s-1]
                      Kcoll_dust_a1, Kcoll_dust_a3)
                      !--wy

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use physconst,     only: rair, cpair, rh2o,    &
                                 rhoh2o, tmelt, pi
      use abortutils,    only: endrun

      implicit none

!     input
      real(r8), intent(in) :: t                ! temperature [K]
      real(r8), intent(in) :: pres             ! pressure [Pa]
      real(r8), intent(in) :: eswtr            ! saturation vapor pressure of water [Pa]
      real(r8), intent(in) :: r3lx             ! volume mean drop radius [m]
      real(r8), intent(in) :: rhwincloud       ! in-cloud relative humidity over water [ ]
      real(r8), intent(in) :: r_bc             ! model radii of BC modes [m]
      real(r8), intent(in) :: r_dust_a1, r_dust_a3 ! model radii of dust modes [m]
!     output
      real(r8), intent(out) :: Kcoll_bc, &      ! collision kernel [cm3 s-1]
                               Kcoll_dust_a1, Kcoll_dust_a3                              
!     local variables
      real(r8) :: a,b,c,a_f,b_f,c_f,f
      real(r8) :: tc          ! temperature [deg C]
      real(r8) :: rho_air     ! air density [kg m-3]
      real(r8) :: viscos_air  ! dynamic viscosity of air [kg m-1 s-1]
      real(r8) :: Ktherm_air  ! thermal conductivity of air [J/(m s K)]
      real(r8) :: lambda      ! mean free path [m]
      real(r8) :: Kn          ! Knudsen number [ ]
      real(r8) :: Re          ! Reynolds number [ ]
      real(r8) :: Pr          ! Prandtl number [ ]
      real(r8) :: Sc          ! Schmidt number [ ]
      real(r8) :: vterm       ! terminal velocity [m s-1]
      real(r8) :: Ktherm      ! thermal conductivity of aerosol [J/(m s K)]
      real(r8) :: Dvap        ! water vapor diffusivity [m2 s-1]
      real(r8) :: Daer        ! aerosol diffusivity [m2 s-1]
      real(r8) :: latvap      ! latent heat of vaporization [J kg-1]
      real(r8) :: kboltz      ! Boltzmann constant [J K-1]
      real(r8) :: G           ! thermodynamic function in Cotton et al. [kg m-1 s-1]
      real(r8) :: r_a         ! aerosol radius [m]
      real(r8) :: f_t         ! factor by Waldmann & Schmidt [ ]
      real(r8) :: Q_heat      ! heat flux [J m-2 s-1]
      real(r8) :: Tdiff_cotton ! temperature difference between droplet and environment [K]
      real(r8) :: K_brownian,K_thermo_cotton,K_diffusio_cotton   ! collision kernels [m3 s-1]
      real(r8) :: K_total     ! total collision kernel [cm3 s-1]
      integer  :: i
        
      Kcoll_bc = 0.
      Kcoll_dust_a1 = 0.
      Kcoll_dust_a3 = 0.

      tc = t-tmelt
      kboltz = 1.38065e-23

      !air viscosity !for tc<0, from depvel_part.F90
      viscos_air = (1.718+0.0049*tc-1.2e-5*tc*tc)*1.e-5
      !air density
      rho_air = pres/(rair*t)
      !mean free path: Seinfeld & Pandis 8.6
      lambda = 2*viscos_air/(pres*SQRT(8/(pi*rair*t)))
      !latent heat of vaporization, varies with T
      latvap = 1000*(-0.0000614342*tc**3 + 0.00158927*tc**2 - 2.36418*tc + 2500.79)
      !droplet terminal velocity after Chen & Liu, QJRMS 2004
      a = 8.8462e2
      b = 9.7593e7
      c = -3.4249e-11
      a_f = 3.1250e-1
      b_f = 1.0552e-3
      c_f = -2.4023
!     f = EXP(EXP(a_f+b_f*(ALOG(r3lx))**3+c_f*rho_air**1.5))
      f = EXP(EXP(a_f+b_f*(LOG(r3lx))**3+c_f*rho_air**1.5))
      vterm = (a+(b+c*r3lx)*r3lx)*r3lx*f

      !Reynolds number
      Re = 2*vterm*r3lx*rho_air/viscos_air
      !thermal conductivity of air: Seinfeld & Pandis eq. 15.75
      Ktherm_air = 1.e-3*(4.39+0.071*t)  !J/(m s K)
      !Prandtl number
      Pr = viscos_air*cpair/Ktherm_air
      !water vapor diffusivity: Pruppacher & Klett 13-3
      Dvap = 0.211e-4*(t/273.15)*(101325./pres) 
      !G-factor = rhoh2o*Xi in Rogers & Yau, p. 104
      G = rhoh2o/((latvap/(rh2o*t)-1)*latvap*rhoh2o/(Ktherm_air*t) &
                 +rhoh2o*rh2o*t/(Dvap*eswtr))
      
      !variables depending on aerosol radius
      do i = 1,3 !loop over 3 aerosol modes
         if (i .eq. 1) r_a = r_bc
         if (i .eq. 2) r_a = r_dust_a1
         if (i .eq. 3) r_a = r_dust_a3
         !Knudsen number (Seinfeld & Pandis 8.1)
         Kn = lambda/r_a
         !aerosol diffusivity
         Daer = kboltz*t*(1+Kn)/(6*pi*r_a*viscos_air)
         !Schmidt number
         Sc = viscos_air/(Daer*rho_air)

         !Young (1974) first equ. on page 771
         K_brownian = 4*pi*r3lx*Daer*(1+0.3*Re**0.5*Sc**0.33)

         !thermal conductivities from Seinfeld & Pandis, Table 8.6
         if (i .eq. 1) Ktherm = 4.2 !Carbon
         if (i .eq. 2 .or. i .eq. 3) Ktherm = 0.72 !clay
         !form factor
         f_t = 0.4*(1.+1.45*Kn+0.4*Kn*EXP(-1./Kn))                   &
                  *(Ktherm_air+2.5*Kn*Ktherm)                        &
               /((1.+3.*Kn)*(2.*Ktherm_air+5.*Kn*Ktherm+Ktherm))
         !calculate T-Tc as in Cotton et al.
         Tdiff_cotton = -G*(rhwincloud-1)*latvap/Ktherm_air
         Q_heat = Ktherm_air/r3lx*(1+0.3*Re**0.5*Pr**0.33)*Tdiff_cotton
         K_thermo_cotton = 4*pi*r3lx*r3lx*f_t*Q_heat/pres
         K_diffusio_cotton = -(1./f_t)*(rh2o*t/latvap)*K_thermo_cotton
         K_total = 1.e6*(K_brownian+K_thermo_cotton+K_diffusio_cotton)  !conversion m3/s -> cm3/s
         !set K to 0 if negative
         if (K_total .lt. 0.) K_total = 0.

         if (i .eq. 1) Kcoll_bc = K_total
         if (i .eq. 2) Kcoll_dust_a1 = K_total
         if (i .eq. 3) Kcoll_dust_a3 = K_total
        
      end do
      

      end subroutine collkernel
