module nucleate_ice

!---------------------------------------------------------------------------------
! Purpose:
!   Ice nucleation code.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use wv_saturation,  only: svp_water, svp_ice
use cam_logfile,    only: iulog
use physconst,      only: pi    ! sxj add  pi

implicit none
private
save

public :: nucleati

!===============================================================================
contains
!===============================================================================

subroutine nucleati(  &
   wbar, tair, relhum, cldn, qc,              &
   nfice, rhoair, so4_num, dst_num, soot_num, &
   nuci, onihf, oniimm, onidep, onimey,       &
   qi_pice, ni_pice, pmid, wpice, weff, fhom)        !sxj 

   !---------------------------------------------------------------
   ! Purpose:
   !  The parameterization of ice nucleation.
   !
   ! Method: The current method is based on Liu & Penner (2005)
   !  It related the ice nucleation with the aerosol number, temperature and the
   !  updraft velocity. It includes homogeneous freezing of sulfate, immersion
   !  freezing of soot, and Meyers et al. (1992) deposition nucleation
   !
   ! Authors: Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010
   !----------------------------------------------------------------

!++ classnuc wy
use classnuc, only: classnuc_in, preexisting_ice
!-- classnuc wy

   ! Input Arguments
   real(r8), intent(in) :: wbar        ! grid cell mean vertical velocity (m/s)
   real(r8), intent(in) :: tair        ! temperature (K)
   real(r8), intent(in) :: relhum      ! relative humidity with respective to liquid
   real(r8), intent(in) :: cldn        ! new value of cloud fraction    (fraction)
   real(r8), intent(in) :: qc          ! liquid water mixing ratio (kg/kg)
   real(r8), intent(in) :: nfice       ! ice mass fraction
   real(r8), intent(in) :: rhoair      ! air density (kg/m3)
   real(r8), intent(in) :: so4_num     ! so4 aerosol number (#/cm^3)
   real(r8), intent(in) :: dst_num     ! total dust aerosol number (#/cm^3)
   real(r8), intent(in) :: soot_num    ! soot (hydrophilic) aerosol number (#/cm^3)
!++sxj
   real(r8), optional, intent(in) :: qi_pice     ! grid-mean preexisting cloud ice mass mixing ratio (kg/kg)
   real(r8), optional, intent(in) :: ni_pice     ! grid-mean preexisting cloud ice number conc (#/kg) 
   real(r8), optional, intent(in) :: pmid        ! pressure at layer midpoints (pa)
!--sxj

   ! Output Arguments
!++sxj
   real(r8), optional, intent(out) :: wpice      ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at Shom
   real(r8), optional, intent(out) :: weff       ! effective Vertical velocity for ice nucleation (m/s); weff=wbar-wpice
   real(r8), optional, intent(out) :: fhom       ! how much fraction of cloud can reach Shom
!--sxj
   real(r8), intent(out) :: nuci       ! ice number nucleated (#/kg)
   real(r8), intent(out) :: onihf      ! nucleated number from homogeneous freezing of so4
   real(r8), intent(out) :: oniimm     ! nucleated number from immersion freezing
   real(r8), intent(out) :: onidep     ! nucleated number from deposition nucleation
   real(r8), intent(out) :: onimey     ! nucleated number from deposition nucleation  (meyers: mixed phase)

   ! Local workspace
!++sxj
   real(r8), parameter :: Shet   = 1.3_r8     ! het freezing threshold
   real(r8) :: wpicehet   ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at shet
   real(r8) :: weffhet    ! effective Vertical velocity for ice nucleation (m/s)  weff=wbar-wpicehet 
   real(r8), parameter :: rhoice = 0.5e3      ! kg/m3, Wpice is not sensitive to rhoice
   real(r8), parameter :: mincld = 0.0001_r8  ! fraction
   real(r8), parameter :: minweff= 0.001_r8   ! m/s
   real(r8), parameter :: gamma1=1.0_r8 
   real(r8), parameter :: gamma2=1.0_r8 
   real(r8), parameter :: gamma3=2.0_r8 
   real(r8), parameter :: gamma4=6.0_r8 
   real(r8), parameter :: ci = rhoice*pi/6._r8
!--sxj
   real(r8) :: nihf                      ! nucleated number from homogeneous freezing of so4
   real(r8) :: niimm                     ! nucleated number from immersion freezing
   real(r8) :: nidep                     ! nucleated number from deposition nucleation
   real(r8) :: nimey                     ! nucleated number from deposition nucleation (meyers)
   real(r8) :: n1, ni                    ! nucleated number
   real(r8) :: tc, A, B, C, regm, RHw    ! work variable
   real(r8) :: esl, esi, deles           ! work variable
   real(r8) :: subgrid
!++sxj  used for SUBROUTINE Vpreice
   real(r8) :: Ni_preice        ! cloud ice number conc (1/m3)   
   real(r8) :: lami,Ri_preice   ! mean cloud ice radius (m)
   real(r8) :: Shom             ! hom freezing saturation ratio (default 1.5)
   real(r8) :: detaT,RHimean    ! temperature standard deviation, mean cloudy RHi
!--sxj
   !-------------------------------------------------------------------------------

!++sxj
  if(preexisting_ice) then
    Ni_preice=ni_pice*rhoair               ! (convert from #/kg -> #/m3)
    Ni_preice=Ni_preice/max(mincld,cldn)   ! in cloud ice number density 
    if (Ni_preice.gt.10.0_r8) then    ! > 0.01/L   
       Shom=-1.5_r8                   ! hom Si, if Shom<1 , Shom will be recalculated in SUBROUTINE Vpreice, based on Shom
       lami=(gamma4*ci*ni_pice/qi_pice)**(1._r8/3._r8)
       Ri_preice=0.5_r8/lami          ! radius
       Ri_preice=max(Ri_preice,1e-8_r8)       ! >0.01micron
       call Vpreice(pmid,tair,Ri_preice,Ni_preice,Shom,wpice)
       call Vpreice(pmid,tair,Ri_preice,Ni_preice,Shet,wpicehet)
    else
       wpice=0.0_r8
       wpicehet=0.0_r8
    endif            
    weff=max(wbar-wpice,minweff)
    wpice=min(wpice,wbar)
    weffhet=max(wbar-wpicehet,minweff)
    wpicehet=min(wpicehet,wbar)

    detaT=wbar/0.23_r8
    RHimean=1.0_r8
    call frachom(tair,RHimean,detaT,fhom)
  end if
!--sxj

   ni = 0._r8
   tc = tair - 273.15_r8

   ! initialize
   niimm = 0._r8
   nidep = 0._r8
   nihf  = 0._r8

!++ sxj
if(preexisting_ice) then
   if(so4_num.ge.1.0e-10_r8 .and. (soot_num+dst_num).ge.1.0e-10_r8 .and. cldn.gt.0._r8) then

      !-----------------------------
      ! RHw parameterization for heterogeneous immersion nucleation
      A = 0.0073_r8
      B = 1.477_r8
      C = 131.74_r8
      RHw=(A*tc*tc+B*tc+C)*0.01_r8   ! RHi ~ 120-130%

      subgrid = 1.2_r8

      if((tc.le.-35.0_r8) .and. ((relhum*svp_water(tair)/svp_ice(tair)*subgrid).ge.1.2_r8)) then ! use higher RHi threshold

         A = -1.4938_r8 * log(soot_num+dst_num) + 12.884_r8
         B = -10.41_r8  * log(soot_num+dst_num) - 67.69_r8
         !regm = A * log(wbar) + B
         regm = A * log(weff) + B   ! sxj, wbar is replaced by weff

         if(tc.gt.regm) then    ! heterogeneous nucleation only
            !if(tc.lt.-40._r8 .and. wbar.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation
            if(tc.lt.-40._r8 .and. weff.gt.1._r8) then !sxj ::  wbar is replaced by weff
               !call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
               call hf(tc,weff,relhum,subgrid,so4_num,nihf) !sxj ::  wbar is replaced by weff
               niimm=0._r8
               nidep=0._r8
               if (nihf.gt.1e-3_r8) then ! sxj ::  hom occur,  add preexisting ice
                  niimm=min(dst_num,Ni_preice*1e-6_r8)       ! assuming dst_num freeze firstly
                  nihf=nihf + Ni_preice*1e-6_r8 - niimm      !   
               !else
               !   nihf  = 0._r8     
               endif
               nihf=nihf*fhom
               !n1=nihf   !sxj, add niimm
               n1=nihf+niimm 
            else
               !call hetero(tc,wbar,soot_num+dst_num,niimm,nidep)
               call hetero(tc,weffhet,soot_num+dst_num,niimm,nidep) !sxj ::  wbar is replaced by weffhet
               if (niimm.gt.1e-6_r8)then ! SXJ:: het freezing occur, add preexisting ice
                  niimm=niimm+Ni_preice*1e-6_r8
                  niimm=min(dst_num,niimm)        ! niimm < dst_num 
               !else
               !   niimm  = 0._r8 
               endif
               nihf=0._r8
               n1=niimm+nidep
            endif
         elseif (tc.lt.regm-5._r8) then ! homogeneous nucleation only
            !call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
            call hf(tc,weff,relhum,subgrid,so4_num,nihf) !sxj ::  wbar is replaced by weff
            niimm=0._r8
            nidep=0._r8
            if (nihf.gt.1e-3_r8) then ! sxj ::  hom occur,  add preexisting ice
               niimm=min(dst_num,Ni_preice*1e-6_r8)       ! assuming dst_num freeze firstly
               nihf=nihf + Ni_preice*1e-6_r8 - niimm      !  
            !else
            !   nihf  = 0._r8 
            endif
            nihf=nihf*fhom
            !n1=nihf   !sxj, add niimm
            n1=nihf+niimm 
         else        ! transition between homogeneous and heterogeneous: interpolate in-between
            !if(tc.lt.-40._r8 .and. wbar.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation
            if(tc.lt.-40._r8 .and. weff.gt.1._r8) then !sxj ::  wbar is replaced by weff
               !call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
               call hf(tc,weff,relhum,subgrid,so4_num,nihf)  !sxj ::  wbar is replaced by weff
               niimm=0._r8
               nidep=0._r8
               if (nihf.gt.1e-3_r8) then ! sxj ::  hom occur,  add preexisting ice
                  niimm=min(dst_num,Ni_preice*1e-6_r8)       ! assuming dst_num freeze firstly
                  nihf=nihf + Ni_preice*1e-6_r8 - niimm      !  
               !else
               !   nihf  = 0._r8
               endif
               nihf=nihf*fhom
               !n1=nihf   !sxj, add niimm
               n1=nihf+niimm 
            else

               !call hf(regm-5._r8,wbar,relhum,subgrid,so4_num,nihf)
               !call hetero(regm,wbar,soot_num+dst_num,niimm,nidep)
               call hf(regm-5._r8,weff,relhum,subgrid,so4_num,nihf)    !sxj ::  wbar is replaced by weff
               call hetero(regm,weffhet,soot_num+dst_num,niimm,nidep)  !sxj ::  wbar is replaced by weffhet
               nihf=nihf*fhom
               if(nihf.le.(niimm+nidep)) then
                  n1=nihf
               else
                  n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc-regm)/5._r8)
               endif
               if (n1.gt.1e-3_r8) then   !sxj::   add preexisting ice
                  n1=n1+Ni_preice*1e-6_r8                
                  niimm=min(dst_num,n1)  ! assuming all dst_num freezing earlier than hom  !!
                  nihf=n1-niimm 
               else
                  n1    = 0._r8
                  niimm = 0._r8
                  nihf  = 0._r8                         
               endif
            endif
         endif

         ni=n1

      endif
   endif
else
   if(so4_num.ge.1.0e-10_r8 .and. (soot_num+dst_num).ge.1.0e-10_r8 .and. cldn.gt.0._r8) then

      !-----------------------------
      ! RHw parameterization for heterogeneous immersion nucleation
      A = 0.0073_r8
      B = 1.477_r8
      C = 131.74_r8
      RHw=(A*tc*tc+B*tc+C)*0.01_r8   ! RHi ~ 120-130%

      subgrid = 1.2_r8

      if((tc.le.-35.0_r8) .and. ((relhum*svp_water(tair)/svp_ice(tair)*subgrid).ge.1.2_r8)) then ! use higher RHi threshold

         A = -1.4938_r8 * log(soot_num+dst_num) + 12.884_r8
         B = -10.41_r8  * log(soot_num+dst_num) - 67.69_r8
         regm = A * log(wbar) + B

         if(tc.gt.regm) then    ! heterogeneous nucleation only
            if(tc.lt.-40._r8 .and. wbar.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation
               call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
               niimm=0._r8
               nidep=0._r8
               n1=nihf
            else
               call hetero(tc,wbar,soot_num+dst_num,niimm,nidep)
               nihf=0._r8
               n1=niimm+nidep
            endif
         elseif (tc.lt.regm-5._r8) then ! homogeneous nucleation only
            call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
            niimm=0._r8
            nidep=0._r8
            n1=nihf
         else        ! transition between homogeneous and heterogeneous: interpolate in-between
            if(tc.lt.-40._r8 .and. wbar.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation
               call hf(tc,wbar,relhum,subgrid,so4_num,nihf)
               niimm=0._r8
               nidep=0._r8
               n1=nihf
            else

               call hf(regm-5._r8,wbar,relhum,subgrid,so4_num,nihf)
               call hetero(regm,wbar,soot_num+dst_num,niimm,nidep)

               if(nihf.le.(niimm+nidep)) then
                  n1=nihf
               else
                  n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc-regm)/5._r8)
               endif
            endif
         endif

         ni=n1

      endif
   endif
end if
!-- sxj

   ! deposition/condensation nucleation in mixed clouds (-40<T<0C) (Meyers, 1992)
   if(tc.lt.0._r8 .and. tc.gt.-37._r8 .and. qc.gt.1.e-12_r8) then
      esl = svp_water(tair)     ! over water in mixed clouds
      esi = svp_ice(tair)     ! over ice
      deles = (esl - esi)
      nimey=1.e-3_r8*exp(12.96_r8*deles/esi - 0.639_r8) 
   else
      nimey=0._r8
   endif
   !++ classnuc wy
   if(classnuc_in) nimey = 0._r8
   !-- classnuc wy
   nuci=ni+nimey
   if(nuci.gt.9999._r8.or.nuci.lt.0._r8) then
      write(iulog, *) 'Warning: incorrect ice nucleation number (nuci reset =0)'
      write(iulog, *) ni, tair, relhum, wbar, nihf, niimm, nidep,deles,esi,dst_num,so4_num
      nuci=0._r8
   endif

   nuci   = nuci*1.e+6_r8/rhoair    ! change unit from #/cm3 to #/kg
   onimey = nimey*1.e+6_r8/rhoair
   onidep = nidep*1.e+6_r8/rhoair
   oniimm = niimm*1.e+6_r8/rhoair
   onihf  = nihf*1.e+6_r8/rhoair

end subroutine nucleati

!===============================================================================

subroutine hetero(T,ww,Ns,Nis,Nid)

    real(r8), intent(in)  :: T, ww, Ns
    real(r8), intent(out) :: Nis, Nid

    real(r8) A11,A12,A21,A22,B11,B12,B21,B22
    real(r8) A,B,C

!---------------------------------------------------------------------
! parameters

      A11 = 0.0263_r8
      A12 = -0.0185_r8
      A21 = 2.758_r8
      A22 = 1.3221_r8
      B11 = -0.008_r8
      B12 = -0.0468_r8
      B21 = -0.2667_r8
      B22 = -1.4588_r8

!     ice from immersion nucleation (cm^-3)

      B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
      C =  A21+B21*log(Ns)

      Nis = exp(A22) * Ns**B22 * exp(B*T) * ww**C
      Nis = min(Nis,Ns)

      Nid = 0.0_r8    ! don't include deposition nucleation for cirrus clouds when T<-37C

end subroutine hetero

!===============================================================================

subroutine hf(T,ww,RH,subgrid,Na,Ni)

      real(r8), intent(in)  :: T, ww, RH, subgrid, Na
      real(r8), intent(out) :: Ni

      real(r8)    A1_fast,A21_fast,A22_fast,B1_fast,B21_fast,B22_fast
      real(r8)    A2_fast,B2_fast
      real(r8)    C1_fast,C2_fast,k1_fast,k2_fast
      real(r8)    A1_slow,A2_slow,B1_slow,B2_slow,B3_slow
      real(r8)    C1_slow,C2_slow,k1_slow,k2_slow
      real(r8)    regm
      real(r8)    A,B,C
      real(r8)    RHw

!---------------------------------------------------------------------
! parameters

      A1_fast  =0.0231_r8
      A21_fast =-1.6387_r8  !(T>-64 deg)
      A22_fast =-6.045_r8   !(T<=-64 deg)
      B1_fast  =-0.008_r8
      B21_fast =-0.042_r8   !(T>-64 deg)
      B22_fast =-0.112_r8   !(T<=-64 deg)
      C1_fast  =0.0739_r8
      C2_fast  =1.2372_r8

      A1_slow  =-0.3949_r8
      A2_slow  =1.282_r8
      B1_slow  =-0.0156_r8
      B2_slow  =0.0111_r8
      B3_slow  =0.0217_r8
      C1_slow  =0.120_r8
      C2_slow  =2.312_r8

      Ni = 0.0_r8

!----------------------------
!RHw parameters
      A = 6.0e-4_r8*log(ww)+6.6e-3_r8
      B = 6.0e-2_r8*log(ww)+1.052_r8
      C = 1.68_r8  *log(ww)+129.35_r8
      RHw=(A*T*T+B*T+C)*0.01_r8

      if((T.le.-37.0_r8) .and. ((RH*subgrid).ge.RHw)) then

        regm = 6.07_r8*log(ww)-55.0_r8

        if(T.ge.regm) then    ! fast-growth regime

          if(T.gt.-64.0_r8) then
            A2_fast=A21_fast
            B2_fast=B21_fast
          else
            A2_fast=A22_fast
            B2_fast=B22_fast
          endif

          k1_fast = exp(A2_fast + B2_fast*T + C2_fast*log(ww))
          k2_fast = A1_fast+B1_fast*T+C1_fast*log(ww)

          Ni = k1_fast*Na**(k2_fast)
          Ni = min(Ni,Na)

        else       ! slow-growth regime

          k1_slow = exp(A2_slow + (B2_slow+B3_slow*log(ww))*T + C2_slow*log(ww))
          k2_slow = A1_slow+B1_slow*T+C1_slow*log(ww)

          Ni = k1_slow*Na**(k2_slow)
          Ni = min(Ni,Na)

        endif

      end if

end subroutine hf

!===============================================================================
!++sxj  based on KL code 2006
!  VERTICAL VELOCITY CALCULATED FROM DEPOSITIONAL LOSS TERM
  SUBROUTINE Vpreice(P_in,T_in,R_in,C_in,S_in,V_out)
     IMPLICIT NONE
     ! SUBROUTINE variables
     REAL(r8), INTENT(in)  :: P_in       ! [Pa],INITIAL AIR pressure 
     REAL(r8), INTENT(in)  :: T_in       ! [K] ,INITIAL AIR temperature 
     REAL(r8), INTENT(in)  :: R_in       ! [m],INITIAL MEAN  ICE CRYSTAL NUMBER RADIUS 
     REAL(r8), INTENT(in)  :: C_in       ! [m-3],INITIAL TOTAL ICE CRYSTAL NUMBER DENSITY, [1/cm3]
     REAL(r8), INTENT(in)  :: S_in       ! [-],INITIAL ICE SATURATION RATIO;; if <1, use hom threshold Si 
     REAL(r8), INTENT(out) :: V_out      ! [m/s],EFFECTIVE VERTICAL VELOCITY 
     ! SUBROUTINE parameters
     REAL(r8), PARAMETER :: ALPHAc  = 0.5  ! density of ice (g/cm3), !!!V is not related to ALPHAc 
     REAL(r8), PARAMETER :: FA1c    = 0.601272523_r8        
     REAL(r8), PARAMETER :: FA2c    = 0.000342181855_r8
     REAL(r8), PARAMETER :: FA3c    = 1.49236645E-12_r8        
     REAL(r8), PARAMETER :: WVP1c   = 3.6E+10_r8   
     REAL(r8), PARAMETER :: WVP2c   = 6145.0_r8
     REAL(r8), PARAMETER :: FVTHc   = 11713803.0_r8
     REAL(r8), PARAMETER :: THOUBKc = 7.24637701E+18_r8
     REAL(r8), PARAMETER :: SVOLc   = 3.23E-23_r8    ! SVOL=XMW/RHOICE
     REAL(r8), PARAMETER :: FDc     = 249.239822    
     REAL(r8), PARAMETER :: FPIVOLc = 3.89051704E+23_r8         
     REAL(r8) :: T,P,S,R,C
     REAL(r8) :: A1,A2,A3,B1,B2
     REAL(r8) :: T_1,PICE,FLUX,ALP4,CISAT,DLOSS,VICE
     T=T_in          ! K  , K
     P=P_in*1e-2_r8  ! Pa , hpa
     IF (S_in.LT.1.0_r8) THEN
        S=2.349_r8-(T/259.0_r8) !threshold Si according to Ren & McKenzie, 2005
     ELSE
        S=S_in                  ! ICE SATURATION RATIO, -,  >1
     ENDIF
     R=R_in*1e2_r8   ! m  , cm
     C=C_in*1e-6_r8  ! m-3, cm-3
     T_1   = 1.0_r8/ T
     PICE  = WVP1c * EXP(-(WVP2c*T_1))
     ALP4  = 0.25_r8 * ALPHAc      
     FLUX  = ALP4 * SQRT(FVTHc*T)
     CISAT = THOUBKc * PICE * T_1   
     A1    = ( FA1c * T_1 - FA2c ) * T_1 
     A2    = 1.0_r8/ CISAT      
     A3    = FA3c * T_1 / P
     B1    = FLUX * SVOLc * CISAT * ( S-1.0_r8 ) 
     B2    = FLUX * FDc * P * T_1**1.94_r8 
     DLOSS = FPIVOLc * C * B1 * R**2 / ( 1.0_r8+ B2 * R )         
     VICE  = ( A2 + A3 * S ) * DLOSS / ( A1 * S )  ! 2006,(19)
     V_out = VICE*1e-2_r8  ! m/s , cm/s
  END SUBROUTINE Vpreice

! How much fraction of cirrus might reach Shom
! base on "A cirrus cloud scheme for general circulation models" B. Ka Ì<88>rcher and U. Burkhardt 2008
        subroutine frachom(Tmean,RHimean,detaT,fhom)
           implicit none
           real(r8),parameter ::  seta = 6132.9_r8  ! K
           integer,parameter :: Nbin=200            ! -3 ~ 3
           real(r8),intent(in)  :: Tmean,RHimean,detaT
           real(r8),intent(out) :: fhom
           real(r8) :: Tbin(Nbin),PDF_T(Nbin),Sbin(Nbin)
           real(r8) :: Sihom,deta
           integer i
           Sihom=2.349_r8-Tmean/259.0_r8
           fhom=0.0_r8
           do i=Nbin,1,-1
              deta=(i-0.5_r8-Nbin/2)*6.0_r8/Nbin   ! -3 ~ 3
              Tbin(i)=Tmean-deta*detaT
              Sbin(i)=RHimean*exp(deta*detaT*seta/Tmean**2.0_r8)
              PDF_T(i)=exp(-deta**2.0_r8/2.0_r8)*6.0_r8/(sqrt(2.0_r8*Pi)*Nbin)
              if (Sbin(i).ge.Sihom) then
                 fhom=fhom+PDF_T(i)
              else
                 exit
              endif
           enddo
           fhom=fhom/0.997_r8   !-3~3
        endsubroutine frachom
!--sxj

end module nucleate_ice

