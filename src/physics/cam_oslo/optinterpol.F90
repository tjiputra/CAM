module optinterpol

! Purpose: To interpolate between look-up table entries for SW optical aerosol properties.

!     Optimized for speed by Arild Burud and Egil Storen (NoSerC), June-July 2002
!--------------------------------------------------------------------------------

! Updated for new kcomp1.out including condensed SOA - Alf Kirkevaag, May 2013.
! Extended for new SOA treatment for  kcomp1-4.out and treating SOA as coagulated OC
! for kcomp5-10 - Alf Kirkevaag, August 2015, and also rewritten to a more generalized
! for for interpolations using common subroutines interpol*dim.

  use shr_kind_mod, only: r8 => shr_kind_r8
  use opttab
  use opttab_lw
  use commondefinitions, only: nmodes, nbmodes
  implicit none

  private 
  save

  public interpol0
  public interpol1
  public interpol2to3
  public interpol4
  public interpol5to10  

 contains

!********************************************************************************************

subroutine interpol0 (lchnk, ncol, daylight, Nnatk, omega, gass, bex, ske, lw_on, kabs)

   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none


!
! Input arguments
!
   integer, intent(in) :: lchnk                       ! chunk identifier
   integer, intent(in) :: ncol                        ! number of atmospheric columns
   logical, intent(in) :: daylight(pcols)             ! calculations also at (polar) night if daylight=.true.
   logical, intent(in) :: lw_on                       ! LW calculations are performed if true
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
!
! Output arguments
!
   real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
   real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
   real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
   real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
   real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absorption coefficient
!
!---------------------------Local variables-----------------------------
!
      integer i, ierr, kcomp, k, icol


      kcomp=0

        do i=1,nbands
          do icol=1,ncol
            do k=1,pver
              omega(icol,k,kcomp,i)=0.0_r8
              gass(icol,k,kcomp,i)=0.0_r8
              bex(icol,k,kcomp,i)=0.0_r8
              ske(icol,k,kcomp,i)=0.0_r8
            end do
          end do
        end do
        do i=1,nlwbands
          do icol=1,ncol
            do k=1,pver
              kabs(icol,k,kcomp,i)=0.0_r8
            end do
          end do
        end do
         
!      SW optical parameters

        do k=1,pver
          do icol=1,ncol

!           if(Nnatk(icol,k,kcomp)>0.0_r8) then
           if(daylight(icol)) then
           do i=1,nbands   ! i = wavelength index
              omega(icol,k,kcomp,i)=om0(i)
              gass(icol,k,kcomp,i)=g0(i) 
              bex(icol,k,kcomp,i)=be0(i)
              ske(icol,k,kcomp,i)=ke0(i)
           end do          ! i
           endif  ! daylight        

          end do ! icol
        end do ! k

       if(lw_on) then

!      LW optical parameters

        do k=1,pver
          do icol=1,ncol

            do i=1,nlwbands   ! i = wavelength index
               kabs(icol,k,kcomp,i)=ka0(i)
            end do            ! i

          end do ! icol
        end do ! k

       endif ! lw_on

      return
end subroutine interpol0


!********************************************************************************************

subroutine interpol1 (lchnk, ncol, daylight, xrh, irh1, irh2, &
                      mplus10, Nnatk, xfombgin, &
                      Camk, xfacin, omega, gass, bex, ske, lw_on, kabs, cxstot)


   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none

!
! Input arguments
!
   integer, intent(in) :: lchnk                     ! chunk identifier
   integer, intent(in) :: ncol                      ! number of atmospheric columns
   integer, intent(in) :: mplus10                   ! mode number (0) or number + 10 (1)
   logical, intent(in) :: daylight(pcols)           ! only daylight calculations if .true.
   logical, intent(in) :: lw_on                       ! LW calculations are performed if true
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes)   ! modal internally mixed (cond+aq) SO4 conc.
   real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
   integer,  intent(in) :: irh1(pcols,pver)
   integer,  intent(in) :: irh2(pcols,pver)
   real(r8), intent(in) :: xfombgin(pcols,pver)     ! SOA/(SOA+H2SO4) for the background mode 
   real(r8), intent(in) :: xfacin(pcols,pver,4)     ! SOA/(SOA+H2SO4) for condensated mass 
!
!
! Input-Output arguments
!
   real(r8), intent(inout) :: cxstot(pcols,pver)    ! excess internally mixed mass ("table overshoot") 
!
!
! Output arguments
!
   real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
   real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
   real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
   real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
   real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absoption coefficient
!
!---------------------------Local variables-----------------------------
!
      integer i, ierr, irelh, ictot, ifac, ifombg, kcomp, k, icol, kc10
!      integer irh1(pcols,pver), irh2(pcols,pver), ict1(pcols,pver), &
      integer ict1(pcols,pver), &
        ict2(pcols,pver), ifac1(pcols,pver), ifac2(pcols,pver), &
        ifombg1(pcols,pver), ifombg2(pcols,pver)
      real(r8) a, b, camdiff
!      real(r8) xrh(pcols,pver), xct(pcols,pver), xfac(pcols,pver), &
      real(r8) xct(pcols,pver), xfac(pcols,pver), &
        xfombg(pcols,pver), cxs(pcols,pver)
 
!     Temporary storage of often used array elements
      integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2, t_ifo1, t_ifo2
      real(r8) t_fac1, t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, &
        t_cat1, t_cat2, t_fombg1, t_fombg2, t_xfombg
      real(r8) d2mx(4), dxm1(4), invd(4)
      real(r8) opt4d(2,2,2,2)
      real(r8) ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2  
      real(r8) kabs1, kabs2

 
      do k=1,pver
        do icol=1,ncol
          xct(icol,k)    = 0.0_r8
          xfac(icol,k)   = 0.0_r8
          xfombg(icol,k) = 0.0_r8
        end do
      end do

!      write(*,*) 'Before kcomp-loop'
        do kcomp=1,1

           if(mplus10==0) then
             kc10=kcomp
           else
             kc10=kcomp+10
           endif 

      do k=1,pver
        do icol=1,ncol

          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
                 /(Nnatk(icol,k,kc10)+eps),cate(kcomp,1)),cate(kcomp,16))
          xfac(icol,k) = min(max(xfacin(icol,k,kcomp),fac(1)),fac(6))
          xfombg(icol,k) = min(max(xfombgin(icol,k),fombg(1)),fombg(6))
          camdiff=Camk(icol,k,kcomp)-xct(icol,k) &
                         *(Nnatk(icol,k,kc10)+eps)
          cxs(icol,k)=max(0.0_r8,camdiff)
          cxstot(icol,k)=cxstot(icol,k)+cxs(icol,k)

!        if(icol.eq.1) then
!          write(*,*) 'dir: k, kc10, cxs =', k, kc10, cxs(icol,k)

!         write(*,*) 'Before cate-loop', kc10
          do ictot=1,15
            if(xct(icol,k)>=cate(kcomp,ictot).and. &
            xct(icol,k) <= cate(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
          end do ! ictot

!         write(*,*) 'Before fac-loop', kcomp
          do ifac=1,5
            if(xfac(icol,k)>=fac(ifac).and. &
             xfac(icol,k) <= fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
          end do ! ifac

!         write(*,*) 'Before fombg-loop', kcomp
          do ifombg=1,5
            if(xfombg(icol,k) >= fombg(ifombg).and. &
            xfombg(icol,k) <= fombg(ifombg+1)) then
             ifombg1(icol,k)=ifombg
             ifombg2(icol,k)=ifombg+1
            endif
          end do ! ifombg

        end do ! icol
      end do ! k


!      write(*,*) 'Before init-loop', kc10
        do i=1,nbands
          do icol=1,ncol
            do k=1,pver
              omega(icol,k,kc10,i)=0.0_r8
              gass(icol,k,kc10,i)=0.0_r8
              bex(icol,k,kc10,i)=0.0_r8
              ske(icol,k,kc10,i)=0.0_r8
            end do
          end do
        end do
        do i=1,nlwbands
          do icol=1,ncol
            do k=1,pver
              kabs(icol,k,kc10,i)=0.0_r8
             end do
          end do
        end do
         
        do k=1,pver
          do icol=1,ncol

!      Collect all the vector elements into temporary storage
!      to avoid cache conflicts and excessive cross-referencing

      t_irh1 = irh1(icol,k)
      t_irh2 = irh2(icol,k)
      t_ict1 = ict1(icol,k)
      t_ict2 = ict2(icol,k)
      t_ifc1 = ifac1(icol,k)
      t_ifc2 = ifac2(icol,k)
      t_ifo1 = ifombg1(icol,k)
      t_ifo2 = ifombg2(icol,k)

!      write(*,*) 't_irh1,t_irh2=',t_irh1,t_irh2
!      write(*,*) 't_ict1,t_ict2=',t_ict1,t_ict2
!      write(*,*) 't_ifc1,t_ifc2=',t_ifc1,t_ifc2
!      write(*,*) 't_ifo1,t_ifo2=',t_ifo1,t_ifo2

      t_rh1  = rh(t_irh1)
      t_rh2  = rh(t_irh2)
      t_cat1 = cate(kcomp,t_ict1)
      t_cat2 = cate(kcomp,t_ict2)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)
      t_fombg1 = fombg(t_ifo1)
      t_fombg2 = fombg(t_ifo2)

!      write(*,*) 't_rh1,t_rh2,t_cat1,t_cat2=',t_rh1,t_rh2,t_cat1,t_cat2
!      write(*,*) 't_fac1,t_fac2=',t_fac1,t_fac2

      t_xrh  = xrh(icol,k)
      t_xct  = xct(icol,k)
      t_xfac = xfac(icol,k)
      t_xfombg = xfombg(icol,k)

!     partial lengths along each dimension (1-4) for interpolation 
      d2mx(1) = (t_rh2-t_xrh)
      dxm1(1) = (t_xrh-t_rh1)
      invd(1) = 1.0_r8/(t_rh2-t_rh1)
      d2mx(2) = (t_fombg2-t_xfombg)
      dxm1(2) = (t_xfombg-t_fombg1)
      invd(2) = 1.0_r8/(t_fombg2-t_fombg1)
      d2mx(3) = (t_cat2-t_xct)
      dxm1(3) = (t_xct-t_cat1)
      invd(3) = 1.0_r8/(t_cat2-t_cat1)
      d2mx(4) = (t_fac2-t_xfac)
      dxm1(4) = (t_xfac-t_fac1)
      invd(4) = 1.0_r8/(t_fac2-t_fac1)


!      SW optical parameters
       if(daylight(icol)) then

      do i=1,nbands            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  single scattering albedo:

!     end points as basis for multidimentional linear interpolation  
      opt4d(1,1,1,1)=om1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
      opt4d(1,1,1,2)=om1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
      opt4d(1,1,2,1)=om1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
      opt4d(1,1,2,2)=om1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
      opt4d(1,2,1,1)=om1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
      opt4d(1,2,1,2)=om1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
      opt4d(1,2,2,1)=om1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
      opt4d(1,2,2,2)=om1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
      opt4d(2,1,1,1)=om1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
      opt4d(2,1,1,2)=om1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
      opt4d(2,1,2,1)=om1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
      opt4d(2,1,2,2)=om1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
      opt4d(2,2,1,1)=om1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
      opt4d(2,2,1,2)=om1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
      opt4d(2,2,2,1)=om1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
      opt4d(2,2,2,2)=om1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

!     interpolation in the fac, cat and fombg dimensions
      call lininterpol4dim (d2mx, dxm1, invd, opt4d, ome1, ome2)

!     finally, interpolation in the rh dimension
!      write(*,*) 'Before omega'
      omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) &
                          /(t_rh2-t_rh1)    
!alt       omega(icol,k,kc10,i)=(d2mx(1)*ome1+dxm1(1)*ome2)*invd(1)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  asymmetry factor   

!     end points as basis for multidimentional linear interpolation  
      opt4d(1,1,1,1)=g1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
      opt4d(1,1,1,2)=g1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
      opt4d(1,1,2,1)=g1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
      opt4d(1,1,2,2)=g1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
      opt4d(1,2,1,1)=g1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
      opt4d(1,2,1,2)=g1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
      opt4d(1,2,2,1)=g1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
      opt4d(1,2,2,2)=g1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
      opt4d(2,1,1,1)=g1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
      opt4d(2,1,1,2)=g1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
      opt4d(2,1,2,1)=g1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
      opt4d(2,1,2,2)=g1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
      opt4d(2,2,1,1)=g1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
      opt4d(2,2,1,2)=g1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
      opt4d(2,2,2,1)=g1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
      opt4d(2,2,2,2)=g1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

!     interpolation in the fac, cat and fombg dimensions
      call lininterpol4dim (d2mx, dxm1, invd, opt4d, ge1, ge2)

!     finally, interpolation in the rh dimension (dim. 1)
!      write(*,*) 'Before gass'
      gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                   /(t_rh2-t_rh1)    
!alt      gass(icol,k,kc10,i)=(d2mx(1)*ge1+dxm1(1)*ge2)*invd(1)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol extinction   

!     end points as basis for multidimentional linear interpolation  
      opt4d(1,1,1,1)=be1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
      opt4d(1,1,1,2)=be1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
      opt4d(1,1,2,1)=be1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
      opt4d(1,1,2,2)=be1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
      opt4d(1,2,1,1)=be1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
      opt4d(1,2,1,2)=be1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
      opt4d(1,2,2,1)=be1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
      opt4d(1,2,2,2)=be1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
      opt4d(2,1,1,1)=be1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
      opt4d(2,1,1,2)=be1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
      opt4d(2,1,2,1)=be1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
      opt4d(2,1,2,2)=be1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
      opt4d(2,2,1,1)=be1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
      opt4d(2,2,1,2)=be1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
      opt4d(2,2,2,1)=be1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
      opt4d(2,2,2,2)=be1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

!     interpolation in the fac, cat and fombg dimensions
      call lininterpol4dim (d2mx, dxm1, invd, opt4d, bex1, bex2)

      bex1=max(bex1,1.e-30_r8)
      bex2=max(bex2,1.e-30_r8)

!     finally, interpolation in the rh dimension
!      write(*,*) 'Before bex'
      if(t_xrh <= 0.37_r8) then
        bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
            /(t_rh2-t_rh1)    
!alt        bex(icol,k,kc10,i)=(d2mx(1)*bex1+dxm1(1)*bex2)*invd(1)
      else
        a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
        b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
        bex(icol,k,kc10,i)=e**(a*t_xrh+b)
!alt        a=(log(bex2)-log(bex1))*invd(1)
!alt        b=(t_rh2*log(bex1)-t_rh1*log(bex2))*invd(1)
!alt        bex(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i

!      if(bex(icol,k,kc10,8)<1.e-20_r8) then
!        write(*,995) 'bex(8)=', kc10, t_xrh, t_xct, t_xfac, t_xfombg, bex(icol,k,kc10,8)
!      endif


      do i=4,4            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol specific extinction   

!     end points as basis for multidimentional linear interpolation  
      opt4d(1,1,1,1)=ke1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
      opt4d(1,1,1,2)=ke1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
      opt4d(1,1,2,1)=ke1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
      opt4d(1,1,2,2)=ke1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
      opt4d(1,2,1,1)=ke1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
      opt4d(1,2,1,2)=ke1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
      opt4d(1,2,2,1)=ke1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
      opt4d(1,2,2,2)=ke1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
      opt4d(2,1,1,1)=ke1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
      opt4d(2,1,1,2)=ke1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
      opt4d(2,1,2,1)=ke1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
      opt4d(2,1,2,2)=ke1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
      opt4d(2,2,1,1)=ke1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
      opt4d(2,2,1,2)=ke1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
      opt4d(2,2,2,1)=ke1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
      opt4d(2,2,2,2)=ke1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

!     interpolation in the fac, cat and fombg dimensions
      call lininterpol4dim (d2mx, dxm1, invd, opt4d, ske1, ske2)

      ske1=max(ske1,1.e-30_r8)
      ske2=max(ske2,1.e-30_r8)

!     finally, interpolation in the rh dimension
!      write(*,*) 'Before ske'
      if(t_xrh <= 0.37_r8) then
        ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
            /(t_rh2-t_rh1)    
!alt        ske(icol,k,kc10,i)=(d2mx(1)*ske1+dxm1(1)*ske2)*invd(1)
      else
        a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
        b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
        ske(icol,k,kc10,i)=e**(a*t_xrh+b)
!alt        a=(log(ske2)-log(ske1))*invd(1)
!alt        b=(t_rh2*log(ske1)-t_rh1*log(ske2))*invd(1)
!alt        ske(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i

       endif  ! daylight

      if (lw_on) then
 
!      LW optical parameters
      do i=1,nlwbands            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol specific absorption in LW   

!     end points as basis for multidimentional linear interpolation  
      opt4d(1,1,1,1)=ka1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
      opt4d(1,1,1,2)=ka1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
      opt4d(1,1,2,1)=ka1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
      opt4d(1,1,2,2)=ka1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
      opt4d(1,2,1,1)=ka1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
      opt4d(1,2,1,2)=ka1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
      opt4d(1,2,2,1)=ka1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
      opt4d(1,2,2,2)=ka1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
      opt4d(2,1,1,1)=ka1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
      opt4d(2,1,1,2)=ka1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
      opt4d(2,1,2,1)=ka1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
      opt4d(2,1,2,2)=ka1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
      opt4d(2,2,1,1)=ka1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
      opt4d(2,2,1,2)=ka1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
      opt4d(2,2,2,1)=ka1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
      opt4d(2,2,2,2)=ka1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

!     interpolation in the fac, cat and fombg dimensions
      call lininterpol4dim (d2mx, dxm1, invd, opt4d, kabs1, kabs2)

      kabs1=max(kabs1,1.e-30)
      kabs2=max(kabs2,1.e-30)

!      write(*,*) 'Before kabs'
      if(t_xrh <= 0.37) then
        kabs(icol,k,kc10,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
        b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
        kabs(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i

      endif ! lw_on
         
          end do ! icol
        end do ! k

!       write(*,*) 'kcomp, omega(1,26,kcomp,4)=', kcomp, omega(1,26,kcomp,4)
!       write(*,*) 'kcomp, gass(1,26,kcomp,4)=', kcomp, gass(1,26,kcomp,4)
!       write(*,*) 'kcomp, bex(1,26,kcomp,4)=', kcomp, bex(1,26,kcomp,4)
!       write(*,*) 'kcomp, ske(1,26,kcomp,4)=', kcomp, ske(1,26,kcomp,4)

        end do  ! kcomp

      return
end subroutine interpol1


!********************************************************************************************

subroutine interpol2to3 (lchnk, ncol, daylight, xrh, irh1, irh2, &
                         mplus10, Nnatk, Camk, xfacsoain, &
                         omega, gass, bex, ske, lw_on, kabs, cxstot)


   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
!
! Input arguments
!
   integer, intent(in) :: lchnk                     ! chunk identifier
   integer, intent(in) :: ncol                      ! number of atmospheric columns
   integer, intent(in) :: mplus10                   ! mode number (0) or number + 10 (1)
   logical, intent(in) :: daylight(pcols)           ! only daylight calculations if .true.
   logical, intent(in) :: lw_on                       ! LW calculations are performed if true
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes) ! modal internally mixed SO4+BC+OC conc.
   real(r8), intent(in) :: xrh(pcols,pver)          ! level relative humidity (fraction)
   integer,  intent(in) :: irh1(pcols,pver)
   integer,  intent(in) :: irh2(pcols,pver)
   real(r8), intent(in) :: xfacsoain(pcols,pver,4)  ! mass fraction SOA/(H2SO4+SOA) for condensate
!
!
! Input-Output arguments
!
   real(r8), intent(inout) :: cxstot(pcols,pver)            ! excess internally mixed mass ("table overshoot") 
!
!
! Output arguments
!
   real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
   real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
   real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
   real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
   real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absorption coefficient
!
!---------------------------Local variables-----------------------------
!
      integer i, ierr, irelh, ictot, ifac, kcomp, k, icol, kc10
      integer ict1(pcols,pver), &
        ict2(pcols,pver), ifac1(pcols,pver), ifac2(pcols,pver)
      real(r8) a, b, camdiff
      real(r8) xct(pcols,pver), xfac(pcols,pver), &
        cxs(pcols,pver)
 
!     Temporary storage of often used array elements
      integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2
      real(r8) t_fac1, t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, &
        t_cat1, t_cat2
      real(r8) d2mx(3), dxm1(3), invd(3)
      real(r8) opt3d(2,2,2)
      real(r8) ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2  
      real(r8) kabs1, kabs2


      do k=1,pver
        do icol=1,ncol
          xct(icol,k)    = 0.0_r8
          xfac(icol,k)   = 0.0_r8
        end do
      end do

!      write(*,*) 'Before kcomp-loop'
!      do kcomp=2,3
      do kcomp=2,2

           if(mplus10==0) then
             kc10=kcomp
           else
             kc10=kcomp+10
           endif 

        do k=1,pver
         do icol=1,ncol

          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
                 /(Nnatk(icol,k,kc10)+eps),cate(kcomp,1)),cate(kcomp,16))
          xfac(icol,k) = min(max(xfacsoain(icol,k,kcomp),fac(1)),fac(6))
          camdiff=Camk(icol,k,kcomp)-xct(icol,k) &
                         *(Nnatk(icol,k,kc10)+eps)
          cxs(icol,k)=max(0.0_r8,camdiff)
          cxstot(icol,k)=cxstot(icol,k)+cxs(icol,k)

!        if(icol.eq.1) then
!          write(*,*) 'dir: k, kcomp, cxs =', k, kc10, cxs(icol,k)
!          write(*,*) 'k, kcomp, xct =', k, kc10, xct(icol,k)
!        endif

!         write(*,*) 'Before cate-loop', kc10
          do ictot=1,15
            if(xct(icol,k)>=cate(kcomp,ictot).and. &
            xct(icol,k) <= cate(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
          end do ! ictot

!         write(*,*) 'Before fac-loop', kcomp
          do ifac=1,5
            if(xfac(icol,k)>=fac(ifac).and. &
             xfac(icol,k) <= fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
          end do ! ifac

         end do ! icol
        end do ! k

!      write(*,*) 'Before init-loop', kc10
        do i=1,nbands
          do icol=1,ncol
            do k=1,pver
              omega(icol,k,kc10,i)=0.0_r8
              gass(icol,k,kc10,i)=0.0_r8
              bex(icol,k,kc10,i)=0.0_r8
              ske(icol,k,kc10,i)=0.0_r8
             end do
          end do
        end do
        do i=1,nlwbands
          do icol=1,ncol
            do k=1,pver
              kabs(icol,k,kc10,i)=0.0_r8
             end do
          end do
        end do
          
        do k=1,pver
          do icol=1,ncol

!      Collect all the vector elements into temporary storage
!      to avoid cache conflicts and excessive cross-referencing

      t_irh1 = irh1(icol,k)
      t_irh2 = irh2(icol,k)
      t_ict1 = ict1(icol,k)
      t_ict2 = ict2(icol,k)
      t_ifc1 = ifac1(icol,k)
      t_ifc2 = ifac2(icol,k)

!      write(*,*) 't_irh1,t_irh2=',t_irh1,t_irh2
!      write(*,*) 't_ict1,t_ict2=',t_ict1,t_ict2
!      write(*,*) 't_ifc1,t_ifc2=',t_ifc1,t_ifc2
!      write(*,*) 't_ifa1,t_ifa2=',t_ifa1,t_ifa2

      t_rh1  = rh(t_irh1)
      t_rh2  = rh(t_irh2)
      t_cat1 = cate(kcomp,t_ict1)
      t_cat2 = cate(kcomp,t_ict2)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)

!      write(*,*) 't_rh1,t_rh2,t_cat1,t_cat2=',t_rh1,t_rh2,t_cat1,t_cat2
!      write(*,*) 't_fac1,t_fac2=',t_fac1,t_fac2

      t_xrh  = xrh(icol,k)
      t_xct  = xct(icol,k)
      t_xfac = xfac(icol,k)

!     partial lengths along each dimension (1-4) for interpolation 
      d2mx(1) = (t_rh2-t_xrh)
      dxm1(1) = (t_xrh-t_rh1)
      invd(1) = 1.0_r8/(t_rh2-t_rh1)
      d2mx(2) = (t_cat2-t_xct)
      dxm1(2) = (t_xct-t_cat1)
      invd(2) = 1.0_r8/(t_cat2-t_cat1)
      d2mx(3) = (t_fac2-t_xfac)
      dxm1(3) = (t_xfac-t_fac1)
      invd(3) = 1.0_r8/(t_fac2-t_fac1)


!      SW optical parameters
       if(daylight(icol)) then

       do i=1,nbands            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  single scattering albedo:

!     end points as basis for multidimentional linear interpolation  
      opt3d(1,1,1)=om2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
      opt3d(1,1,2)=om2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
      opt3d(1,2,1)=om2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
      opt3d(1,2,2)=om2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
      opt3d(2,1,1)=om2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
      opt3d(2,1,2)=om2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
      opt3d(2,2,1)=om2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
      opt3d(2,2,2)=om2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

!     interpolation in the (fac and) cat dimension
      call lininterpol3dim (d2mx, dxm1, invd, opt3d, ome1, ome2)

!     finally, interpolation in the rh dimension
!      write(*,*) 'Before omega'
      omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) &
                          /(t_rh2-t_rh1)    
!      write(*,*) omega(icol,k,kc10,i)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  asymmetry factor   

!     end points as basis for multidimentional linear interpolation  
      opt3d(1,1,1)=g2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
      opt3d(1,1,2)=g2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
      opt3d(1,2,1)=g2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
      opt3d(1,2,2)=g2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
      opt3d(2,1,1)=g2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
      opt3d(2,1,2)=g2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
      opt3d(2,2,1)=g2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
      opt3d(2,2,2)=g2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

!     interpolation in the (fac and) cat dimension
      call lininterpol3dim (d2mx, dxm1, invd, opt3d, ge1, ge2)

!     finally, interpolation in the rh dimension
!      write(*,*) 'Before gass'
      gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                          /(t_rh2-t_rh1)    
!      write(*,*) gass(icol,k,kc10,i)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol extinction   

!     end points as basis for multidimentional linear interpolation  
      opt3d(1,1,1)=be2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
      opt3d(1,1,2)=be2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
      opt3d(1,2,1)=be2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
      opt3d(1,2,2)=be2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
      opt3d(2,1,1)=be2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
      opt3d(2,1,2)=be2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
      opt3d(2,2,1)=be2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
      opt3d(2,2,2)=be2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

!     interpolation in the (fac and) cat dimension
      call lininterpol3dim (d2mx, dxm1, invd, opt3d, bex1, bex2)

      bex1=max(bex1,1.e-30)
      bex2=max(bex2,1.e-30)

!     finally, interpolation in the rh dimension
!      write(*,*) 'Before bex'
      if(t_xrh <= 0.37) then
       bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
        b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
        bex(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i


      do i=4,4            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol specific extinction   

!     end points as basis for multidimentional linear interpolation  
      opt3d(1,1,1)=ke2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
      opt3d(1,1,2)=ke2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
      opt3d(1,2,1)=ke2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
      opt3d(1,2,2)=ke2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
      opt3d(2,1,1)=ke2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
      opt3d(2,1,2)=ke2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
      opt3d(2,2,1)=ke2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
      opt3d(2,2,2)=ke2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

!     interpolation in the (fac and) cat dimension
      call lininterpol3dim (d2mx, dxm1, invd, opt3d, ske1, ske2)

      ske1=max(ske1,1.e-30)
      ske2=max(ske2,1.e-30)

!     finally, interpolation in the rh dimension
!      write(*,*) 'Before ske'
      if(t_xrh <= 0.37) then
        ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
       b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
        ske(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i

       endif  ! daylight

      if (lw_on) then

!      LW optical parameters
      do i=1,nlwbands            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol specific absorption in LW   

!     end points as basis for multidimentional linear interpolation  
      opt3d(1,1,1)=ka2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
      opt3d(1,1,2)=ka2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
      opt3d(1,2,1)=ka2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
      opt3d(1,2,2)=ka2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
      opt3d(2,1,1)=ka2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
      opt3d(2,1,2)=ka2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
      opt3d(2,2,1)=ka2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
      opt3d(2,2,2)=ka2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

!     interpolation in the (fac and) cat dimension
      call lininterpol3dim (d2mx, dxm1, invd, opt3d, kabs1, kabs2)

      kabs1=max(kabs1,1.e-30_r8)
      kabs2=max(kabs2,1.e-30_r8)

!      write(*,*) 'Before kabs'
      if(t_xrh <= 0.37_r8) then
        kabs(icol,k,kc10,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
       b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
        kabs(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i

      endif ! lw_on

          end do ! icol
        end do ! k

!       write(*,*) 'kcomp, omega(1,26,kcomp,4)=', kcomp, omega(1,26,kcomp,4)
!       write(*,*) 'kcomp, gass(1,26,kcomp,4)=', kcomp, gass(1,26,kcomp,4)
!       write(*,*) 'kcomp, bex(1,26,kcomp,4)=', kcomp, bex(1,26,kcomp,4)
!       write(*,*) 'kcomp, ske(1,26,kcomp,4)=', kcomp, ske(1,26,kcomp,4)

      end do  ! kcomp

      return
end subroutine interpol2to3

!********************************************************************************************

subroutine interpol4 (lchnk, ncol, daylight, xrh, irh1, irh2, &
                      mplus10, Nnatk, xfbcbgin, &
                      Camk, xfacin, xfaqin, omega, gass, bex, ske, lw_on, kabs, cxstot)


   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none


!
! Input arguments
!
   integer, intent(in) :: lchnk                       ! chunk identifier
   integer, intent(in) :: ncol                        ! number of atmospheric columns
   integer, intent(in) :: mplus10                     ! mode number (0) or number + 10 (1)
   logical, intent(in) :: daylight(pcols)             ! only daylight calculations if .true.
   logical, intent(in) :: lw_on                       ! LW calculations are performed if true
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes)   ! modal internally mixed SO4+BC+OC conc.
   real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
   integer,  intent(in) :: irh1(pcols,pver)
   integer,  intent(in) :: irh2(pcols,pver)
   real(r8), intent(in) :: xfacin(pcols,pver,4)       ! modal (OC+BC)/(SO4+BC+OC)=SOA/(Sulfate+SOA)
   real(r8), intent(in) :: xfaqin(pcols,pver)         ! modal SO4(aq)/SO4
   real(r8), intent(in) :: xfbcbgin(pcols,pver)       ! mass fraction BC/(BC+OC) for the background mode 
!
!
! Input-Output arguments
!
   real(r8), intent(inout) :: cxstot(pcols,pver)    ! excess internally mixed mass ("table overshoot") 
!
!
! Output arguments
!
   real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
   real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
   real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
   real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
   real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absorption coefficient
!
!---------------------------Local variables-----------------------------
!
      integer i, ierr, irelh, ictot, ifac, ifbcbg, ifaq, kcomp, k, kc10, icol
      integer ict1(pcols,pver), &
        ict2(pcols,pver), ifac1(pcols,pver), ifac2(pcols,pver),     &
        ifbcbg1(pcols,pver), ifbcbg2(pcols,pver), ifaq1(pcols,pver),    &
        ifaq2(pcols,pver)
      real(r8) a, b, camdiff
      real(r8) xct(pcols,pver), xfac(pcols,pver), &
        xfbcbg(pcols,pver), xfaq(pcols,pver), cxs(pcols,pver)

!     Temporary storage of often used array elements
      integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifa1, t_ifa2,         &
       t_ifb1, t_ifb2, t_ifc1, t_ifc2
      real(r8) t_faq1, t_faq2, t_xfaq, t_fbcbg1, t_fbcbg2, t_xfbcbg, t_fac1, &
       t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, t_cat1, t_cat2

      real(r8) d2mx(5), dxm1(5), invd(5)
      real(r8) opt5d(2,2,2,2,2)
      real(r8) ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2  
      real(r8) kabs1, kabs2

      do k=1,pver
        do icol=1,ncol
          xct(icol,k)  = 0.0_r8
        end do
      end do

!      write(*,*) 'Before kcomp-loop'
        do kcomp=4,4

           if(mplus10==0) then
             kc10=kcomp
           else
             kc10=kcomp+10
           endif 

      do k=1,pver
        do icol=1,ncol

            xfac(icol,k) = 0.0_r8
            xfbcbg(icol,k) = 0.0_r8
            xfaq(icol,k) = 0.0_r8

!            if(Nnatk(icol,k,kcomp)>0.0_r8) then

          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
                 /(Nnatk(icol,k,kc10)+eps),cate(kcomp,1)),cate(kcomp,16))
          xfac(icol,k) = min(max(xfacin(icol,k,kcomp),fac(1)),fac(6))
          xfbcbg(icol,k) = min(max(xfbcbgin(icol,k),fbcbg(1)),fbcbg(6))
          xfaq(icol,k) = min(max(xfaqin(icol,k),faq(1)),faq(6))
          camdiff=Camk(icol,k,kcomp)-xct(icol,k) &
                         *(Nnatk(icol,k,kcomp)+eps)
          cxs(icol,k)=max(0.0_r8,camdiff)
          cxstot(icol,k)=cxstot(icol,k)+cxs(icol,k)

!        if(icol.eq.1) then
!          write(*,*) 'dir: k, kcomp, cxs =', k, kcomp, cxs(icol,k)
!        endif

!        write(*,*) 'Before cate-loop', kc10
         do ictot=1,15
            if(xct(icol,k)>=cate(kcomp,ictot).and. &
            xct(icol,k) <= cate(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
         end do ! ictot

!        write(*,*) 'Before fac-loop', kc10
         do ifac=1,5
            if(xfac(icol,k)>=fac(ifac).and. &
             xfac(icol,k) <= fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
         end do ! ifac

!        write(*,*) 'Before fbcbg-loop', kc10
         do ifbcbg=1,5
            if(xfbcbg(icol,k) >= fbcbg(ifbcbg).and. &
             xfbcbg(icol,k) <= fbcbg(ifbcbg+1)) then
             ifbcbg1(icol,k)=ifbcbg
             ifbcbg2(icol,k)=ifbcbg+1
            endif
         end do ! ifbcbg

!        write(*,*) 'Before faq-loop', kc10
         do ifaq=1,5
            if(xfaq(icol,k) >= faq(ifaq).and. &
            xfaq(icol,k) <= faq(ifaq+1)) then
             ifaq1(icol,k)=ifaq
             ifaq2(icol,k)=ifaq+1
            endif
         end do ! ifaq

        end do ! icol
      end do ! k

!      write(*,*) 'Before init-loop', kc10
        do i=1,nbands
          do icol=1,ncol
            do k=1,pver
              omega(icol,k,kc10,i)=0.0_r8
              gass(icol,k,kc10,i)=0.0_r8
              bex(icol,k,kc10,i)=0.0_r8
              ske(icol,k,kc10,i)=0.0_r8
            end do
          end do
        end do
        do i=1,nlwbands
          do icol=1,ncol
            do k=1,pver
              kabs(icol,k,kc10,i)=0.0_r8
             end do
          end do
        end do
         
        do k=1,pver
          do icol=1,ncol

!      Collect all the vector elements into temporary storage
!      to avoid cache conflicts and excessive cross-referencing

      t_irh1 = irh1(icol,k)
      t_irh2 = irh2(icol,k)
      t_ict1 = ict1(icol,k)
      t_ict2 = ict2(icol,k)
      t_ifc1 = ifac1(icol,k)
      t_ifc2 = ifac2(icol,k)
      t_ifb1 = ifbcbg1(icol,k)
      t_ifb2 = ifbcbg2(icol,k)
      t_ifa1 = ifaq1(icol,k)
      t_ifa2 = ifaq2(icol,k)

!      write(*,*) 't_irh1,t_irh2=',t_irh1,t_irh2
!      write(*,*) 't_ict1,t_ict2=',t_ict1,t_ict2
!      write(*,*) 't_ifc1,t_ifc2=',t_ifc1,t_ifc2
!      write(*,*) 't_ifb1,t_ifb2=',t_ifb1,t_ifb2
!      write(*,*) 't_ifa1,t_ifa2=',t_ifa1,t_ifa2

      t_rh1  = rh(t_irh1)
      t_rh2  = rh(t_irh2)
      t_cat1 = cate(kcomp,t_ict1)
      t_cat2 = cate(kcomp,t_ict2)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)
      t_fbcbg1 = fbcbg(t_ifb1)
      t_fbcbg2 = fbcbg(t_ifb2)
      t_faq1 = faq(t_ifa1)
      t_faq2 = faq(t_ifa2)

!      write(*,*) 't_rh1,t_rh2,t_cat1,t_cat2=',t_rh1,t_rh2,t_cat1,t_cat2
!      write(*,*) 't_fac1,t_fac2,t_fbcbg1,t_fbcbg2=',t_fac1,t_fac2,t_fbcbg1,t_fbcbg2
!      write(*,*) 't_faq1,t_faq2=',t_faq1,t_faq2

      t_xrh  = xrh(icol,k)
      t_xct  = xct(icol,k)
      t_xfac = xfac(icol,k)
      t_xfbcbg = xfbcbg(icol,k)
      t_xfaq = xfaq(icol,k)

!     partial lengths along each dimension (1-5) for interpolation 
      d2mx(1) = (t_rh2-t_xrh)
      dxm1(1) = (t_xrh-t_rh1)
      invd(1) = 1.0_r8/(t_rh2-t_rh1)
      d2mx(2) = (t_fbcbg2-t_xfbcbg)
      dxm1(2) = (t_xfbcbg-t_fbcbg1)
      invd(2) = 1.0_r8/(t_fbcbg2-t_fbcbg1)
      d2mx(3) = (t_cat2-t_xct)
      dxm1(3) = (t_xct-t_cat1)
      invd(3) = 1.0_r8/(t_cat2-t_cat1)
      d2mx(4) = (t_fac2-t_xfac)
      dxm1(4) = (t_xfac-t_fac1)
      invd(4) = 1.0_r8/(t_fac2-t_fac1)
      d2mx(5) = (t_faq2-t_xfaq)
      dxm1(5) = (t_xfaq-t_faq1)
      invd(5) = 1.0_r8/(t_faq2-t_faq1)

!      SW optical parameters
       if(daylight(icol)) then

      do i=1,nbands            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  single scattering albedo:

      opt5d(1,1,1,1,1)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
      opt5d(1,1,1,1,2)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
      opt5d(1,1,1,2,1)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
      opt5d(1,1,1,2,2)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
      opt5d(1,1,2,1,1)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
      opt5d(1,1,2,1,2)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
      opt5d(1,1,2,2,1)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
      opt5d(1,1,2,2,2)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
      opt5d(1,2,1,1,1)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
      opt5d(1,2,1,1,2)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
      opt5d(1,2,1,2,1)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
      opt5d(1,2,1,2,2)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
      opt5d(1,2,2,1,1)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
      opt5d(1,2,2,1,2)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
      opt5d(1,2,2,2,1)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
      opt5d(1,2,2,2,2)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
      opt5d(2,1,1,1,1)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
      opt5d(2,1,1,1,2)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
      opt5d(2,1,1,2,1)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
      opt5d(2,1,1,2,2)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
      opt5d(2,1,2,1,1)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
      opt5d(2,1,2,1,2)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
      opt5d(2,1,2,2,1)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
      opt5d(2,1,2,2,2)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
      opt5d(2,2,1,1,1)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
      opt5d(2,2,1,1,2)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
      opt5d(2,2,1,2,1)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
      opt5d(2,2,1,2,2)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
      opt5d(2,2,2,1,1)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
      opt5d(2,2,2,1,2)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
      opt5d(2,2,2,2,1)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
      opt5d(2,2,2,2,2)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

!     interpolation in the faq, fac, cat and fbcbg dimensions
      call lininterpol5dim (d2mx, dxm1, invd, opt5d, ome1, ome2)

!     finally, interpolation in the rh dimension 
!      write(*,*) 'Before omega'
      omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) &
                          /(t_rh2-t_rh1)    
!      write(*,*) omega(icol,k,kc10,i)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  asymmetry factor   

      opt5d(1,1,1,1,1)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
      opt5d(1,1,1,1,2)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
      opt5d(1,1,1,2,1)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
      opt5d(1,1,1,2,2)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
      opt5d(1,1,2,1,1)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
      opt5d(1,1,2,1,2)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
      opt5d(1,1,2,2,1)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
      opt5d(1,1,2,2,2)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
      opt5d(1,2,1,1,1)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
      opt5d(1,2,1,1,2)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
      opt5d(1,2,1,2,1)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
      opt5d(1,2,1,2,2)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
      opt5d(1,2,2,1,1)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
      opt5d(1,2,2,1,2)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
      opt5d(1,2,2,2,1)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
      opt5d(1,2,2,2,2)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
      opt5d(2,1,1,1,1)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
      opt5d(2,1,1,1,2)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
      opt5d(2,1,1,2,1)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
      opt5d(2,1,1,2,2)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
      opt5d(2,1,2,1,1)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
      opt5d(2,1,2,1,2)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
      opt5d(2,1,2,2,1)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
      opt5d(2,1,2,2,2)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
      opt5d(2,2,1,1,1)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
      opt5d(2,2,1,1,2)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
      opt5d(2,2,1,2,1)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
      opt5d(2,2,1,2,2)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
      opt5d(2,2,2,1,1)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
      opt5d(2,2,2,1,2)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
      opt5d(2,2,2,2,1)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
      opt5d(2,2,2,2,2)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

!     interpolation in the faq, fac, cat and fbcbg dimensions
      call lininterpol5dim (d2mx, dxm1, invd, opt5d, ge1, ge2)

!     finally, interpolation in the rh dimension 
!      write(*,*) 'Before gass'
      gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                   /(t_rh2-t_rh1)    

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol extinction   

      opt5d(1,1,1,1,1)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
      opt5d(1,1,1,1,2)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
      opt5d(1,1,1,2,1)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
      opt5d(1,1,1,2,2)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
      opt5d(1,1,2,1,1)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
      opt5d(1,1,2,1,2)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
      opt5d(1,1,2,2,1)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
      opt5d(1,1,2,2,2)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
      opt5d(1,2,1,1,1)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
      opt5d(1,2,1,1,2)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
      opt5d(1,2,1,2,1)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
      opt5d(1,2,1,2,2)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
      opt5d(1,2,2,1,1)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
      opt5d(1,2,2,1,2)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
      opt5d(1,2,2,2,1)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
      opt5d(1,2,2,2,2)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
      opt5d(2,1,1,1,1)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
      opt5d(2,1,1,1,2)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
      opt5d(2,1,1,2,1)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
      opt5d(2,1,1,2,2)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
      opt5d(2,1,2,1,1)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
      opt5d(2,1,2,1,2)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
      opt5d(2,1,2,2,1)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
      opt5d(2,1,2,2,2)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
      opt5d(2,2,1,1,1)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
      opt5d(2,2,1,1,2)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
      opt5d(2,2,1,2,1)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
      opt5d(2,2,1,2,2)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
      opt5d(2,2,2,1,1)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
      opt5d(2,2,2,1,2)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
      opt5d(2,2,2,2,1)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
      opt5d(2,2,2,2,2)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

!     interpolation in the faq, fac, cat and fbcbg dimensions
      call lininterpol5dim (d2mx, dxm1, invd, opt5d, bex1, bex2)

      bex1=max(bex1,1.e-30_r8)
      bex2=max(bex2,1.e-30_r8)

!     finally, interpolation in the rh dimension 
!      write(*,*) 'Before bex'
      if(t_xrh <= 0.37_r8) then
       bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
        b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
        bex(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i


      do i=4,4            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol specific extinction   

      opt5d(1,1,1,1,1)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
      opt5d(1,1,1,1,2)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
      opt5d(1,1,1,2,1)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
      opt5d(1,1,1,2,2)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
      opt5d(1,1,2,1,1)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
      opt5d(1,1,2,1,2)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
      opt5d(1,1,2,2,1)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
      opt5d(1,1,2,2,2)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
      opt5d(1,2,1,1,1)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
      opt5d(1,2,1,1,2)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
      opt5d(1,2,1,2,1)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
      opt5d(1,2,1,2,2)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
      opt5d(1,2,2,1,1)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
      opt5d(1,2,2,1,2)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
      opt5d(1,2,2,2,1)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
      opt5d(1,2,2,2,2)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
      opt5d(2,1,1,1,1)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
      opt5d(2,1,1,1,2)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
      opt5d(2,1,1,2,1)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
      opt5d(2,1,1,2,2)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
      opt5d(2,1,2,1,1)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
      opt5d(2,1,2,1,2)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
      opt5d(2,1,2,2,1)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
      opt5d(2,1,2,2,2)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
      opt5d(2,2,1,1,1)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
      opt5d(2,2,1,1,2)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
      opt5d(2,2,1,2,1)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
      opt5d(2,2,1,2,2)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
      opt5d(2,2,2,1,1)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
      opt5d(2,2,2,1,2)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
      opt5d(2,2,2,2,1)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
      opt5d(2,2,2,2,2)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

!     interpolation in the faq, fac, cat and fbcbg dimensions
      call lininterpol5dim (d2mx, dxm1, invd, opt5d, ske1, ske2)

      ske1=max(ske1,1.e-30_r8)
      ske2=max(ske2,1.e-30_r8)

!     finally, interpolation in the rh dimension 
!      write(*,*) 'Before ske'
      if(t_xrh <= 0.37_r8) then
        ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
       b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
        ske(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i

       endif  ! lw

      if (lw_on) then

!      LW optical parameters

      do i=1,nlwbands            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol specific absorption   

      opt5d(1,1,1,1,1)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
      opt5d(1,1,1,1,2)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
      opt5d(1,1,1,2,1)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
      opt5d(1,1,1,2,2)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
      opt5d(1,1,2,1,1)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
      opt5d(1,1,2,1,2)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
      opt5d(1,1,2,2,1)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
      opt5d(1,1,2,2,2)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
      opt5d(1,2,1,1,1)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
      opt5d(1,2,1,1,2)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
      opt5d(1,2,1,2,1)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
      opt5d(1,2,1,2,2)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
      opt5d(1,2,2,1,1)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
      opt5d(1,2,2,1,2)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
      opt5d(1,2,2,2,1)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
      opt5d(1,2,2,2,2)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
      opt5d(2,1,1,1,1)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
      opt5d(2,1,1,1,2)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
      opt5d(2,1,1,2,1)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
      opt5d(2,1,1,2,2)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
      opt5d(2,1,2,1,1)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
      opt5d(2,1,2,1,2)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
      opt5d(2,1,2,2,1)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
      opt5d(2,1,2,2,2)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
      opt5d(2,2,1,1,1)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
      opt5d(2,2,1,1,2)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
      opt5d(2,2,1,2,1)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
      opt5d(2,2,1,2,2)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
      opt5d(2,2,2,1,1)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
      opt5d(2,2,2,1,2)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
      opt5d(2,2,2,2,1)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
      opt5d(2,2,2,2,2)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

!     interpolation in the faq, fac, cat and fbcbg dimensions
      call lininterpol5dim (d2mx, dxm1, invd, opt5d, kabs1, kabs2)

      kabs1=max(kabs1,1.e-30_r8)
      kabs2=max(kabs2,1.e-30_r8)

!      write(*,*) 'Before kabs'
      if(t_xrh <= 0.37_r8) then
        kabs(icol,k,kc10,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
       b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
        kabs(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i

      endif ! lw_on
         
          end do ! icol
        end do ! k

!       write(*,*) 'kcomp, omega(1,26,kc10,4)=', kcomp, omega(1,26,kc10,4)
!       write(*,*) 'kcomp, gass(1,26,kc10,4)=', kcomp, gass(1,26,kc10,4)
!       write(*,*) 'kcomp, bex(1,26,kc10,4)=', kcomp, bex(1,26,kc10,4)
!       write(*,*) 'kcomp, ske(1,26,kc10,4)=', kcomp, ske(1,26,kc10,4)

      end do  ! kcomp

      return
end subroutine interpol4


!********************************************************************************************

subroutine interpol5to10 (lchnk, ncol, daylight, xrh, irh1, irh2, &
                          Nnatk, Camk, &
                          xfacin, xfbcin, xfaqin, omega, gass, bex, ske, lw_on, kabs, cxstot)


   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none


!
! Input arguments
!
   integer, intent(in) :: lchnk                       ! chunk identifier
   integer, intent(in) :: ncol                        ! number of atmospheric columns
   logical, intent(in) :: daylight(pcols)             ! only daylight calculations if .true.
   logical, intent(in) :: lw_on                       ! LW calculations are performed if true
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes)   ! modal internally mixed SO4+BC+OC conc.
   real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
   integer,  intent(in) :: irh1(pcols,pver)
   integer,  intent(in) :: irh2(pcols,pver)
   real(r8), intent(in) :: xfacin(pcols,pver,nbmodes) ! modal (OC+BC)/(SO4+BC+OC)
   real(r8), intent(in) :: xfbcin(pcols,pver,nbmodes) ! modal BC/(OC+BC)
   real(r8), intent(in) :: xfaqin(pcols,pver,nbmodes) ! modal SO4(aq)/SO4
!
!
! Input-Output arguments
!
   real(r8), intent(inout) :: cxstot(pcols,pver)    ! excess internally mixed mass ("table overshoot") 
!
!
! Output arguments
!
   real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
   real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
   real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
   real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
   real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absorption coefficient
!
!---------------------------Local variables-----------------------------
!
      integer i, ierr, irelh, ictot, ifac, ifbc, ifaq, kcomp, k, icol
      integer ict1(pcols,pver), &
        ict2(pcols,pver), ifac1(pcols,pver), ifac2(pcols,pver),     &
        ifbc1(pcols,pver), ifbc2(pcols,pver), ifaq1(pcols,pver),    &
        ifaq2(pcols,pver)
      real(r8) a, b, camdiff
      real(r8) xct(pcols,pver), xfac(pcols,pver,nbmodes), &
        xfbc(pcols,pver,nbmodes), xfaq(pcols,pver,nbmodes), cxs(pcols,pver)

!     Temporary storage of often used array elements
      integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifa1, t_ifa2,         &
       t_ifb1, t_ifb2, t_ifc1, t_ifc2
      real(r8) t_faq1, t_faq2, t_xfaq, t_fbc1, t_fbc2, t_xfbc, t_fac1, &
       t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, t_cat1, t_cat2
      real(r8) d2mx(5), dxm1(5), invd(5)
      real(r8) opt5d(2,2,2,2,2)
      real(r8) ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2  
      real(r8) kabs1, kabs2


      do k=1,pver
        do icol=1,ncol
          xct(icol,k)  = 0.0_r8
        end do
      end do

!      write(*,*) 'Before kcomp-loop'
        do kcomp=5,10

      do k=1,pver
        do icol=1,ncol

            xfac(icol,k,kcomp) = 0.0_r8
            xfbc(icol,k,kcomp) = 0.0_r8
            xfaq(icol,k,kcomp) = 0.0_r8

          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
                 /(Nnatk(icol,k,kcomp)+eps),cat(kcomp,1)),cat(kcomp,6))
          xfac(icol,k,kcomp) = min(max(xfacin(icol,k,kcomp),fac(1)),fac(6))
          xfbc(icol,k,kcomp) = min(max(xfbcin(icol,k,kcomp),fbc(1)),fbc(6))
          xfaq(icol,k,kcomp) = min(max(xfaqin(icol,k,kcomp),faq(1)),faq(6))
          camdiff=Camk(icol,k,kcomp)-xct(icol,k) &
                         *(Nnatk(icol,k,kcomp)+eps)
          cxs(icol,k)=max(0.0_r8,camdiff)
          cxstot(icol,k)=cxstot(icol,k)+cxs(icol,k)

!        if(icol.eq.1) then
!          write(*,*) 'dir: k, kcomp, cxs =', k, kcomp, cxs(icol,k)
!        endif

!        write(*,*) 'Before cat-loop', kcomp
         do ictot=1,5
            if(xct(icol,k)>=cat(kcomp,ictot).and. &
            xct(icol,k) <= cat(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
         end do ! ictot

!        write(*,*) 'Before fac-loop', kcomp
         do ifac=1,5
            if(xfac(icol,k,kcomp)>=fac(ifac).and. &
             xfac(icol,k,kcomp) <= fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
         end do ! ifac

!        write(*,*) 'Before fbc-loop', kcomp
         do ifbc=1,5
            if(xfbc(icol,k,kcomp) >= fbc(ifbc).and. &
             xfbc(icol,k,kcomp) <= fbc(ifbc+1)) then
             ifbc1(icol,k)=ifbc
             ifbc2(icol,k)=ifbc+1
            endif
         end do ! ifbc

!        write(*,*) 'Before faq-loop', kcomp
         do ifaq=1,5
            if(xfaq(icol,k,kcomp) >= faq(ifaq).and. &
            xfaq(icol,k,kcomp) <= faq(ifaq+1)) then
             ifaq1(icol,k)=ifaq
             ifaq2(icol,k)=ifaq+1
            endif
         end do ! ifaq

        end do ! icol
      end do ! k


!      write(*,*) 'Before init-loop', kcomp
        do i=1,nbands
          do icol=1,ncol
            do k=1,pver
              omega(icol,k,kcomp,i)=0.0_r8
              gass(icol,k,kcomp,i)=0.0_r8
              bex(icol,k,kcomp,i)=0.0_r8
              ske(icol,k,kcomp,i)=0.0_r8
            end do
          end do
        end do
        do i=1,nlwbands
          do icol=1,ncol
            do k=1,pver
              kabs(icol,k,kcomp,i)=0.0_r8
            end do
          end do
        end do
         
        do k=1,pver
          do icol=1,ncol

!      Collect all the vector elements into temporary storage
!      to avoid cache conflicts and excessive cross-referencing

      t_irh1 = irh1(icol,k)
      t_irh2 = irh2(icol,k)
      t_ict1 = ict1(icol,k)
      t_ict2 = ict2(icol,k)
      t_ifc1 = ifac1(icol,k)
      t_ifc2 = ifac2(icol,k)
      t_ifb1 = ifbc1(icol,k)
      t_ifb2 = ifbc2(icol,k)
      t_ifa1 = ifaq1(icol,k)
      t_ifa2 = ifaq2(icol,k)

!      write(*,*) 't_irh1,t_irh2=',t_irh1,t_irh2
!      write(*,*) 't_ict1,t_ict2=',t_ict1,t_ict2
!      write(*,*) 't_ifc1,t_ifc2=',t_ifc1,t_ifc2
!      write(*,*) 't_ifb1,t_ifb2=',t_ifb1,t_ifb2
!      write(*,*) 't_ifa1,t_ifa2=',t_ifa1,t_ifa2

      t_rh1  = rh(t_irh1)
      t_rh2  = rh(t_irh2)
      t_cat1 = cat(kcomp,t_ict1)
      t_cat2 = cat(kcomp,t_ict2)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)
      t_fbc1 = fbc(t_ifb1)
      t_fbc2 = fbc(t_ifb2)
      t_faq1 = faq(t_ifa1)
      t_faq2 = faq(t_ifa2)

!      write(*,*) 't_rh1,t_rh2,t_cat1,t_cat2=',t_rh1,t_rh2,t_cat1,t_cat2
!      write(*,*) 't_fac1,t_fac2,t_fbc1,t_fbc2=',t_fac1,t_fac2,t_fbc1,t_fbc2
!      write(*,*) 't_faq1,t_faq2=',t_faq1,t_faq2

      t_xrh  = xrh(icol,k)
      t_xct  = xct(icol,k)
      t_xfac = xfac(icol,k,kcomp)
      t_xfbc = xfbc(icol,k,kcomp)
      t_xfaq = xfaq(icol,k,kcomp)

!     partial lengths along each dimension (1-5) for interpolation 
      d2mx(1) = (t_rh2-t_xrh)
      dxm1(1) = (t_xrh-t_rh1)
      invd(1) = 1.0_r8/(t_rh2-t_rh1)
      d2mx(2) = (t_cat2-t_xct)
      dxm1(2) = (t_xct-t_cat1)
      invd(2) = 1.0_r8/(t_cat2-t_cat1)
      d2mx(3) = (t_fac2-t_xfac)
      dxm1(3) = (t_xfac-t_fac1)
      invd(3) = 1.0_r8/(t_fac2-t_fac1)
      d2mx(4) = (t_fbc2-t_xfbc)
      dxm1(4) = (t_xfbc-t_fbc1)
      invd(4) = 1.0_r8/(t_fbc2-t_fbc1)
      d2mx(5) = (t_faq2-t_xfaq)
      dxm1(5) = (t_xfaq-t_faq1)
      invd(5) = 1.0_r8/(t_faq2-t_faq1)


!      SW optical parameters
       if(daylight(icol)) then

      do i=1,nbands            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  single scattering albedo:

      opt5d(1,1,1,1,1)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,1,1,2)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,1,2,1)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,1,2,2)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,1,2,1,1)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,2,1,2)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,2,2,1)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,2,2,2)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,1,1,1)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,1,1,2)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,1,2,1)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,1,2,2)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,2,1,1)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,2,1,2)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,2,2,1)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,2,2,2)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,1,1,1)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,1,1,2)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,1,2,1)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,1,2,2)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,2,1,1)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,2,1,2)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,2,2,1)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,2,2,2)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,1,1,1)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,1,1,2)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,1,2,1)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,1,2,2)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,2,1,1)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,2,1,2)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,2,2,1)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,2,2,2)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

!     interpolation in the faq, fbc, fac and cat dimensions
      call lininterpol5dim (d2mx, dxm1, invd, opt5d, ome1, ome2)

!     finally, interpolation in the rh dimension 
!      write(*,*) 'Before omega'
      omega(icol,k,kcomp,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) &
                          /(t_rh2-t_rh1)    
!      write(*,*) omega(icol,k,kcomp,i)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  asymmetry factor   

      opt5d(1,1,1,1,1)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,1,1,2)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,1,2,1)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,1,2,2)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,1,2,1,1)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,2,1,2)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,2,2,1)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,2,2,2)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,1,1,1)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,1,1,2)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,1,2,1)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,1,2,2)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,2,1,1)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,2,1,2)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,2,2,1)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,2,2,2)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,1,1,1)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,1,1,2)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,1,2,1)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,1,2,2)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,2,1,1)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,2,1,2)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,2,2,1)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,2,2,2)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,1,1,1)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,1,1,2)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,1,2,1)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,1,2,2)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,2,1,1)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,2,1,2)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,2,2,1)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,2,2,2)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

!     interpolation in the faq, fbc, fac and cat dimensions
      call lininterpol5dim (d2mx, dxm1, invd, opt5d, ge1, ge2)

!     finally, interpolation in the rh dimension 
!      write(*,*) 'Before gass'
      gass(icol,k,kcomp,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                   /(t_rh2-t_rh1)    

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol extinction   

      opt5d(1,1,1,1,1)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,1,1,2)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,1,2,1)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,1,2,2)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,1,2,1,1)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,2,1,2)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,2,2,1)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,2,2,2)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,1,1,1)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,1,1,2)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,1,2,1)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,1,2,2)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,2,1,1)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,2,1,2)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,2,2,1)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,2,2,2)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,1,1,1)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,1,1,2)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,1,2,1)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,1,2,2)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,2,1,1)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,2,1,2)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,2,2,1)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,2,2,2)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,1,1,1)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,1,1,2)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,1,2,1)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,1,2,2)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,2,1,1)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,2,1,2)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,2,2,1)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,2,2,2)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

!     interpolation in the faq, fbc, fac and cat dimensions
      call lininterpol5dim (d2mx, dxm1, invd, opt5d, bex1, bex2)

      bex1=max(bex1,1.e-30_r8)
      bex2=max(bex2,1.e-30_r8)

!     finally, interpolation in the rh dimension 
!      write(*,*) 'Before bex'
      if(t_xrh <= 0.37_r8) then
       bex(icol,k,kcomp,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
        b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
        bex(icol,k,kcomp,i)=e**(a*t_xrh+b)
      endif

      end do ! i


      do i=4,4            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol specific extinction   

      opt5d(1,1,1,1,1)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,1,1,2)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,1,2,1)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,1,2,2)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,1,2,1,1)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,2,1,2)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,2,2,1)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,2,2,2)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,1,1,1)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,1,1,2)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,1,2,1)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,1,2,2)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,2,1,1)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,2,1,2)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,2,2,1)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,2,2,2)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,1,1,1)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,1,1,2)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,1,2,1)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,1,2,2)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,2,1,1)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,2,1,2)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,2,2,1)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,2,2,2)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,1,1,1)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,1,1,2)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,1,2,1)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,1,2,2)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,2,1,1)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,2,1,2)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,2,2,1)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,2,2,2)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

!     interpolation in the faq, fbc, fac and cat dimensions
      call lininterpol5dim (d2mx, dxm1, invd, opt5d, ske1, ske2)

      ske1=max(ske1,1.e-30_r8)
      ske2=max(ske2,1.e-30_r8)

!     finally, interpolation in the rh dimension 
!      write(*,*) 'Before ske'
      if(t_xrh <= 0.37_r8) then
        ske(icol,k,kcomp,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
       b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
        ske(icol,k,kcomp,i)=e**(a*t_xrh+b)
      endif

      end do ! i

       endif  ! daylight
      
      if (lw_on) then

!      LW optical parameters

      do i=1,nlwbands            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol specific absorption   

      opt5d(1,1,1,1,1)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,1,1,2)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,1,2,1)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,1,2,2)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,1,2,1,1)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,2,1,2)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,2,2,1)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,2,2,2)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,1,1,1)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,1,1,2)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,1,2,1)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,1,2,2)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,2,1,1)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,2,1,2)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,2,2,1)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,2,2,2)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,1,1,1)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,1,1,2)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,1,2,1)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,1,2,2)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,2,1,1)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,2,1,2)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,2,2,1)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,2,2,2)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,1,1,1)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,1,1,2)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,1,2,1)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,1,2,2)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,2,1,1)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,2,1,2)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,2,2,1)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,2,2,2)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

!     interpolation in the faq, fbc, fac and cat dimensions
      call lininterpol5dim (d2mx, dxm1, invd, opt5d, kabs1, kabs2)

      kabs1=max(kabs1,1.e-30_r8)
      kabs2=max(kabs2,1.e-30_r8)

!      write(*,*) 'Before kabs'
      if(t_xrh <= 0.37_r8) then
        kabs(icol,k,kcomp,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
       b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
        kabs(icol,k,kcomp,i)=e**(a*t_xrh+b)
      endif

      end do ! i

      endif ! lw_on
   
          end do ! icol
        end do ! k

!       write(*,*) 'kcomp, omega(1,26,kcomp,4)=', kcomp, omega(1,26,kcomp,4)
!       write(*,*) 'kcomp, gass(1,26,kcomp,4)=', kcomp, gass(1,26,kcomp,4)
!       write(*,*) 'kcomp, bex(1,26,kcomp,4)=', kcomp, bex(1,26,kcomp,4)
!       write(*,*) 'kcomp, ske(1,26,kcomp,4)=', kcomp, ske(1,26,kcomp,4)

      end do  ! kcomp

      return
end subroutine interpol5to10


!********************************************************************************************


end module optinterpol
