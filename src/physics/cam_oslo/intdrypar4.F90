subroutine intdrypar4 (lchnk, ncol, Nnatk, xfbcbgnin, xfbcbgin, Camk, xfacin, xfaqin,  & 
           cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125, & 
           cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125, &
           cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol, &
           aaerosn,aaeroln,vaerosn,vaeroln,cknorm,cknlt05,ckngt125)

   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8
   use opttab, only: fbcbg, cate, cat, fac, faq, fbc, nbmp1
   use commondefinitions, only: nmodes, nbmodes

   implicit none

#include <aerodry.h>
!
! Input arguments
!
   integer, intent(in) :: lchnk                     ! chunk identifier
   integer, intent(in) :: ncol                      ! number of atmospheric columns
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes) ! modal internally mixed (cond+aq) SO4 conc.
   real(r8), intent(in) :: xfacin(pcols,pver,4)     ! modal (OC+BC)/(SO4+BC+OC)=SOA/(Sulfate+SOA)
   real(r8), intent(in) :: xfaqin(pcols,pver)       ! SO4(aq)/SO4
   real(r8), intent(in) :: xfbcbgin(pcols,pver)     ! mass fraction BC/(BC+OC) for the background mode (4)
   real(r8), intent(in) :: xfbcbgnin(pcols,pver)    ! mass fraction BC/(BC+OC) for the background mode (14)
!
! Input-Output arguments
!
   real(r8), intent(inout) :: &
     aaerosn(pcols,pver,nbmp1:nmodes), aaeroln(pcols,pver,nbmp1:nmodes), &
     vaerosn(pcols,pver,nbmp1:nmodes), vaeroln(pcols,pver,nbmp1:nmodes), &
     cknorm(pcols,pver,0:nmodes), cknlt05(pcols,pver,0:nmodes), ckngt125(pcols,pver,0:nmodes)
!
!
! Output arguments: Modal mass concentrations (cint), area (aaero) and volume (vaero)
! (for AeroCom determination of particle effective radii) of each constituent. cint*05 
! and cint*125 are  for r<0.5um and r>1.25um, respectively. aaeros and vaeros are
! integrated over r<0.5um, and aaerol and vaerol over r>0.5um.  
!
   real(r8), intent(out) :: &
     cintbg(pcols,pver,0:nbmodes), cintbg05(pcols,pver,0:nbmodes), cintbg125(pcols,pver,0:nbmodes), & 
     cintbc(pcols,pver,0:nbmodes), cintbc05(pcols,pver,0:nbmodes), cintbc125(pcols,pver,0:nbmodes), & 
     cintoc(pcols,pver,0:nbmodes), cintoc05(pcols,pver,0:nbmodes), cintoc125(pcols,pver,0:nbmodes), &
     cintsc(pcols,pver,0:nbmodes), cintsc05(pcols,pver,0:nbmodes), cintsc125(pcols,pver,0:nbmodes), &
     cintsa(pcols,pver,0:nbmodes), cintsa05(pcols,pver,0:nbmodes), cintsa125(pcols,pver,0:nbmodes), &
     aaeros(pcols,pver,0:nbmodes), aaerol(pcols,pver,0:nbmodes),                                    &
     vaeros(pcols,pver,0:nbmodes), vaerol(pcols,pver,0:nbmodes)
!
!---------------------------Local variables-----------------------------
!
      real(r8) a, b, e, eps
      real(r8) xct(pcols,pver), xfac(pcols,pver), xfaq(pcols,pver), &
        xfbcbgn(pcols,pver), xfbcbg(pcols,pver)

      integer iv, ierr, ifbcbg, ictot, ifac, ifaq, kcomp, k, icol
      integer ifbcbg1(pcols,pver), ifbcbg2(pcols,pver), &
       ict1(pcols,pver), ict2(pcols,pver), ifac1(pcols,pver), &
       ifac2(pcols,pver), ifaq1(pcols,pver), ifaq2(pcols,pver)
      integer ifbcbgn1(pcols,pver), ifbcbgn2(pcols,pver) 

!      Temporary storage of often used array elements
      integer t_ifb1, t_ifb2
      integer t_ict1, t_ict2, t_ifc1, t_ifc2, t_ifa1, t_ifa2
      real(r8) t_fbcbg1, t_fbcbg2  
      real(r8)    t_faq1, t_faq2, t_xfaq
      real(r8)    t_fac1, t_fac2, t_xfac
      real(r8)    t_xct,  t_cat1, t_cat2

      real(r8) t_xfbcbg
      real(r8) d2mx(4), dxm1(4), invd(4)
      real(r8) opt4d(2,2,2,2)
      real(r8) opt1, opt2, opt

      parameter (e=2.718281828_r8, eps=1.0e-60_r8)


!      write(*,*) 'Before kcomp-loop'

!       Mode 4, BC&OC(Ait):

          kcomp=4

!      initialize output fields
      do k=1,pver
         do icol=1,ncol
        cintbg(icol,k,kcomp)=0.0_r8
        cintbg05(icol,k,kcomp)=0.0_r8
        cintbg125(icol,k,kcomp)=0.0_r8
        cintbc(icol,k,kcomp)=0.0_r8
        cintbc05(icol,k,kcomp)=0.0_r8
        cintbc125(icol,k,kcomp)=0.0_r8
        cintoc(icol,k,kcomp)=0.0_r8
        cintoc05(icol,k,kcomp)=0.0_r8
        cintoc125(icol,k,kcomp)=0.0_r8
        cintsc(icol,k,kcomp)=0.0_r8
        cintsc05(icol,k,kcomp)=0.0_r8
        cintsc125(icol,k,kcomp)=0.0_r8
        cintsa(icol,k,kcomp)=0.0_r8
        cintsa05(icol,k,kcomp)=0.0_r8
        cintsa125(icol,k,kcomp)=0.0_r8
        aaeros(icol,k,kcomp)=0.0_r8
        aaerol(icol,k,kcomp)=0.0_r8
        vaeros(icol,k,kcomp)=0.0_r8
        vaerol(icol,k,kcomp)=0.0_r8
         end do
       end do

!      write(*,*) 'Before x-loop'
      do k=1,pver
         do icol=1,ncol

          if(Nnatk(icol,k,kcomp)>0.0_r8) then
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
          xfbcbg(icol,k) = min(max(xfbcbgin(icol,k),fbcbg(1)),fbcbg(6))
          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
                 /(Nnatk(icol,k,kcomp)+eps),cate(kcomp,1)),cate(kcomp,16))
          xfac(icol,k) = min(max(xfacin(icol,k,kcomp),fac(1)),fac(6))
          xfaq(icol,k) = min(max(xfaqin(icol,k),faq(1)),faq(6))

!   write(*,*) 'xct  =', icol, k, xct(icol,k)
!   write(*,*) 'xfac =', icol, k, xfac(icol,k)
!   write(*,*) 'xfaq =', icol, k, xfaq(icol,k)

!      write(*,*) 'Before fbcbg-loop', kcomp   (mode 4)
     do ifbcbg=1,5
            if(xfbcbg(icol,k) >= fbcbg(ifbcbg).and. &
             xfbcbg(icol,k) <= fbcbg(ifbcbg+1)) then
             ifbcbg1(icol,k)=ifbcbg
             ifbcbg2(icol,k)=ifbcbg+1
            endif
      end do ! ifbcbg

!      write(*,*) 'Before cat-loop', kcomp
      do ictot=1,15
            if(xct(icol,k)>=cate(kcomp,ictot).and. &
            xct(icol,k)<=cate(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
      end do ! ictot

!      write(*,*) 'Before fac-loop', kcomp
      do ifac=1,5
            if(xfac(icol,k)>=fac(ifac).and. &
             xfac(icol,k)<=fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
      end do ! ifac

!      write(*,*) 'Before faq-loop', kcomp
      do ifaq=1,5
            if(xfaq(icol,k)>=faq(ifaq).and. &
            xfaq(icol,k)<=faq(ifaq+1)) then
             ifaq1(icol,k)=ifaq
             ifaq2(icol,k)=ifaq+1
            endif
      end do ! ifaq

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
           endif

          end do ! icol
        end do ! k


        do k=1,pver 
          do icol=1,ncol
         
           if(Nnatk(icol,k,kcomp)>0.0_r8) then

!      Collect all the vector elements into temporary storage
!      to avoid cache conflicts and excessive cross-referencing
      t_ifb1 = ifbcbg1(icol,k)
      t_ifb2 = ifbcbg2(icol,k)
      t_ict1 = ict1(icol,k)
      t_ict2 = ict2(icol,k)
      t_ifc1 = ifac1(icol,k)
      t_ifc2 = ifac2(icol,k)
      t_ifa1 = ifaq1(icol,k)
      t_ifa2 = ifaq2(icol,k)
      t_fbcbg1 = fbcbg(t_ifb1)
      t_fbcbg2 = fbcbg(t_ifb2)
      t_cat1 = cate(kcomp,t_ict1)
      t_cat2 = cate(kcomp,t_ict2)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)
      t_faq1 = faq(t_ifa1)
      t_faq2 = faq(t_ifa2)
      t_xfbcbg = xfbcbg(icol,k)
      t_xct  = xct(icol,k)
      t_xfac = xfac(icol,k)
      t_xfaq = xfaq(icol,k)

!     partial lengths along each dimension (1-5) for interpolation 
      d2mx(1) = (t_fbcbg2-t_xfbcbg)
      dxm1(1) = (t_xfbcbg-t_fbcbg1)
      invd(1) = 1.0_r8/(t_fbcbg2-t_fbcbg1)
      d2mx(2) = (t_cat2-t_xct)
      dxm1(2) = (t_xct-t_cat1)
      invd(2) = 1.0_r8/(t_cat2-t_cat1)
      d2mx(3) = (t_fac2-t_xfac)
      dxm1(3) = (t_xfac-t_fac1)
      invd(3) = 1.0_r8/(t_fac2-t_fac1)
      d2mx(4) = (t_faq2-t_xfaq)
      dxm1(4) = (t_xfaq-t_faq1)
      invd(4) = 1.0_r8/(t_faq2-t_faq1)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

         do iv=1,19  ! variable number

!     end points as basis for multidimentional linear interpolation  
      opt4d(1,1,1,1)=a4var(iv,t_ifb1,t_ict1,t_ifc1,t_ifa1)
      opt4d(1,1,1,2)=a4var(iv,t_ifb1,t_ict1,t_ifc1,t_ifa2)
      opt4d(1,1,2,1)=a4var(iv,t_ifb1,t_ict1,t_ifc2,t_ifa1)
      opt4d(1,1,2,2)=a4var(iv,t_ifb1,t_ict1,t_ifc2,t_ifa2)
      opt4d(1,2,1,1)=a4var(iv,t_ifb1,t_ict2,t_ifc1,t_ifa1)
      opt4d(1,2,1,2)=a4var(iv,t_ifb1,t_ict2,t_ifc1,t_ifa2)
      opt4d(1,2,2,1)=a4var(iv,t_ifb1,t_ict2,t_ifc2,t_ifa1)
      opt4d(1,2,2,2)=a4var(iv,t_ifb1,t_ict2,t_ifc2,t_ifa2)
      opt4d(2,1,1,1)=a4var(iv,t_ifb2,t_ict1,t_ifc1,t_ifa1)
      opt4d(2,1,1,2)=a4var(iv,t_ifb2,t_ict1,t_ifc1,t_ifa2)
      opt4d(2,1,2,1)=a4var(iv,t_ifb2,t_ict1,t_ifc2,t_ifa1)
      opt4d(2,1,2,2)=a4var(iv,t_ifb2,t_ict1,t_ifc2,t_ifa2)
      opt4d(2,2,1,1)=a4var(iv,t_ifb2,t_ict2,t_ifc1,t_ifa1)
      opt4d(2,2,1,2)=a4var(iv,t_ifb2,t_ict2,t_ifc1,t_ifa2)
      opt4d(2,2,2,1)=a4var(iv,t_ifb2,t_ict2,t_ifc2,t_ifa1)
      opt4d(2,2,2,2)=a4var(iv,t_ifb2,t_ict2,t_ifc2,t_ifa2)

!     interpolation in the faq, fac and cat dimensions
      call lininterpol4dim (d2mx, dxm1, invd, opt4d, opt1, opt2)

!     finally, interpolation in the fbcbg dimension 
      opt = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

!      if(k.eq.1) write(*,*) 'opt4 =', opt

!      write(*,*) 'Before array'

       if(iv==1) then
         cintbg(icol,k,kcomp)=opt
       elseif(iv==2) then
         cintbg05(icol,k,kcomp)=opt
       elseif(iv==3) then
         cintbg125(icol,k,kcomp)=opt
       elseif(iv==4) then
        cintbc(icol,k,kcomp)=opt
       elseif(iv==5) then
        cintbc05(icol,k,kcomp)=opt
       elseif(iv==6) then
        cintbc125(icol,k,kcomp)=opt
       elseif(iv==7) then
        cintoc(icol,k,kcomp)=opt
       elseif(iv==8) then
        cintoc05(icol,k,kcomp)=opt
       elseif(iv==9) then
        cintoc125(icol,k,kcomp)=opt
       elseif(iv==10) then
        cintsc(icol,k,kcomp)=opt
       elseif(iv==11) then
        cintsc05(icol,k,kcomp)=opt
       elseif(iv==12) then
        cintsc125(icol,k,kcomp)=opt
       elseif(iv==13) then
        cintsa(icol,k,kcomp)=opt
       elseif(iv==14) then
        cintsa05(icol,k,kcomp)=opt
       elseif(iv==15) then
        cintsa125(icol,k,kcomp)=opt
       elseif(iv==16) then
        aaeros(icol,k,kcomp)=opt
       elseif(iv==17) then
        aaerol(icol,k,kcomp)=opt
       elseif(iv==18) then
        vaeros(icol,k,kcomp)=opt
       elseif(iv==19) then
        vaerol(icol,k,kcomp)=opt
       endif

         end do ! iv=1,19 

           endif
         
       end do ! icol
      end do ! k

      kcomp=14
      do k=1,pver 
        do icol=1,ncol

          xfbcbgn(icol,k) = min(max(xfbcbgnin(icol,k),fbcbg(1)),fbcbg(6))
!         write(*,*) 'Before fbcbgn-loop', kcomp   (mode 14)
          do ifbcbg=1,5
            if(xfbcbgn(icol,k) >= fbcbg(ifbcbg).and. &
             xfbcbgn(icol,k) <= fbcbg(ifbcbg+1)) then
             ifbcbgn1(icol,k)=ifbcbg
             ifbcbgn2(icol,k)=ifbcbg+1
            endif
          end do ! ifbcbg

         t_ifb1 = ifbcbgn1(icol,k)
         t_ifb2 = ifbcbgn2(icol,k)
         t_fbcbg1 = fbcbg(t_ifb1)
         t_fbcbg2 = fbcbg(t_ifb2)
         t_xfbcbg = xfbcbgn(icol,k)

         d2mx(1) = (t_fbcbg2-t_xfbcbg)
         dxm1(1) = (t_xfbcbg-t_fbcbg1)
         invd(1) = 1.0_r8/(t_fbcbg2-t_fbcbg1)

!     Only interpolation in the fbcbg dimension for mode 14
         opt1 = a4var(1,1,1,1,1)
         opt2 = a4var(1,2,1,1,1)
         cknorm(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)
         opt1 = a4var(2,1,1,1,1)
         opt2 = a4var(2,2,1,1,1)
         cknlt05(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)
         opt1 = a4var(3,1,1,1,1)
         opt2 = a4var(3,2,1,1,1)
         ckngt125(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)
         opt1 = a4var(16,1,1,1,1)
         opt2 = a4var(16,2,1,1,1)
!        (The remaining variables are actually independent of fbcbg,
!        but we follow the same procedure anyway:) 
         aaerosn(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)
         opt1 = a4var(17,1,1,1,1)
         opt2 = a4var(17,2,1,1,1)
         aaeroln(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)
         opt1 = a4var(18,1,1,1,1)
         opt2 = a4var(18,2,1,1,1)
         vaerosn(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)
         opt1 = a4var(19,1,1,1,1)
         opt2 = a4var(19,2,1,1,1)
         vaeroln(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

       end do ! icol
      end do ! k

      return
end subroutine intdrypar4




