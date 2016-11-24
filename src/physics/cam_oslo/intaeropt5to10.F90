subroutine intaeropt5to10 (lchnk, ncol, xrh, irh1, irh2, Nnatk, Camk, xfacin, xfbcin, xfaqin, &
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1,             &
           beoc550lt1, beoc550gt1, besu550lt1, besu550gt1,             &
           backsc550, babg550, babc550, baoc550, basu550) 


   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8
   use opttab, only: cate, cat, fac, faq, fbc, rh 
   use commondefinitions, only: nmodes, nbmodes

   implicit none

#include <aerocopt.h>
!
! Input arguments
!
   integer, intent(in) :: lchnk                     ! chunk identifier
   integer, intent(in) :: ncol                      ! number of atmospheric columns
   real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
   integer,  intent(in) :: irh1(pcols,pver)
   integer,  intent(in) :: irh2(pcols,pver)
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes) ! modal internally mixed SO4+BC+OC conc.
   real(r8), intent(in) :: xfacin(pcols,pver,nbmodes) ! modal (OC+BC)/(SO4+BC+OC)
   real(r8), intent(in) :: xfbcin(pcols,pver,nbmodes) ! modal BC/(OC+BC)
   real(r8), intent(in) :: xfaqin(pcols,pver,nbmodes) ! modal SO4(aq)/SO4
!
! Output arguments: Modal total and absorption extiction coefficients (for AeroCom)
! for 550nm (1) and 865nm (2), and for r<1um (lt1) and r>1um (gt1).
! March 2009: + backscatter coefficient, backsc550 (km-1 sr-1).
! Rewritten by Alf Kirkevaag September 2015 to a more generalized for for 
! interpolations using common subroutines interpol*dim.
!
   real(r8), intent(out) :: &
     bext440(pcols,pver,0:nbmodes), babs440(pcols,pver,0:nbmodes), &
     bext500(pcols,pver,0:nbmodes), babs500(pcols,pver,0:nbmodes), &
     bext550(pcols,pver,0:nbmodes), babs550(pcols,pver,0:nbmodes), &
     bext670(pcols,pver,0:nbmodes), babs670(pcols,pver,0:nbmodes), &
     bext870(pcols,pver,0:nbmodes), babs870(pcols,pver,0:nbmodes), &
     bebg440(pcols,pver,0:nbmodes), & ! babg440(pcols,pver,0:nbmodes), &
     bebg500(pcols,pver,0:nbmodes), & ! babg500(pcols,pver,0:nbmodes), &
     bebg550(pcols,pver,0:nbmodes), babg550(pcols,pver,0:nbmodes), &
     bebg670(pcols,pver,0:nbmodes), & ! babg670(pcols,pver,0:nbmodes), &
     bebg870(pcols,pver,0:nbmodes), & ! babg870(pcols,pver,0:nbmodes), &
     bebc440(pcols,pver,0:nbmodes), & ! babc440(pcols,pver,0:nbmodes), &
     bebc500(pcols,pver,0:nbmodes), & ! babc500(pcols,pver,0:nbmodes), &
     bebc550(pcols,pver,0:nbmodes), babc550(pcols,pver,0:nbmodes), &
     bebc670(pcols,pver,0:nbmodes), & ! babc670(pcols,pver,0:nbmodes), &
     bebc870(pcols,pver,0:nbmodes), & ! babc870(pcols,pver,0:nbmodes), &
     beoc440(pcols,pver,0:nbmodes), & ! baoc440(pcols,pver,0:nbmodes), &
     beoc500(pcols,pver,0:nbmodes), & ! baoc500(pcols,pver,0:nbmodes), &
     beoc550(pcols,pver,0:nbmodes), baoc550(pcols,pver,0:nbmodes), &
     beoc670(pcols,pver,0:nbmodes), & ! baoc670(pcols,pver,0:nbmodes), &
     beoc870(pcols,pver,0:nbmodes), & ! baoc870(pcols,pver,0:nbmodes), &
     besu440(pcols,pver,0:nbmodes), & ! basu440(pcols,pver,0:nbmodes), &
     besu500(pcols,pver,0:nbmodes), & ! basu500(pcols,pver,0:nbmodes), &
     besu550(pcols,pver,0:nbmodes), basu550(pcols,pver,0:nbmodes), &
     besu670(pcols,pver,0:nbmodes), & ! basu670(pcols,pver,0:nbmodes), &
     besu870(pcols,pver,0:nbmodes), & ! basu870(pcols,pver,0:nbmodes), &
     bebg550lt1(pcols,pver,0:nbmodes), bebg550gt1(pcols,pver,0:nbmodes), &
     bebc550lt1(pcols,pver,0:nbmodes), bebc550gt1(pcols,pver,0:nbmodes), &
     beoc550lt1(pcols,pver,0:nbmodes), beoc550gt1(pcols,pver,0:nbmodes), &
     besu550lt1(pcols,pver,0:nbmodes), besu550gt1(pcols,pver,0:nbmodes), &
     backsc550(pcols,pver,0:nbmodes)
!
!---------------------------Local variables-----------------------------
!
      real(r8) a, b, e, eps
      real(r8) xct(pcols,pver), xfac(pcols,pver,nbmodes), &
        xfbc(pcols,pver,nbmodes), xfaq(pcols,pver,nbmodes)

      integer i, iv, ierr, irelh, ictot, ifac, ifbc, ifaq, kcomp, k, icol
      integer ict1(pcols,pver),&
       ict2(pcols,pver), ifac1(pcols,pver), ifac2(pcols,pver),     &
       ifbc1(pcols,pver), ifbc2(pcols,pver), ifaq1(pcols,pver),    &
       ifaq2(pcols,pver)

!      Temporary storage of often used array elements
      integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifa1, t_ifa2
      integer t_ifb1, t_ifb2, t_ifc1, t_ifc2
      real(r8)    t_faq1, t_faq2, t_xfaq
      real(r8)    t_fbc1, t_fbc2, t_xfbc
      real(r8)    t_fac1, t_fac2, t_xfac
      real(r8)    t_xrh, t_xct, t_rh1, t_rh2
      real(r8)    t_cat1, t_cat2
      real(r8) d2mx(5), dxm1(5), invd(5)
      real(r8) opt5d(2,2,2,2,2)
      real(r8) opt1, opt2, opt(38)

      parameter (e=2.718281828_r8, eps=1.0e-60_r8)


!      write(*,*) 'Before kcomp-loop'

!       Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.):

        do kcomp=5,10        

!        write(*,*) 'kcomp = ', kcomp 

!     initialize all output fields
      do k=1,pver
        do icol=1,ncol
         bext440(icol,k,kcomp)=0.0_r8 
         babs440(icol,k,kcomp)=0.0_r8 
         bext500(icol,k,kcomp)=0.0_r8 
         babs500(icol,k,kcomp)=0.0_r8 
         bext550(icol,k,kcomp)=0.0_r8 
         babs550(icol,k,kcomp)=0.0_r8 
         bext670(icol,k,kcomp)=0.0_r8 
         babs670(icol,k,kcomp)=0.0_r8 
         bext870(icol,k,kcomp)=0.0_r8 
         babs870(icol,k,kcomp)=0.0_r8 
         bebg440(icol,k,kcomp)=0.0_r8 
!         babg440(icol,k,kcomp)=0.0_r8 
         bebg500(icol,k,kcomp)=0.0_r8 
!         babg500(icol,k,kcomp)=0.0_r8 
         bebg550(icol,k,kcomp)=0.0_r8 
         babg550(icol,k,kcomp)=0.0_r8 
         bebg670(icol,k,kcomp)=0.0_r8 
!         babg670(icol,k,kcomp)=0.0_r8 
         bebg870(icol,k,kcomp)=0.0_r8 
!         babg870(icol,k,kcomp)=0.0_r8 
         bebc440(icol,k,kcomp)=0.0_r8 
!         babc440(icol,k,kcomp)=0.0_r8 
         bebc500(icol,k,kcomp)=0.0_r8 
!         babc500(icol,k,kcomp)=0.0_r8 
         bebc550(icol,k,kcomp)=0.0_r8 
         babc550(icol,k,kcomp)=0.0_r8 
         bebc670(icol,k,kcomp)=0.0_r8 
!         babc670(icol,k,kcomp)=0.0_r8 
         bebc870(icol,k,kcomp)=0.0_r8 
!         babc870(icol,k,kcomp)=0.0_r8 
         beoc440(icol,k,kcomp)=0.0_r8 
!         baoc440(icol,k,kcomp)=0.0_r8 
         beoc500(icol,k,kcomp)=0.0_r8 
!         baoc500(icol,k,kcomp)=0.0_r8 
         beoc550(icol,k,kcomp)=0.0_r8 
         baoc550(icol,k,kcomp)=0.0_r8 
         beoc670(icol,k,kcomp)=0.0_r8 
!         baoc670(icol,k,kcomp)=0.0_r8 
         beoc870(icol,k,kcomp)=0.0_r8 
!         baoc870(icol,k,kcomp)=0.0_r8 
         besu440(icol,k,kcomp)=0.0_r8 
!         basu440(icol,k,kcomp)=0.0_r8 
         besu500(icol,k,kcomp)=0.0_r8 
!         basu500(icol,k,kcomp)=0.0_r8 
         besu550(icol,k,kcomp)=0.0_r8 
         basu550(icol,k,kcomp)=0.0_r8 
         besu670(icol,k,kcomp)=0.0_r8 
!         basu670(icol,k,kcomp)=0.0_r8 
         besu870(icol,k,kcomp)=0.0_r8 
!         basu870(icol,k,kcomp)=0.0_r8 
         bebg550lt1(icol,k,kcomp)=0.0_r8 
         bebg550gt1(icol,k,kcomp)=0.0_r8 
         bebc550lt1(icol,k,kcomp)=0.0_r8 
         bebc550gt1(icol,k,kcomp)=0.0_r8 
         beoc550lt1(icol,k,kcomp)=0.0_r8 
         beoc550gt1(icol,k,kcomp)=0.0_r8 
         besu550lt1(icol,k,kcomp)=0.0_r8 
         besu550gt1(icol,k,kcomp)=0.0_r8 
         backsc550(icol,k,kcomp)=0.0_r8 
        end do
      end do

!      write(*,*) 'Before x-loop'
      do k=1,pver
         do icol=1,ncol

          if(Nnatk(icol,k,kcomp).gt.0) then

          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
                 /(Nnatk(icol,k,kcomp)+eps),cat(kcomp,1)),cat(kcomp,6))
          xfac(icol,k,kcomp) = min(max(xfacin(icol,k,kcomp),fac(1)),fac(6))
          xfbc(icol,k,kcomp) = min(max(xfbcin(icol,k,kcomp),fbc(1)),fbc(6))
          xfaq(icol,k,kcomp) = min(max(xfaqin(icol,k,kcomp),faq(1)),faq(6))

      do ictot=1,5
            if(xct(icol,k).ge.cat(kcomp,ictot).and. &
            xct(icol,k).le.cat(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
      end do ! ictot

      do ifac=1,5
            if(xfac(icol,k,kcomp).ge.fac(ifac).and. &
             xfac(icol,k,kcomp).le.fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
      end do ! ifac

      do ifbc=1,5
            if(xfbc(icol,k,kcomp).ge.fbc(ifbc).and. &
             xfbc(icol,k,kcomp).le.fbc(ifbc+1)) then
             ifbc1(icol,k)=ifbc
             ifbc2(icol,k)=ifbc+1
            endif
      end do ! ifbc

      do ifaq=1,5
            if(xfaq(icol,k,kcomp).ge.faq(ifaq).and. &
            xfaq(icol,k,kcomp).le.faq(ifaq+1)) then
             ifaq1(icol,k)=ifaq
             ifaq2(icol,k)=ifaq+1
            endif
      end do ! ifaq
           endif

          end do ! icol
        end do ! k


        do k=1,pver 
          do icol=1,ncol
         
           if(Nnatk(icol,k,kcomp).gt.0) then

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


         do iv=1,38  ! variable number

      opt5d(1,1,1,1,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,1,1,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,1,2,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,1,2,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,1,2,1,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,1,2,1,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,1,2,2,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,1,2,2,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,1,1,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,1,1,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,1,2,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,1,2,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(1,2,2,1,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(1,2,2,1,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(1,2,2,2,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(1,2,2,2,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,1,1,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,1,1,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,1,2,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,1,2,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,1,2,1,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,1,2,1,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,1,2,2,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,1,2,2,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,1,1,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,1,1,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,1,2,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,1,2,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      opt5d(2,2,2,1,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      opt5d(2,2,2,1,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      opt5d(2,2,2,2,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      opt5d(2,2,2,2,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

!     interpolation in the faq, fbc, fac and cat dimensions
      call lininterpol5dim (d2mx, dxm1, invd, opt5d, opt1, opt2)

!     finally, interpolation in the rh dimension 
!      write(*,*) 'Before opt'

      opt(iv)=((t_rh2-t_xrh)*opt1+(t_xrh-t_rh1)*opt2) &
                          /(t_rh2-t_rh1)    

!      write(*,*) opt(iv)

         end do ! iv=1,38 
 
         bext440(icol,k,kcomp)=opt(1)
         bext500(icol,k,kcomp)=opt(2)
         bext670(icol,k,kcomp)=opt(3)
         bext870(icol,k,kcomp)=opt(4)
         bebg440(icol,k,kcomp)=opt(5)
         bebg500(icol,k,kcomp)=opt(6)
         bebg670(icol,k,kcomp)=opt(7)
         bebg870(icol,k,kcomp)=opt(8)
         bebc440(icol,k,kcomp)=opt(9)
         bebc500(icol,k,kcomp)=opt(10)
         bebc670(icol,k,kcomp)=opt(11)
         bebc870(icol,k,kcomp)=opt(12)
         beoc440(icol,k,kcomp)=opt(13)
         beoc500(icol,k,kcomp)=opt(14)
         beoc670(icol,k,kcomp)=opt(15)
         beoc870(icol,k,kcomp)=opt(16)
         besu440(icol,k,kcomp)=opt(17)
         besu500(icol,k,kcomp)=opt(18)
         besu670(icol,k,kcomp)=opt(19)
         besu870(icol,k,kcomp)=opt(20)
         babs440(icol,k,kcomp)=opt(21)
         babs500(icol,k,kcomp)=opt(22)
         babs550(icol,k,kcomp)=opt(23)
         babs670(icol,k,kcomp)=opt(24)
         babs870(icol,k,kcomp)=opt(25)
         bebg550lt1(icol,k,kcomp)=opt(26)
         bebg550gt1(icol,k,kcomp)=opt(27)
         bebc550lt1(icol,k,kcomp)=opt(28)
         bebc550gt1(icol,k,kcomp)=opt(29)
         beoc550lt1(icol,k,kcomp)=opt(30)
         beoc550gt1(icol,k,kcomp)=opt(31)
         besu550lt1(icol,k,kcomp)=opt(32)
         besu550gt1(icol,k,kcomp)=opt(33)
         backsc550(icol,k,kcomp)=opt(34)
         babg550(icol,k,kcomp)=opt(35)
         babc550(icol,k,kcomp)=opt(36)
         baoc550(icol,k,kcomp)=opt(37)
         basu550(icol,k,kcomp)=opt(38)
         bebg550(icol,k,kcomp)=opt(26)+opt(27)
         bebc550(icol,k,kcomp)=opt(28)+opt(29)
         beoc550(icol,k,kcomp)=opt(30)+opt(31)
         besu550(icol,k,kcomp)=opt(32)+opt(33)
         bext550(icol,k,kcomp)=bebg550(icol,k,kcomp)+bebc550(icol,k,kcomp) &
                              +beoc550(icol,k,kcomp)+besu550(icol,k,kcomp)

           endif
         
!     if(beocgt1(icol,k,kcomp)>beoc1(icol,k,kcomp)) then
!       write(*,*) '5to10,kcomp,beocgt1,beoc1=', kcomp, beocgt1(icol,k,kcomp), beoc1(icol,k,kcomp)  
!     endif

       end do ! icol
      end do ! k

        end do  ! kcomp

      return

end subroutine intaeropt5to10




