
 subroutine opticsAtConstRh (lchnk, ncol, pint, rhoda, Nnatk, xrh, irh1, irh2, &
                             Cam, camnull, fcm, fbcm, faqm, &
                             faqm4, fnbc, faitbc, focm, f_soana,  &
                             vnbc, vaitbc, v_soana)

!     Extra AeroCom diagnostics requiring table look-ups with RH=0%.

   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_history,  only: outfld
   use constituents, only: pcnst
   use opttab
   use const
   use aerosoldef
   use commondefinitions
!soa   use optinterpol,     only: interpol0,interpol1,interpol2to3,interpol4,interpol5to10
   use physics_types,   only: physics_state

   implicit none

!
! Input arguments
!
   integer, intent(in) :: lchnk                   ! chunk identifier
   integer, intent(in) :: ncol                    ! number of atmospheric columns
   real(r8), intent(in) :: pint(pcols,pverp)      ! Model interface pressures (10*Pa)
   real(r8), intent(in) :: rhoda(pcols,pver)      ! Density of dry air (kg/m^3)
   real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
   integer,  intent(in) :: irh1(pcols,pver)
   integer,  intent(in) :: irh2(pcols,pver)
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! aerosol mode number concentration  
   real(r8), intent(in) :: fnbc(pcols,pver) 
   real(r8), intent(in) :: faitbc(pcols,pver)
   real(r8), intent(in) :: Cam(pcols,pver,nbmodes)
   real(r8), intent(in) :: fbcm(pcols,pver,nbmodes)
   real(r8), intent(in) :: fcm(pcols,pver,nbmodes)
   real(r8), intent(in) :: faqm(pcols,pver,nbmodes)
   real(r8), intent(in) :: camnull(pcols,pver,nbmodes)
   real(r8), intent(in) :: faqm4(pcols,pver) 
   real(r8), intent(in) :: f_soana(pcols,pver)
   real(r8), intent(in) :: focm(pcols,pver,4)     ! = fraction of added mass which is either SOA cond. or OC coag.
   real(r8), intent(in) :: vnbc(pcols,pver)
   real(r8), intent(in) :: vaitbc(pcols,pver)
   real(r8), intent(in) :: v_soana(pcols,pver)
!
!
!---------------------------Local variables-----------------------------
!
   integer  i, k, icol, mplus10
   integer  iloop

   real(r8) deltah  ! soa: vnbc, vaitbc
   real(r8) dod550dry(pcols), abs550dry(pcols)
!
   real(r8) bext440(pcols,pver,0:nbmodes), babs440(pcols,pver,0:nbmodes), &
            bext500(pcols,pver,0:nbmodes), babs500(pcols,pver,0:nbmodes), &
            bext550(pcols,pver,0:nbmodes), babs550(pcols,pver,0:nbmodes), &
            bext670(pcols,pver,0:nbmodes), babs670(pcols,pver,0:nbmodes), &
            bext870(pcols,pver,0:nbmodes), babs870(pcols,pver,0:nbmodes), &
            bebg440(pcols,pver,0:nbmodes), babg440(pcols,pver,0:nbmodes), &
            bebg500(pcols,pver,0:nbmodes), babg500(pcols,pver,0:nbmodes), &
            bebg550(pcols,pver,0:nbmodes), babg550(pcols,pver,0:nbmodes), &
            bebg670(pcols,pver,0:nbmodes), babg670(pcols,pver,0:nbmodes), &
            bebg870(pcols,pver,0:nbmodes), babg870(pcols,pver,0:nbmodes), &
            bebc440(pcols,pver,0:nbmodes), babc440(pcols,pver,0:nbmodes), &
            bebc500(pcols,pver,0:nbmodes), babc500(pcols,pver,0:nbmodes), &
            bebc550(pcols,pver,0:nbmodes), babc550(pcols,pver,0:nbmodes), &
            bebc670(pcols,pver,0:nbmodes), babc670(pcols,pver,0:nbmodes), &
            bebc870(pcols,pver,0:nbmodes), babc870(pcols,pver,0:nbmodes), &
            beoc440(pcols,pver,0:nbmodes), baoc440(pcols,pver,0:nbmodes), &
            beoc500(pcols,pver,0:nbmodes), baoc500(pcols,pver,0:nbmodes), &
            beoc550(pcols,pver,0:nbmodes), baoc550(pcols,pver,0:nbmodes), &
            beoc670(pcols,pver,0:nbmodes), baoc670(pcols,pver,0:nbmodes), &
            beoc870(pcols,pver,0:nbmodes), baoc870(pcols,pver,0:nbmodes), &
            besu440(pcols,pver,0:nbmodes), basu440(pcols,pver,0:nbmodes), &
            besu500(pcols,pver,0:nbmodes), basu500(pcols,pver,0:nbmodes), &
            besu550(pcols,pver,0:nbmodes), basu550(pcols,pver,0:nbmodes), &
            besu670(pcols,pver,0:nbmodes), basu670(pcols,pver,0:nbmodes), &
            besu870(pcols,pver,0:nbmodes), basu870(pcols,pver,0:nbmodes)
   real(r8) bebglt1(pcols,pver,0:nbmodes), bebggt1(pcols,pver,0:nbmodes), &
            bebclt1(pcols,pver,0:nbmodes), bebcgt1(pcols,pver,0:nbmodes), & 
            beoclt1(pcols,pver,0:nbmodes), beocgt1(pcols,pver,0:nbmodes), & 
            bes4lt1(pcols,pver,0:nbmodes), bes4gt1(pcols,pver,0:nbmodes), &
            backsc550(pcols,pver,0:nbmodes), backsc550x(pcols,pver,nbmp1:nmodes), &
            ec550dry_aer(pcols,pver), abs550dry_aer(pcols,pver)
   real(r8) bebglt1t(pcols,pver), bebclt1t(pcols,pver), &
            beoclt1t(pcols,pver), bes4lt1t(pcols,pver)
   real(r8) basu550tot(pcols,pver), babc550tot(pcols,pver), baoc550tot(pcols,pver), &
            babc550xt(pcols,pver), baoc550xt(pcols,pver), &
            ba550x(pcols,pver,nbmp1:nmodes), belt1x(pcols,pver,nbmp1:nmodes)
!           Additional AeroCom Phase III output:   
   real(r8) ec440dry_aer(pcols,pver), abs440dry_aer(pcols,pver), &
            ec870dry_aer(pcols,pver), abs870dry_aer(pcols,pver), &
            be550lt1_aer(pcols,pver,0:nbmodes), ec550drylt1_aer(pcols,pver), &
            abs550dry_bc(pcols,pver), abs550dry_oc(pcols,pver), &
            abs550dry_su(pcols,pver), abs550dry_ss(pcols,pver), &
            abs550dry_du(pcols,pver), ec550drylt1_bc(pcols,pver), &
            ec550drylt1_oc(pcols,pver), ec550drylt1_su(pcols,pver), &
            ec550drylt1_ss(pcols,pver), ec550drylt1_du(pcols,pver) 
!  
   real(r8) bext440n(pcols,pver,0:nbmodes), babs440n(pcols,pver,0:nbmodes), &
            bext500n(pcols,pver,0:nbmodes), babs500n(pcols,pver,0:nbmodes), &
            bext550n(pcols,pver,0:nbmodes), babs550n(pcols,pver,0:nbmodes), &
            bext670n(pcols,pver,0:nbmodes), babs670n(pcols,pver,0:nbmodes), &
            bext870n(pcols,pver,0:nbmodes), babs870n(pcols,pver,0:nbmodes), &
            bebg440n(pcols,pver,0:nbmodes), babg440n(pcols,pver,0:nbmodes), &
            bebg500n(pcols,pver,0:nbmodes), babg500n(pcols,pver,0:nbmodes), &
            bebg550n(pcols,pver,0:nbmodes), babg550n(pcols,pver,0:nbmodes), &
            bebg670n(pcols,pver,0:nbmodes), babg670n(pcols,pver,0:nbmodes), &
            bebg870n(pcols,pver,0:nbmodes), babg870n(pcols,pver,0:nbmodes), &
            bebc440n(pcols,pver,0:nbmodes), babc440n(pcols,pver,0:nbmodes), &
            bebc500n(pcols,pver,0:nbmodes), babc500n(pcols,pver,0:nbmodes), &
            bebc550n(pcols,pver,0:nbmodes), babc550n(pcols,pver,0:nbmodes), &
            bebc670n(pcols,pver,0:nbmodes), babc670n(pcols,pver,0:nbmodes), &
            bebc870n(pcols,pver,0:nbmodes), babc870n(pcols,pver,0:nbmodes), &
            beoc440n(pcols,pver,0:nbmodes), baoc440n(pcols,pver,0:nbmodes), &
            beoc500n(pcols,pver,0:nbmodes), baoc500n(pcols,pver,0:nbmodes), &
            beoc550n(pcols,pver,0:nbmodes), baoc550n(pcols,pver,0:nbmodes), &
            beoc670n(pcols,pver,0:nbmodes), baoc670n(pcols,pver,0:nbmodes), &
            beoc870n(pcols,pver,0:nbmodes), baoc870n(pcols,pver,0:nbmodes), &
            besu440n(pcols,pver,0:nbmodes), basu440n(pcols,pver,0:nbmodes), &
            besu500n(pcols,pver,0:nbmodes), basu500n(pcols,pver,0:nbmodes), &
            besu550n(pcols,pver,0:nbmodes), basu550n(pcols,pver,0:nbmodes), &
            besu670n(pcols,pver,0:nbmodes), basu670n(pcols,pver,0:nbmodes), &
            besu870n(pcols,pver,0:nbmodes), basu870n(pcols,pver,0:nbmodes)
   real(r8) bebglt1n(pcols,pver,0:nbmodes), bebggt1n(pcols,pver,0:nbmodes), &
            bebclt1n(pcols,pver,0:nbmodes), bebcgt1n(pcols,pver,0:nbmodes), & 
            beoclt1n(pcols,pver,0:nbmodes), beocgt1n(pcols,pver,0:nbmodes), & 
            bes4lt1n(pcols,pver,0:nbmodes), bes4gt1n(pcols,pver,0:nbmodes), &
            backsc550n(pcols,pver,0:nbmodes)
   real(r8) bedustlt1(pcols,pver), bedustgt1(pcols,pver), &
            besslt1(pcols,pver), bessgt1(pcols,pver)
   real(r8) bbclt1xt(pcols,pver), &
            boclt1xt(pcols,pver), bocgt1xt(pcols,pver)


!000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

      do k=1,pver
        do icol=1,ncol
          ec550dry_aer(icol,k)=0.0_r8 
          abs550dry_aer(icol,k)=0.0_r8 
          ec550drylt1_aer(icol,k)=0.0_r8 
          abs550dry_bc(icol,k)=0.0_r8
          abs550dry_oc(icol,k)=0.0_r8
          abs550dry_su(icol,k)=0.0_r8
          abs550dry_ss(icol,k)=0.0_r8
          abs550dry_du(icol,k)=0.0_r8
          ec440dry_aer(icol,k)=0.0_r8 
          abs440dry_aer(icol,k)=0.0_r8 
          ec870dry_aer(icol,k)=0.0_r8 
          abs870dry_aer(icol,k)=0.0_r8 
          basu550tot(icol,k)=0.0_r8 
          babc550tot(icol,k)=0.0_r8 
          baoc550tot(icol,k)=0.0_r8 
          bebglt1t(icol,k)=0.0_r8
          bebclt1t(icol,k)=0.0_r8
          beoclt1t(icol,k)=0.0_r8
          bes4lt1t(icol,k)=0.0_r8
          bedustlt1(icol,k)=0.0_r8
          besslt1(icol,k)=0.0_r8
        end do
      end do
      do icol=1,ncol
        dod550dry(icol)=0.0_r8 
        abs550dry(icol)=0.0_r8 
      end do

      do iloop=1,1

!     BC(ax) mode (hydrophobic, so no rhum needed here):
        call intaeropt0(lchnk, ncol, Nnatk, &
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         &
           backsc550, babg550, babc550, baoc550, basu550)


!     SO4(Ait), BC(Ait) and OC(Ait) modes:
      mplus10=0
        call intaeropt2to3(lchnk, ncol, xrh, irh1, irh2, mplus10, Nnatk, Cam, focm,&
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         &
           backsc550, babg550, babc550, baoc550, basu550)

      mplus10=0
        call intaeropt1(lchnk, ncol, xrh, irh1, irh2, mplus10,         &
           Nnatk, f_soana, Cam, focm,                                  &
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         &
           backsc550, babg550, babc550, baoc550, basu550)

!     BC&OC(Ait) mode:    ------ fcm invalid here (=0). Using faitbc instead
      mplus10=0
        call intaeropt4(lchnk, ncol, xrh, irh1, irh2, mplus10, Nnatk, faitbc, Cam, focm, faqm4, &
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         &
           backsc550, babg550, babc550, baoc550, basu550)
  
!     SO4(Ait75) (5), Mineral (6-7) and Sea-salt (8-10) modes:
        call intaeropt5to10(lchnk, ncol, xrh, irh1, irh2, Nnatk, Cam, fcm, fbcm, faqm, &
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebglt1, bebggt1, bebclt1, bebcgt1,                         &
           beoclt1, beocgt1, bes4lt1, bes4gt1,                         &
           backsc550, babg550, babc550, baoc550, basu550)

!     then to the externally mixed SO4(n), BC(n) and OC(n) modes:
      mplus10=1
        call intaeropt2to3(lchnk, ncol, xrh, irh1, irh2, mplus10, Nnatk, camnull, focm,&
           bext440n, bext500n, bext550n, bext670n, bext870n,           &
           bebg440n, bebg500n, bebg550n, bebg670n, bebg870n,           &
           bebc440n, bebc500n, bebc550n, bebc670n, bebc870n,           &
           beoc440n, beoc500n, beoc550n, beoc670n, beoc870n,           &
           besu440n, besu500n, besu550n, besu670n, besu870n,           &
           babs440n, babs500n, babs550n, babs670n, babs870n,           &
           bebglt1n, bebggt1n, bebclt1n, bebcgt1n,                     &
           beoclt1n, beocgt1n, bes4lt1n, bes4gt1n,                     &
           backsc550n, babg550n, babc550n, baoc550n, basu550n)

!     and finally the BC&OC(n) mode:
      mplus10=1
        call intaeropt4(lchnk, ncol, xrh, irh1, irh2, mplus10, Nnatk, fnbc, camnull, focm, faqm4, &
           bext440n, bext500n, bext550n, bext670n, bext870n,             &
           bebg440n, bebg500n, bebg550n, bebg670n, bebg870n,             &
           bebc440n, bebc500n, bebc550n, bebc670n, bebc870n,             &
           beoc440n, beoc500n, beoc550n, beoc670n, beoc870n,             &
           besu440n, besu500n, besu550n, besu670n, besu870n,             &
           babs440n, babs500n, babs550n, babs670n, babs870n,             &
           bebglt1n, bebggt1n, bebclt1n, bebcgt1n,                       &
           beoclt1n, beocgt1n, bes4lt1n, bes4gt1n,                       &
           backsc550n, babg550n, babc550n, baoc550n, basu550n)

      end do ! iloop

!       Calculation of dry extinction and absorption for all r and for r<0.5um
        do k=1,pver
          do icol=1,ncol

            do i=0,10
              ec550dry_aer(icol,k)  = ec550dry_aer(icol,k)+Nnatk(icol,k,i)*bext550(icol,k,i)
              abs550dry_aer(icol,k) = abs550dry_aer(icol,k)+Nnatk(icol,k,i)*babs550(icol,k,i)
              ec440dry_aer(icol,k)  = ec440dry_aer(icol,k)+Nnatk(icol,k,i)*bext440(icol,k,i)
              abs440dry_aer(icol,k) = abs440dry_aer(icol,k)+Nnatk(icol,k,i)*babs440(icol,k,i)
              ec870dry_aer(icol,k)  = ec870dry_aer(icol,k)+Nnatk(icol,k,i)*bext870(icol,k,i)
              abs870dry_aer(icol,k) = abs870dry_aer(icol,k)+Nnatk(icol,k,i)*babs870(icol,k,i)
              basu550tot(icol,k) = basu550tot(icol,k)+Nnatk(icol,k,i)*basu550(icol,k,i)
              babc550tot(icol,k) = babc550tot(icol,k)+Nnatk(icol,k,i)*babc550(icol,k,i)
              baoc550tot(icol,k) = baoc550tot(icol,k)+Nnatk(icol,k,i)*baoc550(icol,k,i)
              bes4lt1t(icol,k) = bes4lt1t(icol,k)+Nnatk(icol,k,i)*bes4lt1(icol,k,i)
              bebclt1t(icol,k) = bebclt1t(icol,k)+Nnatk(icol,k,i)*bebclt1(icol,k,i)
              beoclt1t(icol,k) = beoclt1t(icol,k)+Nnatk(icol,k,i)*beoclt1(icol,k,i)
            enddo
            do i=11,14
              ec550dry_aer(icol,k)  = ec550dry_aer(icol,k)+Nnatk(icol,k,i)*bext550n(icol,k,i-10)
              abs550dry_aer(icol,k) = abs550dry_aer(icol,k)+Nnatk(icol,k,i)*babs550n(icol,k,i-10)
              ec440dry_aer(icol,k)  = ec440dry_aer(icol,k)+Nnatk(icol,k,i)*bext440n(icol,k,i-10)
              abs440dry_aer(icol,k) = abs440dry_aer(icol,k)+Nnatk(icol,k,i)*babs440n(icol,k,i-10)
              ec870dry_aer(icol,k)  = ec870dry_aer(icol,k)+Nnatk(icol,k,i)*bext870n(icol,k,i-10)
              abs870dry_aer(icol,k) = abs870dry_aer(icol,k)+Nnatk(icol,k,i)*babs870n(icol,k,i-10)
              ba550x(icol,k,i)=babs550n(icol,k,i-10)
              belt1x(icol,k,i)=bebglt1n(icol,k,i-10)
            enddo

!lt1+
            do i=6,7
              bedustlt1(icol,k) = bedustlt1(icol,k) + Nnatk(icol,k,i)*bebglt1(icol,k,i)
            enddo
            do i=8,10
              besslt1(icol,k) = besslt1(icol,k) + Nnatk(icol,k,i)*bebglt1(icol,k,i)
            enddo             
            ec550drylt1_du(icol,k) = bedustlt1(icol,k)
            ec550drylt1_ss(icol,k) = besslt1(icol,k)

!soa: *(1-v_soan) for the sulfate volume fraction of mode 11
            bbclt1xt(icol,k) = Nnatk(icol,k,12)*belt1x(icol,k,12) &
                             + Nnatk(icol,k,14)*belt1x(icol,k,14)*vnbc(icol,k)
!soa + v_soan part of mode 11 for the OC volume fraction of that mode
            boclt1xt(icol,k) = Nnatk(icol,k,13)*belt1x(icol,k,13) &
                             + Nnatk(icol,k,14)*belt1x(icol,k,14)*(1.0_r8-vnbc(icol,k)) 

!soa: *(1-v_soana) for the sulfate volume fraction of mode 1
            ec550drylt1_su(icol,k) = bes4lt1t(icol,k)                         &  ! condensate
                  + Nnatk(icol,k,1)*bebglt1(icol,k,1)*(1.0_r8-v_soana(icol,k))&  ! background, SO4(Ait) mode (1)
                  + Nnatk(icol,k,5)*bebglt1(icol,k,5)                            ! background, SO4(Ait75) mode (5)
            ec550drylt1_bc(icol,k) = bebclt1t(icol,k)+bbclt1xt(icol,k)        &  ! coagulated + n-mode BC (12)
                   + Nnatk(icol,k,2)*bebglt1(icol,k,2)                        &  ! background, BC(Ait) mode (2)
                   + Nnatk(icol,k,4)*bebglt1(icol,k,4)*vaitbc(icol,k)         &  ! background in OC&BC(Ait) mode (4)
                   + Nnatk(icol,k,0)*bebglt1(icol,k,0)                           ! background, BC(ax) mode (0)
!soa + v_soan part of mode 11 for the OC volume fraction of that mode
            ec550drylt1_oc(icol,k) = beoclt1t(icol,k)+boclt1xt(icol,k)        &  ! coagulated + n-mode OC (13)
                   + Nnatk(icol,k,3)*bebglt1(icol,k,3)                        &  ! background, OC(Ait) mode (3)
                   + Nnatk(icol,k,4)*bebglt1(icol,k,4)*(1.0_r8-vaitbc(icol,k))&  ! background in OC&BC(Ait) mode (4)
                   + Nnatk(icol,k,1)*bebglt1(icol,k,1)*v_soana(icol,k)

            ec550drylt1_aer(icol,k) = ec550drylt1_su(icol,k)+ec550drylt1_bc(icol,k) &
              + ec550drylt1_oc(icol,k) + ec550drylt1_ss(icol,k)+ec550drylt1_du(icol,k)
            ec550drylt1_aer(icol,k) = 1.e-3_r8*ec550drylt1_aer(icol,k)
!lt1-

            abs550dry_du(icol,k) = Nnatk(icol,k,6)*babg550(icol,k,6) &
                                 + Nnatk(icol,k,7)*babg550(icol,k,7)
            abs550dry_ss(icol,k) = Nnatk(icol,k,8)*babg550(icol,k,8) &
                                 + Nnatk(icol,k,9)*babg550(icol,k,9) &
                                 + Nnatk(icol,k,10)*babg550(icol,k,10)
!soa: *(1-v_soana) for the sulfate volume fraction of mode 1
            abs550dry_su(icol,k) = basu550tot(icol,k)                   &  ! condensate:w

           + (1.0_r8-v_soana(icol,k))*Nnatk(icol,k,1)*babg550(icol,k,1) &  ! background, SO4(Ait) mode (1)
                                    + Nnatk(icol,k,5)*babg550(icol,k,5)    ! background, SO4(Ait75) mode (5)

!soa: *(1-v_soan) for the sulfate volume fraction
            babc550xt(icol,k) = Nnatk(icol,k,12)*ba550x(icol,k,12)  &
                              + Nnatk(icol,k,14)*ba550x(icol,k,14)*vnbc(icol,k)
!soa 
            baoc550xt(icol,k) = Nnatk(icol,k,13)*ba550x(icol,k,13) &
                              + Nnatk(icol,k,14)*ba550x(icol,k,14)*(1.0_r8-vnbc(icol,k)) 

            abs550dry_bc(icol,k) = babc550tot(icol,k)+babc550xt(icol,k) &     ! coagulated + n-mode BC (12)
                                 + Nnatk(icol,k,2)*babg550(icol,k,2) &        ! background, BC(Ait) mode (2)
                  + vaitbc(icol,k)*Nnatk(icol,k,4)*babg550(icol,k,4) &        ! background in OC&BC(Ait) mode (4)
                                 + Nnatk(icol,k,0)*babg550(icol,k,0)          ! background, BC(ax) mode (0)
!soa 
            abs550dry_oc(icol,k) = baoc550tot(icol,k)+baoc550xt(icol,k) &     ! coagulated + n-mode OC (13)
                  + v_soana(icol,k)*Nnatk(icol,k,1)*babg550(icol,k,1) &       ! SOA fraction of mode 1
                                  + Nnatk(icol,k,3)*babg550(icol,k,3) &       ! background, OC(Ait) mode (3)
          + (1.0_r8-vaitbc(icol,k))*Nnatk(icol,k,4)*babg550(icol,k,4)         ! background in OC&BC(Ait) mode (4)

            deltah=1.e-4_r8*(pint(icol,k+1)-pint(icol,k))/(rhoda(icol,k)*9.8_r8)
            dod550dry(icol) = dod550dry(icol)+ec550dry_aer(icol,k)*deltah
            abs550dry(icol) = abs550dry(icol)+abs550dry_aer(icol,k)*deltah

            ec550dry_aer(icol,k)  = 1.e-3_r8*ec550dry_aer(icol,k)
            abs550dry_aer(icol,k) = 1.e-3_r8*abs550dry_aer(icol,k)
            ec440dry_aer(icol,k)  = 1.e-3_r8*ec440dry_aer(icol,k)
            abs440dry_aer(icol,k) = 1.e-3_r8*abs440dry_aer(icol,k)
            ec870dry_aer(icol,k)  = 1.e-3_r8*ec870dry_aer(icol,k)
            abs870dry_aer(icol,k) = 1.e-3_r8*abs870dry_aer(icol,k)

            abs550dry_bc(icol,k)  = 1.e-3_r8*abs550dry_bc(icol,k)
            abs550dry_oc(icol,k)  = 1.e-3_r8*abs550dry_oc(icol,k)
            abs550dry_su(icol,k)  = 1.e-3_r8*abs550dry_su(icol,k)
            abs550dry_ss(icol,k)  = 1.e-3_r8*abs550dry_ss(icol,k)
            abs550dry_du(icol,k)  = 1.e-3_r8*abs550dry_du(icol,k)

          enddo
        enddo 

        call outfld('ECDRYAER',ec550dry_aer,pcols,lchnk)
        call outfld('ABSDRYAE',abs550dry_aer,pcols,lchnk)
        call outfld('OD550DRY',dod550dry,pcols,lchnk)       ! 2D variable
        call outfld('AB550DRY',abs550dry,pcols,lchnk)       ! 2D variable

        call outfld('ECDRY440',ec440dry_aer,pcols,lchnk)
        call outfld('ABSDR440',abs440dry_aer,pcols,lchnk)
        call outfld('ECDRY870',ec870dry_aer,pcols,lchnk)
        call outfld('ABSDR870',abs870dry_aer,pcols,lchnk)

        call outfld('ECDRYLT1',ec550drylt1_aer,pcols,lchnk)
!       Since we do not have enough look-up table info to take abs550drylt1_aer,
!       instead take out abs550dry for each constituent:
        call outfld('ABSDRYBC',abs550dry_bc,pcols,lchnk)
        call outfld('ABSDRYOC',abs550dry_oc,pcols,lchnk)
        call outfld('ABSDRYSU',abs550dry_su,pcols,lchnk)
        call outfld('ABSDRYSS',abs550dry_ss,pcols,lchnk)
        call outfld('ABSDRYDU',abs550dry_du,pcols,lchnk)

!000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000


      return
end subroutine opticsAtConstRh

