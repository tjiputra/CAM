module opttab_lw

! Purpose: To read in LW look-up tables and calculate optical properties for the aerosols
! for the subroutine interpol*_lw.

!   Based on opttab.F90 and modified for new wavelength bands and look-up tables 
!   by Alf Kirkevaag in January 2014, and extended to include SOA in August 2015.


  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile,  only: iulog
  use opttab
  implicit none

  private 
  save


  ! Interfaces
  public initopt_lw


!     Array bounds in the tabulated optical parameters  
   integer, public, parameter :: nlwbands=16    ! number of aerosol spectral bands in LW
   
  real(r8), public :: ka0(nlwbands)
  real(r8), public :: ka1(nlwbands,10,6,16,6)
  real(r8), public :: ka2to3(nlwbands,10,16,6,2:3)
  real(r8), public :: ka4(nlwbands,10,6,16,6,6)
  real(r8), public :: ka5to10(nlwbands,10,6,6,6,6,5:10)


 contains

subroutine initopt_lw

!---------------------------------------------------------------
!   Modified by Egil Storen/NoSerC July 2002.
!   The sequence of the indices in arrays om1, g1, be1 and ke1
!   (common block /tab1/) has been rearranged to avoid cache
!   problems while running subroutine interpol1. Files also 
!   involved by this modification: interpol1.F and opttab.h.
!   Modified for new aerosol schemes by Alf Kirkevaag in January 
!   2006. Based on opttab.F90 and modified for new wavelength 
!   bands and look-up tables by Alf Kirkevaag in January 2014, 
!   and for SOA in August 2015.
!---------------------------------------------------------------

    use oslo_control, only: oslo_getopts, dir_string_length
!   use shr_kind_mod, only: r8 => shr_kind_r8
!   use inpgraer


!   implicit none

!     Tabulating the 'kcomp'-files to save computing time.

      integer kcomp, iwl, irelh, ictot, ifac, ifbc, ifaq
      integer ifombg, ifbcbg                                  ! soa
      integer ic, ifil, lin, linmax
      real(r8) catot, relh, frac, fabc, fraq, frombg, frbcbg
      real(r8) spabs
      real(r8) rh2(10)
      real(r8) :: eps2 = 1.e-2_r8
      real(r8) :: eps4 = 1.e-4_r8
      real(r8) :: eps6 = 1.e-6_r8
      real(r8) :: eps7 = 1.e-7_r8
      character(len=dir_string_length) :: aerotab_table_dir

      call oslo_getopts(aerotab_table_dir_out = aerotab_table_dir)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      open(40,file=trim(aerotab_table_dir)//'/lwkcomp1.out' &
             ,form="formatted",status="old")
      open(41,file=trim(aerotab_table_dir)//'/lwkcomp2.out' &
             ,form="formatted",status="old")
      open(42,file=trim(aerotab_table_dir)//'/lwkcomp3.out' &
             ,form="formatted",status="old")
      open(43,file=trim(aerotab_table_dir)//'/lwkcomp4.out' &
             ,form="formatted",status="old")
      open(44,file=trim(aerotab_table_dir)//'/lwkcomp5.out' &
             ,form="formatted",status="old")
      open(45,file=trim(aerotab_table_dir)//'/lwkcomp6.out' &
             ,form="formatted",status="old")
      open(46,file=trim(aerotab_table_dir)//'/lwkcomp7.out' &
             ,form="formatted",status="old")
      open(47,file=trim(aerotab_table_dir)//'/lwkcomp8.out' &
             ,form="formatted",status="old")
      open(48,file=trim(aerotab_table_dir)//'/lwkcomp9.out' &
             ,form="formatted",status="old")
      open(49,file=trim(aerotab_table_dir)//'/lwkcomp10.out'& 
             ,form="formatted",status="old")
      open(50,file=trim(aerotab_table_dir)//'/lwkcomp0.out'& 
             ,form="formatted",status="old")
 
!     Skipping the header-text in all input files (Later: use it to check AeroTab - CAM5-Oslo consistency!)
      do ifil = 40,50
        call checkTableHeader (ifil)
      enddo

!     Then reading in the look-up table entries for each file (lwkcomp*.out)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!       Mode 0, BC(ax)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

        ifil = 11
        linmax=nlwbands
        do lin = 1,linmax

          read(39+ifil,996) kcomp, iwl, relh, spabs

          ka0(iwl)=spabs  ! unit m^2/g

!      write(*,*) 'kcomp, ka =', kcomp, ka0(iwl)

        end do

    do iwl=1,nlwbands
     if(ka0(iwl)<=0.0_r8) then
      write(iulog,*) 'ka0 =', iwl, ka0(iwl)
      write(iulog,*) 'Error in initialization of ka0'
      stop
     endif
    enddo

        write(iulog,*)'lw mode 0 ok' 


!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!soa       New mode 1 (SO4 + condensate from H2SO4 and SOA)
!       Mode 1 (H2SO4 + condesate from H2SO4 and SOA)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

        ifil = 1
!soa        linmax=nlwbands*10*16*6
        linmax=nlwbands*10*6*16*6
        do lin = 1,linmax

          read(39+ifil,997) kcomp, iwl, relh, frombg, catot, frac, spabs

       	  do ic=1,10
	   if(abs(relh-rh(ic))<eps4) then
	    irelh=ic
	    goto 121
	   endif
	  end do
  121     continue

 	  do ic=1,6
	   if(abs(frombg-fombg(ic))<eps4) then
	    ifombg=ic
	    goto 122
	   endif
	  end do
  122     continue

 	  do ic=1,16
	   if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
	    ictot=ic
	    goto 131
	   endif
	  end do
  131     continue

 	  do ic=1,6
	   if(abs(frac-fac(ic))<eps4) then
	    ifac=ic
	    goto 141
	   endif
	  end do
  141     continue

          ka1(iwl,irelh,ifombg,ictot,ifac)=spabs  ! unit m^2/g

!      write(*,*) 'kcomp, ka =', kcomp, ka1(iwl,irelh,ifombg,ictot,ifac)
!      if(ifil==1) write(iulog,*) 'iwl,irelh,ifombg,ictot,ifac,ka =', &
!                  iwl,irelh,ictot,ifac,ka1(iwl,irelh,ifombg,ictot,ifac)

        end do  ! lin

    do iwl=1,nlwbands
    do irelh=1,10
    do ictot=1,16
    do ifac=1,6
     if(ka1(iwl,irelh,ifombg,ictot,ifac)<=0.0_r8) then
      write(iulog,*) 'ka1 =', iwl, irelh, ifombg, ictot, ifac, ka1(iwl,irelh,ifombg,ictot,ifac)
      write(iulog,*) 'Error in initialization of ka1'
      stop
     endif
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'lw new mode 1 ok' 


!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!soa       Modes 2 to 3 (BC/OC + condensate from H2SO4)
!       Modes 2 to 3 (BC or OC + condensate from H2SO4 and SOA)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

!soa      linmax = nlwbands*10*16 
      linmax = nlwbands*10*16*6 
      do ifil = 2,3
        do lin = 1,linmax 

          read(39+ifil,994) kcomp, iwl, relh, catot, frac, spabs

       	  do ic=1,10
	   if(abs(relh-rh(ic))<eps4) then
	    irelh=ic
	    goto 61
	   endif
	  end do
   61     continue

 	  do ic=1,16
	   if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
	    ictot=ic
	    goto 71
	   endif
	  end do
   71     continue

 	  do ic=1,6
	   if(abs(frac-fac(ic))<eps4) then
	    ifac=ic
	    goto 72
	   endif
	  end do
   72     continue

          ka2to3(iwl,irelh,ictot,ifac,kcomp)=spabs  ! unit m^2/g

!      write(*,*) 'kcomp, ka =', kcomp, ka2to3(iwl,irelh,ictot,ifac,kcomp)
!      if(ifil==2) write(iulog,*) 'iwl,irelh,ictot,ifac,kcomp,ka =', &
!                  iwl,irelh,ictot,kcomp,ka2to3(iwl,irelh,ictot,ifac,kcomp)

        end do  ! lin
      end do    ! ifil

    do kcomp=2,3
    do iwl=1,nlwbands
    do irelh=1,10
    do ictot=1,16
     if(ka2to3(iwl,irelh,ictot,ifac,kcomp)<=0.0_r8) then
      write(iulog,*) 'ka2to3 =', iwl, irelh, ictot, ifac, ka2to3(iwl,irelh,ictot,ifac,kcomp)
      write(iulog,*) 'Error in initialization of ka2to3'
      stop
     endif
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'lw mode 2-3 ok' 


!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!soa       Mode 4 (BC&OC + condesate from H2SO4 + wetphase (NH4)2SO4)
!       Mode 4 (BC&OC + condesate from H2SO4 and SOA + wetphase (NH4)2SO4)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

        ifil = 4
!soa        linmax = nlwbands*10*16*6*6 
        linmax = nlwbands*10*6*16*6*6 
        do lin = 1,linmax

          read(39+ifil,995) kcomp, iwl, relh, frbcbg, catot, frac, fraq, spabs

       	  do ic=1,10
	   if(abs(relh-rh(ic))<eps4) then
	    irelh=ic
	    goto 81
	   endif
	  end do
   81     continue

 	  do ic=1,6
	   if(abs(frbcbg-fbcbg(ic))<eps4) then
	    ifbcbg=ic
	    goto 92
	   endif
	  end do
   92     continue

 	  do ic=1,16
	   if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
	    ictot=ic
	    goto 91
	   endif
	  end do
   91     continue

 	  do ic=1,6
	   if(abs(frac-fac(ic))<eps4) then
	    ifac=ic
	    goto 101
	   endif
	  end do
  101     continue

	  do ic=1,6
	   if(abs(fraq-faq(ic))<eps4) then
	    ifaq=ic
	    goto 111
	   endif
	  end do
  111     continue

          ka4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)=spabs  ! unit m^2/g

!      write(*,*) 'kcomp, ka =', kcomp, ka4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
        end do

    do iwl=1,nlwbands
    do irelh=1,10
    do ictot=1,16
    do ifac=1,6
    do ifaq=1,6
     if(ka4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)<=0.0_r8) then
      write(iulog,*) 'ka4 =', iwl, irelh, ifbcbg, ictot, ifac, ifaq, ka4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
      write(iulog,*) 'Error in initialization of ka4'
      stop
     endif
    enddo
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'lw mode 4 ok' 


!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!       Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      linmax = nlwbands*10*6*6*6*6 
      do ifil = 5,10
        do lin = 1,linmax   

          read(39+ifil,993) kcomp, iwl, relh, catot, frac, fabc, fraq, spabs

       	  do ic=1,10
	   if(abs(relh-rh(ic))<eps4) then
	    irelh=ic
	    goto 11
	   endif
	  end do
   11     continue

 	  do ic=1,6
	   if(abs(catot-cat(kcomp,ic))<eps6) then
	    ictot=ic
	    goto 21
	   endif
	  end do
   21     continue

 	  do ic=1,6
	   if(abs(frac-fac(ic))<eps4) then
	    ifac=ic
	    goto 31
	   endif
	  end do
   31     continue

 	  do ic=1,6
	   if(abs(fabc-fbc(ic))<eps4) then
	    ifbc=ic
	    goto 41
	   endif
	  end do
   41     continue

	  do ic=1,6
	   if(abs(fraq-faq(ic))<eps4) then
	    ifaq=ic
	    goto 51
	   endif
	  end do
   51     continue

          ka5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=spabs  ! unit m^2/g

!      write(*,*) 'kcomp, ka =', kcomp, ka5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp) 
        end do
      end do


    do kcomp=5,10
    do iwl=1,nlwbands
    do irelh=1,10
    do ictot=1,6
    do ifac=1,6
    do ifaq=1,6
     if(ka5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)<=0.0_r8) then
      write(iulog,*) 'ka5to10 =', iwl, irelh, ictot, ifac, ifbc, ifaq, ka5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)
      write(iulog,*) 'Error in initialization of ka5to10'
      stop
     endif
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'lw mode 5-10 ok' 

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

 993     format(2I3,f8.3,3(x,e10.3),f7.2,x,e12.5)      ! 5-10
 994     format(2I3,f8.3,2(x,e10.3),x,e12.5)           ! 2-3
 995     format(2I3,f8.3,3(x,e10.3),f7.2,x,e12.5)      ! 4
 996     format(2I3,f8.3,x,e12.5)                      ! 0   
 997     format(2I3,f8.3,3(x,e10.3),x,e12.5)           ! 1


      do ifil=40,50
        close (ifil)
      end do 
      return
end subroutine initopt_lw


end module opttab_lw
