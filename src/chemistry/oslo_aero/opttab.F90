module opttab

! Purpose: To read in SW look-up tables and calculate optical properties for the aerosols
!   For subroutine interpol.

!   Modified for new wavelength bands and look-up tables 
!   by Alf Kirkevaag in December 2013.
!   Updated for reading inout files with extra header info - Alf Kirkevaag, May 2015
!   Extended for new SOA treatment - Alf Kirkevaag, August 2015


  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile,  only: iulog
  implicit none

  private 
  save


  ! Interfaces
  public initopt


!     Array bounds in the tabulated optical parameters  
!CAM4-Oslo   integer, public, parameter :: nbands=12    ! number of aerosol spectral bands
   integer, public, parameter :: nbands=14    ! number of aerosol spectral bands in SW
   integer, public, parameter :: nbmp1=11     ! number of first non-background mode

!soa
   real(r8), public, dimension(6) :: fombg = (/ 0.0_r8, 0.1_r8,  0.3_r8, 0.5_r8, 0.7_r8, 0.999_r8    /)
   real(r8), public, dimension(6) :: fbcbg = (/ 0.0_r8, 0.1_r8,  0.3_r8, 0.5_r8, 0.7_r8, 0.999_r8    /)
!soa
   real(r8), public, dimension(6) :: fac = (/ 0.0_r8, 0.1_r8,  0.3_r8, 0.5_r8, 0.7_r8, 0.999_r8    /)
   real(r8), public, dimension(6) :: fbc = (/ 0.0_r8, 0.01_r8, 0.1_r8, 0.3_r8, 0.7_r8, 0.999_r8    /)
   real(r8), public, dimension(6) :: faq = (/ 0.0_r8, 0.25_r8, 0.5_r8, 0.75_r8,0.85_r8,1.0_r8      /)
   real(r8), public, dimension(10) :: rh = (/ 0.0_r8, 0.37_r8, 0.47_r8,0.65_r8,0.75_r8,            &
!test_old                                      0.8_r8, 0.85_r8, 0.9_r8, 0.95_r8,0.98_r8             /)
                                      0.8_r8, 0.85_r8, 0.9_r8, 0.95_r8,0.995_r8             /)
   real(r8), public, dimension(9) :: S   = (/ 1.0005_r8, 1.001_r8,  1.0015_r8, 1.002_r8,           &
                                      1.0025_r8, 1.003_r8,  1.005_r8,  1.008_r8,  1.01_r8  /)  


!        real(r8), public, dimension(5:10,6) :: cat = reshape ( (/ &
!   1.e-10_r8, 1.e-10_r8, 1.e-10_r8, 1.e-10_r8, 1.e-10_r8, 1.e-10_r8, & 
!   5.e-4_r8 , 0.01_r8  , 0.02_r8  , 1.e-4_r8 , 0.005_r8 , 0.02_r8  , &
!   2.e-3_r8 , 0.05_r8  , 0.1_r8   , 6.e-4_r8 , 0.025_r8 , 0.1_r8   , &
!   0.01_r8  , 0.2_r8   , 0.5_r8   , 2.5e-3_r8, 0.1_r8   , 0.5_r8   , &
!   0.04_r8  , 0.8_r8   , 2.0_r8   , 1.e-2_r8 , 0.4_r8   , 2.0_r8   , &
!   0.15_r8  , 4.0_r8   , 8.0_r8   , 3.5e-2_r8, 2.0_r8   , 8.0_r8     &
!                                                      /), (/6,6/) )

!  With look-up tables for the Salter et al. 2015) sea-salt parametrisation: 
        real(r8), public, dimension(5:10,6) :: cat = reshape ( (/ &
   1.e-10_r8, 1.e-10_r8, 1.e-10_r8, 1.e-10_r8, 1.e-10_r8, 1.e-10_r8, & 
   5.e-4_r8 , 0.01_r8  , 0.02_r8  , 5.e-4_r8 , 0.01_r8  , 0.02_r8  , &
   2.e-3_r8 , 0.05_r8  , 0.1_r8   , 2.e-3_r8 , 0.05_r8  , 0.1_r8   , &
   0.01_r8  , 0.2_r8   , 0.5_r8   , 0.01_r8  , 0.2_r8   , 0.5_r8   , &
   0.04_r8  , 0.8_r8   , 2.0_r8   , 0.04_r8  , 0.8_r8   , 2.0_r8   , &
   0.15_r8  , 4.0_r8   , 8.0_r8   , 0.15_r8  , 4.0_r8   , 8.0_r8     &
                                                      /), (/6,6/) )

        real(r8), public, dimension(4,16) :: cate = reshape ( (/ &
   1.e-10_r8, 1.e-10_r8, 1.e-10_r8, 1.e-10_r8*1.904e-3_r8, & 
   1.e-5_r8 , 1.e-5_r8 , 1.e-4_r8 , 0.01_r8*1.904e-3_r8  , &
   2.e-5_r8 , 2.e-5_r8 , 2.e-4_r8 , 0.05_r8*1.904e-3_r8  , &
   4.e-5_r8 , 4.e-5_r8 , 4.e-4_r8 , 0.1_r8*1.904e-3_r8   , &
   8.e-5_r8 , 8.e-5_r8 , 8.e-4_r8 , 0.2_r8*1.904e-3_r8   , &
   1.5e-4_r8, 1.5e-4_r8, 1.5e-3_r8, 0.4_r8*1.904e-3_r8   , &
   3.e-4_r8 , 3.e-4_r8 , 3.e-3_r8 , 0.7_r8*1.904e-3_r8   , &
   6.e-4_r8 , 6.e-4_r8 , 6.e-3_r8 , 1.0_r8*1.904e-3_r8   , &
   1.2e-3_r8, 1.2e-3_r8, 1.2e-2_r8, 1.5_r8*1.904e-3_r8   , &
   2.5e-3_r8, 2.5e-3_r8, 2.5e-2_r8, 2.5_r8*1.904e-3_r8   , &
   5.e-3_r8 , 5.e-3_r8 , 5.e-2_r8 , 5.0_r8*1.904e-3_r8   , &
   1.e-2_r8 , 1.e-2_r8 , 0.1_r8   , 10.0_r8*1.904e-3_r8  , &
   2.e-2_r8 , 2.e-2_r8 , 0.2_r8   , 25.0_r8*1.904e-3_r8  , &
   4.e-2_r8 , 4.e-2_r8 , 0.4_r8   , 50.0_r8*1.904e-3_r8  , &
   8.e-2_r8 , 8.e-2_r8 , 0.8_r8   , 100.0_r8*1.904e-3_r8 , &
   0.15_r8  , 0.15_r8  , 1.5_r8   , 500.0_r8*1.904e-3_r8 /), (/4,16/) )
   
  real(r8), public :: om1(nbands,10,6,16,6)
  real(r8), public :: g1 (nbands,10,6,16,6)
  real(r8), public :: be1(nbands,10,6,16,6)
  real(r8), public :: ke1(nbands,10,6,16,6)
!
  real(r8), public :: om2to3(nbands,10,16,6,2:3)
  real(r8), public :: g2to3 (nbands,10,16,6,2:3)
  real(r8), public :: be2to3(nbands,10,16,6,2:3)
  real(r8), public :: ke2to3(nbands,10,16,6,2:3)

  real(r8), public :: om4(nbands,10,6,16,6,6)
  real(r8), public :: g4 (nbands,10,6,16,6,6)
  real(r8), public :: be4(nbands,10,6,16,6,6)
  real(r8), public :: ke4(nbands,10,6,16,6,6)

  real(r8), public :: om0(nbands)
  real(r8), public :: g0(nbands)
  real(r8), public :: be0(nbands)
  real(r8), public :: ke0(nbands)

  real(r8), public :: om5to10(nbands,10,6,6,6,6,5:10)
  real(r8), public :: g5to10(nbands,10,6,6,6,6,5:10)
  real(r8), public :: be5to10(nbands,10,6,6,6,6,5:10)
  real(r8), public :: ke5to10(nbands,10,6,6,6,6,5:10)

  real(r8), public :: e,eps
  parameter (e=2.718281828_r8, eps=1.0e-30_r8)
  
 contains

subroutine initopt

!---------------------------------------------------------------
!   Modified by Egil Storen/NoSerC July 2002.
!   The sequence of the indices in arrays om1, g1, be1 and ke1
!   (common block /tab1/) has been rearranged to avoid cache
!   problems while running subroutine interpol1. Files also 
!   involved by this modification: interpol1.F and opttab.h.
!   Modified for new aerosol schemes by Alf Kirkevaag in January 
!   2006. Modified for new wavelength bands and look-up tables 
!   by Alf Kirkevaag in December 2013, and for SOA in August 2015.
!---------------------------------------------------------------

      use oslo_control, only : oslo_getopts, dir_string_length
!   use shr_kind_mod, only: r8 => shr_kind_r8
!   use inpgraer


!   implicit none

!     Tabulating the 'kcomp'-files to save computing time.


      integer kcomp, iwl, irelh, ictot, ifac, ifbc, ifaq
      integer ifombg, ifbcbg                                  ! soa
      integer ic, ifil, lin, linmax
      real(r8) catot, relh, frac, fabc, fraq, frombg, frbcbg
      real(r8) ssa, ass, ext, spext
      real(r8) :: eps2 = 1.e-2_r8
      real(r8) :: eps4 = 1.e-4_r8
      real(r8) :: eps6 = 1.e-6_r8
      real(r8) :: eps7 = 1.e-7_r8
      character(len=dir_string_length) :: aerotab_table_dir

      call oslo_getopts(aerotab_table_dir_out= aerotab_table_dir)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      open(40,file=trim(aerotab_table_dir)//'/kcomp1.out' &
             ,form='formatted',status='old')
      open(41,file=trim(aerotab_table_dir)//'/kcomp2.out' &
             ,form='formatted',status='old')
      open(42,file=trim(aerotab_table_dir)//'/kcomp3.out' &
             ,form='formatted',status='old')
      open(43,file=trim(aerotab_table_dir)//'/kcomp4.out' &
             ,form='formatted',status='old')
      open(44,file=trim(aerotab_table_dir)//'/kcomp5.out' &
             ,form='formatted',status='old')
      open(45,file=trim(aerotab_table_dir)//'/kcomp6.out' &
             ,form='formatted',status='old')
      open(46,file=trim(aerotab_table_dir)//'/kcomp7.out' &
             ,form='formatted',status='old')
      open(47,file=trim(aerotab_table_dir)//'/kcomp8.out' &
             ,form='formatted',status='old')
      open(48,file=trim(aerotab_table_dir)//'/kcomp9.out' &
             ,form='formatted',status='old')
      open(49,file=trim(aerotab_table_dir)//'/kcomp10.out'& 
             ,form='formatted',status='old')
      open(50,file=trim(aerotab_table_dir)//'/kcomp0.out'& 
             ,form='formatted',status='old')
 
!     Skipping the header-text in all input files (Later: use it to check AeroTab - CAM5-Oslo consistency!)
      do ifil = 40,50
        call checkTableHeader (ifil)
      enddo

!     Then reading in the look-up table entries for each file (kcomp*.out)
     
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!       Mode 0, BC(ax)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

        ifil = 11
        linmax=nbands
        do lin = 1,linmax

          read(39+ifil,996) kcomp, iwl, relh, &
                          ssa, ass, ext, spext
          om0(iwl)=ssa    
          g0 (iwl)=ass
          be0(iwl)=ext    ! unit km^-1
          ke0(iwl)=spext  ! unit m^2/g

!      write(iulog,*) 'kcomp, om =', kcomp, om0(iwl)
!      write(iulog,*) 'kcomp, g  =', kcomp, g0(iwl)
!      write(iulog,*) 'kcomp, be =', kcomp, be0(iwl)
!      write(iulog,*) 'kcomp, ke =', kcomp, ke0(iwl)

        end do

    do iwl=1,nbands
     if(be0(iwl)<=0.0_r8) then
      write(iulog,*) 'be0 =', iwl, be0(iwl)
      write(iulog,*) 'Error in initialization of be0'
      stop
     endif
    enddo

        write(iulog,*)'mode 0 ok' 


!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!       Mode 1 (H2SO4 and SOA + condesate from H2SO4 and SOA)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

!soa      linmax = nbands*10*16   ! 14*10*16
      linmax = nbands*10*6*16*6   ! 14*10*6*16*6
      do ifil = 1,1
        do lin = 1,linmax 

          read(39+ifil,995) kcomp, iwl, relh, frombg, catot, frac, &
                          ssa, ass, ext, spext

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
	   if(abs(frombg-fombg(ic))<eps4) then
	    ifombg=ic
	    goto 72
	   endif
	  end do
   72     continue

 	  do ic=1,6
	   if(abs(frac-fac(ic))<eps4) then
	    ifac=ic
	    goto 73
	   endif
	  end do
   73     continue

          om1(iwl,irelh,ifombg,ictot,ifac)=ssa    
          g1 (iwl,irelh,ifombg,ictot,ifac)=ass
          be1(iwl,irelh,ifombg,ictot,ifac)=ext    ! unit km^-1
          ke1(iwl,irelh,ifombg,ictot,ifac)=spext  ! unit m^2/g

!      write(iulog,*) 'kcomp, om =', kcomp, om1(iwl,irelh,ifombg,ictot,ifac)
!      write(iulog,*) 'kcomp, g  =', kcomp, g1(iwl,irelh,ifombg,ictot,ifac)
!      write(iulog,*) 'kcomp, be =', kcomp, be1(iwl,irelh,ifombg,ictot,ifac)
!      write(iulog,*) 'kcomp, ke =', kcomp, ke1(iwl,irelh,ifombg,ictot,ifac)

        end do  ! lin
      end do    ! ifil

    do kcomp=1,1
    do iwl=1,nbands
    do irelh=1,10
    do ifombg=1,6   !soa
    do ictot=1,16
    do ifac=1,6     !soa
     if(be1(iwl,irelh,ifombg,ictot,ifac)<=0.0_r8) then
      write(iulog,*) 'be1 =', iwl, irelh, ifombg, ictot, be1(iwl,irelh,ifombg,ictot,ifac)
      write(iulog,*) 'Error in initialization of be1'
      stop
     endif
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'mode 1 ok' 

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!       Modes 2 to 3 (BC/OC + condensate from H2SO4 and SOA)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

!soa        linmax=nbands*10*16*6
        linmax=nbands*10*16*6
      do ifil = 2,3
        do lin = 1,linmax

          read(39+ifil,997) kcomp, iwl, relh, catot, frac, &
                          ssa, ass, ext, spext

       	  do ic=1,10
	   if(abs(relh-rh(ic))<eps4) then
	    irelh=ic
	    goto 121
	   endif
	  end do
  121     continue

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

          om2to3(iwl,irelh,ictot,ifac,kcomp)=ssa    
          g2to3 (iwl,irelh,ictot,ifac,kcomp)=ass
          be2to3(iwl,irelh,ictot,ifac,kcomp)=ext    ! unit km^-1
          ke2to3(iwl,irelh,ictot,ifac,kcomp)=spext  ! unit m^2/g

!      write(iulog,*) 'kcomp, om =', kcomp, om2to3(iwl,irelh,ictot,ifac,kcomp)
!      write(iulog,*) 'kcomp, g  =', kcomp, g2to3(iwl,irelh,ictot,ifac,kcomp)
!      write(iulog,*) 'kcomp, be =', kcomp, be2to3(iwl,irelh,ictot,ifac,kcomp)
!      write(iulog,*) 'kcomp, ke =', kcomp, ke2to3(iwl,irelh,ictot,ifac,kcomp)

        end do  ! lin
      enddo     ! ifil

    do kcomp=2,3
    do iwl=1,nbands
    do irelh=1,10
    do ictot=1,16
    do ifac=1,6
     if(be2to3(iwl,irelh,ictot,ifac,kcomp)<=0.0_r8) then
      write(iulog,*) 'be2to3 =', iwl, irelh, ictot, ifac, be2to3(iwl,irelh,ictot,ifac,kcomp)
      write(iulog,*) 'Error in initialization of be2to3'
      stop
     endif
    enddo
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'modes 2-3 ok' 

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!soa       Mode 4 (BC&OC + condesate from H2SO4 + wet phase (NH4)2SO4)
!       Mode 4 (BC&OC + condensate from H2SO4 and SOA + wet phase (NH4)2SO4)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

        ifil = 4
!soa        linmax = nbands*10*16*6*6 
        linmax = nbands*10*6*16*6*6 
        do lin = 1,linmax

          read(39+ifil,993) kcomp, iwl, relh, frbcbg, catot, frac, fraq, &
                          ssa, ass, ext, spext

       	  do ic=1,10
	   if(abs(relh-rh(ic))<eps4) then
	    irelh=ic
	    goto 81
	   endif
	  end do
   81     continue

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

 	  do ic=1,6
	   if(abs(frbcbg-fbcbg(ic))<eps4) then
	    ifbcbg=ic
	    goto 112
	   endif
	  end do
  112     continue

          om4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)=ssa    
          g4 (iwl,irelh,ifbcbg,ictot,ifac,ifaq)=ass
          be4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)=ext    ! unit km^-1
          ke4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)=spext  ! unit m^2/g

!      write(iulog,*) 'kcomp, om =', kcomp, om4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
!      write(iulog,*) 'kcomp, g  =', kcomp, g4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
!      write(iulog,*) 'kcomp, be =', kcomp, be4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
!      write(iulog,*) 'kcomp, ke =', kcomp, ke4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
        end do

    do iwl=1,nbands
    do irelh=1,10
    do ifbcbg=1,6   !soa
    do ictot=1,16
    do ifac=1,6
    do ifaq=1,6
     if(be4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)<=0.0_r8) then
      write(iulog,*) 'be4 =', iwl, irelh, ifbcbg, ictot, ifac, ifaq, be4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
      write(iulog,*) 'Error in initialization of be4'
      stop
     endif
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'mode 4 ok' 


!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!       Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      linmax = nbands*10*6*6*6*6     ! 14*10*6*6*6*6
      do ifil = 5,10
        do lin = 1,linmax   

          read(39+ifil,993) kcomp, iwl, relh, catot, frac, fabc, fraq, &
                          ssa, ass, ext, spext

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

          om5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=ssa    
          g5to10 (iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=ass
          be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=ext    ! unit km^-1
          ke5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=spext  ! unit m^2/g

!      write(iulog,*) 'kcomp, om =', kcomp, om5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp) 
!      write(iulog,*) 'kcomp, g  =', kcomp, g5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp) 
!      write(iulog,*) 'kcomp, be =', kcomp, be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp) 
!      write(iulog,*) 'kcomp, ke =', kcomp, ke5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp) 
        end do
      end do


    do kcomp=5,10
    do iwl=1,nbands
    do irelh=1,10
    do ictot=1,6
    do ifac=1,6
    do ifaq=1,6
     if(be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)<=0.0_r8) then
      write(iulog,*) 'be5to10 =', iwl, irelh, ictot, ifac, ifbc, ifaq, be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)
      write(iulog,*) 'Error in initialization of be5to10'
      stop
     endif
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'modes 5-10 ok' 


!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc


  993 format(2I3,f8.3,3(x,e10.3),f7.2,4(x,e12.5))
  995 format(2I3,f8.3,3(x,e10.3),4(x,e12.5))
  996 format(2I3,f8.3,4(x,e12.5))
  997 format(2I3,f8.3,2(x,e10.3),4(x,e12.5))


      do ifil=40,50
        close (ifil)
      end do 
      return
end subroutine initopt


end module opttab

