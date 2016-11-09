
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol )
 
      use ppgrid, only : pcols, pver


      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol*pver)
      real(r8), intent(inout) :: rate(ncol*pver,max(1,rxntot))

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      integer   ::  offset
      real(r8)  :: itemp(ncol*pver)
      real(r8)  :: exp_fac(ncol*pver)
      real(r8)  :: ko(ncol*pver)
      real(r8)  :: kinf(ncol*pver)

      rate(:,156) = 8.00e-14_r8
      rate(:,157) = 3.90e-17_r8
      rate(:,160) = 4.20e-13_r8
      rate(:,161) = 8.50e-2_r8
      rate(:,162) = 1.30e-16_r8
      rate(:,164) = 1.00e-20_r8
      rate(:,165) = 2.58e-04_r8
      rate(:,173) = 1.20e-10_r8
      rate(:,174) = 2.07e-10_r8
      rate(:,175) = 1.204e-10_r8
      rate(:,176) = 2.088e-10_r8
      rate(:,177) = 1.17e-10_r8
      rate(:,178) = 4.644e-11_r8
      rate(:,179) = 7.65e-11_r8
      rate(:,180) = 1.794e-10_r8
      rate(:,181) = 1.30e-10_r8
      rate(:,182) = 2.607e-10_r8
      rate(:,183) = 1.80e-10_r8
      rate(:,184) = 9.75e-11_r8
      rate(:,185) = 4.50e-11_r8
      rate(:,186) = 9.90e-11_r8
      rate(:,187) = 1.20e-10_r8
      rate(:,188) = 4.62e-10_r8
      rate(:,189) = 2.57e-10_r8
      rate(:,190) = 2.14e-11_r8
      rate(:,191) = 1.90e-10_r8
      rate(:,192) = 1.31e-10_r8
      rate(:,193) = 3.50e-11_r8
      rate(:,194) = 9.00e-12_r8
      rate(:,195) = 1.20e-10_r8
      rate(:,196) = 9.90e-11_r8
      rate(:,197) = 3.30e-12_r8
      rate(:,198) = 9.00e-11_r8
      rate(:,199) = 3.00e-11_r8
      rate(:,201) = 1.20e-10_r8
      rate(:,204) = 7.20e-11_r8
      rate(:,205) = 6.90e-12_r8
      rate(:,206) = 1.60e-12_r8
      rate(:,210) = 1.80e-12_r8
      rate(:,213) = 1.80e-12_r8
      rate(:,224) = 5.00e-12_r8
      rate(:,225) = 7.00e-13_r8
      rate(:,226) = 5.00e-11_r8
      rate(:,243) = 1.00e-11_r8
      rate(:,244) = 2.20e-11_r8
      rate(:,245) = 3.50e-12_r8
      rate(:,272) = 1.70e-13_r8
      rate(:,323) = 4.00e-13_r8
      rate(:,335) = 1.00e-14_r8
      rate(:,338) = 7.00e-13_r8
      rate(:,341) = 2.00e-13_r8
      rate(:,342) = 6.80e-14_r8
      rate(:,351) = 1.00e-12_r8
      rate(:,352) = 1.00e-11_r8
      rate(:,353) = 1.15e-11_r8
      rate(:,356) = 4.00e-14_r8
      rate(:,374) = 3.00e-12_r8
      rate(:,377) = 6.7e-13_r8
      rate(:,379) = 5.40e-11_r8
      rate(:,382) = 3.5e-13_r8
      rate(:,393) = 2.40e-12_r8
      rate(:,396) = 1.40e-11_r8
      rate(:,399) = 5.00e-12_r8
      rate(:,407) = 2.0e-12_r8
      rate(:,415) = 4.e-11_r8
      rate(:,416) = 4.e-11_r8
      rate(:,417) = 2.40e-12_r8
      rate(:,418) = 2.40e-12_r8
      rate(:,422) = 1.3e-11_r8
      rate(:,425) = 1.40e-11_r8
      rate(:,426) = 1.40e-11_r8
      rate(:,430) = 2.40e-12_r8
      rate(:,432) = 1.4e-11_r8
      rate(:,434) = 4.0e-11_r8
      rate(:,435) = 7.0e-11_r8
      rate(:,436) = 1.0e-10_r8
      rate(:,439) = 2.40e-12_r8
      rate(:,445) = 3.50e-12_r8
      rate(:,446) = 6.7e-12_r8
      rate(:,448) = 1.6e-12_r8
      rate(:,457) = 2.1e-12_r8
      rate(:,458) = 2.8e-13_r8
      rate(:,469) = 4.7e-11_r8
      rate(:,487) = 1.7e-11_r8
      rate(:,488) = 8.4e-11_r8
      rate(:,495) = 2.8e-13_r8
      rate(:,497) = 2.0e-10_r8
      rate(:,499) = 1.2e-14_r8
      rate(:,501) = 1.9e-11_r8
      rate(:,505) = 3.3e-11_r8
      rate(:,506) = 2.3e-11_r8
      rate(:,507) = 5.7e-11_r8
      rate(:,508) = 1e-12_r8
      rate(:,512) = 3.4e-11_r8
      rate(:,516) = 2.4e-12_r8
      rate(:,517) = 2.0e-11_r8
      rate(:,518) = 2.0e-11_r8
      rate(:,520) = 1.2e-14_r8
      rate(:,521) = 1.34e-11_r8
      rate(:,522) = 1.34e-11_r8
      rate(:,531) = 6.60E-11_r8
      rate(:,532) = 2.30E-12_r8
      rate(:,533) = 1.20E-11_r8
      rate(:,537) = 1.40E-11_r8
      rate(:,538) = 2.80E-11_r8
      rate(:,539) = 5.70E-11_r8
      rate(:,540) = 1.90E-12_r8
      rate(:,546) = 6.34e-8_r8
      rate(:,547) = 6.34e-8_r8
      rate(:,569) = 9.0e-10_r8
      rate(:,570) = 1.0e-10_r8
      rate(:,571) = 4.4e-10_r8
      rate(:,572) = 4.0e-10_r8
      rate(:,573) = 2.0e-10_r8
      rate(:,574) = 1.0e-12_r8
      rate(:,575) = 6.0e-11_r8
      rate(:,576) = 5.0e-16_r8
 
      do n = 1,pver
        offset = (n-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,n)
      end do

      exp_fac(:) = exp( -2060._r8 * itemp(:) )
      rate(:,154) = 8.00e-12_r8 * exp_fac(:)
      rate(:,166) = 8.00e-12_r8 * exp_fac(:)
      rate(:,158) = 1.80e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:,159) = 3.50e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:,163) = 3.60e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:,167) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:,168) = 2.64e-11_r8 * exp_fac(:)
      rate(:,169) = 6.60e-12_r8 * exp_fac(:)
      rate(:,170) = 1.63e-10_r8 * exp( 60._r8 * itemp(:) )
      exp_fac(:) = exp( 20._r8 * itemp(:) )
      rate(:,171) = 7.26e-11_r8 * exp_fac(:)
      rate(:,172) = 4.64e-11_r8 * exp_fac(:)
      rate(:,200) = 1.08e-10_r8 * exp( 105._r8 * itemp(:) )
      exp_fac(:) = exp( -470._r8 * itemp(:) )
      rate(:,203) = 1.40e-10_r8 * exp_fac(:)
      rate(:,221) = 1.40e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( 180._r8 * itemp(:) )
      rate(:,207) = 1.80e-11_r8 * exp_fac(:)
      rate(:,333) = 4.20e-12_r8 * exp_fac(:)
      rate(:,361) = 4.20e-12_r8 * exp_fac(:)
      rate(:,366) = 4.20e-12_r8 * exp_fac(:)
      rate(:,386) = 4.20e-12_r8 * exp_fac(:)
      rate(:,413) = 4.40e-12_r8 * exp_fac(:)
      rate(:,414) = 4.40e-12_r8 * exp_fac(:)
      rate(:,502) = 4.2e-12_r8 * exp_fac(:)
      rate(:,509) = 4.2e-12_r8 * exp_fac(:)
      rate(:,513) = 4.2e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -940._r8 * itemp(:) )
      rate(:,208) = 1.70e-12_r8 * exp_fac(:)
      rate(:,222) = 1.70e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 250._r8 * itemp(:) )
      rate(:,209) = 4.80e-11_r8 * exp_fac(:)
      rate(:,282) = 1.70e-11_r8 * exp_fac(:)
      rate(:,212) = 2.80e-12_r8 * exp( -1800._r8 * itemp(:) )
      rate(:,214) = 1.60e-11_r8 * exp( -4570._r8 * itemp(:) )
      exp_fac(:) = exp( 200._r8 * itemp(:) )
      rate(:,215) = 3.00e-11_r8 * exp_fac(:)
      rate(:,290) = 5.50e-12_r8 * exp_fac(:)
      rate(:,322) = 3.80e-12_r8 * exp_fac(:)
      rate(:,343) = 3.80e-12_r8 * exp_fac(:)
      rate(:,364) = 3.80e-12_r8 * exp_fac(:)
      rate(:,368) = 3.80e-12_r8 * exp_fac(:)
      rate(:,373) = 3.80e-12_r8 * exp_fac(:)
      rate(:,388) = 3.80e-12_r8 * exp_fac(:)
      rate(:,397) = 2.30e-11_r8 * exp_fac(:)
      rate(:,421) = 1.52e-11_r8 * exp_fac(:)
      rate(:,450) = 3.80e-12_r8 * exp_fac(:)
      rate(:,456) = 3.8e-12_r8 * exp_fac(:)
      rate(:,461) = 3.8e-12_r8 * exp_fac(:)
      rate(:,464) = 3.8e-12_r8 * exp_fac(:)
      rate(:,471) = 3.8e-12_r8 * exp_fac(:)
      rate(:,479) = 3.80e-12_r8 * exp_fac(:)
      rate(:,491) = 3.8e-12_r8 * exp_fac(:)
      rate(:,493) = 3.8e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -490._r8 * itemp(:) )
      rate(:,216) = 1.00e-14_r8 * exp_fac(:)
      rate(:,223) = 1.00e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( -2000._r8 * itemp(:) )
      rate(:,218) = 1.40e-12_r8 * exp_fac(:)
      rate(:,411) = 1.05e-14_r8 * exp_fac(:)
      rate(:,451) = 1.05e-14_r8 * exp_fac(:)
      rate(:,220) = 7.80e-13_r8 * exp( -1050._r8 * itemp(:) )
      rate(:,227) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:,228) = 2.10e-11_r8 * exp( 100._r8 * itemp(:) )
      exp_fac(:) = exp( 220._r8 * itemp(:) )
      rate(:,229) = 2.90e-12_r8 * exp_fac(:)
      rate(:,230) = 1.45e-12_r8 * exp_fac(:)
      rate(:,231) = 1.45e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 270._r8 * itemp(:) )
      rate(:,233) = 3.30e-12_r8 * exp_fac(:)
      rate(:,254) = 1.40e-11_r8 * exp_fac(:)
      rate(:,259) = 7.40e-12_r8 * exp_fac(:)
      rate(:,346) = 8.10e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( -1500._r8 * itemp(:) )
      rate(:,234) = 3.00e-12_r8 * exp_fac(:)
      rate(:,250) = 3.00e-12_r8 * exp_fac(:)
      rate(:,291) = 5.80e-12_r8 * exp_fac(:)
      rate(:,235) = 5.10e-12_r8 * exp( 210._r8 * itemp(:) )
      exp_fac(:) = exp( -2450._r8 * itemp(:) )
      rate(:,237) = 1.20e-13_r8 * exp_fac(:)
      rate(:,249) = 1.20e-13_r8 * exp_fac(:)
      rate(:,265) = 3.00e-11_r8 * exp_fac(:)
      rate(:,242) = 1.50e-11_r8 * exp( 170._r8 * itemp(:) )
      rate(:,247) = 1.30e-12_r8 * exp( 380._r8 * itemp(:) )
      rate(:,251) = 2.30e-11_r8 * exp( -200._r8 * itemp(:) )
      rate(:,252) = 3.05e-11_r8 * exp( -2270._r8 * itemp(:) )
      rate(:,253) = 1.10e-11_r8 * exp( -980._r8 * itemp(:) )
      rate(:,255) = 3.60e-11_r8 * exp( -375._r8 * itemp(:) )
      rate(:,256) = 8.10e-11_r8 * exp( -30._r8 * itemp(:) )
      rate(:,257) = 7.10e-12_r8 * exp( -1270._r8 * itemp(:) )
      rate(:,258) = 2.80e-11_r8 * exp( 85._r8 * itemp(:) )
      exp_fac(:) = exp( 230._r8 * itemp(:) )
      rate(:,260) = 6.00e-13_r8 * exp_fac(:)
      rate(:,281) = 1.90e-11_r8 * exp_fac(:)
      rate(:,289) = 1.50e-12_r8 * exp_fac(:)
      exp_fac(:) = exp( 290._r8 * itemp(:) )
      rate(:,261) = 2.60e-12_r8 * exp_fac(:)
      rate(:,263) = 6.40e-12_r8 * exp_fac(:)
      rate(:,288) = 4.10e-13_r8 * exp_fac(:)
      rate(:,466) = 7.5e-12_r8 * exp_fac(:)
      rate(:,476) = 7.5e-12_r8 * exp_fac(:)
      rate(:,482) = 7.5e-12_r8 * exp_fac(:)
      rate(:,484) = 7.5e-12_r8 * exp_fac(:)
      rate(:,262) = 3.3e-12_r8 * exp( -115._r8 * itemp(:) )
      rate(:,266) = 1.00e-12_r8 * exp( -1590._r8 * itemp(:) )
      rate(:,267) = 3.50e-13_r8 * exp( -1370._r8 * itemp(:) )
      rate(:,270) = 1.80e-12_r8 * exp( -250._r8 * itemp(:) )
      rate(:,271) = 1.00e-11_r8 * exp( -3300._r8 * itemp(:) )
      rate(:,273) = 3.40e-12_r8 * exp( -130._r8 * itemp(:) )
      exp_fac(:) = exp( -500._r8 * itemp(:) )
      rate(:,274) = 3.00e-12_r8 * exp_fac(:)
      rate(:,295) = 1.40e-10_r8 * exp_fac(:)
      exp_fac(:) = exp( -840._r8 * itemp(:) )
      rate(:,275) = 3.60e-12_r8 * exp_fac(:)
      rate(:,306) = 2.00e-12_r8 * exp_fac(:)
      rate(:,276) = 1.20e-12_r8 * exp( -330._r8 * itemp(:) )
      rate(:,277) = 6.50e-12_r8 * exp( 135._r8 * itemp(:) )
      rate(:,278) = 1.60e-11_r8 * exp( -780._r8 * itemp(:) )
      rate(:,279) = 4.80e-12_r8 * exp( -310._r8 * itemp(:) )
      exp_fac(:) = exp( -800._r8 * itemp(:) )
      rate(:,280) = 1.70e-11_r8 * exp_fac(:)
      rate(:,308) = 6.30e-12_r8 * exp_fac(:)
      rate(:,283) = 4.50e-12_r8 * exp( 460._r8 * itemp(:) )
      exp_fac(:) = exp( 260._r8 * itemp(:) )
      rate(:,284) = 8.80e-12_r8 * exp_fac(:)
      rate(:,287) = 2.30e-12_r8 * exp_fac(:)
      rate(:,286) = 9.50e-13_r8 * exp( 550._r8 * itemp(:) )
      rate(:,292) = 1.20e-10_r8 * exp( -430._r8 * itemp(:) )
      rate(:,293) = 1.90e-11_r8 * exp( 215._r8 * itemp(:) )
      exp_fac(:) = exp( 0._r8 * itemp(:) )
      rate(:,294) = 1.40e-11_r8 * exp_fac(:)
      rate(:,323) = 4.00e-13_r8 * exp_fac(:)
      rate(:,335) = 1.00e-14_r8 * exp_fac(:)
      rate(:,338) = 7.00e-13_r8 * exp_fac(:)
      rate(:,341) = 2.00e-13_r8 * exp_fac(:)
      rate(:,342) = 6.80e-14_r8 * exp_fac(:)
      rate(:,351) = 1.00e-12_r8 * exp_fac(:)
      rate(:,352) = 1.00e-11_r8 * exp_fac(:)
      rate(:,353) = 1.15e-11_r8 * exp_fac(:)
      rate(:,356) = 4.00e-14_r8 * exp_fac(:)
      rate(:,374) = 3.00e-12_r8 * exp_fac(:)
      rate(:,377) = 6.7e-13_r8 * exp_fac(:)
      rate(:,379) = 5.40e-11_r8 * exp_fac(:)
      rate(:,382) = 3.5e-13_r8 * exp_fac(:)
      rate(:,393) = 2.40e-12_r8 * exp_fac(:)
      rate(:,396) = 1.40e-11_r8 * exp_fac(:)
      rate(:,399) = 5.00e-12_r8 * exp_fac(:)
      rate(:,407) = 2.0e-12_r8 * exp_fac(:)
      rate(:,415) = 4.e-11_r8 * exp_fac(:)
      rate(:,416) = 4.e-11_r8 * exp_fac(:)
      rate(:,417) = 2.40e-12_r8 * exp_fac(:)
      rate(:,418) = 2.40e-12_r8 * exp_fac(:)
      rate(:,422) = 1.3e-11_r8 * exp_fac(:)
      rate(:,425) = 1.40e-11_r8 * exp_fac(:)
      rate(:,426) = 1.40e-11_r8 * exp_fac(:)
      rate(:,430) = 2.40e-12_r8 * exp_fac(:)
      rate(:,432) = 1.4e-11_r8 * exp_fac(:)
      rate(:,434) = 4.0e-11_r8 * exp_fac(:)
      rate(:,435) = 7.0e-11_r8 * exp_fac(:)
      rate(:,436) = 1.0e-10_r8 * exp_fac(:)
      rate(:,439) = 2.40e-12_r8 * exp_fac(:)
      rate(:,445) = 3.50e-12_r8 * exp_fac(:)
      rate(:,446) = 6.7e-12_r8 * exp_fac(:)
      rate(:,448) = 1.6e-12_r8 * exp_fac(:)
      rate(:,457) = 2.1e-12_r8 * exp_fac(:)
      rate(:,458) = 2.8e-13_r8 * exp_fac(:)
      rate(:,469) = 4.7e-11_r8 * exp_fac(:)
      rate(:,487) = 1.7e-11_r8 * exp_fac(:)
      rate(:,488) = 8.4e-11_r8 * exp_fac(:)
      rate(:,495) = 2.8e-13_r8 * exp_fac(:)
      rate(:,497) = 2.0e-10_r8 * exp_fac(:)
      rate(:,499) = 1.2e-14_r8 * exp_fac(:)
      rate(:,501) = 1.9e-11_r8 * exp_fac(:)
      rate(:,505) = 3.3e-11_r8 * exp_fac(:)
      rate(:,506) = 2.3e-11_r8 * exp_fac(:)
      rate(:,507) = 5.7e-11_r8 * exp_fac(:)
      rate(:,508) = 1e-12_r8 * exp_fac(:)
      rate(:,512) = 3.4e-11_r8 * exp_fac(:)
      rate(:,516) = 2.4e-12_r8 * exp_fac(:)
      rate(:,517) = 2.0e-11_r8 * exp_fac(:)
      rate(:,518) = 2.0e-11_r8 * exp_fac(:)
      rate(:,520) = 1.2e-14_r8 * exp_fac(:)
      rate(:,521) = 1.34e-11_r8 * exp_fac(:)
      rate(:,522) = 1.34e-11_r8 * exp_fac(:)
      rate(:,531) = 6.60E-11_r8 * exp_fac(:)
      rate(:,532) = 2.30E-12_r8 * exp_fac(:)
      rate(:,533) = 1.20E-11_r8 * exp_fac(:)
      rate(:,537) = 1.40E-11_r8 * exp_fac(:)
      rate(:,538) = 2.80E-11_r8 * exp_fac(:)
      rate(:,539) = 5.70E-11_r8 * exp_fac(:)
      rate(:,540) = 1.90E-12_r8 * exp_fac(:)
      rate(:,546) = 6.34e-8_r8 * exp_fac(:)
      rate(:,547) = 6.34e-8_r8 * exp_fac(:)
      rate(:,569) = 9.0e-10_r8 * exp_fac(:)
      rate(:,570) = 1.0e-10_r8 * exp_fac(:)
      rate(:,571) = 4.4e-10_r8 * exp_fac(:)
      rate(:,572) = 4.0e-10_r8 * exp_fac(:)
      rate(:,573) = 2.0e-10_r8 * exp_fac(:)
      rate(:,574) = 1.0e-12_r8 * exp_fac(:)
      rate(:,575) = 6.0e-11_r8 * exp_fac(:)
      rate(:,576) = 5.0e-16_r8 * exp_fac(:)
      rate(:,296) = 1.60e-10_r8 * exp( -260._r8 * itemp(:) )
      exp_fac(:) = exp( 400._r8 * itemp(:) )
      rate(:,297) = 6.00e-12_r8 * exp_fac(:)
      rate(:,395) = 5.00e-13_r8 * exp_fac(:)
      rate(:,423) = 5.00e-13_r8 * exp_fac(:)
      rate(:,424) = 5.00e-13_r8 * exp_fac(:)
      rate(:,433) = 5.e-13_r8 * exp_fac(:)
      rate(:,441) = 5.e-13_r8 * exp_fac(:)
      exp_fac(:) = exp( -1100._r8 * itemp(:) )
      rate(:,298) = 2.03e-11_r8 * exp_fac(:)
      rate(:,536) = 3.40E-12_r8 * exp_fac(:)
      rate(:,299) = 1.96e-12_r8 * exp( -1200._r8 * itemp(:) )
      exp_fac(:) = exp( -1520._r8 * itemp(:) )
      rate(:,300) = 1.64e-12_r8 * exp_fac(:)
      rate(:,384) = 8.50e-16_r8 * exp_fac(:)
      rate(:,301) = 9.20e-13_r8 * exp( -1560._r8 * itemp(:) )
      rate(:,302) = 1.42e-12_r8 * exp( -1150._r8 * itemp(:) )
      rate(:,303) = 1.46e-11_r8 * exp( -1040._r8 * itemp(:) )
      exp_fac(:) = exp( -1600._r8 * itemp(:) )
      rate(:,304) = 1.25e-12_r8 * exp_fac(:)
      rate(:,315) = 3.40e-11_r8 * exp_fac(:)
      rate(:,305) = 1.30e-12_r8 * exp( -1770._r8 * itemp(:) )
      rate(:,307) = 9.00e-13_r8 * exp( -360._r8 * itemp(:) )
      rate(:,309) = 4.85e-12_r8 * exp( -850._r8 * itemp(:) )
      rate(:,310) = 2.45e-12_r8 * exp( -1775._r8 * itemp(:) )
      rate(:,313) = 6.00e-13_r8 * exp( -2058._r8 * itemp(:) )
      rate(:,314) = 5.50e-12_r8 * exp( 125._r8 * itemp(:) )
      rate(:,316) = 9.7e-15_r8 * exp( 625._r8 * itemp(:) )
      exp_fac(:) = exp( 300._r8 * itemp(:) )
      rate(:,317) = 2.80e-12_r8 * exp_fac(:)
      rate(:,370) = 2.90e-12_r8 * exp_fac(:)
      rate(:,318) = 4.10e-13_r8 * exp( 750._r8 * itemp(:) )
      rate(:,319) = 5.00e-13_r8 * exp( -424._r8 * itemp(:) )
      rate(:,320) = 1.90e-14_r8 * exp( 706._r8 * itemp(:) )
      rate(:,321) = 2.90e-12_r8 * exp( -345._r8 * itemp(:) )
      rate(:,324) = 2.40e12_r8 * exp( -7000._r8 * itemp(:) )
      rate(:,325) = 2.60e-12_r8 * exp( 265._r8 * itemp(:) )
      exp_fac(:) = exp( 700._r8 * itemp(:) )
      rate(:,326) = 7.50e-13_r8 * exp_fac(:)
      rate(:,334) = 7.50e-13_r8 * exp_fac(:)
      rate(:,340) = 7.50e-13_r8 * exp_fac(:)
      rate(:,362) = 7.50e-13_r8 * exp_fac(:)
      rate(:,367) = 7.50e-13_r8 * exp_fac(:)
      rate(:,371) = 8.60e-13_r8 * exp_fac(:)
      rate(:,387) = 7.50e-13_r8 * exp_fac(:)
      rate(:,394) = 8.00e-13_r8 * exp_fac(:)
      rate(:,419) = 8.00e-13_r8 * exp_fac(:)
      rate(:,420) = 8.00e-13_r8 * exp_fac(:)
      rate(:,431) = 8.00e-13_r8 * exp_fac(:)
      rate(:,440) = 8.00e-13_r8 * exp_fac(:)
      rate(:,449) = 7.50e-13_r8 * exp_fac(:)
      rate(:,455) = 7.5e-13_r8 * exp_fac(:)
      rate(:,460) = 7.5e-13_r8 * exp_fac(:)
      rate(:,463) = 7.5e-13_r8 * exp_fac(:)
      rate(:,470) = 7.5e-13_r8 * exp_fac(:)
      rate(:,478) = 7.50e-13_r8 * exp_fac(:)
      rate(:,490) = 7.5e-13_r8 * exp_fac(:)
      rate(:,492) = 7.5e-13_r8 * exp_fac(:)
      rate(:,503) = 7.5e-13_r8 * exp_fac(:)
      rate(:,510) = 7.5e-13_r8 * exp_fac(:)
      rate(:,514) = 7.5e-13_r8 * exp_fac(:)
      rate(:,329) = 7.20e-11_r8 * exp( -70._r8 * itemp(:) )
      rate(:,331) = 7.66e-12_r8 * exp( -1020._r8 * itemp(:) )
      rate(:,336) = 1.60e11_r8 * exp( -4150._r8 * itemp(:) )
      exp_fac(:) = exp( -2630._r8 * itemp(:) )
      rate(:,337) = 1.2e-14_r8 * exp_fac(:)
      rate(:,357) = 1.2e-14_r8 * exp_fac(:)
      exp_fac(:) = exp( 365._r8 * itemp(:) )
      rate(:,339) = 2.60e-12_r8 * exp_fac(:)
      rate(:,454) = 2.6e-12_r8 * exp_fac(:)
      rate(:,459) = 2.6e-12_r8 * exp_fac(:)
      rate(:,462) = 2.6e-12_r8 * exp_fac(:)
      rate(:,472) = 2.6e-12_r8 * exp_fac(:)
      rate(:,480) = 2.6e-12_r8 * exp_fac(:)
      rate(:,489) = 2.6e-12_r8 * exp_fac(:)
      rate(:,494) = 2.6e-12_r8 * exp_fac(:)
      rate(:,344) = 4.63e-12_r8 * exp( 350._r8 * itemp(:) )
      exp_fac(:) = exp( -1900._r8 * itemp(:) )
      rate(:,345) = 1.40e-12_r8 * exp_fac(:)
      rate(:,359) = 6.50e-15_r8 * exp_fac(:)
      rate(:,378) = 6.50e-15_r8 * exp_fac(:)
      exp_fac(:) = exp( 1040._r8 * itemp(:) )
      rate(:,348) = 4.30e-13_r8 * exp_fac(:)
      rate(:,400) = 4.30e-13_r8 * exp_fac(:)
      rate(:,467) = 4.3E-13_r8 * exp_fac(:)
      rate(:,477) = 4.3E-13_r8 * exp_fac(:)
      rate(:,481) = 4.3E-13_r8 * exp_fac(:)
      rate(:,483) = 4.3E-13_r8 * exp_fac(:)
      exp_fac(:) = exp( 500._r8 * itemp(:) )
      rate(:,349) = 2.00e-12_r8 * exp_fac(:)
      rate(:,350) = 2.90e-12_r8 * exp_fac(:)
      rate(:,372) = 7.10e-13_r8 * exp_fac(:)
      rate(:,401) = 2.00e-12_r8 * exp_fac(:)
      rate(:,504) = 2e-12_r8 * exp_fac(:)
      rate(:,511) = 2e-12_r8 * exp_fac(:)
      rate(:,515) = 2e-12_r8 * exp_fac(:)
      rate(:,354) = 6.90e-12_r8 * exp( -230._r8 * itemp(:) )
      rate(:,360) = 4.60e-13_r8 * exp( -1156._r8 * itemp(:) )
      rate(:,363) = 3.75e-13_r8 * exp( -40._r8 * itemp(:) )
      rate(:,365) = 8.70e-12_r8 * exp( -615._r8 * itemp(:) )
      rate(:,375) = 8.40e-13_r8 * exp( 830._r8 * itemp(:) )
      rate(:,376) = 1.40e-12_r8 * exp( -1860._r8 * itemp(:) )
      rate(:,380) = 4.8e-12_r8 * exp( 120._r8 * itemp(:) )
      rate(:,381) = 5.1e-14_r8 * exp( 693._r8 * itemp(:) )
      rate(:,383) = 4.13e-12_r8 * exp( 452._r8 * itemp(:) )
      rate(:,385) = 2.30e-12_r8 * exp( -170._r8 * itemp(:) )
      exp_fac(:) = exp( 360._r8 * itemp(:) )
      rate(:,389) = 9.60e-12_r8 * exp_fac(:)
      rate(:,391) = 2.70e-12_r8 * exp_fac(:)
      rate(:,392) = 1.30e-13_r8 * exp_fac(:)
      rate(:,398) = 5.30e-12_r8 * exp_fac(:)
      rate(:,429) = 2.70e-12_r8 * exp_fac(:)
      rate(:,438) = 2.7e-12_r8 * exp_fac(:)
      rate(:,390) = 1.50e-15_r8 * exp( -2100._r8 * itemp(:) )
      exp_fac(:) = exp( 530._r8 * itemp(:) )
      rate(:,402) = 4.60e-12_r8 * exp_fac(:)
      rate(:,403) = 2.30e-12_r8 * exp_fac(:)
      rate(:,408) = 4.40e-15_r8 * exp( -2500._r8 * itemp(:) )
      rate(:,409) = 7.52e-16_r8 * exp( -1521._r8 * itemp(:) )
      rate(:,410) = 2.54e-11_r8 * exp( 410._r8 * itemp(:) )
      rate(:,412) = 3.03e-12_r8 * exp( -446._r8 * itemp(:) )
      rate(:,427) = 1.6e9_r8 * exp( -8300._r8 * itemp(:) )
      exp_fac(:) = exp( 175._r8 * itemp(:) )
      rate(:,428) = 1.86e-11_r8 * exp_fac(:)
      rate(:,437) = 1.86e-11_r8 * exp_fac(:)
      rate(:,442) = 1.3e-12_r8 * exp( 640._r8 * itemp(:) )
      rate(:,443) = 1.90e-12_r8 * exp( 190._r8 * itemp(:) )
      rate(:,447) = 5.4e-14_r8 * exp( 870._r8 * itemp(:) )
      rate(:,452) = 2.3e-12_r8 * exp( -193._r8 * itemp(:) )
      rate(:,453) = 4.7e-13_r8 * exp( 1220._r8 * itemp(:) )
      rate(:,468) = 1.7e-12_r8 * exp( 352._r8 * itemp(:) )
      rate(:,473) = 5.9e-12_r8 * exp( 225._r8 * itemp(:) )
      rate(:,496) = 1.2e-11_r8 * exp( 440._r8 * itemp(:) )
      exp_fac(:) = exp( -580._r8 * itemp(:) )
      rate(:,498) = 6.3e-16_r8 * exp_fac(:)
      rate(:,519) = 6.3e-16_r8 * exp_fac(:)
      rate(:,500) = 1.2e-12_r8 * exp( 490._r8 * itemp(:) )
      rate(:,527) = 1.70e-12_r8 * exp( -710._r8 * itemp(:) )
      rate(:,529) = 2.10E-11_r8 * exp( -2200.0_r8 * itemp(:) )
      rate(:,530) = 7.20E-14_r8 * exp( -1070.0_r8 * itemp(:) )
      rate(:,534) = 2.70E-11_r8 * exp( 335._r8 * itemp(:) )
      rate(:,535) = 1.60E-13_r8 * exp( -2280.0_r8 * itemp(:) )
      rate(:,543) = 9.60e-12_r8 * exp( -234._r8 * itemp(:) )
      rate(:,545) = 1.90e-13_r8 * exp( 520._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)
 
      n = ncol*pver

      ko(:) = 4.40e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( rate(:,202), m, 0.6_r8, ko, kinf, n )

      ko(:) = 6.90e-31_r8 * itemp(:)**1.0_r8
      kinf(:) = 2.60e-11_r8
      call jpl( rate(:,211), m, 0.6_r8, ko, kinf, n )

      ko(:) = 4.28e-33_r8
      kinf(:) = 9.30e-15_r8 * itemp(:)**(-4.42_r8)
      call jpl( rate(:,219), m, 0.8_r8, ko, kinf, n )

      ko(:) = 9.00e-32_r8 * itemp(:)**1.5_r8
      kinf(:) = 3.0e-11_r8
      call jpl( rate(:,232), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.50e-31_r8 * itemp(:)**1.8_r8
      kinf(:) = 2.2e-11_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,236), m, 0.6_r8, ko, kinf, n )

      ko(:) = 2.40e-30_r8 * itemp(:)**3.0_r8
      kinf(:) = 1.6e-12_r8 * itemp(:)**(-0.1_r8)
      call jpl( rate(:,238), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.80e-30_r8 * itemp(:)**3.0_r8
      kinf(:) = 2.8e-11_r8
      call jpl( rate(:,240), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.90e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 4.0e-12_r8 * itemp(:)**0.3_r8
      call jpl( rate(:,246), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.80e-31_r8 * itemp(:)**3.4_r8
      kinf(:) = 1.5e-11_r8 * itemp(:)**1.9_r8
      call jpl( rate(:,264), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.90e-32_r8 * itemp(:)**3.6_r8
      kinf(:) = 3.7e-12_r8 * itemp(:)**1.6_r8
      call jpl( rate(:,268), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.20e-31_r8 * itemp(:)**3.2_r8
      kinf(:) = 6.9e-12_r8 * itemp(:)**2.9_r8
      call jpl( rate(:,285), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.90e-33_r8 * itemp(:)**1.0_r8
      kinf(:) = 1.10e-12_r8 * itemp(:)**(-1.3_r8)
      call jpl( rate(:,312), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.20e-30_r8 * itemp(:)**2.4_r8
      kinf(:) = 2.2e-10_r8 * itemp(:)**0.7_r8
      call jpl( rate(:,327), m, 0.6_r8, ko, kinf, n )

      ko(:) = 1.60e-29_r8 * itemp(:)**3.3_r8
      kinf(:) = 3.1e-10_r8 * itemp(:)
      call jpl( rate(:,328), m, 0.6_r8, ko, kinf, n )

      ko(:) = 5.50e-30_r8
      kinf(:) = 8.3e-13_r8 * itemp(:)**(-2.0_r8)
      call jpl( rate(:,330), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.60e-29_r8 * itemp(:)**3.1_r8
      kinf(:) = 9.00e-12_r8 * itemp(:)**0.85_r8
      call jpl( rate(:,332), m, 0.48_r8, ko, kinf, n )

      ko(:) = 9.70e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.30e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,347), m, 0.6_r8, ko, kinf, n )

      ko(:) = 8.00e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3.00e-11_r8
      call jpl( rate(:,358), m, 0.5_r8, ko, kinf, n )

      ko(:) = 8.00e-27_r8 * itemp(:)**3.5_r8
      kinf(:) = 3.00e-11_r8
      call jpl( rate(:,406), m, 0.5_r8, ko, kinf, n )

      ko(:) = 9.70e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.30e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,465), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.70e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.30e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,474), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.70e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.30e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,485), m, 0.6_r8, ko, kinf, n )

      ko(:) = 9.70e-29_r8 * itemp(:)**5.6_r8
      kinf(:) = 9.30e-12_r8 * itemp(:)**1.5_r8
      call jpl( rate(:,486), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )
 
      use ppgrid, only : pcols, pver


      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol*pver)
      real(r8), intent(inout) :: rate(ncol*pver,max(1,rxntot))

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      integer   ::  offset
      integer   ::  k
      real(r8)  :: itemp(ncol*kbot)
      real(r8)  :: exp_fac(ncol*kbot)
      real(r8)  :: ko(ncol*kbot)
      real(r8)  :: kinf(ncol*kbot)
      real(r8)  :: wrk(ncol*kbot)
 
      n = ncol*kbot

      rate(:n,156) = 8.00e-14_r8
      rate(:n,157) = 3.90e-17_r8
      rate(:n,162) = 1.30e-16_r8
      rate(:n,164) = 1.00e-20_r8
      rate(:n,205) = 6.90e-12_r8
      rate(:n,224) = 5.00e-12_r8
      rate(:n,225) = 7.00e-13_r8
      rate(:n,570) = 1.0e-10_r8
      rate(:n,571) = 4.4e-10_r8
      rate(:n,572) = 4.0e-10_r8
      rate(:n,573) = 2.0e-10_r8
      rate(:n,574) = 1.0e-12_r8
      rate(:n,575) = 6.0e-11_r8
 
      do k = 1,kbot
        offset = (k-1)*ncol
        itemp(offset+1:offset+ncol) = 1._r8 / temp(:ncol,k)
      end do

      rate(:n,154) = 8.00e-12_r8 * exp( -2060._r8 * itemp(:) )
      rate(:n,158) = 1.80e-15_r8 * exp( 45._r8 * itemp(:) )
      rate(:n,159) = 3.50e-11_r8 * exp( -135._r8 * itemp(:) )
      rate(:n,163) = 3.60e-18_r8 * exp( -220._r8 * itemp(:) )
      rate(:n,167) = 2.15e-11_r8 * exp( 110._r8 * itemp(:) )
      exp_fac(:) = exp( 55._r8 * itemp(:) )
      rate(:n,168) = 2.64e-11_r8 * exp_fac(:)
      rate(:n,169) = 6.60e-12_r8 * exp_fac(:)
      rate(:n,203) = 1.40e-10_r8 * exp( -470._r8 * itemp(:) )
      rate(:n,207) = 1.80e-11_r8 * exp( 180._r8 * itemp(:) )
      rate(:n,208) = 1.70e-12_r8 * exp( -940._r8 * itemp(:) )
      rate(:n,209) = 4.80e-11_r8 * exp( 250._r8 * itemp(:) )
      rate(:n,215) = 3.00e-11_r8 * exp( 200._r8 * itemp(:) )
      rate(:n,216) = 1.00e-14_r8 * exp( -490._r8 * itemp(:) )
      rate(:n,227) = 1.50e-11_r8 * exp( -3600._r8 * itemp(:) )
      rate(:n,228) = 2.10e-11_r8 * exp( 100._r8 * itemp(:) )
      rate(:n,233) = 3.30e-12_r8 * exp( 270._r8 * itemp(:) )
      rate(:n,234) = 3.00e-12_r8 * exp( -1500._r8 * itemp(:) )
      rate(:n,235) = 5.10e-12_r8 * exp( 210._r8 * itemp(:) )

      itemp(:) = 300._r8 * itemp(:)

      ko(:) = 4.40e-32_r8 * itemp(:)**1.3_r8
      kinf(:) = 7.5e-11_r8 * itemp(:)**(-0.2_r8)
      call jpl( wrk, m, 0.6_r8, ko, kinf, n )
      rate(:n,202) = wrk(:)























      end subroutine setrxt_hrates

      end module mo_setrxt
