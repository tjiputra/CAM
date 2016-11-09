      module mo_indprd
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: indprd
      contains
      subroutine indprd( class, prod, nprod, y, extfrc, rxt, chnkpnts )
      use chem_mods, only : gas_pcnst, extcnt, rxntot
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: chnkpnts
      integer, intent(in) :: nprod
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(in) :: extfrc(chnkpnts,extcnt)
      real(r8), intent(inout) :: prod(chnkpnts,nprod)
!--------------------------------------------------------------------
! ... "independent" production for Explicit species
!--------------------------------------------------------------------
      if( class == 1 ) then
         prod(:,1) = 0._r8
         prod(:,2) =.100_r8*rxt(:,359)*y(:,47)*y(:,1)
         prod(:,3) =rxt(:,229)*y(:,7)*y(:,5)
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,8) = 0._r8
         prod(:,9) = 0._r8
         prod(:,10) = 0._r8
         prod(:,11) = 0._r8
         prod(:,12) = 0._r8
         prod(:,13) = 0._r8
         prod(:,14) = 0._r8
         prod(:,15) = 0._r8
         prod(:,16) = 0._r8
         prod(:,17) = 0._r8
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = 0._r8
         prod(:,21) = 0._r8
         prod(:,22) = (rxt(:,311)*y(:,16) +rxt(:,312)*y(:,16) +rxt(:,323)*y(:,132) + &
                 rxt(:,338)*y(:,39) +.500_r8*rxt(:,351)*y(:,44) + &
                 .800_r8*rxt(:,352)*y(:,42) +rxt(:,353)*y(:,43) + &
                 .500_r8*rxt(:,406)*y(:,62) +1.800_r8*rxt(:,512)*y(:,109))*y(:,185) &
                  + (rxt(:,346)*y(:,6) +.900_r8*rxt(:,349)*y(:,187) + &
                 2.000_r8*rxt(:,350)*y(:,190) +2.000_r8*rxt(:,402)*y(:,198) + &
                 rxt(:,425)*y(:,201) +rxt(:,442)*y(:,203))*y(:,190) &
                  + (.200_r8*rxt(:,359)*y(:,47) +.100_r8*rxt(:,384)*y(:,59) + &
                 .270_r8*rxt(:,498)*y(:,105) +.270_r8*rxt(:,499)*y(:,106))*y(:,1) &
                  + (.450_r8*rxt(:,400)*y(:,186) +rxt(:,401)*y(:,187) + &
                 2.000_r8*rxt(:,403)*y(:,198))*y(:,198) + (.900_r8*rxt(:,509)*y(:,6) + &
                 .500_r8*rxt(:,511)*y(:,187))*y(:,216) +rxt(:,65)*y(:,44) &
                  +.400_r8*rxt(:,66)*y(:,46) +.800_r8*rxt(:,98)*y(:,109) +rxt(:,96) &
                 *y(:,110)
         prod(:,23) = 0._r8
         prod(:,24) = 0._r8
         prod(:,25) = (rxt(:,457)*y(:,89) +rxt(:,465)*y(:,211) +rxt(:,485)*y(:,214) + &
                 rxt(:,486)*y(:,215))*y(:,7) +.200_r8*rxt(:,513)*y(:,217)*y(:,6) &
                  +.500_r8*rxt(:,508)*y(:,108)*y(:,8) +.500_r8*rxt(:,406)*y(:,185) &
                 *y(:,62) +rxt(:,547)*y(:,141)
         prod(:,26) =rxt(:,527)*y(:,185)*y(:,139) +rxt(:,546)*y(:,140)
!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      else if( class == 4 ) then
         prod(:,193) = 0._r8
         prod(:,196) = (rxt(:,60) +rxt(:,152))*y(:,130) +.180_r8*rxt(:,62)*y(:,12)
         prod(:,199) =rxt(:,7)*y(:,4)
         prod(:,172) =1.440_r8*rxt(:,62)*y(:,12)
         prod(:,168) = (rxt(:,60) +rxt(:,152))*y(:,130) +.380_r8*rxt(:,62)*y(:,12) &
                  + extfrc(:,3)
         prod(:,194) =.330_r8*rxt(:,62)*y(:,12) + extfrc(:,22)
         prod(:,145) = (rxt(:,136) +.800_r8*rxt(:,139) +rxt(:,148) + &
                 .800_r8*rxt(:,151)) + extfrc(:,21)
         prod(:,201) = + extfrc(:,1)
         prod(:,197) = + extfrc(:,2)
         prod(:,190) = 0._r8
         prod(:,200) = 0._r8
         prod(:,101) = 0._r8
         prod(:,73) = 0._r8
         prod(:,186) =rxt(:,61)*y(:,12) +rxt(:,39)*y(:,111) +rxt(:,50)*y(:,112)
         prod(:,94) = 0._r8
         prod(:,58) = 0._r8
         prod(:,38) = 0._r8
         prod(:,192) =.180_r8*rxt(:,62)*y(:,12)
         prod(:,202) = (rxt(:,61) +.330_r8*rxt(:,62))*y(:,12)
         prod(:,191) = 0._r8
         prod(:,120) = 0._r8
         prod(:,198) =.050_r8*rxt(:,62)*y(:,12)
         prod(:,188) =rxt(:,39)*y(:,111) +2.000_r8*rxt(:,42)*y(:,113) &
                  +2.000_r8*rxt(:,43)*y(:,114) +2.000_r8*rxt(:,44)*y(:,115) +rxt(:,47) &
                 *y(:,116) +4.000_r8*rxt(:,40)*y(:,117) +3.000_r8*rxt(:,41)*y(:,118) &
                  +rxt(:,52)*y(:,120) +rxt(:,48)*y(:,121) +rxt(:,49)*y(:,122) &
                  +2.000_r8*rxt(:,45)*y(:,123) +rxt(:,46)*y(:,124)
         prod(:,52) = 0._r8
         prod(:,195) = 0._r8
         prod(:,82) = 0._r8
         prod(:,39) = 0._r8
         prod(:,184) = 0._r8
         prod(:,146) = 0._r8
         prod(:,156) = 0._r8
         prod(:,64) = 0._r8
         prod(:,185) =rxt(:,50)*y(:,112) +rxt(:,51)*y(:,119) +rxt(:,52)*y(:,120) &
                  +2.000_r8*rxt(:,55)*y(:,125) +2.000_r8*rxt(:,56)*y(:,126) &
                  +3.000_r8*rxt(:,53)*y(:,127) +2.000_r8*rxt(:,54)*y(:,128)
         prod(:,189) = 0._r8
         prod(:,148) = 0._r8
         prod(:,137) = 0._r8
         prod(:,117) = 0._r8
         prod(:,46) =rxt(:,43)*y(:,114) +rxt(:,44)*y(:,115) +rxt(:,47)*y(:,116) &
                  +rxt(:,51)*y(:,119) +rxt(:,52)*y(:,120) +rxt(:,49)*y(:,122) &
                  +2.000_r8*rxt(:,45)*y(:,123) +2.000_r8*rxt(:,46)*y(:,124) +rxt(:,55) &
                 *y(:,125) +2.000_r8*rxt(:,56)*y(:,126)
         prod(:,59) =rxt(:,42)*y(:,113) +rxt(:,44)*y(:,115) +rxt(:,48)*y(:,121)
         prod(:,63) = 0._r8
         prod(:,136) =rxt(:,51)*y(:,119) +rxt(:,46)*y(:,124)
         prod(:,161) = 0._r8
         prod(:,149) = 0._r8
         prod(:,143) = 0._r8
         prod(:,171) = 0._r8
         prod(:,114) = 0._r8
         prod(:,115) = 0._r8
         prod(:,182) = 0._r8
         prod(:,109) = 0._r8
         prod(:,99) = 0._r8
         prod(:,69) = 0._r8
         prod(:,112) = 0._r8
         prod(:,62) = 0._r8
         prod(:,113) = 0._r8
         prod(:,124) = 0._r8
         prod(:,150) = 0._r8
         prod(:,167) = 0._r8
         prod(:,116) = 0._r8
         prod(:,102) = 0._r8
         prod(:,95) = 0._r8
         prod(:,147) = 0._r8
         prod(:,78) = 0._r8
         prod(:,108) = 0._r8
         prod(:,83) = 0._r8
         prod(:,85) = 0._r8
         prod(:,111) = 0._r8
         prod(:,160) = 0._r8
         prod(:,123) = 0._r8
         prod(:,107) = 0._r8
         prod(:,125) = 0._r8
         prod(:,75) = 0._r8
         prod(:,53) = 0._r8
         prod(:,54) = 0._r8
         prod(:,131) = 0._r8
         prod(:,127) = 0._r8
         prod(:,152) = 0._r8
         prod(:,100) = 0._r8
         prod(:,88) = 0._r8
         prod(:,151) = 0._r8
         prod(:,49) = 0._r8
         prod(:,55) = 0._r8
         prod(:,51) = 0._r8
         prod(:,50) = 0._r8
         prod(:,103) = 0._r8
         prod(:,93) = 0._r8
         prod(:,104) = 0._r8
         prod(:,74) = 0._r8
         prod(:,122) = 0._r8
         prod(:,65) = 0._r8
         prod(:,86) = 0._r8
         prod(:,96) = 0._r8
         prod(:,70) = 0._r8
         prod(:,118) = 0._r8
         prod(:,76) = 0._r8
         prod(:,132) = 0._r8
         prod(:,56) = 0._r8
         prod(:,87) = 0._r8
         prod(:,77) = 0._r8
         prod(:,60) = 0._r8
         prod(:,106) = 0._r8
         prod(:,134) = 0._r8
         prod(:,164) = 0._r8
         prod(:,41) = 0._r8
         prod(:,57) = 0._r8
         prod(:,98) = 0._r8
         prod(:,89) = 0._r8
         prod(:,133) = 0._r8
         prod(:,130) = 0._r8
         prod(:,155) = 0._r8
         prod(:,157) = 0._r8
         prod(:,8) = 0._r8
         prod(:,9) = + extfrc(:,5)
         prod(:,158) = 0._r8
         prod(:,166) = 0._r8
         prod(:,163) = 0._r8
         prod(:,119) = 0._r8
         prod(:,162) = 0._r8
         prod(:,177) = 0._r8
         prod(:,178) = 0._r8
         prod(:,61) = 0._r8
         prod(:,42) = 0._r8
         prod(:,180) = 0._r8
         prod(:,174) = 0._r8
         prod(:,179) = 0._r8
         prod(:,79) = 0._r8
         prod(:,181) = 0._r8
         prod(:,153) = 0._r8
         prod(:,80) = 0._r8
         prod(:,40) = 0._r8
         prod(:,144) = 0._r8
         prod(:,90) = 0._r8
         prod(:,159) = 0._r8
         prod(:,91) = 0._r8
         prod(:,142) = 0._r8
         prod(:,66) = 0._r8
         prod(:,165) = 0._r8
         prod(:,169) = 0._r8
         prod(:,135) = 0._r8
         prod(:,97) = 0._r8
         prod(:,43) = 0._r8
         prod(:,81) = 0._r8
         prod(:,170) = 0._r8
         prod(:,175) = 0._r8
         prod(:,173) = 0._r8
         prod(:,44) = 0._r8
         prod(:,176) = 0._r8
         prod(:,67) = 0._r8
         prod(:,128) = 0._r8
         prod(:,71) = 0._r8
         prod(:,141) = 0._r8
         prod(:,92) = 0._r8
         prod(:,154) = + extfrc(:,4)
         prod(:,72) = 0._r8
         prod(:,84) = 0._r8
         prod(:,129) = 0._r8
         prod(:,183) = 0._r8
         prod(:,68) = 0._r8
         prod(:,45) = 0._r8
         prod(:,37) = 0._r8
         prod(:,1) = 0._r8
         prod(:,2) = 0._r8
         prod(:,3) = 0._r8
         prod(:,4) = 0._r8
         prod(:,5) = 0._r8
         prod(:,6) = 0._r8
         prod(:,7) = 0._r8
         prod(:,10) = 0._r8
         prod(:,11) = 0._r8
         prod(:,12) = 0._r8
         prod(:,13) = 0._r8
         prod(:,14) = 0._r8
         prod(:,15) = 0._r8
         prod(:,16) = 0._r8
         prod(:,17) = 0._r8
         prod(:,18) = 0._r8
         prod(:,19) = 0._r8
         prod(:,20) = + extfrc(:,6)
         prod(:,21) = + extfrc(:,8)
         prod(:,22) = + extfrc(:,10)
         prod(:,23) = 0._r8
         prod(:,24) = 0._r8
         prod(:,25) = + extfrc(:,12)
         prod(:,26) = + extfrc(:,7)
         prod(:,27) = 0._r8
         prod(:,28) = 0._r8
         prod(:,29) = + extfrc(:,13)
         prod(:,30) = 0._r8
         prod(:,31) = 0._r8
         prod(:,32) = 0._r8
         prod(:,33) = 0._r8
         prod(:,34) = + extfrc(:,9)
         prod(:,35) = + extfrc(:,11)
         prod(:,36) = + extfrc(:,14)
         prod(:,187) = 0._r8
         prod(:,48) = 0._r8
         prod(:,47) = 0._r8
         prod(:,126) = (rxt(:,132) +rxt(:,144)) + extfrc(:,18)
         prod(:,138) = + extfrc(:,16)
         prod(:,105) = (rxt(:,136) +rxt(:,137) +rxt(:,148) +rxt(:,149)) + extfrc(:,17)
         prod(:,121) = + extfrc(:,15)
         prod(:,139) = 0._r8
         prod(:,110) = (rxt(:,137) +1.200_r8*rxt(:,139) +rxt(:,149) + &
                 1.200_r8*rxt(:,151)) + extfrc(:,19)
         prod(:,140) = (rxt(:,132) +rxt(:,136) +rxt(:,137) +rxt(:,144) +rxt(:,148) + &
                 rxt(:,149)) + extfrc(:,20)
      end if
      end subroutine indprd
      end module mo_indprd
