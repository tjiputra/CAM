      module mo_nln_matrix
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: nlnmat
      contains
      subroutine nlnmat01( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,1711) = -(rxt(k,154)*y(k,2) + rxt(k,173)*y(k,227) + rxt(k,203)*y(k,20) &
                      + rxt(k,208)*y(k,185) + rxt(k,216)*y(k,186) + rxt(k,234)*y(k,6) &
                      + rxt(k,237)*y(k,7) + rxt(k,251)*y(k,183) + rxt(k,278)*y(k,184) &
                      + rxt(k,337)*y(k,36) + rxt(k,359)*y(k,47) + rxt(k,384)*y(k,59) &
                      + rxt(k,390)*y(k,60) + rxt(k,411)*y(k,73) + rxt(k,458)*y(k,89) &
                      + rxt(k,498)*y(k,105) + rxt(k,499)*y(k,106) + rxt(k,533) &
                      *y(k,143) + rxt(k,536)*y(k,144))
         mat(k,1927) = -rxt(k,154)*y(k,1)
         mat(k,2021) = -rxt(k,173)*y(k,1)
         mat(k,2144) = -rxt(k,203)*y(k,1)
         mat(k,1859) = -rxt(k,208)*y(k,1)
         mat(k,1627) = -rxt(k,216)*y(k,1)
         mat(k,2124) = -rxt(k,234)*y(k,1)
         mat(k,1969) = -rxt(k,237)*y(k,1)
         mat(k,1452) = -rxt(k,251)*y(k,1)
         mat(k,1342) = -rxt(k,278)*y(k,1)
         mat(k,445) = -rxt(k,337)*y(k,1)
         mat(k,955) = -rxt(k,359)*y(k,1)
         mat(k,1242) = -rxt(k,384)*y(k,1)
         mat(k,1111) = -rxt(k,390)*y(k,1)
         mat(k,793) = -rxt(k,411)*y(k,1)
         mat(k,394) = -rxt(k,458)*y(k,1)
         mat(k,862) = -rxt(k,498)*y(k,1)
         mat(k,903) = -rxt(k,499)*y(k,1)
         mat(k,584) = -rxt(k,533)*y(k,1)
         mat(k,1309) = -rxt(k,536)*y(k,1)
         mat(k,1927) = mat(k,1927) + rxt(k,153)*y(k,3)
         mat(k,1627) = mat(k,1627) + .150_r8*rxt(k,348)*y(k,190) + .150_r8*rxt(k,400) &
                      *y(k,198)
         mat(k,1293) = .150_r8*rxt(k,348)*y(k,186)
         mat(k,1260) = .150_r8*rxt(k,400)*y(k,186)
         mat(k,1418) = rxt(k,153)*y(k,2)
         mat(k,1930) = -(rxt(k,153)*y(k,3) + rxt(k,154)*y(k,1) + 4._r8*rxt(k,155) &
                      *y(k,2) + rxt(k,207)*y(k,185) + rxt(k,214)*y(k,19) + rxt(k,215) &
                      *y(k,186) + rxt(k,218)*y(k,21) + rxt(k,232)*y(k,6) + (rxt(k,235) &
                      + rxt(k,236)) * y(k,7) + rxt(k,243)*y(k,8) + rxt(k,258)*y(k,25) &
                      + rxt(k,271)*y(k,28) + rxt(k,272)*y(k,29) + rxt(k,275)*y(k,30) &
                      + rxt(k,281)*y(k,32) + rxt(k,291)*y(k,33) + rxt(k,292)*y(k,34) &
                      + rxt(k,293)*y(k,35) + rxt(k,315)*y(k,15) + rxt(k,529)*y(k,142) &
                      + (rxt(k,567) + rxt(k,568)) * y(k,222) + rxt(k,574)*y(k,221))
         mat(k,1421) = -rxt(k,153)*y(k,2)
         mat(k,1714) = -rxt(k,154)*y(k,2)
         mat(k,1862) = -rxt(k,207)*y(k,2)
         mat(k,1073) = -rxt(k,214)*y(k,2)
         mat(k,1630) = -rxt(k,215)*y(k,2)
         mat(k,501) = -rxt(k,218)*y(k,2)
         mat(k,2127) = -rxt(k,232)*y(k,2)
         mat(k,1972) = -(rxt(k,235) + rxt(k,236)) * y(k,2)
         mat(k,1536) = -rxt(k,243)*y(k,2)
         mat(k,1889) = -rxt(k,258)*y(k,2)
         mat(k,1329) = -rxt(k,271)*y(k,2)
         mat(k,755) = -rxt(k,272)*y(k,2)
         mat(k,875) = -rxt(k,275)*y(k,2)
         mat(k,1479) = -rxt(k,281)*y(k,2)
         mat(k,770) = -rxt(k,291)*y(k,2)
         mat(k,675) = -rxt(k,292)*y(k,2)
         mat(k,480) = -rxt(k,293)*y(k,2)
         mat(k,1654) = -rxt(k,315)*y(k,2)
         mat(k,266) = -rxt(k,529)*y(k,2)
         mat(k,557) = -(rxt(k,567) + rxt(k,568)) * y(k,2)
         mat(k,402) = -rxt(k,574)*y(k,2)
         mat(k,2024) = (rxt(k,168)+rxt(k,169))*y(k,3)
         mat(k,1862) = mat(k,1862) + 2.000_r8*rxt(k,210)*y(k,185)
         mat(k,744) = rxt(k,228)*y(k,6) + rxt(k,229)*y(k,7) + rxt(k,227)*y(k,3) &
                      + rxt(k,570)*y(k,219)
         mat(k,2127) = mat(k,2127) + rxt(k,228)*y(k,5)
         mat(k,1972) = mat(k,1972) + rxt(k,229)*y(k,5)
         mat(k,2147) = rxt(k,206)*y(k,186)
         mat(k,1630) = mat(k,1630) + rxt(k,206)*y(k,20)
         mat(k,586) = rxt(k,532)*y(k,3)
         mat(k,1312) = rxt(k,535)*y(k,3)
         mat(k,1421) = mat(k,1421) + (rxt(k,168)+rxt(k,169))*y(k,227) + rxt(k,227) &
                      *y(k,5) + rxt(k,532)*y(k,143) + rxt(k,535)*y(k,144) + rxt(k,573) &
                      *y(k,221) + rxt(k,565)*y(k,218)
         mat(k,682) = rxt(k,570)*y(k,5) + 1.150_r8*rxt(k,578)*y(k,223)
         mat(k,402) = mat(k,402) + rxt(k,573)*y(k,3)
         mat(k,509) = rxt(k,565)*y(k,3)
         mat(k,690) = rxt(k,577)*y(k,223)
         mat(k,701) = 1.150_r8*rxt(k,578)*y(k,219) + rxt(k,577)*y(k,220)
         mat(k,2027) = -((rxt(k,168) + rxt(k,169)) * y(k,3) + rxt(k,170)*y(k,228) &
                      + rxt(k,173)*y(k,1) + rxt(k,190)*y(k,133) + rxt(k,191)*y(k,134) &
                      + rxt(k,195)*y(k,19) + (rxt(k,196) + rxt(k,197)) * y(k,28) &
                      + (rxt(k,198) + rxt(k,199)) * y(k,33) + rxt(k,200)*y(k,17))
         mat(k,1424) = -(rxt(k,168) + rxt(k,169)) * y(k,227)
         mat(k,2001) = -rxt(k,170)*y(k,227)
         mat(k,1717) = -rxt(k,173)*y(k,227)
         mat(k,76) = -rxt(k,190)*y(k,227)
         mat(k,149) = -rxt(k,191)*y(k,227)
         mat(k,1075) = -rxt(k,195)*y(k,227)
         mat(k,1332) = -(rxt(k,196) + rxt(k,197)) * y(k,227)
         mat(k,772) = -(rxt(k,198) + rxt(k,199)) * y(k,227)
         mat(k,145) = -rxt(k,200)*y(k,227)
         mat(k,1424) = mat(k,1424) + rxt(k,224)*y(k,226)
         mat(k,683) = .850_r8*rxt(k,578)*y(k,223)
         mat(k,429) = rxt(k,224)*y(k,3)
         mat(k,702) = .850_r8*rxt(k,578)*y(k,219)
         mat(k,1067) = -(rxt(k,195)*y(k,227) + rxt(k,212)*y(k,185) + rxt(k,214)*y(k,2) &
                      + rxt(k,252)*y(k,183) + rxt(k,295)*y(k,136))
         mat(k,2011) = -rxt(k,195)*y(k,19)
         mat(k,1838) = -rxt(k,212)*y(k,19)
         mat(k,1916) = -rxt(k,214)*y(k,19)
         mat(k,1441) = -rxt(k,252)*y(k,19)
         mat(k,661) = -rxt(k,295)*y(k,19)
         mat(k,2134) = rxt(k,205)*y(k,186)
         mat(k,1607) = rxt(k,205)*y(k,20)
         mat(k,1032) = -((rxt(k,311) + rxt(k,312)) * y(k,185))
         mat(k,1834) = -(rxt(k,311) + rxt(k,312)) * y(k,16)
         mat(k,1688) = .560_r8*rxt(k,359)*y(k,47) + .620_r8*rxt(k,411)*y(k,73) &
                      + .630_r8*rxt(k,337)*y(k,36) + .230_r8*rxt(k,498)*y(k,105) &
                      + .230_r8*rxt(k,499)*y(k,106) + .560_r8*rxt(k,384)*y(k,59) &
                      + .650_r8*rxt(k,390)*y(k,60)
         mat(k,1915) = rxt(k,315)*y(k,15) + rxt(k,529)*y(k,142)
         mat(k,1834) = mat(k,1834) + rxt(k,314)*y(k,15) + rxt(k,353)*y(k,43) &
                      + .700_r8*rxt(k,512)*y(k,109) + rxt(k,375)*y(k,53) &
                      + .350_r8*rxt(k,330)*y(k,131) + rxt(k,530)*y(k,142)
         mat(k,2102) = .400_r8*rxt(k,466)*y(k,211) + .170_r8*rxt(k,482)*y(k,214) &
                      + .350_r8*rxt(k,484)*y(k,215) + .225_r8*rxt(k,509)*y(k,216) &
                      + .220_r8*rxt(k,391)*y(k,199) + .250_r8*rxt(k,438)*y(k,203)
         mat(k,1510) = rxt(k,313)*y(k,15) + .220_r8*rxt(k,393)*y(k,199) + rxt(k,376) &
                      *y(k,53) + .500_r8*rxt(k,439)*y(k,203)
         mat(k,1368) = .125_r8*rxt(k,511)*y(k,216) + .110_r8*rxt(k,395)*y(k,199) &
                      + .200_r8*rxt(k,441)*y(k,203)
         mat(k,1640) = rxt(k,315)*y(k,2) + rxt(k,314)*y(k,185) + rxt(k,313)*y(k,8) &
                      + rxt(k,256)*y(k,183) + rxt(k,280)*y(k,184)
         mat(k,1603) = .160_r8*rxt(k,467)*y(k,211) + .070_r8*rxt(k,481)*y(k,214) &
                      + .140_r8*rxt(k,483)*y(k,215)
         mat(k,1439) = rxt(k,256)*y(k,15)
         mat(k,1336) = rxt(k,280)*y(k,15)
         mat(k,945) = .560_r8*rxt(k,359)*y(k,1)
         mat(k,782) = .620_r8*rxt(k,411)*y(k,1)
         mat(k,1276) = .220_r8*rxt(k,396)*y(k,199) + .250_r8*rxt(k,442)*y(k,203)
         mat(k,442) = .630_r8*rxt(k,337)*y(k,1)
         mat(k,804) = rxt(k,353)*y(k,185)
         mat(k,620) = .400_r8*rxt(k,466)*y(k,6) + .160_r8*rxt(k,467)*y(k,186)
         mat(k,642) = .170_r8*rxt(k,482)*y(k,6) + .070_r8*rxt(k,481)*y(k,186)
         mat(k,987) = .350_r8*rxt(k,484)*y(k,6) + .140_r8*rxt(k,483)*y(k,186)
         mat(k,854) = .230_r8*rxt(k,498)*y(k,1)
         mat(k,895) = .230_r8*rxt(k,499)*y(k,1)
         mat(k,1007) = .225_r8*rxt(k,509)*y(k,6) + .125_r8*rxt(k,511)*y(k,187)
         mat(k,977) = .700_r8*rxt(k,512)*y(k,185)
         mat(k,1229) = .560_r8*rxt(k,384)*y(k,1)
         mat(k,1103) = .650_r8*rxt(k,390)*y(k,1)
         mat(k,1209) = .220_r8*rxt(k,391)*y(k,6) + .220_r8*rxt(k,393)*y(k,8) &
                      + .110_r8*rxt(k,395)*y(k,187) + .220_r8*rxt(k,396)*y(k,190)
         mat(k,1115) = rxt(k,375)*y(k,185) + rxt(k,376)*y(k,8)
         mat(k,1128) = .250_r8*rxt(k,438)*y(k,6) + .500_r8*rxt(k,439)*y(k,8) &
                      + .200_r8*rxt(k,441)*y(k,187) + .250_r8*rxt(k,442)*y(k,190)
         mat(k,196) = .350_r8*rxt(k,330)*y(k,185)
         mat(k,263) = rxt(k,529)*y(k,2) + rxt(k,530)*y(k,185)
         mat(k,1860) = -(rxt(k,207)*y(k,2) + rxt(k,208)*y(k,1) + rxt(k,209)*y(k,186) &
                      + (4._r8*rxt(k,210) + 4._r8*rxt(k,211)) * y(k,185) + rxt(k,212) &
                      *y(k,19) + rxt(k,213)*y(k,21) + rxt(k,219)*y(k,17) + rxt(k,220) &
                      *y(k,18) + rxt(k,226)*y(k,5) + rxt(k,240)*y(k,7) + rxt(k,241) &
                      *y(k,9) + rxt(k,244)*y(k,8) + rxt(k,247)*y(k,10) + (rxt(k,259) &
                      + rxt(k,260)) * y(k,25) + rxt(k,270)*y(k,28) + rxt(k,274) &
                      *y(k,29) + rxt(k,276)*y(k,30) + rxt(k,282)*y(k,32) + rxt(k,290) &
                      *y(k,33) + (rxt(k,311) + rxt(k,312)) * y(k,16) + rxt(k,314) &
                      *y(k,15) + rxt(k,321)*y(k,14) + rxt(k,322)*y(k,13) + rxt(k,323) &
                      *y(k,132) + rxt(k,330)*y(k,131) + rxt(k,331)*y(k,37) + rxt(k,332) &
                      *y(k,36) + rxt(k,338)*y(k,39) + rxt(k,343)*y(k,38) + rxt(k,344) &
                      *y(k,40) + rxt(k,351)*y(k,44) + rxt(k,352)*y(k,42) + rxt(k,353) &
                      *y(k,43) + rxt(k,354)*y(k,41) + rxt(k,356)*y(k,46) + rxt(k,358) &
                      *y(k,47) + rxt(k,364)*y(k,49) + rxt(k,365)*y(k,48) + rxt(k,368) &
                      *y(k,51) + rxt(k,369)*y(k,50) + rxt(k,373)*y(k,54) + rxt(k,374) &
                      *y(k,52) + rxt(k,375)*y(k,53) + rxt(k,377)*y(k,64) + rxt(k,379) &
                      *y(k,55) + rxt(k,383)*y(k,59) + rxt(k,385)*y(k,57) + rxt(k,388) &
                      *y(k,58) + rxt(k,389)*y(k,60) + rxt(k,397)*y(k,61) + rxt(k,406) &
                      *y(k,62) + rxt(k,407)*y(k,65) + rxt(k,410)*y(k,73) + rxt(k,415) &
                      *y(k,66) + rxt(k,416)*y(k,67) + rxt(k,421)*y(k,82) + rxt(k,422) &
                      *y(k,78) + rxt(k,428)*y(k,77) + rxt(k,434)*y(k,68) + rxt(k,435) &
                      *y(k,70) + rxt(k,436)*y(k,69) + rxt(k,437)*y(k,76) + rxt(k,443) &
                      *y(k,81) + rxt(k,445)*y(k,56) + rxt(k,448)*y(k,63) + rxt(k,450) &
                      *y(k,74) + rxt(k,452)*y(k,86) + rxt(k,453)*y(k,87) + rxt(k,456) &
                      *y(k,90) + rxt(k,461)*y(k,91) + rxt(k,464)*y(k,92) + rxt(k,468) &
                      *y(k,83) + rxt(k,469)*y(k,84) + rxt(k,471)*y(k,98) + rxt(k,473) &
                      *y(k,99) + rxt(k,479)*y(k,85) + rxt(k,487)*y(k,101) + rxt(k,488) &
                      *y(k,102) + rxt(k,491)*y(k,103) + rxt(k,493)*y(k,104) + rxt(k,496) &
                      *y(k,105) + rxt(k,497)*y(k,106) + rxt(k,505)*y(k,107) + rxt(k,506) &
                      *y(k,110) + rxt(k,507)*y(k,108) + rxt(k,512)*y(k,109) + rxt(k,517) &
                      *y(k,71) + rxt(k,518)*y(k,72) + rxt(k,527)*y(k,139) + rxt(k,530) &
                      *y(k,142) + rxt(k,531)*y(k,143) + rxt(k,534)*y(k,144) + rxt(k,541) &
                      *y(k,137) + (rxt(k,543) + rxt(k,544)) * y(k,138))
         mat(k,1928) = -rxt(k,207)*y(k,185)
         mat(k,1712) = -rxt(k,208)*y(k,185)
         mat(k,1628) = -rxt(k,209)*y(k,185)
         mat(k,1072) = -rxt(k,212)*y(k,185)
         mat(k,500) = -rxt(k,213)*y(k,185)
         mat(k,144) = -rxt(k,219)*y(k,185)
         mat(k,52) = -rxt(k,220)*y(k,185)
         mat(k,743) = -rxt(k,226)*y(k,185)
         mat(k,1970) = -rxt(k,240)*y(k,185)
         mat(k,2044) = -rxt(k,241)*y(k,185)
         mat(k,1534) = -rxt(k,244)*y(k,185)
         mat(k,373) = -rxt(k,247)*y(k,185)
         mat(k,1887) = -(rxt(k,259) + rxt(k,260)) * y(k,185)
         mat(k,1327) = -rxt(k,270)*y(k,185)
         mat(k,753) = -rxt(k,274)*y(k,185)
         mat(k,873) = -rxt(k,276)*y(k,185)
         mat(k,1477) = -rxt(k,282)*y(k,185)
         mat(k,769) = -rxt(k,290)*y(k,185)
         mat(k,1034) = -(rxt(k,311) + rxt(k,312)) * y(k,185)
         mat(k,1652) = -rxt(k,314)*y(k,185)
         mat(k,711) = -rxt(k,321)*y(k,185)
         mat(k,330) = -rxt(k,322)*y(k,185)
         mat(k,706) = -rxt(k,323)*y(k,185)
         mat(k,198) = -rxt(k,330)*y(k,185)
         mat(k,188) = -rxt(k,331)*y(k,185)
         mat(k,446) = -rxt(k,332)*y(k,185)
         mat(k,457) = -rxt(k,338)*y(k,185)
         mat(k,246) = -rxt(k,343)*y(k,185)
         mat(k,1061) = -rxt(k,344)*y(k,185)
         mat(k,424) = -rxt(k,351)*y(k,185)
         mat(k,999) = -rxt(k,352)*y(k,185)
         mat(k,806) = -rxt(k,353)*y(k,185)
         mat(k,174) = -rxt(k,354)*y(k,185)
         mat(k,360) = -rxt(k,356)*y(k,185)
         mat(k,956) = -rxt(k,358)*y(k,185)
         mat(k,306) = -rxt(k,364)*y(k,185)
         mat(k,58) = -rxt(k,365)*y(k,185)
         mat(k,465) = -rxt(k,368)*y(k,185)
         mat(k,920) = -rxt(k,369)*y(k,185)
         mat(k,312) = -rxt(k,373)*y(k,185)
         mat(k,1041) = -rxt(k,374)*y(k,185)
         mat(k,1120) = -rxt(k,375)*y(k,185)
         mat(k,800) = -rxt(k,377)*y(k,185)
         mat(k,274) = -rxt(k,379)*y(k,185)
         mat(k,1243) = -rxt(k,383)*y(k,185)
         mat(k,414) = -rxt(k,385)*y(k,185)
         mat(k,220) = -rxt(k,388)*y(k,185)
         mat(k,1112) = -rxt(k,389)*y(k,185)
         mat(k,241) = -rxt(k,397)*y(k,185)
         mat(k,453) = -rxt(k,406)*y(k,185)
         mat(k,1029) = -rxt(k,407)*y(k,185)
         mat(k,794) = -rxt(k,410)*y(k,185)
         mat(k,473) = -rxt(k,415)*y(k,185)
         mat(k,382) = -rxt(k,416)*y(k,185)
         mat(k,580) = -rxt(k,421)*y(k,185)
         mat(k,65) = -rxt(k,422)*y(k,185)
         mat(k,157) = -rxt(k,428)*y(k,185)
         mat(k,337) = -rxt(k,434)*y(k,185)
         mat(k,236) = -rxt(k,435)*y(k,185)
         mat(k,764) = -rxt(k,436)*y(k,185)
         mat(k,249) = -rxt(k,437)*y(k,185)
         mat(k,177) = -rxt(k,443)*y(k,185)
         mat(k,160) = -rxt(k,445)*y(k,185)
         mat(k,538) = -rxt(k,448)*y(k,185)
         mat(k,528) = -rxt(k,450)*y(k,185)
         mat(k,92) = -rxt(k,452)*y(k,185)
         mat(k,101) = -rxt(k,453)*y(k,185)
         mat(k,215) = -rxt(k,456)*y(k,185)
         mat(k,170) = -rxt(k,461)*y(k,185)
         mat(k,281) = -rxt(k,464)*y(k,185)
         mat(k,115) = -rxt(k,468)*y(k,185)
         mat(k,120) = -rxt(k,469)*y(k,185)
         mat(k,231) = -rxt(k,471)*y(k,185)
         mat(k,152) = -rxt(k,473)*y(k,185)
         mat(k,572) = -rxt(k,479)*y(k,185)
         mat(k,132) = -rxt(k,487)*y(k,185)
         mat(k,141) = -rxt(k,488)*y(k,185)
         mat(k,301) = -rxt(k,491)*y(k,185)
         mat(k,602) = -rxt(k,493)*y(k,185)
         mat(k,863) = -rxt(k,496)*y(k,185)
         mat(k,904) = -rxt(k,497)*y(k,185)
         mat(k,368) = -rxt(k,505)*y(k,185)
         mat(k,495) = -rxt(k,506)*y(k,185)
         mat(k,915) = -rxt(k,507)*y(k,185)
         mat(k,983) = -rxt(k,512)*y(k,185)
         mat(k,418) = -rxt(k,517)*y(k,185)
         mat(k,258) = -rxt(k,518)*y(k,185)
         mat(k,48) = -rxt(k,527)*y(k,185)
         mat(k,265) = -rxt(k,530)*y(k,185)
         mat(k,585) = -rxt(k,531)*y(k,185)
         mat(k,1310) = -rxt(k,534)*y(k,185)
         mat(k,835) = -rxt(k,541)*y(k,185)
         mat(k,203) = -(rxt(k,543) + rxt(k,544)) * y(k,185)
         mat(k,1712) = mat(k,1712) + rxt(k,203)*y(k,20) + rxt(k,216)*y(k,186) &
                      + .360_r8*rxt(k,359)*y(k,47) + .320_r8*rxt(k,411)*y(k,73) &
                      + .130_r8*rxt(k,337)*y(k,36) + .630_r8*rxt(k,498)*y(k,105) &
                      + .630_r8*rxt(k,499)*y(k,106) + .360_r8*rxt(k,384)*y(k,59) &
                      + .240_r8*rxt(k,390)*y(k,60)
         mat(k,1928) = mat(k,1928) + rxt(k,214)*y(k,19) + rxt(k,315)*y(k,15) &
                      + rxt(k,215)*y(k,186) + rxt(k,218)*y(k,21) + rxt(k,271)*y(k,28) &
                      + rxt(k,272)*y(k,29) + rxt(k,291)*y(k,33) + rxt(k,292)*y(k,34)
         mat(k,2022) = rxt(k,195)*y(k,19) + rxt(k,200)*y(k,17) + 2.000_r8*rxt(k,170) &
                      *y(k,228) + rxt(k,196)*y(k,28) + rxt(k,198)*y(k,33)
         mat(k,1072) = mat(k,1072) + rxt(k,214)*y(k,2) + rxt(k,195)*y(k,227)
         mat(k,1860) = mat(k,1860) + .300_r8*rxt(k,322)*y(k,13) + .500_r8*rxt(k,368) &
                      *y(k,51) + .100_r8*rxt(k,397)*y(k,61) + .500_r8*rxt(k,343) &
                      *y(k,38) + .600_r8*rxt(k,421)*y(k,82) + .650_r8*rxt(k,330) &
                      *y(k,131)
         mat(k,2125) = rxt(k,233)*y(k,186)
         mat(k,1534) = mat(k,1534) + rxt(k,245)*y(k,186)
         mat(k,330) = mat(k,330) + .300_r8*rxt(k,322)*y(k,185)
         mat(k,144) = mat(k,144) + rxt(k,200)*y(k,227)
         mat(k,1652) = mat(k,1652) + rxt(k,315)*y(k,2)
         mat(k,2145) = rxt(k,203)*y(k,1) + 2.000_r8*rxt(k,204)*y(k,186)
         mat(k,1628) = mat(k,1628) + rxt(k,216)*y(k,1) + rxt(k,215)*y(k,2) &
                      + rxt(k,233)*y(k,6) + rxt(k,245)*y(k,8) + 2.000_r8*rxt(k,204) &
                      *y(k,20) + rxt(k,255)*y(k,183) + .450_r8*rxt(k,348)*y(k,190) &
                      + .200_r8*rxt(k,387)*y(k,197) + .400_r8*rxt(k,477)*y(k,213) &
                      + .400_r8*rxt(k,481)*y(k,214) + .400_r8*rxt(k,483)*y(k,215) &
                      + .450_r8*rxt(k,400)*y(k,198) + .150_r8*rxt(k,371)*y(k,195)
         mat(k,500) = mat(k,500) + rxt(k,218)*y(k,2)
         mat(k,1996) = 2.000_r8*rxt(k,170)*y(k,227) + rxt(k,294)*y(k,136)
         mat(k,1453) = rxt(k,255)*y(k,186)
         mat(k,1327) = mat(k,1327) + rxt(k,271)*y(k,2) + rxt(k,196)*y(k,227)
         mat(k,753) = mat(k,753) + rxt(k,272)*y(k,2)
         mat(k,769) = mat(k,769) + rxt(k,291)*y(k,2) + rxt(k,198)*y(k,227)
         mat(k,674) = rxt(k,292)*y(k,2)
         mat(k,664) = rxt(k,294)*y(k,228)
         mat(k,956) = mat(k,956) + .360_r8*rxt(k,359)*y(k,1)
         mat(k,794) = mat(k,794) + .320_r8*rxt(k,411)*y(k,1)
         mat(k,465) = mat(k,465) + .500_r8*rxt(k,368)*y(k,185)
         mat(k,1294) = .450_r8*rxt(k,348)*y(k,186)
         mat(k,446) = mat(k,446) + .130_r8*rxt(k,337)*y(k,1)
         mat(k,545) = .200_r8*rxt(k,387)*y(k,186)
         mat(k,407) = .400_r8*rxt(k,477)*y(k,186)
         mat(k,646) = .400_r8*rxt(k,481)*y(k,186)
         mat(k,991) = .400_r8*rxt(k,483)*y(k,186)
         mat(k,863) = mat(k,863) + .630_r8*rxt(k,498)*y(k,1)
         mat(k,904) = mat(k,904) + .630_r8*rxt(k,499)*y(k,1)
         mat(k,1243) = mat(k,1243) + .360_r8*rxt(k,384)*y(k,1)
         mat(k,1112) = mat(k,1112) + .240_r8*rxt(k,390)*y(k,1)
         mat(k,241) = mat(k,241) + .100_r8*rxt(k,397)*y(k,185)
         mat(k,1261) = .450_r8*rxt(k,400)*y(k,186)
         mat(k,246) = mat(k,246) + .500_r8*rxt(k,343)*y(k,185)
         mat(k,1051) = .150_r8*rxt(k,371)*y(k,186)
         mat(k,580) = mat(k,580) + .600_r8*rxt(k,421)*y(k,185)
         mat(k,198) = mat(k,198) + .650_r8*rxt(k,330)*y(k,185)
      end do
      end subroutine nlnmat01
      subroutine nlnmat02( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,740) = -(rxt(k,226)*y(k,185) + rxt(k,227)*y(k,3) + rxt(k,228)*y(k,6) &
                      + (rxt(k,229) + rxt(k,230) + rxt(k,231)) * y(k,7) + rxt(k,570) &
                      *y(k,219))
         mat(k,1811) = -rxt(k,226)*y(k,5)
         mat(k,1407) = -rxt(k,227)*y(k,5)
         mat(k,2087) = -rxt(k,228)*y(k,5)
         mat(k,1949) = -(rxt(k,229) + rxt(k,230) + rxt(k,231)) * y(k,5)
         mat(k,680) = -rxt(k,570)*y(k,5)
         mat(k,1910) = rxt(k,574)*y(k,221) + rxt(k,225)*y(k,226)
         mat(k,1407) = mat(k,1407) + rxt(k,572)*y(k,221)
         mat(k,554) = 1.100_r8*rxt(k,579)*y(k,223)
         mat(k,400) = rxt(k,574)*y(k,2) + rxt(k,572)*y(k,3)
         mat(k,688) = .200_r8*rxt(k,577)*y(k,223)
         mat(k,427) = rxt(k,225)*y(k,2)
         mat(k,698) = 1.100_r8*rxt(k,579)*y(k,222) + .200_r8*rxt(k,577)*y(k,220)
         mat(k,2132) = -(rxt(k,228)*y(k,5) + rxt(k,232)*y(k,2) + rxt(k,233)*y(k,186) &
                      + rxt(k,234)*y(k,1) + rxt(k,242)*y(k,8) + rxt(k,263)*y(k,25) &
                      + rxt(k,284)*y(k,32) + rxt(k,317)*y(k,187) + rxt(k,325)*y(k,188) &
                      + rxt(k,333)*y(k,191) + rxt(k,339)*y(k,189) + rxt(k,346) &
                      *y(k,190) + rxt(k,361)*y(k,193) + rxt(k,366)*y(k,194) + rxt(k,370) &
                      *y(k,195) + (rxt(k,380) + rxt(k,381)) * y(k,196) + rxt(k,386) &
                      *y(k,197) + (rxt(k,391) + rxt(k,392)) * y(k,199) + rxt(k,398) &
                      *y(k,198) + rxt(k,413)*y(k,201) + rxt(k,414)*y(k,202) + rxt(k,429) &
                      *y(k,79) + rxt(k,438)*y(k,203) + (rxt(k,446) + rxt(k,447) &
                      ) * y(k,200) + rxt(k,454)*y(k,207) + rxt(k,459)*y(k,208) &
                      + rxt(k,462)*y(k,206) + rxt(k,466)*y(k,211) + rxt(k,472) &
                      *y(k,212) + rxt(k,476)*y(k,213) + rxt(k,480)*y(k,204) + rxt(k,482) &
                      *y(k,214) + rxt(k,484)*y(k,215) + rxt(k,489)*y(k,209) + rxt(k,494) &
                      *y(k,210) + rxt(k,502)*y(k,205) + rxt(k,509)*y(k,216) + rxt(k,513) &
                      *y(k,217) + rxt(k,571)*y(k,219))
         mat(k,747) = -rxt(k,228)*y(k,6)
         mat(k,1935) = -rxt(k,232)*y(k,6)
         mat(k,1635) = -rxt(k,233)*y(k,6)
         mat(k,1719) = -rxt(k,234)*y(k,6)
         mat(k,1541) = -rxt(k,242)*y(k,6)
         mat(k,1894) = -rxt(k,263)*y(k,6)
         mat(k,1484) = -rxt(k,284)*y(k,6)
         mat(k,1394) = -rxt(k,317)*y(k,6)
         mat(k,319) = -rxt(k,325)*y(k,6)
         mat(k,658) = -rxt(k,333)*y(k,6)
         mat(k,830) = -rxt(k,339)*y(k,6)
         mat(k,1298) = -rxt(k,346)*y(k,6)
         mat(k,734) = -rxt(k,361)*y(k,6)
         mat(k,722) = -rxt(k,366)*y(k,6)
         mat(k,1054) = -rxt(k,370)*y(k,6)
         mat(k,438) = -(rxt(k,380) + rxt(k,381)) * y(k,6)
         mat(k,547) = -rxt(k,386)*y(k,6)
         mat(k,1225) = -(rxt(k,391) + rxt(k,392)) * y(k,6)
         mat(k,1265) = -rxt(k,398)*y(k,6)
         mat(k,1169) = -rxt(k,413)*y(k,6)
         mat(k,1203) = -rxt(k,414)*y(k,6)
         mat(k,1100) = -rxt(k,429)*y(k,6)
         mat(k,1141) = -rxt(k,438)*y(k,6)
         mat(k,939) = -(rxt(k,446) + rxt(k,447)) * y(k,6)
         mat(k,326) = -rxt(k,454)*y(k,6)
         mat(k,517) = -rxt(k,459)*y(k,6)
         mat(k,391) = -rxt(k,462)*y(k,6)
         mat(k,623) = -rxt(k,466)*y(k,6)
         mat(k,289) = -rxt(k,472)*y(k,6)
         mat(k,409) = -rxt(k,476)*y(k,6)
         mat(k,617) = -rxt(k,480)*y(k,6)
         mat(k,648) = -rxt(k,482)*y(k,6)
         mat(k,993) = -rxt(k,484)*y(k,6)
         mat(k,354) = -rxt(k,489)*y(k,6)
         mat(k,640) = -rxt(k,494)*y(k,6)
         mat(k,818) = -rxt(k,502)*y(k,6)
         mat(k,1016) = -rxt(k,509)*y(k,6)
         mat(k,974) = -rxt(k,513)*y(k,6)
         mat(k,684) = -rxt(k,571)*y(k,6)
         mat(k,1935) = mat(k,1935) + rxt(k,235)*y(k,7)
         mat(k,1867) = rxt(k,226)*y(k,5)
         mat(k,747) = mat(k,747) + rxt(k,226)*y(k,185) + 2.000_r8*rxt(k,230)*y(k,7) &
                      + rxt(k,227)*y(k,3)
         mat(k,1977) = rxt(k,235)*y(k,2) + 2.000_r8*rxt(k,230)*y(k,5) + rxt(k,537) &
                      *y(k,144)
         mat(k,1315) = rxt(k,537)*y(k,7)
         mat(k,1425) = rxt(k,227)*y(k,5) + rxt(k,224)*y(k,226)
         mat(k,430) = rxt(k,224)*y(k,3)
         mat(k,1973) = -((rxt(k,229) + rxt(k,230) + rxt(k,231)) * y(k,5) + (rxt(k,235) &
                      + rxt(k,236)) * y(k,2) + rxt(k,237)*y(k,1) + rxt(k,238)*y(k,8) &
                      + rxt(k,240)*y(k,185) + rxt(k,246)*y(k,186) + rxt(k,264)*y(k,25) &
                      + rxt(k,285)*y(k,32) + rxt(k,347)*y(k,190) + rxt(k,404)*y(k,198) &
                      + rxt(k,457)*y(k,89) + rxt(k,465)*y(k,211) + rxt(k,474)*y(k,213) &
                      + rxt(k,485)*y(k,214) + rxt(k,486)*y(k,215) + rxt(k,537) &
                      *y(k,144))
         mat(k,745) = -(rxt(k,229) + rxt(k,230) + rxt(k,231)) * y(k,7)
         mat(k,1931) = -(rxt(k,235) + rxt(k,236)) * y(k,7)
         mat(k,1715) = -rxt(k,237)*y(k,7)
         mat(k,1537) = -rxt(k,238)*y(k,7)
         mat(k,1863) = -rxt(k,240)*y(k,7)
         mat(k,1631) = -rxt(k,246)*y(k,7)
         mat(k,1890) = -rxt(k,264)*y(k,7)
         mat(k,1480) = -rxt(k,285)*y(k,7)
         mat(k,1295) = -rxt(k,347)*y(k,7)
         mat(k,1262) = -rxt(k,404)*y(k,7)
         mat(k,395) = -rxt(k,457)*y(k,7)
         mat(k,622) = -rxt(k,465)*y(k,7)
         mat(k,408) = -rxt(k,474)*y(k,7)
         mat(k,647) = -rxt(k,485)*y(k,7)
         mat(k,992) = -rxt(k,486)*y(k,7)
         mat(k,1313) = -rxt(k,537)*y(k,7)
         mat(k,1715) = mat(k,1715) + rxt(k,234)*y(k,6)
         mat(k,1931) = mat(k,1931) + rxt(k,232)*y(k,6) + rxt(k,243)*y(k,8)
         mat(k,1863) = mat(k,1863) + rxt(k,244)*y(k,8) + rxt(k,247)*y(k,10) &
                      + rxt(k,448)*y(k,63) + rxt(k,377)*y(k,64) + .700_r8*rxt(k,415) &
                      *y(k,66) + rxt(k,517)*y(k,71)
         mat(k,2128) = rxt(k,234)*y(k,1) + rxt(k,232)*y(k,2) + 2.000_r8*rxt(k,242) &
                      *y(k,8) + rxt(k,317)*y(k,187) + rxt(k,233)*y(k,186) + rxt(k,263) &
                      *y(k,25) + rxt(k,284)*y(k,32) + rxt(k,366)*y(k,194) + rxt(k,346) &
                      *y(k,190) + rxt(k,380)*y(k,196) + rxt(k,446)*y(k,200) &
                      + rxt(k,386)*y(k,197) + rxt(k,480)*y(k,204) + .800_r8*rxt(k,502) &
                      *y(k,205) + rxt(k,462)*y(k,206) + rxt(k,454)*y(k,207) &
                      + rxt(k,459)*y(k,208) + rxt(k,466)*y(k,211) + rxt(k,472) &
                      *y(k,212) + rxt(k,476)*y(k,213) + rxt(k,482)*y(k,214) &
                      + rxt(k,484)*y(k,215) + rxt(k,489)*y(k,209) + rxt(k,494) &
                      *y(k,210) + .900_r8*rxt(k,509)*y(k,216) + 1.600_r8*rxt(k,513) &
                      *y(k,217) + .920_r8*rxt(k,413)*y(k,201) + .920_r8*rxt(k,414) &
                      *y(k,202) + rxt(k,391)*y(k,199) + rxt(k,398)*y(k,198) &
                      + rxt(k,339)*y(k,189) + rxt(k,361)*y(k,193) + rxt(k,333) &
                      *y(k,191) + rxt(k,370)*y(k,195) + rxt(k,429)*y(k,79) &
                      + rxt(k,438)*y(k,203) + rxt(k,325)*y(k,188)
         mat(k,1537) = mat(k,1537) + rxt(k,243)*y(k,2) + rxt(k,244)*y(k,185) &
                      + 2.000_r8*rxt(k,242)*y(k,6) + rxt(k,245)*y(k,186) + rxt(k,382) &
                      *y(k,55) + 2.000_r8*rxt(k,516)*y(k,217) + rxt(k,417)*y(k,201) &
                      + rxt(k,418)*y(k,202) + rxt(k,393)*y(k,199) + rxt(k,399) &
                      *y(k,198) + rxt(k,430)*y(k,79) + rxt(k,439)*y(k,203)
         mat(k,374) = rxt(k,247)*y(k,185)
         mat(k,1391) = rxt(k,317)*y(k,6) + .500_r8*rxt(k,515)*y(k,217)
         mat(k,1631) = mat(k,1631) + rxt(k,233)*y(k,6) + rxt(k,245)*y(k,8)
         mat(k,1890) = mat(k,1890) + rxt(k,263)*y(k,6)
         mat(k,1480) = mat(k,1480) + rxt(k,284)*y(k,6)
         mat(k,720) = rxt(k,366)*y(k,6)
         mat(k,1295) = mat(k,1295) + rxt(k,346)*y(k,6)
         mat(k,539) = rxt(k,448)*y(k,185)
         mat(k,801) = rxt(k,377)*y(k,185)
         mat(k,474) = .700_r8*rxt(k,415)*y(k,185)
         mat(k,419) = rxt(k,517)*y(k,185)
         mat(k,275) = rxt(k,382)*y(k,8)
         mat(k,437) = rxt(k,380)*y(k,6)
         mat(k,937) = rxt(k,446)*y(k,6)
         mat(k,546) = rxt(k,386)*y(k,6)
         mat(k,616) = rxt(k,480)*y(k,6)
         mat(k,817) = .800_r8*rxt(k,502)*y(k,6)
         mat(k,390) = rxt(k,462)*y(k,6)
         mat(k,325) = rxt(k,454)*y(k,6)
         mat(k,516) = rxt(k,459)*y(k,6)
         mat(k,622) = mat(k,622) + rxt(k,466)*y(k,6)
         mat(k,288) = rxt(k,472)*y(k,6)
         mat(k,408) = mat(k,408) + rxt(k,476)*y(k,6)
         mat(k,647) = mat(k,647) + rxt(k,482)*y(k,6)
         mat(k,992) = mat(k,992) + rxt(k,484)*y(k,6)
         mat(k,353) = rxt(k,489)*y(k,6)
         mat(k,639) = rxt(k,494)*y(k,6)
         mat(k,1014) = .900_r8*rxt(k,509)*y(k,6)
         mat(k,973) = 1.600_r8*rxt(k,513)*y(k,6) + 2.000_r8*rxt(k,516)*y(k,8) &
                      + .500_r8*rxt(k,515)*y(k,187)
         mat(k,1166) = .920_r8*rxt(k,413)*y(k,6) + rxt(k,417)*y(k,8)
         mat(k,1200) = .920_r8*rxt(k,414)*y(k,6) + rxt(k,418)*y(k,8)
         mat(k,1222) = rxt(k,391)*y(k,6) + rxt(k,393)*y(k,8)
         mat(k,1262) = mat(k,1262) + rxt(k,398)*y(k,6) + rxt(k,399)*y(k,8)
         mat(k,829) = rxt(k,339)*y(k,6)
         mat(k,732) = rxt(k,361)*y(k,6)
         mat(k,657) = rxt(k,333)*y(k,6)
         mat(k,1052) = rxt(k,370)*y(k,6)
         mat(k,1099) = rxt(k,429)*y(k,6) + rxt(k,430)*y(k,8)
         mat(k,1138) = rxt(k,438)*y(k,6) + rxt(k,439)*y(k,8)
         mat(k,318) = rxt(k,325)*y(k,6)
         mat(k,1530) = -(rxt(k,238)*y(k,7) + rxt(k,242)*y(k,6) + rxt(k,243)*y(k,2) &
                      + rxt(k,244)*y(k,185) + rxt(k,245)*y(k,186) + rxt(k,313)*y(k,15) &
                      + rxt(k,345)*y(k,40) + rxt(k,360)*y(k,47) + rxt(k,376)*y(k,53) &
                      + rxt(k,382)*y(k,55) + rxt(k,393)*y(k,199) + rxt(k,399)*y(k,198) &
                      + rxt(k,412)*y(k,73) + rxt(k,417)*y(k,201) + rxt(k,418)*y(k,202) &
                      + rxt(k,430)*y(k,79) + rxt(k,439)*y(k,203) + rxt(k,500)*y(k,105) &
                      + rxt(k,501)*y(k,106) + rxt(k,508)*y(k,108) + rxt(k,516) &
                      *y(k,217) + rxt(k,545)*y(k,138))
         mat(k,1966) = -rxt(k,238)*y(k,8)
         mat(k,2121) = -rxt(k,242)*y(k,8)
         mat(k,1924) = -rxt(k,243)*y(k,8)
         mat(k,1856) = -rxt(k,244)*y(k,8)
         mat(k,1624) = -rxt(k,245)*y(k,8)
         mat(k,1648) = -rxt(k,313)*y(k,8)
         mat(k,1059) = -rxt(k,345)*y(k,8)
         mat(k,952) = -rxt(k,360)*y(k,8)
         mat(k,1118) = -rxt(k,376)*y(k,8)
         mat(k,272) = -rxt(k,382)*y(k,8)
         mat(k,1218) = -rxt(k,393)*y(k,8)
         mat(k,1257) = -rxt(k,399)*y(k,8)
         mat(k,790) = -rxt(k,412)*y(k,8)
         mat(k,1161) = -rxt(k,417)*y(k,8)
         mat(k,1195) = -rxt(k,418)*y(k,8)
         mat(k,1095) = -rxt(k,430)*y(k,8)
         mat(k,1134) = -rxt(k,439)*y(k,8)
         mat(k,859) = -rxt(k,500)*y(k,8)
         mat(k,900) = -rxt(k,501)*y(k,8)
         mat(k,913) = -rxt(k,508)*y(k,8)
         mat(k,969) = -rxt(k,516)*y(k,8)
         mat(k,201) = -rxt(k,545)*y(k,8)
         mat(k,1708) = rxt(k,237)*y(k,7)
         mat(k,1924) = mat(k,1924) + rxt(k,236)*y(k,7) + rxt(k,275)*y(k,30) &
                      + rxt(k,293)*y(k,35)
         mat(k,1856) = mat(k,1856) + rxt(k,241)*y(k,9) + rxt(k,276)*y(k,30) &
                      + rxt(k,356)*y(k,46) + .500_r8*rxt(k,406)*y(k,62)
         mat(k,1966) = mat(k,1966) + rxt(k,237)*y(k,1) + rxt(k,236)*y(k,2)
         mat(k,2040) = rxt(k,241)*y(k,185) + rxt(k,297)*y(k,136)
         mat(k,1449) = rxt(k,277)*y(k,30)
         mat(k,872) = rxt(k,275)*y(k,2) + rxt(k,276)*y(k,185) + rxt(k,277)*y(k,183)
         mat(k,479) = rxt(k,293)*y(k,2)
         mat(k,663) = rxt(k,297)*y(k,9)
         mat(k,358) = rxt(k,356)*y(k,185)
         mat(k,450) = .500_r8*rxt(k,406)*y(k,185)
      end do
      end subroutine nlnmat02
      subroutine nlnmat03( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,2050) = -(rxt(k,241)*y(k,185) + rxt(k,297)*y(k,136))
         mat(k,1866) = -rxt(k,241)*y(k,9)
         mat(k,666) = -rxt(k,297)*y(k,9)
         mat(k,1866) = mat(k,1866) + rxt(k,240)*y(k,7)
         mat(k,1976) = rxt(k,240)*y(k,185)
         mat(k,1540) = rxt(k,313)*y(k,15) + rxt(k,345)*y(k,40) + rxt(k,376)*y(k,53) &
                      + rxt(k,545)*y(k,138)
         mat(k,1658) = rxt(k,313)*y(k,8)
         mat(k,1333) = (rxt(k,551)+rxt(k,556)+rxt(k,562))*y(k,30)
         mat(k,878) = (rxt(k,551)+rxt(k,556)+rxt(k,562))*y(k,28)
         mat(k,1063) = rxt(k,345)*y(k,8)
         mat(k,1122) = rxt(k,376)*y(k,8)
         mat(k,204) = rxt(k,545)*y(k,8)
         mat(k,369) = -(rxt(k,247)*y(k,185))
         mat(k,1776) = -rxt(k,247)*y(k,10)
         mat(k,1940) = rxt(k,246)*y(k,186)
         mat(k,1565) = rxt(k,246)*y(k,7)
         mat(k,1938) = rxt(k,238)*y(k,8)
         mat(k,1489) = rxt(k,238)*y(k,7)
         mat(k,1382) = -(rxt(k,262)*y(k,25) + rxt(k,317)*y(k,6) + rxt(k,318)*y(k,186) &
                      + (4._r8*rxt(k,319) + 4._r8*rxt(k,320)) * y(k,187) + rxt(k,341) &
                      *y(k,189) + rxt(k,349)*y(k,190) + rxt(k,363)*y(k,193) + rxt(k,372) &
                      *y(k,195) + rxt(k,395)*y(k,199) + rxt(k,401)*y(k,198) + rxt(k,423) &
                      *y(k,201) + rxt(k,424)*y(k,202) + rxt(k,433)*y(k,79) + rxt(k,441) &
                      *y(k,203) + rxt(k,504)*y(k,205) + rxt(k,511)*y(k,216) + rxt(k,515) &
                      *y(k,217))
         mat(k,1879) = -rxt(k,262)*y(k,187)
         mat(k,2117) = -rxt(k,317)*y(k,187)
         mat(k,1620) = -rxt(k,318)*y(k,187)
         mat(k,824) = -rxt(k,341)*y(k,187)
         mat(k,1288) = -rxt(k,349)*y(k,187)
         mat(k,727) = -rxt(k,363)*y(k,187)
         mat(k,1048) = -rxt(k,372)*y(k,187)
         mat(k,1217) = -rxt(k,395)*y(k,187)
         mat(k,1256) = -rxt(k,401)*y(k,187)
         mat(k,1160) = -rxt(k,423)*y(k,187)
         mat(k,1194) = -rxt(k,424)*y(k,187)
         mat(k,1094) = -rxt(k,433)*y(k,187)
         mat(k,1133) = -rxt(k,441)*y(k,187)
         mat(k,813) = -rxt(k,504)*y(k,187)
         mat(k,1010) = -rxt(k,511)*y(k,187)
         mat(k,968) = -rxt(k,515)*y(k,187)
         mat(k,1704) = .280_r8*rxt(k,359)*y(k,47) + .050_r8*rxt(k,411)*y(k,73)
         mat(k,1852) = .700_r8*rxt(k,322)*y(k,13) + rxt(k,338)*y(k,39)
         mat(k,2117) = mat(k,2117) + rxt(k,346)*y(k,190) + .830_r8*rxt(k,482)*y(k,214) &
                      + .170_r8*rxt(k,484)*y(k,215)
         mat(k,1382) = mat(k,1382) + .900_r8*rxt(k,349)*y(k,190)
         mat(k,328) = .700_r8*rxt(k,322)*y(k,185)
         mat(k,1620) = mat(k,1620) + .450_r8*rxt(k,348)*y(k,190) + .330_r8*rxt(k,481) &
                      *y(k,214) + .070_r8*rxt(k,483)*y(k,215)
         mat(k,950) = .280_r8*rxt(k,359)*y(k,1)
         mat(k,789) = .050_r8*rxt(k,411)*y(k,1)
         mat(k,456) = rxt(k,338)*y(k,185)
         mat(k,1288) = mat(k,1288) + rxt(k,346)*y(k,6) + .900_r8*rxt(k,349)*y(k,187) &
                      + .450_r8*rxt(k,348)*y(k,186) + 4.000_r8*rxt(k,350)*y(k,190) &
                      + rxt(k,425)*y(k,201) + rxt(k,426)*y(k,202) + rxt(k,396) &
                      *y(k,199) + rxt(k,402)*y(k,198) + rxt(k,432)*y(k,79) &
                      + rxt(k,442)*y(k,203)
         mat(k,644) = .830_r8*rxt(k,482)*y(k,6) + .330_r8*rxt(k,481)*y(k,186)
         mat(k,989) = .170_r8*rxt(k,484)*y(k,6) + .070_r8*rxt(k,483)*y(k,186)
         mat(k,1160) = mat(k,1160) + rxt(k,425)*y(k,190)
         mat(k,1194) = mat(k,1194) + rxt(k,426)*y(k,190)
         mat(k,1217) = mat(k,1217) + rxt(k,396)*y(k,190)
         mat(k,1256) = mat(k,1256) + rxt(k,402)*y(k,190)
         mat(k,1094) = mat(k,1094) + rxt(k,432)*y(k,190)
         mat(k,1133) = mat(k,1133) + rxt(k,442)*y(k,190)
         mat(k,327) = -(rxt(k,322)*y(k,185))
         mat(k,1770) = -rxt(k,322)*y(k,13)
         mat(k,1351) = rxt(k,318)*y(k,186)
         mat(k,1559) = rxt(k,318)*y(k,187)
         mat(k,142) = -(rxt(k,200)*y(k,227) + rxt(k,219)*y(k,185))
         mat(k,2007) = -rxt(k,200)*y(k,17)
         mat(k,1742) = -rxt(k,219)*y(k,17)
         mat(k,50) = -(rxt(k,220)*y(k,185))
         mat(k,1729) = -rxt(k,220)*y(k,18)
         mat(k,1650) = -(rxt(k,256)*y(k,183) + rxt(k,280)*y(k,184) + rxt(k,313)*y(k,8) &
                      + rxt(k,314)*y(k,185) + rxt(k,315)*y(k,2) + rxt(k,316)*y(k,186))
         mat(k,1451) = -rxt(k,256)*y(k,15)
         mat(k,1341) = -rxt(k,280)*y(k,15)
         mat(k,1532) = -rxt(k,313)*y(k,15)
         mat(k,1858) = -rxt(k,314)*y(k,15)
         mat(k,1926) = -rxt(k,315)*y(k,15)
         mat(k,1626) = -rxt(k,316)*y(k,15)
         mat(k,1710) = .500_r8*rxt(k,359)*y(k,47) + .910_r8*rxt(k,411)*y(k,73) &
                      + rxt(k,337)*y(k,36) + .340_r8*rxt(k,498)*y(k,105) &
                      + .340_r8*rxt(k,499)*y(k,106) + .600_r8*rxt(k,384)*y(k,59) &
                      + .120_r8*rxt(k,390)*y(k,60)
         mat(k,1858) = mat(k,1858) + .300_r8*rxt(k,322)*y(k,13) + .500_r8*rxt(k,351) &
                      *y(k,44) + rxt(k,356)*y(k,46) + .500_r8*rxt(k,406)*y(k,62) &
                      + .400_r8*rxt(k,448)*y(k,63) + .300_r8*rxt(k,415)*y(k,66) &
                      + .680_r8*rxt(k,512)*y(k,109) + rxt(k,321)*y(k,14) &
                      + .800_r8*rxt(k,352)*y(k,42)
         mat(k,2123) = rxt(k,317)*y(k,187) + rxt(k,366)*y(k,194) + .500_r8*rxt(k,380) &
                      *y(k,196) + .100_r8*rxt(k,446)*y(k,200) + .320_r8*rxt(k,502) &
                      *y(k,205) + .340_r8*rxt(k,509)*y(k,216) + .920_r8*rxt(k,413) &
                      *y(k,201) + .250_r8*rxt(k,391)*y(k,199) + rxt(k,398)*y(k,198) &
                      + .500_r8*rxt(k,333)*y(k,191) + rxt(k,370)*y(k,195) &
                      + .250_r8*rxt(k,438)*y(k,203)
         mat(k,1532) = mat(k,1532) + .500_r8*rxt(k,382)*y(k,55) + rxt(k,417)*y(k,201) &
                      + .250_r8*rxt(k,393)*y(k,199) + rxt(k,399)*y(k,198)
         mat(k,1387) = rxt(k,317)*y(k,6) + (4.000_r8*rxt(k,319)+2.000_r8*rxt(k,320)) &
                      *y(k,187) + rxt(k,262)*y(k,25) + rxt(k,349)*y(k,190) &
                      + .950_r8*rxt(k,504)*y(k,205) + .930_r8*rxt(k,511)*y(k,216) &
                      + .750_r8*rxt(k,515)*y(k,217) + 1.500_r8*rxt(k,423)*y(k,201) &
                      + .750_r8*rxt(k,424)*y(k,202) + .880_r8*rxt(k,395)*y(k,199) &
                      + 2.000_r8*rxt(k,401)*y(k,198) + .700_r8*rxt(k,341)*y(k,189) &
                      + rxt(k,363)*y(k,193) + .800_r8*rxt(k,372)*y(k,195) &
                      + .800_r8*rxt(k,433)*y(k,79) + .800_r8*rxt(k,441)*y(k,203)
         mat(k,329) = .300_r8*rxt(k,322)*y(k,185)
         mat(k,1626) = mat(k,1626) + .450_r8*rxt(k,400)*y(k,198) + .150_r8*rxt(k,371) &
                      *y(k,195)
         mat(k,1885) = rxt(k,262)*y(k,187)
         mat(k,954) = .500_r8*rxt(k,359)*y(k,1)
         mat(k,792) = .910_r8*rxt(k,411)*y(k,1)
         mat(k,718) = rxt(k,366)*y(k,6)
         mat(k,1292) = rxt(k,349)*y(k,187) + rxt(k,425)*y(k,201) + .250_r8*rxt(k,396) &
                      *y(k,199) + rxt(k,402)*y(k,198) + .250_r8*rxt(k,442)*y(k,203)
         mat(k,423) = .500_r8*rxt(k,351)*y(k,185)
         mat(k,359) = rxt(k,356)*y(k,185)
         mat(k,444) = rxt(k,337)*y(k,1)
         mat(k,452) = .500_r8*rxt(k,406)*y(k,185)
         mat(k,537) = .400_r8*rxt(k,448)*y(k,185)
         mat(k,472) = .300_r8*rxt(k,415)*y(k,185)
         mat(k,273) = .500_r8*rxt(k,382)*y(k,8)
         mat(k,436) = .500_r8*rxt(k,380)*y(k,6)
         mat(k,935) = .100_r8*rxt(k,446)*y(k,6)
         mat(k,815) = .320_r8*rxt(k,502)*y(k,6) + .950_r8*rxt(k,504)*y(k,187)
         mat(k,861) = .340_r8*rxt(k,498)*y(k,1)
         mat(k,902) = .340_r8*rxt(k,499)*y(k,1)
         mat(k,1012) = .340_r8*rxt(k,509)*y(k,6) + .930_r8*rxt(k,511)*y(k,187)
         mat(k,982) = .680_r8*rxt(k,512)*y(k,185)
         mat(k,971) = .750_r8*rxt(k,515)*y(k,187)
         mat(k,1163) = .920_r8*rxt(k,413)*y(k,6) + rxt(k,417)*y(k,8) &
                      + 1.500_r8*rxt(k,423)*y(k,187) + rxt(k,425)*y(k,190)
         mat(k,1197) = .750_r8*rxt(k,424)*y(k,187)
         mat(k,1241) = .600_r8*rxt(k,384)*y(k,1)
         mat(k,1110) = .120_r8*rxt(k,390)*y(k,1)
         mat(k,1220) = .250_r8*rxt(k,391)*y(k,6) + .250_r8*rxt(k,393)*y(k,8) &
                      + .880_r8*rxt(k,395)*y(k,187) + .250_r8*rxt(k,396)*y(k,190)
         mat(k,1259) = rxt(k,398)*y(k,6) + rxt(k,399)*y(k,8) + 2.000_r8*rxt(k,401) &
                      *y(k,187) + .450_r8*rxt(k,400)*y(k,186) + rxt(k,402)*y(k,190) &
                      + 4.000_r8*rxt(k,403)*y(k,198)
         mat(k,827) = .700_r8*rxt(k,341)*y(k,187)
         mat(k,730) = rxt(k,363)*y(k,187)
         mat(k,710) = rxt(k,321)*y(k,185)
         mat(k,998) = .800_r8*rxt(k,352)*y(k,185)
         mat(k,655) = .500_r8*rxt(k,333)*y(k,6)
         mat(k,1050) = rxt(k,370)*y(k,6) + .800_r8*rxt(k,372)*y(k,187) &
                      + .150_r8*rxt(k,371)*y(k,186)
         mat(k,1097) = .800_r8*rxt(k,433)*y(k,187)
         mat(k,1136) = .250_r8*rxt(k,438)*y(k,6) + .800_r8*rxt(k,441)*y(k,187) &
                      + .250_r8*rxt(k,442)*y(k,190)
         mat(k,2153) = -(rxt(k,202)*y(k,3) + rxt(k,203)*y(k,1) + (rxt(k,204) + rxt(k,205) &
                      + rxt(k,206)) * y(k,186))
         mat(k,1426) = -rxt(k,202)*y(k,20)
         mat(k,1720) = -rxt(k,203)*y(k,20)
         mat(k,1636) = -(rxt(k,204) + rxt(k,205) + rxt(k,206)) * y(k,20)
         mat(k,1936) = rxt(k,214)*y(k,19) + rxt(k,207)*y(k,185)
         mat(k,2030) = rxt(k,195)*y(k,19) + rxt(k,197)*y(k,28) + rxt(k,199)*y(k,33)
         mat(k,1077) = rxt(k,214)*y(k,2) + rxt(k,195)*y(k,227) + rxt(k,212)*y(k,185) &
                      + rxt(k,252)*y(k,183) + rxt(k,295)*y(k,136)
         mat(k,1035) = rxt(k,311)*y(k,185)
         mat(k,1868) = rxt(k,207)*y(k,2) + rxt(k,212)*y(k,19) + rxt(k,311)*y(k,16) &
                      + rxt(k,226)*y(k,5) + rxt(k,314)*y(k,15) + rxt(k,530)*y(k,142) &
                      + rxt(k,531)*y(k,143) + rxt(k,534)*y(k,144)
         mat(k,748) = rxt(k,226)*y(k,185)
         mat(k,1660) = rxt(k,314)*y(k,185)
         mat(k,1461) = rxt(k,252)*y(k,19)
         mat(k,1334) = rxt(k,197)*y(k,227)
         mat(k,773) = rxt(k,199)*y(k,227)
         mat(k,667) = rxt(k,295)*y(k,19)
         mat(k,267) = rxt(k,530)*y(k,185)
         mat(k,587) = rxt(k,531)*y(k,185)
         mat(k,1316) = rxt(k,534)*y(k,185)
         mat(k,1625) = -((rxt(k,204) + rxt(k,205) + rxt(k,206)) * y(k,20) + rxt(k,209) &
                      *y(k,185) + rxt(k,215)*y(k,2) + rxt(k,216)*y(k,1) &
                      + 4._r8*rxt(k,217)*y(k,186) + rxt(k,233)*y(k,6) + rxt(k,245) &
                      *y(k,8) + rxt(k,246)*y(k,7) + (rxt(k,254) + rxt(k,255) &
                      ) * y(k,183) + rxt(k,261)*y(k,25) + rxt(k,279)*y(k,184) &
                      + rxt(k,283)*y(k,32) + rxt(k,316)*y(k,15) + rxt(k,318)*y(k,187) &
                      + rxt(k,326)*y(k,188) + rxt(k,334)*y(k,191) + rxt(k,340) &
                      *y(k,189) + rxt(k,348)*y(k,190) + rxt(k,362)*y(k,193) + rxt(k,367) &
                      *y(k,194) + rxt(k,371)*y(k,195) + rxt(k,387)*y(k,197) + rxt(k,394) &
                      *y(k,199) + rxt(k,400)*y(k,198) + rxt(k,419)*y(k,201) + rxt(k,420) &
                      *y(k,202) + rxt(k,431)*y(k,79) + rxt(k,440)*y(k,203) + rxt(k,449) &
                      *y(k,200) + rxt(k,455)*y(k,207) + rxt(k,460)*y(k,208) + rxt(k,463) &
                      *y(k,206) + rxt(k,467)*y(k,211) + rxt(k,470)*y(k,212) + rxt(k,477) &
                      *y(k,213) + rxt(k,478)*y(k,204) + rxt(k,481)*y(k,214) + rxt(k,483) &
                      *y(k,215) + rxt(k,490)*y(k,209) + rxt(k,492)*y(k,210) + rxt(k,503) &
                      *y(k,205) + rxt(k,510)*y(k,216) + rxt(k,514)*y(k,217))
         mat(k,2142) = -(rxt(k,204) + rxt(k,205) + rxt(k,206)) * y(k,186)
         mat(k,1857) = -rxt(k,209)*y(k,186)
         mat(k,1925) = -rxt(k,215)*y(k,186)
         mat(k,1709) = -rxt(k,216)*y(k,186)
         mat(k,2122) = -rxt(k,233)*y(k,186)
         mat(k,1531) = -rxt(k,245)*y(k,186)
         mat(k,1967) = -rxt(k,246)*y(k,186)
         mat(k,1450) = -(rxt(k,254) + rxt(k,255)) * y(k,186)
         mat(k,1884) = -rxt(k,261)*y(k,186)
         mat(k,1340) = -rxt(k,279)*y(k,186)
         mat(k,1474) = -rxt(k,283)*y(k,186)
         mat(k,1649) = -rxt(k,316)*y(k,186)
         mat(k,1386) = -rxt(k,318)*y(k,186)
         mat(k,316) = -rxt(k,326)*y(k,186)
         mat(k,654) = -rxt(k,334)*y(k,186)
         mat(k,826) = -rxt(k,340)*y(k,186)
         mat(k,1291) = -rxt(k,348)*y(k,186)
         mat(k,729) = -rxt(k,362)*y(k,186)
         mat(k,717) = -rxt(k,367)*y(k,186)
         mat(k,1049) = -rxt(k,371)*y(k,186)
         mat(k,544) = -rxt(k,387)*y(k,186)
         mat(k,1219) = -rxt(k,394)*y(k,186)
         mat(k,1258) = -rxt(k,400)*y(k,186)
         mat(k,1162) = -rxt(k,419)*y(k,186)
         mat(k,1196) = -rxt(k,420)*y(k,186)
         mat(k,1096) = -rxt(k,431)*y(k,186)
         mat(k,1135) = -rxt(k,440)*y(k,186)
         mat(k,934) = -rxt(k,449)*y(k,186)
         mat(k,323) = -rxt(k,455)*y(k,186)
         mat(k,513) = -rxt(k,460)*y(k,186)
         mat(k,388) = -rxt(k,463)*y(k,186)
         mat(k,621) = -rxt(k,467)*y(k,186)
         mat(k,286) = -rxt(k,470)*y(k,186)
         mat(k,406) = -rxt(k,477)*y(k,186)
         mat(k,614) = -rxt(k,478)*y(k,186)
         mat(k,645) = -rxt(k,481)*y(k,186)
         mat(k,990) = -rxt(k,483)*y(k,186)
         mat(k,351) = -rxt(k,490)*y(k,186)
         mat(k,637) = -rxt(k,492)*y(k,186)
         mat(k,814) = -rxt(k,503)*y(k,186)
         mat(k,1011) = -rxt(k,510)*y(k,186)
         mat(k,970) = -rxt(k,514)*y(k,186)
         mat(k,1709) = mat(k,1709) + rxt(k,208)*y(k,185) + .280_r8*rxt(k,359)*y(k,47) &
                      + .370_r8*rxt(k,411)*y(k,73) + .130_r8*rxt(k,337)*y(k,36) &
                      + .570_r8*rxt(k,498)*y(k,105) + .570_r8*rxt(k,499)*y(k,106) &
                      + .280_r8*rxt(k,384)*y(k,59) + .140_r8*rxt(k,390)*y(k,60)
         mat(k,1925) = mat(k,1925) + rxt(k,315)*y(k,15) + rxt(k,218)*y(k,21)
         mat(k,1033) = rxt(k,312)*y(k,185)
         mat(k,1857) = mat(k,1857) + rxt(k,208)*y(k,1) + rxt(k,312)*y(k,16) &
                      + rxt(k,244)*y(k,8) + rxt(k,219)*y(k,17) + rxt(k,220)*y(k,18) &
                      + rxt(k,213)*y(k,21) + rxt(k,259)*y(k,25) + rxt(k,282)*y(k,32) &
                      + .500_r8*rxt(k,406)*y(k,62) + rxt(k,407)*y(k,65) &
                      + .300_r8*rxt(k,415)*y(k,66) + rxt(k,416)*y(k,67) + rxt(k,434) &
                      *y(k,68) + rxt(k,436)*y(k,69) + rxt(k,435)*y(k,70) &
                      + .280_r8*rxt(k,468)*y(k,83) + .730_r8*rxt(k,469)*y(k,84) &
                      + rxt(k,353)*y(k,43) + .650_r8*rxt(k,452)*y(k,86) &
                      + .380_r8*rxt(k,487)*y(k,101) + .800_r8*rxt(k,453)*y(k,87) &
                      + .630_r8*rxt(k,488)*y(k,102) + .200_r8*rxt(k,512)*y(k,109) &
                      + .200_r8*rxt(k,397)*y(k,61) + rxt(k,321)*y(k,14) + rxt(k,354) &
                      *y(k,41) + rxt(k,352)*y(k,42) + rxt(k,374)*y(k,52) &
                      + .350_r8*rxt(k,330)*y(k,131) + rxt(k,323)*y(k,132) + rxt(k,541) &
                      *y(k,137) + .500_r8*rxt(k,544)*y(k,138)
         mat(k,2122) = mat(k,2122) + rxt(k,317)*y(k,187) + rxt(k,366)*y(k,194) &
                      + rxt(k,380)*y(k,196) + rxt(k,446)*y(k,200) + rxt(k,480) &
                      *y(k,204) + .800_r8*rxt(k,502)*y(k,205) + rxt(k,462)*y(k,206) &
                      + rxt(k,454)*y(k,207) + .400_r8*rxt(k,466)*y(k,211) + rxt(k,472) &
                      *y(k,212) + .170_r8*rxt(k,482)*y(k,214) + .830_r8*rxt(k,484) &
                      *y(k,215) + rxt(k,489)*y(k,209) + rxt(k,494)*y(k,210) &
                      + .900_r8*rxt(k,509)*y(k,216) + .920_r8*rxt(k,413)*y(k,201) &
                      + .920_r8*rxt(k,414)*y(k,202) + .470_r8*rxt(k,391)*y(k,199) &
                      + rxt(k,339)*y(k,189) + rxt(k,361)*y(k,193) + .250_r8*rxt(k,333) &
                      *y(k,191) + rxt(k,429)*y(k,79) + rxt(k,438)*y(k,203) &
                      + rxt(k,325)*y(k,188)
         mat(k,1531) = mat(k,1531) + rxt(k,244)*y(k,185) + rxt(k,313)*y(k,15) &
                      + rxt(k,417)*y(k,201) + rxt(k,418)*y(k,202) + .470_r8*rxt(k,393) &
                      *y(k,199) + rxt(k,430)*y(k,79) + rxt(k,439)*y(k,203)
         mat(k,1386) = mat(k,1386) + rxt(k,317)*y(k,6) + 4.000_r8*rxt(k,319)*y(k,187) &
                      + rxt(k,262)*y(k,25) + .900_r8*rxt(k,349)*y(k,190) + rxt(k,504) &
                      *y(k,205) + rxt(k,511)*y(k,216) + .500_r8*rxt(k,515)*y(k,217) &
                      + rxt(k,423)*y(k,201) + rxt(k,424)*y(k,202) + .730_r8*rxt(k,395) &
                      *y(k,199) + rxt(k,401)*y(k,198) + rxt(k,341)*y(k,189) &
                      + rxt(k,363)*y(k,193) + .300_r8*rxt(k,372)*y(k,195) &
                      + 1.200_r8*rxt(k,433)*y(k,79) + .800_r8*rxt(k,441)*y(k,203)
         mat(k,143) = rxt(k,219)*y(k,185)
         mat(k,51) = rxt(k,220)*y(k,185)
         mat(k,1649) = mat(k,1649) + rxt(k,315)*y(k,2) + rxt(k,313)*y(k,8) &
                      + rxt(k,256)*y(k,183) + rxt(k,280)*y(k,184)
         mat(k,2142) = mat(k,2142) + rxt(k,202)*y(k,3)
         mat(k,1625) = mat(k,1625) + .160_r8*rxt(k,467)*y(k,211) + .070_r8*rxt(k,481) &
                      *y(k,214) + .330_r8*rxt(k,483)*y(k,215)
         mat(k,499) = rxt(k,218)*y(k,2) + rxt(k,213)*y(k,185) + rxt(k,253)*y(k,183)
         mat(k,1450) = mat(k,1450) + rxt(k,256)*y(k,15) + rxt(k,253)*y(k,21)
         mat(k,1884) = mat(k,1884) + rxt(k,259)*y(k,185) + rxt(k,262)*y(k,187)
         mat(k,1340) = mat(k,1340) + rxt(k,280)*y(k,15)
         mat(k,1474) = mat(k,1474) + rxt(k,282)*y(k,185)
         mat(k,953) = .280_r8*rxt(k,359)*y(k,1)
         mat(k,791) = .370_r8*rxt(k,411)*y(k,1)
         mat(k,717) = mat(k,717) + rxt(k,366)*y(k,6)
         mat(k,1291) = mat(k,1291) + .900_r8*rxt(k,349)*y(k,187) + rxt(k,425)*y(k,201) &
                      + rxt(k,426)*y(k,202) + .470_r8*rxt(k,396)*y(k,199) + rxt(k,432) &
                      *y(k,79) + rxt(k,442)*y(k,203)
         mat(k,443) = .130_r8*rxt(k,337)*y(k,1)
         mat(k,451) = .500_r8*rxt(k,406)*y(k,185)
         mat(k,1027) = rxt(k,407)*y(k,185)
         mat(k,471) = .300_r8*rxt(k,415)*y(k,185)
         mat(k,381) = rxt(k,416)*y(k,185)
         mat(k,336) = rxt(k,434)*y(k,185)
         mat(k,763) = rxt(k,436)*y(k,185)
         mat(k,235) = rxt(k,435)*y(k,185)
         mat(k,435) = rxt(k,380)*y(k,6)
         mat(k,934) = mat(k,934) + rxt(k,446)*y(k,6)
         mat(k,114) = .280_r8*rxt(k,468)*y(k,185)
         mat(k,119) = .730_r8*rxt(k,469)*y(k,185)
         mat(k,614) = mat(k,614) + rxt(k,480)*y(k,6)
         mat(k,814) = mat(k,814) + .800_r8*rxt(k,502)*y(k,6) + rxt(k,504)*y(k,187)
         mat(k,805) = rxt(k,353)*y(k,185)
         mat(k,91) = .650_r8*rxt(k,452)*y(k,185)
         mat(k,131) = .380_r8*rxt(k,487)*y(k,185)
         mat(k,100) = .800_r8*rxt(k,453)*y(k,185)
         mat(k,388) = mat(k,388) + rxt(k,462)*y(k,6)
         mat(k,323) = mat(k,323) + rxt(k,454)*y(k,6)
         mat(k,621) = mat(k,621) + .400_r8*rxt(k,466)*y(k,6) + .160_r8*rxt(k,467) &
                      *y(k,186)
         mat(k,286) = mat(k,286) + rxt(k,472)*y(k,6)
         mat(k,645) = mat(k,645) + .170_r8*rxt(k,482)*y(k,6) + .070_r8*rxt(k,481) &
                      *y(k,186)
         mat(k,990) = mat(k,990) + .830_r8*rxt(k,484)*y(k,6) + .330_r8*rxt(k,483) &
                      *y(k,186)
         mat(k,140) = .630_r8*rxt(k,488)*y(k,185)
         mat(k,351) = mat(k,351) + rxt(k,489)*y(k,6)
         mat(k,637) = mat(k,637) + rxt(k,494)*y(k,6)
         mat(k,860) = .570_r8*rxt(k,498)*y(k,1)
         mat(k,901) = .570_r8*rxt(k,499)*y(k,1)
         mat(k,1011) = mat(k,1011) + .900_r8*rxt(k,509)*y(k,6) + rxt(k,511)*y(k,187)
         mat(k,981) = .200_r8*rxt(k,512)*y(k,185)
         mat(k,970) = mat(k,970) + .500_r8*rxt(k,515)*y(k,187)
         mat(k,1162) = mat(k,1162) + .920_r8*rxt(k,413)*y(k,6) + rxt(k,417)*y(k,8) &
                      + rxt(k,423)*y(k,187) + rxt(k,425)*y(k,190)
         mat(k,1196) = mat(k,1196) + .920_r8*rxt(k,414)*y(k,6) + rxt(k,418)*y(k,8) &
                      + rxt(k,424)*y(k,187) + rxt(k,426)*y(k,190)
         mat(k,1240) = .280_r8*rxt(k,384)*y(k,1)
         mat(k,1109) = .140_r8*rxt(k,390)*y(k,1)
         mat(k,1219) = mat(k,1219) + .470_r8*rxt(k,391)*y(k,6) + .470_r8*rxt(k,393) &
                      *y(k,8) + .730_r8*rxt(k,395)*y(k,187) + .470_r8*rxt(k,396) &
                      *y(k,190)
         mat(k,240) = .200_r8*rxt(k,397)*y(k,185)
         mat(k,1258) = mat(k,1258) + rxt(k,401)*y(k,187)
         mat(k,826) = mat(k,826) + rxt(k,339)*y(k,6) + rxt(k,341)*y(k,187) &
                      + 2.400_r8*rxt(k,342)*y(k,189)
         mat(k,729) = mat(k,729) + rxt(k,361)*y(k,6) + rxt(k,363)*y(k,187)
         mat(k,709) = rxt(k,321)*y(k,185)
         mat(k,173) = rxt(k,354)*y(k,185)
         mat(k,997) = rxt(k,352)*y(k,185)
         mat(k,1039) = rxt(k,374)*y(k,185)
         mat(k,654) = mat(k,654) + .250_r8*rxt(k,333)*y(k,6)
         mat(k,345) = rxt(k,335)*y(k,3)
         mat(k,1049) = mat(k,1049) + .300_r8*rxt(k,372)*y(k,187)
         mat(k,1096) = mat(k,1096) + rxt(k,429)*y(k,6) + rxt(k,430)*y(k,8) &
                      + 1.200_r8*rxt(k,433)*y(k,187) + rxt(k,432)*y(k,190)
         mat(k,1135) = mat(k,1135) + rxt(k,438)*y(k,6) + rxt(k,439)*y(k,8) &
                      + .800_r8*rxt(k,441)*y(k,187) + rxt(k,442)*y(k,190)
         mat(k,197) = .350_r8*rxt(k,330)*y(k,185)
         mat(k,705) = rxt(k,323)*y(k,185)
         mat(k,316) = mat(k,316) + rxt(k,325)*y(k,6)
         mat(k,834) = rxt(k,541)*y(k,185)
         mat(k,202) = .500_r8*rxt(k,544)*y(k,185)
         mat(k,1416) = rxt(k,202)*y(k,20) + rxt(k,335)*y(k,192)
      end do
      end subroutine nlnmat03
      subroutine nlnmat04( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,496) = -(rxt(k,213)*y(k,185) + rxt(k,218)*y(k,2) + rxt(k,253)*y(k,183))
         mat(k,1792) = -rxt(k,213)*y(k,21)
         mat(k,1902) = -rxt(k,218)*y(k,21)
         mat(k,1431) = -rxt(k,253)*y(k,21)
         mat(k,1792) = mat(k,1792) + 2.000_r8*rxt(k,211)*y(k,185)
         mat(k,1573) = 2.000_r8*rxt(k,217)*y(k,186)
         mat(k,2000) = -(rxt(k,170)*y(k,227) + rxt(k,294)*y(k,136) + rxt(k,542) &
                      *y(k,145))
         mat(k,2026) = -rxt(k,170)*y(k,228)
         mat(k,665) = -rxt(k,294)*y(k,228)
         mat(k,183) = -rxt(k,542)*y(k,228)
         mat(k,1074) = rxt(k,212)*y(k,185)
         mat(k,1864) = rxt(k,212)*y(k,19) + 2.000_r8*rxt(k,210)*y(k,185) + rxt(k,241) &
                      *y(k,9) + rxt(k,247)*y(k,10) + rxt(k,322)*y(k,13) + rxt(k,314) &
                      *y(k,15) + rxt(k,209)*y(k,186) + rxt(k,213)*y(k,21) + rxt(k,270) &
                      *y(k,28) + rxt(k,274)*y(k,29) + rxt(k,290)*y(k,33) + rxt(k,344) &
                      *y(k,40) + rxt(k,338)*y(k,39) + rxt(k,368)*y(k,51) + rxt(k,351) &
                      *y(k,44) + rxt(k,331)*y(k,37) + .500_r8*rxt(k,389)*y(k,60) &
                      + rxt(k,365)*y(k,48) + rxt(k,364)*y(k,49) + rxt(k,369)*y(k,50) &
                      + rxt(k,373)*y(k,54) + rxt(k,375)*y(k,53) + (rxt(k,443) &
                       +rxt(k,444))*y(k,81) + rxt(k,323)*y(k,132) + rxt(k,527) &
                      *y(k,139)
         mat(k,2048) = rxt(k,241)*y(k,185)
         mat(k,375) = rxt(k,247)*y(k,185)
         mat(k,331) = rxt(k,322)*y(k,185)
         mat(k,1656) = rxt(k,314)*y(k,185)
         mat(k,2149) = rxt(k,206)*y(k,186)
         mat(k,1632) = rxt(k,209)*y(k,185) + rxt(k,206)*y(k,20)
         mat(k,502) = rxt(k,213)*y(k,185)
         mat(k,1331) = rxt(k,270)*y(k,185) + (rxt(k,552)+rxt(k,557)+rxt(k,563)) &
                      *y(k,29) + (rxt(k,553)+rxt(k,564))*y(k,34)
         mat(k,756) = rxt(k,274)*y(k,185) + (rxt(k,552)+rxt(k,557)+rxt(k,563))*y(k,28)
         mat(k,771) = rxt(k,290)*y(k,185)
         mat(k,676) = (rxt(k,553)+rxt(k,564))*y(k,28)
         mat(k,1062) = rxt(k,344)*y(k,185)
         mat(k,458) = rxt(k,338)*y(k,185)
         mat(k,466) = rxt(k,368)*y(k,185)
         mat(k,425) = rxt(k,351)*y(k,185)
         mat(k,189) = rxt(k,331)*y(k,185)
         mat(k,1113) = .500_r8*rxt(k,389)*y(k,185)
         mat(k,59) = rxt(k,365)*y(k,185)
         mat(k,307) = rxt(k,364)*y(k,185)
         mat(k,921) = rxt(k,369)*y(k,185)
         mat(k,313) = rxt(k,373)*y(k,185)
         mat(k,1121) = rxt(k,375)*y(k,185)
         mat(k,178) = (rxt(k,443)+rxt(k,444))*y(k,185)
         mat(k,707) = rxt(k,323)*y(k,185)
         mat(k,49) = rxt(k,527)*y(k,185)
         mat(k,1447) = -(rxt(k,251)*y(k,1) + rxt(k,252)*y(k,19) + rxt(k,253)*y(k,21) &
                      + (rxt(k,254) + rxt(k,255)) * y(k,186) + rxt(k,256)*y(k,15) &
                      + rxt(k,273)*y(k,29) + rxt(k,277)*y(k,30) + rxt(k,329)*y(k,37))
         mat(k,1706) = -rxt(k,251)*y(k,183)
         mat(k,1070) = -rxt(k,252)*y(k,183)
         mat(k,498) = -rxt(k,253)*y(k,183)
         mat(k,1622) = -(rxt(k,254) + rxt(k,255)) * y(k,183)
         mat(k,1646) = -rxt(k,256)*y(k,183)
         mat(k,752) = -rxt(k,273)*y(k,183)
         mat(k,871) = -rxt(k,277)*y(k,183)
         mat(k,187) = -rxt(k,329)*y(k,183)
         mat(k,1922) = rxt(k,258)*y(k,25) + rxt(k,271)*y(k,28)
         mat(k,2016) = rxt(k,196)*y(k,28) + rxt(k,191)*y(k,134)
         mat(k,1854) = rxt(k,259)*y(k,25) + rxt(k,270)*y(k,28)
         mat(k,2119) = rxt(k,263)*y(k,25)
         mat(k,1384) = rxt(k,262)*y(k,25)
         mat(k,1881) = rxt(k,258)*y(k,2) + rxt(k,259)*y(k,185) + rxt(k,263)*y(k,6) &
                      + rxt(k,262)*y(k,187) + (4.000_r8*rxt(k,265)+2.000_r8*rxt(k,267)) &
                      *y(k,25) + rxt(k,287)*y(k,32) + rxt(k,538)*y(k,144)
         mat(k,1324) = rxt(k,271)*y(k,2) + rxt(k,196)*y(k,227) + rxt(k,270)*y(k,185)
         mat(k,1471) = rxt(k,287)*y(k,25)
         mat(k,148) = rxt(k,191)*y(k,227)
         mat(k,1306) = rxt(k,538)*y(k,25)
         mat(k,1427) = rxt(k,277)*y(k,30)
         mat(k,1870) = 2.000_r8*rxt(k,266)*y(k,25)
         mat(k,1317) = (rxt(k,552)+rxt(k,557)+rxt(k,563))*y(k,29) + (rxt(k,551) &
                       +rxt(k,556)+rxt(k,562))*y(k,30)
         mat(k,749) = (rxt(k,552)+rxt(k,557)+rxt(k,563))*y(k,28)
         mat(k,867) = rxt(k,277)*y(k,183) + (rxt(k,551)+rxt(k,556)+rxt(k,562))*y(k,28)
         mat(k,1888) = -(rxt(k,258)*y(k,2) + (rxt(k,259) + rxt(k,260)) * y(k,185) &
                      + rxt(k,261)*y(k,186) + rxt(k,262)*y(k,187) + rxt(k,263)*y(k,6) &
                      + rxt(k,264)*y(k,7) + (4._r8*rxt(k,265) + 4._r8*rxt(k,266) &
                      + 4._r8*rxt(k,267) + 4._r8*rxt(k,268)) * y(k,25) + (rxt(k,286) &
                      + rxt(k,287) + rxt(k,288)) * y(k,32) + rxt(k,538)*y(k,144))
         mat(k,1929) = -rxt(k,258)*y(k,25)
         mat(k,1861) = -(rxt(k,259) + rxt(k,260)) * y(k,25)
         mat(k,1629) = -rxt(k,261)*y(k,25)
         mat(k,1390) = -rxt(k,262)*y(k,25)
         mat(k,2126) = -rxt(k,263)*y(k,25)
         mat(k,1971) = -rxt(k,264)*y(k,25)
         mat(k,1478) = -(rxt(k,286) + rxt(k,287) + rxt(k,288)) * y(k,25)
         mat(k,1311) = -rxt(k,538)*y(k,25)
         mat(k,1713) = rxt(k,251)*y(k,183)
         mat(k,1929) = mat(k,1929) + rxt(k,272)*y(k,29) + rxt(k,275)*y(k,30)
         mat(k,2023) = rxt(k,197)*y(k,28)
         mat(k,1861) = mat(k,1861) + rxt(k,274)*y(k,29)
         mat(k,1629) = mat(k,1629) + rxt(k,255)*y(k,183)
         mat(k,1454) = rxt(k,251)*y(k,1) + rxt(k,255)*y(k,186) + rxt(k,273)*y(k,29)
         mat(k,253) = rxt(k,540)*y(k,144)
         mat(k,1328) = rxt(k,197)*y(k,227)
         mat(k,754) = rxt(k,272)*y(k,2) + rxt(k,274)*y(k,185) + rxt(k,273)*y(k,183)
         mat(k,874) = rxt(k,275)*y(k,2)
         mat(k,1311) = mat(k,1311) + rxt(k,540)*y(k,26)
         mat(k,250) = -(rxt(k,540)*y(k,144))
         mat(k,1300) = -rxt(k,540)*y(k,26)
         mat(k,1872) = 2.000_r8*rxt(k,267)*y(k,25) + rxt(k,286)*y(k,32)
         mat(k,1463) = rxt(k,286)*y(k,25)
         mat(k,1869) = 2.000_r8*rxt(k,268)*y(k,25)
         mat(k,1322) = -((rxt(k,196) + rxt(k,197)) * y(k,227) + rxt(k,270)*y(k,185) &
                      + rxt(k,271)*y(k,2) + (rxt(k,551) + rxt(k,556) + rxt(k,562) &
                      ) * y(k,30) + (rxt(k,552) + rxt(k,557) + rxt(k,563)) * y(k,29) &
                      + (rxt(k,553) + rxt(k,564)) * y(k,34))
         mat(k,2012) = -(rxt(k,196) + rxt(k,197)) * y(k,28)
         mat(k,1850) = -rxt(k,270)*y(k,28)
         mat(k,1918) = -rxt(k,271)*y(k,28)
         mat(k,870) = -(rxt(k,551) + rxt(k,556) + rxt(k,562)) * y(k,28)
         mat(k,751) = -(rxt(k,552) + rxt(k,557) + rxt(k,563)) * y(k,28)
         mat(k,670) = -(rxt(k,553) + rxt(k,564)) * y(k,28)
         mat(k,1068) = rxt(k,252)*y(k,183)
         mat(k,1850) = mat(k,1850) + rxt(k,260)*y(k,25)
         mat(k,1642) = rxt(k,256)*y(k,183)
         mat(k,1618) = rxt(k,254)*y(k,183)
         mat(k,497) = rxt(k,253)*y(k,183)
         mat(k,1443) = rxt(k,252)*y(k,19) + rxt(k,256)*y(k,15) + rxt(k,254)*y(k,186) &
                      + rxt(k,253)*y(k,21) + rxt(k,273)*y(k,29) + rxt(k,329)*y(k,37)
         mat(k,1877) = rxt(k,260)*y(k,185)
         mat(k,751) = mat(k,751) + rxt(k,273)*y(k,183)
         mat(k,186) = rxt(k,329)*y(k,183)
         mat(k,750) = -(rxt(k,272)*y(k,2) + rxt(k,273)*y(k,183) + rxt(k,274)*y(k,185) &
                      + (rxt(k,552) + rxt(k,557) + rxt(k,563)) * y(k,28))
         mat(k,1911) = -rxt(k,272)*y(k,29)
         mat(k,1434) = -rxt(k,273)*y(k,29)
         mat(k,1812) = -rxt(k,274)*y(k,29)
         mat(k,1320) = -(rxt(k,552) + rxt(k,557) + rxt(k,563)) * y(k,29)
         mat(k,1812) = mat(k,1812) + rxt(k,276)*y(k,30)
         mat(k,1589) = rxt(k,261)*y(k,25)
         mat(k,1873) = rxt(k,261)*y(k,186)
         mat(k,868) = rxt(k,276)*y(k,185)
         mat(k,869) = -(rxt(k,275)*y(k,2) + rxt(k,276)*y(k,185) + rxt(k,277)*y(k,183) &
                      + (rxt(k,551) + rxt(k,556) + rxt(k,562)) * y(k,28))
         mat(k,1914) = -rxt(k,275)*y(k,30)
         mat(k,1822) = -rxt(k,276)*y(k,30)
         mat(k,1437) = -rxt(k,277)*y(k,30)
         mat(k,1321) = -(rxt(k,551) + rxt(k,556) + rxt(k,562)) * y(k,30)
         mat(k,1952) = rxt(k,264)*y(k,25)
         mat(k,1875) = rxt(k,264)*y(k,7)
         mat(k,1871) = rxt(k,288)*y(k,32)
         mat(k,1318) = (rxt(k,553)+rxt(k,564))*y(k,34)
         mat(k,1462) = rxt(k,288)*y(k,25)
         mat(k,668) = (rxt(k,553)+rxt(k,564))*y(k,28)
         mat(k,1337) = -(rxt(k,278)*y(k,1) + rxt(k,279)*y(k,186) + rxt(k,280)*y(k,15))
         mat(k,1703) = -rxt(k,278)*y(k,184)
         mat(k,1619) = -rxt(k,279)*y(k,184)
         mat(k,1643) = -rxt(k,280)*y(k,184)
         mat(k,1919) = rxt(k,281)*y(k,32) + rxt(k,291)*y(k,33)
         mat(k,2013) = rxt(k,198)*y(k,33)
         mat(k,1851) = rxt(k,282)*y(k,32) + rxt(k,290)*y(k,33)
         mat(k,2116) = rxt(k,284)*y(k,32)
         mat(k,1878) = (rxt(k,286)+rxt(k,287))*y(k,32)
         mat(k,1469) = rxt(k,281)*y(k,2) + rxt(k,282)*y(k,185) + rxt(k,284)*y(k,6) + ( &
                      + rxt(k,286)+rxt(k,287))*y(k,25) + 4.000_r8*rxt(k,289)*y(k,32) &
                      + rxt(k,539)*y(k,144)
         mat(k,767) = rxt(k,291)*y(k,2) + rxt(k,198)*y(k,227) + rxt(k,290)*y(k,185)
         mat(k,1304) = rxt(k,539)*y(k,32)
         mat(k,1472) = -(rxt(k,281)*y(k,2) + rxt(k,282)*y(k,185) + rxt(k,283)*y(k,186) &
                      + rxt(k,284)*y(k,6) + rxt(k,285)*y(k,7) + (rxt(k,286) + rxt(k,287) &
                      + rxt(k,288)) * y(k,25) + 4._r8*rxt(k,289)*y(k,32) + rxt(k,539) &
                      *y(k,144))
         mat(k,1923) = -rxt(k,281)*y(k,32)
         mat(k,1855) = -rxt(k,282)*y(k,32)
         mat(k,1623) = -rxt(k,283)*y(k,32)
         mat(k,2120) = -rxt(k,284)*y(k,32)
         mat(k,1965) = -rxt(k,285)*y(k,32)
         mat(k,1882) = -(rxt(k,286) + rxt(k,287) + rxt(k,288)) * y(k,32)
         mat(k,1307) = -rxt(k,539)*y(k,32)
         mat(k,1707) = rxt(k,278)*y(k,184)
         mat(k,1923) = mat(k,1923) + rxt(k,292)*y(k,34) + rxt(k,293)*y(k,35)
         mat(k,2017) = rxt(k,199)*y(k,33)
         mat(k,1339) = rxt(k,278)*y(k,1)
         mat(k,768) = rxt(k,199)*y(k,227)
         mat(k,673) = rxt(k,292)*y(k,2)
         mat(k,478) = rxt(k,293)*y(k,2)
         mat(k,766) = -((rxt(k,198) + rxt(k,199)) * y(k,227) + rxt(k,290)*y(k,185) &
                      + rxt(k,291)*y(k,2))
         mat(k,2010) = -(rxt(k,198) + rxt(k,199)) * y(k,33)
         mat(k,1814) = -rxt(k,290)*y(k,33)
         mat(k,1912) = -rxt(k,291)*y(k,33)
         mat(k,1639) = rxt(k,280)*y(k,184)
         mat(k,1590) = rxt(k,279)*y(k,184)
         mat(k,1335) = rxt(k,280)*y(k,15) + rxt(k,279)*y(k,186)
         mat(k,669) = -(rxt(k,292)*y(k,2) + (rxt(k,553) + rxt(k,564)) * y(k,28))
         mat(k,1906) = -rxt(k,292)*y(k,34)
         mat(k,1319) = -(rxt(k,553) + rxt(k,564)) * y(k,34)
         mat(k,1585) = rxt(k,283)*y(k,32)
         mat(k,1465) = rxt(k,283)*y(k,186)
         mat(k,475) = -(rxt(k,293)*y(k,2))
         mat(k,1901) = -rxt(k,293)*y(k,35)
         mat(k,1944) = rxt(k,285)*y(k,32)
         mat(k,1464) = rxt(k,285)*y(k,7)
         mat(k,74) = -(rxt(k,190)*y(k,227))
         mat(k,2005) = -rxt(k,190)*y(k,133)
         mat(k,146) = -(rxt(k,191)*y(k,227))
         mat(k,2008) = -rxt(k,191)*y(k,134)
         mat(k,1065) = rxt(k,295)*y(k,136)
         mat(k,2031) = rxt(k,297)*y(k,136)
         mat(k,1980) = rxt(k,294)*y(k,136)
         mat(k,659) = rxt(k,295)*y(k,19) + rxt(k,297)*y(k,9) + rxt(k,294)*y(k,228)
      end do
      end subroutine nlnmat04
      subroutine nlnmat05( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,660) = -(rxt(k,294)*y(k,228) + rxt(k,295)*y(k,19) + rxt(k,297)*y(k,9))
         mat(k,1982) = -rxt(k,294)*y(k,136)
         mat(k,1066) = -rxt(k,295)*y(k,136)
         mat(k,2032) = -rxt(k,297)*y(k,136)
         mat(k,2009) = 2.000_r8*rxt(k,190)*y(k,133) + rxt(k,191)*y(k,134)
         mat(k,75) = 2.000_r8*rxt(k,190)*y(k,227)
         mat(k,147) = rxt(k,191)*y(k,227)
         mat(k,944) = -(rxt(k,358)*y(k,185) + rxt(k,359)*y(k,1) + rxt(k,360)*y(k,8))
         mat(k,1827) = -rxt(k,358)*y(k,47)
         mat(k,1683) = -rxt(k,359)*y(k,47)
         mat(k,1503) = -rxt(k,360)*y(k,47)
         mat(k,1683) = mat(k,1683) + .130_r8*rxt(k,411)*y(k,73)
         mat(k,781) = .130_r8*rxt(k,411)*y(k,1)
         mat(k,780) = -(rxt(k,410)*y(k,185) + rxt(k,411)*y(k,1) + rxt(k,412)*y(k,8))
         mat(k,1815) = -rxt(k,410)*y(k,73)
         mat(k,1675) = -rxt(k,411)*y(k,73)
         mat(k,1494) = -rxt(k,412)*y(k,73)
         mat(k,713) = -(rxt(k,366)*y(k,6) + rxt(k,367)*y(k,186))
         mat(k,2085) = -rxt(k,366)*y(k,194)
         mat(k,1587) = -rxt(k,367)*y(k,194)
         mat(k,1809) = rxt(k,358)*y(k,47) + .500_r8*rxt(k,368)*y(k,51)
         mat(k,942) = rxt(k,358)*y(k,185)
         mat(k,460) = .500_r8*rxt(k,368)*y(k,185)
         mat(k,1056) = -(rxt(k,344)*y(k,185) + rxt(k,345)*y(k,8))
         mat(k,1837) = -rxt(k,344)*y(k,40)
         mat(k,1513) = -rxt(k,345)*y(k,40)
         mat(k,1691) = .500_r8*rxt(k,359)*y(k,47) + .100_r8*rxt(k,384)*y(k,59)
         mat(k,1837) = mat(k,1837) + .800_r8*rxt(k,448)*y(k,63) + .500_r8*rxt(k,343) &
                      *y(k,38) + rxt(k,354)*y(k,41)
         mat(k,2105) = rxt(k,366)*y(k,194) + rxt(k,380)*y(k,196) + .400_r8*rxt(k,446) &
                      *y(k,200) + rxt(k,386)*y(k,197) + rxt(k,339)*y(k,189) &
                      + .270_r8*rxt(k,361)*y(k,193)
         mat(k,1513) = mat(k,1513) + rxt(k,382)*y(k,55)
         mat(k,1371) = .800_r8*rxt(k,341)*y(k,189)
         mat(k,1606) = .200_r8*rxt(k,387)*y(k,197)
         mat(k,947) = .500_r8*rxt(k,359)*y(k,1)
         mat(k,715) = rxt(k,366)*y(k,6)
         mat(k,534) = .800_r8*rxt(k,448)*y(k,185)
         mat(k,271) = rxt(k,382)*y(k,8)
         mat(k,434) = rxt(k,380)*y(k,6)
         mat(k,930) = .400_r8*rxt(k,446)*y(k,6)
         mat(k,542) = rxt(k,386)*y(k,6) + .200_r8*rxt(k,387)*y(k,186)
         mat(k,1231) = .100_r8*rxt(k,384)*y(k,1)
         mat(k,823) = rxt(k,339)*y(k,6) + .800_r8*rxt(k,341)*y(k,187) &
                      + 3.200_r8*rxt(k,342)*y(k,189)
         mat(k,244) = .500_r8*rxt(k,343)*y(k,185)
         mat(k,726) = .270_r8*rxt(k,361)*y(k,6)
         mat(k,172) = rxt(k,354)*y(k,185)
         mat(k,455) = -(rxt(k,338)*y(k,185))
         mat(k,1787) = -rxt(k,338)*y(k,39)
         mat(k,1670) = .120_r8*rxt(k,359)*y(k,47)
         mat(k,1353) = .100_r8*rxt(k,349)*y(k,190)
         mat(k,1570) = .150_r8*rxt(k,348)*y(k,190) + .150_r8*rxt(k,400)*y(k,198)
         mat(k,940) = .120_r8*rxt(k,359)*y(k,1)
         mat(k,1270) = .100_r8*rxt(k,349)*y(k,187) + .150_r8*rxt(k,348)*y(k,186)
         mat(k,1251) = .150_r8*rxt(k,400)*y(k,186)
         mat(k,459) = -(rxt(k,368)*y(k,185))
         mat(k,1788) = -rxt(k,368)*y(k,51)
         mat(k,1571) = rxt(k,367)*y(k,194)
         mat(k,712) = rxt(k,367)*y(k,186)
         mat(k,1287) = -(rxt(k,346)*y(k,6) + rxt(k,347)*y(k,7) + rxt(k,348)*y(k,186) &
                      + rxt(k,349)*y(k,187) + 4._r8*rxt(k,350)*y(k,190) + rxt(k,396) &
                      *y(k,199) + rxt(k,425)*y(k,201) + rxt(k,426)*y(k,202) + rxt(k,432) &
                      *y(k,79) + rxt(k,442)*y(k,203))
         mat(k,2115) = -rxt(k,346)*y(k,190)
         mat(k,1958) = -rxt(k,347)*y(k,190)
         mat(k,1617) = -rxt(k,348)*y(k,190)
         mat(k,1381) = -rxt(k,349)*y(k,190)
         mat(k,1216) = -rxt(k,396)*y(k,190)
         mat(k,1159) = -rxt(k,425)*y(k,190)
         mat(k,1193) = -rxt(k,426)*y(k,190)
         mat(k,1093) = -rxt(k,432)*y(k,190)
         mat(k,1132) = -rxt(k,442)*y(k,190)
         mat(k,1701) = .080_r8*rxt(k,411)*y(k,73) + .060_r8*rxt(k,498)*y(k,105) &
                      + .060_r8*rxt(k,499)*y(k,106) + .280_r8*rxt(k,384)*y(k,59) &
                      + .100_r8*rxt(k,390)*y(k,60)
         mat(k,1848) = rxt(k,344)*y(k,40) + .500_r8*rxt(k,351)*y(k,44) &
                      + .650_r8*rxt(k,512)*y(k,109) + rxt(k,375)*y(k,53)
         mat(k,2115) = mat(k,2115) + rxt(k,386)*y(k,197) + .530_r8*rxt(k,391)*y(k,199) &
                      + rxt(k,398)*y(k,198) + rxt(k,370)*y(k,195)
         mat(k,1523) = rxt(k,345)*y(k,40) + .530_r8*rxt(k,393)*y(k,199) + rxt(k,399) &
                      *y(k,198) + rxt(k,376)*y(k,53)
         mat(k,1381) = mat(k,1381) + .260_r8*rxt(k,395)*y(k,199) + rxt(k,401)*y(k,198) &
                      + .300_r8*rxt(k,372)*y(k,195)
         mat(k,1617) = mat(k,1617) + .200_r8*rxt(k,387)*y(k,197) + .450_r8*rxt(k,400) &
                      *y(k,198) + .150_r8*rxt(k,371)*y(k,195)
         mat(k,788) = .080_r8*rxt(k,411)*y(k,1)
         mat(k,1057) = rxt(k,344)*y(k,185) + rxt(k,345)*y(k,8)
         mat(k,1287) = mat(k,1287) + .530_r8*rxt(k,396)*y(k,199)
         mat(k,421) = .500_r8*rxt(k,351)*y(k,185)
         mat(k,543) = rxt(k,386)*y(k,6) + .200_r8*rxt(k,387)*y(k,186)
         mat(k,857) = .060_r8*rxt(k,498)*y(k,1)
         mat(k,898) = .060_r8*rxt(k,499)*y(k,1)
         mat(k,979) = .650_r8*rxt(k,512)*y(k,185)
         mat(k,1236) = .280_r8*rxt(k,384)*y(k,1)
         mat(k,1108) = .100_r8*rxt(k,390)*y(k,1)
         mat(k,1216) = mat(k,1216) + .530_r8*rxt(k,391)*y(k,6) + .530_r8*rxt(k,393) &
                      *y(k,8) + .260_r8*rxt(k,395)*y(k,187) + .530_r8*rxt(k,396) &
                      *y(k,190)
         mat(k,1255) = rxt(k,398)*y(k,6) + rxt(k,399)*y(k,8) + rxt(k,401)*y(k,187) &
                      + .450_r8*rxt(k,400)*y(k,186) + 4.000_r8*rxt(k,403)*y(k,198)
         mat(k,1047) = rxt(k,370)*y(k,6) + .300_r8*rxt(k,372)*y(k,187) &
                      + .150_r8*rxt(k,371)*y(k,186)
         mat(k,1117) = rxt(k,375)*y(k,185) + rxt(k,376)*y(k,8)
         mat(k,420) = -(rxt(k,351)*y(k,185))
         mat(k,1783) = -rxt(k,351)*y(k,44)
         mat(k,1569) = .400_r8*rxt(k,348)*y(k,190) + .400_r8*rxt(k,400)*y(k,198)
         mat(k,1269) = .400_r8*rxt(k,348)*y(k,186)
         mat(k,1249) = .400_r8*rxt(k,400)*y(k,186)
         mat(k,355) = -(rxt(k,356)*y(k,185))
         mat(k,1774) = -rxt(k,356)*y(k,46)
         mat(k,1939) = rxt(k,347)*y(k,190)
         mat(k,1268) = rxt(k,347)*y(k,7)
         mat(k,184) = -(rxt(k,329)*y(k,183) + rxt(k,331)*y(k,185))
         mat(k,1428) = -rxt(k,329)*y(k,37)
         mat(k,1751) = -rxt(k,331)*y(k,37)
         mat(k,439) = -(rxt(k,328)*y(k,183) + rxt(k,332)*y(k,185) + rxt(k,337)*y(k,1))
         mat(k,1430) = -rxt(k,328)*y(k,36)
         mat(k,1785) = -rxt(k,332)*y(k,36)
         mat(k,1669) = -rxt(k,337)*y(k,36)
         mat(k,158) = -(rxt(k,445)*y(k,185))
         mat(k,1746) = -rxt(k,445)*y(k,56)
         mat(k,1666) = .050_r8*rxt(k,498)*y(k,105) + .050_r8*rxt(k,499)*y(k,106)
         mat(k,843) = .050_r8*rxt(k,498)*y(k,1)
         mat(k,884) = .050_r8*rxt(k,499)*y(k,1)
         mat(k,447) = -(rxt(k,406)*y(k,185))
         mat(k,1786) = -rxt(k,406)*y(k,62)
         mat(k,1943) = rxt(k,404)*y(k,198)
         mat(k,1250) = rxt(k,404)*y(k,7)
         mat(k,530) = -(rxt(k,448)*y(k,185))
         mat(k,1795) = -rxt(k,448)*y(k,63)
         mat(k,2074) = rxt(k,447)*y(k,200)
         mat(k,924) = rxt(k,447)*y(k,6)
         mat(k,796) = -(rxt(k,377)*y(k,185))
         mat(k,1816) = -rxt(k,377)*y(k,64)
         mat(k,1816) = mat(k,1816) + .500_r8*rxt(k,416)*y(k,67) + rxt(k,434)*y(k,68) &
                      + rxt(k,436)*y(k,69) + rxt(k,435)*y(k,70)
         mat(k,1495) = rxt(k,360)*y(k,47)
         mat(k,943) = rxt(k,360)*y(k,8)
         mat(k,377) = .500_r8*rxt(k,416)*y(k,185)
         mat(k,335) = rxt(k,434)*y(k,185)
         mat(k,759) = rxt(k,436)*y(k,185)
         mat(k,233) = rxt(k,435)*y(k,185)
         mat(k,1020) = -(rxt(k,407)*y(k,185))
         mat(k,1833) = -rxt(k,407)*y(k,65)
         mat(k,1833) = mat(k,1833) + .300_r8*rxt(k,415)*y(k,66) + .500_r8*rxt(k,416) &
                      *y(k,67)
         mat(k,2101) = rxt(k,381)*y(k,196) + rxt(k,392)*y(k,199)
         mat(k,469) = .300_r8*rxt(k,415)*y(k,185)
         mat(k,379) = .500_r8*rxt(k,416)*y(k,185)
         mat(k,433) = rxt(k,381)*y(k,6)
         mat(k,1208) = rxt(k,392)*y(k,6)
         mat(k,467) = -(rxt(k,415)*y(k,185))
         mat(k,1789) = -rxt(k,415)*y(k,66)
         mat(k,2071) = .080_r8*rxt(k,413)*y(k,201)
         mat(k,1143) = .080_r8*rxt(k,413)*y(k,6)
         mat(k,376) = -(rxt(k,416)*y(k,185))
         mat(k,1777) = -rxt(k,416)*y(k,67)
         mat(k,2064) = .080_r8*rxt(k,414)*y(k,202)
         mat(k,1173) = .080_r8*rxt(k,414)*y(k,6)
         mat(k,333) = -(rxt(k,434)*y(k,185))
         mat(k,1771) = -rxt(k,434)*y(k,68)
         mat(k,1560) = rxt(k,431)*y(k,79)
         mat(k,1079) = rxt(k,431)*y(k,186)
         mat(k,758) = -(rxt(k,436)*y(k,185))
         mat(k,1813) = -rxt(k,436)*y(k,69)
         mat(k,2088) = rxt(k,429)*y(k,79)
         mat(k,1493) = rxt(k,430)*y(k,79)
         mat(k,1356) = .800_r8*rxt(k,433)*y(k,79)
         mat(k,1271) = rxt(k,432)*y(k,79)
         mat(k,1082) = rxt(k,429)*y(k,6) + rxt(k,430)*y(k,8) + .800_r8*rxt(k,433) &
                      *y(k,187) + rxt(k,432)*y(k,190)
         mat(k,232) = -(rxt(k,435)*y(k,185))
         mat(k,1757) = -rxt(k,435)*y(k,70)
         mat(k,1349) = .200_r8*rxt(k,433)*y(k,79)
         mat(k,1078) = .200_r8*rxt(k,433)*y(k,187)
         mat(k,415) = -(rxt(k,517)*y(k,185))
         mat(k,1782) = -rxt(k,517)*y(k,71)
         mat(k,2069) = .200_r8*rxt(k,502)*y(k,205) + .200_r8*rxt(k,513)*y(k,217)
         mat(k,1352) = .500_r8*rxt(k,515)*y(k,217)
         mat(k,808) = .200_r8*rxt(k,502)*y(k,6)
         mat(k,961) = .200_r8*rxt(k,513)*y(k,6) + .500_r8*rxt(k,515)*y(k,187)
         mat(k,255) = -(rxt(k,518)*y(k,185))
         mat(k,1761) = -rxt(k,518)*y(k,72)
         mat(k,1551) = rxt(k,514)*y(k,217)
         mat(k,960) = rxt(k,514)*y(k,186)
         mat(k,268) = -(rxt(k,379)*y(k,185) + rxt(k,382)*y(k,8))
         mat(k,1763) = -rxt(k,379)*y(k,55)
         mat(k,1491) = -rxt(k,382)*y(k,55)
         mat(k,431) = -((rxt(k,380) + rxt(k,381)) * y(k,6))
         mat(k,2070) = -(rxt(k,380) + rxt(k,381)) * y(k,196)
         mat(k,1784) = rxt(k,379)*y(k,55)
         mat(k,269) = rxt(k,379)*y(k,185)
         mat(k,928) = -((rxt(k,446) + rxt(k,447)) * y(k,6) + rxt(k,449)*y(k,186))
         mat(k,2095) = -(rxt(k,446) + rxt(k,447)) * y(k,200)
         mat(k,1597) = -rxt(k,449)*y(k,200)
         mat(k,1826) = rxt(k,445)*y(k,56) + rxt(k,450)*y(k,74)
         mat(k,159) = rxt(k,445)*y(k,185)
         mat(k,523) = rxt(k,450)*y(k,185)
         mat(k,519) = -(rxt(k,450)*y(k,185))
         mat(k,1794) = -rxt(k,450)*y(k,74)
         mat(k,1575) = rxt(k,449)*y(k,200)
         mat(k,923) = rxt(k,449)*y(k,186)
      end do
      end subroutine nlnmat05
      subroutine nlnmat06( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,410) = -(rxt(k,385)*y(k,185))
         mat(k,1781) = -rxt(k,385)*y(k,57)
         mat(k,2068) = .800_r8*rxt(k,446)*y(k,200)
         mat(k,922) = .800_r8*rxt(k,446)*y(k,6)
         mat(k,541) = -(rxt(k,386)*y(k,6) + rxt(k,387)*y(k,186))
         mat(k,2075) = -rxt(k,386)*y(k,197)
         mat(k,1576) = -rxt(k,387)*y(k,197)
         mat(k,1796) = rxt(k,385)*y(k,57) + rxt(k,388)*y(k,58)
         mat(k,411) = rxt(k,385)*y(k,185)
         mat(k,217) = rxt(k,388)*y(k,185)
         mat(k,216) = -(rxt(k,388)*y(k,185))
         mat(k,1755) = -rxt(k,388)*y(k,58)
         mat(k,1547) = .800_r8*rxt(k,387)*y(k,197)
         mat(k,540) = .800_r8*rxt(k,387)*y(k,186)
         mat(k,109) = -(rxt(k,468)*y(k,185))
         mat(k,1737) = -rxt(k,468)*y(k,83)
         mat(k,116) = -(rxt(k,469)*y(k,185))
         mat(k,1738) = -rxt(k,469)*y(k,84)
         mat(k,1738) = mat(k,1738) + .180_r8*rxt(k,468)*y(k,83)
         mat(k,110) = .180_r8*rxt(k,468)*y(k,185)
         mat(k,607) = -(rxt(k,478)*y(k,186) + rxt(k,480)*y(k,6))
         mat(k,1580) = -rxt(k,478)*y(k,204)
         mat(k,2076) = -rxt(k,480)*y(k,204)
         mat(k,1801) = .650_r8*rxt(k,468)*y(k,83) + rxt(k,479)*y(k,85)
         mat(k,113) = .650_r8*rxt(k,468)*y(k,185)
         mat(k,564) = rxt(k,479)*y(k,185)
         mat(k,563) = -(rxt(k,479)*y(k,185))
         mat(k,1797) = -rxt(k,479)*y(k,85)
         mat(k,1577) = rxt(k,478)*y(k,204)
         mat(k,606) = rxt(k,478)*y(k,186)
         mat(k,810) = -(rxt(k,502)*y(k,6) + rxt(k,503)*y(k,186) + rxt(k,504)*y(k,187))
         mat(k,2091) = -rxt(k,502)*y(k,205)
         mat(k,1593) = -rxt(k,503)*y(k,205)
         mat(k,1359) = -rxt(k,504)*y(k,205)
         mat(k,1818) = rxt(k,505)*y(k,107) + rxt(k,496)*y(k,105) + rxt(k,497)*y(k,106)
         mat(k,363) = rxt(k,505)*y(k,185)
         mat(k,847) = rxt(k,496)*y(k,185)
         mat(k,888) = rxt(k,497)*y(k,185)
         mat(k,362) = -(rxt(k,505)*y(k,185))
         mat(k,1775) = -rxt(k,505)*y(k,107)
         mat(k,1564) = rxt(k,503)*y(k,205)
         mat(k,807) = rxt(k,503)*y(k,186)
         mat(k,1667) = .100_r8*rxt(k,498)*y(k,105) + .100_r8*rxt(k,499)*y(k,106)
         mat(k,844) = .100_r8*rxt(k,498)*y(k,1)
         mat(k,885) = .100_r8*rxt(k,499)*y(k,1)
         mat(k,803) = -(rxt(k,353)*y(k,185))
         mat(k,1817) = -rxt(k,353)*y(k,43)
         mat(k,1817) = mat(k,1817) + rxt(k,436)*y(k,69) + .200_r8*rxt(k,352)*y(k,42) &
                      + .650_r8*rxt(k,330)*y(k,131)
         mat(k,2090) = .600_r8*rxt(k,480)*y(k,204) + rxt(k,462)*y(k,206) &
                      + .700_r8*rxt(k,454)*y(k,207) + .400_r8*rxt(k,466)*y(k,211) &
                      + .170_r8*rxt(k,484)*y(k,215) + .170_r8*rxt(k,489)*y(k,209) &
                      + .340_r8*rxt(k,494)*y(k,210) + .050_r8*rxt(k,414)*y(k,202) &
                      + .250_r8*rxt(k,438)*y(k,203)
         mat(k,1496) = .050_r8*rxt(k,418)*y(k,202) + .250_r8*rxt(k,439)*y(k,203)
         mat(k,1358) = .100_r8*rxt(k,441)*y(k,203)
         mat(k,1592) = .160_r8*rxt(k,467)*y(k,211) + .070_r8*rxt(k,483)*y(k,215)
         mat(k,1273) = .250_r8*rxt(k,442)*y(k,203)
         mat(k,760) = rxt(k,436)*y(k,185)
         mat(k,610) = .600_r8*rxt(k,480)*y(k,6)
         mat(k,387) = rxt(k,462)*y(k,6)
         mat(k,322) = .700_r8*rxt(k,454)*y(k,6)
         mat(k,619) = .400_r8*rxt(k,466)*y(k,6) + .160_r8*rxt(k,467)*y(k,186)
         mat(k,985) = .170_r8*rxt(k,484)*y(k,6) + .070_r8*rxt(k,483)*y(k,186)
         mat(k,349) = .170_r8*rxt(k,489)*y(k,6)
         mat(k,632) = .340_r8*rxt(k,494)*y(k,6)
         mat(k,1178) = .050_r8*rxt(k,414)*y(k,6) + .050_r8*rxt(k,418)*y(k,8)
         mat(k,994) = .200_r8*rxt(k,352)*y(k,185)
         mat(k,1126) = .250_r8*rxt(k,438)*y(k,6) + .250_r8*rxt(k,439)*y(k,8) &
                      + .100_r8*rxt(k,441)*y(k,187) + .250_r8*rxt(k,442)*y(k,190)
         mat(k,195) = .650_r8*rxt(k,330)*y(k,185)
         mat(k,87) = -(rxt(k,452)*y(k,185))
         mat(k,1734) = -rxt(k,452)*y(k,86)
         mat(k,126) = -(rxt(k,487)*y(k,185))
         mat(k,1739) = -rxt(k,487)*y(k,101)
         mat(k,97) = -(rxt(k,453)*y(k,185))
         mat(k,1736) = -rxt(k,453)*y(k,87)
         mat(k,1736) = mat(k,1736) + .530_r8*rxt(k,452)*y(k,86)
         mat(k,89) = .530_r8*rxt(k,452)*y(k,185)
         mat(k,1735) = .120_r8*rxt(k,452)*y(k,86)
         mat(k,88) = .120_r8*rxt(k,452)*y(k,185)
         mat(k,385) = -(rxt(k,462)*y(k,6) + rxt(k,463)*y(k,186))
         mat(k,2065) = -rxt(k,462)*y(k,206)
         mat(k,1566) = -rxt(k,463)*y(k,206)
         mat(k,1778) = .350_r8*rxt(k,452)*y(k,86) + rxt(k,464)*y(k,92)
         mat(k,90) = .350_r8*rxt(k,452)*y(k,185)
         mat(k,278) = rxt(k,464)*y(k,185)
         mat(k,321) = -(rxt(k,454)*y(k,6) + rxt(k,455)*y(k,186))
         mat(k,2060) = -rxt(k,454)*y(k,207)
         mat(k,1558) = -rxt(k,455)*y(k,207)
         mat(k,1769) = .200_r8*rxt(k,469)*y(k,84) + .140_r8*rxt(k,453)*y(k,87) &
                      + rxt(k,456)*y(k,90)
         mat(k,117) = .200_r8*rxt(k,469)*y(k,185)
         mat(k,98) = .140_r8*rxt(k,453)*y(k,185)
         mat(k,212) = rxt(k,456)*y(k,185)
         mat(k,392) = -(rxt(k,457)*y(k,7) + rxt(k,458)*y(k,1))
         mat(k,1941) = -rxt(k,457)*y(k,89)
         mat(k,1668) = -rxt(k,458)*y(k,89)
         mat(k,1779) = .070_r8*rxt(k,469)*y(k,84) + .060_r8*rxt(k,453)*y(k,87) &
                      + .070_r8*rxt(k,488)*y(k,102)
         mat(k,2066) = rxt(k,459)*y(k,208)
         mat(k,118) = .070_r8*rxt(k,469)*y(k,185)
         mat(k,99) = .060_r8*rxt(k,453)*y(k,185)
         mat(k,511) = rxt(k,459)*y(k,6)
         mat(k,139) = .070_r8*rxt(k,488)*y(k,185)
         mat(k,211) = -(rxt(k,456)*y(k,185))
         mat(k,1754) = -rxt(k,456)*y(k,90)
         mat(k,1546) = rxt(k,455)*y(k,207)
         mat(k,320) = rxt(k,455)*y(k,186)
         mat(k,512) = -(rxt(k,459)*y(k,6) + rxt(k,460)*y(k,186))
         mat(k,2073) = -rxt(k,459)*y(k,208)
         mat(k,1574) = -rxt(k,460)*y(k,208)
         mat(k,1671) = rxt(k,458)*y(k,89)
         mat(k,1793) = rxt(k,461)*y(k,91)
         mat(k,2073) = mat(k,2073) + rxt(k,476)*y(k,213)
         mat(k,1574) = mat(k,1574) + .400_r8*rxt(k,477)*y(k,213)
         mat(k,393) = rxt(k,458)*y(k,1)
         mat(k,169) = rxt(k,461)*y(k,185)
         mat(k,405) = rxt(k,476)*y(k,6) + .400_r8*rxt(k,477)*y(k,186)
         mat(k,167) = -(rxt(k,461)*y(k,185))
         mat(k,1747) = -rxt(k,461)*y(k,91)
         mat(k,1544) = rxt(k,460)*y(k,208)
         mat(k,510) = rxt(k,460)*y(k,186)
         mat(k,276) = -(rxt(k,464)*y(k,185))
         mat(k,1764) = -rxt(k,464)*y(k,92)
         mat(k,1552) = rxt(k,463)*y(k,206)
         mat(k,383) = rxt(k,463)*y(k,186)
         mat(k,2061) = .200_r8*rxt(k,480)*y(k,204) + .500_r8*rxt(k,462)*y(k,206) &
                      + .060_r8*rxt(k,494)*y(k,210)
         mat(k,604) = .200_r8*rxt(k,480)*y(k,6)
         mat(k,384) = .500_r8*rxt(k,462)*y(k,6)
         mat(k,626) = .060_r8*rxt(k,494)*y(k,6)
         mat(k,2055) = .200_r8*rxt(k,480)*y(k,204) + .200_r8*rxt(k,494)*y(k,210)
         mat(k,603) = .200_r8*rxt(k,480)*y(k,6)
         mat(k,624) = .200_r8*rxt(k,494)*y(k,6)
         mat(k,2072) = .200_r8*rxt(k,480)*y(k,204) + .150_r8*rxt(k,494)*y(k,210)
         mat(k,605) = .200_r8*rxt(k,480)*y(k,6)
         mat(k,627) = .150_r8*rxt(k,494)*y(k,6)
         mat(k,2056) = .210_r8*rxt(k,494)*y(k,210)
         mat(k,625) = .210_r8*rxt(k,494)*y(k,6)
         mat(k,618) = -(rxt(k,465)*y(k,7) + rxt(k,466)*y(k,6) + rxt(k,467)*y(k,186))
         mat(k,1946) = -rxt(k,465)*y(k,211)
         mat(k,2077) = -rxt(k,466)*y(k,211)
         mat(k,1581) = -rxt(k,467)*y(k,211)
         mat(k,1740) = .100_r8*rxt(k,468)*y(k,83) + .230_r8*rxt(k,487)*y(k,101)
         mat(k,111) = .100_r8*rxt(k,468)*y(k,185)
         mat(k,127) = .230_r8*rxt(k,487)*y(k,185)
         mat(k,284) = -(rxt(k,470)*y(k,186) + rxt(k,472)*y(k,6))
         mat(k,1553) = -rxt(k,470)*y(k,212)
         mat(k,2058) = -rxt(k,472)*y(k,212)
         mat(k,1765) = .070_r8*rxt(k,468)*y(k,83) + .060_r8*rxt(k,487)*y(k,101) &
                      + rxt(k,471)*y(k,98)
         mat(k,112) = .070_r8*rxt(k,468)*y(k,185)
         mat(k,129) = .060_r8*rxt(k,487)*y(k,185)
         mat(k,228) = rxt(k,471)*y(k,185)
         mat(k,227) = -(rxt(k,471)*y(k,185))
         mat(k,1756) = -rxt(k,471)*y(k,98)
         mat(k,1548) = rxt(k,470)*y(k,212)
         mat(k,283) = rxt(k,470)*y(k,186)
         mat(k,150) = -(rxt(k,473)*y(k,185))
         mat(k,1744) = -rxt(k,473)*y(k,99)
         mat(k,2054) = rxt(k,472)*y(k,212)
         mat(k,282) = rxt(k,472)*y(k,6)
         mat(k,404) = -(rxt(k,474)*y(k,7) + rxt(k,476)*y(k,6) + rxt(k,477)*y(k,186))
         mat(k,1942) = -rxt(k,474)*y(k,213)
         mat(k,2067) = -rxt(k,476)*y(k,213)
         mat(k,1568) = -rxt(k,477)*y(k,213)
         mat(k,1780) = rxt(k,473)*y(k,99)
         mat(k,151) = rxt(k,473)*y(k,185)
         mat(k,641) = -(rxt(k,481)*y(k,186) + rxt(k,482)*y(k,6) + rxt(k,485)*y(k,7))
         mat(k,1583) = -rxt(k,481)*y(k,214)
         mat(k,2079) = -rxt(k,482)*y(k,214)
         mat(k,1947) = -rxt(k,485)*y(k,214)
         mat(k,986) = -(rxt(k,483)*y(k,186) + rxt(k,484)*y(k,6) + rxt(k,486)*y(k,7))
         mat(k,1600) = -rxt(k,483)*y(k,215)
         mat(k,2098) = -rxt(k,484)*y(k,215)
         mat(k,1953) = -rxt(k,486)*y(k,215)
         mat(k,1937) = rxt(k,474)*y(k,213)
         mat(k,403) = rxt(k,474)*y(k,7)
         mat(k,137) = -(rxt(k,488)*y(k,185))
         mat(k,1741) = -rxt(k,488)*y(k,102)
         mat(k,1741) = mat(k,1741) + .150_r8*rxt(k,487)*y(k,101)
         mat(k,128) = .150_r8*rxt(k,487)*y(k,185)
         mat(k,348) = -(rxt(k,489)*y(k,6) + rxt(k,490)*y(k,186))
         mat(k,2063) = -rxt(k,489)*y(k,209)
         mat(k,1563) = -rxt(k,490)*y(k,209)
         mat(k,1773) = .300_r8*rxt(k,488)*y(k,102) + rxt(k,491)*y(k,103)
         mat(k,138) = .300_r8*rxt(k,488)*y(k,185)
         mat(k,297) = rxt(k,491)*y(k,185)
         mat(k,296) = -(rxt(k,491)*y(k,185))
         mat(k,1766) = -rxt(k,491)*y(k,103)
         mat(k,1554) = rxt(k,490)*y(k,209)
         mat(k,347) = rxt(k,490)*y(k,186)
         mat(k,630) = -(rxt(k,492)*y(k,186) + rxt(k,494)*y(k,6))
         mat(k,1582) = -rxt(k,492)*y(k,210)
         mat(k,2078) = -rxt(k,494)*y(k,210)
         mat(k,1803) = .560_r8*rxt(k,487)*y(k,101) + rxt(k,493)*y(k,104)
         mat(k,130) = .560_r8*rxt(k,487)*y(k,185)
         mat(k,594) = rxt(k,493)*y(k,185)
         mat(k,592) = -(rxt(k,493)*y(k,185))
         mat(k,1800) = -rxt(k,493)*y(k,104)
         mat(k,1579) = rxt(k,492)*y(k,210)
         mat(k,628) = rxt(k,492)*y(k,186)
      end do
      end subroutine nlnmat06
      subroutine nlnmat07( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,848) = -(rxt(k,496)*y(k,185) + rxt(k,498)*y(k,1) + rxt(k,500)*y(k,8))
         mat(k,1821) = -rxt(k,496)*y(k,105)
         mat(k,1678) = -rxt(k,498)*y(k,105)
         mat(k,1498) = -rxt(k,500)*y(k,105)
         mat(k,889) = -(rxt(k,497)*y(k,185) + rxt(k,499)*y(k,1) + rxt(k,501)*y(k,8))
         mat(k,1823) = -rxt(k,497)*y(k,106)
         mat(k,1679) = -rxt(k,499)*y(k,106)
         mat(k,1499) = -rxt(k,501)*y(k,106)
         mat(k,13) = -(rxt(k,521)*y(k,185))
         mat(k,1726) = -rxt(k,521)*y(k,147)
         mat(k,19) = -(rxt(k,522)*y(k,185))
         mat(k,1727) = -rxt(k,522)*y(k,148)
         mat(k,908) = -(rxt(k,507)*y(k,185) + rxt(k,508)*y(k,8))
         mat(k,1824) = -rxt(k,507)*y(k,108)
         mat(k,1500) = -rxt(k,508)*y(k,108)
         mat(k,1680) = .330_r8*rxt(k,498)*y(k,105) + .330_r8*rxt(k,499)*y(k,106)
         mat(k,1824) = mat(k,1824) + rxt(k,517)*y(k,71)
         mat(k,2093) = .800_r8*rxt(k,502)*y(k,205) + .800_r8*rxt(k,513)*y(k,217)
         mat(k,1500) = mat(k,1500) + rxt(k,516)*y(k,217)
         mat(k,1361) = rxt(k,504)*y(k,205) + .500_r8*rxt(k,515)*y(k,217)
         mat(k,416) = rxt(k,517)*y(k,185)
         mat(k,811) = .800_r8*rxt(k,502)*y(k,6) + rxt(k,504)*y(k,187)
         mat(k,849) = .330_r8*rxt(k,498)*y(k,1)
         mat(k,890) = .330_r8*rxt(k,499)*y(k,1)
         mat(k,963) = .800_r8*rxt(k,513)*y(k,6) + rxt(k,516)*y(k,8) &
                      + .500_r8*rxt(k,515)*y(k,187)
         mat(k,1006) = -(rxt(k,509)*y(k,6) + rxt(k,510)*y(k,186) + rxt(k,511)*y(k,187))
         mat(k,2100) = -rxt(k,509)*y(k,216)
         mat(k,1602) = -rxt(k,510)*y(k,216)
         mat(k,1367) = -rxt(k,511)*y(k,216)
         mat(k,1832) = rxt(k,507)*y(k,108) + rxt(k,506)*y(k,110)
         mat(k,1508) = .500_r8*rxt(k,508)*y(k,108)
         mat(k,911) = rxt(k,507)*y(k,185) + .500_r8*rxt(k,508)*y(k,8)
         mat(k,491) = rxt(k,506)*y(k,185)
         mat(k,976) = -(rxt(k,512)*y(k,185))
         mat(k,1829) = -rxt(k,512)*y(k,109)
         mat(k,1685) = .300_r8*rxt(k,498)*y(k,105) + .300_r8*rxt(k,499)*y(k,106)
         mat(k,2097) = .900_r8*rxt(k,509)*y(k,216)
         mat(k,1364) = rxt(k,511)*y(k,216)
         mat(k,853) = .300_r8*rxt(k,498)*y(k,1)
         mat(k,894) = .300_r8*rxt(k,499)*y(k,1)
         mat(k,1004) = .900_r8*rxt(k,509)*y(k,6) + rxt(k,511)*y(k,187)
         mat(k,487) = -(rxt(k,506)*y(k,185))
         mat(k,1791) = -rxt(k,506)*y(k,110)
         mat(k,1572) = rxt(k,510)*y(k,216)
         mat(k,1001) = rxt(k,510)*y(k,186)
         mat(k,964) = -(rxt(k,513)*y(k,6) + rxt(k,514)*y(k,186) + rxt(k,515)*y(k,187) &
                      + rxt(k,516)*y(k,8))
         mat(k,2096) = -rxt(k,513)*y(k,217)
         mat(k,1598) = -rxt(k,514)*y(k,217)
         mat(k,1363) = -rxt(k,515)*y(k,217)
         mat(k,1504) = -rxt(k,516)*y(k,217)
         mat(k,1828) = rxt(k,518)*y(k,72)
         mat(k,1504) = mat(k,1504) + rxt(k,500)*y(k,105) + rxt(k,501)*y(k,106) &
                      + .500_r8*rxt(k,508)*y(k,108)
         mat(k,257) = rxt(k,518)*y(k,185)
         mat(k,852) = rxt(k,500)*y(k,8)
         mat(k,893) = rxt(k,501)*y(k,8)
         mat(k,909) = .500_r8*rxt(k,508)*y(k,8)
         mat(k,1155) = -(rxt(k,413)*y(k,6) + rxt(k,417)*y(k,8) + rxt(k,419)*y(k,186) &
                      + rxt(k,423)*y(k,187) + rxt(k,425)*y(k,190))
         mat(k,2110) = -rxt(k,413)*y(k,201)
         mat(k,1518) = -rxt(k,417)*y(k,201)
         mat(k,1612) = -rxt(k,419)*y(k,201)
         mat(k,1376) = -rxt(k,423)*y(k,201)
         mat(k,1282) = -rxt(k,425)*y(k,201)
         mat(k,1843) = .600_r8*rxt(k,410)*y(k,73)
         mat(k,785) = .600_r8*rxt(k,410)*y(k,185)
         mat(k,1189) = -(rxt(k,414)*y(k,6) + rxt(k,418)*y(k,8) + rxt(k,420)*y(k,186) &
                      + rxt(k,424)*y(k,187) + rxt(k,426)*y(k,190))
         mat(k,2111) = -rxt(k,414)*y(k,202)
         mat(k,1519) = -rxt(k,418)*y(k,202)
         mat(k,1613) = -rxt(k,420)*y(k,202)
         mat(k,1377) = -rxt(k,424)*y(k,202)
         mat(k,1283) = -rxt(k,426)*y(k,202)
         mat(k,1844) = .400_r8*rxt(k,410)*y(k,73)
         mat(k,786) = .400_r8*rxt(k,410)*y(k,185)
         mat(k,153) = -(rxt(k,428)*y(k,185))
         mat(k,1745) = -rxt(k,428)*y(k,77)
         mat(k,63) = -(rxt(k,422)*y(k,185))
         mat(k,1731) = -rxt(k,422)*y(k,78)
         mat(k,1731) = mat(k,1731) + .600_r8*rxt(k,421)*y(k,82)
         mat(k,573) = .600_r8*rxt(k,421)*y(k,185)
         mat(k,1234) = -(rxt(k,383)*y(k,185) + rxt(k,384)*y(k,1))
         mat(k,1846) = -rxt(k,383)*y(k,59)
         mat(k,1699) = -rxt(k,384)*y(k,59)
         mat(k,1699) = mat(k,1699) + .200_r8*rxt(k,411)*y(k,73)
         mat(k,2113) = .560_r8*rxt(k,413)*y(k,201)
         mat(k,1521) = .600_r8*rxt(k,417)*y(k,201)
         mat(k,1379) = .440_r8*rxt(k,423)*y(k,201)
         mat(k,787) = .200_r8*rxt(k,411)*y(k,1)
         mat(k,1285) = .610_r8*rxt(k,425)*y(k,201)
         mat(k,1157) = .560_r8*rxt(k,413)*y(k,6) + .600_r8*rxt(k,417)*y(k,8) &
                      + .440_r8*rxt(k,423)*y(k,187) + .610_r8*rxt(k,425)*y(k,190)
         mat(k,1104) = -(rxt(k,389)*y(k,185) + rxt(k,390)*y(k,1))
         mat(k,1840) = -rxt(k,389)*y(k,60)
         mat(k,1693) = -rxt(k,390)*y(k,60)
         mat(k,1693) = mat(k,1693) + .300_r8*rxt(k,411)*y(k,73)
         mat(k,2107) = .360_r8*rxt(k,413)*y(k,201)
         mat(k,1515) = .400_r8*rxt(k,417)*y(k,201)
         mat(k,1373) = .310_r8*rxt(k,423)*y(k,201)
         mat(k,784) = .300_r8*rxt(k,411)*y(k,1)
         mat(k,1279) = .390_r8*rxt(k,425)*y(k,201)
         mat(k,1152) = .360_r8*rxt(k,413)*y(k,6) + .400_r8*rxt(k,417)*y(k,8) &
                      + .310_r8*rxt(k,423)*y(k,187) + .390_r8*rxt(k,425)*y(k,190)
         mat(k,1214) = -((rxt(k,391) + rxt(k,392)) * y(k,6) + rxt(k,393)*y(k,8) &
                      + rxt(k,394)*y(k,186) + rxt(k,395)*y(k,187) + rxt(k,396) &
                      *y(k,190))
         mat(k,2112) = -(rxt(k,391) + rxt(k,392)) * y(k,199)
         mat(k,1520) = -rxt(k,393)*y(k,199)
         mat(k,1614) = -rxt(k,394)*y(k,199)
         mat(k,1378) = -rxt(k,395)*y(k,199)
         mat(k,1284) = -rxt(k,396)*y(k,199)
         mat(k,1845) = rxt(k,383)*y(k,59) + .500_r8*rxt(k,389)*y(k,60) &
                      + .200_r8*rxt(k,397)*y(k,61)
         mat(k,1233) = rxt(k,383)*y(k,185)
         mat(k,1106) = .500_r8*rxt(k,389)*y(k,185)
         mat(k,238) = .200_r8*rxt(k,397)*y(k,185)
         mat(k,237) = -(rxt(k,397)*y(k,185))
         mat(k,1758) = -rxt(k,397)*y(k,61)
         mat(k,1549) = rxt(k,394)*y(k,199)
         mat(k,1205) = rxt(k,394)*y(k,186)
         mat(k,1254) = -(rxt(k,398)*y(k,6) + rxt(k,399)*y(k,8) + rxt(k,400)*y(k,186) &
                      + rxt(k,401)*y(k,187) + rxt(k,402)*y(k,190) + 4._r8*rxt(k,403) &
                      *y(k,198) + rxt(k,404)*y(k,7))
         mat(k,2114) = -rxt(k,398)*y(k,198)
         mat(k,1522) = -rxt(k,399)*y(k,198)
         mat(k,1616) = -rxt(k,400)*y(k,198)
         mat(k,1380) = -rxt(k,401)*y(k,198)
         mat(k,1286) = -rxt(k,402)*y(k,198)
         mat(k,1957) = -rxt(k,404)*y(k,198)
         mat(k,1847) = .500_r8*rxt(k,389)*y(k,60) + .500_r8*rxt(k,397)*y(k,61)
         mat(k,1107) = .500_r8*rxt(k,389)*y(k,185)
         mat(k,239) = .500_r8*rxt(k,397)*y(k,185)
         mat(k,822) = -(rxt(k,339)*y(k,6) + rxt(k,340)*y(k,186) + rxt(k,341)*y(k,187) &
                      + 4._r8*rxt(k,342)*y(k,189))
         mat(k,2092) = -rxt(k,339)*y(k,189)
         mat(k,1594) = -rxt(k,340)*y(k,189)
         mat(k,1360) = -rxt(k,341)*y(k,189)
         mat(k,1819) = rxt(k,331)*y(k,37) + .500_r8*rxt(k,343)*y(k,38)
         mat(k,1436) = rxt(k,329)*y(k,37)
         mat(k,185) = rxt(k,331)*y(k,185) + rxt(k,329)*y(k,183)
         mat(k,243) = .500_r8*rxt(k,343)*y(k,185)
         mat(k,242) = -(rxt(k,343)*y(k,185))
         mat(k,1759) = -rxt(k,343)*y(k,38)
         mat(k,1550) = rxt(k,340)*y(k,189)
         mat(k,820) = rxt(k,340)*y(k,186)
         mat(k,56) = -(rxt(k,365)*y(k,185))
         mat(k,1730) = -rxt(k,365)*y(k,48)
         mat(k,724) = -(rxt(k,361)*y(k,6) + rxt(k,362)*y(k,186) + rxt(k,363)*y(k,187))
         mat(k,2086) = -rxt(k,361)*y(k,193)
         mat(k,1588) = -rxt(k,362)*y(k,193)
         mat(k,1355) = -rxt(k,363)*y(k,193)
         mat(k,1810) = rxt(k,365)*y(k,48) + rxt(k,364)*y(k,49)
         mat(k,57) = rxt(k,365)*y(k,185)
         mat(k,303) = rxt(k,364)*y(k,185)
         mat(k,302) = -(rxt(k,364)*y(k,185))
         mat(k,1767) = -rxt(k,364)*y(k,49)
         mat(k,1555) = rxt(k,362)*y(k,193)
         mat(k,723) = rxt(k,362)*y(k,186)
         mat(k,916) = -(rxt(k,369)*y(k,185))
         mat(k,1825) = -rxt(k,369)*y(k,50)
         mat(k,1681) = .520_r8*rxt(k,498)*y(k,105) + .520_r8*rxt(k,499)*y(k,106)
         mat(k,1825) = mat(k,1825) + .800_r8*rxt(k,448)*y(k,63) + .500_r8*rxt(k,512) &
                      *y(k,109)
         mat(k,2094) = .500_r8*rxt(k,380)*y(k,196) + .250_r8*rxt(k,446)*y(k,200) &
                      + .040_r8*rxt(k,502)*y(k,205) + .270_r8*rxt(k,509)*y(k,216) &
                      + .820_r8*rxt(k,361)*y(k,193)
         mat(k,1501) = .500_r8*rxt(k,382)*y(k,55)
         mat(k,1362) = .025_r8*rxt(k,504)*y(k,205) + .150_r8*rxt(k,511)*y(k,216) &
                      + .820_r8*rxt(k,363)*y(k,193)
         mat(k,533) = .800_r8*rxt(k,448)*y(k,185)
         mat(k,270) = .500_r8*rxt(k,382)*y(k,8)
         mat(k,432) = .500_r8*rxt(k,380)*y(k,6)
         mat(k,927) = .250_r8*rxt(k,446)*y(k,6)
         mat(k,812) = .040_r8*rxt(k,502)*y(k,6) + .025_r8*rxt(k,504)*y(k,187)
         mat(k,850) = .520_r8*rxt(k,498)*y(k,1)
         mat(k,891) = .520_r8*rxt(k,499)*y(k,1)
         mat(k,1003) = .270_r8*rxt(k,509)*y(k,6) + .150_r8*rxt(k,511)*y(k,187)
         mat(k,975) = .500_r8*rxt(k,512)*y(k,185)
         mat(k,725) = .820_r8*rxt(k,361)*y(k,6) + .820_r8*rxt(k,363)*y(k,187)
         mat(k,308) = -(rxt(k,373)*y(k,185))
         mat(k,1768) = -rxt(k,373)*y(k,54)
         mat(k,1556) = .850_r8*rxt(k,371)*y(k,195)
         mat(k,1042) = .850_r8*rxt(k,371)*y(k,186)
         mat(k,708) = -(rxt(k,321)*y(k,185))
         mat(k,1808) = -rxt(k,321)*y(k,14)
         mat(k,1354) = 2.000_r8*rxt(k,320)*y(k,187) + .250_r8*rxt(k,504)*y(k,205) &
                      + .250_r8*rxt(k,511)*y(k,216) + .250_r8*rxt(k,515)*y(k,217) &
                      + .250_r8*rxt(k,423)*y(k,201) + .250_r8*rxt(k,424)*y(k,202) &
                      + .250_r8*rxt(k,395)*y(k,199) + .300_r8*rxt(k,341)*y(k,189) &
                      + .500_r8*rxt(k,372)*y(k,195) + .200_r8*rxt(k,433)*y(k,79) &
                      + .300_r8*rxt(k,441)*y(k,203)
         mat(k,809) = .250_r8*rxt(k,504)*y(k,187)
         mat(k,1002) = .250_r8*rxt(k,511)*y(k,187)
         mat(k,962) = .250_r8*rxt(k,515)*y(k,187)
         mat(k,1145) = .250_r8*rxt(k,423)*y(k,187)
         mat(k,1176) = .250_r8*rxt(k,424)*y(k,187)
         mat(k,1206) = .250_r8*rxt(k,395)*y(k,187)
         mat(k,821) = .300_r8*rxt(k,341)*y(k,187)
         mat(k,1043) = .500_r8*rxt(k,372)*y(k,187)
         mat(k,1081) = .200_r8*rxt(k,433)*y(k,187)
         mat(k,1125) = .300_r8*rxt(k,441)*y(k,187)
      end do
      end subroutine nlnmat07
      subroutine nlnmat08( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,171) = -(rxt(k,354)*y(k,185))
         mat(k,1748) = -rxt(k,354)*y(k,41)
         mat(k,1348) = .200_r8*rxt(k,341)*y(k,189)
         mat(k,819) = .200_r8*rxt(k,341)*y(k,187) + .800_r8*rxt(k,342)*y(k,189)
         mat(k,995) = -(rxt(k,352)*y(k,185))
         mat(k,1831) = -rxt(k,352)*y(k,42)
         mat(k,1831) = mat(k,1831) + .700_r8*rxt(k,415)*y(k,66) + .500_r8*rxt(k,416) &
                      *y(k,67) + rxt(k,435)*y(k,70)
         mat(k,2099) = .225_r8*rxt(k,509)*y(k,216) + .050_r8*rxt(k,414)*y(k,202) &
                      + .530_r8*rxt(k,391)*y(k,199) + .250_r8*rxt(k,438)*y(k,203)
         mat(k,1507) = .050_r8*rxt(k,418)*y(k,202) + .530_r8*rxt(k,393)*y(k,199) &
                      + .250_r8*rxt(k,439)*y(k,203)
         mat(k,1366) = .125_r8*rxt(k,511)*y(k,216) + .260_r8*rxt(k,395)*y(k,199) &
                      + .100_r8*rxt(k,441)*y(k,203)
         mat(k,1275) = .530_r8*rxt(k,396)*y(k,199) + .250_r8*rxt(k,442)*y(k,203)
         mat(k,468) = .700_r8*rxt(k,415)*y(k,185)
         mat(k,378) = .500_r8*rxt(k,416)*y(k,185)
         mat(k,234) = rxt(k,435)*y(k,185)
         mat(k,1005) = .225_r8*rxt(k,509)*y(k,6) + .125_r8*rxt(k,511)*y(k,187)
         mat(k,1180) = .050_r8*rxt(k,414)*y(k,6) + .050_r8*rxt(k,418)*y(k,8)
         mat(k,1207) = .530_r8*rxt(k,391)*y(k,6) + .530_r8*rxt(k,393)*y(k,8) &
                      + .260_r8*rxt(k,395)*y(k,187) + .530_r8*rxt(k,396)*y(k,190)
         mat(k,343) = rxt(k,335)*y(k,3)
         mat(k,1127) = .250_r8*rxt(k,438)*y(k,6) + .250_r8*rxt(k,439)*y(k,8) &
                      + .100_r8*rxt(k,441)*y(k,187) + .250_r8*rxt(k,442)*y(k,190)
         mat(k,1409) = rxt(k,335)*y(k,192)
         mat(k,1036) = -(rxt(k,374)*y(k,185))
         mat(k,1835) = -rxt(k,374)*y(k,52)
         mat(k,1835) = mat(k,1835) + .500_r8*rxt(k,368)*y(k,51) + .500_r8*rxt(k,406) &
                      *y(k,62) + .700_r8*rxt(k,415)*y(k,66) + .500_r8*rxt(k,416) &
                      *y(k,67)
         mat(k,2103) = .050_r8*rxt(k,414)*y(k,202) + .220_r8*rxt(k,391)*y(k,199) &
                      + .250_r8*rxt(k,438)*y(k,203)
         mat(k,1511) = .050_r8*rxt(k,418)*y(k,202) + .220_r8*rxt(k,393)*y(k,199) &
                      + .250_r8*rxt(k,439)*y(k,203)
         mat(k,1369) = .230_r8*rxt(k,395)*y(k,199) + .200_r8*rxt(k,372)*y(k,195) &
                      + .100_r8*rxt(k,441)*y(k,203)
         mat(k,461) = .500_r8*rxt(k,368)*y(k,185)
         mat(k,1277) = .220_r8*rxt(k,396)*y(k,199) + .250_r8*rxt(k,442)*y(k,203)
         mat(k,448) = .500_r8*rxt(k,406)*y(k,185)
         mat(k,470) = .700_r8*rxt(k,415)*y(k,185)
         mat(k,380) = .500_r8*rxt(k,416)*y(k,185)
         mat(k,1183) = .050_r8*rxt(k,414)*y(k,6) + .050_r8*rxt(k,418)*y(k,8)
         mat(k,1210) = .220_r8*rxt(k,391)*y(k,6) + .220_r8*rxt(k,393)*y(k,8) &
                      + .230_r8*rxt(k,395)*y(k,187) + .220_r8*rxt(k,396)*y(k,190)
         mat(k,1044) = .200_r8*rxt(k,372)*y(k,187)
         mat(k,1129) = .250_r8*rxt(k,438)*y(k,6) + .250_r8*rxt(k,439)*y(k,8) &
                      + .100_r8*rxt(k,441)*y(k,187) + .250_r8*rxt(k,442)*y(k,190)
         mat(k,651) = -(rxt(k,333)*y(k,6) + rxt(k,334)*y(k,186))
         mat(k,2080) = -rxt(k,333)*y(k,191)
         mat(k,1584) = -rxt(k,334)*y(k,191)
         mat(k,1805) = rxt(k,332)*y(k,36)
         mat(k,440) = rxt(k,332)*y(k,185)
         mat(k,342) = -(rxt(k,335)*y(k,3))
         mat(k,1398) = -rxt(k,335)*y(k,192)
         mat(k,2062) = .750_r8*rxt(k,333)*y(k,191)
         mat(k,650) = .750_r8*rxt(k,333)*y(k,6)
         mat(k,1543) = rxt(k,334)*y(k,191)
         mat(k,649) = rxt(k,334)*y(k,186)
         mat(k,247) = -(rxt(k,437)*y(k,185))
         mat(k,1760) = -rxt(k,437)*y(k,76)
         mat(k,2057) = .870_r8*rxt(k,414)*y(k,202)
         mat(k,1490) = .950_r8*rxt(k,418)*y(k,202)
         mat(k,1350) = .750_r8*rxt(k,424)*y(k,202)
         mat(k,1267) = rxt(k,426)*y(k,202)
         mat(k,1172) = .870_r8*rxt(k,414)*y(k,6) + .950_r8*rxt(k,418)*y(k,8) &
                      + .750_r8*rxt(k,424)*y(k,187) + rxt(k,426)*y(k,190)
         mat(k,1045) = -(rxt(k,370)*y(k,6) + rxt(k,371)*y(k,186) + rxt(k,372)*y(k,187))
         mat(k,2104) = -rxt(k,370)*y(k,195)
         mat(k,1605) = -rxt(k,371)*y(k,195)
         mat(k,1370) = -rxt(k,372)*y(k,195)
         mat(k,1690) = .060_r8*rxt(k,498)*y(k,105) + .060_r8*rxt(k,499)*y(k,106)
         mat(k,1836) = .150_r8*rxt(k,512)*y(k,109) + rxt(k,369)*y(k,50) + rxt(k,373) &
                      *y(k,54)
         mat(k,855) = .060_r8*rxt(k,498)*y(k,1)
         mat(k,896) = .060_r8*rxt(k,499)*y(k,1)
         mat(k,978) = .150_r8*rxt(k,512)*y(k,185)
         mat(k,917) = rxt(k,369)*y(k,185)
         mat(k,309) = rxt(k,373)*y(k,185)
         mat(k,1116) = -(rxt(k,375)*y(k,185) + rxt(k,376)*y(k,8))
         mat(k,1841) = -rxt(k,375)*y(k,53)
         mat(k,1516) = -rxt(k,376)*y(k,53)
         mat(k,1694) = .500_r8*rxt(k,384)*y(k,59) + .880_r8*rxt(k,390)*y(k,60)
         mat(k,1841) = mat(k,1841) + rxt(k,377)*y(k,64) + rxt(k,374)*y(k,52)
         mat(k,2108) = .400_r8*rxt(k,480)*y(k,204) + .170_r8*rxt(k,482)*y(k,214) &
                      + .170_r8*rxt(k,484)*y(k,215) + .510_r8*rxt(k,489)*y(k,209) &
                      + .540_r8*rxt(k,494)*y(k,210) + .050_r8*rxt(k,414)*y(k,202) &
                      + .250_r8*rxt(k,391)*y(k,199) + .250_r8*rxt(k,438)*y(k,203)
         mat(k,1516) = mat(k,1516) + .050_r8*rxt(k,418)*y(k,202) + .250_r8*rxt(k,393) &
                      *y(k,199) + .250_r8*rxt(k,439)*y(k,203)
         mat(k,1374) = .240_r8*rxt(k,395)*y(k,199) + .500_r8*rxt(k,372)*y(k,195) &
                      + .100_r8*rxt(k,441)*y(k,203)
         mat(k,1610) = .070_r8*rxt(k,481)*y(k,214) + .070_r8*rxt(k,483)*y(k,215)
         mat(k,1280) = .250_r8*rxt(k,396)*y(k,199) + .250_r8*rxt(k,442)*y(k,203)
         mat(k,797) = rxt(k,377)*y(k,185)
         mat(k,613) = .400_r8*rxt(k,480)*y(k,6)
         mat(k,643) = .170_r8*rxt(k,482)*y(k,6) + .070_r8*rxt(k,481)*y(k,186)
         mat(k,988) = .170_r8*rxt(k,484)*y(k,6) + .070_r8*rxt(k,483)*y(k,186)
         mat(k,350) = .510_r8*rxt(k,489)*y(k,6)
         mat(k,635) = .540_r8*rxt(k,494)*y(k,6)
         mat(k,1187) = .050_r8*rxt(k,414)*y(k,6) + .050_r8*rxt(k,418)*y(k,8)
         mat(k,1232) = .500_r8*rxt(k,384)*y(k,1)
         mat(k,1105) = .880_r8*rxt(k,390)*y(k,1)
         mat(k,1213) = .250_r8*rxt(k,391)*y(k,6) + .250_r8*rxt(k,393)*y(k,8) &
                      + .240_r8*rxt(k,395)*y(k,187) + .250_r8*rxt(k,396)*y(k,190)
         mat(k,1037) = rxt(k,374)*y(k,185)
         mat(k,1046) = .500_r8*rxt(k,372)*y(k,187)
         mat(k,1130) = .250_r8*rxt(k,438)*y(k,6) + .250_r8*rxt(k,439)*y(k,8) &
                      + .100_r8*rxt(k,441)*y(k,187) + .250_r8*rxt(k,442)*y(k,190)
         mat(k,1088) = -(rxt(k,429)*y(k,6) + rxt(k,430)*y(k,8) + rxt(k,431)*y(k,186) &
                      + rxt(k,432)*y(k,190) + rxt(k,433)*y(k,187))
         mat(k,2106) = -rxt(k,429)*y(k,79)
         mat(k,1514) = -rxt(k,430)*y(k,79)
         mat(k,1608) = -rxt(k,431)*y(k,79)
         mat(k,1278) = -rxt(k,432)*y(k,79)
         mat(k,1372) = -rxt(k,433)*y(k,79)
         mat(k,1514) = mat(k,1514) + rxt(k,412)*y(k,73)
         mat(k,783) = rxt(k,412)*y(k,8)
         mat(k,1732) = rxt(k,407)*y(k,65)
         mat(k,2053) = .100_r8*rxt(k,509)*y(k,216)
         mat(k,1017) = rxt(k,407)*y(k,185)
         mat(k,1000) = .100_r8*rxt(k,509)*y(k,6)
         mat(k,1131) = -(rxt(k,438)*y(k,6) + rxt(k,439)*y(k,8) + rxt(k,440)*y(k,186) &
                      + rxt(k,441)*y(k,187) + rxt(k,442)*y(k,190))
         mat(k,2109) = -rxt(k,438)*y(k,203)
         mat(k,1517) = -rxt(k,439)*y(k,203)
         mat(k,1611) = -rxt(k,440)*y(k,203)
         mat(k,1375) = -rxt(k,441)*y(k,203)
         mat(k,1281) = -rxt(k,442)*y(k,203)
         mat(k,1842) = rxt(k,428)*y(k,77) + rxt(k,422)*y(k,78) + rxt(k,437)*y(k,76) &
                      + rxt(k,443)*y(k,81) + .400_r8*rxt(k,421)*y(k,82)
         mat(k,155) = rxt(k,428)*y(k,185)
         mat(k,64) = rxt(k,422)*y(k,185)
         mat(k,248) = rxt(k,437)*y(k,185)
         mat(k,176) = rxt(k,443)*y(k,185)
         mat(k,576) = .400_r8*rxt(k,421)*y(k,185)
         mat(k,175) = -((rxt(k,443) + rxt(k,444)) * y(k,185))
         mat(k,1749) = -(rxt(k,443) + rxt(k,444)) * y(k,81)
         mat(k,1545) = rxt(k,440)*y(k,203)
         mat(k,1124) = rxt(k,440)*y(k,186)
         mat(k,574) = -(rxt(k,421)*y(k,185))
         mat(k,1798) = -rxt(k,421)*y(k,82)
         mat(k,1578) = rxt(k,419)*y(k,201) + rxt(k,420)*y(k,202)
         mat(k,1144) = rxt(k,419)*y(k,186)
         mat(k,1175) = rxt(k,420)*y(k,186)
         mat(k,193) = -(rxt(k,327)*y(k,183) + rxt(k,330)*y(k,185))
         mat(k,1429) = -rxt(k,327)*y(k,131)
         mat(k,1752) = -rxt(k,330)*y(k,131)
         mat(k,704) = -(rxt(k,323)*y(k,185))
         mat(k,1807) = -rxt(k,323)*y(k,132)
         mat(k,1674) = .120_r8*rxt(k,359)*y(k,47) + .110_r8*rxt(k,411)*y(k,73) &
                      + .370_r8*rxt(k,337)*y(k,36) + .050_r8*rxt(k,498)*y(k,105) &
                      + .050_r8*rxt(k,499)*y(k,106) + .120_r8*rxt(k,384)*y(k,59) &
                      + .330_r8*rxt(k,390)*y(k,60)
         mat(k,1807) = mat(k,1807) + .350_r8*rxt(k,330)*y(k,131)
         mat(k,2084) = rxt(k,325)*y(k,188)
         mat(k,1586) = rxt(k,326)*y(k,188)
         mat(k,941) = .120_r8*rxt(k,359)*y(k,1)
         mat(k,779) = .110_r8*rxt(k,411)*y(k,1)
         mat(k,441) = .370_r8*rxt(k,337)*y(k,1)
         mat(k,845) = .050_r8*rxt(k,498)*y(k,1)
         mat(k,886) = .050_r8*rxt(k,499)*y(k,1)
         mat(k,1227) = .120_r8*rxt(k,384)*y(k,1)
         mat(k,1102) = .330_r8*rxt(k,390)*y(k,1)
         mat(k,194) = .350_r8*rxt(k,330)*y(k,185)
         mat(k,315) = rxt(k,325)*y(k,6) + rxt(k,326)*y(k,186)
         mat(k,314) = -(rxt(k,325)*y(k,6) + rxt(k,326)*y(k,186))
         mat(k,2059) = -rxt(k,325)*y(k,188)
         mat(k,1557) = -rxt(k,326)*y(k,188)
         mat(k,1637) = rxt(k,316)*y(k,186)
         mat(k,1557) = mat(k,1557) + rxt(k,316)*y(k,15)
         mat(k,832) = -(rxt(k,541)*y(k,185))
         mat(k,1820) = -rxt(k,541)*y(k,137)
         mat(k,1677) = rxt(k,536)*y(k,144)
         mat(k,1820) = mat(k,1820) + (rxt(k,543)+.500_r8*rxt(k,544))*y(k,138) &
                      + rxt(k,530)*y(k,142) + rxt(k,534)*y(k,144)
         mat(k,1951) = rxt(k,537)*y(k,144)
         mat(k,1497) = rxt(k,545)*y(k,138)
         mat(k,1874) = rxt(k,538)*y(k,144)
         mat(k,251) = rxt(k,540)*y(k,144)
         mat(k,1466) = rxt(k,539)*y(k,144)
         mat(k,200) = (rxt(k,543)+.500_r8*rxt(k,544))*y(k,185) + rxt(k,545)*y(k,8)
         mat(k,262) = rxt(k,530)*y(k,185)
         mat(k,1302) = rxt(k,536)*y(k,1) + rxt(k,534)*y(k,185) + rxt(k,537)*y(k,7) &
                      + rxt(k,538)*y(k,25) + rxt(k,540)*y(k,26) + rxt(k,539)*y(k,32) &
                      + rxt(k,535)*y(k,3)
         mat(k,1408) = rxt(k,535)*y(k,144)
         mat(k,199) = -((rxt(k,543) + rxt(k,544)) * y(k,185) + rxt(k,545)*y(k,8))
         mat(k,1753) = -(rxt(k,543) + rxt(k,544)) * y(k,138)
         mat(k,1488) = -rxt(k,545)*y(k,138)
         mat(k,260) = -(rxt(k,529)*y(k,2) + rxt(k,530)*y(k,185))
         mat(k,1898) = -rxt(k,529)*y(k,142)
         mat(k,1762) = -rxt(k,530)*y(k,142)
         mat(k,581) = -(rxt(k,531)*y(k,185) + rxt(k,532)*y(k,3) + rxt(k,533)*y(k,1))
         mat(k,1799) = -rxt(k,531)*y(k,143)
         mat(k,1403) = -rxt(k,532)*y(k,143)
         mat(k,1672) = -rxt(k,533)*y(k,143)
         mat(k,1303) = -(rxt(k,534)*y(k,185) + rxt(k,535)*y(k,3) + rxt(k,536)*y(k,1) &
                      + rxt(k,537)*y(k,7) + rxt(k,538)*y(k,25) + rxt(k,539)*y(k,32) &
                      + rxt(k,540)*y(k,26))
         mat(k,1849) = -rxt(k,534)*y(k,144)
         mat(k,1411) = -rxt(k,535)*y(k,144)
         mat(k,1702) = -rxt(k,536)*y(k,144)
         mat(k,1959) = -rxt(k,537)*y(k,144)
         mat(k,1876) = -rxt(k,538)*y(k,144)
         mat(k,1467) = -rxt(k,539)*y(k,144)
         mat(k,252) = -rxt(k,540)*y(k,144)
         mat(k,1702) = mat(k,1702) + rxt(k,533)*y(k,143)
         mat(k,1917) = rxt(k,529)*y(k,142)
         mat(k,1849) = mat(k,1849) + rxt(k,531)*y(k,143)
         mat(k,264) = rxt(k,529)*y(k,2)
         mat(k,582) = rxt(k,533)*y(k,1) + rxt(k,531)*y(k,185) + rxt(k,532)*y(k,3)
         mat(k,1411) = mat(k,1411) + rxt(k,532)*y(k,143)
      end do
      end subroutine nlnmat08
      subroutine nlnmat09( ofl, ofu, chnkpnts, mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k,180) = -(rxt(k,542)*y(k,228))
         mat(k,1981) = -rxt(k,542)*y(k,145)
         mat(k,1750) = rxt(k,541)*y(k,137)
         mat(k,831) = rxt(k,541)*y(k,185)
         mat(k,1979) = rxt(k,542)*y(k,145)
         mat(k,179) = rxt(k,542)*y(k,228)
         mat(k,47) = -(rxt(k,527)*y(k,185))
         mat(k,1728) = -rxt(k,527)*y(k,139)
         mat(k,1721) = .0031005_r8*rxt(k,410)*y(k,73) + .1364005_r8*rxt(k,468)*y(k,83) &
                      + .0023005_r8*rxt(k,452)*y(k,86) + .1677005_r8*rxt(k,487) &
                      *y(k,101) + .0508005_r8*rxt(k,496)*y(k,105) &
                      + .2202005_r8*rxt(k,497)*y(k,106) + .2381005_r8*rxt(k,521) &
                      *y(k,147) + .5931005_r8*rxt(k,522)*y(k,148)
         mat(k,774) = .0031005_r8*rxt(k,410)*y(k,185)
         mat(k,104) = .1364005_r8*rxt(k,468)*y(k,185)
         mat(k,82) = .0023005_r8*rxt(k,452)*y(k,185)
         mat(k,121) = .1677005_r8*rxt(k,487)*y(k,185)
         mat(k,838) = .0508005_r8*rxt(k,496)*y(k,185)
         mat(k,879) = .2202005_r8*rxt(k,497)*y(k,185)
         mat(k,8) = .2381005_r8*rxt(k,521)*y(k,185)
         mat(k,14) = .5931005_r8*rxt(k,522)*y(k,185)
         mat(k,1722) = .0035005_r8*rxt(k,410)*y(k,73) + .0101005_r8*rxt(k,468)*y(k,83) &
                      + .0008005_r8*rxt(k,452)*y(k,86) + .0174005_r8*rxt(k,487) &
                      *y(k,101) + .1149005_r8*rxt(k,496)*y(k,105) &
                      + .2067005_r8*rxt(k,497)*y(k,106) + .1308005_r8*rxt(k,521) &
                      *y(k,147) + .1534005_r8*rxt(k,522)*y(k,148)
         mat(k,775) = .0035005_r8*rxt(k,410)*y(k,185)
         mat(k,105) = .0101005_r8*rxt(k,468)*y(k,185)
         mat(k,83) = .0008005_r8*rxt(k,452)*y(k,185)
         mat(k,122) = .0174005_r8*rxt(k,487)*y(k,185)
         mat(k,839) = .1149005_r8*rxt(k,496)*y(k,185)
         mat(k,880) = .2067005_r8*rxt(k,497)*y(k,185)
         mat(k,9) = .1308005_r8*rxt(k,521)*y(k,185)
         mat(k,15) = .1534005_r8*rxt(k,522)*y(k,185)
         mat(k,1661) = .0348005_r8*rxt(k,498)*y(k,105) + .0653005_r8*rxt(k,499) &
                      *y(k,106)
         mat(k,1723) = .0003005_r8*rxt(k,410)*y(k,73) + .0763005_r8*rxt(k,468)*y(k,83) &
                      + .0843005_r8*rxt(k,452)*y(k,86) + .086_r8*rxt(k,487)*y(k,101) &
                      + .0348005_r8*rxt(k,496)*y(k,105) + .0653005_r8*rxt(k,497) &
                      *y(k,106) + .0348005_r8*rxt(k,521)*y(k,147) &
                      + .0459005_r8*rxt(k,522)*y(k,148)
         mat(k,776) = .0003005_r8*rxt(k,410)*y(k,185)
         mat(k,106) = .0763005_r8*rxt(k,468)*y(k,185)
         mat(k,84) = .0843005_r8*rxt(k,452)*y(k,185)
         mat(k,123) = .086_r8*rxt(k,487)*y(k,185)
         mat(k,840) = .0348005_r8*rxt(k,498)*y(k,1) + .0348005_r8*rxt(k,496)*y(k,185)
         mat(k,881) = .0653005_r8*rxt(k,499)*y(k,1) + .0653005_r8*rxt(k,497)*y(k,185)
         mat(k,10) = .0348005_r8*rxt(k,521)*y(k,185)
         mat(k,16) = .0459005_r8*rxt(k,522)*y(k,185)
         mat(k,1662) = .0271005_r8*rxt(k,411)*y(k,73) + .0554005_r8*rxt(k,498) &
                      *y(k,105) + .1284005_r8*rxt(k,499)*y(k,106)
         mat(k,1724) = .0271005_r8*rxt(k,410)*y(k,73) + .2157005_r8*rxt(k,468)*y(k,83) &
                      + .0443005_r8*rxt(k,452)*y(k,86) + .0512005_r8*rxt(k,487) &
                      *y(k,101) + .0554005_r8*rxt(k,496)*y(k,105) &
                      + .1284005_r8*rxt(k,497)*y(k,106) + .0076005_r8*rxt(k,521) &
                      *y(k,147) + .0085005_r8*rxt(k,522)*y(k,148)
         mat(k,1486) = .0590245_r8*rxt(k,412)*y(k,73) + .1749305_r8*rxt(k,500) &
                      *y(k,105) + .1749305_r8*rxt(k,501)*y(k,106)
         mat(k,777) = .0271005_r8*rxt(k,411)*y(k,1) + .0271005_r8*rxt(k,410)*y(k,185) &
                      + .0590245_r8*rxt(k,412)*y(k,8)
         mat(k,107) = .2157005_r8*rxt(k,468)*y(k,185)
         mat(k,85) = .0443005_r8*rxt(k,452)*y(k,185)
         mat(k,124) = .0512005_r8*rxt(k,487)*y(k,185)
         mat(k,841) = .0554005_r8*rxt(k,498)*y(k,1) + .0554005_r8*rxt(k,496)*y(k,185) &
                      + .1749305_r8*rxt(k,500)*y(k,8)
         mat(k,882) = .1284005_r8*rxt(k,499)*y(k,1) + .1284005_r8*rxt(k,497)*y(k,185) &
                      + .1749305_r8*rxt(k,501)*y(k,8)
         mat(k,11) = .0076005_r8*rxt(k,521)*y(k,185)
         mat(k,17) = .0085005_r8*rxt(k,522)*y(k,185)
         mat(k,1663) = .1278005_r8*rxt(k,498)*y(k,105) + .114_r8*rxt(k,499)*y(k,106)
         mat(k,1725) = .0474005_r8*rxt(k,410)*y(k,73) + .0232005_r8*rxt(k,468)*y(k,83) &
                      + .1621005_r8*rxt(k,452)*y(k,86) + .1598005_r8*rxt(k,487) &
                      *y(k,101) + .1278005_r8*rxt(k,496)*y(k,105) + .114_r8*rxt(k,497) &
                      *y(k,106) + .0113005_r8*rxt(k,521)*y(k,147) &
                      + .0128005_r8*rxt(k,522)*y(k,148)
         mat(k,1487) = .0250245_r8*rxt(k,412)*y(k,73) + .5901905_r8*rxt(k,500) &
                      *y(k,105) + .5901905_r8*rxt(k,501)*y(k,106)
         mat(k,778) = .0474005_r8*rxt(k,410)*y(k,185) + .0250245_r8*rxt(k,412)*y(k,8)
         mat(k,108) = .0232005_r8*rxt(k,468)*y(k,185)
         mat(k,86) = .1621005_r8*rxt(k,452)*y(k,185)
         mat(k,125) = .1598005_r8*rxt(k,487)*y(k,185)
         mat(k,842) = .1278005_r8*rxt(k,498)*y(k,1) + .1278005_r8*rxt(k,496)*y(k,185) &
                      + .5901905_r8*rxt(k,500)*y(k,8)
         mat(k,883) = .114_r8*rxt(k,499)*y(k,1) + .114_r8*rxt(k,497)*y(k,185) &
                      + .5901905_r8*rxt(k,501)*y(k,8)
         mat(k,12) = .0113005_r8*rxt(k,521)*y(k,185)
         mat(k,18) = .0128005_r8*rxt(k,522)*y(k,185)
         mat(k,1413) = -(rxt(k,153)*y(k,2) + rxt(k,163)*y(k,225) + rxt(k,168)*y(k,227) &
                      + rxt(k,202)*y(k,20) + rxt(k,224)*y(k,226) + rxt(k,227)*y(k,5) &
                      + rxt(k,335)*y(k,192) + rxt(k,532)*y(k,143) + rxt(k,535) &
                      *y(k,144) + rxt(k,565)*y(k,218) + (rxt(k,572) + rxt(k,573) &
                      ) * y(k,221) + rxt(k,575)*y(k,222))
         mat(k,1921) = -rxt(k,153)*y(k,3)
         mat(k,78) = -rxt(k,163)*y(k,3)
         mat(k,2015) = -rxt(k,168)*y(k,3)
         mat(k,2138) = -rxt(k,202)*y(k,3)
         mat(k,428) = -rxt(k,224)*y(k,3)
         mat(k,742) = -rxt(k,227)*y(k,3)
         mat(k,344) = -rxt(k,335)*y(k,3)
         mat(k,583) = -rxt(k,532)*y(k,3)
         mat(k,1305) = -rxt(k,535)*y(k,3)
         mat(k,508) = -rxt(k,565)*y(k,3)
         mat(k,401) = -(rxt(k,572) + rxt(k,573)) * y(k,3)
         mat(k,556) = -rxt(k,575)*y(k,3)
         mat(k,1705) = 2.000_r8*rxt(k,154)*y(k,2) + 2.000_r8*rxt(k,173)*y(k,227) &
                      + rxt(k,208)*y(k,185) + rxt(k,234)*y(k,6) + rxt(k,237)*y(k,7) &
                      + rxt(k,203)*y(k,20) + 2.000_r8*rxt(k,216)*y(k,186) + rxt(k,251) &
                      *y(k,183) + rxt(k,278)*y(k,184) + rxt(k,533)*y(k,143) &
                      + rxt(k,536)*y(k,144)
         mat(k,1921) = mat(k,1921) + 2.000_r8*rxt(k,154)*y(k,1) + 2.000_r8*rxt(k,155) &
                      *y(k,2) + rxt(k,207)*y(k,185) + rxt(k,235)*y(k,7) + rxt(k,243) &
                      *y(k,8) + rxt(k,215)*y(k,186) + rxt(k,258)*y(k,25) + rxt(k,281) &
                      *y(k,32) + rxt(k,162)*y(k,225)
         mat(k,2015) = mat(k,2015) + 2.000_r8*rxt(k,173)*y(k,1)
         mat(k,1853) = rxt(k,208)*y(k,1) + rxt(k,207)*y(k,2) + rxt(k,247)*y(k,10) &
                      + rxt(k,209)*y(k,186) + rxt(k,260)*y(k,25)
         mat(k,742) = mat(k,742) + rxt(k,231)*y(k,7)
         mat(k,2118) = rxt(k,234)*y(k,1) + rxt(k,571)*y(k,219)
         mat(k,1963) = rxt(k,237)*y(k,1) + rxt(k,235)*y(k,2) + rxt(k,231)*y(k,5)
         mat(k,1527) = rxt(k,243)*y(k,2) + rxt(k,245)*y(k,186)
         mat(k,370) = rxt(k,247)*y(k,185)
         mat(k,1383) = rxt(k,318)*y(k,186)
         mat(k,2138) = mat(k,2138) + rxt(k,203)*y(k,1) + rxt(k,205)*y(k,186)
         mat(k,1621) = 2.000_r8*rxt(k,216)*y(k,1) + rxt(k,215)*y(k,2) + rxt(k,209) &
                      *y(k,185) + rxt(k,245)*y(k,8) + rxt(k,318)*y(k,187) + rxt(k,205) &
                      *y(k,20) + 2.000_r8*rxt(k,217)*y(k,186) + rxt(k,254)*y(k,183) &
                      + rxt(k,261)*y(k,25) + rxt(k,279)*y(k,184) + rxt(k,283)*y(k,32) &
                      + rxt(k,367)*y(k,194) + rxt(k,340)*y(k,189) + rxt(k,362) &
                      *y(k,193)
         mat(k,1446) = rxt(k,251)*y(k,1) + rxt(k,254)*y(k,186)
         mat(k,1880) = rxt(k,258)*y(k,2) + rxt(k,260)*y(k,185) + rxt(k,261)*y(k,186) + ( &
                      + 2.000_r8*rxt(k,265)+2.000_r8*rxt(k,266))*y(k,25) + (rxt(k,287) &
                       +rxt(k,288))*y(k,32)
         mat(k,1338) = rxt(k,278)*y(k,1) + rxt(k,279)*y(k,186)
         mat(k,1470) = rxt(k,281)*y(k,2) + rxt(k,283)*y(k,186) + (rxt(k,287) &
                       +rxt(k,288))*y(k,25) + 2.000_r8*rxt(k,289)*y(k,32)
         mat(k,716) = rxt(k,367)*y(k,186)
         mat(k,825) = rxt(k,340)*y(k,186)
         mat(k,728) = rxt(k,362)*y(k,186)
         mat(k,583) = mat(k,583) + rxt(k,533)*y(k,1)
         mat(k,1305) = mat(k,1305) + rxt(k,536)*y(k,1)
         mat(k,1413) = mat(k,1413) + 2.000_r8*rxt(k,163)*y(k,225)
         mat(k,78) = mat(k,78) + rxt(k,162)*y(k,2) + 2.000_r8*rxt(k,163)*y(k,3)
         mat(k,681) = rxt(k,571)*y(k,6)
         mat(k,80) = -(rxt(k,156)*y(k,2) + rxt(k,157)*y(k,3) + rxt(k,159)*y(k,1))
         mat(k,1897) = -rxt(k,156)*y(k,224)
         mat(k,1397) = -rxt(k,157)*y(k,224)
         mat(k,1665) = -rxt(k,159)*y(k,224)
         mat(k,2006) = rxt(k,168)*y(k,3)
         mat(k,1397) = mat(k,1397) + rxt(k,168)*y(k,227)
         mat(k,77) = -(rxt(k,162)*y(k,2) + rxt(k,163)*y(k,3))
         mat(k,1896) = -rxt(k,162)*y(k,225)
         mat(k,1396) = -rxt(k,163)*y(k,225)
         mat(k,1664) = rxt(k,159)*y(k,224)
         mat(k,1896) = mat(k,1896) + rxt(k,156)*y(k,224)
         mat(k,1396) = mat(k,1396) + rxt(k,157)*y(k,224)
         mat(k,79) = rxt(k,159)*y(k,1) + rxt(k,156)*y(k,2) + rxt(k,157)*y(k,3)
         mat(k,550) = -((rxt(k,567) + rxt(k,568)) * y(k,2) + rxt(k,575)*y(k,3) &
                      + rxt(k,579)*y(k,223))
         mat(k,1904) = -(rxt(k,567) + rxt(k,568)) * y(k,222)
         mat(k,1402) = -rxt(k,575)*y(k,222)
         mat(k,694) = -rxt(k,579)*y(k,222)
         mat(k,677) = -(rxt(k,570)*y(k,5) + rxt(k,571)*y(k,6) + rxt(k,578)*y(k,223))
         mat(k,737) = -rxt(k,570)*y(k,219)
         mat(k,2081) = -rxt(k,571)*y(k,219)
         mat(k,695) = -rxt(k,578)*y(k,219)
         mat(k,1404) = rxt(k,575)*y(k,222) + rxt(k,572)*y(k,221) + rxt(k,565)*y(k,218)
         mat(k,551) = rxt(k,575)*y(k,3)
         mat(k,398) = rxt(k,572)*y(k,3)
         mat(k,504) = rxt(k,565)*y(k,3)
         mat(k,396) = -((rxt(k,572) + rxt(k,573)) * y(k,3) + rxt(k,574)*y(k,2))
         mat(k,1399) = -(rxt(k,572) + rxt(k,573)) * y(k,221)
         mat(k,1899) = -rxt(k,574)*y(k,221)
         mat(k,503) = -(rxt(k,565)*y(k,3))
         mat(k,1401) = -rxt(k,565)*y(k,218)
         mat(k,1903) = rxt(k,568)*y(k,222) + rxt(k,574)*y(k,221)
         mat(k,549) = rxt(k,568)*y(k,2)
         mat(k,397) = rxt(k,574)*y(k,2)
         mat(k,686) = -(rxt(k,577)*y(k,223))
         mat(k,696) = -rxt(k,577)*y(k,220)
         mat(k,1908) = rxt(k,567)*y(k,222)
         mat(k,738) = rxt(k,570)*y(k,219)
         mat(k,2082) = rxt(k,571)*y(k,219)
         mat(k,1405) = rxt(k,573)*y(k,221)
         mat(k,552) = rxt(k,567)*y(k,2)
         mat(k,678) = rxt(k,570)*y(k,5) + rxt(k,571)*y(k,6)
         mat(k,399) = rxt(k,573)*y(k,3)
         mat(k,426) = -(rxt(k,224)*y(k,3) + rxt(k,225)*y(k,2))
         mat(k,1400) = -rxt(k,224)*y(k,226)
         mat(k,1900) = -rxt(k,225)*y(k,226)
         mat(k,1900) = mat(k,1900) + rxt(k,567)*y(k,222)
         mat(k,548) = rxt(k,567)*y(k,2) + .900_r8*rxt(k,579)*y(k,223)
         mat(k,685) = .800_r8*rxt(k,577)*y(k,223)
         mat(k,693) = .900_r8*rxt(k,579)*y(k,222) + .800_r8*rxt(k,577)*y(k,220)
         mat(k,697) = -(rxt(k,577)*y(k,220) + rxt(k,578)*y(k,219) + rxt(k,579) &
                      *y(k,222))
         mat(k,687) = -rxt(k,577)*y(k,223)
         mat(k,679) = -rxt(k,578)*y(k,223)
         mat(k,553) = -rxt(k,579)*y(k,223)
      end do
      end subroutine nlnmat09
      subroutine nlnmat_finit( ofl, ofu, chnkpnts, mat, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(chnkpnts,nzcnt)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
      integer :: k
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
      do k = ofl,ofu
         mat(k, 1) = lmat(k, 1)
         mat(k, 2) = lmat(k, 2)
         mat(k, 3) = lmat(k, 3)
         mat(k, 4) = lmat(k, 4)
         mat(k, 5) = lmat(k, 5)
         mat(k, 6) = lmat(k, 6)
         mat(k, 7) = lmat(k, 7)
         mat(k, 13) = mat(k, 13) + lmat(k, 13)
         mat(k, 19) = mat(k, 19) + lmat(k, 19)
         mat(k, 20) = lmat(k, 20)
         mat(k, 21) = lmat(k, 21)
         mat(k, 22) = lmat(k, 22)
         mat(k, 23) = lmat(k, 23)
         mat(k, 24) = lmat(k, 24)
         mat(k, 25) = lmat(k, 25)
         mat(k, 26) = lmat(k, 26)
         mat(k, 27) = lmat(k, 27)
         mat(k, 28) = lmat(k, 28)
         mat(k, 29) = lmat(k, 29)
         mat(k, 30) = lmat(k, 30)
         mat(k, 31) = lmat(k, 31)
         mat(k, 32) = lmat(k, 32)
         mat(k, 33) = lmat(k, 33)
         mat(k, 34) = lmat(k, 34)
         mat(k, 35) = lmat(k, 35)
         mat(k, 36) = lmat(k, 36)
         mat(k, 37) = lmat(k, 37)
         mat(k, 38) = lmat(k, 38)
         mat(k, 39) = lmat(k, 39)
         mat(k, 40) = lmat(k, 40)
         mat(k, 41) = lmat(k, 41)
         mat(k, 42) = lmat(k, 42)
         mat(k, 43) = lmat(k, 43)
         mat(k, 44) = lmat(k, 44)
         mat(k, 45) = lmat(k, 45)
         mat(k, 46) = lmat(k, 46)
         mat(k, 47) = mat(k, 47) + lmat(k, 47)
         mat(k, 50) = mat(k, 50) + lmat(k, 50)
         mat(k, 53) = lmat(k, 53)
         mat(k, 54) = lmat(k, 54)
         mat(k, 55) = lmat(k, 55)
         mat(k, 56) = mat(k, 56) + lmat(k, 56)
         mat(k, 60) = lmat(k, 60)
         mat(k, 61) = lmat(k, 61)
         mat(k, 62) = lmat(k, 62)
         mat(k, 63) = mat(k, 63) + lmat(k, 63)
         mat(k, 66) = lmat(k, 66)
         mat(k, 67) = lmat(k, 67)
         mat(k, 68) = lmat(k, 68)
         mat(k, 69) = lmat(k, 69)
         mat(k, 70) = lmat(k, 70)
         mat(k, 71) = lmat(k, 71)
         mat(k, 72) = lmat(k, 72)
         mat(k, 73) = lmat(k, 73)
         mat(k, 74) = mat(k, 74) + lmat(k, 74)
         mat(k, 75) = mat(k, 75) + lmat(k, 75)
         mat(k, 77) = mat(k, 77) + lmat(k, 77)
         mat(k, 78) = mat(k, 78) + lmat(k, 78)
         mat(k, 79) = mat(k, 79) + lmat(k, 79)
         mat(k, 80) = mat(k, 80) + lmat(k, 80)
         mat(k, 81) = lmat(k, 81)
         mat(k, 87) = mat(k, 87) + lmat(k, 87)
         mat(k, 93) = lmat(k, 93)
         mat(k, 94) = lmat(k, 94)
         mat(k, 95) = lmat(k, 95)
         mat(k, 96) = lmat(k, 96)
         mat(k, 97) = mat(k, 97) + lmat(k, 97)
         mat(k, 102) = lmat(k, 102)
         mat(k, 103) = lmat(k, 103)
         mat(k, 109) = mat(k, 109) + lmat(k, 109)
         mat(k, 116) = mat(k, 116) + lmat(k, 116)
         mat(k, 126) = mat(k, 126) + lmat(k, 126)
         mat(k, 133) = lmat(k, 133)
         mat(k, 134) = lmat(k, 134)
         mat(k, 135) = lmat(k, 135)
         mat(k, 136) = lmat(k, 136)
         mat(k, 137) = mat(k, 137) + lmat(k, 137)
         mat(k, 142) = mat(k, 142) + lmat(k, 142)
         mat(k, 146) = mat(k, 146) + lmat(k, 146)
         mat(k, 147) = mat(k, 147) + lmat(k, 147)
         mat(k, 148) = mat(k, 148) + lmat(k, 148)
         mat(k, 150) = mat(k, 150) + lmat(k, 150)
         mat(k, 153) = mat(k, 153) + lmat(k, 153)
         mat(k, 154) = lmat(k, 154)
         mat(k, 156) = lmat(k, 156)
         mat(k, 157) = mat(k, 157) + lmat(k, 157)
         mat(k, 158) = mat(k, 158) + lmat(k, 158)
         mat(k, 161) = lmat(k, 161)
         mat(k, 162) = lmat(k, 162)
         mat(k, 163) = lmat(k, 163)
         mat(k, 164) = lmat(k, 164)
         mat(k, 165) = lmat(k, 165)
         mat(k, 166) = lmat(k, 166)
         mat(k, 167) = mat(k, 167) + lmat(k, 167)
         mat(k, 168) = lmat(k, 168)
         mat(k, 170) = mat(k, 170) + lmat(k, 170)
         mat(k, 171) = mat(k, 171) + lmat(k, 171)
         mat(k, 175) = mat(k, 175) + lmat(k, 175)
         mat(k, 177) = mat(k, 177) + lmat(k, 177)
         mat(k, 180) = mat(k, 180) + lmat(k, 180)
         mat(k, 181) = lmat(k, 181)
         mat(k, 182) = lmat(k, 182)
         mat(k, 184) = mat(k, 184) + lmat(k, 184)
         mat(k, 190) = lmat(k, 190)
         mat(k, 191) = lmat(k, 191)
         mat(k, 192) = lmat(k, 192)
         mat(k, 193) = mat(k, 193) + lmat(k, 193)
         mat(k, 199) = mat(k, 199) + lmat(k, 199)
         mat(k, 205) = lmat(k, 205)
         mat(k, 206) = lmat(k, 206)
         mat(k, 207) = lmat(k, 207)
         mat(k, 208) = lmat(k, 208)
         mat(k, 209) = lmat(k, 209)
         mat(k, 210) = lmat(k, 210)
         mat(k, 211) = mat(k, 211) + lmat(k, 211)
         mat(k, 213) = lmat(k, 213)
         mat(k, 214) = lmat(k, 214)
         mat(k, 215) = mat(k, 215) + lmat(k, 215)
         mat(k, 216) = mat(k, 216) + lmat(k, 216)
         mat(k, 218) = lmat(k, 218)
         mat(k, 219) = lmat(k, 219)
         mat(k, 220) = mat(k, 220) + lmat(k, 220)
         mat(k, 221) = lmat(k, 221)
         mat(k, 222) = lmat(k, 222)
         mat(k, 223) = lmat(k, 223)
         mat(k, 224) = lmat(k, 224)
         mat(k, 225) = lmat(k, 225)
         mat(k, 226) = lmat(k, 226)
         mat(k, 227) = mat(k, 227) + lmat(k, 227)
         mat(k, 230) = lmat(k, 230)
         mat(k, 231) = mat(k, 231) + lmat(k, 231)
         mat(k, 232) = mat(k, 232) + lmat(k, 232)
         mat(k, 237) = mat(k, 237) + lmat(k, 237)
         mat(k, 242) = mat(k, 242) + lmat(k, 242)
         mat(k, 244) = mat(k, 244) + lmat(k, 244)
         mat(k, 245) = lmat(k, 245)
         mat(k, 246) = mat(k, 246) + lmat(k, 246)
         mat(k, 247) = mat(k, 247) + lmat(k, 247)
         mat(k, 250) = mat(k, 250) + lmat(k, 250)
         mat(k, 253) = mat(k, 253) + lmat(k, 253)
         mat(k, 254) = lmat(k, 254)
         mat(k, 255) = mat(k, 255) + lmat(k, 255)
         mat(k, 256) = lmat(k, 256)
         mat(k, 258) = mat(k, 258) + lmat(k, 258)
         mat(k, 259) = lmat(k, 259)
         mat(k, 260) = mat(k, 260) + lmat(k, 260)
         mat(k, 261) = lmat(k, 261)
         mat(k, 263) = mat(k, 263) + lmat(k, 263)
         mat(k, 268) = mat(k, 268) + lmat(k, 268)
         mat(k, 276) = mat(k, 276) + lmat(k, 276)
         mat(k, 277) = lmat(k, 277)
         mat(k, 279) = lmat(k, 279)
         mat(k, 280) = lmat(k, 280)
         mat(k, 281) = mat(k, 281) + lmat(k, 281)
         mat(k, 284) = mat(k, 284) + lmat(k, 284)
         mat(k, 290) = lmat(k, 290)
         mat(k, 291) = lmat(k, 291)
         mat(k, 292) = lmat(k, 292)
         mat(k, 293) = lmat(k, 293)
         mat(k, 294) = lmat(k, 294)
         mat(k, 295) = lmat(k, 295)
         mat(k, 296) = mat(k, 296) + lmat(k, 296)
         mat(k, 298) = lmat(k, 298)
         mat(k, 299) = lmat(k, 299)
         mat(k, 300) = lmat(k, 300)
         mat(k, 301) = mat(k, 301) + lmat(k, 301)
         mat(k, 302) = mat(k, 302) + lmat(k, 302)
         mat(k, 304) = lmat(k, 304)
         mat(k, 305) = lmat(k, 305)
         mat(k, 306) = mat(k, 306) + lmat(k, 306)
         mat(k, 308) = mat(k, 308) + lmat(k, 308)
         mat(k, 310) = lmat(k, 310)
         mat(k, 311) = lmat(k, 311)
         mat(k, 312) = mat(k, 312) + lmat(k, 312)
         mat(k, 314) = mat(k, 314) + lmat(k, 314)
         mat(k, 316) = mat(k, 316) + lmat(k, 316)
         mat(k, 317) = lmat(k, 317)
         mat(k, 321) = mat(k, 321) + lmat(k, 321)
         mat(k, 327) = mat(k, 327) + lmat(k, 327)
         mat(k, 329) = mat(k, 329) + lmat(k, 329)
         mat(k, 330) = mat(k, 330) + lmat(k, 330)
         mat(k, 332) = lmat(k, 332)
         mat(k, 333) = mat(k, 333) + lmat(k, 333)
         mat(k, 334) = lmat(k, 334)
         mat(k, 336) = mat(k, 336) + lmat(k, 336)
         mat(k, 338) = lmat(k, 338)
         mat(k, 339) = lmat(k, 339)
         mat(k, 340) = lmat(k, 340)
         mat(k, 341) = lmat(k, 341)
         mat(k, 342) = mat(k, 342) + lmat(k, 342)
         mat(k, 345) = mat(k, 345) + lmat(k, 345)
         mat(k, 346) = lmat(k, 346)
         mat(k, 348) = mat(k, 348) + lmat(k, 348)
         mat(k, 355) = mat(k, 355) + lmat(k, 355)
         mat(k, 356) = lmat(k, 356)
         mat(k, 357) = lmat(k, 357)
         mat(k, 358) = mat(k, 358) + lmat(k, 358)
         mat(k, 361) = lmat(k, 361)
         mat(k, 362) = mat(k, 362) + lmat(k, 362)
         mat(k, 364) = lmat(k, 364)
         mat(k, 365) = lmat(k, 365)
         mat(k, 366) = lmat(k, 366)
         mat(k, 367) = lmat(k, 367)
         mat(k, 368) = mat(k, 368) + lmat(k, 368)
         mat(k, 369) = mat(k, 369) + lmat(k, 369)
         mat(k, 371) = lmat(k, 371)
         mat(k, 372) = lmat(k, 372)
         mat(k, 373) = mat(k, 373) + lmat(k, 373)
         mat(k, 374) = mat(k, 374) + lmat(k, 374)
         mat(k, 376) = mat(k, 376) + lmat(k, 376)
         mat(k, 385) = mat(k, 385) + lmat(k, 385)
         mat(k, 392) = mat(k, 392) + lmat(k, 392)
         mat(k, 396) = mat(k, 396) + lmat(k, 396)
         mat(k, 404) = mat(k, 404) + lmat(k, 404)
         mat(k, 410) = mat(k, 410) + lmat(k, 410)
         mat(k, 412) = lmat(k, 412)
         mat(k, 413) = lmat(k, 413)
         mat(k, 415) = mat(k, 415) + lmat(k, 415)
         mat(k, 416) = mat(k, 416) + lmat(k, 416)
         mat(k, 417) = lmat(k, 417)
         mat(k, 419) = mat(k, 419) + lmat(k, 419)
         mat(k, 420) = mat(k, 420) + lmat(k, 420)
         mat(k, 422) = lmat(k, 422)
         mat(k, 424) = mat(k, 424) + lmat(k, 424)
         mat(k, 426) = mat(k, 426) + lmat(k, 426)
         mat(k, 431) = mat(k, 431) + lmat(k, 431)
         mat(k, 439) = mat(k, 439) + lmat(k, 439)
         mat(k, 447) = mat(k, 447) + lmat(k, 447)
         mat(k, 449) = lmat(k, 449)
         mat(k, 454) = lmat(k, 454)
         mat(k, 455) = mat(k, 455) + lmat(k, 455)
         mat(k, 459) = mat(k, 459) + lmat(k, 459)
         mat(k, 462) = lmat(k, 462)
         mat(k, 463) = lmat(k, 463)
         mat(k, 464) = lmat(k, 464)
         mat(k, 465) = mat(k, 465) + lmat(k, 465)
         mat(k, 467) = mat(k, 467) + lmat(k, 467)
         mat(k, 475) = mat(k, 475) + lmat(k, 475)
         mat(k, 476) = lmat(k, 476)
         mat(k, 477) = lmat(k, 477)
         mat(k, 478) = mat(k, 478) + lmat(k, 478)
         mat(k, 479) = mat(k, 479) + lmat(k, 479)
         mat(k, 481) = lmat(k, 481)
         mat(k, 482) = lmat(k, 482)
         mat(k, 483) = lmat(k, 483)
         mat(k, 484) = lmat(k, 484)
         mat(k, 485) = lmat(k, 485)
         mat(k, 486) = lmat(k, 486)
         mat(k, 487) = mat(k, 487) + lmat(k, 487)
         mat(k, 488) = lmat(k, 488)
         mat(k, 489) = lmat(k, 489)
         mat(k, 490) = lmat(k, 490)
         mat(k, 492) = lmat(k, 492)
         mat(k, 493) = lmat(k, 493)
         mat(k, 494) = lmat(k, 494)
         mat(k, 495) = mat(k, 495) + lmat(k, 495)
         mat(k, 496) = mat(k, 496) + lmat(k, 496)
         mat(k, 500) = mat(k, 500) + lmat(k, 500)
         mat(k, 503) = mat(k, 503) + lmat(k, 503)
         mat(k, 504) = mat(k, 504) + lmat(k, 504)
         mat(k, 505) = lmat(k, 505)
         mat(k, 506) = lmat(k, 506)
         mat(k, 507) = lmat(k, 507)
         mat(k, 512) = mat(k, 512) + lmat(k, 512)
         mat(k, 518) = lmat(k, 518)
         mat(k, 519) = mat(k, 519) + lmat(k, 519)
         mat(k, 522) = lmat(k, 522)
         mat(k, 524) = lmat(k, 524)
         mat(k, 526) = lmat(k, 526)
         mat(k, 527) = lmat(k, 527)
         mat(k, 528) = mat(k, 528) + lmat(k, 528)
         mat(k, 529) = lmat(k, 529)
         mat(k, 530) = mat(k, 530) + lmat(k, 530)
         mat(k, 533) = mat(k, 533) + lmat(k, 533)
         mat(k, 534) = mat(k, 534) + lmat(k, 534)
         mat(k, 536) = lmat(k, 536)
         mat(k, 537) = mat(k, 537) + lmat(k, 537)
         mat(k, 539) = mat(k, 539) + lmat(k, 539)
         mat(k, 541) = mat(k, 541) + lmat(k, 541)
         mat(k, 550) = mat(k, 550) + lmat(k, 550)
         mat(k, 560) = lmat(k, 560)
         mat(k, 561) = lmat(k, 561)
         mat(k, 562) = lmat(k, 562)
         mat(k, 563) = mat(k, 563) + lmat(k, 563)
         mat(k, 567) = lmat(k, 567)
         mat(k, 570) = lmat(k, 570)
         mat(k, 571) = lmat(k, 571)
         mat(k, 572) = mat(k, 572) + lmat(k, 572)
         mat(k, 574) = mat(k, 574) + lmat(k, 574)
         mat(k, 575) = lmat(k, 575)
         mat(k, 577) = lmat(k, 577)
         mat(k, 578) = lmat(k, 578)
         mat(k, 579) = lmat(k, 579)
         mat(k, 581) = mat(k, 581) + lmat(k, 581)
         mat(k, 588) = lmat(k, 588)
         mat(k, 589) = lmat(k, 589)
         mat(k, 590) = lmat(k, 590)
         mat(k, 591) = lmat(k, 591)
         mat(k, 592) = mat(k, 592) + lmat(k, 592)
         mat(k, 596) = lmat(k, 596)
         mat(k, 599) = lmat(k, 599)
         mat(k, 601) = lmat(k, 601)
         mat(k, 602) = mat(k, 602) + lmat(k, 602)
         mat(k, 607) = mat(k, 607) + lmat(k, 607)
         mat(k, 618) = mat(k, 618) + lmat(k, 618)
         mat(k, 630) = mat(k, 630) + lmat(k, 630)
         mat(k, 641) = mat(k, 641) + lmat(k, 641)
         mat(k, 651) = mat(k, 651) + lmat(k, 651)
         mat(k, 659) = mat(k, 659) + lmat(k, 659)
         mat(k, 660) = mat(k, 660) + lmat(k, 660)
         mat(k, 662) = lmat(k, 662)
         mat(k, 669) = mat(k, 669) + lmat(k, 669)
         mat(k, 671) = lmat(k, 671)
         mat(k, 674) = mat(k, 674) + lmat(k, 674)
         mat(k, 677) = mat(k, 677) + lmat(k, 677)
         mat(k, 678) = mat(k, 678) + lmat(k, 678)
         mat(k, 684) = mat(k, 684) + lmat(k, 684)
         mat(k, 686) = mat(k, 686) + lmat(k, 686)
         mat(k, 697) = mat(k, 697) + lmat(k, 697)
         mat(k, 704) = mat(k, 704) + lmat(k, 704)
         mat(k, 708) = mat(k, 708) + lmat(k, 708)
         mat(k, 713) = mat(k, 713) + lmat(k, 713)
         mat(k, 724) = mat(k, 724) + lmat(k, 724)
         mat(k, 735) = lmat(k, 735)
         mat(k, 739) = lmat(k, 739)
         mat(k, 740) = mat(k, 740) + lmat(k, 740)
         mat(k, 750) = mat(k, 750) + lmat(k, 750)
         mat(k, 752) = mat(k, 752) + lmat(k, 752)
         mat(k, 753) = mat(k, 753) + lmat(k, 753)
         mat(k, 757) = lmat(k, 757)
         mat(k, 758) = mat(k, 758) + lmat(k, 758)
         mat(k, 763) = mat(k, 763) + lmat(k, 763)
         mat(k, 765) = lmat(k, 765)
         mat(k, 766) = mat(k, 766) + lmat(k, 766)
         mat(k, 767) = mat(k, 767) + lmat(k, 767)
         mat(k, 773) = mat(k, 773) + lmat(k, 773)
         mat(k, 780) = mat(k, 780) + lmat(k, 780)
         mat(k, 796) = mat(k, 796) + lmat(k, 796)
         mat(k, 798) = lmat(k, 798)
         mat(k, 799) = lmat(k, 799)
         mat(k, 801) = mat(k, 801) + lmat(k, 801)
         mat(k, 802) = lmat(k, 802)
         mat(k, 803) = mat(k, 803) + lmat(k, 803)
         mat(k, 804) = mat(k, 804) + lmat(k, 804)
         mat(k, 805) = mat(k, 805) + lmat(k, 805)
         mat(k, 810) = mat(k, 810) + lmat(k, 810)
         mat(k, 822) = mat(k, 822) + lmat(k, 822)
         mat(k, 832) = mat(k, 832) + lmat(k, 832)
         mat(k, 833) = lmat(k, 833)
         mat(k, 836) = lmat(k, 836)
         mat(k, 848) = mat(k, 848) + lmat(k, 848)
         mat(k, 868) = mat(k, 868) + lmat(k, 868)
         mat(k, 869) = mat(k, 869) + lmat(k, 869)
         mat(k, 871) = mat(k, 871) + lmat(k, 871)
         mat(k, 872) = mat(k, 872) + lmat(k, 872)
         mat(k, 874) = mat(k, 874) + lmat(k, 874)
         mat(k, 876) = lmat(k, 876)
         mat(k, 878) = mat(k, 878) + lmat(k, 878)
         mat(k, 889) = mat(k, 889) + lmat(k, 889)
         mat(k, 908) = mat(k, 908) + lmat(k, 908)
         mat(k, 910) = lmat(k, 910)
         mat(k, 912) = lmat(k, 912)
         mat(k, 914) = lmat(k, 914)
         mat(k, 916) = mat(k, 916) + lmat(k, 916)
         mat(k, 918) = lmat(k, 918)
         mat(k, 919) = lmat(k, 919)
         mat(k, 928) = mat(k, 928) + lmat(k, 928)
         mat(k, 944) = mat(k, 944) + lmat(k, 944)
         mat(k, 964) = mat(k, 964) + lmat(k, 964)
         mat(k, 975) = mat(k, 975) + lmat(k, 975)
         mat(k, 976) = mat(k, 976) + lmat(k, 976)
         mat(k, 977) = mat(k, 977) + lmat(k, 977)
         mat(k, 978) = mat(k, 978) + lmat(k, 978)
         mat(k, 979) = mat(k, 979) + lmat(k, 979)
         mat(k, 981) = mat(k, 981) + lmat(k, 981)
         mat(k, 982) = mat(k, 982) + lmat(k, 982)
         mat(k, 986) = mat(k, 986) + lmat(k, 986)
         mat(k, 995) = mat(k, 995) + lmat(k, 995)
         mat(k, 996) = lmat(k, 996)
         mat(k, 997) = mat(k, 997) + lmat(k, 997)
         mat(k, 998) = mat(k, 998) + lmat(k, 998)
         mat(k,1006) = mat(k,1006) + lmat(k,1006)
         mat(k,1018) = lmat(k,1018)
         mat(k,1019) = lmat(k,1019)
         mat(k,1020) = mat(k,1020) + lmat(k,1020)
         mat(k,1021) = lmat(k,1021)
         mat(k,1022) = lmat(k,1022)
         mat(k,1024) = lmat(k,1024)
         mat(k,1025) = lmat(k,1025)
         mat(k,1027) = mat(k,1027) + lmat(k,1027)
         mat(k,1028) = lmat(k,1028)
         mat(k,1030) = lmat(k,1030)
         mat(k,1032) = mat(k,1032) + lmat(k,1032)
         mat(k,1036) = mat(k,1036) + lmat(k,1036)
         mat(k,1038) = lmat(k,1038)
         mat(k,1039) = mat(k,1039) + lmat(k,1039)
         mat(k,1040) = lmat(k,1040)
         mat(k,1045) = mat(k,1045) + lmat(k,1045)
         mat(k,1055) = lmat(k,1055)
         mat(k,1056) = mat(k,1056) + lmat(k,1056)
         mat(k,1058) = lmat(k,1058)
         mat(k,1060) = lmat(k,1060)
         mat(k,1067) = mat(k,1067) + lmat(k,1067)
         mat(k,1088) = mat(k,1088) + lmat(k,1088)
         mat(k,1103) = mat(k,1103) + lmat(k,1103)
         mat(k,1104) = mat(k,1104) + lmat(k,1104)
         mat(k,1107) = mat(k,1107) + lmat(k,1107)
         mat(k,1108) = mat(k,1108) + lmat(k,1108)
         mat(k,1109) = mat(k,1109) + lmat(k,1109)
         mat(k,1110) = mat(k,1110) + lmat(k,1110)
         mat(k,1115) = mat(k,1115) + lmat(k,1115)
         mat(k,1116) = mat(k,1116) + lmat(k,1116)
         mat(k,1117) = mat(k,1117) + lmat(k,1117)
         mat(k,1119) = lmat(k,1119)
         mat(k,1131) = mat(k,1131) + lmat(k,1131)
         mat(k,1155) = mat(k,1155) + lmat(k,1155)
         mat(k,1171) = lmat(k,1171)
         mat(k,1189) = mat(k,1189) + lmat(k,1189)
         mat(k,1196) = mat(k,1196) + lmat(k,1196)
         mat(k,1214) = mat(k,1214) + lmat(k,1214)
         mat(k,1228) = lmat(k,1228)
         mat(k,1229) = mat(k,1229) + lmat(k,1229)
         mat(k,1234) = mat(k,1234) + lmat(k,1234)
         mat(k,1236) = mat(k,1236) + lmat(k,1236)
         mat(k,1237) = lmat(k,1237)
         mat(k,1254) = mat(k,1254) + lmat(k,1254)
         mat(k,1287) = mat(k,1287) + lmat(k,1287)
         mat(k,1301) = lmat(k,1301)
         mat(k,1303) = mat(k,1303) + lmat(k,1303)
         mat(k,1312) = mat(k,1312) + lmat(k,1312)
         mat(k,1322) = mat(k,1322) + lmat(k,1322)
         mat(k,1324) = mat(k,1324) + lmat(k,1324)
         mat(k,1334) = mat(k,1334) + lmat(k,1334)
         mat(k,1337) = mat(k,1337) + lmat(k,1337)
         mat(k,1382) = mat(k,1382) + lmat(k,1382)
         mat(k,1401) = mat(k,1401) + lmat(k,1401)
         mat(k,1404) = mat(k,1404) + lmat(k,1404)
         mat(k,1406) = lmat(k,1406)
         mat(k,1413) = mat(k,1413) + lmat(k,1413)
         mat(k,1421) = mat(k,1421) + lmat(k,1421)
         mat(k,1424) = mat(k,1424) + lmat(k,1424)
         mat(k,1439) = mat(k,1439) + lmat(k,1439)
         mat(k,1443) = mat(k,1443) + lmat(k,1443)
         mat(k,1444) = lmat(k,1444)
         mat(k,1445) = lmat(k,1445)
         mat(k,1447) = mat(k,1447) + lmat(k,1447)
         mat(k,1450) = mat(k,1450) + lmat(k,1450)
         mat(k,1469) = mat(k,1469) + lmat(k,1469)
         mat(k,1472) = mat(k,1472) + lmat(k,1472)
         mat(k,1479) = mat(k,1479) + lmat(k,1479)
         mat(k,1527) = mat(k,1527) + lmat(k,1527)
         mat(k,1530) = mat(k,1530) + lmat(k,1530)
         mat(k,1536) = mat(k,1536) + lmat(k,1536)
         mat(k,1537) = mat(k,1537) + lmat(k,1537)
         mat(k,1540) = mat(k,1540) + lmat(k,1540)
         mat(k,1541) = mat(k,1541) + lmat(k,1541)
         mat(k,1573) = mat(k,1573) + lmat(k,1573)
         mat(k,1625) = mat(k,1625) + lmat(k,1625)
         mat(k,1640) = mat(k,1640) + lmat(k,1640)
         mat(k,1641) = lmat(k,1641)
         mat(k,1650) = mat(k,1650) + lmat(k,1650)
         mat(k,1660) = mat(k,1660) + lmat(k,1660)
         mat(k,1664) = mat(k,1664) + lmat(k,1664)
         mat(k,1705) = mat(k,1705) + lmat(k,1705)
         mat(k,1711) = mat(k,1711) + lmat(k,1711)
         mat(k,1714) = mat(k,1714) + lmat(k,1714)
         mat(k,1717) = mat(k,1717) + lmat(k,1717)
         mat(k,1733) = lmat(k,1733)
         mat(k,1743) = lmat(k,1743)
         mat(k,1851) = mat(k,1851) + lmat(k,1851)
         mat(k,1852) = mat(k,1852) + lmat(k,1852)
         mat(k,1854) = mat(k,1854) + lmat(k,1854)
         mat(k,1857) = mat(k,1857) + lmat(k,1857)
         mat(k,1860) = mat(k,1860) + lmat(k,1860)
         mat(k,1864) = mat(k,1864) + lmat(k,1864)
         mat(k,1881) = mat(k,1881) + lmat(k,1881)
         mat(k,1888) = mat(k,1888) + lmat(k,1888)
         mat(k,1889) = mat(k,1889) + lmat(k,1889)
         mat(k,1903) = mat(k,1903) + lmat(k,1903)
         mat(k,1909) = lmat(k,1909)
         mat(k,1930) = mat(k,1930) + lmat(k,1930)
         mat(k,1970) = mat(k,1970) + lmat(k,1970)
         mat(k,1972) = mat(k,1972) + lmat(k,1972)
         mat(k,1973) = mat(k,1973) + lmat(k,1973)
         mat(k,1976) = mat(k,1976) + lmat(k,1976)
         mat(k,1977) = mat(k,1977) + lmat(k,1977)
         mat(k,1984) = lmat(k,1984)
         mat(k,1996) = mat(k,1996) + lmat(k,1996)
         mat(k,1998) = lmat(k,1998)
         mat(k,2000) = mat(k,2000) + lmat(k,2000)
         mat(k,2001) = mat(k,2001) + lmat(k,2001)
         mat(k,2004) = lmat(k,2004)
         mat(k,2005) = mat(k,2005) + lmat(k,2005)
         mat(k,2008) = mat(k,2008) + lmat(k,2008)
         mat(k,2009) = mat(k,2009) + lmat(k,2009)
         mat(k,2011) = mat(k,2011) + lmat(k,2011)
         mat(k,2013) = mat(k,2013) + lmat(k,2013)
         mat(k,2014) = lmat(k,2014)
         mat(k,2015) = mat(k,2015) + lmat(k,2015)
         mat(k,2016) = mat(k,2016) + lmat(k,2016)
         mat(k,2019) = lmat(k,2019)
         mat(k,2020) = lmat(k,2020)
         mat(k,2022) = mat(k,2022) + lmat(k,2022)
         mat(k,2024) = mat(k,2024) + lmat(k,2024)
         mat(k,2027) = mat(k,2027) + lmat(k,2027)
         mat(k,2029) = lmat(k,2029)
         mat(k,2030) = mat(k,2030) + lmat(k,2030)
         mat(k,2044) = mat(k,2044) + lmat(k,2044)
         mat(k,2047) = lmat(k,2047)
         mat(k,2050) = mat(k,2050) + lmat(k,2050)
         mat(k,2082) = mat(k,2082) + lmat(k,2082)
         mat(k,2083) = lmat(k,2083)
         mat(k,2087) = mat(k,2087) + lmat(k,2087)
         mat(k,2127) = mat(k,2127) + lmat(k,2127)
         mat(k,2132) = mat(k,2132) + lmat(k,2132)
         mat(k,2153) = mat(k,2153) + lmat(k,2153)
         mat(k, 229) = 0._r8
         mat(k, 285) = 0._r8
         mat(k, 287) = 0._r8
         mat(k, 324) = 0._r8
         mat(k, 352) = 0._r8
         mat(k, 386) = 0._r8
         mat(k, 389) = 0._r8
         mat(k, 514) = 0._r8
         mat(k, 515) = 0._r8
         mat(k, 520) = 0._r8
         mat(k, 521) = 0._r8
         mat(k, 525) = 0._r8
         mat(k, 531) = 0._r8
         mat(k, 532) = 0._r8
         mat(k, 535) = 0._r8
         mat(k, 555) = 0._r8
         mat(k, 558) = 0._r8
         mat(k, 559) = 0._r8
         mat(k, 565) = 0._r8
         mat(k, 566) = 0._r8
         mat(k, 568) = 0._r8
         mat(k, 569) = 0._r8
         mat(k, 593) = 0._r8
         mat(k, 595) = 0._r8
         mat(k, 597) = 0._r8
         mat(k, 598) = 0._r8
         mat(k, 600) = 0._r8
         mat(k, 608) = 0._r8
         mat(k, 609) = 0._r8
         mat(k, 611) = 0._r8
         mat(k, 612) = 0._r8
         mat(k, 615) = 0._r8
         mat(k, 629) = 0._r8
         mat(k, 631) = 0._r8
         mat(k, 633) = 0._r8
         mat(k, 634) = 0._r8
         mat(k, 636) = 0._r8
         mat(k, 638) = 0._r8
         mat(k, 652) = 0._r8
         mat(k, 653) = 0._r8
         mat(k, 656) = 0._r8
         mat(k, 672) = 0._r8
         mat(k, 689) = 0._r8
         mat(k, 691) = 0._r8
         mat(k, 692) = 0._r8
         mat(k, 699) = 0._r8
         mat(k, 700) = 0._r8
         mat(k, 703) = 0._r8
         mat(k, 714) = 0._r8
         mat(k, 719) = 0._r8
         mat(k, 721) = 0._r8
         mat(k, 731) = 0._r8
         mat(k, 733) = 0._r8
         mat(k, 736) = 0._r8
         mat(k, 741) = 0._r8
         mat(k, 746) = 0._r8
         mat(k, 761) = 0._r8
         mat(k, 762) = 0._r8
         mat(k, 795) = 0._r8
         mat(k, 816) = 0._r8
         mat(k, 828) = 0._r8
         mat(k, 837) = 0._r8
         mat(k, 846) = 0._r8
         mat(k, 851) = 0._r8
         mat(k, 856) = 0._r8
         mat(k, 858) = 0._r8
         mat(k, 864) = 0._r8
         mat(k, 865) = 0._r8
         mat(k, 866) = 0._r8
         mat(k, 877) = 0._r8
         mat(k, 887) = 0._r8
         mat(k, 892) = 0._r8
         mat(k, 897) = 0._r8
         mat(k, 899) = 0._r8
         mat(k, 905) = 0._r8
         mat(k, 906) = 0._r8
         mat(k, 907) = 0._r8
         mat(k, 925) = 0._r8
         mat(k, 926) = 0._r8
         mat(k, 929) = 0._r8
         mat(k, 931) = 0._r8
         mat(k, 932) = 0._r8
         mat(k, 933) = 0._r8
         mat(k, 936) = 0._r8
         mat(k, 938) = 0._r8
         mat(k, 946) = 0._r8
         mat(k, 948) = 0._r8
         mat(k, 949) = 0._r8
         mat(k, 951) = 0._r8
         mat(k, 957) = 0._r8
         mat(k, 958) = 0._r8
         mat(k, 959) = 0._r8
         mat(k, 965) = 0._r8
         mat(k, 966) = 0._r8
         mat(k, 967) = 0._r8
         mat(k, 972) = 0._r8
         mat(k, 980) = 0._r8
         mat(k, 984) = 0._r8
         mat(k,1008) = 0._r8
         mat(k,1009) = 0._r8
         mat(k,1013) = 0._r8
         mat(k,1015) = 0._r8
         mat(k,1023) = 0._r8
         mat(k,1026) = 0._r8
         mat(k,1031) = 0._r8
         mat(k,1053) = 0._r8
         mat(k,1064) = 0._r8
         mat(k,1069) = 0._r8
         mat(k,1071) = 0._r8
         mat(k,1076) = 0._r8
         mat(k,1080) = 0._r8
         mat(k,1083) = 0._r8
         mat(k,1084) = 0._r8
         mat(k,1085) = 0._r8
         mat(k,1086) = 0._r8
         mat(k,1087) = 0._r8
         mat(k,1089) = 0._r8
         mat(k,1090) = 0._r8
         mat(k,1091) = 0._r8
         mat(k,1092) = 0._r8
         mat(k,1098) = 0._r8
         mat(k,1101) = 0._r8
         mat(k,1114) = 0._r8
         mat(k,1123) = 0._r8
         mat(k,1137) = 0._r8
         mat(k,1139) = 0._r8
         mat(k,1140) = 0._r8
         mat(k,1142) = 0._r8
         mat(k,1146) = 0._r8
         mat(k,1147) = 0._r8
         mat(k,1148) = 0._r8
         mat(k,1149) = 0._r8
         mat(k,1150) = 0._r8
         mat(k,1151) = 0._r8
         mat(k,1153) = 0._r8
         mat(k,1154) = 0._r8
         mat(k,1156) = 0._r8
         mat(k,1158) = 0._r8
         mat(k,1164) = 0._r8
         mat(k,1165) = 0._r8
         mat(k,1167) = 0._r8
         mat(k,1168) = 0._r8
         mat(k,1170) = 0._r8
         mat(k,1174) = 0._r8
         mat(k,1177) = 0._r8
         mat(k,1179) = 0._r8
         mat(k,1181) = 0._r8
         mat(k,1182) = 0._r8
         mat(k,1184) = 0._r8
         mat(k,1185) = 0._r8
         mat(k,1186) = 0._r8
         mat(k,1188) = 0._r8
         mat(k,1190) = 0._r8
         mat(k,1191) = 0._r8
         mat(k,1192) = 0._r8
         mat(k,1198) = 0._r8
         mat(k,1199) = 0._r8
         mat(k,1201) = 0._r8
         mat(k,1202) = 0._r8
         mat(k,1204) = 0._r8
         mat(k,1211) = 0._r8
         mat(k,1212) = 0._r8
         mat(k,1215) = 0._r8
         mat(k,1221) = 0._r8
         mat(k,1223) = 0._r8
         mat(k,1224) = 0._r8
         mat(k,1226) = 0._r8
         mat(k,1230) = 0._r8
         mat(k,1235) = 0._r8
         mat(k,1238) = 0._r8
         mat(k,1239) = 0._r8
         mat(k,1244) = 0._r8
         mat(k,1245) = 0._r8
         mat(k,1246) = 0._r8
         mat(k,1247) = 0._r8
         mat(k,1248) = 0._r8
         mat(k,1252) = 0._r8
         mat(k,1253) = 0._r8
         mat(k,1263) = 0._r8
         mat(k,1264) = 0._r8
         mat(k,1266) = 0._r8
         mat(k,1272) = 0._r8
         mat(k,1274) = 0._r8
         mat(k,1289) = 0._r8
         mat(k,1290) = 0._r8
         mat(k,1296) = 0._r8
         mat(k,1297) = 0._r8
         mat(k,1299) = 0._r8
         mat(k,1308) = 0._r8
         mat(k,1314) = 0._r8
         mat(k,1323) = 0._r8
         mat(k,1325) = 0._r8
         mat(k,1326) = 0._r8
         mat(k,1330) = 0._r8
         mat(k,1343) = 0._r8
         mat(k,1344) = 0._r8
         mat(k,1345) = 0._r8
         mat(k,1346) = 0._r8
         mat(k,1347) = 0._r8
         mat(k,1357) = 0._r8
         mat(k,1365) = 0._r8
         mat(k,1385) = 0._r8
         mat(k,1388) = 0._r8
         mat(k,1389) = 0._r8
         mat(k,1392) = 0._r8
         mat(k,1393) = 0._r8
         mat(k,1395) = 0._r8
         mat(k,1410) = 0._r8
         mat(k,1412) = 0._r8
         mat(k,1414) = 0._r8
         mat(k,1415) = 0._r8
         mat(k,1417) = 0._r8
         mat(k,1419) = 0._r8
         mat(k,1420) = 0._r8
         mat(k,1422) = 0._r8
         mat(k,1423) = 0._r8
         mat(k,1432) = 0._r8
         mat(k,1433) = 0._r8
         mat(k,1435) = 0._r8
         mat(k,1438) = 0._r8
         mat(k,1440) = 0._r8
         mat(k,1442) = 0._r8
         mat(k,1448) = 0._r8
         mat(k,1455) = 0._r8
         mat(k,1456) = 0._r8
         mat(k,1457) = 0._r8
         mat(k,1458) = 0._r8
         mat(k,1459) = 0._r8
         mat(k,1460) = 0._r8
         mat(k,1468) = 0._r8
         mat(k,1473) = 0._r8
         mat(k,1475) = 0._r8
         mat(k,1476) = 0._r8
         mat(k,1481) = 0._r8
         mat(k,1482) = 0._r8
         mat(k,1483) = 0._r8
         mat(k,1485) = 0._r8
         mat(k,1492) = 0._r8
         mat(k,1502) = 0._r8
         mat(k,1505) = 0._r8
         mat(k,1506) = 0._r8
         mat(k,1509) = 0._r8
         mat(k,1512) = 0._r8
         mat(k,1524) = 0._r8
         mat(k,1525) = 0._r8
         mat(k,1526) = 0._r8
         mat(k,1528) = 0._r8
         mat(k,1529) = 0._r8
         mat(k,1533) = 0._r8
         mat(k,1535) = 0._r8
         mat(k,1538) = 0._r8
         mat(k,1539) = 0._r8
         mat(k,1542) = 0._r8
         mat(k,1561) = 0._r8
         mat(k,1562) = 0._r8
         mat(k,1567) = 0._r8
         mat(k,1591) = 0._r8
         mat(k,1595) = 0._r8
         mat(k,1596) = 0._r8
         mat(k,1599) = 0._r8
         mat(k,1601) = 0._r8
         mat(k,1604) = 0._r8
         mat(k,1609) = 0._r8
         mat(k,1615) = 0._r8
         mat(k,1633) = 0._r8
         mat(k,1634) = 0._r8
         mat(k,1638) = 0._r8
         mat(k,1644) = 0._r8
         mat(k,1645) = 0._r8
         mat(k,1647) = 0._r8
         mat(k,1651) = 0._r8
         mat(k,1653) = 0._r8
         mat(k,1655) = 0._r8
         mat(k,1657) = 0._r8
         mat(k,1659) = 0._r8
         mat(k,1673) = 0._r8
         mat(k,1676) = 0._r8
         mat(k,1682) = 0._r8
         mat(k,1684) = 0._r8
         mat(k,1686) = 0._r8
         mat(k,1687) = 0._r8
         mat(k,1689) = 0._r8
         mat(k,1692) = 0._r8
         mat(k,1695) = 0._r8
         mat(k,1696) = 0._r8
         mat(k,1697) = 0._r8
         mat(k,1698) = 0._r8
         mat(k,1700) = 0._r8
         mat(k,1716) = 0._r8
         mat(k,1718) = 0._r8
         mat(k,1772) = 0._r8
         mat(k,1790) = 0._r8
         mat(k,1802) = 0._r8
         mat(k,1804) = 0._r8
         mat(k,1806) = 0._r8
         mat(k,1830) = 0._r8
         mat(k,1839) = 0._r8
         mat(k,1865) = 0._r8
         mat(k,1883) = 0._r8
         mat(k,1886) = 0._r8
         mat(k,1891) = 0._r8
         mat(k,1892) = 0._r8
         mat(k,1893) = 0._r8
         mat(k,1895) = 0._r8
         mat(k,1905) = 0._r8
         mat(k,1907) = 0._r8
         mat(k,1913) = 0._r8
         mat(k,1920) = 0._r8
         mat(k,1932) = 0._r8
         mat(k,1933) = 0._r8
         mat(k,1934) = 0._r8
         mat(k,1945) = 0._r8
         mat(k,1948) = 0._r8
         mat(k,1950) = 0._r8
         mat(k,1954) = 0._r8
         mat(k,1955) = 0._r8
         mat(k,1956) = 0._r8
         mat(k,1960) = 0._r8
         mat(k,1961) = 0._r8
         mat(k,1962) = 0._r8
         mat(k,1964) = 0._r8
         mat(k,1968) = 0._r8
         mat(k,1974) = 0._r8
         mat(k,1975) = 0._r8
         mat(k,1978) = 0._r8
         mat(k,1983) = 0._r8
         mat(k,1985) = 0._r8
         mat(k,1986) = 0._r8
         mat(k,1987) = 0._r8
         mat(k,1988) = 0._r8
         mat(k,1989) = 0._r8
         mat(k,1990) = 0._r8
         mat(k,1991) = 0._r8
         mat(k,1992) = 0._r8
         mat(k,1993) = 0._r8
         mat(k,1994) = 0._r8
         mat(k,1995) = 0._r8
         mat(k,1997) = 0._r8
         mat(k,1999) = 0._r8
         mat(k,2002) = 0._r8
         mat(k,2003) = 0._r8
         mat(k,2018) = 0._r8
         mat(k,2025) = 0._r8
         mat(k,2028) = 0._r8
         mat(k,2033) = 0._r8
         mat(k,2034) = 0._r8
         mat(k,2035) = 0._r8
         mat(k,2036) = 0._r8
         mat(k,2037) = 0._r8
         mat(k,2038) = 0._r8
         mat(k,2039) = 0._r8
         mat(k,2041) = 0._r8
         mat(k,2042) = 0._r8
         mat(k,2043) = 0._r8
         mat(k,2045) = 0._r8
         mat(k,2046) = 0._r8
         mat(k,2049) = 0._r8
         mat(k,2051) = 0._r8
         mat(k,2052) = 0._r8
         mat(k,2089) = 0._r8
         mat(k,2129) = 0._r8
         mat(k,2130) = 0._r8
         mat(k,2131) = 0._r8
         mat(k,2133) = 0._r8
         mat(k,2135) = 0._r8
         mat(k,2136) = 0._r8
         mat(k,2137) = 0._r8
         mat(k,2139) = 0._r8
         mat(k,2140) = 0._r8
         mat(k,2141) = 0._r8
         mat(k,2143) = 0._r8
         mat(k,2146) = 0._r8
         mat(k,2148) = 0._r8
         mat(k,2150) = 0._r8
         mat(k,2151) = 0._r8
         mat(k,2152) = 0._r8
         mat(k, 1) = mat(k, 1) - dti
         mat(k, 2) = mat(k, 2) - dti
         mat(k, 3) = mat(k, 3) - dti
         mat(k, 4) = mat(k, 4) - dti
         mat(k, 5) = mat(k, 5) - dti
         mat(k, 6) = mat(k, 6) - dti
         mat(k, 7) = mat(k, 7) - dti
         mat(k, 13) = mat(k, 13) - dti
         mat(k, 19) = mat(k, 19) - dti
         mat(k, 20) = mat(k, 20) - dti
         mat(k, 21) = mat(k, 21) - dti
         mat(k, 22) = mat(k, 22) - dti
         mat(k, 23) = mat(k, 23) - dti
         mat(k, 24) = mat(k, 24) - dti
         mat(k, 25) = mat(k, 25) - dti
         mat(k, 26) = mat(k, 26) - dti
         mat(k, 27) = mat(k, 27) - dti
         mat(k, 28) = mat(k, 28) - dti
         mat(k, 29) = mat(k, 29) - dti
         mat(k, 30) = mat(k, 30) - dti
         mat(k, 31) = mat(k, 31) - dti
         mat(k, 32) = mat(k, 32) - dti
         mat(k, 33) = mat(k, 33) - dti
         mat(k, 34) = mat(k, 34) - dti
         mat(k, 35) = mat(k, 35) - dti
         mat(k, 36) = mat(k, 36) - dti
         mat(k, 37) = mat(k, 37) - dti
         mat(k, 38) = mat(k, 38) - dti
         mat(k, 39) = mat(k, 39) - dti
         mat(k, 40) = mat(k, 40) - dti
         mat(k, 41) = mat(k, 41) - dti
         mat(k, 42) = mat(k, 42) - dti
         mat(k, 43) = mat(k, 43) - dti
         mat(k, 44) = mat(k, 44) - dti
         mat(k, 45) = mat(k, 45) - dti
         mat(k, 46) = mat(k, 46) - dti
         mat(k, 47) = mat(k, 47) - dti
         mat(k, 50) = mat(k, 50) - dti
         mat(k, 53) = mat(k, 53) - dti
         mat(k, 56) = mat(k, 56) - dti
         mat(k, 60) = mat(k, 60) - dti
         mat(k, 63) = mat(k, 63) - dti
         mat(k, 66) = mat(k, 66) - dti
         mat(k, 69) = mat(k, 69) - dti
         mat(k, 71) = mat(k, 71) - dti
         mat(k, 74) = mat(k, 74) - dti
         mat(k, 77) = mat(k, 77) - dti
         mat(k, 80) = mat(k, 80) - dti
         mat(k, 87) = mat(k, 87) - dti
         mat(k, 93) = mat(k, 93) - dti
         mat(k, 97) = mat(k, 97) - dti
         mat(k, 102) = mat(k, 102) - dti
         mat(k, 109) = mat(k, 109) - dti
         mat(k, 116) = mat(k, 116) - dti
         mat(k, 126) = mat(k, 126) - dti
         mat(k, 133) = mat(k, 133) - dti
         mat(k, 137) = mat(k, 137) - dti
         mat(k, 142) = mat(k, 142) - dti
         mat(k, 146) = mat(k, 146) - dti
         mat(k, 150) = mat(k, 150) - dti
         mat(k, 153) = mat(k, 153) - dti
         mat(k, 158) = mat(k, 158) - dti
         mat(k, 161) = mat(k, 161) - dti
         mat(k, 164) = mat(k, 164) - dti
         mat(k, 167) = mat(k, 167) - dti
         mat(k, 171) = mat(k, 171) - dti
         mat(k, 175) = mat(k, 175) - dti
         mat(k, 180) = mat(k, 180) - dti
         mat(k, 184) = mat(k, 184) - dti
         mat(k, 190) = mat(k, 190) - dti
         mat(k, 193) = mat(k, 193) - dti
         mat(k, 199) = mat(k, 199) - dti
         mat(k, 205) = mat(k, 205) - dti
         mat(k, 211) = mat(k, 211) - dti
         mat(k, 216) = mat(k, 216) - dti
         mat(k, 221) = mat(k, 221) - dti
         mat(k, 227) = mat(k, 227) - dti
         mat(k, 232) = mat(k, 232) - dti
         mat(k, 237) = mat(k, 237) - dti
         mat(k, 242) = mat(k, 242) - dti
         mat(k, 247) = mat(k, 247) - dti
         mat(k, 250) = mat(k, 250) - dti
         mat(k, 255) = mat(k, 255) - dti
         mat(k, 260) = mat(k, 260) - dti
         mat(k, 268) = mat(k, 268) - dti
         mat(k, 276) = mat(k, 276) - dti
         mat(k, 284) = mat(k, 284) - dti
         mat(k, 290) = mat(k, 290) - dti
         mat(k, 296) = mat(k, 296) - dti
         mat(k, 302) = mat(k, 302) - dti
         mat(k, 308) = mat(k, 308) - dti
         mat(k, 314) = mat(k, 314) - dti
         mat(k, 321) = mat(k, 321) - dti
         mat(k, 327) = mat(k, 327) - dti
         mat(k, 333) = mat(k, 333) - dti
         mat(k, 339) = mat(k, 339) - dti
         mat(k, 342) = mat(k, 342) - dti
         mat(k, 348) = mat(k, 348) - dti
         mat(k, 355) = mat(k, 355) - dti
         mat(k, 362) = mat(k, 362) - dti
         mat(k, 369) = mat(k, 369) - dti
         mat(k, 376) = mat(k, 376) - dti
         mat(k, 385) = mat(k, 385) - dti
         mat(k, 392) = mat(k, 392) - dti
         mat(k, 396) = mat(k, 396) - dti
         mat(k, 404) = mat(k, 404) - dti
         mat(k, 410) = mat(k, 410) - dti
         mat(k, 415) = mat(k, 415) - dti
         mat(k, 420) = mat(k, 420) - dti
         mat(k, 426) = mat(k, 426) - dti
         mat(k, 431) = mat(k, 431) - dti
         mat(k, 439) = mat(k, 439) - dti
         mat(k, 447) = mat(k, 447) - dti
         mat(k, 455) = mat(k, 455) - dti
         mat(k, 459) = mat(k, 459) - dti
         mat(k, 467) = mat(k, 467) - dti
         mat(k, 475) = mat(k, 475) - dti
         mat(k, 483) = mat(k, 483) - dti
         mat(k, 487) = mat(k, 487) - dti
         mat(k, 496) = mat(k, 496) - dti
         mat(k, 503) = mat(k, 503) - dti
         mat(k, 512) = mat(k, 512) - dti
         mat(k, 519) = mat(k, 519) - dti
         mat(k, 530) = mat(k, 530) - dti
         mat(k, 541) = mat(k, 541) - dti
         mat(k, 550) = mat(k, 550) - dti
         mat(k, 563) = mat(k, 563) - dti
         mat(k, 574) = mat(k, 574) - dti
         mat(k, 581) = mat(k, 581) - dti
         mat(k, 592) = mat(k, 592) - dti
         mat(k, 607) = mat(k, 607) - dti
         mat(k, 618) = mat(k, 618) - dti
         mat(k, 630) = mat(k, 630) - dti
         mat(k, 641) = mat(k, 641) - dti
         mat(k, 651) = mat(k, 651) - dti
         mat(k, 660) = mat(k, 660) - dti
         mat(k, 669) = mat(k, 669) - dti
         mat(k, 677) = mat(k, 677) - dti
         mat(k, 686) = mat(k, 686) - dti
         mat(k, 697) = mat(k, 697) - dti
         mat(k, 704) = mat(k, 704) - dti
         mat(k, 708) = mat(k, 708) - dti
         mat(k, 713) = mat(k, 713) - dti
         mat(k, 724) = mat(k, 724) - dti
         mat(k, 740) = mat(k, 740) - dti
         mat(k, 750) = mat(k, 750) - dti
         mat(k, 758) = mat(k, 758) - dti
         mat(k, 766) = mat(k, 766) - dti
         mat(k, 780) = mat(k, 780) - dti
         mat(k, 796) = mat(k, 796) - dti
         mat(k, 803) = mat(k, 803) - dti
         mat(k, 810) = mat(k, 810) - dti
         mat(k, 822) = mat(k, 822) - dti
         mat(k, 832) = mat(k, 832) - dti
         mat(k, 848) = mat(k, 848) - dti
         mat(k, 869) = mat(k, 869) - dti
         mat(k, 889) = mat(k, 889) - dti
         mat(k, 908) = mat(k, 908) - dti
         mat(k, 916) = mat(k, 916) - dti
         mat(k, 928) = mat(k, 928) - dti
         mat(k, 944) = mat(k, 944) - dti
         mat(k, 964) = mat(k, 964) - dti
         mat(k, 976) = mat(k, 976) - dti
         mat(k, 986) = mat(k, 986) - dti
         mat(k, 995) = mat(k, 995) - dti
         mat(k,1006) = mat(k,1006) - dti
         mat(k,1020) = mat(k,1020) - dti
         mat(k,1032) = mat(k,1032) - dti
         mat(k,1036) = mat(k,1036) - dti
         mat(k,1045) = mat(k,1045) - dti
         mat(k,1056) = mat(k,1056) - dti
         mat(k,1067) = mat(k,1067) - dti
         mat(k,1088) = mat(k,1088) - dti
         mat(k,1104) = mat(k,1104) - dti
         mat(k,1116) = mat(k,1116) - dti
         mat(k,1131) = mat(k,1131) - dti
         mat(k,1155) = mat(k,1155) - dti
         mat(k,1189) = mat(k,1189) - dti
         mat(k,1214) = mat(k,1214) - dti
         mat(k,1234) = mat(k,1234) - dti
         mat(k,1254) = mat(k,1254) - dti
         mat(k,1287) = mat(k,1287) - dti
         mat(k,1303) = mat(k,1303) - dti
         mat(k,1322) = mat(k,1322) - dti
         mat(k,1337) = mat(k,1337) - dti
         mat(k,1382) = mat(k,1382) - dti
         mat(k,1413) = mat(k,1413) - dti
         mat(k,1447) = mat(k,1447) - dti
         mat(k,1472) = mat(k,1472) - dti
         mat(k,1530) = mat(k,1530) - dti
         mat(k,1625) = mat(k,1625) - dti
         mat(k,1650) = mat(k,1650) - dti
         mat(k,1711) = mat(k,1711) - dti
         mat(k,1860) = mat(k,1860) - dti
         mat(k,1888) = mat(k,1888) - dti
         mat(k,1930) = mat(k,1930) - dti
         mat(k,1973) = mat(k,1973) - dti
         mat(k,2000) = mat(k,2000) - dti
         mat(k,2027) = mat(k,2027) - dti
         mat(k,2050) = mat(k,2050) - dti
         mat(k,2132) = mat(k,2132) - dti
         mat(k,2153) = mat(k,2153) - dti
      end do
      end subroutine nlnmat_finit
      subroutine nlnmat( ofl, ofu, chnkpnts, mat, y, rxt, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      integer, intent(in) :: ofl
      integer, intent(in) :: ofu
      integer, intent(in) :: chnkpnts
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(chnkpnts,nzcnt)
      real(r8), intent(in) :: y(chnkpnts,gas_pcnst)
      real(r8), intent(in) :: rxt(chnkpnts,rxntot)
      real(r8), intent(inout) :: mat(chnkpnts,nzcnt)
      call nlnmat01( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat02( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat03( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat04( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat05( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat06( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat07( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat08( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat09( ofl, ofu, chnkpnts, mat, y, rxt )
      call nlnmat_finit( ofl, ofu, chnkpnts, mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
