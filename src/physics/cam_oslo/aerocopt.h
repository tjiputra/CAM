!    For subroutines initaeropt and intaeropt1to3,4,6to10:

!SOA       common /aerocopt1/ bep1to3, bep4, bep5to10
!soa       common /aerocopt1/ bep1wsoa, bep1to3, bep4, bep5to10
!soa       common /aerocopt1/ bep1wsoa, bep2to3, bep4, bep5to10
       common /aerocopt1/ bep1, bep2to3, bep4, bep5to10

!soa       real(r8) bep1wsoa(38,10,16,6)
       real(r8) bep1(38,10,6,16,6)
!SOA
!soa       real(r8) bep1to3(38,10,16,3)
       real(r8) bep2to3(38,10,16,6,2:3)
!soa       real(r8) bep4(38,10,16,6,6)
       real(r8) bep4(38,10,6,16,6,6)
!soa
       real(r8) bep5to10(38,10,6,6,6,6,5:10)
