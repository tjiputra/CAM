module classnuc

  implicit none
  private
  save
 
  public classnuc_in, preexisting_ice

! Note: classnuc_in can only be set to true for MODAL_AERO, 
! Classnuc_in =.false. in a CAM5 compset gives back the original CAM5
!#ifdef MODAL_AERO
# if (defined MODAL_AERO) || (defined OSLO_AERO)
  logical, parameter :: classnuc_in = .true.
#else
  logical, parameter :: classnuc_in = .false.
#endif     

  logical, parameter :: preexisting_ice = .false.
    
    
end module classnuc
