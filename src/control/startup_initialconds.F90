module startup_initialconds
!----------------------------------------------------------------------- 
! 
! Wrapper for calls to initialize buffers and read initial/topo files
! 
!-----------------------------------------------------------------------

implicit none
private
save

public :: initial_conds ! Read in initial conditions (dycore dependent)

!======================================================================= 
contains
!======================================================================= 

subroutine initial_conds(dyn_in)

   ! This routine calls the dycore dependent
   ! routine read_inidat.

#if (defined BFB_CAM_SCAM_IOP )
   use history_defaults, only: initialize_iop_history
#endif
   
   use inidat,        only: read_inidat
   use dyn_comp,      only: dyn_import_t
   use cam_pio_utils, only: clean_iodesc_list

   type(dyn_import_t),  intent(inout) :: dyn_in
   !-----------------------------------------------------------------------

#if (defined BFB_CAM_SCAM_IOP )
   call initialize_iop_history
#endif

   ! Read in initial data
   call read_inidat(dyn_in)

   call clean_iodesc_list()

end subroutine initial_conds

!======================================================================= 

end module startup_initialconds
