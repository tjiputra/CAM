module fv_control_mod

use shr_kind_mod, only: r8=>shr_kind_r8

implicit none
public
save

real(r8) :: tmass0
real(r8) :: zgsint

integer :: nsplit                  ! Lagrangian time splits
integer :: nspltrac                ! Tracer time splits
integer :: nspltvrm                ! Vertical re-mapping time splits

! _ord = 1: first order upwind
! _ord = 2: 2nd order van Leer (Lin et al 1994)
! _ord = 3: standard PPM 
! _ord = 4: enhanced PPM (default)
integer :: iord                     ! scheme to be used in E-W direction
integer :: jord                     ! scheme to be used in N-S direction
integer :: kord                     ! scheme to be used for vertical mapping

logical :: dyn_conservative = .false.  ! Flag indicating whether the dynamics is conservative
integer :: filtcw                   ! flag for filtering c-grid winds
integer :: ct_overlap               ! nonzero for overlap of cd_core and trac2d, 0 otherwise
integer :: trac_decomp              ! size of tracer domain decomposition for trac2d
integer :: fft_flt                  ! 0 => FFT/algebraic filter; 1 => FFT filter
integer :: div24del2flag            ! 2 for 2nd order div damping, 4 for 4th order div damping,
                                    ! 42 for 4th order div damping plus 2nd order velocity damping
real(r8):: del2coef                 ! strength of 2nd order velocity damping

end module fv_control_mod
