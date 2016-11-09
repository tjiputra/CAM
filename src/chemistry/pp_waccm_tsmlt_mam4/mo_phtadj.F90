      module mo_phtadj
      private
      public :: phtadj
      contains
      subroutine phtadj( p_rate, inv, m, ncol, nlev )
      use chem_mods, only : nfs, phtcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol, nlev
      real(r8), intent(in) :: inv(ncol,nlev,max(1,nfs))
      real(r8), intent(in) :: m(ncol,nlev)
      real(r8), intent(inout) :: p_rate(ncol,nlev,max(1,phtcnt))
!--------------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------------
      integer :: k
      real(r8) :: im(ncol,nlev)
      do k = 1,nlev
         im(:ncol,k) = 1._r8 / m(:ncol,k)
         p_rate(:,k,132) = p_rate(:,k,132) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,136) = p_rate(:,k,136) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,137) = p_rate(:,k,137) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,139) = p_rate(:,k,139) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,144) = p_rate(:,k,144) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,148) = p_rate(:,k,148) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,149) = p_rate(:,k,149) * inv(:,k, 2) * im(:,k)
         p_rate(:,k,151) = p_rate(:,k,151) * inv(:,k, 2) * im(:,k)
      end do
      end subroutine phtadj
      end module mo_phtadj
