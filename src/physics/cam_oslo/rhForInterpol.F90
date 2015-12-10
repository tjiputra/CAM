
   subroutine rhForInterpol (lchnk, ncol, rhum, xrh, irh1, irh2)


   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8
   use opttab, only: rh

   implicit none

!
! Input arguments
!
      integer, intent(in) :: lchnk              ! chunk identifier
      integer, intent(in) :: ncol               ! number of atmospheric columns
      real(r8), intent(in) :: rhum(pcols,pver)  ! level relative humidity (fraction)
!
! Output arguments
!
      real(r8), intent(out) :: xrh(pcols,pver)
      integer, intent(out) :: irh1(pcols,pver)
      integer, intent(out) :: irh2(pcols,pver)
!
!---------------------------Local variables-----------------------------
!
      integer k, icol, irelh 
!
!------------------------------------------------------------------------
!

!      write(*,*) 'Before xrh-loop'
      do k=1,pver
        do icol=1,ncol
          xrh(icol,k)  = min(max(rhum(icol,k),rh(1)),rh(10))
        end do 
      end do

!      write(*,*) 'Before rh-loop'
      do irelh=1,9
       do k=1,pver
        do icol=1,ncol
           if(xrh(icol,k) >= rh(irelh).and. &
             xrh(icol,k)<=rh(irelh+1)) then
             irh1(icol,k)=irelh
             irh2(icol,k)=irelh+1
           endif
         end do
       end do
      end do

!      write(*,*) 'xrh, irh1, irh2 =', xrh(1,26), irh1(1,26), irh2(1,26)
 
      return

end subroutine rhForInterpol
