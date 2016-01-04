subroutine status( x_pt, ist_pt )

  use mod_cnst, only: npt, ndim
  use mod_set , only: x_fld, nx1, nx2, nx3

  implicit none

  !..io
  real(8), intent(in)   :: x_pt(ndim,npt)
  integer         , intent(inout):: ist_pt(npt)

  !..tmp
  real(8):: bnd_in(1:ndim), bnd_out(1:ndim)
  integer:: i, j


  do j = 1, ndim
     bnd_in(j)  = x_fld(j,1,1,1)
     bnd_out(j) = x_fld(j,nx1,nx2,nx3)
  end do

  !! status now
  do i = 1, npt
     if( ist_pt(i) /= 0 ) cycle

     do j = 1, ndim
        if ( x_pt(j,i) < bnd_in(j) .or. x_pt(j,i) > bnd_out(j) ) then
           ist_pt(i) =  1
        end if
     end do

  end do


  return

end subroutine status
