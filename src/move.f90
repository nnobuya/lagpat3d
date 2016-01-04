subroutine move( dt, dt0, ist_pt, v_pt, v_pt_p, x_pt )

  use mod_cnst, only: npt, ndim
  use mod_set , only: nx1, nx2, nx3, int_t, x_fld

  implicit none

  integer, intent(in):: ist_pt(1:npt)
  real(8), intent(in)   :: dt, dt0
  real(8), intent(in)   :: v_pt(ndim,npt)
  real(8), intent(in)   :: v_pt_p(1:ndim,0:4,1:npt)
  real(8), intent(inout):: x_pt(ndim,npt)

  real(8):: c1, c2, c3, c4
  real(8):: v_comb(1:ndim)
  integer:: i


  do i = 1, npt

     if ( ist_pt(i) /= 0 ) cycle

     !..time evol.x
     if      ( int_t == 0 ) then
        v_comb(1:ndim) = v_pt(1:ndim,i)
     else if ( int_t == 1 ) then
        v_comb(1:ndim) = 0.5d0 *( v_pt(1:ndim,i) + v_pt_p(1:ndim,0,i) )
     else if ( int_t == 2 ) then

        c1 = dt /dt0 - 0.d5
        c2 = 0.d5
        c3 = 1.d0 - dt /dt0
        c4 = 0.d0

        v_comb(1:ndim) &
             & = c1 *v_pt_p(1:ndim,1,i) + c2 *v_pt_p(1:ndim,2,i) &
               + c3 *v_pt_p(1:ndim,3,i) + c4 *v_pt_p(1:ndim,4,i)
     else
        stop 'error: no method'
     end if

     x_pt(1,i) = x_pt(1,i) + v_comb(1) *dt
     if ( nx2 >= 2 ) x_pt(2,i) = x_pt(2,i) + v_comb(2) *dt
     if ( nx3 >= 2 ) x_pt(3,i) = x_pt(3,i) + v_comb(3) *dt
     !..avoid negative z
     x_pt(3,i) = max(0.d0, x_pt(3,i))

  end do


  return

end subroutine move
