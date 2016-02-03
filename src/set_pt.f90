subroutine set_pt(istg, ti, ist_pt, id, dma, x_pt, v_pt, d_fld)

  use mod_cnst, only: npt, npt_x1, npt_x2, npt_x3, ndim, rm_sol
  use mod_set , only: k_zoku, nx1, nx2, nx3, x_fld, dx_fld

  implicit none

  !..io
  real(8), intent(in) :: d_fld(1:nx1,1:nx2,1:nx3)
  integer, intent(out):: istg, id(1:ndim,1:npt), ist_pt(npt)
  real(8), intent(out):: ti, dma(1:npt), x_pt(1:ndim,1:npt), v_pt(1:ndim,1:npt)

  !..local
  integer:: i, j, k, i_tmp, j_tmp, ipt


  if      ( k_zoku == 0 ) then

     !..set x_pt and dma
      call pt_set_3d(ndim, nx1, nx2, nx3, npt, npt_x1, npt_x2, npt_x3, &
          & x_fld(1:ndim,1:nx1,1:nx2,1:nx3), d_fld(1:nx1,1:nx2,1:nx3), &
          & x_pt(1:ndim,1:npt), dma(1:npt))
    !  in: others
     ! out: x_pt, dma


     !..other values
     istg               = 0
     ist_pt(1:npt)      = 0
     v_pt(1:ndim,1:npt) = 0.d0

     do ipt = 1, npt
        if (dma(ipt) <= 0.d0) ist_pt(ipt) = -1
     end do


     write(70,*) '#      npt,      nx1,      nx2,      nx3,'
     write(70,'(4i10)') npt, nx1, nx2, nx3
     write(70,*) '#'

     ipt = 0
     do k = 1, npt_x3
        do j = 1, npt_x2
           do i = 1, npt_x1
              ipt = ipt + 1

              id(1,ipt) = i
              id(2,ipt) = j
              id(3,ipt) = k
              write(70,'(4i10, 1p, 10e14.5)') &
                   & ipt, i, j, k, dma(ipt), x_pt(1:ndim,ipt)
           end do
        end do
     end do
     close(70)

  else if ( k_zoku == 1 ) then

     read(91,*)
     read(91,*) ti, istg
     do i = 1, npt
        read(91,*) &
             & i_tmp, j_tmp, &
             & dma(i), x_pt(1:2,i), v_pt(i,1:2), ist_pt(i)
     end do

  else
     !! error
     stop '### Error: "k_zoku" isn''t 0 or 1.  ###'

  end if


  return

end subroutine set_pt
