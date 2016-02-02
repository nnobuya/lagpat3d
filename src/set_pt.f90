subroutine set_pt(istg, ti, ist_pt, id, dma, x_pt, v_pt, d_fld)

  use mod_cnst, only: npt, npt_x1, npt_x2, npt_x3, ndim
  use mod_set , only: k_zoku, nx1, nx2, nx3, x_fld, dx_fld

  implicit none

  !..io
  integer, intent(out):: istg, id(1:ndim,1:npt), ist_pt(npt)
  real(8), intent(out):: ti, dma(1:npt), x_pt(1:ndim,1:npt), v_pt(1:ndim,1:npt)

  !..local
  real(8):: d_fld(1:nx1,1:nx2,1:nx3)
  integer:: i, j, k, i_tmp, j_tmp, ipt, nskip1, nskip2, nskip3


  if      ( k_zoku == 0 ) then

     !..tmp
     !write(100) npt_x1, npt_x2, npt_x3
     !write(100) nx1, nx2, nx3
     !write(100) x_fld(1:ndim,1:nx1,1:nx2,1:nx3)
     !write(100) d_fld(1:nx1,1:nx2,1:nx3)
     !write(100) dx_fld(1,1:nx1), dx_fld(2,1:nx2), dx_fld(3,1:nx3)
     !--------------------------------


     x_pt   = 0.d0
     nskip1 = nx1 /npt_x1
     nskip2 = nx2 /npt_x2
     nskip3 = nx3 /npt_x3

     ipt = 0
     dma = 0.d0
     do k = 1, nx3, nskip3
        do j = 1, nx2, nskip2
           do i = 1, nx1, nskip1
            ipt = ipt + 1
            if( ipt < npt )then
            x_pt(1:ndim,ipt) = x_fld(1:ndim,i,j,k)
            dma(ipt) = d_fld(i,j,k) &
                 & *dx_fld(1,i) *dble(nskip1) &
                 & *dx_fld(2,j) *dble(nskip2) &
                 & *dx_fld(3,k) *dble(nskip3)
!            write(10000,'(99es12.3)')x_pt(1,ipt),x_pt(2,ipt),x_pt(3,ipt)
            end if
           end do
        end do
     end do

     istg               = 0
     ist_pt(1:npt)      = 0
     v_pt(1:ndim,1:npt) = 0.d0

     if(ipt < npt)ist_pt(ipt+1:npt) = - 1


     write(60,*)
     write(60,'(4i10)') npt, nx1, nx2, nx3
     write(60,*)
     write(60,*)
     write(60,*)

     do ipt = 1, npt
        write(60,'(4i10, 1p, 10e14.5)') &
             & ipt, 1, 1, 1, dma(ipt), x_pt(1:ndim,ipt)
     end do

     close(60)

     write(*,*) nx1, nx2, nx3

     !stop 'db: set_pt'

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


  id(1:ndim,1:npt) = 0

  return

end subroutine set_pt
