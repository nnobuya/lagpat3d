subroutine pt_set_3d(ndim, nx1, nx2, nx3, npt, npt_x1, npt_x2, npt_x3, &
     & x_fld, d_fld, x_pt, dma)

  implicit none

  !..param.
  real(8), parameter:: rm_sol = 1.989d33
  real(8), parameter:: de_high = 1.e10
  logical, parameter:: debug = .false.

  !..io
  integer, intent(in) :: ndim, nx1, nx2, nx3, npt, npt_x1, npt_x2, npt_x3
  real(8), intent(in) :: x_fld(1:ndim,1:nx1,1:nx2,1:nx3), &
       & d_fld(1:nx1,1:nx2,1:nx3)
  real(8), intent(out):: x_pt(1:ndim,1:npt), dma(1:npt)

  !..local
  integer:: nd_x1, nd_x2, nd_x3
  real(8):: grid_mass, total
  real(8):: x1(1:nx1+1), x2(1:nx2+1), x3(1:nx3+1)
  real(8), dimension(1:nx1,1:nx2,1:nx3):: dv, dma_fld, fac
  integer:: i1_in(1:npt_x1), i2_in(1:npt_x2), i3_in(1:npt_x3), &
       & i1_ou(1:npt_x1), i2_ou(1:npt_x2), i3_ou(1:npt_x3)

  integer:: i11, i12, i21, i22, i31, i32
  integer:: i, j, k, l, ipt


  write(*,*)
  write(*,'(" ---------------", "pt_set_3d(): setting tracer particles", &
       & " ---------------------")')
  write(*,*)


  ! ------------------------------------------------------------------ !
  !     make grid                                                      !
  ! ------------------------------------------------------------------ !

  !..x1-grid
  x1(2:nx1+1) = x_fld(1,1:nx1,1,1)
  x2(2:nx2+1) = x_fld(2,1,1:nx2,1)
  x3(1:nx3)   = x_fld(3,1,1,1:nx3)

  x1(1) = x1(2) - (x1(3)- x1(2))
  x2(1) = x2(2) - (x2(3)- x2(2))
  x3(nx3+1) = x3(nx3) + (x3(nx3)- x3(nx3-1))

  do k = 1, nx3
     do j = 1, nx2
        do i = 1, nx1
           dv(i,j,k) = (x1(i+1) - x1(i)) *(x2(j+1) - x2(j)) *(x3(k+1) - x3(k))
        end do
     end do
  end do

  dma_fld(1:nx1,1:nx2,1:nx3) = &
       & d_fld(1:nx1,1:nx2,1:nx3) *dv(1:nx1,1:nx2,1:nx3)/rm_sol
  fac(1:nx1,1:nx2,1:nx3) = &
       & dma_fld(1:nx1,1:nx2,1:nx3) /sum(dma_fld(1:nx1,1:nx2,1:nx3))

  grid_mass = 0.d0
  do k = 1, nx3
     do j = 1, nx2
        do i = 1, nx1
           if (d_fld(i,j,k) < de_high) then
              grid_mass = grid_mass + dma_fld(i,j,k)
           else
              write(*,'(5x, "avoid NS region: density >", 1pe10.2, " g/cc")') &
                   & de_high
           end if
        end do
     end do
  end do

  !     make grid                                                      !
  ! ------------------------------------------------------------------ !



  ! ------------------------------------------------------------------ !
  !     grid data check                                                !
  ! ------------------------------------------------------------------ !
  if (debug) then
     do k = 1, 1
        do j = 1, nx2
           do i = 1, nx1
              write(61,'(1p, *(e14.5))') &
              & x_fld(1,i,j,k), x_fld(2,i,j,k), log10(d_fld(i,j,k))
           end do
           write(61,*)
        end do
     end do
     close(61)

     do k = 1, nx3
        do j = nx2/2, nx2/2
           do i = 1, nx1
              write(62,'(1p, *(e14.5))') &
              & x_fld(1,i,j,k), x_fld(3,i,j,k), log10(d_fld(i,j,k))
           end do
           write(62,*)
        end do
     end do
     close(62)

     do k = 1, nx3
        do j = 1, nx2
           do i = nx1/2, nx1/2
              write(63,'(1p, *(e14.5))') &
              & x_fld(2,i,j,k), x_fld(3,i,j,k), log10(d_fld(i,j,k))
           end do
        end do
        write(63,*)
     end do
     close(63)
  end if
  !     grid data check                                                !
  ! ------------------------------------------------------------------ !



  ! ------------------------------------------------------------------ !
  !     set particle position and mass                                 !
  ! ------------------------------------------------------------------ !

  nd_x1 = nx1 /npt_x1
  nd_x2 = nx2 /npt_x2
  nd_x3 = nx3 /npt_x3

  write(70,'("# pt-grid x1")')
  do i = 1, npt_x1
     if (i == 1) then
        i1_in(i) = 1
     else
        i1_in(i) = i1_ou(i-1) + 1
     end if

     if (i <= mod(nx1,npt_x1))then
        i1_ou(i) = i1_in(i) + nd_x1
     else
        i1_ou(i) = i1_in(i) + nd_x1 - 1
     end if

     write(70,'(4i10)') i, i1_in(i), i1_ou(i), i1_ou(i) - i1_in(i) +1
  end do


  write(70,'("# pt-grid x2")')
  do i = 1, npt_x2
     if (i == 1) then
        i2_in(i) = 1
     else
        i2_in(i) = i2_ou(i-1) + 1
     end if

     if (i <= mod(nx2,npt_x2))then
        i2_ou(i) = i2_in(i) + nd_x2
     else
        i2_ou(i) = i2_in(i) + nd_x2 - 1
     end if

     write(70,'(4i10)') i, i2_in(i), i2_ou(i), i2_ou(i) - i2_in(i) +1
  end do


  write(70,'("# pt-grid x3")')
  do i = 1, npt_x3
     if (i == 1) then
        i3_in(i) = 1
     else
        i3_in(i) = i3_ou(i-1) + 1
     end if

     if (i <= mod(nx3,npt_x3))then
        i3_ou(i) = i3_in(i) + nd_x3
     else
        i3_ou(i) = i3_in(i) + nd_x3 - 1
     end if

     write(70,'(4i10)') i, i3_in(i), i3_ou(i), i3_ou(i) - i3_in(i) +1
  end do



  do k = 1, npt_x3
     do j = 1, npt_x2
        do i = 1, npt_x1

           i11 = i1_in(i)
           i12 = i1_ou(i)
           i21 = i2_in(j)
           i22 = i2_ou(j)
           i31 = i3_in(k)
           i32 = i3_ou(k)

           total = sum(fac(i11:i12, i21:i22, i31:i32))


           do l = 1, ndim
              x_pt(l,ipt) = sum(fac(i11:i12, i21:i22, i31:i32)&
                   & *x_fld(l, i11:i12, i21:i22, i31:i32)) /total
           end do

           if (maxval(d_fld(i11:i12, i21:i22, i31:i32)) > de_high) then
              dma(ipt) = 0.d0
           else
              dma(ipt) = sum(fac(i11:i12, i21:i22, i31:i32) &
                   & *dma_fld(i11:i12, i21:i22, i31:i32)) /total
           end if

        end do
     end do
  end do

  !     set particle position                                          !
  ! ------------------------------------------------------------------ !


  write(*,'(" ------------------------------------------------------", &
       & "--------------------")')

  return

end subroutine pt_set_3d
